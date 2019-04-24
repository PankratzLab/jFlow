package org.genvisis.seq.manage.mosdepth;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;

import org.genvisis.cnv.Launch;
import org.genvisis.cnv.filesys.Compression;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.filesys.SampleList;
import org.genvisis.cnv.manage.MarkerDataLoader;
import org.genvisis.cnv.manage.TempFileTranspose;
import org.genvisis.cnv.manage.TransposeData;
import org.pankratzlab.common.CmdLine;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.GenomeBuild;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.qsub.Qsub;

/**
 * Framework / common-code for marker-dominant parsing of data files into Genvisis files.
 */
public abstract class AbstractParsingPipeline {

  protected String projDir;
  protected String propFileDir;
  protected String projName;
  // PBS job ID for identifying job-specific files
  protected String jobID;
  protected Project proj;
  // Initalize new logger - after project creation, this will be reassigned to the project's logger
  protected Logger log = new Logger();
  // default X/Y scale factor
  protected double scaleFactor;
  // ideally this will be scaled per-project based on number of markers and samples
  protected int numMarkersPerFile;
  protected long fingerprintForMarkerFiles = 0L;

  public AbstractParsingPipeline(int scale, int mkrsPerFile) {
    this.scaleFactor = scale;
    this.numMarkersPerFile = mkrsPerFile;
  }

  public void setProjectName(String projName2) {
    projName = projName2;
  }

  public void setProjectPropertiesDir(String propFileDir2) {
    propFileDir = propFileDir2;
  }

  public void setProjectDir(String projDir2) {
    projDir = projDir2;
  }

  public void setMarkersPerMDRAF(int mkrsPerFile) {
    numMarkersPerFile = mkrsPerFile;
  }

  public void setJobID(String jobId) {
    this.jobID = jobId;
  }

  public void setScaleFactor(double scale) {
    scaleFactor = scale;
  }

  protected abstract byte getNullStatus();

  /**
   * Create a Genvisis project, or load it if one with the specified name already exists.
   */
  public void createProject() {
    String propFile = ext.verifyDirFormat(propFileDir)
                      + ext.replaceWithLinuxSafeCharacters(projName) + ".properties";
    if (!Files.exists(propFile)) {
      Files.write((new Project()).PROJECT_NAME.getName() + "=" + projName, propFile);
      proj = new Project(propFile);
      log = proj.getLog();
      proj.PROJECT_NAME.setValue(projName);
      proj.PROJECT_DIRECTORY.setValue(projDir);
      proj.XY_SCALE_FACTOR.setValue(scaleFactor);
      proj.TARGET_MARKERS_FILENAMES.setValue(new String[] {});
      proj.ID_HEADER.setValue("NULL");

      setAdditionalProjectProperties();

      proj.saveProperties();
      log.reportTime("Created project properties file: " + propFile);
    } else {
      log.reportTime("Project properties file already exists at " + propFile
                     + "; skipping creation.");
      proj = new Project(propFile);
      log = proj.getLog();
    }
  }

  /**
   * Used in {@link AbstractParsingPipeline#createProject()} to set additional properties such as
   * {@link ARRAY} or {@link GenomeBuild}.
   */
  protected abstract void setAdditionalProjectProperties();

  public void writeLookup() {
    if (!Files.exists(proj.MARKERLOOKUP_FILENAME.getValue())) {
      doWriteLookup();
      log.reportTime("Created marker lookup file: " + proj.MARKERLOOKUP_FILENAME.getValue());
    } else {
      log.reportTime("Project marker lookup file already exists; skipping creation.");
    }
  }

  /**
   * Called by {@link AbstractParsingPipeline#writeLookup()} if the marker lookup file doesn't
   * already exist.
   */
  protected abstract void doWriteLookup();

  public void createSampleList() {
    String[] samples = parseSamples();
    fingerprintForMarkerFiles = MarkerSet.fingerprint(samples);
    if (!Files.exists(proj.SAMPLELIST_FILENAME.getValue())) {
      SampleList sl = new SampleList(samples);
      sl.serialize(proj.SAMPLELIST_FILENAME.getValue());
      sl = null;
      samples = null;
      log.reportTime("Created sample list file: " + proj.SAMPLELIST_FILENAME.getValue());
    } else {
      log.reportTime("Project sample list file already exists; skipping creation.");
    }
  }

  /**
   * Loads sample names
   * 
   * @return a write-safe array of sample names
   */
  protected abstract String[] parseSamples();

  protected abstract int getNumSamples();

  protected abstract int getNumMarkers();

  protected RandomAccessFile openMDRAF(String filename, int nInd, byte nullStatus,
                                       String[] mkrNames) throws IOException {
    byte[] mkrBytes = Compression.objToBytes(mkrNames);
    byte[] mdRAFHeader = TransposeData.getParameterSectionForMdRaf(nInd, mkrNames.length,
                                                                   nullStatus,
                                                                   fingerprintForMarkerFiles,
                                                                   mkrBytes);
    mkrBytes = null;

    String file = proj.MARKER_DATA_DIRECTORY.getValue(true, true) + filename;
    if (Files.exists(file)) {
      new File(file).delete();
    }

    RandomAccessFile mdRAF = new RandomAccessFile(file, "rw");
    mdRAF.write(mdRAFHeader);
    mdRAFHeader = null;

    return mdRAF;
  }

  protected void buildOutliers() {
    try {
      MarkerDataLoader.buildOutliersFromMDRAFs(proj);
    } catch (ClassNotFoundException | IOException e) {
      log.reportError("Problem occurred while loading outliers from marker files. Attempting to continue...");
      log.reportException(e);
    }
  }

  protected void createSampRAFsFromMDRAFs() {
    TempFileTranspose tft = new TempFileTranspose(proj, proj.PROJECT_DIRECTORY.getValue() + "temp/",
                                                  jobID);
    tft.setupMarkerListFile();
    tft.setupSampleListFile();

    int gb = 240;
    int wall = 240;
    int proc = 12;

    String file = setupTransposeScripts(proj.PROJECT_DIRECTORY.getValue(), getNumSamples(),
                                        getNumMarkers(), gb, wall, proc);
    CmdLine.run("qsub " + file, proj.PROJECT_DIRECTORY.getValue());
  }

  private String setupTransposeScripts(String runDir, int numSamples, int numMarkers, int qGBLim,
                                       int qWallLim, int qProcLim) {
    String jar = Launch.getJarLocation();
    String jobName1 = "tempTransposeFirst.qsub";
    String jobName2 = "tempTransposeSecond.qsub";

    long bytesPerMkrF = Sample.getNBytesPerSampleMarker(getNullStatus()) * numSamples
                        * numMarkersPerFile + TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN;
    long bytesPerSmpF = Sample.getNBytesPerSampleMarker(getNullStatus()) * numMarkers
                        + Sample.PARAMETER_SECTION_BYTES;

    long bLim = ((long) qGBLim) * 1024 * 1024 * 1024;

    int numF = 1;
    while (((numF + 1) * bytesPerMkrF) < (bLim * .8)) {
      numF++;
    }
    numF = Math.min(numF, qProcLim);

    String jobCmd1 = "cd " + runDir + "\njava -jar " + jar + " " + TempFileTranspose.class.getName()
                     + " proj=" + proj.getPropertyFilename() + " jobID=$PBS_JOBID type=M qsub="
                     + runDir + jobName1;

    Qsub.qsub(runDir + jobName1, jobCmd1, qGBLim * 1024, qWallLim, numF);

    numF = 1;
    while (((numF + 1) * bytesPerSmpF) < (bLim * .8)) {
      numF++;
    }
    numF = Math.min(numF, qProcLim);

    String jobCmd2 = "cd " + runDir + "\njava -jar " + jar + " " + TempFileTranspose.class.getName()
                     + " proj=" + proj.getPropertyFilename() + " jobID=$PBS_JOBID type=S qsub="
                     + runDir + jobName2;
    Qsub.qsub(runDir + jobName2, jobCmd2, qGBLim * 1024, qWallLim, numF);

    return jobName1;
  }

  protected abstract void createMarkerPositions();

}
