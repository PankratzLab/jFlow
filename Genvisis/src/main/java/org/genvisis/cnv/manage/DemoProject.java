package org.genvisis.cnv.manage;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.SampleList;
import org.genvisis.cnv.manage.ExtProjectDataParser.ProjectDataParserBuilder;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.ext;

public class DemoProject extends Project {
  public static enum DEMO_TYPE {
                                MARKER_FOCUS, SAMPLE_FOCUS
  }

  private static void copyGeneTrack(Project original, Project demo) {
    String geneTrack = original.getGeneTrackFilename(true);
    if (geneTrack != null && Files.exists(geneTrack)) {
      demo.setProperty(demo.GENETRACK_FILENAME, ext.removeDirectoryInfo(geneTrack));
      original.getLog()
              .reportTimeInfo("Detected " + geneTrack + ", copying to "
                              + demo.GENETRACK_FILENAME.getValue(false, false)
                              + "\n\t (this takes a while due to byte by byte copying)");
      if (!Files.exists(demo.GENETRACK_FILENAME.getValue(false, false))) {
        Files.copyFileExactlyByteByByte(geneTrack, demo.GENETRACK_FILENAME.getValue(false, false));
      }
    }
  }

  private static void copyStratResults(Project original, Project demo) {
    String[] strats = Array.toStringArray(original.getStratResults());
    if (strats != null && strats.length > 0) {
      for (String strat : strats) {
        Files.copyFile(strat, demo.PROJECT_DIRECTORY.getValue() + ext.removeDirectoryInfo(strat));
      }
    } else {
      original.getLog().reportTimeWarning("Did not find any stratification results");
    }

  }

  private static ExtProjectDataParser getParserFor(Project proj, String file, boolean sampleBased) {
    ExtProjectDataParser.ProjectDataParserBuilder builder = new ProjectDataParserBuilder();
    builder.dataKeyColumnIndex(0);
    builder.stringDataIndices(new int[][] {{0}});
    builder.treatAllNumeric(false);
    builder.hasHeader(false);
    builder.skipUnFound(false);
    builder.sampleBased(sampleBased);
    ExtProjectDataParser sampleParser = null;
    try {
      sampleParser = builder.build(proj, file);
    } catch (FileNotFoundException e) {
      proj.getLog().reportException(e);
      e.printStackTrace();
    }
    sampleParser.loadData();
    return sampleParser;
  }

  private static String[] loadMarkers(Project proj, String markersFile, DEMO_TYPE dType) {
    String[] markersToUse = null;
    String targetMarkerFile = proj.TARGET_MARKERS_FILENAMES.getDefaultValueString();
    if (markersFile != null) {
      markersToUse = getParserFor(proj, markersFile, false).getStringDataAt(0, true);
    } else if (Files.exists(targetMarkerFile) && dType == DEMO_TYPE.MARKER_FOCUS) {
      markersToUse = getParserFor(proj, targetMarkerFile, false).getStringDataAt(0, true);
    } else {
      if (dType == DEMO_TYPE.MARKER_FOCUS) {
        proj.getLog().reportTimeWarning("A marker subset file was not provided and " + markersFile
                                        + " did not exist, exporting all markers");
      } else {
        proj.getLog().reportTimeInfo(
                                     "A marker subset file was not provided, exporting all markers for demotype "
                                     + DEMO_TYPE.SAMPLE_FOCUS);

      }
    }
    return markersToUse;
  }

  private static boolean[] loadSamples(Project proj, String samplesFile, DEMO_TYPE dType) {
    boolean[] samplesToUse = new boolean[proj.getSamples().length];
    Arrays.fill(samplesToUse, true);// defualt to all samples;
    // String sampleSubsetFile = proj.getFilename(proj.SAMPLE_SUBSET_FILENAME);
    String sampleSubsetFile = proj.SAMPLE_SUBSET_FILENAME.getValue();
    if (samplesFile != null) {
      samplesToUse = getParserFor(proj, samplesFile, true).getDataPresent();
    } else if (Files.exists(sampleSubsetFile) && dType == DEMO_TYPE.SAMPLE_FOCUS) {
      samplesToUse = getParserFor(proj, sampleSubsetFile, true).getDataPresent();
    } else {
      if (dType == DEMO_TYPE.SAMPLE_FOCUS) {
        proj.getLog()
            .reportTimeWarning("A sample subset file was not provided and " + sampleSubsetFile
                               + " did not exist, exporting all samples");
      } else {
        proj.getLog().reportTimeInfo(
                                     "A sample subset file was not provided, exporting all samples for demotype "
                                     + DEMO_TYPE.MARKER_FOCUS);
      }
    }
    return samplesToUse;
  }

  public static void main(String[] args) {
    test();
  }

  public static void test() {
    Project proj = new Project(null, false);
    DemoProject demoProject = new DemoProject(proj, proj.PROJECT_DIRECTORY.getValue() + "Demo/",
                                              true, DEMO_TYPE.MARKER_FOCUS);
    String[] adfa = null;
    demoProject.createProjectDemo(adfa, null, 8);
  }

  private final Project proj;

  /**
   * where the new project directory will be placed
   */
  private final String demoDirectory;

  private final boolean overwriteExisting;

  private boolean fail;

  private final DEMO_TYPE dType;

  public DemoProject(Project proj, String demoDirectory, boolean overwriteExisting,
                     DEMO_TYPE dType) {
    super();
    this.proj = proj;// (Project) proj.clone();
    this.demoDirectory = demoDirectory;
    this.overwriteExisting = overwriteExisting;
    this.dType = dType;
    fail = false;
    init();
  }

  // private void copyFileIfExists(String property) {
  private void copyFileIfExists(FileProperty property, FileProperty other) {
    String file = property.getValue(false, false);
    if (file != null && Files.exists(file)) {
      System.out.println("Copying " + file + " to " + other.getValue(false, false));
      Files.copyFile(file, other.getValue(false, false));
    } else {
      proj.getLog().reportTimeWarning("Did not find file " + file + " cannot copy to demo");
    }
  }

  private void copyFiles(String[] filesOriginal, String[] filesNew) {
    if (filesOriginal != null && filesOriginal.length > 0) {
      for (int i = 0; i < filesOriginal.length; i++) {
        if (filesOriginal[i] != null && Files.exists(filesOriginal[i])) {
          Files.copyFile(filesOriginal[i], filesNew[i]);
        } else {
          proj.getLog()
              .reportTimeWarning("Did not find file " + filesOriginal[i] + " cannot copy to demo");
        }
      }
    }
  }

  private void copyIfFilesExists(String propertyWithMultipleFiles) {
    // String propertyValue = proj.getProperty(propertyWithMultipleFiles);
    String propertyValue = proj.getProperty(propertyWithMultipleFiles).getValueString();
    String propertyValueNew = propertyValue;
    propertyValueNew =
        propertyValueNew.replace(proj.PROJECT_DIRECTORY.getValue(), PROJECT_DIRECTORY.getValue());
    setProperty(propertyWithMultipleFiles, propertyValueNew);

    // String[] filesOriginal = proj.getFilenames(propertyWithMultipleFiles);
    // String[] filesOriginal =
    // proj.<MultiFileProperty>getProperty(propertyWithMultipleFiles).getValue();
    // String[] filesNew =
    // this.<MultiFileProperty>getProperty(propertyWithMultipleFiles).getValue(true);
    String[] filesOriginal =
        proj.<StringListProperty>getProperty(propertyWithMultipleFiles).getValue();
    String[] filesNew = this.<StringListProperty>getProperty(propertyWithMultipleFiles).getValue();
    copyFiles(filesOriginal, filesNew);
  }

  public boolean createProjectDemo(String markersFile, String samplesFile, int numThreads) {
    if (!fail) {
      System.out.println(proj.MARKERSET_FILENAME.getValue());
      return createProjectDemo(loadMarkers(proj, markersFile, dType),
                               loadSamples(proj, samplesFile, dType), numThreads);
    } else {
      return false;
    }
  }

  public boolean createProjectDemo(String[] markersToExport, boolean[] samplesToUse,
                                   int numThreads) {

    boolean created = true;
    if (!fail) {
      if (markersToExport == null) {
        markersToExport = proj.getMarkerNames();
        proj.getLog().reportTimeWarning("Since the markers to export were not provided, all "
                                        + markersToExport.length + " markers will be exported");
      }
      if (samplesToUse == null) {
        samplesToUse = new boolean[proj.getSamples().length];
        Arrays.fill(samplesToUse, true);
        proj.getLog().reportTimeWarning("Since the samples to export were not provided, all "
                                        + samplesToUse.length + " samples will be exported");
      }
      if (samplesToUse.length != proj.getSamples().length) {
        proj.getLog().reportTimeError("The array length provided does (" + samplesToUse.length
                                      + ") does not contain boolean values for all samples");
        created = false;
        return created;
      } else {
        // Markers.orderMarkers(markersToExport, getFilename(proj.MARKER_POSITION_FILENAME),
        // getFilename(proj.MARKERSET_FILENAME, true, true), proj.getLog());
        Markers.orderMarkers(markersToExport, proj.MARKER_POSITION_FILENAME.getValue(),
                             MARKERSET_FILENAME.getValue(true, true), proj.getLog());
        long fingerPrint = getMarkerSet().getFingerprint();
        proj.getLog().reportTimeInfo("Attempting to subset the samples...");
        new File(SAMPLELIST_FILENAME.getValue(false, false)).delete();
        created = FocusedSample.focusAProject(proj, this, markersToExport, samplesToUse,
                                              fingerPrint, numThreads, overwriteExisting, getLog());
        if (created) {
          if (dType == DEMO_TYPE.MARKER_FOCUS) {
            proj.getLog()
                .reportTimeInfo("Finished subsetting the samples...Attempting to transpose the data");
            TransposeData.transposeData(this, 2000000000, false);
            // Files.writeList(markersToExport, getFilename(proj.DISPLAY_MARKERS_FILENAME));
            Files.writeList(markersToExport, DISPLAY_MARKERS_FILENAMES.getValue()[0]);
          }
          SampleList.generateSampleList(this)
                    .writeToTextFile(PROJECT_DIRECTORY.getValue() + "ListOfSamples.txt");
          getSamples();
        } else {
          proj.getLog().reportTimeError("Could not create focused samples, halting...");
        }
      }
    } else {
      created = false;
    }
    return created;
  }

  public DEMO_TYPE getdType() {
    return dType;
  }

  private void init() {
    String demoProjectDirectory =
        demoDirectory + proj.PROJECT_NAME.getValue().replaceAll(" ", "") + "_" + dType + "/";
    setProperty(PROJECT_DIRECTORY, demoProjectDirectory);
    if (!Files.exists(PROJECT_DIRECTORY.getValue()) || overwriteExisting) {
      new File(PROJECT_DIRECTORY.getValue()).mkdirs();
      DATA_DIRECTORY.getValue(true, false);
      SAMPLE_DIRECTORY.getValue(true, false);
      RESULTS_DIRECTORY.getValue(true, false);

      if (dType == DEMO_TYPE.MARKER_FOCUS) {
        MARKER_DATA_DIRECTORY.getValue(true, false);
      }

      // single file Copies
      copyFileIfExists(proj.MARKER_POSITION_FILENAME, MARKER_POSITION_FILENAME);
      copyFileIfExists(proj.SAMPLE_DATA_FILENAME, SAMPLE_DATA_FILENAME);
      copyFileIfExists(proj.SAMPLE_QC_FILENAME, SAMPLE_QC_FILENAME);
      copyFileIfExists(proj.MOSAIC_RESULTS_FILENAME, MOSAIC_RESULTS_FILENAME);
      copyFileIfExists(proj.SEXCHECK_RESULTS_FILENAME, SEXCHECK_RESULTS_FILENAME);
      copyFileIfExists(proj.INTENSITY_PC_FILENAME, INTENSITY_PC_FILENAME);
      copyFileIfExists(proj.CLUSTER_FILTER_COLLECTION_FILENAME, CLUSTER_FILTER_COLLECTION_FILENAME);
      copyFileIfExists(proj.ANNOTATION_FILENAME, ANNOTATION_FILENAME);
      copyFileIfExists(proj.BLAST_ANNOTATION_FILENAME, BLAST_ANNOTATION_FILENAME);
      Files.copyFile(proj.BLAST_ANNOTATION_FILENAME.getValue() + ".tbi",
                     BLAST_ANNOTATION_FILENAME.getValue() + ".tbi");
      // copyFileIfExists(proj.MARKER_POSITION_FILENAME.getName());
      // copyFileIfExists(proj.SAMPLE_DATA_FILENAME.getName());
      // copyFileIfExists(proj.SAMPLE_QC_FILENAME.getName());
      // copyFileIfExists(proj.MOSAIC_RESULTS_FILENAME.getName());
      // copyFileIfExists(proj.SEXCHECK_RESULTS_FILENAME.getName());
      // copyFileIfExists(proj.INTENSITY_PC_FILENAME.getName());
      // copyFileIfExists(proj.CLUSTER_FILTER_COLLECTION_FILENAME.getName());
      // copyFileIfExists(proj.ANNOTATION_FILENAME.getName());
      // multi-file Copies
      // copyIfFilesExists(proj.TWOD_LOADED_FILENAMES);
      copyIfFilesExists(proj.TWOD_LOADED_FILENAMES.getName());

      // copyIfFilesExists(proj.INDIVIDUAL_CNV_LIST_FILENAMES);//copying to both marker and sample
      // focus in case of marker focused cnvs
      copyIfFilesExists(proj.INDIVIDUAL_CNV_LIST_FILENAMES.getName());// copying to both marker and
                                                                      // sample focus in case of
                                                                      // marker focused cnvs

      copyStratResults(proj, this);
      copyGeneTrack(proj, this);

    } else {
      proj.getLog().reportTimeError(demoProjectDirectory
                                    + " exists and the overwrite option was not flagged, halting");
      fail = true;
    }
  }

  public boolean isFail() {
    return fail;
  }

}
