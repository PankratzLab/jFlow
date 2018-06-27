package org.genvisis.affy;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import org.genvisis.CLI;
import org.genvisis.cnv.analysis.CentroidCompute;
import org.genvisis.cnv.filesys.Compression;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerLookup;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.filesys.SampleList;
import org.genvisis.cnv.manage.MarkerDataLoader;
import org.genvisis.cnv.manage.Markers;
import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import com.google.common.collect.Sets;
import com.googlecode.charts4j.collect.Lists;

public class AffyParsingPipeline {

  private static final int MAX_MKRS_PER_MDRAF = 1500;

  private Project proj;
  private String callFile;
  private String confFile;
  private String intFile;

  private String delim = "\t";

  private BufferedReader confReader;
  private BufferedReader callReader;
  private BufferedReader sigReader;

  private String[] confHeader;
  private String[] callHeader;
  private String[] sigHeader;

  private String[] samples;
  int numSamples;
  private long fingerprint;

  public void setProject(Project proj) {
    this.proj = proj;
  }

  public void setGenotypeCallFile(String callFile) {
    this.callFile = callFile;
  }

  public void setConfidencesFile(String confFile) {
    this.confFile = confFile;
  }

  public void setNormIntensitiesFile(String quantNormFile) {
    this.intFile = quantNormFile;
  }

  public void run() {
    try {
      initReaders();
    } catch (IOException e) {
      proj.getLog().reportError("Couldn't open files for reading: \n\t" + confFile + "\n\t"
                                + callFile + "\n\t" + intFile);
      try {
        closeReaders();
      } catch (IOException e1) {}
      return;
    }

    boolean canXYBeNegative = true;
    byte nullStatus = Sample.updateNullStatus(new float[0], new float[0], new float[0],
                                              new float[0], new float[0], new byte[0], new byte[0],
                                              canXYBeNegative);
    int bytesPerMarker = numSamples * Sample.getNBytesPerSampleMarker(nullStatus);
    long mem = (long) (Runtime.getRuntime().maxMemory() * 0.8);
    long markersInMemory = mem / (long) bytesPerMarker;
    int numMarkersPerFile = Math.min(MAX_MKRS_PER_MDRAF, (int) markersInMemory);

    Hashtable<String, String> lookup = new Hashtable<>();

    int markerFileIndex = 0;
    String mkrFile = proj.MARKER_DATA_DIRECTORY.getValue() + "markers." + markerFileIndex
                     + ".mdRAF";

    new File(proj.MARKER_DATA_DIRECTORY.getValue()).mkdirs();

    List<String> allMarkers = Lists.newArrayList();
    String[] names = new String[numMarkersPerFile];
    byte[][] mkrBytes = new byte[numMarkersPerFile][];
    int mkrCount = 0;
    Hashtable<String, Float> oorTable = new Hashtable<>();
    MarkerData md = null;

    String[] allSampsInProj = new String[numSamples];
    for (int i = 0; i < numSamples; i++) {
      allSampsInProj[i] = "Sample_" + i;
    }

    int dump = 5;

    try {
      while ((md = parseLine()) != null) {
        names[mkrCount] = md.getMarkerName();
        allMarkers.add(md.getMarkerName());
        try {
          if (dump > 0) {
            md.dump(null, ext.parseDirectoryOfFile(mkrFile) + ext.rootOf(mkrFile) + "_dump_original"
                          + ".xln",
                    allSampsInProj, true, new Logger());
          }
          mkrBytes[mkrCount] = md.compress(mkrCount, nullStatus, oorTable, canXYBeNegative);
          if (dump > 0) {
            MarkerData md1 = MarkerDataLoader.parseMarkerData((byte) -1, -1, nullStatus,
                                                              fingerprint, true, true, true, true,
                                                              true, allSampsInProj, names, mkrCount,
                                                              oorTable, mkrBytes[mkrCount],
                                                              new Logger());
            md1.dump(null, ext.parseDirectoryOfFile(mkrFile) + ext.rootOf(mkrFile)
                           + "_dump_compressed" + ".xln",
                     allSampsInProj, true, new Logger());
            dump--;
          }
        } catch (Elision e1) {
          proj.getLog()
              // TODO FOR-REVIEW remove existing mdRAF files automatically?
              .reportError("Unexpected error while compressing data: " + e1.getMessage()
                           + ". Parsing will stop.  Please remove any existing .mdRAF files and try again.");
          try {
            closeReaders();
          } catch (IOException e2) {}
          return;
        }
        mkrCount++;
        if (mkrCount == numMarkersPerFile) {
          try {
            writeRAF(numMarkersPerFile, lookup, nullStatus, mkrFile, names, mkrBytes, mkrCount,
                     oorTable, bytesPerMarker);
          } catch (IOException e) {
            proj.getLog()
                // TODO FOR-REVIEW remove existing mdRAF files automatically?
                .reportError("Uexpected error occurred while writing marker file: " + e.getMessage()
                             + ". Parsing will stop.  Please remove any existing .mdRAF files and try again.");
            try {
              closeReaders();
            } catch (IOException e2) {}
            return;
          }

          names = new String[numMarkersPerFile];
          mkrBytes = new byte[numMarkersPerFile][];
          oorTable.clear();
          mkrCount = 0;
          markerFileIndex++;
          mkrFile = proj.MARKER_DATA_DIRECTORY.getValue() + "markers." + markerFileIndex + ".mdRAF";
        }
      }
      try {
        writeRAF(numMarkersPerFile, lookup, nullStatus, mkrFile, names, mkrBytes, mkrCount,
                 oorTable, bytesPerMarker);
      } catch (IOException e) {
        proj.getLog()
            // TODO FOR-REVIEW remove existing mdRAF files automatically?
            .reportError("Uexpected error occurred while writing marker file: " + e.getMessage()
                         + ". Parsing will stop.  Please remove any existing .mdRAF files and try again.");
        try {
          closeReaders();
        } catch (IOException e2) {}
        return;
      }

      names = null;
      mkrBytes = null;
      oorTable = null;
      try {
        closeReaders();
      } catch (IOException e) {
        proj.getLog()
            .reportTimeWarning("Exception occurred when closing input files.  This may not be a problem.");
      }
      writeMarkerLookupAndList(allMarkers, lookup);
    } catch (IOException e) {
      proj.getLog()
          // TODO FOR-REVIEW remove existing mdRAF files automatically?
          .reportError("Uexpected error occurred while reading input files: " + e.getMessage()
                       + ". Parsing will stop.  Please remove any existing .mdRAF files and try again.");
      try {
        closeReaders();
      } catch (IOException e2) {}
    }
  }

  private void writeMarkerLookupAndList(List<String> allMarkers, Hashtable<String, String> lookup) {
    new MarkerLookup(lookup).serialize(proj.MARKERLOOKUP_FILENAME.getValue());
    Markers.orderMarkers(allMarkers.toArray(new String[allMarkers.size()]), proj, proj.getLog());
  }

  private void writeRAF(int numMarkersPerFile, Hashtable<String, String> lookup, byte nullStatus,
                        String mkrFile, String[] names, byte[][] mkrBytes, int mkrsToWrite,
                        Hashtable<String, Float> oorTable,
                        int numBytesPerMarker) throws IOException {
    byte[] mkrNmBytes = Compression.objToBytes(names);
    byte[] param = TransposeData.getParameterSectionForMdRaf(numSamples, numMarkersPerFile,
                                                             nullStatus, fingerprint, mkrNmBytes);
    RandomAccessFile raf = new RandomAccessFile(mkrFile, "rw");
    raf.seek(0);
    raf.write(param);
    for (int i = 0; i < mkrsToWrite; i++) {
      long seek = TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN + mkrBytes.length
                  + i * numBytesPerMarker;
      // seek to location of marker in file, as we may be writing out of order
      if (raf.getFilePointer() != seek) {
        raf.seek(seek);
      }
      raf.write(mkrBytes[i]);
      lookup.put(names[i], ext.removeDirectoryInfo(mkrFile) + "\t" + i);
    }

    byte[] oorBytes = Compression.objToBytes(oorTable);
    raf.write(Compression.intToBytes(oorBytes.length));
    raf.write(oorBytes);

    raf.close();
    proj.getLog().reportTime("Wrote marker file " + mkrFile);
  }

  @SuppressWarnings("deprecation")
  void initReaders() throws IOException {
    confReader = Files.getAppropriateReader(confFile);
    callReader = Files.getAppropriateReader(callFile);
    sigReader = Files.getAppropriateReader(intFile);

    String line;

    line = null;
    while ((line = confReader.readLine()) != null && line.charAt(0) == '#') {
      // seek to header
    }
    confHeader = line.trim().split(delim, -1);

    line = null;
    while ((line = callReader.readLine()) != null && line.charAt(0) == '#') {
      // seek to header
    }
    callHeader = line.trim().split(delim, -1);

    line = null;
    while ((line = sigReader.readLine()) != null && line.charAt(0) == '#') {
      // seek to header
    }
    sigHeader = line.trim().split(delim, -1);

    ensureHeaderContentsSame();
    samples = ArrayUtils.subArray(confHeader, 1);
    SampleList sl = new SampleList(samples);
    sl.serialize(proj.SAMPLELIST_FILENAME.getValue());
    Files.writeArray(samples, proj.PROJECT_DIRECTORY.getValue() + "ListOfSamples.txt");
    Files.writeArray(samples, proj.MARKER_DATA_DIRECTORY.getValue(true, false)
                              + TransposeData.TEMP_SAMPLES_FILE);
    numSamples = samples.length;
    fingerprint = org.genvisis.cnv.filesys.MarkerSet.fingerprint(ArrayUtils.subArray(confHeader,
                                                                                     1));
  }

  private void ensureHeaderContentsSame() {
    if (confHeader.length != callHeader.length || confHeader.length != sigHeader.length
        || callHeader.length != sigHeader.length) {
      throw new IllegalStateException("File headers have different lengths: [Conf: "
                                      + confHeader.length + "], [Call: " + callHeader.length
                                      + "], [Sig: " + sigHeader.length + "].");
    }

    Set<String> confs = Sets.newHashSet(confHeader);
    Set<String> calls = Sets.newHashSet(callHeader);
    Set<String> sigs = Sets.newHashSet(sigHeader);

    boolean diff = Sets.symmetricDifference(confs, calls).size() > 0
                   || Sets.symmetricDifference(confs, sigs).size() > 0
                   || Sets.symmetricDifference(calls, sigs).size() > 0;
    if (diff) {
      throw new IllegalStateException("File headers have different values.");
    }
  }

  MarkerData parseLine() throws IOException {
    String confLine = confReader.readLine();
    String callLine = callReader.readLine();
    String sigLineA = sigReader.readLine();
    if (confLine == null && callLine == null && sigLineA != null) {
      return parseCNMarker(sigLineA);
    }
    String sigLineB = sigReader.readLine();

    if (confLine == null && callLine == null && sigLineA == null && sigLineB == null) {
      return null;
    }

    return parseSNPMarker(confLine, callLine, sigLineA, sigLineB);
  }

  private MarkerData parseCNMarker(String sigLine) {
    String[] sigs = sigLine.trim().split(delim, -1);

    String mkr = sigs[0];

    float[] gcs = new float[numSamples];
    float[] xRaws = null;
    float[] yRaws = null;
    float[] xs = new float[numSamples];
    float[] ys = new float[numSamples];
    float[] thetas = null;
    float[] rs = null;
    float[] bafs = new float[numSamples];
    float[] lrrs = new float[numSamples];
    byte[] abGenos = new byte[numSamples];
    byte[] forwardGenos = null;

    double scale = proj.XY_SCALE_FACTOR.getValue();
    for (int i = 0; i < numSamples; i++) {
      xs[i] = ys[i] = (float) (AffySNP6Tables.power2(sigs[i + 1]) / scale);
      bafs[i] = lrrs[i] = gcs[i] = 0;
      abGenos[i] = -1;
    }
    return new MarkerData(mkr, (byte) 0, 0, fingerprint, gcs, xRaws, yRaws, xs, ys, thetas, rs,
                          bafs, lrrs, abGenos, forwardGenos);
  }

  private MarkerData parseSNPMarker(String confLine, String callLine, String sigLineA,
                                    String sigLineB) {
    String[] confs = confLine.trim().split(delim, -1);
    String[] calls = callLine.trim().split(delim, -1);
    String[] sigsA = sigLineA.trim().split(delim, -1);
    String[] sigsB = sigLineB.trim().split(delim, -1);
    ensureSame(confs[0], calls[0], sigsA[0], sigsB[0]);

    String mkr = calls[0];

    float[] gcs = new float[numSamples];
    float[] xRaws = null;
    float[] yRaws = null;
    float[] xs = new float[numSamples];
    float[] ys = new float[numSamples];
    float[] thetas = null;
    float[] rs = null;
    float[] bafs = null;
    float[] lrrs = null;
    byte[] abGenos = new byte[numSamples];
    byte[] forwardGenos = ArrayUtils.byteArray(numSamples, (byte) 0);

    double scale = proj.XY_SCALE_FACTOR.getValue();
    for (int i = 0; i < numSamples; i++) {
      abGenos[i] = (byte) Integer.parseInt(calls[i + 1]);
      gcs[i] = Float.parseFloat(confs[i + 1]);
      xs[i] = (float) (AffySNP6Tables.power2(sigsA[i + 1]) / scale);
      ys[i] = (float) (AffySNP6Tables.power2(sigsB[i + 1]) / scale);
    }
    MarkerData md = new MarkerData(mkr, (byte) 0, 0, fingerprint, gcs, xRaws, yRaws, xs, ys, thetas,
                                   rs, bafs, lrrs, abGenos, forwardGenos);
    CentroidCompute centroid = md.getCentroid(null, null, false, 0, 0, null, true, proj.getLog());
    md = null;
    return new MarkerData(mkr, (byte) 0, 0, fingerprint, gcs, xRaws, yRaws, xs, ys, thetas, rs,
                          centroid.getRecomputedBAF(), centroid.getRecomputedLRR(), abGenos,
                          forwardGenos);
  }

  private void ensureSame(String conf, String call, String sigA, String sigB) {
    if (conf.equals(call) && sigA.startsWith(conf) && sigB.startsWith(conf)) {
      return;
    }
    throw new IllegalStateException("Affymetrix parsing became desynchronised - expected these to all be the same: [Call: "
                                    + call + "], [Conf: " + conf + "], [SigA: " + sigA + ", SigB: "
                                    + sigB + "].");
  }

  void closeReaders() throws IOException {
    confReader.close();
    callReader.close();
    sigReader.close();
  }

  private static final String ARG_CALL_FILE = "call";
  private static final String ARG_CONF_FILE = "conf";
  private static final String ARG_NORM_INT_FILE = "norm";
  private static final String DESC_CALL_FILE = "A file containing genotype calls created by AffyPowerTools";
  private static final String DESC_CONF_FILE = "A file containing confidence values created by AffyPowerTools";
  private static final String DESC_NORM_INT_FILE = "A file containing normalized intensities, created by AffyPowerTools";

  public static void main(String[] args) {
    CLI cli = new CLI(AffyParsingPipeline.class);

    cli.addArg(CLI.ARG_PROJ, CLI.DESC_PROJ);
    cli.addArg(ARG_CALL_FILE, DESC_CALL_FILE);
    cli.addArg(ARG_CONF_FILE, DESC_CONF_FILE);
    cli.addArg(ARG_NORM_INT_FILE, DESC_NORM_INT_FILE);

    cli.parseWithExit(args);

    AffyParsingPipeline app = new AffyParsingPipeline();
    app.setProject(new Project(cli.get(CLI.ARG_PROJ)));
    app.setConfidencesFile(cli.get(ARG_CONF_FILE));
    app.setGenotypeCallFile(cli.get(ARG_CALL_FILE));
    app.setNormIntensitiesFile(cli.get(ARG_NORM_INT_FILE));
    app.run();

  }

}
