package org.genvisis.cnv.affy;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Collection;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import org.genvisis.cnv.analysis.CentroidCompute;
import org.genvisis.cnv.filesys.Compression;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.MarkerLookup;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.filesys.SampleList;
import org.genvisis.cnv.manage.Markers;
import org.genvisis.cnv.manage.TransposeData;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.Elision;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.ext;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public class AffyParsingPipeline {

  private Project proj;
  private String callFile;
  private String confFile;
  private String intFile;

  private String delim = "\t";

  private BufferedReader confReader;
  private BufferedReader callReader;
  private BufferedReader sigReader;

  private String[] samples;
  int numSamples;
  Map<Marker, String> markerFileMap;
  Map<String, String[]> markersInFileMap;
  Map<Marker, Integer> markerIndexInFileMap;
  private long fingerprint;
  private int maxSize = DEFAULT_MDRAF_SIZE;

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

  public void setMaxMDRAFSize(int size) {
    this.maxSize = size;
  }

  private void loadAndSortMarkers() throws IOException {
    Set<String> allMarkers = new HashSet<>();
    BufferedReader reader = Files.getAppropriateReader(proj.PROJECT_DIRECTORY.getValue()
                                                       + proj.PROJECT_NAME.getValue()
                                                       + ".probesetIdsAll.txt");
    String line = null;
    boolean before = true;
    while ((line = reader.readLine()) != null) {
      // skip comment lines
      if (line.startsWith("#%")) continue;
      if (before) {
        // skip header
        before = false;
        continue;
      }
      String snp = line.split("\t")[0];
      if (snp.endsWith("-A") || snp.endsWith("-B")) {
        snp = snp.substring(0, snp.length() - 2);
      }
      allMarkers.add(snp);
    }
    reader.close();

    // create naive MarkerSet file 
    Markers.orderMarkers(allMarkers.toArray(new String[allMarkers.size()]), proj);
  }

  private void binMarkers(long numMkrsPerFile) {
    List<Marker> markers = proj.getMarkerSet().markersAsList();

    ArrayListMultimap<String, String> markerListMap = ArrayListMultimap.create();

    markerFileMap = Maps.newHashMap();
    Hashtable<String, String> lookup = new Hashtable<>();
    markerIndexInFileMap = Maps.newHashMap();

    for (int i = 0; i < markers.size(); i++) {
      long markerFileIndex = i / numMkrsPerFile;
      String mkrFile = proj.MARKER_DATA_DIRECTORY.getValue() + "markers." + markerFileIndex
                       + ".mdRAF";
      int markerIndexInFile = (int) (i % numMkrsPerFile);
      markerFileMap.put(markers.get(i), mkrFile);
      markerIndexInFileMap.put(markers.get(i), markerIndexInFile);
      lookup.put(markers.get(i).getName(),
                 ext.removeDirectoryInfo(mkrFile) + "\t" + markerIndexInFile);
      markerListMap.put(mkrFile, markers.get(i).getName());
    }
    markersInFileMap = Maps.newHashMap();
    for (Entry<String, Collection<String>> m : markerListMap.asMap().entrySet()) {
      markersInFileMap.put(m.getKey(), m.getValue().toArray(new String[m.getValue().size()]));
    }
    new MarkerLookup(lookup).serialize(proj.MARKERLOOKUP_FILENAME.getValue());
  }

  public void run() {
    long startTimeNanos = System.nanoTime();
    // ensure directory exists
    new File(proj.MARKER_DATA_DIRECTORY.getValue()).mkdirs();

    // open files and parse header to determine samples
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

    // load a list of all markers, sort into chr/pos order, and create the MarkerSet file
    try {
      loadAndSortMarkers();
    } catch (IOException e3) {
      proj.getLog()
          .reportError("Encountered problem when reading the confidences file.  Parsing will fail; error message: "
                       + e3.getMessage());
      return;
    }

    boolean canXYBeNegative = false;
    byte nullStatus = Sample.computeNullStatus(new float[0], new float[0], new float[0],
                                               new float[0], new float[0], new byte[0], null,
                                               canXYBeNegative);
    long bytesPerMarker = ((long) numSamples) * Sample.getNBytesPerSampleMarker(nullStatus);
    long maxSizeBytes = ((long) maxSize) * 1024 * 1024;
    long mem = Math.min(maxSizeBytes, (long) (Runtime.getRuntime().maxMemory() * 0.8));
    long markersInMemory = mem / ((long) bytesPerMarker);
    long mkrsInProj = proj.getMarkerSet().markersAsList().size();
    long numMarkersPerFile = Math.min((long) Math.min(Integer.MAX_VALUE * .75, markersInMemory),
                                      mkrsInProj);

    proj.getLog()
        .reportTime("Parsing data into " + numMarkersPerFile
                    + (numMarkersPerFile == mkrsInProj ? " (all)" : "")
                    + " markers per mdRAF file.");

    // assign each marker to a file
    binMarkers(numMarkersPerFile);

    Map<String, Integer> markerNameBytesLengthMap = Maps.newHashMap();
    Map<String, RandomAccessFile> rafMap = Maps.newHashMap();
    Map<String, Long> afterLastMarkerPosition = Maps.newHashMap();

    // setup out of range value tables for each marker file
    Map<String, Hashtable<String, Float>> oorTables = Maps.newHashMap();
    for (String file : markerFileMap.values()) {
      oorTables.put(file, new Hashtable<>());
      afterLastMarkerPosition.put(file, 0L);
    }

    Map<String, Marker> markerNameMap = proj.getMarkerSet().getMarkerNameMap();
    MarkerData md = null;
    byte[] mkrBytes;
    try {
      while ((md = parseLine()) != null) {
        Marker marker = markerNameMap.get(md.getMarkerName());
        if (marker == null) {
          proj.getLog().reportTimeWarning("No Marker object found for " + md.getMarkerName());
          continue;
        }

        String mkrFile = markerFileMap.get(marker);
        if (mkrFile == null) {
          proj.getLog().reportTimeWarning("No file found for " + marker.getName());
          continue;
        }

        RandomAccessFile raf = rafMap.get(mkrFile);
        if (raf == null) {
          raf = new RandomAccessFile(mkrFile, "rw");
          rafMap.put(mkrFile, raf);

          String[] names = markersInFileMap.get(mkrFile);
          byte[] mkrNmBytes = Compression.objToBytes(names);
          byte[] param = TransposeData.getParameterSectionForMdRaf(numSamples, names.length,
                                                                   nullStatus, fingerprint,
                                                                   mkrNmBytes);
          raf.seek(0);
          raf.write(param);
          markerNameBytesLengthMap.put(mkrFile, mkrNmBytes.length);
        }

        int indInFile = markerIndexInFileMap.get(marker);

        long seek = TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN;
        seek += markerNameBytesLengthMap.get(mkrFile);
        seek += ((long) indInFile) * ((long) bytesPerMarker);

        try {
          mkrBytes = md.compress(indInFile, nullStatus, oorTables.get(mkrFile), canXYBeNegative);
        } catch (Elision e1) {
          proj.getLog()
              .reportError("Unexpected error while compressing data: " + e1.getMessage()
                           + ". Parsing will stop.  Please remove any existing .mdRAF files and try again.");
          try {
            closeReaders();
          } catch (IOException e2) {}
          return;
        }

        if (raf.getFilePointer() != seek) {
          raf.seek(seek);
        }
        raf.write(mkrBytes);
        if (seek + mkrBytes.length > afterLastMarkerPosition.get(mkrFile)) {
          afterLastMarkerPosition.put(mkrFile, seek + mkrBytes.length);
        }
      }
      try {
        closeReaders();
      } catch (IOException e) {
        proj.getLog()
            .reportTimeWarning("Exception occurred when closing input files.  This may not be a problem.");
      }

      // use a set so we only process each file once
      HashSet<String> files = new HashSet<>(markerFileMap.values());
      for (String mkrFile : files) {
        Hashtable<String, Float> oorTable = oorTables.get(mkrFile);
        RandomAccessFile raf = rafMap.get(mkrFile);
        if (raf.getFilePointer() != afterLastMarkerPosition.get(mkrFile)) {
          raf.seek(afterLastMarkerPosition.get(mkrFile));
        }
        if (oorTable.isEmpty()) {
          raf.write(Compression.intToBytes(0));
        } else {
          byte[] oorBytes = Compression.objToBytes(oorTable);
          raf.write(Compression.intToBytes(oorBytes.length));
          raf.write(oorBytes);
        }
        raf.close();
      }

    } catch (IOException e) {
      proj.getLog()
          .reportError("Uexpected error occurred while reading input files: " + e.getMessage()
                       + ". Parsing will stop.  Please remove any existing .mdRAF files and try again.");
      proj.getLog().reportException(e);
      try {
        closeReaders();
      } catch (IOException e2) {}
    }

    proj.getLog()
        .reportTime("Parsing finished in " + ext.getTimeElapsedNanos(startTimeNanos) + ".");
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
    String[] confHeader = line.trim().split(delim, -1);

    line = null;
    while ((line = callReader.readLine()) != null && line.charAt(0) == '#') {
      // seek to header
    }
    String[] callHeader = line.trim().split(delim, -1);

    line = null;
    while ((line = sigReader.readLine()) != null && line.charAt(0) == '#') {
      // seek to header
    }
    String[] sigHeader = line.trim().split(delim, -1);

    ensureHeaderContentsSame(confHeader, callHeader, sigHeader);
    samples = ArrayUtils.subArray(confHeader, 1);
    SampleList sl = new SampleList(samples);
    sl.serialize(proj.SAMPLELIST_FILENAME.getValue());
    Files.writeArray(samples, proj.PROJECT_DIRECTORY.getValue() + "ListOfSamples.txt");
    Files.writeArray(samples, proj.MARKER_DATA_DIRECTORY.getValue(true, false)
                              + TransposeData.TEMP_SAMPLES_FILE);
    numSamples = samples.length;
    fingerprint = org.genvisis.cnv.filesys.MarkerSet.fingerprint(samples);
  }

  private void ensureHeaderContentsSame(String[] confHeader, String[] callHeader,
                                        String[] sigHeader) {
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
      xs[i] = ys[i] = (float) (power2(sigs[i + 1]) / scale);
      bafs[i] = lrrs[i] = gcs[i] = 0;
      abGenos[i] = -1;
    }
    MarkerData md = new MarkerData(mkr, (byte) 0, 0, fingerprint, gcs, xRaws, yRaws, xs, ys, thetas,
                                   rs, bafs, lrrs, abGenos, forwardGenos);
    CentroidCompute centroid = md.getCentroid(null, null, true, 0, 0, null, true, proj.getLog());
    md = null;
    return new MarkerData(mkr, (byte) 0, 0, fingerprint, gcs, xRaws, yRaws, xs, ys, thetas, rs,
                          bafs, centroid.getRecomputedLRR(), abGenos, forwardGenos);
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
    byte[] forwardGenos = null;

    double scale = proj.XY_SCALE_FACTOR.getValue();
    for (int i = 0; i < numSamples; i++) {
      abGenos[i] = (byte) Integer.parseInt(calls[i + 1]);
      gcs[i] = Float.parseFloat(confs[i + 1]);
      xs[i] = (float) (power2(sigsA[i + 1]) / scale);
      ys[i] = (float) (power2(sigsB[i + 1]) / scale);
    }
    MarkerData md = new MarkerData(mkr, (byte) 0, 0, fingerprint, gcs, xRaws, yRaws, xs, ys, thetas,
                                   rs, bafs, lrrs, abGenos, forwardGenos);
    CentroidCompute centroid = md.getCentroid(null, null, false, 0, 0, null, true, proj.getLog());
    md = null;
    return new MarkerData(mkr, (byte) 0, 0, fingerprint, gcs, xRaws, yRaws, xs, ys, thetas, rs,
                          centroid.getRecomputedBAF(), centroid.getRecomputedLRR(), abGenos,
                          forwardGenos);
  }

  private static double power2(String signal) {
    return (Math.pow(2, Double.parseDouble(signal)));
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
  private static final String ARG_MDRAF_MAX_SIZE_MB = "size";
  private static final String DESC_CALL_FILE = "A file containing genotype calls created by AffyPowerTools";
  private static final String DESC_CONF_FILE = "A file containing confidence values created by AffyPowerTools";
  private static final String DESC_NORM_INT_FILE = "A file containing normalized intensities, created by AffyPowerTools";
  private static final String DESC_MDRAF_MAX_SIZE_MB = "Maximum size of marker data files in megabytes (default 2048, i.e. 2Gb)";

  private static final int DEFAULT_MDRAF_SIZE = 2048;

  public static void main(String[] args) {
    CLI cli = new CLI(AffyParsingPipeline.class);

    cli.addArg(CLI.ARG_PROJ, CLI.DESC_PROJ);
    cli.addArg(ARG_CALL_FILE, DESC_CALL_FILE);
    cli.addArg(ARG_CONF_FILE, DESC_CONF_FILE);
    cli.addArg(ARG_NORM_INT_FILE, DESC_NORM_INT_FILE);
    cli.addArg(ARG_MDRAF_MAX_SIZE_MB, DESC_MDRAF_MAX_SIZE_MB, false);

    cli.parseWithExit(args);

    AffyParsingPipeline app = new AffyParsingPipeline();
    app.setProject(new Project(cli.get(CLI.ARG_PROJ)));
    app.setConfidencesFile(cli.get(ARG_CONF_FILE));
    app.setGenotypeCallFile(cli.get(ARG_CALL_FILE));
    app.setNormIntensitiesFile(cli.get(ARG_NORM_INT_FILE));
    app.setMaxMDRAFSize(cli.has(ARG_MDRAF_MAX_SIZE_MB) ? cli.getI(ARG_MDRAF_MAX_SIZE_MB)
                                                       : DEFAULT_MDRAF_SIZE);
    app.run();

  }

}
