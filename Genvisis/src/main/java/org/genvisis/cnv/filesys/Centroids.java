package org.genvisis.cnv.filesys;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.Arrays;
import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

import org.genvisis.cnv.analysis.CentroidCompute;
import org.genvisis.cnv.hmm.PFB;
import org.genvisis.cnv.manage.MarkerDataLoader;
import org.genvisis.cnv.manage.TextExport;
import org.genvisis.cnv.prop.Property;
import org.genvisis.cnv.qc.SexChecks;
import org.genvisis.cnv.var.SampleData;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.SerializedFiles;
import org.pankratzlab.common.WorkerTrain;
import org.pankratzlab.common.WorkerTrain.AbstractProducer;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.stats.Maths;

/**
 * @author lane0212
 */
public class Centroids implements Serializable, TextExport {

  public static final long serialVersionUID = 1L;
  public static final String[] ILLUMINA_CENTROID_SUFFIXES = {"Name", "AA T Mean", "AA R Mean",
                                                             "AB T Mean", "AB R Mean", "BB T Mean",
                                                             "BB R Mean"};

  private final float[][][] centroids; // marker, genotype (0=AA, 1=AB, 2=BB), coordinates (0=Mean
                                       // Theta, 1=Mean R) (a.k.a. follows the suffix order above)
  private final long fingerprint;

  public Centroids(float[][][] centroids, long fingerprint) {
    this.centroids = centroids;
    this.fingerprint = fingerprint;
  }

  public float[][][] getCentroids() {
    return centroids;
  }

  public long getFingerprint() {
    return fingerprint;
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename);
  }

  @Override
  public void exportToText(Project proj, String outputFile) {
    PrintWriter writer;
    float[][][] centroids;
    MarkerSetInfo markerSet;
    String[] markerNames;
    Logger log = proj.getLog();

    markerSet = proj.getMarkerSet();
    markerNames = markerSet.getMarkerNames();
    centroids = getCentroids();

    if (markerNames.length != centroids.length) {
      log.reportError("Error - mismatched number of markers in centroid object and the project's marker set; aborting");
      return;
    }

    if (markerSet.getFingerprint() != getFingerprint()) {
      log.reportError("Error - mismatched marker fingerprints in centroid object and the project's marker set ; aborting");
      return;
    }

    try {
      writer = Files.openAppropriateWriter(outputFile);
      writer.println("marker_fingerprint=" + getFingerprint());
      writer.println("MarkerName\tAA_Theta_Mean\tAA_R_Mean\tAB_Theta_Mean\tAB_R_Mean\tBB_Theta_Mean\tBB_R_Mean");
      for (int i = 0; i < markerNames.length; i++) {
        writer.print(markerNames[i]);
        for (int j = 0; j < 3; j++) {
          if (centroids[i][j] == null) {
            writer.print("\t.\t.");
          } else {
            writer.print("\t" + centroids[i][j][0] + "\t" + centroids[i][j][1]);
          }
        }
        writer.println();
      }
      writer.close();
    } catch (Exception e) {
      log.reportException(e);
    }
  }

  public static Centroids load(String filename) {
    return (Centroids) SerializedFiles.readSerial(filename, true);
  }

  public static float calcR(float x, float y) {
    return (float) (Math.max(x, 0.0001) + Math.max(y, 0.0001));
  }

  public static float calcTheta(float x, float y) {
    return (float) (Math.atan(Math.max(y, 0.0001) / Math.max(x, 0.0001)) * 2 / Math.PI);
  }

  public static float calcBAF(float theta, float[][] centroids) {
    if (centroids[0] != null && theta < centroids[0][0]) {
      return 0;
    } else if (centroids[1] != null && theta < centroids[1][0]) {
      if (centroids[0] == null) {
        return 0.50f;
      } else {
        return 0.5f * (theta - centroids[0][0]) / (centroids[1][0] - centroids[0][0]);
      }
    } else if (centroids[2] != null && theta < centroids[2][0]) {
      if (centroids[1] == null) {
        return 1.0f;
      } else {
        return 0.5f + 0.5f * (theta - centroids[1][0]) / (centroids[2][0] - centroids[1][0]);
      }
    } else {
      if (centroids[2] == null) {
        return 0.50f;
      } else {
        return 1;
      }
    }
  }

  public static float calcLRR(float theta, float r, float[][] centroids) {
    float estimatedR;

    // if (centroids[2][1] < 0.0000001) {
    // centroids[2][0] = 1;
    // }

    if (centroids[0] != null && theta < centroids[0][0]
        || (centroids[1] == null && centroids[2] == null)) {
      estimatedR = centroids[0][1];
    } else if (centroids[1] != null && theta < centroids[1][0]) {
      if (centroids[0] == null) {
        estimatedR = centroids[1][1];
      } else {
        estimatedR = centroids[0][1]
                     + (theta - centroids[0][0]) * (centroids[1][1] - centroids[0][1])
                       / (centroids[1][0] - centroids[0][0]);
      }
    } else if (centroids[2] != null && theta < centroids[2][0]) {
      if (centroids[1] == null) {
        estimatedR = centroids[2][1];
      } else {
        estimatedR = centroids[1][1]
                     + (theta - centroids[1][0]) * (centroids[2][1] - centroids[1][1])
                       / (centroids[2][0] - centroids[1][0]);
      }
    } else {
      if (centroids[2] == null) {
        estimatedR = centroids[1][1];
      } else {
        estimatedR = centroids[2][1];
      }
    }

    return (float) Maths.log2(r / estimatedR);
  }

  public static void parseIlluminaCentroidsFromCSV(Project proj, String filename) {
    BufferedReader reader;
    String[] line, header;
    Hashtable<String, String> hash;
    int[] indices;
    float[][][] centroids;
    String[] markerNames;
    int index;
    boolean missing;
    MarkerSetInfo markerSet;

    markerSet = proj.getMarkerSet();
    markerNames = markerSet.getMarkerNames();
    centroids = new float[markerNames.length][][];

    hash = new Hashtable<>();
    for (int i = 0; i < markerNames.length; i++) {
      hash.put(markerNames[i], i + "");
    }

    try {
      reader = new BufferedReader(new FileReader(proj.PROJECT_DIRECTORY.getValue() + filename));
      header = reader.readLine().trim().split(",");
      indices = ArrayUtils.intArray(ILLUMINA_CENTROID_SUFFIXES.length, -1);
      for (int i = 0; i < header.length; i++) {
        index = ext.indexOfEndsWith(header[i], ILLUMINA_CENTROID_SUFFIXES, true);
        if (index >= 0) {
          indices[index] = i;
        }
      }
      if (ArrayUtils.min(indices) == -1) {
        System.err.println("Error - Need a column header ending with the following suffixes; missing at least one");
        System.err.println("        " + ArrayUtils.toStr(ILLUMINA_CENTROID_SUFFIXES, "  "));
      }
      // for (int i = 0; i < markerNames.length; i++) {
      while (reader.ready()) {
        line = reader.readLine().trim().split(",");
        if (!hash.containsKey(line[indices[0]])) {
          System.err.println("Error - marker '" + line[indices[0]]
                             + "' was not found in MarkerSet");
          System.exit(1);
        }
        index = Integer.parseInt(hash.get(line[indices[0]]));
        centroids[index] = new float[3][2];
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 2; k++) {
            centroids[index][j][k] = Float.parseFloat(line[indices[1 + j * 2 + k]]);
          }
        }
      }
      missing = false;
      for (int i = 0; i < centroids.length; i++) { // might want to generate an error log or display
                                                   // the number if greater than, say, 10
        if (centroids[i] == null) {
          if (!missing) {
            System.err.println("Error - did not find a centroid for the following markers:");
            missing = true;
          }
          System.err.println("  " + markerNames[i]);
        }

      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in "
                         + proj.PROJECT_DIRECTORY.getValue());
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }

    // new Centroids(centroids,
    // markerSet.getFingerprint()).serialize(proj.getFilename(proj.ORIGINAL_CENTROIDS_FILENAME));
    // Files.backup(proj.getFilename(proj.CUSTOM_CENTROIDS_FILENAME), "", "");
    // new Centroids(centroids,
    // markerSet.getFingerprint()).serialize(proj.getFilename(proj.CUSTOM_CENTROIDS_FILENAME));
    new Centroids(centroids,
                  markerSet.getFingerprint()).serialize(proj.ORIGINAL_CENTROIDS_FILENAME.getValue());
    Files.backup(proj.CUSTOM_CENTROIDS_FILENAME.getValue(), "", "");
    new Centroids(centroids,
                  markerSet.getFingerprint()).serialize(proj.CUSTOM_CENTROIDS_FILENAME.getValue());
  }

  public static void parseCentroidsFromGenotypes(Project proj, boolean[] samplesToBeUsed,
                                                 double missingnessThreshold) {
    String[] samples, markerNames;
    float[][][] centroids;
    float[][] centroid;
    int count;
    float[] thetas, rs;
    double[] meanThetas, meanRs;
    byte[] genotypes;
    int[] counts;
    SampleList sampleList;
    MarkerDataLoader markerDataLoader;
    MarkerData markerData;
    long time;
    MarkerSetInfo markerSet;

    time = new Date().getTime();
    System.out.println("Computing centroids from genotype means");
    sampleList = proj.getSampleList();
    samples = sampleList.getSamples();
    if (samples.length != samplesToBeUsed.length) {
      System.err.println("Error - mismatched number of samples in project versus sample mask");
      return;
    }
    markerSet = proj.getMarkerSet();
    markerNames = markerSet.getMarkerNames();
    centroids = new float[markerNames.length][][];

    markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames);
    time = new Date().getTime();
    for (int i = 0; i < markerNames.length; i++) {
      markerData = markerDataLoader.requestMarkerData(i);

      genotypes = markerData.getAbGenotypes();
      thetas = markerData.getThetas();
      rs = markerData.getRs();
      meanThetas = new double[5];
      meanRs = new double[5];
      counts = new int[5];

      for (int k = 0; k < samples.length; k++) {
        if (samplesToBeUsed[k] && !Float.isNaN(thetas[k]) && !Float.isNaN(rs[k])) {
          meanThetas[0] += thetas[k];
          meanRs[0] += rs[k];
          counts[0]++;
          meanThetas[genotypes[k] + 2] += thetas[k];
          meanRs[genotypes[k] + 2] += rs[k];
          counts[genotypes[k] + 2]++;
        }
      }
      for (int k = 0; k < 5; k++) {
        meanThetas[k] /= counts[k];
        meanRs[k] /= counts[k];
      }
      centroid = new float[3][];
      if (counts[1] >= counts[0] * missingnessThreshold) {
        for (int k = 0; k < 3; k++) {
          centroid[k] = new float[] {(float) meanThetas[0], (float) meanRs[0]};
        }
      } else {
        for (int k = 0; k < 3; k++) {
          if (counts[k + 2] > 0) {
            centroid[k] = new float[] {(float) meanThetas[k + 2], (float) meanRs[k + 2]};
          } else {
            centroid[k] = null;
          }
        }
      }

      centroids[i] = centroid;
      markerDataLoader.releaseIndex(i);
    }

    count = 0;
    for (int i = 0; i < centroids.length; i++) {
      if (centroids[i] == null) {
        if (count == 0) {
          System.out.println("The following marker(s) could not be computed:");
        }
        System.out.println(markerNames[i]);
        count++;
      }
    }
    if (count > 0) {
      System.out.println("Computed mean genotyped centroids for " + (centroids.length - count)
                         + " of " + centroids.length + " markers, " + count + " missing");
    } else {
      System.out.println("Computed mean genotyped centroids for all " + centroids.length
                         + " markers");
    }
    // new Centroids(centroids,
    // markerSet.getFingerprint()).serialize(proj.getFilename(proj.GENOTYPE_CENTROIDS_FILENAME));
    new Centroids(centroids,
                  markerSet.getFingerprint()).serialize(proj.GENOTYPE_CENTROIDS_FILENAME.getValue(true, false));
    System.out.println("Computation took " + ext.getTimeElapsed(time));
  }

  public static void generateChimeraCentroids(Project proj, String intensityOnlyFlagFile) {
    float[][][] cents, clusteredCents, unclusteredCents;
    Centroids clustered, unclustered;
    Hashtable<String, String> hash;
    String[] markerNames;
    MarkerSetInfo markerSet;
    boolean problem;
    String flag;
    long time;

    time = new Date().getTime();
    markerSet = proj.getMarkerSet();
    markerNames = markerSet.getMarkerNames();
    hash = HashVec.loadFileToHashString(proj.PROJECT_DIRECTORY.getValue() + intensityOnlyFlagFile,
                                        false);
    // if (!Files.exists(proj.getFilename(proj.GENOTYPE_CENTROIDS_FILENAME), jar)) {
    // System.err.println("Error - file '"+proj.getFilename(proj.GENOTYPE_CENTROIDS_FILENAME)+"'
    // does not exist in the project's data directory");
    if (!Files.exists(proj.GENOTYPE_CENTROIDS_FILENAME.getValue())) {
      System.err.println("Error - file '" + proj.GENOTYPE_CENTROIDS_FILENAME.getValue()
                         + "' does not exist in the project's data directory");
      return;
    }
    if (!Files.exists(proj.ORIGINAL_CENTROIDS_FILENAME.getValue())) {
      System.err.println("Error - file '" + proj.ORIGINAL_CENTROIDS_FILENAME.getValue()
                         + "' does not exist in the project's data directory");
      return;
    }
    clustered = Centroids.load(proj.GENOTYPE_CENTROIDS_FILENAME.getValue());
    unclustered = Centroids.load(proj.ORIGINAL_CENTROIDS_FILENAME.getValue());
    if (clustered.getFingerprint() != unclustered.getFingerprint()) {
      System.err.println("Error - the two centroid files cannot be merged, because they do not have the same fingerprint");
      return;
    }
    if (clustered.getFingerprint() != markerSet.getFingerprint()) {
      System.err.println("Error - the centroid files do not match the fingerprint of the current marker set");
      return;
    }

    problem = false;
    cents = new float[markerNames.length][][];
    clusteredCents = clustered.getCentroids();
    unclusteredCents = unclustered.getCentroids();
    for (int i = 0; i < markerNames.length; i++) {
      flag = hash.get(markerNames[i]);
      if (flag == null) {
        System.err.println("Error - no flag for marker '" + markerNames[i] + "'");
        problem = true;
      } else if (flag.equals("1")) {
        cents[i] = unclusteredCents[i];
      } else if (flag.equals("0")) {
        cents[i] = clusteredCents[i];
      } else {
        System.err.println("Error - invalid flag for marker '" + markerNames[i] + "'; found '"
                           + flag + "', expecting '1' or '0'");
        problem = true;
      }
    }

    if (problem) {
      System.err.println("Error - chimera centroids generation failed");
    } else {
      // new Centroids(cents,
      // markerSet.getFingerprint()).serialize(proj.getFilename(proj.CHIMERA_CENTROIDS_FILENAME));
      new Centroids(cents,
                    markerSet.getFingerprint()).serialize(proj.CHIMERA_CENTROIDS_FILENAME.getValue());
      System.out.println("Created chimera centroids in " + ext.getTimeElapsed(time));
    }
  }

  public static void recompute(Project proj, String centroidsFile) {
    recompute(proj, centroidsFile, false, 1);
  }

  /**
   * Recompute thread, applies centroids and saves new sample file
   */
  private static class RecomputeWorker implements Callable<Hashtable<String, Float>> {

    private final Project proj;
    private final String sample;
    private final Centroids centroids;
    private final boolean preserveBafs;

    public RecomputeWorker(Project proj, String sample, Centroids centroids, boolean preserveBafs) {
      super();
      this.proj = proj;
      this.sample = sample;
      this.centroids = centroids;
      this.preserveBafs = preserveBafs;
    }

    @Override
    public Hashtable<String, Float> call() throws Exception {
      Hashtable<String, Float> outliers = new Hashtable<>();
      Sample original = proj.getFullSampleFromRandomAccessFile(sample);
      Sample sample = new Sample(original.getSampleName(), original.getFingerprint(),
                                 original.getGCs(), original.getXs(), original.getYs(),
                                 preserveBafs ? original.getBAFs()
                                              : original.getBAFs(centroids.getCentroids()),
                                 original.getLRRs(centroids.getCentroids()),
                                 original.getForwardGenotypes(), original.getAB_Genotypes(),
                                 original.getCanXYBeNegative());
      sample.saveToRandomAccessFile(proj.SAMPLE_DIRECTORY.getValue(false, true)
                                    + original.getSampleName() + Sample.SAMPLE_FILE_EXTENSION,
                                    outliers, sample.getSampleName());
      return outliers;
    }

  }

  /**
   * Manages centroid application to each sample
   */
  private static class RecomputeProducer extends AbstractProducer<Hashtable<String, Float>> {

    private final Project proj;
    private final String[] samples;
    private final Centroids centroids;
    private final boolean preserveBafs;
    private int index;

    public RecomputeProducer(Project proj, String[] samples, Centroids centroids,
                             boolean preserveBafs) {
      super();
      this.proj = proj;
      this.samples = samples;
      this.centroids = centroids;
      this.preserveBafs = preserveBafs;
      index = 0;
    }

    @Override
    public boolean hasNext() {
      return index < samples.length;
    }

    @Override
    public Callable<Hashtable<String, Float>> next() {
      String currentSample = samples[index];
      RecomputeWorker worker = new RecomputeWorker(proj, currentSample, centroids, preserveBafs);
      index++;
      return worker;
    }
  }

  /**
   * @param proj
   * @param centroidsFile
   * @param preserveBafs bafs will not be recomputed from the centroids, useful if an atypical baf
   *          value is used
   */
  public static void recompute(Project proj, String centroidsFile, boolean preserveBafs,
                               int numThreads) {
    MarkerSetInfo markerSet;
    Centroids centroids;
    // Sample original, sample;
    String[] samples;
    // float[][][] cents;

    markerSet = proj.getMarkerSet();
    centroids = load(centroidsFile);
    if (centroids.getFingerprint() != markerSet.getFingerprint()) {
      System.err.println("Error - fingerprint for Centroids file '" + centroidsFile
                         + "' does not match the fingerprint for the current MarkerSet");
    }

    // cents = centroids.getCentroids();
    samples = proj.getSamples();
    Hashtable<String, Float> outliers = new Hashtable<>();
    RecomputeProducer producer = new RecomputeProducer(proj, samples, centroids, preserveBafs);
    try (WorkerTrain<Hashtable<String, Float>> train = new WorkerTrain<>(producer, numThreads, 10,
                                                                         proj.getLog())) {
      while (train.hasNext()) {
        Hashtable<String, Float> currentOutliers = train.next();
        outliers.putAll(currentOutliers);
      }
    }
    // for (int i = 0; i < samples.length; i++) {
    // original = proj.getFullSampleFromRandomAccessFile(samples[i]);
    // sample = new Sample(original.getSampleName(), original.getFingerprint(), original.getGCs(),
    // original.getXs(), original.getYs(), preserveBafs ? original.getBAFs() :
    // original.getBAFs(cents), original.getLRRs(cents), original.getForwardGenotypes(),
    // original.getAB_Genotypes(), original.getCanXYBeNegative());
    // sample.saveToRandomAccessFile(proj.SAMPLE_DIRECTORY.getValue(false, true) +
    // original.getSampleName() + Sample.SAMPLE_DATA_FILE_EXTENSION, outliers,
    // sample.getSampleName());
    // }
    if (outliers.size() > 0) {
      if (Files.exists(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser")) {
        Files.copyFile(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser",
                       ext.addToRoot(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser",
                                     ext.getTimestampForFilename()));
      }
    }
    SerializedFiles.writeSerial(outliers,
                                proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");

  }

  public static float[][] computeClusterCenters(MarkerData markerData, boolean[] samplesToBeUsed,
                                                double missingnessThreshold) {
    float[][] centers;
    float[] xs, ys;
    double[] meanXs, meanYs;
    byte[] genotypes;
    int[] counts;

    centers = new float[3][2];
    genotypes = markerData.getAbGenotypes();
    xs = markerData.getXs();
    ys = markerData.getYs();
    meanXs = new double[5];
    meanYs = new double[5];
    counts = new int[5];

    for (int k = 0; k < xs.length; k++) {
      if ((samplesToBeUsed == null || samplesToBeUsed[k]) && !Float.isNaN(xs[k])
          && !Float.isNaN(ys[k])) {
        meanXs[0] += xs[k];
        meanYs[0] += ys[k];
        counts[0]++;
        meanXs[genotypes[k] + 2] += xs[k];
        meanYs[genotypes[k] + 2] += ys[k];
        counts[genotypes[k] + 2]++;
      }
    }
    for (int k = 0; k < 5; k++) {
      meanXs[k] /= counts[k];
      meanYs[k] /= counts[k];
    }
    if (counts[1] >= counts[0] * missingnessThreshold) {
      for (int k = 0; k < 3; k++) {
        // centers[k] = new float[] { (float)meanXs[0], (float)meanYs[0] };
        centers[k] = null;
      }
    } else {
      for (int k = 0; k < 3; k++) {
        if (counts[k + 2] > 0) {
          centers[k] = new float[] {(float) meanXs[k + 2], (float) meanYs[k + 2]};
        } else {
          centers[k] = null;
        }
      }
    }

    return centers;
  }

  /**
   * @param proj
   * @param centFilename File path FROM THE PROJECT'S DIRECTORY
   * @param exportFilename File path FROM THE PROJECT'S DIRECTORY
   */
  public static void exportToText(Project proj, String centFilename, String exportFilename) {

    Centroids centObject;
    String dir;

    dir = proj.PROJECT_DIRECTORY.getValue();
    centObject = Centroids.load(dir + centFilename);

    centObject.exportToText(proj, dir + exportFilename);
  }

  public static enum CENTROID_STRATEGY {
    USE_ORIG_LRRS,
    USE_CENT_IF_EXISTS_OTHERWISE_ORIG,
    USE_CENT_IF_EXISTS_OTHERWISE_COMPUTE,
    COMPUTE_CENT;
  }

  @SuppressWarnings("unchecked")
  public static Centroids[] computeSexSpecificCentroids(final Project proj, String[] pfbFiles,
                                                        String[] centFiles, int threads) {
    PrintWriter writerM;
    PrintWriter writerF;
    MarkerSetInfo markerSet;
    String sampleDataFile;
    String[] allMarkers;
    String[] samples;
    String[] header;
    byte[] markerChrs;
    int[] markerPos;
    final boolean[] inclSampAll;
    final boolean[] inclSampFemales;
    final boolean[] inclSampMales;
    final boolean[] chromPlus11Markers = proj.getMarkerForChrsBoolean(new int[] {11, 23, 24, 25,
                                                                                 26});
    final boolean[] chromMarkers = proj.getMarkerForChrsBoolean(new int[] {23, 24, 25, 26});
    final int markerCount = ArrayUtils.booleanArraySum(chromPlus11Markers);
    int[] sampleSex;
    final float[][][] rawCentroidsFemale;
    final float[][][] rawCentroidsMale;
    Vector<String>[] markerLists;
    final Logger log = proj.getLog();
    ExecutorService computeHub;
    final ConcurrentLinkedQueue<Integer>[] markerIndexQueues;
    final Hashtable<Integer, String[][]> pfbInfo;
    final Hashtable<Integer, Integer>[] fullToTruncMarkerIndices;
    Hashtable<String, Vector<String>> sexData;
    final MarkerDataLoader[] markerDataLoaders;
    SampleData sampleData;
    int threadCount = threads == -1 ? Runtime.getRuntime().availableProcessors() : threads;

    markerSet = proj.getMarkerSet();
    sampleData = proj.getSampleData(false);

    inclSampAll = proj.getSamplesToInclude(null);
    if (!sampleData.hasExcludedIndividuals()) {
      log.report("Warning - there is no 'Exclude' column in SampleData.txt; centroids will be determined using all samples.");
    }
    samples = proj.getSamples();
    sampleDataFile = proj.SAMPLE_DATA_FILENAME.getValue(false, false);
    header = Files.getHeaderOfFile(sampleDataFile, proj.getLog());
    int sexInd = -1;
    for (int i = 0; i < header.length; i++) {
      if (("CLASS=" + SexChecks.EST_SEX_HEADER).toUpperCase().equals(header[i].toUpperCase())) {
        sexInd = i;
        break;
      }
    }
    if (sexInd == -1) {
      log.reportError("Error - no estimated sex found in sample data file - please run SexChecks with -check argument to generate the required data");
      return null;
    }
    sexData = HashVec.loadFileToHashVec(sampleDataFile, 0, new int[] {sexInd}, "\t", true, false);

    inclSampMales = Arrays.copyOf(inclSampAll, inclSampAll.length);
    inclSampFemales = Arrays.copyOf(inclSampAll, inclSampAll.length);
    sampleSex = new int[inclSampAll.length];
    for (int i = 0; i < samples.length; i++) {
      int sex = sampleData.getSexForIndividual(samples[i]);
      if (sex == -1) {
        sex = Integer.parseInt(sexData.get(samples[i]).get(0));
      }
      sampleSex[i] = sex;
      if (sex == 1) {
        inclSampFemales[i] = false;
      } else if (sex == 2) {
        inclSampMales[i] = false;
      } else {
        // Leave these for now, but when computing LRRs and BAFs, will need to be crafty....
      }
    }

    markerIndexQueues = new ConcurrentLinkedQueue[threadCount];
    markerLists = new Vector[threadCount];
    fullToTruncMarkerIndices = new Hashtable[threadCount];
    markerDataLoaders = new MarkerDataLoader[threadCount];
    for (int i = 0; i < threadCount; i++) {
      markerLists[i] = new Vector<>();
      markerIndexQueues[i] = new ConcurrentLinkedQueue<>();
      fullToTruncMarkerIndices[i] = new Hashtable<>();
    }

    allMarkers = markerSet.getMarkerNames();
    markerChrs = markerSet.getChrs();
    markerPos = markerSet.getPositions();

    int qInd = 0;
    for (int i = 0; i < markerChrs.length; i++) {
      if (chromPlus11Markers[i]) {
        markerLists[qInd].add(allMarkers[i]);
        fullToTruncMarkerIndices[qInd].put(i, markerLists[qInd].size() - 1);
        markerIndexQueues[qInd].add(Integer.valueOf(i));
        qInd = (qInd + 1) % threadCount;
      }
    }

    for (int i = 0; i < threadCount; i++) {
      markerDataLoaders[i] = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj,
                                                                                     ArrayUtils.toStringArray(markerLists[i]));
    }

    rawCentroidsMale = new float[allMarkers.length][][];
    rawCentroidsFemale = new float[allMarkers.length][][];

    pfbInfo = new Hashtable<>();
    log.report("Computing sex-specific centroids for " + markerCount + " sex-specific markers on "
               + threadCount + " thread(s).");

    PFB pfb = PFB.loadPFB(proj);
    if (pfb == null) {
      String newPFB = PFB.populationBAF(proj);
      pfb = PFB.loadPFB(proj, newPFB);
    }
    final PFB autoPFB = pfb;

    computeHub = Executors.newFixedThreadPool(threadCount);

    // counter for tracking how many centroids have been calculated
    AtomicInteger numCentroidsCalculated = new AtomicInteger(0);

    for (int i = 0; i < threadCount; i++) {
      final int myIndex = i;
      final long myStartTime = System.currentTimeMillis();
      computeHub.execute(new Runnable() {

        @Override
        public void run() {
          int myMarkerCount = 0;
          while (!markerIndexQueues[myIndex].isEmpty()) {
            Integer indexInt = markerIndexQueues[myIndex].poll();
            if (indexInt == null) {
              continue;
            }
            int index = indexInt.intValue();

            if (!chromMarkers[index] && !chromPlus11Markers[index]) {
              continue;
            }
            numCentroidsCalculated.incrementAndGet();
            int markerIndex = fullToTruncMarkerIndices[myIndex].get(index);
            MarkerData markerData = markerDataLoaders[myIndex].requestMarkerData(markerIndex);
            if (!chromMarkers[index]) {
              rawCentroidsMale[index] = new float[][] {{Float.NaN, Float.NaN},
                                                       {Float.NaN, Float.NaN},
                                                       {Float.NaN, Float.NaN}};
              rawCentroidsFemale[index] = new float[][] {{Float.NaN, Float.NaN},
                                                         {Float.NaN, Float.NaN},
                                                         {Float.NaN, Float.NaN}};
              pfbInfo.put(index,
                          new String[][] {{markerData.getMarkerName(),
                                           Integer.toString(markerData.getChr()),
                                           Integer.toString(markerData.getPosition()),
                                           Double.toString(autoPFB.getPfbs()[index])},
                                          {markerData.getMarkerName(),
                                           Integer.toString(markerData.getChr()),
                                           Integer.toString(markerData.getPosition()),
                                           Double.toString(autoPFB.getPfbs()[index])}});
            } else {

              CentroidCompute centCompM = new CentroidCompute(markerData, null, inclSampMales,
                                                              false, // NOT
                                                              // intensity
                                                              // only
                                                              1, // no filtering
                                                              0, // no filtering
                                                              null, // no filtering
                                                              true, // median, not mean
                                                              proj.getLog());

              CentroidCompute centCompF = new CentroidCompute(markerData, null, inclSampFemales,
                                                              false, // NOT
                                                              // intensity
                                                              // only
                                                              1, // no filtering
                                                              0, // no filtering
                                                              null, // no filtering
                                                              true, // median, not mean
                                                              proj.getLog());

              centCompM.computeCentroid(true);
              centCompF.computeCentroid(true);

              rawCentroidsMale[index] = centCompM.getCentroid();
              rawCentroidsFemale[index] = centCompF.getCentroid();

              float[] bafCnt = new float[] {0, 0};
              float[] bafSum = new float[] {0, 0};
              float[] genCnt = new float[] {0, 0};
              float[] bafM = centCompM.getRecomputedBAF();
              byte[] genM = centCompM.getClustGenotypes();
              float[] bafF = centCompF.getRecomputedBAF();
              byte[] genF = centCompF.getClustGenotypes();
              for (int s = 0; s < inclSampAll.length; s++) {
                if (inclSampMales[s]) {
                  if (!Float.isNaN(bafM[s])) {
                    bafSum[0] += bafM[s];
                    bafCnt[0]++;
                    if (genM[s] >= 0) {
                      genCnt[0]++;
                    }
                  }
                }
                if (inclSampFemales[s]) {
                  if (!Float.isNaN(bafF[s])) {
                    bafSum[1] += bafF[s];
                    bafCnt[1]++;
                    if (genF[s] >= 0) {
                      genCnt[1]++;
                    }
                  }
                }
              }

              pfbInfo.put(index, new String[][] {
                                                 {markerData.getMarkerName(),
                                                  Integer.toString(markerData.getChr()),
                                                  Integer.toString(markerData.getPosition()),
                                                  Float.toString(genCnt[0] > 0 ? (bafSum[0]
                                                                                  / bafCnt[0])
                                                                               : 2)},
                                                 {markerData.getMarkerName(),
                                                  Integer.toString(markerData.getChr()),
                                                  Integer.toString(markerData.getPosition()),
                                                  Float.toString(genCnt[1] > 0 ? (bafSum[1]
                                                                                  / bafCnt[1])
                                                                               : 2)}});

            }

            if (markerIndex > 0 && numCentroidsCalculated.get() % 10000 == 0) {
              log.report(ext.getTime() + "\t...sex centroids computed up to marker "
                         + (numCentroidsCalculated.get()) + " of " + markerCount);
            }
            markerDataLoaders[myIndex].releaseIndex(markerIndex);

            myMarkerCount++;
          }

          log.report("Thread " + myIndex + " processed " + myMarkerCount + " markers in "
                     + ext.getTimeElapsed(myStartTime));
        }
      });
    }

    computeHub.shutdown();
    try {
      computeHub.awaitTermination(Long.MAX_VALUE, java.util.concurrent.TimeUnit.NANOSECONDS);
    } catch (InterruptedException e) {
      log.report("Centroid computation was interrupted - .pfb and .cent files may not be complete or correct.");
    }
    computeHub = null;

    int nullCnt = 0;
    for (float[][] element : rawCentroidsFemale) {
      if (element == null) {
        nullCnt++;
      }
    }
    log.report(nullCnt + " null cent entries");

    if (pfbFiles != null) {
      log.report("Writing sex-specific (plus chr11) PFB files");

      try {
        new File(ext.parseDirectoryOfFile(pfbFiles[0])).mkdirs();
        new File(ext.parseDirectoryOfFile(pfbFiles[1])).mkdirs();
        writerM = Files.openAppropriateWriter(pfbFiles[0]);
        writerF = Files.openAppropriateWriter(pfbFiles[1]);

        writerM.println("Name\tChr\tPosition\tPFB");
        writerF.println("Name\tChr\tPosition\tPFB");

        for (int i = 0; i < allMarkers.length; i++) {
          if (!chromPlus11Markers[i]) {
            writerM.println(allMarkers[i] + "\t" + markerChrs[i] + "\t" + markerPos[i] + "\tNaN");
            writerF.println(allMarkers[i] + "\t" + markerChrs[i] + "\t" + markerPos[i] + "\tNaN");
            continue;
          }
          String[][] pfbEntry = pfbInfo.get(Integer.valueOf(i));
          writerM.println(pfbEntry[0][0] + "\t" + pfbEntry[0][1] + "\t" + pfbEntry[0][2] + "\t"
                          + pfbEntry[0][3]);
          writerF.println(pfbEntry[1][0] + "\t" + pfbEntry[1][1] + "\t" + pfbEntry[1][2] + "\t"
                          + pfbEntry[1][3]);
        }

        writerM.flush();
        writerF.flush();

        writerM.close();
        writerF.close();

        proj.PFB_MALE_FILENAME.setValue(pfbFiles[0]);
        proj.PFB_FEMALE_FILENAME.setValue(pfbFiles[1]);
      } catch (IOException e1) {
        log.reportError("Error - problem occured when writing to new sex-specific .pfb files");
        log.reportException(e1);
      }

      writerM = null;
      writerF = null;
    }

    Centroids[] centroids = new Centroids[2];
    centroids[0] = new Centroids(rawCentroidsMale, markerSet.getFingerprint());
    centroids[1] = new Centroids(rawCentroidsFemale, markerSet.getFingerprint());

    if (centFiles != null) {
      log.report("Writing sex-specific Centroid files");

      centroids[0].serialize(centFiles[0]);
      Centroids.exportToText(proj, centFiles[0], centFiles[0] + ".txt", allMarkers);
      proj.SEX_CENTROIDS_MALE_FILENAME.setValue(centFiles[0]);

      centroids[1].serialize(centFiles[1]);
      Centroids.exportToText(proj, centFiles[1], centFiles[1] + ".txt", allMarkers);
      proj.SEX_CENTROIDS_FEMALE_FILENAME.setValue(centFiles[1]);

      proj.saveProperties(new Property[] {proj.SEX_CENTROIDS_MALE_FILENAME,
                                          proj.SEX_CENTROIDS_FEMALE_FILENAME});
    }

    return centroids;
  }

  public static void exportToText(Project proj, String centFilename, String exportFilename,
                                  String[] markerNames) {
    PrintWriter writer;
    Centroids centObject;
    float[][][] centroids;
    String dir;

    dir = proj.PROJECT_DIRECTORY.getValue();
    String file = centFilename.startsWith(dir) || centFilename.contains(":")
                  || centFilename.startsWith("/") ? centFilename : dir + centFilename;
    centObject = Centroids.load(file);
    centroids = centObject.getCentroids();

    if (markerNames.length != centroids.length) {
      System.err.println("Error - mismatched number of markers in the project's marker set and the imported centroids file ("
                         + centFilename + "); aborting");
      return;
    }

    if (MarkerSet.fingerprint(markerNames) != centObject.getFingerprint()) {
      System.err.println("Error - mismatched marker fingerprints in the project's marker set and the imported centroids file ("
                         + centFilename + "); aborting");
      return;
    }

    String outFile = exportFilename.startsWith(dir) || exportFilename.contains(":")
                     || exportFilename.startsWith("/") ? exportFilename : dir + exportFilename;
    try {
      writer = Files.openAppropriateWriter(outFile);
      writer.println("marker_fingerprint=" + centObject.getFingerprint());
      writer.println("MarkerName\tAA_Theta_Mean\tAA_R_Mean\tAB_Theta_Mean\tAB_R_Mean\tBB_Theta_Mean\tBB_R_Mean");
      for (int i = 0; i < markerNames.length; i++) {
        writer.print(markerNames[i]);
        for (int j = 0; j < 3; j++) {
          if (centroids[i] == null || centroids[i][j] == null) {
            writer.print("\t.\t.");
          } else {
            writer.print("\t" + centroids[i][j][0] + "\t" + centroids[i][j][1]);
          }
        }
        writer.println();
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + outFile);
      e.printStackTrace();
    }
  }

  public static void importFromText(Project proj, String importFilename, String centFilename) {

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = null;
    // String centFile = "ForNathan_table.csv";
    // String centFile = "SNP Table Myers raw dataset final 022908.csv";
    // String centFile = "Myers_final_042208_ReclusteredCNV_SNP_Table2.csv";
    // String centFile = "Theta_R_mean_dev_550.csv";
    String centFile = "CentroidExample.csv";
    String intensityFlags = "";
    // String intensityFlags = "intesityFlag.txt";
    // String clusteredCentroids = proj.DEFAULT_ORIGINAL_CENTROIDS_FILENAME;
    // String unclusteredCentroids = proj.DEFAULT_GENOTYPE_CENTROIDS_FILENAME;
    boolean fromGenotypes = false;
    boolean sexSpecific = false;
    boolean projComputeDump = false;
    Project proj;
    String compute = "";
    String importFile = null;
    String exportFile = null;
    int numThreads = 1;

    String usage = "\n" + "cnv.filesys.Centroids requires 0-1 arguments\n"
                   + "   (1) project properties filename (i.e. proj="
                   + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
                   + "   (2) filename (i.e. file=" + centFile + " (default))\n" + " OR\n"
                   + "   (2) generate centroids from genotypes (i.e. -fromGenotypes (not the default))\n"
                   + "   (3) compute and dump centroids for project (i.e. -computeDump (not the default))\n"
                   + " OR\n"
                   + "   (2) file with intensity only flags (i.e. flags=intensityFlags.dat (not the default))\n"
                   + "   (3) centroid file for clustered markers (see \"GENOTYPE_CENTROIDS_FILENAME\" in the Project properties file)\n"
                   + "   (4) centroid file for intensity only markers (see \"GENOTYPE_CENTROIDS_FILENAME\" in the Project properties file)\n"
                   + " OR\n"
                   + "   (2) recompute BAF/LRR and generate new Sample files using these centroids (i.e. compute=genotype.cent (not the default))\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("file=")) {
        centFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("import=")) {
        importFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("export=")) {
        exportFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-fromGenotypes")) {
        fromGenotypes = true;
        numArgs--;
      } else if (arg.startsWith("-computeDump")) {
        projComputeDump = true;
        numArgs--;
      } else if (arg.startsWith("threads=")) {
        numThreads = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("flags=")) {
        intensityFlags = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("compute=")) {
        compute = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-sexSpecific")) {
        sexSpecific = true;
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    proj = new Project(filename);
    // fromGenotypes = true;
    // // compute = "genotype.cent";
    //
    // centFile = "data/genotype.cent";
    // exportFile = "data/genotype.cent.xln";
    try {
      if (exportFile != null) {
        exportToText(proj, centFile, exportFile);
      } else if (sexSpecific) {
        String dir, malePFB, femalePFB, centFilePathM, centFilePathF;
        dir = proj.DATA_DIRECTORY.getValue();
        malePFB = dir + "males.pfb";
        femalePFB = dir + "females.pfb";
        centFilePathM = dir + "sexSpecific_Male.cent";
        centFilePathF = dir + "sexSpecific_Female.cent";
        computeSexSpecificCentroids(proj, new String[] {malePFB, femalePFB},
                                    new String[] {centFilePathM, centFilePathF}, numThreads);
      } else if (importFile != null) {
        importFromText(proj, importFile, centFile);
      } else if (fromGenotypes) {
        parseCentroidsFromGenotypes(proj, ArrayUtils.booleanArray(proj.getSamples().length, true),
                                    1);
      } else if (projComputeDump) {
        CentroidCompute.computeAndDumpCentroids(proj);
      } else if (!compute.equals("")) {
        recompute(proj, compute, false, numThreads);
      } else if (!intensityFlags.equals("")) {
        generateChimeraCentroids(proj, intensityFlags);
      } else {
        parseIlluminaCentroidsFromCSV(proj, centFile);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
