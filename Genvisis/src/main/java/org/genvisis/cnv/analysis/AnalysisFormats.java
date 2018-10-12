// quantile normalize

package org.genvisis.cnv.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;
import java.util.Vector;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import org.genvisis.cnv.filesys.Centroids;
import org.genvisis.cnv.filesys.Centroids.CENTROID_STRATEGY;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.manage.MarkerDataLoader;
import org.genvisis.cnv.prop.Property;
import org.genvisis.cnv.qc.SexChecks;
import org.genvisis.cnv.var.SampleData;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.Matrix;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ProgressMonitor;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Segment;

public class AnalysisFormats implements Runnable {

  public static final String[] PROGRAM_OPTIONS = {"PennCNV", "QuantiSNP"};
  public static final String[] OUTPUT_DIRECTORIES = {"PennCNV/", "QuantiSNP/"};
  public static final int PENN_CNV = 1;
  public static final int QUANTISNP = 2;
  // public static final int EM_ITERATIONS = 10;
  public static final int EM_ITERATIONS = 25;
  private final Project proj;
  private final String[] samples;
  private final int program;
  private final HashSet<String> hash;
  private final int threadCount;

  public AnalysisFormats(Project proj, String[] samples, int program, HashSet<String> hash,
                         int threadCount) {
    this.proj = proj;
    this.samples = samples;
    this.program = program;
    this.hash = hash;
    this.threadCount = threadCount;
  }

  @Override
  public void run() {
    switch (program) {
      case PENN_CNV:
        penncnv(proj, samples, hash, null, threadCount);
        break;
      case QUANTISNP:
        quantisnp(proj, samples, hash);
        break;
      default:
        System.err.println("Error - invalid program option: " + program);
        break;
    }

  }

  public static void penncnv(final Project proj, final String[] samples,
                             final HashSet<String> markersToWrite, String subDir, int threadCount) {
    exportPenncnvSamples(proj, samples, markersToWrite, subDir, threadCount);

    // Create the scripts for building cnvs from the penncnv data
    // TODO: sex checks (3rd flag) - should probably determine if they would be
    // appropriate or not.
    PennCNV.doBatch(proj, true, true, true, false, 1, true, null, null, null, false, true, false,
                    threadCount);
  }

  public static void exportPenncnvSamples(final Project proj, final String[] samples,
                                          final Set<String> markersToWrite, String subDir,
                                          int threadCount) {
    final String[] markerNames = proj.getMarkerNames();
    final boolean gzip;
    final String dir;
    final String sampleDir;
    final Logger log = proj.getLog();

    final String PROG_KEY = "GeneratePennCNV";

    dir = proj.PENNCNV_DATA_DIRECTORY.getValue(false, false) + (subDir == null ? "" : subDir);
    sampleDir = proj.SAMPLE_DIRECTORY.getValue(false, true);
    new File(dir).mkdirs();
    // gzip = proj.getBoolean(proj.PENNCNV_GZIP_YESNO);
    gzip = proj.PENNCNV_GZIP_YESNO.getValue();

    // int threadCount = Runtime.getRuntime().availableProcessors();

    final ConcurrentLinkedQueue<Integer>[] sampleIndexQueues = new ConcurrentLinkedQueue[threadCount];
    for (int i = 0; i < threadCount; i++) {
      sampleIndexQueues[i] = new ConcurrentLinkedQueue<>();
    }
    for (int i = 0; i < samples.length; i++) {
      sampleIndexQueues[i % threadCount].add(i);
    }

    proj.initializeProgressMonitor(null);
    ExecutorService computeHub = Executors.newFixedThreadPool(threadCount);
    for (int threadI = 0; threadI < threadCount; threadI++) {
      final int myIndex = threadI;
      final long myStartTime = System.currentTimeMillis();
      final String MY_PROG_KEY = PROG_KEY + "_thread" + myIndex;

      computeHub.execute(new Runnable() {

        @Override
        public void run() {
          PrintWriter writer;
          int mySampleCount = 0;
          String sampleName;
          float[] lrrs, bafs;
          byte[] genotypes;
          Sample mySample;
          int skippedExports = 0;

          proj.getProgressMonitor()
              .beginDeterminateTask(MY_PROG_KEY,
                                    "Generate PennCNV Files in Thread " + (myIndex + 1),
                                    sampleIndexQueues[myIndex].size(),
                                    ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);

          while (!sampleIndexQueues[myIndex].isEmpty()) {
            int sampleIndex = sampleIndexQueues[myIndex].poll();
            sampleName = samples[sampleIndex];
            String exportFileName = dir + sampleName + (gzip ? ".gz" : "");
            if (!Files.exists(exportFileName)) {
              log.report(ext.getTime() + "\tExporting " + (sampleIndex + 1) + " of "
                         + samples.length + "\t" + sampleName);
              if (Files.exists(sampleDir + sampleName + Sample.SAMPLE_FILE_EXTENSION)) {
                mySample = Sample.loadFromRandomAccessFile(sampleDir + sampleName
                                                           + Sample.SAMPLE_FILE_EXTENSION, false,
                                                           false, true, true, true);
              } else {
                log.reportError("Error - the " + sampleName + Sample.SAMPLE_FILE_EXTENSION
                                + " is not found.");
                proj.getProgressMonitor().endTask(MY_PROG_KEY);
                return;
              }
              lrrs = mySample.getLRRs();
              bafs = mySample.getBAFs();
              genotypes = mySample.getAB_Genotypes();

              try {
                writer = Files.getAppropriateWriter(exportFileName);
                writer.println("Name\t" + sampleName + ".GType\t" + sampleName + ".Log R Ratio\t"
                               + sampleName + ".B Allele Freq");
                for (int j = 0; j < markerNames.length; j++) {
                  if (markersToWrite == null || markersToWrite.contains(markerNames[j])) {
                    writer.println(markerNames[j] + "\t"
                                   + (genotypes[j] == -1 ? "NC" : Sample.AB_PAIRS[genotypes[j]])
                                   + "\t" + lrrs[j] + "\t" + bafs[j]);
                  }
                }
                writer.close();
              } catch (Exception e) {
                log.reportError("Error writing PennCNV data for " + sampleName);
                log.reportException(e);
              }
            } else {
              skippedExports++;
            }

            proj.getProgressMonitor().updateTask(MY_PROG_KEY);
            mySampleCount++;
          }

          proj.getProgressMonitor().endTask(MY_PROG_KEY);
          log.report("Thread " + myIndex + " processed " + mySampleCount + " samples in "
                     + ext.getTimeElapsed(myStartTime)
                     + (skippedExports > 0 ? "; skipped " + skippedExports
                                             + " samples that had been exported previously"
                                           : ""));

        }
      });
    }

    computeHub.shutdown();
    try {
      computeHub.awaitTermination(Long.MAX_VALUE, java.util.concurrent.TimeUnit.NANOSECONDS);
    } catch (InterruptedException e) {
      log.report("Sample export was interrupted - exported sample files may not be complete or correct.");
    }
  }

  @SuppressWarnings("unchecked")
  public static String[] pennCNVSexHackMultiThreaded(Project proj, String gcModelFile,
                                                     CENTROID_STRATEGY chr11Strategy,
                                                     boolean useExcluded, int threadCount) {
    String sampleDataFile;
    final String sampleDir;
    String sexDir;
    final String maleDir;
    final String femaleDir;
    String malePFBFile, femalePFBFile, newGCFile, centFilePathM, centFilePathF;
    final String[] allMarkers;
    final String[] allSamples;
    String[] header;
    final SampleData sampleData;
    MarkerSetInfo ms;
    byte[] markerChrs;
    final boolean gzip;
    final boolean[] centroidsMarkersList;
    final boolean[] includeMarkersList;
    boolean[] includeSamplesList;
    final Hashtable<String, Vector<String>> sexData;
    final HashSet<String> chr11Markers = new HashSet<>();
    Centroids[] centroids;

    String centFile = proj.CUSTOM_CENTROIDS_FILENAME.getValue();

    boolean computeCentroids = chr11Strategy == CENTROID_STRATEGY.COMPUTE_CENT;
    Centroids tempCentroids = null;
    if (chr11Strategy == CENTROID_STRATEGY.USE_CENT_IF_EXISTS_OTHERWISE_ORIG
        || chr11Strategy == CENTROID_STRATEGY.USE_CENT_IF_EXISTS_OTHERWISE_COMPUTE) {
      if (Files.exists(centFile)) {
        tempCentroids = Centroids.load(centFile);
      } else if (chr11Strategy == CENTROID_STRATEGY.USE_CENT_IF_EXISTS_OTHERWISE_COMPUTE) {
        computeCentroids = true;
      }
    }
    if (computeCentroids) {
      tempCentroids = CentroidCompute.computeAndDumpCentroids(proj);
    }
    final Centroids autoCentroids = tempCentroids;

    final Logger log = proj.getLog();

    sexDir = proj.PENNCNV_DATA_DIRECTORY.getValue() + "sexSpecific/";

    maleDir = sexDir + "male/";
    femaleDir = sexDir + "female/";

    new File(sexDir).mkdirs();
    new File(maleDir).mkdir();
    new File(femaleDir).mkdir();

    malePFBFile = sexDir + "males.pfb";
    femalePFBFile = sexDir + "females.pfb";
    newGCFile = sexDir + "sexSpecific.gcModel";

    centFilePathM = sexDir + "sexSpecific_Male.cent";
    centFilePathF = sexDir + "sexSpecific_Female.cent";

    ms = proj.getMarkerSet();
    sampleData = proj.getSampleData(false);
    { // block is similar to getChromosomalMarkersOnly, but also puts markers into other
     // datastructures, hence not using the utility method;
      allMarkers = ms.getMarkerNames();
      markerChrs = ms.getChrs();
      centroidsMarkersList = new boolean[allMarkers.length];
      includeMarkersList = new boolean[allMarkers.length];

      for (int i = 0; i < markerChrs.length; i++) {
        switch (markerChrs[i]) {
          case 23:
          case 24:
          case 25:
          case 26:
            centroidsMarkersList[i] = true;
            includeMarkersList[i] = true;
            break;
          case 11:
            includeMarkersList[i] = true;
            chr11Markers.add(allMarkers[i]);
            break;
          default:
            includeMarkersList[i] = false;
            break;
        }
      }
    }

    final int[] markerIndicesToUse = ArrayUtils.booleanArrayToIndices(includeMarkersList);

    if (Files.exists(centFilePathM) && Files.exists(centFilePathF)) {
      centroids = new Centroids[] {Centroids.load(centFilePathM), Centroids.load(centFilePathM)};
    } else {
      centroids = Centroids.computeSexSpecificCentroids(proj,
                                                        new String[] {malePFBFile, femalePFBFile},
                                                        new String[] {centFilePathM, centFilePathF},
                                                        threadCount);
    }
    final float[][][] rawCentroidsMale;
    final float[][][] rawCentroidsFemale;
    rawCentroidsMale = centroids[0].getCentroids();
    rawCentroidsFemale = centroids[1].getCentroids();

    log.report("Exporting sex-specific sample data");

    sampleDir = proj.SAMPLE_DIRECTORY.getValue(false, true);
    gzip = proj.PENNCNV_GZIP_YESNO.getValue();

    if (!useExcluded) {
      includeSamplesList = proj.getSamplesToInclude(null);
      if (!sampleData.hasExcludedIndividuals()) {
        log.report("Warning - there is no 'Exclude' column in SampleData.txt; centroids will be determined using all samples.");
      }
      allSamples = ArrayUtils.subArray(proj.getSamples(), includeSamplesList);
    } else {
      allSamples = proj.getSamples();
    }

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

    if (threadCount == -1) {
      threadCount = Runtime.getRuntime().availableProcessors();
    }

    final ConcurrentLinkedQueue<Integer>[] sampleIndexQueues = new ConcurrentLinkedQueue[threadCount];
    for (int i = 0; i < threadCount; i++) {
      sampleIndexQueues[i] = new ConcurrentLinkedQueue<>();
    }
    for (int i = 0; i < allSamples.length; i++) {
      sampleIndexQueues[i % threadCount].add(i);
    }

    log.report("Exporting " + markerIndicesToUse.length + " markers per sample.");

    ExecutorService computeHub = Executors.newFixedThreadPool(threadCount);
    for (int threadI = 0; threadI < threadCount; threadI++) {
      final int myIndex = threadI;
      final long myStartTime = System.currentTimeMillis();
      computeHub.execute(new Runnable() {

        @Override
        public void run() {
          PrintWriter writer;
          int mySampleCount = 0;
          String sampleName;
          float[] thetas, rs, lrrs, bafs;
          float[][][] cent;
          byte[] genotypes;
          int skippedExports = 0;

          while (!sampleIndexQueues[myIndex].isEmpty()) {
            Sample mySample;
            int sampleIndex = sampleIndexQueues[myIndex].poll();
            sampleName = allSamples[sampleIndex];
            int sex = sampleData.getSexForIndividual(sampleName);
            if (sex == -1) {
              sex = Integer.parseInt(sexData.get(sampleName.toUpperCase()).get(0));
            }
            boolean compFemale = SexChecks.EstimatedSex.values()[sex].getKaryotype().contains("XX");
            String exportFileName = (compFemale ? femaleDir : maleDir) + sampleName
                                    + (gzip ? ".gz" : "");
            if (!Files.exists(exportFileName)) {
              log.report(ext.getTime() + "\tExporting " + (sampleIndex + 1) + " of "
                         + allSamples.length + ".");
              if (Files.exists(sampleDir + sampleName + Sample.SAMPLE_FILE_EXTENSION)) {
                mySample = Sample.loadFromRandomAccessFile(sampleDir + sampleName
                                                           + Sample.SAMPLE_FILE_EXTENSION, false,
                                                           true, false, false, true);
              } else {
                log.reportError("Error - the " + sampleName + Sample.SAMPLE_FILE_EXTENSION
                                + " is not found.");
                // TODO okay to just skip this sample instead of halting entirely?
                continue;
              }

              thetas = mySample.getThetas();
              rs = mySample.getRs();
              genotypes = mySample.getAB_Genotypes();
              cent = compFemale ? rawCentroidsFemale : rawCentroidsMale;
              lrrs = autoCentroids == null ? mySample.getLRRs()
                                           : mySample.getLRRs(autoCentroids.getCentroids());
              bafs = autoCentroids == null ? mySample.getBAFs()
                                           : mySample.getBAFs(autoCentroids.getCentroids());

              if (lrrs == null) {
                if (autoCentroids != null) {
                  log.reportError("LRRs for sample " + sampleName
                                  + " were null after centroid transformation.");
                }
                if (mySample.getXs() == null) {
                  log.reportError("Xs for sample " + sampleName + " were null.");
                }
              }

              try {
                writer = Files.getAppropriateWriter(exportFileName);
                writer.println("Name\t" + sampleName + ".GType\t" + sampleName + ".Log R Ratio\t"
                               + sampleName + ".B Allele Freq");

                int skip = 0;
                for (int m = 0; m < markerIndicesToUse.length; m++) {
                  int j = markerIndicesToUse[m];
                  if (null == cent[j]) {
                    skip++;
                    continue;
                  }

                  float lrr;
                  float baf;

                  if (chr11Markers.contains(allMarkers[j])) {
                    lrr = lrrs[j];
                    baf = bafs[j];
                  } else {
                    lrr = Centroids.calcLRR(thetas[j], rs[j], cent[j]);
                    baf = Centroids.calcBAF(thetas[j], cent[j]);
                  }

                  writer.println(allMarkers[j] + "\t"
                                 + (genotypes[j] == -1 ? "NC" : Sample.AB_PAIRS[genotypes[j]])
                                 + "\t" + lrr + "\t" + baf);
                }
                if (skip > 0) {
                  System.out.println("Skipped " + skip + " markers due to null centroid entries.");
                }
                writer.flush();
                writer.close();
              } catch (Exception e) {
                log.reportError("Error writing sex-specific (" + (compFemale ? "female" : "male")
                                + ") PennCNV data for " + sampleName);
                log.reportException(e);
              }
            } else {
              skippedExports++;
            }
            mySampleCount++;
          }

          log.report("Thread " + myIndex + " processed " + mySampleCount + " samples in "
                     + ext.getTimeElapsed(myStartTime)
                     + (skippedExports > 0 ? "; skipped " + skippedExports
                                             + " samples that had been exported previously"
                                           : ""));
        }
      });
    }

    computeHub.shutdown();
    try {
      computeHub.awaitTermination(Long.MAX_VALUE, java.util.concurrent.TimeUnit.NANOSECONDS);
    } catch (InterruptedException e) {
      log.report("Sample export was interrupted - exported sample files may not be complete or correct.");
    }
    computeHub = null;

    filterSexSpecificGCModel(proj, gcModelFile, newGCFile);

    return new String[] {malePFBFile, femalePFBFile, newGCFile};
  }

  public static String[] pennCNVSexHackSingleThreaded(Project proj, String gcModelFile) {
    // exports data for chr23-chr26 and recodes them as chr1-chr4 in a new subdirectory
    // ~/penndata/sexSpecific/
    boolean gzip, writeNewPFBs, writeCentroids, writeGCFile;
    boolean[] inclSampAll, inclSampMales, inclSampFemales;
    int[] sampleSex;
    byte[] markerChrs, genM, genF, genotypes;
    float[] bafCnt, bafSum, genCnt, bafM, bafF, thetas, rs;
    String[] allMarkers, sexMarkers, samples, header, centFilePathM, centFilePathF;
    Logger log;
    MarkerSetInfo ms;
    MarkerDataLoader markerDataLoader;
    SampleData sampleData;
    PrintWriter writer;
    Sample samp;
    String sampleDataFile, sampleDir, sexDir, pennDir, pennData, maleDir, femaleDir, malePFBFile,
        femalePFBFile, newGCFile;
    float[][][] rawCentroidsMale, rawCentroidsFemale;

    Vector<String[]> malePFBs, femalePFBs;
    Vector<String> markerList;
    Hashtable<String, Vector<String>> sexData;
    Hashtable<String, Integer> sexMarkerToIndex = new Hashtable<>();

    Centroids centroidsMale, centroidsFemale;

    log = proj.getLog();

    pennDir = proj.getProperty(proj.PENNCNV_RESULTS_DIRECTORY);
    pennData = proj.getProperty(proj.PENNCNV_DATA_DIRECTORY);
    sexDir = proj.PROJECT_DIRECTORY.getValue() + pennDir + pennData + "sexSpecific/";

    maleDir = sexDir + "male/";
    femaleDir = sexDir + "female/";

    new File(sexDir).mkdirs();
    new File(maleDir).mkdir();
    new File(femaleDir).mkdir();

    malePFBFile = sexDir + "males.pfb";
    femalePFBFile = sexDir + "females.pfb";
    newGCFile = sexDir + "sexSpecific.gcModel";

    // centFilePathM = sexDir + "sexSpecific_Male.cent";
    // centFilePathF = sexDir + "sexSpecific_Female.cent";

    centFilePathM = new String[] {pennDir + pennData + "sexSpecific/sexSpecific_Male.cent", ""};
    centFilePathF = new String[] {pennDir + pennData + "sexSpecific/sexSpecific_Female.cent", ""};
    centFilePathM[1] = proj.PROJECT_DIRECTORY.getValue() + centFilePathM[0];
    centFilePathF[1] = proj.PROJECT_DIRECTORY.getValue() + centFilePathF[0];

    writeCentroids = !Files.exists(centFilePathM[1]) || !Files.exists(centFilePathF[1]);
    writeNewPFBs = !Files.exists(malePFBFile) && !Files.exists(femalePFBFile);
    writeGCFile = !Files.exists(newGCFile);

    if (!writeCentroids && !writeNewPFBs && !writeGCFile) {
      return new String[] {malePFBFile, femalePFBFile, newGCFile};
    }

    ms = proj.getMarkerSet();
    sampleData = proj.getSampleData(false);

    allMarkers = ms.getMarkerNames();
    markerChrs = ms.getChrs();
    markerList = new Vector<>();

    for (int i = 0; i < markerChrs.length; i++) {
      switch (markerChrs[i]) {
        case 23:
        case 24:
        case 25:
        case 26:
          markerList.add(allMarkers[i]);
          sexMarkerToIndex.put(allMarkers[i], Integer.valueOf(i));
          break;
      }
    }
    sexMarkers = ArrayUtils.toStringArray(markerList);

    markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, sexMarkers);

    inclSampAll = proj.getSamplesToInclude(null);
    if (!sampleData.hasExcludedIndividuals()) {
      log.report("Warning - there is no 'Exclude' column in SampleData.txt; centroids will be determined using all samples.");
    }
    samples = proj.getSamples();// Array.subArray(proj.getSamples(), inclSampAll);
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

    inclSampMales = new boolean[inclSampAll.length];
    inclSampFemales = new boolean[inclSampAll.length];
    sampleSex = new int[inclSampAll.length];
    for (int i = 0; i < inclSampAll.length; i++) {
      int sex = sampleData.getSexForIndividual(samples[i]);
      if (sex == -1) {
        sex = Integer.parseInt(sexData.get(samples[i].toUpperCase()).get(0));
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

    // TODO should we write these out as they're computed instead of storing and writing later?
    malePFBs = new Vector<>();
    femalePFBs = new Vector<>();

    rawCentroidsMale = new float[sexMarkers.length][][];
    rawCentroidsFemale = new float[sexMarkers.length][][];

    log.report("Computing sex-specific centroids for " + sexMarkers.length
               + " sex-specific markers on one thread.");
    CentroidCompute centCompM;
    CentroidCompute centCompF;
    for (int i = 0; i < sexMarkers.length; i++) {
      MarkerData markerData = markerDataLoader.requestMarkerData(i);

      centCompM = new CentroidCompute(markerData, null, inclSampMales, false, // NOT intensity only
                                      1, // no filtering
                                      0, // no filtering
                                      null, // no filtering
                                      true, // median, not mean
                                      proj.getLog());

      centCompF = new CentroidCompute(markerData, null, inclSampFemales, false, // NOT intensity
                                      // only
                                      1, // no filtering
                                      0, // no filtering
                                      null, // no filtering
                                      true, // median, not mean
                                      proj.getLog());

      centCompM.computeCentroid(true);
      centCompF.computeCentroid(true);

      rawCentroidsMale[i] = centCompM.getCentroid();
      rawCentroidsFemale[i] = centCompF.getCentroid();

      bafCnt = new float[] {0, 0};
      bafSum = new float[] {0, 0};
      genCnt = new float[] {0, 0};
      bafM = centCompM.getRecomputedBAF();
      genM = centCompM.getClustGenotypes();
      bafF = centCompF.getRecomputedBAF();
      genF = centCompF.getClustGenotypes();
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

      malePFBs.add(new String[] {markerData.getMarkerName(), "" + (markerData.getChr() - 22),
                                 "" + markerData.getPosition(),
                                 "" + (genCnt[0] > 0 ? (bafSum[0] / bafCnt[0]) : 2)});
      femalePFBs.add(new String[] {markerData.getMarkerName(), "" + (markerData.getChr() - 22),
                                   "" + markerData.getPosition(),
                                   "" + (genCnt[1] > 0 ? (bafSum[1] / bafCnt[1]) : 2)});
      if (i > 0 && i % 10000 == 0) {
        log.report(ext.getTime() + "\t...sex centroids computed up to marker " + i + " of "
                   + sexMarkers.length);
      }

      markerDataLoader.releaseIndex(i);
      centCompM = null;
      centCompF = null;
    }

    if (writeNewPFBs) {
      log.report("Writing sex-specific PFB files");

      try {
        writer = Files.openAppropriateWriter(malePFBFile);
        writer.println("Name\tChr\tPosition\tPFB");
        for (String[] male : malePFBs) {
          writer.println(male[0] + "\t" + male[1] + "\t" + male[2] + "\t" + male[3]);
        }
        writer.close();
      } catch (IOException e1) {
        log.reportError("Error - problem occured when writing to new male-only .pfb file");
        log.reportException(e1);
      }

      try {
        writer = Files.openAppropriateWriter(femalePFBFile);
        writer.println("Name\tChr\tPosition\tPFB");
        for (String[] female : femalePFBs) {
          writer.println(female[0] + "\t" + female[1] + "\t" + female[2] + "\t" + female[3]);
        }
        writer.close();
      } catch (IOException e1) {
        log.reportError("Error - problem occured when writing to new female-only .pfb file");
        log.reportException(e1);
      }
    }
    malePFBs = null;
    femalePFBs = null;

    if (writeCentroids) {
      log.report("Writing sex-specific Centroid files");

      centroidsMale = new Centroids(rawCentroidsMale, MarkerSet.fingerprint(sexMarkers));
      centroidsMale.serialize(centFilePathM[1]);
      Centroids.exportToText(proj, centFilePathM[0], centFilePathM[0] + ".txt", sexMarkers);
      proj.SEX_CENTROIDS_MALE_FILENAME.setValue(centFilePathM[1]);

      centroidsFemale = new Centroids(rawCentroidsFemale, MarkerSet.fingerprint(sexMarkers));
      centroidsFemale.serialize(centFilePathF[1]);
      Centroids.exportToText(proj, centFilePathF[0], centFilePathF[0] + ".txt", sexMarkers);
      proj.SEX_CENTROIDS_FEMALE_FILENAME.setValue(centFilePathF[1]);

      proj.saveProperties(new Property[] {proj.SEX_CENTROIDS_MALE_FILENAME,
                                          proj.SEX_CENTROIDS_FEMALE_FILENAME});
    }
    centroidsMale = null;
    centroidsFemale = null;

    log.report("Exporting sex-specific sample data");

    sampleDir = proj.SAMPLE_DIRECTORY.getValue(false, true);
    // gzip = proj.getBoolean(proj.PENNCNV_GZIP_YESNO);
    gzip = proj.PENNCNV_GZIP_YESNO.getValue();

    int skippedExports = 0;
    for (int i = 0; i < samples.length; i++) {
      int sex = sampleData.getSexForIndividual(samples[i]);
      if (sex == -1) {
        sex = Integer.parseInt(sexData.get(samples[i].toUpperCase()).get(0));
      }
      boolean compFemale = SexChecks.EstimatedSex.values()[sex].getKaryotype().contains("XX");

      String exportFileName = (compFemale ? femaleDir : maleDir) + samples[i] + (gzip ? ".gz" : "");
      if (!Files.exists(exportFileName)) {
        log.report(ext.getTime() + "\tTransforming " + (i + 1) + " of " + samples.length);
        if (Files.exists(sampleDir + samples[i] + Sample.SAMPLE_FILE_EXTENSION)) {
          samp = Sample.loadFromRandomAccessFile(sampleDir + samples[i]
                                                 + Sample.SAMPLE_FILE_EXTENSION, false, true, false,
                                                 false, true);
        } else {
          log.reportError("Error - the " + samples[i] + Sample.SAMPLE_FILE_EXTENSION
                          + " is not found.");
          // TODO okay to just skip this sample instead of halting entirely?
          continue;
        }

        thetas = samp.getThetas();
        rs = samp.getRs();
        genotypes = samp.getAB_Genotypes();

        try {
          writer = Files.getAppropriateWriter(exportFileName);
          writer.println("Name\t" + samples[i] + ".GType\t" + samples[i] + ".Log R Ratio\t"
                         + samples[i] + ".B Allele Freq");
          for (int j = 0; j < sexMarkers.length; j++) {
            int markerIndex = sexMarkerToIndex.get(sexMarkers[j]).intValue();

            float lrr = Centroids.calcLRR(thetas[markerIndex], rs[markerIndex],
                                          (compFemale ? rawCentroidsFemale[j]
                                                      : rawCentroidsMale[j]));
            float baf = Centroids.calcBAF(thetas[markerIndex], (compFemale ? rawCentroidsFemale[j]
                                                                           : rawCentroidsMale[j]));

            writer.println(sexMarkers[j] + "\t"
                           + (genotypes[markerIndex] == -1 ? "NC"
                                                           : Sample.AB_PAIRS[genotypes[markerIndex]])
                           + "\t" + lrr + "\t" + baf);
          }
          writer.close();
        } catch (Exception e) {
          log.reportError("Error writing sex-specific (" + (compFemale ? "female" : "male")
                          + ") PennCNV data for " + samples[i]);
          log.reportException(e);
        }
      } else {
        skippedExports++;
      }
    }

    log.report(skippedExports > 0 ? "Skipped " + skippedExports + " of " + samples.length
                                    + " samples that had been exported previously"
                                  : "");

    if (writeGCFile) {
      filterSexSpecificGCModel(proj, gcModelFile, newGCFile);
    }

    return new String[] {malePFBFile, femalePFBFile, newGCFile};
  }

  public static String filterSexSpecificGCModel(Project proj, String gcModelFile,
                                                String newGCFile) {
    // TODO combine method with filter methods in PennCNV - only difference is changing chr #
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] chrs = new String[] {"11", "23", "X", "24", "Y", "25", "XY", "26", "M"};

    try {
      (new File(ext.parseDirectoryOfFile(newGCFile))).mkdirs();
      reader = new BufferedReader(new FileReader(gcModelFile));
      writer = Files.openAppropriateWriter(newGCFile);

      String temp;
      String[] line;
      writer.println(reader.readLine());
      gc: while ((temp = reader.readLine()) != null) {
        line = temp.trim().split(PSF.Regex.GREEDY_WHITESPACE);
        for (String chr : chrs) {
          if (line[1].equals(chr)) {
            byte chrVal = 0;
            if ("11".equals(chr)) {
              chrVal = 11;
            } else if ("23".equals(chr) || "X".equals(chr)) {
              chrVal = 23;
            } else if ("24".equals(chr) || "Y".equals(chr)) {
              chrVal = 24;
            } else if ("25".equals(chr) || "XY".equals(chr)) {
              chrVal = 25;
            } else if ("26".equals(chr) || "M".equals(chr)) {
              chrVal = 26;
            }

            writer.println(line[0] + "\t" + (chrVal == 11 ? chrVal : (chrVal - 22)) + "\t" + line[2]
                           + "\t" + line[3]);
            continue gc;
          }
        }
      }
    } catch (IOException e) {
      proj.getLog().reportError("Error - transforming sex-specific gcModel failed");
      proj.getLog().reportException(e);
      return gcModelFile;
    } finally {
      if (reader != null) {
        try {
          reader.close();
        } catch (IOException e) {
          proj.getLog()
              .reportError("Error - couldn't properly close file reader for " + gcModelFile);
          proj.getLog().reportException(e);
        }
        reader = null;
      }
      if (writer != null) {
        writer.close();
        writer = null;
      }
    }

    return newGCFile;
  }

  public static void quantisnp(Project proj, String[] samples, HashSet<String> hash) {
    PrintWriter writer;
    MarkerSetInfo set;
    String[] markerNames = proj.getMarkerNames();
    Sample samp;
    float[] lrrs, bafs;
    byte[] chrs;
    int[] positions;

    set = proj.getMarkerSet();
    chrs = set.getChrs();
    positions = set.getPositions();

    new File(proj.PROJECT_DIRECTORY.getValue() + "quanti_data/").mkdirs();
    for (int i = 0; i < samples.length; i++) {
      System.out.println(ext.getTime() + "\tTransforming " + (i + 1) + " of " + samples.length);
      samp = proj.getPartialSampleFromRandomAccessFile(samples[i]);
      set.checkFingerprint(samp);
      lrrs = samp.getLRRs();
      bafs = samp.getBAFs();

      try {
        writer = Files.openAppropriateWriter(proj.PROJECT_DIRECTORY.getValue() + "quanti_data/"
                                             + samples[i]);
        writer.println("Name\tChr\tPosition\t" + samples[i] + ".Log R Ratio\t" + samples[i]
                       + ".B Allele Freq");
        for (int j = 0; j < markerNames.length; j++) {
          if (hash == null || hash.contains(markerNames[j])) {
            writer.println(markerNames[j] + "\t" + chrs[j] + "\t" + positions[j] + "\t" + lrrs[j]
                           + "\t" + bafs[j]);
          }
        }
        writer.close();
      } catch (Exception e) {
        System.err.println("Error writing QuantiSNP data for " + samples[i]);
        e.printStackTrace();
      }
    }
  }

  public static void batchQuantiSNP(Project proj, int numBatches) {
    Vector<String[]> v = new Vector<>();
    Hashtable<String, String> genders;
    String[] inputs, outputs;
    String commands, gender;

    // genders = HashVec.loadFileToHashString(proj.getFilename(proj.SAMPLE_DATA_FILENAME), "DNA",
    // new String[] {"CLASS=Gender"}, "");
    genders = HashVec.loadFileToHashString(proj.SAMPLE_DATA_FILENAME.getValue(), "DNA",
                                           new String[] {"CLASS=Gender"}, "");

    inputs = new File(proj.PROJECT_DIRECTORY.getValue()
                      + "quanti_data/").list(new FilenameFilter() {

                        @Override
                        public boolean accept(File file, String filename) {
                          return filename.endsWith(".qs");
                        }
                      });

    if (inputs == null) {
      System.err.println("Error - QuantiSNP inputs files have not yet been created");
      return;
    }

    outputs = new File(proj.RESULTS_DIRECTORY.getValue(false, true)
                       + "QuantiSNP/").list(new FilenameFilter() {

                         @Override
                         public boolean accept(File file, String filename) {
                           return filename.endsWith("_output.out");
                         }
                       });

    if (outputs == null) {
      System.out.println("Found " + inputs.length
                         + " samples; creating output directory for QuantiSNP in "
                         + proj.RESULTS_DIRECTORY.getValue(false, true) + "QuantiSNP/");
      outputs = new String[0];
    } else {
      System.out.println("Found " + inputs.length + " samples, as well as results for "
                         + outputs.length + " that have been done (not necessarily the same ones)");
    }

    for (String input : inputs) {
      if (ext.indexOfStr(ext.rootOf(input) + "_output.out", outputs) == -1) {
        if (genders.containsKey(ext.rootOf(input))) {
          gender = genders.get(ext.rootOf(input));
          if (gender.equals("M") || gender.equals("1")) {
            gender = "male";
          } else if (gender.equals("F") || gender.equals("2")) {
            gender = "female";
          } else {
            System.err.println("Error - '" + gender
                               + "' is not a valid gender (expecting M/F or 1/2)");
          }
        } else {
          System.err.println("Error - no gender found for subject '" + ext.rootOf(input) + "'");
          gender = null;
        }

        v.add(new String[] {ext.rootOf(input), gender});
      }
    }

    System.out.println("Made " + numBatches + " batch files that will take care of the " + v.size()
                       + " files yet to parse");

    // commands = "quantisnp.exe --config ../windows/config.dat --emiters "+EM_ITERATIONS+"
    // --Lsetting 2000000 --maxcopy 3 --printRS --doGCcorrect --gcdir ../gc/b36/ --output
    // "+OUTPUT_DIRECTORIES[1]+"[%0].out --gender [%1]--input-files ../source/[%0].qs 300\n\n";
    commands = "quantisnp --output " + proj.RESULTS_DIRECTORY.getValue(false, true)
               + OUTPUT_DIRECTORIES[1]
               + "[%0].out --gender [%1] --input-files ../source/[%0].qs 300\n\n";
    Files.batchIt("batch", null, numBatches, commands, Matrix.toStringArrays(v));
  }

  // TODO convert this to an Executor
  public static void launch(Project proj, int program, String markers, int numThreads) {
    Vector<Vector<String>> sampleLists = new Vector<>();
    String[] samples = proj.getSamples();
    Thread[] threads;
    HashSet<String> hash;

    for (int i = 0; i < numThreads; i++) {
      sampleLists.add(new Vector<String>());
    }
    for (int i = 0; i < samples.length; i++) {
      sampleLists.elementAt(i % numThreads).add(samples[i]);
    }
    if (markers == null) {
      hash = null;
    } else {
      hash = HashVec.loadFileToHashSet(proj.PROJECT_DIRECTORY.getValue() + markers, false);
    }
    threads = new Thread[numThreads];
    for (int i = 0; i < numThreads; i++) {
      threads[i] = new Thread(new AnalysisFormats(proj,
                                                  ArrayUtils.toStringArray(sampleLists.elementAt(i)),
                                                  program, hash, 1));
      threads[i].start();
      try {
        Thread.sleep(100L);
      } catch (InterruptedException ex) {}
    }
  }

  public static void filter(Project proj, String regions, String list, String outfile) {
    PrintWriter writer;
    MarkerSetInfo markers;
    Segment[] segs;
    String[] markerNames;
    byte[] chrs;
    int[] positions;
    int countFromList, countInRegions, countOverlap;
    HashSet<String> hash;

    if (outfile == null) {
      System.err.println("Error - outfile is defined as null; need to provide a filename before results can be filtered");
      return;
    }

    if (regions.equals("")) {
      segs = new Segment[0];
    } else {
      segs = Segment.loadUCSCregions(proj.PROJECT_DIRECTORY.getValue() + regions, false);
    }

    if (list.equals("")) {
      hash = new HashSet<>();
    } else {
      hash = HashVec.loadFileToHashSet(proj.PROJECT_DIRECTORY.getValue() + list, false);
    }

    markers = proj.getMarkerSet();
    markerNames = markers.getMarkerNames();
    chrs = markers.getChrs();
    positions = markers.getPositions();

    try {
      writer = Files.openAppropriateWriter(proj.PROJECT_DIRECTORY.getValue() + outfile);
      countFromList = countInRegions = countOverlap = 0;
      for (int i = 0; i < markerNames.length; i++) {
        if (segs.length > 0
            && Segment.overlapsAny(new Segment(chrs[i], positions[i], positions[i]), segs)) {
          countInRegions++;
          if (hash.contains(markerNames[i])) {
            countFromList++;
            countOverlap++;
          }
        } else if (hash.contains(markerNames[i])) {
          countFromList++;
        } else {
          writer.println(markerNames[i]);
        }
      }
      System.out.println("Started off with " + chrs.length + " markers in the dataset");
      System.out.println("   " + countFromList + " of " + hash.size()
                         + " markers on the list were removed ("
                         + ext.formDeci((countFromList - countOverlap) / (double) chrs.length * 100,
                                        2, true)
                         + "% of total)");
      System.out.println("   " + countInRegions
                         + " were found within the list of regions and removed ("
                         + ext.formDeci((countInRegions - countOverlap) / (double) chrs.length
                                        * 100, 2, true)
                         + "% of total)");
      System.out.println("   " + countOverlap + " overlap in filtering criteria ("
                         + ext.formDeci(countOverlap / (double) chrs.length * 100, 2, true)
                         + "% of total)");
      System.out.println("Leaving behind "
                         + (chrs.length - countFromList - countInRegions + countOverlap)
                         + " in final marker list ("
                         + ext.formDeci((chrs.length - countFromList - countInRegions
                                         + countOverlap)
                                        / (double) chrs.length * 100, 2, true)
                         + "% of total)");
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + outfile);
      e.printStackTrace();
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = null;
    int numThreads = 6;
    int program = PENN_CNV;
    String filterRegions = "";
    String filterList = "";
    String markers = null;
    String gcmodel = null;
    Project proj;

    String usage = "\n" + "filesys.AnalysisFormats requires 0-1 arguments\n"
                   + "   (1) project properties filename (i.e. proj="
                   + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
                   + "   (2) number of threads to use (i.e. threads=" + numThreads + " (default))\n"
                   + "   (3) filter markers out within specified regions (i.e. filterRegions=problematicRegions.dat (not the default))\n"
                   + "   (4) filter markers out from list (i.e. filterList=drops.dat (not the default))\n"
                   + "   (5) input/output file of final list of markers to use (all markers if null) (i.e. markers="
                   + markers + " (default))\n" + "   (6) program option (i.e. program=" + program
                   + " (default))\n" + " OR \n" + "   (1) Project properties file (i.e. proj= )\n"
                   + "   (2) GCMODEL File (i.e. gcmodel= )\n" + "" + "";
    for (int i = 0; i < PROGRAM_OPTIONS.length; i++) {
      usage += "           " + (i + 1) + " = " + PROGRAM_OPTIONS[i] + "\n";
    }

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("threads=")) {
        numThreads = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("program=")) {
        program = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("filterRegions=")) {
        filterRegions = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("filterList=")) {
        filterList = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("markers=")) {
        markers = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("gcmodel=")) {
        gcmodel = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    // filterRegions = "data/problematicRegions.dat";
    // filterList = "data/drops.dat";
    // markers = "finalMarkerList.dat";
    try {
      proj = new Project(filename);
      if (gcmodel != null) {
        String pennData = proj.getProperty(proj.PENNCNV_DATA_DIRECTORY);
        String sexDir = pennData + "sexSpecific/";
        String newGCFile = sexDir + "sexSpecific.gcModel";
        filterSexSpecificGCModel(proj, gcmodel, newGCFile);
      } else if (!filterRegions.equals("") || !filterList.equals("")) {
        filter(proj, filterRegions, filterList, markers);
      } else {
        launch(proj, program, markers, numThreads);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
