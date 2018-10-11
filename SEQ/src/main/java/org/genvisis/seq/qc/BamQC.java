package org.genvisis.seq.qc;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import org.genvisis.seq.analysis.GATK_Genotyper.ANNOTATION_BUILD;
import org.genvisis.seq.analysis.SNPEFF;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.SerializedFiles;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Positions;
import org.pankratzlab.common.filesys.Segment;
import org.pankratzlab.common.stats.Histogram;
import org.pankratzlab.common.stats.Histogram.DynamicHistogram;
import com.google.common.primitives.Doubles;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class BamQC {

  public static final String[] QC_HEADER = {"Input File", "numTotal", "numUnique", "numDuplicated",
                                            "numUnMapped", "numOnTarget", "NumSecondaryAlignments",
                                            "NumInvalidAlignments", "PercentDuplicated",
                                            "PercentUnMapped", "PercentOnTarget",
                                            "AverageInsertSize", "AverageOnTargetInsertSize",
                                            "OnTargetInsertSizeStdev", "MappingQualityFilter",
                                            "PhreadScoreFilter", "Total Base Pairs Targeted"};
  public static final String[] HIST_HEADER = {"Bin"};
  public static int NUM_GC_HISTOGRAMS = 2;
  public static int NO_QC_GC_HISTOGRAM = 0;
  public static int QC_GC_HISTOGRAM = 1;
  public static int NUM_INSERT_HISTOGRAMS = 2;
  public static int NO_QC_INSERT_HISTOGRAM = 0;
  public static int QC_INSERT_HISTOGRAM = 1;
  public static final String READ_DEPTH_AT = "Percent Coverage At";
  // private static final int DEFUALT_GC_SIG_FIGS = 2;
  // private static final double DEFUALT_MIN_GC = 0.0D;
  // private static final double DEFUALT_MAX_GC = 1.0D;
  // private static final double NORM_UNITS = 1000000.0D;
  // private static final int DEFUALT_INSERT_SIG_FIGS = 0;
  // private static final double DEFUALT_MIN_INSERT = 0.0D;
  // private static final double DEFUALT_MAX_INSERT = 2000.0D;
  private final String inputSamOrBamFile;
  private String libraryReadDepthResultsFile;
  private int numDuplicated;
  private int numUnMapped;
  private int numOnTarget;
  private int numTotal;
  private int numUnique;
  private int totalTargetedBasePairs;
  private int numSecondaryAlignments;
  private int numInvalidAligments;
  private double percentDuplicated;
  private double percentUnMapped;
  private double percentOnTarget;
  private double[] percentCoverageAtDepth;
  private double[] totalPercentGCAtDepth;
  private final double[] averageInsertSize;
  private double[] stDevInsertSize;

  private final Histogram.DynamicHistogram[] gHistograms;
  private final Histogram.DynamicHistogram[] insertHistograms;
  private final FilterNGS filterNGS;
  public static final String COMMAND_DIR = "dir=";
  public static final String COMMAND_COMMON_EXT = "commonExt=";
  public static final String COMMAND_INPUT_FILE = "fileOfinputSamOrBams=";
  public static final String COMMAND_TARGET_LIBRARY = "target=";
  public static final String COMMAND_OUTPUT = "output=";
  public static final String COMMAND_SKIP_NUM_LINES = "skip=";
  public static final String COMMAND_NUM_THREADS = "numThreads=";
  public static final String COMMAND_LOG_FILE = "log=";
  public static final String COMMAND_MAPPING_SCORE = "mapQ=";
  public static final String COMMAND_PHREAD_SCORE = "phread=";
  public static final String COMMAND_READ_DEPTH = "readDepth=";
  public static final String COMMAND_BAITS_LIBRARY = "baits=";
  public static final String COMMAND_OUTPUT_DIRECTORY = "outDir=";
  public static final String COMMAND_BAITS_AS_TARGET = "-baits";
  public static final String COMMAND_NORM_READS = "normReads=";

  public BamQC(String inputSamOrBamFile, String outputDir, FilterNGS filterNGS) {
    this.inputSamOrBamFile = inputSamOrBamFile;
    libraryReadDepthResultsFile = SerializedFiles.getSerializedFileName(outputDir,
                                                                        inputSamOrBamFile);
    numTotal = 0;
    numUnique = 0;
    numDuplicated = 0;
    numUnMapped = 0;
    numOnTarget = 0;
    numSecondaryAlignments = 0;
    numInvalidAligments = 0;
    averageInsertSize = new double[NUM_INSERT_HISTOGRAMS];
    totalTargetedBasePairs = 0;
    percentDuplicated = 0;
    percentUnMapped = 0;
    percentOnTarget = 0;
    percentCoverageAtDepth = new double[filterNGS.getReadDepthFilter().length];
    totalPercentGCAtDepth = new double[filterNGS.getReadDepthFilter().length];
    gHistograms = Histogram.DynamicHistogram.initHistograms(NUM_GC_HISTOGRAMS, 0, 1.0, 2);
    insertHistograms = Histogram.DynamicHistogram.initHistograms(NUM_INSERT_HISTOGRAMS, 0.0, 2000.0,
                                                                 0);
    this.filterNGS = filterNGS;
  }

  public void addToGHistoGram(double gcContent, int which) {
    gHistograms[which].addDataPointToHistogram(gcContent);
  }

  public void addToInsertHistoGram(double insertSize, int which) {
    insertHistograms[which].addDataPointToHistogram(insertSize);
  }

  public Histogram.DynamicHistogram getGcHistogram(int which) {
    return gHistograms[which];
  }

  public Histogram.DynamicHistogram getInsertHistogram(int which) {
    return insertHistograms[which];
  }

  public void computePercents() {
    percentDuplicated = computePercent(numDuplicated, numTotal);
    percentUnMapped = computePercent(numUnMapped, numUnique);
    percentOnTarget = computePercent(numOnTarget, numUnique);
    averageInsertSize[NO_QC_INSERT_HISTOGRAM] /= ArrayUtils.sum(getInsertHistogram(NO_QC_INSERT_HISTOGRAM).getCounts());
    averageInsertSize[QC_INSERT_HISTOGRAM] /= ArrayUtils.sum(getInsertHistogram(QC_INSERT_HISTOGRAM).getCounts());

    stDevInsertSize = new double[2];
    stDevInsertSize[NO_QC_INSERT_HISTOGRAM] = getStDevInsertSize(getInsertHistogram(NO_QC_INSERT_HISTOGRAM));
    stDevInsertSize[QC_INSERT_HISTOGRAM] = getStDevInsertSize(getInsertHistogram(QC_INSERT_HISTOGRAM));

  }

  private double getStDevInsertSize(DynamicHistogram dynamicHistogram) {
    ArrayList<Double> inserts = new ArrayList<>();
    double[] bins = dynamicHistogram.getBins();
    for (int i = 0; i < bins.length; i++) {
      int count = dynamicHistogram.getCounts()[i];
      for (int j = 0; j < count; j++) {
        inserts.add(bins[i]);
      }

    }
    double[] finals = Doubles.toArray(inserts);
    return ArrayUtils.stdev(finals);
  }

  public void setPercentCoverageAtDepth(double[] percentCoverageAtDepth) {
    this.percentCoverageAtDepth = percentCoverageAtDepth;
  }

  public double[] getTotalPercentGCAtDepth() {
    return totalPercentGCAtDepth;
  }

  public void setTotalPercentGCAtDepth(double[] totalPercentGCAtDepth) {
    this.totalPercentGCAtDepth = totalPercentGCAtDepth;
  }

  public void addNumTotal() {
    numTotal += 1;
  }

  public double addInsertSize(SAMRecord samRecord, int which) {
    double insert = (0.0D / 0.0D);
    if (samRecord.getProperPairFlag()) {
      insert = Math.abs(samRecord.getInferredInsertSize());
      averageInsertSize[which] += insert;
      addToInsertHistoGram(insert, which);
    }
    return insert;
  }

  public String getLibraryReadDepthResultsFile() {
    return libraryReadDepthResultsFile;
  }

  public void setLibraryReadDepthResultsFile(String libraryReadDepthResultsFile) {
    this.libraryReadDepthResultsFile = libraryReadDepthResultsFile;
  }

  public void setTotalTargetedBasePairs(int totalTargetedBasePairs) {
    this.totalTargetedBasePairs = totalTargetedBasePairs;
  }

  public void addNumUnique(SAMRecord samRecord) {
    if (!samRecord.getDuplicateReadFlag()) {
      numUnique += 1;
    }
  }

  public void addNumUnMapped(SAMRecord samRecord) {
    if (samRecord.getMateUnmappedFlag()) {
      numUnMapped += 1;
    }
  }

  public void addNumDuplicate(SAMRecord samRecord) {
    if (samRecord.getDuplicateReadFlag()) {
      numDuplicated += 1;
    }
  }

  public void addNumOnTarget() {
    numOnTarget += 1;
  }

  public String getInputSamOrBamFile() {
    return inputSamOrBamFile;
  }

  public int getNumDuplicated() {
    return numDuplicated;
  }

  public int getNumUnMapped() {
    return numUnMapped;
  }

  public int getNumOnTarget() {
    return numOnTarget;
  }

  public int getNumTotalReads() {
    return numTotal;
  }

  public double getPercentDuplicated() {
    return percentDuplicated;
  }

  public double getPercentUnMapped() {
    return percentUnMapped;
  }

  public double getPercentOnTarget() {
    return percentOnTarget;
  }

  public int getNumSecondaryAlignments() {
    return numSecondaryAlignments;
  }

  public int getNumInvalidAligments() {
    return numInvalidAligments;
  }

  public boolean secondaryAlignment(SAMRecord samRecord) {
    if (samRecord.isSecondaryOrSupplementary()) {
      numSecondaryAlignments += 1;
      return true;
    }
    return false;
  }

  public boolean invalidAlignment(SAMRecord samRecord) {
    if (samRecord.isValid() != null) {
      numInvalidAligments += 1;
      return true;
    }
    return false;
  }

  public String getSummary() {
    String summary = "";
    summary = summary + inputSamOrBamFile + "\t" + numTotal + "\t" + numUnique + "\t"
              + numDuplicated + "\t" + numUnMapped + "\t" + numOnTarget + "\t"
              + numSecondaryAlignments + "\t" + numInvalidAligments + "\t" + percentDuplicated
              + "\t" + percentUnMapped + "\t" + percentOnTarget + "\t"
              + averageInsertSize[NO_QC_INSERT_HISTOGRAM] + "\t" + "\t"
              + stDevInsertSize[QC_INSERT_HISTOGRAM] + "\t" + filterNGS.getMappingQualityFilter()
              + "\t" + filterNGS.getPhreadScoreFilter() + "\t" + totalTargetedBasePairs + "\t"
              + ArrayUtils.toStr(percentCoverageAtDepth) + "\t"
              + ArrayUtils.toStr(totalPercentGCAtDepth);
    return summary;
  }

  public static BamQC[] qcBams(String dir, String outputDir, String commonExt,
                               String fileOfinputSamOrBams, String targetLibraryFile,
                               String baitLibraryFile, int skipNumLines, FilterNGS filterNGS,
                               int numThreads, String output, String snpEffLocation,
                               boolean baitsAsTarget, double normalizeDepthsTo,
                               boolean summarizeLib, Logger log) {
    long time = System.currentTimeMillis();
    log.report(ext.getTime() + " Info - beginning qc");
    BamQC[] bamQCs = null;
    if ((dir == null) && (fileOfinputSamOrBams == null)) {
      log.reportError("Error - must provide a directory to search or a file listing files to qc");
    } else {
      String[] inputbams = determineInputFiles(dir, commonExt, fileOfinputSamOrBams, log);
      if ((inputbams == null) || (inputbams.length < 1)) {
        log.reportError("Error - did not find any files to qc");
      } else {
        LibraryNGS libraryNGS = null;
        LibraryNGS.BaitsLibrary lBaitsLibrary = null;
        if ((baitsAsTarget) && ((baitLibraryFile == null) || (!Files.exists(baitLibraryFile)))) {
          log.reportError("Error - baits a library was flagged, but could not find the baits file"
                          + (baitLibraryFile == null ? "" : baitLibraryFile));
          return null;
        }
        if ((targetLibraryFile == null) && (!baitsAsTarget)) {
          log.report("Warning - a target library file was not provided, on target percentages will not be computed");
        } else {
          if (baitLibraryFile != null) {
            lBaitsLibrary = LibraryNGS.BaitsLibrary.loadBaitLibrary(baitLibraryFile, log);
          }
          if (baitsAsTarget) {
            libraryNGS = new LibraryNGS(lBaitsLibrary.getBaits(), log);
          } else {
            libraryNGS = LibraryNGS.getLibraryNGS(targetLibraryFile, skipNumLines, log);
          }
          if (lBaitsLibrary != null) {
            log.report(ext.getTime() + " Info - mapping baits library...");
            libraryNGS.mapBaits(lBaitsLibrary, baitsAsTarget);
          }
        }
        bamQCs = qcBams(inputbams, outputDir, libraryNGS, filterNGS, normalizeDepthsTo, numThreads,
                        log);
        summarize(bamQCs, outputDir, filterNGS, output, log);
        String librarySummary = outputDir
                                + ext.addToRoot(output,
                                                baitsAsTarget ? ".libraryBaitsResults.summary"
                                                              : ".libraryResults.summary");
        if (summarizeLib) {
          LibraryNGS.summarizeLibraries(libraryNGS, getLibraryReadDepthResultsFiles(bamQCs),
                                        librarySummary, filterNGS, log);
          SNPEFF snpeff = new SNPEFF(snpEffLocation, true, true, log);
          snpeff.runSnpEFFCount(inputbams,
                                outputDir + ext.addToRoot(output,
                                                          new StringBuilder(String.valueOf(baitsAsTarget ? ".libraryBaitsResults.summary"
                                                                                                         : ".libraryResults.summary")).append("count")
                                                                                                                                      .toString()),
                                ANNOTATION_BUILD.HG19.getSnpEffBuild(),
                                baitsAsTarget ? null : targetLibraryFile, numThreads);
        }
      }
    }
    log.report(ext.getTime() + " Info - finished qc in " + ext.getTimeElapsed(time));
    return bamQCs;
  }

  public static BamQC[] qcBams(String[] inputbams, String outputDir, LibraryNGS libraryNGS,
                               FilterNGS filterNGS, double normalizeDepthsTo, int numThreads,
                               Logger log) {
    BamQC[] bamQCs = new BamQC[inputbams.length];
    ExecutorService executor = Executors.newFixedThreadPool(numThreads);
    Hashtable<String, Future<BamQC>> tmpResults = new Hashtable<>();
    for (int i = 0; i < bamQCs.length; i++) {
      tmpResults.put(i + "", executor.submit(new WorkerBamQC(inputbams[i], outputDir, libraryNGS,
                                                             filterNGS, normalizeDepthsTo, log)));
    }
    for (int i = 0; i < bamQCs.length; i++) {
      try {
        try {
          bamQCs[i] = tmpResults.get(i + "").get();
        } catch (NullPointerException npe) {
          log.reportError("Could not get qc for " + inputbams[i]);
        }
      } catch (InterruptedException e) {
        log.reportError("Error - in file " + inputbams[i]);
        log.reportException(e);
      } catch (ExecutionException e) {
        log.reportError("Error - in file " + inputbams[i]);
        log.reportException(e);
      }
    }
    executor.shutdown();
    try {
      executor.awaitTermination(10L, TimeUnit.DAYS);
    } catch (InterruptedException e) {
      log.reportException(e);
    }
    return bamQCs;
  }

  private static String[] getLibraryReadDepthResultsFiles(BamQC[] bamQCs) {
    String[] libraryReadDepthResultsFiles = new String[bamQCs.length];
    for (int i = 0; i < bamQCs.length; i++) {
      libraryReadDepthResultsFiles[i] = bamQCs[i].getLibraryReadDepthResultsFile();
    }
    return libraryReadDepthResultsFiles;
  }

  public static BamQC qcBam(String inputSamOrBamFile, String outputDir, LibraryNGS libraryNGS,
                            FilterNGS filterNGS, double normalizeDepthsTo, Logger log) {

    LibraryNGS.ReadDepth readDepth = null;
    if (libraryNGS != null) {
      readDepth = new LibraryNGS.ReadDepth(libraryNGS, filterNGS, log);
    }
    BamQC bamQC = initBamQC(inputSamOrBamFile, outputDir, filterNGS);
    System.out.println(bamQC.getLibraryReadDepthResultsFile());
    // if (!Files.exists(bamQC.getLibraryReadDepthResultsFile())) {

    SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
    samReaderFactory.validationStringency(ValidationStringency.LENIENT);
    SamReader reader = samReaderFactory.open(new File(inputSamOrBamFile));
    log.report(ext.getTime() + " Info - beginning processing of " + inputSamOrBamFile);
    int numTotalReads = 0;
    for (SAMRecord samRecord : reader) {
      if (!bamQC.invalidAlignment(samRecord)) {
        bamQC.addInsertSize(samRecord, NO_QC_INSERT_HISTOGRAM);
        if ((samRecord.getReadUnmappedFlag()) || (!bamQC.secondaryAlignment(samRecord))) {
          double currentGC = getGCContent(samRecord);
          bamQC.addToGHistoGram(currentGC, NO_QC_GC_HISTOGRAM);
          numTotalReads++;
          if (passesMapQ(samRecord, filterNGS)) {
            bamQC.addNumTotal();
            bamQC.addNumDuplicate(samRecord);
            bamQC.addNumUnique(samRecord);
            bamQC.addNumUnMapped(samRecord);
            if ((!samRecord.getReferenceName().equals("*")) && (libraryNGS != null)
                && (readDepth != null)) {
              Segment segment = new Segment(Positions.chromosomeNumber(samRecord.getReferenceName(),
                                                                       log),
                                            samRecord.getAlignmentStart(),
                                            samRecord.getAlignmentEnd());
              if (segment.getChr() >= 0) {
                int[] libraryIndices = libraryNGS.indicesInLibrary(segment);
                if ((!samRecord.getDuplicateReadFlag()) && (!samRecord.getReadUnmappedFlag())
                    && (libraryIndices != null)) {
                  bamQC.addInsertSize(samRecord, QC_INSERT_HISTOGRAM);
                  bamQC.addToGHistoGram(currentGC, QC_GC_HISTOGRAM);
                  for (int libraryIndice : libraryIndices) {
                    readDepth.addCounts(libraryIndice,
                                        libraryNGS.getTargetSegments()[libraryIndice], samRecord,
                                        filterNGS);
                  }
                  bamQC.addNumOnTarget();
                }
              }
            }
            if ((numTotalReads != 0) && (numTotalReads % 1000000 == 0)) {
              float usedMemory = Runtime.getRuntime().totalMemory()
                                 - Runtime.getRuntime().freeMemory();
              float freeMemory = Runtime.getRuntime().maxMemory() - usedMemory;
              float maxMemory = Runtime.getRuntime().maxMemory();
              log.report(ext.getTime() + "Info - processed " + numTotalReads + " total reads from "
                         + inputSamOrBamFile + "\nFree memory: "
                         + Math.round(freeMemory / maxMemory * 100.0D) + "%");
            }
          }
        }
      }
    }
    try {
      reader.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
    bamQC.computePercents();
    if ((libraryNGS != null) && (readDepth != null)) {
      double normalizeFactor = 1.0D;
      if (normalizeDepthsTo != 0.0D) {
        normalizeFactor = normalizeDepthsTo * 1000000.0D / bamQC.getNumOnTarget();
        log.report("Info -normalizing depth results from " + bamQC.getInputSamOrBamFile() + " with "
                   + bamQC.getNumOnTarget() + " reads to " + normalizeDepthsTo * 1000000.0D
                   + " reads");
        log.report("Info -normalizing depth results from " + bamQC.getInputSamOrBamFile()
                   + " using a factor of (depth)*" + normalizeFactor);
      }
      LibraryNGS.LibraryReadDepthResults lDepthResults = readDepth.getDepthResults(libraryNGS,
                                                                                   filterNGS,
                                                                                   normalizeFactor);
      bamQC.setPercentCoverageAtDepth(lDepthResults.getTotalPercentCoveredAtDepth());
      bamQC.setTotalTargetedBasePairs(lDepthResults.getTotalBasePairsTargeted());
      bamQC.setTotalPercentGCAtDepth(lDepthResults.getTotalPercentCoveredAtDepth());
      lDepthResults.dump(ext.rootOf(bamQC.getLibraryReadDepthResultsFile(), false) + ".txt", log);
      lDepthResults.serialize(bamQC.getLibraryReadDepthResultsFile());
    } else {
      double[] nans = new double[filterNGS.getReadDepthFilter().length];
      Arrays.fill(nans, Double.NaN);
      bamQC.setPercentCoverageAtDepth(nans);
      bamQC.setTotalTargetedBasePairs(0);
    }
    log.report(ext.getTime() + " Summary for: " + inputSamOrBamFile + "\n");
    log.report(ext.getTime() + " Summary:" + bamQC.getSummary());
    if (bamQC.getNumInvalidAligments() > 0) {
      log.reportError("Warning - " + bamQC.getNumInvalidAligments()
                      + " read(s) were invalid and were skipped in file " + inputSamOrBamFile);
    }
    if (bamQC.getNumSecondaryAlignments() > 0) {
      log.report("Info - " + bamQC.getNumSecondaryAlignments() + " secondary alignment(s) in "
                 + inputSamOrBamFile + " were skipped for computing QC metrics");
    }
    // } else {
    // log.reportTimeInfo("loading pre-computed library results");
    // LibraryNGS.LibraryReadDepthResults lDepthResults =
    // LibraryNGS.LibraryReadDepthResults.load(bamQC.getLibraryReadDepthResultsFile(), false);
    // bamQC.setPercentCoverageAtDepth(lDepthResults.getTotalPercentCoveredAtDepth());
    // bamQC.setTotalTargetedBasePairs(lDepthResults.getTotalBasePairsTargeted());
    // }
    return bamQC;
  }

  private static double getGCContent(SAMRecord samRecord) {
    String sr = samRecord.getReadString();
    int gscs = 0;
    for (int i = 0; i < sr.length(); i++) {
      String base = sr.charAt(i) + "";
      if ((base.equals("G")) || (base.equals("C"))) {
        gscs++;
      }
    }
    return gscs / sr.length();
  }

  private static class WorkerBamQC implements Callable<BamQC> {

    private final String inputSamOrBamFile;
    private final String outputDir;
    private final LibraryNGS libraryNGS;
    private final FilterNGS filterNGS;
    private final double normalizeDepthsTo;
    private final Logger log;

    public WorkerBamQC(String inputSamOrBamFile, String outputDir, LibraryNGS libraryNGS,
                       FilterNGS filterNGS, double normalizeDepthsTo, Logger log) {
      this.inputSamOrBamFile = inputSamOrBamFile;
      this.outputDir = outputDir;
      this.libraryNGS = libraryNGS;
      this.filterNGS = filterNGS;
      this.normalizeDepthsTo = normalizeDepthsTo;
      this.log = log;
    }

    @Override
    public BamQC call() {
      return BamQC.qcBam(inputSamOrBamFile, outputDir, libraryNGS, filterNGS, normalizeDepthsTo,
                         log);
    }
  }

  private static String[] determineInputFiles(String dir, String commonExt,
                                              String fileOfinputSamOrBams, Logger log) {
    String[] inputbams = null;
    if (fileOfinputSamOrBams != null) {
      if (Files.exists(fileOfinputSamOrBams)) {
        log.report(ext.getTime() + " Info - loading files to qc from " + fileOfinputSamOrBams);
        inputbams = HashVec.loadFileToStringArray(fileOfinputSamOrBams, false, new int[] {0},
                                                  false);
        if (!Files.exists("", inputbams)) {
          log.reportError("Error - found files that do not exist in " + fileOfinputSamOrBams);
          inputbams = null;
        }
      } else {
        log.reportError("Error - could not find file " + fileOfinputSamOrBams);
      }
    } else {
      log.report(ext.getTime() + " Info - loading all files from dir " + dir + " with extension "
                 + commonExt);
      inputbams = Files.toFullPaths(Files.list(dir, "", commonExt, true), dir);
    }
    return inputbams;
  }

  private static double computePercent(int top, int bottom) {
    return bottom > 0 ? (double) top / bottom : 0;
  }

  private static BamQC initBamQC(String inputSamOrBamFile, String outputDir, FilterNGS filterNGS) {
    if (filterNGS == null) {
      return new BamQC(inputSamOrBamFile, outputDir, new FilterNGS(0.0, 0.0, new int[1]));
    }
    BamQC bamQC = new BamQC(inputSamOrBamFile, outputDir, filterNGS);
    return bamQC;
  }

  private static void summarize(BamQC[] bamQCs, String outputDir, FilterNGS filterNGS,
                                String output, Logger log) {
    if (outputDir == null) {
      for (BamQC bamQC : bamQCs) {
        if (bamQC != null) {
          output = ext.parseDirectoryOfFile(bamQC.getInputSamOrBamFile())
                   + ext.removeDirectoryInfo(output);
          break;
        }
      }
    } else {
      output = outputDir + output;
    }
    summarizeQC(bamQCs, filterNGS, output, log);
    summarizeGCContent(bamQCs, ext.addToRoot(output, "raw"), NO_QC_GC_HISTOGRAM, log);
    summarizeGCContent(bamQCs, ext.addToRoot(output, "onTarget"), QC_GC_HISTOGRAM, log);
    summarizeInsertSize(bamQCs, ext.addToRoot(output, "raw"), NO_QC_INSERT_HISTOGRAM, log);
    summarizeInsertSize(bamQCs, ext.addToRoot(output, "onTarget"), QC_INSERT_HISTOGRAM, log);
  }

  private static void summarizeQC(BamQC[] bamQCs, FilterNGS filterNGS, String output, Logger log) {
    try {
      PrintWriter writer = Files.openAppropriateWriter(output);
      writer.print(ArrayUtils.toStr(QC_HEADER));
      for (int i = 0; i < filterNGS.getReadDepthFilter().length; i++) {
        writer.print("\tPercent Coverage At " + filterNGS.getReadDepthFilter()[i]);
      }
      for (int i = 0; i < filterNGS.getReadDepthFilter().length; i++) {
        writer.print("\tPercent GC At coverage " + filterNGS.getReadDepthFilter()[i]);
      }
      writer.println();
      for (BamQC bamQC : bamQCs) {
        if (bamQC != null) {
          writer.println(bamQC.getSummary());
        }
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + output);
      log.reportException(e);
    }
  }

  private static void summarizeGCContent(BamQC[] bamQCs, String output, int which, Logger log) {
    try {
      PrintWriter writer = Files.openAppropriateWriter(ext.rootOf(output, false) + ".gcs");
      double[] bins = null;
      writer.print(ArrayUtils.toStr(HIST_HEADER));
      for (BamQC bamQC : bamQCs) {
        if (bamQC != null) {
          writer.print("\t" + ext.rootOf(bamQC.getInputSamOrBamFile()));
          if (bins == null) {
            bins = bamQC.getGcHistogram(which).getBins();
          }
        }
      }
      writer.println();
      for (int i = 0; i < bins.length; i++) {
        writer.print(bins[i]);
        for (BamQC bamQC : bamQCs) {
          if (bamQC != null) {
            writer.print("\t" + bamQC.getGcHistogram(which).getCounts()[i]);
          }
        }
        writer.println();
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + output);
      log.reportException(e);
    }
  }

  private static void summarizeInsertSize(BamQC[] bamQCs, String output, int which, Logger log) {
    try {
      PrintWriter writer = Files.openAppropriateWriter(ext.rootOf(output, false) + ".insert");
      double[] bins = null;
      writer.print(ArrayUtils.toStr(HIST_HEADER));
      for (BamQC bamQC : bamQCs) {
        if (bamQC != null) {
          writer.print("\t" + ext.rootOf(bamQC.getInputSamOrBamFile()));
          if (bins == null) {
            bins = bamQC.getInsertHistogram(which).getBins();
          }
        }
      }
      writer.println();
      for (int i = 0; i < bins.length; i++) {
        writer.print(bins[i]);
        for (BamQC bamQC : bamQCs) {
          if (bamQC != null) {
            writer.print("\t" + bamQC.getInsertHistogram(which).getCounts()[i]);
          }
        }
        writer.println();
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + output);
      log.reportException(e);
    }
  }

  private static boolean passesMapQ(SAMRecord samRecord, FilterNGS filterNGS) {
    boolean passesMapQ = false;
    if (filterNGS.getMappingQualityFilter() == 0.0) {
      passesMapQ = true;
    } else if (((filterNGS.getMappingQualityFilter() == 0.0)
                || (samRecord.getMappingQuality() != 255))
               && (samRecord.getMappingQuality() >= filterNGS.getMappingQualityFilter())) {
      passesMapQ = true;
    }
    return passesMapQ;
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String dir = null;
    String commonExt = ".bam";
    String fileOfinputSamOrBams = null;
    String targetLibraryFile = null;
    String baitLibraryFile = null;
    String output = "bamQC.txt";
    String outputDir = null;
    int skipNumLines = 2;
    int numThreads = 4;
    double mappingQuality = 0;
    double phreadScore = 0.0;
    int[] readDepth = {0, 1, 2, 3, 4, 10, 20, 30, 40};
    boolean baitsAsTarget = false;
    double normalizeDepthsTo = 0.0D;
    String snpEffLocation = null;
    String logfile = "bamQC.log";

    String usage = "\nseq.BamQC requires 0-1 arguments\n";
    usage = usage + "    NOTE: all bam files are assumed to be sorted:\n";
    usage = usage + "   (1) directory containing bam or sam files to qc (i.e. dir=" + dir
            + " (no default))\n";
    usage = usage
            + "   (2) full path to a file listing bam or sam files (one per line, in first column, no header) to qc (i.e. fileOfinputSamOrBams="
            + fileOfinputSamOrBams + " (no default))\n";
    usage = usage
            + "   (3) full path to a target library file to compute on target percentage (i.e. target="
            + targetLibraryFile + " (no default))\n";
    usage = usage
            + "   (3) full path to a bait library file to map to the target library (i.e. baits="
            + baitLibraryFile + " (no default))\n";
    usage = usage + "   (3) full path to output directory (i.e. outDir=" + outputDir
            + " (defualts to directory of first file qc -ed))\n";

    usage = usage + "   (4) output summary file name (i.e. output=" + output
            + " (default - defaults to directory of first file qc -ed))\n";
    usage = usage + "   (5) file extension to search for in a directory (i.e. commonExt="
            + commonExt + " (default))\n";
    usage = usage
            + "   (6) if using a target library file, number of lines to skip before targets are listed (i.e. skip="
            + skipNumLines + " (default))\n";
    usage = usage + "   (7) number of threads  (i.e. numThreads=" + numThreads + " (default))\n";
    usage = usage + "   (8) log file  (i.e. log=" + logfile + " (default))\n";
    usage = usage + "   (9) mapping quality score filter to compute qc metrics with (i.e. mapQ="
            + mappingQuality + " (default))\n";
    usage = usage + "   (10) phread score filter to compute qc metrics with (i.e. phread="
            + phreadScore + " (default))\n";
    usage = usage
            + "   (11) read depth filter (comma-delimited if multiple) to compute qc metrics with (i.e. readDepth="
            + ArrayUtils.toStr(ArrayUtils.toStringArray(readDepth), ",") + " (default))\n";
    usage = usage + "   (12) compute QC against the baits library (i.e. -baits" + baitsAsTarget
            + " (default))\n";
    usage = usage + "   (13) normalize coverage to number of Reads (in millions)  (i.e. normReads="
            + normalizeDepthsTo + " (default, no normalizing))\n";
    usage = usage + "   (14) full path to the SNP EFF directory (i.e. snpEff= ( no default))\n";
    for (String arg : args) {
      if ((arg.equals("-h")) || (arg.equals("-help")) || (arg.equals("/h"))
          || (arg.equals("/help"))) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("fileOfinputSamOrBams=")) {
        fileOfinputSamOrBams = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("target=")) {
        targetLibraryFile = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("-baits")) {
        baitsAsTarget = true;
        numArgs--;
      } else if (arg.startsWith("outDir=")) {
        outputDir = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("output=")) {
        output = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("commonExt=")) {
        commonExt = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("baits=")) {
        baitLibraryFile = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("skip=")) {
        skipNumLines = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("numThreads=")) {
        numThreads = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("mapQ=")) {
        mappingQuality = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("phread=")) {
        phreadScore = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("phread=")) {
        phreadScore = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("normReads=")) {
        normalizeDepthsTo = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("snpEff=")) {
        snpEffLocation = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("readDepth=")) {
        readDepth = ArrayUtils.toIntArray(ext.parseStringArg(arg, "").split(","));
        phreadScore = ext.parseDoubleArg(arg);
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      Logger log = new Logger(logfile);
      FilterNGS filterNGS = new FilterNGS(mappingQuality, phreadScore, readDepth);
      qcBams(dir, outputDir, commonExt, fileOfinputSamOrBams, targetLibraryFile, baitLibraryFile,
             skipNumLines, filterNGS, numThreads, output, snpEffLocation, baitsAsTarget,
             normalizeDepthsTo, true, log);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}