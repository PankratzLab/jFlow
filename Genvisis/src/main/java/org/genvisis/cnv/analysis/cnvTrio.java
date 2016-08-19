package org.genvisis.cnv.analysis;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Vector;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.qc.CNVTrioFilter;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.CNVFilter;
import org.genvisis.common.CNVFilter.CNVFilterPass;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.Segment;

import com.google.common.primitives.Bytes;
import com.google.common.primitives.Doubles;

// TODO, split the centromeres or remove?
/**
 * Class for filtering denovo calls in offspring by BEAST SCORE and LRR_SD, and a few other metrics
 * Filtering at the default metrics here seems to work OK
 */
public class cnvTrio extends CNVariant {
  private static final long serialVersionUID = 1L;
  private static final String[] SPLITS = {"\t"};
  private static final String[] COMBINED_TRIOS = {"denovo_joint.cnv", "denovo_trio.cnv",
                                                  ".filtered", ".cnv"};
  private static final String[] TRIOS = {".jointcnv", ".triocnv"};
  private static final String[] BEAST_OUTPUT = {".trio.cnv", ".beast.cnv", ".filtered", ".qc",
                                                "parental.cnv"};
  public static final String[] BEAST_HEADER = {"NUM_MARKERS", "IDNA", "ILRR_SD", "IBAF1585_SD",
                                               "IHEIGHT", "INumbProbes", "IBEAST", "FADNA",
                                               "FALRR_SD", "FABAF1585_SD", "FAHEIGHT",
                                               "FANumProbes", "FABEAST", "MODNA", "MOLRR_SD",
                                               "MOBAF1585_SD", "MOHEIGHT", "MONumProbes", "MOBEAST",
                                               "MinBeastDiff", "MaxTrioLRR_SD", "MaxTrioBAF1585_SD",
                                               "NumCalls", "DenovoLength(BP)", "DenovoMedianLRR",
                                               "NumRawOverlappingCNVs"};
  public static final String[] OUTPUT_HEADER = {"IDNA", "FADNA", "MODNA", "ILRR_SD", "FALRR_SD",
                                                "MOLRR_SD", "IBAF1585_SD", "FABAF1585_SD",
                                                "MOBAF1585_SD", "IHEIGHT", "FAHEIGHT", "MOHEIGHT",
                                                "IMedianLRR", "FAMedianLRR", "MOMedianLRR",
                                                "NumCalls", "INumMarkers", "FANumMarkers",
                                                "MONumMarkers", "IBEAST", "FABEAST", "MOBEAST",
                                                "NumRawOverlappingCNVs"};
  // TODO better map for OUT and Filt
  public static final String[] FILTERED_OUTPUT_HEADER = {"IDNA", "FADNA", "MODNA", "ILRR_SD",
                                                         "FALRR_SD", "MOLRR_SD", "IBAF1585_SD",
                                                         "FABAF1585_SD", "MOBAF1585_SD", "IHEIGHT",
                                                         "FAHEIGHT", "MOHEIGHT", "IMedianLRR",
                                                         "FAMedianLRR", "MOMedianLRR", "NumCalls",
                                                         "MinBeastDiff", "UCSC", "UCSC_LINK"};

  public static final String COMMAND_PARSE = "-parse";
  public static final String COMMAND_PROJECT = "proj=";
  public static final String COMMAND_TRIO_RESULTS = "trioResults=";
  public static final String COMMAND_FILTER_FILE = "filterFile=";
  public static final String COMMAND_BUILD = "build=";
  public static final String COMMAND_EXCLUDE = "-exclude";
  public static final String COMMAND_OUTPUT = "output=";

  public static final int DEFAULT_BUILD = 36;

  private final String IDNA, FADNA, MODNA;
  private final double ILRR_SD, FALRR_SD, MOLRR_SD;
  private final double IBAF1585_SD, FABAF1585_SD, MOBAF1585_SD;
  private final double IHEIGHT, FAHEIGHT, MOHEIGHT;
  private final double iMedianLRR, fAMedianLRR, mOMedianLRR;
  private double minBeastHeightDifference;
  private final int numCalls;

  public cnvTrio(String[] plinkLine, String iDNA, String fADNA, String mODNA, double iLRR_SD,
                 double fALRR_SD, double mOLRR_SD, double iBAF1585_SD, double fABAF1585_SD,
                 double mOBAF1585_SD, double iHEIGHT, double fAHEIGHT, double mOHEIGHT,
                 double iMedianLRR, double fAMedianLRR, double mOMedianLRR, int numCalls) {
    super(plinkLine);
    IDNA = iDNA;
    FADNA = fADNA;
    MODNA = mODNA;
    ILRR_SD = iLRR_SD;
    FALRR_SD = fALRR_SD;
    MOLRR_SD = mOLRR_SD;
    IBAF1585_SD = iBAF1585_SD;
    FABAF1585_SD = fABAF1585_SD;
    MOBAF1585_SD = mOBAF1585_SD;
    IHEIGHT = iHEIGHT;
    FAHEIGHT = fAHEIGHT;
    MOHEIGHT = mOHEIGHT;
    this.numCalls = numCalls;
    this.iMedianLRR = iMedianLRR;
    this.fAMedianLRR = fAMedianLRR;
    this.mOMedianLRR = mOMedianLRR;
  }

  public void computeMinHeightDist(double maxBeastHeightParents) {
    if (Math.abs(FAHEIGHT) > maxBeastHeightParents || Math.abs(MOHEIGHT) > maxBeastHeightParents
        || Math.abs(IHEIGHT) < maxBeastHeightParents) {
      setMinBeastHeightDifference(0);
    }
    setMinBeastHeightDifference(Math.min(Math.abs(IHEIGHT
                                                  - checkParentOppositeHeight(IHEIGHT, FAHEIGHT)),
                                         Math.abs(IHEIGHT
                                                  - checkParentOppositeHeight(IHEIGHT, MOHEIGHT))));
  }

  public double getMinBeastHeightDifference() {
    return minBeastHeightDifference;
  }

  public String getFullSummary(String build) {
    String summary = "";
    summary += toPlinkFormat();
    summary += "\t" + IDNA;
    summary += "\t" + FADNA;
    summary += "\t" + MODNA;
    summary += "\t" + ILRR_SD;
    summary += "\t" + FALRR_SD;
    summary += "\t" + MOLRR_SD;
    summary += "\t" + IBAF1585_SD;
    summary += "\t" + FABAF1585_SD;
    summary += "\t" + MOBAF1585_SD;
    summary += "\t" + IHEIGHT;
    summary += "\t" + FAHEIGHT;
    summary += "\t" + MOHEIGHT;
    summary += "\t" + iMedianLRR;
    summary += "\t" + fAMedianLRR;
    summary += "\t" + mOMedianLRR;
    summary += "\t" + numCalls;
    summary += "\t" + minBeastHeightDifference;
    summary += "\t" + getUCSClocation();
    summary += "\t" + getUCSCLink(build);
    return summary;
  }

  public String getCNV(String familyID_individualID) {
    return getCNV(familyID_individualID.split("\t")[0],
                  familyID_individualID.split("\t")[1]).toPlinkFormat();
  }

  public CNVariant getCNV(String familyID, String individualID) {
    return new CNVariant(familyID, individualID, getChr(), getStart(), getStop(), getCN(),
                         getScore(), getNumMarkers(), getSource());
  }

  public String getListSummary() {
    String summary = "";
    summary += IDNA + "\t" + getUCSClocation() + "\tHeightDiff="
               + ext.formDeci(minBeastHeightDifference, 2) + ";LRR_SD=" + ext.formDeci(ILRR_SD, 2)
               + ";MedianLRR=" + ext.formDeci(iMedianLRR, 2) + ";" + "CHILD\n";
    summary += FADNA + "\t" + getUCSClocation() + "\tHeightDiff="
               + ext.formDeci(minBeastHeightDifference, 2) + ";LRR_SD=" + ext.formDeci(FALRR_SD, 2)
               + ";MedianLRR=" + ext.formDeci(fAMedianLRR, 2) + ";" + "FATHER\n";
    summary += MODNA + "\t" + getUCSClocation() + "\tHeightDiff="
               + ext.formDeci(minBeastHeightDifference, 2) + ";LRR_SD=" + ext.formDeci(MOLRR_SD, 2)
               + ";MedianLRR=" + ext.formDeci(mOMedianLRR, 2) + ";" + "MOTHER";
    return summary;
  }

  public void setMinBeastHeightDifference(double minBeastHeightDifference) {
    this.minBeastHeightDifference = minBeastHeightDifference;
  }

  public int getNumCalls() {
    return numCalls;
  }

  public String getIDNA() {
    return IDNA;
  }

  public String getFADNA() {
    return FADNA;
  }

  public String getMODNA() {
    return MODNA;
  }

  public double getILRR_SD() {
    return ILRR_SD;
  }

  public double getFALRR_SD() {
    return FALRR_SD;
  }

  public double getMOLRR_SD() {
    return MOLRR_SD;
  }

  public double getIBAF1585_SD() {
    return IBAF1585_SD;
  }

  public double getFABAF1585_SD() {
    return FABAF1585_SD;
  }

  public double getMOBAF1585_SD() {
    return MOBAF1585_SD;
  }

  public double getIHEIGHT() {
    return IHEIGHT;
  }

  public double getFAHEIGHT() {
    return FAHEIGHT;
  }

  public double getMOHEIGHT() {
    return MOHEIGHT;
  }

  private static cnvTrio[] loadTrios(String fullPathToTrioFile, Logger log) {
    ArrayList<cnvTrio> tmpTrios = new ArrayList<cnvTrio>();
    try {
      BufferedReader reader = Files.getAppropriateReader(fullPathToTrioFile);
      String[] header = reader.readLine().trim().split("\t");
      int[] cnvIndices = ext.indexFactors(PLINK_CNV_HEADER, header, true, log, true, false);
      int[] trioIndics = ext.indexFactors(OUTPUT_HEADER, header, true, false);
      boolean allThere = true;
      for (int i = 0; i < PLINK_CNV_HEADER.length; i++) {
        if (cnvIndices[i] < 0) {
          log.reportError("Error - could not find column " + PLINK_CNV_HEADER[i] + " in file "
                          + fullPathToTrioFile);
          allThere = false;
        }
      }
      for (int i = 0; i < OUTPUT_HEADER.length; i++) {
        if (trioIndics[i] < 0) {
          log.reportError("Error - could not find column " + OUTPUT_HEADER[i] + " in file "
                          + fullPathToTrioFile);
          allThere = false;
        }
      }
      if (!allThere) {
        return null;
      }

      while (reader.ready()) {
        String[] line = reader.readLine().trim().split("\t");

        String[] plinkLine = Array.subArray(line, cnvIndices);
        String[] trio = Array.subArray(line, trioIndics);
        String iDNA = trio[0];
        String fADNA = trio[1];
        String mODNA = trio[2];

        try {
          double iLRR_SD = Double.parseDouble(trio[3]);
          double fALRR_SD = Double.parseDouble(trio[4]);
          double mOLRR_SD = Double.parseDouble(trio[5]);
          double iBAF1585_SD = Double.parseDouble(trio[6]);
          double fABAF1585_SD = Double.parseDouble(trio[7]);
          double mOBAF1585_SD = Double.parseDouble(trio[8]);
          double iHEIGHT = Double.parseDouble(trio[9]);
          double fAHEIGHT = Double.parseDouble(trio[10]);
          double mOHEIGHT = Double.parseDouble(trio[11]);
          double iMedianLRR = Double.parseDouble(trio[12]);
          double fAMedianLRR = Double.parseDouble(trio[13]);
          double mOMedianLRR = Double.parseDouble(trio[14]);
          int numCalls = Integer.parseInt(trio[15]);
          tmpTrios.add(new cnvTrio(plinkLine, iDNA, fADNA, mODNA, iLRR_SD, fALRR_SD, mOLRR_SD,
                                   iBAF1585_SD, fABAF1585_SD, mOBAF1585_SD, iHEIGHT, fAHEIGHT,
                                   mOHEIGHT, iMedianLRR, fAMedianLRR, mOMedianLRR, numCalls));

        } catch (NumberFormatException nfe) {
          log.reportError("Error - improper number on line " + line);
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + fullPathToTrioFile + "\" not found in current directory");
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + fullPathToTrioFile + "\"");
    }
    return tmpTrios.toArray(new cnvTrio[tmpTrios.size()]);
  }

  /**
   * Steps are:
   * <p>
   * 1. load the trios to Trio[]
   * <p>
   * 2. Parse the penncnv trio/joint results and collect in one .cnv file using
   * DeNovoCNV.parsePennCnvResult, or supply existing cnv file
   * <p>
   * 3. assign variants to each trio in Trio[]
   * <p>
   * 4. process trios (compute beast and other qc metrics), multithreaded
   * <p>
   * 5. summarize the trio calls into output file
   * <p>
   */
  public static void computeMetrics(Project proj, String trioFile, int fileType, String cnvFile,
                                    String output, int numThreads) {
    Logger log = proj.getLog();
    TrioQC[] trios = TrioQC.loadTrios(proj, trioFile);
    CNVariant[] cnVariants;
    Hashtable<String, Integer> rawRegionFrequency = new Hashtable<String, Integer>();
    if (cnvFile == null) {
      DeNovoCNV.parsePennCnvResult(proj, proj.PENNCNV_RESULTS_DIRECTORY.getValue(false, true),
                                   proj.DATA_DIRECTORY.getValue(false, true) + trioFile,
                                   TRIOS[fileType]);
      rawRegionFrequency = determinRegionFrequencies(CNVariant.loadPlinkFile(
                                                                             proj.DATA_DIRECTORY.getValue(false,
                                                                                                          true)
                                                                             + COMBINED_TRIOS[fileType],
                                                                             false),
                                                     proj.getLog());
      cnVariants = CNVariant.loadPlinkFile(proj.DATA_DIRECTORY.getValue(false, true)
                                           + COMBINED_TRIOS[fileType], false);
    } else {
      cnVariants = CNVariant.loadPlinkFile(proj.PROJECT_DIRECTORY.getValue() + cnvFile, false);
      // rawRegionFrequency = determinRegionFrequencies(cnVariants, proj.getLog());
    }
    parseTrios(cnVariants, trios, proj.getSampleData(0, false),
               (cnvFile == null ? proj.PROJECT_DIRECTORY.getValue() + output
                                  + COMBINED_TRIOS[fileType]
                                : proj.PROJECT_DIRECTORY.getValue() + output + ext.rootOf(cnvFile)),
               log);
    processTrios(proj, trios, rawRegionFrequency, numThreads);
    summarizeTrios(proj, trios, (cnvFile == null ? output + COMBINED_TRIOS[fileType]
                                                 : output + ext.rootOf(cnvFile)));

  }

  private static Hashtable<String, Integer> determinRegionFrequencies(CNVariant[] cnVariants,
                                                                      Logger log) {
    Hashtable<String, Integer> regionFrequency = new Hashtable<String, Integer>();
    ArrayList<String> uniqstmp = new ArrayList<String>();
    log.report(ext.getTime() + " Info - building population cnv frequency map for "
               + cnVariants.length + " cnvs, this may take some time...");
    for (int i = 0; i < cnVariants.length; i++) {
      if (!regionFrequency.containsKey(cnVariants[i].getUCSClocation())) {
        regionFrequency.put(cnVariants[i].getUCSClocation(), 0);
        uniqstmp.add(cnVariants[i].getUCSClocation());
      }
    }
    String[] uniq = uniqstmp.toArray(new String[uniqstmp.size()]);
    for (String element : uniq) {
      for (CNVariant cnVariant : cnVariants) {
        Segment curSegment = new Segment(element);
        if (closeAndOverlaps(curSegment, cnVariant)) {
          regionFrequency.put(element, regionFrequency.get(element) + 1);
        }
      }
    }
    for (CNVariant cnVariant : cnVariants) {
      if (regionFrequency.containsKey(cnVariant.getUCSClocation())
          && regionFrequency.get(cnVariant.getUCSClocation()) == 0) {
        regionFrequency.remove(cnVariant.getUCSClocation());// shrink it
      }
    }
    log.report(ext.getTime() + " Info - finished building population cnv frequency map");

    return regionFrequency;

  }

  private static boolean closeAndOverlaps(Segment seg1, Segment seg2) {
    boolean closeAndOverlaps = false;
    double maxRatio = (double) Math.max(seg1.getSize(), seg2.getSize())
                      / Math.min(seg1.getSize(), seg2.getSize());
    if (seg1.significantOverlap(seg2) && maxRatio < 10) {
      closeAndOverlaps = true;
    }
    return closeAndOverlaps;
  }

  /**
   * Assigns offspring CNV calls to trios by matching through sampleData
   * <p>
   * Essentially parses the array of CNVariants and assigns to the appropriate trio
   * <p>
   * If FID/IIDs are not unique (correspond to more than one DNA, the CNV will be dropped due to the
   * ambiguity)
   *
   * @param cnVariants
   * @param trios
   * @param sampleData
   * @param log
   */

  // TODO, resolve the issue of non-unique FID/IIDs (-> might have to be done manually), this is a
  // problem if there are identical FID/IID combos, and the DNA happens to match. We will not catch
  // it.
  private static void parseTrios(CNVariant[] cnVariants, TrioQC[] trios, SampleData sampleData,
                                 String parentOutput, Logger log) {
    Hashtable<String, Integer> trioIndex = hashOffspringTrios(trios, sampleData);
    Hashtable<String, Integer> parentalIndex = hashParentalTrios(trios, sampleData);
    ArrayList<CNVariant> parentalCNVs = new ArrayList<CNVariant>();
    int count = 0;
    for (CNVariant cnVariant : cnVariants) {
      String key = cnVariant.getFamilyID() + "\t" + cnVariant.getIndividualID();
      if (trioIndex.containsKey(key)) {
        if (trios[trioIndex.get(key)].addCNVTrio(sampleData, cnVariant, log)) {
          count++;
        }
      } else if (parentalIndex.containsKey(key)) {
        parentalCNVs.add(cnVariant);
        log.reportError("Warning - the cnv " + cnVariant.toPlinkFormat()
                        + " was not defined as an offspring in the trio file \n Reporting parental cnvs to "
                        + parentOutput);
      } else {
        log.reportError("Warning - the cnv " + cnVariant.toPlinkFormat()
                        + " was not defined as an offspring or parent file in the trio file \n Skipping");

      }
    }
    log.report(ext.getTime() + " Found " + count + " offspring cnvs");
    if (count != cnVariants.length) {
      log.reportError("Warning - Due to duplicate FID/IID combos, mismatched ids, or cnvs that were not for an offspring, only "
                      + count + "/" + cnVariants.length + " cnvs are being analyzed");
    }
    try {
      PrintWriter writer = new PrintWriter(new FileWriter(parentOutput + BEAST_OUTPUT[4]));
      writer.println(Array.toStr(PLINK_CNV_HEADER));
      for (int i = 0; i < parentalCNVs.size(); i++) {
        writer.println(parentalCNVs.get(i).toPlinkFormat());
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + parentOutput + BEAST_OUTPUT[4]);
      log.reportException(e);
    }
    // BEAST_OUTPUT[4]
  }

  /**
   * To speed things up a bit, we send out the heavy beast lifting to other threads
   */
  private static void processTrios(Project proj, TrioQC[] trios,
                                   Hashtable<String, Integer> rawRegionFrequency, int numThreads) {
    MarkerSet markerSet = proj.getMarkerSet();
    int[][] indi = markerSet.getIndicesByChr();
    ExecutorService executor = Executors.newFixedThreadPool(numThreads);
    for (int i = 0; i < trios.length; i++) {
      Runnable worker = new WorkerThread(proj, proj.getSampleData(0, false), trios[i],
                                         markerSet.getChrs(), markerSet.getPositions(),
                                         Matrix.clone(indi), rawRegionFrequency, i);
      executor.execute(worker);
    }
    executor.shutdown();
    try {
      executor.awaitTermination(7, TimeUnit.DAYS);
    } catch (InterruptedException e) {
    }
    proj.getLog().report(ext.getTime() + "Info - completed beast score computation");
  }

  /**
   * WorkerThreads which processes each trio to compute beast scores and qc metrics
   */

  private static class WorkerThread implements Runnable {
    private final Project proj;
    private final TrioQC trio;
    private final byte[] chrs;
    private final int[] positions;
    private final int[][] indicesByChr;
    Hashtable<String, Integer> rawRegionFrequency;
    private final SampleData sampleData;
    // private Logger log;
    private final int threadID;

    public WorkerThread(Project proj, SampleData sampleData, TrioQC trio, byte[] chr,
                        int[] positions, int[][] indicesByChr,
                        Hashtable<String, Integer> rawRegionFrequency, int threadID) {
      super();
      this.proj = proj;
      this.trio = trio;
      chrs = chr;
      this.positions = positions;
      this.indicesByChr = indicesByChr;
      this.sampleData = sampleData;
      this.rawRegionFrequency = rawRegionFrequency;
      this.threadID = threadID;
    }

    @Override
    public void run() {
      Logger log = proj.getLog();

      log.report(ext.getTime() + " Info - trio (" + (threadID + 1) + ") " + trio.getTrio());
      trio.computeTrioBeast(proj, sampleData, chrs, positions, indicesByChr, rawRegionFrequency);
    }
  }

  /**
   * Creates two output files, one with all the qc metrics, and another with just the list of
   * filtered trio calls
   *
   * @param proj
   * @param trios
   * @param fileType
   * @param log
   */
  private static void summarizeTrios(Project proj, TrioQC[] trios, String fileType) {
    PrintWriter writerFullSummary = Files.getAppropriateWriter(proj.PROJECT_DIRECTORY.getValue()
                                                               + fileType + BEAST_OUTPUT[0]);
    writerFullSummary.print(Array.toStr(CNVariant.PLINK_CNV_HEADER) + "\t"
                            + Array.toStr(OUTPUT_HEADER) + "\n");
    for (TrioQC trio : trios) {
      if (trio.hasSummary()) {
        writerFullSummary.println(Array.toStr(trio.getFullSummary(), "\n"));
      } else {
        proj.getLog().report("Info - not reporting " + trio.getTrio()
                             + ", did not have any offspring cnvs...");
      }
    }
    writerFullSummary.close();
  }

  /**
   * Assign the offspring FID/IIDs to a hash so we know which trio to add a particular CNVariant to
   *
   * @param trios
   * @param sampleData
   * @return
   */
  private static Hashtable<String, Integer> hashOffspringTrios(TrioQC[] trios,
                                                               SampleData sampleData) {
    Hashtable<String, Integer> trioIndex = new Hashtable<String, Integer>();
    for (int i = 0; i < trios.length; i++) {
      String key = trios[i].getOffspringFIDIID();
      trioIndex.put(key, i);
    }
    return trioIndex;
  }

  private static Hashtable<String, Integer> hashParentalTrios(TrioQC[] trios,
                                                              SampleData sampleData) {
    Hashtable<String, Integer> trioIndex = new Hashtable<String, Integer>();
    for (int i = 0; i < trios.length; i++) {
      trioIndex.put(sampleData.lookup(trios[i].getFADNA())[1], i);
      trioIndex.put(sampleData.lookup(trios[i].getMODNA())[1], i);
    }
    return trioIndex;
  }

  /**
   * Helper class to hold trio related information, and handle the computation/QC/summarization
   *
   */
  public static class TrioQC {
    private static final String[] PED_TRIO_HEADER = {"fId", "iId", "faId", "moId", "iDna", "faDna",
                                                     "moDna"};
    private final String fID;
    private final String iID;
    private final String faID;
    private final String moID;
    private final String iDNA;
    private final String faDNA;
    private final String moDNA;
    private final ArrayList<CNVariant> ICNV;
    private final ArrayList<String> fullSummary;
    private final ArrayList<String> filteredCNVS;

    public TrioQC(String fID, String iID, String faID, String moID, String iDNA, String faDNA,
                  String moDNA) {
      super();
      this.fID = fID;
      this.iID = iID;
      this.faID = faID;
      this.moID = moID;
      this.iDNA = iDNA;
      this.faDNA = faDNA;
      this.moDNA = moDNA;
      ICNV = new ArrayList<CNVariant>();
      fullSummary = new ArrayList<String>();
      filteredCNVS = new ArrayList<String>();
    }

    public String getFaID() {
      return faID;
    }

    public String getMoID() {
      return moID;
    }

    public String getIDNA() {
      return iDNA;
    }

    public String getFADNA() {
      return faDNA;
    }

    public String getMODNA() {
      return moDNA;
    }

    public String getOffspringFIDIID() {
      return fID + "\t" + iID;
    }

    public String getTrio() {
      return getIDNA() + "\t" + getFADNA() + "\t" + getMODNA();
    }

    /**
     * @return the full summary for all calls passed to the trio
     */
    public String[] getFullSummary() {
      return fullSummary.toArray(new String[fullSummary.size()]);
    }

    /**
     * @return the filtered cnv calls
     */
    public String[] getFilteredCNVS() {
      return filteredCNVS.toArray(new String[filteredCNVS.size()]);
    }

    public boolean hasSummary() {
      return fullSummary.size() > 0;
    }

    public boolean hasFilteredCNVS() {
      return filteredCNVS.size() > 0;
    }

    /***
     * The main event, grabs LRR_SD and BeastScore and few others. Summarizes QC across trio... max
     * LRR_SD across the trio and minimum beast score (height) distance from offspring to parent
     * Writes to two files: writerFullSummary the full summary of the beast-> writerFullSummary
     * writerFiltered a new filtered cnv file with passing calls, and calls assigned to each of the
     * parents for ease of vis in comp plot->trailer
     */
    public void computeTrioBeast(Project proj, SampleData sampleData, byte[] chrs, int[] positions,
                                 int[][] indicesByChr,
                                 Hashtable<String, Integer> rawRegionFrequency) {
      CNVariant[] analyzeCNV = ICNV.toArray(new CNVariant[ICNV.size()]);
      Logger log = proj.getLog();

      if (analyzeCNV.length == 0) {
        log.reportError("Warning - trio " + getTrio() + " does not have any CNVs defined");
      } else {
        int[][] targetIndices = getCNVIndices(chrs, positions, analyzeCNV, log);
        Sample iSamp = proj.getFullSampleFromRandomAccessFile(iDNA);
        Sample faSamp = proj.getFullSampleFromRandomAccessFile(faDNA);
        Sample moSamp = proj.getFullSampleFromRandomAccessFile(moDNA);

        float iBAF1585SD = getBAF1585SD(iSamp.getBAFs(), log);
        float faBAF1585SD = getBAF1585SD(faSamp.getBAFs(), log);
        float moBAF1585SD = getBAF1585SD(moSamp.getBAFs(), log);

        float[] iLrr = iSamp.getLRRs();
        float[] faLrr = faSamp.getLRRs();
        float[] moLrr = moSamp.getLRRs();

        final int end = Bytes.indexOf(chrs, (byte) 23);
        float iStDev = Array.stdev(Array.subArray(iLrr, 0, end), true);
        float faStDev = Array.stdev(Array.subArray(faLrr, 0, end), true);
        float moStDev = Array.stdev(Array.subArray(moLrr, 0, end), true);

        BeastScore iBeast = getBeast(iLrr, indicesByChr, targetIndices, log);
        BeastScore faBeast = getBeast(faLrr, indicesByChr, targetIndices, log);
        BeastScore moBeast = getBeast(moLrr, indicesByChr, targetIndices, log);

        for (int i = 0; i < analyzeCNV.length; i++) {
          // { "NUM_MARKERS", "IDNA", "ILRR_SD", "IBAF1585_SD", "IHEIGHT", "INumbProbes", "IBEAST",
          // "FADNA", "FALRR_SD", "FABAF1585_SD", "FAHEIGHT", "FANumProbes", "FABEAST", "MODNA",
          // "MOLRR_SD", "MOBAF1585_SD", "MOHEIGHT", "MONumProbes", "MOBEAST", "MinBeastDiff",
          // "MaxTrioLRR_SD", "MaxTrioBAF1585_SD", "NumCalls", "DenovoLength(BP)",
          // "DenovoMedianLRR", "NumRawOverlappingCNVs" };
          String fSum = "";
          fSum += analyzeCNV[i].toPlinkFormat();
          fSum += "\t" + iDNA + "\t" + faDNA + "\t" + moDNA;
          fSum += "\t" + iStDev + "\t" + faStDev + "\t" + moStDev;
          fSum += "\t" + iBAF1585SD + "\t" + faBAF1585SD + "\t" + moBAF1585SD;
          fSum += "\t" + iBeast.getBeastHeights()[i] + "\t" + faBeast.getBeastHeights()[i] + "\t"
                  + moBeast.getBeastHeights()[i];
          fSum += "\t" + getregionMedianLRR(iLrr, targetIndices[i]) + "\t"
                  + getregionMedianLRR(faLrr, targetIndices[i]) + "\t"
                  + getregionMedianLRR(moLrr, targetIndices[i]) + "\t" + analyzeCNV.length;
          fSum += "\t" + iBeast.getBeastLengths()[i] + "\t" + faBeast.getBeastLengths()[i] + "\t"
                  + moBeast.getBeastLengths()[i];
          fSum += "\t" + iBeast.getBeastScores()[i] + "\t" + faBeast.getBeastScores()[i] + "\t"
                  + moBeast.getBeastScores()[i];
          fSum += "\t"
                  + (rawRegionFrequency.containsKey(analyzeCNV[i].getUCSClocation()) ? rawRegionFrequency.get(analyzeCNV[i].getUCSClocation())
                                                                                     : "0");
          fullSummary.add(fSum);
        }
      }
    }

    private static float getregionMedianLRR(float[] lrrs, int[] indices) {
      ArrayList<Float> tmp = new ArrayList<Float>(indices.length);
      float median = Float.NaN;
      if (indices != null) {
        for (int i = 0; i < indices.length; i++) {
          if (!Float.isNaN(lrrs[indices[i]])) {
            tmp.add(lrrs[indices[i]]);
          }
        }
        if (tmp.size() > 0) {
          median = (float) Array.median(Doubles.toArray(tmp));
        }
      }
      return median;
    }

    private BeastScore getBeast(float[] lrrs, int[][] indicesByChr, int[][] targetIndices,
                                Logger log) {
      BeastScore iBeast = new BeastScore(lrrs, indicesByChr, targetIndices, log);
      iBeast.computeBeastScores(BeastScore.DEFAULT_ALPHA);
      return iBeast;
    }

    private static float getBAF1585SD(float[] bafs, Logger log) {
      for (int j = 0; j < bafs.length; j++) {
        if (bafs[j] < 0.15 || bafs[j] > 0.85) {
          bafs[j] = Float.NaN;
        }
      }
      return Array.stdev(bafs, true);
    }

    /**
     * Add the variant if lookup matches
     */
    private boolean addCNVTrio(SampleData sampleData, CNVariant cnVariant, Logger log) {
      if (checkValid(sampleData, cnVariant)) {
        ICNV.add(cnVariant);
        return true;
      } else {
        log.reportError("\nError - The current offspring has DNA: " + iDNA + ", FID: " + fID
                        + ", and IID: " + iID + ", the current cnv call has FID: "
                        + cnVariant.getFamilyID() + " and IID: " + cnVariant.getIndividualID());
        log.reportError("The lookup of the cnv FID/IID combo returned a mismatched DNA: "
                        + sampleData.lookup(cnVariant.getFamilyID() + "\t"
                                            + cnVariant.getIndividualID())[0]);
        log.reportError("This can happen if there are identical FID/IID combos in the file coorresponding to different DNA IDs. If so, please make the FID/IID combos unique");
        log.reportError("Skipping CNV " + cnVariant.toPlinkFormat()
                        + " due to unknown source DNA\n");
        return false;
      }
    }

    /**
     * Make sure the FID/IID lookup returns the correct DNA
     *
     * @param sampleData
     * @param cnVariant
     * @return
     */
    private boolean checkValid(SampleData sampleData, CNVariant cnVariant) {
      String lookup = sampleData.lookup(cnVariant.getFamilyID() + "\t"
                                        + cnVariant.getIndividualID())[0];
      String FIDIID = fID + "\t" + iID;
      return iDNA.equals(lookup)
             && FIDIID.equals(cnVariant.getFamilyID() + "\t" + cnVariant.getIndividualID());
    }

    /**
     * @param proj
     * @param trioFile should have header cnvTrio.PED_TRIO_HEADER
     * @param log
     * @return
     */
    private static TrioQC[] loadTrios(Project proj, String trioFile) {
      ArrayList<TrioQC> alTrio = new ArrayList<TrioQC>();
      int trioIndex = 0;
      Logger log = proj.getLog();

      try {
        BufferedReader reader =
                              Files.getReader(proj.DATA_DIRECTORY.getValue(false, true) + trioFile,
                                              false, true, false);
        String[] line = reader.readLine().trim().split(SPLITS[0]);
        int[] indices = ext.indexFactors(PED_TRIO_HEADER, line, true, true);
        while (reader.ready()) {
          line = reader.readLine().trim().split(SPLITS[0]);
          alTrio.add(new TrioQC(line[indices[0]], line[indices[1]], line[indices[2]],
                                line[indices[3]], line[indices[4]], line[indices[5]],
                                line[indices[6]]));
          trioIndex++;
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        log.reportError("Error: file \"" + proj.PROJECT_DIRECTORY.getValue() + trioFile
                        + "\" not found in current directory");
      } catch (IOException ioe) {
        log.reportError("Error reading file \"" + proj.PROJECT_DIRECTORY.getValue() + trioFile
                        + "\"");
      }
      if (trioIndex == 0) {
        log.reportError("Error - did not find any trios in "
                        + proj.DATA_DIRECTORY.getValue(false, true) + trioFile);
      }
      log.report(ext.getTime() + " Info - found "
                 + (trioIndex > 1 ? trioIndex + " trios" : trioIndex + " trio"));
      return alTrio.toArray(new TrioQC[alTrio.size()]);
    }
  }

  /**
   * Get the indices for the markers contained in the cnv
   *
   * @param chr
   * @param positions
   * @param cnVariants
   * @param log
   * @return
   */
  private static int[][] getCNVIndices(byte[] chr, int[] positions, CNVariant[] cnVariants,
                                       Logger log) {
    int[][] indices = new int[cnVariants.length][];
    for (int i = 0; i < cnVariants.length; i++) {
      indices[i] = BeastScore.getCNVMarkerIndices(chr, positions, cnVariants[i], log);
    }
    return indices;
  }

  /**
   * if offSpringHeight and parentHeight are different signs, return 0 for parentHeight Note, checks
   * for a real cnv should be done prior to using this function
   *
   * @param offSpringHeight
   * @param parentHeight
   * @return
   */
  private static double checkParentOppositeHeight(double offSpringHeight, double parentHeight) {
    if (offSpringHeight * parentHeight >= 0) {
      return parentHeight;
    } else {
      return 0;
    }
  }

  public static void parse(Project proj, String trioResultsFile, CNVTrioFilter trioFilter,
                           boolean excludeFromSampleData, String ouput) {
    SampleData sampleData = proj.getSampleData(0, false);
    cnvTrio[] rawCNVTrios = loadTrios(proj.PROJECT_DIRECTORY.getValue() + trioResultsFile,
                                      proj.getLog());
    ArrayList<Double> tmpMinBestDiffs = new ArrayList<Double>();
    ArrayList<cnvTrio> tmpFilteredCnvTrios = new ArrayList<cnvTrio>();
    for (cnvTrio rawCNVTrio : rawCNVTrios) {
      CNVFilterPass filterPass = trioFilter.getCNVTrioFilterPass(rawCNVTrio);
      proj.getLog().report(filterPass.getReasonNotPassing() + "\t"
                           + rawCNVTrio.getFullSummary(trioFilter.getHGBuild()));
      if (filterPass.passedFilter()) {

        tmpFilteredCnvTrios.add(rawCNVTrio);
        tmpMinBestDiffs.add(rawCNVTrio.getMinBeastHeightDifference());
      }

    }
    double[] beastDiffsToSort = Doubles.toArray(tmpMinBestDiffs);
    cnvTrio[] filteredCNVTrios =
                               tmpFilteredCnvTrios.toArray(new cnvTrio[tmpFilteredCnvTrios.size()]);
    int[] sorted = org.genvisis.common.Sort.quicksort(beastDiffsToSort, 1);
    try {
      PrintWriter writerSummary = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()
                                                                 + ouput + COMBINED_TRIOS[2]));
      PrintWriter writerCNV = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()
                                                             + ouput + COMBINED_TRIOS[2]
                                                             + COMBINED_TRIOS[3]));
      // TODO potential bug? gets first CNV filename and writes to file
      // PrintWriter writerList = new PrintWriter(new
      // FileWriter(proj.getFilename(proj.INDIVIDUAL_CNV_LIST_FILENAMES)));
      PrintWriter writerList =
                             new PrintWriter(new FileWriter(proj.INDIVIDUAL_CNV_LIST_FILENAMES.getValue()[0]));
      // PrintWriter writerRegion = new PrintWriter(new
      // FileWriter(proj.getFilename(proj.REGION_LIST_FILENAMES)));
      PrintWriter writerRegion =
                               new PrintWriter(new FileWriter(proj.REGION_LIST_FILENAMES.getValue()[0]));
      writerSummary.println(Array.toStr(PLINK_CNV_HEADER) + "\t"
                            + Array.toStr(FILTERED_OUTPUT_HEADER));
      writerCNV.println(Array.toStr(PLINK_CNV_HEADER));
      Hashtable<String, String> track = new Hashtable<String, String>();

      for (int i = 0; i < filteredCNVTrios.length; i++) {
        cnvTrio currentTrio = filteredCNVTrios[sorted[i]];
        if (!excludeFromSampleData || !trioHasExcluded(currentTrio, sampleData)) {
          writerSummary.println(currentTrio.getFullSummary(trioFilter.getHGBuild()));
          writerCNV.println(currentTrio.toPlinkFormat());
          writerCNV.println(currentTrio.getCNV(sampleData.lookup(currentTrio.getFADNA())[1]));
          writerCNV.println(currentTrio.getCNV(sampleData.lookup(currentTrio.getMODNA())[1]));
          writerList.println(currentTrio.getListSummary());
          if (!track.containsKey(currentTrio.getUCSClocation())) {
            track.put(currentTrio.getUCSClocation(), "has");
            writerRegion.println(currentTrio.getUCSClocation());
          }
        }
      }
      writerSummary.close();
      writerCNV.close();
      writerList.close();
      writerRegion.close();
    } catch (Exception e) {
      proj.getLog().reportError("Error writing to " + proj.PROJECT_DIRECTORY.getValue() + ouput
                                + COMBINED_TRIOS[2]);
      proj.getLog().reportException(e);
    }

  }

  private static boolean trioHasExcluded(cnvTrio trio, SampleData sampleData) {
    if (!sampleData.hasExcludedIndividuals()) {
      return false;
    } else {
      return sampleData.individualShouldBeExcluded(trio.getIDNA())
             || sampleData.individualShouldBeExcluded(trio.getFADNA())
             || sampleData.individualShouldBeExcluded(trio.getMODNA());
    }
  }

  public static String[] getParserParams() {
    String[] params = new String[5];
    params[0] = "#To intialize the trio parser, provide the following arguments";
    params[1] = "#the full path to a project properties file";
    params[2] = COMMAND_PROJECT;
    params[3] = "# the path (relative to the project directory) for the trio results file";
    params[4] = COMMAND_TRIO_RESULTS;
    params = Array.concatAll(params, CNVTrioFilter.getDefaultCNVTrioParams());
    return params;
  }

  public static void fromParameters(String filename, Logger log) {
    Vector<String> params;
    params = Files.parseControlFile(filename, CNVTrioFilter.COMMAND_CNV_TRIO_CRF, getParserParams(),
                                    log);
    if (params != null) {
      params.add(COMMAND_PARSE);
      main(Array.toStringArray(params));
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = null;
    String trioFile = "PedigreeOfTrios.txt";
    String logFile = null;
    String filenameOfProblematicRegions = null;
    String cnvFile = null;
    String output = "trio_";
    String trioResultsFile = null;
    int fileType = 1;
    boolean parse = false;
    boolean excludeFromSampleData = true;
    int numThreads = 7;
    Project proj;
    CNVTrioFilter trioFilter;
    String usage = "\n";
    usage += "jlDev.cnvTrio requires 1 arguments\n";
    usage += "   (1) project file (i.e. " + COMMAND_PROJECT + filename + " (no default))\n";
    usage += "    TO compute QC metrics for trio calls:";
    usage += "   (1) a file listing pedgiree information for the trios (i.e. trioFile=" + trioFile
             + " (default))\n";
    usage += "   (2) a pre-defined cnvFile (i.e. cnvFile=" + cnvFile + " (no default))\n";
    usage += "   (3) output base name (i.e. " + COMMAND_OUTPUT + output + " (default))\n";
    usage += "   (4) number of threads to use for the analysis (i.e. numThreads=" + numThreads
             + " ( default))\n";
    usage += "   (5) log filename (i.e. logFile=" + logFile + " (default))\n";

    usage += "   OR: filter a results file";

    usage += "   (1) parse the results of a trio run (i.e " + COMMAND_PARSE + " (not the default))";
    usage += "   (2) the trio results file to parse (i.e " + COMMAND_TRIO_RESULTS + trioResultsFile
             + " (no default))";
    usage += "   (3) filter out cnvs in known problematicRegions (i.e. " + COMMAND_FILTER_FILE
             + filenameOfProblematicRegions + " (no default))\n";
    usage += "   (4) log filename (i.e. logFile=" + logFile + " (default))\n";

    usage += "   (5) the minimum height score difference between offspring and parents (i.e. "
             + CNVTrioFilter.COMMAND_MIN_BEAST_HEIGHT_DIFFERENCE
             + CNVTrioFilter.DEFAULT_MIN_BEAST_HEIGHT_DIFFERENCE + " (default))\n";
    usage += "   (6) the maximum height score for parents  (i.e. "
             + CNVTrioFilter.COMMAND_MAX_BEAST_HEIGHT_PARENTS
             + CNVTrioFilter.DEFAULT_MAX_BEAST_HEIGHT_PARENTS + " (default))\n";
    usage += "   (7) maximum log R Ratio standard deviation across the trio (i.e."
             + CNVTrioFilter.COMMAND_MAX_TRIO_LRR_SD + CNVTrioFilter.DEFAULT_MAX_TRIO_LRR_SD
             + " (default))\n";
    usage += "   (8) maximum BAF1585 SD (i.e. " + CNVTrioFilter.COMMAND_MAX_TRIO_1585_SD
             + CNVTrioFilter.DEFAULT_MAX_TRIO_1585_SD + " (default))\n";
    usage += "   (9) maximum number of denovo cnvs in a sample (i.e. "
             + CNVTrioFilter.COMMAND_MAX_NUM_CALLS + CNVTrioFilter.DEFAULT_MAX_NUM_CALLS
             + " (default))\n";
    usage += "   (10) minimum number of probes for a cnv call;  (i.e. "
             + CNVFilter.COMMAND_MIN_NUM_MARKERS + CNVFilter.DEFAULT_MIN_NUM_MARKERS
             + " (default))\n";
    usage += "   (11) minimum score for a cnv (i.e. " + CNVFilter.COMMAND_MIN_SCORE
             + CNVFilter.DEFAULT_MIN_SCORE_THRESHOLD + " (default))\n";
    usage += "   (12) the genome build (can be \"36\" or \"37\") (i.e. " + CNVFilter.COMMAND_BUILD
             + CNVFilter.DEFAULT_BUILD + " (default))\n";
    usage += "   (13) the minimum length in base pairs of a denovo cnv (i.e. "
             + CNVFilter.COMMAND_MIN_SIZE + CNVFilter.NO_FILTER_MIN_SIZE + " (default))\n";
    usage += "   (14) output base name (i.e." + COMMAND_OUTPUT + output + " (default))\n";
    usage += "   (15) exclude samples from sample data (i.e. -exclude (default))\n";
    usage += "   (16) if a common reference file is provided, keep common variants in (i.e. "
             + CNVFilter.COMMAND_COMMON_IN + CNVFilter.DEFAULT_COMMON_IN + "  (default))\n";
    usage += "   (17) break up centromeres (i.e. " + CNVFilter.COMMAND_BREAK_UP_CENTROMERES
             + CNVFilter.DEFAULT_BREAK_UP_CENTROMERES + " (default))\n";
    usage += "";

    if (ext.indexOfStr(COMMAND_PROJECT, args, true, false) >= 0) {
      proj =
           new Project(ext.parseStringArg(args[ext.indexOfStr(COMMAND_PROJECT, args, true, false)],
                                          ""),
                       logFile, false);
    } else {
      proj = new Project(filename, logFile, false);
    }
    if (ext.indexOfStr(COMMAND_PARSE, args, true, false) >= 0) {
      trioFilter = CNVTrioFilter.setupCNVTrioFilterFromArgs(proj, args, true, proj.getLog());
    } else {
      trioFilter = null;// don't need it to compute metrics
    }

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith(COMMAND_PROJECT)) {
        filename = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("trioFile=")) {
        trioFile = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("fileType=")) {
        fileType = ext.indexOfStr(arg.split("=")[1], TRIOS);
        numArgs--;
      } else if (arg.startsWith("logFile=")) {
        logFile = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("cnvFile=")) {
        cnvFile = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("numThreads=")) {
        numThreads = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith(COMMAND_OUTPUT)) {
        output = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith(COMMAND_TRIO_RESULTS)) {
        trioResultsFile = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith(COMMAND_PARSE)) {
        parse = true;
        numArgs--;
      } else if (trioFilter.isCommandLineFilterInEffect(arg)) {
        numArgs--;

      } else {
        proj.getLog().report("Invalid argument " + arg);
      }
    }

    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    if (parse) {
      parse(proj, trioResultsFile, trioFilter, excludeFromSampleData, output);
      System.exit(1);
    } else {
      computeMetrics(proj, trioFile, fileType, cnvFile, output, numThreads);

    }
  }
}

//
// fSum += "\t" + iBeast.getBeastHeights()[i] + "\t" + faBeast.getBeastHeights()[i] + "\t" +
// moBeast.getBeastHeights()[i];
//
// fSum += fSum += fSum +=
//
// fSum += analyzeCNV[i].toPlinkFormat() + "\t" + targetIndices[i].length;
// fSum += "\t" + iDNA + "\t" + iStDev + "\t" + iBAF1585SD + "\t" + iBeast.getSummaryAt(i);
// fSum += "\t" + faDNA + "\t" + faStDev + "\t" + faBAF1585SD + "\t" + faBeast.getSummaryAt(i);
// fSum += "\t" + moDNA + "\t" + moStDev + "\t" + moBAF1585SD + "\t" + moBeast.getSummaryAt(i);
// fSum += "\t" + minBeastDiff + "\t" + maxStDev + "\t" + maxBAF1585SD + "\t" + analyzeCNV.length +
// "\t" + analyzeCNV[i].getSize() + "\t" + getDenovoMedianLRR(iLrr, targetIndices[i]);
// fSum += "\t" + getQCString(faBeast.getBeastHeights()[i], moBeast.getBeastHeights()[i],
// minBeastDiff, maxStDev, maxBAF1585SD, analyzeCNV.length, analyzeCNV[i].getSize(), filter);
// fSum += "\t" + analyzeCNV[i].getUCSClocation() + "\t" +
// analyzeCNV[i].getUCSCLink(filter.getBuild() == 36 ? "hg18" : "hg19") + "\t" +
// filter.getFrequencyForRegion(analyzeCNV[i].getUCSClocation());
// fullSummary.add(fSum);
//
// if (checkAll(faBeast.getBeastHeights()[i], moBeast.getBeastHeights()[i], analyzeCNV.length,
// maxStDev, maxBAF1585SD, minBeastDiff, analyzeCNV[i].getSize(), filter)) {
// filteredCNVS.add(analyzeCNV[i].toPlinkFormat());
// filteredCNVS.add(getVariant(sampleData.lookup(faDNA)[1].split("\t")[0], faID,
// analyzeCNV[i]).toPlinkFormat());
// filteredCNVS.add(getVariant(sampleData.lookup(moDNA)[1].split("\t")[0], moID,
// analyzeCNV[i]).toPlinkFormat());
// }
// }
// }

// /**
// *
// * @return a string with 1s for good and 0s for bad
// */
// private static String getQCString(float faHeight, float moHeight, float minBeastDiff, float
// maxStDev, float maxBAF1585SD, int numCalls, int sizeBP, Filter filter) {
// String qcString = "";
// if (minBeastDiff >= filter.getMinBeastHeightDiffThreshold()) {
// qcString += "1\t";
// } else {
// qcString += "0\t";
// }
// if (faHeight <= filter.getMaxBeastHeightParents() && moHeight <=
// filter.getMaxBeastHeightParents()) {
// qcString += "1\t";
// } else {
// qcString += "0\t";
// }
// if (checkLrrSdThreshold(maxStDev, filter)) {
// qcString += "1\t";
// } else {
// qcString += "0\t";
// }
// if (checkBAF1585SD(maxStDev, filter)) {
// qcString += "1\t";
// } else {
// qcString += "0\t";
// }
// if (numCalls <= filter.getMaxNumCallsThreshold()) {
// qcString += "1\t";
// } else {
// qcString += "0\t";
// }
// if (checkSize(filter, sizeBP)) {
// qcString += "1\t";
// } else {
// qcString += "0\t";
// }
// if (checkAll(faHeight, moHeight, numCalls, maxStDev, maxBAF1585SD, minBeastDiff, sizeBP, filter))
// {
// qcString += "1";
// } else {
// qcString += "0";
// }
// return qcString;
// }
//
//
// private static boolean checkAll(float faHeight, float moHeight, int numCalls, float maxStDev,
// float maxBAF1585SD, float minBeastDiff, int sizeBP, Filter filter) {
// return faHeight <= filter.getMaxBeastHeightParents() && moHeight <=
// filter.getMaxBeastHeightParents() && minBeastDiff >= filter.getMinBeastHeightDiffThreshold() &&
// checkLrrSdThreshold(maxStDev, filter) && numCalls <= filter.getMaxNumCallsThreshold() &&
// checkBAF1585SD(maxBAF1585SD, filter) && checkSize(filter, sizeBP);
// }
//
//
//
//
// private static boolean checkSize(Filter filter, int sizeBP) {
// return sizeBP >= filter.getDenovoMinLength() && sizeBP <= filter.getDenovoMaxLength();
// }
//
// private static boolean checkBAF1585SD(float maxStDev, Filter filter) {
// return maxStDev <= filter.getMaxBAF1585SD();
// }
//
// private static boolean checkLrrSdThreshold(float maxStDev, Filter filter) {
// return maxStDev <= filter.getMaxLrrSdThreshold();
// }
//
// /**
// * If the parents have a height that exceeds the realCNVHeight (MaxBeastHeightParents), we defualt
// to 0; If the offspring's height does not surpass the realCNVHeight (MaxBeastHeightParents), we
// defualt to 0;
// * <p>
// * <p>
// * Also check if the parents have in opposite signed height as the offspring, if so, we set the
// parent's height to 0. Otherwise, we return the minimum distance between parents and offspring
// *
// * @param offSpringHeight
// * beast height of the offspring
// * @param faHeigt
// * father beast height
// * @param moHeight
// * mother beast height
// * @param realCNVHeight
// * definition of height for a cnv
// * @return
// */
// private static float minHeightDist(float offSpringHeight, float faHeigt, float moHeight, Filter
// filter) {
// float minHeightDist = 0;
// if (Math.abs(faHeigt) > filter.getMaxBeastHeightParents() || Math.abs(moHeight) >
// filter.getMaxBeastHeightParents() || Math.abs(offSpringHeight) <
// filter.getMaxBeastHeightParents()) {
// return minHeightDist;
// }
// faHeigt = checkParentOppositeHeight(offSpringHeight, faHeigt);
// moHeight = checkParentOppositeHeight(offSpringHeight, moHeight);
// return Math.min(Math.abs(offSpringHeight - faHeigt), Math.abs(offSpringHeight - moHeight));
// //}
//
// /**
// * A small helper class to facilitate passing thresholds around to anyone who wants em
// *
// */
// private static class Filter {
// private float minBeastHeightDiffThreshold;
// private float maxLrrSdThreshold;
// private float maxBeastHeightParents;
// private int maxNumCallsThreshold, numProbes, build, denovoMinLength, denovoMaxLength;
// private float maxBAF1585SD;
// private float score;
// private String filenameOfProblematicRegions;
// private Hashtable<String, Integer> rawRegionFrequency;
//
// /**
// * @param minBeastHeightDiffThreshold
// * minimum height difference between parents and offspring (also defines what height is a cnv)
// * @param maxLrrSdThreshold
// * max LRR_SD across the trio
// * @param maxBeastHeightParents
// * maximum height of the parents (i.e height of a real cnv)
// * @param maxNumCallsThreshold
// * maximum denovo calls allowed
// * @param numProbes
// * minimum number of probes for a call
// * @param maxBAF1585SD
// * maximum BAF_1585_SD across the trio
// */
// public Filter(float minBeastHeightDiffThreshold, float maxLrrSdThreshold, float
// maxBeastHeightParents, int maxNumCallsThreshold, int numProbes, float maxBAF1585SD, float score,
// int build, String filenameOfProblematicRegions, int denovoMinLength, int denovoMaxLength) {
// super();
// this.minBeastHeightDiffThreshold = minBeastHeightDiffThreshold;
// this.maxLrrSdThreshold = maxLrrSdThreshold;
// this.maxBeastHeightParents = maxBeastHeightParents;
// this.maxNumCallsThreshold = maxNumCallsThreshold;
// this.numProbes = numProbes;
// this.maxBAF1585SD = maxBAF1585SD;
// this.score = score;
// this.build = build;
// this.filenameOfProblematicRegions = filenameOfProblematicRegions;
// this.denovoMinLength = denovoMinLength;
// this.denovoMaxLength = denovoMaxLength;
//
// }
//
// public String getParamsUsed() {
// String params = "";
// params += "Minimum beast height difference= " + minBeastHeightDiffThreshold + "\n";
// params += "Maximum trio LRR SD = " + maxLrrSdThreshold + "\n";
// params += "Maximumn beast height for parents = " + maxBeastHeightParents + "\n";
// params += "Maximumn number of calls = " + maxNumCallsThreshold + "\n";
// params += "Minumum number of probes = " + numProbes + "\n";
// params += "Maximum trio BAF 15/85 SD = " + maxBAF1585SD + "\n";
// params += "Minumum Score = " + score + "\n";
// params += "Genome Build = " + build + "\n";
// params += "File name of probelamtic regions = " + filenameOfProblematicRegions + "\n";
// params += "Minimum denovo length = " + denovoMinLength + "\n";
// params += "Maximum denovo length = " + denovoMaxLength + "\n";
// return params;
// }
//
// public void setRawRegionFrequency(Hashtable<String, Integer> rawRegionFrequency) {
// this.rawRegionFrequency = rawRegionFrequency;
// }
//
// public int getBuild() {
// return build;
// }
//
// public String getFilenameOfProblematicRegions() {
// return filenameOfProblematicRegions;
// }
//
// public float getScore() {
// return score;
// }
//
//
//
//
//
//
// public int getNumProbes() {
// return numProbes;
// }
//
//
// // }
//
//
// private static boolean trioContainsExclude(SampleData sampleData, Trio trio) {
// boolean hasExcluded = false;
// if (sampleData.hasExcludedIndividuals()) {
// if (sampleData.individualShouldBeExcluded(trio.getFADNA())) {
// hasExcluded = true;
// } else if (sampleData.individualShouldBeExcluded(trio.getMODNA())) {
// hasExcluded = true;
// } else if (sampleData.individualShouldBeExcluded(trio.getIDNA())) {
// hasExcluded = true;
// }
// }
// return hasExcluded;
// }
//
// /**
// * @param proj
// * @param cnvFile
// * converts this cnv file to Region list, Individual CNV list, and a bed type file
// * @param log
// */
// private static void writeRegionList(Project proj, String cnvFile) {
// CNVariant[] cnvs = CNVariant.toCNVariantArray(CNVariant.loadPlinkFile(cnvFile, null, false));
// proj.getLog().report("Info - updating regions and list using " + cnvs.length + " cnvs");
// String list = proj.getFilename(proj.INDIVIDUAL_CNV_LIST_FILENAMES).replaceAll(";", "");
// String regionList = proj.getFilename(proj.REGION_LIST_FILENAMES).replaceAll(";", "");
// ArrayList<String> lists = new ArrayList<String>();
// ArrayList<String> bed = new ArrayList<String>();
// ArrayList<String> regions = new ArrayList<String>();
// SampleData sampleData = proj.getSampleData(0, false);
// Hashtable<String, String> track = new Hashtable<String, String>();
// for (int i = 0; i < cnvs.length; i++) {
// String[] id = sampleData.lookup(cnvs[i].getFamilyID() + "\t" + cnvs[i].getIndividualID());
// if (id != null) {
// String UCSC = cnvs[i].getUCSClocation();
// if (!track.containsKey(UCSC)) {
// bed.add("chr" + cnvs[i].getChr() + "\t" + cnvs[i].getStart() + "\t" + cnvs[i].getStop() + "\t" +
// cnvs[i].getUCSClocation());
// regions.add(UCSC);
// track.put(UCSC, UCSC);
// }
// lists.add(sampleData.lookup(cnvs[i].getFamilyID() + "\t" + cnvs[i].getIndividualID())[0] + "\t" +
// UCSC);
// } else {
// proj.getLog().report("Warning - could not look up the sample corresponding to " +
// cnvs[i].toPlinkFormat() + " in sample data, and will not be present in " + list + " or " +
// regionList);
// }
// }
// Files.writeList(lists.toArray(new String[lists.size()]), list);
// Files.writeList(regions.toArray(new String[regions.size()]), regionList);
// Files.writeList(bed.toArray(new String[bed.size()]), regionList + ".bed");
// }
