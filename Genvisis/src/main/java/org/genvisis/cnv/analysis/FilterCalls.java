package org.genvisis.cnv.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.UCSCtrack;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Aliases;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CNVFilter;
import org.genvisis.common.CNVFilter.CNVFilterPass;
import org.genvisis.common.Files;
import org.genvisis.common.GenomicPosition;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.GeneSet;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.filesys.SegmentLists;
import org.genvisis.filesys.SnpMarkerSet;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;
import com.google.common.collect.SetMultimap;
import com.google.common.collect.Sets;
import com.google.common.primitives.Doubles;

public class FilterCalls {

  public static final int DEFAULT_MIN_SIZE_KB = 0;
  public static final int DEFAULT_MIN_NUM_SNPS = 1;
  public static final double DEFAULT_MIN_SCORE = 10.0;
  public static final boolean DEFAULT_FILTER_REGIONS_FLAG = false;
  public static final String DEFAULT_PROBLEMATIC_REGIONS = "data/problematicRegions.dat";
  // public static final String DEFAULT_COMMON_CNP_REFERENCE = "data/polymorphic_CNPs.txt";
  // public static final String DEFAULT_COMMON_CNP_REFERENCE = "data/common_CNPs.txt";
  public static final String DEFAULT_COMMON_CNP_REFERENCE = "data/all_CNPs.txt";
  public static final String[] DEFAULT_REGION_DIRECTORIES = {"C:\\Documents and Settings\\npankrat\\My Documents\\CNV_PD\\",
                                                             "/home/npankrat/", "P:\\",
                                                             "/home/pankrat2/"};
  public static final int COMMON_IN = 1;
  public static final int COMMON_OUT = 2;
  public static final int COMMON_IGNORED = 3;
  public static final int DEFAULT_COMMON_IN_OUT_OR_IGNORED = COMMON_IGNORED;
  public static final boolean DEFAULT_BREAK_CENTROMERE = false;
  /**
   * Score/Probe thresholds for CNVStats, altering these will alter the number of columns in the
   * outputted stats file
   */
  private static final double[][] CNV_STATS_THRESHOLDS = new double[][] {{10, 10}, {10, 20}};

  public static final float DEFAULT_CLEAN_FACTOR = 0.2f;

  private static final class CNVComparator implements Comparator<CNVariant> {

    @Override
    public int compare(CNVariant o1, CNVariant o2) {
      if (o1.getStart() < o2.getStart()) {
        return -1;
      } else if (o1.getStart() > o2.getStart()) {
        return 1;
      } else {
        if (o1.getStop() < o2.getStop()) {
          return -1;
        } else if (o1.getStop() > o2.getStop()) {
          return 1;
        } else {
          return 0;
        }
      }
    }
  }

  final static CNVComparator cnvComparator = new CNVComparator();

  private static class CNVFilterNode {

    public static final int HET_DEL = 0;
    public static final int HOM_DEL = 1;
    public static final int HET_DUP = 2;
    public static final int HOM_DUP = 3;

    final CNVariant cnv;
    final int popCnt;
    final ArrayList<CNVariant> major = new ArrayList<>();
    final ArrayList<CNVariant> minor = new ArrayList<>();

    public CNVFilterNode(CNVariant myCNV, int pop) {
      cnv = myCNV;
      popCnt = pop;
    }

    public void addMajor(CNVariant cnv) {
      major.add(cnv);
    }

    public void addMinor(CNVariant cnv) {
      minor.add(cnv);
    }

    public int countType(int type, boolean lookMajor) {
      int lookingFor = -1;
      switch (type) {
        case HET_DEL:
          lookingFor = 1;
          break;
        case HOM_DEL:
          lookingFor = 0;
          break;
        case HET_DUP:
          lookingFor = 3;
          break;
        case HOM_DUP:
          lookingFor = 4;
          break;
        default:
          // TODO show error message
          return 0;
      }

      int cnt = 0;
      int ind = 0;
      ArrayList<CNVariant> pop = lookMajor ? major : minor;
      for (CNVariant cnv : pop) {
        System.out.println((ind++) + " " + cnv.getCN());
        // TODO output warning if CN > 4
        if (cnv.getCN() == lookingFor || (lookingFor == 4 && cnv.getCN() > lookingFor)) {
          cnt++;
        }
      }

      return cnt;
    }

    @Override
    public String toString() {
      // System.out.println("Major: " + major.size() + " of " + popCnt + " - " +
      // ext.formDeci((((double)(major.size()) / ((double)popCnt)) * 100.0), 3, true) + "%");
      // System.out.println("Minor: " + minor.size() + " of " + popCnt + " - " +
      // ext.formDeci((((double)(minor.size()) / ((double)popCnt)) * 100.0), 3, true) + "%");

      return cnv.toPlinkFormat() + "\t"
             + ext.formDeci((((double) (major.size()) / ((double) popCnt)) * 100.0), 3, true) + "\t"
             + ext.formDeci((((double) (minor.size()) / ((double) popCnt)) * 100.0), 3, true) + "\t"
             + "(" + countType(HOM_DEL, true) + "," + countType(HET_DEL, true) + ","
             + countType(HET_DUP, true) + "," + countType(HOM_DUP, true) + ")\t" + "("
             + countType(HOM_DEL, false) + "," + countType(HET_DEL, false) + ","
             + countType(HET_DUP, false) + "," + countType(HOM_DUP, false) + ")";
    }

  }

  /**
   * Write a file about the contents of a given CNV file.<br />
   * Output format:<br />
   * <br />
   * <code>|	SAMPLE/DNA	|	FID	|	IID	|	Exclude	|	LRRSD	|	#CNVs	|	#CNVs_c10p10	|	#CNVs_c20p10	|</code>
   * <br />
   * <br />
   * Columns can change depending on an internal array, CNV_STATS_THRESHOLDS, which define the
   * thresholds for the last few columns
   *
   * @param proj Project
   * @param dir Path to directory
   * @param filenameNoExt Filename (without extension) of CNV file
   * @throws IOException
   */
  public static void CNVStats(Project proj, String dir, String filenameNoExt) throws IOException {
    String qcFile, cnvFile, outputFile;
    PrintWriter writer;
    BufferedReader reader;
    SampleData sampleData;

    // find .cnv and .fam file from fileroot
    qcFile = (proj == null ? dir : proj.PROJECT_DIRECTORY.getValue()) + "Sample_QC.xln";
    cnvFile = dir + filenameNoExt + ".cnv";
    // famFile = dir + filenameNoExt + ".fam";

    outputFile = dir + filenameNoExt + "_CNVStats.xln";

    sampleData = proj == null ? null : proj.getSampleData(false);

    List<CNVariant> cnvList = CNVariant.loadPlinkFile(cnvFile, null, true);
    HashMap<String, ArrayList<CNVariant>[]> cnvMap = new HashMap<>();
    for (CNVariant cnv : cnvList) {
      ArrayList<CNVariant>[] indivLists = cnvMap.get(cnv.getFamilyID() + "\t"
                                                     + cnv.getIndividualID());
      if (indivLists == null) {
        indivLists = new ArrayList[CNV_STATS_THRESHOLDS.length + 1];
        for (int i = 0; i < CNV_STATS_THRESHOLDS.length + 1; i++) {
          indivLists[i] = new ArrayList<>();
        }
        cnvMap.put(cnv.getFamilyID() + "\t" + cnv.getIndividualID(), indivLists);
      }
      indivLists[0].add(cnv);
      for (int i = 0; i < CNV_STATS_THRESHOLDS.length; i++) {
        if (cnv.getScore() > CNV_STATS_THRESHOLDS[i][0]
            && cnv.getNumMarkers() > CNV_STATS_THRESHOLDS[i][1]) {
          indivLists[i + 1].add(cnv);
        }
      }
    }

    /*
     * OUTPUT | SAMPLE/DNA | FID | IID | Excluded | LRRSD | CNV COUNTS .... | | |
     */
    String header = "SAMPLE/DNA\tFID\tIID\tExclude\tLRRSD\t#CNVs";// +"#CNVs_c10p10" + "\t" +
                                                                  // "#CNVs_c20p10";
    for (double[] element : CNV_STATS_THRESHOLDS) {
      header = header + "\t#CNVs_c" + element[0] + "p" + element[1];
    }
    writer = new PrintWriter(outputFile);
    writer.println(header);
    reader = new BufferedReader(new FileReader(qcFile));
    String line, SID, FID, IID, LRRSD;
    boolean excluded;
    String[] cnts;
    String[] data;
    reader.readLine();
    while (reader.ready()) {
      line = reader.readLine();
      data = line.split("\t");
      SID = data[0];
      FID = data[1];
      IID = data[2];

      LRRSD = data[11];
      excluded = sampleData == null ? false : sampleData.individualShouldBeExcluded(SID);

      ArrayList<CNVariant>[] indivLists = cnvMap.get(FID + "\t" + IID);
      if (indivLists == null) {
        cnts = new String[] {".", ".", "."};
      } else {
        cnts = new String[indivLists.length];
        for (int i = 0; i < cnts.length; i++) {
          cnts[i] = indivLists[i].size() + "";
        }
      }

      writer.println(SID + "\t" + FID + "\t" + IID + "\t" + (excluded ? "1" : "0") + "\t" + LRRSD
                     + "\t" + cnts[0] + "\t" + cnts[1] + "\t" + cnts[2]);
    }
    reader.close();
    writer.flush();
    writer.close();

  }

  /**
   * Take a list of CNVs and search through other lists of CNVs, counting number of overlapping
   * (both major and minor overlap) CNVs with a minimum score and probe count
   *
   * @param cnvList List of CNVs for which to compile stats
   * @param cnvFiles Lists of CNV filenames in which to look for overlapping CNVs
   * @param outputFile Name of output file
   * @param score Minimum score threshold for comparison
   * @param probes Minimum probe-count threshold for comparison
   */
  public static void cnvStats(String cnvList, String[] cnvFiles, String outputFile, double score,
                              int probes, double overlapThreshold) {
    System.out.println(ext.getTime() + "] Loading Plink-formatted CNV files...");
    CNVariant[] srcCNVs = CNVariant.loadPlinkFile(cnvList);
    ArrayList<CNVariant> compCNVs = new ArrayList<>();
    HashSet<String> ids = new HashSet<>();
    for (String file : cnvFiles) {
      CNVariant[] cnvs = CNVariant.loadPlinkFile(file);
      for (CNVariant cnv : cnvs) {
        compCNVs.add(cnv);
        ids.add(cnv.getFamilyID() + "\t" + cnv.getIndividualID());
      }
    }

    ArrayList<CNVFilterNode> outputNodes = new ArrayList<>();
    ArrayList<String> majorMatches = new ArrayList<>();
    ArrayList<String> minorMatches = new ArrayList<>();

    System.out.println(ext.getTime() + "] Analyzing CNV overlap...");
    for (CNVariant cnv : srcCNVs) {
      CNVFilterNode cnvNode = new CNVFilterNode(cnv, ids.size());
      outputNodes.add(cnvNode);

      for (CNVariant comp : compCNVs) {
        if (cnv.getFamilyID().equals(comp.getFamilyID())
            && cnv.getIndividualID().equals(comp.getIndividualID())) {
          continue;
        }
        if (comp.getScore() < score) {
          continue;
        }
        if (/* (comp.getCN() == 0 && comp.getNumMarkers() < 3) || (comp.getCN() != 0 && */comp.getNumMarkers() < probes/*
                                                                                                                        * )
                                                                                                                        */) {
          continue;
        }
        if (comp.getCN() == 2) {
          continue;
        }
        int overlap = cnv.amountOfOverlapInBasepairs(comp);
        if (overlap == -1) {
          continue;
        }
        if (overlap >= ((cnv.getSize()) * overlapThreshold)) {
          cnvNode.addMajor(comp);
          majorMatches.add(comp.getFamilyID() + "\t" + comp.getIndividualID() + "\t" + cnv.getChr()
                           + "\t" + cnv.getStart() + "\t" + cnv.getStop() + "\t" + comp.getCN());
        } else {
          cnvNode.addMinor(comp);
          minorMatches.add(comp.getFamilyID() + "\t" + comp.getIndividualID() + "\t" + cnv.getChr()
                           + "\t" + cnv.getStart() + "\t" + cnv.getStop() + "\t" + comp.getCN());
        }
      }

    }

    System.out.println(ext.getTime() + "] Writing results to \"" + outputFile + "\"...");

    String header = ArrayUtils.toStr(CNVariant.PLINK_CNV_HEADER, "\t")
                    + "\t%major\t%minor\tStats(M)\tStats(m)";

    PrintWriter writer;
    try {
      writer = Files.openAppropriateWriter(outputFile);
      writer.println(header);
      for (CNVFilterNode node : outputNodes) {
        writer.println(node.toString());
      }
      writer.flush();
      writer.close();

      writer = Files.openAppropriateWriter(ext.rootOf(outputFile, false) + ".major");
      writer.println("FID\tIID\tCHR\tBP1\tBP2\tCN");
      for (String str : majorMatches) {
        writer.println(str);
      }
      writer.flush();
      writer.close();

      writer = Files.openAppropriateWriter(ext.rootOf(outputFile, false) + ".minor");
      writer.println("FID\tIID\tCHR\tBP1\tBP2\tCN");
      for (String str : minorMatches) {
        writer.println(str);
      }
      writer.flush();
      writer.close();

    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }

  private static class MergedCNVariant {

    int markerStart;
    int markerStop;
    CNVariant cnv;
    final ArrayList<MergedCNVariant> originalCNVs = new ArrayList<>();
    double medianLRR;
    double stdevLRR;

    public MergedCNVariant(int start, int stop, CNVariant cnvariant) {
      markerStart = start;
      markerStop = stop;
      cnv = cnvariant;
    }

    public MergedCNVariant addOriginalCNV(MergedCNVariant cnv) {
      if (cnv != null) {
        originalCNVs.add(cnv);
      }
      return this;
    }

  }

  public static CNVariant[] mergeCNVsInMemory(Project proj, CNVariant[] inCNVs,
                                              double distanceQuotient, boolean cnvsAsPositions) {
    Logger log = proj.getLog();
    int initialCount = inCNVs.length;
    ArrayList<CNVariant> mergedCNVs = new ArrayList<>();

    HashMap<String, HashMap<Byte, ArrayList<CNVariant>>> indivChrCNVMap = new HashMap<>();

    log.report(ext.getTime() + "] Organizing CNVs in memory...");
    for (CNVariant cnv : inCNVs) {
      HashMap<Byte, ArrayList<CNVariant>> chrCNVMap = indivChrCNVMap.get(cnv.getFamilyID() + "\t"
                                                                         + cnv.getIndividualID());
      if (chrCNVMap == null) {
        chrCNVMap = new HashMap<>();
        indivChrCNVMap.put(cnv.getFamilyID() + "\t" + cnv.getIndividualID(), chrCNVMap);
      }
      ArrayList<CNVariant> cnvList = chrCNVMap.get(cnv.getChr());
      if (cnvList == null) {
        cnvList = new ArrayList<>();
        chrCNVMap.put(cnv.getChr(), cnvList);
      }
      cnvList.add(cnv);
    }

    MarkerSetInfo markerSet = proj.getMarkerSet();
    int[][] positions = null;
    if (cnvsAsPositions) {
      LocusSet<CNVariant> tmp = new LocusSet<CNVariant>(inCNVs, true, proj.getLog()) {

        /**
         *
         */
        private static final long serialVersionUID = 1L;

      };
      positions = tmp.getStartsAndStopsByChromosome();
    } else {

      positions = markerSet.getPositionsByChr();
    }
    Hashtable<String, String> droppedMarkerNames = proj.getFilteredHash();
    String[] markerNames = markerSet.getMarkerNames();
    SampleData sampleData = proj.getSampleData(false);

    log.report(ext.getTime() + "] Cleaning CNVs...");
    int size = indivChrCNVMap.size();
    int cnt = 0;
    int step = size / 10;
    int cnvCount = 0;
    for (java.util.Map.Entry<String, HashMap<Byte, ArrayList<CNVariant>>> entry : indivChrCNVMap.entrySet()) {
      if (cnt > 0 && cnt % step == 0) {
        log.report(ext.getTime() + "] \t" + ((cnt / step) * 10) + "% complete");
      }
      String fidiid = entry.getKey();

      String ind;
      float[] lrrs;
      try {
        ind = sampleData.lookup(fidiid)[0];
      } catch (NullPointerException npe) {
        log.reportError("Error - could not look up the sample " + fidiid
                        + " in the sample data file " + proj.SAMPLE_DATA_FILENAME.getValue()
                        + ", cannot load sample to compute beast score");
        log.reportError("Error - please ensure that the sample names correspond to the varaints being processed with FID="
                        + fidiid.split("\t")[0] + " and IID=" + fidiid.split("\t")[1]);
        continue;
      }
      try {
        lrrs = proj.getFullSampleFromRandomAccessFile(ind).getLRRs();
      } catch (NullPointerException npe) {
        log.reportError("Error - could not load data for the sample " + ind + "\t" + fidiid
                        + ", please ensure samples have been parsed prior to computing beast score");
        log.report("Skipping beast score for sample " + ind + "\t" + fidiid);
        continue;
      }

      for (java.util.Map.Entry<Byte, ArrayList<CNVariant>> cnvLists : entry.getValue().entrySet()) {
        byte chr = cnvLists.getKey();

        if (cnvLists.getValue().size() == 0) {
          continue;
        } else if (cnvLists.getValue().size() == 1) {
          mergedCNVs.add(cnvLists.getValue().get(0));
          continue;
        }

        ArrayList<CNVariant> cnvList = cnvLists.getValue();
        Collections.sort(cnvList, cnvComparator);

        LinkedList<MergedCNVariant> tempChromo = new LinkedList<>();

        // create objects for all CNVs while also setting start/end marker indices, accounting for
        // dropped markers
        for (int i = 0; i < cnvList.size(); i++) {
          CNVariant cnv = cnvList.get(i);
          int firstSNPIndex = Arrays.binarySearch(positions[chr], cnv.getStart());
          int lastSNPIndex = Arrays.binarySearch(positions[chr], cnv.getStop());

          if (firstSNPIndex < 0 || lastSNPIndex < 0) {
            log.reportTimeWarning("Cannot merge CNV with bounds outside known start/stop: "
                                  + cnv.getStart() + "-" + cnv.getStop());
            continue;
          }

          // account for dropped markers
          if (droppedMarkerNames.contains(markerNames[firstSNPIndex])) {
            int f1 = firstSNPIndex, f2 = firstSNPIndex;
            while (f1 > (i == 0 ? 0 : tempChromo.get(i - 1).markerStop)
                   && droppedMarkerNames.contains(markerNames[f1])) {
              f1--;
            }
            while (f2 < lastSNPIndex - 1 && droppedMarkerNames.contains(markerNames[f2])) {
              f2++;
            }
            firstSNPIndex = ((firstSNPIndex - f1) < (f2 - firstSNPIndex) ? f1 : f2);
          }
          if (droppedMarkerNames.contains(markerNames[lastSNPIndex])) {
            int s1 = firstSNPIndex, s2 = lastSNPIndex;
            while (s1 > firstSNPIndex && droppedMarkerNames.contains(markerNames[s1])) {
              s1--;
            }
            while (s2 < markerNames.length - 1 && droppedMarkerNames.contains(markerNames[s2])) {
              s2++;
            }

            lastSNPIndex = ((lastSNPIndex - s1) < (s2 - lastSNPIndex) ? s1 : s2);
          }

          MergedCNVariant orig = new MergedCNVariant(firstSNPIndex, lastSNPIndex, cnvList.get(i));
          orig.addOriginalCNV(orig);
          tempChromo.add(orig);
        }

        LinkedList<MergedCNVariant> chromo = new LinkedList<>();

        // create full list: CNVs and objects for the spaces between CNVs, also accounting for
        // dropped markers
        chromo.add(tempChromo.get(0));
        for (int i = 1; i < tempChromo.size(); i++) {
          MergedCNVariant cnv1 = tempChromo.get(i - 1);
          MergedCNVariant cnv2 = tempChromo.get(i);

          if (cnv2.markerStart - cnv1.markerStop <= 1) {
            chromo.add(cnv2);
            continue;
          }

          int m1 = cnv1.markerStop + 1;
          int m2 = cnv2.markerStart - 1;

          if (droppedMarkerNames.contains(markerNames[m1])) {
            int f1 = m1;
            while (f1 < m2 - 1 && droppedMarkerNames.contains(markerNames[f1])) {
              f1++;
            }
            m1 = f1;
          }
          if (droppedMarkerNames.contains(markerNames[m2])) {
            int s1 = m2;
            while (s1 > m1 + 1 && droppedMarkerNames.contains(markerNames[s1])) {
              s1--;
            }
            m2 = s1;
          }

          chromo.add(new MergedCNVariant(m1, m2, null));
          chromo.add(cnv2);
        }
        tempChromo = null;

        // iterate through full array and compute median and stdev of LRRs
        for (int i = 0; i < chromo.size(); i++) {
          MergedCNVariant cnv = chromo.get(i);
          setLRRMedStdDev(cnv, droppedMarkerNames, markerNames, lrrs);
        }

        HashMap<Integer, MergedCNVariant> removed = new HashMap<>();
        TreeSet<Integer> indexList = new TreeSet<>();
        boolean done = false;
        while (!done) {
          int index = 0;
          MergedCNVariant cnv = chromo.get(0);
          // find index of first ICS (inter-CNV space); should never be index 0 or max
          for (index = 1; index < chromo.size() - 1 && cnv.cnv != null; index++) {
            cnv = chromo.get(index);
          }

          if (cnv.cnv != null) {
            // iterated through entire list and no ICSs left!
            done = true;
            break;
          }
          index--;

          MergedCNVariant actualCNV1 = chromo.get(index - 1);
          MergedCNVariant actualCNV2 = chromo.get(index + 1);

          int bpSize = countBP(actualCNV1) + countBP(actualCNV2);

          // same CN
          boolean sameCN = actualCNV1.cnv.getCN() == actualCNV2.cnv.getCN();
          // less than distanceQuotient percent space vs # of markers
          float mkQ = (float) ((cnv.markerStop - cnv.markerStart + 1))
                      / (float) ((actualCNV1.cnv.getNumMarkers() + actualCNV2.cnv.getNumMarkers()));
          boolean markerQuotient = mkQ < distanceQuotient;
          // less than 100% of total called base pairs
          float bpQ = (actualCNV2.cnv.getStart() - actualCNV1.cnv.getStop() + 1) / (float) bpSize;
          boolean basepairQuotient = bpQ < 1;

          if (!sameCN || !markerQuotient || !basepairQuotient) {
            // remove ICS placeholder
            removed.put(index, chromo.remove(index));
            indexList.add(index);
            continue;
          }

          double[] medSDLimits = {Math.max(actualCNV1.medianLRR - actualCNV1.stdevLRR,
                                           actualCNV2.medianLRR - actualCNV2.stdevLRR),
                                  Math.min(actualCNV1.medianLRR + actualCNV1.stdevLRR,
                                           actualCNV2.medianLRR + actualCNV2.stdevLRR)};

          if (cnv.medianLRR > medSDLimits[0] && cnv.medianLRR < medSDLimits[1]) {
            MergedCNVariant newCNV = new MergedCNVariant(actualCNV1.markerStart,
                                                         actualCNV2.markerStop,
                                                         new CNVariant(actualCNV1.cnv.getFamilyID(),
                                                                       actualCNV1.cnv.getIndividualID(),
                                                                       chr,
                                                                       actualCNV1.cnv.getStart(),
                                                                       actualCNV2.cnv.getStop(),
                                                                       actualCNV1.cnv.getCN(),
                                                                       Math.max(actualCNV1.cnv.getScore(),
                                                                                actualCNV2.cnv.getScore()),
                                                                       actualCNV2.markerStop - actualCNV1.markerStart + 1,
                                                                       actualCNV1.cnv.getSource()));

            setLRRMedStdDev(newCNV, actualCNV1, actualCNV2, droppedMarkerNames, markerNames, lrrs);

            // remove ICS
            chromo.remove(index);
            // remove CNV2
            chromo.remove(index);
            // remove CNV1
            chromo.remove(index - 1);

            chromo.add(index - 1, newCNV);

            // add removed ICSs, in case we've altered things enough to now combine more
            Integer addBack = indexList.pollLast();
            while (addBack != null) {
              chromo.add(addBack, removed.remove(addBack));
              addBack = indexList.pollLast();
            }
          } else {
            // not close enough to combine, so remove ICS placeholder
            removed.put(index, chromo.remove(index));
            indexList.add(index);
            continue;
          }

        }

        cnvCount += chromo.size();

        for (int i = 0; i < chromo.size(); i++) {
          // shouldn't have any ICSs remaining... but check anyway, just to be safe?
          if (chromo.get(i).cnv != null) {
            mergedCNVs.add(chromo.get(i).cnv);
          }
        }
      }
      cnt++;
    }

    log.report(ext.getTime() + "] CNV cleaning complete!  Started with " + initialCount
               + " CNVs, reduced to " + cnvCount + " CNVs");

    return mergedCNVs.toArray(new CNVariant[mergedCNVs.size()]);
  }

  /**
   * @param proj
   * @param in
   * @param out
   * @param distanceQuotient
   */
  public static void mergeCNVs(Project proj, String in, String out, double distanceQuotient) {
    Logger log = proj.getLog();
    log.report(ext.getTime() + "] Loading CNV file...");
    CNVariant[] srcCNVs = CNVariant.loadPlinkFile(in);
    int initialCount = srcCNVs.length;

    try {
      PrintWriter writer = Files.openAppropriateWriter(out);
      writer.println(ArrayUtils.toStr(CNVariant.PLINK_CNV_HEADER, "\t"));

      HashMap<String, HashMap<Byte, ArrayList<CNVariant>>> indivChrCNVMap = new HashMap<>();

      log.report(ext.getTime() + "] Organizing CNVs in memory...");
      for (CNVariant cnv : srcCNVs) {
        HashMap<Byte, ArrayList<CNVariant>> chrCNVMap = indivChrCNVMap.get(cnv.getFamilyID() + "\t"
                                                                           + cnv.getIndividualID());
        if (chrCNVMap == null) {
          chrCNVMap = new HashMap<>();
          indivChrCNVMap.put(cnv.getFamilyID() + "\t" + cnv.getIndividualID(), chrCNVMap);
        }
        ArrayList<CNVariant> cnvList = chrCNVMap.get(cnv.getChr());
        if (cnvList == null) {
          cnvList = new ArrayList<>();
          chrCNVMap.put(cnv.getChr(), cnvList);
        }
        cnvList.add(cnv);
      }

      srcCNVs = null;

      MarkerSetInfo markerSet = proj.getMarkerSet();
      int[][] positions = markerSet.getPositionsByChr();
      Hashtable<String, String> droppedMarkerNames = proj.getFilteredHash();
      String[] markerNames = markerSet.getMarkerNames();
      SampleData sampleData = proj.getSampleData(false);

      log.report(ext.getTime() + "] Cleaning CNVs...");
      int size = indivChrCNVMap.size();
      int cnt = 0;
      int step = size / 10;
      int cnvCount = 0;
      for (java.util.Map.Entry<String, HashMap<Byte, ArrayList<CNVariant>>> entry : indivChrCNVMap.entrySet()) {
        if (cnt > 0 && cnt % step == 0) {
          log.report(ext.getTime() + "] \t" + ((cnt / step) * 10) + "% complete");
        }
        // log.report(ext.getTime() + "] Analyzing CNVs for {" + entry.getKey() + "}");
        String fidiid = entry.getKey();

        String ind;
        float[] lrrs;
        try {
          ind = sampleData.lookup(fidiid)[0];
        } catch (NullPointerException npe) {
          // log.reportError("Error - could not look up the sample " + fidiid + " in the sample data
          // file " + proj.getFilename(proj.SAMPLE_DATA_FILENAME) + ", cannot load sample to compute
          // beast score");
          log.reportError("Error - could not look up the sample " + fidiid
                          + " in the sample data file " + proj.SAMPLE_DATA_FILENAME.getValue()
                          + ", cannot load sample to compute beast score");
          log.reportError("Error - please ensure that the sample names correspond to the varaints being processed with FID="
                          + fidiid.split("\t")[0] + " and IID=" + fidiid.split("\t")[1]);
          continue;
        }
        try {
          lrrs = proj.getFullSampleFromRandomAccessFile(ind).getLRRs();
        } catch (NullPointerException npe) {
          log.reportError("Error - could not load data for the sample " + ind + "\t" + fidiid
                          + ", please ensure samples have been parsed prior to computing beast score");
          log.report("Skipping beast score for sample " + ind + "\t" + fidiid);
          continue;
        }

        for (java.util.Map.Entry<Byte, ArrayList<CNVariant>> cnvLists : entry.getValue()
                                                                             .entrySet()) {
          byte chr = cnvLists.getKey();

          if (cnvLists.getValue().size() == 0) {
            continue;
          } else if (cnvLists.getValue().size() == 1) {
            writer.println(cnvLists.getValue().get(0).toPlinkFormat());
            continue;
          }

          ArrayList<CNVariant> cnvList = cnvLists.getValue();
          Collections.sort(cnvList, cnvComparator);

          LinkedList<MergedCNVariant> tempChromo = new LinkedList<>();

          // create objects for all CNVs while also setting start/end marker indices, accounting for
          // dropped markers
          for (int i = 0; i < cnvList.size(); i++) {
            int firstSNPIndex = ArrayUtils.binarySearch(positions[chr], cnvList.get(i).getStart(),
                                                        true);
            int lastSNPIndex = ArrayUtils.binarySearch(positions[chr], cnvList.get(i).getStop(),
                                                       true);

            // account for dropped markers
            if (droppedMarkerNames.contains(markerNames[firstSNPIndex])) {
              int f1 = firstSNPIndex, f2 = firstSNPIndex;
              while (f1 > (i == 0 ? 0 : tempChromo.get(i - 1).markerStop)
                     && droppedMarkerNames.contains(markerNames[f1])) {
                f1--;
              }
              while (f2 < lastSNPIndex - 1 && droppedMarkerNames.contains(markerNames[f2])) {
                f2++;
              }
              firstSNPIndex = ((firstSNPIndex - f1) < (f2 - firstSNPIndex) ? f1 : f2);
            }
            if (droppedMarkerNames.contains(markerNames[lastSNPIndex])) {
              int s1 = firstSNPIndex, s2 = lastSNPIndex;
              while (s1 > firstSNPIndex && droppedMarkerNames.contains(markerNames[s1])) {
                s1--;
              }
              while (s2 < markerNames.length - 1 && droppedMarkerNames.contains(markerNames[s2])) {
                s2++;
              }

              lastSNPIndex = ((lastSNPIndex - s1) < (s2 - lastSNPIndex) ? s1 : s2);
            }

            MergedCNVariant orig = new MergedCNVariant(firstSNPIndex, lastSNPIndex, cnvList.get(i));
            orig.addOriginalCNV(orig);
            tempChromo.add(orig);
          }

          LinkedList<MergedCNVariant> chromo = new LinkedList<>();

          // create full list: CNVs and objects for the spaces between CNVs, also accounting for
          // dropped markers
          chromo.add(tempChromo.get(0));
          for (int i = 1; i < tempChromo.size(); i++) {
            MergedCNVariant cnv1 = tempChromo.get(i - 1);
            MergedCNVariant cnv2 = tempChromo.get(i);

            if (cnv2.markerStart - cnv1.markerStop <= 1) {
              chromo.add(cnv2);
              continue;
            }

            int m1 = cnv1.markerStop + 1;
            int m2 = cnv2.markerStart - 1;

            if (droppedMarkerNames.contains(markerNames[m1])) {
              int f1 = m1;
              while (f1 < m2 - 1 && droppedMarkerNames.contains(markerNames[f1])) {
                f1++;
              }
              m1 = f1;
            }
            if (droppedMarkerNames.contains(markerNames[m2])) {
              int s1 = m2;
              while (s1 > m1 + 1 && droppedMarkerNames.contains(markerNames[s1])) {
                s1--;
              }
              m2 = s1;
            }

            chromo.add(new MergedCNVariant(m1, m2, null));
            chromo.add(cnv2);
          }
          tempChromo = null;

          // iterate through full array and compute median and stdev of LRRs
          for (int i = 0; i < chromo.size(); i++) {
            MergedCNVariant cnv = chromo.get(i);
            setLRRMedStdDev(cnv, droppedMarkerNames, markerNames, lrrs);
          }

          HashMap<Integer, MergedCNVariant> removed = new HashMap<>();
          TreeSet<Integer> indexList = new TreeSet<>();
          boolean done = false;
          while (!done) {
            int index = 0;
            MergedCNVariant cnv = chromo.get(0);
            // find index of first ICS (inter-CNV space); should never be index 0 or max
            for (index = 1; index < chromo.size() - 1 && cnv.cnv != null; index++) {
              cnv = chromo.get(index);
            }

            if (cnv.cnv != null) {
              // iterated through entire list and no ICSs left!
              done = true;
              break;
            }
            index--;

            MergedCNVariant actualCNV1 = chromo.get(index - 1);
            MergedCNVariant actualCNV2 = chromo.get(index + 1);

            int bpSize = countBP(actualCNV1) + countBP(actualCNV2);

            // same CN
            boolean sameCN = actualCNV1.cnv.getCN() == actualCNV2.cnv.getCN();
            // less than distanceQuotient percent space vs # of markers
            float mkQ = (float) ((cnv.markerStop - cnv.markerStart + 1))
                        / (float) ((actualCNV1.cnv.getNumMarkers()
                                    + actualCNV2.cnv.getNumMarkers()));
            boolean markerQuotient = mkQ < distanceQuotient;
            // less than 100% of total called base pairs
            float bpQ = (actualCNV2.cnv.getStart() - actualCNV1.cnv.getStop() + 1) / (float) bpSize;
            boolean basepairQuotient = bpQ < 1;

            if (!sameCN || !markerQuotient || !basepairQuotient) {
              // remove ICS placeholder
              removed.put(index, chromo.remove(index));
              indexList.add(index);
              continue;
            }

            double[] medSDLimits = {Math.max(actualCNV1.medianLRR - actualCNV1.stdevLRR,
                                             actualCNV2.medianLRR - actualCNV2.stdevLRR),
                                    Math.min(actualCNV1.medianLRR + actualCNV1.stdevLRR,
                                             actualCNV2.medianLRR + actualCNV2.stdevLRR)};

            if (cnv.medianLRR > medSDLimits[0] && cnv.medianLRR < medSDLimits[1]) {
              MergedCNVariant newCNV = new MergedCNVariant(actualCNV1.markerStart,
                                                           actualCNV2.markerStop,
                                                           new CNVariant(actualCNV1.cnv.getFamilyID(),
                                                                         actualCNV1.cnv.getIndividualID(),
                                                                         chr,
                                                                         actualCNV1.cnv.getStart(),
                                                                         actualCNV2.cnv.getStop(),
                                                                         actualCNV1.cnv.getCN(),
                                                                         Math.max(actualCNV1.cnv.getScore(),
                                                                                  actualCNV2.cnv.getScore()),
                                                                         actualCNV2.markerStop - actualCNV1.markerStart
                                                                                                              + 1,
                                                                         actualCNV1.cnv.getSource()));

              setLRRMedStdDev(newCNV, actualCNV1, actualCNV2, droppedMarkerNames, markerNames,
                              lrrs);

              // remove ICS
              chromo.remove(index);
              // remove CNV2
              chromo.remove(index);
              // remove CNV1
              chromo.remove(index - 1);

              chromo.add(index - 1, newCNV);

              // add removed ICSs, in case we've altered things enough to now combine more
              Integer addBack = indexList.pollLast();
              while (addBack != null) {
                chromo.add(addBack, removed.remove(addBack));
                addBack = indexList.pollLast();
              }
            } else {
              // not close enough to combine, so remove ICS placeholder
              removed.put(index, chromo.remove(index));
              indexList.add(index);
              continue;
            }

          }

          cnvCount += chromo.size();

          for (int i = 0; i < chromo.size(); i++) {
            // shouldn't have any ICSs remaining... but check anyway, just to be safe?
            if (chromo.get(i).cnv != null) {
              writer.println(chromo.get(i).cnv.toPlinkFormat());
              // if (newLargeCNV(chromo.get(i), positions[chr])) {
              // verboseWriter.println(chromo.get(i).cnv.toPlinkFormat());
              // }
            }
          }
        }
        cnt++;
      }

      writer.flush();
      writer.close();
      // verboseWriter.flush();
      // verboseWriter.close();

      log.report(ext.getTime() + "] CNV cleaning complete!  Started with " + initialCount
                 + " CNVs, reduced to " + cnvCount + " CNVs");
    } catch (IOException e) {
      proj.getLog().reportException(e);
    }
  }

  //
  // private static boolean newLargeCNV(MergedCNVariant cleanCNVariant, int[] positions) {
  // int newSize = positions[cleanCNVariant.markerStop] - positions[cleanCNVariant.markerStart] + 1;
  // int sz = 1000;
  // if (cleanCNVariant.originalCNVs.get(0).cnv.getCN() < 0) {
  // sz = 250;
  // }
  // return newSize > sz * 1000;
  // }

  private static void setLRRMedStdDev(MergedCNVariant cnv,
                                      Hashtable<String, String> droppedMarkerNames,
                                      String[] markerNames, float[] lrrs) {
    ArrayList<Float> lrr = new ArrayList<>();
    for (int m = cnv.markerStart; m <= cnv.markerStop; m++) {
      if (droppedMarkerNames.contains(markerNames[m])) {
        continue;
      }
      if (!Double.isNaN(lrrs[m])) {
        lrr.add(lrrs[m]);
      }
    }
    double[] lrrsDubs = Doubles.toArray(lrr);
    cnv.medianLRR = ArrayUtils.median(lrrsDubs);
    cnv.stdevLRR = ArrayUtils.stdev(lrrsDubs);
  }

  private static void setLRRMedStdDev(MergedCNVariant newCNV, MergedCNVariant oldCNV1,
                                      MergedCNVariant oldCNV2,
                                      Hashtable<String, String> droppedMarkerNames,
                                      String[] markerNames, float[] lrrs) {
    ArrayList<MergedCNVariant> allOriginalCNVs = new ArrayList<>();
    allOriginalCNVs.addAll(oldCNV1.originalCNVs);
    allOriginalCNVs.addAll(oldCNV2.originalCNVs);

    ArrayList<Float> lrr = new ArrayList<>();
    for (MergedCNVariant cnv : allOriginalCNVs) {
      for (int m = cnv.markerStart; m <= cnv.markerStop; m++) {
        if (droppedMarkerNames.contains(markerNames[m])) {
          continue;
        }
        if (!Double.isNaN(lrrs[m])) {
          lrr.add(lrrs[m]);
        }
      }
    }

    double[] lrrsDubs = Doubles.toArray(lrr);
    newCNV.medianLRR = ArrayUtils.median(lrrsDubs);
    newCNV.stdevLRR = ArrayUtils.stdev(lrrsDubs);

    newCNV.originalCNVs.addAll(allOriginalCNVs);
  }

  private static int countBP(MergedCNVariant actualCNV1) {
    int bpCount = 0;
    for (MergedCNVariant cnv : actualCNV1.originalCNVs) {
      bpCount += cnv.cnv.getSize();
    }
    return bpCount;
  }

  public static void mergeCNVs(String in, String out, float distanceQuotient, String bimFile) {
    List<CNVariant> inputCNVs = CNVariant.loadPlinkFile(in, null, true);

    int[][] positions = null;
    if (bimFile != null) {
      MarkerSetInfo markerSet = MarkerSet.load(bimFile);
      positions = markerSet.getPositionsByChr();
    }

    ArrayList<CNVariant> newCNVs = new ArrayList<>();
    newCNVs.addAll(inputCNVs);
    inputCNVs = null;
    int cnt = 1;
    int startCnt = newCNVs.size();
    int actStart = startCnt;
    do {
      startCnt = newCNVs.size();
      newCNVs = getMergedCNVs(newCNVs, distanceQuotient, positions);
      cnt++;
    } while (startCnt > newCNVs.size());

    try {
      PrintWriter writer = Files.openAppropriateWriter(out);
      writer.println(ArrayUtils.toStr(CNVariant.PLINK_CNV_HEADER, "\t"));
      for (CNVariant cnv : newCNVs) {
        writer.println(cnv.toPlinkFormat());
      }
      writer.flush();
      writer.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
    System.out.println(cnt + " iterations; from " + actStart + " to " + newCNVs.size() + " CNVs");
  }

  private static ArrayList<CNVariant> getMergedCNVs(List<CNVariant> inputCNVs, float distFactor,
                                                    int[][] positions) {

    HashMap<String, HashMap<Byte, ArrayList<CNVariant>>> indivChrCNVMap = new HashMap<>();

    for (CNVariant cnv : inputCNVs) {
      HashMap<Byte, ArrayList<CNVariant>> chrCNVMap = indivChrCNVMap.get(cnv.getFamilyID() + "\t"
                                                                         + cnv.getIndividualID());
      if (chrCNVMap == null) {
        chrCNVMap = new HashMap<>();
        indivChrCNVMap.put(cnv.getFamilyID() + "\t" + cnv.getIndividualID(), chrCNVMap);
      }
      ArrayList<CNVariant> cnvList = chrCNVMap.get(cnv.getChr());
      if (cnvList == null) {
        cnvList = new ArrayList<>();
        chrCNVMap.put(cnv.getChr(), cnvList);
      }
      cnvList.add(cnv);
    }

    ArrayList<CNVariant> newCNVs = new ArrayList<>();

    StringBuilder status = new StringBuilder();

    for (java.util.Map.Entry<String, HashMap<Byte, ArrayList<CNVariant>>> entry : indivChrCNVMap.entrySet()) {
      String fidiid = entry.getKey();

      for (java.util.Map.Entry<Byte, ArrayList<CNVariant>> chrEntry : entry.getValue().entrySet()) {
        byte chr = chrEntry.getKey();
        ArrayList<CNVariant> cnvList = chrEntry.getValue();

        Collections.sort(cnvList, cnvComparator);
        if (cnvList.size() == 1) {
          newCNVs.add(cnvList.get(0));
        } else {
          for (int i = 0; i < cnvList.size() - 1; i++) {
            CNVariant curr = cnvList.get(i);
            CNVariant next = cnvList.get(i + 1);
            if (curr.getCN() != next.getCN()) {
              newCNVs.add(curr);
              if (i == cnvList.size() - 2) {
                newCNVs.add(next);
              }
              continue;
            }

            float[] szDiff = getCombinedSizeAndMarkerCountBetween(curr, next, positions);
            float sz = szDiff[0];
            float diff = szDiff[1];

            if (diff / sz <= distFactor) {
              newCNVs.add(new CNVariant(curr.getFamilyID(), curr.getIndividualID(), chr,
                                        curr.getStart(), next.getStop(), curr.getCN(),
                                        Math.max(curr.getScore(), next.getScore()),
                                        (int) (curr.getNumMarkers() + next.getNumMarkers()
                                               + (null == positions ? 0 : szDiff[1])),
                                        0));
              status.append(fidiid).append(" > ").append(curr.getChr()).append("{")
                    .append(curr.getStart()).append(", ").append(next.getStop()).append("}")
                    .append("\n");
              // skip combined CNVs
              i++;
              if (i == cnvList.size() - 2) {
                newCNVs.add(cnvList.get(i + 1));
              }
            } else {
              newCNVs.add(curr);
              if (i == cnvList.size() - 2) {
                newCNVs.add(next);
              }
            }

          }
        }

      }

    }

    indivChrCNVMap = null;

    System.out.println("Created " + (inputCNVs.size() - newCNVs.size()) + " combined CNVs");
    System.out.println(status.toString());
    System.out.println("------------------------------------------------");

    return newCNVs;
  }

  private static float[] getCombinedSizeAndMarkerCountBetween(CNVariant curr, CNVariant next,
                                                              int[][] positions) {
    if (null == positions) {
      return new float[] {curr.getSize() + next.getSize(), next.getStart() - curr.getStop()};
    }

    int lastSNPIndex1 = ArrayUtils.binarySearch(positions[curr.getChr()], curr.getStop(), true);
    int firstSNPIndex2 = ArrayUtils.binarySearch(positions[next.getChr()], next.getStart(), true);

    float sz = (curr.getNumMarkers() + next.getNumMarkers());
    float dist = firstSNPIndex2 - lastSNPIndex1 - 2;

    return new float[] {sz, dist};
  }

  // /**
  // * Given a list of CNVs, and a list of sample IDs, find all CNVs that overlap at least one CNV
  // belonging to a sample ID on the given list.
  // *
  // * Useful for data without controls, but with known affected samples.
  // *
  // * @param dir location of in / out files
  // * @param in Plink-formatted CNV file for input
  // * @param out desired name of filtered CNV file
  // * @param listFile file name (full path) of the list of sample IDs
  // */
  // public static void filterForAllCNVsSharedInGroup(String dir, String in, String out, String
  // listFile) {
  // PrintWriter writer;
  // Vector<CNVariant> cnvs = CNVariant.loadPlinkFile(dir + in, null, false);
  // boolean[] remove = new boolean[cnvs.size()];
  // HashSet<String> indivList = HashVec.loadFileToHashSet(listFile, false);
  //
  // for (int i = 0; i < remove.length; i++) {
  // CNVariant examine = cnvs.get(i);
  //
  // boolean mark = true;
  // for (CNVariant groupCNV : cnvs) {
  // if (!indivList.contains(groupCNV.getIndividualID())) { // loop through only those belonging to
  // listed sample IDs
  // continue;
  // }
  //
  // if (examine.overlaps(groupCNV)) {
  // mark = false;
  // break;
  // }
  // }
  // remove[i] = mark;
  // }
  //
  // try {
  // writer = Files.openAppropriateWriter(dir + out);
  // writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER, "\t"));
  // for (int i = 0; i < remove.length; i++) {
  // if (!remove[i]) writer.println(cnvs.get(i).toPlinkFormat());
  // }
  // writer.flush();
  // writer.close();
  // } catch (IOException e) {
  // e.printStackTrace();
  // }
  // }
  //
  //
  // public static void filterForGroupCNVs(String dir, String in, String out, String listFile,
  // String excludeFile, boolean excludeCommon) {
  // PrintWriter writer;
  // Vector<CNVariant> cnvs = CNVariant.loadPlinkFile(dir + in, null, false);
  // boolean[] remove = new boolean[cnvs.size()];
  // HashSet<String> indivList = HashVec.loadFileToHashSet(listFile, false);
  // HashSet<String> excludeList = excludeFile == null ? new HashSet<String>() :
  // HashVec.loadFileToHashSet(excludeFile, false);
  //
  // for (int i = 0; i < remove.length; i++) {
  // CNVariant examine = cnvs.get(i);
  // if (!indivList.contains(examine.getIndividualID()) ||
  // excludeList.contains(examine.getIndividualID())) { // remove all CNVs not belonging to a listed
  // sample ID
  // remove[i] = true;
  // continue;
  // }
  //
  // if (excludeCommon) {
  // for (int j = 0; j < remove.length; j++) {
  // if (i == j) continue;
  // CNVariant check = cnvs.get(j);
  // if (examine.significantOverlap(check, true) && !indivList.contains(check.getIndividualID())) {
  // // if the sample ID of the second CNV isn't on the list, and the current cnv is on the list,
  // remove both
  // remove[i] = true;
  // remove[j] = true;
  // }
  // }
  // }
  //
  // }
  //
  // try {
  // writer = Files.openAppropriateWriter(dir + out);
  // writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER, "\t"));
  // for (int i = 0; i < remove.length; i++) {
  // if (!remove[i]) writer.println(cnvs.get(i).toPlinkFormat());
  // }
  // writer.flush();
  // writer.close();
  // } catch (IOException e) {
  // e.printStackTrace();
  // }
  //
  // }
  //
  //
  // public static void filterOutCommonCNVs(String dir, String in, String out, double pct, boolean
  // checkLarger) {
  // PrintWriter writer;
  // Vector<CNVariant> cnvs = CNVariant.loadPlinkFile(dir + in, null, false);
  // boolean[] remove = new boolean[cnvs.size()];
  //
  // for (int i = 0; i < remove.length; i++) {
  // CNVariant examine = cnvs.get(i);
  // int overlapCnt = 0;
  //
  // for (int j = 0; j < remove.length; j++) {
  // if (i == j) continue;
  //
  // if (examine.significantOverlap(cnvs.get(j), checkLarger)) {
  // overlapCnt++;
  // }
  // }
  //
  // if (((double)overlapCnt) / ((double)cnvs.size()) > pct) {
  // remove[i] = true;
  // } else {
  // remove[i] = false;
  // }
  // }
  //
  // try {
  // writer = Files.openAppropriateWriter(dir + out);
  // writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER, "\t"));
  // for (int i = 0; i < remove.length; i++) {
  // if (!remove[i]) writer.println(cnvs.get(i).toPlinkFormat());
  // }
  // writer.flush();
  // writer.close();
  // } catch (IOException e) {
  // e.printStackTrace();
  // }
  // }
  //
  // public static void filterCommonCNVs(String dir, String in, String out, double pct) {
  // PrintWriter writer;
  // Vector<CNVariant> cnvs = CNVariant.loadPlinkFile(dir + in, null, false);
  // boolean[] remove = new boolean[cnvs.size()];
  //
  // HashMap<Integer, ArrayList<CNVariant>> chrToCNVs = new HashMap<Integer,
  // ArrayList<CNVariant>>();
  // HashMap<CNVariant, Integer> indexMap = new HashMap<CNVariant, Integer>();
  // for (int index = 0; index < cnvs.size(); index++) {
  // ArrayList<CNVariant> chrList = chrToCNVs.get((int) cnvs.get(index).getChr());
  // if (chrList == null) {
  // chrList = new ArrayList<CNVariant>();
  // chrToCNVs.put((int) cnvs.get(index).getChr(), chrList);
  // }
  // chrList.add(cnvs.get(index));
  // indexMap.put(cnvs.get(index), index);
  // }
  //
  // int hundredthStep = (cnvs.size() / 100) + 1;
  //
  // for (int i = 0; i < remove.length; i++) {
  // CNVariant examine = cnvs.get(i);
  // int overlapCnt = 0;
  //
  // ArrayList<CNVariant> compCNVs = chrToCNVs.get((int) examine.getChr());
  // for (CNVariant comp : compCNVs) {
  // if (i == indexMap.get(comp).intValue()) continue;
  //
  // if (examine.significantOverlap(comp)) {
  // overlapCnt++;
  // }
  // }
  //
  // if (((double)overlapCnt) / ((double)cnvs.size()) > pct) {
  // remove[i] = true;
  // } else {
  // remove[i] = false;
  // }
  //
  // if (i > 0 && i % hundredthStep == 0) {
  // System.out.println(ext.getTime() + "] Examined " + i + " of " + cnvs.size() + " CNVs...");
  // }
  // }
  //
  // try {
  // writer = Files.openAppropriateWriter(dir + out);
  // writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER, "\t"));
  // for (int i = 0; i < remove.length; i++) {
  // if (!remove[i]) writer.println(cnvs.get(i).toPlinkFormat());
  // }
  // writer.flush();
  // writer.close();
  // } catch (IOException e) {
  // e.printStackTrace();
  // }
  // }
  //

  public static void filterExclusions(Project proj, String cnvFile) {
    PrintWriter writer;
    BufferedReader reader;
    String path = cnvFile.substring(0, cnvFile.lastIndexOf('/'));
    String fileext = cnvFile.substring(cnvFile.lastIndexOf('/'));
    String filenm = fileext.substring(0, fileext.lastIndexOf('.'));
    String newFile = path + filenm + ".excluded.cnv";

    String[] excludes = ArrayUtils.subArray(proj.getSamples(), proj.getSamplesToExclude());
    HashSet<String> excludeSet = new HashSet<>();
    for (String exclude : excludes) {
      excludeSet.add(exclude);
    }

    try {
      reader = new BufferedReader(new FileReader(cnvFile));
      writer = Files.openAppropriateWriter(newFile);

      while (reader.ready()) {
        String line = reader.readLine();
        if (excludeSet.contains(line.split("\t")[0])) {
          continue;
        }
        writer.println(line);
      }

      writer.flush();

      reader.close();
      writer.close();
    } catch (FileNotFoundException e) {
      proj.getLog().reportException(e);
    } catch (IOException e) {
      proj.getLog().reportException(e);
    }
  }

  public static void filterExclusions(String dir, String in, String out, String indivFile,
                                      boolean exclude) {
    PrintWriter writer;
    List<CNVariant> cnvs = CNVariant.loadPlinkFile(dir + in, null, true);
    boolean[] remove = ArrayUtils.booleanArray(cnvs.size(), !exclude);
    HashSet<String> indivList = HashVec.convertHashNullToHashSet(HashVec.loadFileToHashString(indivFile,
                                                                                              new int[] {0,
                                                                                                         1},
                                                                                              null,
                                                                                              false,
                                                                                              "\t",
                                                                                              false,
                                                                                              false));

    int cnt = 0;
    for (int i = 0; i < remove.length; i++) {
      CNVariant examine = cnvs.get(i);
      if (indivList.contains(examine.getFamilyID() + "\t" + examine.getIndividualID())) {
        cnt++;
        remove[i] = exclude;
        continue;
      }
    }

    System.out.println((exclude ? "Removed " : "Matched ") + cnt + " of " + cnvs.size() + " CNVs");

    try {
      writer = Files.openAppropriateWriter(dir + out);
      writer.println(ArrayUtils.toStr(CNVariant.PLINK_CNV_HEADER, "\t"));
      for (int i = 0; i < remove.length; i++) {
        if (!remove[i]) {
          writer.println(cnvs.get(i).toPlinkFormat());
        }
      }
      writer.flush();
      writer.close();
    } catch (IOException e) {
      e.printStackTrace();
    }

  }

  private static void configureFilterBasics(CNVFilter filter, String filenameOfProblematicRegions,
                                            int commonInOutOrIgnore, String individualsToKeepFile,
                                            int[][] pos, int[][] bnds, boolean makeUCSCtrack,
                                            int build, Logger log) {
    filter.setBuild(build);
    if (pos != null && bnds != null) {
      filter.setPositions(pos);
      filter.setCentromereBoundaries(bnds);
      filter.computeCentromereMidPoints();
      filter.setBreakupCentromeres(true);
    }
    if (commonInOutOrIgnore != COMMON_IGNORED) {
      filter.setCommonIn(commonInOutOrIgnore == COMMON_IN);
      filter.setCommonReference(Segment.loadUCSCregions(Files.firstDirectoryThatExists(DEFAULT_REGION_DIRECTORIES,
                                                                                       true, true,
                                                                                       log)
                                                        + DEFAULT_COMMON_CNP_REFERENCE, false));
    }
    filter.setIndividualsToKeepFromFile(individualsToKeepFile);
    if (filenameOfProblematicRegions != null && !"".equals(filenameOfProblematicRegions)) {
      filter.setProblemRegions(Segment.loadUCSCregions(filenameOfProblematicRegions, 0, false,
                                                       log));
    }
  }

  public static void filterCNVs(String dir, String in, String out, int[] delSize, int[] dupSize,
                                int[] number, double score, String filenameOfProblematicRegions,
                                int commonInOutOrIgnore, String individualsToKeepFile,
                                boolean breakupCentromeres,
                                String markerSetFilenameToBreakUpCentromeres, boolean makeUCSCtrack,
                                int build, Logger log) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    CNVariant cnv;

    int[][] pos = null, bnds = null;
    if (breakupCentromeres && markerSetFilenameToBreakUpCentromeres != null
        && !markerSetFilenameToBreakUpCentromeres.equals("")
        && Files.exists(markerSetFilenameToBreakUpCentromeres)) {
      if (markerSetFilenameToBreakUpCentromeres.endsWith(".bim")
          || markerSetFilenameToBreakUpCentromeres.endsWith(".map")
          || markerSetFilenameToBreakUpCentromeres.endsWith(".txt")) {
        SnpMarkerSet markerSet = new SnpMarkerSet(markerSetFilenameToBreakUpCentromeres, false,
                                                  log);
        markerSet.sortMarkers();
        pos = markerSet.getPositionsByChr();
        bnds = Positions.determineCentromereBoundariesFromMarkerSet(markerSet.getChrs(),
                                                                    markerSet.getPositions(), build,
                                                                    log);
      } else if (markerSetFilenameToBreakUpCentromeres.endsWith(".ser")) {
        SnpMarkerSet markerSet = SnpMarkerSet.load(markerSetFilenameToBreakUpCentromeres);
        markerSet.sortMarkers();
        pos = markerSet.getPositionsByChr();
        bnds = Positions.determineCentromereBoundariesFromMarkerSet(markerSet.getChrs(),
                                                                    markerSet.getPositions(), build,
                                                                    log);
      }
    }

    // first filter
    CNVFilter filter = new CNVFilter(log);
    configureFilterBasics(filter, filenameOfProblematicRegions, commonInOutOrIgnore,
                          individualsToKeepFile, pos, bnds, makeUCSCtrack, build, log);
    filter.setMinScore(score);

    CNVFilter delFilter0 = new CNVFilter(log);
    configureFilterBasics(delFilter0, filenameOfProblematicRegions, commonInOutOrIgnore,
                          individualsToKeepFile, pos, bnds, makeUCSCtrack, build, log);
    delFilter0.setMinNumMarkers(number[1]);
    delFilter0.setMinSize(delSize[1]);

    CNVFilter delFilter1 = new CNVFilter(log);
    configureFilterBasics(delFilter1, filenameOfProblematicRegions, commonInOutOrIgnore,
                          individualsToKeepFile, pos, bnds, makeUCSCtrack, build, log);
    delFilter1.setMinNumMarkers(number[0]);
    delFilter1.setMinSize(delSize[0]);

    CNVFilter dupFilter3 = new CNVFilter(log);
    configureFilterBasics(dupFilter3, filenameOfProblematicRegions, commonInOutOrIgnore,
                          individualsToKeepFile, pos, bnds, makeUCSCtrack, build, log);
    dupFilter3.setMinNumMarkers(number[0]);
    dupFilter3.setMinSize(dupSize[0]);

    CNVFilter dupFilter4 = new CNVFilter(log);
    configureFilterBasics(dupFilter4, filenameOfProblematicRegions, commonInOutOrIgnore,
                          individualsToKeepFile, pos, bnds, makeUCSCtrack, build, log);
    dupFilter4.setMinNumMarkers(number[0]);
    dupFilter4.setMinSize(dupSize[0]);

    int lineNumber = 1;
    try {
      reader = new BufferedReader(new FileReader(dir + in));
      writer = Files.openAppropriateWriter(dir + out);
      // header
      writer.println(reader.readLine());
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        if (line.length < 8) {
          log.reportError("Error - line " + lineNumber
                          + " is invalid (does not contain 8+ tokens, i.e. PLINK format): " + line);
          continue;
        }
        cnv = new CNVariant(line);

        CNVFilterPass filterResult = filter.getCNVFilterPass(cnv);
        // first pass filter:
        if (filterResult.passedFilter()) {
          // TODO Currently filters out centromeric CNVs!
          if (filterResult.isCentromeric()) {
            CNVariant[] newCNVs = filter.breakUpCentromere(filterResult, cnv);
            for (CNVariant newcnv : newCNVs) {
              // re-filter newly-split CNVs
              if (filter.getCNVFilterPass(newcnv).passedFilter()) {
                boolean passes = false;
                // CN-specific filters
                switch (newcnv.getCN()) {
                  case 0:
                    passes = delFilter0.getCNVFilterPass(newcnv).passedFilter();
                    break;
                  case 1:
                    passes = delFilter1.getCNVFilterPass(newcnv).passedFilter();
                    break;
                  case 3:
                    passes = dupFilter3.getCNVFilterPass(newcnv).passedFilter();
                    break;
                  case 4:
                  case 5:
                  case 6:
                    passes = dupFilter4.getCNVFilterPass(newcnv).passedFilter();
                    break;
                }
                if (passes) {
                  writer.println(newcnv.toPlinkFormat());
                }
              }
            }
          } else {
            boolean passes = false;
            // CN-specific filters
            switch (cnv.getCN()) {
              case 0:
                passes = delFilter0.getCNVFilterPass(cnv).passedFilter();
                break;
              case 1:
                passes = delFilter1.getCNVFilterPass(cnv).passedFilter();
                break;
              case 3:
                passes = dupFilter3.getCNVFilterPass(cnv).passedFilter();
                break;
              case 4:
              case 5:
              case 6:
                passes = dupFilter4.getCNVFilterPass(cnv).passedFilter();
                break;
            }
            if (passes) {
              writer.println(cnv.toPlinkFormat());
            }
          }
        }
        lineNumber++;
      }

      reader.close();
      writer.flush();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + in + "\" not found in current directory");
      return;
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + in + "\"");
      return;
    }
  }

  @Deprecated
  public static void filter(String dir, String in, String out, int[] delSize, int[] dupSize,
                            int[] number, double score, String filenameOfProblematicRegions,
                            int commonInOutOrIgnore, String individualsToKeepFile,
                            boolean breakupCentromeres,
                            String markerSetFilenameToBreakUpCentromeres, boolean makeUCSCtrack,
                            int build, Logger log) {
    String[] individualsToKeepList;
    if (individualsToKeepFile != null && !new File(individualsToKeepFile).exists()) {
      System.err.println("Error - could not find \"" + individualsToKeepFile
                         + "\" in directory; will not be able to filter by indiviudals");
      individualsToKeepFile = null;
    }
    individualsToKeepList = individualsToKeepFile == null ? null
                                                          : HashVec.loadFileToStringArray(individualsToKeepFile,
                                                                                          false,
                                                                                          new int[] {0,
                                                                                                     1},
                                                                                          true,
                                                                                          false,
                                                                                          "\t");

    filter(dir, in, out, delSize, dupSize, number, score, filenameOfProblematicRegions,
           commonInOutOrIgnore, individualsToKeepList, breakupCentromeres,
           markerSetFilenameToBreakUpCentromeres, makeUCSCtrack, build, log);
  }

  @Deprecated
  public static void filter(String dir, String in, String out, int[] delSize, int[] dupSize,
                            int[] number, double score, String filenameOfProblematicRegions,
                            int commonInOutOrIgnore, String[] individualsToKeepList,
                            boolean breakupCentromeres,
                            String markerSetFilenameToBreakUpCentromeres, boolean makeUCSCtrack,
                            int build, Logger log) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    CNVariant cnv;
    Segment[] problemRegions, centromereMidpoints, commonReference;
    HashSet<String> indHash;
    int countGiant, countCentromeric, countGiantCentromeric;
    int[][] centromereBoundaries;

    problemRegions = filenameOfProblematicRegions == null ? new Segment[0]
                                                          : Segment.loadUCSCregions(filenameOfProblematicRegions,
                                                                                    0, false, log);
    centromereBoundaries = Positions.determineCentromereBoundariesFromMarkerSet(markerSetFilenameToBreakUpCentromeres,
                                                                                build, log);
    centromereMidpoints = Positions.computeCentromereMidpoints(centromereBoundaries);
    commonReference = commonInOutOrIgnore != COMMON_IGNORED ? Segment.loadUCSCregions(Files.firstDirectoryThatExists(DEFAULT_REGION_DIRECTORIES,
                                                                                                                     true,
                                                                                                                     true,
                                                                                                                     log)
                                                                                      + DEFAULT_COMMON_CNP_REFERENCE,
                                                                                      false)
                                                            : new Segment[0];
    indHash = individualsToKeepList == null ? null : HashVec.loadToHashSet(individualsToKeepList);

    try {
      reader = new BufferedReader(new FileReader(dir + in));
      writer = Files.openAppropriateWriter(dir + out);
      System.out.println("Writing to '" + dir + out + "'");
      writer.println(reader.readLine());
      countGiant = 0;
      countCentromeric = 0;
      countGiantCentromeric = 0;
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        cnv = new CNVariant(line);
        if (((cnv.getCN() == 1 && cnv.getSize() >= delSize[0] * 1000) || // heterozygous deletion
             (cnv.getCN() == 0 && cnv.getSize() >= delSize[1] * 1000) || // homozygous deletion
             (cnv.getCN() > 2 && cnv.getSize() >= dupSize[0] * 1000) // duplications
        // ignoring homozygotic duplications
        ) && (((cnv.getCN() == 1 || cnv.getCN() == 3 || cnv.getCN() == 4)
               && cnv.getNumMarkers() >= number[0])
              || (cnv.getCN() == 0 && cnv.getNumMarkers() >= number[1]))
            && cnv.getScore() > score && !inOneOfTheseRegions(cnv, problemRegions)) {

          if ((commonInOutOrIgnore == COMMON_IGNORED
               || (commonInOutOrIgnore == COMMON_IN && inOneOfTheseRegions(cnv, commonReference))
               || (commonInOutOrIgnore == COMMON_OUT && !inOneOfTheseRegions(cnv, commonReference)))
              && (indHash == null || indHash.contains(line[0] + "\t" + line[1]))) {
            if (cnv.overlaps(centromereMidpoints[cnv.getChr()])) {
              if (breakupCentromeres) {
                // System.out.println("Splitting "+cnv.getUCSClocation()+" due to overlap with
                // "+centromereMidpoints[cnv.getChr()].getUCSClocation()+" using boundaries
                // "+Array.toStr(centromereBoundaries[cnv.getChr()], ", "));

                // TODO update marker count for newly-broken CNV
                line[3] = cnv.getStart() + "";
                line[4] = centromereBoundaries[cnv.getChr()][0] + "";
                writer.println(ArrayUtils.toStr(line));
                line[3] = centromereBoundaries[cnv.getChr()][1] + "";
                line[4] = cnv.getStop() + "";
                writer.println(ArrayUtils.toStr(line));
                // return;
              }
              countCentromeric++;

              if (cnv.getSize() > 10000000) {
                countGiantCentromeric++;
              }

              // System.err.println("Warning - a CNV for
              // "+cnv.getFamilyID()+","+cnv.getIndividualID()+" spans a centromere
              // ("+cnv.getUCSClocation()+") with "+cnv.getNumMarkers()+" markers");
            } else {
              writer.println(ArrayUtils.toStr(line));
            }
          }
          if (cnv.getSize() > 10000000 || cnv.getNumMarkers() > 500) {
            // System.err.println("Warning - "+cnv.getFamilyID()+","+cnv.getIndividualID()+" has a
            // gigantic CNV spanning "+ext.prettyUpDistance(cnv.getSize(), 0)+" and
            // "+cnv.getNumMarkers()+" markers ("+cnv.getUCSClocation()+")");
            countGiant++;
          }
        }
      }
      System.err.println("Identified " + countCentromeric
                         + " CNVs that spanned centromeres; these were "
                         + (breakupCentromeres ? "broken up into two CNVs, one on each side of the centromere"
                                               : "retained as is"));
      System.err.println("Identified " + countGiant
                         + " gigantic CNVs ( 10+ Mb or 500+ probes ), of which "
                         + countGiantCentromeric + " spanned a centromere");
      reader.close();
      writer.close();
      if (makeUCSCtrack) {
        UCSCtrack.makeTrack(dir + out, dir + ext.rootOf(out) + ".bed.gz", log);
      }
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + in + "\" not found in current directory");
      return;
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + in + "\"");
      return;
    }
  }

  public static void filterOnSegments(String dir, String filein, String fileout, String segmentFile,
                                      boolean excludeInsteadOfInclude) {
    BufferedReader reader;
    PrintWriter writer;
    Segment[][] genesByChr;
    CNVariant cnv;
    SegmentLists segList;

    if (segmentFile.endsWith(".segs")) {
      segList = SegmentLists.load(segmentFile);
    } else if (new File(segmentFile + ".segs").exists()) {
      segList = SegmentLists.load(segmentFile + ".segs");
    } else {
      segList = SegmentLists.parseUCSCSegmentList(segmentFile, false);
      segList.serialize(segmentFile + ".segs");
    }

    genesByChr = segList.getLists();

    try {
      reader = new BufferedReader(new FileReader(dir + filein));
      writer = Files.openAppropriateWriter(dir + fileout);
      writer.println(reader.readLine());
      while (reader.ready()) {
        cnv = new CNVariant(reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE));
        if (genesByChr[cnv.getChr()] != null
            && Segment.overlapsAny(new Segment(cnv.getChr(), cnv.getStart(), cnv.getStop()),
                                   genesByChr[cnv.getChr()])) {
          writer.println(cnv.toPlinkFormat());
        }
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + filein + "\" not found in current directory");
      return;
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + filein + "\"");
      return;
    }
  }

  /**
   * @param cnvsAsPositions get the start and stop positions from the {@link CNVariant}s themselves
   *          so that a {@link MarkerSet} is not needed
   */
  public static CNVariant[] filterBasedOnNumberOfCNVsAtLocusInMemory(Project proj, CNVariant[] cnvs,
                                                                     int totalRequired,
                                                                     int delRequired,
                                                                     int dupRequired,
                                                                     int totalLimitedTo,
                                                                     int delLimitedTo,
                                                                     int dupLimitedTo,
                                                                     double proportionOfProbesThatNeedToPassForFinalInclusion,
                                                                     boolean cnvsAsPositions) {
    MarkerSetInfo markerSet;
    int[][] positions;
    int[][][] counts;
    int firstSNP, lastSNP, indel;
    ;
    int index;
    boolean[][] acceptableSNPs;
    boolean accepted;
    int dels, dups;
    int countAcceptable;
    long time;

    time = new Date().getTime();

    markerSet = proj.getMarkerSet();
    if (cnvsAsPositions) {
      LocusSet<CNVariant> tmp = new LocusSet<CNVariant>(cnvs, true, proj.getLog()) {

        /**
         *
         */
        private static final long serialVersionUID = 1L;

      };
      positions = tmp.getStartsAndStopsByChromosome();
    } else {
      positions = markerSet.getPositionsByChr();
    }
    counts = new int[positions.length][][];
    acceptableSNPs = new boolean[positions.length][];
    for (int i = 0; i < positions.length; i++) {
      counts[i] = new int[positions[i].length][2];
      acceptableSNPs[i] = new boolean[positions[i].length];
    }

    System.out.println(ext.getTime() + "\tDetermining acceptability...");
    for (int i = 0; i < cnvs.length; i++) {
      firstSNP = ArrayUtils.binarySearch(positions[cnvs[i].getChr()], cnvs[i].getStart(), true);
      lastSNP = ArrayUtils.binarySearch(positions[cnvs[i].getChr()], cnvs[i].getStop(), true);
      if (firstSNP == -1 || lastSNP == -1) {
        System.err.println("Error - could not locate start or stop position for "
                           + cnvs[i].getUCSClocation());
      } else {
        indel = cnvs[i].getCN() < 2 ? 0 : 1;
        for (int j = firstSNP; j <= lastSNP; j++) {
          counts[cnvs[i].getChr()][j][indel]++;
        }
      }
    }

    for (int i = 0; i < positions.length; i++) {
      for (int j = 0; j < positions[i].length; j++) {
        dels = counts[i][j][0];
        dups = counts[i][j][1];
        acceptableSNPs[i][j] = dels + dups >= totalRequired && dels >= delRequired
                               && dups >= dupRequired && dels + dups <= totalLimitedTo
                               && dels <= delLimitedTo && dups <= dupLimitedTo;
      }
    }

    ArrayList<CNVariant> acceptable = new ArrayList<>();

    System.out.println(ext.getTime() + "\tFiltering CNVs...");
    for (CNVariant cnv : cnvs) {
      firstSNP = ArrayUtils.binarySearch(positions[cnv.getChr()], cnv.getStart(), true);
      lastSNP = ArrayUtils.binarySearch(positions[cnv.getChr()], cnv.getStop(), true);
      indel = cnv.getCN() < 2 ? 0 : 1;

      if (firstSNP == -1 || lastSNP == -1) {
        accepted = false;
      } else {
        if (proportionOfProbesThatNeedToPassForFinalInclusion < 1.0) {
          countAcceptable = 0;
          for (int j = firstSNP; j <= lastSNP; j++) {
            if (acceptableSNPs[cnv.getChr()][j]) {
              countAcceptable++;
            }
          }
          accepted = (double) countAcceptable
                     / (double) (lastSNP - firstSNP
                                 + 1) > proportionOfProbesThatNeedToPassForFinalInclusion;
        } else {
          index = firstSNP;
          accepted = false;
          while (!accepted && index <= lastSNP) {
            if (acceptableSNPs[cnv.getChr()][index]) {
              accepted = true;
            }
            index++;
          }
        }
      }

      if (accepted) {
        acceptable.add(cnv);
      }
    }

    System.out.println("Finished in " + ext.getTimeElapsed(time));
    return acceptable.toArray(new CNVariant[acceptable.size()]);
  }

  public static void filterBasedOnNumberOfCNVsAtLocus(Project proj, String filein, String fileout,
                                                      int totalRequired, int delRequired,
                                                      int dupRequired, int totalLimitedTo,
                                                      int delLimitedTo, int dupLimitedTo,
                                                      double proportionOfProbesThatNeedToPassForFinalInclusion) {

    long time = new Date().getTime();

    MarkerDetailSet markerSet = proj.getMarkerSet();
    Multiset<Marker> delCounts = HashMultiset.create();
    Multiset<Marker> dupCounts = HashMultiset.create();

    System.out.println(ext.getTime() + "\tLoading plink file...");
    CNVariant[] cnvs = CNVariant.loadPlinkFile(filein);

    System.out.println(ext.getTime() + "\tDetermining acceptability...");
    SetMultimap<GenomicPosition, Marker> genomicPositionMap = markerSet.getGenomicPositionMap();
    for (int i = 0; i < cnvs.length; i++) {
      byte chr = cnvs[i].getChr();
      int start = cnvs[i].getStart();
      int stop = cnvs[i].getStop();
      Set<Marker> firstSnpMatches = genomicPositionMap.get(new GenomicPosition(chr, start));
      Set<Marker> lastSnpMatches = genomicPositionMap.get(new GenomicPosition(chr, stop));
      if (firstSnpMatches.isEmpty() || lastSnpMatches.isEmpty()) {
        System.err.println("Error - could not locate start or stop position for "
                           + cnvs[i].getUCSClocation());
      } else {
        Set<Marker> markersInCNV = markerSet.getMarkersInSeg(cnvs[i]);
        if (cnvs[i].getCN() < 2) {
          delCounts.addAll(markersInCNV);
        } else {
          dupCounts.addAll(markersInCNV);
        }
      }
    }
    Set<Marker> acceptables = Sets.newHashSet();
    for (Marker marker : markerSet.getMarkers()) {
      int dels = delCounts.count(marker);
      int dups = dupCounts.count(marker);
      if (dels + dups >= totalRequired && dels >= delRequired && dups >= dupRequired
          && dels + dups <= totalLimitedTo && dels <= delLimitedTo && dups <= dupLimitedTo) {
        acceptables.add(marker);
      }
    }

    System.out.println(ext.getTime() + "\tFiltering CNVs...");
    try (PrintWriter writer = Files.openAppropriateWriter(fileout)) {
      writer.println(ArrayUtils.toStr(CNVariant.PLINK_CNV_HEADER));
      for (CNVariant cnv : cnvs) {
        byte chr = cnv.getChr();
        int start = cnv.getStart();
        int stop = cnv.getStop();
        Set<Marker> firstSnpMatches = genomicPositionMap.get(new GenomicPosition(chr, start));
        Set<Marker> lastSnpMatches = genomicPositionMap.get(new GenomicPosition(chr, stop));
        final boolean accepted;
        if (firstSnpMatches.isEmpty() || lastSnpMatches.isEmpty()) {
          accepted = false;
        } else {
          Set<Marker> markersInCNV = markerSet.getMarkersInSeg(cnv);
          if (proportionOfProbesThatNeedToPassForFinalInclusion < 1.0) {
            int countAcceptable = Sets.intersection(acceptables, markersInCNV).size();
            accepted = countAcceptable
                       / (double) (markersInCNV.size()) > proportionOfProbesThatNeedToPassForFinalInclusion;
          } else {
            accepted = acceptables.containsAll(markersInCNV);
          }
        }

        if (accepted) {
          writer.println(cnv.toPlinkFormat());
        }
      }
    } catch (IOException e) {
      System.err.println("Error writing to " + fileout);
      e.printStackTrace();
    }

    System.out.println("Finished in " + ext.getTimeElapsed(time));
  }

  public static boolean inOneOfTheseRegions(CNVariant cnv, Segment[] regions) {
    for (Segment region : regions) {
      if (cnv.significantOverlap(region)) {
        return true;
      }
    }
    return false;
  }

  // public static boolean spansCentromereMidpoint(CNVariant cnv, Segment[] midpoints) {
  // for (int i = 0; i<midpoints.length; i++) {
  // if (cnv.overlaps(midpoints[i])) {
  // return true;
  // }
  // }
  // return false;
  // }

  public static String getFilename(String root, int delSize, int dupSize, int number, double score,
                                   int commonInOutOrIgnore) {
    return root + "_" + (delSize == dupSize ? delSize + "kb" : delSize + "," + dupSize + "kb") + "_"
           + number + "SNP_" + score + "_"
           + (commonInOutOrIgnore == COMMON_IN ? "isCNP"
                                               : (commonInOutOrIgnore == COMMON_OUT ? "notCNP"
                                                                                    : "CNPstatusIgnored"))
           + ".cnv";

  }

  public static void union(String firstCNVfile, String secondCNVfile, String outputfile) {
    PrintWriter writer;
    CNVariant[] list1, list2;
    int count;
    boolean unique;

    list1 = CNVariant.loadPlinkFile(firstCNVfile);
    list2 = CNVariant.loadPlinkFile(secondCNVfile);

    try {
      writer = Files.openAppropriateWriter(outputfile);
      for (CNVariant element : list1) {
        writer.println(element.toPlinkFormat());
      }
      for (CNVariant element : list2) {
        count = 0;
        unique = true;
        while (unique && count < list1.length) {
          if (list1[count].equalsIncludingIndividual(element)) {
            unique = false;
          }
          count++;
        }
        if (unique) {
          writer.println(element.toPlinkFormat());
        }
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + outputfile);
      e.printStackTrace();
    }
  }

  public static void stdFilters(String dir, String filename, boolean makeUCSCtracks, String pedfile,
                                int build) {
    String root;
    Logger log;

    log = new Logger();
    root = ext.rootOf(filename);
    FilterCalls.filter(dir, filename, root + "_allAbove10.0_unfiltered.cnv", new int[] {1, 1},
                       new int[] {1, 1}, new int[] {1, 1}, 10, null, COMMON_IGNORED, pedfile, true,
                       null, makeUCSCtracks, build, log);
    FilterCalls.filter(dir, filename, root + "_allAbove10.0_filtered.cnv", new int[] {1, 1},
                       new int[] {1, 1}, new int[] {1, 1}, 10, DEFAULT_PROBLEMATIC_REGIONS,
                       COMMON_IGNORED, pedfile, true, null, makeUCSCtracks, build, log);
    FilterCalls.filter(dir, filename, root + "_ConservativeCalls.cnv", new int[] {100, 100},
                       new int[] {100, 100}, new int[] {20, 20}, 10, DEFAULT_PROBLEMATIC_REGIONS,
                       COMMON_IGNORED, pedfile, true, null, makeUCSCtracks, build, log);

    FilterCalls.filterOnSegments(dir, root + "_allAbove10.0_filtered.cnv",
                                 root + "_allAbove10.0_filtered_inGenes.cnv",
                                 Aliases.getPathToFileInReferenceDirectory(GeneSet.REFSEQ_SEGS,
                                                                           true, log),
                                 false);
    FilterCalls.filterOnSegments(dir, root + "_allAbove10.0_filtered.cnv",
                                 root + "_allAbove10.0_filtered_inExons.cnv",
                                 Aliases.getPathToFileInReferenceDirectory(GeneSet.REFSEQ_EXONS,
                                                                           true, log),
                                 false);

    // FilterCalls.filterBasedOnNumberOfCNVsAtLocus(new
    // Project(cnv.Launch.getDefaultDebugProjectFile(), false),
    // dir+root+"_0kb_5SNP_10.0_CNPstatusIgnored.cnv", dir+root+"_0kb_5SNP_10.0_3dels_LT2dups.cnv",
    // 0, 3, 0, Integer.MAX_VALUE, 1);
    // FilterCalls.filterBasedOnNumberOfCNVsAtLocus(new
    // Project(cnv.Launch.getDefaultDebugProjectFile(), false),
    // dir+root+"_0kb_5SNP_10.0_CNPstatusIgnored.cnv", dir+root+"_0kb_5SNP_10.0_3dups_LT2dels.cnv",
    // 0, 0, 3, 1, Integer.MAX_VALUE);
    // FilterCalls.filterBasedOnNumberOfCNVsAtLocus(new
    // Project(cnv.Launch.getDefaultDebugProjectFile(), false),
    // dir+root+"_0kb_5SNP_10.0_CNPstatusIgnored.cnv", dir+root+"_0kb_5SNP_10.0_2dels_2dups.cnv", 0,
    // 2, 2, Integer.MAX_VALUE, Integer.MAX_VALUE);
    // FilterCalls.filterBasedOnNumberOfCNVsAtLocus(new
    // Project(cnv.Launch.getDefaultDebugProjectFile(), false),
    // dir+root+"_0kb_5SNP_10.0_CNPstatusIgnored.cnv", dir+root+"_0kb_5SNP_10.0_3anythings.cnv", 3,
    // 0, 0, Integer.MAX_VALUE, Integer.MAX_VALUE);
    // FilterCalls.filterBasedOnNumberOfCNVsAtLocus(new
    // Project(cnv.Launch.getDefaultDebugProjectFile(), false),
    // dir+root+"_0kb_5SNP_10.0_CNPstatusIgnored.cnv", dir+root+"_0kb_5SNP_10.0_5anythings.cnv", 5,
    // 0, 0, Integer.MAX_VALUE, Integer.MAX_VALUE);
    // FilterCalls.filterOnSegments(dir+root+"_0kb_5SNP_10.0_3anythings.cnv",
    // dir+root+"_CommonInGenes.cnv", GeneSet.DIRECTORY+GeneSet.REFSEQ_SEGS, true);

    // FilterCalls.union(dir+root+"_0kb_5SNP_10.0_3anythings.cnv",
    // dir+root+"_100kb_20SNP_10.0_CNPstatusIgnored.cnv", dir+"unionOfConservativeAndCommon.cnv");
    // FilterCalls.filterOnSegments(dir+"unionOfConservativeAndCommon.cnv",
    // dir+"unionOfConservativeAndCommonInGenes.cnv", GeneSet.DIRECTORY+GeneSet.REFSEQ_SEGS, true);
  }

  public static void fromParameters(String filename, Logger log) {
    List<String> params;
    String problematicRegionsLocation = "";

    problematicRegionsLocation = Aliases.getPathToFileInReferenceDirectory("problematicRegions_hg19.dat",
                                                                           false, log);
    if (problematicRegionsLocation == null) {
      problematicRegionsLocation = "problematicRegions_hg19.dat";
    }

    params = Files.parseControlFile(filename, "filterCNVs",
                                    new String[] {"dir=", "in=penncnv.cnv", "out=conf15used.cnv",
                                                  "# minimum size of heterozygous deletions / duplications (in kb):",
                                                  "delSize=0", "dupSize=0",
                                                  "# minimum size of homozygous deletions / duplications (in kb):",
                                                  "hDelSize=0", "hDupSize=0",
                                                  "# minimum number of heterozygous SNPs:",
                                                  "number=15",
                                                  "# minimum number of probes for homozygous deletions:",
                                                  "hNumber=15", "minScore=10.0",
                                                  "filterFile=" + problematicRegionsLocation,
                                                  "# pedfile to be used as a filter:",
                                                  "ped=plink.fam",
                                                  "# if CNV spans centromere, break into two spanning actual markers",
                                                  "breakCentromere=true",
                                                  "# if breakCentromere is set to true, custom marker set to determine the last and first marker of the centromeres",
                                                  "markerFile=",
                                                  "# make a UCSC track (.bed file) as well",
                                                  "ucsc=true", "",
                                                  "# ALTERNATIVELY, in addition to dir/in/out and ignoring all other filters, you can",
                                                  "# keep only CNVs overlapping these segments (simply uncomment the following argument):",
                                                  "#segs=gene_region.dat",
                                                  "# exclude instead of include:",
                                                  "#excludeSegsInstead=true", "",
                                                  "# ALTERNATIVELY, in addition to dir/in/out, and ignoring all other filters, you can",
                                                  "#  filter CNVs by frequency of occurrence:",
                                                  "# Project Properties File",
                                                  "# proj=project.properties",
                                                  "# FAM file, or list or IDs (used to count number of individuals)",
                                                  "# famFile=plink.fam",
                                                  "# Percentage of total required before reporting CNVs at a locus",
                                                  "# totalRequired=0",
                                                  "# Percentage of deletions required before reporting CNVs at a locus",
                                                  "# delRequired=0",
                                                  "# Percentage of duplications required before reporting CNVs at a locus",
                                                  "# dupRequired=0",
                                                  "# Percentage limit of total before excluding CNVs at a locus",
                                                  "# totalLimitedTo=1",
                                                  "# Percentage limit of deletions before excluding CNVs at a locus",
                                                  "# delLimitedTo=1",
                                                  "# Percentage limit of duplications before excluding CNVs at a locus",
                                                  "# dupLimitedTo=1",
                                                  "# Proportion of probes that need to pass for final inclusion ",
                                                  "# proportion=0.5", "",

                                    }, log);

    if (params != null) {
      params.add("log=" + log.getFilename());
      main(ArrayUtils.toStringArray(params));
    }
  }

  // dir="D:/SIDS and IQ/" in="penncnv.cnv" out="IQ/c10_p3,15.cnv" delSize=15 hDelSize=4 dupSize=15
  // hDupSize=15 filterFile="N:/statgen/NCBI/problematicRegions_hg19.dat" breakCentromere=true
  // in=D:/data/gedi_gwas/regions.cnv list=D:/data/gedi_gwas/data/merged_split.cnv
  // out=D:/data/gedi_gwas/regions_major.cnv minScore=10 number=15 -stats
  // famFile=D:/data/gedi_gwas/gedi_gwas_plink.fam
  // in=D:/data/ny_registry/new_york/data/cnvlist_puv.cnv
  // list=D:/data/ny_registry/new_york/penncnvShadow/combinedAX.cnv
  // out=D:/data/ny_registry/new_york/penncnvShadow/cnvstats_puv_repeat.cnv minScore=10 number=15
  // -stats

  public static void main(String[] args) {
    int numArgs = args.length;
    int delSize = DEFAULT_MIN_SIZE_KB;
    int dupSize = DEFAULT_MIN_SIZE_KB;
    int hDelSize = DEFAULT_MIN_SIZE_KB;
    int hDupSize = DEFAULT_MIN_SIZE_KB;
    int number = DEFAULT_MIN_NUM_SNPS;
    int hNumber = DEFAULT_MIN_NUM_SNPS;
    int build = 37;
    int inOutIgnore = COMMON_IGNORED;
    double score = DEFAULT_MIN_SCORE;
    double overlap = 0.5;
    String dir = "";
    String in = "conf.cnv";
    String filenameOfProblematicRegions = null;
    String out = getFilename("conf", delSize, dupSize, number, score, inOutIgnore);
    String pedfile = null;
    String segs = "";
    String markerSetFilenameToBreakUpCentromeres = null;
    String listFile = null;
    String excludeFile = null;
    String projName = null;
    String logfile = null;
    String[] listFiles = null;
    boolean excludeSegs = false;
    boolean standards = false;
    boolean tracks = false;
    boolean breakCent = DEFAULT_BREAK_CENTROMERE;
    boolean exclude = false;
    boolean stats = false;
    boolean merge = false;
    float mergeFactor = DEFAULT_CLEAN_FACTOR;

    double totalRequired, delRequired, dupRequired, totalLimitedTo, delLimitedTo, dupLimitedTo,
        proportionOfProbesThatNeedToPassForFinalInclusion;
    totalRequired = delRequired = dupRequired = totalLimitedTo = delLimitedTo = dupLimitedTo = proportionOfProbesThatNeedToPassForFinalInclusion = 0.0;
    String famFile = null;

    String usage = "cnv.analysis.FilterCalls requires 2+ arguments\n"
                   + "   (1) directory (i.e. dir=" + dir + " (default))\n"
                   + "   (2) file in (i.e. in=" + in + " (default))\n"
                   + "   (3) file out (i.e. out=" + out + " (default))\n"
                   + "   (4) minimum size of a deletion (in kb) (i.e. delSize=" + delSize
                   + " (default))\n" + "   (5) minimum size of a duplication (in kb) (i.e. dupSize="
                   + dupSize + " (default))\n"
                   + "   (6) (Optional) minimum size of a homozygous deletion (in kb) (i.e. hDelSize="
                   + hDelSize + " (default))\n"
                   + "   (7) (Optional) minimum size of a homozygous duplication (in kb) (i.e. hDupSize="
                   + hDelSize + " (default))\n"
                   + "   (8) minimum number of heterozygous SNPs (i.e. number=" + number
                   + " (default))\n"
                   + "   (9) (Optional) minimum number of homozygous SNPs (i.e. hNumber=" + number
                   + " (default))\n" + "   (9) minimum score (i.e. minScore=" + score
                   + " (default))\n"
                   + "   (10) filter out cnvs in known problematicRegions (i.e. filterFile="
                   + filenameOfProblematicRegions + " (default))\n"
                   + "   (11) pedfile to use as a filter (i.e. ped=" + pedfile + " (default))\n"
                   + "   (12) if CNV spans centromere, break into two spanning actual markers (i.e. breakCentromere="
                   + breakCent + " (default))\n"
                   + "   (13) build of the genome to use for centromeres (i.e. build=" + build
                   + " (default))\n"
                   + "   (14) custom marker set to determine the last and first marker of the centromeres (i.e. markerFile=plink.bim (not the default))\n"
                   + "   (15) make UCSC track as well (i.e. ucsc=true (default))\n" + "  OR\n"
                   + "   (1) keep only CNVs overlapping these segments (i.e. segs=gene_region.dat (not the default))\n"
                   + "   (2) exclude instead of include (i.e. excludeSegsInstead=false (default))\n"
                   + "  OR\n"
                   + "   (1) perform all standard filters (i.e. -std (not the default))\n"
                   + "   (2) make UCSC tracks as well (i.e. ucsc=false (default))\n" + "  OR\n"
                   + "   (1) directory (i.e. dir=" + dir + " (default))\n"
                   + "   (2) file in (i.e. in=" + in + " (default))\n"
                   + "   (3) file out (i.e. out=" + out + " (default))\n"
                   + "   (4) Project Properties File (i.e. proj=project.properties (not the default))\n"
                   + "   (5) FAM file, or list of IDs (used to count number of individuals) (i.e. famFile=plink.fam (not the default))\n"
                   + "   (6) Percentage of total required before reporting CNVs at a locus (i.e. totalRequired=0 (default))\n"
                   + "   (7) Percentage of deletions required before reporting CNVs at a locus (i.e. delRequired=0 (default))\n"
                   + "   (8) Percentage of duplications required before reporting CNVs at a locus (i.e. dupRequired=0 (default))\n"
                   + "   (9) Percentage limit of total before excluding CNVs at a locus (i.e. totalLimitedTo=1 (default))\n"
                   + "  (10) Percentage limit of deletions before excluding CNVs at a locus (i.e. delLimitedTo=1 (default))\n"
                   + "  (11) Percentage limit of duplications before excluding CNVs at a locus (i.e. dupLimitedTo=1 (default))\n"
                   + "  (12) Proportion of probes that need to pass for final inclusion (i.e. proportion=0.5 (default))\n"
                   + "  OR\n" + "   (1) directory (i.e. dir=" + dir + " (default))\n"
                   + "   (2) file in (i.e. in=" + in + " (default))\n"
                   + "   (3) file out (i.e. out=" + out + " (default))\n"
                   + "   (4) remove all CNVs not belonging to a group of sample IDs (i.e. list=/path/to/list.txt (not the default))\n"
                   + "   (5) (Optional) remove all CNVs that overlap other CNVs not belonging to the given list of sample IDs (i.e. -exclude (not the default))"
                   + "";

    System.out.println();
    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("in=")) {
        in = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("out=")) {
        out = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("delSize=")) {
        delSize = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("dupSize=")) {
        dupSize = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("hDelSize=")) {
        hDelSize = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("hDupSize=")) {
        hDupSize = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("number=")) {
        number = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("hNumber=")) {
        hNumber = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("minScore=")) {
        score = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("cnps=")) {
        inOutIgnore = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("filterFile=")) {
        filenameOfProblematicRegions = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("-std")) {
        standards = true;
        numArgs--;
      } else if (arg.startsWith("ucsc=")) {
        tracks = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("segs=")) {
        segs = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("breakCentromere=")) {
        breakCent = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("markerFile=")) {
        markerSetFilenameToBreakUpCentromeres = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("excludeSegsInstead=")) {
        excludeSegs = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("ped=")) {
        pedfile = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("build=")) {
        build = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("overlap=")) {
        overlap = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("-exclude")) {
        exclude = true;
        numArgs--;
      } else if (arg.startsWith("-stats")) {
        stats = true;
        numArgs--;
      } else if (arg.startsWith("-merge")) {
        merge = true;
        numArgs--;
      } else if (arg.startsWith("mergeFactor=")) {
        mergeFactor = ext.parseFloatArg(arg);
        numArgs--;
      } else if (arg.startsWith("list=")) {
        listFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("lists=")) {
        listFiles = arg.split("=")[1].split(",");
        numArgs--;
      } else if (arg.startsWith("excludeFile=")) {
        excludeFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("proj=")) {
        projName = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("famFile=")) {
        famFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("totalRequired=")) {
        totalRequired = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("delRequired=")) {
        delRequired = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("dupRequired=")) {
        dupRequired = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("totalLimitedTo=")) {
        totalLimitedTo = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("delLimitedTo=")) {
        delLimitedTo = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("dupLimitedTo=")) {
        dupLimitedTo = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("proportion=")) {
        proportionOfProbesThatNeedToPassForFinalInclusion = ext.parseDoubleArg(arg);
        numArgs--;
      } else {
        System.err.println("Error - don't know what to do with argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      return;
    }

    // dir = "D:/data/GEDI/penn_results/custom_gediBoth/";
    // filterBasedOnNumberOfCNVsAtLocus(new Project(cnv.Launch.getDefaultDebugProjectFile(), false),
    // dir+"conf15_usedFiltered.cnv", dir+"conf15_usedFilteredRare.cnv", 0, 0, 0, 275, 275, 275,
    // 0.50);
    // System.exit(1);

    // filter("D:/data/GEDI/penn_results/custom_gediBoth/", "penncnv.cnv", "conf1checkers.cnv", 0,
    // 0, 1, 0, null, -1, "plink.fam", true, null, true, 37, new Logger());
    // System.exit(1);

    // FilterCalls.filterOnSegments("D:/data/GEDI/global/homoDelsOverlappingGenesOnly/", "conf.cnv",
    // "conf_overlappingGenes.cnv", GeneSet.DIRECTORY+GeneSet.REFSEQ_SEGS, false);
    // FilterCalls.filterOnSegments("D:/data/GEDI/penn_results/custom_gediBoth/conf15_usedFilteredRare/homoDels/",
    // "conf.cnv", "conf_overlappingGenes.cnv", GeneSet.DIRECTORY+GeneSet.REFSEQ_SEGS, false);
    // System.exit(1);

    // breakCent = true;
    // out = "noCentromeric.cnv";

    /*
     * dir, in, out, segs, excludeSegs --> filterOnSegments listFile, group, dir, in, out -->
     * filterForAllCNVsSharedInGroup listFile, stats, in, out, score, number --> filterList
     * listFile, dir, in, out, excludeFile, exclude --> filterForGroupCNVs listFiles, in, out,
     * score, number --> filterLists projName, dir, in --> CNVStats excludeFile, dir, in, out -->
     * filterExclusions proj, in, out, req[], lim[], proportion (dir, in, out) --> (filter) delSize
     * hDelSize dupSize hDupSize number hNumber score filenameOfProblematicRegions pedfile breakCent
     * markerSetFilenameToBreakUpCentromeres tracks build logfile --> (filterOutCommonCNVs) -common
     * pct --> (filterForAllCNVsSharedInGroup) -group listFile --> (filterForGroupCNVs) listFile
     * excludeFile -exclude --> (filterExclusions) excludeFile (dir, in) --> (stdFilters) -std
     * tracks pedfile build --> (CNVStats) projName (in, out) --> (filterList) -stats listFile score
     * number --> (filterLists) -stats listFiles score number --> (cleanCNVs) -clean
     */
    try {
      // FilterCalls.filterOnSegments(dir+"conf_100kb_20SNP_10.0_CNPstatusIgnored.cnv",
      // dir+"ConservativeGeneCentric.cnv", GeneSet.DIRECTORY+GeneSet.REFSEQ_SEGS, true);
      // MakeUCSCtrack.makeTrack(dir, "ConservativeGeneCentric.cnv");
      if (standards) {
        stdFilters(dir, in, tracks, pedfile, build);
      } else if (!segs.equals("")) {
        filterOnSegments(dir, in, out, segs, excludeSegs);
      } else if (listFile != null && stats) {
        cnvStats(in, new String[] {listFile}, /* famFile, */out, score, number, overlap);
      } else if (listFiles != null && stats) {
        // in=D:/data/ny_registry/new_york/data/cnvlist.cnv
        // list=D:/data/ny_registry/new_york/penncnvShadow/penncnv.cnv
        // out=D:/data/ny_registry/new_york/penncnvShadow/cnvstats_auto.cnv minScore=10 number=15
        // -stats
        // -stats in= out= list= famFile= minScore=10 number=15 overlap=.50
        cnvStats(in, listFiles, out, score, number, overlap);
      } else if (projName != null) {
        Project proj = new Project(projName);
        if (dir != null && !"".equals(dir)) {
          CNVStats(proj, dir, in);
        } else if (merge) {
          mergeCNVs(proj, in, out, mergeFactor);
        } else {
          int famCnt = Files.countLines(famFile, 0);
          filterBasedOnNumberOfCNVsAtLocus(proj, in, out, (int) (totalRequired * famCnt * 2.0),
                                           (int) (delRequired * famCnt * 2.0),
                                           (int) (dupRequired * famCnt * 2.0),
                                           (int) (totalLimitedTo * famCnt * 2.0),
                                           (int) (delLimitedTo * famCnt * 2.0),
                                           (int) (dupLimitedTo * famCnt * 2.0),
                                           proportionOfProbesThatNeedToPassForFinalInclusion);
        }
        // proj=D:/projects/NY_Registry_Combo_Data.properties
        // dir=D:/data/ny_registry/new_york/penncnvShadow/ in=penncnv
      } else if (excludeFile != null) {
        filterExclusions(dir, in, out, excludeFile, exclude);
      } else if (merge) {
        mergeCNVs(in, out, mergeFactor, null);
      } else {
        filterCNVs(dir, in, out, new int[] {delSize, hDelSize}, new int[] {dupSize, hDupSize},
                   new int[] {number, hNumber}, score, filenameOfProblematicRegions,
                   DEFAULT_COMMON_IN_OUT_OR_IGNORED, pedfile, breakCent,
                   markerSetFilenameToBreakUpCentromeres, tracks, build, new Logger(logfile));
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}

/*
 * in=D:/data/ny_registry/new_york/data/cnvlist_puv.cnv
 * lists=D:/data/ny_registry/new_york/penncnvShadow/sexSpecific/male/recodedM.cnv
 * out=D:/data/ny_registry/new_york/penncnvShadow/cnvstats_puv_repeat.cnv minScore=10 number=15
 * -stats in=D:/data/ny_registry/new_york/data/cnvlist_puv.cnv
 * lists=D:/data/ny_registry/new_york/penncnvShadow/sexSpecific/male/recodedM.cnv,D:/data/
 * ny_registry/new_york/penncnvShadow/sexSpecific/female/recodedF.cnv
 * out=D:/data/ny_registry/new_york/penncnvShadow/cnvstats_puv_repeat.cnv minScore=0 number=0 -stats
 * ---------------------- in=D:/data/ny_registry/new_york/data/cnvlist_puv.cnv
 * lists=D:/data/ny_registry/new_york/penncnvShadow/sexSpecific/male/recodedM.cnv,D:/data/
 * ny_registry/new_york/penncnvShadow/sexSpecific/female/recodedF.cnv
 * out=D:/data/ny_registry/new_york/penncnvShadow/cnvstats_puv_all.cnv minScore=0 number=0 -stats
 * in=D:/data/ny_registry/new_york/data/cnvlist_puv.cnv
 * lists=D:/data/ny_registry/new_york/penncnvShadow/sexSpecific/male/recodedM.cnv,D:/data/
 * ny_registry/new_york/penncnvShadow/sexSpecific/female/recodedF.cnv
 * out=D:/data/ny_registry/new_york/penncnvShadow/cnvstats_puv_filtered.cnv minScore=10 number=15
 * -stats
 */
