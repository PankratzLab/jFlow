package org.genvisis.cnv.park;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.Hashtable;
import java.util.List;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.stats.Maths;

public class CompareCalls_dev {
  public static final String DEFAULT_ROOT = "C:/data/pennComp/";
  public static final String DEFAULT_MARKER_FILE = "affygw6.hg18.pfb";
  public static final String DEFAULT_CNV_FILE = "all_gw6.cnv";
  public static final String[] DEFAULT_COMP_LISTS = {"roots1.txt", "dups1.txt"};
  public static final String DEFAULT_COMP_LIST = "rootsANDdubs.comp.txt";
  public static final String DEFAULT_LRR_LOOKUP = "LRR_full_samps.txt";
  public static final Double DEFAULT_LRR_FILTER = .30;
  public static final Double DEFAULT_CONF_FILTER = 10.0;
  public static final Integer DEFAULT_NUM_MARKERS_FILTER = 20;
  public static final Double DEFAULT_BP_FRACTION = 0.3;
  public static final Double DEFAULT_NUM_MARKER_FRACTION = 0.3;

  private static boolean alreadyDefinedID(Hashtable<String, Hashtable<String, Integer>> defineCompHash,
                                          String[] line) {
    boolean alreadyDefined = false;
    for (String element : line) {
      // same id not allowed in same line
      if (defineCompHash.containsKey(element)) {
        System.err.println("Warning - duplicate IDs were detectected in the comparision file: "
                           + element + " was seen twice, only the first comparison will be used ");
        alreadyDefined = true;
      }
    }
    return alreadyDefined;
  }

  private static void assignMax(Hashtable<String, Hashtable<String, Integer>> defineCompHash,
                                int maxNumComparisions) {
    defineCompHash.put("MAX_COMP_NUMBER", new Hashtable<String, Integer>());
    defineCompHash.get("MAX_COMP_NUMBER").put("MAX_COMP_NUMBER", maxNumComparisions);
  }

  // computes average percent match for each copy number, and a total for all copy numbers
  private static double[] averageCNPercents(double[] goodCalls, int[] cnNumbers) {
    double[] averageCNPercent = new double[6];
    for (int i = 0; i < goodCalls.length; i++) {
      averageCNPercent[i] = goodCalls[i] / cnNumbers[i];
    }
    return averageCNPercent;
  }

  // return is used for parameter iteration
  private static double[] binIt(double startVal, double stopVal, int numBins) {
    double[] values = new double[numBins + 1];
    double inc = getIncrement(startVal, stopVal, numBins);
    for (int i = 0; i < numBins + 1; i++) {
      values[i] = (inc * i) + startVal;
    }
    return values;
  }

  // takes gap divided by proposed new length and compares to fraction desired
  private static boolean bpGapCheck(CNVariant thisCNV, CNVariant nextCNV, double bpFraction) {
    int gapLengthBP = nextCNV.getStart() - thisCNV.getStop() - 1;
    int proposedNewBPLength = nextCNV.getStop() - thisCNV.getStart() + 1;
    return testFrac(bpFraction, (gapLengthBP / proposedNewBPLength));
  }

  private static boolean checkAll(double bpFraction, double markerFraction,
                                  Hashtable<Byte, List<Integer>> markers, CNVariant thisCNV,
                                  CNVariant nextCNV) {
    return checkChr(thisCNV, nextCNV) && bpGapCheck(thisCNV, nextCNV, bpFraction)
           && checkCN(thisCNV, nextCNV)
           && markerGapCheck(thisCNV, nextCNV, markers, markerFraction);
  }

  private static boolean checkChr(CNVariant thisCNV, CNVariant nextCNV) {
    return thisCNV.getChr() == nextCNV.getChr();
  }

  private static boolean checkCN(CNVariant thisCNV, CNVariant nextCNV) {
    return thisCNV.getCN() == nextCNV.getCN();
  }

  // merges CNVs based on gap percentages in base pairs and number of markers, copy number aware
  public static String cleanCNVs(String rootDir, String cnvfile, String markerFile,
                                 double bpFraction, double markerFraction) {
    PrintWriter writer;
    String[] inds;
    String file;
    CNVariant[] cnvs;
    CNVariant[] fileCNVs;
    Hashtable<String, ArrayList<CNVariant>> allIndCNVs;
    int changecount = 0;
    int nochange = 0;
    int totalCount = 0;
    int newCalls = 0;
    fileCNVs = CNVariant.loadPlinkFile(cnvfile, false);
    inds = toStringArray(getIDList(fileCNVs));
    allIndCNVs = getindCNVsfromFile(fileCNVs);
    Hashtable<Byte, List<Integer>> markers = getMarkerLookup(rootDir, markerFile);

    file = cnvfile.replace(".cnv", ".clean.cnv");
    System.out.println("cleaning " + cnvfile);
    writer = Files.getAppropriateWriter(file);
    writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));

    for (String ind : inds) {
      CNVariant mergedCNV;
      cnvs = toCNVArray(allIndCNVs.get(ind));
      // TODO add number of CNVs as param?
      // if (cnvs.length > 100) {
      // continue;
      // }
      int startCombine = 0;
      boolean combining = false;
      for (int j = 0; j < cnvs.length; j++) {
        totalCount++;
        if (j + 1 < cnvs.length) {
          CNVariant thisCNV = cnvs[j];
          CNVariant nextCNV = cnvs[j + 1];
          // main check for proper gap length and CN
          if (checkAll(bpFraction, markerFraction, markers, thisCNV, nextCNV)) {
            if (!combining) {
              startCombine = j;
              combining = true;
              changecount++;
            }
            changecount++;
          } else {
            // merge any previous calls that can be combined, else print the current cnv call
            if (combining) {
              newCalls++;
              mergedCNV = mergeCNVs(cnvs, startCombine, j, markers);
              writer.println(mergedCNV.toPlinkFormat());
              combining = false;
            } else {
              nochange++;
              writer.println(thisCNV.toPlinkFormat());
            }
          }
        } else {
          // in case the last call for an individual was merged
          if (combining) {
            newCalls++;
            mergedCNV = mergeCNVs(cnvs, startCombine, j, markers);
            writer.println(mergedCNV.toPlinkFormat());
            combining = false;
          } else {
            nochange++;
            writer.println(cnvs[j].toPlinkFormat());
          }
        }
      }
    }
    System.out.println("Info - " + changecount + " individual calls were combined into " + newCalls
                       + " new calls,\t" + nochange + " Calls were Unchanged,\tA total "
                       + totalCount + " of Calls Examined from " + inds.length + " individuals");
    writer.close();
    return file;
  }

  // compare all combinations of files, where Ids are matched as defined in the comp file, copy
  // number aware comparisons
  public static CompareCalls_dev[] compare(String rootDir, String[] files, String compFile) {
    String[] inds1;
    String[] inds2;
    int numPassingAndPresent = 0;
    int[][] counts;
    String[][] comparedIDs;
    int[] exactMatches, sigOlapMatches, olapMatches;
    int[][] indCNVCalls;
    double[] indPercentMatches;
    double[] goodCalls;
    double averageCNPercent[];
    int[] cnNumbers;
    CNVariant[][] cnvs;
    CNVariant[][] fileCNVs;

    int[][] allPossibleCombinations = Maths.getIndicesForAllCombinations(files.length, 2);
    CompareCalls_dev[] comparedCalls = new CompareCalls_dev[allPossibleCombinations.length];
    Hashtable<Integer, ArrayList<String>> fileIDs;
    Hashtable<Integer, Hashtable<String, ArrayList<CNVariant>>> indFileCNVs;
    Hashtable<String, Hashtable<String, Integer>> defineCompHash =
        new Hashtable<String, Hashtable<String, Integer>>();
    defineCompHash = defineCompLists(rootDir, compFile);
    fileCNVs = getFileCNVs(files);
    fileIDs = getFileIDs(fileCNVs);
    indFileCNVs = getindCNVsfromFiles(fileCNVs);

    for (int i = 0; i < allPossibleCombinations.length; i++) {
      inds1 = toStringArray(fileIDs.get(allPossibleCombinations[i][0]));
      inds2 = toStringArray(fileIDs.get(allPossibleCombinations[i][1]));
      if (inds1.length != inds2.length || inds1.length == 0 || inds2.length == 0) {
        System.err.println("Info - there are " + inds1.length + " ids  in "
                           + files[allPossibleCombinations[i][0]] + " and " + inds2.length
                           + " ids in " + files[allPossibleCombinations[i][1]]);
        if (inds1.length == 0 || inds2.length == 0) {
          System.err.println("Info - skipping the comparison between "
                             + files[allPossibleCombinations[i][0]] + " and "
                             + files[allPossibleCombinations[i][1]]);
        }
      }
      try {
        int init = getMaxNum(inds1, inds2);
        exactMatches = new int[init];
        sigOlapMatches = new int[init];
        olapMatches = new int[init];
        comparedIDs = new String[2][init];
        indPercentMatches = new double[init];
        indCNVCalls = new int[2][init];
        cnNumbers = new int[6];
        goodCalls = new double[6];
        numPassingAndPresent = 0;

        for (int j = 0; j < inds1.length; j++) {
          for (String element : inds2) {
            if (defineCompHash.containsKey(inds1[j])
                && defineCompHash.get(inds1[j]).containsKey(element)) {
              cnvs = getCompareCNVs(inds1[j], element, allPossibleCombinations[i], indFileCNVs);
              counts = countMatches(cnvs);

              // checkCNVs_b_a(cnvs, counts);
              olapMatches[j] = counts[2][2];
              exactMatches[j] = counts[2][3];
              sigOlapMatches[j] = counts[2][4];
              comparedIDs[0][j] = inds1[j];
              comparedIDs[1][j] = element;
              indCNVCalls[0][j] = cnvs[0].length;
              indCNVCalls[1][j] = cnvs[1].length;
              // TODO john do copynumber aware call percent. OK I did it. now calculates copy number
              // specific
              // percent matches
              if (cnvs[0].length > 0 || cnvs[1].length > 0) {
                cnNumbers = getCNNumbers(cnvs[0], cnvs[1], cnNumbers);
                goodCalls = goodCallCN(cnvs[0], cnvs[1], counts, goodCalls);
                double callPercent = ((double) ((2 * exactMatches[j]) + sigOlapMatches[j])
                                      / (cnvs[0].length + cnvs[1].length));
                // TODO tabulate with CN presence only
                numPassingAndPresent++;
                indPercentMatches[j] = callPercent;

              } else {
                indPercentMatches[j] = Double.NaN;
              }
            }
          }
        }
        averageCNPercent = averageCNPercents(goodCalls, cnNumbers);
        System.err.println("Info - " + numPassingAndPresent
                           + " with valid comparison matched Ids were found and compared");
        comparedCalls[i] = new CompareCalls_dev(numPassingAndPresent, averageCNPercent, cnNumbers,
                                                files, exactMatches, sigOlapMatches, olapMatches,
                                                comparedIDs, indCNVCalls, indPercentMatches);
      } catch (Exception e) {
        System.err.println("Error comparing " + files[allPossibleCombinations[i][0]] + " and "
                           + files[allPossibleCombinations[i][1]]);
        e.printStackTrace();
      }
    }
    return comparedCalls;
  }

  // copy number aware, compare all v all
  private static int[][] countMatches(CNVariant[][] cnvs) {
    // CN states as defined here are ,0,1,2,3,4,5...where 2 is a total non-copy number aware count,
    // and 5 is for non matching (overlap,exact,sigoverlap)
    PrintWriter misswriter, makewriter;
    int[][] counts = new int[6][5];
    try {
      misswriter = new PrintWriter(new FileWriter("C:/data/pennComp/misses.cnv", true));
      makewriter = new PrintWriter(new FileWriter("C:/data/pennComp/makes.cnv", true));
      int match;
      int CN = 5;
      for (int a = 0; a < cnvs[0].length; a++) {
        match = 0;
        for (int b = 0; b < cnvs[1].length; b++) {
          if (cnvs[1][b].getCN() == cnvs[0][a].getCN()) {
            if (cnvs[0][a].equals(cnvs[1][b])) {
              match = 3;
              cnvs[1][b].setSource(99);
              CN = cnvs[1][b].getCN();
            } else if (match < 2 && cnvs[0][a].significantOverlap(cnvs[1][b])) {
              match = 4;
              CN = cnvs[1][b].getCN();
              // an overlap is not assumed in both directions
              // cnvs[1][b].setSource(99);
            } else if (match < 2 && cnvs[0][a].overlaps(cnvs[1][b])) {
              match = 2;
              CN = cnvs[1][b].getCN();
              cnvs[1][b].setSource(99);
            }
          }
        }
        if (match == 0) {
          misswriter.println(cnvs[0][a].toPlinkFormat());
        }
        if (match != 0) {
          makewriter.println(cnvs[0][a].toPlinkFormat());
        }
        counts[CN][match]++;
        counts[2][match]++;
      }
      CN = 5;
      for (int b = 0; b < cnvs[1].length; b++) {
        match = 1;
        for (int a = 0; a < cnvs[0].length; a++) {
          if (cnvs[1][b].getCN() == cnvs[0][a].getCN()) {
            if (cnvs[1][b].getSource() != 99 && cnvs[1][b].equals(cnvs[0][a])) {
              match = 3;
              CN = cnvs[1][b].getCN();
            } else if (match < 2 && cnvs[1][b].getSource() != 99
                       && cnvs[1][b].significantOverlap(cnvs[0][a])) {
              match = 4;
              CN = cnvs[1][b].getCN();
            } else if (match < 2 && cnvs[1][b].getSource() != 99
                       && cnvs[1][b].overlaps(cnvs[0][a])) {
              match = 2;
              CN = cnvs[1][b].getCN();
            }
          }
        }
        if (match == 1 && cnvs[1][b].getSource() != 99) {
          misswriter.println(cnvs[1][b].toPlinkFormat());
        }
        if (match != 1) {
          makewriter.println(cnvs[1][b].toPlinkFormat());
        }
        if (cnvs[1][b].getSource() != 99) {
          counts[CN][match]++;
          counts[2][match]++;
        }
      }
      misswriter.close();
      makewriter.close();
    } catch (IOException e) {
      System.err.println("Error creating " + "C:/data/pennComp/misses.cnv");
      e.printStackTrace();
    }
    return counts;
  }

  private static void defineComparisons(Hashtable<String, Hashtable<String, Integer>> defineCompHash,
                                        String[] line) {
    Hashtable<String, Integer> defined = new Hashtable<String, Integer>();
    for (int i = 0; i < line.length; i++) {
      for (int j = 0; j < line.length; j++) {
        if (i != j && line[i] == line[j]) {
          System.err.println("Warning - duplicate IDs were detectected in the comparision file: "
                             + line[i] + " in column " + i + "and " + line[j] + " in column " + j);
        }

        defined.put(line[j], j);
        defineCompHash.put(line[i], defined);
      }
    }
  }

  // defines comparisons for appropriate parsing and comparing
  private static Hashtable<String, Hashtable<String, Integer>> defineCompLists(String rootDir,
                                                                               String compFile) {
    Hashtable<String, Hashtable<String, Integer>> defineCompHash =
        new Hashtable<String, Hashtable<String, Integer>>();
    BufferedReader reader;
    String[] line;
    int maxNumComparisions = 0;
    System.out.println("HI");
    try {
      reader = new BufferedReader(new FileReader(rootDir + compFile));
      while (reader.ready()) {
        line = reader.readLine().trim().split("\t");
        if (line.length > maxNumComparisions) {
          maxNumComparisions = line.length;
        }
        // same id in multiple lines not allowed
        if (alreadyDefinedID(defineCompHash, line)) {
          continue;
        }

        defineComparisons(defineCompHash, line);
      }
      reader.close();
      assignMax(defineCompHash, maxNumComparisions);
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + rootDir + compFile
                         + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + rootDir + compFile + "\"");
      System.exit(2);
    }
    return defineCompHash;
  }

  // do both filtering and cleaning
  public static String[] filterAndClean(String rootDir, String cnvFile, String LRR_lookup,
                                        String compFile, double lrrFilter, double confFilter,
                                        int numMarkers, String markerFile, double bpFraction,
                                        double markerFraction) {
    String[] filteredFiles, cleanedFiles;
    filteredFiles =
        filterCNVs(rootDir, cnvFile, LRR_lookup, compFile, lrrFilter, confFilter, numMarkers);
    cleanedFiles = new String[filteredFiles.length];
    for (int i = 0; i < filteredFiles.length; i++) {
      cleanedFiles[i] =
          cleanCNVs(rootDir, filteredFiles[i], markerFile, bpFraction, markerFraction);
    }
    return cleanedFiles;
  }

  public static String[] filterCNVs(String rootDir, String cnvFile, String LRR_lookup,
                                    String compFile, double lrrFilter, double confFilter,
                                    int numMarkers) {
    Hashtable<String, Double> lrrs = new Hashtable<String, Double>();
    Hashtable<String, Hashtable<String, Integer>> defineCompHash =
        new Hashtable<String, Hashtable<String, Integer>>();
    BufferedReader reader;
    String[] line;

    defineCompHash = defineCompLists(rootDir, compFile);
    int maxRepNumber = defineCompHash.get("MAX_COMP_NUMBER").get("MAX_COMP_NUMBER");
    PrintWriter[] writers = new PrintWriter[maxRepNumber];
    String[] files = new String[maxRepNumber];
    int[] filteredCounts = new int[maxRepNumber];
    int[] includeCounts = new int[maxRepNumber];
    int totalCounts = 0;
    try {

      lrrs = getLRRs(rootDir, LRR_lookup);
      reader = new BufferedReader(new FileReader(rootDir + cnvFile));
      if (!ext.checkHeader(reader.readLine().trim().split("[\\s]+"), CNVariant.PLINK_CNV_HEADER,
                           false)) {
        reader.close();
        System.err.println("quitting comparison");
      } else {
        for (int i = 0; i < maxRepNumber; i++) {
          files[i] = rootDir + "intermediateFiles/Filtered_LRR_" + lrrFilter + "_conf_" + confFilter
                     + "_numMarkers_" + numMarkers + "_rep" + i + ".cnv";
          writers[i] = new PrintWriter(new FileWriter(files[i]));
          writers[i].println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
          filteredCounts[i] = 0;
        }
        while (reader.ready()) {
          line = reader.readLine().trim().split("\t");
          totalCounts++;
          if (defineCompHash.containsKey(line[0])) {
            int fileNumber = defineCompHash.get(line[0]).get(line[0]);
            // main check for filter parameters
            if (lrrs.get(line[0]) <= lrrFilter && Double.parseDouble(line[6]) >= confFilter
                && Integer.parseInt(line[7]) >= numMarkers) {
              writers[fileNumber].println(Array.toStr(line));
              includeCounts[fileNumber]++;
            } else {
              filteredCounts[fileNumber]++;
            }
          }
        }
        for (int i = 0; i < maxRepNumber; i++) {
          writers[i].close();
          int totalInCompFile = includeCounts[i] + filteredCounts[i];
          System.out.println("Info - Out of a total of " + totalCounts + " cnv calls in " + rootDir
                             + cnvFile + ", a total of " + totalInCompFile
                             + " calls had matched ids listed in " + compFile + " column " + i
                             + ": " + filteredCounts[i] + " were filtered out, " + includeCounts[i]
                             + " were included for comparison in " + files[i]);
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + rootDir + cnvFile
                         + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + rootDir + cnvFile + "\"");
      System.exit(2);
    }
    return files;

  }

  private static int[] getCNNumbers(CNVariant[] cnvs0, CNVariant[] cnvs1, int[] cnNumbers) {
    for (int i = 0; i < cnvs0.length; i++) {
      cnNumbers[cnvs0[i].getCN()]++;
      cnNumbers[2]++;
    }
    for (int i = 0; i < cnvs1.length; i++) {
      cnNumbers[cnvs1[i].getCN()]++;
      cnNumbers[2]++;
    }
    return cnNumbers;

  }

  private static CNVariant[][] getCompareCNVs(String ind1, String ind2, int[] combinations,
                                              Hashtable<Integer, Hashtable<String, ArrayList<CNVariant>>> indFileCNVs) {
    return new CNVariant[][] {toCNVArray(indFileCNVs.get(combinations[0]).get(ind1)),
                              toCNVArray(indFileCNVs.get(combinations[1]).get(ind2))};
  }

  private static String[] getComparisons(CompareCalls_dev comparedCall) {
    String[] comparisons = new String[2];
    comparisons[0] = comparedCall.getCnvFiles()[0].replaceFirst(".*rep", "rep");
    comparisons[1] = comparedCall.getCnvFiles()[1].replaceFirst(".*rep", "rep");
    return comparisons;
  }

  private static CNVariant[][] getFileCNVs(String[] files) {
    CNVariant[][] fileCNVs = new CNVariant[files.length][];
    for (int i = 0; i < files.length; i++) {
      fileCNVs[i] = CNVariant.loadPlinkFile(files[i], false);
    }
    return fileCNVs;
  }

  private static Hashtable<Integer, ArrayList<String>> getFileIDs(CNVariant[][] filesCNVs) {
    Hashtable<Integer, ArrayList<String>> fileIds = new Hashtable<Integer, ArrayList<String>>();
    for (int i = 0; i < filesCNVs.length; i++) {
      fileIds.put(i, getIDList(filesCNVs[i]));
    }
    return fileIds;
  }

  private static String getFileName(String rootDir, String compareType, double start, double stop,
                                    String[] comparisons) {
    String file;
    if (compareType.equals("lrr")) {
      file = rootDir + "test" + comparisons[0] + "_vs_" + comparisons[1] + "_" + compareType
             + "_Concordance_" + start + "_" + stop + "_numMarkers_" + DEFAULT_NUM_MARKERS_FILTER
             + "_Conf_" + DEFAULT_CONF_FILTER + ".concord";
    } else if (compareType.equals("conf")) {
      file = rootDir + "test" + comparisons[0] + "_vs_" + comparisons[1] + "_" + compareType
             + "_Concordance_" + start + "_" + stop + "_lrr_" + DEFAULT_LRR_FILTER + "_numMarkers_"
             + DEFAULT_NUM_MARKERS_FILTER + ".concord";
    } else {
      System.err.println("invalid invalid Comparison type");
      file = "NA";
    }
    return file;
  }

  private static String getFileName(String rootDir, String compareType, int start, int stop,
                                    String[] comparisons) {
    String file;
    if (compareType.equals("numMarkers")) {
      file = rootDir + "test" + comparisons[0] + "_vs_" + comparisons[1] + "_" + compareType
             + "_Concordance_" + start + "_" + stop + "_LRR_" + DEFAULT_LRR_FILTER + "Conf_"
             + DEFAULT_CONF_FILTER + ".concord";
    } else {
      System.err.println("invalid invalid Comparison type");
      file = "NA";
    }
    return file;
  }

  private static String getHeader(String iterType) {
    return iterType
           + "\tCN=any;avgPercent\tCN=any;numCallsAnalyzed\tCN=0;avgPercent\tCN=0;numCallsAnalyzed\tCN=1;avgPercent\tCN=1;numCallsAnalyzed\tCN=3;avgPercent\tCN=3;numCallsAnalyzed\tCN=4;avgPercent\tCN=4;numCallsAnalyzed\tPairs of Individuals PassingFilters";
  }

  public static ArrayList<String> getIDList(CNVariant[] afileCNVs) {
    ArrayList<String> al = new ArrayList<String>();
    Hashtable<String, Integer> tracker = new Hashtable<String, Integer>();
    for (int k = 0; k < afileCNVs.length; k++) {
      // TODO
      // change to ind id when i update input file
      if (!tracker.containsKey(afileCNVs[k].getFamilyID())) {
        tracker.put(afileCNVs[k].getFamilyID(), 1);
        al.add(afileCNVs[k].getFamilyID());
      }
    }
    return al;
  }

  private static double getIncrement(double startVal, double stopVal, int numBins) {
    return (stopVal - startVal) / numBins;
  }

  public static Hashtable<String, ArrayList<CNVariant>> getindCNVsfromFile(CNVariant[] fileCNVs) {
    Hashtable<String, ArrayList<CNVariant>> indFileCNVs;
    indFileCNVs = new Hashtable<String, ArrayList<CNVariant>>();
    for (int k = 0; k < fileCNVs.length; k++) {
      if (!indFileCNVs.containsKey(fileCNVs[k].getFamilyID())) {
        // TODO
        // change to ind id when i update input file
        indFileCNVs.put(fileCNVs[k].getFamilyID(), new ArrayList<CNVariant>());
      }
      indFileCNVs.get(fileCNVs[k].getFamilyID()).add(fileCNVs[k]);
    }
    return indFileCNVs;
  }

  private static Hashtable<Integer, Hashtable<String, ArrayList<CNVariant>>> getindCNVsfromFiles(CNVariant[][] fileCNVs) {
    Hashtable<Integer, Hashtable<String, ArrayList<CNVariant>>> indCNVs =
        new Hashtable<Integer, Hashtable<String, ArrayList<CNVariant>>>();
    for (int i = 0; i < fileCNVs.length; i++) {
      indCNVs.put(i, getindCNVsfromFile(fileCNVs[i]));
    }
    return indCNVs;
  }

  // used in filtering
  private static Hashtable<String, Double> getLRRs(String rootDir, String LRR_lookup) {
    Hashtable<String, Double> lrrs = new Hashtable<String, Double>();
    BufferedReader reader;
    String[] line;
    try {
      reader = new BufferedReader(new FileReader(rootDir + LRR_lookup));
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        if (lrrs.containsKey(line[0])) {
          System.err.println("Warning - duplicate samples found in " + rootDir + LRR_lookup
                             + " only using one of them");
        } else {
          lrrs.put(line[0], Double.parseDouble(line[1]));
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + rootDir + LRR_lookup
                         + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + rootDir + LRR_lookup + "\"");
      System.exit(2);
    }
    return lrrs;
  }

  private static Hashtable<Byte, List<Integer>> getMarkerLookup(String rootDir,
                                                                String markerPositionsFile) {
    Hashtable<Byte, List<Integer>> markers = new Hashtable<Byte, List<Integer>>();
    BufferedReader reader;
    byte chr;
    String[] line;
    try {
      reader = Files.getAppropriateReader(rootDir + markerPositionsFile);
      do {
        line = reader.readLine().trim().split("\t", -1);
      } while (reader.ready() && !line[0].startsWith("Name"));
      // System.out.println(line[0]);
      if (!reader.ready()) {
        // TODO add index factors
        System.err.println("DID not FIND THE INDEX FACTORS THAT NEED TO BE ADDED");
      }
      while (reader.ready()) {
        line = reader.readLine().trim().split("\t", -1);
        if (line[1].startsWith("X")) {
          chr = (byte) 23;
        } else if (line[1].equals("Y")) {
          chr = (byte) 24;
        } else if (line[1].equals("MT") || line[1].equals("0") || line[1].equals("0")) {
          continue;
        } else {
          chr = Byte.valueOf(line[1]);
          if (!markers.containsKey(chr)) {
            List<Integer> bpList = new ArrayList<Integer>();
            markers.put(chr, bpList);
          }
          markers.get(chr).add(Integer.parseInt(line[2]));
        }
      }
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + markerPositionsFile
                         + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + rootDir + markerPositionsFile + "\"");
      System.exit(2);
    }
    sortBPbyChromosome(markers);
    return markers;
  }

  private static int getMaxNum(String[] inds1, String[] inds2) {
    return Math.max(inds1.length, inds2.length);
  }

  private static String getPrintString(CompareCalls_dev comparedCall) {
    String printString = "NA";

    double[] avePercents = comparedCall.getAveragePercent();
    int[] totalCallsAnalyzed = comparedCall.getTotalCallsAnalyzed();
    printString = avePercents[2] + "\t" + totalCallsAnalyzed[2] + "\t" + avePercents[0] + "\t"
                  + totalCallsAnalyzed[0] + "\t" + avePercents[1] + "\t" + totalCallsAnalyzed[1]
                  + "\t" + avePercents[3] + "\t" + totalCallsAnalyzed[3] + "\t" + avePercents[4]
                  + "\t" + totalCallsAnalyzed[4] + "\t" + comparedCall.getnumPassingAndPresent();

    return printString;
  }

  private static int getStartIndex(CNVariant CNV, Hashtable<Byte, List<Integer>> markers) {
    return markers.get(CNV.getChr()).indexOf(CNV.getStart());
  }

  private static int getStopIndex(CNVariant CNV, Hashtable<Byte, List<Integer>> markers) {
    return markers.get(CNV.getChr()).indexOf(CNV.getStop());
  }

  private static PrintWriter[] getWriters(String rootDir, CompareCalls_dev[] comparedCalls,
                                          String compareType, double start, double stop,
                                          boolean append) {
    PrintWriter[] writers = new PrintWriter[comparedCalls.length];
    for (int i = 0; i < comparedCalls.length; i++) {
      String[] comparisons = getComparisons(comparedCalls[i]);
      String file = getFileName(rootDir, compareType, start, stop, comparisons);
      try {
        if (comparedCalls[i].getnumPassingAndPresent() != 0) {
          writers[i] = new PrintWriter(new FileWriter(file, append));
        }
      } catch (IOException e) {
        System.err.println("Error creating " + file);
        e.printStackTrace();
      }
    }
    return writers;
  }

  private static PrintWriter[] getWriters(String rootDir, CompareCalls_dev[] comparedCalls,
                                          String compareType, int start, int stop, boolean append) {
    PrintWriter[] writers = new PrintWriter[comparedCalls.length];
    for (int i = 0; i < comparedCalls.length; i++) {
      String[] comparisons = getComparisons(comparedCalls[i]);
      String file = getFileName(rootDir, compareType, start, stop, comparisons);
      try {
        if (comparedCalls[i].getnumPassingAndPresent() != 0) {
          writers[i] = new PrintWriter(new FileWriter(file, append));
        }
      } catch (IOException e) {
        System.err.println("Error creating " + file);
        e.printStackTrace();
      }
    }
    return writers;
  }

  // computes average percent match for each copy number, and a total for all copy numbers
  private static double[] goodCallCN(CNVariant[] cnvs0, CNVariant[] cnvs1, int[][] counts,
                                     double[] goodCalls) {
    // CN state percents as defined here are ,0,1,2,3,4,5...where 2 is a total non-copy number aware
    // count, and 5 is for non matching (overlap,exact,sigoverlap)
    int[] indCNNumbers = new int[6];
    indCNNumbers = getCNNumbers(cnvs0, cnvs1, indCNNumbers);
    for (int i = 0; i < counts.length; i++) {
      // System.out.println(i);
      if (indCNNumbers[i] > 0) {
        // System.out.println(cnvs0[0].getFamilyID() + "\t" + cnvs1[0].getFamilyID() + "\t" +
        // counts[i][3] + "\t" + counts[i][4] + "\t" + indCNNumbers[i] + "\t" + i);
        goodCalls[i] += 2 * counts[i][3] + counts[i][4];
        // System.out.println(cnvs0[0].getFamilyID() + "\t" + cnvs1[0].getFamilyID() + "\t" +
        // counts[i][3] + "\t" + counts[i][4] + "\t" + indCNNumbers[i] + "\t" + i + "\t" + ((double)
        // ((2 * counts[i][3]) + counts[i][4]) / (indCNNumbers[i])));

      }
    }
    return goodCalls;

  }

  // iterate(0.2, 0.6, 20, 0.0, 100, 10, 0, 20, rootDir, cnvFile, LRR_lookup,
  // "rootsANDdubs.compsmall.txt", markerFile);
  // iterate over parameters(lrr, conf, numMarkers) and compare results among duplicates
  public static void iterate(double lrrstart, double lrrstop, int numLrrBins, double confstart,
                             double confstop, int numConfBins, int probestart, int probeStop,
                             String rootDir, String cnvFile, String LRR_lookup, String compFile,
                             String markerFile) {
    PrintWriter[] probeWriters, lrrWriters, confWriters;
    CompareCalls_dev[] comparedCalls;
    double[] lrrValues = new double[numLrrBins];
    lrrValues = binIt(lrrstart, lrrstop, numLrrBins);
    double[] confValues = new double[numConfBins];
    confValues = binIt(confstart, confstop, numConfBins);
    long time;
    time = new Date().getTime();
    try {

      // for (int i = probestart; i < probeStop + 1; i += 5) {
      // for (int k = 0; k < lrrValues.length; k++) {
      // // for (int l = 0; l < confValues.length; l++) {
      // String iterType = "numMarkers";
      // comparedCalls = compare(rootDir, filterAndClean(rootDir, cnvFile, LRR_lookup, compFile,
      // lrrValues[k], DEFAULT_CONF_FILTER, i, markerFile, DEFAULT_BP_FRACTION,
      // DEFAULT_NUM_MARKER_FRACTION), compFile);
      // if (i == probestart && k == 0) {
      // probeWriters = getWriters(rootDir, comparedCalls, iterType, probestart, probeStop, false);
      // } else {
      // probeWriters = getWriters(rootDir, comparedCalls, iterType, probestart, probeStop, true);
      // }
      // for (int j = 0; j < comparedCalls.length; j++) {
      // if (comparedCalls[j].getnumPassingAndPresent() == 0) {
      // continue;
      // } else {
      // if (i == probestart && k == 0) {
      // probeWriters[j].println(getHeader("lrr\t" + iterType));
      // }
      // probeWriters[j].println(lrrValues[k] + "\t" + i + "\t" + getPrintString(comparedCalls[j]));
      // probeWriters[j].close();
      // }
      // // }
      // }
      // }
      // }

      for (int i = probestart; i < probeStop + 1; i += 5) {
        String iterType = "numMarkers";
        comparedCalls =
            compare(rootDir,
                    filterAndClean(rootDir, cnvFile, LRR_lookup, compFile, DEFAULT_LRR_FILTER,
                                   DEFAULT_CONF_FILTER, i, markerFile, DEFAULT_BP_FRACTION,
                                   DEFAULT_NUM_MARKER_FRACTION),
                    DEFAULT_COMP_LIST);
        if (i == probestart) {
          probeWriters = getWriters(rootDir, comparedCalls, iterType, probestart, probeStop, false);
        } else {
          probeWriters = getWriters(rootDir, comparedCalls, iterType, probestart, probeStop, true);
        }
        for (int j = 0; j < comparedCalls.length; j++) {
          if (comparedCalls[j].getnumPassingAndPresent() == 0) {
            continue;
          } else {
            if (i == probestart) {
              probeWriters[j].println(getHeader(iterType));
            }
            probeWriters[j].println(i + "\t" + getPrintString(comparedCalls[j]));
            probeWriters[j].close();
          }
        }
      }

      for (int i = 0; i < lrrValues.length; i++) {
        String iterType = "lrr";
        comparedCalls =
            compare(rootDir,
                    filterAndClean(rootDir, cnvFile, LRR_lookup, compFile, lrrValues[i],
                                   DEFAULT_CONF_FILTER, DEFAULT_NUM_MARKERS_FILTER, markerFile,
                                   DEFAULT_BP_FRACTION, DEFAULT_NUM_MARKER_FRACTION),
                    DEFAULT_COMP_LIST);
        if (i == 0) {
          lrrWriters = getWriters(rootDir, comparedCalls, iterType, lrrstart, lrrstop, false);
        } else {
          lrrWriters = getWriters(rootDir, comparedCalls, iterType, lrrstart, lrrstop, true);
        }
        for (int j = 0; j < comparedCalls.length; j++) {
          if (comparedCalls[j].getnumPassingAndPresent() == 0) {
            continue;
          } else {
            if (i == 0) {
              lrrWriters[j].println(getHeader(iterType));
            }
            lrrWriters[j].println(lrrValues[i] + "\t" + getPrintString(comparedCalls[j]));
            lrrWriters[j].close();
          }
        }
      }

      for (int i = 0; i < confValues.length; i++) {
        String iterType = "conf";
        comparedCalls =
            compare(rootDir,
                    filterAndClean(rootDir, cnvFile, LRR_lookup, compFile, DEFAULT_LRR_FILTER,
                                   confValues[i], DEFAULT_NUM_MARKERS_FILTER, markerFile,
                                   DEFAULT_BP_FRACTION, DEFAULT_NUM_MARKER_FRACTION),
                    DEFAULT_COMP_LIST);
        if (i == 0) {
          confWriters = getWriters(rootDir, comparedCalls, iterType, confstart, confstop, false);
        } else {
          confWriters = getWriters(rootDir, comparedCalls, iterType, confstart, confstop, true);
        }
        for (int j = 0; j < comparedCalls.length; j++) {
          System.out.println(j + "\t" + comparedCalls[j].getAveragePercent()[2]);
          if (comparedCalls[j].getnumPassingAndPresent() == 0) {
            continue;
          } else {
            if (i == 0) {
              confWriters[j].println(getHeader(iterType));
            }
            confWriters[j].println(confValues[i] + "\t" + getPrintString(comparedCalls[j]));
            confWriters[j].close();
          }
        }
      }
      System.out.println("Iterations took " + ext.getTimeElapsed(time));
    } catch (Exception e) {
      System.err.println("Error comparing ");
      e.printStackTrace();
    }

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String rootDir = DEFAULT_ROOT;
    String cnvFile = DEFAULT_CNV_FILE;
    // String[] compFiles = DEFAULT_COMP_LISTS;
    String LRR_lookup = DEFAULT_LRR_LOOKUP;
    double lrrFilter = DEFAULT_LRR_FILTER;
    double confFilter = DEFAULT_CONF_FILTER;
    int numMarkers = DEFAULT_NUM_MARKERS_FILTER;
    String markerFile = DEFAULT_MARKER_FILE;
    // String usage = "\n"+
    // "cnv.park.CompareCalls requires 0-1 arguments\n"+"" +
    // " (1) root directory where files are located\n"+
    // " (2) cnv filename (i.e. cnv=all_gw6.cnv)\n"+
    // " (3) marker filename with chromosome and base pair locations (i.e markers=hg18markers.txtn"+
    // " (4) lrr lookup filename with sample ids and LRR_SD (i.e. lrrLookup=lrrLookup.txt)\n"+
    // " (5) comparison file defining duplicates (comp=comp.txt)\n"+
    // " ADD the following if you want to dump the data to a text file\n"+
    // " (6) dump all regions for a particular transformation (i.e. dumpAll=(default))\n"+
    // " (7) the index of the transformation to export or -1 for all (i.e. transIndex=
    // (default))\n"+
    // " OR:\n"+
    // " (6) dump all transformations for a particular region (i.e. dump=chr8:25129632-25130278 (not
    // the default))\n"+
    // "";
    // options, clean file, filter file, iterate
    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        // System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        rootDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("files=")) {
        // files = args[i].split("=")[1].split(",");
        numArgs--;
      }
    }
    if (numArgs != 0) {
      // System.err.println(usage);
      System.exit(1);
    }
    try {
      // public static void iterate(double lrrstart, double lrrstop, int numLrrBins, double
      // confstart, double confstop, int numConfBins, int probestart, int probeStop, String rootDir,
      // String cnvFile, String LRR_lookup, String compFile, String markerFile) {

      // getMarkerLookup(rootDir, markerFile);
      // public static String[] filterAndClean(String rootDir, String cnvFile, String LRR_lookup,
      // String compFile, double lrrFilter, double confFilter, int numMarkers, String markerFile,
      // double bpFraction, double markerFraction) {
      filterAndClean(rootDir, cnvFile, LRR_lookup, "Full_SAMP.txt", lrrFilter, confFilter,
                     numMarkers, markerFile, 0.3, 0.3);
      // // cleanCNVs(rootDir, filterCNVs(rootDir, cnvFile, LRR_lookup, DEFAULT_COMP_LIST,
      // lrrFilter, confFilter, numMarkers)[0], markerFile, 0.2, 0.2);
      // compare(rootDir, filterAndClean(rootDir, cnvFile, LRR_lookup, DEFAULT_COMP_LIST, lrrFilter,
      // confFilter, numMarkers, markerFile, 0.2, 0.2), DEFAULT_COMP_LIST);
      // defineCompLists(rootDir, DEFAULT_COMP_LIST);
      // getLRRs(rootDir, LRR_lookup);
      // iterate(0.2, 0.6, 30, 0, 100, 10, 0, 100, rootDir, cnvFile, LRR_lookup, DEFAULT_COMP_LIST,
      // markerFile);
      // filterCNVs(rootDir, cnvFile, LRR_lookup, DEFAULT_COMP_LIST, lrrFilter, confFilter,
      // numMarkers);
      // compare(rootDir, filterCNVs(rootDir, cnvFile, LRR_lookup, compFiles, lrrFilter, confFilter,
      // numMarkers));
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  private static boolean markerGapCheck(CNVariant thisCNV, CNVariant nextCNV,
                                        Hashtable<Byte, List<Integer>> markers,
                                        double markerFraction) {
    int gapLengthMarkers = getStartIndex(nextCNV, markers) - getStopIndex(thisCNV, markers) - 1;
    int proposedNewMarkerLength =
        getStopIndex(nextCNV, markers) - getStartIndex(thisCNV, markers) + 1;
    return testFrac(markerFraction, ((double) gapLengthMarkers / proposedNewMarkerLength));
  }

  // new confidence value = old+old+..., as per PennCNV folks
  private static CNVariant mergeCNVs(CNVariant[] cnvs, int startCombine, int stopCombine,
                                     Hashtable<Byte, List<Integer>> markers) {
    CNVariant mergedCNV;
    double newConf = 0.0;
    int newNumMarkers = 0;
    for (int j = startCombine; j <= stopCombine; j++) {
      newConf += cnvs[j].getScore();
      if (j + 1 <= stopCombine) {
        newNumMarkers += getStartIndex(cnvs[j + 1], markers) - getStopIndex(cnvs[j], markers) - 1;
      }
      newNumMarkers += cnvs[j].getNumMarkers();
    }
    mergedCNV =
        new CNVariant(cnvs[startCombine].getFamilyID(), cnvs[startCombine].getIndividualID(),
                      cnvs[startCombine].getChr(), cnvs[startCombine].getStart(),
                      cnvs[stopCombine].getStop(), cnvs[startCombine].getCN(), newConf,
                      newNumMarkers, cnvs[startCombine].getSource());
    return mergedCNV;
  }

  // so that an index represents an absolute position within chromosome markers
  private static void sortBPbyChromosome(Hashtable<Byte, List<Integer>> markers) {
    for (int i = 1; i <= markers.keySet().size(); i++) {
      if (markers.containsKey((byte) i)) {
        Collections.sort(markers.get((byte) i));
      } else {
        System.err.println("chromosome " + i + " is undefined");
      }
    }
  }

  private static boolean testFrac(double Fraction, double FractionToTest) {
    if (FractionToTest <= Fraction) {
      return true;
    } else {
      return false;
    }
  }

  public static CNVariant[] toCNVArray(ArrayList<CNVariant> cnvVariants) {
    return cnvVariants.toArray(new CNVariant[cnvVariants.size()]);
  }

  public static String[] toStringArray(ArrayList<String> stringList) {
    return stringList.toArray(new String[stringList.size()]);
  }

  private final int numPassingAndPresent;

  private final double[] averagePercent;

  private final int[] totalCallsAnalyzed;

  private final String[] cnvFiles;

  private final int[] exactMatches;

  private final int[] olapMatches;

  private final int[] sigOlapMatches;

  private final double[] indPercentMatches;

  private final int[][] indCNVCalls;

  private final String[][] comparedIDs;

  public CompareCalls_dev(int numPassingAndPresent, double averagePercent[],
                          int[] totalCallsAnalyzed, String[] cnvFiles, int[] exactMatches,
                          int[] sigOlapMatches, int[] olapMatches, String[][] comparedIDs,
                          int[][] indCNVCalls, double[] indPercentMatches) {
    this.numPassingAndPresent = numPassingAndPresent;
    this.averagePercent = averagePercent;
    this.totalCallsAnalyzed = totalCallsAnalyzed;
    this.cnvFiles = cnvFiles;
    this.exactMatches = exactMatches;
    this.sigOlapMatches = sigOlapMatches;
    this.olapMatches = olapMatches;
    this.comparedIDs = comparedIDs;
    this.indCNVCalls = indCNVCalls;
    this.indPercentMatches = indPercentMatches;
  }

  public double[] getAveragePercent() {
    return averagePercent;
  }

  public String[] getCnvFiles() {
    return cnvFiles;
  }

  public String[][] getcomparedIDs() {
    return comparedIDs;
  }

  public int[] getExactMatches() {
    return exactMatches;
  }

  public int[][] getIndCNVCalls() {
    return indCNVCalls;
  }

  public double[] getindPercentMatches() {
    return indPercentMatches;
  }

  public int getnumPassingAndPresent() {
    return numPassingAndPresent;
  }

  public int[] getOlapMatches() {
    return olapMatches;
  }

  public int[] getSigOlapMatches() {
    return sigOlapMatches;
  }

  public int[] getTotalCallsAnalyzed() {
    return totalCallsAnalyzed;
  }
}

// private static Hashtable<Integer, Hashtable<String, ArrayList<CNVariant>>> getCNVHash(String[]
// files) {
// BufferedReader reader;
// String[] line;
// ArrayList<CNVariant> v;
// Hashtable<Integer, Hashtable<String, ArrayList<CNVariant>>> hash = new Hashtable<Integer,
// Hashtable<String, ArrayList<CNVariant>>>();
// Hashtable<String, ArrayList<CNVariant>> source = new Hashtable<String, ArrayList<CNVariant>>();
// for (int i = 0; i < files.length; i++) {
// try {
// reader = new BufferedReader(new FileReader(files[i]));
// if (!ext.checkHeader(reader.readLine().trim().split("[\\s]+"), CNVariant.PLINK_CNV_HEADER,
// false)) {
// reader.close();
// System.err.println("quitting comparison");
// }
// while (reader.ready()) {
// line = reader.readLine().trim().split("[\\s]+");
// if (hash.containsKey(i)) {
// source = hash.get(i);
// } else {
// hash.put(i, source = new Hashtable<String, ArrayList<CNVariant>>());
// }
// if (source.containsKey(line[0])) {
// v = source.get(line[0]);
// } else {
// source.put(line[0], v = new ArrayList<CNVariant>());
// }
// v.add(new CNVariant(line, i));
// }
// reader.close();
// } catch (FileNotFoundException fnfe) {
// System.err.println("Error: file \"" + files[i] + "\" not found in current directory");
// System.exit(1);
// } catch (IOException ioe) {
// System.err.println("Error reading file \"" + files[i] + "\"");
// System.exit(2);
// }
// }
// return hash;
// }

// if (testFrac(markerFraction, ((double) gapLengthMarkers / proposedNewMarkerLength))) {
// System.out.println(nextCNV.getIndividualID() + "\t" + thisCNV.getIndividualID() + "\t" +
// markers.get(nextCNV.getChr()).indexOf(nextCNV.getStart()) + "\t" +
// (markers.get(thisCNV.getChr()).indexOf(thisCNV.getStop()) + 1));
// System.out.println(gapLengthMarkers + "\t" + proposedNewMarkerLength + "\t" +
// markerFractionToTest + "\t" + testFrac(markerFraction, markerFractionToTest));
//
// }

// System.out.println(nextCNV.getIndividualID() + "\t" + thisCNV.getIndividualID() + "\t" +
// markers.get(nextCNV.getChr()).indexOf(nextCNV.getStart()) + "\t" +
// (markers.get(thisCNV.getChr()).indexOf(thisCNV.getStop()) + 1));
// System.out.println(nextCNV.getChr() + "\t" + nextCNV.getStart() + "\t" + thisCNV.getChr() + "\t"
// + thisCNV.getStart());

// for (int i = 0; i < lrrValues.length; i++) {
// comparedCalls = compare(rootDir, filterCNVs(rootDir, cnvFile, LRR_lookup, compFiles,
// lrrValues[i], DEFAULT_CONF_FILTER, DEFAULT_NUM_MARKERS_FILTER));
// for (int j = 0; j < comparedCalls.getInds().length; j++) {
// if (comparedCalls.getIndCNVCalls()[0][j] > 0 && comparedCalls.getIndCNVCalls()[1][j] > 0) {
// lrrWriter.println(comparedCalls.getInds()[j] + "\t" + lrrValues[i] + "\t" +
// comparedCalls.getIndCNVCalls()[0][j] + "\t" + comparedCalls.getIndCNVCalls()[1][j] + "\t" +
// comparedCalls.getExactMatches()[j] + "\t" + comparedCalls.getSigOlapMatches()[j] + "\t" +
// comparedCalls.getindPercentMatches()[j] + "\t" + comparedCalls.getAveragePercent() + "\t" +
// comparedCalls.getnumPassingAndPresent() + "\t" + comparedCalls.getTotalCallsAnalyzed() + "\t" +
// Array.toStr(comparedCalls.getCnvFiles()) + "\t");
// }
// }
// // lrrWriter.println(lrrValues[i] + "\t" + comparedCalls.getAveragePercent() + "\t" +
// comparedCalls.getnumPassingAndPresent() + "\t" + comparedCalls.getTotalCallsAnalyzed() + "\t" +
// Array.toStr(comparedCalls.getCnvFiles()));
// }
// for (int i = 0; i < confValues.length; i++) {
// comparedCalls = compare(rootDir, filterCNVs(rootDir, cnvFile, LRR_lookup, compFiles,
// DEFAULT_LRR_FILTER, confValues[i], DEFAULT_NUM_MARKERS_FILTER));
// for (int j = 0; j < comparedCalls.getInds().length; j++) {
// if (comparedCalls.getIndCNVCalls()[0][j] > 0 && comparedCalls.getIndCNVCalls()[1][j] > 0) {
// confWriter.println(comparedCalls.getInds()[j] + "\t" + confValues[i] + "\t" +
// comparedCalls.getIndCNVCalls()[0][j] + "\t" + comparedCalls.getIndCNVCalls()[1][j] + "\t" +
// comparedCalls.getExactMatches()[j] + "\t" + comparedCalls.getSigOlapMatches()[j] + "\t" +
// comparedCalls.getindPercentMatches()[j] + "\t" + comparedCalls.getAveragePercent() + "\t" +
// comparedCalls.getnumPassingAndPresent() + "\t" + comparedCalls.getTotalCallsAnalyzed() + "\t" +
// Array.toStr(comparedCalls.getCnvFiles()) + "\t");
// }
// }
//
// // confWriter.println(confValues[i] + "\t" + comparedCalls.getAveragePercent() + "\t" +
// comparedCalls.getnumPassingAndPresent() + "\t" + comparedCalls.getTotalCallsAnalyzed() + "\t" +
// Array.toStr(comparedCalls.getCnvFiles()));
// }
// // for (int j = 0; j < lrrValues.length; j++) {
// // for (int i = probestart; i < probeStop + 1; i++) {
// // comparedCalls = compare(rootDir, filterCNVs(rootDir, cnvFile, LRR_lookup, compFiles,
// lrrValues[j], DEFAULT_CONF_FILTER, i));
// // lrrProbeWriter.println(lrrValues[j] + "\t" + i + "\t" + comparedCalls.getAveragePercent() +
// "\t" + comparedCalls.getnumPassingAndPresent() + "\t" + comparedCalls.getTotalCallsAnalyzed() +
// "\t" + Array.toStr(comparedCalls.getCnvFiles()));
// // }
// // }

// confWriter = new PrintWriter(new FileWriter(rootDir + "confConcordance_" + confstart + "_" +
// confstop + "LRR" + DEFAULT_LRR_FILTER + "numMarkers" + DEFAULT_NUM_MARKERS_FILTER + ".concord"));
// lrrWriter = new PrintWriter(new FileWriter(rootDir + "lrrConcordance_" + lrrstart + "_" + lrrstop
// + "numMarkers" + DEFAULT_NUM_MARKERS_FILTER + "Conf" + DEFAULT_CONF_FILTER + ".concord"));

// if (isDefined(inds, allPossibleCombinations, defineCompHash, hash, allPossibleCombinations[i][0],
// allPossibleCombinations[i][1], j, k)) {
// // System.out.println("hhragdsafdsacomparing " + inds[j] + " and " + inds[k]);
// cnvs = new CNVariant[][] { CNVariant.toCNVariantArray(),
// CNVariant.toCNVariantArray(hash.get(inds[j]).get(allPossibleCombinations[i][1] + "")) };
// if (hash.get(inds[j]).containsKey(allPossibleCombinations[i][0] + "")) {
// cnvs = new CNVariant[][] {
// CNVariant.toCNVariantArray(hash.get(inds[j]).get(allPossibleCombinations[i][0] + "")),
// CNVariant.toCNVariantArray(hash.get(inds[k]).get(allPossibleCombinations[i][1] + "")) };
// }
//
// System.out.println(cnvs[0].length);
//
// counts = new int[5];
// // // we handle samples missing in one of the files later
// checkCNVs_a_b(cnvs, counts);
// checkCNVs_b_a(cnvs, counts);
//
// exactMatches[j] = counts[3];
//
// sigOlapMatches[j] = counts[4];
// comparedIDs[0][j] = inds[j].replaceAll("\t.*", "");
// indCNVCalls[0][j] = cnvs[0].length;
// indCNVCalls[1][j] = cnvs[1].length;
// totalCallsAnalyzed += cnvs[0].length + cnvs[1].length;
//
// if (cnvs[0].length > 0 && cnvs[1].length > 0) {
// System.out.println(numPassingAndPresent);
// numPassingAndPresent++;
// double callPercent = ((double) ((2 * counts[3]) + counts[4]) / (cnvs[0].length +
// cnvs[1].length));
// goodCallPPercent += callPercent;
// indPercentMatches[j] = callPercent;
// } else {
// indPercentMatches[j] = Double.NaN;
// }
// // public CompareCalls_dev(int numPassingAndPresent, double averagePercent, int
// totalCallsAnalyzed, String[] cnvFiles, int[] exactMatches, int[] sigOlapMatches, String[] inds,
// int[][] indCNVCalls, double[] indPercentMatches) {
//
// }

// private static Hashtable<String, Hashtable<String, ArrayList<CNVariant>>> getCNVs(String file) {
// ArrayList<CNVariant> v = new ArrayList<CNVariant>();
// BufferedReader reader;
// String[] line;
// Hashtable<String, Hashtable<String, ArrayList<CNVariant>>> hash = new Hashtable<String,
// Hashtable<String, ArrayList<CNVariant>>>();
// Hashtable<String, ArrayList<CNVariant>> source = new Hashtable<String, ArrayList<CNVariant>>();
// try {
// reader = new BufferedReader(new FileReader(file));
// if (!ext.checkHeader(reader.readLine().trim().split("[\\s]+"), CNVariant.PLINK_CNV_HEADER,
// false)) {
// reader.close();
// System.err.println("quitting comparison");
// }
// while (reader.ready()) {
// line = reader.readLine().trim().split("[\\s]+");
// if (hash.containsKey(line[0] + "\t" + line[1])) {
// source = hash.get(line[0] + "\t" + line[1]);
// } else {
// hash.put(line[0] + "\t" + line[1], source = new Hashtable<String, ArrayList<CNVariant>>());
// }
// if (source.containsKey(0 + "")) {
// v = source.get(0 + "");
// } else {
// source.put(0 + "", v = new ArrayList<CNVariant>());
// }
// v.add(new CNVariant(line, 0));
// }
// reader.close();
// } catch (FileNotFoundException fnfe) {
// System.err.println("Error: file \"" + file + "\" not found in current directory");
// System.exit(1);
// } catch (IOException ioe) {
// System.err.println("Error reading file \"" + file + "\"");
// System.exit(2);
// }
// return hash;
// }

// private static boolean isDefined(String[] inds, int[][] allPossibleCombinations,
// Hashtable<String, Hashtable<String, Integer>> defineCompHash, Hashtable<Integer,
// Hashtable<String, ArrayList<CNVariant>>> hash, int combo1, int combo2, int j, int k) {
// if (defineCompHash.containsKey(inds[j]) && defineCompHash.get(inds[j]).containsKey(inds[k])) {
// if (hash.containsKey(combo1) && hash.containsKey(combo1)) {
// if (hash.get(combo1).containsKey(inds[j])) {
// if (hash.get(combo2).containsKey(inds[k])) {
// return true;
// } else {
// return false;
// }
// } else {
// return false;
// }
// } else {
// return false;
// }
// } else {
// return false;
// }
//
// // return ((defineCompHash.containsKey(inds[j]) &&
// defineCompHash.get(inds[j]).containsKey(inds[k]))||defineCompHash.containsKey(inds[j]) &&
// defineCompHash.get(inds[j]).containsKey(inds[k])) && hash.containsKey(inds[j])&&
// hash.containsKey(inds[k]) && (hash.get(inds[j]).containsKey(allPossibleCombinations[i][0] + "") )
// && hash.get(inds[k]).containsKey(allPossibleCombinations[i][1] +
// "")||hash.get(inds[j]).containsKey(allPossibleCombinations[i][0] + "") ) &&
// hash.get(inds[k]).containsKey(allPossibleCombinations[i][1] + "");
// }
// for (int i = probestart; i < probeStop + 1; i++) {
// comparedCalls = compare(rootDir, filterAndClean(rootDir, cnvFile, LRR_lookup, compFile,
// DEFAULT_LRR_FILTER, DEFAULT_CONF_FILTER, i, markerFile, DEFAULT_BP_FRACTION,
// DEFAULT_NUM_MARKER_FRACTION), DEFAULT_COMP_LIST);
// if (i == probestart) {
// probeWriters = getWriters(rootDir, comparedCalls, "probe", probestart, probeStop, false);
// } else {
// probeWriters = getWriters(rootDir, comparedCalls, "probe", probestart, probeStop, true);
// }
// for (int j = 0; j < comparedCalls.length; j++) {
// if (comparedCalls[j].getTotalCallsAnalyzed() == 0) {
// continue;
// } else {
// probeWriters[j].println(i + "\t" + comparedCalls[j].getAveragePercent());
// probeWriters[j].close();
// System.out.println(i + "\t" + comparedCalls[j].getAveragePercent());
// }
// }
// }
// for (int j = 0; j < comparedCalls.length; j++) {
// System.out.println(j);
// if (comparedCalls[j] != null) {
// probeWriter.println(comparedCalls[j].getcomparedIDs()[0][0]);
// // + "\t" + 20 + "\t" + comparedCalls[j].getIndCNVCalls()[0][j] + "\t" +
// comparedCalls[j].getIndCNVCalls()[1][j] + "\t" + comparedCalls[j].getExactMatches()[j] + "\t" +
// comparedCalls[j].getSigOlapMatches()[j] + "\t" + comparedCalls[j].getindPercentMatches()[j] +
// "\t" + comparedCalls[j].getAveragePercent() + "\t" + comparedCalls[j].getnumPassingAndPresent() +
// "\t" + comparedCalls[j].getTotalCallsAnalyzed() + "\t" +
// Array.toStr(comparedCalls[j].getCnvFiles()) + "\t");
// }
// // && comparedCalls[j].getIndCNVCalls()[0][j] > 0 && comparedCalls[j].getIndCNVCalls()[1][j] > 0
// }

// probeWriter.close();
// }
// lrrProbeWriter.close();
