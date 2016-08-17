package org.genvisis.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Map;
import java.util.Vector;

import org.genvisis.bioinformatics.Alleles;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Numbers;
import org.genvisis.common.ext;
import org.genvisis.stats.ProbDist;

public class ResultsPackager {
  public static final String[] IBC_OUTPUT_FORMAT1 =
      {"CHR", "POS", "SNP", "STRAND (Illumina)", "STRAND (HapMap)", "N", "EFFECT_ALLELE1",
       "NON_EFFECT_ALLELE", "EA_FREQ", "BETA", "SE", "P_VAL"};
  public static final String[] TRADITIONAL_OUTPUT_FORMAT =
      {"Chr", "Position", "MarkerName", "Strand", "HapMapStrand", "N", "Effect_allele",
       "Reference_allele", "Freq1", "BETA", "SE", "P-value"};
  public static final String[] STANDARD_OUTPUT_FORMAT =
      {"MarkerName", "Chr", "Position", "Effect_allele", "Reference_allele",
       "Effect_allele_frequency", "N", "BETA", "SE", "P-value"};
  // public static final String[] ABSOLUTE_MINIMUM_OUTPUT_FORMAT = {"MarkerName", "Effect_allele",
  // "BETA", "SE", "P-value"}; // ChiSquare value
  public static final String[] EMIM_OUTPUT_FORMAT_SNPS =
      {"Chr", "Pos", "MarkerName", "allele_A", "allele_B"};
  public static final String[] EMIM_OUTPUT_FORMAT_MENDEL_ERRORS = {"Mendel_Errors"};
  public static final String[] EMIM_OUTPUT_FORMAT_HWE = {"HWE_GENO", "HWE_P", "SigHWE"};
  public static final String[] EMIM_OUTPUT_FORMAT_TDT =
      {"tdt_T", "tdt_U", "tdt_OR", "tdt_P", "tdt_L95", "tdt_U95", "tdt_BETA", "tdt_SE"};
  public static final String[] EMIM_OUTPUT_FORMAT_EMIM_RESULTS =
      {"freq", "C_lnR1", "C_se_lnR1", "C_lnR2", "C_se_lnR2", "C_lnS1", "C_se_lnS1", "C_lnS2",
       "C_se_lnS2", "CM_lnR1", "CM_se_lnR1", "CM_lnR2", "CM_se_lnR2", "CM_lnS1", "CM_se_lnS1",
       "CM_lnS2", "CM_se_lnS2", "M_lnR1", "M_se_lnR1", "M_lnR2", "M_se_lnR2", "M_lnS1", "M_se_lnS1",
       "M_lnS2", "M_se_lnS2"};
  public static final String[] EMIM_OUTPUT_FORMAT_EMIM_PVALS =
      {"pVal_C_df2", "pVal_C_df2_Excel", "pVal_C_df1", "pVal_C_df1_Excel", "pVal_M_df2",
       "pVal_M_df2_Excel", "pVal_M_df1", "pVal_M_df1_Excel", "pVal_CM-C_df2", "pVal_CM-C_df2_Excel",
       "pVal_CM-C_df1", "pVal_CM-C_df1_Excel", "pVal_CM-M_df2", "pVal_CM-M_df2_Excel",
       "pVal_CM-M_df1", "pVal_CM-M_df1_Excel"};

  public static final String[] PLINK_REQS = {"SNP", "A1", "TEST", "NMISS", "OR", "BETA", "SE", "P"};
  public static final String[] SOL_REQS =
      {"Variant_ID", "Beta", "Se", "Pvalue", "CAF", "CAC", "N0", "N1", "N2", "NMISS"};
  // Standard Errors in emimsummary.out are currently labeled "sd"
  public static final String[] EMIM_REQS =
      {"snpID", "freq", "lnR1", "sd_lnR1", "lnR2", "sd_lnR2", "lnS1", "sd_lnS1", "lnS2", "sd_lnS2",
       "lnliknull", "lnlikfull"};
  public static final String[] TDT_REQS = {"T", "U", "OR", "P", "L95", "U95"};

  public static String[] addPrefixToArray(String prefix, String[] array, Logger log) {
    String[] result;

    result = new String[array.length];
    for (int i = 0; i < array.length; i++) {
      result[i] = prefix + array[i];
    }

    return result;
  }

  public static void createFromParameters(String filename, Logger log) {
    Vector<String> params;

    params = Files.parseControlFile(filename, "results",
                                    new String[] {"dir=",
                                                  "results=plink.assoc.linear or plink.assoc.logistic",
                                                  "type=plink, gwaf, sol, emim, or  forest",
                                                  "map=plink.map", "freq=plink.frq", "customFreq=",
                                                  "list=specificSNPs.txt # leave blank for all",
                                                  "filter=1.0 # set to 0.001 to report only those variants with a p<=0.001",
                                                  "out=finalProduct.out"},
                                    log);

    if (params != null) {
      params.add("log=" + log.getFilename());
      main(Array.toStringArray(params));
    }
  }

  /**
   * 
   * @param fullPathStatResults {{name of analysis ("." for main/overall), path of results
   *        file},...}
   * @param fullPathMarkerList Path to list of markers to use
   * @param markerColumnName name to use for marker column
   * @param analyses names of analyses
   * @param columnNamesOfAnalyses String array for each analysis (i) where [i][0] = beta (OR for
   *        tdt), [i][1] = SE (U95 for tdt), [i][2] = primary result, [i][3+] = secondary results
   * @param columnDisplayNames String array for each analysis corresponding to columns in
   *        columnNamesOfAnalyses with desired display text (ignored for [i][0] and [i][1], null to
   *        use column names
   * @param fullPathOutFile Path to generate output file to
   * @param log a log
   * 
   *        Examples of using this method:
   * 
   *        getForestPlotParameterFile(new String[][] {{"extragonadal",
   *        "D:/temp/Poynter_emim/testing/allFinalPoynter_results_pvals_1.xln"}, {"germinoma",
   *        "D:/temp/Poynter_emim/testing/allFinalPoynter_results_pvals_2.xln"}},
   *        "D:/temp/Poynter_emim/testing/markerList.txt", "MarkerName", new String[] {"tdt",
   *        "emim"}, new String[][] {{"tdt_OR", "tdt_U95", "tdt_P"}, {"C_lnR1", "C_se_lnR1",
   *        "pVal_C_df1"}}, null "D:/temp/Poynter_emim/testing/forestplot.xln", null);
   * 
   *        getForestPlotParameterFile(HashVec.loadFileToStringMatrix("/home/pankrat2/shared/Poynter_emim/allFinalPoynter/fileList_allFinalPoynter.txt",
   *        false, null, false), "/home/pankrat2/shared/Poynter_emim/markerList.txt", "MarkerName",
   *        new String[] {"TDT", "EMIM Maternal Effect"}, new String[][] {{"tdt_OR", "tdt_U95",
   *        "tdt_P"}, {"M_lnS1", "M_se_lnS1", "pVal_M_df1", "pVal_CM-C_df1"}}, new String[][]
   *        {{"","","p"},{"","","p","removing Child Effect p"}},
   *        "/home/pankrat2/shared/Poynter_emim/allFinalPoynter/allFinalPoynter_forestplot.xln",
   *        null);
   */
  public static void getForestPlotParameterFile(String[][] fullPathStatResults,
                                                String fullPathMarkerList, String markerColumnName,
                                                String[] analyses, String[][] columnNamesOfAnalyses,
                                                String[][] columnDisplayNames,
                                                String fullPathOutFile, Logger log) {
    Double beta, se;
    int markerIndex, mainIndex, columnIndex;
    String filename;
    String[][] markerList;
    String[] markers, displayMarkers, columnNamesToLoad, out1, out2;
    String[][][] statResults;
    Vector<String> columnsTmp;
    boolean isMain;

    if (columnDisplayNames == null) {
      columnDisplayNames = columnNamesOfAnalyses;
    }

    columnsTmp = new Vector<String>();
    for (String[] columnNamesOfAnalyse : columnNamesOfAnalyses) {
      for (int j = 0; j < columnNamesOfAnalyse.length; j++) {
        columnsTmp.add(columnNamesOfAnalyse[j]);
      }
    }
    columnNamesToLoad = columnsTmp.toArray(new String[0]);

    markerList = HashVec.loadFileToStringMatrix(fullPathMarkerList, false, null, false);
    markers = new String[markerList.length];
    displayMarkers = new String[markerList.length];
    for (int i = 0; i < markerList.length; i++) {
      markers[i] = markerList[i][0];
      displayMarkers[i] = markerList[i][1];
    }
    statResults = new String[fullPathStatResults.length][][];
    for (int i = 0; i < statResults.length; i++) {
      System.out.println("Loading " + fullPathStatResults[i][1]);
      statResults[i] =
          loadFile(fullPathStatResults[i][1], markerColumnName, markers, columnNamesToLoad, log);
    }

    out2 = new String[analyses.length * markers.length];
    mainIndex = -1;
    for (int i = 0; i < analyses.length; i++) {
      filename =
          ext.rootOf(fullPathOutFile) + "_" + (analyses[i].equalsIgnoreCase("tdt") ? "_" : "")
                 + ext.replaceWithLinuxSafeCharacters(analyses[i], true)
                 + fullPathOutFile.substring(fullPathOutFile.lastIndexOf("."));
      out1 = new String[markers.length + 1];
      out1[0] = "Gene\t" + markerColumnName;

      // Header of the output file
      for (int j = 0; j < fullPathStatResults.length; j++) {
        // if (analyses[i].equalsIgnoreCase("tdt")) {
        isMain = (fullPathStatResults[j][0].equals(".") || fullPathStatResults[j][0].equals(""));
        if (isMain && mainIndex == -1) {
          mainIndex = j;
        }
        out1[0] += "\tbeta" + (isMain ? "" : "." + fullPathStatResults[j][0]) + "\tse"
                   + (isMain ? "" : "." + fullPathStatResults[j][0]) + "\tOR"
                   + (isMain ? "" : "." + fullPathStatResults[j][0]) + "\tLowerCI"
                   + (isMain ? "" : "." + fullPathStatResults[j][0]) + "\tUpperCI"
                   + (isMain ? "" : "." + fullPathStatResults[j][0]);
        for (int k = 2; k < columnNamesOfAnalyses[i].length; k++) {
          out1[0] +=
              "\t" + columnNamesOfAnalyses[i][k] + (isMain ? "" : "." + fullPathStatResults[j][0]);
        }
        // } else {
        // for (int k = 0; k < columnNamesOfAnalyses[i].length; k++) {
        // line[0] += "\t" + fullPathStatResults[j][0] + "_" + columnNamesOfAnalyses[i][k];
        // }
        // }
      }

      // content of the output file
      for (int j = 1; j < out1.length; j++) {
        markerIndex = j - 1;
        int out2index = markerIndex * analyses.length + i;
        out1[j] = "Gene\t" + markers[markerIndex];
        for (int k = 0; k < statResults.length; k++) {
          if (k == mainIndex) {
            out2[out2index] = markers[markerIndex] + "\t" + filename + "\t" + analyses[i];
          }

          if (statResults[k][markerIndex][0] == null) {
            out1[j] += "\t.\t.\t.\t.\t.";
            for (int l = 2; l < columnNamesOfAnalyses[i].length; l++) {
              out1[j] += "\t.";
            }
          } else {
            if (analyses[i].equalsIgnoreCase("tdt")) {
              int betaInd = ext.indexOfStr(columnNamesOfAnalyses[i][0], columnNamesToLoad);
              int seInd = ext.indexOfStr(columnNamesOfAnalyses[i][1], columnNamesToLoad);
              beta = Math.log(Double.parseDouble(statResults[k][markerIndex][betaInd]));
              se = ((Math.log(Double.parseDouble(statResults[k][markerIndex][seInd])) - beta)
                    / 1.96);
              out1[j] += "\t" + beta + "\t" + se + "\t" + Math.exp(beta) + "\t"
                         + Math.exp(beta - 1.96 * se) + "\t" + Math.exp(beta + 1.96 * se);
            } else {
              int betaInd = ext.indexOfStr(columnNamesOfAnalyses[i][0], columnNamesToLoad);
              int seInd = ext.indexOfStr(columnNamesOfAnalyses[i][1], columnNamesToLoad);
              beta = Double.parseDouble(statResults[k][markerIndex][betaInd]);
              se = Double.parseDouble(statResults[k][markerIndex][seInd]);
              out1[j] += "\t" + statResults[k][markerIndex][betaInd] + "\t"
                         + statResults[k][markerIndex][seInd] + "\t" + Math.exp(beta) + "\t"
                         + Math.exp(beta - 1.96 * se) + "\t" + Math.exp(beta + 1.96 * se);
            }

            for (int l = 2; l < columnNamesOfAnalyses[i].length; l++) {
              columnIndex = 0;
              for (int m = 0; m < i; m++) {
                columnIndex += columnNamesOfAnalyses[m].length;
              }
              columnIndex += l;
              out1[j] += "\t" + statResults[k][markerIndex][columnIndex];

              if (k == mainIndex) {
                if (l == 2) {
                  out2[out2index] += " " + columnDisplayNames[i][l] + "="
                                     + ext.prettyP(Double.parseDouble(statResults[mainIndex][markerIndex][columnIndex]));
                } else {
                  out2[out2index] += " (" + columnDisplayNames[i][l] + "="
                                     + ext.prettyP(Double.parseDouble(statResults[mainIndex][markerIndex][columnIndex]))
                                     + ")";
                }

              }
            }
          }
          if (k == mainIndex) {
            out2[out2index] += "\t" + displayMarkers[markerIndex];
          }
        }
      }
      Files.writeList(out1, ext.parseDirectoryOfFile(fullPathOutFile) + filename);

      for (int j = 0; j < markers.length; j++) {
        if (statResults[mainIndex][j][0] == null) {
        } else {
          if (analyses[i].equalsIgnoreCase("tdt")) {
          } else {
          }
        }
        filename = ext.rootOf(filename) + filename.substring(filename.lastIndexOf("."));
      }
    }

    Files.writeList(out2, ext.parseDirectoryOfFile(fullPathOutFile) + ext.rootOf(fullPathOutFile)
                          + ".input");
  }

  private static String getOutputString(Hashtable<Long, String> snpList,
                                        Hashtable<String, String[]> mendelErrors,
                                        Hashtable<String, String[]> hwe, double hweThreshold,
                                        String freq, String[] lineC, Map<String, Integer> indicesC,
                                        String[] lineM, Map<String, Integer> indicesM,
                                        String[] lineCM, Map<String, Integer> indicesCM,
                                        Hashtable<String, String[]> tdtResults, double[] pvals,
                                        String[] pvalEquations, Logger log) {
    String result = null, tmp, snp;
    String[] tmp2;
    long index;

    index =
        Long.parseLong(lineC[indicesC.get("snpID")].substring(0,
                                                              lineC[indicesC.get("snpID")].indexOf(".")));
    if (snpList.containsKey(index)) {
      tmp = "";
      snp = snpList.get(index).split("\t")[2];
      if (mendelErrors != null) {
        tmp2 = mendelErrors.get(snp);
        for (String element : tmp2) {
          tmp += ("\t" + element);
        }
      }
      if (hwe != null) {
        tmp2 = hwe.get(snp);
        for (String element : tmp2) {
          tmp += ("\t" + element);
        }
        if ((!tmp2[1].equals("NA")) && Double.parseDouble(tmp2[1]) < hweThreshold) {
          tmp += ("\t1");
        } else {
          tmp += ("\t0");
        }
      }
      if (tdtResults != null) {
        tmp2 = tdtResults.get(snp);
        for (String element : tmp2) {
          tmp += ("\t" + element);
        }
        double beta, se;
        try {
          beta = Math.log(Double.parseDouble(tmp2[ext.indexOfStr("OR", TDT_REQS)]));
          try {
            se = ((Math.abs(Math.log(Double.parseDouble(tmp2[ext.indexOfStr("L95", TDT_REQS)]))
                            - beta)
                   + Math.abs(Math.log(Double.parseDouble(tmp2[ext.indexOfStr("U95", TDT_REQS)]))
                              - beta))
                  / 2)
                 / 1.96;
          } catch (NumberFormatException nfe) {
            se = Double.NaN;
          }
        } catch (NumberFormatException nfe) {
          beta = Double.NaN;
          se = Double.NaN;
        }
        tmp += "\t" + (!Numbers.isFinite(beta) ? "." : beta);
        tmp += "\t" + (!Numbers.isFinite(se) ? "." : se);
      }
      result = snpList.get(index) + tmp + "\t" + freq + "\t" + lineC[indicesC.get("lnR1")] + "\t"
               + lineC[indicesC.get("sd_lnR1")] + "\t" + lineC[indicesC.get("lnR2")] + "\t"
               + lineC[indicesC.get("sd_lnR2")] + "\t" + lineC[indicesC.get("lnS1")] + "\t"
               + lineC[indicesC.get("sd_lnS1")] + "\t" + lineC[indicesC.get("lnS2")] + "\t"
               + lineC[indicesC.get("sd_lnS2")] + "\t" + lineCM[indicesCM.get("lnR1")] + "\t"
               + lineCM[indicesCM.get("sd_lnR1")] + "\t" + lineCM[indicesCM.get("lnR2")] + "\t"
               + lineCM[indicesCM.get("sd_lnR2")] + "\t" + lineCM[indicesCM.get("lnS1")] + "\t"
               + lineCM[indicesCM.get("sd_lnS1")] + "\t" + lineCM[indicesCM.get("lnS2")] + "\t"
               + lineCM[indicesCM.get("sd_lnS2")] + "\t" + lineM[indicesM.get("lnR1")] + "\t"
               + lineM[indicesM.get("sd_lnR1")] + "\t" + lineM[indicesM.get("lnR2")] + "\t"
               + lineM[indicesM.get("sd_lnR2")] + "\t" + lineM[indicesM.get("lnS1")] + "\t"
               + lineM[indicesM.get("sd_lnS1")] + "\t" + lineM[indicesM.get("lnS2")] + "\t"
               + lineM[indicesM.get("sd_lnS2")];
      // result = snpList.get(index) + "\t" + lineC[indicesC[1]];
      for (int i = 0; i < pvals.length; i++) {
        result += "\t" + pvals[i] + "\t" + pvalEquations[i];
      }
    } else {
      log.reportError("Error - no map position for " + index);
    }

    return result;
  }

  private static String[] getPvalEquations(String logLikilihood_null_C, String logLikilihood_full_C,
                                           String logLikilihood_null_M, String logLikilihood_full_M,
                                           String logLikilihood_full_CM, Logger log) {

    // Lines up with odd entries in EMIM_OUTPUT_FORMAT_EMIM_PVALS
    return new String[] {"=1-CHISQ.DIST(2 * (" + logLikilihood_full_C + "-"
                         + logLikilihood_null_C + "),2,TRUE)",
                         "=1-CHISQ.DIST(2 * (" + logLikilihood_full_C + "-" + logLikilihood_null_C
                                                               + "),1,TRUE)",
                         "=1-CHISQ.DIST(2 * (" + logLikilihood_full_M + "-" + logLikilihood_null_M + "),2,TRUE)",
                         "=1-CHISQ.DIST(2 * (" + logLikilihood_full_M + "-" + logLikilihood_null_M + "),1,TRUE)",
                         "=1-CHISQ.DIST(2 * (" + logLikilihood_full_CM + "-" + logLikilihood_full_C + "),2,TRUE)",
                         "=1-CHISQ.DIST(2 * (" + logLikilihood_full_CM + "-" + logLikilihood_full_C + "),1,TRUE)",
                         "=1-CHISQ.DIST(2 * (" + logLikilihood_full_CM + "-" + logLikilihood_full_M + "),2,TRUE)",
                         "=1-CHISQ.DIST(2 * (" + logLikilihood_full_CM + "-" + logLikilihood_full_M + "),1,TRUE)"};
  }

  private static double getPvalue(double diffLogLikilihood, int degreeFreedom, Logger log) {
    if (diffLogLikilihood < 0) {
      return 1;
    } else {
      return ProbDist.ChiDist(diffLogLikilihood, degreeFreedom);
    }
  }



  private static double[] getPvalues(double logLikilihood_null_C, double logLikilihood_full_C,
                                     double logLikilihood_null_M, double logLikilihood_full_M,
                                     double logLikilihood_full_CM, Logger log) {

    // Lines up with even entries in EMIM_OUTPUT_FORMAT_EMIM_PVALS
    return new double[] {getPvalue(2 * (logLikilihood_full_C - logLikilihood_null_C), 2, log),
                         getPvalue(2 * (logLikilihood_full_C - logLikilihood_null_C), 1, log),
                         getPvalue(2 * (logLikilihood_full_M - logLikilihood_null_M), 2, log),
                         getPvalue(2 * (logLikilihood_full_M - logLikilihood_null_M), 1, log),
                         getPvalue(2 * (logLikilihood_full_CM - logLikilihood_full_C), 2, log),
                         getPvalue(2 * (logLikilihood_full_CM - logLikilihood_full_C), 1, log),
                         getPvalue(2 * (logLikilihood_full_CM - logLikilihood_full_M), 2, log),
                         getPvalue(2 * (logLikilihood_full_CM - logLikilihood_full_M), 1, log)};
  }

  public static String[][] loadFile(String fullPathStatResults, String nameOfMarkerColumn,
                                    String[] markersToBeLoaded, String[] columnNamesToBeLoaded,
                                    Logger log) {
    int markerColumnIndex;
    int[] indices;
    String[] line;
    String[][] result;
    BufferedReader reader;

    result = new String[markersToBeLoaded.length][columnNamesToBeLoaded.length];
    try {
      reader = new BufferedReader(new FileReader(fullPathStatResults));
      line = reader.readLine().split("\t");
      markerColumnIndex = ext.indexFactors(new String[] {nameOfMarkerColumn}, line, false, true)[0];
      indices = ext.indexFactors(columnNamesToBeLoaded, line, false, true);
      while (reader.ready()) {
        line = reader.readLine().split("\t");
        for (int i = 0; i < markersToBeLoaded.length; i++) {
          if (markersToBeLoaded[i].equalsIgnoreCase(line[markerColumnIndex])) {
            for (int j = 0; j < indices.length; j++) {
              result[i][j] = line[indices[j]];
            }
            break;
          }
        }
      }
      reader.close();
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    }
    return result;
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String resultsFile = null;
    String mapFile = "plink.bim";
    String mendelErrorFile;
    String hweFile;
    String freqFile = "plink.frq";
    String customFreqFile = null;
    String markersToReport = null;
    String outfile = null;
    String dir = "";
    String type = "plink";
    String resultsFileChild = null;
    String resultsFileMom = null;
    String resultsFileChildMom = null;
    String resultsFileTdt = null;
    Logger log;
    String logfile = null;
    double filter = 1;
    double callRateThreshold = 0;
    double pThreshold = .000001;

    String usage = "\n" + "gwas.ResultsPackager requires 0-1 arguments\n"
                   + "   (0) name of directory of all other files (i.e. dir=" + dir
                   + " (default))\n"
                   + "   (1) name of results file (i.e. results=all_results.csv (not the default))\n"
                   + "   (2) type of gwas file (i.e. type=plink (default) other options include =gwaf, =forest, =sol, and =emim)\n"
                   + "   (3) name of map file (i.e. map=" + mapFile + " (default))\n"
                   + "   (4) name of original freq file used to generate the gwaf files (i.e. freq="
                   + freqFile
                   + " (default; needs to be original freq file from when the gwaf files were made))\n"
                   + "   (5) (optional) name of freq file limited to those indiviudals used in anlayses (i.e. customFreq=femalesOnly.frq (not the default; only used in gwaf))\n"
                   + "   (6) (optional) list of markers to include (i.e. list=list.txt (not the default))\n"
                   + "   (7) (optional) name of output file (i.e. out=[results file]_out.csv (default; when I get around to coding it, it will be comma instead of tab delimited if ending in .csv))\n"
                   + "   (8) (optional) limit to those results with a p-value less than the specified filter (i.e. filter=0.001 (not the default))\n"
                   + "   (9) (optional) minimum call rate threshold (currently for SOL parser only) (i.e. callRateThreshold="
                   + callRateThreshold + " (default))\n" +

                   "";

    type = null;
    mendelErrorFile = null;
    hweFile = null;

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("results=")) {
        resultsFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("resultschild=")) {
        resultsFileChild = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("resultsmom=")) {
        resultsFileMom = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("resultschildmom=")) {
        resultsFileChildMom = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("resultstdt=")) {
        resultsFileTdt = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("resultshwe=")) {
        hweFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("resultsmendel=")) {
        mendelErrorFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("type=")) {
        type = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("map=")) {
        mapFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("freq=")) {
        freqFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("customFreq=")) {
        customFreqFile = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("list=")) {
        markersToReport = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("filter=")) {
        filter = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("callRateThreshold=")) {
        callRateThreshold = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("pthreshold=")) {
        pThreshold = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("out=")) {
        outfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    // getForestPlotParameterFile(new String[][] {{".",
    // "D:/temp/Poynter_emim/testing/allFinalWhitePoynter_results_pvals_1.xln"}, {"extragonadal",
    // "D:/temp/Poynter_emim/testing/allFinalWhitePoynter_results_pvals_2.xln"}, {"intracranial",
    // "D:/temp/Poynter_emim/testing/allFinalWhitePoynter_results_pvals_3.xln"}},
    // "D:/temp/Poynter_emim/testing/markerList.txt", "MarkerName", new String[] {"tdt", "emim"},
    // new String[][] {{"tdt_OR", "tdt_U95", "tdt_P"}, {"C_lnR1", "C_sd_lnR1", "pVal_C_df1",
    // "pVal_C_df2"}}, "D:/temp/Poynter_emim/testing/testing_forestplot.xln", null);
    // getForestPlotParameterFile(new String[][] {{".",
    // "D:/temp/Poynter_emim/testing/testing_results_pVals_df1_df2.xln"}},
    // "D:/temp/Poynter_emim/testing/markerList.txt", "MarkerName", new String[] {"tdt", "emim"},
    // new String[][] {{"tdt_OR", "tdt_U95", "tdt_P"}, {"C_lnR1", "C_sd_lnR1", "pVal_C_df1",
    // "pVal_C_df2"}}, "D:/temp/Poynter_emim/testing/testing_forestplot.xln", null);
    // System.exit(0);

    // getForestPlotParameterFile(HashVec.loadFileToStringMatrix("/home/pankrat2/shared/Poynter_emim/allFinalPoynter/fileList_allFinalPoynter.txt",
    // false, null, false), "/home/pankrat2/shared/Poynter_emim/markerList.txt", "MarkerName", new
    // String[] {"tdt", "emim"}, new String[][] {{"tdt_OR", "tdt_U95", "tdt_P"}, {"C_lnR1",
    // "C_sd_lnR1", "pVal_C_df1"}},
    // "/home/pankrat2/shared/Poynter_emim/allFinalPoynter/allFinalPoynter_forestplot.xln", null);
    // getForestPlotParameterFile(HashVec.loadFileToStringMatrix("/home/pankrat2/shared/Poynter_emim/allFinalWhitePoynter/fileList_allFinalWhitePoynter.txt",
    // false, null, false), "/home/pankrat2/shared/Poynter_emim/markerList.txt", "MarkerName", new
    // String[] {"tdt", "emim"}, new String[][] {{"tdt_OR", "tdt_U95", "tdt_P"}, {"C_lnR1",
    // "C_sd_lnR1", "pVal_C_df1"}},
    // "/home/pankrat2/shared/Poynter_emim/allFinalWhitePoynter/allFinalWhitePoynter_forestplot.xln",
    // null);
    // getForestPlotParameterFile(HashVec.loadFileToStringMatrix("/home/pankrat2/shared/Poynter_emim/completeTriosPoynter/fileList_completeTriosPoynter.txt",
    // false, null, false), "/home/pankrat2/shared/Poynter_emim/markerList.txt", "MarkerName", new
    // String[] {"tdt", "emim"}, new String[][] {{"tdt_OR", "tdt_U95", "tdt_P"}, {"C_lnR1",
    // "C_sd_lnR1", "pVal_C_df1"}},
    // "/home/pankrat2/shared/Poynter_emim/completeTriosPoynter/completeTriosPoynter_forestplot.xln",
    // null);
    // getForestPlotParameterFile(HashVec.loadFileToStringMatrix("/home/pankrat2/shared/Poynter_emim/completeWhiteTriosPoynter/fileList_completeWhiteTriosPoynter.txt",
    // false, null, false), "/home/pankrat2/shared/Poynter_emim/markerList.txt", "MarkerName", new
    // String[] {"tdt", "emim"}, new String[][] {{"tdt_OR", "tdt_U95", "tdt_P"}, {"C_lnR1",
    // "C_sd_lnR1", "pVal_C_df1"}},
    // "/home/pankrat2/shared/Poynter_emim/completeWhiteTriosPoynter/completeWhiteTriosPoynter_forestplot.xln",
    // null);
    // System.exit(0);

    // type = "emim";

    // pThreshold = .05;
    // resultsFileChild = "D:/temp/Poynter_emim/testing/emimsummary_C_allelic.out";
    // resultsFileMom = "D:/temp/Poynter_emim/testing/emimsummary_M_allelic.out";
    // resultsFileChildMom = "D:/temp/Poynter_emim/testing/emimsummary_CM_allelic.out";
    // resultsFileTdt = "D:/temp/Poynter_emim/testing/plink.tdt";
    // mapFile = "D:/temp/Poynter_emim/testing/plink_noChr23_24_25_26.bim";
    // mendelErrorFile = "D:/temp/Poynter_emim/testing/plink.lmendel";
    // hweFile = "D:/temp/Poynter_emim/testing/hardy.hwe";
    // outfile = "D:/temp/Poynter_emim/testing/testing_results_pVals_df1_df2.xln";

    // resultsFileChild = "D:/logan/emim/emim_516/emimsummary_C_1equals2.out";
    // resultsFileMom = "D:/logan/emim/emim_516/emimsummary_M_1equals2.out";
    // resultsFileChildMom = "D:/logan/emim/emim_516/emimsummary_CM_1equals2.out";
    // mapFile = "D:/logan/emim/emim_516/plink.bim";
    // outfile = "D:/logan/emim/emim_516/results_pVals_1equals2.xln";
    // resultsFileTdt = "D:/logan/emim/emim_516/plink.tdt";
    // mendelErrorFile = "D:/logan/emim/emim_516/plink.lmendel";
    // hweFile = "D:/logan/emim/emim_516/plink.hwe";

    // resultsFileChild = "D:/logan/emim/emim_276/emimsummary_C_1equals2.out";
    // resultsFileMom = "D:/logan/emim/emim_276/emimsummary_M_1equals2.out";
    // resultsFileChildMom = "D:/logan/emim/emim_276/emimsummary_CM_1equals2.out";
    // mapFile = "D:/logan/emim/emim_276/plink.bim";
    // outfile = "D:/logan/emim/emim_276/results_pVals_1equals2.xln";
    // resultsFileTdt = "D:/logan/emim/emim_276/plink.tdt";
    // mendelErrorFile = "D:/logan/emim/emim_276/plink.lmendel";
    // hweFile = "D:/logan/emim/emim_276/plink.hwe";

    // resultsFileChild = "C:/projects/Poynter_emim/allFinalPoynter/emimsummary_C_1=2.out";
    // resultsFileMom = "C:/projects/Poynter_emim/allFinalPoynter/emimsummary_M_1=2.out";
    // resultsFileChildMom = "C:/projects/Poynter_emim/allFinalPoynter/emimsummary_CM_1=2.out";
    // resultsFileTdt = "C:/projects/Poynter_emim/allFinalPoynter/plink.tdt";
    // mapFile = "C:/projects/Poynter_emim/allFinalPoynter/allFinalPoynter_noChr23_24_25_26.bim";
    // mendelErrorFile = "C:/projects/Poynter_emim/allFinalPoynter/plink.lmendel";
    // hweFile = "C:/projects/Poynter_emim/hardy.hwe";
    // outfile = "C:/projects/Poynter_emim/allFinalPoynter/allFinalPoynter_results_pVals_1=2.xln";

    // resultsFileChild =
    // "C:/projects/Poynter_emim/completeWhiteTriosPoynter/emimsummary_C_1=2.out";
    // resultsFileMom = "C:/projects/Poynter_emim/completeWhiteTriosPoynter/emimsummary_M_1=2.out";
    // resultsFileChildMom =
    // "C:/projects/Poynter_emim/completeWhiteTriosPoynter/emimsummary_CM_1=2.out";
    // resultsFileTdt = "C:/projects/Poynter_emim/completeWhiteTriosPoynter/plink.tdt";
    // mapFile =
    // "C:/projects/Poynter_emim/completeWhiteTriosPoynter/completeWhiteTriosPoynter_noChr23_24_25_26.bim";
    // mendelErrorFile = "C:/projects/Poynter_emim/completeWhiteTriosPoynter/plink.lmendel";
    // hweFile = "C:/projects/Poynter_emim/hardy.hwe";
    // outfile =
    // "C:/projects/Poynter_emim/completeWhiteTriosPoynter/completeWhiteTriosPoynter_results_pVals_1=2.xln";

    // resultsFileChild = "C:/projects/Poynter_emim/allFinalPoynter/emimsummary_C.out";
    // resultsFileMom = "C:/projects/Poynter_emim/allFinalPoynter//emimsummary_M.out";
    // resultsFileChildMom = "C:/projects/Poynter_emim/allFinalPoynter/emimsummary_CM.out";
    // resultsFileTdt = "C:/projects/Poynter_emim/allFinalPoynter/plink.tdt";
    // mapFile = "C:/projects/Poynter_emim/allFinalPoynter/allFinalPoynter_noChr23_24_25_26.bim";
    // mendelErrorFile = "C:/projects/Poynter_emim/allFinalPoynter/plink.lmendel";
    // hweFile = "C:/projects/Poynter_emim/hardy.hwe";
    // outfile = "C:/projects/Poynter_emim/allFinalPoynter/allFinalPoynter_results_pVals.xln";

    // resultsFileChild = "C:/projects/Poynter_emim/allFinalWhitePoynter/emimsummary_C.out";
    // resultsFileMom = "C:/projects/Poynter_emim/allFinalWhitePoynter/emimsummary_M.out";
    // resultsFileChildMom = "C:/projects/Poynter_emim/allFinalWhitePoynter/emimsummary_CM.out";
    // resultsFileTdt = "C:/projects/Poynter_emim/allFinalWhitePoynter/plink.tdt";
    // mapFile =
    // "C:/projects/Poynter_emim/allFinalWhitePoynter/allFinalWhitePoynter_noChr23_24_25_26.bim";
    // mendelErrorFile = "C:/projects/Poynter_emim/allFinalWhitePoynter/plink.lmendel";
    // hweFile = "C:/projects/Poynter_emim/hardy.hwe";
    // outfile =
    // "C:/projects/Poynter_emim/allFinalWhitePoynter/allFinalWhitePoynter_results_pVals.xln";

    // resultsFileChild = "C:/projects/Poynter_emim/completeTriosPoynter/emimsummary_C.out";
    // resultsFileMom = "C:/projects/Poynter_emim/completeTriosPoynter/emimsummary_M.out";
    // resultsFileChildMom = "C:/projects/Poynter_emim/completeTriosPoynter/emimsummary_CM.out";
    // resultsFileTdt = "C:/projects/Poynter_emim/completeTriosPoynter/plink.tdt";
    // mapFile =
    // "C:/projects/Poynter_emim/completeTriosPoynter/completeTriosPoynter_noChr23_24_25_26.bim";
    // mendelErrorFile = "C:/projects/Poynter_emim/completeTriosPoynter/plink.lmendel";
    // hweFile = "C:/projects/Poynter_emim/hardy.hwe";
    // outfile =
    // "C:/projects/Poynter_emim/completeTriosPoynter/completeTriosPoynter_results_pVals.xln";

    // resultsFileChild = "C:/projects/Poynter_emim/completeWhiteTriosPoynter/emimsummary_C.out";
    // resultsFileMom = "C:/projects/Poynter_emim/completeWhiteTriosPoynter/emimsummary_M.out";
    // resultsFileChildMom =
    // "C:/projects/Poynter_emim/completeWhiteTriosPoynter/emimsummary_CM.out";
    // resultsFileTdt = "C:/projects/Poynter_emim/completeWhiteTriosPoynter/plink.tdt";
    // mapFile =
    // "C:/projects/Poynter_emim/completeWhiteTriosPoynter/completeWhiteTriosPoynter_noChr23_24_25_26.bim";
    // mendelErrorFile = "C:/projects/Poynter_emim/completeWhiteTriosPoynter/plink.lmendel";
    // hweFile = "C:/projects/Poynter_emim/hardy.hwe";
    // outfile =
    // "C:/projects/Poynter_emim/completeWhiteTriosPoynter/completeWhiteTriosPoynter_results_pVals.xln";

    // resultsFileChild = "C:/projects/Poynter_emim/allFinalPoynterNoSibs/emimsummary_C.out";
    // resultsFileMom = "C:/projects/Poynter_emim/allFinalPoynterNoSibs/emimsummary_M.out";
    // resultsFileChildMom = "C:/projects/Poynter_emim/allFinalPoynterNoSibs/emimsummary_CM.out";
    // resultsFileTdt = "C:/projects/Poynter_emim/allFinalPoynterNoSibs/plink.tdt";
    // mapFile =
    // "C:/projects/Poynter_emim/allFinalPoynterNoSibs/allFinalPoynterNoSibs_noChr23_24_25_26.bim";
    // mendelErrorFile = "C:/projects/Poynter_emim/allFinalPoynterNoSibs/plink.lmendel";
    // hweFile = "C:/projects/Poynter_emim/hardy.hwe";
    // outfile =
    // "C:/projects/Poynter_emim/allFinalPoynterNoSibs/allFinalPoynterNoSibs_results_pVals.xln";

    try {
      // parseBP();
      // parseIBCFormatFromGWAF("D:/BOSS/GWAF/reclustered/PTA/", "pta_HmPCs_results.csv",
      // "plink.map", "plink.frq", null, "list.txt", "PTA_results.txt");
      // parseIBCFormatFromGWAF("D:/CARe/CARe_geno_data_and_misc/IBC/FHS/iSELECT/gwaf/",
      // "pta_fhs_parsedResults.csv", "grk100.map", "grk100.frq", null, "list.txt",
      // "PTA_fhs_results.txt");

      // parseStdFormatFromPlink("D:/Myron/CALICO/T2DM/", "t2dm.assoc.logistic", "ADD", "plink.bim",
      // "t2dm.frq", null, "cardia_page_t2dm_results.txt", new Logger());
      // parseStdFormatFromPlink("D:/Myron/CALICO/AgeMenarche/", "menarche1.assoc.linear", "ADD",
      // "plink.bim", "menarche.frq", null, "cardia_page_menarche1_results.txt", new Logger());
      // parseStdFormatFromPlink("D:/Myron/CALICO/AgeMenarche/", "menarche2.assoc.linear", "ADD",
      // "plink.bim", "menarche.frq", null, "cardia_page_menarche2_results.txt", new Logger());
      // System.exit(1);

      // parseSOLformat("D:/data/SOL/", "MODEL3slim_b.out",
      // "N:/statgen/CALICo_SOL/SOL-2013-04-05_Metabochip-mappingfile.txt", null, null, 1.0, 0.95,
      // "MODEL3slim_b_test.out", new Logger());
      // System.exit(1);

      log = new Logger(logfile);
      if (type.equalsIgnoreCase("gwaf")) {
        parseIBCFormatFromGWAF(dir, resultsFile, mapFile, freqFile, customFreqFile, markersToReport,
                               filter, outfile, log);
      } else if (type.equalsIgnoreCase("plink")) {
        parseStdFormatFromPlink(dir, resultsFile, "Add", mapFile, freqFile, markersToReport, filter,
                                outfile, log);
      } else if (type.equalsIgnoreCase("sol")) {
        parseSOLformat(dir, resultsFile,
                       "N:/statgen/CALICo_SOL/SOL-2013-04-05_Metabochip-mappingfile.txt", freqFile,
                       markersToReport, filter, callRateThreshold, outfile, log);
      } else if (type.equalsIgnoreCase("emim")) {
        parseEmimFormat(resultsFileChild, resultsFileMom, resultsFileChildMom, resultsFileTdt,
                        mapFile, mendelErrorFile, hweFile, freqFile, pThreshold, outfile, log);
      } else if (type.equalsIgnoreCase("forest")) {
        String mkrFile = "/home/pankrat2/shared/Poynter_emim/gwasHits.txt";
        String mkrColNm = "markerName";
        String[] analyses = {"tdt", "emim_child", "emim_maternal"};
        String[][] analysisNms =
            {{"tdt_OR", "tdt_U95", "tdt_P"}, {"C_lnR1", "C_se_lnR1", "pVal_C_df1"},
             {"CM_lnS1", "CM_se_lnS1", "pVal_CM-C_df1"}};

        String[][] files =
            {{"/home/pankrat2/shared/Poynter_emim/allFinalWhitePoynter/fileList_allFinalWhitePoynter.txt",
              "/home/pankrat2/shared/Poynter_emim/allFinalWhitePoynter/allFinalWhitePoynter_forestplot.xln"},
             {"/home/pankrat2/shared/Poynter_emim/completeTriosPoynter/fileList_completeTriosPoynter.txt",
              "/home/pankrat2/shared/Poynter_emim/completeTriosPoynter/completeTriosPoynter_forestplot.xln"},
             {"/home/pankrat2/shared/Poynter_emim/allFinalPoynter/fileList_allFinalPoynter.txt",
              "/home/pankrat2/shared/Poynter_emim/allFinalPoynter/allFinalPoynter_forestplot.xln"},
             {"/home/pankrat2/shared/Poynter_emim/completeWhiteTriosPoynter/fileList_completeWhiteTriosPoynter.txt",
              "/home/pankrat2/shared/Poynter_emim/completeWhiteTriosPoynter/completeWhiteTriosPoynter_forestplot.xln"},};

        boolean oddsRatio = true;
        String sortFileName = "/home/pankrat2/shared/Poynter_emim/forestPlotDisplayOrder.txt";
        for (String[] fileSet : files) {
          getForestPlotParameterFile(HashVec.loadFileToStringMatrix(fileSet[0], false, null, false),
                                     mkrFile, mkrColNm, analyses, analysisNms, null, fileSet[1],
                                     null);
          org.genvisis.cnv.plots.ForestPlot fp =
              new org.genvisis.cnv.plots.ForestPlot(ext.rootOf(fileSet[1], false) + ".input", null);
          fp.setOddsRatioDisplay(oddsRatio);
          fp.loadOrderFile(sortFileName, true);
          fp.screenCapAll("forestPlots", oddsRatio, false);
        }
        String[] dirs = {"/home/pankrat2/shared/Poynter_emim/allFinalWhitePoynter/",
                         "/home/pankrat2/shared/Poynter_emim/completeTriosPoynter/",
                         "/home/pankrat2/shared/Poynter_emim/allFinalPoynter/",
                         "/home/pankrat2/shared/Poynter_emim/completeWhiteTriosPoynter/"};
        for (String rundir : dirs) {
          String fileOut = ext.rootOf(mkrFile, true);
          String temp = rundir.substring(0, rundir.length() - 1);
          temp = temp.substring(temp.lastIndexOf('/') + 1);
          fileOut += "_" + temp + "_plink_tdtWithCi.tdt";
          Files.filterByKeys(mkrFile, rundir + "plink_tdtWithCi.tdt", rundir + fileOut, 2, true);
        }

      } else {
        System.err.println("Error - unknown results type: '" + type + "'");
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void parseBP() {
    String[] phenos = {"DBP", "HTN", "MAP", "PP", "SBP"};
    String[] types = {"M", "F"};
    String dir;

    dir = "C:/Users/npankrat/Documents/BloodPressure/";
    dir = "D:/BOSS/IBC_meta_analyses/BloodPressure/";
    for (String type : types) {
      for (String pheno : phenos) {
        // parseIBCFormatFromGWAF(dir, phenos[j]+"."+types[i]+".results.csv", "plink.map",
        // "PLINK.FRQ", null, null, "BOSS-EHLS."+phenos[j]+"."+types[i]+".txt");
        parseIBCFormatFromGWAF(dir, pheno + "." + type + ".results.csv", "plink.map", "plink.frq",
                               null, "list.txt", 1,
                               "BOSS-EHLS_" + pheno + "_" + type + "_072512.txt", new Logger());
        System.out.print(";BOSS-EHLS_" + pheno + "_" + type + "_072512.txt,11");
      }
    }
  }

  public static void parseEmimFormat(String childResultsFile, String momResultsFile,
                                     String childMomResultsFile, String tdtResultsFile,
                                     String mapFile, String mendelErrorFile, String hweFile,
                                     String frqFile, double pValueThreshold, String outfile,
                                     Logger log) {
    BufferedReader reader1, reader2, reader3;
    PrintWriter writer;
    String[] lineC, lineM, lineCM, pvalEquations;
    String temp;
    Hashtable<Long, String> snpList; // , freqHash; // , customFreqHash;
    String delimiter1, delimiter2, delimiter3;
    Map<String, Integer> indicesC, indicesM, indicesCM;
    double hweThreshold;
    double[] pvals;
    long index;
    Hashtable<String, String[]> mendelErrors = null, hwe = null, tdtResults = null, freqs = null;

    if (outfile == null) {
      outfile = ext.rootOf(childResultsFile, false) + "parsedResults.txt";
    }

    if (!Files.exists(childResultsFile)) {
      log.reportError("Error: could not find results file '" + childResultsFile + "'");
      return;
    }
    try {
      snpList = new Hashtable<Long, String>();
      index = 1;
      reader1 = Files.getAppropriateReader(mapFile);
      while (reader1.ready()) {
        lineC = reader1.readLine().split("\t");
        snpList.put(index, lineC[0] + "\t" + lineC[3] + "\t" + lineC[1] + "\t" + lineC[4] + "\t"
                           + lineC[5]);
        index++;
      }
      reader1.close();

      if (mendelErrorFile != null) {
        mendelErrors = SkatMeta2.loadFile(mendelErrorFile, null, new String[] {"SNP"},
                                          new String[] {"N"}, null, null);
      }
      if (hweFile != null) {
        hwe = SkatMeta2.loadFile(hweFile, null, new String[] {"SNP"}, new String[] {"GENO", "p"},
                                 new String[] {"TEST==UNAFF"}, log);
      }
      if (frqFile != null) {
        freqs = SkatMeta2.loadFile(frqFile, null, new String[] {"SNP"}, new String[] {"MAF"}, null,
                                   log);
      }
      if (tdtResultsFile != null) {
        // tdtResults = one.SkatMeta.loadFile(tdtResultsFile, null, new String[] {"SNP"}, new
        // String[] {"T", "U", "OR", "P"}, null, null);
        tdtResults =
            SkatMeta2.loadFile(tdtResultsFile, null, new String[] {"SNP"}, TDT_REQS, null, null);
      }

      hweThreshold = 0.05 / Files.countLines(childResultsFile, 1);
      log.report("");
      reader1 = Files.getAppropriateReader(childResultsFile);
      reader3 = Files.getAppropriateReader(childMomResultsFile);
      reader2 = Files.getAppropriateReader(momResultsFile);
      writer = Files.getAppropriateWriter(outfile);
      temp = reader1.readLine().trim();
      delimiter1 = ext.determineDelimiter(temp);
      lineC = temp.split(delimiter1);
      indicesC = ext.indexMap(EMIM_REQS, lineC, false, log, false, false);
      temp = reader2.readLine().trim();
      delimiter2 = ext.determineDelimiter(temp);
      lineM = temp.split(delimiter2);
      indicesM = ext.indexMap(EMIM_REQS, lineM, false, log, false, false); // TODO EMIM_REQS
      temp = reader3.readLine().trim();
      delimiter3 = ext.determineDelimiter(temp);
      lineCM = temp.split(delimiter3);
      indicesCM = ext.indexMap(EMIM_REQS, lineCM, false, log, false, false); // TODO EMIM_REQS

      writer.println(Array.toStr(EMIM_OUTPUT_FORMAT_SNPS)
                     + (mendelErrorFile == null ? ""
                                                : ("\t"
                                                   + Array.toStr(EMIM_OUTPUT_FORMAT_MENDEL_ERRORS)))
                     + (hweFile == null ? "" : ("\t" + Array.toStr(EMIM_OUTPUT_FORMAT_HWE)))
                     + (tdtResultsFile == null ? "" : ("\t" + Array.toStr(EMIM_OUTPUT_FORMAT_TDT)))
                     + "\t" + Array.toStr(EMIM_OUTPUT_FORMAT_EMIM_RESULTS) + "\t"
                     + Array.toStr(EMIM_OUTPUT_FORMAT_EMIM_PVALS));
      while (reader1.ready()) {
        // lineC = reader1.readLine().replaceAll("\\*", " ").trim().split(delimiter1, -1);
        lineC =
            reader1.readLine()
                   .replaceAll("\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*", "            NaN")
                   .trim().split(delimiter1, -1);
        String snp = snpList
                            .get(Long.parseLong(lineC[indicesC.get("snpID")].substring(0,
                                                                                       lineC[indicesC.get("snpID")].indexOf("."))))
                            .split("\t")[2];
        String freq;
        if (freqs == null) {
          freq = lineC[indicesC.get("freq")];
        } else {
          freq = freqs.get(snp)[0];
        }
        // lineM = reader2.readLine().replaceAll("\\*", " ").trim().split(delimiter2, -1);
        // lineCM = reader3.readLine().replaceAll("\\*", " ").trim().split(delimiter3, -1);
        lineM =
            reader2.readLine()
                   .replaceAll("\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*", "            NaN")
                   .trim().split(delimiter2, -1);
        lineCM =
            reader3.readLine()
                   .replaceAll("\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*", "            NaN")
                   .trim().split(delimiter3, -1);

        if (!lineC[indicesC.get("snpID")].equals(lineM[indicesM.get("snpID")])
            || !lineC[indicesC.get("snpID")].equals(lineCM[indicesCM.get("snpID")])) {
          log.reportError("Error - SNP ID in the files are not lined up. Child SNP ID: "
                          + lineC[indicesC.get("snpID")] + "; Mom SNP ID: "
                          + lineM[indicesM.get("snpID")] + "; ChildMom SNP ID: "
                          + lineCM[indicesCM.get("snpID")] + ".");
          return;
        }
        pvals = getPvalues(Double.parseDouble(lineC[indicesC.get("lnliknull")]),
                           Double.parseDouble(lineC[indicesC.get("lnlikfull")]),
                           Double.parseDouble(lineM[indicesM.get("lnliknull")]),
                           Double.parseDouble(lineM[indicesM.get("lnlikfull")]),
                           Double.parseDouble(lineCM[indicesCM.get("lnlikfull")]), log);
        temp = tdtResults.get(snp)[3];
        // if (tdtResults.get(temp)[3].equals("NA")) {
        // System.out.println(temp + "\t" + tdtResults.get(test)[3]);
        //// System.out.println(temp + "\t" + Double.parseDouble(tdtResults.get(test)[3]));
        // }
        if (pValueThreshold >= 1
            || ((!temp.equals("NA")
                 && pvalsPass(Array.addDoubleToArray(Double.parseDouble(temp), pvals),
                              pValueThreshold))
                && Double.parseDouble(freq) >= .01)) {
          pvalEquations =
              getPvalEquations(lineC[indicesC.get("lnliknull")], lineC[indicesC.get("lnlikfull")],
                               lineM[indicesM.get("lnliknull")], lineM[indicesM.get("lnlikfull")],
                               lineCM[indicesCM.get("lnlikfull")], log);
          // writer.println(getOutputString(snpList, lineC, indicesC, lineM, indicesM, lineCM,
          // indicesCM, log) + "\t" + pvals[0] + "\t" + pvalEquations[0] + "\t" + pvals[1] + "\t" +
          // pvalEquations[1] + "\t" + pvals[2] + "\t" + pvalEquations[2]);
          // writer.println(getOutputString(snpList, mendelErrors, hwe, hweThreshold, lineC,
          // indicesC, lineM, indicesM, lineCM, indicesCM, tdtResults, log) + "\t" + pvals[0] + "\t"
          // + pvalEquations[0] + "\t" + pvals[1] + "\t" + pvalEquations[1] + "\t" + pvals[2] + "\t"
          // + pvalEquations[2]);
          writer.println(getOutputString(snpList, mendelErrors, hwe, hweThreshold, freq, lineC,
                                         indicesC, lineM, indicesM, lineCM, indicesCM, tdtResults,
                                         pvals, pvalEquations, log));
        }
      }
      reader1.close();
      reader3.close();
      reader2.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: one of the files was not found within the directory");
      fnfe.printStackTrace();
      return;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + childResultsFile + "\"");
      return;
    }

    log.report("Parsed Emim result is ready at: " + outfile);
  }

  public static void parseIBCFormatFromGWAF(String dir, String resultsFile, String mapFile,
                                            String originalFrqFile, String customFrqFile,
                                            String markersToReport, double filter, String outfile,
                                            Logger log) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String temp, trav;
    HashSet<String> markerHash;
    Hashtable<String, String> mapHash, originalFreqHash, customFreqHash;
    String delimiter;
    String[] alleles;
    String freq;

    if (outfile == null) {
      outfile = ext.rootOf(resultsFile, false) + "_out.csv";
    }

    if (markersToReport != null) {
      markerHash = HashVec.loadFileToHashSet(dir + markersToReport, false);
    } else {
      markerHash = null;
    }

    mapHash = HashVec.loadFileToHashString(dir + mapFile, new int[] {1}, new int[] {0, 3}, false,
                                           "\t", false, false, false);
    originalFreqHash =
        HashVec.loadFileToHashString(dir + originalFrqFile, new int[] {1}, new int[] {2, 3}, false,
                                     "\t", false, false, false); // add 4 if you want global
                                                                 // frequency

    if (customFrqFile != null) {
      System.err.println("Warning - use of custom freq file has not been tested properly; if it works then remove this warning");
      customFreqHash =
          HashVec.loadFileToHashString(originalFrqFile, new int[] {1}, new int[] {2, 3, 4}, false,
                                       "\t", false, false, false); // add 4 if you want global
                                                                   // frequency instead of custom
                                                                   // Freq
    } else {
      customFreqHash = null;
    }

    try {
      reader = Files.getAppropriateReader(dir + resultsFile);
      writer = Files.getAppropriateWriter(dir + outfile);
      temp = reader.readLine().trim();
      delimiter = ext.determineDelimiter(temp);
      line = temp.split(delimiter);
      ext.checkHeader(line, GWAF.HEADER_SUMMARY, true);
      // writer.println(Array.toStr(IBC_OUTPUT_FORMAT));
      writer.println(Array.toStr(TRADITIONAL_OUTPUT_FORMAT));
      while (reader.ready()) {
        line = reader.readLine().trim().split(delimiter);
        trav = line[0];
        if ((markerHash == null || markerHash.contains(trav)) && !line[3].equals("")
            && (filter >= 1
                || (!ext.isMissingValue(line[3]) && Double.parseDouble(line[3]) <= filter))) {
          if (mapHash.containsKey(trav)) {
            writer.print(mapHash.get(trav));
          } else if (mapHash.containsKey(ext.replaceAllWith(trav, ".", "-"))) {
            writer.print(mapHash.get(ext.replaceAllWith(trav, ".", "-")));
          } else {
            log.reportError("Error - no map position for " + trav);
            writer.print(".\t.");
          }
          writer.print("\t" + trav + "\tNA\t+\t" + line[4]);

          if (originalFreqHash.containsKey(trav)) {
            alleles = originalFreqHash.get(trav).split("[\\s]+");
          } else if (originalFreqHash.containsKey(ext.replaceAllWith(trav, ".", "-"))) {
            alleles = originalFreqHash.get(ext.replaceAllWith(trav, ".", "-")).split("[\\s]+");
          } else {
            alleles = new String[] {".", "."};
            log.reportError("Error - no alleles from original .frq file for " + trav);
          }
          writer.print("\t" + Array.toStr(alleles));

          freq = line[5];
          if (customFreqHash != null) {
            if (customFreqHash.containsKey(trav)) {
              freq =
                  Alleles.getAlleleFreqForA1(alleles[0], customFreqHash.get(trav).split("\t"))[2];
            } else if (customFreqHash.containsKey(ext.replaceAllWith(trav, ".", "-"))) {
              freq =
                  Alleles.getAlleleFreqForA1(alleles[0],
                                             customFreqHash.get(ext.replaceAllWith(trav, ".", "-"))
                                                           .split("\t"))[2];
            } else {
              log.reportError("Error - no alleles from custom .frq file for " + trav);
              freq = ".";
            }
          }
          writer.println("\t" + freq + "\t" + line[1] + "\t" + line[2] + "\t" + line[3]);
        }
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + resultsFile + "\" not found in current directory");
      return;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + resultsFile + "\"");
      return;
    }
  }

  public static void parseSOLformat(String dir, String resultsFile, String mapFile, String freqFile,
                                    String markersToReport, double filter, double callRateThreshold,
                                    String outfile, Logger log) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String temp, trav;
    HashSet<String> markerHash;
    Hashtable<String, String> mapHash; // , freqHash; // , customFreqHash;
    String delimiter;
    int[] indices;
    boolean logistic;
    int nnmiss, mac;
    double callrate;
    Vector<String> lowCallrateMarkers;

    dir = ext.verifyDirFormat(dir);

    if (outfile == null) {
      outfile = ext.rootOf(resultsFile, false) + ".out";
    }

    if (markersToReport != null) {
      markerHash = HashVec.loadFileToHashSet(dir + markersToReport, false);
    } else {
      markerHash = null;
    }

    if (Files.exists(dir + mapFile)) {
      mapFile = dir + mapFile;
    } else if (!Files.exists(mapFile)) {
      log.reportError("Error: could not find map file '" + mapFile
                      + "' in its absolute path or in '" + dir + "'");
      return;
    }
    // if (!Files.exists(dir+freqFile)) {
    // log.reportError("Error: could not find freq file '"+dir+freqFile+"'");
    // return;
    // }
    if (!Files.exists(dir + resultsFile)) {
      log.reportError("Error: could not find results file '" + dir + resultsFile + "'");
      return;
    }
    mapHash = HashVec.loadFileToHashString(mapFile, new int[] {1}, new int[] {0, 2, 3, 4}, false,
                                           "\t", false, false, false);

    // TODO
    // freqHash = HashVec.loadFileToHashString(dir+freqFile, new int[] {1}, new int[] {2,3,4},
    // false, "\t", false, false, false); // 4 gets global frequency

    try {
      reader = Files.getAppropriateReader(dir + resultsFile);
      new File(ext.parseDirectoryOfFile(dir + outfile)).mkdirs();
      writer = Files.getAppropriateWriter(dir + outfile);
      temp = reader.readLine().trim();
      delimiter = ext.determineDelimiter(temp);
      line = temp.split(delimiter);
      indices = ext.indexFactors(SOL_REQS, line, false, log, false, false);

      // TODO: Need to revisit after logistic has been implemented
      if (indices[4] == -1 && indices[5] == -1) {
        log.reportError("Error - results file did not contain a column for 'OR' (logistic) or 'BETA' (linear); aborting");
        return;
        // } else if (indices[4] != -1 && indices[5] != -1) {
        // log.reportError("Error - results file contain a column for both 'OR' (logistic) and
        // 'BETA' (linear); aborting");
        // return;
        // } else if (indices[4] != -1) {
        // logistic = true;
      } else {
        logistic = false;
      }

      if (logistic) {
        log.reportError("Logistic parsing has not yet been implemented; check the column order etc, parse it properly and re-run");
        return;
      }

      if (indices[6] == -1) {
        log.reportError("Warning - results file did not contain a column for 'StdErr/SE'; values will be set to NA");
      }

      lowCallrateMarkers = new Vector<String>();
      writer.println(Array.toStr(STANDARD_OUTPUT_FORMAT));
      while (reader.ready()) {
        line = reader.readLine().trim().split(delimiter);
        trav = line[indices[0]];
        nnmiss = Integer.parseInt(line[indices[6]]) + Integer.parseInt(line[indices[7]])
                 + Integer.parseInt(line[indices[8]]);
        callrate = (double) nnmiss / (double) (nnmiss + Integer.parseInt(line[indices[9]]));
        if (callrate <= callRateThreshold) {
          lowCallrateMarkers.add(trav);
        }

        if ((markerHash == null || markerHash.contains(trav)) && (callrate > callRateThreshold)
            && (filter >= 1 || (!ext.isMissingValue(line[indices[3]])
                                && Double.parseDouble(line[indices[3]]) <= filter))) {
          writer.print(trav); // MarkerName
          if (mapHash.containsKey(trav)) {
            writer.print("\t" + mapHash.get(trav)); // chr, pos, A1, A2
          } else {
            log.reportError("Error - no map position for " + trav);
            writer.print("\t.\t.\t.\t."); // null, null
          }
          // if (line[indices[1]].equals("NA")) {
          // writer.print("\t"+line[indices[1]]+"\t0\t0"); // missing data, null, null
          // } else {
          //
          // if (freqHash.containsKey(trav)) {
          // writer.print("\t"+Array.toStr(Alleles.getAlleleFreqForA1(line[indices[1]],
          // freqHash.get(trav).split("\t")))); // a1, a2, a1_freq
          // } else if (freqHash.containsKey(ext.replaceAllWith(trav, ".", "-"))) {
          // writer.print("\t"+Array.toStr(Alleles.getAlleleFreqForA1(line[indices[1]],
          // freqHash.get(ext.replaceAllWith(trav, ".", "-")).split("\t")))); // a1, a2, a1_freq
          // } else {
          // log.reportError("Error - no frequency for "+trav);
          // writer.print("\t"+line[indices[1]]+"\t.\t."); // a1, null, null
          // }
          //
          // }

          // writer.print("\t"+line[indices[9]]); // NMISS here means number missing not, "non
          // missing"

          nnmiss = Integer.parseInt(line[indices[6]]) + Integer.parseInt(line[indices[7]])
                   + Integer.parseInt(line[indices[8]]);
          mac = Integer.parseInt(line[indices[7]]) + Integer.parseInt(line[indices[8]]) * 2;

          writer.print("\t" + ext.formDeci((double) mac / (double) nnmiss / 2.0, 6));
          writer.print("\t" + nnmiss);

          // TODO logistic has not been implemented yet, so this needs to be revisited when it is
          if (logistic) {
            log.reportError("Logistic parsing has not yet been implemented; check the column order etc, parse it properly and re-run");
            if (line[indices[1]].equals("NA")) {
              writer.print("\t" + line[indices[1]]); // NA OR -> NA beta
            } else {
              writer.print("\t"
                           + ext.formDeci(Math.log(Double.parseDouble(line[indices[1]])), 5, true)); // OR
                                                                                                     // ->
                                                                                                     // beta
            }
          } else {
            writer.print("\t" + line[indices[1]]); // beta
          }
          writer.println("\t" + (indices[2] == -1 ? "NA" : line[indices[2]]) + "\t"
                         + line[indices[3]]); // se,
                                              // pval
        }
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + resultsFile + "\" not found in current directory");
      return;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + resultsFile + "\"");
      return;
    }

    Files.writeList(Array.toStringArray(lowCallrateMarkers),
                    dir + outfile + "lowCallRateMarkers.out");
  }

  public static void parseStdFormatFromPlink(String dir, String resultsFile, String test,
                                             String mapFile, String freqFile,
                                             String markersToReport, double filter, String outfile,
                                             Logger log) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String temp, trav;
    HashSet<String> markerHash;
    Hashtable<String, String> mapHash, freqHash; // , customFreqHash;
    String delimiter;
    int[] indices;
    boolean logistic;

    if (outfile == null) {
      outfile = ext.rootOf(resultsFile, false) + ".out";
    }

    if (markersToReport != null) {
      if (!Files.exists(dir + markersToReport)) {
        log.reportError("Error: could not find markers to report file '" + dir + markersToReport
                        + "'");
        return;
      }
      markerHash = HashVec.loadFileToHashSet(dir + markersToReport, false);
    } else {
      markerHash = null;
    }

    if (!Files.exists(dir + mapFile)) {
      log.reportError("Error: could not find map file '" + dir + mapFile + "'");
      return;
    }
    if (!Files.exists(dir + freqFile)) {
      log.reportError("Error: could not find freq file '" + dir + freqFile + "'");
      return;
    }
    if (!Files.exists(dir + resultsFile)) {
      log.reportError("Error: could not find results file '" + dir + resultsFile + "'");
      return;
    }
    mapHash = HashVec.loadFileToHashString(dir + mapFile, new int[] {1}, new int[] {0, 3}, false,
                                           "\t", false, false, false);
    freqHash = HashVec.loadFileToHashString(dir + freqFile, new int[] {1}, new int[] {2, 3, 4},
                                            false, "\t", false, false, false); // 4 gets global
                                                                               // frequency

    try {
      reader = Files.getAppropriateReader(dir + resultsFile);
      new File(ext.parseDirectoryOfFile(dir + outfile)).mkdirs();
      writer = Files.getAppropriateWriter(dir + outfile);
      temp = reader.readLine().trim();
      delimiter = ext.determineDelimiter(temp);
      line = temp.split(delimiter);
      indices = ext.indexFactors(PLINK_REQS, line, false, log, false, false);
      if (indices[4] == -1 && indices[5] == -1) {
        log.reportError("Error - results file did not contain a column for 'OR' (logistic) or 'BETA' (linear); aborting");
        return;
      } else if (indices[4] != -1 && indices[5] != -1) {
        log.reportError("Error - results file contain a column for both 'OR' (logistic) and 'BETA' (linear); aborting");
        return;
      } else if (indices[4] != -1) {
        logistic = true;
      } else {
        logistic = false;
      }

      if (indices[6] == -1) {
        log.reportError("Warning - results file did not contain a column for 'StdErr/SE'; values will be set to NA");
      }


      writer.println(Array.toStr(STANDARD_OUTPUT_FORMAT));
      while (reader.ready()) {
        line = reader.readLine().trim().split(delimiter);
        trav = line[indices[0]];
        if ((markerHash == null || markerHash.contains(trav)) && !line[3].equals("")
            && line[indices[2]].equalsIgnoreCase(test)
            && (filter >= 1 || (!ext.isMissingValue(line[indices[7]])
                                && Double.parseDouble(line[indices[7]]) <= filter))) {
          writer.print(trav); // MarkerName
          if (mapHash.containsKey(trav)) {
            writer.print("\t" + mapHash.get(trav)); // chr, pos
          } else if (mapHash.containsKey(ext.replaceAllWith(trav, ".", "-"))) {
            writer.print("\t" + mapHash.get(ext.replaceAllWith(trav, ".", "-"))); // chr, pos
          } else {
            log.reportError("Error - no map position for " + trav);
            writer.print("\t.\t."); // null, null
          }
          if (line[indices[1]].equals("NA")) {
            writer.print("\t" + line[indices[1]] + "\t0\t0"); // missing data, null, null
          } else {
            // try {
            if (freqHash.containsKey(trav)) {
              writer.print("\t" + Array.toStr(Alleles.getAlleleFreqForA1(line[indices[1]],
                                                                         freqHash.get(trav)
                                                                                 .split("\t")))); // a1,
                                                                                                  // a2,
                                                                                                  // a1_freq
            } else if (freqHash.containsKey(ext.replaceAllWith(trav, ".", "-"))) {
              writer.print("\t" + Array.toStr(Alleles.getAlleleFreqForA1(line[indices[1]],
                                                                         freqHash.get(ext.replaceAllWith(trav,
                                                                                                         ".",
                                                                                                         "-"))
                                                                                 .split("\t")))); // a1,
                                                                                                  // a2,
                                                                                                  // a1_freq
            } else {
              log.reportError("Error - no frequency for " + trav);
              writer.print("\t" + line[indices[1]] + "\t.\t."); // a1, null, null
            }
            // } catch (Exception e) {
            // System.out.println("hola");
            // }
          }
          writer.print("\t" + line[indices[3]]); // NMISS
          if (logistic) {
            if (line[indices[4]].equals("NA")) {
              writer.print("\t" + line[indices[4]]); // NA OR -> NA beta
            } else {
              writer.print("\t"
                           + ext.formDeci(Math.log(Double.parseDouble(line[indices[4]])), 5, true)); // OR
                                                                                                     // ->
                                                                                                     // beta
            }
          } else {
            writer.print("\t" + line[indices[5]]); // beta
          }
          // writer.println("\t"+line[indices[6]]+"\t"+line[indices[7]]); // se, pval
          writer.println("\t" + (indices[6] == -1 ? "NA" : line[indices[6]]) + "\t"
                         + line[indices[7]]); // se,
                                              // pval
        }
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + resultsFile + "\" not found in current directory");
      return;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + resultsFile + "\"");
      return;
    }
  }

  private static boolean pvalsPass(double[] pvals, double thresh) {
    for (double pval : pvals) {
      if (pval <= thresh) {
        return true;
      }
    }
    return false;
  }
}
