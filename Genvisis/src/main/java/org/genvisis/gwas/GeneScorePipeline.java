package org.genvisis.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.stream.Collectors;
import org.genvisis.bioinformatics.MapSNPsAndGenes;
import org.genvisis.bioinformatics.Sequence;
import org.genvisis.common.Aliases;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.DosageData;
import org.genvisis.filesys.SnpMarkerSet;
import org.genvisis.gwas.MergeExtractPipeline.DataSource;
import org.genvisis.gwas.parsing.AbstractColumnFilter;
import org.genvisis.gwas.parsing.ColumnFilter;
import org.genvisis.gwas.parsing.DataLine;
import org.genvisis.gwas.parsing.FileColumn;
import org.genvisis.gwas.parsing.FileParser;
import org.genvisis.gwas.parsing.FileParserFactory;
import org.genvisis.gwas.parsing.StandardFileColumns;
import org.genvisis.seq.manage.StrandOps;
import org.genvisis.seq.manage.StrandOps.AlleleOrder;
import org.genvisis.stats.ProbDist;
import org.genvisis.stats.RegressionModel;
import com.google.common.base.Joiner;
import com.google.common.primitives.Doubles;

public class GeneScorePipeline {

  private static final String ARG_WORKING_DIR = "workDir=";
  private static final String ARG_INDEX_THRESH = "indexThresh=";
  private static final String ARG_WINDOW_SIZE = "minWinSize=";
  private static final String ARG_WINDOW_EXT = "winThresh=";
  private static final String ARG_MISS_THRESH = "missThresh=";

  private static float DEFAULT_INDEX_THRESHOLD = (float) 0.00000005;
  private static int DEFAULT_WINDOW_MIN_SIZE_PER_SIDE = 500000;// 500kb each side is technically a
                                                               // 1M window until the next hit
                                                               // region, but we now take this into
                                                               // consideration in the main
                                                               // algorithm
  private static float DEFAULT_WINDOW_EXTENSION_THRESHOLD = (float) 0.000005; // (float)0.00001;
  private static String[] DEFAULT_ADDL_ANNOT_VAR_NAMES = new String[0];
  private static double DEFAULT_MIN_MISS_THRESH = 0.5;

  private static final String PLINK_FRQ_DIR = "N:/statgen/1000G_work/Frequencies/";
  private static final String TAG = "##";
  private static final String PLINK_FRQ_FILE_PATTERN = "chr" + TAG + "_eu_unrel.frq.xln";
  private static final String CROSS_FILTERED_DATAFILE = "bimData.xln";
  private static final String DATA_SOURCE_FILENAME = "data.txt";

  private static final FileColumn<String> MARKER_COL = StandardFileColumns.snp("Marker");
  private static final FileColumn<String> EFFECT_ALLELE_COL = StandardFileColumns.a1("EffectAllele");
  private static final FileColumn<String> NON_EFFECT_ALLELE_COL = StandardFileColumns.a2("NonEffectAllele");
  private static final FileColumn<Double> BETA_COL = StandardFileColumns.beta("beta");

  private static final String[][] LINKERS = {Aliases.MARKER_NAMES, Aliases.ALLELES[0],
                                             Aliases.EFFECTS};

  private static final String REGRESSION_HEADER = "STUDY\tDATAFILE\tINDEX-THRESHOLD\tFACTOR\tBASE-R-SQR\tR-SQR\tR-DIFF\tP-VALUE\tBETA\tSE\tNUM\t#DATASNPs\t#PLINKSNPs\t#HITSNPs\tB-F-SCORE\tINVCHI-SCORE\tEXCEL-SIG";
  private final String metaDir;

  private float[] indexThresholds = new float[] {DEFAULT_INDEX_THRESHOLD};
  private int[] windowMinSizePerSides = new int[] {DEFAULT_WINDOW_MIN_SIZE_PER_SIDE};
  private float[] windowExtensionThresholds = new float[] {DEFAULT_WINDOW_EXTENSION_THRESHOLD};
  private double minMissThresh = DEFAULT_MIN_MISS_THRESH;

  // private int numThreads = 1;
  // private boolean runPlink = false;
  // private boolean runRegression = false;
  // private boolean writeHist = false;

  private final ArrayList<String> metaFiles = new ArrayList<String>();
  private final ArrayList<Study> studies = new ArrayList<GeneScorePipeline.Study>();
  private final HashMap<String, Constraint> analysisConstraints = new HashMap<String, GeneScorePipeline.Constraint>();
  private final HashMap<String, HashMap<String, Integer>> dataCounts = new HashMap<String, HashMap<String, Integer>>();

  // private int bimChrIndex = 0;
  // private int bimMkrIndex = 1;
  // private int bimPosIndex = 3;
  // private int bimA1Index = 4;
  // private int bimA2Index = 5;

  private final int hitsMkrIndex = 1;

  private final Logger log;

  private static class HitMarker {

    private final String effectAllele;
    private final String nonEffectAllele;
    private final Double effect;

    /**
     * @param effectAllele
     * @param nonEffectAllele
     * @param effect
     */
    private HitMarker(String effectAllele, String nonEffectAllele, Double effect) {
      super();
      this.effectAllele = effectAllele;
      this.nonEffectAllele = nonEffectAllele;
      this.effect = effect;
    }

    public String getEffectAllele() {
      return effectAllele;
    }

    public String getNonEffectAllele() {
      return nonEffectAllele;
    }

    public Double getEffect() {
      return effect;
    }

  }

  private class Study {

    String studyDir;
    String studyName;

    String dataSource;
    ArrayList<DataSource> dataSources;
    HashMap<String, DosageData> data = new HashMap<String, DosageData>();

    ArrayList<String> phenoFiles = new ArrayList<String>();

    HashMap<String, PhenoData> phenoData = new HashMap<String, GeneScorePipeline.PhenoData>();

    // constraint -> datafile -> phenofile
    HashMap<String, HashMap<String, HashMap<String, RegressionResult>>> regressions = new HashMap<String, HashMap<String, HashMap<String, RegressionResult>>>();
    // constraint -> datafile
    HashMap<String, HashMap<String, double[]>> scores = new HashMap<String, HashMap<String, double[]>>();
    // constraint -> datafile
    HashMap<String, HashMap<String, Integer>> hitSnpCounts = new HashMap<String, HashMap<String, Integer>>();
    // constraint -> datafile
    HashMap<String, HashMap<String, Integer>> hitWindowCnts = new HashMap<String, HashMap<String, Integer>>();
    // constraint -> datafile
    HashMap<String, HashMap<String, Integer>> dataCounts = new HashMap<String, HashMap<String, Integer>>();
    Map<String, Map<String, HitMarker>> markerData = new HashMap<>();

    /**
     * Returns a hashSet containing all markers present in the data files that are present in
     * hitMkrSet.
     *
     * @param hitMkrSet
     * @return
     */
    public HashSet<String> retrieveMarkers(String dataKey, HashSet<String> hitMkrSet) {
      HashSet<String> returnMarkers = new HashSet<String>();
      if (data.get(dataKey).isEmpty()) {
        return returnMarkers;
      }
      String[] mkrs = data.get(dataKey).getMarkerSet().getMarkerNames();
      for (String mkr : mkrs) {
        if (hitMkrSet.contains(mkr)) {
          returnMarkers.add(mkr);
        }
      }
      return returnMarkers;
    }

    public boolean isPlinkData() {
      if (dataSources != null) {
        for (DataSource ds : dataSources) {
          if (!(ds.dataFile.endsWith(".bed") || ds.dataFile.endsWith(".ped"))) {
            return false;
          }
        }
        return true;
      }
      return false;
    }

    public void loadDataSources(String dataKey, String[] hitMkrs) {
      SnpMarkerSet markerSet = new SnpMarkerSet(hitMkrs, false, log);
      markerSet.parseSNPlocations(log);
      int[][] markerLocations = markerSet.getChrAndPositionsAsInts();
      dataSources = MergeExtractPipeline.parseDataFile(null, markerLocations, null, dataSource, 0,
                                                       log);
      if (dataSources.size() == 0) {
        // error
        log.reportError("Error - no data sources loaded from file: " + dataSource);
      } else {
        log.reportTime("Loading data file " + dataSources.get(0).dataFile);
        DosageData d0 = new DosageData(dataSources.get(0).dataFile, dataSources.get(0).idFile,
                                       dataSources.get(0).mapFile, null, hitMkrs, true, log);
        if (dataSources.size() > 1) {
          for (int i = 1; i < dataSources.size(); i++) {
            log.reportTime("Loading data file " + dataSources.get(i).dataFile);
            DosageData d1 = new DosageData(dataSources.get(i).dataFile, dataSources.get(i).idFile,
                                           dataSources.get(i).mapFile, null, hitMkrs, true, log);
            d0 = DosageData.combine(d0, d1, DosageData.COMBINE_OP.EITHER_IF_OTHER_MISSING, false, 0,
                                    log);
            System.gc();
          }
        }
        if (d0.isEmpty()) {
          log.reportError("no data for key: " + dataKey);
        }
        data.put(dataKey, d0);
        System.gc();
      }
    }
  }

  private class RegressionResult {

    boolean logistic;
    double rsq;
    double baseRSq;
    double pval;
    double beta;
    double se;
    int num;
    public double stats;

    private void dummy() {
      logistic = true;
      rsq = Double.NaN;
      baseRSq = Double.NaN;
      pval = Double.NaN;
      beta = Double.NaN;
      se = Double.NaN;
      num = 0;
      stats = Double.NaN;
    }

  }

  private class Constraint {

    final float indexThreshold;
    final int windowMinSizePerSide;
    final float windowExtensionThreshold;

    public Constraint(float i, int m, float w) {
      indexThreshold = i;
      windowMinSizePerSide = m;
      windowExtensionThreshold = w;
    }
  }

  private class PhenoData {

    String phenoName;
    HashMap<String, PhenoIndiv> indivs = new HashMap<String, GeneScorePipeline.PhenoIndiv>();
    ArrayList<String> covars = new ArrayList<String>();
  }

  private class PhenoIndiv {

    String fid, iid;
    double depvar;
    HashMap<String, Double> covars = new HashMap<String, Double>();
  }

  private static String[] readPlinkFile(int chr, String file, HashMap<String, String> rsToFull) {
    String[] rsPos = new String[Positions.CHROMOSOME_LENGTHS_MAX[chr]];
    try {
      BufferedReader reader = Files.getAppropriateReader(file);
      String temp = reader.readLine();
      String[] line;
      while ((temp = reader.readLine()) != null) {
        line = temp.split(PSF.Regex.GREEDY_WHITESPACE);
        String[] parts = line[1].split(":");
        String rsOrChr = parts[0];
        String posStr = parts[1];
        // TODO incorporate allele examination
        // String a1 = parts[2];
        // String a2 = parts[3];
        String freqStr = line[4];
        int pos = Integer.parseInt(posStr);
        rsPos[pos] = (rsOrChr.startsWith("rs")) ? rsOrChr : line[1];
        rsToFull.put((rsOrChr.startsWith("rs")) ? rsOrChr : line[1], line[1] + "\t" + freqStr);
      }
    } catch (NumberFormatException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    return rsPos;
  }

  private static HashMap<String, Double> get1000GFreq(HashMap<String, int[]> markerMap) {
    // BufferedReader reader = null;
    String currFile, /* temp, */plinkFile;
    // String[] line;
    HashMap<Integer, ArrayList<String>> mkrsByChr, nonRSByChr;
    ArrayList<String> chrMkrs, nonRS;
    HashMap<String, Double> plinkFreqs = new HashMap<String, Double>();

    mkrsByChr = new HashMap<Integer, ArrayList<String>>();
    nonRSByChr = new HashMap<Integer, ArrayList<String>>();
    for (java.util.Map.Entry<String, int[]> entry : markerMap.entrySet()) {
      chrMkrs = mkrsByChr.get(entry.getValue()[0]);
      nonRS = nonRSByChr.get(entry.getValue()[0]);
      if (chrMkrs == null) {
        chrMkrs = new ArrayList<String>();
        mkrsByChr.put(entry.getValue()[0], chrMkrs);
      }
      if (nonRS == null) {
        nonRS = new ArrayList<String>();
        nonRSByChr.put(entry.getValue()[0], nonRS);
      }
      chrMkrs.add(entry.getKey());
      if (!entry.getKey().startsWith("rs")) {
        nonRS.add(entry.getKey());
      }
    }

    plinkFile = PLINK_FRQ_DIR + PLINK_FRQ_FILE_PATTERN;
    for (Integer chr : mkrsByChr.keySet()) {
      // skip this chromosome if we don't have any markers for which we need freqencies
      if (mkrsByChr.get(chr) == null || mkrsByChr.get(chr).isEmpty()) {
        continue;
      }
      currFile = plinkFile.replace(TAG, chr.toString());

      HashMap<String, String> rsToFullData = new HashMap<String, String>();
      String[] posRS = readPlinkFile(chr, currFile, rsToFullData);

      chrMkrs = mkrsByChr.get(chr);
      nonRS = nonRSByChr.get(chr);

      ArrayList<String> rsNotFound = new ArrayList<String>();

      for (String mkr : chrMkrs) {
        if (rsToFullData.containsKey(mkr)) {
          plinkFreqs.put(mkr, Double.valueOf(rsToFullData.get(mkr).split("\t")[1]));
        } else {
          rsNotFound.add(mkr);
          continue;
        }
      }

      nonRS.addAll(rsNotFound);

      for (String nonRSMkr : nonRS) {
        int pos = markerMap.get(nonRSMkr)[1];
        String mkrNm = pos > 0 && pos < posRS.length && null != posRS[pos]
                       && !"".equals(posRS[pos]) ? posRS[pos] : "";
        plinkFreqs.put(nonRSMkr,
                       "".equals(mkrNm) ? 0.0
                                        : Double.valueOf(rsToFullData.get(mkrNm).split("\t")[1]));
      }

      // Hashtable<String, Vector<String>> plinkData = HashVec.loadFileToHashVec(currFile, 1, new
      // int[]{4}, "\t", true, false);
      //

      //
      // ArrayList<String> found = new ArrayList<String>();
      //
      // for (String mkr : chrMkrs) {
      // if (plinkData.containsKey(mkr))
      // }
      //
      // try {
      // reader = Files.getAppropriateReader(currFile);
      // temp = reader.readLine();
      // while((temp = reader.readLine()) != null) {
      // line = temp.split(PSF.Regex.GREEDY_WHITESPACE);
      // String rsOrPos = line[1];
      // String a1 = line[2];
      // String a2 = line[3];
      // String freqStr = line[4];
      //
      // if (rsOrPos.startsWith("rs")) {
      // if (mkrsByChr.get(chr).contains(rsOrPos.split(":")[0])) {
      // plinkFreqs.put(rsOrPos.split(":")[0], Double.valueOf(freqStr));
      // mkrsByChr.get(chr).remove(rsOrPos.split(":")[0]);
      // }
      // } else if (!nonRSByChr.get(chr).isEmpty()){
      // // check positions
      // int plinkPos = Integer.parseInt(rsOrPos.split(":")[1]);
      // for (String nonRSMkr : nonRSByChr.get(chr)) {
      // int mkrPos = markerMap.get(nonRSMkr)[1];
      // if (Math.abs(plinkPos - mkrPos) <= ACCEPTABLE_SEPARATION) {
      // plinkFreqs.put(nonRSMkr, Double.valueOf(freqStr));
      // }
      // }
      // }
      // }
      // reader.close();
      // reader = null;
      // } catch (IOException e) {
      // // TODO Auto-generated catch block
      // e.printStackTrace();
      // if (reader != null) {
      // try {
      // reader.close();
      // } catch (IOException e1) {
      // // TODO Auto-generated catch block
      // // e1.printStackTrace();
      // }
      // reader = null;
      // }
      // }
    }

    return plinkFreqs;
  }

  public static void preprocessDataFiles(String[] files, Logger log) {
    BufferedReader reader;
    String temp, delimiter;
    String[] header, snps = null, data, line;
    String[][] factors;
    int[] indices;
    HashMap<String, int[]> markerMap;
    Hashtable<String, Vector<String>> fileData;
    HashMap<String, Double> freqs;

    factors = new String[][] {Aliases.MARKER_NAMES, Aliases.CHRS, Aliases.POSITIONS,
                              Aliases.PVALUES, Aliases.ALLELE_FREQS, Aliases.EFFECTS};
    for (String filename : files) {
      try {
        reader = Files.getAppropriateReader(filename);
        temp = reader.readLine();
        delimiter = ext.determineDelimiter(temp);
        header = temp.trim().split(delimiter);
        indices = ext.indexFactors(factors, header, false, false, true, true, log);
        markerMap = new HashMap<String, int[]>();
        String errorMsg = "";
        if (indices[0] == -1) {
          errorMsg = "ERROR - no MarkerName column found";
          // ERROR - couldn't find MarkerName column! COMPLETE FAIL
        }
        if (indices[5] == -1) {
          // NO BETAS! COMPLETE FAIL
          errorMsg = errorMsg.equals("") ? "ERROR - no Beta/Effect column found"
                                         : errorMsg + "; no Beta/Effect column found";
        }
        if (indices[3] == -1) {
          // NO PVALUES! COMPLETE FAIL
          errorMsg = errorMsg.equals("") ? "ERROR - no P-Value column found"
                                         : errorMsg + "; no P-Value column found";
        }
        if (errorMsg.equals("")) {
          if (indices[1] == -1 || indices[2] == -1) {
            // No chromosomes/positions
            snps = HashVec.loadFileToStringArray(filename, true, new int[] {indices[0]}, false);// fileData.keySet().toArray(new
            // String[fileData.size()]);
            Files.writeArray(snps, ext.rootOf(filename, false) + ".snps");
            MapSNPsAndGenes.procSNPsToGenes(ext.parseDirectoryOfFile(filename),
                                            ext.rootOf(filename, true) + ".snps",
                                            MapSNPsAndGenes.DEFAULT_WIGGLE_ROOM, (byte) 37, log,
                                            true, false, false, null, null, null, false); // TODO
                                                                                                                                                                                                                                                        // should
                                                                                                                                                                                                                                                        // run
                                                                                                                                                                                                                                                        // SnpEff
                                                                                                                                                                                                                                                        // too?
            data = ArrayUtils.toStringArray(HashVec.loadFileToVec(ext.rootOf(filename, false)
                                                                  + "_positions.xln", true, false,
                                                                  false));
            for (String element : data) {
              line = element.trim().split(PSF.Regex.GREEDY_WHITESPACE);
              markerMap.put(line[0],
                            new int[] {Positions.chromosomeNumber(line[1]),
                                       ext.isMissingValue(line[2]) ? -1
                                                                   : Integer.parseInt(line[2])});
            }
          } else {
            // fileData = HashVec.loadFileToHashVec(filename, indices[0], new int[]{indices[1],
            // indices[2]}, "\t", true, false);
            // for (String key : fileData.keySet()) {
            // markerMap.put(key, new int[]{Positions.chromosomeNumber(fileData.get(key).get(0)),
            // Integer.parseInt(fileData.get(key).get(1))});
            // }
            // snps = fileData.keySet().toArray(new String[fileData.size()]);
            // fileData = null;
          }
          if (indices[4] == -1) {
            // no frequencies
            if (markerMap.isEmpty()) {
              fileData = HashVec.loadFileToHashVec(filename, indices[0],
                                                   new int[] {indices[1], indices[2]}, "\t", true,
                                                   false);
              for (String key : fileData.keySet()) {
                markerMap.put(key,
                              new int[] {Positions.chromosomeNumber(fileData.get(key).get(0)
                                                                            .split("\t")[0]),
                                         new BigDecimal(fileData.get(key).get(0)
                                                                .split("\t")[1]).intValueExact()});
              }
              fileData = null;
            }
            freqs = get1000GFreq(markerMap);
          } else {
            freqs = null;// new HashMap<String, Double>();
            // while ((temp = reader.readLine()) != null) {
            // line = temp.trim().split(PSF.Regex.GREEDY_WHITESPACE);
            // freqs.put(line[indices[0]], Double.valueOf(line[indices[4]]));
            // }
          }

          StringBuilder newHeaderSB = new StringBuilder("SNP\tChr\tPos\tFreq\tP\tBeta");
          for (int i = 0; i < header.length; i++) {
            if (i != indices[0] && i != indices[1] && i != indices[2] && i != indices[3]
                && i != indices[4] && i != indices[5]) {
              newHeaderSB.append("\t").append(header[i]);
            }
          }

          PrintWriter metaWriter = Files.getAppropriateWriter(ext.rootOf(filename, false)
                                                              + ".meta");
          metaWriter.println(newHeaderSB.toString());
          while ((temp = reader.readLine()) != null) {
            line = temp.trim().split(PSF.Regex.GREEDY_WHITESPACE);
            String snp = line[indices[0]];
            String chr = indices[1] == -1 ? "" + markerMap.get(snp)[0] : line[indices[1]];
            String pos = indices[2] == -1 ? "" + markerMap.get(snp)[1] : line[indices[2]];
            String pval = line[indices[3]];
            String freq = indices[4] == -1 ? "" + (freqs == null || freqs.isEmpty()
                                                   || freqs.get(snp) == null ? 0.0 : freqs.get(snp))
                                           : line[indices[4]];
            String beta = line[indices[5]];
            StringBuilder writeLineSB = new StringBuilder();
            writeLineSB.append(snp).append("\t").append(chr).append("\t").append(pos).append("\t")
                       .append(freq).append("\t").append(pval).append("\t").append(beta);
            for (int i = 0; i < line.length; i++) {
              if (i != indices[0] && i != indices[1] && i != indices[2] && i != indices[3]
                  && i != indices[4] && i != indices[5]) {
                writeLineSB.append("\t").append(line[i]);
              }
            }
            metaWriter.println(writeLineSB.toString());
          }
          metaWriter.flush();
          metaWriter.close();
        } else {
          log.reportError(errorMsg);
          reader.close();
          continue;
        }
      } catch (IOException e) {
        e.printStackTrace();
      }
    }
  }

  public GeneScorePipeline(String metaDir, float[] indexThresholds, int[] windowMins,
                           float[] windowExtThresholds, double missThresh, Logger log) {
    this.log = log;
    this.metaDir = metaDir;
    // this.numThreads = numThreads;
    // this.runPlink = plink;
    // this.runRegression = runPlink && regression;
    // this.writeHist = runPlink && histogram;
    this.minMissThresh = missThresh;
    this.indexThresholds = indexThresholds;
    windowMinSizePerSides = windowMins;
    windowExtensionThresholds = windowExtThresholds;
    setFilePrefices();
    loadStudyFolders();
    // instantiate inner hashmaps:
    for (Study study : studies) {
      for (String pref : analysisConstraints.keySet()) {
        HashMap<String, HashMap<String, RegressionResult>> res = new HashMap<String, HashMap<String, RegressionResult>>();
        for (String dFile : metaFiles) {
          String dataFile = ext.rootOf(dFile, false);
          HashMap<String, RegressionResult> res2 = new HashMap<String, RegressionResult>();
          res.put(dataFile, res2);
        }
        study.regressions.put(pref, res);

        HashMap<String, Integer> cntMap = new HashMap<String, Integer>();
        HashMap<String, Integer> hitMap = new HashMap<String, Integer>();
        HashMap<String, Integer> cntMap2 = new HashMap<String, Integer>();
        study.hitWindowCnts.put(pref, cntMap);
        study.hitSnpCounts.put(pref, hitMap);
        study.dataCounts.put(pref, cntMap2);
      }
    }
    loadDataCounts();
    runMetaHitWindowsAndLoadData();
  }

  private void setFilePrefices() {
    for (float i : indexThresholds) {
      for (int m : windowMinSizePerSides) {
        for (float w : windowExtensionThresholds) {
          StringBuilder prefixSB = new StringBuilder();
          prefixSB.append(ext.formSciNot(i, 4, false)).append("_")
                  .append(ext.formSciNot(m, 4, false)).append("_")
                  .append(ext.formSciNot(w, 4, false));
          analysisConstraints.put(prefixSB.toString(), new Constraint(i, m, w));
        }
      }
    }
  }

  private void loadStudyFolders() {
    File dir = new File(metaDir);
    File[] fs = dir.listFiles();
    for (File f : fs) {
      if (f.isDirectory()) {
        Study study = new Study();
        study.studyName = ext.rootOf(f.getAbsolutePath(), true);
        study.studyDir = f.getAbsolutePath() + "/";
        for (File f1 : f.listFiles()) {
          if (f1.getName().endsWith(".pheno")) {
            study.phenoFiles.add(f1.getName());
          }
          if (f1.getName().equals(DATA_SOURCE_FILENAME)) {
            study.dataSource = f1.getAbsolutePath();
          }
        }
        if (study.dataSource == null) {
          log.reportError("Error - data source file {" + DATA_SOURCE_FILENAME
                          + "} missing for study " + study.studyName);
        }
        studies.add(study);
      } else if (f.getAbsolutePath().endsWith(".meta")) {
        metaFiles.add(f.getName());
      }
    }
  }

  private void loadDataCounts() {
    String countsFile = metaDir + "data.cnt";

    if ((new File(countsFile).exists())) {
      try {
        BufferedReader reader = Files.getAppropriateReader(countsFile);
        String line = null;
        while ((line = reader.readLine()) != null && !"".equals(line)) {
          String[] temp = line.split("\t");
          HashMap<String, Integer> dFileCnts = dataCounts.get(temp[0]);
          if (dFileCnts == null) {
            dFileCnts = new HashMap<String, Integer>();
            dataCounts.put(temp[0], dFileCnts);
          }
          dFileCnts.put(temp[1], Integer.parseInt(temp[2]));
        }
      } catch (NumberFormatException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      } catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
    }

    for (String dFile : metaFiles) {
      HashMap<String, Constraint> threshNeed = new HashMap<String, GeneScorePipeline.Constraint>();

      HashMap<String, Integer> cnts = dataCounts.get(dFile);
      if (cnts == null) {
        threshNeed.putAll(analysisConstraints);
        cnts = new HashMap<String, Integer>();
        dataCounts.put(dFile, cnts);
      } else {
        for (String pref : analysisConstraints.keySet()) {
          if (!cnts.containsKey(pref)) {
            threshNeed.put(pref, analysisConstraints.get(pref));
          }
        }
      }

      if (threshNeed.size() > 0) {
        try {
          BufferedReader reader = Files.getAppropriateReader(metaDir + dFile);
          String line = reader.readLine();
          String[] dataHdrs = line.split(PSF.Regex.GREEDY_WHITESPACE);
          int[] indices = ext.indexFactors(Aliases.PVALUES, dataHdrs, false);
          int ind = -1;
          for (int i : indices) {
            if (i > 0) {
              ind = i;
              break;
            }
          }
          if (ind > 0) {
            while ((line = reader.readLine()) != null) {
              String[] parts = line.split("\t");
              if (!ext.isMissingValue(parts[ind]) && ext.isValidDouble(parts[ind])) {
                double pval = Double.parseDouble(parts[ind]);
                for (java.util.Map.Entry<String, Constraint> constraint : threshNeed.entrySet()) {
                  if (pval < constraint.getValue().indexThreshold) {
                    Integer cnt = dataCounts.get(dFile).get(constraint.getKey());
                    if (cnt == null) {
                      cnt = 0;
                    }
                    cnt = cnt + 1;
                    dataCounts.get(dFile).put(constraint.getKey(), cnt);
                  }
                }
              }
            }
          }
        } catch (NumberFormatException e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
        } catch (IOException e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
      }
    }

    StringBuilder output = new StringBuilder();
    for (java.util.Map.Entry<String, HashMap<String, Integer>> entry : dataCounts.entrySet()) {
      for (java.util.Map.Entry<String, Integer> subEntry : entry.getValue().entrySet()) {
        output.append(entry.getKey()).append("\t").append(subEntry.getKey()).append("\t")
              .append(subEntry.getValue()).append("\n");
      }
    }
    Files.write(output.toString(), countsFile);

  }

  private void runMetaHitWindowsAndLoadData() {
    String[][] factors = new String[][] {Aliases.MARKER_NAMES, Aliases.EFFECTS,
                                         Aliases.ALLELE_FREQS, Aliases.PVALUES};

    for (String dFile : metaFiles) {
      String dataFile = ext.rootOf(dFile, false);
      HashSet<String> hitMkrSet = new HashSet<String>();
      for (java.util.Map.Entry<String, Constraint> filePrefix : analysisConstraints.entrySet()) {

        int metaCount = Files.countLines(metaDir + dFile, 0);

        if (metaCount > 1000) {
          String[][] results = HitWindows.determine(metaDir + dFile,
                                                    filePrefix.getValue().indexThreshold,
                                                    filePrefix.getValue().windowMinSizePerSide,
                                                    filePrefix.getValue().windowExtensionThreshold,
                                                    DEFAULT_ADDL_ANNOT_VAR_NAMES, new Logger());
          if (results == null) {
            log.reportError("HitWindows result was null for " + dFile + ". Using all SNPs");
          } else {
            log.report(ext.getTime() + "]\tFound " + results.length + " hit windows");
            for (String[] result : results) {
              hitMkrSet.add(result[1]);
            }
          }
        }
        if (hitMkrSet.isEmpty()) {
          String[] mkrs = HashVec.loadFileToStringArray(metaDir + dFile, true, new int[] {0},
                                                        false);
          if (mkrs == null) {
            log.reportError(".meta file was empty for " + dFile);
          } else {
            log.report(ext.getTime() + "]\tUsing all " + mkrs.length + " SNPs in .meta file");
            for (String mkr : mkrs) {
              hitMkrSet.add(mkr);
            }
          }
        }
        // uncomment to use all markers in dataFile

        if (!hitMkrSet.isEmpty()) {
          // read betas and freqs for hitwindow markers
          HashMap<String, double[]> dataMarkers = new HashMap<String, double[]>();
          try {
            BufferedReader reader = Files.getAppropriateReader(metaDir + dFile);
            String line = reader.readLine();
            String[] temp = line.split(PSF.Regex.GREEDY_WHITESPACE);
            int[] indices = ext.indexFactors(factors, temp, false, false, true, true, new Logger());
            while ((line = reader.readLine()) != null) {
              String mkr = line.split(PSF.Regex.GREEDY_WHITESPACE)[indices[0]];
              if (hitMkrSet.contains(mkr)) {
                if ((indices[1] != -1
                     && ext.isMissingValue(line.split(PSF.Regex.GREEDY_WHITESPACE)[indices[1]]))
                    || ext.isMissingValue(line.split(PSF.Regex.GREEDY_WHITESPACE)[indices[2]])
                    || ext.isMissingValue(line.split(PSF.Regex.GREEDY_WHITESPACE)[indices[3]])) {
                  hitMkrSet.remove(mkr);
                  continue;
                }
                dataMarkers.put(mkr,
                                new double[] {indices[1] == -1 ? Double.NaN
                                                               : Double.parseDouble(line.split(PSF.Regex.GREEDY_WHITESPACE)[indices[1]]),
                                              Double.parseDouble(line.split(PSF.Regex.GREEDY_WHITESPACE)[indices[2]]),
                                              Double.parseDouble(line.split(PSF.Regex.GREEDY_WHITESPACE)[indices[3]])});
              }
            }
            reader.close();

            double dataScore1 = getBetaFreqScore(dataMarkers);
            double dataScore2 = getChiDistRevScore(dataMarkers);

            String[] hitMkrs = hitMkrSet.toArray(new String[hitMkrSet.size()]);
            // cross-ref PLINK markers
            for (Study study : studies) {
              study.loadDataSources(dataFile + "\t" + filePrefix.getKey(), hitMkrs);

              HashMap<String, double[]> bimSubsetMarkers = new HashMap<String, double[]>();
              HashSet<String> bimMkrSet = study.retrieveMarkers(dataFile + "\t"
                                                                + filePrefix.getKey(), hitMkrSet);

              // pull betas and freqs for union markers
              for (String mkr : bimMkrSet) {
                if (hitMkrSet.contains(mkr)) {
                  bimSubsetMarkers.put(mkr, dataMarkers.get(mkr));
                }
              }

              // apply equation and set value for overall results
              double bimScore1 = getBetaFreqScore(bimSubsetMarkers);
              double bimScore2 = getChiDistRevScore(bimSubsetMarkers);

              HashMap<String, double[]> fileMap = study.scores.get(filePrefix.getKey());
              if (fileMap == null) {
                fileMap = new HashMap<String, double[]>();
                study.scores.put(filePrefix.getKey(), fileMap);
              }
              fileMap.put(dataFile, new double[] {bimScore1 / dataScore1, bimScore2 / dataScore2});
            }
          } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
          }
        }
      }
    }
  }

  private void createAffectedPhenoFiles(Study study) {
    // String famFile = study.studyDir + "plink.fam";
    String affFile = study.studyDir + "AFFECTED.pheno";

    // skip if we don't have PHENO data
    if (!study.isPlinkData()) {
      return;
    }
    // if affected.pheno file already exists, skip
    if (study.phenoFiles.contains("AFFECTED.pheno")) {
      return;
    }
    // if any of the data folders have been created, skip creating affected.pheno file
    for (String dFile : metaFiles) {
      String dataFile = ext.rootOf(dFile, false);
      if ((new File(study.studyDir + dataFile + "/")).exists()) {
        return;
      }
    }

    BufferedReader reader;
    // 0 - fid
    // 1 - iid
    // 2
    // 3
    // 4 - sex
    // 5 - pheno
    ArrayList<String> fam = new ArrayList<String>();
    ArrayList<String> pheno = new ArrayList<String>();
    String temp;
    try {
      for (int i = 0; i < study.dataSources.size(); i++) {
        reader = Files.getAppropriateReader(study.dataSources.get(i).idFile); // only satisfies
                                                                             // 'isPlinkData()'
                                                                             // if only one data
                                                                             // source is present
        while ((temp = reader.readLine()) != null) {
          String[] line = temp.split(PSF.Regex.GREEDY_WHITESPACE);
          String affLine = line[0] + "\t" + line[1] + "\t"
                           + (ext.isMissingValue(line[5]) ? "."
                                                          : -1 * (Integer.parseInt(line[5]) - 2))
                           + "\t" + line[4];
          if (!ext.isMissingValue(line[5])) {
            if (ext.isValidDouble(line[5]) && Double.parseDouble(line[5]) != 0.0) {
              fam.add(affLine);
            }
            pheno.add(line[5]);
          }
        }

        String[] unique = ArrayUtils.unique(pheno.toArray(new String[] {}));
        if (unique.length == 1) {
          log.report("Error - no variance in pheno data from .fam file for study '"
                     + study.studyName + "'");
          continue;
        }
      }
    } catch (NumberFormatException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

    String header = "FID\tIID\tPHENO\tMALE";
    PrintWriter writer = Files.getAppropriateWriter(affFile);
    writer.println(header);
    for (String line : fam) {
      writer.println(line);
    }
    writer.flush();
    writer.close();

    study.phenoFiles.add("AFFECTED.pheno");
  }

  private void loadPhenoFiles(Study study) {
    for (String pheno : study.phenoFiles) {
      PhenoData pd = new PhenoData();
      pd.phenoName = pheno;

      try {
        BufferedReader reader = Files.getAppropriateReader(study.studyDir + pheno);
        String[] header = reader.readLine().split("\t");
        // fid == header[0]
        // iid == header[1]
        // depvar = header[2];
        ArrayList<String> covars = new ArrayList<String>();
        for (int i = 3; i < header.length; i++) {
          covars.add(header[i]);
        }
        pd.covars.addAll(covars);
        String temp = reader.readLine();
        indiv: do {
          String[] line = temp.split("\t");

          if (!ext.isMissingValue(line[2])) {
            PhenoIndiv pi = new PhenoIndiv();
            pi.fid = line[0];
            pi.iid = line[1];
            pi.depvar = Double.parseDouble(line[2]);
            for (int i = 3; i < line.length; i++) {
              if (ext.isMissingValue(line[i])) {
                continue indiv;
              }
              pi.covars.put(header[i], Double.parseDouble(line[i]));
            }
            pd.indivs.put(pi.fid + "\t" + pi.iid, pi);
          }

        } while ((temp = reader.readLine()) != null);

        reader.close();
      } catch (IOException e) {
        e.printStackTrace();
      }

      study.phenoData.put(pheno, pd);
    }
  }

  private double getBetaFreqScore(HashMap<String, double[]> markerMap) {
    double sum = 0.0;
    // (Beta^2 * 2 * MAF (1-MAF))
    for (java.util.Map.Entry<String, double[]> entry : markerMap.entrySet()) {
      double score = 2 * (entry.getValue()[0] * entry.getValue()[0]) * entry.getValue()[1]
                     * (1 - entry.getValue()[1]);
      sum += score;
    }
    return sum;
  }

  private double getChiDistRevScore(HashMap<String, double[]> markerMap) {
    double sum = 0.0;
    for (java.util.Map.Entry<String, double[]> entry : markerMap.entrySet()) {
      sum += ProbDist.ChiDistReverse(entry.getValue()[2], 1);
    }
    return sum;
  }

  public void runPipeline() {
    log.report(ext.getTime() + "]\tProcessing study data [" + studies.size() + " total]:");
    // if (numThreads == 1) {
    // for (String studyDir : studyFolders) {
    // processStudy(studyDir);
    // }
    for (Study study : studies) {
      createAffectedPhenoFiles(study);
      loadPhenoFiles(study);
      processStudy(study);
    }
    writeResults();
    // } else {
    // ExecutorService server = Executors.newFixedThreadPool(numThreads);
    //
    // }
    log.report(ext.getTime() + "]\tProcessing Complete!");
  }

  private void processStudy(Study study) {
    try {
      createFolders(study);
      crossFilterMarkerData(study);
      runHitWindows(study);
      extractHitMarkerData(study);
      runScore(study);
      runRegression(study);
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  private void createFolders(Study study) {
    for (String dataFile : metaFiles) {
      String dataFolder = study.studyDir + ext.rootOf(dataFile, true) + "/";
      for (String constraints : analysisConstraints.keySet()) {
        String constraintFolder = dataFolder + constraints + "/";
        File f = new File(constraintFolder);
        if (!(f.exists())) {
          f.mkdirs();
        }
      }
    }
  }

  private void crossFilterMarkerData(Study study) throws IOException {
    for (String dFile : metaFiles) {
      String dataFile = ext.rootOf(dFile, false);
      for (java.util.Map.Entry<String, Constraint> constraintEntry : analysisConstraints.entrySet()) {
        String crossFilterFile = study.studyDir + dataFile + "/" + constraintEntry.getKey() + "/"
                                 + CROSS_FILTERED_DATAFILE;
        if ((new File(crossFilterFile).exists())) {
          log.report(ext.getTime() + "]\tCross-filtered data file already exists! [ --> '"
                     + crossFilterFile + "']");
          study.hitSnpCounts.get(constraintEntry.getKey())
                            .put(dataFile, Files.countLines(crossFilterFile, 1));
          continue;
        }
        if (study.data.get(dataFile + "\t" + constraintEntry.getKey()).isEmpty()) {
          study.hitSnpCounts.get(constraintEntry.getKey()).put(dataFile, 0);
          continue;
        }

        log.report(ext.getTime() + "]\tCross-filtering data and .BIM files [ --> '"
                   + crossFilterFile + "']");
        BufferedReader dataReader;
        PrintWriter dataWriter;
        HashMap<String, int[]> mkrsBim;

        mkrsBim = new HashMap<String, int[]>();
        SnpMarkerSet markerSet = study.data.get(dataFile + "\t" + constraintEntry.getKey())
                                           .getMarkerSet();
        int cntAmbig = 0;

        String[] mkrNames = markerSet.getMarkerNames();
        String[][] alleles = markerSet.getAlleles();
        int[][] chrPos = markerSet.getChrAndPositionsAsInts();

        for (int i = 0; i < mkrNames.length; i++) {
          String a1 = (alleles[i][0]).toUpperCase();
          String a2 = (alleles[i][1]).toUpperCase();
          if (Sequence.validAllele(a1) && Sequence.validAllele(a2) // TODO validAllele requires a
                                                                  // 1 character allele
              && !a1.equals(Sequence.flip(a2))) {
            mkrsBim.put(mkrNames[i], chrPos[i]);
          } else {
            cntAmbig++;
          }
        }

        log.report(ext.getTime() + "]\tFound " + cntAmbig
                   + " ambiguous markers (will be excluded)");
        dataReader = Files.getAppropriateReader(metaDir + dFile);
        dataWriter = new PrintWriter(crossFilterFile);

        int cnt = 0;
        String dataHdr = dataReader.readLine();
        String[] dataHdrs = dataHdr.split(PSF.Regex.GREEDY_WHITESPACE);
        int[] indices = ext.indexFactors(Aliases.PVALUES, dataHdrs, false);
        int ind = -1;
        for (int i : indices) {
          if (i > 0) {
            ind = i;
            break;
          }
        }
        String line;
        dataWriter.println("MarkerName\tChr\tPosition\t"
                           + ArrayUtils.toStr(ArrayUtils.subArray(dataHdrs, 1))); // Allele1\tAllele2\tFreq.Allele1.HapMapCEU\tb\tSE\tp\tN
        while ((line = dataReader.readLine()) != null) {
          String[] parts = line.split(PSF.Regex.GREEDY_WHITESPACE);
          if (mkrsBim.containsKey(parts[0])
              && (ind == -1
                  || (ext.isValidDouble(parts[ind])
                      && Double.parseDouble(parts[ind]) < constraintEntry.getValue().indexThreshold))) {
            int[] chrPosBim = mkrsBim.get(parts[0]);
            dataWriter.print(parts[0]);
            dataWriter.print("\t");
            dataWriter.print(chrPosBim[0]);
            dataWriter.print("\t");
            dataWriter.print(chrPosBim[1]);
            dataWriter.print("\t");
            dataWriter.println(ArrayUtils.toStr(ArrayUtils.subArray(parts, 1)));
            cnt++;
          }
        }
        dataWriter.flush();
        dataWriter.close();
        dataReader.close();

        study.hitSnpCounts.get(constraintEntry.getKey()).put(dataFile, cnt);
      }
    }
  }

  private void runHitWindows(Study study) {
    for (String dFile : metaFiles) {
      String dataFile = ext.rootOf(dFile, false);

      for (java.util.Map.Entry<String, Constraint> filePrefix : analysisConstraints.entrySet()) {
        File prefDir = new File(study.studyDir + dataFile + "/" + filePrefix.getKey() + "/");
        String crossFilterFile = prefDir + "/" + CROSS_FILTERED_DATAFILE;
        String hitsFile = prefDir + "/hits_" + filePrefix.getKey() + ".out";
        if ((new File(hitsFile)).exists()) {
          log.report(ext.getTime() + "]\tHit window analysis file already exists! [ --> '"
                     + hitsFile + "']");
          study.hitWindowCnts.get(filePrefix.getKey()).put(dataFile, Files.countLines(hitsFile, 1));
          continue;
        }
        if (study.data.get(dataFile + "\t" + filePrefix.getKey()).isEmpty()) {
          continue;
        }

        log.report(ext.getTime() + "]\tRunning hit window analysis [ --> '" + hitsFile + "']");
        String[][] results = HitWindows.determine(crossFilterFile,
                                                  filePrefix.getValue().indexThreshold,
                                                  filePrefix.getValue().windowMinSizePerSide,
                                                  filePrefix.getValue().windowExtensionThreshold,
                                                  DEFAULT_ADDL_ANNOT_VAR_NAMES, log);
        if (results == null) {
          log.reportError("Error - HitWindows result from " + crossFilterFile + " was null");
        } else {
          int hits = results.length - 1; // Don't count header
          log.report(ext.getTime() + "]\tFound " + hits + " hit windows");
          Files.writeMatrix(results, hitsFile, "\t");
          study.hitWindowCnts.get(filePrefix.getKey()).put(dataFile, hits);
        }
      }
    }
  }

  private void extractHitMarkerData(Study study) {
    for (String dFile : metaFiles) {
      String dataFile = ext.rootOf(dFile, false);

      for (java.util.Map.Entry<String, Constraint> filePrefix : analysisConstraints.entrySet()) {
        if (study.data.get(dataFile + "\t" + filePrefix.getKey()).isEmpty()) {
          continue;
        }

        File prefDir = new File(study.studyDir + dataFile + "/" + filePrefix.getKey() + "/");

        String crossFilterFile = prefDir + "/" + CROSS_FILTERED_DATAFILE;
        String hitsFile = prefDir + "/hits_" + filePrefix.getKey() + ".out";
        String[] hitMarkers = HashVec.loadFileToStringArray(hitsFile, true,
                                                            new int[] {hitsMkrIndex}, false);
        HashSet<String> hitMrkSet = new HashSet<String>();
        for (String mkr : hitMarkers) {
          hitMrkSet.add(mkr);
        }

        ColumnFilter hitMarkerFilter = new AbstractColumnFilter(MARKER_COL) {

          @Override
          public boolean filter(DataLine values) {
            return hitMrkSet.contains(values.get(MARKER_COL, null));
          }
        };
        try (FileParser crossFilterParser = FileParserFactory.setup(crossFilterFile, MARKER_COL)
                                                        .optionalColumns(EFFECT_ALLELE_COL,
                                                                         NON_EFFECT_ALLELE_COL,
                                                                         BETA_COL)
                                                             .filter(hitMarkerFilter).build()) {
          String subsetFile = prefDir + "/subsetData_" + filePrefix.getKey() + ".xln";
          Map<String, HitMarker> dataList = crossFilterParser.parseToFileAndLoad(subsetFile, "\t",
                                                                                 true, MARKER_COL)
                                                             .entrySet().stream()
                                                             .collect(Collectors.toMap(Map.Entry::getKey,
                                                                                       e -> formHitMarker(e.getValue())));
          study.markerData.put(dFile + "\t" + filePrefix.getKey(), dataList);
        } catch (IOException e) {
          log.reportIOException(crossFilterFile);
        }
      }
    }
  }

  private static final HitMarker formHitMarker(DataLine hitMarkerDataLine) {
    return new HitMarker(hitMarkerDataLine.get(EFFECT_ALLELE_COL, null),
                         hitMarkerDataLine.get(NON_EFFECT_ALLELE_COL, null),
                         hitMarkerDataLine.get(BETA_COL, null));
  }

  private static final String SCORE_FILE = "score.profile";

  private void runScore(Study study) {
    for (String dFile : metaFiles) {
      String dataFile = ext.rootOf(dFile, false);
      for (java.util.Map.Entry<String, Constraint> filePrefix : analysisConstraints.entrySet()) {
        File prefDir = new File(study.studyDir + dataFile + "/" + filePrefix.getKey() + "/");
        if (!prefDir.exists()) {
          log.report(ext.getTime() + "]\tError - no subfolder for '" + filePrefix.getKey()
                     + "' analysis");
          continue;
        }
        if ((new File(prefDir + "/" + SCORE_FILE)).exists()) {
          log.report(ext.getTime() + "]\tPlink analysis results file already exists! [ --> '"
                     + prefDir + "/" + SCORE_FILE + "']");
          continue;
        }
        DosageData data = study.data.get(dataFile + "\t" + filePrefix.getKey());
        if (data.isEmpty()) {
          continue;
        }

        Map<String, HitMarker> hitMarkerData = study.markerData.get(dFile + "\t"
                                                                    + filePrefix.getKey());
        PrintWriter scoreWriter;
        scoreWriter = Files.getAppropriateWriter(prefDir + "/" + SCORE_FILE);

        String[][] ids = data.getIds();
        float[][] dose = data.getDosageValues();
        if (dose == null) {
          data.computeDosageValues(log);
          dose = data.getDosageValues();
        }
        if (dose == null) {
          log.reportError("No dosage data available for {" + dataFile + "\t" + filePrefix.getKey()
                          + "}");
          scoreWriter.close();
          continue;
        }
        String[] markers = data.getMarkerSet().getMarkerNames();
        String[][] alleles = data.getMarkerSet().getAlleles();
        Map<String, Integer> matchedMarkerIndices = new HashMap<>();
        Map<String, Float> matchedMarkerFreqs = new HashMap<>();
        Map<String, AlleleOrder> matchedMarkerAlleleOrders = new HashMap<>();
        for (int m = 0; m < markers.length; m++) {
          String mkr = markers[m];
          HitMarker hitMarker = hitMarkerData.get(mkr);
          if (hitMarker == null) continue;
          AlleleOrder alleleOrder = determineAlleleOrder(alleles[m], hitMarker);
          if (!(alleleOrder.equals(AlleleOrder.SAME) || alleleOrder.equals(AlleleOrder.OPPOSITE))) {
            Joiner alleleJoiner = Joiner.on('/');
            log.reportError("Alleles in study (" + alleleJoiner.join(alleles[m])
                            + ") do not match source alleles ("
                            + alleleJoiner.join(hitMarker.getEffectAllele(),
                                                hitMarker.getNonEffectAllele())
                            + ") for " + mkr);
            continue;
          }
          matchedMarkerAlleleOrders.put(mkr, alleleOrder);
          matchedMarkerIndices.put(mkr, m);
          int cnt = 0;
          float tot = 0;
          for (int i = 0; i < ids.length; i++) {
            if (!Float.isNaN(dose[m][i])) {
              tot += dose[m][i];
              cnt++;
            }
          }
          matchedMarkerFreqs.put(mkr, tot / cnt);
        }
        for (int i = 0; i < ids.length; i++) {
          float scoreSum = 0;
          float cnt2 = 0;
          int cnt = 0;

          for (Map.Entry<String, Integer> mkrIndexEntry : matchedMarkerIndices.entrySet()) {
            String mkr = mkrIndexEntry.getKey();
            int mkrIndex = mkrIndexEntry.getValue();
            float dosage = dose[mkrIndex][i];
            boolean isNaN = Float.isNaN(dosage);
            HitMarker hitMarker = hitMarkerData.get(mkr);
            float beta = hitMarker.getEffect().floatValue();
            AlleleOrder alleleOrder = matchedMarkerAlleleOrders.get(mkr);
            if (alleleOrder.equals(StrandOps.AlleleOrder.SAME)) {
              cnt += isNaN ? 0 : 1;
              cnt2 += isNaN ? 0 : (2.0 - dosage);
              scoreSum += (2.0 - (isNaN ? matchedMarkerFreqs.get(mkr) : dosage)) * beta;
            } else if (alleleOrder.equals(StrandOps.AlleleOrder.OPPOSITE)) {
              cnt += isNaN ? 0 : 1;
              cnt2 += isNaN ? 0 : dosage;
              scoreSum += (isNaN ? matchedMarkerFreqs.get(mkr) : dosage) * beta;
            } else {
              throw new IllegalStateException("Mismatched alleles were not caught when mapping");
            }
            // int code = Metal.determineStrandConfig(new String[]{}, new String[]{});
            // if (code == Metal.STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND || code ==
            // Metal.STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND) {
            // cnt += isNaN ? 0 : 1;
            // cnt2 += isNaN ? 0 : (2.0 - dosage);
            // scoreSum += (2.0 - (isNaN ? freqs.get(mkr) : dosage)) * beta;
            // } else if (code == Metal.STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND || code ==
            // Metal.STRAND_CONFIG_SAME_ORDER_SAME_STRAND) {
            // cnt += isNaN ? 0 : 1;
            // cnt2 += isNaN ? 0 : dosage;
            // scoreSum += (isNaN ? freqs.get(mkr) : dosage) * beta;
            // } else {
            // System.err.println("Error - METAL CODE " + code + " - " + mkr + " == mkrAlleles: {" +
            // mkrAllele1 + ", " + mkrAllele2 + "}; scoringAlleles: {" + scoringAllele1 + ", " +
            // scoringAllele2 + "}");
            // }
          }

          double mkrRatio = markers.length / (double) cnt;
          if (mkrRatio > minMissThresh) {
            scoreWriter.println(ids[i][0] + "\t" + ids[i][1] + "\t" + markers.length + "\t"
                                + (2 * cnt) + "\t" + cnt2 + "\t" + ext.formDeci(scoreSum, 3));
          }
        }

        scoreWriter.flush();
        scoreWriter.close();
      }
    }
  }

  private static StrandOps.AlleleOrder determineAlleleOrder(String[] alleles, HitMarker hitMarker) {
    return StrandOps.determineStrandConfig(alleles,
                                           new String[] {hitMarker.getEffectAllele(),
                                                         hitMarker.getNonEffectAllele()})
                    .getAlleleOrder();
  }

  // private void runPlink(Study study) {
  // for (String dFile : dataFiles) {
  // String dataFile = ext.rootOf(dFile, false);
  //
  // for (java.util.Map.Entry<String, Constraint> filePrefix : analysisConstraints.entrySet()) {
  // File prefDir = new File(study.studyDir + dataFile + "/" + filePrefix.getKey() + "/");
  // if (!prefDir.exists()) {
  // log.report(ext.getTime()+"]\tError - no subfolder for '" + filePrefix.getKey() + "' analysis");
  // continue;
  // }
  // if ((new File(prefDir + "/plink.profile")).exists()) {
  // log.report(ext.getTime()+"]\tPlink analysis results file already exists! [ --> '" + prefDir +
  // "/plink.profile" + "']");
  // continue;
  // }
  // String mkrDataFile = prefDir + "/subsetData_" + filePrefix.getKey() + ".xln";
  // log.report(ext.getTime()+"]\tRunning plink command [ --> '");
  // String cmd = "plink" + /*(plink2 ? "2" : "") +*/ " --noweb --bfile ../../" + study.plinkPref +
  // " --score " + mkrDataFile;
  // log.report(cmd + "']");
  // /*boolean results = */CmdLine.run(cmd, prefDir.getAbsolutePath());
  // }
  // }
  // }

  private void runRegression(Study study) {
    for (String dFile : metaFiles) {
      String dataFile = ext.rootOf(dFile, false);

      for (java.util.Map.Entry<String, Constraint> filePrefix : analysisConstraints.entrySet()) {
        if (study.data.get(dataFile + "\t" + filePrefix.getKey()).isEmpty()) {
          continue;
        }
        try {
          File prefDir = new File(study.studyDir + dataFile + "/" + filePrefix.getKey() + "/");

          String scoreFile = prefDir + "/" + SCORE_FILE;

          HashMap<String, Double> scoreData = new HashMap<String, Double>();

          BufferedReader scoreReader = Files.getAppropriateReader(scoreFile);
          String line = scoreReader.readLine();
          while ((line = scoreReader.readLine()) != null) {
            String[] parts = line.split("\t");
            String score = parts[parts.length - 1];
            scoreData.put(parts[0] + "\t" + parts[1], Double.parseDouble(score));
          }
          scoreReader.close();

          // if (!(new File(prefDir + "/scores.hist").exists())) {
          // double[] scores = new double[scoreData.size()];
          // int ind = 0;
          // for (double data : scoreData.values()) {
          // scores[ind] = data;
          // ind++;
          // }
          // String[] unq = Array.unique(Array.toStringArray(scores));
          // if (unq.length == 1) {
          // log.report(ext.getTime()+"]\tError - no variance in scores for " + dataFile + " / " +
          // filePrefix.getKey() + " -- no .hist file created");
          // } else {
          // Files.write((new Histogram(scores)).getSummary().trim(), prefDir + "/scores.hist");
          // }
          // }

          for (int i = 0; i < study.phenoFiles.size(); i++) {
            PhenoData pd = study.phenoData.get(study.phenoFiles.get(i));
            ArrayList<Double> depData = new ArrayList<Double>();
            ArrayList<double[]> baselineIndeps = new ArrayList<double[]>();
            ArrayList<double[]> indepData = new ArrayList<double[]>();

            for (java.util.Map.Entry<String, PhenoIndiv> indiv : pd.indivs.entrySet()) {
              if (scoreData.containsKey(indiv.getKey())) {
                PhenoIndiv pdi = pd.indivs.get(indiv.getKey());
                depData.add(pdi.depvar);
                double[] baseData = new double[pd.covars.size()];
                double[] covarData = new double[pd.covars.size() + 1];
                covarData[0] = scoreData.get(indiv.getKey());
                for (int k = 1; k < pd.covars.size() + 1; k++) {
                  Double d;
                  d = pdi.covars.get(pd.covars.get(k - 1));
                  if (d == null) {
                    log.reportError("Covar value missing for individual: " + indiv.getKey() + " | "
                                    + pd.covars.get(k - 1));
                  } else {
                    baseData[k - 1] = d.doubleValue();
                  }
                  d = pdi.covars.get(pd.covars.get(k - 1));
                  if (d == null) {
                    log.reportError("Covar value missing for individual: " + indiv.getKey() + " | "
                                    + pd.covars.get(k - 1));
                  } else {
                    covarData[k] = pdi.covars.get(pd.covars.get(k - 1));
                  }
                }
                baselineIndeps.add(baseData);
                indepData.add(covarData);
              }
            }

            double[][] baseCovars = new double[baselineIndeps.size()][];
            double[][] covars = new double[indepData.size()][];
            for (int k = 0; k < covars.length; k++) {
              covars[k] = indepData.get(k);
              baseCovars[k] = baselineIndeps.get(k);
            }
            RegressionModel baseModel = RegressionModel.determineAppropriate(Doubles.toArray(depData),
                                                                             baseCovars, false,
                                                                             true);
            RegressionModel model = RegressionModel.determineAppropriate(Doubles.toArray(depData),
                                                                         covars, false, true);

            RegressionResult rr = new RegressionResult();
            if (model.analysisFailed()) {
              rr.baseRSq = baseModel.analysisFailed() ? Double.NaN : baseModel.getRsquare();
              rr.beta = Double.NaN;
              rr.se = Double.NaN;
              rr.rsq = Double.NaN;
              rr.pval = Double.NaN;
              rr.num = depData.size();
              rr.logistic = model.isLogistic();
            } else {
              int ind = -1;
              for (int l = 0; l < model.getVarNames().length; l++) {
                if ("Indep 1".equals(model.getVarNames()[l])) {
                  ind = l;
                  break;
                }
              }
              if (ind == -1) {
                rr.baseRSq = baseModel.analysisFailed() ? Double.NaN : baseModel.getRsquare();
                rr.beta = Double.NaN;
                rr.se = Double.NaN;
                rr.rsq = model.getRsquare();
                rr.pval = Double.NaN;
                rr.num = depData.size();
                rr.logistic = model.isLogistic();
                rr.stats = Double.NaN;
              } else {
                rr.baseRSq = baseModel.analysisFailed() ? Double.NaN : baseModel.getRsquare();
                rr.beta = model.getBetas()[ind];
                rr.se = model.getSEofBs()[ind];
                rr.rsq = model.getRsquare();
                rr.pval = model.getSigs()[ind];
                rr.num = depData.size();
                rr.logistic = model.isLogistic();
                rr.stats = model.getStats()[ind];
              }
            }
            study.regressions.get(filePrefix.getKey()).get(dataFile).put(pd.phenoName, rr);

            baseModel = null;
            model = null;
          }
        } catch (IOException e) {
          e.printStackTrace();
        }
      }
    }
  }

  private void writeResults() {
    PrintWriter writer;

    String resFile = metaDir + "results.xln";
    log.report(ext.getTime() + "]\tWriting regression results... [ --> " + resFile + "]");
    writer = Files.getAppropriateWriter(resFile);
    writer.println(REGRESSION_HEADER);

    for (Study study : studies) {
      for (String dFile : metaFiles) {
        String dataFile = ext.rootOf(dFile, false);
        for (java.util.Map.Entry<String, Constraint> filePrefix : analysisConstraints.entrySet()) {
          HashMap<String, RegressionResult> phenoResults = study.regressions.get(filePrefix.getKey())
                                                                            .get(dataFile);

          String resultPrefix = study.studyName + "\t" + dataFile + "\t"
                                + ext.formSciNot(filePrefix.getValue().indexThreshold, 5, false)
                                + "\t";

          String middle = (new StringBuilder()).append(dataCounts.get(dFile)
                                                                 .get(filePrefix.getKey()))
                                               .append("\t")
                                               .append(study.hitSnpCounts.get(filePrefix.getKey())
                                                                         .get(dataFile))
                                               .append("\t")
                                               .append(study.hitWindowCnts.get(filePrefix.getKey())
                                                                          .get(dataFile))
                                               .append("\t")
                                               .append(study.scores.get(filePrefix.getKey())
                                                                   .get(dataFile)[0])
                                               .append("\t")
                                               .append(study.scores.get(filePrefix.getKey())
                                                                   .get(dataFile)[1])
                                               .append("\t").toString();
          if (study.phenoFiles.isEmpty()) {
            RegressionResult rr = new RegressionResult();
            rr.dummy();

            String pvalExcl = rr.num == 0 ? "."
                                          : (rr.logistic ? "=NORMSDIST(" + Math.sqrt(rr.stats) + ")"
                                                         : "=TDIST(" + Math.abs(rr.stats) + ","
                                                           + rr.num + ",2)");

            StringBuilder sb = new StringBuilder(resultPrefix).append("--").append("\t")
                                                              .append(rr.baseRSq).append("\t")
                                                              .append(rr.rsq).append("\t")
                                                              .append((Double.isNaN(rr.rsq) ? Double.NaN
                                                                                            : (Double.isNaN(rr.baseRSq) ? rr.rsq
                                                                                                                        : (new BigDecimal(rr.rsq
                                                                                                                                          + "")).subtract(new BigDecimal(rr.baseRSq
                                                                                                                                                                         + "")))))
                                                              .append("\t").append(rr.pval)
                                                              .append("\t").append(rr.beta)
                                                              .append("\t").append(rr.se)
                                                              .append("\t").append(rr.num)
                                                              .append("\t").append(middle)
                                                              .append(pvalExcl);

            writer.println(sb.toString());
          } else {
            for (String pheno : study.phenoFiles) {
              RegressionResult rr = phenoResults.get(pheno);

              if (rr == null) {
                rr = new RegressionResult();
                rr.dummy();
              }
              String pvalExcl = rr.num == 0 ? "."
                                            : (rr.logistic ? "=NORMSDIST(" + Math.sqrt(rr.stats)
                                                             + ")"
                                                           : "=TDIST(" + Math.abs(rr.stats) + ","
                                                             + rr.num + ",2)");

              StringBuilder sb = new StringBuilder(resultPrefix).append(pheno).append("\t")
                                                                .append(rr.baseRSq).append("\t")
                                                                .append(rr.rsq).append("\t")
                                                                .append((Double.isNaN(rr.rsq) ? Double.NaN
                                                                                              : (Double.isNaN(rr.baseRSq) ? rr.rsq
                                                                                                                          : (new BigDecimal(rr.rsq
                                                                                                                                            + "")).subtract(new BigDecimal(rr.baseRSq
                                                                                                                                                                           + "")))))
                                                                .append("\t").append(rr.pval)
                                                                .append("\t").append(rr.beta)
                                                                .append("\t").append(rr.se)
                                                                .append("\t").append(rr.num)
                                                                .append("\t").append(middle)
                                                                .append(pvalExcl);

              writer.println(sb.toString());
            }
          }
        }
      }
    }

    writer.flush();
    writer.close();
  }

  public static final String COMMAND_GENESCORE = "geneScore";

  public static void fromParameters(String filename, Logger log) {
    String dir = ext.verifyDirFormat(new File(ext.parseDirectoryOfFile(filename)).getAbsolutePath());
    List<String> params;
    params = Files.parseControlFile(filename, COMMAND_GENESCORE,
                                    new String[] {ARG_WORKING_DIR + dir,
                                                  ARG_INDEX_THRESH + DEFAULT_INDEX_THRESHOLD,
                                                  ARG_WINDOW_SIZE + DEFAULT_WINDOW_MIN_SIZE_PER_SIDE,
                                                  ARG_WINDOW_EXT + DEFAULT_WINDOW_EXTENSION_THRESHOLD},
                                    log);
    if (params == null) {
      setupCRF(filename, log);
      return;
    }
    if (params != null) {
      params.add("log=" + log.getFilename());
      main(ArrayUtils.toStringArray(params));
    }
  }

  private static void setupCRF(String crfFile, Logger log) {
    String dir = ext.parseDirectoryOfFile(crfFile);
    String[] files = (new File(dir)).list();
    HashSet<String> potentialRoots = new HashSet<>();
    HashSet<String> preprocess = new HashSet<>();
    HashSet<String> metas = new HashSet<>();
    HashSet<String> phenos = new HashSet<>();
    for (String f : files) {
      if (f.endsWith(".bed") || f.endsWith(".bim") || f.endsWith(".fam")) {
        potentialRoots.add(ext.rootOf(f));
      } else if (f.endsWith(".result") || f.endsWith(".results")) {
        preprocess.add(f);
      } else if (f.endsWith(".meta")) {
        metas.add(f);
      } else if (f.endsWith(".pheno")) {
        phenos.add(f);
      }
    }

    HashSet<String> validData = new HashSet<>();
    for (String s : potentialRoots) {
      if (PSF.Plink.allFilesExist(dir + s, true)) {
        validData.add(s);
      }
    }

    createDataFile(dir, validData, phenos);

    log.reportTime("Processing " + preprocess.size() + " results files...");
    preprocessDataFiles(preprocess.toArray(new String[preprocess.size()]), log);
    log.reportTime("Done!");
  }

  private static void createDataFile(String baseDir, HashSet<String> plinkRoots,
                                     HashSet<String> phenoFiles) {
    String dir = baseDir + "plink/";
    new File(dir).mkdir();
    PrintWriter writer = Files.getAppropriateWriter(dir + "data.txt");
    for (String plinkRoot : plinkRoots) {
      writer.println(plinkRoot + "\t" + getFull(baseDir + PSF.Plink.getBED(plinkRoot)) + "\t"
                     + getFull(baseDir + PSF.Plink.getBIM(plinkRoot)) + "\t"
                     + getFull(baseDir + PSF.Plink.getFAM(plinkRoot)));
    }
    writer.flush();
    writer.close();

    for (String s : phenoFiles) {
      Files.copyFile(baseDir + s, dir + s);
    }
  }

  private static String getFull(String file) {
    return (new File(file)).getAbsolutePath();
  }

  // private void writeToForestInput() {
  // String resultsFile = metaDir + "regressions.out";
  // PrintWriter writer = Files.getAppropriateWriter(resultsFile);
  // StringBuilder header = new StringBuilder("Name\tbeta\tse");
  // for (Study study : studies) {
  // header.append("\tbeta.").append(study.studyName).append("\tse.").append(study.studyName);
  // }
  // writer.println(header);
  //
  // StringBuilder dataSB;
  // for (String constraint : analysisConstraints.keySet()) {
  // dataSB = new StringBuilder("GeneScore_").append(constraint);
  // dataSB.append("\t");
  // // TODO meta beta
  // dataSB.append("\t");
  // // TODO meta se
  // for (Study study : studies) {
  // dataSB.append("\t").append(study.regressionResults.get(constraint)[0]).append("\t").append(study.regressionResults.get(constraint)[1]);
  // }
  // writer.println(dataSB.toString());
  // }
  //
  // writer.flush();
  // writer.close();
  // }

  public static void main(String[] args) {
    int numArgs = args.length;

    String workDir = null;
    String logFile = null;

    float[] iT = new float[] {DEFAULT_INDEX_THRESHOLD};
    int[] mZ = new int[] {DEFAULT_WINDOW_MIN_SIZE_PER_SIDE};
    float[] wT = new float[] {DEFAULT_WINDOW_EXTENSION_THRESHOLD};
    double mT = DEFAULT_MIN_MISS_THRESH;

    String[] processList = null;
    boolean process = false;
    // boolean test = true;
    // if (test) {
    // preprocessDataFiles(new String[]{
    // "D:/GeneScorePipe/Telomere/telo.xln"
    // });
    // "D:/GeneScorePipe/Cancer/InputCancer.xln",
    // "D:/GeneScorePipe/height/GeneScorePipeline/metas/ExtremeHeight.xln",
    // "D:/GeneScorePipe/height/GeneScorePipeline/metas/height_full.xln",
    // "D:/height/GeneScorePipeline/HeightScoring/TannerSexCombined.xln"
    // });
    // return;
    // }

    String usage = "\n"
                   + "GeneScorePipeline is a convention-driven submodule.  It relies on a standard folder structure and file naming scheme:\n"
                   + "\tThe directory and file structure must conform to the following:\n"
                   + "\t\t>Root Directory ['workDir' argument]\n" + "\t\t\t>SNP Effect files:\n"
                   + "\t\t\t\t-Effect files must end with '.meta'.\n"
                   + "\t\t\t\t-Effect files may be hand-constructed, or may be generated with the 'preprocess' command from a .xln file\n"
                   + "\t\t\t\t-Effect files contain, at minimum, SNP, Freq, P-value, and Beta/Effect, and, if created with the preprocessor, will include any additional information present in the .xln file\n"
                   + "\t\t\t\t-HitWindows analysis will be run on SNP Effect files, with results being used in regression analysis; the\n"
                   + "\t\t\t\tadditional arguments to GeneScorePipeline affect only the HitWindows processing.\n"
                   + "\t\t\t>Data Source Directory 1\n"
                   + "\t\t\t\t>data.txt file [defines location of data, which may be in an arbitrary location in the filesystem]\n"
                   + "\t\t\t\t\tExample:\n"
                   + "\t\t\t\t\tdataLabel1\tfullPathDataFile1\tFullPathMapFile1\tFullPathIdFile1\n"
                   + "\t\t\t\t\tdataLabel2\tfullPathDataFile2\tFullPathMapFile2\tFullPathIdFile2\n"
                   + "\t\t\t\t\tdataLabel3\tdir1\tdataFileExt1\tmapFileExt1\tidFile3\n"
                   + "\t\t\t\t\tdataLabel4\tdir2\tdataFileExt2\tmapFileExt2\tidFile4\n"
                   + "\t\t\t\t\t>For data files that contain map or ID info, 'FullPathMapFile' / 'FullPathIdFile' / 'mapFileExt1' / 'idFile' can be replaced with a period ('.')."
                   + "\t\t\t\t>PhenoData1.pheno file [Any '.pheno' files will be included as covariate data in the regression analysis]\n"
                   + "\t\t\t\t\t[Note: if data is in PLINK format and contains valid affected status information, an AFFECTED.PHENO file will be created]\n"
                   +

                   "\t\t\t\t>Pheno2.pheno file\n" + "\t\t\t>Data Source Directory 2\n"
                   + "\t\t\t\t>data.txt file\n" + "\t\t\t\t>Pheno3.pheno file\n" + "\t\t\t>...\n"
                   + "\n" + "\n" + "gwas.GeneScorePipeline requires 1+ arguments\n"
                   + "   (1) Pre-process data files (i.e. process=path/to/file1.xln,path/to/file2.xln (not the default)) \n"
                   + "  OR\n"
                   + "   (1) Metastudy directory root, containing subdirectories for each study (i.e. workDir=C:/ (not the default))\n"
                   + "       OPTIONAL:\n"
                   + "   (2) p-value threshold (or comma-delimited list) for index SNPs (i.e. "
                   + ARG_INDEX_THRESH + DEFAULT_INDEX_THRESHOLD + " (default))\n"
                   + "   (3) minimum num bp per side of window (or comma delimited list) (i.e. "
                   + ARG_WINDOW_SIZE + DEFAULT_WINDOW_MIN_SIZE_PER_SIDE + " (default))\n"
                   + "   (4) p-value threshold to extend the window (or comma delimited list) (i.e. "
                   + ARG_WINDOW_EXT + DEFAULT_WINDOW_EXTENSION_THRESHOLD + " (default))\n"
                   + "   (5) minimum ratio of missing data for an individual's gene loading score to be included in the final analysis (i.e. "
                   + ARG_MISS_THRESH + DEFAULT_MIN_MISS_THRESH + " (default))\n" +
                   // " (8) Number of threads to use for computation (i.e. threads=" + threads + "
                   // (default))\n" +
                   "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith(ARG_WORKING_DIR)) {
        workDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("process=")) {
        processList = arg.split("=")[1].split(",");
        process = true;
      } else if (arg.startsWith(ARG_INDEX_THRESH)) {
        String[] lst = arg.split("=")[1].split(",");
        int cntValid = 0;
        for (String poss : lst) {
          if (ext.isValidDouble(poss)) {
            cntValid++;
          }
        }
        iT = new float[cntValid];
        int ind = 0;
        for (String poss : lst) {
          if (ext.isValidDouble(poss)) {
            iT[ind] = Float.parseFloat(poss);
            ind++;
          }
        }
        numArgs--;
      } else if (arg.startsWith(ARG_WINDOW_SIZE)) {
        String[] lst = arg.split("=")[1].split(",");
        int cntValid = 0;
        for (String poss : lst) {
          if (ext.isValidDouble(poss)) {
            cntValid++;
          }
        }
        mZ = new int[cntValid];
        int ind = 0;
        for (String poss : lst) {
          if (ext.isValidInteger(poss)) {
            mZ[ind] = Integer.parseInt(poss);
            ind++;
          }
        }
        numArgs--;
      } else if (arg.startsWith(ARG_WINDOW_EXT)) {
        String[] lst = arg.split("=")[1].split(",");
        int cntValid = 0;
        for (String poss : lst) {
          if (ext.isValidDouble(poss)) {
            cntValid++;
          }
        }
        wT = new float[cntValid];
        int ind = 0;
        for (String poss : lst) {
          if (ext.isValidDouble(poss)) {
            wT[ind] = Float.parseFloat(poss);
            ind++;
          }
        }
        numArgs--;
      } else if (arg.startsWith(ARG_MISS_THRESH)) {
        mT = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logFile = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0 || args.length == 0) {
      System.err.println(usage);
      System.exit(1);
    }

    Logger log = new Logger(logFile, true);

    if (process) {
      preprocessDataFiles(processList, log);
      return;
    }

    File dir = new File(workDir);
    if (!dir.isDirectory()) {
      System.err.println("Error - argument 'workDir' must be a valid directory");
      System.exit(1);
    }
    // if (regress && !runPlink) {
    // System.err.println("Error - '-runPlink' option is required for '-regress' option");
    // System.exit(1);
    // }
    // if (writeHist && !runPlink) {
    // System.err.println("Error - '-runPlink' option is required for '-writeHist' option");
    // System.exit(1);
    // }
    GeneScorePipeline gsp = new GeneScorePipeline(workDir, iT, mZ, wT, mT, log);
    gsp.runPipeline();
  }

}
