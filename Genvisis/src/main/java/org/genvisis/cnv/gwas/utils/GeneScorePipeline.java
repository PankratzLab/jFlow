package org.genvisis.cnv.gwas.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.Vector;
import java.util.stream.Collectors;
import org.apache.commons.math3.stat.inference.TTest;
import org.apache.commons.math3.stat.inference.WilcoxonSignedRankTest;
import org.genvisis.cnv.Resources;
import org.genvisis.cnv.Resources.CHROMOSOME;
import org.genvisis.cnv.gwas.utils.MergeExtractPipeline.DataSource;
import org.genvisis.cnv.plots.AFPlot;
import org.genvisis.cnv.plots.AFPlot.POPULATION;
import org.genvisis.seq.GenomeBuild;
import org.genvisis.seq.manage.StrandOps;
import org.genvisis.seq.manage.StrandOps.AlleleOrder;
import org.genvisis.seq.manage.StrandOps.CONFIG;
import org.pankratzlab.common.Aliases;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.GenomicPosition;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.bioinformatics.Sequence;
import org.pankratzlab.common.filesys.Positions;
import org.pankratzlab.common.parsing.AbstractColumnFilter;
import org.pankratzlab.common.parsing.AbstractFileColumn;
import org.pankratzlab.common.parsing.AliasedFileColumn;
import org.pankratzlab.common.parsing.ColumnFilter;
import org.pankratzlab.common.parsing.ColumnFilters;
import org.pankratzlab.common.parsing.DataLine;
import org.pankratzlab.common.parsing.DoubleFilter;
import org.pankratzlab.common.parsing.DoubleWrapperColumn;
import org.pankratzlab.common.parsing.ExplicitIndexedFileColumn;
import org.pankratzlab.common.parsing.FileColumn;
import org.pankratzlab.common.parsing.FileParser;
import org.pankratzlab.common.parsing.FileParserFactory;
import org.pankratzlab.common.parsing.IndexedFileColumn;
import org.pankratzlab.common.parsing.ParseFailureException;
import org.pankratzlab.common.parsing.StandardFileColumns;
import org.pankratzlab.common.stats.Maths.COMPARISON;
import org.pankratzlab.common.stats.ProbDist;
import org.pankratzlab.common.stats.RegressionModel;
import org.pankratzlab.utils.bioinformatics.MapSNPsAndGenes;
import org.pankratzlab.utils.filesys.SnpMarkerSet;
import org.pankratzlab.utils.gwas.DosageData;
import org.pankratzlab.utils.gwas.DosageData.Trio;
import org.pankratzlab.utils.gwas.windows.HitWindows;
import com.google.common.base.Joiner;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Maps;
import com.google.common.collect.Multiset;
import com.google.common.collect.Table;
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

  private static final String CROSS_FILTERED_DATAFILE = "bimData.xln";
  private static final String DATA_SOURCE_FILENAME = "data.txt";

  private static final String MARKER_COL_NAME = "Marker";
  private static final String EFFECT_ALLELE_COL_NAME = "EffectAllele";
  private static final String NON_EFFECT_ALLELE_COL_NAME = "NonEffectAllele";
  private static final String BETA_COL_NAME = "beta";
  private static final String REGRESSION_HEADER = new StringJoiner("\t").add("STUDY")
                                                                        .add("DATAFILE")
                                                                        .add("INDEX-THRESHOLD")
                                                                        .add("FACTOR")
                                                                        .add("BASE-R-SQR")
                                                                        .add("R-SQR").add("R-DIFF")
                                                                        .add("P-VALUE")
                                                                        .add("EXCEL-SIG")
                                                                        .add("BETA").add("SE")
                                                                        .add("NUM").add("CASES")
                                                                        .add("CONTROLS")
                                                                        .add("PAIRED-T-P-VALUE")
                                                                        .add("WILCOXON-SIGNED-RANK-P-VALUE")
                                                                        .add("PAIRED-STAT-NUM-TRIOS")
                                                                        .add("#sigInMeta")
                                                                        .add("#indexVariantsInMeta")
                                                                        .add("#indexVariantsInDataset")
                                                                        .add("B-F-SCORE")
                                                                        .add("INVCHI-SCORE")
                                                                        .toString();
  private final String metaDir;

  private float[] indexThresholds = new float[] {DEFAULT_INDEX_THRESHOLD};
  private int[] windowMinSizePerSides = new int[] {DEFAULT_WINDOW_MIN_SIZE_PER_SIDE};
  private float[] windowExtensionThresholds = new float[] {DEFAULT_WINDOW_EXTENSION_THRESHOLD};
  private double minMissThresh = DEFAULT_MIN_MISS_THRESH;

  // private int numThreads = 1;
  // private boolean runPlink = false;
  // private boolean runRegression = false;
  // private boolean writeHist = false;

  private final ArrayList<String> metaFiles = new ArrayList<>();
  private final ArrayList<Study> studies = new ArrayList<>();
  private final HashMap<String, Constraint> analysisConstraints = new HashMap<>();
  private final HashMap<String, HashMap<String, Integer>> dataCounts = new HashMap<>();

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
    HashMap<String, DosageData> data = new HashMap<>();

    ArrayList<String> phenoFiles = new ArrayList<>();

    HashMap<String, PhenoData> phenoData = new HashMap<>();

    // constraint -> datafile
    Table<String, String, TrioScoreTest.Results> trioScoreTests = HashBasedTable.create();
    // constraint -> datafile -> phenofile
    HashMap<String, HashMap<String, HashMap<String, RegressionResult>>> regressions = new HashMap<>();
    // constraint -> datafile
    HashMap<String, HashMap<String, double[]>> scores = new HashMap<>();
    // constraint -> datafile
    HashMap<String, HashMap<String, Integer>> hitSnpCounts = new HashMap<>();
    // constraint -> datafile
    HashMap<String, HashMap<String, Integer>> hitWindowCnts = new HashMap<>();
    // constraint -> datafile
    HashMap<String, HashMap<String, Integer>> dataCounts = new HashMap<>();
    Map<String, Map<String, HitMarker>> markerData = new HashMap<>();

    /**
     * Returns a hashSet containing all markers present in the data files that are present in
     * hitMkrSet.
     *
     * @param hitMkrSet
     * @return
     */
    public HashSet<String> retrieveMarkers(String dataKey, Set<String> hitMkrSet) {
      HashSet<String> returnMarkers = new HashSet<>();
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

    public void loadDataSources(String dataKey, Map<String, int[]> markerLocationMap) {
      String[] mkrsArray = markerLocationMap.keySet().toArray(new String[0]);
      SnpMarkerSet markerSet = new SnpMarkerSet(mkrsArray, false, log);
      markerSet.parseSNPlocations(log);
      int[][] markerLocations = markerSet.getChrAndPositionsAsInts();
      dataSources = MergeExtractPipeline.parseDataFile(null, markerLocations, null, dataSource, 0,
                                                       log);
      if (dataSources.size() == 0) {
        // error
        log.reportError("Error - no data sources loaded from file: " + dataSource + "; expected "
                        + markerLocations.length);
      } else {
        final String serOutput = studyDir + ext.replaceWithLinuxSafeCharacters(dataKey)
                                 + "_dosage.ser";
        if (!Files.exists(serOutput)) {
          log.reportTime("Loading data file " + dataSources.get(0).dataFile);
          DosageData d0 = new DosageData(dataSources.get(0).dataFile, dataSources.get(0).idFile,
                                         dataSources.get(0).mapFile, null, mkrsArray, true, log);
          if (dataSources.size() > 1) {
            for (int i = 1; i < dataSources.size(); i++) {
              log.reportTime("Loading data file " + dataSources.get(i).dataFile);
              new DosageData(dataSources.get(i).dataFile, dataSources.get(i).idFile,
                             dataSources.get(i).mapFile, null, mkrsArray, markerLocationMap, null,
                             true, log);
              DosageData d1 = new DosageData(dataSources.get(i).dataFile, dataSources.get(i).idFile,
                                             dataSources.get(i).mapFile, null, mkrsArray, true,
                                             log);
              d0 = DosageData.combine(d0, d1, DosageData.COMBINE_OP.EITHER_IF_OTHER_MISSING, false,
                                      0, log);
              System.gc();
            }
          }
          if (d0.isEmpty()) {
            log.reportError("no data for key: " + dataKey);
          }
          data.put(dataKey, d0);
          d0.serialize(serOutput);
          System.gc();
        } else {
          log.reportTime(serOutput + " already exists, loading previously loaded data for "
                         + dataSource);
          data.put(dataKey, DosageData.load(serOutput));
        }
      }
    }
  }

  private static class RegressionResult {

    private static class Builder {

      private boolean logistic;
      private double rsq;
      private double baseRSq;
      private double pval;
      private double beta;
      private double se;
      private int num;
      private int nCases;
      private int nControls;
      private double stats;

      private Builder() {

      }

      private RegressionResult build() {
        return new RegressionResult(this);
      }

      /**
       * @param logistic the logistic to set
       */
      private Builder setLogistic(boolean logistic) {
        this.logistic = logistic;
        return this;
      }

      /**
       * @param rsq the rsq to set
       */
      private Builder setRsq(double rsq) {
        this.rsq = rsq;
        return this;
      }

      /**
       * @param baseRSq the baseRSq to set
       */
      private Builder setBaseRSq(double baseRSq) {
        this.baseRSq = baseRSq;
        return this;
      }

      /**
       * @param pval the pval to set
       */
      private Builder setPval(double pval) {
        this.pval = pval;
        return this;
      }

      /**
       * @param beta the beta to set
       */
      private Builder setBeta(double beta) {
        this.beta = beta;
        return this;
      }

      /**
       * @param se the se to set
       */
      private Builder setSe(double se) {
        this.se = se;
        return this;
      }

      /**
       * @param num the num to set
       */
      private Builder setNum(int num) {
        this.num = num;
        return this;
      }

      /**
       * @param nCases the nCases to set
       */
      private Builder setnCases(int nCases) {
        this.nCases = nCases;
        return this;
      }

      /**
       * @param nControls the nControls to set
       */
      private Builder setnControls(int nControls) {
        this.nControls = nControls;
        return this;
      }

      /**
       * @param stats the stats to set
       */
      private Builder setStats(double stats) {
        this.stats = stats;
        return this;
      }

    }

    private final boolean logistic;
    private final double rsq;
    private final double baseRSq;
    private final double pval;
    private final double beta;
    private final double se;
    private final int num;
    private final int nCases;
    private final int nControls;
    private final double stats;

    private RegressionResult(Builder builder) {
      super();
      logistic = builder.logistic;
      rsq = builder.rsq;
      baseRSq = builder.baseRSq;
      pval = builder.pval;
      beta = builder.beta;
      se = builder.se;
      num = builder.num;
      nCases = builder.nCases;
      nControls = builder.nControls;
      stats = builder.stats;
    }

    /**
     * @return the logistic
     */
    boolean isLogistic() {
      return logistic;
    }

    /**
     * @return the rsq
     */
    double getRsq() {
      return rsq;
    }

    /**
     * @return the baseRSq
     */
    double getBaseRSq() {
      return baseRSq;
    }

    /**
     * @return the pval
     */
    double getPval() {
      return pval;
    }

    /**
     * @return the beta
     */
    double getBeta() {
      return beta;
    }

    /**
     * @return the se
     */
    double getSe() {
      return se;
    }

    /**
     * @return the num
     */
    int getNum() {
      return num;
    }

    /**
     * @return the nCases
     */
    int getnCases() {
      return nCases;
    }

    /**
     * @return the nControls
     */
    int getnControls() {
      return nControls;
    }

    /**
     * @return the stats
     */
    double getStats() {
      return stats;
    }

    private static Builder builder() {
      return new Builder();
    }

    private static RegressionResult dummy() {
      return builder().setLogistic(true).setRsq(Double.NaN).setBaseRSq(Double.NaN)
                      .setPval(Double.NaN).setBeta(Double.NaN).setSe(Double.NaN).setNum(0)
                      .setnCases(0).setnControls(0).setStats(Double.NaN).build();
    }

  }

  private static class TrioScoreTest {

    private static class Results {

      private final int trioCount;
      private final double wilcoxonPVal;
      private final double pairedTPVal;

      private Results(TrioScoreTest trioScoreTest) {
        super();
        double[] caseScoresArray = Doubles.toArray(trioScoreTest.caseScores);
        double[] parentMeansArray = Doubles.toArray(trioScoreTest.parentScores);
        if (caseScoresArray.length == parentMeansArray.length) {
          trioCount = caseScoresArray.length;
        } else {
          throw new IllegalStateException("Case and mean parent score arrays must be parallel");
        }
        this.wilcoxonPVal = new WilcoxonSignedRankTest().wilcoxonSignedRankTest(caseScoresArray,
                                                                                parentMeansArray,
                                                                                false);
        this.pairedTPVal = new TTest().pairedTTest(caseScoresArray, parentMeansArray);
      }

      /**
       * @return the number of trios used in the paired stats calculated
       */
      public int getTrioCount() {
        return trioCount;
      }

      /**
       * @return the wilcoxonPVal
       */
      public double getWilcoxonPVal() {
        return wilcoxonPVal;
      }

      /**
       * @return the pairedTPVal
       */
      public double getPairedTPVal() {
        return pairedTPVal;
      }

    }

    private final List<Double> caseScores;
    private final List<Double> parentScores;

    public TrioScoreTest() {
      caseScores = new ArrayList<>();
      parentScores = new ArrayList<>();
    }

    public void addScore(double caseScore, double fatherScore, double motherScore) {
      caseScores.add(caseScore);
      parentScores.add((fatherScore + motherScore) / 2.0);
    }

    public Results runTests() {
      if (caseScores.size() < 2) return null;
      return new Results(this);
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
    HashMap<String, PhenoIndiv> indivs = new HashMap<>();
    ArrayList<String> covars = new ArrayList<>();
  }

  private class PhenoIndiv {

    private final String fid;
    private final String iid;
    private final double depvar;
    private final Map<String, Double> covars;

    /**
     * @param fid
     * @param iid
     * @param depvar
     */
    public PhenoIndiv(String fid, String iid, double depvar) {
      super();
      this.fid = fid;
      this.iid = iid;
      this.depvar = depvar;
      this.covars = new HashMap<>();
    }

    /**
     * @return the fid
     */
    public String getFid() {
      return fid;
    }

    /**
     * @return the iid
     */
    public String getIid() {
      return iid;
    }

    /**
     * @return the depvar
     */
    public double getDepvar() {
      return depvar;
    }

    /**
     * @return an unmodifiable view of the covars
     */
    public Map<String, Double> getCovars() {
      return Collections.unmodifiableMap(covars);
    }

    public void addCovar(String name, double value) {
      if (covars.putIfAbsent(name, value) != null) {
        log.reportTimeWarning("Duplicate covars named '" + name + "', only one will be used");
      }
    }

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

  private static HashMap<String, Double> get1000GFreq(HashMap<String, int[]> markerMap,
                                                      AFPlot.POPULATION population,
                                                      GenomeBuild build, Logger log) {
    AFPlot.POPULATION pop = population == null ? POPULATION.ALL : population;
    HashMap<Integer, HashSet<String>> mkrsByChr, nonRSByChr;
    HashSet<String> chrMkrs, nonRS;
    HashMap<String, Double> markerFreqs = new HashMap<>();

    mkrsByChr = new HashMap<>();
    nonRSByChr = new HashMap<>();
    for (java.util.Map.Entry<String, int[]> entry : markerMap.entrySet()) {
      chrMkrs = mkrsByChr.get(entry.getValue()[0]);
      nonRS = nonRSByChr.get(entry.getValue()[0]);
      if (chrMkrs == null) {
        chrMkrs = new HashSet<>();
        mkrsByChr.put(entry.getValue()[0], chrMkrs);
      }
      if (nonRS == null) {
        nonRS = new HashSet<>();
        nonRSByChr.put(entry.getValue()[0], nonRS);
      }
      chrMkrs.add(entry.getKey());
      if (!entry.getKey().startsWith("rs")) {
        nonRS.add(entry.getKey());
      }
    }

    for (Integer chr : mkrsByChr.keySet()) {
      // skip this chromosome if we don't have any markers for which we need freqencies
      if (mkrsByChr.get(chr) == null || mkrsByChr.get(chr).isEmpty()) {
        continue;
      }

      String afFile = Resources.genome(build, log)
                               .chr(CHROMOSOME.valueOf("C" + Integer.toString(chr)))
                               .getG1Kphase3v5AlleleFreq().get();

      FileColumn<String> snpCol = new AliasedFileColumn("SNP", "ID");
      FileColumn<Byte> chrCol = StandardFileColumns.chr("CHROM");
      FileColumn<Integer> posCol = StandardFileColumns.pos("POS");
      FileColumn<Double> afCol = new DoubleWrapperColumn(new AliasedFileColumn(pop.name(),
                                                                               pop.colName));

      chrMkrs = mkrsByChr.get(chr);
      nonRS = nonRSByChr.get(chr);
      long parseFails = 0;

      try (FileParser parser = FileParserFactory.setup(afFile, snpCol, posCol, chrCol, afCol)
                                                .build()) {
        for (DataLine line : parser) {
          String snp = line.getString(snpCol);
          try {
            if (snp.startsWith("rs")) {
              if (chrMkrs.contains(snp)) {
                markerFreqs.put(snp, line.get(afCol));
                chrMkrs.remove(snp);
              }
            } else if (nonRS.contains(snp)) {
              markerFreqs.put(snp, line.get(afCol));
              nonRS.remove(snp);
            }
          } catch (ParseFailureException e) {
            parseFails++;
          }
        }
      } catch (IOException e) {
        e.printStackTrace();
      }

      log.reportTime("Chromosome " + chr + " -- Couldn't find " + chrMkrs.size() + " RS-ids and "
                     + nonRS.size() + " non-RS snps."
                     + (parseFails > 0 ? "  Failed to parse allele frequencies for " + parseFails
                                         + " markers."
                                       : ""));
    }

    return markerFreqs;
  }

  public static void preprocessDataFiles(String[] files, AFPlot.POPULATION pop, GenomeBuild build,
                                         Logger log) {
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
        markerMap = new HashMap<>();
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
            freqs = get1000GFreq(markerMap, pop, build, log);
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
    this.metaDir = ext.verifyDirFormat(metaDir);
    // this.numThreads = numThreads;
    // this.runPlink = plink;
    // this.runRegression = runPlink && regression;
    // this.writeHist = runPlink && histogram;
    this.minMissThresh = missThresh;
    this.indexThresholds = indexThresholds;
    windowMinSizePerSides = windowMins;
    windowExtensionThresholds = windowExtThresholds;
    if (indexThresholds.length != windowExtThresholds.length) {
      throw new IllegalArgumentException();
    }
    setFilePrefices();
    loadStudyFolders();
    // instantiate inner hashmaps:
    for (Study study : studies) {
      for (String pref : analysisConstraints.keySet()) {
        HashMap<String, HashMap<String, RegressionResult>> res = new HashMap<>();
        for (String dFile : metaFiles) {
          String dataFile = ext.rootOf(dFile, false);
          HashMap<String, RegressionResult> res2 = new HashMap<>();
          res.put(dataFile, res2);
        }
        study.regressions.put(pref, res);

        HashMap<String, Integer> cntMap = new HashMap<>();
        HashMap<String, Integer> hitMap = new HashMap<>();
        HashMap<String, Integer> cntMap2 = new HashMap<>();
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
          final float indexThresh = i;
          final int minSize = m;
          final float windowThresh;
          if (w < indexThresh) {
            log.reportError("Window extension threshold (" + w
                            + ") is more stringent than index threshold (" + indexThresh
                            + "), index threshold will be used as window extension threshold");
            windowThresh = indexThresh;
          } else windowThresh = w;
          StringBuilder prefixSB = new StringBuilder();
          prefixSB.append(ext.formSciNot(indexThresh, 4, false)).append("_")
                  .append(ext.formSciNot(minSize, 4, false)).append("_")
                  .append(ext.formSciNot(windowThresh, 4, false));
          analysisConstraints.put(prefixSB.toString(),
                                  new Constraint(indexThresh, minSize, windowThresh));
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
            dFileCnts = new HashMap<>();
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
    } else {

      for (String dFile : metaFiles) {
        HashMap<String, Constraint> threshNeed = new HashMap<>();

        HashMap<String, Integer> cnts = dataCounts.get(dFile);
        if (cnts == null) {
          threshNeed.putAll(analysisConstraints);
          cnts = new HashMap<>();
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
            int[] indices = ext.indexFactors(Aliases.PVALUES, dataHdrs, false, log, false);
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
  }

  private void runMetaHitWindowsAndLoadData() {
    String[][] factors = new String[][] {Aliases.MARKER_NAMES, Aliases.EFFECTS,
                                         Aliases.ALLELE_FREQS, Aliases.PVALUES};

    for (String dFile : metaFiles) {
      String dataFile = ext.rootOf(dFile, false);
      Map<String, int[]> hitMkrLocations = new HashMap<>();
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
              try {
                hitMkrLocations.put(result[1], new int[] {Integer.parseInt(result[2]),
                                                          Integer.parseInt(result[3])});
              } catch (NumberFormatException nfe) {
                log.reportError("Failed to parse position for result " + ArrayUtils.toStr(result));
                hitMkrLocations.put(result[1], new int[] {-1, -1});
              }
            }
          }
        }
        if (hitMkrLocations.isEmpty()) {
          FileColumn<String> mkrColumn = new ExplicitIndexedFileColumn("mkr", 0);
          FileColumn<Byte> chrColumn = StandardFileColumns.chr("chr");
          FileColumn<Integer> posColumn = StandardFileColumns.pos("pos");
          try {
            hitMkrLocations = Maps.transformValues(FileParserFactory.setup(metaDir + dFile,
                                                                           mkrColumn)
                                                                    .optionalColumns(chrColumn,
                                                                                     posColumn)
                                                                    .build().load(false, mkrColumn),
                                                   d -> new int[] {d.get(chrColumn, (byte) -1),
                                                                   d.get(posColumn, -1)});
          } catch (IOException e) {
            log.reportIOException(metaDir + dFile);
          }

          if (hitMkrLocations.isEmpty()) {
            log.reportError(".meta file was empty for " + dFile);
          } else {
            log.report(ext.getTime() + "]\tUsing all " + hitMkrLocations.size()
                       + " SNPs in .meta file");
          }
        }
        // uncomment to use all markers in dataFile

        if (!hitMkrLocations.isEmpty()) {
          // read betas and freqs for hitwindow markers
          HashMap<String, double[]> dataMarkers = new HashMap<>();
          try {
            BufferedReader reader = Files.getAppropriateReader(metaDir + dFile);
            String line = reader.readLine();
            String[] temp = line.split(PSF.Regex.GREEDY_WHITESPACE);
            int[] indices = ext.indexFactors(factors, temp, false, false, true, true, new Logger());
            while ((line = reader.readLine()) != null) {
              String mkr = line.split(PSF.Regex.GREEDY_WHITESPACE)[indices[0]];
              if (hitMkrLocations.containsKey(mkr)) {
                if ((indices[1] != -1
                     && ext.isMissingValue(line.split(PSF.Regex.GREEDY_WHITESPACE)[indices[1]]))
                    || ext.isMissingValue(line.split(PSF.Regex.GREEDY_WHITESPACE)[indices[2]])
                    || ext.isMissingValue(line.split(PSF.Regex.GREEDY_WHITESPACE)[indices[3]])) {
                  hitMkrLocations.remove(mkr);
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

            // cross-ref PLINK markers
            for (Study study : studies) {
              study.loadDataSources(dataFile + "\t" + filePrefix.getKey(), hitMkrLocations);

              HashMap<String, double[]> bimSubsetMarkers = new HashMap<>();
              HashSet<String> bimMkrSet = study.retrieveMarkers(dataFile + "\t"
                                                                + filePrefix.getKey(),
                                                                hitMkrLocations.keySet());

              // pull betas and freqs for union markers
              for (String mkr : bimMkrSet) {
                if (hitMkrLocations.containsKey(mkr)) {
                  bimSubsetMarkers.put(mkr, dataMarkers.get(mkr));
                }
              }

              // apply equation and set value for overall results
              double bimScore1 = getBetaFreqScore(bimSubsetMarkers);
              double bimScore2 = getChiDistRevScore(bimSubsetMarkers);

              HashMap<String, double[]> fileMap = study.scores.get(filePrefix.getKey());
              if (fileMap == null) {
                fileMap = new HashMap<>();
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
    ArrayList<String> fam = new ArrayList<>();
    ArrayList<String> pheno = new ArrayList<>();
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
        ArrayList<String> covars = new ArrayList<>();
        for (int i = 3; i < header.length; i++) {
          covars.add(header[i]);
        }
        pd.covars.addAll(covars);
        String temp = reader.readLine();
        do {
          String[] line = temp.split("\t");

          if (!ext.isMissingValue(line[2])) {
            String fid = line[0];
            String iid = line[1];
            double depvar = Double.parseDouble(line[2]);
            PhenoIndiv pi = new PhenoIndiv(fid, iid, depvar);
            for (int i = 3; i < line.length; i++) {
              if (!ext.isMissingValue(line[i])) {
                pi.addCovar(header[i], Double.parseDouble(line[i]));
              }
            }
            pd.indivs.put(pi.getFid() + "\t" + pi.getIid(), pi);
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

        log.report(ext.getTime() + "]\tCross-filtering data and .meta files [ --> '"
                   + crossFilterFile + "']");
        HashMap<String, GenomicPosition> mkrsMeta = new HashMap<>();
        Map<String, String[]> mkrAlleles = new HashMap<>();
        SnpMarkerSet markerSet = study.data.get(dataFile + "\t" + constraintEntry.getKey())
                                           .getMarkerSet();
        List<Integer> inval = new ArrayList<>();
        List<Integer> ambig = new ArrayList<>();

        String[] mkrNames = markerSet.getMarkerNames();
        String[][] alleles = markerSet.getAlleles();
        int[][] chrPos = markerSet.getChrAndPositionsAsInts();

        for (int i = 0; i < mkrNames.length; i++) {
          String a1 = (alleles[i][0]).toUpperCase();
          String a2 = (alleles[i][1]).toUpperCase();
          boolean validAlleles = (a1.length() > 1 || a2.length() > 1)
                                 || (Sequence.validBase(a1) && Sequence.validBase(a2));
          if (validAlleles) {
            if (!a1.equals(Sequence.flip(a2))) {
              mkrsMeta.put(mkrNames[i], new GenomicPosition((byte) chrPos[i][0], chrPos[i][1]));
              mkrAlleles.put(mkrNames[i], new String[] {a1, a2});
            } else {
              mkrsMeta.put(mkrNames[i], new GenomicPosition((byte) chrPos[i][0], chrPos[i][1]));
              mkrAlleles.put(mkrNames[i], new String[] {a1, a2});
              ambig.add(i);
            }
          } else {
            inval.add(i);
          }
        }

        if (ambig.size() > 0) {
          log.report(ext.getTime() + "]\tFound " + ambig.size() + " ambiguous markers out of "
                     + mkrNames.length
                     + " (either A/T or G/C allele pairs); make sure you have everything on the same strand, or you many run into problems!");
          if (ambig.size() < 10) {
            StringBuilder builder = new StringBuilder(mkrNames[ambig.get(0)]);
            for (int i = 1; i < ambig.size(); i++) {
              builder.append(", ").append(mkrNames[ambig.get(i)]);
            }
            log.report(ext.getTime() + "]\tAmbiguous markers: " + builder.toString());
          } else {
            PrintWriter writer = Files.getAppropriateWriter(metaDir + "ambiguous.txt");
            for (Integer i : ambig) {
              writer.println(mkrNames[i] + "\t" + chrPos[i][0] + "\t" + chrPos[i][1] + "\t"
                             + alleles[i][0] + "\t" + alleles[i][1]);
            }
            writer.close();
          }
        }
        if (inval.size() > 0) {
          log.report(ext.getTime() + "]\tFound " + inval.size() + " invalid markers out of "
                     + mkrNames.length + ".");
          if (inval.size() < 10) {
            StringBuilder builder = new StringBuilder(inval.get(0));
            for (int i = 1; i < inval.size(); i++) {
              builder.append(", ").append(inval.get(i));
            }
            log.report(ext.getTime() + "]\tInvalid markers: " + builder.toString());
          } else {
            PrintWriter writer = Files.getAppropriateWriter(metaDir + "invalid.txt");
            for (Integer i : inval) {
              writer.println(mkrNames[i] + "\t" + chrPos[i][0] + "\t" + chrPos[i][1] + "\t"
                             + alleles[i][0] + "\t" + alleles[i][1]);
            }
            writer.close();
          }
        }

        final String outDelim = "\t";

        final AliasedFileColumn markerCol = StandardFileColumns.snp("MarkerName");
        final FileColumn<Byte> chrLinkedColumn = new AbstractFileColumn<Byte>("Chr") {

          @Override
          public Byte getValue(String[] line) throws ParseFailureException {
            String mkr = markerCol.getValue(line);
            if (mkr == null) {
              System.err.println("Error - couldn't identify marker in line: "
                                 + ArrayUtils.toStr(line));
              return null;
            }
            GenomicPosition gp = mkrsMeta.get(mkr);
            if (gp == null) {
              System.err.println("Error - GenomicPosition not found for marker: " + mkr);
              return null;
            }
            return gp.getChr();
          }
        };
        final FileColumn<Integer> posLinkedColumn = new AbstractFileColumn<Integer>("Position") {

          @Override
          public Integer getValue(String[] line) throws ParseFailureException {
            return mkrsMeta.get(markerCol.getValue(line)).getPosition();
          }
        };

        final AliasedFileColumn a1Column = StandardFileColumns.a1("a1");
        final AliasedFileColumn a2Column = StandardFileColumns.a2("a2");
        final IndexedFileColumn<Double> pColumn = StandardFileColumns.pVal("p");
        final FileColumn<String> remainingColumns = StandardFileColumns.allExcept(outDelim,
                                                                                  markerCol,
                                                                                  pColumn,
                                                                                  StandardFileColumns.chr("excludedChr"),
                                                                                  StandardFileColumns.pos("excludedPos"));
        final ColumnFilter mkrsBimFilter = new AbstractColumnFilter(markerCol) {

          @Override
          public boolean filter(DataLine values) {
            return mkrsMeta.containsKey(values.get(markerCol, null));
          }
        };
        final ColumnFilter mkrsAlleleFilter = new AbstractColumnFilter(markerCol, a1Column,
                                                                       a2Column) {

          @Override
          public boolean filter(DataLine values) {
            String mkr = values.get(markerCol, null);
            String[] dataAlleles = mkrAlleles.get(mkr);
            String[] metaAlleles = new String[] {values.get(a1Column, "NA"),
                                                 values.get(a2Column, "NA")};
            CONFIG config = StrandOps.determineStrandConfig(dataAlleles, metaAlleles);
            AlleleOrder alleleOrder = config.getAlleleOrder();
            if (config == CONFIG.AMBIGUOUS) {
              boolean a1 = dataAlleles[0].equalsIgnoreCase(metaAlleles[0])
                           || dataAlleles[0].equalsIgnoreCase(metaAlleles[1]);
              boolean a2 = dataAlleles[1].equalsIgnoreCase(metaAlleles[0])
                           || dataAlleles[1].equalsIgnoreCase(metaAlleles[1]);
              return a1 && a2;
            } else if (alleleOrder == AlleleOrder.SAME
                       || alleleOrder == AlleleOrder.OPPOSITE) return true;
            else {
              Joiner alleleJoiner = Joiner.on('/');
              log.reportError("Alleles in study (" + alleleJoiner.join(dataAlleles)
                              + ") do not match source alleles (" + alleleJoiner.join(metaAlleles)
                              + ") for " + mkr);
              return false;
            }
          }
        };
        final ColumnFilter pValThreshold = new DoubleFilter(pColumn, COMPARISON.LT,
                                                            constraintEntry.getValue().indexThreshold);
        final ColumnFilter pValMissing = new AbstractColumnFilter(pColumn) {

          @Override
          public boolean filter(DataLine values) {
            return !values.has(pColumn);
          }
        };
        try (FileParser crossFilterParser = FileParserFactory.setup(metaDir + dFile, markerCol,
                                                                    chrLinkedColumn,
                                                                    posLinkedColumn,
                                                                    remainingColumns)
                                                             .optionalColumns(pColumn)
                                                             .filter(mkrsBimFilter)
                                                             .filter(mkrsAlleleFilter)
                                                             .filter(ColumnFilters.or(pValMissing,
                                                                                      pValThreshold))
                                                             .build()) {
          int hitCount = crossFilterParser.parseToFile(crossFilterFile, outDelim);
          log.report(ext.getTime() + "]\tDropped "
                     + crossFilterParser.getFilteredCount(pValThreshold) + " snps for p-value < "
                     + constraintEntry.getValue().indexThreshold);
          log.report(ext.getTime() + "]\tDropped "
                     + crossFilterParser.getFilteredCount(mkrsAlleleFilter)
                     + " snps for mismatched alleles.");
          log.report(ext.getTime() + "]\tDropped "
                     + crossFilterParser.getFilteredCount(mkrsBimFilter)
                     + " snps for lacking data.");

          study.hitSnpCounts.get(constraintEntry.getKey()).put(dataFile, hitCount);
        }
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
        HashSet<String> hitMrkSet = new HashSet<>();
        for (String mkr : hitMarkers) {
          hitMrkSet.add(mkr);
        }

        final FileColumn<String> markerCol = StandardFileColumns.snp(MARKER_COL_NAME);
        final FileColumn<?> effectAlleleCol = StandardFileColumns.a1(EFFECT_ALLELE_COL_NAME);
        final FileColumn<?> nonEffectAlleleCol = StandardFileColumns.a2(NON_EFFECT_ALLELE_COL_NAME);
        final FileColumn<?> betaCol = StandardFileColumns.beta(BETA_COL_NAME);
        ColumnFilter hitMarkerFilter = new AbstractColumnFilter(markerCol) {

          @Override
          public boolean filter(DataLine values) {
            return hitMrkSet.contains(values.get(markerCol, null));
          }
        };
        try (FileParser crossFilterParser = FileParserFactory.setup(crossFilterFile, markerCol)
                                                             .optionalColumns(effectAlleleCol,
                                                                              nonEffectAlleleCol,
                                                                              betaCol)
                                                             .filter(hitMarkerFilter).build()) {
          String subsetFile = prefDir + "/subsetData_" + filePrefix.getKey() + ".xln";
          Map<String, HitMarker> dataList = crossFilterParser.parseToFileAndLoad(subsetFile, "\t",
                                                                                 true, markerCol)
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
    return new HitMarker(hitMarkerDataLine.get(StandardFileColumns.a1(EFFECT_ALLELE_COL_NAME),
                                               "NA"),
                         hitMarkerDataLine.get(StandardFileColumns.a2(NON_EFFECT_ALLELE_COL_NAME),
                                               "NA"),
                         hitMarkerDataLine.get(StandardFileColumns.beta(BETA_COL_NAME), null));
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
          if (hitMarker == null) {
            log.report(ext.getTime() + "]\tNo HitMarker data available for " + mkr);
            continue;
          }
          CONFIG config = determineAlleleConfig(alleles[m], hitMarker);
          AlleleOrder alleleOrder = config.getAlleleOrder();
          if (config == CONFIG.AMBIGUOUS) {
            alleleOrder = hitMarker.getEffectAllele()
                                   .compareToIgnoreCase(alleles[m][0]) == 0 ? StrandOps.AlleleOrder.SAME
                                                                            : StrandOps.AlleleOrder.OPPOSITE;
          } else if (!(alleleOrder.equals(AlleleOrder.SAME)
                       || alleleOrder.equals(AlleleOrder.OPPOSITE))) {
            Joiner alleleJoiner = Joiner.on('/');
            log.reportError("Alleles in study (" + alleleJoiner.join(alleles[m])
                            + ") do not match source alleles ("
                            + alleleJoiner.join(hitMarker.getEffectAllele() == null ? "NULL"
                                                                                    : hitMarker.getEffectAllele(),
                                                hitMarker.getNonEffectAllele() == null ? "NULL"
                                                                                       : hitMarker.getNonEffectAllele())
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
            if (alleleOrder.equals(StrandOps.AlleleOrder.OPPOSITE)) {
              cnt += isNaN ? 0 : 1;
              cnt2 += isNaN ? 0 : dosage;
              scoreSum += (isNaN ? matchedMarkerFreqs.get(mkr) : dosage) * beta;
            } else if (alleleOrder.equals(StrandOps.AlleleOrder.SAME)) {
              cnt += isNaN ? 0 : 1;
              cnt2 += isNaN ? 0 : (2.0 - dosage);
              scoreSum += (2.0 - (isNaN ? matchedMarkerFreqs.get(mkr) : dosage)) * beta;
            } else {
              throw new IllegalStateException("Mismatched alleles were not caught when cross-filtering");
            }
          }

          double mkrRatio = cnt / (double) markers.length;
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

  private static CONFIG determineAlleleConfig(String[] alleles, HitMarker hitMarker) {
    CONFIG config = StrandOps.determineStrandConfig(new String[] {hitMarker.getEffectAllele(),
                                                                  hitMarker.getNonEffectAllele()},
                                                    alleles);
    return config;
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

          HashMap<String, Double> scoreData = new HashMap<>();

          try (BufferedReader scoreReader = Files.getAppropriateReader(scoreFile)) {
            String line;
            while ((line = scoreReader.readLine()) != null) {
              String[] parts = line.split("\t");
              String score = parts[parts.length - 1];
              scoreData.put(parts[0] + "\t" + parts[1], Double.parseDouble(score));
            }
          }

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

          DosageData data = study.data.get(dataFile + "\t" + filePrefix.getKey());
          List<Trio> trios = data.getTrios();
          String[][] ids = data.getIds();
          TrioScoreTest trioTests = new TrioScoreTest();

          int sibs = 0;

          for (Trio trio : data.getTrios()) {
            String caseID = ids[trio.getChildIndex()][PSF.Plink.FAM_FID_INDEX] + "\t"
                            + ids[trio.getChildIndex()][PSF.Plink.FAM_IID_INDEX];
            if (ids[trio.getChildIndex()].length > PSF.Plink.FAM_AFF_INDEX
                && !PSF.Plink.affIsCase(ids[trio.getChildIndex()][PSF.Plink.FAM_AFF_INDEX])) {
              sibs++;
              continue;
            }
            String fatherID = ids[trio.getFatherIndex()][PSF.Plink.FAM_FID_INDEX] + "\t"
                              + ids[trio.getFatherIndex()][PSF.Plink.FAM_IID_INDEX];
            String motherID = ids[trio.getMotherIndex()][PSF.Plink.FAM_FID_INDEX] + "\t"
                              + ids[trio.getMotherIndex()][PSF.Plink.FAM_IID_INDEX];
            if (scoreData.containsKey(caseID) && scoreData.containsKey(fatherID)
                && scoreData.containsKey(motherID)) {
              double caseScore = scoreData.get(caseID);
              double fatherScore = scoreData.get(fatherID);
              double motherScore = scoreData.get(motherID);
              trioTests.addScore(caseScore, fatherScore, motherScore);
            }
            TrioScoreTest.Results trioTestResults = trioTests.runTests();
            if (trioTestResults != null) study.trioScoreTests.put(filePrefix.getKey(), dataFile,
                                                                  trioTestResults);
          }

          if (sibs > 0) {
            log.reportTimeWarning(sibs + " trios (of " + data.getTrios().size()
                                  + " total trios) were excluded from paired score analyses because the children were not coded as cases in the fam data.");
          }

          for (int i = 0; i < study.phenoFiles.size(); i++) {
            PhenoData pd = study.phenoData.get(study.phenoFiles.get(i));
            ArrayList<Double> depData = new ArrayList<>();
            ArrayList<double[]> baselineIndeps = new ArrayList<>();
            ArrayList<double[]> indepData = new ArrayList<>();

            for (java.util.Map.Entry<String, PhenoIndiv> indiv : pd.indivs.entrySet()) {
              if (scoreData.containsKey(indiv.getKey())) {
                PhenoIndiv pdi = pd.indivs.get(indiv.getKey());
                double[] baseData = new double[pd.covars.size()];
                double[] covarData = new double[pd.covars.size() + 1];
                covarData[0] = scoreData.get(indiv.getKey());
                boolean validCovars = true;
                for (int k = 1; k < pd.covars.size() + 1; k++) {
                  Double d = pdi.getCovars().get(pd.covars.get(k - 1));
                  if (d == null) {
                    log.reportError("Covar value missing for individual: " + indiv.getKey() + " | "
                                    + pd.covars.get(k - 1));
                    validCovars = false;
                  } else {
                    baseData[k - 1] = d.doubleValue();
                    covarData[k] = d.doubleValue();
                  }
                }
                if (validCovars) {
                  depData.add(pdi.getDepvar());
                  baselineIndeps.add(baseData);
                  indepData.add(covarData);
                }
              }
            }

            double[][] baseCovars = new double[baselineIndeps.size()][];
            double[][] covars = new double[indepData.size()][];
            for (int k = 0; k < covars.length; k++) {
              covars[k] = indepData.get(k);
              baseCovars[k] = baselineIndeps.get(k);
            }
            int cases = 0;
            int controls = 0;
            Multiset<Double> phenos = ImmutableMultiset.copyOf(depData);
            if (phenos.entrySet().size() == 2) {
              if (phenos.elementSet()
                        .equals(ImmutableSet.of(Double.valueOf(0.0), Double.valueOf(1.0)))) {
                cases = phenos.count(Double.valueOf(1.0));
                controls = phenos.count(Double.valueOf(0.0));
              } else if (phenos.elementSet()
                               .equals(ImmutableSet.of(Double.valueOf(1.0), Double.valueOf(2.0)))) {
                cases = phenos.count(Double.valueOf(2.0));
                controls = phenos.count(Double.valueOf(1.0));
              } else {
                log.reportError("Unrecognized case/control designations, cannot count cases and controls");
              }
            }
            RegressionModel baseModel = RegressionModel.determineAppropriate(Doubles.toArray(depData),
                                                                             baseCovars, false,
                                                                             true);
            RegressionModel model = RegressionModel.determineAppropriate(Doubles.toArray(depData),
                                                                         covars, false, true);

            RegressionResult.Builder rr = RegressionResult.builder();
            rr.setNum(depData.size());
            rr.setnCases(cases);
            rr.setnControls(controls);
            if (model.analysisFailed()) {
              rr.setBaseRSq(baseModel.analysisFailed() ? Double.NaN : baseModel.getRsquare());
              rr.setBeta(Double.NaN);
              rr.setSe(Double.NaN);
              rr.setRsq(Double.NaN);
              rr.setPval(Double.NaN);
              rr.setLogistic(model.isLogistic());
            } else {
              int ind = -1;
              for (int l = 0; l < model.getVarNames().length; l++) {
                if ("Indep 1".equals(model.getVarNames()[l])) {
                  ind = l;
                  break;
                }
              }
              if (ind == -1) {
                rr.setBaseRSq(baseModel.analysisFailed() ? Double.NaN : baseModel.getRsquare());
                rr.setBeta(Double.NaN);
                rr.setSe(Double.NaN);
                rr.setRsq(model.getRsquare());
                rr.setPval(Double.NaN);
                rr.setLogistic(model.isLogistic());
                rr.setStats(Double.NaN);
              } else {
                rr.setBaseRSq(baseModel.analysisFailed() ? Double.NaN : baseModel.getRsquare());
                rr.setBeta(model.getBetas()[ind]);
                rr.setSe(model.getSEofBs()[ind]);
                rr.setRsq(model.getRsquare());
                rr.setPval(model.getSigs()[ind]);
                rr.setLogistic(model.isLogistic());
                rr.setStats(model.getStats()[ind]);
              }
            }
            study.regressions.get(filePrefix.getKey()).get(dataFile).put(pd.phenoName, rr.build());

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
    String resFile = metaDir + "results.xln";
    log.report(ext.getTime() + "]\tWriting regression results... [ --> " + resFile + "]");
    try (PrintWriter writer = Files.getAppropriateWriter(resFile)) {
      writer.println(REGRESSION_HEADER);

      for (Study study : studies) {
        for (String dFile : metaFiles) {
          String dataFile = ext.rootOf(dFile, false);
          for (java.util.Map.Entry<String, Constraint> filePrefix : analysisConstraints.entrySet()) {
            HashMap<String, RegressionResult> phenoResults = study.regressions.get(filePrefix.getKey())
                                                                              .get(dataFile);
            TrioScoreTest.Results trioTestResults = study.trioScoreTests.get(filePrefix.getKey(),
                                                                             dataFile);
            final String pairedStatTrios;
            final String pairedT;
            final String wilcoxon;
            if (trioTestResults == null) {
              pairedStatTrios = "0";
              pairedT = "N/A";
              wilcoxon = "N/A";
            } else {
              pairedStatTrios = String.valueOf(trioTestResults.getTrioCount());
              pairedT = String.valueOf(trioTestResults.pairedTPVal);
              wilcoxon = String.valueOf(trioTestResults.wilcoxonPVal);
            }

            String resultPrefix = new StringJoiner("\t").add(study.studyName).add(dataFile)
                                                        .add(ext.formSciNot(filePrefix.getValue().indexThreshold,
                                                                            5, false))
                                                        .toString();

            String middle = new StringJoiner("\t").add(pairedT).add(wilcoxon).add(pairedStatTrios)
                                                  .add(String.valueOf(dataCounts.get(dFile)
                                                                                .get(filePrefix.getKey())))
                                                  .add(String.valueOf(study.hitWindowCnts.get(filePrefix.getKey())
                                                                                         .get(dataFile)))
                                                  .add(String.valueOf(study.hitSnpCounts.get(filePrefix.getKey())
                                                                                        .get(dataFile)))
                                                  .add(String.valueOf(study.scores.get(filePrefix.getKey())
                                                                                  .get(dataFile)[0]))
                                                  .add(String.valueOf(study.scores.get(filePrefix.getKey())
                                                                                  .get(dataFile)[1]))
                                                  .toString();
            if (study.phenoFiles.isEmpty()) {
              writeSingleResult("--", RegressionResult.dummy(), resultPrefix, middle, writer);
            } else {
              for (String pheno : study.phenoFiles) {
                RegressionResult rr = phenoResults.get(pheno);
                if (rr == null) {
                  rr = RegressionResult.dummy();
                }
                writeSingleResult(pheno, rr, resultPrefix, middle, writer);
              }
            }
          }
        }
      }
    }
  }

  private void writeSingleResult(String pheno, RegressionResult rr, String resultPrefix,
                                 String middle, PrintWriter writer) {
    String pvalExcl = rr.getNum() == 0 ? "."
                                       : (rr.isLogistic() ? "=(1-(NORM.S.DIST(ABS(" + rr.getBeta()
                                                            + "/" + rr.getSe() + "),TRUE)))*2"
                                                          : "=TDIST(" + Math.abs(rr.getStats())
                                                            + "," + rr.getNum() + ",2)");

    String line = Joiner.on("\t")
                        .join(ImmutableList.builder().add(resultPrefix).add(pheno)
                                           .add(rr.getBaseRSq()).add(rr.getRsq())
                                           .add((Double.isNaN(rr.getRsq()) ? Double.NaN
                                                                           : (Double.isNaN(rr.getBaseRSq()) ? rr.getRsq()
                                                                                                            : (new BigDecimal(rr.getRsq()
                                                                                                                              + "")).subtract(new BigDecimal(rr.getBaseRSq()
                                                                                                                                                             + "")))))
                                           .add(rr.getPval()).add(pvalExcl).add(rr.getBeta())
                                           .add(rr.getSe()).add(rr.getNum()).add(rr.getnCases())
                                           .add(rr.getnControls()).add(middle).build());

    writer.println(line);

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
    preprocessDataFiles(preprocess.toArray(new String[preprocess.size()]), POPULATION.ALL,
                        GenomeBuild.HG19, log);
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

    String workDir = null;
    String logFile = "GeneScorePipeline.log";

    float[] iT = new float[] {DEFAULT_INDEX_THRESHOLD};
    int[] mZ = new int[] {DEFAULT_WINDOW_MIN_SIZE_PER_SIDE};
    float[] wT = new float[] {DEFAULT_WINDOW_EXTENSION_THRESHOLD};
    double mT = DEFAULT_MIN_MISS_THRESH;

    String[] processList = null;
    boolean process = false;
    POPULATION pop = POPULATION.ALL;
    GenomeBuild build = GenomeBuild.HG19;
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

    boolean fail = false;
    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith(ARG_WORKING_DIR)) {
        workDir = arg.split("=")[1];
      } else if (arg.startsWith("process=")) {
        processList = arg.split("=")[1].split(",");
        process = true;
      } else if (arg.startsWith("pop=")) {
        pop = POPULATION.valueOf(arg.split("=")[1]);
      } else if (arg.startsWith("build=")) {
        build = GenomeBuild.valueOf(arg.split("=")[1]);
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
      } else if (arg.startsWith(ARG_MISS_THRESH)) {
        mT = ext.parseDoubleArg(arg);
      } else if (arg.startsWith("log=")) {
        logFile = arg.split("=")[1];
      } else {
        fail = true;
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (fail || args.length == 0) {
      System.err.println(usage);
      System.exit(1);
    }

    Logger log = new Logger(logFile, true);

    if (process) {
      preprocessDataFiles(processList, pop, build, log);
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
