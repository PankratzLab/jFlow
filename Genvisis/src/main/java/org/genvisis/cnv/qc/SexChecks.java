// -Xms1024M -Xmx1024M
package org.genvisis.cnv.qc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NavigableMap;
import java.util.NavigableSet;
import java.util.Set;
import java.util.StringJoiner;
import java.util.Vector;
import org.apache.commons.math3.stat.inference.TTest;
import org.genvisis.cnv.analysis.MosaicismDetect;
import org.genvisis.cnv.analysis.MosaicismDetect.MosaicBuilder;
import org.genvisis.cnv.filesys.ClusterFilterCollection;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.MDL;
import org.genvisis.cnv.manage.MarkerDataLoader;
import org.genvisis.cnv.var.MosaicRegion;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.DoubleVector;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.PSF;
import org.genvisis.common.ProgressMonitor.DISPLAY_MODE;
import org.genvisis.common.ext;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.stats.LogisticRegression;
import org.genvisis.stats.Ttest;
import com.google.common.base.Functions;
import com.google.common.base.Joiner;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ArrayTable;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.google.common.collect.Table;
import com.google.common.math.Quantiles;
import com.google.common.primitives.Doubles;

public class SexChecks {

  public enum BinarySex {
    UNKNOWN, MALE, FEMALE
  }

  public enum EstimatedSex {
    UNKNOWN(BinarySex.UNKNOWN, 0, "Unknown", ""),
    MALE(BinarySex.MALE, 1, "Male", "XY"),
    FEMALE(BinarySex.FEMALE, 2, "Female", "XX"),
    KLINEFELTER(BinarySex.MALE, 2, "Klinefelter", "XXY"),
    UPD_KLINEFELTER(BinarySex.MALE, 2, "UPD Klinefelter", "XXY"),
    MOSAIC_KLINEFELTER(BinarySex.MALE, 2, "Mosaic Klinefelter", "XXY"),
    TRIPLE_X(BinarySex.FEMALE, 3, "Triple X", "XXX"),
    MOSAIC_TRIPLE_X(BinarySex.FEMALE, 3, "Mosaic Triple X", "XXX"),
    TURNER(BinarySex.FEMALE, 1, "Turner", "X"),
    MOSAIC_TURNER(BinarySex.FEMALE, 1, "Mosaic Turner", "X");

    private static final Map<String, EstimatedSex> DESC_MAPPING = Arrays.stream(values())
                                                                        .collect(ImmutableMap.toImmutableMap(EstimatedSex::toString,
                                                                                                             Functions.identity()));

    private static final Map<String, EstimatedSex> HEADER_MAPPING = Arrays.stream(values())
                                                                          .collect(ImmutableMap.toImmutableMap(EstimatedSex::headerString,
                                                                                                               Functions.identity()));

    private final BinarySex binarySex;
    private final int chrXCount;
    private final String desc;
    private final String karyotype;

    private EstimatedSex(BinarySex baseSex, int chrXCount, String desc, String karyotype) {
      this.binarySex = baseSex;
      this.chrXCount = chrXCount;
      this.desc = desc;
      this.karyotype = karyotype;
    }

    /**
     * @return the binarySex
     */
    public BinarySex getBinarySex() {
      return binarySex;
    }

    /**
     * @return the chrXCount
     */
    public int getChrXCount() {
      return chrXCount;
    }

    /**
     * @return the karyotype
     */
    public String getKaryotype() {
      return karyotype;
    }

    @Override
    public String toString() {
      return desc;
    }

    public String headerString() {
      return ordinal() + "=" + toString();
    }

    public static EstimatedSex parse(String desc) {
      return DESC_MAPPING.get(desc);
    }

    public static EstimatedSex parseHeaderString(String headerString) {
      return HEADER_MAPPING.get(headerString);
    }

  }

  public static final String EST_SEX_HEADER = generateEstSexHeader();

  public static final String[] SEX_HEADER = {"Sample", "FID", "IID", "Sex", EST_SEX_HEADER, "Note",
                                             "Check", "Excluded", "Median X e^LRR",
                                             "Median Y e^LRR", "e^LRR Ratio Y:X",
                                             "% X Heterozygote Calls", "% X BAF 0.15-0.85",
                                             "Median X LRR", "Median Y LRR"};

  private static final float XY_ELRR_RATIO_MIN_SEED_MALE = 1.0f;
  private static final float XY_ELRR_RATIO_MAX_SEED_MALE = 1.5f;
  private static final float XY_ELRR_RATIO_MAX_SEED_FEMALE = 0.5f;
  private static final float NUM_SD_FOR_HET_OUTLIERS = 4.0f;
  private static final float NUM_SD_FOR_MALE_X_OUTLIERS = 1.5f;
  private static final float NUM_SD_FOR_MALE_X_FULL_ANEUPLOIDY = 4.0f;
  private static final float NUM_SD_FOR_FEMALE_X_OUTLIERS = 1.5f;
  private static final float NUM_SD_FOR_FEMALE_X_FULL_ANEUPLOIDY = 4.0f;
  private static final float MAX_SD_FOR_Y_OUTLIERS = 5.0f;
  private static final double SEX_DISCRIMINATING_BASE_P_THRESHOLD = 0.001; // This will be
                                                                          // bonferroni corrected
                                                                          // for number of markers
                                                                          // checked
  private static final double MOSAIC_F_CERTAINTY_THRESHOLD = 0.2;
  private static final double MOSAIC_COVERAGE_CERTAINTY_THRESHOLD = 0.8;
  private static final double MOSAIC_COVERAGE_ABSOLUTE_THRESHOLD = 0.5;

  private final Project proj;
  private final Logger log;
  private final MarkerDetailSet markerSet;
  private final List<String> sampleNames;
  private Set<String> qcPassedSamples;

  private Set<Marker> xMarkers;
  private Set<Marker> yMarkers;

  private final Set<Marker> xUseMarkers;
  private final Set<Marker> yUseMarkers;

  // TODO Talk to Nathan about switching from r ratios to e^lrr ratios

  private Map<String, Double> elrrMedX;
  private Map<String, Double> elrrMedY;

  private Table<Marker, String, Double> lrrsX;
  private Table<Marker, String, Double> lrrsY;

  private Table<Marker, String, Byte> genotypesX;
  private Table<Marker, String, Double> bafsX;

  private Set<String> seedMales;
  private Set<String> seedFemales;

  private final Map<String, Double> lrrMedX;
  private final Map<String, Double> lrrMedY;

  private final Map<String, Double> pctXHets;
  private final Map<String, Double> pctXBaf15_85;

  private Map<String, EstimatedSex> sexes;
  private Set<String> uncertains;
  private Multimap<String, String> notes;

  private SexChecks(Project proj, boolean appendToSampleData,
                    String nonCrossHybridizingMarkersFile) {
    long startTime = new Date().getTime();
    this.proj = proj;
    log = proj.getLog();

    markerSet = proj.getMarkerSet();
    sampleNames = Collections.unmodifiableList(Arrays.asList(proj.getSamples()));

    qcPassedSamples = generateQCPassedSamples();

    PSF.checkInterrupted();

    log.report("Finding Sex Chromosome Markers...");
    generateMarkerLists(nonCrossHybridizingMarkersFile);

    log.report("Loading Sex Chromosome Marker Data...");
    gatherMarkerStats();

    log.report("Determining Samples to seed Sex Checks...");
    generateSeedSexLists();
    log.report("Found " + seedMales.size() + " obvious males");
    log.report("Found " + seedFemales.size() + " obvious females");
    log.report("Seeding sex checks using these " + (seedMales.size() + seedFemales.size())
               + " samples (of " + qcPassedSamples.size() + " QC passed samples)");

    log.report("Scanning for markers that express differently by sex...");
    PSF.checkInterrupted();

    xUseMarkers = sexDiscriminatingXMarkers();
    yUseMarkers = sexDiscriminatingYMarkers();

    log.report("Found " + xUseMarkers.size() + " sex differentiating markers out of "
               + xMarkers.size() + " X chromosome markers");
    log.report("Found " + yUseMarkers.size() + " sex differentiating markers out of "
               + yMarkers.size() + " Y chromosome markers");
    PSF.checkInterrupted();

    log.report("Calculating median sample LRR for identified X and Y chromosome markers");
    lrrMedX = calcMedianLRRs(lrrsX, xUseMarkers);
    lrrMedY = calcMedianLRRs(lrrsY, yUseMarkers);
    PSF.checkInterrupted();

    log.report("Calculating sample counts of heterozygote calls for identified X chromosome markers...");
    pctXHets = calcPctHets(genotypesX, xUseMarkers);
    pctXBaf15_85 = calcPctBaf15_85(bafsX, xUseMarkers);
    PSF.checkInterrupted();

    log.report("Estimating sex for each sample...");
    estimateSexes();
    PSF.checkInterrupted();

    log.report("Writing outputs...");
    writeToFile(appendToSampleData);
    log.report("Finished estimating sample sexes in " + ext.getTimeElapsed(startTime));
  }

  private static String generateEstSexHeader() {
    String header = "Estimated Sex";
    for (EstimatedSex sex : EstimatedSex.values()) {
      header += ";" + sex.headerString();
    }
    return header;
  }

  private Set<String> generateQCPassedSamples() {
    if (!Files.exists(proj.SAMPLE_QC_FILENAME.getValue())) {
      log.reportTimeWarning("LRR SD values are required for filtering samples to use in Sex Checks. Sample QC will be run now...");
      int numThreads = proj.NUM_THREADS.getValue();
      LrrSd.init(proj, null, null, numThreads, false);
    }

    Set<String> passedSamples = LrrSd.qcPassedSampleSet(proj);
    if (passedSamples.isEmpty()) {
      throw new IllegalStateException("Sex Checks failed: no LRR SD-filtered samples were found. Please verify that Sample QC was run successfully.");
    }
    return passedSamples;
  }

  private void generateMarkerLists(String nonCrossHybridizingMarkersFile) {
    NavigableMap<Byte, NavigableSet<Marker>> chrMap = markerSet.getChrMap();
    xMarkers = chrMap.get((byte) 23);
    yMarkers = chrMap.get((byte) 24);
    if (nonCrossHybridizingMarkersFile == null) {
      log.reportError("No file of markers that do not cross hybridize was provided, all X and Y chromosome markers will be used to determine sex baselines");
    } else {
      log.report("Using " + nonCrossHybridizingMarkersFile
                 + " to identify markers that do not cross hybridize");
      Map<String, Marker> markerNameMap = markerSet.getMarkerNameMap();
      Set<Marker> nonCrossHybridizingMarkers = HashVec.loadFileToHashSet(nonCrossHybridizingMarkersFile,
                                                                         false)
                                                      .stream().map(markerNameMap::get)
                                                      .collect(ImmutableSet.toImmutableSet());
      xMarkers = Sets.intersection(xMarkers, nonCrossHybridizingMarkers);
      yMarkers = Sets.intersection(yMarkers, nonCrossHybridizingMarkers);
    }
  }

  private void gatherMarkerStats() {
    gatherXMarkerStats();
    gatherYMarkerStats();

  }

  private void gatherXMarkerStats() {
    ClusterFilterCollection clusterFilters = proj.getClusterFilterCollection();
    float gcThreshold = proj.GC_THRESHOLD.getValue().floatValue();
    List<Marker> xMarkerList = ImmutableList.copyOf(xMarkers);
    MDL mdl = new MDL(proj, xMarkerList, Math.max(proj.NUM_THREADS.getValue() - 1, 1), 100);

    lrrsX = ArrayTable.create(xMarkerList, sampleNames);
    genotypesX = ArrayTable.create(xMarkerList, sampleNames);
    bafsX = ArrayTable.create(xMarkerList, sampleNames);
    Table<String, Marker, Double> elrrs = HashBasedTable.create(sampleNames.size(),
                                                                xMarkerList.size());

    for (int m = 0; mdl.hasNext(); m++) {
      MarkerData markerData = mdl.next();
      Marker marker = xMarkerList.get(m);
      if (!markerData.getMarkerName().equals(marker.getName())) {
        throw new IllegalStateException("MarkerData name does not match Marker name at index " + m
                                        + ". Expected " + marker.getName() + ", found "
                                        + markerData.getMarkerName());
      }
      float[] lrrs = markerData.getLRRs();
      byte[] genos = markerData.getAbGenotypesAfterFilters(clusterFilters,
                                                           markerData.getMarkerName(), gcThreshold,
                                                           log);
      float[] bafs = markerData.getBAFs();
      for (int s = 0; s < sampleNames.size(); s++) {
        String sample = sampleNames.get(s);
        double lrr = lrrs[s];
        double baf = bafs[s];
        lrrsX.put(marker, sample, lrr);
        genotypesX.put(marker, sample, genos[s]);
        bafsX.put(marker, sample, baf);
        double elrr = Math.pow(Math.E, lrr);
        if (Double.isFinite(elrr)) elrrs.put(sample, marker, elrr);
      }
    }
    mdl.shutdown();
    elrrMedX = elrrs.rowMap().entrySet().stream()
                    .collect(ImmutableMap.toImmutableMap(Entry::getKey,
                                                         e -> Quantiles.median()
                                                                       .compute(e.getValue()
                                                                                 .values())));
  }

  private void gatherYMarkerStats() {
    ImmutableList<Marker> yMarkerList = ImmutableList.copyOf(yMarkers);
    MDL mdl = new MDL(proj, yMarkerList, 1, 100);

    lrrsY = ArrayTable.create(yMarkerList, sampleNames);
    Table<String, Marker, Double> elrrs = HashBasedTable.create(sampleNames.size(),
                                                                yMarkerList.size());
    for (int m = 0; mdl.hasNext(); m++) {
      MarkerData markerData = mdl.next();
      Marker marker = yMarkerList.get(m);
      if (!markerData.getMarkerName().equals(marker.getName())) {
        throw new IllegalStateException("MarkerData name does not match Marker name at index " + m
                                        + ". Expected " + marker.getName() + ", found "
                                        + markerData.getMarkerName());
      }
      float[] lrrs = markerData.getLRRs();
      for (int s = 0; s < sampleNames.size(); s++) {
        String sample = sampleNames.get(s);
        double lrr = lrrs[s];
        lrrsY.put(marker, sample, lrr);
        double elrr = Math.pow(Math.E, lrr);
        if (Double.isFinite(elrr)) elrrs.put(sample, marker, elrr);
      }
    }
    mdl.shutdown();
    elrrMedY = elrrs.rowMap().entrySet().stream()
                    .collect(ImmutableMap.toImmutableMap(Entry::getKey,
                                                         e -> Quantiles.median()
                                                                       .compute(e.getValue()
                                                                                 .values())));
  }

  private void generateSeedSexLists() {
    ImmutableSet.Builder<String> maleBuilder = ImmutableSet.builder();
    ImmutableSet.Builder<String> femaleBuilder = ImmutableSet.builder();

    for (String sample : qcPassedSamples) {
      double elrrRatio = elrrMedY.get(sample) / elrrMedX.get(sample);
      if (elrrRatio > XY_ELRR_RATIO_MIN_SEED_MALE && elrrRatio < XY_ELRR_RATIO_MAX_SEED_MALE) {
        maleBuilder.add(sample);
      } else if (elrrRatio < XY_ELRR_RATIO_MAX_SEED_FEMALE) {
        femaleBuilder.add(sample);
      }
    }

    seedMales = maleBuilder.build();
    seedFemales = femaleBuilder.build();
  }

  private Set<Marker> sexDiscriminatingXMarkers() {
    ImmutableSet.Builder<Marker> discriminatingMarkers = ImmutableSet.builder();
    TTest tTest = new TTest();
    for (Marker marker : xMarkers) {
      Map<String, Double> maleLrrs = Maps.filterValues(Maps.filterKeys(lrrsX.row(marker),
                                                                       seedMales::contains),
                                                       Double::isFinite);
      Map<String, Double> femaleLrrs = Maps.filterValues(Maps.filterKeys(lrrsX.row(marker),
                                                                         seedFemales::contains),
                                                         Double::isFinite);
      if (maleLrrs.size() > 1 && femaleLrrs.size() > 2
          && ArrayUtils.mean(femaleLrrs.values()) > ArrayUtils.mean(maleLrrs.values())) {
        double pVal = tTest.tTest(Doubles.toArray(maleLrrs.values()),
                                  Doubles.toArray(femaleLrrs.values()));
        if (pVal < SEX_DISCRIMINATING_BASE_P_THRESHOLD / xMarkers.size()) {
          discriminatingMarkers.add(marker);
        }
      }
    }
    return discriminatingMarkers.build();
  }

  private Set<Marker> sexDiscriminatingYMarkers() {
    ImmutableSet.Builder<Marker> discriminatingMarkers = ImmutableSet.builder();
    TTest tTest = new TTest();
    for (Marker marker : yMarkers) {
      Map<String, Double> maleLrrs = Maps.filterValues(Maps.filterKeys(lrrsY.row(marker),
                                                                       seedMales::contains),
                                                       Double::isFinite);
      Map<String, Double> femaleLrrs = Maps.filterValues(Maps.filterKeys(lrrsY.row(marker),
                                                                         seedFemales::contains),
                                                         Double::isFinite);
      if (maleLrrs.size() > 1 && femaleLrrs.size() > 1
          && ArrayUtils.mean(maleLrrs.values()) > ArrayUtils.mean(femaleLrrs.values())) {
        double pVal = tTest.tTest(Doubles.toArray(maleLrrs.values()),
                                  Doubles.toArray(femaleLrrs.values()));
        if (pVal < SEX_DISCRIMINATING_BASE_P_THRESHOLD / yMarkers.size()) {
          discriminatingMarkers.add(marker);
        }
      }
    }
    return discriminatingMarkers.build();
  }

  private Map<String, Double> calcMedianLRRs(Table<Marker, String, Double> lrrs,
                                             Set<Marker> useMarkers) {
    ImmutableMap.Builder<String, Double> medianLRRs = ImmutableMap.builderWithExpectedSize(sampleNames.size());
    for (String sample : sampleNames) {
      medianLRRs.put(sample,
                     Quantiles.median()
                              .compute(Maps.filterEntries(lrrs.column(sample),
                                                          e -> useMarkers.contains(e.getKey())
                                                               && Double.isFinite(e.getValue()))
                                           .values()));
    }
    return medianLRRs.build();
  }

  private Map<String, Double> calcPctHets(Table<Marker, String, Byte> genotypes,
                                          Set<Marker> useMarkers) {
    return genotypes.columnMap().entrySet().stream()
                    .collect(ImmutableMap.toImmutableMap(Entry::getKey,
                                                         e -> pctHet(Maps.filterKeys(e.getValue(),
                                                                                     useMarkers::contains)
                                                                         .values())));

  }

  private static double pctHet(Collection<Byte> genos) {
    return genos.stream().filter(g -> g == (byte) 1).count() / (double) genos.size();
  }

  private Map<String, Double> calcPctBaf15_85(Table<Marker, String, Double> bafs,
                                              Set<Marker> useMarkers) {
    return bafs.columnMap().entrySet().stream()
               .collect(ImmutableMap.toImmutableMap(Entry::getKey,
                                                    e -> pctBaf15_85(Maps.filterKeys(e.getValue(),
                                                                                     useMarkers::contains)
                                                                         .values())));

  }

  private static double pctBaf15_85(Collection<Double> bafs) {
    return bafs.stream().mapToDouble(Double::doubleValue).filter(baf -> baf > 0.15 && baf < 0.85)
               .count()
           / (double) bafs.size();
  }

  private void estimateSexes() {
    Collection<Double> maleMedLRRsX = Maps.filterEntries(lrrMedX,
                                                         e -> seedMales.contains(e.getKey())
                                                              && Double.isFinite(e.getValue()))
                                          .values();
    Collection<Double> femaleMedLRRsX = Maps.filterEntries(lrrMedX,
                                                           e -> seedFemales.contains(e.getKey())
                                                                && Double.isFinite(e.getValue()))
                                            .values();

    double maleMeanX = ArrayUtils.mean(maleMedLRRsX);
    double maleStdDevX = ArrayUtils.stdev(maleMedLRRsX);
    double femaleMeanX = ArrayUtils.mean(femaleMedLRRsX);
    double femaleStdDevX = ArrayUtils.stdev(femaleMedLRRsX);

    Collection<Double> malePctXHets = Maps.filterEntries(pctXHets,
                                                         e -> seedMales.contains(e.getKey())
                                                              && Double.isFinite(e.getValue()))
                                          .values();

    double maleMeanPctXHets = ArrayUtils.mean(malePctXHets);
    double maleStdDevPctXHets = ArrayUtils.stdev(malePctXHets);

    Collection<Double> femalePctXBaf15_85 = Maps.filterEntries(pctXBaf15_85,
                                                               e -> seedFemales.contains(e.getKey())
                                                                    && Double.isFinite(e.getValue()))
                                                .values();

    double femaleMeanPctXBaf15_85 = ArrayUtils.mean(femalePctXBaf15_85);
    double femaleStdDevPctXBaf15_85 = ArrayUtils.stdev(femalePctXBaf15_85);

    Collection<Double> maleMedLRRsY = Maps.filterEntries(lrrMedY,
                                                         e -> seedMales.contains(e.getKey())
                                                              && Double.isFinite(e.getValue()))
                                          .values();
    Collection<Double> femaleMedLRRsY = Maps.filterEntries(lrrMedY,
                                                           e -> seedFemales.contains(e.getKey())
                                                                && Double.isFinite(e.getValue()))
                                            .values();

    double maleMeanY = ArrayUtils.mean(maleMedLRRsY);
    double maleStdDevY = ArrayUtils.stdev(maleMedLRRsY);
    double femaleMeanY = ArrayUtils.mean(femaleMedLRRsY);
    double femaleStdDevY = ArrayUtils.stdev(femaleMedLRRsY);

    sexes = Maps.newHashMapWithExpectedSize(sampleNames.size());
    uncertains = Sets.newHashSet();
    notes = ArrayListMultimap.create();

    Set<Marker> mosaicismCheckUse = mosaicismUse();

    int sdForYOutliers = 0;
    while (sdForYOutliers < MAX_SD_FOR_Y_OUTLIERS) {
      if ((maleMeanY
           - (sdForYOutliers + 1) * maleStdDevY) > (femaleMeanY
                                                    + (sdForYOutliers + 1) * femaleStdDevY)) {
        sdForYOutliers++;
      } else {
        break;
      }
    }
    log.report("Using " + sdForYOutliers
               + " standard deviations from mean male and female Y LRRs to define sex clusters");
    double maleFloorY = maleMeanY - sdForYOutliers * maleStdDevY;
    log.report("Male mean Y LRR:    " + maleMeanY);
    log.report("Male Std Dev Y LRR: " + maleStdDevY);
    log.report("Male Y LRR Floor:   " + maleFloorY);
    double femaleCeilingY = femaleMeanY + sdForYOutliers * femaleStdDevY;
    log.report("Female mean Y LRR:    " + femaleMeanY);
    log.report("Female SD Y LRR:      " + femaleStdDevY);
    log.report("Female Y LRR Ceiling: " + femaleCeilingY);

    String taskName = "SexEstimation";

    proj.getProgressMonitor().beginDeterminateTask(taskName, "Estimating Sexes", sampleNames.size(),
                                                   DISPLAY_MODE.GUI_AND_CONSOLE);

    for (String sample : sampleNames) {
      boolean male = false;
      boolean female = false;
      double lrrY = lrrMedY.get(sample);
      if (lrrY > maleFloorY) {
        male = true;
      } else if (lrrY < femaleCeilingY) {
        female = true;
      } else {
        uncertains.add(sample);
        notes.put(sample, "Median Y LRR (" + ext.formDeci(lrrY, 4)
                          + ") is outside of both male and female acceptance intervals");
        if (seedMales.contains(sample)) {
          male = true;
        } else if (seedFemales.contains(sample)) {
          female = true;
        }
      }

      if (male) {
        if (seedFemales.contains(sample)) {
          uncertains.add(sample);
          notes.put(sample,
                    "Ratio of Median X e^LRR to Median Y e^LRR ("
                            + ext.formDeci(elrrMedY.get(sample) / elrrMedX.get(sample), 4)
                            + ") indicated female");
        } else if (!seedMales.contains(sample)) {
          notes.put(sample,
                    "Ratio of Median X e^LRR to Median Y e^LRR ("
                            + ext.formDeci(elrrMedY.get(sample) / elrrMedX.get(sample), 4)
                            + ") outlier");
        }
        if (pctXHets.get(sample) > (maleMeanPctXHets + NUM_SD_FOR_HET_OUTLIERS * maleStdDevPctXHets)
            && lrrMedX.get(sample) > (maleMeanX + NUM_SD_FOR_MALE_X_OUTLIERS * maleStdDevX)) {
          if (lrrMedX.get(sample) < (maleMeanX + NUM_SD_FOR_MALE_X_FULL_ANEUPLOIDY * maleStdDevX)) {
            uncertains.add(sample);
            notes.put(sample,
                      "Median X LRR (" + ext.formDeci(lrrMedX.get(sample), 4)
                              + ") not elevated enough to call Klinefelter without X heterozygosity ("
                              + ext.formPercent(pctXHets.get(sample), 4) + ")");
          }
          if (checkXMosaicism(sample, mosaicismCheckUse)) {
            sexes.put(sample, EstimatedSex.MOSAIC_KLINEFELTER);
          } else {
            sexes.put(sample, EstimatedSex.KLINEFELTER);
          }
        } else if (lrrMedX.get(sample) > (maleMeanX
                                          + NUM_SD_FOR_MALE_X_FULL_ANEUPLOIDY * maleStdDevX)) {
          sexes.put(sample, EstimatedSex.UPD_KLINEFELTER);
        } else {
          if (pctXHets.get(sample) > (maleMeanPctXHets
                                      + NUM_SD_FOR_HET_OUTLIERS * maleStdDevPctXHets)) {
            uncertains.add(sample);
            notes.put(sample,
                      "X heterozygosity (" + ext.formPercent(pctXHets.get(sample), 4)
                              + ") suggests Klinefelter but Median X LRR ("
                              + ext.formDeci(lrrMedX.get(sample), 4) + ") is not elevated");
          }
          sexes.put(sample, EstimatedSex.MALE);
        }
      } else if (female) {
        if (seedMales.contains(sample)) {
          uncertains.add(sample);
          notes.put(sample,
                    "Ratio of Median X e^LRR to Median Y e^LRR ("
                            + ext.formDeci(elrrMedY.get(sample) / elrrMedX.get(sample), 4)
                            + ") indicated male");
        } else if (!seedFemales.contains(sample)) {
          notes.put(sample,
                    "Ratio of Median X e^LRR to Median Y e^LRR ("
                            + ext.formDeci(elrrMedY.get(sample) / elrrMedX.get(sample), 4)
                            + ") outlier");
        }

        if (lrrMedX.get(sample) > (femaleMeanX + NUM_SD_FOR_FEMALE_X_OUTLIERS * femaleStdDevX)
            && checkXMosaicism(sample, mosaicismCheckUse)) {
          if (lrrMedX.get(sample) > (femaleMeanX
                                     + NUM_SD_FOR_FEMALE_X_FULL_ANEUPLOIDY * femaleStdDevX)) {
            sexes.put(sample, EstimatedSex.TRIPLE_X);
          } else {
            sexes.put(sample, EstimatedSex.MOSAIC_TRIPLE_X);
          }
        } else if (lrrMedX.get(sample) < (femaleMeanX
                                          - NUM_SD_FOR_FEMALE_X_OUTLIERS * femaleStdDevX)
                   && checkXMosaicism(sample, mosaicismCheckUse)) {
          if (lrrMedX.get(sample) < (femaleMeanX
                                     - NUM_SD_FOR_FEMALE_X_FULL_ANEUPLOIDY * femaleStdDevX)
              && pctXBaf15_85.get(sample) < (femaleMeanPctXBaf15_85 - NUM_SD_FOR_HET_OUTLIERS
                                                                      * femaleStdDevPctXBaf15_85)) {
            sexes.put(sample, EstimatedSex.TURNER);
          } else {
            sexes.put(sample, EstimatedSex.MOSAIC_TURNER);
          }
        } else {
          sexes.put(sample, EstimatedSex.FEMALE);
        }
      } else {
        sexes.put(sample, EstimatedSex.UNKNOWN);
        uncertains.add(sample);
      }
      proj.getProgressMonitor().updateTask(taskName);
    }
    proj.getProgressMonitor().endTask(taskName);
  }

  /*
   * Generates a set of every marker except for X chromosome markers that were excluded. Autosomal
   * markers are needed for the mosaicism checker but we want to only check mosaicism on "good" X
   * chromosome markers
   */
  private Set<Marker> mosaicismUse() {
    Set<Marker> allMarkers = ImmutableSet.copyOf(markerSet.getMarkers());
    Set<Marker> excludedXMarkers = Sets.difference(xMarkers, xUseMarkers);
    return Sets.difference(allMarkers, excludedXMarkers);
  }

  private boolean checkXMosaicism(String sample, Set<Marker> use) {
    Map<Marker, Double> bafs = proj.getPartialSampleFromRandomAccessFile(sample, false, false, true,
                                                                         false, false)
                                   .markerBAFMap(proj.getMarkerSet());
    MosaicBuilder mosaicBuilder = new MosaicBuilder();
    mosaicBuilder.use(use);
    MosaicismDetect mosaicismDetect = mosaicBuilder.build(proj, sample, bafs);
    NavigableSet<Marker> chrXMarkers = markerSet.getChrMap().get((byte) 23);
    int xStart = chrXMarkers.first().getPosition();
    int xStop = chrXMarkers.last().getPosition();
    Segment xSegment = new Segment((byte) 23, xStart, xStop);
    LocusSet<MosaicRegion> xMosaic = mosaicismDetect.callMosaic(xSegment, false);
    int numRegions = xMosaic.getLoci().length;
    if (numRegions == 0) {
      return false;
    }
    String notesAdd = (numRegions == 1 ? "Region" : (numRegions + " regions"))
                      + " of X chromosome mocaicism identified: ";
    double totalCoverage = 0.0;
    double weightedSumF = 0.0;
    for (MosaicRegion mr : xMosaic.getLoci()) {
      double regionCoverage = (double) mr.getSize() / xSegment.getSize();
      totalCoverage += regionCoverage;
      weightedSumF += mr.getCustomF() * regionCoverage;
      notesAdd += "F=" + ext.formDeci(mr.getCustomF(), 4) + ", "
                  + ext.formPercent(regionCoverage, 4) + " coverage (" + mr.getStart() + " - "
                  + mr.getStop() + "); ";
    }
    notesAdd += "Total Mosaic Coverage: " + ext.formPercent(totalCoverage, 4) + "; ";
    if (totalCoverage < MOSAIC_COVERAGE_CERTAINTY_THRESHOLD) {
      if (totalCoverage > MOSAIC_COVERAGE_ABSOLUTE_THRESHOLD) {
        uncertains.add(sample);
      } else {
        return false;
      }
    }
    double weightedAverageF = weightedSumF / totalCoverage;
    if (weightedAverageF < MOSAIC_F_CERTAINTY_THRESHOLD) {
      uncertains.add(sample);
    }
    notes.put(sample, notesAdd);
    return true;
  }

  private void writeToFile(boolean appendToSampleData) {
    Multimap<Integer, String> regionLists = ArrayListMultimap.create();

    SampleData sampleData = proj.getSampleData(false);
    String resultsDir = new File(proj.SEXCHECK_RESULTS_FILENAME.getValue(true, false)).getParent()
                        + "/";
    Map<String, String> pedigreeMap = null;
    final String pedFile = proj.PEDIGREE_FILENAME.getValue();
    if (Files.exists(pedFile)) {
      log.report("Loading Pedigree file, assuming standard pedigree.dat file format (FID, IID, FA, MO, SEX, PHENO, DNA)");
      pedigreeMap = HashVec.loadFileToHashString(pedFile, 6, new int[] {4}, "\t", false);
    }

    try (PrintWriter writer = Files.openAppropriateWriter(proj.SEXCHECK_RESULTS_FILENAME.getValue(true,
                                                                                                  false))) {
      writer.println(ArrayUtils.toStr(SEX_HEADER));

      for (String sample : sampleNames) {
        String fid = sampleData.lookupFID(sample);
        String iid = sampleData.lookupIID(sample);
        EstimatedSex sex = sexes.get(sample);
        BinarySex binSex = sex.getBinarySex();

        if (pedigreeMap != null && binSex == BinarySex.UNKNOWN && pedigreeMap.containsKey(sample)) {
          binSex = BinarySex.values()[Integer.parseInt(pedigreeMap.get(sample))];
        }

        StringJoiner lineBuilder = new StringJoiner("\t");

        lineBuilder.add(sample);
        if (fid == null) {
          log.reportError("Error - no data for sample '" + sample + "'");
          lineBuilder.add(".").add(".").add("-9");
        } else {
          lineBuilder.add(fid).add(iid).add(String.valueOf(binSex.ordinal()));
        }
        String note = notes.get(sample).isEmpty() ? "." : Joiner.on("; ").join(notes.get(sample));
        String uncertain = uncertains.contains(sample) ? "1" : "0";
        String qcPassed = qcPassedSamples.contains(sample) ? "0" : "1";
        double elrrX = elrrMedX.get(sample);
        double elrrY = elrrMedY.get(sample);
        lineBuilder.add(String.valueOf(sex.ordinal())).add(note).add(uncertain).add(qcPassed)
                   .add(String.valueOf(elrrX)).add(String.valueOf(elrrY))
                   .add(String.valueOf(elrrY / elrrX)).add(String.valueOf(pctXHets.get(sample)))
                   .add(String.valueOf(pctXBaf15_85.get(sample)))
                   .add(String.valueOf(lrrMedX.get(sample)))
                   .add(String.valueOf(lrrMedY.get(sample)));
        writer.println(lineBuilder.toString());
        // Create sex-specific region files for any "unusual" sex call to allow easy review
        // TODO it would be nice to add these to Trailer automatically but not clear if that
        // requires a UI, and public API
        String dna = sampleData.lookupDNA(sample);
        if (dna == null) dna = sample;
        if (!qcPassedSamples.contains(sample)) {
          addRegion(regionLists, EstimatedSex.values().length, dna, "chr1", "(excluded ) " + note);
        } else if (uncertains.contains(sample)) {
          addRegion(regionLists, EstimatedSex.values().length + 1, dna, "chr1",
                    "(uncertain ) " + note);
        } else if (sex != EstimatedSex.MALE && sex != EstimatedSex.FEMALE) {
          addRegion(regionLists, sex.ordinal(), dna, sex == EstimatedSex.UNKNOWN ? "chr1" : "chr23",
                    "(" + sex.toString().replaceAll("\\s", "") + " ) " + note);
        }
      }
    } catch (Exception e) {
      log.reportError("Error writing to " + proj.SEXCHECK_RESULTS_FILENAME.getValue());
      log.reportException(e);
    }
    writeSexRegions(regionLists, resultsDir + "sexCheck_regions.txt");
    if (appendToSampleData) {
      if (!sampleData.addData(sexes.entrySet().stream()
                                   .collect(ImmutableMap.toImmutableMap(Entry::getKey,
                                                                        e -> String.valueOf(e.getValue()
                                                                                             .getBinarySex()
                                                                                             .ordinal()))),
                              "DNA", new String[] {"CLASS=Sex"}, ".", "", log)) {
        log.reportError("Error - failed to write Binarized Sex to sample data file");
      }
      if (!sampleData.addData(sexes.entrySet().stream()
                                   .collect(ImmutableMap.toImmutableMap(Entry::getKey,
                                                                        e -> String.valueOf(e.getValue()
                                                                                             .ordinal()))),
                              "DNA", new String[] {"CLASS=" + EST_SEX_HEADER}, ".", "", log)) {
        log.reportError("Error - failed to write Estimated Sex to sample data file");
      }
    }
  }

  /**
   * Adds a string of the format "dna\tchr\tnote" to the list for the specified sex value. Lists are
   * created if they do not already exist.
   */
  private void addRegion(Multimap<Integer, String> regionLists, int sex, String dna, String chr,
                         String note) {
    regionLists.put(sex, dna + "\t" + chr + "\t" + note);
  }

  /**
   * Write all samples of all regions to the specified path
   */
  private void writeSexRegions(Multimap<Integer, String> regions, String path) {
    try (PrintWriter out = Files.openAppropriateWriter(path)) {
      log.report("SexChecks -- Creating sex-specific region file: " + path);
      for (Collection<String> samples : regions.asMap().values()) {
        if (samples == null) {
          continue;
        }
        // Insert count information to the region comment
        String suffix = " of " + samples.size();
        int count = 0;
        for (String line : samples) {
          int parIndex = line.indexOf(')');
          out.println(line.substring(0, parIndex) + count++ + suffix + line.substring(parIndex));
        }
      }
    } catch (IOException e) {
      log.reportError("Error writing to " + path);
    }
  }

  /**
   * Maps complete code (e.g. "1=Male") to definitive sex
   */
  public static int mapEstimatedSexToSex(String estCode) {
    for (EstimatedSex sex : EstimatedSex.values()) {
      if (sex.headerString().startsWith(estCode)) {
        return sex.ordinal();
      }
    }
    return 0;
  }

  /**
   * returns the trinary (0=unknown, 1=male, 2=female) sex for the given code value (e.g. "1")
   */
  public static int getMappedSex(String estimatedValue) {
    return EstimatedSex.values()[Integer.parseInt(estimatedValue)].getBinarySex().ordinal();
  }

  /**
   * returns the number of X chromosomes present for a given code value (e.g female=2, Triple X =3)
   */
  public static int getNumXChrSex(String estimatedValue) {
    return EstimatedSex.values()[Integer.parseInt(estimatedValue)].getChrXCount();
  }

  public static void sexCheck(Project proj, boolean appendToSampleData) {
    sexCheck(proj, appendToSampleData, null);
  }

  public static void sexCheck(Project proj, boolean appendToSampleData,
                              String nonCrossHybridizingMarkersFile) {
    new SexChecks(proj, appendToSampleData, nonCrossHybridizingMarkersFile);
  }

  public static void markerByMarker(Project proj) {
    PrintWriter writer;
    String[] samples;
    Vector<double[]> xys, baflrrs;
    float[] xs, ys, lrrs, bafs;
    Vector<String> intensityDeps;
    LogisticRegression lr;
    String output;
    SampleData sampleData;
    int[] sexes;
    MarkerDataLoader markerDataLoader;
    MarkerData markerData;
    String[] markerNames;
    long time;
    Logger log;

    log = proj.getLog();
    sampleData = proj.getSampleData(false);
    samples = proj.getSamples();
    sexes = new int[samples.length];
    for (int i = 0; i < samples.length; i++) {
      sexes[i] = sampleData.getSexForIndividual(samples[i]);
    }

    try {
      writer = Files.openAppropriateWriter(proj.RESULTS_DIRECTORY.getValue(false, true)
                                           + "markerGenderChecks.xln");
      writer.println("SNP\tX abs(T)\tY abs(T)\tBAF abs(T)\tLRR abs(T)\tX p\tY p\tXY r2\tBAF p\tLRR p\tBAF/LRR r2");

      time = new Date().getTime();
      markerNames = proj.getMarkerNames();
      markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames);
      for (int i = 0; i < markerNames.length; i++) {
        markerData = markerDataLoader.requestMarkerData(i);

        output = markerData.getMarkerName();

        xs = markerData.getXs();
        ys = markerData.getYs();
        bafs = markerData.getBAFs();
        lrrs = markerData.getLRRs();

        intensityDeps = new Vector<>();
        xys = new Vector<>();
        baflrrs = new Vector<>();
        for (int s = 0; s < samples.length; s++) {
          if (ext.isValidDouble(lrrs[s] + "")) {
            intensityDeps.add(sexes[s] + "");
            xys.add(new double[] {xs[s], ys[s]});
            baflrrs.add(new double[] {bafs[s], lrrs[s]});
          }
        }

        if (intensityDeps.size() == 0) {
          log.reportError("Warning - no data for marker " + markerData.getMarkerName());
          output += "\t.\t.\t.\t.\t.\t.\t.\t.";
        } else {
          output += "\t"
                    + Math.abs(new Ttest(ArrayUtils.toIntArray(ArrayUtils.toStringArray(intensityDeps)),
                                         Matrix.extractColumn(Matrix.toDoubleArrays(xys), 0))
                                                                                             .getPvalue());
          output += "\t"
                    + Math.abs(new Ttest(ArrayUtils.toIntArray(ArrayUtils.toStringArray(intensityDeps)),
                                         Matrix.extractColumn(Matrix.toDoubleArrays(xys), 1))
                                                                                             .getPvalue());
          output += "\t" + Math
                               .abs(new Ttest(ArrayUtils.toIntArray(ArrayUtils.toStringArray(intensityDeps)),
                                              Matrix.extractColumn(Matrix.toDoubleArrays(baflrrs),
                                                                   0)).getPvalue());
          output += "\t" + Math
                               .abs(new Ttest(ArrayUtils.toIntArray(ArrayUtils.toStringArray(intensityDeps)),
                                              Matrix.extractColumn(Matrix.toDoubleArrays(baflrrs),
                                                                   1)).getPvalue());
        }

        lr = null;
        try {
          lr = new LogisticRegression(intensityDeps, xys);
          output += "\t" + lr.getSigs()[1] + "\t" + lr.getSigs()[2] + "\t"
                    + (lr.getRsquare() < 0 ? "." : lr.getRsquare());
        } catch (Exception e) {
          output += "\t.\t.\t.";
        }
        try {
          lr = new LogisticRegression(intensityDeps, baflrrs);
          output += "\t" + lr.getSigs()[1] + "\t" + lr.getSigs()[2] + "\t"
                    + (lr.getRsquare() < 0 ? "." : lr.getRsquare());
        } catch (Exception e) {
          output += "\t.\t.\t.";
        }

        writer.println(output);
        markerDataLoader.releaseIndex(i);
      }
      log.reportError("Finished in " + ext.getTimeElapsed(time));
      writer.flush();
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing results");
      log.reportException(e);
    }

  }

  public static void dropMarkers(String allMarkers, String markersToDrop) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    HashSet<String> hashSet;

    hashSet = HashVec.loadFileToHashSet(markersToDrop, false);

    try {
      reader = new BufferedReader(new FileReader(allMarkers));
      writer = Files.openAppropriateWriter(ext.rootOf(allMarkers) + "_dropped.out");
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        if (!hashSet.contains(line[0])) {
          writer.println(ArrayUtils.toStr(line));
        }
      }
      writer.flush();
      writer.close();
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + allMarkers + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + allMarkers + "\"");
      System.exit(2);
    }
  }

  public static void identifyPseudoautosomalBreakpoints(Project proj) {
    PrintWriter writer;
    String[] samples;
    float[] lrrs;
    MarkerData markerData;
    SampleData sampleData;
    int[] sexes;
    byte[] abGenotypes;
    String markerName;
    ClusterFilterCollection clusterFilterCollection;
    float gcThreshold;
    long time;
    DoubleVector[] values; // sex
    MarkerDataLoader markerDataLoader;
    String[] markerList;
    String line, eol;
    MarkerSetInfo markerSet;
    String[] markerNames;
    boolean[] sexChrs;
    byte[] chrs;
    int[][] genotypeCounts;
    boolean[] samplesToExclude;
    Logger log;

    if (Files.isWindows()) {
      eol = "\r\n";
    } else {
      eol = "\n";
    }

    log = proj.getLog();
    sampleData = proj.getSampleData(false);
    samplesToExclude = proj.getSamplesToExclude();
    samples = proj.getSamples();
    sexes = new int[samples.length];
    for (int i = 0; i < samples.length; i++) {
      sexes[i] = Math.max(0, sampleData.getSexForIndividual(samples[i]));
    }

    markerSet = proj.getMarkerSet();
    markerNames = markerSet.getMarkerNames();
    chrs = markerSet.getChrs();
    sexChrs = new boolean[chrs.length];
    for (int i = 0; i < chrs.length; i++) {
      sexChrs[i] = chrs[i] >= 23;
    }
    markerList = ArrayUtils.subArray(markerNames, sexChrs);

    clusterFilterCollection = proj.getClusterFilterCollection();
    // gcThreshold = Float.parseFloat(proj.getProperty(Project.GC_THRESHOLD));
    gcThreshold = proj.getProperty(proj.GC_THRESHOLD).floatValue();

    try {
      writer = Files.openAppropriateWriter(proj.RESULTS_DIRECTORY.getValue(true, false)
                                           + "pseudoautosomalSearch.xln");
      writer.println("SNP\tChr\tPosition\tmLRR_M\tmLRR_F\thet_M\thet_F\tmiss_M\tmiss_F");

      markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerList);
      time = new Date().getTime();
      line = "";
      for (int i = 0; i < markerList.length; i++) {
        markerData = markerDataLoader.requestMarkerData(i);

        markerName = markerData.getMarkerName();
        lrrs = markerData.getLRRs();
        abGenotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerName,
                                                            gcThreshold, log);

        genotypeCounts = new int[2][4]; // sex, genotype
        values = new DoubleVector[2]; // sex
        values[0] = new DoubleVector();
        values[1] = new DoubleVector();
        for (int s = 0; s < samples.length; s++) {
          if (ext.isValidDouble(lrrs[s] + "") && !samplesToExclude[s]) {
            if (sexes[s] == 1 || sexes[s] == 2) {
              values[sexes[s] - 1].add((double) lrrs[s]);
              genotypeCounts[sexes[s] - 1][abGenotypes[s] + 1]++;
            }
          }
        }

        line += markerName + "\t" + markerData.getChr() + "\t" + markerData.getPosition();
        if (values[0].size() > 0) {
          line += "\t" + ArrayUtils.mean(Doubles.toArray(values[0]));
        } else {
          line += "\t.";
        }
        if (values[1].size() > 0) {
          line += "\t" + ArrayUtils.mean(Doubles.toArray(values[1]));
        } else {
          line += "\t.";
        }
        if (genotypeCounts[0][1] + genotypeCounts[0][2] + genotypeCounts[0][3] > 0) {
          line += "\t"
                  + (double) genotypeCounts[0][2]
                    / (double) (genotypeCounts[0][1] + genotypeCounts[0][2] + genotypeCounts[0][3]);
        } else {
          line += "\t.";
        }
        if (genotypeCounts[1][1] + genotypeCounts[1][2] + genotypeCounts[1][3] > 0) {
          line += "\t"
                  + (double) genotypeCounts[1][2]
                    / (double) (genotypeCounts[1][1] + genotypeCounts[1][2] + genotypeCounts[1][3]);
        } else {
          line += "\t.";
        }
        line += "\t" + genotypeCounts[0][0];
        line += "\t" + genotypeCounts[1][0];
        line += eol;

        if (line.length() > 25000) {
          writer.print(line);
          writer.flush();
          line = "";
        }
        markerDataLoader.releaseIndex(i);
      }
      writer.print(line);
      log.report("Identified pseudo-autosomal breakpoints from " + markerList.length
                 + " markers in " + ext.getTimeElapsed(time));
      writer.flush();
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing results");
      log.reportException(e);
    }

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    boolean check = false;
    boolean skipSampleData = false;
    String useMarkers = null;
    String markersToDrop = "data/drops.dat";
    String allMarkers = "data/markerListWithIndices.dat";
    boolean drop = false;
    Project proj;
    String filename = null;
    boolean par = false;

    String usage = "\\n" + "qc.SexChecks requires 0-1 arguments\n"
                   + "   (1) project properties filename (i.e. proj="
                   + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
                   + " AND\n" + "   (2) check sex of indiviudals (i.e. -check (not the default))\n"
                   + "   (3) skip adding estimated sex to Sample Data (i.e. -skipSampleData (not the default))\n"
                   + "   (4) filename of list of markers that do not cross hybridize to use for sex determination (i.e. useMarkers=oneHitWonders.txt (not the default))\n"
                   + " OR\n" + "   (2) drop markers (i.e. -drop (not the default))\n"
                   + "   (3) file with all markers (i.e. all=" + allMarkers + " (default file))\n"
                   + "   (4) list of bad markers (i.e. drop=" + markersToDrop + " (default file))\n"
                   + " OR\n"
                   + "   (2) check sex chromosomes for pseudoautosomal regions (i.e. -PARcheck (not the default))\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-check")) {
        check = true;
        numArgs--;
      } else if (arg.startsWith("-skipSampleData")) {
        skipSampleData = true;
        numArgs--;
      } else if (arg.startsWith("useMarkers=")) {
        useMarkers = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-drop")) {
        drop = true;
        numArgs--;
      } else if (arg.startsWith("all=")) {
        allMarkers = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("drop=")) {
        markersToDrop = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("-PARcheck")) {
        par = true;
        numArgs--;
      }
    }

    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      proj = new Project(filename);

      if (check) {
        sexCheck(proj, !skipSampleData, useMarkers);
      } else if (par) {
        identifyPseudoautosomalBreakpoints(proj);
      } else if (drop) {
        dropMarkers(allMarkers, markersToDrop);
      } else {
        markerByMarker(proj);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
