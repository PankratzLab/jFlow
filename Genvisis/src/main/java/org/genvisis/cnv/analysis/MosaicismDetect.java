package org.genvisis.cnv.analysis;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.concurrent.Callable;
import javax.annotation.Nullable;
import org.genvisis.CLI;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.hmm.PennHmm.ViterbiResult;
import org.genvisis.cnv.var.MosaicRegion;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.CNVariant.CNVBuilder;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.stats.Maths;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import be.ac.ulg.montefiore.run.distributions.GaussianMixtureDistribution;

/**
 * @author lane0212 Class that attempts to detect and determine mosiac break points...currently
 *         models bafs as gaussian mixture (n=3) to determine likely candidates
 */
public class MosaicismDetect {

  private static final int DEFAULT_MOVING_FACTOR = 25;
  private static final double DEFAULT_NULL_SIGMA = 2;
  private static final double DEFAULT_BASELINE = -1;

  private final Project proj;
  private final String sample;
  private final int movingFactor;
  private final Map<Marker, ? extends Number> bafs;
  private GaussianMixtureDistribution gd;
  private final double nullSigma;
  private final boolean verbose;
  private double baseLine;
  private double[] means;
  private double[] variances;
  @Nullable
  private final Set<Marker> use;

  public int getMovingFactor() {

    return movingFactor;
  }

  public @Nullable Set<Marker> getUse() {
    return use;
  }

  /**
   * @param seg call on this segment
   * @param force if true, the region is assumed to be mosaic, and is essentially just scored
   * @return
   */
  public <T extends Segment> LocusSet<MosaicRegion> callMosaic(T seg, boolean force) {
    if (seg.getStop() < seg.getStart()) {
      throw new IllegalArgumentException("Segment must have stop that is gte to start");
    }
    Set<Marker> markersInSeg = applyUse(proj.getMarkerSet().getMarkersInSeg(seg));
    int totalMarkers = markersInSeg.size();
    Set<Marker> evalMarkers = new HashSet<>(totalMarkers / 2);

    LocusSet<MosaicRegion> mSet = new LocusSet<>(new MosaicRegion[0], true, proj.getLog());
    Map<Marker, Double> pDensities = Maps.newHashMapWithExpectedSize(markersInSeg.size());
    Map<Marker, Double> nearestNs = Maps.newHashMapWithExpectedSize(markersInSeg.size());
    for (Marker marker : markersInSeg) {
      double baf = bafs.get(marker).doubleValue();
      double pDensity = 0.0;
      double nearestN = 0.0;
      for (int j = 0; j < gd.distributions().length; j++) {
        if (j == 0 || j == 2 || Double.isNaN(baf)) {
          if (j == 0 && Math.abs(baf - means[j]) < nullSigma * Math.sqrt(variances[j])) {
            pDensity = Double.NaN;
            nearestN = -1.0;
          } else if (Math.abs(baf - means[j]) < nullSigma * Math.sqrt(variances[j])) {
            pDensity = Double.NaN;
            nearestN = -1.0;
          }
        }
        double test = !Double.isFinite(baf) ? 0 : (double) (baf);
        double tmp = gd.distributions()[j].probability(test) * Math.sqrt(variances[j]);
        if (tmp > pDensity && !Double.isNaN(pDensity)) {
          pDensity = tmp;
          if (Double.isFinite(baf)) {
            if (j == 0 || j == 2) {
              nearestN = baf < gd.distributions()[1].mean() ? Math.max(baf
                                                                       - gd.distributions()[j].mean(),
                                                                       0)
                                                            : gd.distributions()[j].mean() - baf;
            } else {
              nearestN = baf < gd.distributions()[1].mean() ? gd.distributions()[1].mean() - baf
                                                            : baf - gd.distributions()[1].mean();
            }
          } else {
            nearestN = -1;
          }
        }
      }
      if (!Double.isNaN(pDensity)) {
        evalMarkers.add(marker);
        pDensities.put(marker, pDensity);
        nearestNs.put(marker, nearestN);
      }
    }
    ImmutableList<Marker> mosMarkers = ImmutableList.copyOf(Sets.intersection(markersInSeg,
                                                                              evalMarkers));
    if (!mosMarkers.isEmpty()) {
      Map<Marker, Double> pDensityMovingAverage = Maps.newHashMapWithExpectedSize(mosMarkers.size());
      Maths.MovingAverage movingAverage = new Maths.MovingAverage(movingFactor);
      for (Marker marker : mosMarkers) {
        pDensityMovingAverage.put(marker, movingAverage.add(pDensities.get(marker)));
      }

      Map<Marker, Double> pDensityMovingAverageReverse = Maps.newHashMapWithExpectedSize(mosMarkers.size());
      movingAverage = new Maths.MovingAverage(movingFactor);
      for (Marker marker : mosMarkers.reverse()) {
        pDensityMovingAverageReverse.put(marker, movingAverage.add(pDensities.get(marker)));
      }
      Map<Marker, Integer> states = Maps.newHashMapWithExpectedSize(mosMarkers.size());
      Map<Marker, Double> pDensityScored = Maps.newHashMapWithExpectedSize(mosMarkers.size());
      for (Marker marker : mosMarkers) {
        double[] tD = ArrayUtils.removeNaN(new double[] {pDensityMovingAverage.get(marker),
                                                         pDensityMovingAverageReverse.get(marker)});
        double d = tD.length > 0 ? ArrayUtils.mean(tD) : Double.NaN;
        pDensityScored.put(marker, d);
        if (Double.isFinite(d)) {
          if (d <= baseLine || force || marker.getChr() == 21) {
            states.put(marker, 0);
          } else {
            states.put(marker, 2);
          }
        } else {
          throw new IllegalStateException("Currently NaNs should have been removed");
        }
      }

      ViterbiResult vtr = new ViterbiResult(mosMarkers.stream().mapToInt(states::get).toArray(),
                                            null);
      LocusSet<CNVariant> dud = vtr.analyzeStateSequence(proj, sample, sample, seg.getChr(),
                                                         mosMarkers.stream()
                                                                   .mapToInt(Marker::getPosition)
                                                                   .toArray(),
                                                         mosMarkers.stream().map(Marker::getName)
                                                                   .toArray(String[]::new),
                                                         2, false, verbose);
      MosaicRegion[] tmp = new MosaicRegion[dud.getLoci().length];
      proj.getLog().reportTimeInfo("Scoring mosaic regions");
      for (int i = 0; i < dud.getLoci().length; i++) {
        CNVBuilder builder = new CNVBuilder(dud.getLoci()[i]);
        int numFMarkers = dud.getLoci()[i].getNumMarkers();
        builder.numMarkers((int) proj.getMarkerSet().viewMarkersInSeg(dud.getLoci()[i]).count());
        if (force) {
          builder.chr(seg.getChr());
          builder.start(seg.getStart());
          builder.stop(seg.getStop());
        }

        int[] scoreStopStart = vtr.getIndexStateChange().get(i);
        List<Marker> mosMarkersScored = mosMarkers.subList(scoreStopStart[0],
                                                           scoreStopStart[1] + 1);
        double pdfScore = baseLine
                          - ArrayUtils.mean(mosMarkersScored.stream()
                                                            .mapToDouble(pDensityScored::get)
                                                            .toArray());
        double delta = ArrayUtils.median(mosMarkersScored.stream().map(bafs::get)
                                                         .mapToDouble(Number::doubleValue)
                                                         .map(baf -> Math.abs(baf
                                                                              - gd.distributions()[1].mean()))
                                                         .filter(d -> !Double.isNaN(d)).toArray());
        double factor = dud.getLoci()[i].getSize(); // factor = factor * (double)
                                                   // dud.getLoci()[i].getNumMarkers() /
                                                   // states.length;
        double customF = MosaicismQuant.getDisomyF(delta);
        builder.score(customF);
        double nearestStateScore = ArrayUtils.mean(mosMarkersScored.stream()
                                                                   .mapToDouble(nearestNs::get)
                                                                   .toArray());
        tmp[i] = new MosaicRegion(builder.build(), Math.log10(Math.pow(factor, 2)),
                                  nearestStateScore, pdfScore, delta, Double.NaN, customF);
        tmp[i].setNumFMarkers(numFMarkers);
      }

      mSet = new LocusSet<MosaicRegion>(tmp, true, proj.getLog()) {

        /**
        	 *
        	 */
        private static final long serialVersionUID = 1L;

      };

    } else if (force) {// no markers met criteria, set a blank
      CNVariant blank = new CNVariant(sample, sample, seg.getChr(), seg.getStart(), seg.getStop(),
                                      2, Double.NaN, 0, 99);
      MosaicRegion blankMr = new MosaicRegion(blank, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
                                              Double.NaN, Double.NaN);
      mSet = new LocusSet<MosaicRegion>(new MosaicRegion[] {blankMr}, true, proj.getLog()) {

        /**
         *
         */
        private static final long serialVersionUID = 1L;

      };

    }
    return mSet;
  }

  private void prep() {
    Set<Marker> autosomalMarkers = applyUse(ImmutableSet.copyOf(proj.getMarkerSet()
                                                                    .getAutosomalMarkers()));
    gd = prepareGaussMixture(autosomalMarkers.stream().map(bafs::get)
                                             .mapToDouble(Number::doubleValue).toArray(),
                             .33, .66);
    if (baseLine < 0) {
      baseLine = 0;
      for (int j = 0; j < gd.distributions().length; j++) {
        baseLine += gd.distributions()[j].probability(means[j]
                                                      + nullSigma * Math.sqrt(variances[j]))
                    * Math.sqrt(variances[j]);
      }
    }
    reportDynRange();
  }

  private Set<Marker> applyUse(Set<Marker> inputMarkers) {
    if (use != null) {
      return Sets.intersection(inputMarkers, use);
    }
    return inputMarkers;
  }

  private void reportDynRange() {
    if (verbose) {
      double minDelta = Math.sqrt(gd.distributions()[1].variance());
      proj.getLog().reportTimeInfo("Min proportion Disomy detection ~="
                                   + MosaicismQuant.getDisomyF(minDelta));
      double maxDelta = .5 - nullSigma * Math.sqrt(gd.distributions()[2].variance());// B allele
                                                                                     // generally
                                                                                     // greater
                                                                                     // variance
      proj.getLog().reportTimeInfo("Max proportion Disomy detection ~="
                                   + MosaicismQuant.getDisomyF(maxDelta));
      proj.getLog().reportTimeInfo("Null conf ~=" + baseLine);

    }
  }

  private GaussianMixtureDistribution prepareGaussMixture(double[] autosomalBafs, double r1,
                                                          double r2) {
    if (means == null) {

      means = new double[3];
      variances = new double[3];// defualt to 1 if actual variance is not finite (since baf 0->1
                                // anyway)
      double[] zero_tsMeanVar = getMeanVar(autosomalBafs, 0, r1);
      double[] t_tsMeanVar = getMeanVar(autosomalBafs, r1, r2);
      double[] t_sMeanVar = getMeanVar(autosomalBafs, r2, 1);
      means[0] = zero_tsMeanVar[0];
      variances[0] = Double.isFinite(zero_tsMeanVar[1]) && zero_tsMeanVar[1] > 0 ? zero_tsMeanVar[1]
                                                                                 : 1;
      means[1] = t_tsMeanVar[0];
      variances[1] = Double.isFinite(t_tsMeanVar[1]) && t_tsMeanVar[1] > 0 ? t_tsMeanVar[1] : 1;
      means[2] = t_sMeanVar[0];
      variances[2] = Double.isFinite(t_sMeanVar[1]) && t_sMeanVar[1] > 0 ? t_sMeanVar[1] : 1;

      if (!Double.isFinite(zero_tsMeanVar[1] + t_tsMeanVar[1] + t_sMeanVar[1])
          || zero_tsMeanVar[1] <= 0 || t_tsMeanVar[1] <= 0 || t_sMeanVar[1] <= 0) {
        proj.getLog().reportTimeWarning("Sample " + sample
                                        + " had non-finite or 0 baf variance, setting to 1");
      }
    }

    double[] props = new double[3];
    Arrays.fill(props, (double) 1 / 3);

    return new GaussianMixtureDistribution(means, variances, props);
  }

  private static double[] getMeanVar(double[] autosomalBafs, double r1, double r2) {
    double[] sub = ArrayUtils.getValuesBetween(autosomalBafs, r1, r2, true);
    double[] meanVar = new double[2];
    meanVar[0] = ArrayUtils.mean(sub);
    meanVar[1] = Math.pow(ArrayUtils.stdev(sub), 2);
    return meanVar;
  }

  /**
   * builder for {@link MosaicismDetect}
   */
  public static class MosaicBuilder {

    // init to default params...
    private int movingFactor = DEFAULT_MOVING_FACTOR;
    private double nullSigma = DEFAULT_NULL_SIGMA;
    private boolean verbose = false;
    private double baseLine = DEFAULT_BASELINE;
    private double[] means = null;
    private double[] variances = null;
    private Set<Marker> use = null;// only these markers will be used for the computation

    /**
     * @param use mask markers from computation
     * @return
     */
    public MosaicBuilder use(Set<Marker> use) {
      this.use = use;
      return this;
    }

    /**
     * @param movingFactor set the moving average length for the smoothing step
     * @return
     */
    public MosaicBuilder movingFactor(int movingFactor) {
      this.movingFactor = movingFactor;
      if (movingFactor <= 0) {
        throw new IllegalArgumentException("movingFactor must be positive");
      }
      return this;
    }

    /**
     * @param nullSigma points within this many standard deviations of either of the three
     *          distributions will not be counted
     * @return
     */
    public MosaicBuilder nullSigma(double nullSigma) {
      this.nullSigma = nullSigma;
      if (nullSigma <= 0) {
        throw new IllegalArgumentException("nullSigma must be positive");
      }
      return this;
    }

    /**
     * @param verbose verbose output
     * @return
     */
    public MosaicBuilder verbose(boolean verbose) {
      this.verbose = verbose;
      return this;
    }

    /**
     * @param baseLine threshold for calling mosiac regions
     * @return
     */
    public MosaicBuilder baseLine(double baseLine) {
      this.baseLine = baseLine;
      if (baseLine <= 0) {
        throw new IllegalArgumentException("Baseline must be positive");
      }
      return this;
    }

    /**
     * @param means means of the 3 distributions
     * @param variances variances of the 3 distributions;
     * @return
     */
    public MosaicBuilder meansAndVariances(double[] means, double[] variances) {
      this.means = means;
      this.variances = variances;

      if (means == null || means.length != 3) {
        throw new IllegalArgumentException("Internal error, means must be length 3");
      }
      if (variances == null || variances.length != 3) {
        throw new IllegalArgumentException("Internal error, variances must be length 3");
      }
      return this;
    }

    /**
     * @param proj
     * @param sample
     * @param bafs
     * @return a Mosaic detector
     */
    public MosaicismDetect build(Project proj, String sample, Map<Marker, ? extends Number> bafs) {
      return new MosaicismDetect(this, proj, sample, bafs);
    }
  }

  private static class MosaicWorker implements Callable<LocusSet<MosaicRegion>> {

    private final Project proj;
    private final MosaicBuilder builder;
    private final LocusSet<Segment> segs;
    private final String sample;

    public MosaicWorker(Project proj, MosaicBuilder builder, LocusSet<Segment> segs,
                        String sample) {
      super();
      this.proj = proj;
      this.builder = builder;
      this.segs = segs;
      this.sample = sample;
    }

    @Override
    public LocusSet<MosaicRegion> call() throws Exception {
      Sample samp = proj.getFullSampleFromRandomAccessFile(sample);
      MarkerDetailSet markerSet = proj.getMarkerSet();
      Map<Marker, Double> bafs = samp.markerBAFMap(markerSet);
      MosaicismDetect md = builder.build(proj, sample, bafs);
      ArrayList<MosaicRegion> all = new ArrayList<>();
      for (int i = 0; i < segs.getLoci().length; i++) {
        LocusSet<MosaicRegion> tmp = md.callMosaic(segs.getLoci()[i], false);
        tmp.addAll(all);
      }
      LocusSet<MosaicRegion> allCalls = new LocusSet<>(all.toArray(new MosaicRegion[all.size()]),
                                                       true, proj.getLog());
      BeastScore beastScore = BeastScore.beastInd(proj, null, samp.getLRRs(), allCalls.getLoci());
      for (int i = 0; i < allCalls.getLoci().length; i++) {
        allCalls.getLoci()[i].setBeastScore(beastScore.getBeastScores()[i]);
        allCalls.getLoci()[i].setBeastHeight(beastScore.getBeastHeights()[i]);
        allCalls.getLoci()[i].setBeastLength(beastScore.getBeastLengths()[i]);
      }

      return allCalls;
    }
  }

  /**
   * @author Kitty
   */
  public static class MosaicProducer extends AbstractProducer<LocusSet<MosaicRegion>> {

    private final Project proj;
    private final String[] samples;
    private final MosaicBuilder builder;
    private final LocusSet<Segment> segs;
    private int index;

    /**
     * @param proj
     * @param builder builds the detectors
     * @param samples which samples to use
     * @param segs segments to call (like a chromosome)
     */
    public MosaicProducer(Project proj, MosaicBuilder builder, String[] samples,
                          LocusSet<Segment> segs) {
      super();
      this.proj = proj;
      this.samples = samples;
      this.builder = builder;
      this.segs = segs;
      index = 0;
    }

    @Override
    public boolean hasNext() {
      return index < samples.length;
    }

    @Override
    public Callable<LocusSet<MosaicRegion>> next() {
      final String sample = samples[index];
      MosaicWorker worker = new MosaicWorker(proj, builder, segs, sample);
      index++;
      return worker;
    }
  }

  /**
   * @param builder
   * @param proj
   * @param sample
   * @param bafs
   */
  private MosaicismDetect(MosaicBuilder builder, Project proj, String sample,
                          Map<Marker, ? extends Number> bafs) {
    this.proj = proj;
    this.sample = sample;
    this.bafs = bafs;
    movingFactor = builder.movingFactor;
    nullSigma = builder.nullSigma;
    verbose = builder.verbose;
    baseLine = builder.baseLine;
    means = builder.means;
    variances = builder.variances;
    use = builder.use;
    if (bafs.keySet().size() != proj.getMarkerSet().markersAsList().size()) {
      throw new IllegalArgumentException("Internal error, bafs must be present for entire array, fill with NaN if neccesary");
    }
    prep();
  }

  /**
   * @param proj call mosaicism for this project
   * @param output the output file name
   * @param numThreads number of threads to use
   */
  public static void callMosaicRegions(Project proj, String output, int numThreads) {

    MarkerDetailSet markerSet = proj.getMarkerSet();
    MosaicBuilder builder = new MosaicBuilder();// most customizing can be done in the builder if
                                                // needed
    Set<Marker> use = markerSet.markersAsList().stream()
                               .filter(m -> !proj.ARRAY_TYPE.getValue().isCNOnly(m.getName()))
                               .collect(ImmutableSet.toImmutableSet());
    builder.use(use);
    proj.getLog()
        .reportTimeWarning("Skipping " + (markerSet.markersAsList().size() - use.size())
                           + " copy number only markers for mosaic regions and array type "
                           + proj.ARRAY_TYPE.getValue().toString());
    // develop segments to call on
    Segment[] callSegs = markerSet.getChrMap().tailMap((byte) 1).entrySet().stream()
                                  .filter(e -> !e.getValue().isEmpty()).map(Entry::getKey)
                                  .map(chr -> new Segment(chr, 0, Integer.MAX_VALUE))
                                  .toArray(Segment[]::new);
    LocusSet<Segment> segs = new LocusSet<>(callSegs, true, proj.getLog());

    MosaicProducer producer = new MosaicProducer(proj, builder, proj.getSamples(), segs);
    try (WorkerTrain<LocusSet<MosaicRegion>> train = new WorkerTrain<>(producer, numThreads, 10,
                                                                       proj.getLog())) {
      int numCalled = 0;
      boolean wroteHeader = false;

      PrintWriter writer = Files.getAppropriateWriter(output);
      while (train.hasNext()) {
        LocusSet<MosaicRegion> current = train.next();
        numCalled++;
        if (numCalled % 100 == 0) {
          proj.getLog().reportTimeInfo("Called " + numCalled + " samples for mosaicism");
        }
        for (MosaicRegion region : current.getLoci()) {
          if (!wroteHeader) {
            writer.println(ArrayUtils.toStr(region.getHeader()));
            wroteHeader = true;
          }
          // could definitely put a num marker filter here..
          writer.println(region.toAnalysisString());
        }
      }
      writer.close();
    }
  }

  public static void main(String[] args) {
    CLI c = new CLI(MosaicismDetect.class);

    c.addArgWithDefault(CLI.ARG_PROJ, CLI.DESC_PROJ, "proj.properties");
    c.addArgWithDefault(CLI.ARG_OUTFILE, CLI.DESC_OUTFILE, "out.mos");
    c.addArgWithDefault(CLI.ARG_THREADS, CLI.DESC_THREADS, CLI.EXAMPLE_THREADS);
    c.parseWithExit(args);

    Project proj = new Project(c.get(CLI.ARG_PROJ));

    callMosaicRegions(proj, c.get(CLI.ARG_OUTFILE), c.getI(CLI.ARG_THREADS));

  }

}
