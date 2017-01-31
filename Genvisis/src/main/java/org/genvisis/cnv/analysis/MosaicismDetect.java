package org.genvisis.cnv.analysis;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.concurrent.Callable;

import org.genvisis.CLI;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.hmm.PennHmm.ViterbiResult;
import org.genvisis.cnv.var.MosaicRegion;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Numbers;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.CNVariant.CNVBuilder;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;

import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

import be.ac.ulg.montefiore.run.distributions.GaussianMixtureDistribution;

/**
 * @author lane0212
 *
 *         Class that attempts to detect and determine mosiac break points...currently models bafs
 *         as gaussian mixture (n=3) to determine likely candidates
 */
public class MosaicismDetect {

	private static final int DEFAULT_MOVING_FACTOR = 25;
	private static final double DEFAULT_NULL_SIGMA = 2;
	private static final double DEFAULT_BASELINE = -1;

	private final Project proj;
	private final String sample;
	private final MarkerSet markerSet;
	private final int movingFactor;
	private final double[] bafs;
	private GaussianMixtureDistribution gd;
	private final double nullSigma;
	private final boolean verbose;
	private final int[][] indicesByChr;
	private double baseLine;
	private double[] means;
	private double[] variances;
	private final Hashtable<String, Integer> markerIndices;// can provide a speedup if provided
	private final boolean[] use;

	public int getMovingFactor() {

		return movingFactor;
	}

	public int[][] getIndicesByChr() {
		return indicesByChr;
	}

	public MarkerSet getMarkerSet() {
		return markerSet;
	}



	public boolean[] getUse() {
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
		int[] segIndices = null;
		if (markerIndices == null) {
			segIndices = markerSet.getIndicesOfMarkersIn(seg, indicesByChr, proj.getLog());
		} else {
			String[] tmp = markerSet.getMarkersIn(seg, indicesByChr);
			segIndices = new int[tmp.length];
			for (int i = 0; i < tmp.length; i++) {
				segIndices[i] = markerIndices.get(tmp[i]);
			}
		}
		if (use != null) {
			segIndices = removeMasked(segIndices);
		}
		int totalIndices = segIndices.length;
		ArrayList<Integer> evalIndicestmp = new ArrayList<Integer>(totalIndices / 2);
		LocusSet<CNVariant> dud = new LocusSet<CNVariant>(new CNVariant[0], true, proj.getLog());

		LocusSet<MosaicRegion> mSet = new LocusSet<MosaicRegion>(	new MosaicRegion[0], true,
																															proj.getLog()) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};
		double[] p_density = new double[segIndices.length];
		double[] nearestN = new double[segIndices.length];
		for (int i = 0; i < segIndices.length; i++) {
			double baf = bafs[segIndices[i]];
			p_density[i] = 0;
			for (int j = 0; j < gd.distributions().length; j++) {
				if (j == 0 || j == 2 || Double.isNaN(baf)) {
					if (j == 0 && Math.abs(baf - means[j]) < nullSigma * Math.sqrt(variances[j])) {
						p_density[i] = Double.NaN;
						nearestN[i] = -1;
					} else if (Math.abs(baf - means[j]) < nullSigma * Math.sqrt(variances[j])) {
						p_density[i] = Double.NaN;
						nearestN[i] = -1;
					}
				}
				double test = !Numbers.isFinite(baf) ? 0 : (double) (baf);
				double tmp = gd.distributions()[j].probability(test) * Math.sqrt(variances[j]);
				if (tmp > p_density[i] && !Double.isNaN(p_density[i])) {
					p_density[i] = tmp;
					if (Numbers.isFinite(baf)) {
						if (j == 0 || j == 2) {
							nearestN[i] = baf < gd.distributions()[1].mean()
																																? Math.max(baf
																																						- gd.distributions()[j].mean(),
																																						0)
																																: gd.distributions()[j].mean()
																																	- baf;
						} else {
							nearestN[i] =
													baf < gd.distributions()[1].mean()	? gd.distributions()[1].mean() - baf
																															: baf - gd.distributions()[1].mean();
						}
					} else {
						nearestN[i] = -1;
					}
				}
			}
			if (!Double.isNaN(p_density[i])) {
				evalIndicestmp.add(i);
			}

		}
		int[] evalIndices = Ints.toArray(evalIndicestmp);
		segIndices = ArrayUtils.subArray(segIndices, evalIndices);
		if (segIndices.length > 0) {

			p_density = ArrayUtils.subArray(p_density, evalIndices);
			nearestN = ArrayUtils.subArray(nearestN, evalIndices);
			double[] p_densityMA = ArrayUtils.movingAverageForward(movingFactor, p_density, true);
			double[] p_densityMAReverse = ArrayUtils.reverse(ArrayUtils.movingAverageForward(	movingFactor,
																																												ArrayUtils.reverse(p_density),
																																												true));
			int[] states = new int[p_densityMA.length];
			ArrayList<Double> p_densityScored = new ArrayList<Double>(p_densityMA.length);
			ArrayList<Integer> mosIndicesTmp = new ArrayList<Integer>(p_densityMA.length);
			for (int i = 0; i < p_densityMA.length; i++) {
				double[] tD = ArrayUtils.removeNaN(new double[] {p_densityMA[i], p_densityMAReverse[i]});
				double d = tD.length > 0 ? ArrayUtils.mean(tD) : Double.NaN;
				p_densityScored.add(d);
				if (Numbers.isFinite(d)) {
					if (d <= baseLine || force) {
						states[i] = 0;
						mosIndicesTmp.add(i);
					} else {
						states[i] = 2;
						mosIndicesTmp.add(i);

					}
				} else {
					throw new IllegalStateException("Currently NaNs should have been removed");
				}
			}
			int[] mosIndices = Ints.toArray(mosIndicesTmp);
			int[] positions =
											ArrayUtils.subArray(ArrayUtils.subArray(markerSet.getPositions(), segIndices),
																					mosIndices);
			String[] names = ArrayUtils.subArray(	ArrayUtils.subArray(markerSet.getMarkerNames(),
																																segIndices),
																						mosIndices);
			double[] bafsSub = ArrayUtils.subArray(ArrayUtils.subArray(bafs, segIndices), mosIndices);
			ViterbiResult vtr = new ViterbiResult(ArrayUtils.subArray(states, mosIndices), null);
			dud = vtr.analyzeStateSequence(	proj, sample, sample, seg.getChr(), positions, names, 2, false,
																			verbose);
			MosaicRegion[] tmp = new MosaicRegion[dud.getLoci().length];
			double[] finalPDensit = Doubles.toArray(p_densityScored);
			for (int i = 0; i < dud.getLoci().length; i++) {
				CNVBuilder builder = new CNVBuilder(dud.getLoci()[i]);
				int numFMarkers = dud.getLoci()[i].getNumMarkers();
				builder.numMarkers(markerSet.getIndicesOfMarkersIn(	dud.getLoci()[i], indicesByChr,
																														proj.getLog()).length);
				if (force) {
					builder.chr(seg.getChr());
					builder.start(seg.getStart());
					builder.stop(seg.getStop());
				}

				int[] scoreStopStart = vtr.getIndexStateChange().get(i);
				double[] scored =
												ArrayUtils.subArray(finalPDensit, scoreStopStart[0], scoreStopStart[1] + 1);
				double pdfScore = baseLine - ArrayUtils.mean(scored);
				double delta = ArrayUtils.median(ArrayUtils.removeNaN(ArrayUtils.distFrom(ArrayUtils.subArray(bafsSub,
																																																			scoreStopStart[0],
																																																			scoreStopStart[1] + 1),
																																									gd.distributions()[1].mean())));
				double factor = dud.getLoci()[i].getSize(); // factor = factor * (double)
																										// dud.getLoci()[i].getNumMarkers() /
																										// states.length;
				double customF = MosaicismQuant.getDisomyF(delta);
				builder.score(customF);
				double nearestStateScore = ArrayUtils.mean(ArrayUtils.subArray(	nearestN, scoreStopStart[0],
																																				scoreStopStart[1] + 1));
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
		int[] autosomalIndices = proj.getAutosomalMarkerIndices();
		if (use != null) {
			autosomalIndices = removeMasked(autosomalIndices);
		}
		double[] autosomalBafs = ArrayUtils.subArray(bafs, autosomalIndices);
		gd = prepareGaussMixture(autosomalBafs, .33, .66);
		if (baseLine < 0) {
			baseLine = 0;
			for (int j = 0; j < gd.distributions().length; j++) {
				baseLine +=
									gd.distributions()[j].probability(means[j] + nullSigma * Math.sqrt(variances[j]))
										* Math.sqrt(variances[j]);
			}
		}
		reportDynRange();
	}

	private int[] removeMasked(int[] autosomalIndices) {
		ArrayList<Integer> notMasked = new ArrayList<Integer>();
		for (int autosomalIndice : autosomalIndices) {
			if (use[autosomalIndice]) {
				notMasked.add(autosomalIndice);
			}
		}
		return Ints.toArray(notMasked);
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
			variances[0] = Numbers.isFinite(zero_tsMeanVar[1]) && zero_tsMeanVar[1] > 0
																																									? zero_tsMeanVar[1]
																																									: 1;
			means[1] = t_tsMeanVar[0];
			variances[1] = Numbers.isFinite(t_tsMeanVar[1]) && t_tsMeanVar[1] > 0 ? t_tsMeanVar[1] : 1;
			means[2] = t_sMeanVar[0];
			variances[2] = Numbers.isFinite(t_sMeanVar[1]) && t_sMeanVar[1] > 0 ? t_sMeanVar[1] : 1;

			if (!Numbers.isFinite(zero_tsMeanVar[1] + t_tsMeanVar[1] + t_sMeanVar[1])
					|| zero_tsMeanVar[1] <= 0 || t_tsMeanVar[1] <= 0 || t_sMeanVar[1] <= 0) {
				proj.getLog().reportTimeWarning("Sample "+ sample
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
	 *
	 */
	public static class MosaicBuilder {
		// init to default params...
		private int movingFactor = DEFAULT_MOVING_FACTOR;
		private double nullSigma = DEFAULT_NULL_SIGMA;
		private boolean verbose = false;
		private int[][] indicesByChr = null;
		private double baseLine = DEFAULT_BASELINE;
		private double[] means = null;
		private double[] variances = null;
		private Hashtable<String, Integer> markerIndices = null;// can provide a speedup if provided
		private boolean[] use = null;// only these markers will be used for the computation

		/**
		 * @param use mask markers from computation
		 * @return
		 */
		public MosaicBuilder use(boolean[] use) {
			this.use = use;
			return this;
		}

		/**
		 * @param markerIndices speeds up extraction of which markers to use for a given seg
		 * @return
		 */
		public MosaicBuilder markerIndices(Hashtable<String, Integer> markerIndices) {
			this.markerIndices = markerIndices;
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
		 *        distributions will not be counted
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
		 * @param indicesByChr project marker indices by chromosome
		 * @return
		 */
		public MosaicBuilder indicesByChr(int[][] indicesByChr) {
			this.indicesByChr = indicesByChr;
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
		 * @param markerSet
		 * @param bafs
		 * @return a Mosaic detector
		 */
		public MosaicismDetect build(Project proj, String sample, MarkerSet markerSet, double[] bafs) {
			return new MosaicismDetect(this, proj, sample, markerSet, bafs);
		}
	}

	private static class MosaicWorker implements Callable<LocusSet<MosaicRegion>> {
		private final Project proj;
		private final MosaicBuilder builder;
		private final MarkerSet markerSet;
		private final LocusSet<Segment> segs;
		private final String sample;
		private double[] bafs;
		private double[] lrrs;
		private final int[][] indicesByChr;

		public MosaicWorker(Project proj, MosaicBuilder builder, LocusSet<Segment> segs,
												MarkerSet markerSet, String sample) {
			super();
			this.proj = proj;
			this.builder = builder;
			this.markerSet = markerSet;
			this.segs = segs;
			this.sample = sample;
			indicesByChr = markerSet.getIndicesByChr();

		}

		@Override
		public LocusSet<MosaicRegion> call() throws Exception {
			Sample samp = proj.getFullSampleFromRandomAccessFile(sample);
			bafs = ArrayUtils.toDoubleArray(samp.getBAFs());
			lrrs = ArrayUtils.toDoubleArray(samp.getLRRs());
			builder.indicesByChr(indicesByChr);
			MosaicismDetect md = builder.build(proj, sample, markerSet, bafs);
			ArrayList<MosaicRegion> all = new ArrayList<MosaicRegion>();
			for (int i = 0; i < segs.getLoci().length; i++) {
				LocusSet<MosaicRegion> tmp = md.callMosaic(segs.getLoci()[i], false);
				tmp.addAll(all);
			}
			LocusSet<MosaicRegion> allCalls =
																			new LocusSet<MosaicRegion>(	all.toArray(new MosaicRegion[all.size()]),
																																	true, proj.getLog()) {

																				/**
																				 * 
																				 */
																				private static final long serialVersionUID = 1L;

																			};
			BeastScore beastScore = BeastScore.beastInd(proj, null, ArrayUtils.toFloatArray(lrrs),
																									allCalls.getLoci(), md.getMarkerSet().getChrs(),
																									md.getMarkerSet().getPositions(),
																									md.getIndicesByChr());
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
	 *
	 */
	public static class MosaicProducer extends AbstractProducer<LocusSet<MosaicRegion>> {
		private final Project proj;
		private final String[] samples;
		private final MosaicBuilder builder;
		private final LocusSet<Segment> segs;
		private final MarkerSet markerSet;
		private int index;

		/**
		 * @param proj
		 * @param builder builds the detectors
		 * @param samples which samples to use
		 * @param markerSet
		 * @param segs segments to call (like a chromosome)
		 */
		public MosaicProducer(Project proj, MosaicBuilder builder, String[] samples,
													MarkerSet markerSet, LocusSet<Segment> segs) {
			super();
			this.proj = proj;
			this.samples = samples;
			this.markerSet = markerSet;
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
			MosaicWorker worker = new MosaicWorker(proj, builder, segs, markerSet, sample);
			index++;
			return worker;
		}
	}

	/**
	 * @param builder
	 * @param proj
	 * @param sample
	 * @param markerSet
	 * @param bafs
	 */
	private MosaicismDetect(MosaicBuilder builder, Project proj, String sample, MarkerSet markerSet,
													double[] bafs) {
		this.proj = proj;
		this.sample = sample;
		this.markerSet = markerSet;
		this.bafs = bafs;
		movingFactor = builder.movingFactor;
		nullSigma = builder.nullSigma;
		verbose = builder.verbose;
		indicesByChr = builder.indicesByChr == null	? markerSet.getIndicesByChr()
																								: builder.indicesByChr;
		baseLine = builder.baseLine;
		means = builder.means;
		variances = builder.variances;
		markerIndices = builder.markerIndices;
		use = builder.use;
		if (bafs.length != markerSet.getMarkerNames().length) {
			throw new IllegalArgumentException("Internal error, bafs must be present for entire array, fill with NaN if neccesary");
		}
		if (use != null && use.length != markerSet.getMarkerNames().length) {
			throw new IllegalArgumentException("Internal error, use definition must be present for entire array if not null, fill with NaN if neccesary");
		}
		prep();
	}


	/**
	 * @param proj call mosaicism for this project
	 * @param output the output file name
	 * @param numThreads number of threads to use
	 */
	public static void callMosaicRegions(Project proj, String output, int numThreads) {

		MarkerSet markerSet = proj.getMarkerSet();
		int[][] indicesByChr = markerSet.getIndicesByChr();
		MosaicBuilder builder = new MosaicBuilder();// most customizing can be done in the builder if
																								// needed
		builder.indicesByChr(indicesByChr);

		// develop segments to call on
		ArrayList<Segment> callSegs = new ArrayList<Segment>();

		for (int i = 0; i < indicesByChr.length; i++) {
			if (i > 0 && i < 23 && indicesByChr[i].length > 0) {// weird things happen on non-auto
				callSegs.add(new Segment((byte) i, 0, Integer.MAX_VALUE));
			}
		}
		LocusSet<Segment> segs = new LocusSet<Segment>(	callSegs.toArray(new Segment[callSegs.size()]),
																										true, proj.getLog()) {

			/**
			 *
			 */
			private static final long serialVersionUID = 1L;

		};

		MosaicProducer producer = new MosaicProducer(proj, builder, proj.getSamples(), markerSet, segs);
		WorkerTrain<LocusSet<MosaicRegion>> train =
																							new WorkerTrain<LocusSet<MosaicRegion>>(producer,
																																											numThreads,
																																											10,
																																											proj.getLog());
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



	public static void main(String[] args) {
		CLI c = new CLI(MosaicismDetect.class);

		c.addArgWithDefault(CLI.ARG_PROJ, CLI.DESC_PROJ, "proj.properties");
		c.addArgWithDefault(CLI.ARG_OUTFILE, CLI.DESC_OUTFILE, "out.mos");
		c.addArgWithDefault(CLI.ARG_THREADS, CLI.DESC_THREADS, "24");
		c.parseWithExit(args);

		Project proj = new Project(c.get(CLI.ARG_PROJ), false);

		callMosaicRegions(proj, c.get(CLI.ARG_OUTFILE), c.getI(CLI.ARG_THREADS));

	}


}
