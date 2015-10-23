package cnv.analysis;

import java.util.ArrayList;
import java.util.Arrays;

import be.ac.ulg.montefiore.run.distributions.GaussianMixtureDistribution;
import common.Array;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.hmm.PennHmm.ViterbiResult;
import cnv.var.CNVariant;
import cnv.var.LocusSet;
import cnv.var.LocusSet.TO_STRING_TYPE;
import filesys.Segment;
import cnv.var.CNVariant.CNVBuilder;

/**
 * @author lane0212
 *
 *         Class that attempts to detect and determine mosiac break points...currently models bafs as gaussian mixture (n=3) to determine likely candidates
 */
public class MosaicismDetect {

	private static final int DEFAULT_MOVING_FACTOR = 25;
	private static final double DEFAULT_NULL_SIGMA = 2;
	private static final double DEFAULT_BASELINE = -1;
	private static final double DEFAULT_MIN_PERCENT_STATES = 0.05;

	private Project proj;
	private String sample;
	private MarkerSet markerSet;
	private int movingFactor;
	private double[] bafs;
	private GaussianMixtureDistribution gd;
	private double nullSigma;
	private boolean verbose;
	private int[][] indicesByChr;
	private double baseLine;
	private double[] means;
	private double[] variances;
	private double minPercentStates;

	public int getMovingFactor() {

		return movingFactor;
	}

	public <T extends Segment> LocusSet<CNVariant> callMosaic(T seg) {
		if (seg.getStop() < seg.getStart()) {
			throw new IllegalArgumentException("Segment must have stop that is gte to start");
		}
		int[] segIndices = markerSet.getIndicesOfMarkersIn(seg, indicesByChr, proj.getLog());
		int totalIndices = segIndices.length;
		ArrayList<Integer> evalIndicestmp = new ArrayList<Integer>(totalIndices / 2);
		LocusSet<CNVariant> dud = new LocusSet<CNVariant>(new CNVariant[0], true, proj.getLog()) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};
		double[] p_density = new double[segIndices.length];
		for (int i = 0; i < segIndices.length; i++) {
			double baf = bafs[segIndices[i]];
			p_density[i] = 0;
			for (int j = 0; j < gd.distributions().length; j++) {
				if (j == 0 || j == 2) {
					if (j == 0 && Math.abs(baf - means[j]) < nullSigma * Math.sqrt(variances[j])) {
						p_density[i] = Double.NaN;

					} else if (Math.abs(baf - means[j]) < nullSigma * Math.sqrt(variances[j])) {
						p_density[i] = Double.NaN;

					}
				}
				double test = Double.isFinite(baf) ? (double) (baf) : 0;
				double tmp = (double) gd.distributions()[j].probability(test) * Math.sqrt(variances[j]);
				if (tmp > p_density[i] && !Double.isNaN(p_density[i])) {
					p_density[i] = tmp;
				}
			}
			if (!Double.isNaN(p_density[i])) {
				evalIndicestmp.add(i);
			}

		}
		int[] evalIndices = Array.toIntArray(evalIndicestmp);
		segIndices = Array.subArray(segIndices, evalIndices);
		if (segIndices.length > 0) {

			p_density = Array.subArray(p_density, evalIndices);
			double[] p_densityMA = Array.movingAverageForward(movingFactor, p_density, true);
			double[] p_densityMAReverse = Array.reverse(Array.movingAverageForward(movingFactor, Array.reverse(p_density), true));
			int[] states = new int[p_densityMA.length];
			ArrayList<Double> p_densityScored = new ArrayList<Double>(p_densityMA.length);
			ArrayList<Integer> mosIndicesTmp = new ArrayList<Integer>(p_densityMA.length);
			for (int i = 0; i < p_densityMA.length; i++) {
				double[] tD = Array.removeNaN(new double[] { p_densityMA[i], p_densityMAReverse[i] });
				double d = tD.length > 0 ? Array.mean(tD) : Double.NaN;
				p_densityScored.add(d);
				if (Double.isFinite(d)) {
					if (d <= baseLine) {
						states[i] = 0;
						mosIndicesTmp.add(i);
					} else {
						states[i] = 2;
						mosIndicesTmp.add(i);

					}
				} else {
					states[i] = 2;
				}
			}
			int[] mosIndices = Array.toIntArray(mosIndicesTmp);
			double percentState = (double) evalIndices.length / totalIndices;
			// System.out.println(percentState+"\t"+evalIndices.length+"\t"+totalIndices);

			// if (percentState > minPercentStates) {
			int[] positions = Array.subArray(Array.subArray(markerSet.getPositions(), segIndices), mosIndices);
			String[] names = Array.subArray(Array.subArray(markerSet.getMarkerNames(), segIndices), mosIndices);
			ViterbiResult vtr = new ViterbiResult(Array.subArray(states, mosIndices), null);
			dud = vtr.analyzeStateSequence(proj, sample, sample, seg.getChr(), positions, names, 2, false, verbose);
			CNVariant[] tmp = new CNVariant[dud.getLoci().length];
			double[] finalPDensit = Array.toDoubleArray(p_densityScored);
			for (int i = 0; i < dud.getLoci().length; i++) {
				CNVBuilder builder = new CNVBuilder(dud.getLoci()[i]);
				int[] scoreStopStart = vtr.getIndexStateChange().get(i);
				double[] scored = Array.subArray(finalPDensit, scoreStopStart[0], scoreStopStart[1] + 1);
				double score = baseLine - Array.mean(scored);
				double factor = (double) dud.getLoci()[i].getSize(); // factor = factor * (double) dud.getLoci()[i].getNumMarkers() / states.length;
				score = Math.log10(score * factor);
				builder.score(score);
				tmp[i] = builder.build();
			}

			dud = new LocusSet<CNVariant>(tmp, true, proj.getLog()) {

				/**
					 * 
					 */
				private static final long serialVersionUID = 1L;

			};
			// }

		}
		return dud;
	}

	private void prep() {
		int[] autosomalIndices = proj.getAutosomalMarkerIndices();
		double[] autosomalBafs = Array.subArray(bafs, autosomalIndices);
		this.gd = prepareGaussMixture(autosomalBafs, .33, .66);
		if (baseLine < 0) {
			this.baseLine = 0;
			for (int j = 0; j < gd.distributions().length; j++) {
				baseLine += (double) gd.distributions()[j].probability(means[j] + nullSigma * Math.sqrt(variances[j])) * Math.sqrt(variances[j]);
			}
		}
		reportDynRange();
	}

	private void reportDynRange() {
		if (verbose) {
			double minDelta = Math.sqrt(gd.distributions()[1].variance());
			proj.getLog().reportTimeInfo("Min proportion Disomy detection ~=" + MosaicismQuant.getDisomyF(minDelta));
			double maxDelta = .5 - nullSigma * Math.sqrt(gd.distributions()[2].variance());// B allele generally greater variance
			proj.getLog().reportTimeInfo("Max proportion Disomy detection ~=" + MosaicismQuant.getDisomyF(maxDelta));
			proj.getLog().reportTimeInfo("Null conf ~=" + baseLine);

		}
	}

	private GaussianMixtureDistribution prepareGaussMixture(double[] autosomalBafs, double r1, double r2) {
		if (means == null) {

			this.means = new double[3];
			this.variances = new double[3];
			double[] zero_tsMeanVar = getMeanVar(autosomalBafs, 0, r1);
			double[] t_tsMeanVar = getMeanVar(autosomalBafs, r1, r2);
			double[] t_sMeanVar = getMeanVar(autosomalBafs, r2, 1);
			means[0] = zero_tsMeanVar[0];
			variances[0] = zero_tsMeanVar[1];
			means[1] = t_tsMeanVar[0];
			variances[1] = t_tsMeanVar[1];
			means[2] = t_sMeanVar[0];
			variances[2] = t_sMeanVar[1];
		}

		double[] props = new double[3];
		Arrays.fill(props, (double) 1 / 3);
		GaussianMixtureDistribution gd = new GaussianMixtureDistribution(means, variances, props);
		return gd;
	}

	private static double[] getMeanVar(double[] autosomalBafs, double r1, double r2) {
		double[] sub = Array.getValuesBetween(autosomalBafs, r1, r2, true);
		double[] meanVar = new double[2];
		meanVar[0] = Array.mean(sub);
		meanVar[1] = Math.pow(Array.stdev(sub), 2);
		return meanVar;
	}

	private static void test() {
		Project proj = new Project("C:/workspace/Genvisis/projects/OSv2_hg19.properties", false);

		ArrayList<CNVariant> all = new ArrayList<CNVariant>();
		int movingFactor = 50;
		// String[] samples = new String[] { "7355066051_R03C01", "7330686030_R02C01", "7159911135_R01C02" };
		String[] samples = new String[] { "7355066051_R03C01" };

		for (int i = 0; i < samples.length; i++) {

			String sample = samples[i];
			Sample samp = proj.getFullSampleFromRandomAccessFile(sample);
			MarkerSet markerSet = proj.getMarkerSet();
			MosaicBuilder mosBuilder = new MosaicBuilder();

			MosaicismDetect md = mosBuilder.build(proj, sample, markerSet, Array.toDoubleArray(samp.getBAFs()));
			int[][] te = markerSet.getIndicesByChr();
			for (int j = 0; j < te.length; j++) {
				if (te[j].length > 0 && j < 23) {
					proj.getLog().reportTimeInfo("Calling chr " + j + " for sample " + i);
					LocusSet<CNVariant> hi = md.callMosaic(new Segment((byte) j, 0, markerSet.getPositions()[te[j][te[j].length - 1]] + 10));
					for (int k = 0; k < hi.getLoci().length; k++) {
						System.out.println(hi.getLoci()[k].toPlinkFormat());
					}
					hi.addAll(all);
					ArrayList<CNVariant> tmp = new ArrayList<CNVariant>();
					for (int k = 0; k < all.size(); k++) {
						if (all.get(k).getNumMarkers() > movingFactor) {
							tmp.add(all.get(k));
						}
					}
					all = tmp;

				}

			}
			// }
		}

		LocusSet<CNVariant> mos = new LocusSet<CNVariant>(all.toArray(new CNVariant[all.size()]), true, proj.getLog()) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};

		mos.writeRegions(proj.PROJECT_DIRECTORY.getValue() + "TestMosaic/mos.cnvs", TO_STRING_TYPE.REGULAR, true, proj.getLog());
		// proj.CNV_FILENAMES.addValue(proj.PROJECT_DIRECTORY.getValue() + "TestMosaic/mos.cnvs");

		proj.saveProperties();

	}

	public static void main(String[] args) {

		test();
	}

	public static class MosaicBuilder {
		// init to default params...
		private int movingFactor = DEFAULT_MOVING_FACTOR;
		private double nullSigma = DEFAULT_NULL_SIGMA;
		private boolean verbose = false;
		private int[][] indicesByChr = null;
		private double baseLine = DEFAULT_BASELINE;
		private double[] means = null;
		private double[] variances = null;
		private double minPercentStates = DEFAULT_MIN_PERCENT_STATES;

		/**
		 * @param movingFactor
		 *            set the moving average length for the smoothing step
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
		 * @param nullSigma
		 *            points within this many standard deviations of either of the three distributions will not be counted
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
		 * @param verbose
		 *            verbose output
		 * @return
		 */
		public MosaicBuilder verbose(boolean verbose) {
			this.verbose = verbose;
			return this;
		}

		/**
		 * @param indicesByChr
		 *            project marker indices by chromosome
		 * @return
		 */
		public MosaicBuilder indicesByChr(int[][] indicesByChr) {
			this.indicesByChr = indicesByChr;
			return this;
		}

		/**
		 * @param baseLine
		 *            threshold for calling mosiac regions
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
		 * @param means
		 *            means of the 3 distributions
		 * @param variances
		 *            variances of the 3 distributions;
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
		 * @param minPercentStates
		 *            this percent of all markers in a region of interst must be outside of nullsigma
		 * @return
		 */
		public MosaicBuilder minPercentStates(double minPercentStates) {
			this.minPercentStates = minPercentStates;
			if (minPercentStates <= 0) {
				throw new IllegalArgumentException("minPercentStates must be positive");
			}
			return this;
		}

		public MosaicismDetect build(Project proj, String sample, MarkerSet markerSet, double[] bafs) {
			return new MosaicismDetect(this, proj, sample, markerSet, bafs);
		}
	}

	private MosaicismDetect(MosaicBuilder builder, Project proj, String sample, MarkerSet markerSet, double[] bafs) {
		this.proj = proj;
		this.sample = sample;
		this.markerSet = markerSet;
		this.bafs = bafs;
		this.movingFactor = builder.movingFactor;
		this.nullSigma = builder.nullSigma;
		this.verbose = builder.verbose;
		this.indicesByChr = builder.indicesByChr == null ? markerSet.getIndicesByChr() : builder.indicesByChr;
		this.baseLine = builder.baseLine;
		this.means = builder.means;
		this.variances = builder.variances;
		this.minPercentStates = builder.minPercentStates;
		if (bafs.length != markerSet.getMarkerNames().length) {
			throw new IllegalArgumentException("Internal error, bafs must be present for entire array, fill with NaN if neccesary");
		}
		prep();
	}
}
