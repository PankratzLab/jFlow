package cnv.analysis;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;

import be.ac.ulg.montefiore.run.distributions.GaussianMixtureDistribution;
import common.Array;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.hmm.PennHmm.ViterbiResult;
import cnv.var.CNVariant;
import cnv.var.LocusSet;
import cnv.var.LocusSet.TO_STRING_TYPE;
import cnv.var.MosaicRegion;
import filesys.Segment;
import cnv.var.CNVariant.CNVBuilder;

public class MosaicismDetect implements Iterator<MosaicRegion> {

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

	public MosaicismDetect(Project proj, String sample,int[][] indicesByChr, double[] bafs, MarkerSet markerSet, int movingFactor, double nullSigma, boolean verbose) {
		super();
		this.proj = proj;
		this.sample = sample;
		this.markerSet = markerSet;
		this.indicesByChr =indicesByChr==null? markerSet.getIndicesByChr():indicesByChr;
		this.sample = sample;
		this.bafs = bafs;
		this.movingFactor = movingFactor;
		this.nullSigma = nullSigma;
		this.verbose = verbose;
		if (bafs.length != markerSet.getMarkerNames().length) {
			throw new IllegalArgumentException("Internal error, bafs must be present for entire array, fill with NaN if neccesary");
		}
		prep();
	}

	public int getMovingFactor() {
		return movingFactor;
	}

	public <T extends Segment> LocusSet<CNVariant> callMosaic(T seg) {
		int[] segIndices = markerSet.getIndicesOfMarkersIn(seg, indicesByChr, proj.getLog());
		ArrayList<Integer> evalIndicestmp = new ArrayList<Integer>();
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
					if (j == 0 && Math.abs(baf - means[j]) < 2 * Math.sqrt(variances[j])) {
						p_density[i] = Double.NaN;

					} else if (Math.abs(baf - means[j]) < 2 * Math.sqrt(variances[j])) {
						p_density[i] = Double.NaN;

					}
				}
				double test = Double.isFinite(baf) ? (double) (baf) : 0;
				// System.out.println(test);
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
			ArrayList<Double> p_densityScored = new ArrayList<Double>();
			ArrayList<Integer> mosIndicesTmp = new ArrayList<Integer>();
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
				builder.score(baseLine - Array.mean(scored));
				tmp[i] = builder.build();
			}
			dud = new LocusSet<CNVariant>(tmp, true, proj.getLog()) {

				/**
				 * 
				 */
				private static final long serialVersionUID = 1L;

			};
 
		}
		return dud;
	}

	private void prep() {
		int[] autosomalIndices = proj.getAutosomalMarkerIndices();
		double[] autosomalBafs = Array.subArray(bafs, autosomalIndices);
		this.gd = prepareGaussMixture(autosomalBafs, .33, .66);
		this.baseLine = 0;
		for (int j = 0; j < gd.distributions().length; j++) {
			baseLine += (double) gd.distributions()[j].probability(means[j] + nullSigma * Math.sqrt(variances[j])) * Math.sqrt(variances[j]);
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

		double[] props = new double[3];
		Arrays.fill(props, (double) 1 / 3);
		GaussianMixtureDistribution gd = new GaussianMixtureDistribution(means, variances, props);
		return gd;
	}

	private double[] getMeanVar(double[] autosomalBafs, double r1, double r2) {
		double[] sub = Array.getValuesBetween(autosomalBafs, r1, r2, true);
		double[] meanVar = new double[2];
		meanVar[0] = Array.mean(sub);
		meanVar[1] = Math.pow(Array.stdev(sub), 2);
		return meanVar;
	}

	@Override
	public boolean hasNext() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public MosaicRegion next() {
		// TODO Auto-generated method stub
		return null;
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
			MosaicismDetect md = new MosaicismDetect(proj, sample,null, Array.toDoubleArray(samp.getBAFs()), markerSet, movingFactor, 2, true);
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

}
