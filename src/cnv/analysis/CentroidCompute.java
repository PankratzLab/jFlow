package cnv.analysis;

import java.util.ArrayList;

import stats.Maths;
import common.Array;
import common.Files;
import common.Logger;
import cnv.filesys.Centroids;
import cnv.filesys.ClusterFilterCollection;
import cnv.filesys.MarkerData;
import cnv.filesys.Project;
import cnv.manage.MarkerDataLoader;

/**
 * A class for centroid related computations for a single {@link MarkerData}
 * <p>
 * Currently handles sex chromosomes, intensity only cluster, genotype based clustering, etc...
 * <p>
 * Warning - this stores the markerData, any use of this class should be weary of this fact
 * <p>
 * 
 * TODO, this class should eventually replace the computation portions found in Centroids, we can leave Centroids as a container and storage.
 */
public class CentroidCompute {
	private MarkerData markerData;
	private boolean[] samplesToUse;
	private boolean intensityOnly, sexSpecific, failed, hasCentroid, medianCenter;
	private double[] centerThetas, centerRs;
	private int[] counts, sampleSex;
	private double missingnessThreshold, gcThreshold;
	private byte[] alternateGenotypes;
	private float[][] centroid; // centroid[genotype][centerTheta,centerR]
	private ClusterFilterCollection clusterFilterCollection;
	private Logger log;

	/**
	 * 
	 * @param markerData
	 *            the markerData to use for clustering
	 * @param sampleSex
	 *            (optional, can be null)an array representing sample sex, 1=male,2=female. Only utilized for clustering on sex chromosomes
	 * @param samplesToUse
	 *            (optional, can be null)an array representing samples to use for the clustering
	 * @param intensityOnly
	 *            perform an intensity only clustering, useful for cnv-type probesets, we do not auto-detect cnv-probesets here
	 * @param missingnessThreshold
	 *            threshold for marker call rate: set to 1 for no filtering
	 * @param confThreshold
	 *            threshold for call confidence: set to 0 for no filtering
	 * @param clusterFilterCollection
	 *            (optional, can be null) if provided, these cluster filters will be applied
	 * @param medianCenter
	 *            computes the centroids using median values as opposed to average values.
	 * @param log
	 */

	public CentroidCompute(MarkerData markerData, int[] sampleSex, boolean[] samplesToUse, boolean intensityOnly, double missingnessThreshold, double confThreshold, ClusterFilterCollection clusterFilterCollection, boolean medianCenter, Logger log) {
		super();
		this.failed = false;
		this.hasCentroid = false;
		this.markerData = markerData;
		this.sexSpecific = sexSpecific(markerData.getChr());
		this.samplesToUse = samplesToUse;
		this.sampleSex = sampleSex;
		this.clusterFilterCollection = clusterFilterCollection;
		this.intensityOnly = intensityOnly;
		this.missingnessThreshold = missingnessThreshold;
		this.gcThreshold = confThreshold;
		this.medianCenter = medianCenter;
		this.centerThetas = new double[5];// 0=all,1=missing,2=AA,3=AB,4=BB
		this.centerRs = new double[5];
		this.counts = new int[5];
		this.centroid = null;
		this.alternateGenotypes = null;
		this.log = log;
		checkData();// fail early if mismatched array sizes
	}

	public void setCentroid(float[][] centroid) {
		this.centroid = centroid;
		hasCentroid = true;
	}

	public float[][] getCentroid() {
		return centroid;
	}

	public MarkerData getMarkerData() {
		return markerData;
	}

	/**
	 * Nice if you want to quickly recompute with new Xs and Ys, but all other parameters the same
	 */
	public void setMarkerData(MarkerData markerData) {
		this.markerData = markerData;
	}

	public boolean[] getSamplesToUse() {
		return samplesToUse;
	}

	public void setIntensityOnly(boolean intensityOnly) {
		this.intensityOnly = intensityOnly;
	}

	/**
	 * We favor returning {@link #alternateGenotypes} if they have been explicitly set.
	 * <p>
	 * Next we favor the {@link #clusterFilterCollection} genotypes.
	 * <P>
	 * Lastly we take the {@link #markerData} original genotypes
	 */
	public byte[] getClustGenotypes() {
		if (alternateGenotypes != null) {
			return alternateGenotypes;
		}
		if (clusterFilterCollection == null) {
			return markerData.getAbGenotypes();
		} else {
			return markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerData.getMarkerName(), (float) gcThreshold);
		}
	}

	public byte[] getAlternateGenotypes() {
		return alternateGenotypes;
	}

	public void setAlternateGenotypes(byte[] alternateGenotypes) {
		this.alternateGenotypes = alternateGenotypes;
	}

	public boolean failed() {
		return failed;
	}

	public float[] getRecomputedLRR() {
		if (!hasCentroid) {// we can also set a new centroid
			computeCentroid();// hasCentroid turns true if success
		}
		if (hasCentroid && !failed) {
			return computeLRR();
		} else {
			return null;
		}
	}

	public float[] getRecomputedBAF() {
		if (!hasCentroid) {
			computeCentroid();
		}
		if (hasCentroid && !failed) {
			return computeBAF();
		} else {
			return null;
		}

	}

	public void init() {
		this.centerThetas = new double[5];// 0=all,1=missing,2=AA,3=AB,4=BB
		this.centerRs = new double[5];
		this.counts = new int[5];
		hasCentroid = false;// just to make sure its not passed on
	}

	public void computeCentroid() {
		computeCentroid(false);
	}

	public void computeCentroid(boolean forceAuto) {
		init();
		computeCenters(forceAuto);
		if (!failed) {
			if (intensityOnly) {
				assignIntensityCentroid();
			} else {
				assignGenotypeCentroid();
			}
			if (!failed) {
				hasCentroid = true;
			}
		} else {
			log.reportError("Error - cannot compute centroid for " + markerData.getMarkerName());
		}
	}

	/**
	 * @return recomputed lrrs, or null if they could not be computed
	 */
	private float[] computeLRR() {
		float[] lrrs = null;
		if (checkCanCompute()) {
			lrrs = new float[markerData.getXs().length];
			if (!intensityOnly) {
				float[] thetas = markerData.getThetas();
				float[] rs = markerData.getRs();
				for (int i = 0; i < thetas.length; i++) {
					lrrs[i] = Centroids.calcLRR(thetas[i], rs[i], centroid);
				}
			} else {
				float[] xs = markerData.getXs();
				for (int i = 0; i < xs.length; i++) {
					lrrs[i] = (float) (Maths.log2(xs[i]) - centroid[0][1]);
				}
			}
		}
		return lrrs;
	}

	/**
	 * @return recomputed bafs, or null if they could not be computed
	 */
	private float[] computeBAF() {
		float[] baf = null;
		if (checkCanCompute()) {
			baf = new float[markerData.getXs().length];
			if (!intensityOnly) {
				float[] thetas = markerData.getThetas();
				for (int i = 0; i < thetas.length; i++) {
					baf[i] = Centroids.calcBAF(thetas[i], centroid);
				}
			} else {
				float[] xs = markerData.getXs();
				for (int i = 0; i < xs.length; i++) {
					baf[i] = 1;
				}
			}
		}
		return baf;
	}

	private void checkData() {
		int check = markerData.getThetas().length;
		byte chr = markerData.getChr();
		if (sampleSex == null && sexSpecific) {
			// log.report("Warning - marker " + markerData.getMarkerName() + " is on chromosome " + chr + " and sample sex was not provided, cluster centers may be inaccurate");
		} else if (sampleSex != null && sampleSex.length != check) {
			failed = true;
			log.reportError("Error - mismatched number of samples for data's length");
		}
		if (samplesToUse != null && samplesToUse.length != check) {
			failed = true;
			log.reportError("Error - mismatched number of samples to use for data's length");
		}
		if (sampleSex != null && sexSpecific && chr == 23 && Array.countIf(sampleSex, 2) == 0) {
			log.reportError("Error - sample sex was defined, but no females were found; cannot properly cluster chr X ");
			sampleSex = null;
		}
		if (sampleSex != null && sexSpecific && chr == 24 && Array.countIf(sampleSex, 1) == 0) {
			log.reportError("Error - sample sex was defined, but no males were found; cannot properly cluster chr Y");
			sampleSex = null;
		}
	}

	private boolean checkCanCompute() {
		if (!hasCentroid || centroid == null) {
			log.reportError("Error - do not hava a centroid for " + markerData.getMarkerName() + ", cannot compute LRR/BAF");
			return false;
		} else {
			return true;
		}
	}

	private void assignGenotypeCentroid() {
		this.centroid = assignRegularCentroid(centerThetas, centerRs, counts, missingnessThreshold);// we first assign all the centroid values, which may have NaNs present from 0 counts
		if (counts[2] == 0) {
			// Estimating the AA cluster if needed and possible (from the AB cluster)
			if (counts[3] > 0) {
				centroid[0] = new float[] { (float) (centerThetas[3] - 0.3), (float) (centerRs[3]) };
			} else {
				centroid[0] = null;
			}
		}
		if (counts[3] == 0) {
			// Estimating the AB cluster if needed and possible (from the AA/BB cluster)
			if (counts[2] > 0 && counts[4] > 0) {
				centroid[1] = new float[] { (float) ((centerThetas[2] + centerThetas[4]) / 2), getAltRs(centerRs, counts, 2, 4) };
			} else if (counts[2] > 0 || counts[4] > 0) {
				centroid[1] = new float[] { (float) (0.5), getAltRs(centerRs, counts, 2, 4) };// default to middle, and take any r value available
			} else {
				centroid[1] = null;
			}
		}
		if (counts[4] == 0) {
			// Estimating the BB cluster if needed and possible (from the AB cluster)
			if (counts[3] > 0) {
				centroid[2] = new float[] { (float) (centerThetas[3] + 0.3), (float) (centerRs[3]) };
			} else {
				centroid[2] = null;
			}
		}
		centroid = assignForMissingnessThreshold(centroid, centerThetas, centerRs, counts, missingnessThreshold);// if missingnessThreshold is exceeded, we set a modified centroid
	}

	private void assignIntensityCentroid() {
		this.centroid = new float[3][2];
		centroid[0] = new float[] { (float) centerThetas[2], (float) centerRs[2] };// median/mean of the log base 2 values
		centroid[1] = new float[] { (float) centerThetas[3], (float) centerRs[3] };// mean/median of the log base 2 values
		centroid[2] = new float[] { (float) centerThetas[4], (float) centerRs[4] };// median/mean of the actual X values
	}

	private void computeCenters(boolean forceAuto) {
		if (!failed) {
			if (intensityOnly) {
				computeIntensityCenters(forceAuto);
			} else {
				computeGenotypeCenters(forceAuto);
			}
		}
	}

	private void computeIntensityCenters(boolean forceAuto) {
		byte chr = markerData.getChr();
		float[] xs = markerData.getXs();
		if (xs == null) {
			log.reportError("Error - no X values found, cannot compute genotype centers");
			failed = true;
		} else {
			ArrayList<Double> filteredXs = new ArrayList<Double>(xs.length);
			for (int i = 0; i < xs.length; i++) {
				if ((forceAuto || checkSex(chr, i)) && (samplesToUse == null || samplesToUse[i]) && !Float.isNaN(xs[i])) {
					filteredXs.add((double) xs[i]);
				}
			}

			double[] xsFilt = Array.toDoubleArray(filteredXs);
			double[] xsFiltLog2 = Array.log2(xsFilt);

			centerRs[0] = Double.NaN;
			centerRs[1] = Double.NaN;
			centerRs[2] = (medianCenter ? Array.median(xsFiltLog2) : Array.mean(xsFiltLog2));
			centerRs[3] = (medianCenter ? Array.mean(xsFiltLog2) : Array.median(xsFiltLog2));
			centerRs[4] = (medianCenter ? Array.median(xsFilt) : Array.mean(xsFilt));

			centerThetas[0] = 0.5;
			centerThetas[1] = 0.5;
			centerThetas[2] = 0.5;
			centerThetas[3] = 0.5;
			centerThetas[4] = 0.5;

			counts[0]++;
			// counts[1] intentionally left to 0, can never really have missing data
			counts[2]++;
			counts[3]++;
			counts[4]++;
		}

	}

	/**
	 * We compute the center for thetas and Rs according to samplesToUse (if available), sex (if available and required), confidence threshold, and clusterfilters.
	 */
	private void computeGenotypeCenters(boolean forceAuto) {

		byte chr = markerData.getChr();
		float[] thetas = markerData.getThetas();
		float[] rs = markerData.getRs();
		float[] confs = markerData.getGCs();
		byte[] genotypes;
		ArraySpecialLists slThetas = new ArraySpecialLists(centerThetas.length, thetas.length); // only used for median
		ArraySpecialLists slRs = new ArraySpecialLists(centerRs.length, rs.length);// only used for median

		genotypes = getClustGenotypes();
		failed = checkGenoClusterMarkerData(thetas, rs, confs, genotypes, gcThreshold, log);
		if (!failed) {
			for (int i = 0; i < genotypes.length; i++) {
				boolean check1 = forceAuto || checkSex(chr, i);
				boolean check2 = (samplesToUse == null || samplesToUse[i]);
				boolean check3 = useMarker(thetas[i], rs[i], confs == null ? 0 : confs[i], gcThreshold);

				// System.out.println(i + " - s:" + sampleSex[i] + "->" + check1 + " ; " + samplesToUse[i]);

				if (check1 && check2 && check3) {
					counts[0]++;
					counts[genotypes[i] + 2]++;
					if (medianCenter) {// a bit slower
						slThetas.addTo(0, thetas[i]);
						slThetas.addTo(genotypes[i] + 2, thetas[i]);
						slRs.addTo(0, rs[i]);
						slRs.addTo(genotypes[i] + 2, rs[i]);
					} else {
						centerThetas[0] += thetas[i];
						centerRs[0] += rs[i];
						centerThetas[genotypes[i] + 2] += thetas[i];
						centerRs[genotypes[i] + 2] += rs[i];
					}
				}
			}
			for (int k = 0; k < 5; k++) {
				if (counts[k] > 0) {
					if (medianCenter) {
						centerThetas[k] = slThetas.getMedianAt(k);
						centerRs[k] = slRs.getMedianAt(k);

					} else {
						centerThetas[k] /= counts[k];
						centerRs[k] /= counts[k];
					}
				} else {
					centerThetas[k] = Double.NaN;
					centerRs[k] = Double.NaN;
				}
			}
		} else {
			log.reportError("Error - cannot compute cluster for " + markerData.getMarkerName());
		}
	}

	/**
	 * A data check for genotype based clustering
	 */
	private static boolean checkGenoClusterMarkerData(float[] thetas, float[] rs, float[] confs, byte[] genotypes, double gcThreshold, Logger log) {
		boolean failed = false;
		if (thetas == null) {
			log.reportError("Error - no theta values found, cannot compute genotype centers");
			failed = true;
		}
		if (rs == null) {
			log.reportError("Error - no R values found, cannot compute genotype centers");
			failed = true;
		}
		if (confs == null && gcThreshold > 0) {
			log.reportError("Error - no confidence values found and a confidence threshold was specified greater than 0 (" + gcThreshold + ") cannot compute genotype centers");
			failed = true;
		}
		if (genotypes == null) {
			log.reportError("Error - no genotype values found, cannot compute genotype centers");
			failed = true;
		}
		return failed;
	}

	/**
	 * Essentially just an array of ArrayLists
	 *
	 */
	private static class ArraySpecialLists {
		private ArraySpecialList[] arraySpecialLists;

		public ArraySpecialLists(int num, int initCapacity) {
			this.arraySpecialLists = new ArraySpecialList[num];
			init(initCapacity);
		}

		private void init(int initCapacity) {
			for (int i = 0; i < arraySpecialLists.length; i++) {
				arraySpecialLists[i] = new ArraySpecialList(initCapacity);
			}
		}

		public void addTo(int index, double d) {
			arraySpecialLists[index].add(d);
		}

		public double[] getAt(int index) {
			return Array.toDoubleArray(arraySpecialLists[index]);
		}

		public double getMedianAt(int index) {
			return Array.median(getAt(index));
		}
	}

	/**
	 * Just an ArrayList that can be placed in an array
	 *
	 */
	private static class ArraySpecialList extends ArrayList<Double> {
		private static final long serialVersionUID = 1L;

		public ArraySpecialList(int initCapacity) {
			super(initCapacity);
		}
	}

	private boolean checkSex(byte chr, int i) {
		return sampleSex == null || !sexSpecific || useSampleSex(sampleSex[i], chr);
	}

	private static boolean useMarker(float theta, float r, float conf, double confThreshold) {
		return !Float.isNaN(theta) && !Float.isNaN(r) && !Float.isNaN(conf) && (confThreshold == 0 || conf > confThreshold);
	}

	/**
	 * @param sampleSex
	 *            int 1=male, 2=female
	 * @param chr
	 * @return whether the sex is appropriate for the chromosome
	 */
	private static boolean useSampleSex(int sampleSex, byte chr) {
		boolean use = true;
		if (chr == 23 && sampleSex != 2) {// only cluster females on chromosome 23
			use = false;
		}
		if (chr == 24 && sampleSex != 1) {// only cluster males on chromosome 24
			use = false;
		}
		return use;
	}

	/**
	 * @param chr
	 * @return true if chr is a sex chromosome
	 */
	private static boolean sexSpecific(byte chr) {
		return (chr == 23 || chr == 24);
	}

	/**
	 * In the case where one of the genotype clusters has 0 calls in the samples used to compute them, we attempt to assign an alternate R value
	 * <p>
	 * This could be important if the number of samples used is small and applying to a large subset
	 */
	private static float getAltRs(double[] meanRs, int[] counts, int primaryAlt, int secondaryAlt) {
		float altR;
		if (counts[primaryAlt] > 0) {
			altR = (float) meanRs[primaryAlt];
		} else if (counts[secondaryAlt] > 0) {
			altR = (float) meanRs[secondaryAlt];
		} else {
			altR = Float.NaN;
		}
		return altR;
	}

	private static float[][] assignRegularCentroid(double[] meanThetas, double[] meanRs, int[] counts, double missingnessThreshold) {
		float[][] centroid = new float[3][2];
		for (int k = 0; k < 3; k++) {
			centroid[k] = new float[] { (float) meanThetas[k + 2], (float) meanRs[k + 2] };
		}
		return centroid;
	}

	private static float[][] assignForMissingnessThreshold(float[][] centroid, double[] meanThetas, double[] meanRs, int[] counts, double missingnessThreshold) {
		if (counts[1] >= counts[0] * missingnessThreshold) {
			for (int k = 0; k < 3; k++) {
				centroid[k] = new float[] { (float) meanThetas[0], (float) meanRs[0] };
			}
		}
		return centroid;
	}

	/**
	 * Gives a printout of the correlation to existing LRRs and BAFs
	 */
	public static void test(Project proj) {
		String[] markers = proj.getTargetMarkers();
		MarkerDataLoader markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markers);
		System.out.println(proj.getPropertyFilename() + "\t" + markers.length);
		String min = "";
		double mincor = .90;
		for (int i = 0; i < markers.length; i++) {
			MarkerData markerData = markerDataLoader.requestMarkerData(i);
			CentroidCompute cent = new CentroidCompute(markerData, null, null, false, 1, 0, null, true, proj.getLog());
			cent.computeCentroid();
			double cor = markerData.compareLRRs(cent.getCentroid())[0];
			if (markerData.getChr() < 23 && cor < mincor) {
				// mincor = markerData.compareLRRs(cent.getCentroid())[0];
				min += markers[i] + "\t" + markerData.getChr() + "\t" + markerData.getPosition() + "\t" + cor + "\n";
			}
			if (i % 1000 == 0) {
				System.out.println((i + 1) + "/" + markers.length);
			}
			// System.out.println("LRR:" + Array.toStr(markerData.compareLRRs(cent.getCentroid())) + "\t" + markers[i]);
			// System.out.println(cent.computeBAF()[0] + "\tBAF:\t" + Array.toStr(stats.Correlation.Pearson(Array.toDoubleArray(cent.computeBAF()), Array.toDoubleArray(markerData.getBAFs()))));
			markerDataLoader.releaseIndex(i);
		}
		Files.write(min, proj.getProjectDir() + "badcorrel.txt");
		System.out.println("MINMarker:\t" + min + "\tMINCorrel\t" + mincor);
	}

	public static void main(String[] args) {
		try {
			Project proj = new Project(null, false);
			test(proj);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
