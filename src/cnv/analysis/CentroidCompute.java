package cnv.analysis;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.Callable;

import stats.Maths;
import common.Array;
import common.Files;
import common.Logger;
import common.WorkerTrain;
import common.WorkerTrain.Producer;
import common.ext;
import cnv.filesys.Centroids;
import cnv.filesys.ClusterFilterCollection;
import cnv.filesys.MarkerData;
import cnv.filesys.Project;
import cnv.filesys.Project.ARRAY;
import cnv.manage.MDL;
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

	public ClusterFilterCollection getClusterFilterCollection() {
		return clusterFilterCollection;
	}

	public float[][] getCentroid() {
		return centroid;
	}

	public double getGcThreshold() {
		return gcThreshold;
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
			return markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerData.getMarkerName(), (float) gcThreshold, log);
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
		Logger log;

		log = proj.getLog();
		log.report(proj.getPropertyFilename() + "\t" + markers.length);
		String min = "";
		double mincor = .90;
		for (int i = 0; i < markers.length; i++) {
			MarkerData markerData = markerDataLoader.requestMarkerData(i);
			CentroidCompute cent = new CentroidCompute(markerData, null, null, false, 1, 0, null, true, proj.getLog());
			cent.computeCentroid();
			double cor = markerData.compareLRRs(cent.getCentroid(), proj.getLog())[0];
			if (markerData.getChr() < 23 && cor < mincor) {
				// mincor = markerData.compareLRRs(cent.getCentroid())[0];
				min += markers[i] + "\t" + markerData.getChr() + "\t" + markerData.getPosition() + "\t" + cor + "\n";
			}
			if (i % 1000 == 0) {
				log.report((i + 1) + "/" + markers.length);
			}
			// System.out.println("LRR:" + Array.toStr(markerData.compareLRRs(cent.getCentroid())) + "\t" + markers[i]);
			// System.out.println(cent.computeBAF()[0] + "\tBAF:\t" + Array.toStr(stats.Correlation.Pearson(Array.toDoubleArray(cent.computeBAF()), Array.toDoubleArray(markerData.getBAFs()))));
			markerDataLoader.releaseIndex(i);
		}
		Files.write(min, proj.PROJECT_DIRECTORY.getValue() + "badcorrel.txt");
		log.report("MINMarker:\t" + min + "\tMINCorrel\t" + mincor);
	}

	public static void main(String[] args) {
		try {
			Project proj = new Project(null, false);
			// test(proj);
			test2(proj);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * @author Builder for {@link CentroidCompute}
	 *
	 */
	public static class CentroidBuilder {
		private boolean[] samplesToUse = null;
		private boolean intensityOnly = false;
		private boolean medianCenter = true;
		private int[] sampleSex = null;
		private double missingnessThreshold = 1;
		private double gcThreshold = 0;
		private ClusterFilterCollection clusterFilterCollection = null;

		/**
		 * @param samplesToUse
		 *            only cluster on these samples
		 * @return
		 */
		public CentroidBuilder samplesToUse(boolean[] samplesToUse) {
			this.samplesToUse = samplesToUse;
			return this;
		}

		/**
		 * @param intensityOnly
		 *            this is an intensity only cluster
		 * @return
		 */
		public CentroidBuilder intensityOnly(boolean intensityOnly) {
			this.intensityOnly = intensityOnly;
			return this;
		}

		/**
		 * @param medianCenter
		 *            cluster about the median value
		 * @return
		 */
		public CentroidBuilder medianCenter(boolean medianCenter) {
			this.medianCenter = medianCenter;
			return this;
		}

		/**
		 * @param sampleSex
		 *            sex of each individual
		 * @return
		 */
		public CentroidBuilder sampleSex(int[] sampleSex) {
			this.sampleSex = sampleSex;
			return this;
		}

		/**
		 * @param missingnessThreshold
		 *            call rate for clustering
		 * @return
		 */
		public CentroidBuilder missingnessThreshold(double missingnessThreshold) {
			this.missingnessThreshold = missingnessThreshold;
			return this;
		}

		/**
		 * @param gcThreshold
		 *            gc threshold for clustering
		 * @return
		 */
		public CentroidBuilder gcThreshold(double gcThreshold) {
			this.gcThreshold = gcThreshold;
			return this;
		}

		/**
		 * @param clusterFilterCollection
		 * @return
		 */
		public CentroidBuilder clusterFilterCollection(ClusterFilterCollection clusterFilterCollection) {
			this.clusterFilterCollection = clusterFilterCollection;
			return this;
		}

		public CentroidCompute build(MarkerData markerData, Logger log) {
			return new CentroidCompute(this, markerData, log);
		}
	}

	private CentroidCompute(CentroidBuilder builder, MarkerData markerData, Logger log) {
		this(markerData, builder.sampleSex, builder.samplesToUse, builder.intensityOnly, builder.missingnessThreshold, builder.gcThreshold, builder.clusterFilterCollection, builder.medianCenter, log);
	}

	public static Centroids computeAndDumpCentroids(Project proj) {
		CentroidBuilder builder = new CentroidBuilder();
		return CentroidCompute.computeAndDumpCentroids(proj, proj.CUSTOM_CENTROIDS_FILENAME.getValue(true, false), builder, proj.NUM_THREADS.getValue(), 1);
	}

	public static Centroids computeAndDumpCentroids(Project proj, String fullpathToCentFile, CentroidBuilder builder, int numCentThreads, int numDecompressThreads) {
		return computeAndDumpCentroids(proj, new String[] { fullpathToCentFile }, new CentroidBuilder[] { builder }, numCentThreads, numDecompressThreads)[0];
	}

	/**
	 * @param proj
	 *            Current project
	 * @param fullpathToCentFile
	 *            where the centroids will be stored
	 * @param builder
	 *            a centroid builder to use
	 * @param numCentThreads
	 *            number of threads for computing centroids
	 * @param numDecompressThreads
	 *            number of threads for decompressing marker data
	 */
	public static Centroids[] computeAndDumpCentroids(Project proj, String[] fullpathToCentFiles, CentroidBuilder[] builders, int numCentThreads, int numDecompressThreads) {
		if (builders.length != fullpathToCentFiles.length) {
			throw new IllegalArgumentException("Each centroid builder must have a file name associated");
		}
		String[] markers = proj.getMarkerNames();
		float[][][][] centroids = new float[builders.length][markers.length][][];
		CentroidProducer producer = new CentroidProducer(proj, markers, builders, numDecompressThreads);
		WorkerTrain<CentroidCompute[]> train = new WorkerTrain<CentroidCompute[]>(producer, numCentThreads, 10, proj.getLog());
		int index = 0;
		while (train.hasNext()) {
			CentroidCompute[] centroidCompute = train.next();
			if (index % 10000 == 0) {
				proj.getLog().reportTimeInfo(index + " of " + markers.length);
			}
			for (int i = 0; i < centroidCompute.length; i++) {
				centroids[i][index] = centroidCompute[i].getCentroid();

			}
			index++;
		}
		train.shutdown();
		Centroids[] computed = new Centroids[builders.length];
		for (int i = 0; i < centroids.length; i++) {
			Centroids cent = new Centroids(centroids[i], proj.getMarkerSet().getFingerprint());
			cent.serialize(fullpathToCentFiles[i]);
		}
		return computed;
	}

	private static class CentroidProducer implements Producer<CentroidCompute[]> {
		private final Project proj;
		private String[] markers;
		private final CentroidBuilder[] builders;
		private final MDL mdl;
		private int count;

		public CentroidProducer(Project proj, String[] markers, CentroidBuilder[] builders, int numDecompressThreads) {
			super();
			this.proj = proj;
			this.markers = markers;
			this.builders = builders;
			this.mdl = new MDL(proj, proj.getMarkerSet(), markers, numDecompressThreads, 100);

			this.count = 0;
		}

		public void shutdown() {
			mdl.shutdown();

		}

		@Override
		public boolean hasNext() {
			boolean hasNext = count < markers.length;
			if (!hasNext) {
				mdl.shutdown();
			}
			return hasNext;
		}

		@Override
		public Callable<CentroidCompute[]> next() {
			final MarkerData markerData = mdl.next();
			Callable<CentroidCompute[]> compute = new Callable<CentroidCompute[]>() {

				@Override
				public CentroidCompute[] call() throws Exception {
					CentroidCompute[] computed = new CentroidCompute[builders.length];
					for (int i = 0; i < computed.length; i++) {

						CentroidCompute centroidCompute = builders[i].build(markerData, proj.getLog());
						ARRAY array = proj.getArrayType();
						if ((array == ARRAY.ILLUMINA || array == ARRAY.NGS) && array.isCNOnly(markerData.getMarkerName())) {
							setFakeAB(markerData, centroidCompute, centroidCompute.getClusterFilterCollection(), 0, proj.getLog());
						} else if (array.isCNOnly(markerData.getMarkerName())) {
							if (array != ARRAY.AFFY_GW6_CN && array != ARRAY.AFFY_GW6) {
								proj.getLog().reportTimeError("Intenstity only centroid designed for Affymetrix only");
								return null;
							}
						}
						centroidCompute.computeCentroid();
						computed[i] = centroidCompute;
					}

					return computed;
				}
			};
			count++;
			return compute;
		}

		@Override
		public void remove() {
		}
	}

	public static void test2(Project proj) {
		String[] markers = proj.getMarkerNames();
		long time = System.currentTimeMillis();
		CentroidBuilder builder = new CentroidBuilder();
		CentroidProducer producer = new CentroidProducer(proj, markers, new CentroidBuilder[] { builder }, 2);
		WorkerTrain<CentroidCompute[]> train = new WorkerTrain<CentroidCompute[]>(producer, 6, 100, proj.getLog());
		int index = 0;
		while (train.hasNext()) {
			CentroidCompute[] centroidCompute = train.next();
			if (!centroidCompute[0].getMarkerData().getMarkerName().equals(markers[index])) {
				proj.getLog().reportError(centroidCompute[0].getMarkerData().getMarkerName() + "\t" + markers[index]);
			}
			index++;
		}
		train.shutdown();
		proj.getLog().reportTimeInfo("Took" + ext.getTimeElapsed(time));
	}

	/**
	 * We use this in the case of intensity only probesets...We assign all genotypes to be the same
	 * 
	 */
	public static void setFakeAB(MarkerData markerData, CentroidCompute centroid, ClusterFilterCollection clusterFilterCollection, float gcThreshold, Logger log) {
		byte[] fakeAB = new byte[centroid.getClustGenotypes().length];
		byte[] clustAB = markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerData.getMarkerName(), gcThreshold, log);
		int[] counts = new int[4];
		for (int i = 0; i < clustAB.length; i++) {
			counts[clustAB[i] + 1]++;
		}
		if (counts[0] == clustAB.length) {
			byte tmpCluster = Array.mean(markerData.getXs()) > Array.mean(markerData.getYs()) ? (byte) 0 : (byte) 1;
			Arrays.fill(fakeAB, (byte) tmpCluster);
			centroid.setAlternateGenotypes(fakeAB);

		} else {
			centroid.setAlternateGenotypes(clustAB);
		}
		// assign all to this cluster
	}
}

// private static class CentroidComputeWorker implements Callable<CentroidCompute> {
// private Builder builder;
// private MarkerData markerData;
// private Logger log;
//
// public CentroidComputeWorker(Builder builder, MarkerData markerData, Logger log) {
// super();
// this.builder = builder;
// this.markerData = markerData;
// this.log = log;
// }
//
// @Override
// public CentroidCompute call() throws Exception {
// CentroidCompute centroidCompute = builder.build(markerData, log);
// log.reportTimeInfo("Computing a centroid on " + Thread.currentThread().getName());
// centroidCompute.computeCentroid();
// return centroidCompute;
// }
//
// }
