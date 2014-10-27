package cnv.analysis.pca;

import java.util.Arrays;
import java.util.Hashtable;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import stats.CrossValidation;
import common.Array;
import common.Logger;
import cnv.analysis.CentroidCompute;
import cnv.filesys.Centroids;
import cnv.filesys.ClusterFilterCollection;
import cnv.filesys.MarkerData;
import cnv.filesys.Project;
import cnv.plots.ScatterPlot;

/**
 * Class to facilitate correcting X/Y intensity data; TODO, intensity only correction, should be easy -> just like a single genotype
 *
 */
public class PrincipalComponentsIntensity extends PrincipalComponentsResiduals {
	private static final String[] CORRECTION_METHODS = { "Assign Missing Genotypes to Clusters (by nearest theta)", "Two stage correction (use existing genotypes to find most predictive model)", "N stage correction (use existing genotypes to find most predictive model)", "N stage correction with residual filtering (use existing genotypes to find most predictive model)" };
	public static final String XY_RETURN = "X_Y";
	public static final String THETA_R_RETURN = "THETA_R";
	public static final String BAF_LRR_RETURN = "BAF_LRR";

	private static final int[] CORRECTION_INTS = { 0, 1, 2, 3 };
	private static final int KILL_INT = -99;
	private static final String[] AFFY_INTENSITY_ONLY_FLAG = { "CN_" };
	private static final String[] ILLUMINA_INTENSITY_ONLY_FLAG = { "cnvi" };
	private static final double MIN_CLUSTER_PERCENT = 0.05;
	private CentroidCompute centroid;
	private boolean fail, svdRegression, verbose;
	private int correctionMethod, nStage, numThreads;
	private String[] samples;
	private float[] genotypeThetaCenters, correctedXFull, correctedYFull;
	private float[][] correctedXCluster, correctedYCluster;
	private int[] genoClusterCounts;
	private int numTotalSamples;
	private double residStandardDeviationFilter;
	private boolean[][] genoSampleClusters;// genotype, sample
	private boolean[] forceThisCluster;

	public PrincipalComponentsIntensity(PrincipalComponentsResiduals principalComponentsResiduals, MarkerData markerData, boolean recomputeLRR, int[] sampleSex, boolean[] samplesToUseCluster, double missingnessThreshold, double confThreshold, ClusterFilterCollection clusterFilterCollection, boolean medianCenter, boolean svdRegression, int correctionMethod, int nStage, double residStandardDeviationFilter, int numThreads, boolean verbose, String output) {
		super(principalComponentsResiduals);// we hijack the loading of the PC file and tracking of samples etc ...
		this.samples = proj.getSamples();
		this.svdRegression = svdRegression;
		this.correctionMethod = correctionMethod;
		this.residStandardDeviationFilter = residStandardDeviationFilter;
		this.verbose = verbose;
		this.nStage = nStage;
		this.fail = (markerData.getXs() == null || markerData.getYs() == null || markerData.getAbGenotypes() == null);
		this.centroid = prepareProperCentroid(markerData, sampleSex, samplesToUseCluster, missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, proj.getLog());
		this.numTotalSamples = markerData.getXs().length;
		this.correctedXFull = new float[numTotalSamples];
		this.correctedYFull = new float[numTotalSamples];
		this.forceThisCluster = new boolean[3];
		this.numThreads = numThreads;
		if (centroid.failed()) {
			fail = true;
		}
		if (invalidMethod(correctionMethod, this.getProj().getLog())) {
			fail = true;
		}
	}

	public MarkerData getMarkerDataUsed() {
		return centroid.getMarkerData();
	}

	public CentroidCompute getCentroidCompute() {
		return centroid;
	}

	public byte[] getGenotypesUsed() {
		return centroid.getClustGenotypes();
	}

	public boolean shouldforceOriginalGenotypes() {// to switch genotypes, we demand that we have at least two valid clusters above min size
		int goodClust = 0;
		for (int i = 0; i < genoClusterCounts.length; i++) {
			if (isClusterBigEnough(genoClusterCounts[i], numTotalSamples)) {
				goodClust++;
				// System.out.println("good Enough geno" + i);
			} else {
				forceThisCluster[i] = true;
			}
		}
		if (goodClust >= 1) {
			// System.out.println("good Enough");
			return false;

		} else {
			return true;
		}

	}

	public void correctXYAt(int atComponent) {
		if (!fail) {
			assignGenotypeClusterStatus(centroid.getClustGenotypes());
			if (hasSomeGenotypes()) {
				this.correctedXCluster = new float[genoSampleClusters.length][];
				this.correctedYCluster = new float[genoSampleClusters.length][];
				if (correctionMethod == CORRECTION_INTS[2] && shouldforceOriginalGenotypes()) {
					correctionMethod = KILL_INT;// we kill it right here
				}
				if ((correctionMethod == CORRECTION_INTS[0] && !fail) || correctionMethod == KILL_INT) {
					correctXYRegression(atComponent);
				} else if (correctionMethod == CORRECTION_INTS[1] && !fail) {
					correctionMethod = KILL_INT;
					estimateNewGenotypes(atComponent);// applies new genotypes to clusters
					correctXYAt(atComponent);// estimate with new genotypes
				} else if (correctionMethod == CORRECTION_INTS[2] && !fail) {
					if (verbose) {
						log.report("Correcting at stage " + nStage);
					}
					nStage--;

					if (nStage >= 1) {// we have more stages to go
						estimateNewGenotypes(atComponent);
					}
					if (nStage == 1 && residStandardDeviationFilter != 0) {// get final genotypes (there will be no more missing genotypes)
						residStandardDeviationFilter = 0;
						estimateNewGenotypes(atComponent);
					}
					if (nStage <= 1) {// final regression, can also be killed if gentypes are identical between iterations
						correctionMethod = KILL_INT;
					}
					// applies new genotypes to clusters
					correctXYAt(atComponent);// estimate with new genotypes

				}
				mergeCorrectedClusters();
			} else {
				fail = true;
				proj.getLog().reportError("Error - correction has failed for marker due to not having any genotypes " + centroid.getMarkerData().getMarkerName() + ", returning orgiginal x, y, and genotype values");
				setOriginal();
			}
		} else {
			proj.getLog().reportError("Error - correction has failed during the centroid computation for marker " + centroid.getMarkerData().getMarkerName() + ", returning orgiginal x, y, and genotype values");
			setOriginal();
		}
	}

	public boolean isFail() {
		return fail;
	}

	private boolean hasSomeGenotypes() {
		for (int i = 0; i < genoClusterCounts.length; i++) {
			if (genoClusterCounts[i] > 0) {
				return true;
			}
		}
		return false;
	}

	private void setOriginal() {
		correctedXFull = centroid.getMarkerData().getXs();
		correctedYFull = centroid.getMarkerData().getYs();
		centroid.setAlternateGenotypes(centroid.getMarkerData().getAbGenotypes());
	}

	public float[] getCorrectedXFull() {
		return correctedXFull;
	}

	public float[] getCorrectedYFull() {
		return correctedYFull;
	}

	private float[][] getCorrectedXY(boolean original) {
		if (original || fail || correctedXFull == null || correctedYFull == null) {
			return new float[][] { centroid.getMarkerData().getXs(), centroid.getMarkerData().getYs() };
		} else {
			return new float[][] { getCorrectedXFull(), getCorrectedYFull() };
		}
	}

	public float[][] getCorrectedIntensity(String type, boolean originalGenotypes) {
		boolean original = false;
		if (fail || correctedXFull == null || correctedYFull == null) {
			original = true;
		}
		if (type.equals(XY_RETURN)) {
			return getCorrectedXY(original);
		}
		byte[] genotypesToUse = originalGenotypes ? centroid.getMarkerData().getAbGenotypes() : getGenotypesUsed();
		MarkerData tmpMarkerData = new MarkerData(centroid.getMarkerData().getMarkerName(), centroid.getMarkerData().getChr(), centroid.getMarkerData().getPosition(), centroid.getMarkerData().getFingerprint(), centroid.getMarkerData().getGCs(), null, null, correctedXFull, correctedYFull, null, null, null, null, genotypesToUse, genotypesToUse);

		if (type.equals(THETA_R_RETURN)) {
			return getCorrectedThetaR(original, tmpMarkerData);
		} else if (type.equals(BAF_LRR_RETURN)) {
			return getCorrectedBAFLRR(original, tmpMarkerData);
		} else {
			log.reportError("Error - invalid request for corrected data");
			return null;
		}

	}

	private float[][] getCorrectedBAFLRR(boolean original, MarkerData tmpMarkerData) {
		if (original) {
			return new float[][] { centroid.getMarkerData().getBAFs(), centroid.getMarkerData().getLRRs() };

		} else {
			centroid.setMarkerData(tmpMarkerData);
			centroid.computeCentroid();
			return new float[][] { centroid.getRecomputedBAF(), centroid.getRecomputedLRR() };
		}
	}

	private float[][] getCorrectedThetaR(boolean original, MarkerData tmpMarkerData) {
		if (original) {
			return new float[][] { centroid.getMarkerData().getThetas(), centroid.getMarkerData().getRs() };
		} else {
			return new float[][] { tmpMarkerData.getThetas(), tmpMarkerData.getRs() };
		}
	}

	private void correctXYRegression(int atComponent) {
		CrossValidation[][] cvals = threadIt(atComponent, Array.toDoubleArray(centroid.getMarkerData().getXs()), Array.toDoubleArray(centroid.getMarkerData().getYs()));
		for (int i = 0; i < genoSampleClusters.length; i++) {
			int clusterComponent = getProperComponent(atComponent, genoClusterCounts[i], i, centroid.getMarkerData().getMarkerName(), verbose, proj.getLog());
			if (validClusterComponent(genoClusterCounts[i], clusterComponent, numTotalSamples)) {
				correctedXCluster[i] = Array.toFloatArray(Array.subArray(cvals[i][0].getResiduals(), genoSampleClusters[i]));
				correctedYCluster[i] = Array.toFloatArray(Array.subArray(cvals[i][1].getResiduals(), genoSampleClusters[i]));
			} else {
				genoSampleClusters[i] = new boolean[genoSampleClusters[i].length];
				Arrays.fill(genoSampleClusters[i], false); // not enough individuals, no corrections
			}
		}
	}

	/**
	 * 
	 * Estimate the genotype based on the nearest cluster to the predicted values (across all regression models tested)
	 */
	private void estimateNewGenotypes(int atComponent) {
		double[] Xs = Array.toDoubleArray(centroid.getMarkerData().getXs());
		double[] Ys = Array.toDoubleArray(centroid.getMarkerData().getYs());
		byte[] estimatedGenotypes = new byte[Xs.length];
		double[][] fullXPredicteds = new double[3][];
		double[][] fullYPredicteds = new double[3][];
		double[][] residX = new double[3][];
		double[][] residY = new double[3][];
		double[] residstdevX = new double[3];
		double[] residstdevY = new double[3];
		CrossValidation[][] cvals = threadIt(atComponent, Xs, Ys);// organaized as genotype,x/y

		for (int i = 0; i < genoSampleClusters.length; i++) {
			int clusterComponent = getProperComponent(atComponent, genoClusterCounts[i], i, centroid.getMarkerData().getMarkerName(), verbose, proj.getLog());
			if (validClusterComponent(genoClusterCounts[i], clusterComponent, numTotalSamples)) {

				CrossValidation cvalX = cvals[i][0];
				CrossValidation cvalY = cvals[i][1];

				if (!cvalX.analysisFailed()) {
					fullXPredicteds[i] = cvalX.getPredicteds();
					if (residStandardDeviationFilter != 0) {
						residX[i] = cvalX.getResiduals();
						residstdevX[i] = Array.stdev(Array.subArray(residX[i], genoSampleClusters[i]), true);// compute standard deviation only from members of that cluster
					} else {
						residX[i] = null;
						residstdevX[i] = Double.NaN;
					}
				} else {
					fullXPredicteds[i] = null;
					residX[i] = null;
					residstdevX[i] = Double.NaN;
				}
				if (!cvalY.analysisFailed()) {
					fullYPredicteds[i] = cvalY.getPredicteds();
					if (residStandardDeviationFilter != 0) {
						residY[i] = cvalY.getResiduals();
						residstdevY[i] = Array.stdev(Array.subArray(residY[i], genoSampleClusters[i]), true);// compute standard deviation only from members of that cluster
					} else {
						residY[i] = null;
						residstdevY[i] = Double.NaN;
					}

				} else {
					fullYPredicteds[i] = null;
					residY[i] = null;
					residstdevY[i] = Double.NaN;
				}
			} else {
				fullXPredicteds[i] = null;
				fullYPredicteds[i] = null;
			}
		}

		for (int i = 0; i < samplesToUse.length; i++) {
			byte tmpGenotype = (byte) -1;// default to not using it next round
			if (samplesToUse[i]) {// not in the pc file
				double minDist = Double.MAX_VALUE;
				for (int j = 0; j < genoSampleClusters.length; j++) {
					if (!forceThisCluster[j]) {// this cluster has two few inds, we won't use it for new genotypes
						if (fullXPredicteds[j] != null && fullYPredicteds[j] != null && genoClusterCounts[j] > 1) {// must have predicteds, and at more than one individual
							if (residStandardDeviationFilter == 0 || (Math.abs(residY[j][i]) < residStandardDeviationFilter * residstdevY[j] && Math.abs(residX[j][i]) < residStandardDeviationFilter * residstdevX[j])) {
								double tmpDist = cartDistance(Xs[i], Ys[i], fullXPredicteds[j][i], fullYPredicteds[j][i]);
								if (tmpDist < minDist) {
									minDist = tmpDist;
									tmpGenotype = (byte) j;
								}
							} else {

							}
						}
					}
				}
			}
			estimatedGenotypes[i] = tmpGenotype;
		}
		for (int i = 0; i < genoSampleClusters.length; i++) {// for each genotype

			for (int j = 0; j < genoSampleClusters[i].length; j++) {// sample boolean map
				if (forceThisCluster[i]) {
					if (genoSampleClusters[i][j]) {// if a sample belongs to this cluster
						estimatedGenotypes[j] = centroid.getMarkerData().getAbGenotypes()[j]; // which will always be the original at this point
					}
				} else if (estimatedGenotypes[j] == -1) {
					estimatedGenotypes[j] = centroid.getMarkerData().getAbGenotypes()[j]; // didnt get a new cluster, set to original
				} else {
					// new genotype!
				}
			}

		}
		int[] counts = new int[4];
		for (int i = 0; i < estimatedGenotypes.length; i++) {
			counts[estimatedGenotypes[i] + 1]++;
		}
		// log.report("Counts" + Array.toStr(counts));

		if (Arrays.equals(centroid.getAlternateGenotypes(), estimatedGenotypes)) {
			correctionMethod = KILL_INT;// no change with the new stage, kill after this
		}
		centroid.setAlternateGenotypes(estimatedGenotypes);
		centroid.init();// necessary to clear previous counts
		centroid.computeCentroid();// recompute centroid with these estimated genotypes
	}

	private CrossValidation[][] threadIt(int atComponent, double[] Xs, double[] Ys) {
		CrossValidation[][] cvals = new CrossValidation[3][2];
		ExecutorService executor = Executors.newFixedThreadPool(Math.min(6, numThreads));// da pool of threads, maximum of 6 needed
		Hashtable<String, Future<CrossValidation>> tmpResults = new Hashtable<String, Future<CrossValidation>>();
		for (int i = 0; i < genoSampleClusters.length; i++) {
			int clusterComponent = getProperComponent(atComponent, genoClusterCounts[i], i, centroid.getMarkerData().getMarkerName(), verbose, proj.getLog());
			if (validClusterComponent(genoClusterCounts[i], clusterComponent, numTotalSamples)) {
				String keyX = i + "_" + clusterComponent + "_X";
				String keyY = i + "_" + clusterComponent + "_Y";
				WorkerRegression workerX = new WorkerRegression(this, Xs, genoSampleClusters[i], clusterComponent, svdRegression, "Genotype cluster: " + i + " X values", verbose, log);
				tmpResults.put(keyX, executor.submit(workerX));
				WorkerRegression workerY = new WorkerRegression(this, Ys, genoSampleClusters[i], clusterComponent, svdRegression, "Genotype cluster: " + i + " Y values", verbose, log);
				tmpResults.put(keyY, executor.submit(workerY));
			} else {
				cvals[i][0] = null;
				cvals[i][1] = null;
			}
		}

		for (int i = 0; i < genoSampleClusters.length; i++) {
			int clusterComponent = getProperComponent(atComponent, genoClusterCounts[i], i, centroid.getMarkerData().getMarkerName(), verbose, proj.getLog());
			if (validClusterComponent(genoClusterCounts[i], clusterComponent, numTotalSamples)) {
				String keyX = i + "_" + clusterComponent + "_X";
				String keyY = i + "_" + clusterComponent + "_Y";
				try {
					cvals[i][0] = tmpResults.get(keyX).get();
					cvals[i][1] = tmpResults.get(keyY).get();// get is only applied after the job has finished
				} catch (InterruptedException e) {
					proj.getLog().reportError("Error - could not correct " + centroid.getMarkerData().getMarkerName() + " regression has failed");
					proj.getLog().reportException(e);
					fail = true;
				} catch (ExecutionException e) {
					proj.getLog().reportError("Error - could not correct " + centroid.getMarkerData().getMarkerName() + " regression has failed");
					log.reportException(e);
					fail = true;
				}
			}
		}
		executor.shutdown();
		try {
			executor.awaitTermination(7, TimeUnit.DAYS);
		} catch (InterruptedException e) {
			proj.getLog().reportException(e);
		}
		return cvals;
	}

	private double cartDistance(double x1, double y1, double x2, double y2) {
		return Math.sqrt(Math.pow((x1 - x2), 2) + Math.pow((y1 - y2), 2));
	}

	/**
	 * Consolidates the results of the clustered regression back into single float[] for X and Y
	 * 
	 */
	private void mergeCorrectedClusters() {
		int[] genoIndices = new int[genoSampleClusters.length];
		for (int i = 0; i < numTotalSamples; i++) {
			boolean hasCorrection = false;
			for (int j = 0; j < genoSampleClusters.length; j++) {
				if (genoSampleClusters[j][i] && centroid.getCentroid()[j] != null) {
					// correctedXFull[i] = correctedXCluster[j][genoIndices[j]] + centroid.getMarkerData().getXs()[j];
					// correctedYFull[i] = correctedYCluster[j][genoIndices[j]] + centroid.getMarkerData().getYs()[j];
					correctedXFull[i] = Math.max(0, correctedXCluster[j][genoIndices[j]] + (float) toCartesianX(centroid.getCentroid()[j]));
					correctedYFull[i] = Math.max(0, correctedYCluster[j][genoIndices[j]] + (float) toCartesianY(centroid.getCentroid()[j]));
					genoIndices[j]++;
					hasCorrection = true;
				}
			}
			if (!hasCorrection) {
				if (verbose) {
					proj.getLog().report("Warning - sample " + samples[i] + " could not be corrected for marker " + centroid.getMarkerData().getMarkerName());
				}
				correctedXFull[i] = centroid.getMarkerData().getXs()[i];
				correctedYFull[i] = centroid.getMarkerData().getYs()[i];

			}
		}
	}

	/**
	 * Here we assign samples to each of the three genotypes
	 * <p>
	 * We assign missing genotypes according to {@link PrincipalComponentsIntensity#CORRECTION_METHODS}
	 * <p>
	 * Warning - only samples with principal components will be analyzed
	 * 
	 */
	private void assignGenotypeClusterStatus(byte[] abGenotypes) {
		extractThetaCenters();
		if (!fail) {
			if (centroid.getMarkerData().getAbGenotypes() != null) {

				float[] xs = centroid.getMarkerData().getXs();
				float[] ys = centroid.getMarkerData().getYs();
				this.genoClusterCounts = new int[3];
				this.genoSampleClusters = new boolean[3][abGenotypes.length];
				for (int i = 0; i < abGenotypes.length; i++) {
					if (!samplesToUse[i] || centroid.getSamplesToUse() != null && !centroid.getSamplesToUse()[i]) {// samplesToUse represents inds with a pc those without will be skipped, and the centroid represents those used to cluster, those who were not used to cluster will not be included in the regression model
						assignAllFalseAt(i, genoSampleClusters);
					} else if (Float.isNaN(xs[i]) || Float.isNaN(ys[i])) {
						assignAllFalseAt(i, genoSampleClusters);
					} else if (abGenotypes[i] >= 0) {
						assignTrueAt(i, abGenotypes[i], genoSampleClusters);
						genoClusterCounts[abGenotypes[i]]++;
					} else if (correctionMethod == CORRECTION_INTS[0]) {// assign missing genotype,
						byte nearestThetaGenotype = findNearestTheta(genotypeThetaCenters, Centroids.calcTheta(centroid.getMarkerData().getXs()[i], centroid.getMarkerData().getYs()[i]));
						if (nearestThetaGenotype >= 0) {
							assignTrueAt(i, nearestThetaGenotype, genoSampleClusters);
							genoClusterCounts[nearestThetaGenotype]++;
						} else {
							assignAllFalseAt(i, genoSampleClusters);
							proj.getLog().reportError("Warning - marker " + centroid.getMarkerData().getMarkerName() + " for sample " + samples[i] + " had a missing genotype and could not be assigned to a theta cluster, skipping");
						}
					} else {
						assignAllFalseAt(i, genoSampleClusters);// simply skip missing markers
					}
				}
			} else {
				proj.getLog().reportError("Error - marker " + centroid.getMarkerData().getMarkerName() + " did not have ab genotypes, cannot correct intensity from genotype clusters for " + centroid.getMarkerData().getMarkerName());
			}
		}
	}

	/**
	 * Just retrieves the theta centers from the newly computed centroid. These are currently only used to assign missing genotype values...in {@link PrincipalComponentsIntensity#findNearestTheta(float[], float)}
	 */
	private void extractThetaCenters() {
		float[][] cent = centroid.getCentroid();
		this.genotypeThetaCenters = new float[3];
		boolean atLeastoneCenter = false;
		for (int i = 0; i < cent.length; i++) {
			if (cent[i] != null) {
				genotypeThetaCenters[i] = cent[i][0];
				atLeastoneCenter = true;
			} else {
				genotypeThetaCenters[i] = Float.NaN;
			}
		}
		if (!atLeastoneCenter) {
			proj.getLog().reportError("Error - marker " + centroid.getMarkerData().getMarkerName() + " did not have any cluster centers, cannot correct intensity for " + centroid.getMarkerData().getMarkerName());
			fail = !atLeastoneCenter;
		}
	}

	private static boolean validClusterComponent(int genoClusterCounts, int clusterComponent, int numSamples) {
		return isClusterBigEnough(genoClusterCounts, numSamples) && clusterComponent > 0;
	}

	public static boolean isClusterBigEnough(int count, int numSamples) {
		return (double) count / numSamples > MIN_CLUSTER_PERCENT;
	}

	/**
	 * @param i
	 *            assigns false to all clusters at this index
	 */
	private static void assignAllFalseAt(int i, boolean[][] genoSampleClusters) {
		for (int j = 0; j < genoSampleClusters.length; j++) {
			genoSampleClusters[j][i] = false; // individual was not in PC file so we cannot adjust
		}
	}

	/**
	 * @param i
	 *            assings true to cluster at index geno
	 * @param geno
	 */
	private static void assignTrueAt(int i, byte geno, boolean[][] genoSampleClusters) {
		for (int j = 0; j < genoSampleClusters.length; j++) {
			if (j == geno) {
				genoSampleClusters[j][i] = true;
			} else {
				genoSampleClusters[j][i] = false;
			}
		}
	}

	/**
	 * We check to make sure this cluster has enough individuals to use for a given number of components. We try to use the max number we can
	 */
	private static int getProperComponent(int atComponent, int numGenoCounts, int geno, String markerName, boolean verbose, Logger log) {
		int clusterComponent = atComponent;
		if (numGenoCounts <= clusterComponent) {// number of inds must be greater than the number of PCs being regressed
			clusterComponent = Math.max(numGenoCounts - 1, 0);
			if (clusterComponent > 0) {
				if (verbose) {
					log.report("Warning - marker " + markerName + " only has " + numGenoCounts + (numGenoCounts == 1 ? " individual " : " individuals") + " for the " + ScatterPlot.GENOTYPE_OPTIONS[geno + 1] + " genotype, can only correct for " + clusterComponent + (clusterComponent + clusterComponent == 1 ? " component " : " components"));
				}
			}
		}
		return clusterComponent;
	}

	/**
	 * @param genotypeThetaCenters
	 *            a theta cluster center for each genotype
	 * @param theta
	 *            the theta value to assign
	 * @return -1 if the cluster could not be assigned, else the closest gentype clustering
	 */
	private static byte findNearestTheta(float[] genotypeThetaCenters, float theta) {
		int minGenoDistanceCluster = -1;
		float curMinDist = Float.MAX_VALUE;
		for (int i = 0; i < genotypeThetaCenters.length; i++) {
			if (!Float.isNaN(genotypeThetaCenters[i]) && Math.abs(genotypeThetaCenters[i] - theta) < curMinDist) {
				curMinDist = Math.abs(genotypeThetaCenters[i] - theta);
				minGenoDistanceCluster = i;
			}
		}
		return (byte) minGenoDistanceCluster;
	}

	private static boolean invalidMethod(int correctionMethod, Logger log) {
		if (correctionMethod >= 0 && correctionMethod < CORRECTION_METHODS.length) {
			// log.report("Info - using correction Method: \"" + CORRECTION_METHODS[correctionMethod] + "\"");
			return false;
		} else {
			return true;
		}
	}

	// 0=theta,1=r
	private static float toCartesianX(float[] centroid) {
		return (float) (centroid[1] / (1 + Math.sin(centroid[0] * Math.PI / 2) / Math.cos(centroid[0] * Math.PI / 2)));
	}

	private static float toCartesianY(float[] centroid) {
		return (float) (centroid[1] / (1 + Math.cos(centroid[0] * Math.PI / 2) / Math.sin(centroid[0] * Math.PI / 2)));
	}

	private static CentroidCompute prepareProperCentroid(MarkerData markerData, int[] sampleSex, boolean[] samplesToUseCluster, double missingnessThreshold, double confThreshold, ClusterFilterCollection clusterFilterCollection, boolean medianCenter, Logger log) {
		CentroidCompute centroid;
		if (markerData.getMarkerName().startsWith(AFFY_INTENSITY_ONLY_FLAG[0])) {
			// centroid = markerData.getCentroid(sampleSex, samplesToUseCluster, true, missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, log);
			centroid = markerData.getCentroid(sampleSex, samplesToUseCluster, false, missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, log);
			setFakeAB(markerData, centroid, clusterFilterCollection, .1f);
		} else {
			centroid = markerData.getCentroid(sampleSex, samplesToUseCluster, false, missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, log);
			if (markerData.getMarkerName().startsWith(ILLUMINA_INTENSITY_ONLY_FLAG[0])) {
				setFakeAB(markerData, centroid, clusterFilterCollection, 0);
			}
		}
		centroid.computeCentroid();
		return centroid;
	}

	/**
	 * We use this in the case of intensity only probesets...We assign all genotypes to be the same
	 * 
	 */
	private static void setFakeAB(MarkerData markerData, CentroidCompute centroid, ClusterFilterCollection clusterFilterCollection, float gcThreshold) {
		byte[] fakeAB = new byte[centroid.getClustGenotypes().length];
		byte[] clustAB = markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerData.getMarkerName(), gcThreshold);
		int[] counts = new int[4];
		for (int i = 0; i < clustAB.length; i++) {
			counts[clustAB[i] + 1]++;
		}
		if (counts[0] == clustAB.length) {
			byte tmpCluster = Array.mean(markerData.getXs()) >= Array.mean(markerData.getYs()) ? (byte) 0 : (byte) 1;
			Arrays.fill(fakeAB, (byte) tmpCluster);
			centroid.setAlternateGenotypes(fakeAB);

		} else {
			centroid.setAlternateGenotypes(clustAB);
		}
		// assign all to this cluster
	}

	private static class WorkerRegression implements Callable<CrossValidation> {
		private PrincipalComponentsResiduals principalComponentsResiduals;
		private double[] data;
		private boolean[] samplesTobuildModel;
		private int clusterComponent;
		private boolean svdRegression;
		private String title;
		boolean verbose;

		// private Logger log;

		public WorkerRegression(PrincipalComponentsResiduals principalComponentsResiduals, double[] data, boolean[] samplesTobuildModel, int clusterComponent, boolean svdRegression, String title, boolean verbose, Logger log) {
			super();
			this.principalComponentsResiduals = principalComponentsResiduals;
			this.data = data;
			this.samplesTobuildModel = samplesTobuildModel;
			this.clusterComponent = clusterComponent;
			this.svdRegression = svdRegression;
			this.title = title;
			this.verbose = verbose;
			// this.log = log;
		}

		@Override
		public CrossValidation call() {// acts like run
			return principalComponentsResiduals.getCorrectedDataAt(data, samplesTobuildModel, clusterComponent, svdRegression, title, verbose);

		}
	}

	// private CrossValidation getRegression(double[] data, boolean[] samplesTobuildModel, int clusterComponent,boolean sv) {
	// CrossValidation cvalX = getCorrectedDataAt(data, samplesTobuildModel, clusterComponent, svdRegression, "Genotype cluster: " + i + " X values");
	// return cvalX;
	// }
	public static void main(String[] args) {
		// int numArgs = args.length;
		// String filename = "PrincipalComponentsIntensity.dat";
		// String pcFile = "PCA_final_ExomeMarkers_3050_UnrelatedW_NEWPCLRRMedian_BOTH.PCs.extrapolated.txt";
		String pcFile = "PCA_final_ExomeMarkers_3050_UnrelatedW.PCs.txt";
		String output = "test_corrected";
		// boolean svdRegression = true;
		int numcomponents = 2;
		Project proj = new Project(null, false);
		test(proj, pcFile, numcomponents, true, null, null, 1, 0, null, true, true, proj.getProjectDir() + output);
	}

	// if(Array.booleanArraySum(this.samplesToUse)!=markerData.getXs())
	public static void test(Project proj, String pcFile, int numComponents, boolean recomputeLRR, int[] sampleSex, boolean[] samplesToUseCluster, double missingnessThreshold, double confThreshold, ClusterFilterCollection clusterFilterCollection, boolean medianCenter, boolean svdRegression, String output) {

		// PrincipalComponentsResiduals testResids = new PrincipalComponentsResiduals(proj, pcFile, null, numComponents, false, (float) confThreshold, false, false, null);
		// String[] markers = HashVec.loadFileToStringArray(proj.getFilename(Project.DISPLAY_MARKERS_FILENAME), false, new int[] { 0 }, false);
		// MarkerDataLoader markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markers);
		// for (int i = 0; i < markers.length; i++) {
		// MarkerData markerData = markerDataLoader.requestMarkerData(i);
		// PrincipalComponentsIntensity pcIntensity = new PrincipalComponentsIntensity(testResids, markerData, recomputeLRR, sampleSex, samplesToUseCluster, missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, svdRegression, 1, 5, 2, true, output);
		// pcIntensity.correctXYRegression(100);
		// String corrected = "OrigX\tOrigY\tCorrectedX\tCorrectedY\n";
		// for (int j = 0; j < markerData.getXs().length; j++) {
		// corrected += markerData.getXs()[j] + "\t" + markerData.getYs()[j] + "\t" + pcIntensity.getCorrectedXFull()[j] + "\t" + pcIntensity.getCorrectedYFull()[j] + "\n";
		// }
		//
		// Files.write(corrected, proj.getProjectDir() + markers[i] + "_CompareCorrectedXY.txt");
		// if (i > 5) {
		// System.exit(1);
		// System.out.println((i + 1) + "/" + markers.length);
		//
		// }
		// markerDataLoader.releaseIndex(i);
		// }
	}
}

//
// *
// * Estimate the genotype based on the nearest cluster to the predicted values (across all regression models tested)
// */
// private void estimateNewGenotypes(int atComponent) {
// double[] xMatchedPCs = Array.toDoubleArray(Array.subArray(centroid.getMarkerData().getXs(), samplesToUse));
// double[] yMatchedPCs = Array.toDoubleArray(Array.subArray(centroid.getMarkerData().getYs(), samplesToUse));
// double[][] fullXPredicteds = new double[3][];
// double[][] fullYPredicteds = new double[3][];
// double[][] residX = new double[3][];
// double[][] residY = new double[3][];
// double[] residstdevX = new double[3];
// double[] residstdevY = new double[3];
// byte[] estimatedGenotypes = new byte[xMatchedPCs.length];
// double[] testXs = Array.toDoubleArray(centroid.getMarkerData().getXs());
// double[] testYs = Array.toDoubleArray(centroid.getMarkerData().getYs());
// for (int i = 0; i < genoSampleClusters.length; i++) {
// int clusterComponent = getProperComponent(atComponent, genoClusterCounts[i], i, centroid.getMarkerData().getMarkerName(), verbose, proj.getLog());
// if (validClusterComponent(genoClusterCounts[i], clusterComponent, numTotalSamples)) {
// double[][] train_indeps = genoPCsClusters[i].returnUpToNumComponents(clusterComponent);// train on pcs within a cluster
// long time =System.currentTimeMillis();
// CrossValidation cvalXTest = getCorrectedDataAt(testXs, genoSampleClusters[i], clusterComponent, svdRegression, "Genotype cluster: " +i+" X values");
// CrossValidation cvalYTest = getCorrectedDataAt(testYs, genoSampleClusters[i], clusterComponent, svdRegression, "Genotype cluster: " +i+" X values");
// System.out.println(ext.getTimeElapsed(time) + " for New method");
// time=System.currentTimeMillis();
// CrossValidation cvalX = getPredictedValues(genoXClusters[i], train_indeps, xMatchedPCs, getTrimmedPreppedPCs(clusterComponent,false), svdRegression, proj.getLog());
// CrossValidation cvalY = getPredictedValues(genoYClusters[i], train_indeps, yMatchedPCs, getTrimmedPreppedPCs(clusterComponent,false), svdRegression, proj.getLog());
// System.out.println(ext.getTimeElapsed(time) + " for OLD method");
//
// System.out.println("Are Xs equal" + Arrays.equals(cvalX.getResiduals(), cvalXTest.getResiduals()));
// System.out.println("Are Ys equal" + Arrays.equals(cvalY.getResiduals(), cvalYTest.getResiduals()));
//
//
// if (!cvalX.analysisFailed()) {
// fullXPredicteds[i] = cvalX.getPredicteds();
// if (residStandardDeviationFilter != 0) {
// residX[i] = cvalX.getResiduals();
// residstdevX[i] = Array.stdev(Array.subArray(residX[i], genoSampleClusters[i]), true);// compute standard deviation only from members of that cluster
// } else {
// residX[i] = null;
// residstdevX[i] = Double.NaN;
// }
// } else {
// fullXPredicteds[i] = null;
// residX[i] = null;
// residstdevX[i] = Double.NaN;
// }
// if (!cvalY.analysisFailed()) {
// fullYPredicteds[i] = cvalY.getPredicteds();
// if (residStandardDeviationFilter != 0) {
// residY[i] = cvalY.getResiduals();
// residstdevY[i] = Array.stdev(Array.subArray(residY[i], genoSampleClusters[i]), true);// compute standard deviation only from members of that cluster
// } else {
// residY[i] = null;
// residstdevY[i] = Double.NaN;
// }
//
// } else {
// fullYPredicteds[i] = null;
// residY[i] = null;
// residstdevY[i] = Double.NaN;
// }
// } else {
// fullXPredicteds[i] = null;
// fullYPredicteds[i] = null;
// }
// }
//
// for (int i = 0; i < samplesToUse.length; i++) {
// byte tmpGenotype = (byte) -1;// default to not using it next round
// if (samplesToUse[i]) {// not in the pc file
// double minDist = Double.MAX_VALUE;
// for (int j = 0; j < genoSampleClusters.length; j++) {
// if (!forceThisCluster[j]) {// this cluster has two few inds, we won't use it for new genotypes
// if (fullXPredicteds[j] != null && fullYPredicteds[j] != null && genoPCsClusters[j].getNumStored() > 1) {// must have predicteds, and at more than one individual
// if (residStandardDeviationFilter == 0 || (Math.abs(residY[j][i]) < residStandardDeviationFilter * residstdevY[j] && Math.abs(residX[j][i]) < residStandardDeviationFilter * residstdevX[j])) {
// double tmpDist = cartDistance(xMatchedPCs[i], yMatchedPCs[i], fullXPredicteds[j][i], fullYPredicteds[j][i]);
// if (tmpDist < minDist) {
// minDist = tmpDist;
// tmpGenotype = (byte) j;
// }
// } else {
//
// }
// }
// }
// }
// }
// estimatedGenotypes[i] = tmpGenotype;
// }
// for (int i = 0; i < genoSampleClusters.length; i++) {// for each genotype
//
// for (int j = 0; j < genoSampleClusters[i].length; j++) {// sample boolean map
// if (forceThisCluster[i]) {
// if (genoSampleClusters[i][j]) {// if a sample belongs to this cluster
// estimatedGenotypes[j] = centroid.getMarkerData().getAbGenotypes()[j]; // which will always be the original at this point
// }
// } else if (estimatedGenotypes[j] == -1) {
// estimatedGenotypes[j] = centroid.getMarkerData().getAbGenotypes()[j]; // didnt get a new cluster, set to original
// } else {
// // new genotype!
// }
// }
//
// }
// int[] counts = new int[4];
// for (int i = 0; i < estimatedGenotypes.length; i++) {
// counts[estimatedGenotypes[i] + 1]++;
// }
// log.report("Counts" + Array.toStr(counts));
//
// if (Arrays.equals(centroid.getAlternateGenotypes(), estimatedGenotypes)) {
// correctionMethod = KILL_INT;// no change with the new stage, kill after this
// }
// centroid.setAlternateGenotypes(estimatedGenotypes);
// centroid.init();// necessary to clear previous counts
// centroid.computeCentroid();// recompute centroid with these estimated genotypes
// }
//
//
// private static CrossValidation getPredictedValues(double[] train_deps, double[][] train_indeps, double[] val_deps, double[][] val_indeps, boolean svdRegression, Logger log) {
// CrossValidation cval = new CrossValidation(train_deps, train_indeps, val_deps, val_indeps, true, svdRegression, log);
// cval.train();
// cval.computePredictedValues();
// cval.computeResiduals();
// return cval;
// }
//
//
// private static class ClusterPC {
// private double[][] clusteredPCs;
// private int currentIndex;
//
// public ClusterPC(int numSamples) {
// this.clusteredPCs = new double[numSamples][];
// this.currentIndex = 0;
// }
//
// public void addToNextIndex(double[] data) {
// clusteredPCs[currentIndex] = data;
// currentIndex++;
// }
//
// public boolean allFull() {
// return currentIndex == clusteredPCs.length;
// }
//
// public int getNumStored() {
// return currentIndex;
// }
//
// private double[][] returnUpToNumComponents(int numComponent) {
// double[][] trimmed = new double[clusteredPCs.length][];
// if (numComponent == clusteredPCs[0].length) {
// trimmed = clusteredPCs;
// } else {
// for (int i = 0; i < clusteredPCs.length; i++) {
// trimmed[i] = Array.subArray(clusteredPCs[i], 0, numComponent);
// }
// }
// return trimmed;
// }
// }
//
// private void assignPCClusters() {
// if (!fail) {
// double[][] preppedPCs = this.getPreppedPCs();// sample dominant
// for (int i = 0; i < genoSampleClusters.length; i++) {// initialize underlying arrays
// genoPCsClusters[i] = new ClusterPC(Array.booleanArraySum(genoSampleClusters[i]));
// }
// for (int i = 0; i < genoSampleClusters.length; i++) {// for each genotype
// for (int j = 0; j < genoSampleClusters[i].length; j++) {// for each sample
// if (genoSampleClusters[i][j]) {// if the sample belongs to a genotype
// genoPCsClusters[i].addToNextIndex(preppedPCs[samplesInPc.get(samples[j])]);// add all the components, matched by the index in the pc file
// }
// }
// }
// for (int i = 0; i < genoPCsClusters.length; i++) {
// if (!genoPCsClusters[i].allFull()) {
// proj.getLog().reportError("Error - could not match all individuals with PCs to their respective genotypes");
// }
// }
// }
// }