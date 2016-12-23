package org.genvisis.cnv.analysis.pca;

import java.util.Arrays;
import java.util.HashSet;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.genvisis.cnv.analysis.CentroidCompute;
import org.genvisis.cnv.filesys.Centroids;
import org.genvisis.cnv.filesys.ClusterFilterCollection;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.manage.MDL;
import org.genvisis.cnv.manage.SexOps;
import org.genvisis.cnv.manage.SexOps.SEX_LOAD_TYPE;
import org.genvisis.cnv.plots.ScatterPlot;
import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.stats.CrossValidation;
import org.genvisis.stats.LeastSquares.LS_TYPE;

/**
 * Class to facilitate correcting X/Y intensity data; TODO, intensity only correction, should be
 * easy -> just like a single genotype
 *
 */
public class PrincipalComponentsIntensity extends PrincipalComponentsResiduals {
	private static final long serialVersionUID = 1L;

	private static final String[] CORRECTION_METHODS = {"Assign Missing Genotypes to Clusters (by nearest theta)",
																											"Two stage correction (use existing genotypes to find most predictive model)",
																											"N stage correction (use existing genotypes to find most predictive model)",
																											"N stage correction with residual filtering (use existing genotypes to find most predictive model)"};
	public static final String XY_RETURN = "X_Y";
	public static final String THETA_R_RETURN = "THETA_R";
	public static final String BAF_LRR_RETURN = "BAF_LRR";
	public static final double DEFAULT_CORRECTION_RATIO = 0.1;
	public static final double DEFAULT_RESID_STDV_FILTER = 0.0;

	private static final int[] CORRECTION_INTS = {0, 1, 2, 3};
	private static final int KILL_INT = -99;
	// private static final String[] AFFY_INTENSITY_ONLY_FLAG = { "CN_" };
	// private static final String[] ILLUMINA_INTENSITY_ONLY_FLAG = { "cnvi" };
	private static final double MIN_CLUSTER_PERCENT = 0.001;//
	private static final int MIN_CLUSTER_COUNT = 15;//

	private final CentroidCompute centroid;
	private boolean fail;

	private final boolean verbose;
	private final LS_TYPE lType;
	private int correctionMethod, nStage;

	private final int numThreads;
	private final String[] samples;
	private double[][] sexSpecificChrCovariates;// typically sex, and used
	private float[] genotypeThetaCenters, correctedXFull, correctedYFull;
	private float[][] correctedXCluster, correctedYCluster;
	private int[] genoClusterCounts;
	private final int numTotalSamples;
	private double residStandardDeviationFilter;

	private final double correctionRatio;
	private boolean[][] genoSampleClusters;// genotype, sample
	private final boolean[] forceThisCluster;
	private float[] correctedLRR;
	private final CHROMOSOME_X_STRATEGY sexStrategy;

	/**
	 * The two types of intensity correction we support
	 *
	 */
	public enum CORRECTION_TYPE {
																/**
																 * Multi-stage correction of X and Y intensities within a genotype
																 * cluster
																 */
																XY,
																/**
																 * A simple LRR only correction
																 */
																LRR_ONLY;
	}

	public enum CHROMOSOME_X_STRATEGY {
																				/**
																				 * Females should have ~2x intensity on chrX, etc
																				 */
																				BIOLOGICAL("Reflects reality, ~2x signal for females (vs. males) on chrX"),
																				/**
																				 * sex is explicitly regressed out along with PCs
																				 */
																				ARTIFICIAL("Sex is regressed out of intensity values");

		private String toolTip;

		private CHROMOSOME_X_STRATEGY(String toolTip) {
			this.toolTip = toolTip;

		}

		public String getToolTip() {
			return toolTip;
		}

	}

	public PrincipalComponentsIntensity(final PrincipalComponentsResiduals principalComponentsResiduals,
																			final MarkerData markerData, boolean recomputeLRR,
																			int[] sampleSex, boolean[] samplesToUseCluster,
																			double missingnessThreshold, double confThreshold,
																			ClusterFilterCollection clusterFilterCollection,
																			boolean medianCenter, LS_TYPE lType, int correctionMethod,
																			int nStage, double residStandardDeviationFilter,
																			double correctionRatio, int numThreads, boolean verbose,
																			String output, CHROMOSOME_X_STRATEGY sexStrategy) {
		super(principalComponentsResiduals);// we hijack the loading of the PC file and tracking of
																				// samples etc ...
		samples = proj.getSamples();
		this.lType = lType;
		this.correctionMethod = correctionMethod;
		this.residStandardDeviationFilter = residStandardDeviationFilter;
		this.verbose = verbose;
		this.nStage = nStage;
		fail = (markerData.getXs() == null|| markerData.getYs() == null
						|| markerData.getAbGenotypes() == null);
		centroid = prepareProperCentroid(	getProj().getArrayType(), markerData, sampleSex,
																			samplesToUseCluster, missingnessThreshold, confThreshold,
																			clusterFilterCollection, medianCenter, proj.getLog());
		numTotalSamples = markerData.getXs().length;
		correctedXFull = new float[numTotalSamples];
		correctedYFull = new float[numTotalSamples];
		forceThisCluster = new boolean[3];
		this.numThreads = numThreads;
		this.correctionRatio = correctionRatio;
		this.sexStrategy = sexStrategy;
		this.sexSpecificChrCovariates = setUpSexSpecificStrategy(sexStrategy);
		if (centroid.failed()) {
			fail = true;
		}
		if (invalidMethod(correctionMethod, getProj().getLog())) {
			fail = true;
		}
	}

	private double[][] setUpSexSpecificStrategy(CHROMOSOME_X_STRATEGY sexStrategy) {
		double[][] sexCovariates =null;
		if (centroid.isSexSpecific()) {

			if (sexStrategy == CHROMOSOME_X_STRATEGY.ARTIFICIAL
					&& centroid.getMarkerData().getChr() == 23) {
			int[] sex = SexOps.getSampleSex(proj, SEX_LOAD_TYPE.NUM_X_SEX);

				sexCovariates = new double[sex.length][];
			for (int i = 0; i < sexCovariates.length; i++) {
				if (sex[i] > 0) {
					sexCovariates[i] = new double[] {sex[i]};
				} else {
					sexCovariates[i] = new double[] {Double.NaN};
				}
			}
			return sexCovariates;
			} else if (centroid.getMarkerData().getChr() == 24) {
				int[] sex = SexOps.getSampleSex(proj, SEX_LOAD_TYPE.MAPPED_SEX);
				for (int i = 0; i < sex.length; i++) {
					if (sex[i] != 1) {
						// samplesToUse[i] = false;
						centroid.getSamplesToUse()[i]=false;
					}
				}
			}
		}
		return sexCovariates;
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

	/**
	 * TODO, do we want to force, and keep the correction ratio as well? This method is now a place
	 * holder if we want to revert to up-front forcing,
	 * {@link PrincipalComponentsIntensity#MIN_CLUSTER_PERCENT} must be altered if you ever want to
	 * force clusters.
	 */
	public boolean shouldforceOriginalGenotypes() {// to switch genotypes, we demand that we have at
																									// least two valid clusters above min size
		int goodClust = 0;
		for (int i = 0; i < genoClusterCounts.length; i++) {
			if (isClusterBigEnough(genoClusterCounts[i], numTotalSamples)) {
				goodClust++;
				// System.out.println("good Enough geno" + i);
			} else {
				if (verbose) {
					proj.getLog().report("Info - forcing genotype cluster " + i);
				}
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

	public void correctLRRAt(int atComponent) {
		setOriginal();
		this.correctedLRR = Array.toFloatArray(getCorrectedDataAt(
																															Array.toDoubleArray(centroid.getMarkerData()
																																													.getLRRs()),
																															sexSpecificChrCovariates,
																															centroid.getSamplesToUse(),
																															atComponent, lType,
																															"LRR correction at PC " + atComponent,
																															verbose).getResiduals());
	}



	public float[] getCorrectedLRR() {
		return correctedLRR;
	}

	public void correctXYAt(int atComponent) {
		if (!fail) {
			if (atComponent > getNumComponents()) {
				if (verbose) {
					log.report("Warning - cannot correct at "+ atComponent + ", only " + getNumComponents()
											+ " are available. Setting correction to " + getNumComponents());
				}
				atComponent = getNumComponents();
			}
			assignGenotypeClusterStatus(centroid.getClustGenotypes());
			if (hasSomeGenotypes()) {
				correctedXCluster = new float[genoSampleClusters.length][];
				correctedYCluster = new float[genoSampleClusters.length][];
				if (correctionMethod == CORRECTION_INTS[2] && shouldforceOriginalGenotypes()) {
					correctionMethod = KILL_INT;// we kill it right here
				}
				if ((correctionMethod == CORRECTION_INTS[0] && !fail) || correctionMethod == KILL_INT) {
					correctXYRegression(atComponent, sexSpecificChrCovariates);
				} else if (correctionMethod == CORRECTION_INTS[1] && !fail) {
					correctionMethod = KILL_INT;
					estimateNewGenotypes(atComponent, sexSpecificChrCovariates);// applies new genotypes to
																																			// clusters
					correctXYAt(atComponent);// estimate with new genotypes
				} else if (correctionMethod == CORRECTION_INTS[2] && !fail) {
					if (verbose) {
						log.report("Correcting at stage " + nStage);
					}
					nStage--;

					if (nStage >= 1) {// we have more stages to go
						estimateNewGenotypes(atComponent, sexSpecificChrCovariates);
					}
					if (nStage == 1 && residStandardDeviationFilter != 0) {// get final genotypes (there will
																																	// be no more missing genotypes)
						residStandardDeviationFilter = 0;
						estimateNewGenotypes(atComponent, sexSpecificChrCovariates);
					}
					if (nStage <= 1) {// final regression, can also be killed if gentypes are identical
														// between iterations
						correctionMethod = KILL_INT;
					}
					// applies new genotypes to clusters
					correctXYAt(atComponent);// estimate with new genotypes

				}
				mergeCorrectedClusters();
			} else {
				fail = true;
				// proj.getLog()
				// .reportError("Error - correction has failed for marker due to not having any genotypes "
				// + centroid.getMarkerData().getMarkerName()
				// + ", returning original x, y, and genotype values");
				setOriginal();
			}
		} else {
			proj.getLog()
					.reportError("Error - correction has failed during the centroid computation for marker "
												+ centroid.getMarkerData().getMarkerName()
												+ ", returning original x, y, and genotype values");
			setOriginal();
		}
	}

	public boolean isFail() {
		return fail;
	}

	private boolean hasSomeGenotypes() {
		for (int genoClusterCount : genoClusterCounts) {
			if (genoClusterCount > 0) {
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
			return new float[][] {centroid.getMarkerData().getXs(), centroid.getMarkerData().getYs()};
		} else {
			return new float[][] {getCorrectedXFull(), getCorrectedYFull()};
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
		byte[] genotypesToUse = originalGenotypes	? centroid.getMarkerData().getAbGenotypes()
																							: getGenotypesUsed();
		MarkerData tmpMarkerData = new MarkerData(centroid.getMarkerData().getMarkerName(),
																							centroid.getMarkerData().getChr(),
																							centroid.getMarkerData().getPosition(),
																							centroid.getMarkerData().getFingerprint(),
																							centroid.getMarkerData().getGCs(), null, null,
																							correctedXFull, correctedYFull, null, null, null,
																							null, genotypesToUse, genotypesToUse);

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
			return new float[][] {centroid.getMarkerData().getBAFs(), centroid.getMarkerData().getLRRs()};
		} else {

			centroid.setMarkerData(tmpMarkerData);
			centroid.computeCentroid(sexStrategy == CHROMOSOME_X_STRATEGY.ARTIFICIAL);
			float[] lrrs = centroid.getRecomputedLRR();
			if (isAffyIntensityOnly(getProj().getArrayType(), tmpMarkerData)) {
				centroid.setIntensityOnly(true); // this is to get proper LRRs /BAFs for affy CN_ probes.
																					// The correction process does not need this flag however
			}
			float[] bafs = centroid.getRecomputedBAF();
			return new float[][] {bafs, lrrs};
		}
	}

	private float[][] getCorrectedThetaR(boolean original, MarkerData tmpMarkerData) {
		if (original) {
			return new float[][] {centroid.getMarkerData().getThetas(), centroid.getMarkerData().getRs()};
		} else {
			return new float[][] {tmpMarkerData.getThetas(), tmpMarkerData.getRs()};
		}
	}

	private void correctXYRegression(int atComponent, double[][] extraIndeps) {
		CrossValidation[][] cvals = threadIt(	atComponent,
																					Array.toDoubleArray(centroid.getMarkerData().getXs()),
																					Array.toDoubleArray(centroid.getMarkerData().getYs()),
																					extraIndeps);
		for (int i = 0; i < genoSampleClusters.length; i++) {
			int clusterComponent = getProperComponent(atComponent, genoClusterCounts[i], i,
																								correctionRatio,
																								centroid.getMarkerData().getMarkerName(), verbose,
																								proj.getLog());
			if (cvals[i] != null&& cvals[i][1] != null && cvals[i][0] != null
					&& !cvals[i][0].analysisFailed() && !cvals[i][1].analysisFailed()
					&& validClusterComponent(genoClusterCounts[i], clusterComponent, numTotalSamples)) {
				correctedXCluster[i] = Array.toFloatArray(Array.subArray(	cvals[i][0].getResiduals(),
																																	genoSampleClusters[i]));
				correctedYCluster[i] = Array.toFloatArray(Array.subArray(	cvals[i][1].getResiduals(),
																																	genoSampleClusters[i]));
			} else {
				if (cvals[i] == null|| cvals[i][1] == null || cvals[i][0] == null
						|| cvals[i][0].analysisFailed() || cvals[i][1].analysisFailed()) {
					if (verbose) {
						proj.getLog().reportError("Analysis failed for "+ "Genotype cluster: " + i + " Marker "
																			+ centroid.getMarkerData().getMarkerName());
					}
				}
				genoSampleClusters[i] = new boolean[genoSampleClusters[i].length];
				Arrays.fill(genoSampleClusters[i], false); // not enough individuals, no corrections
			}
		}
	}

	/**
	 *
	 * Estimate the genotype based on the nearest cluster to the predicted values (across all
	 * regression models tested)
	 */
	private void estimateNewGenotypes(int atComponent, double[][] extraIndeps) {
		double[] Xs = Array.toDoubleArray(centroid.getMarkerData().getXs());
		double[] Ys = Array.toDoubleArray(centroid.getMarkerData().getYs());
		byte[] estimatedGenotypes = new byte[Xs.length];
		double[][] fullXPredicteds = new double[3][];
		double[][] fullYPredicteds = new double[3][];
		double[][] residX = new double[3][];
		double[][] residY = new double[3][];
		double[] residstdevX = new double[3];
		double[] residstdevY = new double[3];
		CrossValidation[][] cvals = threadIt(atComponent, Xs, Ys, extraIndeps);// organaized as
																																						// genotype,x/y

		for (int i = 0; i < genoSampleClusters.length; i++) {
			int clusterComponent = getProperComponent(atComponent, genoClusterCounts[i], i,
																								correctionRatio,
																								centroid.getMarkerData().getMarkerName(), verbose,
																								proj.getLog());
			if (validClusterComponent(genoClusterCounts[i], clusterComponent, numTotalSamples)) {

				CrossValidation cvalX = cvals[i][0];
				CrossValidation cvalY = cvals[i][1];

				if (cvalX != null && !cvalX.analysisFailed()) {
					fullXPredicteds[i] = cvalX.getPredicteds();
					if (residStandardDeviationFilter != 0) {
						residX[i] = cvalX.getResiduals();
						residstdevX[i] = Array.stdev(Array.subArray(residX[i], genoSampleClusters[i]), true);// compute
																																																	// standard
																																																	// deviation
																																																	// only
																																																	// from
																																																	// members
																																																	// of
																																																	// that
																																																	// cluster
					} else {
						residX[i] = null;
						residstdevX[i] = Double.NaN;
					}
				} else {
					fullXPredicteds[i] = null;
					residX[i] = null;
					residstdevX[i] = Double.NaN;
				}
				if (cvalY != null && !cvalY.analysisFailed()) {
					fullYPredicteds[i] = cvalY.getPredicteds();
					if (residStandardDeviationFilter != 0) {
						residY[i] = cvalY.getResiduals();
						residstdevY[i] = Array.stdev(Array.subArray(residY[i], genoSampleClusters[i]), true);// compute
																																																	// standard
																																																	// deviation
																																																	// only
																																																	// from
																																																	// members
																																																	// of
																																																	// that
																																																	// cluster
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
					if (!forceThisCluster[j]) {// this cluster has two few inds, we won't use it for new
																			// genotypes
						if (fullXPredicteds[j] != null&& fullYPredicteds[j] != null
								&& genoClusterCounts[j] > 1) {// must have predicteds, and at more than one
																							// individual
							if (!Double.isNaN(fullXPredicteds[j][i]) && !Double.isNaN(fullYPredicteds[j][i])) {
								if (Double.compare(residStandardDeviationFilter, 0) == 0
										|| (Math.abs(residY[j][i]) < residStandardDeviationFilter * residstdevY[j]
												&& Math.abs(residX[j][i]) < residStandardDeviationFilter
																										* residstdevX[j])) {
									double tmpDist = cartDistance(Xs[i], Ys[i], fullXPredicteds[j][i],
																								fullYPredicteds[j][i]);
									if (tmpDist < minDist) {
										minDist = tmpDist;
										tmpGenotype = (byte) j;
									}
								}
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
						estimatedGenotypes[j] = centroid.getMarkerData().getAbGenotypes()[j]; // which will
																																									// always be the
																																									// original at
																																									// this point
					}
				} else if (estimatedGenotypes[j] == -1) {
					estimatedGenotypes[j] = centroid.getMarkerData().getAbGenotypes()[j]; // didnt get a new
																																								// cluster, set to
																																								// original
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

	private CrossValidation[][] threadIt(	int atComponent, double[] Xs, double[] Ys,
																				final double[][] extraIndeps) {
		CrossValidation[][] cvals = new CrossValidation[3][2];
		ExecutorService executor = Executors.newFixedThreadPool(Math.min(6, numThreads));// da pool of
																																											// threads,
																																											// maximum of
																																											// 6
																																											// needed
		for (int i = 0; i < genoSampleClusters.length; i++) {
			int clusterComponent = getProperComponent(atComponent, genoClusterCounts[i], i,
																								correctionRatio,
																								centroid.getMarkerData().getMarkerName(), verbose,
																								proj.getLog());
			if (verbose) {
				proj.getLog().reportTimeInfo("Attempting to enter regression model with "+ clusterComponent
																			+ " principal components for genotype cluster " + i);
			}
			if (validClusterComponent(genoClusterCounts[i], clusterComponent, numTotalSamples)) {
				double[][] extraIndepsToUse = validateExtraIndeps(extraIndeps, genoSampleClusters[i]);
				// if(extraIndeps==null&&e)
				WorkerRegression workerX = new WorkerRegression(this, Xs, genoSampleClusters[i],
																												clusterComponent, lType,
																												"Genotype cluster: "			+ i + " X values for Marker "
																																									+ centroid.getMarkerData()
																																														.getMarkerName(),
																												cvals[i], extraIndepsToUse, 0, verbose,
																												log);

				executor.submit(workerX);
				WorkerRegression workerY = new WorkerRegression(this, Ys, genoSampleClusters[i],
																												clusterComponent, lType,
																												"Genotype cluster: "			+ i + " Y values for Marker "
																																									+ centroid.getMarkerData()
																																														.getMarkerName(),
																												cvals[i], extraIndepsToUse, 1, verbose,
																												log);
				executor.submit(workerY);
			} else {
				cvals[i][0] = null;
				cvals[i][1] = null;
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
		int count = 0;
		boolean hasAnyCorrection = false;
		for (int i = 0; i < numTotalSamples; i++) {
			boolean hasCorrection = false;
			for (int j = 0; j < genoSampleClusters.length; j++) {
				if (genoSampleClusters[j][i] && centroid.getCentroid()[j] != null) {
					// correctedXFull[i] = correctedXCluster[j][genoIndices[j]] +
					// centroid.getMarkerData().getXs()[j];
					// correctedYFull[i] = correctedYCluster[j][genoIndices[j]] +
					// centroid.getMarkerData().getYs()[j];
					correctedXFull[i] = Math.max(0, correctedXCluster[j][genoIndices[j]]
																					+ toCartesianX(centroid.getCentroid()[j]));
					correctedYFull[i] = Math.max(0, correctedYCluster[j][genoIndices[j]]
																					+ toCartesianY(centroid.getCentroid()[j]));
					genoIndices[j]++;
					hasCorrection = true;
					hasAnyCorrection = true;
					count++;
				}
			}
			if (!hasCorrection) {
				// if (verbose) {
				// proj.getLog().report("Warning - sample " + samples[i] + " could not be corrected for
				// marker " + centroid.getMarkerData().getMarkerName());
				// }
				correctedXFull[i] = centroid.getMarkerData().getXs()[i];
				correctedYFull[i] = centroid.getMarkerData().getYs()[i];

			}
		}
		if (verbose) {
			proj.getLog().report("Info - corrected a total of "+ count + " of " + numTotalSamples + " "
														+ (count == 1 ? "sample" : "samples"));
		}
		if (!hasAnyCorrection) {
			proj.getLog().reportTimeWarning("Marker "+ centroid.getMarkerData().getMarkerName()
																			+ " was unable to be corrected...");
		}
	}

	/**
	 * Here we assign samples to each of the three genotypes
	 * <p>
	 * We assign missing genotypes according to
	 * {@link PrincipalComponentsIntensity#CORRECTION_METHODS}
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
				genoClusterCounts = new int[3];
				genoSampleClusters = new boolean[3][abGenotypes.length];
				for (int i = 0; i < abGenotypes.length; i++) {
					if (!samplesToUse[i]
							|| centroid.getSamplesToUse() != null && !centroid.getSamplesToUse()[i]) {// samplesToUse
																																												// represents
																																												// inds with
																																												// a pc
																																												// those
																																												// without
																																												// will be
																																												// skipped,
																																												// and the
																																												// centroid
																																												// represents
																																												// those
																																												// used to
																																												// cluster,
																																												// those who
																																												// were not
																																												// used to
																																												// cluster
																																												// will not
																																												// be
																																												// included
																																												// in the
																																												// regression
																																												// model
						assignAllFalseAt(i, genoSampleClusters);
					} else if (Float.isNaN(xs[i]) || Float.isNaN(ys[i])) {
						assignAllFalseAt(i, genoSampleClusters);
					} else if (abGenotypes[i] >= 0) {
						assignTrueAt(i, abGenotypes[i], genoSampleClusters);
						genoClusterCounts[abGenotypes[i]]++;
					} else if (correctionMethod == CORRECTION_INTS[0]) {// assign missing genotype,
						byte nearestThetaGenotype = findNearestTheta(	genotypeThetaCenters,
																													Centroids.calcTheta(centroid.getMarkerData()
																																											.getXs()[i],
																																							centroid.getMarkerData()
																																											.getYs()[i]));
						if (nearestThetaGenotype >= 0) {
							assignTrueAt(i, nearestThetaGenotype, genoSampleClusters);
							genoClusterCounts[nearestThetaGenotype]++;
						} else {
							assignAllFalseAt(i, genoSampleClusters);
							proj.getLog()
									.reportError("Warning - marker "+ centroid.getMarkerData().getMarkerName()
																+ " for sample " + samples[i]
																+ " had a missing genotype and could not be assigned to a theta cluster, skipping");
						}
					} else {
						assignAllFalseAt(i, genoSampleClusters);// simply skip missing markers
					}
				}
				if (verbose) {
					for (int i = 0; i < genoClusterCounts.length; i++) {
						proj.getLog().reportTimeInfo("Genotype Cluster: " + i + ", n=" + genoClusterCounts[i]);
					}
				}
			} else {
				proj.getLog()
						.reportError("Error - marker "+ centroid.getMarkerData().getMarkerName()
													+ " did not have ab genotypes, cannot correct intensity from genotype clusters for "
													+ centroid.getMarkerData().getMarkerName());
			}
		}
	}

	/**
	 * Just retrieves the theta centers from the newly computed centroid. These are currently only
	 * used to assign missing genotype values...in
	 * {@link PrincipalComponentsIntensity#findNearestTheta(float[], float)}
	 */
	private void extractThetaCenters() {
		float[][] cent = centroid.getCentroid();
		genotypeThetaCenters = new float[3];
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
			proj.getLog()
					.reportError("Error - marker "+ centroid.getMarkerData().getMarkerName()
												+ " did not have any cluster centers, cannot correct intensity for "
												+ centroid.getMarkerData().getMarkerName());
			fail = !atLeastoneCenter;
		}
	}

	private static boolean validClusterComponent(	int genoClusterCounts, int clusterComponent,
																								int numSamples) {
		return isClusterBigEnough(genoClusterCounts, numSamples) && clusterComponent > 0;
	}

	/**
	 * Make sure there is variance in the indeps for this cluster
	 */
	private static double[][] validateExtraIndeps(double[][] extraIndeps,
																								boolean[] genoSampleClusters) {
		if (extraIndeps != null) {
			for (int i = 0; i < extraIndeps[0].length; i++) {
				HashSet<String> variance = new HashSet<String>();

				for (int j = 0; j < genoSampleClusters.length; j++) {
					if (genoSampleClusters[j]) {
						if (!Double.isNaN(extraIndeps[j][i])) {
							variance.add(extraIndeps[j][i] + "");
							if (variance.size() > 1) {
								break;
							}
						}
					}
				}
				if (variance.size() < 2) {
					return null;
				}

			}
			return extraIndeps;
		}

		return null;

	}

	private static boolean isClusterBigEnough(int count, int numSamples) {
		return (double) count / numSamples > MIN_CLUSTER_PERCENT && count >= MIN_CLUSTER_COUNT;
	}

	/**
	 * @param i assigns false to all clusters at this index
	 */
	private static void assignAllFalseAt(int i, boolean[][] genoSampleClusters) {
		for (int j = 0; j < genoSampleClusters.length; j++) {
			genoSampleClusters[j][i] = false; // individual was not in PC file so we cannot adjust
		}
	}

	/**
	 * @param i assings true to cluster at index geno
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
	 * We check to make sure this cluster has enough individuals to use for a given number of
	 * components. We now use the max number we can if correction ratio is 0
	 */
	private static int getProperComponent(int atComponent, int numGenoCounts, int geno,
																				double correctionRatio, String markerName, boolean verbose,
																				Logger log) {
		int clusterComponent = atComponent;
		double currentCorrectionRatio = (double) clusterComponent / numGenoCounts;
		if (currentCorrectionRatio > correctionRatio) {// need adjust the ratio
			clusterComponent = Math.max(Math.round((float) correctionRatio * numGenoCounts - 1), 0);
			clusterComponent = Math.min(clusterComponent, numGenoCounts - 1);// number of inds must be
																																				// greater than the number
																																				// of
																																				// PCs being regressed

			if (clusterComponent > 0) {
				if (verbose) {
					log.report("Warning - marker "+ markerName + " only has " + numGenoCounts
											+ (numGenoCounts == 1 ? " individual " : " individuals") + " for the "
											+ ScatterPlot.GENOTYPE_OPTIONS[geno + 1] + " genotype, can only correct for "
											+ clusterComponent
											+ (clusterComponent + clusterComponent == 1 ? " component " : " components"));
				}
			}
		}
		if (correctionRatio == 1) {
			clusterComponent = Math.max(numGenoCounts - 1, 0);
		}
		clusterComponent = Math.min(atComponent, clusterComponent);
		if (verbose) {
			log.report("Clustering at component "+ clusterComponent + "\tCurrent correction Ratio:"
									+ currentCorrectionRatio + "\tMax correction Ratio:" + correctionRatio
									+ "\tNumber in cluster " + geno + ":" + numGenoCounts);
		}
		return clusterComponent;
	}

	/**
	 * @param genotypeThetaCenters a theta cluster center for each genotype
	 * @param theta the theta value to assign
	 * @return -1 if the cluster could not be assigned, else the closest gentype clustering
	 */
	private static byte findNearestTheta(float[] genotypeThetaCenters, float theta) {
		int minGenoDistanceCluster = -1;
		float curMinDist = Float.MAX_VALUE;
		for (int i = 0; i < genotypeThetaCenters.length; i++) {
			if (!Float.isNaN(genotypeThetaCenters[i])
					&& Math.abs(genotypeThetaCenters[i] - theta) < curMinDist) {
				curMinDist = Math.abs(genotypeThetaCenters[i] - theta);
				minGenoDistanceCluster = i;
			}
		}
		return (byte) minGenoDistanceCluster;
	}

	private static boolean invalidMethod(int correctionMethod, Logger log) {
		if (correctionMethod >= 0 && correctionMethod < CORRECTION_METHODS.length) {
			// log.report("Info - using correction Method: \"" + CORRECTION_METHODS[correctionMethod] +
			// "\"");
			return false;
		} else {
			return true;
		}
	}

	// 0=theta,1=r
	private static float toCartesianX(float[] centroid) {
		return (float) (centroid[1] / (1 + Math.sin(centroid[0] * Math.PI / 2)
																				/ Math.cos(centroid[0] * Math.PI / 2)));
	}

	private static float toCartesianY(float[] centroid) {
		return (float) (centroid[1] / (1 + Math.cos(centroid[0] * Math.PI / 2)
																				/ Math.sin(centroid[0] * Math.PI / 2)));
	}

	private static CentroidCompute prepareProperCentroid(	ARRAY array, MarkerData markerData,
																												int[] sampleSex,
																												boolean[] samplesToUseCluster,
																												double missingnessThreshold,
																												double confThreshold,
																												ClusterFilterCollection clusterFilterCollection,
																												boolean medianCenter, Logger log) {
		CentroidCompute centroid;
		if (isAffyIntensityOnly(array, markerData)) {
			// centroid = markerData.getCentroid(sampleSex, samplesToUseCluster, true,
			// missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, log);
			centroid = markerData.getCentroid(sampleSex, samplesToUseCluster, false, missingnessThreshold,
																				confThreshold, clusterFilterCollection, medianCenter, log);
			CentroidCompute.setFakeAB(markerData, centroid, clusterFilterCollection, .1f, log);
		} else {
			centroid = markerData.getCentroid(sampleSex, samplesToUseCluster, false, missingnessThreshold,
																				confThreshold, clusterFilterCollection, medianCenter, log);
			if (array.isCNOnly(markerData.getMarkerName())) {
				CentroidCompute.setFakeAB(markerData, centroid, clusterFilterCollection, 0, log);
			}
		}
		centroid.computeCentroid();
		return centroid;
	}

	private static boolean isAffyIntensityOnly(ARRAY array, MarkerData markerData) {
		return (array == ARRAY.AFFY_GW6 || array == ARRAY.AFFY_GW6_CN)
						&& array.isCNOnly(markerData.getMarkerName());
	}


	private class WorkerRegression implements Runnable {
		private final PrincipalComponentsResiduals principalComponentsResiduals;
		private final double[] data;
		private final boolean[] samplesTobuildModel;
		private final double[][] extraIndeps;
		private final int clusterComponent;
		private final LS_TYPE lType;
		private final String title;
		private final CrossValidation[] out;
		private final int outIndex;
		private final Logger log;
		boolean verbose;

		public WorkerRegression(PrincipalComponentsResiduals principalComponentsResiduals,
														double[] data, boolean[] samplesTobuildModel, int clusterComponent,
														LS_TYPE lType, String title, CrossValidation[] output,
														double[][] extraIndeps, int outIndex, boolean verbose, Logger log) {
			super();
			this.principalComponentsResiduals = principalComponentsResiduals;
			this.data = data;
			this.samplesTobuildModel = samplesTobuildModel;
			this.clusterComponent = clusterComponent;
			this.lType = lType;
			this.title = title;
			this.verbose = verbose;
			this.out = output;
			this.outIndex = outIndex;
			this.extraIndeps = extraIndeps;
			this.log = log;
		}

		@Override
		public void run() {
			try {
				out[outIndex] = principalComponentsResiduals.getCorrectedDataAt(data, extraIndeps,
																																				samplesTobuildModel,
																																				clusterComponent, lType,
																																				title, verbose);
			} catch (Exception e) {
				log.reportError("Error - could not correct " + title + " regression has failed");
				log.reportException(e);
				fail = true;
			}
		}
	}

	/**
	 * corrects on multiple threads
	 *
	 *
	 */
	public static class PcCorrectionProducer extends AbstractProducer<PrincipalComponentsIntensity> {

		private final PrincipalComponentsResiduals pcResiduals;
		private final MDL mdl;
		private final int[] sampleSex;
		private final boolean[] samplesToUseCluster;
		private final LS_TYPE lType;
		private final int numCorrectionThreads;
		// private final int numDecompressThreads;
		// private final String[] markersToCorrect;
		private final int correctAt;
		private final CORRECTION_TYPE correctionType;
		private final CHROMOSOME_X_STRATEGY sexStrategy;

		public PcCorrectionProducer(PrincipalComponentsResiduals pcResiduals, int correctAt,
																int[] sampleSex, boolean[] samplesToUseCluster, LS_TYPE lType,
																int numCorrectionThreads, int numDecompressThreads,
																String[] markersToCorrect, CORRECTION_TYPE correctionType,
																CHROMOSOME_X_STRATEGY sexStrategy) {
			super();
			this.pcResiduals = pcResiduals;
			this.sampleSex = sampleSex;
			this.samplesToUseCluster = samplesToUseCluster;
			this.lType = lType;
			this.numCorrectionThreads = numCorrectionThreads;
			this.correctionType = correctionType;
			this.sexStrategy = sexStrategy;
			// this.numDecompressThreads = numDecompressThreads;
			// this.markersToCorrect = markersToCorrect;
			this.correctAt = correctAt;
			mdl = new MDL(pcResiduals.getProj(), pcResiduals.getProj().getMarkerSet(), markersToCorrect,
										numDecompressThreads, 0);
		}

		@Override
		public boolean hasNext() {
			return mdl.hasNext();
		}

		@Override
		public Callable<PrincipalComponentsIntensity> next() {
			final MarkerData markerData = mdl.next();
			return new Callable<PrincipalComponentsIntensity>() {
				@Override
				public PrincipalComponentsIntensity call() throws Exception {
					PrincipalComponentsIntensity principalComponentsIntensity =
																																		new PrincipalComponentsIntensity(	pcResiduals,
																																																			markerData,
																																																			true,
																																																			sampleSex,
																																																			samplesToUseCluster,
																																																			1,
																																																			0,
																																																			null,
																																																			true,
																																																			lType,
																																																			2,
																																																			5,
																																																			DEFAULT_RESID_STDV_FILTER,
																																																			DEFAULT_CORRECTION_RATIO,
																																																			numCorrectionThreads,
																																																			false,
																																																			null,
																																																			sexStrategy);
					switch (correctionType) {
						case LRR_ONLY:
							principalComponentsIntensity.correctLRRAt(correctAt);
							break;
						case XY:
							principalComponentsIntensity.correctXYAt(correctAt);
							break;
						default:
							break;

					}
					return principalComponentsIntensity;
				}

			};
		}

		@Override
		public void shutdown() {
			mdl.shutdown();
		}
	}

	public static void main(String[] args) {
		// int numArgs = args.length;
		// String filename = "PrincipalComponentsIntensity.dat";
		// String pcFile =
		// "PCA_final_ExomeMarkers_3050_UnrelatedW_NEWPCLRRMedian_BOTH.PCs.extrapolated.txt";
		String pcFile = "PCA_final_ExomeMarkers_3050_UnrelatedW.PCs.txt";
		String output = "test_corrected";
		// boolean svdRegression = true;
		int numcomponents = 2;
		Project proj = new Project(null, false);
		test(	proj, pcFile, numcomponents, true, null, null, 1, 0, null, true, true,
					proj.PROJECT_DIRECTORY.getValue() + output);
	}

	public static void test(Project proj, String pcFile, int numComponents, boolean recomputeLRR,
													int[] sampleSex, boolean[] samplesToUseCluster,
													double missingnessThreshold, double confThreshold,
													ClusterFilterCollection clusterFilterCollection, boolean medianCenter,
													boolean svdRegression, String output) {


	}
}

