package cnv.qc;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;

import seq.manage.ReferenceGenome;
import stats.CrossValidation;
import cnv.filesys.MarkerSet;
import cnv.filesys.MarkerSet.PreparedMarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.var.LocusSet;
import common.Array;
import common.Files;
import common.Logger;
import common.Positions;
import common.ext;
import filesys.Segment;

/**
 * Class for correcting intensity using gc content, see http://nar.oxfordjournals.org/content/36/19/e126.full.pdf
 * <p>
 * We deviate slightly by determining median - binned gc content from a provided gc file, instead of using the PennCNV reference below...however the refernence can be flagged
 * <p>
 * It is important to note that regression indices are determined by distance between markers, and qc (WF,GCWF) are base-pair bins covering the genome
 * <p>
 * The current version takes ~ 0.5 seconds for a sample with ~600K markers. From PennCNV
 * <p>
 * #ref_median can be GC percentage in 135 1MB non-overlapping sliding windows in chr11 (2004 or 2006 human genome assembly)
 * 
 * <p>
 * my @ref_median = qw/54.8207535282258 56.8381472081218 53.1218950320513 46.9484174679487 39.9367227359694 38.3365384615385 41.9867788461538 40.4431401466837 44.5320512820513 42.1979166666667 41.6984215561224 43.1598557692308 43.4388020833333 40.8104967948718 39.8444475446429 41.5357572115385 38.7496995192308 45.0213249362245 42.3251201923077 43.5287459935897 40.7440808354592 37.0492788461538 36.5006009615385 35.8518016581633 35.2767427884615 35.1972155448718 36.5286192602041 39.4890825320513 36.5779246794872 36.7275641025641 38.3256935586735 37.791266025641 41.1777844551282 41.950534119898 42.3639823717949 41.9208733974359 41.2061543367347 35.4974959935897 35.2123397435897 36.5101841517857 36.7135416666667 36.8268229166667 37.6945153061224 40.7453926282051 47.7049278846154
 * 47.3233173076923 44.7361288265306 46.6585536858974 39.1593549679487 36.5684789540816 38.2718466806667 37.184425 37.184425 37.184425 37.184425 35.9227764423077 41.1157852564103 41.6662348533163 39.7402844551282 40.0149238782051 46.6417211415816 49.9136618589744 45.2016225961538 51.3019172512755 52.0818309294872 51.1320112179487 49.9807185102302 49.9807185102302 49.5874473187766 50.547349024718 50.7186498397436 45.6435347576531 46.3352363782051 42.4091546474359 46.6399274553571 43.7746394230769 45.0160256410256 41.8526642628205 43.8899075255102 38.5112179487179 36.1038661858974 36.1689851721939 39.8506610576923 37.0439703525641 36.8012595663265 40.2521033653846 39.661858974359 37.5013769093564 35.5448717948718 36.9039979272959 35.2046274038462 38.2195512820513 40.074537627551
 * 40.7097355769231 40.5470753205128 38.4104380072343 36.131109775641 35.3915264423077 34.9693080357143 36.2953725961538 37.9602363782051 39.1942362882653 37.4464142628205 36.8879206730769 35.7242588141026 36.7556202168367 37.0639022435897 40.6929086538462 38.385084502551 39.4121594551282 40.2410857371795 42.0772879464286 43.2935697115385 43.2345753205128 40.9113919005102 44.9575320512821 46.2513020833333 46.4753069196429 48.3886217948718 47.8520633012821 43.8001802884615 39.808274872449 44.5042067307692 38.3835136217949 44.9097177933673 45.5366586538462 41.7346754807692 39.2198461415816 41.9489182692308 44.3351362179487 42.7910754145408 42.3190104166667 42.0425681089744 47.0514787946429 45.3482603740699/;
 * 
 * 
 */
public class GcAdjustor {

	public static final String[] GC_ADJUSTOR_TITLE = { "Adjust by GC content" };
	private static final int DEFUALT_GC_MODEL_WINDOW_SNP = 100 * 5120;
	private static final int DEFUALT_GC_MODEL_WINDOW_GC = 5120;

	/**
	 * These defaults are also used by PennCNV - detect_cnv.pl, but might be worth adjusting sometime, somewhere, in a galaxy called msi
	 */
	public static final double DEFAULT_MINAUTOSOMALGC = 15d;
	public static final double DEFAULT_MAXAUTOSOMALGC = 80d;
	public static final double DEFAULT_MIN_DATA_VALUE = -1d;
	public static final double DEFAULT_MAX_DATA_VALUE = 1d;
	public static final int[] DEFAULT_REGRESSION_DISTANCE = { 1000000 };
	public static final int[] DEFAULT_SKIP_PER_CHR = { 0 };
	public static final int DEFUALT_NUM_SNP_MAD = 10;// not available as command line option currently
	public static final double[] DEFUALT_PENNCNV_CHR11_GC_BINS = { 54.8207535282258, 56.8381472081218, 53.1218950320513, 46.9484174679487, 39.9367227359694, 38.3365384615385, 41.9867788461538, 40.4431401466837, 44.5320512820513, 42.1979166666667, 41.6984215561224, 43.1598557692308, 43.4388020833333, 40.8104967948718, 39.8444475446429, 41.5357572115385, 38.7496995192308, 45.0213249362245, 42.3251201923077, 43.5287459935897, 40.7440808354592, 37.0492788461538, 36.5006009615385, 35.8518016581633, 35.2767427884615, 35.1972155448718, 36.5286192602041, 39.4890825320513, 36.5779246794872, 36.7275641025641, 38.3256935586735, 37.791266025641, 41.1777844551282, 41.950534119898, 42.3639823717949, 41.9208733974359, 41.2061543367347, 35.4974959935897, 35.2123397435897, 36.5101841517857, 36.7135416666667, 36.8268229166667, 37.6945153061224, 40.7453926282051, 47.7049278846154, 47.3233173076923, 44.7361288265306, 46.6585536858974, 39.1593549679487, 36.5684789540816, 38.2718466806667, 37.184425, 37.184425, 37.184425, 37.184425, 35.9227764423077, 41.1157852564103, 41.6662348533163, 39.7402844551282, 40.0149238782051, 46.6417211415816, 49.9136618589744, 45.2016225961538, 51.3019172512755, 52.0818309294872, 51.1320112179487, 49.9807185102302, 49.9807185102302, 49.5874473187766, 50.547349024718, 50.7186498397436, 45.6435347576531, 46.3352363782051, 42.4091546474359, 46.6399274553571, 43.7746394230769, 45.0160256410256, 41.8526642628205, 43.8899075255102, 38.5112179487179, 36.1038661858974, 36.1689851721939, 39.8506610576923, 37.0439703525641, 36.8012595663265, 40.2521033653846, 39.661858974359, 37.5013769093564, 35.5448717948718, 36.9039979272959, 35.2046274038462, 38.2195512820513, 40.074537627551, 40.7097355769231, 40.5470753205128, 38.4104380072343, 36.131109775641, 35.3915264423077, 34.9693080357143, 36.2953725961538, 37.9602363782051, 39.1942362882653, 37.4464142628205, 36.8879206730769, 35.7242588141026, 36.7556202168367, 37.0639022435897, 40.6929086538462, 38.385084502551, 39.4121594551282, 40.2410857371795, 42.0772879464286, 43.2935697115385, 43.2345753205128, 40.9113919005102, 44.9575320512821, 46.2513020833333, 46.4753069196429, 48.3886217948718, 47.8520633012821, 43.8001802884615, 39.808274872449, 44.5042067307692, 38.3835136217949, 44.9097177933673, 45.5366586538462, 41.7346754807692, 39.2198461415816, 41.9489182692308, 44.3351362179487, 42.7910754145408, 42.3190104166667, 42.0425681089744, 47.0514787946429, 45.3482603740699 };

	private Project proj;
	private GcModel gcModel;
	private CrossValidation crossValidation;
	private int[] correctedIndices;// the indices we were able to correct, as ordered by the input data
	private int[][] qcIndices;// the indices used to compute qc metrics, as ordered by the input-gc matched arrays
	private int[][] chr11qcIndices;// only used if defualy penncnv GCWF calculation is desired
	private double minimumAutosomalGC, maximimumAutosomalGC, minIntensity, maxIntensity, wfPrior, wfPost, gcwfPrior, gcwfPost;
	private double[] regressionGcs, fullGcs;
	private double[] fullIntensity, regressionIntensity, markerIntensities, correctedIntensities;
	private int regressionDistance, skipPerChr;
	private boolean[] markerMask;
	private boolean fail, verbose;
	private GC_CORRECTION_METHOD correctionMethod;
	private PreparedMarkerSet preparedMarkerSet;

	public enum GC_CORRECTION_METHOD {
		/**
		 * Perform gc correction using the penncnv defaults (chr 11, etc)
		 */
		PENNCNV_GC, /**
		 * Perform gc correction using genvisis defaults (all available markers)
		 */
		GENVISIS_GC;
	}

	// /**
	// * Constructor for PennCNV defaults, minus the chromosome 11 GC content bins for the GCWF computation
	// */
	// public GcAdjustor(Project proj, PreparedMarkerSet preparedMarkerSet, GcModel gcModel, double[] markerIntensities, boolean[] markerMask, GC_CORRECTION_METHOD correctionMethod, boolean verbose) {
	// this(proj, preparedMarkerSet, gcModel, markerIntensities, DEFAULT_MINAUTOSOMALGC, DEFAULT_MAXAUTOSOMALGC, DEFAULT_MIN_DATA_VALUE, DEFAULT_MAX_DATA_VALUE, DEFAULT_REGRESSION_DISTANCE[0], DEFAULT_SKIP_PER_CHR[0], markerMask, correctionMethod, verbose);
	// }
	//
	// /**
	// * @param proj
	// * current project
	// * @param gcModel
	// * valid {@link GcModel} object to use for correction
	// * @param markerIntensities
	// * array of intensities to correct
	// * @param minimumAutosomalGC
	// * gc contents less than this will not be included in the regression model
	// * @param maximimumAutosomalGC
	// * gc contents greater than this will not be included in the regression model
	// * @param minIntensity
	// * intensities less than this will not be included in the regression model
	// * @param maxIntensity
	// * intensities greater than this will not be included in the regression model
	// * @param regressionDistance
	// * the bin distance in base pairs, sliding window for regression, straight bins for qc
	// * @param skipPerChr
	// * skip this many markers at the beginning and end of each chromosome
	// * @param markerMask
	// * use at own peril, mainly a place holder currently
	// * @param PennCNVGCWF
	// * use the fixed PennCNV GC content array, based on chr 11
	// * @param verbose
	// * report
	// */
	// public GcAdjustor(Project proj, PreparedMarkerSet preparedMarkerSet, GcModel gcModel, double[] markerIntensities, double minimumAutosomalGC, double maximimumAutosomalGC, double minIntensity, double maxIntensity, int regressionDistance, int skipPerChr, boolean[] markerMask, GC_CORRECTION_METHOD correctionMethod, boolean verbose) {
	// super();
	// this.proj = proj;
	// this.gcModel = gcModel;
	// this.markerIntensities = markerIntensities;
	// this.minimumAutosomalGC = minimumAutosomalGC;
	// this.maximimumAutosomalGC = maximimumAutosomalGC;
	// this.minIntensity = minIntensity;
	// this.maxIntensity = maxIntensity;
	// this.regressionDistance = regressionDistance;
	// this.skipPerChr = skipPerChr;
	// this.verbose = verbose;
	// this.fail = false;
	// this.wfPrior = Double.NaN;
	// this.wfPost = Double.NaN;
	// this.gcwfPrior = Double.NaN;
	// this.gcwfPost = Double.NaN;
	// this.correctionMethod = correctionMethod;
	// this.markerMask = markerMask;
	// this.preparedMarkerSet = preparedMarkerSet;
	// populateCurrentData();// Initialize everything needed
	// }

	public CrossValidation getCrossValidation() {
		return crossValidation;
	}

	public void correctIntensities() {
		correctIntensities(null, null);
	}

	public double[] getFullGcs() {
		return fullGcs;
	}

	public double[] getFullIntensity() {
		return fullIntensity;
	}

	/**
	 * Builds a regression model from markers (usually ~3K) selected by {@link GcAdjustor#populateCurrentData()}, the corrected intensities are the residuals after applying the model to all intensities with a gc value
	 */
	public void correctIntensities(String sample, GcAdjustorParameter gcParameters) {
		if (!fail) {
			if (gcParameters != null) {
				this.crossValidation = gcParameters.adjust(correctionMethod, sample, fullIntensity, prepForRegression(fullGcs), verbose, proj.getLog());
			} else {
				this.crossValidation = new CrossValidation(regressionIntensity, prepForRegression(regressionGcs), fullIntensity, prepForRegression(fullGcs), true, false, proj.getLog());
				crossValidation.train();
				crossValidation.computePredictedValues();
				crossValidation.computeResiduals();
			}
			if (crossValidation.analysisFailed()) {
				proj.getLog().reportError("Error - the regression model has failed,  reverting to original intensity values");
				assignOriginalIntensities();
			} else {
				if (verbose) {
					proj.getLog().report("Info - regression model determined to be " + crossValidation.getBetas()[0] + " + " + crossValidation.getBetas()[1] + " * GC content");
				}
				assignCorrectedIntensities();
			}
		} else {
			proj.getLog().reportError("Error - cannot train gc model, reverting to original intensity values");
			assignOriginalIntensities();
		}
	}

	/**
	 * @param computePrior
	 *            compute wf and GCWF from original data
	 * @param computePost
	 *            compute wf and GCWF from corrected data, must call {@link GcAdjustor#correctIntensities()} first if this is flagged
	 */

	public void computeQCMetrics(boolean computePrior, boolean computePost) {

		if (!fail) {
			if (qcIndices == null || qcIndices.length == 0) {
				proj.getLog().reportError("Error - can not compute qc metrics, not enough qc markers were found");
				fail = true;
			}
			if (correctionMethod == GC_CORRECTION_METHOD.PENNCNV_GC && (chr11qcIndices == null || chr11qcIndices.length == 0)) {
				proj.getLog().reportError("Error - can not compute qc metrics using chromosome 11, not enough qc markers were found");
				fail = true;
			}
			if (!fail) {
				if (computePost && (crossValidation == null || crossValidation.analysisFailed())) {
					if (verbose) {
						proj.getLog().report("Warning - intensity correction has not been performed or has failed,  corrected qc metrics will not be computed");
					}
				} else {
					if (computePost) {
						double[] wavesCorrected = getWave(crossValidation.getResiduals(), fullGcs, qcIndices, correctionMethod == GC_CORRECTION_METHOD.PENNCNV_GC ? chr11qcIndices : null, verbose, proj.getLog());
						this.wfPost = wavesCorrected[0];
						this.gcwfPost = wavesCorrected[1];
					}
				}
				if (computePrior) {
					double[] wavesOriginal = getWave(fullIntensity, fullGcs, qcIndices, correctionMethod == GC_CORRECTION_METHOD.PENNCNV_GC ? chr11qcIndices : null, verbose, proj.getLog());
					this.wfPrior = wavesOriginal[0];
					this.gcwfPrior = wavesOriginal[1];
				}
			} else {
				proj.getLog().reportError("Error - cannot compute qc metrics");
			}
		}
	}

	public boolean isFail() {
		return fail;
	}

	public String getQCString() {
		return wfPrior + "\t" + wfPost + "\t" + gcwfPrior + "\t" + gcwfPost;
	}

	public String getAnnotatedQCString() {
		return "WF_PRIOR: " + wfPrior + "\tWF_POST: " + wfPost + "\tGCWF_PRIOR" + gcwfPrior + "\tGCWF_POST" + gcwfPost;
	}

	public double[] getCorrectedIntensities() {
		if (fail) {
			return markerIntensities;
		} else {
			return correctedIntensities;
		}
	}

	public double getWfPrior() {
		return wfPrior;
	}

	public double getWfPost() {
		return wfPost;
	}

	public double getGcwfPrior() {
		return gcwfPrior;
	}

	public double getGcwfPost() {
		return gcwfPost;
	}

	/**
	 * We try to gather all we can in one pass, a long loop.... The loop iterates over markers by chromosome, and gathers info for correcting and qc
	 */
	private void populateCurrentData() {
		if (gcModel == null) {
			fail = true;
			proj.getLog().reportError("Error - a gcModel was not supplied");
		}
		if (markerIntensities.length != proj.getMarkerNames().length) {
			fail = true;
			proj.getLog().reportError("Error - the intensity data array must represent every marker in the project");
		} else {

			preparedMarkerSet = preparedMarkerSet == null ? new PreparedMarkerSet(proj.getMarkerSet()) : preparedMarkerSet;
			String[] markers = preparedMarkerSet.getMarkerNames();
			int[][] indicesByChr = preparedMarkerSet.getIndicesByChr();
			byte[] chrs = preparedMarkerSet.getChrs();
			int[] positions = preparedMarkerSet.getPositions();

			ArrayList<Double> tmpRegressGcs = new ArrayList<Double>(3000);// usually around 3K, stores info for regression model
			ArrayList<Double> tmpRegressIntensity = new ArrayList<Double>(3000);// usually around 3K, stores info for regression model
			ArrayList<Double> tmpFullGcs = new ArrayList<Double>(markers.length);// can't be more than this, stores data for full correction
			ArrayList<Double> tmpFullIntensity = new ArrayList<Double>(markers.length);// can't be more than this, stores data for full correction
			ArrayList<int[]> tmpQcIndices = new ArrayList<int[]>(3000);// refers to the index matched gc/LRR values
			ArrayList<int[]> tmpchr11qcIndices = new ArrayList<int[]>(3000);// refers to the index matched gc/LRR values on chromosome 11 only
			ArrayList<Integer> tmpCorrectedIndices = new ArrayList<Integer>(markers.length);// refers to all indices we will correct

			int qcIndex = -1;
			int numPossibleBins = 0;

			for (int i = 0; i < indicesByChr.length; i++) {
				int currentRegressDistance = 0;
				int currentBin = 0;
				ArrayList<Integer> tmpCurrentBin = new ArrayList<Integer>(1000);
				ArrayList<Integer> tmpCurrentBinChr11 = new ArrayList<Integer>(1000);// only used for chromosome 11...for defualt pennCNV behavior, and track even if not needed

				for (int j = 0; j < indicesByChr[i].length; j++) {
					double gc = gcModel.getGcFor(markers[indicesByChr[i][j]]);// returns NaN if not present in the gcModel file
					if (!Double.isNaN(gc) && !Double.isNaN(markerIntensities[indicesByChr[i][j]])) {// we take everything we can for correcting after the model is built
						tmpFullGcs.add(gc);
						tmpFullIntensity.add(markerIntensities[indicesByChr[i][j]]);
						tmpCorrectedIndices.add(indicesByChr[i][j]);
						qcIndex++;// qc indices refer to every marker that has a valid LRR AND GC value
					}
					if (chrs[indicesByChr[i][j]] > 0 && chrs[indicesByChr[i][j]] < 23) {// build from autosomal only
						if (j >= skipPerChr && j <= indicesByChr[i].length - skipPerChr) {// allows to skip first few and last few markers per chromosome
							if (markerMask == null || markerMask[indicesByChr[i][j]]) {// so we can skip poor quality markers/customize the qc
								// note that this binning procedure is slightly different than the regression indices
								if (positions[indicesByChr[i][j]] >= regressionDistance * currentBin && positions[indicesByChr[i][j]] < regressionDistance * (1 + currentBin)) {
									if (!Double.isNaN(gc) && !Double.isNaN(markerIntensities[indicesByChr[i][j]])) {
										tmpCurrentBin.add(qcIndex);
										if (chrs[indicesByChr[i][j]] == 11) {
											tmpCurrentBinChr11.add(qcIndex);
										}
									}
								} else {// On to the next bin for qc, add indices if appropriate

									if (chrs[indicesByChr[i][j]] == 11) {
										if (tmpCurrentBinChr11.size() > DEFUALT_NUM_SNP_MAD) {
											tmpchr11qcIndices.add(Array.toIntegerArray(tmpCurrentBinChr11));
										} else {
											tmpchr11qcIndices.add(null);
										}
									}
									if (tmpCurrentBin.size() > DEFUALT_NUM_SNP_MAD) {// skip small sized bins
										tmpQcIndices.add(Array.toIntegerArray(tmpCurrentBin));

									}

									while (positions[indicesByChr[i][j]] < regressionDistance * (currentBin + 1)) {// scan to the next bin to use, may be necesarry on less dense arrays
										numPossibleBins++;
										currentBin++;
										proj.getLog().report("Warning - found a " + ext.prettyUpDistance(regressionDistance, 2) + " sliding window without any markers");
									}
									numPossibleBins++;
									currentBin++;

									tmpCurrentBin = new ArrayList<Integer>(1000);
									tmpCurrentBinChr11 = new ArrayList<Integer>(1000);
									// TODO, take this out if we start getting too far from PennCNV, should be more accurate though
									// if (!Double.isNaN(gc) && !Double.isNaN(markerIntensities[indicesByChr[i][j]])) {// Note, PennCNV will not use this marker, they skip instead. We'll take it
									// tmpCurrentBin.add(qcIndex);
									// if (chrs[indicesByChr[i][j]] == 11) {
									// tmpCurrentBinChr11.add(qcIndex);
									// }
									// }
								}
								if (!Double.isNaN(gc) && !Double.isNaN(markerIntensities[indicesByChr[i][j]])) {// for populating regression model
									if (positions[indicesByChr[i][j]] - currentRegressDistance > regressionDistance) {
										if (useMarker(markerIntensities[indicesByChr[i][j]], gc, minimumAutosomalGC, maximimumAutosomalGC, minIntensity, maxIntensity)) {
											// The gc and intensity will be used for the correction
											// if(!PennCNVGCWF||chrs[indicesByChr[i][j]] == 11){
											tmpRegressGcs.add(gc);
											tmpRegressIntensity.add(markerIntensities[indicesByChr[i][j]]);
											// }
										}
										currentRegressDistance = positions[indicesByChr[i][j]];
									}
								}
							}
						}
					}
				}
				if (tmpCurrentBin.size() > 0) {
					numPossibleBins++;
				}
				if (tmpCurrentBin.size() > DEFUALT_NUM_SNP_MAD) {// add any leftovers
					tmpQcIndices.add(Array.toIntegerArray(tmpCurrentBin));
					if (chrs[indicesByChr[i][0]] == 11) {
						tmpchr11qcIndices.add(Array.toIntegerArray(tmpCurrentBinChr11));
					}
				}
				tmpCurrentBin = new ArrayList<Integer>(1000);// reset the bins
				tmpCurrentBinChr11 = new ArrayList<Integer>(1000);// reset the bins
			}
			if (tmpRegressGcs.size() == 0) {
				fail = true;
				proj.getLog().reportError("Error - could not find any autosomal markers to train the regression model");
			} else {
				this.regressionGcs = Array.toDoubleArray(tmpRegressGcs);
				this.regressionIntensity = Array.toDoubleArray(tmpRegressIntensity);
				this.fullGcs = Array.toDoubleArray(tmpFullGcs);
				this.fullIntensity = Array.toDoubleArray(tmpFullIntensity);
				this.correctedIndices = Array.toIntegerArray(tmpCorrectedIndices);
				if (verbose) {
					proj.getLog().report("Info - using " + regressionIntensity.length + " of " + markers.length + " markers for regression model");
					proj.getLog().report("Info - " + fullIntensity.length + " of " + markers.length + " markers had a valid gc and a valid LRR for correction ");
				}
				if (regressionIntensity.length == 0) {
					proj.getLog().reportError("Error - did not find enough valid markers for regression");
					fail = true;
				}

				this.qcIndices = new int[tmpQcIndices.size()][];
				this.chr11qcIndices = new int[tmpchr11qcIndices.size()][];
				for (int i = 0; i < qcIndices.length; i++) {
					qcIndices[i] = tmpQcIndices.get(i);
				}
				for (int i = 0; i < chr11qcIndices.length; i++) {
					chr11qcIndices[i] = tmpchr11qcIndices.get(i);
				}
				if (verbose) {
					proj.getLog().report("Info - detected " + qcIndices.length + " " + ext.prettyUpDistance(regressionDistance, 2) + " sliding windows with >" + DEFUALT_NUM_SNP_MAD + " markers out of a possible " + numPossibleBins + " windows to compute WF" + (correctionMethod == GC_CORRECTION_METHOD.PENNCNV_GC ? "" : " and GCWF"));
				}
			}
		}
	}

	/**
	 * @param intensities
	 *            gc matched intensities
	 * @param gcs
	 *            intensity matched gcs
	 * @param WFbins
	 *            bins for WF , and GCWF if pennCNVGCBins is not supplied should reflect all autosomal marker bins with valid gc/intensity
	 * @param pennCNVGCBins
	 *            (optional, can be null) if provided, this should reflect chromosome 11 bins
	 * @param verbose
	 *            report some stuff
	 * @param log
	 * @return
	 */
	private static double[] getWave(double[] intensities, double[] gcs, int[][] WFbins, int[][] pennCNVGCBins, boolean verbose, Logger log) {
		double[] waves = new double[2];// Organizes as WF, GCWF
		Arrays.fill(waves, Double.NaN);
		ArrayList<Double> medianIntensity = new ArrayList<Double>();
		ArrayList<Double> medianGc = new ArrayList<Double>();
		for (int i = 0; i < WFbins.length; i++) {
			medianIntensity.add(Array.median(Array.subArray(intensities, WFbins[i])));
			medianGc.add(Array.median(Array.subArray(gcs, WFbins[i])));
		}

		double wf = Array.mad(Array.toDoubleArray(medianIntensity));
		if (pennCNVGCBins != null) {// Used for PennCNV bins, if not supplied we use what we found above
			if (pennCNVGCBins.length != DEFUALT_PENNCNV_CHR11_GC_BINS.length) {
				log.reportError("Error - default PennCNV GC bins and current data do not match up, computing using full autosomal bins instead");
			} else {
				medianIntensity = new ArrayList<Double>();
				medianGc = new ArrayList<Double>();
				for (int i = 0; i < pennCNVGCBins.length; i++) {
					if (pennCNVGCBins[i] != null) {
						medianGc.add(DEFUALT_PENNCNV_CHR11_GC_BINS[i]);
						medianIntensity.add(Array.median(Array.subArray(intensities, pennCNVGCBins[i])));
					}
				}
			}
			if (verbose) {
				log.report("Info - computing GCWF using " + medianGc.size() + " elements");
			}
		}
		double cc = stats.Correlation.Pearson(Array.toDoubleArray(medianIntensity), Array.toDoubleArray(medianGc))[0];
		waves[0] = cc > 0 ? -1 * wf : wf;
		waves[1] = waves[0] * Math.abs(cc);
		return waves;
	}

	/**
	 * When correction has failed, set to original
	 */
	private void assignOriginalIntensities() {
		correctedIntensities = markerIntensities;
	}

	/**
	 * Assign indices we have corrected to their new values
	 */
	private void assignCorrectedIntensities() {
		assignOriginalIntensities();
		for (int i = 0; i < correctedIndices.length; i++) {
			correctedIntensities[correctedIndices[i]] = crossValidation.getResiduals()[i];
		}
	}

	private static double[][] prepForRegression(double[] toPrep) {
		double[][] prepped = new double[toPrep.length][1];
		for (int i = 0; i < toPrep.length; i++) {
			prepped[i][0] = toPrep[i];
		}
		return prepped;
	}

	private static boolean useMarker(double intensity, double gc, double minimumAutosomalGC, double maximimumAutosomalGC, double minIntensity, double maxIntensity) {
		if (Double.isNaN(intensity)) {
			return false;
		}
		if (Double.isNaN(gc)) {
			return false;
		}
		if (gc <= minimumAutosomalGC || gc >= maximimumAutosomalGC) {
			return false;
		}
		if (intensity <= minIntensity || intensity >= maxIntensity) {
			return false;
		}
		return true;
	}

	/**
	 * Put this in a separate class so it can be passed between samples
	 *
	 */
	public static class GcModel implements Serializable {
		private static final long serialVersionUID = 1L;
		public static final String[] GC_HEADER = { "Name", "Chr", "Position", "GC" };
		private String[] markers;
		private byte[] chrs;
		private int[] positions;
		private double[] gcs;
		private Hashtable<String, Integer> index = new Hashtable<String, Integer>();
		private Logger log;

		public GcModel(GcModel gcmodel) {
			this.markers = gcmodel.markers;
			this.chrs = gcmodel.chrs;
			this.positions = gcmodel.positions;
			this.gcs = gcmodel.gcs;
			this.index = gcmodel.index;
			this.log = gcmodel.log;
		}

		public GcModel(String[] markers, byte[] chrs, int[] positions, double[] gcs, Hashtable<String, Integer> index, Logger log) {
			super();
			this.markers = markers;
			this.chrs = chrs;
			this.positions = positions;
			this.gcs = gcs;
			this.index = index;
			this.log = log;
		}

		public byte[] getChrs() {
			return chrs;
		}

		public int[] getPositions() {
			return positions;
		}

		public Logger getLog() {
			return log;
		}

		public double[] getGcsFor(String[] markers) {
			double[] gcs = new double[markers.length];
			for (int i = 0; i < gcs.length; i++) {
				gcs[i] = getGcFor(markers[i]);
			}
			return gcs;
		}

		public double getGcFor(String marker) {
			if (hasMarker(marker)) {
				return gcs[index.get(marker)];
			} else {
				return Double.NaN;
			}
		}

		public boolean hasMarker(String marker) {
			return index.containsKey(marker);
		}

		public String[] getMarkers() {
			return markers;
		}

		public void Serialize(String fullPathToGcSer) {
			Files.writeSerial(this, fullPathToGcSer);
		}

		public static GcModel loadSerial(String fullPathToGcSer) {
			return (GcModel) Files.readSerial(fullPathToGcSer);
		}

		/**
		 * This generates a gcModel file that takes the bin size to be explicitly the window around the snp,
		 */
		public static GcModel generateSnpWindowModel(Project proj, int snpWindow) {
			String refGenome = proj.REFERENCE_GENOME_FASTA_FILENAME.getValue();
			if (!Files.exists(refGenome)) {
				proj.getLog().reportFileNotFound(refGenome);
				return null;
			} else {
				ReferenceGenome referenceGenome = new ReferenceGenome(refGenome, proj.getLog());
				Hashtable<String, Integer> indices = proj.getMarkerIndices();
				MarkerSet markerSet = proj.getMarkerSet();
				double[] gcs = new double[markerSet.getMarkerNames().length];
				for (int i = 0; i < markerSet.getMarkerNames().length; i++) {
					Segment seg = new Segment(markerSet.getChrs()[i], markerSet.getPositions()[i], markerSet.getPositions()[i]).getBufferedSegment(snpWindow);
					double gc = referenceGenome.getGCContentFor(seg);
					gcs[i] = gc;
				}

				return new GcModel(markerSet.getMarkerNames(), markerSet.getChrs(), markerSet.getPositions(), gcs, indices, proj.getLog());
			}
		}

		/**
		 * uses {@link GcModel#generateFromReferenceGenome(Project, String, String, int)} but with the default gc model window size of 100*5120
		 */
		public static boolean generateFromReferenceGenome(Project proj, String fullPathToReferenceGenome, String fullPathToOutputModel) {
			return generateFromReferenceGenome(proj, fullPathToReferenceGenome, fullPathToOutputModel, DEFUALT_GC_MODEL_WINDOW_SNP, DEFUALT_GC_MODEL_WINDOW_GC);
		}

		/**
		 * @param proj
		 * @param fullPathToReferenceGenome
		 *            reference fasta
		 * @param fullPathToOutputModel
		 *            the output gc model file
		 * @param windowPerSNP
		 *            bp window around each marker to compute gc content in
		 * @return
		 */
		public static boolean generateFromReferenceGenome(Project proj, String fullPathToReferenceGenome, String fullPathToOutputModel, int windowPerSNP, int windowPerGC) {
			if (!Files.exists(fullPathToReferenceGenome)) {
				proj.getLog().reportFileNotFound(fullPathToReferenceGenome);
				return false;
			} else if (Files.exists(fullPathToOutputModel)) {
				proj.getLog().reportTimeWarning(fullPathToOutputModel + " exists, will not create again...");
				return false;
			} else {
				proj.getLog().reportTimeInfo("Generating gc model file at " + fullPathToOutputModel);
				ReferenceGenome referenceGenome = new ReferenceGenome(fullPathToReferenceGenome, proj.getLog());
				LocusSet<Segment> bins = referenceGenome.getBins(windowPerGC);
				proj.getLog().reportTimeInfo("Computing gc content for " + bins.getLoci().length + " bins of " + windowPerGC + " bp");
				double[] gcContents = new double[bins.getLoci().length];
				for (int i = 0; i < bins.getLoci().length; i++) {
					if (i % 1000 == 0) {
						proj.getLog().reportTimeInfo("Queried " + (i + 1) + " bins for gc content");
					}
					gcContents[i] = referenceGenome.getGCContentFor(bins.getLoci()[i]);
				}

				MarkerSet markerSet = proj.getMarkerSet();
				String[] markerNames = markerSet.getMarkerNames();
				int[] positions = markerSet.getPositions();
				byte[] chrs = markerSet.getChrs();
				try {
					PrintWriter writer = new PrintWriter(new FileWriter(fullPathToOutputModel));
					writer.println("Name\tChr\tPosition\tGC");
					for (int i = 0; i < markerNames.length; i++) {
						Segment bufferedMarkerSeg = new Segment(chrs[i], positions[i], positions[i]).getBufferedSegment(windowPerSNP);
						String chr = chrs[i] + "";
						if (chrs[i] == 23) {
							chr = "X";
						} else if (chrs[i] == 24) {
							chr = "Y";
						} else if (chrs[i] == 25) {
							chr = "XY";
						} else if (chrs[i] == 26) {
							chr = "Un";
						}
						double gcContent = Double.NaN;
						int[] overlapping = bins.getOverlappingIndices(bufferedMarkerSeg);

						if (overlapping != null && overlapping.length > 0) {
							Segment[] binOverlap = Array.subArray(bins.getLoci(), overlapping);
							double gcs = 0;
							int bps = 0;
							for (int j = 0; j < overlapping.length; j++) {
								bps += binOverlap[j].getSize();
								gcs += gcContents[overlapping[j]] * binOverlap[j].getSize();
							}
							gcContent = (double) gcs / bps;
						} else if (chrs[i] > 0) {
							String error = "BUG: Did not find any overlapping bins for marker " + markerNames[i] + " for window search " + bufferedMarkerSeg.getUCSClocation();
							proj.getLog().reportTimeError(error);
							writer.close();
							throw new IllegalStateException(error);
						}

						writer.println(markerNames[i] + "\t" + chr + "\t" + positions[i] + "\t" + gcContent);
					}
					writer.close();
					return true;

				} catch (Exception e) {
					proj.getLog().reportError("Error writing to " + fullPathToOutputModel);
					proj.getLog().reportException(e);
					return false;
				}

			}
		}

		public static GcModel populateFromFile(String fullPathToGcModel, boolean verbose, Logger log) {
			ArrayList<String> markers = new ArrayList<String>();
			ArrayList<Byte> chrs = new ArrayList<Byte>();
			ArrayList<Integer> positions = new ArrayList<Integer>();
			ArrayList<Double> gcs = new ArrayList<Double>();
			Hashtable<String, Integer> index = new Hashtable<String, Integer>();
			String fullPathToGcSer = ext.rootOf(fullPathToGcModel, false) + ".gcmodel.ser";
			if (Files.exists(fullPathToGcSer)) {
				log.report("Info - loading gc model file " + fullPathToGcSer);
				return loadSerial(fullPathToGcSer);
			}
			if (!Files.exists(fullPathToGcModel)) {
				log.reportError("Error - could not find gc model file " + fullPathToGcModel);
				return null;
			}
			log.report("Info - parsing gc model file " + fullPathToGcModel);

			int[] indices = ext.indexFactors(Files.getHeaderOfFile(fullPathToGcModel, log), GC_HEADER, true, false);
			if (Array.countIf(indices, -1) > 0) {
				log.reportError("Error - could not find correct header for gc model file " + fullPathToGcModel);
				log.reportError("		 - header must be:" + Array.toStr(GC_HEADER));
				return null;
			} else {
				BufferedReader reader;
				try {
					reader = Files.getAppropriateReader(fullPathToGcModel);
					String[] line = reader.readLine().trim().split("[\\s]+");// header
					int lineNum = 0;
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						try {
							String marker = line[indices[0]];
							byte chr = Positions.chromosomeNumber(line[indices[1]], log);
							int pos = Integer.parseInt(line[indices[2]]);
							double gc = Double.parseDouble(line[indices[3]]);
							markers.add(marker);
							chrs.add(chr);
							positions.add(pos);
							gcs.add(gc);
							index.put(marker, lineNum);
							lineNum++;
						} catch (NumberFormatException nfe) {
							if (verbose) {
								log.reportError("Error - found invalid number format on line " + Array.toStr(line) + " , skipping");
							}
						}
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					log.reportError("Error: file \"" + fullPathToGcModel + "\" not found in current directory");
					return null;
				} catch (IOException ioe) {
					log.reportError("Error reading file \"" + fullPathToGcModel + "\"");
					return null;
				}
				if (markers.size() == 0) {
					log.reportError("Error - did not find any valid markers in gc model file " + fullPathToGcModel);
					return null;
				} else {
					log.report("Info - loaded " + markers.size() + " markers from gc model file " + fullPathToGcModel);
					GcModel gcModel = new GcModel(Array.toStringArray(markers), Array.toByteArray(chrs), Array.toIntegerArray(positions), Array.toDoubleArray(gcs), index, log);
					gcModel.Serialize(fullPathToGcSer);
					return gcModel;
				}
			}
		}
	}

	/**
	 * Similar to {@link GcAdjustor#getComputedAdjustor()} except it takes a {@link Sample} object as input
	 */
	public static GcAdjustor getComputedAdjustor(Project proj, Sample sample, PreparedMarkerSet preparedMarkerSet, GcModel gcModel, GC_CORRECTION_METHOD correctionMethod, boolean computePrior, boolean computePost, boolean verbose) {
		return getComputedAdjustor(proj, preparedMarkerSet, sample.getLRRs(), gcModel, correctionMethod, computePrior, computePost, verbose);
	}

	public static GcAdjustor getComputedAdjustor(Project proj, PreparedMarkerSet preparedMarkerSet, float[] markerIntensities, GcModel gcModel, GC_CORRECTION_METHOD correctionMethod, boolean computePrior, boolean computePost, boolean verbose) {
		return getComputedAdjustor(proj, null, null, preparedMarkerSet, markerIntensities, gcModel, correctionMethod, computePrior, computePost, verbose);
	}

	/**
	 * @param proj
	 * @param sample
	 *            the current sample, only neccesary if a {@link GcAdjustorParameter} is provided
	 * @param gcParameters
	 *            if a valid {@link GcAdjustorParameter} is provided, the regression model will be skipped and the betas from this object will be used
	 * @param preparedMarkerSet
	 *            can be null and will be autogenerated if is
	 * @param markerIntensities
	 *            float[] that is converted to double[] for correction
	 *
	 * @param gcModel
	 *            a valid model
	 * @param pennCNVGCWF
	 *            use the PennCNV defualt chr 11 bins for GCWF calculation
	 * @param computePrior
	 *            compute WF and GCWF prior to gc correction
	 * @param computePost
	 *            compute WF and GCWF post gc correction
	 * @param verbose
	 *            reports things akin to PennCNV
	 * @return a {@link GcAdjustor} object with corrected intensities and qc metrics
	 */
	public static GcAdjustor getComputedAdjustor(Project proj, String sample, GcAdjustorParameter gcParameters, PreparedMarkerSet preparedMarkerSet, float[] markerIntensities, GcModel gcModel, GC_CORRECTION_METHOD correctionMethod, boolean computePrior, boolean computePost, boolean verbose) {
		GCAdjustorBuilder builder = new GCAdjustorBuilder();
		builder.correctionMethod(correctionMethod);
		builder.verbose(verbose);
		GcAdjustor gcAdjustor = builder.build(proj, preparedMarkerSet, gcModel, Array.toDoubleArray(markerIntensities));
		gcAdjustor.correctIntensities(sample, gcParameters);
		gcAdjustor.computeQCMetrics(computePrior, computePost);
		return gcAdjustor;
	}

	public static GcAdjustor getComputedAdjustor(Project proj, GCAdjustorBuilder builder, String sample, GcAdjustorParameter gcParameters, PreparedMarkerSet preparedMarkerSet, float[] markerIntensities, GcModel gcModel, boolean computePrior, boolean computePost) {
		GcAdjustor gcAdjustor = builder.build(proj, preparedMarkerSet, gcModel, Array.toDoubleArray(markerIntensities));
		gcAdjustor.correctIntensities(sample, gcParameters);
		gcAdjustor.computeQCMetrics(computePrior, computePost);
		return gcAdjustor;
	}

	public static void test(Project proj, String fullPathToGcModel, String fullPathToFileOfTestSamples) {
		GcModel gcModel = GcModel.populateFromFile(fullPathToGcModel, false, proj.getLog());
		String[] samplesToTest = Array.subArray(proj.getSamples(), proj.getSamplesToInclude(fullPathToFileOfTestSamples));
		String fileTest = proj.PROJECT_DIRECTORY.getValue() + "testGCWF.txt";
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(fileTest));
			for (int i = 0; i < samplesToTest.length; i++) {
				long time = System.currentTimeMillis();
				Sample samp = proj.getFullSampleFromRandomAccessFile(samplesToTest[i]);
				proj.getLog().report("Testing sample " + samp.getSampleName());
				GCAdjustorBuilder builder = new GCAdjustorBuilder();
				GcAdjustor gcAdjusterNew = builder.build(proj, null, gcModel, Array.toDoubleArray(samp.getLRRs()));
				gcAdjusterNew.correctIntensities();
				gcAdjusterNew.computeQCMetrics(true, true);
				writer.println(samp.getSampleName() + "\t" + gcAdjusterNew.getQCString() + "\t" + ext.getTimeElapsed(time));
			}
			writer.close();
		} catch (Exception e) {
			proj.getLog().reportError("Error writing to " + fileTest);
			proj.getLog().reportException(e);
		}
	}

	public static void main(String[] args) {
		Project proj = new Project(null, null, false);
		String fullPathToGcModel = "D:/data/gedi_gwas/CompareGC_CORRECTION/custom.gcmodel";
		String fullPathToFileOfTestSamples = "D:/data/gedi_gwas/CompareGC_CORRECTION/testGCCorrection.txt";
		test(proj, fullPathToGcModel, fullPathToFileOfTestSamples);
	}

	public static class GCAdjustorBuilder {
		private double minimumAutosomalGC = DEFAULT_MINAUTOSOMALGC;
		private double maximimumAutosomalGC = DEFAULT_MAXAUTOSOMALGC;
		private double minIntensity = DEFAULT_MIN_DATA_VALUE;
		private double maxIntensity = DEFAULT_MAX_DATA_VALUE;
		private int regressionDistance = DEFAULT_REGRESSION_DISTANCE[0];
		private int skipPerChr = DEFAULT_SKIP_PER_CHR[0];
		private boolean[] markerMask = null;
		private boolean verbose = false;
		private GC_CORRECTION_METHOD correctionMethod = GC_CORRECTION_METHOD.GENVISIS_GC;

		public GCAdjustorBuilder minimumAutosomalGC(double minimumAutosomalGC) {
			this.minimumAutosomalGC = minimumAutosomalGC;
			return this;
		}

		public GCAdjustorBuilder maximimumAutosomalGC(double maximimumAutosomalGC) {
			this.maximimumAutosomalGC = maximimumAutosomalGC;
			return this;
		}

		public GCAdjustorBuilder minIntensity(double minIntensity) {
			this.minIntensity = minIntensity;
			return this;
		}

		public GCAdjustorBuilder maxIntensity(double maxIntensity) {
			this.maxIntensity = maxIntensity;
			return this;
		}

		public GCAdjustorBuilder regressionDistance(int regressionDistance) {
			this.regressionDistance = regressionDistance;
			return this;
		}

		public GCAdjustorBuilder skipPerChr(int skipPerChr) {
			this.skipPerChr = skipPerChr;
			return this;
		}

		public GCAdjustorBuilder markerMask(boolean[] markerMask) {
			this.markerMask = markerMask;
			return this;
		}

		public GCAdjustorBuilder verbose(boolean verbose) {
			this.verbose = verbose;
			return this;
		}

		public GCAdjustorBuilder correctionMethod(GC_CORRECTION_METHOD correctionMethod) {
			this.correctionMethod = correctionMethod;
			return this;
		}

		public GcAdjustor build(Project proj, PreparedMarkerSet markerSet, GcModel gcModel, double[] markerIntensities) {
			return new GcAdjustor(this, proj, markerSet, gcModel, markerIntensities);
		}
	}

	private GcAdjustor(GCAdjustorBuilder builder, Project proj, PreparedMarkerSet markerSet, GcModel gcModel, double[] markerIntensities) {
		this.minimumAutosomalGC = builder.minimumAutosomalGC;
		this.maximimumAutosomalGC = builder.maximimumAutosomalGC;
		this.minIntensity = builder.minIntensity;
		this.maxIntensity = builder.maxIntensity;
		this.regressionDistance = builder.regressionDistance;
		this.skipPerChr = builder.skipPerChr;
		this.markerMask = builder.markerMask;
		this.verbose = builder.verbose;
		this.correctionMethod = builder.correctionMethod;
		this.proj = proj;
		this.gcModel = gcModel;
		this.markerIntensities = markerIntensities;
		this.fail = false;
		this.wfPrior = Double.NaN;
		this.wfPost = Double.NaN;
		this.gcwfPrior = Double.NaN;
		this.gcwfPost = Double.NaN;
		this.preparedMarkerSet = markerSet;
		populateCurrentData();// Initialize everything needed
	}
}

// /**
// * @return an array with WF and GCWF prior to and post gc correction
// */
// public static double[] getQCMetrics(Project proj, Sample sample, GcModel gcModel, boolean pennCNVGCWF, boolean verbose) {
// double[] waves = new double[4];// stores WF prior, WF post, GCWF prior, GCWF post
// Arrays.fill(waves, Double.NaN);
// if (gcModel != null) {
// GcAdjustor gcAdjustor = getComputedAdjustor(proj, sample, gcModel, pennCNVGCWF, true, true, verbose);
// waves[0] = gcAdjustor.getWfPrior();
// waves[1] = gcAdjustor.getWfPost();
// waves[2] = gcAdjustor.getGcwfPrior();
// waves[3] = gcAdjustor.getGcwfPost();
// }
// return waves;
// }