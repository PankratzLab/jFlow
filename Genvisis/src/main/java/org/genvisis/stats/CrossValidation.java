package org.genvisis.stats;

import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.stats.LeastSquares.LS_TYPE;

// TODO, verify logistic validation....
public class CrossValidation {

	private double Rsquare, SSerr, avgSEofBs, fullModelR2, fullModelSSerr;
	private double[] SEofBs, stats, sigs;
	private double[] train_deps;
	private double[][] train_indeps;
	private double[] val_deps;
	private double[][] val_indeps;
	private double[] betas, predicteds, residuals;
	private final boolean logistic;
	private boolean analysisFailed;
	private final boolean verbose;
	private final LS_TYPE lType;
	private RegressionModel model;
	private final Logger log;

	/**
	 * @param train_deps dependent variables for training data
	 * @param train_indeps independent variables for training data
	 * @param validation_deps dependent variables for validation data
	 * @param validation_indeps independent variables for validation data
	 * @param verbose verbose reporting of the training regression
	 * @param svdRegression use an SVD based regression
	 * @param log Note training and validation dependent variables can have different lengths,
	 *        however, the number of independent variables must be the same
	 */
	public CrossValidation(	double[] train_deps, double[][] train_indeps, double[] validation_deps,
													double[][] validation_indeps, boolean verbose, LS_TYPE lType,
													Logger log) {
		super();
		this.train_deps = train_deps;
		this.train_indeps = train_indeps;
		val_deps = validation_deps;
		val_indeps = validation_indeps;
		predicteds = new double[validation_deps.length];
		residuals = new double[validation_deps.length];
		SSerr = 0;
		Rsquare = Double.NaN;
		analysisFailed = false;
		this.verbose = verbose;
		this.lType = lType;
		logistic = isLogistic(train_deps == null ? val_deps : train_deps);// in case we are just
																																			// applying betas or something
																																			// like that
		this.log = log;
	}

	/**
	 * We only need to obtain betas from the training data using the appropriate regression method
	 */
	public void train() {
		if (verifyData()) {
			if (logistic) {
				// TODO, logistic validation can be done in a similar manner? , please check before removing
				// this?
				model = new LogisticRegression(train_deps, train_indeps, false, verbose);
				log.reportError("Error - currently can only handle linear cross-validations, I think");
				analysisFailed = true;
			} else {
				model = new LeastSquares(train_deps, train_indeps, null, false, verbose, lType);
			}
			if (model.analysisFailed()) {
				analysisFailed = true;
				log.reportError("Error - regression model has failed");
			} else {
				betas = model.getBetas();
				SEofBs = model.getSEofBs();
				avgSEofBs = Array.mean(SEofBs);
				Rsquare = model.getRsquare();
				stats = model.getStats();
				sigs = model.getSigs();
			}
		} else {
			analysisFailed = true;
		}
	}

	/**
	 * Computes the predicted values for the validation set from the betas obtained in the training
	 * set
	 */
	public void computePredictedValues() {
		if (!analysisFailed) {
			int N = val_deps.length;
			int M = val_indeps[0].length;
			double[][] X = new double[M + 1][N];
			for (int i = 0; i < N; i++) {
				X[0][i] = 1;// constant term
				for (int j = 1; j <= M; j++) {
					X[j][i] = val_indeps[i][j - 1];
				}
			}
			for (int i = 0; i < N; i++) {
				predicteds[i] = betas[0];
				for (int j = 1; j <= M; j++) {
					predicteds[i] += betas[j] * X[j][i];// if X[j][i] is NaN, predicteds[i] and residuals[i]
																							// will be NaN as well
				}
				if (logistic) {
					predicteds[i] = Math.exp(predicteds[i]) / (1 + Math.exp(predicteds[i]));// TODO
				}
			}
		} else {
			log.reportError("Error - could not train data set, cannot compute predicted values");
		}

	}

	/**
	 * Computes the residuals and residual sum of squares for the validation set
	 */
	public void computeResiduals() {
		if (!analysisFailed) {
			int droppedNaN = 0;
			for (int i = 0; i < val_deps.length; i++) {
				residuals[i] = val_deps[i] - predicteds[i];
				if (Double.isNaN(val_deps[i]) || Double.isNaN(predicteds[i])) {
					droppedNaN++;
					residuals[i] = Double.NaN;
				} else {
					SSerr += Math.pow(residuals[i], 2);// TODO is this valid for logistic?
				}
			}
			if (droppedNaN > 0 && verbose) {
				log.report("Warning - "	+ droppedNaN + " "
										+ (droppedNaN > 1 ? "individuals were" : "individual was")
										+ " not included in the residual sum of squares calculation due to missing independent or dependent variables");
			}
		} else {
			log.reportError("Error - could not train data set, cannot compute residuals values");
			SSerr = Double.NaN;
		}
	}

	public static boolean isLogistic(double[] test) {
		Hashtable<Double, Double> tmp = new Hashtable<Double, Double>();
		for (double element : test) {
			tmp.put(element, element);
			if (tmp.size() > 2) {
				return false;
			}
		}
		return true;
	}

	/**
	 * @return true if the training data set and validation data set contain the same (gt 0) number of
	 *         independent variables
	 *         <p>
	 *         Note: add any more data checks needed here
	 */
	public boolean verifyData() {
		boolean verify = true;
		int numIndeps = (train_indeps[0] == null ? 0 : train_indeps[0].length);
		if (numIndeps == 0) {
			verify = false;
			log.reportError("Error - no independant variables were provided");
		}
		if (!verifyEqualLength(train_indeps, numIndeps)) {
			verify = false;
			log.reportError("Error - the number of independent variables to train on must be equal size");
		}
		if (!verifyEqualLength(val_indeps, numIndeps)) {
			verify = false;
			log.reportError("Error - the number of independent variables to validate must be the same as the training set");
		}
		if (train_deps.length < train_indeps[0].length + 1) {
			verify = false;
			log.reportError("Error - must have more observations than individuals in the regression model");
		}
		return verify;
	}

	public double getSSerr() {
		return SSerr;
	}

	public void setBetas(double[] betas) {
		this.betas = betas;
	}

	public double[] getResiduals() {
		return residuals;
	}

	public RegressionModel getModel() {
		return model;
	}

	public boolean analysisFailed() {
		return analysisFailed;
	}

	public double getRsquare() {
		return Rsquare;
	}

	public double[] getSEofBs() {
		return SEofBs;
	}

	public double[] getBetas() {
		return betas;
	}

	public double getAvgSEofBs() {
		return avgSEofBs;
	}

	public double[] getPredicteds() {
		return predicteds;
	}

	// Full model variables must be set from outside, only for storage..
	public double getFullModelR2() {
		return fullModelR2;
	}

	public void setFullModelR2(double fullModelR2) {
		this.fullModelR2 = fullModelR2;
	}

	public double getFullModelSSerr() {
		return fullModelSSerr;
	}

	public void setFullModelSSerr(double fullModelSSerr) {
		this.fullModelSSerr = fullModelSSerr;
	}

	public void setAnalysisFailed(boolean analysisFailed) {
		this.analysisFailed = analysisFailed;
	}

	/**
	 * stats of the training model
	 */
	public double[] getStats() {
		return stats;
	}

	/**
	 * sigs of the training model
	 */
	public double[] getSigs() {
		return sigs;
	}

	/**
	 * So we can save some memory if we return a bunch of validations
	 */
	public void clearInputData(boolean hard) {
		train_deps = new double[0];
		train_indeps = new double[0][0];
		val_deps = new double[0];
		val_indeps = new double[0][0];
		model = null;
		if (hard) {
			predicteds = new double[0];
			residuals = new double[0];
			betas = null;
			SEofBs = null;
		}
	}

	/**
	 * This will compute the residuals and SSerr within a CrossValidation for a single set of training
	 * and validation data
	 */
	public static CrossValidation crossValidate(double[] train_deps, double[][] train_indeps,
																							double[] validation_deps,
																							double[][] validation_indeps, boolean verbose,
																							LS_TYPE lType, Logger log) {
		CrossValidation crossValidation = new CrossValidation(train_deps, train_indeps, validation_deps,
																													validation_indeps, verbose, lType, log);
		crossValidation.train();
		if (!crossValidation.analysisFailed) {
			crossValidation.computePredictedValues();
			crossValidation.computeResiduals();
		}
		return crossValidation;
	}

	/**
	 * In k-fold cross-validation, the original sample is randomly partitioned into k equal size
	 * subsamples.
	 * <p>
	 * Of the k subsamples, a single subsample is retained as the validation data for testing the
	 * model,
	 * <p>
	 * and the remaining k- 1 subsamples are used as training data. -
	 * http://en.wikipedia.org/wiki/Cross-validation_(statistics)#k-fold_cross-validation
	 * <p>
	 * Note: we do not randomly partition the data, we split it up as it comes. If random is desired,
	 * please shuffle before this method
	 * <p>
	 * Note: this is an in-sample cross-validation, the training and test data is created on the fly
	 */
	public static CrossValidation[] kFoldCrossValidate(	double[] deps, double[][] indeps, int kFolds,
																											boolean verbose, LS_TYPE lType, Logger log) {
		if (!foldCheck(deps, kFolds, log)) {
			return new CrossValidation[0];
		}
		int[] chunks = Array.splitUp(deps.length, kFolds);
		CrossValidation[] crossValidations = new CrossValidation[chunks.length];
		boolean[][] folds = getFolds(deps, chunks, log);
		for (int i = 0; i < chunks.length; i++) {
			double[] train_deps = extractDeps(deps, folds[i], true, log);
			double[] val_deps = extractDeps(deps, folds[i], false, log);
			double[][] train_indeps = extractIndeps(indeps, folds[i], true, log);
			double[][] val_indeps = extractIndeps(indeps, folds[i], false, log);
			crossValidations[i] = crossValidate(train_deps, train_indeps, val_deps, val_indeps, verbose,
																					lType, log);
			crossValidations[i].clearInputData(true);// save some memory
		}
		return crossValidations;
	}

	/**
	 * Similar to {@link CrossValidation#kFoldCrossValidate()}, but we keep the validation data
	 * constant
	 */
	public static CrossValidation[] kFoldCrossValidateOutSample(double[] deps, double[][] indeps,
																															double[] val_deps,
																															double[][] val_indeps, int kFolds,
																															boolean verbose, LS_TYPE lType,
																															Logger log) {
		if (!foldCheck(deps, kFolds, log)) {
			return new CrossValidation[0];
		}
		int[] chunks = Array.splitUp(deps.length, kFolds);
		CrossValidation[] crossValidations = new CrossValidation[chunks.length];
		boolean[][] folds = getFolds(deps, chunks, log);
		for (int i = 0; i < chunks.length; i++) {
			double[] train_deps = extractDeps(deps, folds[i], true, log);
			double[][] train_indeps = extractIndeps(indeps, folds[i], true, log);
			crossValidations[i] = crossValidate(train_deps, train_indeps, val_deps, val_indeps, verbose,
																					lType, log);
			crossValidations[i].clearInputData(true);
		}
		return crossValidations;
	}

	/**
	 * If you only care about the average error
	 */
	public static double kfoldAverageSSerr(	double[] deps, double[][] indeps, int kFolds,
																					boolean verbose, LS_TYPE lType, Logger log) {
		CrossValidation[] crossValidations = kFoldCrossValidate(deps, indeps, kFolds, verbose, lType,
																														log);
		return getEstimateError(crossValidations);
	}

	/**
	 * @param crossValidation
	 * @return the average SSerr of non-failing CrossValidations
	 */

	public static double getEstimateError(CrossValidation[] crossValidation) {
		int count = 0;
		double sum = 0;
		for (int i = 0; i < crossValidation.length; i++) {
			if (!crossValidation[i].analysisFailed()) {
				sum += crossValidation[i].getSSerr();
				count++;
			}
		}
		return (count > 0 ? (sum / (count)) : Double.NaN);
	}

	/**
	 * @param crossValidation
	 * @return the average R-squared of non-failing CrossValidations
	 */

	public static double getAverageR2(CrossValidation[] crossValidation) {
		int count = 0;
		double sum = 0;
		for (int i = 0; i < crossValidation.length; i++) {
			if (!crossValidation[i].analysisFailed()) {
				sum += crossValidation[i].getRsquare();
				count++;
			}
		}
		return (count > 0 ? (sum / (count)) : Double.NaN);
	}

	/**
	 * @param crossValidation
	 * @return the average Standared error of all betas for non-failing CrossValidations
	 */

	public static double getAverageSEbetas(CrossValidation[] crossValidation) {
		int count = 0;
		double sum = 0;
		for (int i = 0; i < crossValidation.length; i++) {
			if (!crossValidation[i].analysisFailed()) {
				sum += crossValidation[i].getAvgSEofBs();
				count++;
			}
		}
		return (count > 0 ? (sum / (count)) : Double.NaN);
	}

	// TODO, merge extractDeps and extract Indeps, or use an ext method
	/**
	 * @param deps dependent variables
	 * @param fold a boolean array defining the comparisons
	 * @param train if train, true entries an array of deps corresponding to true entries in fold will
	 *        be returned. Otherwise the false entries will be returned
	 * @param log
	 * @return
	 */
	private static double[] extractDeps(double[] deps, boolean[] fold, boolean train, Logger log) {
		int num = Array.booleanArraySum(fold);
		num = train ? num : fold.length - num;
		double[] extract = new double[num];
		int depIndex = 0;
		for (int i = 0; i < fold.length; i++) {
			if (train && fold[i]) {
				extract[depIndex] = deps[i];
				depIndex++;
			} else if (!train && !fold[i]) {
				extract[depIndex] = deps[i];
				depIndex++;
			}
		}
		return extract;
	}

	// TODO, merge extractDeps and extract Indeps, or use an ext method

	/**
	 * Similar to {@link CrossValidation#extractDeps()}, but for the independent data
	 */
	private static double[][] extractIndeps(double[][] indeps, boolean[] fold, boolean train,
																					Logger log) {
		int num = Array.booleanArraySum(fold);
		num = train ? num : fold.length - num;
		double[][] extract = new double[num][];
		int depIndex = 0;
		for (int i = 0; i < fold.length; i++) {
			if (train && fold[i]) {
				extract[depIndex] = indeps[i];
				depIndex++;
			} else if (!train && !fold[i]) {
				extract[depIndex] = indeps[i];
				depIndex++;
			}
		}
		return extract;
	}

	/**
	 * @param deps
	 * @param chunks we leave the number specified at chunks[i] out of each analysis
	 * @param log
	 *
	 */
	private static boolean[][] getFolds(double[] deps, int[] chunks, Logger log) {
		boolean[][] folds = new boolean[chunks.length][deps.length];
		int depIndex = 0;
		for (int i = 0; i < chunks.length; i++) {
			for (int j = 0; j < chunks[i]; j++) {
				for (int j2 = 0; j2 < chunks.length; j2++) {
					if (j2 != i) {// we leave a chunk out
						folds[j2][depIndex] = true;
					} else {
						folds[j2][depIndex] = false;
					}
				}
				depIndex++;
			}
		}
		return folds;
	}

	private static boolean verifyEqualLength(double[][] data, int length) {
		for (double[] element : data) {
			if (element == null || element.length != length) {
				return false;
			}
		}
		return true;
	}

	private static boolean foldCheck(double[] deps, int kFolds, Logger log) {
		boolean pass = true;
		if (kFolds > deps.length) {
			log.reportError("Error - cannot have number of folds greater than number of dependent variables");
			pass = false;
		}
		if (kFolds < 2) {
			log.reportError("Error - must have at least two folds");
			pass = false;
		}
		return pass;
	}

	/**
	 * Just a test, if we use the same training and validation data, we should always end up with
	 * identical residuals and thus identical SSerr
	 */
	public static void testingStufff(double[] deps, double[][] indeps, String testNum, Logger log) {
		CrossValidation crossValidation = crossValidate(deps, indeps, deps, indeps, false,
																										LS_TYPE.REGULAR, log);
		double[] modelResiduals = crossValidation.getModel().getResiduals();
		double[] extrapResiduals = crossValidation.getResiduals();
		String[] toPrint = new String[modelResiduals.length];
		int modelIndex = 0;
		boolean allEqual = true;
		for (int i = 0; i < extrapResiduals.length; i++) {
			if (!Double.isNaN(extrapResiduals[i])) {
				toPrint[modelIndex] = modelResiduals[modelIndex] + "\t" + extrapResiduals[i];
				if (modelResiduals[modelIndex] != extrapResiduals[i]) {
					allEqual = false;
				}
				modelIndex++;
			}
		}
		if (allEqual) {
			log.report("The extrapolated residuals are identical!!!");
		}
	}
}
