package org.genvisis.cnv.analysis.pca;

import java.io.FileNotFoundException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import org.genvisis.cnv.analysis.pca.PrincipalComponentsResiduals.PrincipalComponentsIterator;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.ExtProjectDataParser;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.Array.BooleanClassifier;
import org.genvisis.common.WorkerTrain.Producer;
import org.genvisis.stats.Correlation;
import org.genvisis.stats.CrossValidation;
import org.genvisis.stats.ICC;
import org.genvisis.stats.LeastSquares.LS_TYPE;

public class CorrectionEvaluator implements Producer<EvaluationResult>, Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private static final String[] STRING_DATA = new String[] { "CENTER" };// CLASS Sex
	private static final String[][] EVAL_MASKS = new String[][] { { "0", "-1", "NaN" }, { "0", "-1" } };
	public static final String[] DOUBLE_DATA = new String[] { "AGE" };
	public static final String[] DOUBLE_DATA_PATTERN = new String[] { "EVAL_DATA", "AGE" };// For Correlation(Spearman by ICC)
	public static final String[] STRING_DATA_PATTERN = new String[] { "EVAL_CLASS", "AGE" };// For ICC
	public static final String[] STRAT_STRING_PATTERN = new String[] { "STRAT_CLASS" };// For ICC

	public static final String[] INDEPS = new String[] { "CLASS=SEX", "AGE" };
	public static final String[] INDEPS_CATS = new String[] { "CENTER" };
	private static final String NO_STRAT = "NO_STRAT";
	// private static final String STRAT_BY = "STRAT_";
	public static final int NUM_PC_SVD_OVERIDE = 160;
	private Project proj;
	private PrincipalComponentsIterator iterator;
	private EvaluationResult[] precomputeds;
	private boolean[][] samplesToInclude;
	private String[] matchDouble, matchString, stratString;
	private double[][] extraIndeps;
	private ExtProjectDataParser parser;
	private int index = 0;
	private Logger log;
	private LS_TYPE lType;

	public CorrectionEvaluator(Project proj, PrincipalComponentsResiduals pcResiduals, EvaluationResult[] precomputeds, int[] order, boolean[][] samplesToInclude, double[][] extraIndeps, LS_TYPE lType) {
		super();
		this.proj = proj;
		this.samplesToInclude = samplesToInclude;
		this.log = proj.getLog();
		this.matchDouble = gatherPatternTitles(proj.SAMPLE_DATA_FILENAME.getValue(), DOUBLE_DATA_PATTERN, log);
		this.matchString = gatherPatternTitles(proj.SAMPLE_DATA_FILENAME.getValue(), STRING_DATA_PATTERN, log);
		this.stratString = gatherPatternTitles(proj.SAMPLE_DATA_FILENAME.getValue(), Array.concatAll(STRAT_STRING_PATTERN, INDEPS_CATS), log);
		loadSampleData();
		this.extraIndeps = extraIndeps;
		this.iterator = precomputeds == null ? new PrincipalComponentsIterator(pcResiduals, order) : null;
		this.precomputeds = precomputeds;
		if (precomputeds != null) {
			log.reportTimeWarning("Using pre-computed estimates to analyze");
		}
		this.lType = lType;

	}

	public static void serialize(CorrectionEvaluator cEvaluator, String fileName) {
		SerializedFiles.writeSerial(cEvaluator, fileName, true);
	}

	public static CorrectionEvaluator readSerial(String fileName, Logger log) {
		return (CorrectionEvaluator) SerializedFiles.readSerial(fileName, false, log, false, true);
	}

	public ExtProjectDataParser getParser() {
		return parser;
	}

	public boolean[][] getSamplesToInclude() {
		return samplesToInclude;
	}

	@Override
	public boolean hasNext() {
		boolean next = false;
		if (precomputeds != null) {
			next = index < precomputeds.length;
		}
		else {
			next = iterator.hasNext();
		}
		return next;
	}

	@Override
	public Callable<EvaluationResult> next() {
		EvaluationWorker worker = new EvaluationWorker(precomputeds == null ? iterator.next() : null, precomputeds == null ? null : precomputeds[index], extraIndeps, matchString, matchDouble, stratString, samplesToInclude, parser, lType, log);
		index++;
		return worker;
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub

	}

	@Override
	public void shutdown() {
		// TODO Auto-generated method stub
	}

	private static class EvaluationWorker implements Callable<EvaluationResult> {
		private PrincipalComponentsResiduals tmpResiduals;
		private double[][] extraIndeps;
		private String[] matchString;
		private String[] matchDouble;
		private String[] stratString;
		private boolean[][] samplesToInclude;
		private EvaluationResult precomputed;
		private LS_TYPE lType;
		private ExtProjectDataParser parser;
		private Logger log;

		public EvaluationWorker(PrincipalComponentsResiduals tmpResiduals, EvaluationResult precomputed, double[][] extraIndeps, String[] matchString, String[] matchDouble, String[] stratString, boolean[][] samplesToInclude, ExtProjectDataParser parser, LS_TYPE lType, Logger log) {
			super();
			this.tmpResiduals = tmpResiduals;
			this.extraIndeps = extraIndeps;
			this.matchString = matchString;
			this.matchDouble = matchDouble;
			this.stratString = Array.concatAll(new String[] { NO_STRAT }, stratString);
			this.samplesToInclude = samplesToInclude;
			this.parser = parser;
			this.lType = lType;
			this.precomputed=precomputed;
			this.log = log;
		}

		@Override
		public EvaluationResult call() throws Exception {
			return evaluate(tmpResiduals, precomputed, extraIndeps, matchString, matchDouble, stratString, samplesToInclude, parser, lType, log);
		}

	}

	private static EvaluationResult evaluate(PrincipalComponentsResiduals tmpResiduals, EvaluationResult precomputed, double[][] extraIndeps, String[] matchString, String[] matchDouble, String[] stratString, boolean[][] samplesToInclude, ExtProjectDataParser parser, LS_TYPE lType, Logger log) {
		String baseTitle = "";
		double[] estimate = null;
		double rsquare = Double.NaN;
		// tmpResiduals.crossValidate(kFolds, numComponentsIter, numThreads, svdRegression, tmpOutput, val_pcs)
		if (precomputed != null) {
			estimate = precomputed.getEstimateData();
			rsquare = precomputed.getrSquared();
			baseTitle = precomputed.getTitle();
		} else {
			baseTitle = "" + tmpResiduals.getNumComponents();
			if (tmpResiduals.getNumComponents() > 0 || extraIndeps != null) {

				CrossValidation cValidation = tmpResiduals.getCorrectedDataAt(tmpResiduals.getMedians(), extraIndeps, samplesToInclude[0], tmpResiduals.getNumComponents(), lType, "HFDS", true);
				estimate = cValidation.getResiduals();
				rsquare = cValidation.getFullModelR2();

			} else {
				estimate = tmpResiduals.getMedians();
			}
			if (estimate.length != tmpResiduals.getProj().getSamples().length) {
				throw new IllegalStateException("Could not obtain estimate for all samples in project");
			}
		}

		EvaluationResult evaluationResult = new EvaluationResult(baseTitle, estimate, rsquare);

		for (int i = 0; i < stratString.length; i++) {
			boolean[][] strat = new boolean[][] { Array.booleanArray(samplesToInclude[1].length, true) };
			String[] stratTitles = new String[] { NO_STRAT };
			if (!stratString[i].equals(NO_STRAT)) {
				BooleanClassifier bClassifier = Array.classifyStringsToBoolean(parser.getStringDataForTitle(stratString[i]), new String[] { "NaN" });
				stratTitles = bClassifier.getTitles();
				strat = bClassifier.getClassified();
			}
			for (int j = 0; j < stratTitles.length; j++) {
				boolean[] finalEval = Array.booleanArray(strat[j].length, false);
				for (int k = 0; k < strat[j].length; k++) {
					finalEval[k] = strat[j][k] && samplesToInclude[1][k];
				}

				for (int k = 0; k < matchString.length; k++) {
					String[] response = Array.subArray(parser.getStringDataForTitle(matchString[k]), finalEval);
					double[] data = Array.subArray(estimate, finalEval);
					ICC icc = new ICC(data, response, EVAL_MASKS[0], null, false, log);
					icc.computeICC();
					evaluationResult.getIccs().add(icc);
					evaluationResult.getNumIndsIcc().add(Array.booleanArraySum(finalEval));
					evaluationResult.getIccTitles().add(matchString[k] + "_" + stratTitles[j]);
					if (Integer.parseInt(baseTitle) % 1 == 0) {
						log.reportTimeInfo("ICC: " + matchString[k] + "_" + stratTitles[j] + " -> " + icc.getICC() + " NumComps = " + baseTitle);
					}
				}
				for (int k = 0; k < matchDouble.length; k++) {
				
					StatPrep result = prepData(estimate, parser.getNumericDataForTitle(matchDouble[k]), finalEval, matchDouble[k], true, log);
					ICC icc = new ICC(result.getFinalData(), result.getFinalResponse(), null, null, false, log);
					icc.computeICC();
					evaluationResult.getIccs().add(icc);
					evaluationResult.getNumIndsIcc().add(Array.booleanArraySum(finalEval));
					evaluationResult.getIccTitles().add(matchDouble[k] + "_" + stratTitles[j]);
					double[][] correlData = new double[][] { result.getInternalEstimate(), result.getExternalEstimate() };
					double[] pearson = Correlation.Pearson(correlData);
					double[] spearman = Correlation.Spearman(correlData);
					evaluationResult.getPearsonCorrels().add(pearson);
					evaluationResult.getNumIndsCorrel().add(result.getInternalEstimate().length);
					evaluationResult.getSpearmanCorrel().add(spearman);
					evaluationResult.getCorrelTitles().add(matchDouble[k] + "_" + stratTitles[j]);
					if (Integer.parseInt(baseTitle) % 1 == 0) {
						log.reportTimeInfo("Spearman: " + matchDouble[k] + "_" + stratTitles[j] + " -> " + Array.toStr(spearman) + " NumComps = " + baseTitle);
					}
				}
			}
		}
		return evaluationResult;
	}

	private static StatPrep prepData(double[] internalEstimate, double[] externalEstimate, boolean[] samplesToInclude, String title, boolean normalize, Logger log) {
		StatPrep result = null;
		if (internalEstimate.length != externalEstimate.length) {
			log.reportTimeError("For " + title + ", internal n=" + internalEstimate.length + " data points do not match external n=" + externalEstimate.length + " datapoints");

		} else {
			ArrayList<Double> tmpInternals = new ArrayList<Double>();
			ArrayList<Double> tmpExternals = new ArrayList<Double>();
			ArrayList<String> tmpResponseInternal = new ArrayList<String>();
			ArrayList<String> tmpResponseExternal = new ArrayList<String>();

			for (int i = 0; i < externalEstimate.length; i++) {
				if (samplesToInclude[i] && !Double.isNaN(internalEstimate[i]) && !Double.isNaN(externalEstimate[i])) {
					tmpResponseInternal.add(i + "");
					tmpInternals.add(internalEstimate[i]);
					tmpResponseExternal.add(i + "");
					tmpExternals.add(externalEstimate[i]);
				}
			}

			double[] internals = Array.toDoubleArray(tmpInternals);
			double[] externals = Array.toDoubleArray(tmpExternals);
			double[] internalNotNorm = internals;
			double[] externalNotNorm = externals;

			double[] finalData = new double[externals.length * 2];
			String[] finalResponse = new String[externals.length * 2];
			int index = 0;
			if (normalize) {
				internals = Array.normalize(internals);
				externals = Array.normalize(externals);
			}
			for (int i = 0; i < externals.length; i++) {
				finalData[index] = internals[i];
				finalResponse[index] = i + "";
				index++;
				finalData[index] = externals[i];
				finalResponse[index] = i + "";
				index++;
			}
			result = new StatPrep(finalData, finalResponse, internalNotNorm, externalNotNorm);
		}
		return result;
	}

	private static class StatPrep {
		private double[] finalData;
		private String[] finalResponse;
		private double[] internalEstimate;
		private double[] externalEstimate;

		public StatPrep(double[] finalData, String[] finalResponse, double[] internalEstimate, double[] externalEstimate) {
			super();
			this.finalData = finalData;
			this.finalResponse = finalResponse;
			this.internalEstimate = internalEstimate;
			this.externalEstimate = externalEstimate;
		}

		public double[] getFinalData() {
			return finalData;
		}

		public String[] getFinalResponse() {
			return finalResponse;
		}

		public double[] getInternalEstimate() {
			return internalEstimate;
		}

		public double[] getExternalEstimate() {
			return externalEstimate;
		}

	}

	private void loadSampleData() {
		log.reportTimeInfo("Found " + matchDouble.length + "(" + Array.toStr(matchDouble) + ") data columns to load matching the patterns defined by " + Array.toStr(DOUBLE_DATA_PATTERN));
		log.reportTimeInfo("Found " + matchString.length + "(" + Array.toStr(matchString) + ") String columns to load matching the patterns defined by " + Array.toStr(STRING_DATA_PATTERN));
		ExtProjectDataParser.ProjectDataParserBuilder builder = new ExtProjectDataParser.ProjectDataParserBuilder();
		builder.sampleBased(true);
		builder.treatAllNumeric(false);
		builder.requireAll(true);
		builder.verbose(true);
		builder.dataKeyColumnName("DNA");
		// System.out.println(Array.toStr(matchDouble));
		// System.exit(1);
		builder.stringDataTitles(Array.concatAll(matchString, stratString));
		builder.numericDataTitles(Array.concatAll(matchDouble));
		try {
			log.reportTimeInfo("Loading " + proj.SAMPLE_DATA_FILENAME.getValue());
			this.parser = builder.build(proj, proj.SAMPLE_DATA_FILENAME.getValue());
			parser.determineIndicesFromTitles();
			parser.loadData();
			log.reportTimeInfo("Finished loading " + proj.SAMPLE_DATA_FILENAME.getValue());

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	private String[] gatherPatternTitles(String dataFile, String[] patterns, Logger log) {

		String[] header = Files.getHeaderOfFile(dataFile, log);
		ArrayList<String> matches = new ArrayList<String>();
		for (int i = 0; i < header.length; i++) {
			for (int j = 0; j < patterns.length; j++) {
				if (header[i].startsWith(patterns[j])) {
					matches.add(header[i]);
				}
			}
		}
		return matches.toArray(new String[matches.size()]);
	}
}
// if (cds != null) {
// icc = new ICC(estimate, cds, EVAL_MASKS[0], null, false, log);
// icc.computeICC();
// log.reportTimeInfo("DFDSFDSFICC: " + matchString[i] + " -> " + icc.getICC() + " NumComps = " + tmpResiduals.getNumComponents());
// log.reportTimeInfo("OTHERMETHOD"+others[ic]);
// }
// String[] samples = Array.subArray(tmpResiduals.getProj().getSamples(), samplesToInclude[1]);
// String[] cds = null;
// for (int j = 0; j < classDefinitions.length; j++) {
// if (classDefinitions[j].getClassTitle().equals("DuplicateNotNA") && matchString[i].equals("EVAL_CLASS_DUPLICATE")) {
// cds =classDefinitions[j].getClassDefs();
// ic=j;
// }
// }
//
// if (tmpResiduals.getNumComponents() == 100) {
// String[] wtf = new String[tmpResiduals.getProj().getSamples().length];
// for (int j = 0; j < tmpResiduals.getProj().getSamples().length; j++) {
// wtf[j] = tmpResiduals.getProj().getSamples()[j] + "\t" + tmpResiduals.getMedians()[j] + "\t" + estimate[j];
// }
// Files.writeList(wtf, tmpResiduals.getProj().PROJECT_DIRECTORY.getValue() + "WTF.resid");
// //System.exit(1);
//
// }
// // for (int j = 0; j < data.length; j++) {
// // ArrayList<String> dupResponse = new ArrayList<String>();
// // ArrayList<Double> dupData = new ArrayList<Double>();
//
// // if (matchString[i].equals("EVAL_CLASS_DUPLICATE") && ext.indexOfStr(response[j], EVAL_MASKS[0]) < 0) {
// // dupResponse.add(response[j]);
// // dupData.add(data[j]);
// // System.out.println(samples[j] + "\t" + matchString[i] + "\t" + data[j] + "\t" + response[j]);
// // if (tmpResiduals.getNumComponents() < 3) {
// // String name = tmpResiduals.getProj().PROJECT_DIRECTORY.getValue() + "ICC/icc." + tmpResiduals.getNumComponents() + ".txt";
// //
// // try {
// // PrintWriter writer = new PrintWriter(new FileWriter(name));
// // for (int k = 0; k < dupResponse.size(); k++) {
// // writer.println(dupResponse.get(k) + "\t" + dupData.get(k));
// // }
// // writer.close();
// // } catch (Exception e) {
// // log.reportError("Error writing to " + name);
// // log.reportException(e);
// // }
// // }
// // }
// // }
//
// ClassDefinition[] classDefinitions = ClassDefinition.getClassDefinitionsFromSampleData(tmpResiduals.getProj());
// double[] others = IntensityCorrectionQC.computeAt(tmpResiduals.getMedians(), svd, tmpResiduals, classDefinitions, samplesToInclude[0], tmpResiduals.getNumComponents(), log);

