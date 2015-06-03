package cnv.analysis.pca;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import stats.Correlation;
import stats.CrossValidation;
import stats.ICC;

import common.Array;
import common.Files;
import common.Logger;
import common.WorkerTrain.Producer;
import cnv.analysis.pca.PrincipalComponentsResiduals.PrincipalComponentsIterator;
import cnv.filesys.Project;
import cnv.manage.ExtProjectDataParser;

class CorrectionEvaluator implements Producer<EvaluationResult> {
	private static final String[] STRING_DATA = new String[] {};// CLASS Sex
	private static final String[][] EVAL_MASKS = new String[][] { { "0" }, { "0", "-1" } };
	private static final String[] DOUBLE_DATA = new String[] { "AGE" };
	private static final String[] DOUBLE_DATA_PATTERN = new String[] { "EVAL_DATA", "AGE" };// For Correlation(Spearman by ICC)
	private static final String[] STRING_DATA_PATTERN = new String[] { "EVAL_CLASS", "AGE" };// For ICC
	public static final String[] INDEPS = new String[] { "CLASS=SEX", "AGE" };

	// private static final String GEN_TAG = "GEN_ESTIMATE";
	private Project proj;
	private PrincipalComponentsIterator iterator;
	private boolean[] samplesToInclude;
	private String[] matchDouble, matchString;
	private double[][] extraIndeps;
	private ExtProjectDataParser parser;
	private Logger log;

	public CorrectionEvaluator(Project proj, PrincipalComponentsResiduals pcResiduals, int[] order, boolean[] samplesToExclude, double[][] extraIndeps) {
		super();
		this.proj = proj;
		this.samplesToInclude = samplesToExclude;
		this.log = proj.getLog();
		this.matchDouble = gatherPatternTitles(proj.SAMPLE_DATA_FILENAME.getValue(), DOUBLE_DATA_PATTERN, log);
		this.matchString = gatherPatternTitles(proj.SAMPLE_DATA_FILENAME.getValue(), STRING_DATA_PATTERN, log);
		loadSampleData();
		this.extraIndeps = extraIndeps;
		this.iterator = new PrincipalComponentsIterator(pcResiduals, order);

	}

	public ExtProjectDataParser getParser() {
		return parser;
	}

	@Override
	public boolean hasNext() {
		return iterator.hasNext();
	}

	@Override
	public Callable<EvaluationResult> next() {
		return new EvaluationWorker(iterator.next(), extraIndeps, matchString, matchDouble, samplesToInclude, parser, log);
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
		private boolean[] samplesToInclude;
		private ExtProjectDataParser parser;
		private Logger log;

		public EvaluationWorker(PrincipalComponentsResiduals tmpResiduals, double[][] extraIndeps, String[] matchString, String[] matchDouble, boolean[] samplesToInclude, ExtProjectDataParser parser, Logger log) {
			super();
			this.tmpResiduals = tmpResiduals;
			this.extraIndeps = extraIndeps;
			this.matchString = matchString;
			this.matchDouble = matchDouble;
			this.samplesToInclude = samplesToInclude;
			this.parser = parser;
			this.log = log;
		}

		@Override
		public EvaluationResult call() throws Exception {
			return evaluate(tmpResiduals, extraIndeps, matchString, matchDouble, samplesToInclude, parser, log);
		}

	}

	private static EvaluationResult evaluate(PrincipalComponentsResiduals tmpResiduals, double[][] extraIndeps, String[] matchString, String[] matchDouble, boolean[] samplesToInclude, ExtProjectDataParser parser, Logger log) {
		String baseTitle = "" + tmpResiduals.getNumComponents();
		EvaluationResult evaluationResult = new EvaluationResult(baseTitle);
		double[] estimate = null;

		if (tmpResiduals.getNumComponents() > 0) {
			CrossValidation cValidation = tmpResiduals.getCorrectedDataAt(tmpResiduals.getMedians(), extraIndeps, samplesToInclude, tmpResiduals.getNumComponents(), false, "HFDS", true);
			estimate = cValidation.getResiduals();
		} else {
			estimate = tmpResiduals.getMedians();
		}
		for (int i = 0; i < matchString.length; i++) {
			String[] response = parser.getStringDataForTitle(matchString[i]);
			ICC icc = new ICC(estimate, response, EVAL_MASKS[0], null, false, log);
			icc.computeICC();
			evaluationResult.getIccs().add(icc);
			evaluationResult.getIccTitles().add(matchString[i]);
		}
		for (int i = 0; i < matchDouble.length; i++) {
			StatPrep result = prepData(estimate, parser.getNumericDataForTitle(matchDouble[i]), matchDouble[i], true, log);
			ICC icc = new ICC(result.getFinalData(), result.getFinalResponse(), null, null, false, log);
			icc.computeICC();
			evaluationResult.getIccs().add(icc);
			evaluationResult.getIccTitles().add(matchDouble[i]);
			double[][] correlData = new double[][] { result.getInternalEstimate(), result.getExternalEstimate() };
			double[] pearson = Correlation.Pearson(correlData);
			double[] spearman = Correlation.Spearman(correlData);
			evaluationResult.getPearsonCorrels().add(pearson);
			evaluationResult.getSpearmanCorrel().add(spearman);
			evaluationResult.getCorrelTitles().add(matchDouble[i]);
		}
		return evaluationResult;
	}

	private static StatPrep prepData(double[] internalEstimate, double[] externalEstimate, String title, boolean normalize, Logger log) {
		StatPrep result = null;
		if (internalEstimate.length != externalEstimate.length) {
			log.reportTimeError("For " + title + ", internal n=" + internalEstimate.length + " data points do not match external n=" + externalEstimate.length + " datapoints");

		} else {
			ArrayList<Double> tmpInternals = new ArrayList<Double>();
			ArrayList<Double> tmpExternals = new ArrayList<Double>();
			ArrayList<String> tmpResponseInternal = new ArrayList<String>();
			ArrayList<String> tmpResponseExternal = new ArrayList<String>();

			for (int i = 0; i < externalEstimate.length; i++) {
				if (!Double.isNaN(internalEstimate[i]) && !Double.isNaN(externalEstimate[i])) {
					tmpResponseInternal.add(i + "");
					tmpInternals.add(internalEstimate[i]);
					tmpResponseExternal.add(i + "");
					tmpExternals.add(externalEstimate[i]);
				}
			}

			double[] internals = Array.toDoubleArray(tmpInternals);
			double[] externals = Array.toDoubleArray(tmpExternals);
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
			result = new StatPrep(finalData, finalResponse, internals, externals);
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
		ExtProjectDataParser.Builder builder = new ExtProjectDataParser.Builder();
		builder.sampleBased(true);
		builder.treatAllNumeric(false);
		builder.requireAll(true);
		builder.verbose(true);
		builder.dataKeyColumnName("DNA");
		builder.stringDataTitles(Array.concatAll(STRING_DATA, matchString));
		builder.numericDataTitles(Array.concatAll(DOUBLE_DATA, matchDouble));
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

	// public static void run(Project proj, String markesToEvaluate, String outputDir) {
	// proj.getLog().reportTimeInfo("Loading " + proj.INTENSITY_PC_FILENAME.getValue());
	// if (outputDir == null) {
	// outputDir = proj.PROJECT_DIRECTORY.getValue() + "mitoEval/";
	//
	// }
	// new File(outputDir).mkdirs();
	//
	// PrincipalComponentsResiduals pcResiduals = proj.loadPcResids();
	// pcResiduals.setMarkersToAssessFile(markesToEvaluate);
	// pcResiduals.computeAssessmentDataMedians();
	// StatsCrossTabRank sTabRank = pcResiduals.getStatRankFor(pcResiduals.getMedians(), null, "RAW_MEDIANS", STAT_TYPE.LIN_REGRESSION, VALUE_TYPE.STAT, proj.getLog());
	// int[] pcs = new int[sTabRank.getOrder().length];
	//
	// for (int i = 0; i < sTabRank.getOrder().length; i++) {
	// pcs[i] = sTabRank.getOrder()[i] + 1;
	// }
	// proj.getLog().reportTimeInfo("Finished loading " + proj.INTENSITY_PC_FILENAME.getValue());
	// boolean[] samplesToInclude = proj.getSamplesToInclude(null);
	//
	// CorrectionEvaluator cEvaluator = new CorrectionEvaluator(proj, pcResiduals, pcs, samplesToInclude, null);
	// String output = outputDir + "mitos.eval.txt";
	// // sTabRank.dump(ext.addToRoot(output, ".rank"), proj.getLog());
	// try {
	// PrintWriter writer = new PrintWriter(new FileWriter(output));
	// WorkerTrain<EvaluationResult> train = new WorkerTrain<EvaluationResult>(cEvaluator, 8, 3, proj.getLog());
	// int index = 0;
	// while (train.hasNext()) {
	// EvaluationResult result = train.next();
	// if (index == 0) {
	// writer.println(Array.toStr(result.getHeader()));
	// }
	// writer.println(Array.toStr(result.getData()));
	// index++;
	// }
	// index++;
	//
	// writer.close();
	// } catch (Exception e) {
	// proj.getLog().reportError("Error writing to " + output);
	// proj.getLog().reportException(e);
	// }
	//
	// // cEvaluator.evaluate();
	// }
	//
	// public static void main(String[] args) {
	// int numArgs = args.length;
	// String proj = null;
	// String markers = null;
	// String defaultDir = null;
	//
	// String usage = "\n" + "cnv.analysis.pca.CorrectionEvaluator requires 0-1 arguments\n";
	// usage += "   (1) project filename (i.e. proj=" + proj + " (default))\n" + "";
	// usage += "   (2) markers to Evaluate (i.e. markers=" + markers + " (default))\n" + "";
	// usage += PSF.Ext.getOutputDirCommand(3, defaultDir);
	// for (int i = 0; i < args.length; i++) {
	// if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
	// System.err.println(usage);
	// System.exit(1);
	// } else if (args[i].startsWith("proj=")) {
	// proj = args[i].split("=")[1];
	// numArgs--;
	// } else if (args[i].startsWith("markers=")) {
	// markers = args[i].split("=")[1];
	// numArgs--;
	// } else if (args[i].startsWith(PSF.Ext.OUTPUT_DIR_COMMAND)) {
	// defaultDir = args[i].split("=")[1];
	// numArgs--;
	// } else if (args[i].startsWith("log=")) {
	// logfile = args[i].split("=")[1];
	// numArgs--;
	// } else {
	// System.err.println("Error - invalid argument: " + args[i]);
	// }
	// }
	// if (numArgs != 0) {
	// System.err.println(usage);
	// System.exit(1);
	// }
	// try {
	// run(new Project(proj, false), markers, defaultDir);
	// } catch (Exception e) {
	// e.printStackTrace();
	// }
	// }

}
