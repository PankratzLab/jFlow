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
	private static final String[][] EVAL_MASKS = new String[][] { { "0", "-1", "NaN" }, { "0", "-1" } };
	private static final String[] DOUBLE_DATA = new String[] { "AGE" };
	private static final String[] DOUBLE_DATA_PATTERN = new String[] { "EVAL_DATA", "AGE" };// For Correlation(Spearman by ICC)
	private static final String[] STRING_DATA_PATTERN = new String[] { "EVAL_CLASS", "AGE" };// For ICC
	public static final String[] INDEPS = new String[] { "CLASS=SEX", "AGE" };
	public static final int NUM_PC_SVD_OVERIDE = 160;
	private Project proj;
	private PrincipalComponentsIterator iterator;
	private boolean[][] samplesToInclude;
	private String[] matchDouble, matchString;
	private double[][] extraIndeps;
	private ExtProjectDataParser parser;
	private Logger log;
	private boolean svd;

	public CorrectionEvaluator(Project proj, PrincipalComponentsResiduals pcResiduals, int[] order, boolean[][] samplesToExclude, double[][] extraIndeps, boolean svd) {
		super();
		this.proj = proj;
		this.samplesToInclude = samplesToExclude;
		this.log = proj.getLog();
		this.matchDouble = gatherPatternTitles(proj.SAMPLE_DATA_FILENAME.getValue(), DOUBLE_DATA_PATTERN, log);
		this.matchString = gatherPatternTitles(proj.SAMPLE_DATA_FILENAME.getValue(), STRING_DATA_PATTERN, log);
		loadSampleData();
		this.extraIndeps = extraIndeps;
		this.iterator = new PrincipalComponentsIterator(pcResiduals, order);
		this.svd = svd;

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
		return new EvaluationWorker(iterator.next(), extraIndeps, matchString, matchDouble, samplesToInclude, parser, svd, log);
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
		private boolean[][] samplesToInclude;
		private boolean svd;
		private ExtProjectDataParser parser;
		private Logger log;

		public EvaluationWorker(PrincipalComponentsResiduals tmpResiduals, double[][] extraIndeps, String[] matchString, String[] matchDouble, boolean[][] samplesToInclude, ExtProjectDataParser parser, boolean svd, Logger log) {
			super();
			this.tmpResiduals = tmpResiduals;
			this.extraIndeps = extraIndeps;
			this.matchString = matchString;
			this.matchDouble = matchDouble;
			this.samplesToInclude = samplesToInclude;
			this.parser = parser;
			this.svd = svd;
			this.log = log;
		}

		@Override
		public EvaluationResult call() throws Exception {
			return evaluate(tmpResiduals, extraIndeps, matchString, matchDouble, samplesToInclude, parser, svd, log);
		}

	}

	private static EvaluationResult evaluate(PrincipalComponentsResiduals tmpResiduals, double[][] extraIndeps, String[] matchString, String[] matchDouble, boolean[][] samplesToInclude, ExtProjectDataParser parser, boolean svd, Logger log) {
		String baseTitle = "" + tmpResiduals.getNumComponents();
		double[] estimate = null;
		CrossValidation cValidation = null;
		if (tmpResiduals.getNumComponents() > 0 || extraIndeps != null) {
			
			cValidation = tmpResiduals.getCorrectedDataAt(tmpResiduals.getMedians(), extraIndeps, samplesToInclude[0], tmpResiduals.getNumComponents(), svd, "HFDS", true);
			estimate = cValidation.getResiduals();
		} else {
			estimate = tmpResiduals.getMedians();
		}
		if (estimate.length != tmpResiduals.getProj().getSamples().length) {
			throw new IllegalStateException("Could not obtain estimate for all samples in project");
		}

		EvaluationResult evaluationResult = new EvaluationResult(baseTitle, estimate, cValidation == null ? 0 : cValidation.getRsquare());
		for (int i = 0; i < matchString.length; i++) {
			String[] response = Array.subArray(parser.getStringDataForTitle(matchString[i]), samplesToInclude[1]);
			double[] data = Array.subArray(estimate, samplesToInclude[1]);
			ICC icc = new ICC(data, response, EVAL_MASKS[0], null, false, log);
			icc.computeICC();
			evaluationResult.getIccs().add(icc);
			evaluationResult.getIccTitles().add(matchString[i]);
			log.reportTimeInfo("ICC: " + matchString[i] + " -> " + icc.getICC() + " NumComps = " + tmpResiduals.getNumComponents());
		}
		for (int i = 0; i < matchDouble.length; i++) {
			StatPrep result = prepData(estimate, parser.getNumericDataForTitle(matchDouble[i]), samplesToInclude[1], matchDouble[i], true, log);
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
		ExtProjectDataParser.Builder builder = new ExtProjectDataParser.Builder();
		builder.sampleBased(true);
		builder.treatAllNumeric(false);
		builder.requireAll(true);
		builder.verbose(true);
		builder.dataKeyColumnName("DNA");
		builder.stringDataTitles(Array.concatAll(STRING_DATA, matchString));
		builder.numericDataTitles(Array.concatAll(DOUBLE_DATA, matchDouble, INDEPS));
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

