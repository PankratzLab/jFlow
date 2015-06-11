package cnv.analysis.pca;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;

import stats.Rscript.COLUMNS_MULTIPLOT;
import stats.Rscript.GEOM_POINT_SIZE;
import stats.Rscript.PLOT_DEVICE;
import stats.Rscript.RScatter;
import stats.Rscript.RScatters;
import stats.Rscript.SCATTER_TYPE;
import stats.StatsCrossTabs.STAT_TYPE;
import stats.StatsCrossTabs.StatsCrossTabRank;
import stats.StatsCrossTabs.VALUE_TYPE;
import common.Array;
import common.Files;
import common.Logger;
import common.PSF;
import common.WorkerTrain;
import common.ext;
import cnv.analysis.pca.PCSelector.SelectionResult;
import cnv.filesys.Project;
import cnv.manage.ExtProjectDataParser;
import cnv.manage.TransposeData;
import cnv.qc.LrrSd;

class CorrectionIterator {
	private Project proj;
	private String markesToEvaluate;
	private String samplesToBuildModels;
	private ITERATION_TYPE iType;
	private ORDER_TYPE oType;
	private MODEL_BUILDER_TYPE bType;
	private String outputDir;
	private boolean svd;
	private int numthreads;
	private IterationResult iterationResult;

	public CorrectionIterator(Project proj, String markesToEvaluate, String samplesToBuildModels, ITERATION_TYPE iType, ORDER_TYPE oType, MODEL_BUILDER_TYPE bType, String outputDir, boolean svd, int numthreads) {
		super();
		this.proj = proj;
		this.markesToEvaluate = markesToEvaluate;
		this.samplesToBuildModels = samplesToBuildModels;
		this.iType = iType;
		this.oType = oType;
		this.bType = bType;
		this.outputDir = outputDir;
		this.svd = svd;
		this.numthreads = numthreads;
	}

	public void run() {
		this.iterationResult = run(proj, markesToEvaluate, samplesToBuildModels, iType, oType, bType, outputDir, svd, numthreads);
	}

	public IterationResult getIterationResult() {
		return iterationResult;
	}

	public enum ITERATION_TYPE {
		/**
		 * The evaluation happens with the addition of other independent variables in addition to PCs
		 */
		WITHOUT_INDEPS, /**
		 * other independent variables are not added
		 */
		WITH_INDEPS;

	}

	public enum MODEL_BUILDER_TYPE {
		/**
		 * We load the additional sample (union with not excluded) file for model building
		 */
		WITH_BUILDERS, /**
		 * We build models with everyone, except excluded individuals
		 */
		WITHOUT_BUILDERS;

	}

	public enum ORDER_TYPE {
		/**
		 * PCs are regressed by their natural ordering (PC1,2,3)
		 */
		NATURAL, /**
		 * PCS are ranked by the amount of variance explained in the estimate of interest
		 */
		RANK_R2,

		/**
		 * PCS are ranked by spearman abs r
		 */
		// RANK_R,
		/**
		 * PCs are filtered for an association with known QC metrics from {@link LrrSd}
		 */
		QC_ASSOCIATION;

	}

	private IterationResult run(Project proj, String markesToEvaluate, String samplesToBuildModels, ITERATION_TYPE iType, ORDER_TYPE oType, MODEL_BUILDER_TYPE bType, String outputDir, boolean svd, int numthreads) {
		Logger log = proj.getLog();

		String output = outputDir + "correctionEval_" + iType + "_" + oType + "_" + bType;
		IterationResult iterationResult = new IterationResult(output, iType, oType, bType);
		if (!Files.exists(iterationResult.getOutputSer()) || oType == ORDER_TYPE.QC_ASSOCIATION) {

			//
			// iterationResult.plotRank(log);
			// iterationResult.plotSummary(new String[] { "Rsquare_correction", "ICC_EVAL_CLASS_DUPLICATE_ALL", "PEARSON_CORREL_AGE", "PEARSON_CORREL_EVAL_DATA_SEX","PEARSON_CORREL_EVAL_DATA_resid.mtDNaN.qPCR.MT001","PEARSON_CORREL_EVAL_DATA_resid.mtDNA.qPCR" }, log);
			// // rScatter.setxLabel("Principal Component ("+oType+")");
			//
			// // public RScatter(String dataFile, String rSriptFile, String output, String dataXvalueColumn, String[] dataYvalueColumns, SCATTER_TYPE sType, Logger log) {
			//
			// System.exit(1);

			proj.getLog().reportTimeInfo("Loading " + proj.INTENSITY_PC_FILENAME.getValue());

			new File(outputDir).mkdirs();
			log.reportTimeInfo("Beginning iteration evaluation:");
			log.reportTimeInfo("PC file: " + proj.INTENSITY_PC_FILENAME.getValue());
			log.reportTimeInfo("Iteration type : " + iType);
			log.reportTimeInfo("Order type : " + oType);
			log.reportTimeInfo("Model building type: " + bType);

			PrincipalComponentsResiduals pcResiduals = proj.loadPcResids();
			pcResiduals.setMarkersToAssessFile(markesToEvaluate);
			pcResiduals.setHomozygousOnly(true);
			pcResiduals.computeAssessmentDataMedians();
			boolean[] samplesForModels = null;
			boolean valid = true;

			switch (bType) {
			case WITHOUT_BUILDERS:
				samplesForModels = proj.getSamplesToInclude(null, true, true);
				break;
			case WITH_BUILDERS:
				if (!Files.exists(samplesToBuildModels)) {
					log.reportTimeError("Model building type was set to " + bType + " but the sample file " + samplesToBuildModels + " did not exist");
					valid = false;
				} else {
					samplesForModels = proj.getSamplesToInclude(samplesToBuildModels, true, true);
				}
				break;
			default:
				break;

			}

			CorrectionEvaluator cEvaluator = new CorrectionEvaluator(proj, pcResiduals, null, null, null, svd);
			int[] order = null;
			double[][] extraIndeps = null;
			StatsCrossTabRank sTabRank = null;
			switch (iType) {
			case WITHOUT_INDEPS:
				log.reportTimeInfo("Evaluating with " + Array.booleanArraySum(samplesForModels) + " samples, no additional independent variables");
				break;
			case WITH_INDEPS:
				extraIndeps = loadIndeps(cEvaluator, CorrectionEvaluator.INDEPS, new double[][] { { 0, 3, 4, 5, 6, 7, 8, 9, 10, Double.NaN }, { -1, Double.NaN } }, log);
				if (extraIndeps == null) {
					log.reportTimeError("type = " + iType + " and were missing some of the following " + Array.toStr(CorrectionEvaluator.INDEPS));
					log.reportTimeError("Available = " + Array.toStr(cEvaluator.getParser().getNumericDataTitles()));
					valid = false;
				} else {
					boolean[] tmpInclude = new boolean[samplesForModels.length];
					Arrays.fill(tmpInclude, false);
					for (int i = 0; i < extraIndeps.length; i++) {

						if (samplesForModels[i]) {
							boolean hasNan = false;
							for (int j = 0; j < extraIndeps[i].length; j++) {
								if (Double.isNaN(extraIndeps[i][j])) {
									hasNan = true;
								}
							}
							if (!hasNan) {
								tmpInclude[i] = true;
							}
						}
					}
					log.reportTimeInfo("Original number of samples: " + Array.booleanArraySum(samplesForModels));
					log.reportTimeInfo("Number of samples with valid independant variables, final evaluation set: " + Array.booleanArraySum(tmpInclude));
					samplesForModels = tmpInclude;
				}
				break;
			default:
				break;

			}
			iterationResult.setValid(valid);
			if (valid) {
				sTabRank = pcResiduals.getStatRankFor(pcResiduals.getMedians(), extraIndeps, samplesForModels, "RAW_MEDIANS", STAT_TYPE.LIN_REGRESSION, VALUE_TYPE.STAT, proj.getLog());
				sTabRank.dump(iterationResult.getOutputRank(), oType != ORDER_TYPE.NATURAL, log);

				switch (oType) {
				case NATURAL:
					order = null;
					break;
				case RANK_R2:
					order = new int[sTabRank.getOrder().length];
					for (int i = 0; i < sTabRank.getOrder().length; i++) {
						order[i] = sTabRank.getOrder()[i] + 1;// one based for pcs
					}
					break;
				case QC_ASSOCIATION:

					SelectionResult result = PCSelector.select(proj, 0.10, STAT_TYPE.SPEARMAN_CORREL);
					order = result.getOrder();
					if (result == null || order.length < 1) {
						log.reportTimeError("Could not select PCs from QC metrics, trying again");
						Files.copyFile(proj.SAMPLE_QC_FILENAME.getValue(), proj.SAMPLE_QC_FILENAME.getValue() + ext.getTimestampForFilename());
						new File(proj.SAMPLE_QC_FILENAME.getValue()).delete();
						LrrSd.init(proj, null, null, numthreads);
						result = PCSelector.select(proj, 0.10, STAT_TYPE.SPEARMAN_CORREL);
						order = result.getOrder();
						if (order.length < 1) {
							return null;
						}
					}
				default:
					break;
				}
				if (!Files.exists(iterationResult.getOutputSer())) {
					boolean[] samplesToEvaluate = proj.getSamplesToInclude(null);
					log.reportTimeInfo(Array.booleanArraySum(samplesForModels) + " samples for models");
					log.reportTimeInfo(Array.booleanArraySum(samplesToEvaluate) + " samples for evaluation");

					cEvaluator = new CorrectionEvaluator(proj, pcResiduals, order, new boolean[][] { samplesForModels, samplesToEvaluate }, extraIndeps, svd);
					ArrayList<EvaluationResult> store = new ArrayList<EvaluationResult>();

					try {
						PrintWriter writer = new PrintWriter(new FileWriter(iterationResult.getOutputSummary()));
						WorkerTrain<EvaluationResult> train = new WorkerTrain<EvaluationResult>(cEvaluator, numthreads, numthreads, proj.getLog());
						int index = 0;
						while (train.hasNext()) {
							EvaluationResult result = train.next();
							result.setItType(iType);
							result.setOrType(oType);
							result.setbType(bType);
							if (index == 0) {
								writer.println(Array.toStr(result.getHeader()));
							}
							writer.println(Array.toStr(result.getData()));
							index++;
							result.shrink();
							store.add(result);
						}
						writer.close();
					} catch (Exception e) {
						proj.getLog().reportError("Error writing to " + iterationResult.getOutputSummary());
						proj.getLog().reportException(e);
					}
					EvaluationResult.serialize(store.toArray(new EvaluationResult[store.size()]), iterationResult.getOutputSer());
				}
			}
		}
		return iterationResult;
	}

	private static class IterationResult {
		private String outputRoot;
		private String outputSer;
		private String outputRank;
		private String outputSummary;
		private String rankRscript;
		private String evalRscript;
		private ITERATION_TYPE iType;
		private MODEL_BUILDER_TYPE bType;
		private ORDER_TYPE oType;
		private String rankplot;
		private String evalPlot;
		private boolean valid;

		public IterationResult(String outputRoot, ITERATION_TYPE iType, ORDER_TYPE oType, MODEL_BUILDER_TYPE bType) {
			super();
			this.outputRoot = outputRoot;
			this.outputSer = outputRoot + ".summary.ser";
			this.outputSummary = outputRoot + ".summary.txt";
			this.outputRank = outputRoot + ".rank.txt";
			this.rankRscript = outputRoot + ".rank.rscript";
			this.rankplot = outputRoot + ".rank.Plot.pdf";
			this.evalRscript = outputRoot + ".eval.rscript";
			this.evalPlot = outputRoot + ".eval.Plot.pdf";
			this.iType = iType;
			this.oType = oType;
			this.bType = bType;
		}

		public RScatter plotRank(Logger log) {
			RScatter rScatter = new RScatter(outputRank, rankRscript, ext.rootOf(outputRank), rankplot, "OriginalOrder", new String[] { "Stat" }, SCATTER_TYPE.POINT, log);
			rScatter.setyRange(new double[] { 0, 1 });
			rScatter.setxLabel("PC (" + oType + " - sorted)");
			rScatter.setyLabel("Rsq");
			rScatter.setTitle(iType + " " + bType);
			rScatter.setgPoint_SIZE(GEOM_POINT_SIZE.GEOM_POINT);
			rScatter.execute();
			return rScatter;
		}

		public RScatter plotSummary(String[] dataColumns, Logger log) {
			RScatter rScatter = new RScatter(outputSummary, evalRscript, ext.rootOf(outputSummary), evalPlot, "Evaluated", dataColumns, SCATTER_TYPE.POINT, log);
			rScatter.setyRange(new double[] { -1, 1 });
			rScatter.setxLabel("PC (" + oType + " - sorted)");
			rScatter.setTitle(iType + " " + bType);
			rScatter.setgPoint_SIZE(GEOM_POINT_SIZE.GEOM_POINT);
			rScatter.execute();
			rScatter.setOutput(ext.addToRoot(evalPlot, ".trim"));
			rScatter.setxRange(new double[] { 0, 50 });
			rScatter.execute();
			rScatter.setOutput(evalPlot);
			rScatter.setxRange(null);
			return rScatter;

		}

		public void setValid(boolean valid) {
			this.valid = valid;
		}

		public String getOutputSer() {
			return outputSer;
		}

		public String getOutputRank() {
			return outputRank;
		}

		public String getOutputSummary() {
			return outputSummary;
		}

	}

	private static double[][] loadIndeps(CorrectionEvaluator cEvaluator, String[] indepHeaders, double[][] indepMasks, Logger log) {
		ExtProjectDataParser parser = cEvaluator.getParser();

		double[][] extraIndeps = new double[Array.booleanArraySum(parser.getDataPresent())][indepHeaders.length];
		for (int i = 0; i < indepHeaders.length; i++) {
			int curSum = 0;
			if (parser.getNumericDataForTitle(indepHeaders[i]).length <= 0) {
				log.reportTimeError("Did not find " + indepHeaders[i] + " to load for independent variable");
				return null;
			} else {
				double[] data = parser.getNumericDataForTitle(indepHeaders[i]);
				for (int j = 0; j < data.length; j++) {
					boolean mask = false;
					for (int k = 0; k < indepMasks[i].length; k++) {
						if (data[j] == indepMasks[i][k]) {
							mask = true;
						}
					}
					if (mask) {
						extraIndeps[j][i] = Double.NaN;
					} else {
						extraIndeps[j][i] = data[j];
						curSum++;
					}
				}
			}
			log.reportTimeInfo(curSum + " samples for independant variable " + indepHeaders[i]);
		}
		return extraIndeps;
	}

	private static CorrectionIterator[] getIterations(Project proj, String markesToEvaluate, String samplesToBuildModels, String outputDir, boolean svd, int numthreads) {
		ArrayList<CorrectionIterator> cIterators = new ArrayList<CorrectionIterator>();
		System.out.println("JDOFJSDF remember the pcs");
		for (int i = 0; i < ITERATION_TYPE.values().length; i++) {
			for (int j = 0; j < ORDER_TYPE.values().length; j++) {
				for (int j2 = 0; j2 < MODEL_BUILDER_TYPE.values().length; j2++) {
					cIterators.add(new CorrectionIterator(proj, markesToEvaluate, samplesToBuildModels, ITERATION_TYPE.values()[i], ORDER_TYPE.values()[j], MODEL_BUILDER_TYPE.values()[j2], outputDir, svd, numthreads));
				}
			}
		}
		return cIterators.toArray(new CorrectionIterator[cIterators.size()]);
	}

	public static CorrectionIterator[] runAll(Project proj, String markesToEvaluate, String samplesToBuildModels, String outputDir, boolean svd, int numthreads) {
		proj.INTENSITY_PC_NUM_COMPONENTS.setValue(400);
		System.out.println("JDOFJSDF remember the pcs");
		if (outputDir == null) {
			outputDir = proj.PROJECT_DIRECTORY.getValue() + ext.rootOf(proj.INTENSITY_PC_FILENAME.getValue()) + "_eval/";
		}
		if (!Files.exists(proj.MARKER_DATA_DIRECTORY.getValue() + "markers.0.mdRAF")) {
			proj.getLog().reportTimeWarning("Did not see " + proj.MARKER_DATA_DIRECTORY.getValue() + "markers.0.mdRAF, attempting to transpose now");
			TransposeData.transposeData(proj, 2000000000, false);
		}
		CorrectionIterator[] cIterators = getIterations(proj, markesToEvaluate, samplesToBuildModels, outputDir, svd, numthreads);
		ArrayList<RScatter> rScatters = new ArrayList<RScatter>();

		for (int i = 0; i < cIterators.length; i++) {
			cIterators[i].run();
			IterationResult iterationResult = cIterators[i].getIterationResult();
			iterationResult.plotRank(proj.getLog());
			RScatter rScatter = iterationResult.plotSummary(new String[] { "Rsquare_correction", "ICC_EVAL_CLASS_DUPLICATE_ALL", "ICC_EVAL_CLASS_DUPLICATE_SAME_VISIT", "ICC_EVAL_CLASS_FC", "SPEARMAN_CORREL_AGE", "SPEARMAN_CORREL_EVAL_DATA_SEX", "SPEARMAN_CORREL_EVAL_DATA_resid.mtDNaN.qPCR.MT001", "SPEARMAN_CORREL_EVAL_DATA_resid.mtDNA.qPCR", "SPEARMAN_CORREL_EVAL_DATA_Mt_DNA_relative_copy_number" }, proj.getLog());
			rScatters.add(rScatter);
			rScatters.add(iterationResult.plotRank(proj.getLog()));
		}
		String outputRoot = outputDir + "finalSummary";
		RScatters finalScatters = new RScatters(rScatters.toArray(new RScatter[rScatters.size()]), outputRoot + ".rscript", outputRoot + ".pdf", COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_2, PLOT_DEVICE.PDF, proj.getLog());
		finalScatters.execute();
		return cIterators;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String proj = null;
		String markers = null;
		String defaultDir = null;
		boolean svd = false;
		int numThreads = 3;
		String samplesToBuildModels = null;
		String usage = "\n" + "cnv.analysis.pca.CorrectionEvaluator requires 0-1 arguments\n";
		usage += "   (1) project filename (i.e. proj=" + proj + " (default))\n" + "";
		usage += "   (2) markers to Evaluate (i.e. markers=" + markers + " (default))\n" + "";
		usage += PSF.Ext.getOutputDirCommand(3, defaultDir);
		usage += PSF.Ext.getNumThreadsCommand(4, numThreads);
		usage += "   (5) svd regression (i.e.-svd (not default))\n" + "";
		usage += "   (6) samples to generate models (i.e.samples= (no default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				proj = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("markers=")) {
				markers = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("samples=")) {
				samplesToBuildModels = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-svd")) {
				svd = true;
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.OUTPUT_DIR_COMMAND)) {
				defaultDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			runAll(new Project(proj, false), markers, samplesToBuildModels, defaultDir, svd, numThreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
