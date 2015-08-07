package cnv.analysis.pca;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.Callable;

import stats.CategoricalPredictor;
import stats.CategoricalPredictor.DummyCoding;
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
import common.WorkerTrain.Producer;
import common.ext;
import cnv.analysis.pca.PCSelector.SELECTION_TYPE;
import cnv.analysis.pca.PCSelector.SelectionResult;
import cnv.filesys.Project;
import cnv.manage.ExtProjectDataParser;
import cnv.manage.TransposeData;
import cnv.qc.LrrSd;

class CorrectionIterator implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
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
		 * PCS are ranked by the amount of variance explained within a stepwise regression
		 */
		STEPWISE_RANK_R2,

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
		//QC_ASSOCIATION;

	}

	private IterationResult run(Project proj, String markesToEvaluate, String samplesToBuildModels, ITERATION_TYPE iType, ORDER_TYPE oType, MODEL_BUILDER_TYPE bType, String outputDir, boolean svd, int numthreads) {
		Logger log = proj.getLog();
		CorrectionEvaluator cEvaluator = null;
		String output = outputDir + "correctionEval_" + iType + "_" + oType + "_" + bType;
		IterationResult iterationResult = new IterationResult(output, iType, oType, bType);
		if (!Files.exists(iterationResult.getOutputSer())) {

			// || !Files.exists(iterationResult.getBasePrep())
			// || oType == ORDER_TYPE.QC_ASSOCIATION
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
			pcResiduals.fillInMissing();
			pcResiduals.setMarkersToAssessFile(markesToEvaluate);

			pcResiduals.setHomozygousOnly(true);
			pcResiduals.computeAssessmentDataMedians();
			boolean[] samplesForModels = null;
			boolean valid = true;
			if (!Files.exists(iterationResult.getBasePrep())) {
				 //|| oType == ORDER_TYPE.QC_ASSOCIATION
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

				cEvaluator = new CorrectionEvaluator(proj, pcResiduals, null, null, null, svd);
				int[] order = null;
				double[][] extraIndeps = null;
				StatsCrossTabRank sTabRank = null;
				switch (iType) {
				case WITHOUT_INDEPS:
					log.reportTimeInfo("Evaluating with " + Array.booleanArraySum(samplesForModels) + " samples, no additional independent variables");
					break;
				case WITH_INDEPS:
					extraIndeps = loadIndeps(cEvaluator, CorrectionEvaluator.INDEPS, new double[][] { { 0, 3, 4, 5, 6, 7, 8, 9, 10, Double.NaN }, { -1, Double.NaN } }, CorrectionEvaluator.INDEPS_CATS, new String[] { "NaN" }, log);

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
					sTabRank = pcResiduals.getStatRankFor(pcResiduals.getMedians(), extraIndeps, samplesForModels, "RAW_MEDIANS", STAT_TYPE.LIN_REGRESSION, VALUE_TYPE.STAT, oType == ORDER_TYPE.STEPWISE_RANK_R2, numthreads, proj.getLog());
					sTabRank.dump(iterationResult.getOutputRank(), oType != ORDER_TYPE.NATURAL&&oType !=ORDER_TYPE.STEPWISE_RANK_R2, log);
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
					case STEPWISE_RANK_R2:
						order = new int[sTabRank.getOrder().length];
						for (int i = 0; i < sTabRank.getOrder().length; i++) {
							order[i] = sTabRank.getOrder()[i] + 1;// one based for pcs
						}
						// log.reportTimeInfo("PC steArray.toStr(order));
						break;
//					case QC_ASSOCIATION:
//
//						SelectionResult result = PCSelector.select(proj, 0.05, STAT_TYPE.SPEARMAN_CORREL, SELECTION_TYPE.EFFECTIVE_M_CORRECTED);
//						order = result.getOrder();
//						if (result == null || order.length < 1) {
//							log.reportTimeError("Could not select PCs from QC metrics, trying again");
//							Files.copyFile(proj.SAMPLE_QC_FILENAME.getValue(), proj.SAMPLE_QC_FILENAME.getValue() + ext.getTimestampForFilename());
//							new File(proj.SAMPLE_QC_FILENAME.getValue()).delete();
//							LrrSd.init(proj, null, null, numthreads);
//							result = PCSelector.select(proj, 0.05, STAT_TYPE.SPEARMAN_CORREL, SELECTION_TYPE.EFFECTIVE_M_CORRECTED);
//							order = result.getOrder();
//							if (order.length < 1) {
//								return null;
//							}
//						}
					default:
						break;
					}
					boolean[] samplesToEvaluate = proj.getSamplesToInclude(null);
					cEvaluator = new CorrectionEvaluator(proj, pcResiduals, order, new boolean[][] { samplesForModels, samplesToEvaluate }, extraIndeps, svd);

					BasicPrep basicPrep = new BasicPrep(cEvaluator.getParser().getNumericData(), cEvaluator.getParser().getNumericDataTitles(), samplesToEvaluate, samplesForModels);
					BasicPrep.serialize(basicPrep, iterationResult.getBasePrep());
					iterationResult.setBasicPrep(basicPrep);
					if (!Files.exists(iterationResult.getOutputSer())) {
						log.reportTimeInfo(Array.booleanArraySum(samplesForModels) + " samples for models");
						log.reportTimeInfo(Array.booleanArraySum(samplesToEvaluate) + " samples for evaluation");

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
		} else {
			BasicPrep basicPrep = BasicPrep.readSerial(iterationResult.getBasePrep(), log);
			iterationResult.setBasicPrep(basicPrep);
		}
		return iterationResult;
	}

	private static class BasicPrep implements Serializable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		private double[][] numericData;
		private boolean[] samplesForModels;
		private boolean[] samplesToEvaluate;
		private String[] numericTitles;

		public BasicPrep(double[][] numericData, String[] numericTitles, boolean[] samplesToEvaluate, boolean[] samplesForModels) {
			super();
			this.numericData = numericData;
			this.samplesToEvaluate = samplesToEvaluate;
			this.samplesForModels = samplesForModels;
			this.numericTitles = numericTitles;
		}

		public boolean[] getSamplesForModels() {
			return samplesForModels;
		}

		public boolean[] getSamplesToEvaluate() {
			return samplesToEvaluate;
		}

		public double[][] getNumericData() {
			return numericData;
		}

		public String[] getNumericTitles() {
			return numericTitles;
		}

		public static void serialize(BasicPrep basicPrep, String fileName) {
			Files.writeSerial(basicPrep, fileName, true);
		}

		public static BasicPrep readSerial(String fileName, Logger log) {
			return (BasicPrep) Files.readSerial(fileName, false, log, false, true);
		}

	}

	private static class IterationResult implements Serializable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		private static final String[] MERLIN_ADDITIONS = new String[] { "NUM_PC", "PERCENT_HERITABLITY" };
		private String outputRoot;
		private String outputSer;
		private String outputRank;
		private String outputSummary;
		private String rankRscript;
		private String evalRscript;
		private String heritRscript;
		private String heritSummary;
		private BasicPrep basicPrep;
		private ITERATION_TYPE iType;
		private MODEL_BUILDER_TYPE bType;
		private ORDER_TYPE oType;
		private String rankplot;
		private String evalPlot;
		private String heritPlot;
		private String basePrep;
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
			this.basePrep = outputRoot + "basePrep.ser";
			this.iType = iType;
			this.oType = oType;
			this.bType = bType;
		}

		public String getBasePrep() {
			return basePrep;
		}

		public void setBasicPrep(BasicPrep basicPrep) {
			this.basicPrep = basicPrep;
		}

		public BasicPrep getBasicPrep() {
			return basicPrep;
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

		public RScatter plotSummary(String[] dataColumns,int index, Logger log) {
			RScatter rScatter = new RScatter(outputSummary, evalRscript, ext.rootOf(outputSummary)+"_"+index, evalPlot, "Evaluated", dataColumns, SCATTER_TYPE.POINT, log);
			rScatter.setyRange(new double[] { -1, 1 });
			rScatter.setxLabel("PC (" + oType + " - sorted)");
			rScatter.setTitle(iType + " " + bType);
			// rScatter.setgPoint_SIZE(GEOM_POINT_SIZE.GEOM_POINT);
			rScatter.execute();
			rScatter.setOutput(ext.addToRoot(evalPlot, ".trim"));
			rScatter.setxRange(new double[] { 0, 50 });
			rScatter.execute();
			rScatter.setOutput(evalPlot);
			rScatter.setxRange(null);
			return rScatter;
		}

		public RScatter plotHeritability(Project proj, String pedFile, boolean[] samplesToEvaluate, Logger log) {
			this.heritRscript = outputRoot + ".summary.heritability.rscript";
			this.heritPlot = outputRoot + ".summary.heritability.pdf";
			this.heritSummary = outputRoot + ".summary.heritability_summary.parsed.xln";
			String tmpHerit = outputRoot + ".summary.heritability_summary.xln";

			if (pedFile != null) {

				if (!Files.exists(heritSummary) || Files.countLines(heritSummary, false) != Files.countLines(outputSummary, false)) {
					if (!Files.exists(tmpHerit) || Files.countLines(tmpHerit, false) != Files.countLines(outputSummary, false)) {
						// System.out.println(Files.exists(tmpHerit) + "\t" + tmpHerit);
						// System.out.println(outputSer);
						// System.out.println(outputRoot);

						EvaluationResult.prepareHeritability(proj, pedFile, samplesToEvaluate, outputSer);
					}
					try {
						BufferedReader reader = Files.getAppropriateReader(tmpHerit);
						int merlinIndex = ext.indexOfStr("Merlin_est.", Files.getHeaderOfFile(tmpHerit, log));

						if (merlinIndex < 0) {
							log.reportTimeError("Could not find Merlin_est. in " + tmpHerit);
							return null;
						}
						PrintWriter writer = new PrintWriter(new FileWriter(heritSummary));
						writer.println(Array.toStr(MERLIN_ADDITIONS) + "\t" + reader.readLine().trim());
						int index = 0;
						while (reader.ready()) {
							String[] line = reader.readLine().trim().split("\t");
							try {
								double est = Double.parseDouble(line[merlinIndex].replaceAll("%", ""));

								writer.println(index + "\t" + est + "\t" + Array.toStr(line));
							} catch (NumberFormatException nfe) {
								log.reportTimeWarning("Skipping line " + Array.toStr(line) + " , invalid estimate");
							}
							index++;
						}
						reader.close();
						writer.close();
					} catch (FileNotFoundException fnfe) {
						log.reportError("Error: file \"" + tmpHerit + "\" not found in current directory");
						return null;
					} catch (IOException ioe) {
						log.reportError("Error reading file \"" + tmpHerit + "\"");
						return null;
					} catch (Exception e) {
						log.reportError("Error writing to " + heritSummary);
						log.reportException(e);
					}
				} else {
					log.reportFileExists(heritSummary);
				}
				int numFamIndex = ext.indexOfStr("n_Families_size>1", Files.getHeaderOfFile(tmpHerit, log));
				int numSampsIndex = ext.indexOfStr("n_Samples", Files.getHeaderOfFile(tmpHerit, log));
				String[] line = Files.getHeaderOfFile(tmpHerit, "\t", new String[] { "Model" }, log);
				int numFam = Integer.parseInt(line[numFamIndex]);
				int numSamps = Integer.parseInt(line[numSampsIndex]);

				RScatter rScatter = new RScatter(heritSummary, heritRscript, ext.rootOf(heritSummary), heritPlot, MERLIN_ADDITIONS[0], new String[] { MERLIN_ADDITIONS[1] }, SCATTER_TYPE.POINT, log);
				rScatter.setyRange(new double[] { 0, 100 });
				rScatter.setxLabel("PC (" + oType + " - sorted)");
				rScatter.setTitle(iType + " " + bType + "; nInd=" + numSamps + ", nFam=" + numFam);
				rScatter.setyLabel("Percent Heritability");
				rScatter.setgPoint_SIZE(GEOM_POINT_SIZE.GEOM_POINT);
				rScatter.setOverWriteExisting(true);
				rScatter.execute();
				return rScatter;
			}

			return null;
		}

		public void setValid(boolean valid) {
			this.valid = valid;
		}

		public String getOutputSer() {
			return outputSer;
		}

		public void setOutputSer(String outputSer) {
			this.outputSer = outputSer;
		}

		public void setOutputRoot(String outputRoot) {
			this.outputRoot = outputRoot;
		}

		public String getOutputRank() {
			return outputRank;
		}

		public String getOutputSummary() {
			return outputSummary;
		}

	}

	private static double[][] loadIndeps(CorrectionEvaluator cEvaluator, String[] indepHeaders, double[][] indepMasks, String[] catHeaders, String[] catHeaderMask, Logger log) {
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
		for (int i = 0; i < catHeaders.length; i++) {

			if (parser.getStringDataForTitle(catHeaders[i]).length <= 0) {
				log.reportTimeError("Did not find " + catHeaders[i] + " to load for independent variable");
				return null;
			} else {
				String[] data = parser.getStringDataForTitle(catHeaders[i]);
				CategoricalPredictor predictor = new CategoricalPredictor(data, catHeaderMask, log);
				DummyCoding dummyCoding = predictor.createDummyCoding(false);
				for (int j = 0; j < dummyCoding.getDummyBoolean().length; j++) {
					double[] dummyData = Array.doubleArray(dummyCoding.getTitles().length, Double.NaN);
					if (dummyCoding.getDummyBoolean()[j]) {
						for (int k = 0; k < dummyData.length; k++) {
							dummyData[k] = dummyCoding.getDummyData()[k][j];
						}
					}
					extraIndeps[j] = Array.concatDubs(extraIndeps[j], dummyData);
				}
			}
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

	public static CorrectionIterator[] runAll(Project proj, String markesToEvaluate, String samplesToBuildModels, String outputDir, String pcFile, String pedFile, boolean svd, int numthreads) {
		proj.INTENSITY_PC_NUM_COMPONENTS.setValue(400);
		if (pcFile != null) {
			proj.INTENSITY_PC_FILENAME.setValue(pcFile);
		}
		System.out.println("JDOFJSDF remember the pcs");
		if (outputDir == null) {
			outputDir = proj.PROJECT_DIRECTORY.getValue() + ext.rootOf(proj.INTENSITY_PC_FILENAME.getValue()) + "_eval/";
		}
		if (!Files.exists(proj.MARKER_DATA_DIRECTORY.getValue() + "markers.0.mdRAF")) {
			proj.getLog().reportTimeWarning("Did not see " + proj.MARKER_DATA_DIRECTORY.getValue() + "markers.0.mdRAF, attempting to transpose now");
			TransposeData.transposeData(proj, 2000000000, false);
		}
//		if (!Files.exists(proj.SAMPLE_QC_FILENAME.getValue())) {
//			LrrSd.init(proj, null, null, null, null, numthreads);
//		}
		CorrectionIterator[] cIterators = getIterations(proj, markesToEvaluate, samplesToBuildModels, outputDir, svd, numthreads);
		ArrayList<RScatter> rScatters = new ArrayList<RScatter>();

		for (int i = 0; i < cIterators.length; i++) {
			cIterators[i].run();
			// IterationResult iterationResult = cIterators[i].getIterationResult();
			// RScatter rScatter = iterationResult.plotSummary(new String[] { "Rsquare_correction", "ICC_EVAL_CLASS_DUPLICATE_ALL", "ICC_EVAL_CLASS_DUPLICATE_SAME_VISIT", "ICC_EVAL_CLASS_FC", "SPEARMAN_CORREL_AGE", "SPEARMAN_CORREL_EVAL_DATA_SEX", "SPEARMAN_CORREL_EVAL_DATA_resid.mtDNaN.qPCR.MT001", "SPEARMAN_CORREL_EVAL_DATA_resid.mtDNA.qPCR", "SPEARMAN_CORREL_EVAL_DATA_Mt_DNA_relative_copy_number" }, proj.getLog());
			// rScatters.add(rScatter);
			// rScatters.add(iterationResult.plotRank(proj.getLog()));
			// if (pedFile != null) {
			// proj.getLog().reportTimeInfo("Since a ped file " + pedFile + " was provided, we will generate heritability estimates");
			// EvaluationResult.prepareHeritability(proj, pedFile, iterationResult.getOutputSer());
			// }

		}
		String[] plotTitlesForMain = new String[] { "Rsquare_correction", "ICC_EVAL_CLASS_DUPLICATE_ALL", "ICC_EVAL_CLASS_DUPLICATE_SAME_VISIT", "ICC_EVAL_CLASS_FC", "SPEARMAN_CORREL_AGE", "SPEARMAN_CORREL_EVAL_DATA_SEX", "SPEARMAN_CORREL_EVAL_DATA_resid.mtDNaN.qPCR.MT001", "SPEARMAN_CORREL_EVAL_DATA_resid.mtDNA.qPCR", "SPEARMAN_CORREL_EVAL_DATA_Mt_DNA_relative_copy_number", "SPEARMAN_CORREL_EVAL_DATA_Ratio.ND1", "SPEARMAN_CORREL_EVAL_DATA_qpcr.qnorm.exprs" };
		String[] plotTitlesForMito = new String[] { "Rsquare_correction", "SPEARMAN_CORREL_EVAL_DATA_Mt_DNA_relative_copy_number_NO_STRAT", "ICC_EVAL_CLASS_FC_NO_STRAT","SPEARMAN_CORREL_AGE_NO_STRAT", "SPEARMAN_CORREL_EVAL_DATA_SEX_NO_STRAT" };
		String[] plotTitlesForMitoFC = new String[] { "Rsquare_correction","ICC_EVAL_CLASS_FC_NO_STRAT" ,
				"SPEARMAN_CORREL_EVAL_DATA_Mt_DNA_relative_copy_number_NO_STRAT",
				"SPEARMAN_CORREL_EVAL_DATA_Mt_DNA_relative_copy_number_PT",
				"SPEARMAN_CORREL_EVAL_DATA_Mt_DNA_relative_copy_number_BU",
				"SPEARMAN_CORREL_EVAL_DATA_Mt_DNA_relative_copy_number_NY",
				"SPEARMAN_CORREL_EVAL_DATA_Mt_DNA_relative_copy_number_DK" };

		String[][] plotTitlesForSummary = new String[][] { plotTitlesForMain,plotTitlesForMito,plotTitlesForMitoFC };
		String[] subsetDataHeritability = new String[] { "EVAL_DATA_Mt_DNA_relative_copy_number" };
		IterSummaryProducer producer = new IterSummaryProducer(proj, cIterators, plotTitlesForSummary, pedFile);
		if (pedFile != null) {
			producer.setSubsetDataHeritability(subsetDataHeritability);
		}

		WorkerTrain<RScatter[]> summaryTrain = new WorkerTrain<RScatter[]>(producer, numthreads, numthreads, proj.getLog());
		while (summaryTrain.hasNext()) {
			RScatter[] rScattersTmp = summaryTrain.next();
			for (int i = 0; i < rScattersTmp.length; i++) {
				rScatters.add(rScattersTmp[i]);
			}
		}

		String outputRoot = outputDir + "finalSummary";
		RScatters finalScatters = new RScatters(rScatters.toArray(new RScatter[rScatters.size()]), outputRoot + ".rscript", outputRoot + ".pdf", COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_2, PLOT_DEVICE.PDF, proj.getLog());
		finalScatters.execute();
		return cIterators;
	}

	public static class IterSummaryProducer implements Producer<RScatter[]> {
		private Project proj;
		private CorrectionIterator[] cIterators;
		private String[][] plotTitlesForSummary;
		private String[] subsetDataHeritability;
		private String pedFile;
		private int index;

		public IterSummaryProducer(Project proj, CorrectionIterator[] cIterators, String[][] plotTitlesForSummary, String pedFile) {
			super();
			this.proj = proj;
			this.cIterators = cIterators;
			this.pedFile = pedFile;
			this.plotTitlesForSummary = plotTitlesForSummary;
			this.index = 0;
		}

		public void setSubsetDataHeritability(String[] subsetDataHeritability) {
			this.subsetDataHeritability = subsetDataHeritability;
		}

		@Override
		public boolean hasNext() {
			// TODO Auto-generated method stub
			return index < cIterators.length;
		}

		@Override
		public Callable<RScatter[]> next() {
			final CorrectionIterator tmp = cIterators[index];
			final Logger log = proj.getLog();
			Callable<RScatter[]> callable = new Callable<RScatter[]>() {
				@Override
				public RScatter[] call() throws Exception {
					IterationResult iterationResult = tmp.getIterationResult();
					ArrayList<RScatter> scatters = new ArrayList<RScatter>();
					// Add multiplots here
					for (int i = 0; i < plotTitlesForSummary.length; i++) {
						RScatter rScatterSummary = iterationResult.plotSummary(plotTitlesForSummary[i], i,proj.getLog());
						scatters.add(rScatterSummary);
					}
					if (pedFile != null) {

						if (subsetDataHeritability == null) {
							scatters.add(iterationResult.plotHeritability(proj, pedFile, iterationResult.getBasicPrep().getSamplesToEvaluate(), log));
							scatters.add(iterationResult.plotRank(log));
						} else {
							if (iterationResult.getBasicPrep() == null) {
								log.reportTimeError("must have data basic prep object");
								return null;
							} else {

								scatters.add(iterationResult.plotHeritability(proj, pedFile, iterationResult.getBasicPrep().getSamplesToEvaluate(), log));

								for (int i = 0; i < subsetDataHeritability.length; i++) {
									IterationResult tmpR = tmp.getIterationResult();

									int index = ext.indexOfStr(subsetDataHeritability[i], iterationResult.getBasicPrep().getNumericTitles());
									if (index >= 0) {
										String newSer = ext.addToRoot(iterationResult.getOutputSer(), "." + subsetDataHeritability[i] + ".summary");
										Files.copyFileUsingFileChannels(new File(iterationResult.getOutputSer()), new File(newSer), log);
										double[] data = iterationResult.getBasicPrep().getNumericData()[index];
										tmpR.setOutputRoot(ext.rootOf(ext.rootOf(newSer, false), false));
										tmpR.setOutputSer(newSer);

										boolean[] currentModel = new boolean[iterationResult.getBasicPrep().getSamplesToEvaluate().length];
										for (int j = 0; j < currentModel.length; j++) {
											if (iterationResult.getBasicPrep().getSamplesToEvaluate()[j] && !Double.isNaN(data[j])) {
												currentModel[j] = true;
											} else {
												currentModel[j] = false;
											}
										}
										RScatter tmpherit = tmpR.plotHeritability(proj, pedFile, currentModel, log);
										tmpherit.setTitle(subsetDataHeritability[i] + " - matched samples n=" + Array.booleanArraySum(currentModel));
										scatters.add(tmpR.plotHeritability(proj, pedFile, currentModel, log));
									} else {
										log.reportTimeWarning("Skipping " + subsetDataHeritability[i]);
									}
								}
								scatters.add(iterationResult.plotRank(log));
							}
						}
					} else {
						scatters.add(iterationResult.plotRank(log));
					}
					return scatters.toArray(new RScatter[scatters.size()]);
				}
			};
			index++;
			return callable;
		}

		@Override
		public void remove() {
			// TODO Auto-generated method stub

		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String proj = null;
		String markers = null;
		String defaultDir = null;
		boolean svd = false;
		int numThreads = 3;
		String samplesToBuildModels = null;
		String pcFile = null;
		String pedFile = null;
		String usage = "\n" + "cnv.analysis.pca.CorrectionEvaluator requires 0-1 arguments\n";
		usage += "   (1) project filename (i.e. proj=" + proj + " (default))\n" + "";
		usage += "   (2) markers to Evaluate (i.e. markers=" + markers + " (default))\n" + "";
		usage += PSF.Ext.getOutputDirCommand(3, defaultDir);
		usage += PSF.Ext.getNumThreadsCommand(4, numThreads);
		usage += "   (5) svd regression (i.e.-svd (not default))\n" + "";
		usage += "   (6) samples to generate models (i.e.samples= (no default))\n" + "";
		usage += "   (7) ped file to generate heritability (i.e.ped= (no default))\n" + "";
		usage += "   (8) alternate pc file to use (i.e.pcFile= (no default))\n" + "";

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
			} else if (args[i].startsWith("ped=")) {
				pedFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-svd")) {
				svd = true;
				numArgs--;
			} else if (args[i].startsWith("pcFile=")) {
				pcFile = ext.parseStringArg(args[i], "");
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
			runAll(new Project(proj, false), markers, samplesToBuildModels, defaultDir, pcFile, pedFile, svd, numThreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
