package cnv.analysis.pca;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;

import stats.StatsCrossTabs.STAT_TYPE;
import stats.StatsCrossTabs.StatsCrossTabRank;
import stats.StatsCrossTabs.VALUE_TYPE;
import common.Array;
import common.Files;
import common.Logger;
import common.PSF;
import common.WorkerTrain;
import common.ext;
import cnv.filesys.Project;
import cnv.manage.ExtProjectDataParser;
import cnv.qc.LrrSd;

class CorrectionIterator {

	public enum ITERATION_TYPE {
		/**
		 * The evaluation happens with the addition of other independent variables in addition to PCs
		 */
		WITHOUT_INDEPS, WITH_INDEPS;

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
		 * PCs are filtered for an association with known QC metrics from {@link LrrSd}
		 */
		QC_ASSOCIATION;

	}

	private static String run(Project proj, String markesToEvaluate, String samplesToBuildModels, ITERATION_TYPE iType, ORDER_TYPE oType, String outputDir, boolean svd, int numthreads) {
		Logger log = proj.getLog();
		if (outputDir == null) {
			outputDir = proj.PROJECT_DIRECTORY.getValue() + ext.rootOf(proj.INTENSITY_PC_FILENAME.getValue()) + "_eval/";

		}
		String output = outputDir + "correctionEval_" + iType + "_" + oType;
		String ser = output + ".summary.ser";
		if (!Files.exists(ser)) {
			proj.getLog().reportTimeInfo("Loading " + proj.INTENSITY_PC_FILENAME.getValue());

			new File(outputDir).mkdirs();
			log.reportTimeInfo("Beginning iteration evaluation:");
			log.reportTimeInfo("PC file: " + proj.INTENSITY_PC_FILENAME.getValue());
			log.reportTimeInfo("Iteration type : " + iType);
			log.reportTimeInfo("Order type : " + oType);

			PrincipalComponentsResiduals pcResiduals = proj.loadPcResids();
			pcResiduals.setMarkersToAssessFile(markesToEvaluate);
			pcResiduals.setHomozygousOnly(true);
			pcResiduals.computeAssessmentDataMedians();

			boolean[] samplesForModels = proj.getSamplesToInclude(samplesToBuildModels, true, true);

			CorrectionEvaluator cEvaluator = new CorrectionEvaluator(proj, pcResiduals, null, null, null, svd);
			int[] order = null;
			double[][] extraIndeps = null;
			StatsCrossTabRank sTabRank = null;
			boolean valid = true;
			switch (iType) {
			case WITHOUT_INDEPS:
				log.reportTimeInfo("Evaluating with " + Array.booleanArraySum(samplesForModels) + " samples, no additional independent variables");
				break;
			case WITH_INDEPS:
				extraIndeps = loadIndeps(cEvaluator, CorrectionEvaluator.INDEPS, new double[][] { { 0, 3, 4, 5, 6, 7, 8, 9, 10 }, { -1 } }, log);
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

			if (valid) {
				sTabRank = pcResiduals.getStatRankFor(pcResiduals.getMedians(), extraIndeps, samplesForModels, "RAW_MEDIANS", STAT_TYPE.LIN_REGRESSION, VALUE_TYPE.STAT, proj.getLog());
				sTabRank.dump(output + sTabRank.getRankedTo() + ".rank.txt", oType != ORDER_TYPE.NATURAL, log);
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
					order = PCSelector.select(proj, 0.10);
					if (order.length < 1) {
						log.reportTimeError("Could not select PCs from QC metrics");
						return null;
					}
				default:
					break;
				}
				boolean[] samplesToEvaluate = proj.getSamplesToInclude(null);
				log.reportTimeInfo(Array.booleanArraySum(samplesForModels) + " samples for models");
				log.reportTimeInfo(Array.booleanArraySum(samplesToEvaluate) + " samples for evaluation");

				cEvaluator = new CorrectionEvaluator(proj, pcResiduals, order, new boolean[][] { samplesForModels, samplesToEvaluate }, extraIndeps, svd);
				String summary = output + ".summary.txt";
				ArrayList<EvaluationResult> store = new ArrayList<EvaluationResult>();

				try {
					PrintWriter writer = new PrintWriter(new FileWriter(summary));
					WorkerTrain<EvaluationResult> train = new WorkerTrain<EvaluationResult>(cEvaluator, numthreads, numthreads, proj.getLog());
					int index = 0;
					while (train.hasNext()) {
						EvaluationResult result = train.next();
						result.setItType(iType);
						result.setOrType(oType);
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
					proj.getLog().reportError("Error writing to " + summary);
					proj.getLog().reportException(e);
				}
				EvaluationResult.serialize(store.toArray(new EvaluationResult[store.size()]), ser);
			}
		} else {
			log.reportFileExists(ser);

		}
		return ser;
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

	public static void runAll(Project proj, String markesToEvaluate, String samplesToBuildModels, String outputDir, boolean svd, int numthreads) {
		proj.INTENSITY_PC_NUM_COMPONENTS.setValue(900);
		System.out.println("JDOFJSDF remember the pcs");
		for (int i = 0; i < ITERATION_TYPE.values().length; i++) {
			for (int j = 0; j < ORDER_TYPE.values().length; j++) {
				String ser = run(proj, markesToEvaluate, samplesToBuildModels, ITERATION_TYPE.values()[i], ORDER_TYPE.values()[j], outputDir, svd, numthreads);
				EvaluationResult[] results = EvaluationResult.readSerial(ser, proj.getLog());
				for (int k = 0; k < results.length; k++) {
					// System.out.println(Array.toStr(results[i].getData()));
				}
			}
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
