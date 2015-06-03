package cnv.analysis.pca;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

import stats.StatsCrossTabs.STAT_TYPE;
import stats.StatsCrossTabs.StatsCrossTabRank;
import stats.StatsCrossTabs.VALUE_TYPE;
import common.Array;
import common.Logger;
import common.PSF;
import common.WorkerTrain;
import common.ext;
import cnv.filesys.Project;
import cnv.manage.ExtProjectDataParser;

class CorrectionIterator {

	public enum ITERATION_TYPE {
		/**
		 * The evaluation happens with the addition of variables in addition to PCs
		 */
		WITH_INDEPS, WITHOUT_INDEPS;

	}

	public enum ORDER_TYPE {
		/**
		 * PCs are regressed by their natural ordering (PC1,2,3)
		 */
		NATURAL, /**
		 * PCS are ranked by the amount of variance in the estimate of interest
		 */
		RANK;

	}

	private static void run(Project proj, String markesToEvaluate, ITERATION_TYPE iType, ORDER_TYPE oType, String outputDir, int numthreads) {
		proj.getLog().reportTimeInfo("Loading " + proj.INTENSITY_PC_FILENAME.getValue());
		if (outputDir == null) {
			outputDir = proj.PROJECT_DIRECTORY.getValue() + "mitoEval/";

		}
		new File(outputDir).mkdirs();
		Logger log = proj.getLog();
		log.reportTimeInfo("Beginning iteration evaluation:");
		log.reportTimeInfo("PC file: " + proj.INTENSITY_PC_FILENAME.getValue());
		log.reportTimeInfo("Iteration type : " + iType);
		log.reportTimeInfo("Order type : " + oType);

		PrincipalComponentsResiduals pcResiduals = proj.loadPcResids();
		pcResiduals.setMarkersToAssessFile(markesToEvaluate);
		pcResiduals.computeAssessmentDataMedians();
		String output = outputDir + "correctionEval_" + iType + "_" + oType;
		boolean[] samplesToInclude = proj.getSamplesToInclude(null);

		CorrectionEvaluator cEvaluator = new CorrectionEvaluator(proj, pcResiduals, null, null, null);
		int[] order = null;
		double[][] extraIndeps = null;
		StatsCrossTabRank sTabRank = null;
		boolean valid = true;
		switch (iType) {
		case WITHOUT_INDEPS:
			log.reportTimeInfo("Evaluating with " + Array.booleanArraySum(samplesToInclude) + " samples");
			break;
		case WITH_INDEPS:
			extraIndeps = loadIndeps(cEvaluator, CorrectionEvaluator.INDEPS, new double[][] { { 0, 3, 4, 5, 6, 7, 8, 9, 10 }, {} }, log);
			if (extraIndeps == null) {
				log.reportTimeError("type = " + iType + " and were missing some of the following " + Array.toStr(CorrectionEvaluator.INDEPS));
				valid = false;
			} else {
				boolean[] tmpInclude = new boolean[samplesToInclude.length];
				for (int i = 0; i < extraIndeps.length; i++) {

					if (samplesToInclude[i]) {
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
				log.reportTimeInfo("Original number of samples: " + Array.booleanArraySum(samplesToInclude));
				log.reportTimeInfo("Number of samples with valid independant variables: " + Array.booleanArraySum(tmpInclude));
				samplesToInclude = tmpInclude;
			}
			break;
		default:
			break;

		}

		if (valid) {
			sTabRank = pcResiduals.getStatRankFor(pcResiduals.getMedians(), null, "RAW_MEDIANS", STAT_TYPE.LIN_REGRESSION, VALUE_TYPE.STAT, proj.getLog());

			sTabRank.dump(output + ".rank.txt", oType != ORDER_TYPE.NATURAL, log);

			switch (oType) {
			case NATURAL:
				order = null;
				break;
			case RANK:
				order = new int[sTabRank.getOrder().length];
				for (int i = 0; i < sTabRank.getOrder().length; i++) {
					order[i] = sTabRank.getOrder()[i] + 1;
				}
				break;
			default:
				break;
			}
			cEvaluator = new CorrectionEvaluator(proj, pcResiduals, order, samplesToInclude, extraIndeps);
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(output));
				WorkerTrain<EvaluationResult> train = new WorkerTrain<EvaluationResult>(cEvaluator, numthreads, numthreads, proj.getLog());
				int index = 0;
				while (train.hasNext()) {
					EvaluationResult result = train.next();
					if (index == 0) {
						writer.println(Array.toStr(result.getHeader()));
					}
					writer.println(Array.toStr(result.getData()));
					index++;
				}
				index++;

				writer.close();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + output);
				proj.getLog().reportException(e);
			}
		}
	}

	private static double[][] loadIndeps(CorrectionEvaluator cEvaluator, String[] indepHeaders, double[][] indepMasks, Logger log) {
		ExtProjectDataParser parser = cEvaluator.getParser();

		double[][] extraIndeps = new double[Array.booleanArraySum(parser.getDataPresent())][];
		for (int i = 0; i < indepHeaders.length; i++) {
			if (parser.getStringDataForTitle(indepHeaders[i]).length <= 0) {
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
						extraIndeps[i][j] = Double.NaN;
					} else {
						extraIndeps[i][j] = data[j];

					}
				}
			}
			log.reportTimeInfo(Array.removeNaN(extraIndeps[i]).length + " samples for independant variable " + indepHeaders[i]);
		}
		return extraIndeps;
	}

	public static void runAll(Project proj, String markesToEvaluate, String outputDir, int numthreads) {
		for (int i = 0; i < ITERATION_TYPE.values().length; i++) {
			for (int j = 0; j < ORDER_TYPE.values().length; j++) {
				run(proj, markesToEvaluate, ITERATION_TYPE.values()[i], ORDER_TYPE.values()[j], outputDir, numthreads);
			}
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String proj = null;
		String markers = null;
		String defaultDir = null;
		int numThreads = 3;

		String usage = "\n" + "cnv.analysis.pca.CorrectionEvaluator requires 0-1 arguments\n";
		usage += "   (1) project filename (i.e. proj=" + proj + " (default))\n" + "";
		usage += "   (2) markers to Evaluate (i.e. markers=" + markers + " (default))\n" + "";
		usage += PSF.Ext.getOutputDirCommand(3, defaultDir);
		usage += PSF.Ext.getNumThreadsCommand(4, numThreads);

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
			runAll(new Project(proj, false), markers, defaultDir, numThreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
