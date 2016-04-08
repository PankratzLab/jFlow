package one.JL;

import java.io.File;

import common.Array;
import common.Files;
import common.HashVec;
import common.PSF;
import common.ext;
import cnv.filesys.Project;
import cnv.filesys.Project.ARRAY;
import cnv.manage.ExtProjectDataParser;
import cnv.manage.Stats;
import cnv.qc.MarkerMetrics;
import cnv.qc.SexChecks;
import cnv.var.SampleData;

/**
 * testing some auto-marker selection for PCA
 *
 */
public class AutoMito {

	private static void run(Project proj, String name, String mitoMarks, int maxNumMarkers, double callRateSamp, double callRateMarker, double hetExMarker, double sexZscore, double lrrSDSamp, int numThreads) {
		long arrayLength = maxNumMarkers * proj.getSamples().length;
		if (arrayLength >= Integer.MAX_VALUE) {
			proj.getLog().reportTimeWarning("Maximum number of markers set to: " + maxNumMarkers);
			proj.getLog().reportTimeWarning("Number of samples : " + proj.getSamples().length);
			proj.getLog().reportTimeWarning(proj.getSamples().length + " X " + maxNumMarkers + " = " + arrayLength + " , which is greater than max java integer");
			maxNumMarkers = Math.round(Integer.MAX_VALUE / proj.getSamples().length + 1);
			proj.getLog().reportTimeWarning("Updated max num markers to " + maxNumMarkers);
		}

		String qcDir = proj.PROJECT_DIRECTORY.getValue() + name + "_PCA_INPUT_QC/";
		new File(qcDir).mkdirs();
		String temporarySampleData = qcDir + "sampleInfo.txt";
		proj.SAMPLE_DATA_FILENAME.setValue(temporarySampleData);
		SampleData.createMinimalSampleData(proj);

		String sexCheck = qcDir + ext.rootOf(proj.SEXCHECK_RESULTS_FILENAME.getValue());
		proj.SEXCHECK_RESULTS_FILENAME.setValue(sexCheck);
		if (!Files.exists(sexCheck)) {
			SexChecks.sexCheck(proj);
		}
		int sexIndex = ext.indexOfStr(SexChecks.EST_SEX_HEADER, Files.getHeaderOfFile(proj.SAMPLE_DATA_FILENAME.getValue(), proj.getLog()));
		if (sexIndex < 0) {
			throw new IllegalArgumentException("Could not find " + SexChecks.EST_SEX_HEADER + " in " + proj.SAMPLE_DATA_FILENAME.getValue());
		}
		String temporarySampleDataWithSex = qcDir + ext.addToRoot(temporarySampleData, ".withSex.txt");

		if (ext.indexOfStr(SexChecks.EST_SEX_HEADER, Files.getHeaderOfFile(proj.SAMPLE_DATA_FILENAME.getValue(), proj.getLog())) < 0) {
			String[][] matrix = HashVec.loadFileToStringMatrix(proj.SAMPLE_DATA_FILENAME.getValue(), false, null, false);
			String[][] withSex = new String[matrix.length][];
			withSex[0] = Array.concatAll(matrix[0], new String[] { "CLASS=" + SampleData.EUPHEMISMS[1] });
			for (int i = 1; i < withSex.length; i++) {
				int sex = Integer.parseInt(matrix[i][sexIndex]);
				withSex[i] = Array.concatAll(matrix[i], new String[] { (sex <= 2 ? sex + "" : "-1") });
				Files.writeMatrix(withSex, temporarySampleDataWithSex, "\t");
			}
		}
		proj.SAMPLE_DATA_FILENAME.setValue(temporarySampleDataWithSex);

		new File(qcDir).mkdirs();
		String baseMarkerQC = qcDir + name + "_base_markerQC.txt";
		if (!Files.exists(proj.MARKER_METRICS_FILENAME.getValue())) {
			MarkerMetrics.fullQC(proj, null, null, false, numThreads);
		}
		proj.getLog().reportTimeInfo("Loading " + baseMarkerQC);
		ExtProjectDataParser parser;
		parser = MarkerMetrics.developParser(proj, baseMarkerQC);

		proj = new Project(proj.getPropertyFilename(), false);// reset everything

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "AutoMito.dat";
		String mitoMarks = null;
		String name = "autoMarkerSelectionTester";
		int numThreads = 24;
		int maxNumMarkers = 1000000;
		double callRateSamp = .96;
		double callRateMarker = .96;
		double hetExMarker = .1;
		double sexZscore = 2.5;
		double lrrSDSamp = Double.NaN;

		String usage = "\n" +
				"one.JL.AutoMito requires 0-1 arguments\n" +
				"   (1) an existing project filename (i.e. proj=" + filename + " (default))\n" +
				"   (2) full path to mitochondrial markers (i.e. mitoMarks=" + mitoMarks + " (default))\n" +
				"   (3) analysis name (i.e. name=" + name + " (default))\n" +
				PSF.Ext.getNumThreadsCommand(4, numThreads) +
				"   (5) maximum Number of markers(i.e. maxNumMarkers=" + maxNumMarkers + " (default))\n" +
				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("mitoMarks=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("name=")) {
				name = args[i].split("=")[1];
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
			Project proj = new Project(filename, false);
			proj.setPropertyFilename(filename);
			if (Double.isNaN(lrrSDSamp)) {
				if (proj.getArrayType() == ARRAY.ILLUMINA) {
					lrrSDSamp = 0.30;
				} else if (proj.getArrayType() == ARRAY.AFFY_GW6 || proj.getArrayType() == ARRAY.AFFY_GW6_CN) {
					lrrSDSamp = 0.35;
				} else {
					throw new IllegalArgumentException("did not expect to see array type " + proj.getArrayType() + " here");
				}
			}
			run(proj, name, mitoMarks, maxNumMarkers, callRateSamp, callRateMarker, hetExMarker, sexZscore, lrrSDSamp, numThreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
