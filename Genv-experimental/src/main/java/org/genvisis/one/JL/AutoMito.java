package org.genvisis.one.JL;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.manage.ExtProjectDataParser;
import org.genvisis.cnv.manage.ExtProjectDataParser.ProjectDataParserBuilder;
import org.genvisis.cnv.qc.LrrSd;
import org.genvisis.cnv.qc.MarkerMetrics;
import org.genvisis.cnv.qc.SampleQC;
import org.genvisis.cnv.qc.SexChecks;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;

/**
 * testing some auto-marker selection for PCA
 *
 */
public class AutoMito {

	private static void run(Project proj, String name, String mitoMarks, int maxNumMarkers,
													double callRateSamp, double callRateMarker, double hetExMarker,
													double sexZscore, double lrrSDSamp,
													int numThreads) throws FileNotFoundException {
		long arrayLength = maxNumMarkers * proj.getSamples().length;
		if (arrayLength >= Integer.MAX_VALUE) {
			proj.getLog().reportTimeWarning("Maximum number of markers set to: " + maxNumMarkers);
			proj.getLog().reportTimeWarning("Number of samples : " + proj.getSamples().length);
			proj.getLog().reportTimeWarning(proj.getSamples().length + " X " + maxNumMarkers + " = "
																			+ arrayLength + " , which is greater than max java integer");
			maxNumMarkers = Math.round(Integer.MAX_VALUE / (proj.getSamples().length + 1f));
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
			SexChecks.sexCheck(proj, true);
		}

		String temporarySampleDataWithSex = qcDir + ext.addToRoot(temporarySampleData, ".withSex.txt");

		if (ext.indexOfStr(SexChecks.EST_SEX_HEADER,
											 Files.getHeaderOfFile(proj.SAMPLE_DATA_FILENAME.getValue(),
																						 proj.getLog())) < 0) {
			ProjectDataParserBuilder builder = new ProjectDataParserBuilder();
			builder.sampleBased(true);
			builder.dataKeyColumnName("Sample");
			builder.treatAllNumeric(false);
			builder.stringDataTitles(new String[] {SexChecks.EST_SEX_HEADER});
			builder.requireAll(true);
			System.out.println(ext.indexOfStr(SexChecks.EST_SEX_HEADER,
																				Files.getHeaderOfFile(proj.SEXCHECK_RESULTS_FILENAME.getValue(),
																															proj.getLog())));
			ExtProjectDataParser parser = builder.build(proj, proj.SEXCHECK_RESULTS_FILENAME.getValue());

			String[][] matrix = HashVec.loadFileToStringMatrix(proj.SAMPLE_DATA_FILENAME.getValue(),
																												 false, null, false);
			int sexIndex = ext.indexOfStr("CLASS=" + SampleData.EUPHEMISMS[1], matrix[0]);

			for (int i = 1; i < matrix.length; i++) {
				matrix[i][sexIndex] = parser.getStringDataForTitle(SexChecks.EST_SEX_HEADER)[i - 1];
			}
			Files.writeMatrix(matrix, temporarySampleDataWithSex, "\t");

		}
		proj.SAMPLE_DATA_FILENAME.setValue(temporarySampleDataWithSex);

		String baseSampleQC = qcDir + name + "_baseSampleQC.txt";
		proj.SAMPLE_QC_FILENAME.setValue(baseSampleQC);

		new File(qcDir).mkdirs();
		// MarkerSet markerSet = proj.getMarkerSet();// remove CN only

		String baseMarkerQC = qcDir + name + "_base_markerQC.txt";
		String baseMarkers = ext.addToRoot(baseMarkerQC, ".nonCNOnlyMarkers");
		proj.MARKER_METRICS_FILENAME.setValue(baseMarkerQC);
		ArrayList<String> nonCN_Only = new ArrayList<String>();
		String[] autos = proj.getAutosomalMarkers();
		for (int i = 0; i < autos.length; i++) {
			if (!proj.getArrayType().isCNOnly(autos[i])) {
				nonCN_Only.add(autos[i]);
			}
		}
		Files.writeIterable(nonCN_Only, baseMarkers);
		if (!Files.exists(proj.SAMPLE_QC_FILENAME.getValue())) {
			LrrSd.init(proj, null, baseMarkers, baseMarkers, numThreads, null, false);
		}

		SampleQC sampleQC = SampleQC.loadSampleQC(proj);
		boolean[] samplesPassing = ArrayUtils.booleanArray(proj.getSamples().length, false);
		String sampleFiltRound1 = qcDir + name + "_base_sampleFilter.txt";
		double[] lrr = sampleQC.getDataFor("LRR_SD");
		double[] callRate = sampleQC.getDataFor("Genotype_callrate");
		for (int i = 0; i < callRate.length; i++) {
			if (callRate[i] > callRateSamp && lrr[i] < lrrSDSamp) {
				samplesPassing[i] = true;
			}
		}
		proj.getLog()
				.reportTimeInfo(ArrayUtils.booleanArraySum(samplesPassing) + " samples passed initial QC");
		Files.writeArray(ArrayUtils.subArray(proj.getSamples(), samplesPassing), sampleFiltRound1);

		if (!Files.exists(proj.MARKER_METRICS_FILENAME.getValue())) {
			MarkerMetrics.fullQC(proj, samplesPassing, baseMarkers, false, numThreads);
		}

		proj.getLog().reportTimeInfo("Loading " + baseMarkerQC);
		ExtProjectDataParser parser = MarkerMetrics.developParser(proj, baseMarkerQC);
		parser.getNumericDataForTitle("CallRate");
		parser.getNumericDataForTitle("HetEx");
		parser.getNumericDataForTitle("LRR_SEX_z");

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

		String usage = "\n" + "one.JL.AutoMito requires 0-1 arguments\n"
									 + "   (1) an existing project filename (i.e. proj=" + filename + " (default))\n"
									 + "   (2) full path to mitochondrial markers (i.e. mitoMarks=" + mitoMarks
									 + " (default))\n" + "   (3) analysis name (i.e. name=" + name + " (default))\n"
									 + PSF.Ext.getNumThreadsCommand(4, numThreads)
									 + "   (5) maximum Number of markers(i.e. maxNumMarkers=" + maxNumMarkers
									 + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("proj=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("mitoMarks=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("name=")) {
				name = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(arg);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
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
				} else if (proj.getArrayType() == ARRAY.AFFY_GW6
									 || proj.getArrayType() == ARRAY.AFFY_GW6_CN) {
					lrrSDSamp = 0.35;
				} else {
					throw new IllegalArgumentException("did not expect to see array type "
																						 + proj.getArrayType() + " here");
				}
			}
			run(proj, name, mitoMarks, maxNumMarkers, callRateSamp, callRateMarker, hetExMarker,
					sexZscore, lrrSDSamp, numThreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
