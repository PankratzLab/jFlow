package org.genvisis.cnv.qc;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;

import javax.jms.IllegalStateException;

import org.genvisis.cnv.analysis.pca.PrincipalComponentsCompute;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.ExtProjectDataParser;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.PSF;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

/**
 * @author lane0212 Detecting batch effects similar to http://www.ncbi.nlm.nih.gov/pubmed/23958724
 */
public class BatchEffectDetector {

	private static String[] determineMarkers(Project proj, String outDir, int numMarks, int threads) throws FileNotFoundException {
		if (numMarks < 0) {
			return proj.getAutosomalMarkers();
		} else {
			String qcFile = outDir + ext.removeDirectoryInfo(proj.MARKER_METRICS_FILENAME.getValue());
			if (!Files.exists(qcFile)) {
				String[] autosomal = proj.getAutosomalNonCNMarkers();
				Files.writeList(autosomal, outDir + "detectMarkers.txt");

				proj.MARKER_METRICS_FILENAME.setValue(qcFile);
				MarkerMetrics.fullQC(proj, null, outDir + "detectMarkers.txt", false, threads);
				// TODO, throw command to run from workflow
				// throw new FileNotFoundException("Sorry, marker metrics required: Perform step X of the workflow to generate -> " + proj.MARKER_METRICS_FILENAME.getValue());
			}
			ExtProjectDataParser parser = MarkerMetrics.developParser(proj, qcFile);
			double[] lrrSD = parser.getNumericDataForTitle("SD_LRR");

			String[] markers = Array.subArray(proj.getMarkerNames(), parser.getDataPresent());
			int[] order = Sort.trickSort(lrrSD, markers);

			ArrayList<String> markersToUse = new ArrayList<String>();

			String summary = outDir + "markerSDSummary";
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(summary));
				writer.println("MarkerName\tSD_RANK\tSD");
				for (int i = 0; i < order.length; i++) {
					if (markersToUse.size() < numMarks) {
						String mark = markers[order[i]];
						markersToUse.add(mark);
						writer.println(mark + "\t" + (i + 1) + "\t" + lrrSD[order[i]]);
					} else {
						break;
					}
				}
				writer.close();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + summary);
				proj.getLog().reportException(e);
			}

			return null;

		}
	}

	private static void loadData(Project proj, String[] markers, String outputBase) {
		double[][] data = PrincipalComponentsCompute.getData(proj, markers, null, false, true, true, true, outputBase, null);
	}

	private static void run(Project proj, String outDir, int numMarks, int threads) throws FileNotFoundException {
		outDir = proj.PROJECT_DIRECTORY.getValue() + outDir;
		new File(outDir).mkdirs();

		String[] markers = determineMarkers(proj, outDir, numMarks, threads);

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		int numthreads = 6;
		String outDir = "BatchEffects/";
		String batchFile = "batches.txt";
		int numSDMarkers = 5000;
		String logfile = null;

		String usage = "\n" +
				"cnv.qc.BatchEffectDetector requires 0-1 arguments\n" +
				"   (1) project (i.e. proj=" + filename + " (default))\n" +
				"   (2) output directory relative to project (i.e. outDir=" + outDir + " (default))\n" +

				"   (3) batch definition file (i.e. outDir=" + batchFile + " (default))\n" +

				"   (4) num markers (SD ranked) (i.e. numSDMarkers=" + numSDMarkers + " (default))\n" +
				PSF.Ext.getNumThreadsCommand(5, numthreads) +

				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("outDir=")) {
				outDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("batchFile=")) {
				batchFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("numSDMarkers=")) {
				numSDMarkers = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numthreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
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
			run(proj, outDir, numSDMarkers, numthreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
