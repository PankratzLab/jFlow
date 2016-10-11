package org.genvisis.cnv.manage;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.stats.ContingencyTable;
import org.genvisis.stats.ProbDist;
import org.genvisis.stats.Stats;

public class temp {
	public static void yap(Project proj) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;

		String filename =
										"C:\\Documents and Settings\\npankrat\\My Documents\\CNV_PD\\results\\paperComp\\finalLargeRare\\nums.txt";
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(filename + "_chi.xln"));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				writer.print(Array.toStr(line)	+ "\t"
											+ Stats.FishersExact(	Double.parseDouble(line[0]), Double.parseDouble(line[1]),
																						Double.parseDouble(line[2]),
																						Double.parseDouble(line[3]), true));
				writer.println("\t" + ProbDist.ChiDist(	ContingencyTable.ChiSquare(
																																					new double[][] {{	Double.parseDouble(line[0]),
																																														Double.parseDouble(line[1])},
																																													{	Double.parseDouble(line[2]),
																																														Double.parseDouble(line[3])}},
																																					true, true),
																								1));
			}
			writer.close();
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}

		System.exit(1);

		boolean wassup = Boolean.valueOf("true");
		System.out.println(wassup);

		System.exit(1);

		Files.list("", "", true);

		System.exit(1);

		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV_PD\\data\\";
		String[] peeps = HashVec.loadFileToStringArray(dir	+ "peeps.txt", false, false,
																										new int[] {0, 1}, true);
		String[] cnvs = HashVec.loadFileToStringArray(dir	+ "parkin_introns.cnv", false, false,
																									new int[] {2, 3, 4, 5, 6, 7}, true);

		try {
			writer = new PrintWriter(new FileWriter(dir + "real_introns.cnv"));
			writer.println("FID\tIID\t" + cnvs[0]);
			for (String peep : peeps) {
				for (int j = 1; j < cnvs.length; j++) {
					writer.println(peep + "\t" + cnvs[j]);
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to file");
			e.printStackTrace();
		}

		System.exit(1);

		System.out.println(proj.getSamples().length + " samples");
		System.out.println(proj.getMarkerNames().length + " markers");

		System.exit(1);

		// PrintWriter writer;
		String[] markers = new String[] {"cnvi0017132", "cnvi0021043", "rs7435827", "rs6552182"};
		MarkerData[] markerData = MarkerSet.loadFromList(proj, markers);
		String[] samples = proj.getSamples();
		double[] sums = new double[samples.length];
		int[] counts = new int[samples.length];
		float[] lrrs;

		for (MarkerData element : markerData) {
			lrrs = element.getLRRs();
			for (int j = 0; j < samples.length; j++) {
				if (!Double.isNaN(lrrs[j])) {
					sums[j] += lrrs[j];
					counts[j]++;
				}
			}
		}

		try {
			writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + "LRRmeans.xln"));
			for (int i = 0; i < samples.length; i++) {
				writer.println(samples[i] + "\t" + (sums[i] / counts[i]));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing LRR means");
			e.printStackTrace();
		}

		System.exit(1);

		// String[] data = new String[] {"B", "C", "A", "D", "E"};
		//
		// int[] keys = Sort.quicksort(data);
		// int[] again = Sort.quicksort(keys);
		// for (int i = 0; i < keys.length; i++) {
		// System.out.println(data[i]+"\t"+keys[i]+"\t"+again[i]);
		// }
		//
		//
		// System.exit(1);
		//

		// String[] markerNames = proj.getMarkerNames();
		// Sample.load(proj.getDir(Project.IND_DIRECTORY)+"96M5656.samp",
		// false).writeToFile(markerNames, "96M5656.xln");
		// FullSample.load(proj.getDir(Project.SAMPLE_DIRECTORY)+"96M5656.fsamp",
		// false).writeToFile(markerNames, "96M5656f.xln");
		// proj.getMarkerSet().writeToFile("MarkerSet.xln");
		// String[] line =
		// proj.getMarkerLookup().get("rs2627690").split("[\\s]+");
		// MarkerData[] markerData =
		// MarkerDataCollection.load(proj.getDir(Project.PLOT_DIRECTORY)+line[0],
		// false).getCollection();
		// markerData[Integer.parseInt(line[1])].writeToFile(proj.getSamples(),
		// "rs2627690.xln");
	}

	public static void main(String[] args) throws IOException {
		try {
			yap(new Project(org.genvisis.cnv.Launch.getDefaultDebugProjectFile(true), false));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
