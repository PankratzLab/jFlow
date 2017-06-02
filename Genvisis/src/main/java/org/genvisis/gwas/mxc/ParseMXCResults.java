package org.genvisis.gwas.mxc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;

import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;
import org.genvisis.gwas.HitWindows;

public class ParseMXCResults {
	private static String resultsGeneCol = "gene_name";
	private static String genesGeneCol = "GRCh37.p5-Primary Assembly_name";
	private static String geneStartCol = "GRCh37.p5-Primary Assembly_start";
	private static String geneEndCol = "GRCh37.p5-Primary Assembly_stop";
	private static String chrCol = "GRCh37.p5-Primary Assembly_chr";
	private static String genesFile = "/home/pankrat2/public/bin/NCBI/genes37.xln";

	private static String usage = "";

	private static void cleanResultsFile(String inputFile, String outputFile,
																			 Logger log) throws IOException {
		String[][] mxc = mxcToMatrix(inputFile);

		String[] keys = getKeys(mxc, resultsGeneCol, log);

		String[][] results = findGenes(mxc, keys, genesFile, log);


		// output cleaned stuff to a file
		writeToFile(results, outputFile, log);
	}

	private static void writeToFile(String[][] matrix, String filename, Logger log) {
		PrintWriter writer;

		try {
			writer = new PrintWriter(new FileWriter(filename));
			for (String[] element : matrix) {
				for (int j = 0; j < element.length; j++) {
					writer.print((j == 0 ? "" : ",") + element[j]);
				}
				writer.println();
			}

			writer.close();
		} catch (Exception e) {
			log.report("Error writing to " + filename);
			e.printStackTrace();
		}
	}

	private static String[] getKeys(String[][] mxc, String colname, Logger log) {
		String[] headers = mxc[0];
		int index = -1;
		for (int i = 0; i < headers.length; i++)
			if (colname.equals(headers[i]))
				index = i;


		if (index == -1) {
			log.report("Could not find column: " + colname);
			System.exit(0);
		}

		return Matrix.extractColumn(mxc, index);
	}

	private static String[][] findGenes(String[][] mxc, String[] keys, String filename,
																			Logger log) throws IOException {
		String[][] results = null;
		String[] headers = {"gene_name", "start_pos",
												"stop_pos", "chr"};

		String[][] genes = loadGenes(filename, log);

		// get the results of the comparison
		try {
			results = Files.combineInMemory(keys, genes, "NA", true, false, log);
		} catch (Elision e) {
			log.reportError(e.getMessage());
		}

		// append results to the original matrix
		results = append(mxc, results);

		// headers are marked as NA because the column names don't match
		// put in R-safe headers in place of NA values
		for (int j = 0; j < results[0].length; j++) {
			results[0][j] = (j < mxc[0].length) ? mxc[0][j] : headers[j - mxc[0].length];
		}

		// filter out rows where no gene can be found
		results = removeNAs(results, genes, keys, results[0].length - 4);

		return results;
	}

	private static String[][] loadGenes(String filename, Logger log) {
		int[] cols = new int[4];

		try {
			BufferedReader b = new BufferedReader(new FileReader(filename));
			String[] line = b.readLine().split(Files.determineDelimiter(filename, log));

			for (int i = 0; i < line.length; i++) {
				if (genesGeneCol.equals(line[i]))
					cols[0] = i;
				else if (geneStartCol.equals(line[i]))
					cols[1] = i;
				else if (geneEndCol.equals(line[i]))
					cols[2] = i;
				else if (chrCol.equals(line[i]))
					cols[3] = i;
			}

			b.close();
		} catch (FileNotFoundException e) {
			log.report("ERROR - Could not find file " + filename);
			e.printStackTrace();
		} catch (IOException e) {
			log.report("ERROR - Problem reading file " + filename);
			e.printStackTrace();
		}

		return HashVec.loadFileToStringMatrix(filename, false, cols, false);
	}

	private static String[][] append(String[][] start, String[][] toAppend) {
		String[][] m = new String[start.length][start[0].length + toAppend[0].length];

		for (int i = 0; i < m.length; i++) {
			for (int j = 0; j < m[0].length; j++) {
				if (start[0].length <= j)
					m[i][j] = toAppend[i][j - start[0].length];
				else
					m[i][j] = start[i][j];
			}
		}

		return m;
	}

	private static String[][] removeNAs(String[][] mxc, String[][] genes, String[] keys, int col) {
		boolean[] rowsToKeep = new boolean[mxc.length];
		if (col < 0)
			return mxc;

		for (int i = 0; i < mxc.length; i++) {
			if (mxc[i][col].equals("NA")) {
				rowsToKeep[i] = false;
				if (keys != null) {
					// check for synonyms and fix results if possible
					String[] synInfo = getSynonym(keys[i]);
					if (synInfo != null) {
						mxc[i][col] = synInfo[0];
						mxc[i][col + 1] = synInfo[1];
						mxc[i][col + 2] = synInfo[2];
						mxc[i][col + 3] = synInfo[3];
						rowsToKeep[i] = true;
					}
				}
			} else {
				rowsToKeep[i] = true;
			}
		}

		mxc = Matrix.subset(mxc, rowsToKeep);
		return mxc;
	}

	private static String[] getSynonym(String key) {
		return null;
	}


	// TODO: Use genvisis qq and manhattan plots instead of R's
	private static void generateRScript(String outputFile) {
		ArrayList<String> r = new ArrayList<String>();
		outputFile = new File(outputFile).getAbsolutePath();

		String filename = outputFile.split("\\.")[0];

		r.add("library(\"qqman\", lib.loc=\"/panfs/roc/groups/5/pankrat2/cole0482/R/x86_64-pc-linux-gnu-library/3.2\")");
		r.add("data=read.csv(\"" + filename + ".csv" + "\")");
		r.add("png(\"" + filename + "_manhattan.png\")");
		r.add("manhattan(data, chr=\"chr\", bp=\"start_pos\", p=\"pvalue\")");
		r.add("dev.off()");


		// write rscript to a file
		Files.writeIterable(r, filename + "_manhattan.R");

		r.clear();
		r.add("library(\"qqman\", lib.loc=\"/panfs/roc/groups/5/pankrat2/cole0482/R/x86_64-pc-linux-gnu-library/3.2\")");
		r.add("data=read.csv(\"" + filename + ".csv" + "\")");
		r.add("png(\"" + filename + "_qq.png\")");
		r.add("qq(data$pvalue)");
		r.add("dev.off()");

		Files.writeIterable(r, filename + "_qq.R");
	}


	public static void markersUsed(String mapfile, String posfile, String mxcfile, String metalfile) {
		Hashtable<String, String> mxc = HashVec.loadFileToHashString(mxcfile, "gene",
																																 new String[] {"gene_name"}, ",");
		Hashtable<String, String> metal = HashVec.loadFileToHashString(metalfile, "MarkerName",
																																	 new String[] {"P-value"}, " ");
		String[][] map = HashVec.loadFileToStringMatrix(mapfile, true, new int[] {0, 1}, false);
		String[][] pos = HashVec.loadFileToStringMatrix(posfile, true, null, false);

		float sig = (float) (0.05 / metal.keySet().size());
		float sugg = sig * 100;

		String root = ext.parseDirectoryOfFile(metalfile);

		// maps gene to related rsids
		Hashtable<String, HashSet<String>> markers = new Hashtable<String, HashSet<String>>();

		for (int j = 0; j < map.length; j++) {
			String key = map[j][0].split("\\.")[0];
			if (!markers.containsKey(key)) {
				markers.put(key, new HashSet<String>());
			}
			markers.get(key).add(map[j][1]);
		}


		// maps rsid to a chr and position
		Hashtable<String, String> posmap = new Hashtable<String, String>();
		for (int i = 0; i < pos.length; i++) {
			posmap.put(pos[i][0], pos[i][1] + "," + pos[i][2]);
		}

		ArrayList<String> results = new ArrayList<String>();
		results.add("MarkerName,p,sig,Chr,Position");
		ArrayList<String> mxc_results = new ArrayList<String>();
		mxc_results.add("MarkerName,Gene,p,sig,Chr,Position");
		for (String s : metal.keySet()) {
			if (!metal.get(s).equals("NA")) {
				double p = Double.parseDouble(metal.get(s));
				String ind = p < sig ? "***" : (p < sugg ? "**" : ".");
				String chrAndPos = posmap.get(s);

				results.add(s + "," + p + "," + ind + "," + chrAndPos);
			}
		}

		for (String g : mxc.keySet()) {
			HashSet<String> rsids = markers.get(g);
			if (rsids != null) {
				for (String s : rsids) {
					if (metal.get(s) != null && !metal.get(s).equals("NA")) {
						double p = Double.parseDouble(metal.get(s));
						String ind = p < sig ? "***" : (p < sugg ? "**" : ".");
						String chrAndPos = posmap.get(s);
						String geneName = mxc.get(g);
						mxc_results.add(s + "," + geneName + "," + p + "," + ind + "," + chrAndPos);
					}
				}
			}
		}

		Files.writeIterable(results, root + "relatedMarkers_all.csv");
		Files.writeIterable(mxc_results, root + "relatedMarkers_mxc.csv");

		String[][] regions = HitWindows.determine(root + "relatedMarkers_all.csv",
																							sig, 500000, sugg, new String[] {},
																							new Logger());
		String[][] mxc_regions = HitWindows.determine(root + "relatedMarkers_mxc.csv",
																									sig, 500000, sugg, new String[] {"Gene"},
																									new Logger());

		Files.writeMatrix(regions, root + "metal_regions.csv", ",");
		Files.writeMatrix(mxc_regions, root + "mxc_regions.csv", ",");
	}

	private static String[][] mxcToMatrix(String filename) {
		String[][] mxc = HashVec.loadFileToStringMatrix(filename, false, null, false);
		int p = -1;
		for (int i = 0; i < mxc[0].length; i++)
			if (mxc[0][i].equals("pvalue"))
				p = i;
		mxc = removeNAs(mxc, null, null, p);
		return mxc;
	}

	public static void main(String[] args) throws IOException {
		if (args.length == 0) {
			System.out.println(usage);
			System.exit(1);
		}
		Logger log = new Logger();

		String mxcFile = "";
		String outputFile = "";
		String metal = null;
		String map = null;
		String posmap = null;

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("help")) {
				System.out.println(usage);
				System.exit(1);
			} else if (arg.startsWith("mxc="))
				mxcFile = arg.split("=")[1];
			else if (arg.startsWith("genes="))
				genesFile = arg.split("=")[1];
			else if (arg.startsWith("out="))
				outputFile = arg.split("=")[1];
			else if (arg.startsWith("chr="))
				chrCol = arg.split("=")[1];
			else if (arg.startsWith("mGeneCol="))
				resultsGeneCol = arg.split("=")[1];
			else if (arg.startsWith("gGeneCol="))
				genesGeneCol = arg.split("=")[1];
			else if (arg.startsWith("startPos="))
				geneStartCol = arg.split("=")[1];
			else if (arg.startsWith("endPos="))
				geneEndCol = arg.split("=")[1];
			else if (arg.startsWith("metal="))
				metal = arg.split("=")[1];
			else if (arg.startsWith("map="))
				map = arg.split("=")[1];
			else if (arg.startsWith("posmap="))
				posmap = arg.split("=")[1];
			else
				log.report("Invalid argument: " + arg);
		}

		if (metal != null && map != null) {
			markersUsed(map, posmap, mxcFile, metal);
		} else {
			cleanResultsFile(mxcFile, outputFile, log);
			generateRScript(outputFile);
		}
	}
}
