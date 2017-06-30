package org.genvisis.gwas.mxc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

import org.genvisis.common.Aliases;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;

public class ParseMXCResults {
	private static String resultsGeneCol = "gene_name";
	private static String genesGeneCol = "GRCh37.p5-Primary Assembly_name";
	private static String geneStartCol = "GRCh37.p5-Primary Assembly_start";
	private static String geneEndCol = "GRCh37.p5-Primary Assembly_stop";
	private static String chrCol = "GRCh37.p5-Primary Assembly_chr";
	private static String genesFile = "/home/pankrat2/public/bin/NCBI/genes37.xln";

	private static String usage = "";

	private static String[] getKeys(String[][] matrix, String colname, Logger log) {
		String[] headers = matrix[0];
		int index = -1;
		for (int i = 0; i < headers.length; i++)
			if (colname.equals(headers[i]))
				index = i;


		if (index == -1) {
			log.report("Could not find column: " + colname);
			System.exit(0);
		}

		return Matrix.extractColumn(matrix, index);
	}


	public static String[][] loadGenes(String filename, boolean omitHeader, Logger log) {
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

		return HashVec.loadFileToStringMatrix(filename, omitHeader, cols, false);
	}

	private static String[][] removeNAs(String[][] matrix, int[] cols) {
		boolean[] rowsToKeep = new boolean[matrix.length];
		if (cols == null || cols.length == 0)
			return matrix;

		for (int c : cols) {
			for (int i = 0; i < matrix.length; i++) {
				rowsToKeep[i] = !(matrix[i][c].equals("NA") || matrix[i][c].equals("."));
			}
		}

		matrix = Matrix.subset(matrix, rowsToKeep);
		return matrix;
	}


	private static void generateRScript(String outputFile) {
		ArrayList<String> r = new ArrayList<String>();
		outputFile = new File(outputFile).getAbsolutePath();

		String filename = ext.rootOf(outputFile, false);

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

	public static void addMetalHits(String posfile, String mxcfile, String mapfile, String metalfile,
																	Logger log) {
		// read in the metal gwas results file
		String[][] metal = null;

		String[] line = Files.getHeaderOfFile(metalfile, log);

		int[] valueIndices = ext.indexFactors(new String[][] {Aliases.MARKER_NAMES, Aliases.PVALUES},
																					line, false,
																					true, false, false);

		metal = HashVec.loadFileToStringMatrix(metalfile, true, valueIndices, false);

		float bf = (float) 0.05 / metal.length;

		String[][] positions = HashVec.loadFileToStringMatrix(posfile, false, null, false);

		String[] keys = Matrix.extractColumn(metal, 0);
		String[][] results = null;

		// combine snps, pvals, and positions;
		try {
			results = Files.combineInMemory(keys, positions, "NA", true, true, log);
		} catch (Elision e) {
			log.reportError(e.getMessage());
		}

		// metal should now be of the form MarkerName, P-value, (marker), Chr, Position
		metal = ArrayUtils.append(metal, results);
		metal = removeNAs(metal, new int[] {3, 4});
		Arrays.sort(metal, new Comparator<String[]>() {
			@Override
			public int compare(final String[] entry1, final String[] entry2) {
				final String chr1 = entry1[3];
				final String chr2 = entry2[3];
				if (chr1.equals(chr2)) {
					final int pos1 = Integer.parseInt(entry1[4]);
					final int pos2 = Integer.parseInt(entry2[4]);
					return pos1 - pos2;
				}
				return chr1.compareTo(chr2);
			}
		});

		String[][] mxc = HashVec.loadFileToStringMatrix(mxcfile, false, null, false);

		// get position ranges for genes
		String[][] genePositions = loadGenes(genesFile, true, log);
		genePositions = removeNAs(genePositions, new int[] {1});

		Arrays.sort(genePositions, new Comparator<String[]>() {
			@Override
			public int compare(final String[] entry1, final String[] entry2) {
				final String chr1 = entry1[3];
				final String chr2 = entry2[3];
				if (chr1.equals(chr2)) {
					final int pos1 = Integer.parseInt(entry1[1]);
					final int pos2 = Integer.parseInt(entry2[1]);
					return pos1 - pos2;
				}
				return chr1.compareTo(chr2);
			}
		});

		// maps gene to related snps
		String[] mxcGenes = getKeys(mxc, resultsGeneCol, new Logger());
		int startPos;
		int endPos;
		int snpPos;
		double pval;
		String chr, snpChr;
		String[] snp, g;
		int start = 1;
		String[][] sig = new String[genePositions.length][];

		for (int gene = 0; gene < genePositions.length; gene++) {
			g = genePositions[gene];
			startPos = Integer.parseInt(g[1]) < 500000 ? 0 : Integer.parseInt(g[1]) - 500000;
			endPos = Integer.parseInt(g[2]) + 500000;
			chr = g[3];

			sig[gene] = new String[] {g[0], "0", "0"};

			while (start < metal.length && !metal[start][3].equals(chr)) {
				start++;
			}

			for (int j = start; j < metal.length; j++) {
				snp = metal[j];
				pval = Double.parseDouble(snp[1]);
				snpChr = snp[3];
				snpPos = Integer.parseInt(snp[4]);

				if (!snpChr.equals(chr))
					break;

				// genes are sorted by start position, so if we're not there yet, we can move it up
				if (snpPos < startPos) {
					start = j;
					continue;
				}

				if (snpPos > endPos)
					break;

				if (pval < bf * 100) {
					int[] numSig = new int[] {Integer.parseInt(sig[gene][1]), Integer.parseInt(sig[gene][2])};
					if (pval < bf)
						sig[gene][1] = numSig[0] + 1 + "";

					sig[gene][2] = numSig[1] + 1 + "";
				}
			}
		}


		// append num of sig snps to mxc file
		try {
			results = Files.combineInMemory(mxcGenes, sig, "NA", true, true, log);
			results[0] = new String[] {"gene", "numSig", "numSug"};
		} catch (Elision e) {
			log.reportError(e.getMessage());
		}

		mxc = ArrayUtils.append(mxc, results);
		// write to mxc file
		Files.writeMatrix(mxc, ext.addToRoot(mxcfile, "_sig"), ",");

		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(ext.addToRoot(mxcfile,
																																							"_no_sig_markers.csv")));
			writer.write(ArrayUtils.toStr(mxc[0], ","));
			for (String[] s : mxc) {
				if (s[s.length - 2].equals("0"))
					writer.write(ArrayUtils.toStr(s, ",") + "\n");
			}
			writer.close();
		} catch (IOException e) {
			log.reportError("Unable to write to " + ext.rootOf(mxcfile) + "_no_sig_markers.csv");
			e.printStackTrace();
		}

		generateRScript(ext.addToRoot(mxcfile, "_plot"));
	}

	public static void main(String[] args) throws IOException {
		if (args.length == 0) {
			System.out.println(usage);
			System.exit(1);
		}
		Logger log = new Logger();

		String mxcFile = "";
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

		addMetalHits(posmap, mxcFile, map, metal, log);

	}
}
