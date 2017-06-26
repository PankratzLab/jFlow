package org.genvisis.gwas.mxc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

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


	private static String[][] loadGenes(String filename, boolean omitHeader, Logger log) {
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

	private static String[][] removeNAs(String[][] matrix, int col) {
		boolean[] rowsToKeep = new boolean[matrix.length];
		if (col < 0)
			return matrix;

		for (int i = 0; i < matrix.length; i++) {
			rowsToKeep[i] = !(matrix[i][col].equals("NA") || matrix[i][col].equals("."));
		}

		matrix = Matrix.subset(matrix, rowsToKeep);
		return matrix;
	}

	public static void addMetalHits(String posfile, String mxcfile, String mapfile, String metalfile,
																	Logger log) {
		// read in the metal gwas results file
		String[][] metal = null;
		try {
			BufferedReader b = new BufferedReader(new FileReader(new File(metalfile)));
			String temp = b.readLine();
			String[] line = null;
			if (temp != null) {
				line = temp.split("\t");
			} else {
				log.reportError("Could not find MarkerName in header for " + metalfile);
				System.exit(1);
			}

			int[] valueIndices = ext.indexFactors(new String[] {"MarkerName", "P-value"}, line, false,
																						true);

			metal = HashVec.loadFileToStringMatrix(metalfile, true, valueIndices, false);
			b.close();
		} catch (Exception e) {
			log.reportError("Unable to load gwas summary file. Aborting.");
			System.exit(1);
		}
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
		metal = append(metal, results);

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
		genePositions = removeNAs(genePositions, 1);

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
		ArrayList<String> snpsOnGene = new ArrayList<String>();
		int startPos;
		int endPos;
		int snpPos;
		double pval;
		String chr, snpChr;
		String[] snp, g;
		int start = 1;
		String[][] sig = new String[genePositions.length][];

		snpsOnGene.add("gene,snp,chr,pos,p,sig");

		for (int gene = 0; gene < genePositions.length; gene++) {
			g = genePositions[gene];
			startPos = Integer.parseInt(g[1]) < 500000 ? 0 : Integer.parseInt(g[1]) - 500000;
			endPos = Integer.parseInt(g[2]) + 500000;
			chr = g[3];

			sig[gene] = new String[] {g[0], "0", "0"};

			while (start < metal.length && !metal[start][3].equals(chr)) {
				start++;
			}

			if (g[0].equals("SNCA")) {
				System.out.println("Starting SNCA, with range " + startPos + " : " + endPos);
				System.out.println("Starting search at position " + metal[start][4]);
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
					if (pval < bf) {
						snpsOnGene.add(g[0] + "," + snp[0] + "," + chr + "," + snpPos + "," + pval + ","
													 + "***");
						sig[gene][1] = numSig[0] + 1 + "";
					} else {
						snpsOnGene.add(g[0] + "," + snp[0] + "," + chr + "," + snpPos + "," + pval + "," + "*");
					}
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

		mxc = append(mxc, results);
		// write to mxc file
		Files.writeMatrix(mxc, ext.addToRoot(mxcfile, "_sig"), ",");
		Files.writeIterable(snpsOnGene, ext.parseDirectoryOfFile(mxcfile) + "sigSNPs.csv");
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
