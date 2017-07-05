package org.genvisis.gwas.mxc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

import org.genvisis.cnv.plots.ManhattanPlot;
import org.genvisis.common.Aliases;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;
import org.genvisis.gwas.AlleleVerification;

public class ParseMXCResults {
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

		String[] line = Files.getHeaderOfFile(filename, log);
		int[] cols = ext.indexFactors(new String[][] {new String[] {"Assembly_name", "Gene",
																																"SKATgene"},
																									Aliases.POSITIONS_START,
																									Aliases.POSITIONS_STOP, Aliases.CHRS},
																	line, true, false, false, false, log, false);


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

	public static void addMetalHits(String posfile, String mxcfile, String metalfile,
																	String genesFile,
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
		String[] mxcGenes = getKeys(mxc, "gene_name", new Logger());
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

			sig[gene] = new String[] {g[0], "0", "0", g[1], g[3]};

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
			results[0] = new String[] {"gene", "numSig", "numSug", "pos", "chr"};
		} catch (Elision e) {
			log.reportError(e.getMessage());
		}

		mxc = ArrayUtils.append(mxc, results);
		// write to mxc file
		Files.writeMatrix(mxc, ext.addToRoot(mxcfile, "_sig"), ",");

		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(ext.addToRoot(mxcfile,
																																							"_no_sig_markers")));
			writer.write(ArrayUtils.toStr(mxc[0], ",") + "\n");
			for (String[] s : mxc) {
				if (s[s.length - 4].equals("0"))
					writer.write(ArrayUtils.toStr(s, ",") + "\n");
			}
			writer.close();
		} catch (IOException e) {
			log.reportError("Unable to write to " + ext.rootOf(mxcfile) + "_no_sig_markers.csv");
			e.printStackTrace();
		}

		plot(ext.addToRoot(mxcfile, "_no_sig_markers"), ext.rootOf(mxcfile, false) + "_manhattan.png",
				 log);

		generateRScript(ext.addToRoot(mxcfile, "_plot"));
	}

	// TODO: make this write the commands to a a file instead?
	// Actually this should be in parse s/t main here can make the script executor
	// and parse's main can call run
	private static void run(String data, String db, String posmap, String covar, String out,
													boolean overwrite,
													String mxcFolder, Logger log) {
		String gwas_folder = ext.parseDirectoryOfFile(new File(data).getAbsolutePath());

		db = new File(db).getAbsolutePath();
		covar = new File(covar).getAbsolutePath();

		out = new File(out).getAbsolutePath();

		String py = "MetaXcan.py";

		String delimiter = Files.determineDelimiter(data, log);
		String[] header = Files.getHeaderOfFile(data, delimiter, log);
		int[] index = ext.indexFactors(new String[][] {Aliases.MARKER_NAMES, Aliases.REF_ALLELES,
																									 Aliases.ALT_ALLELES, Aliases.EFFECTS,
																									 Aliases.STD_ERRS, Aliases.PVALUES},
																	 header, false, true, false, false);

		if (index[0] == -1) {
			log.reportError("Unable to find SNP column. Aborting.");
			System.exit(1);
		} else if (index[1] == -1 || index[2] == -1) {
			log.reportError("Unable to find allele columns. Aborting.");
			System.exit(1);
		} else if ((index[5] == -1 && index[4] == -1) || index[3] == -1) {
			log.reportError("Unable to find necessary pval or effect columns. Beta column and either P-value or SE columns are required. Aborting.");
			System.exit(1);
		}

		String markerColumn = header[index[0]];
		String a1 = header[index[1]];
		String a2 = header[index[2]];
		String effect = " --beta_column " + header[index[3]]
										+ (index[5] == -1 ? " --se_column " + header[index[4]]
																			: " --pvalue_column " + header[index[5]]);

		// build mxc command
		String command = "./" + py
										 + " --model_db_path " + db + " --covariance " + covar + " --gwas_folder "
										 + gwas_folder + " --gwas_file_pattern " + new File(data).getName()
										 + " --output_file " + out + " --effect_allele_column " + a1
										 + " --non_effect_allele_column " + a2 + " --snp_column " + markerColumn
										 + effect
										 + (overwrite ? " --overwrite" : "");

		System.out.println(command);

		// run MetaXcan on the given inputs
		boolean runSuccess = CmdLine.run(command, ext.parseDirectoryOfFile(mxcFolder), null, null,
																		 new Logger("MetaXcan.log"), false);

		if (!runSuccess || !new File(out).exists()) {
			System.out.println("Error running MetaXcan with the given inputs.");
			System.exit(0);
		}
	}

	public static void plot(String filename, String out, Logger log) {
		ManhattanPlot mp = new ManhattanPlot(null);

		try {
			boolean load = mp.loadFileAuto(filename);
			if (!load) {
				log.reportError("Unable to load " + filename + " for manhattan plot.");
				return;
			}
			while (!mp.isDataLoaded()) {
				Thread.sleep(200);
			}
		} catch (InterruptedException e) {
			log.reportError("Problem creating manhattan plot from " + filename);
			e.printStackTrace();
		}

		mp.getManPan().setSize(800, 400);
		mp.screenshot(out);
	}

	public static void main(String[] args) throws IOException {
		if (args.length == 0) {
			System.out.println(usage);
			System.exit(1);
		}
		Logger log = new Logger();

		String data = null;
		String posmap = null;
		String out = "";
		String db = "DGN-HapMap-2015/DGN-WB_0.5.db";
		String covar = "covariance.DGN-WB_0.5.txt.gz";
		String mxcFolder = "MetaXcan/software";
		String freqFile = "freq.tbl";
		String refFile = "1000G.xln";
		String genesFile = "/home/pankrat2/public/bin/NCBI/genes37.xln";

		boolean verify = false;
		boolean overwrite = false;

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("help")) {
				System.out.println(usage);
				System.exit(1);
			} else if (arg.startsWith("genes="))
				genesFile = ext.parseStringArg(arg);
			else if (arg.startsWith("data="))
				data = ext.parseStringArg(arg);
			else if (arg.startsWith("posmap="))
				posmap = ext.parseStringArg(arg);
			else if (arg.startsWith("out="))
				out = ext.parseStringArg(arg);
			else if (arg.startsWith("db="))
				db = ext.parseStringArg(arg);
			else if (arg.startsWith("covar="))
				covar = ext.parseStringArg(arg);
			else if (arg.startsWith("mxc="))
				mxcFolder = ext.parseStringArg(arg);
			else if (arg.startsWith("freq="))
				freqFile = ext.parseStringArg(arg);
			else if (arg.startsWith("ref="))
				refFile = ext.parseStringArg(arg);
			else if (arg.startsWith("-verify"))
				verify = true;
			else if (arg.startsWith("-overwrite"))
				overwrite = true;
			else
				log.report("Invalid argument: " + arg);
		}

		if (verify) {
			AlleleVerification.verifyAlleles(data, refFile, freqFile, posmap, false, log);
			data = ext.rootOf(data, false) + "_allele_verified.txt";
		}

		// run(data, db, posmap, covar, out, overwrite, mxcFolder, log);

		// take the output mxc file and find the number of hits for each gene range
		addMetalHits(posmap, out, data, genesFile, new Logger("parseMXC.log"));

	}
}
