// -Xms1024M -Xmx1024M
package org.genvisis.dead;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.ext;

public class SummaryFor_dbGaP {
	public static final String HWE_FILE = "plink.hwe";
	public static final String CHIP_BATCH_NAME = "IlluminaHumanCNV370Duo";
	public static final String[] FINAL_HEADER = {	"Chip_batch_name",
																								"SNP_ID  # original ID or ss#  used by the chip vender. ",
																								"Original Allele 1  (probe allele)",
																								"Original Allele 2  (probe allele)",
																								"Genotype 1-1 Case", "Genotype 1-2 Case",
																								"Genotype 2-2 Case", "Genotype 1-1 Control",
																								"Genotype 1-2 Control", "Genotype 2-2 Control",
																								"HWE Cases", "HWE Controls",
																								"Odds ratio (for allele 1)", "Test Statistic",
																								"P value"};

	public static void summarize(String dir, String results, String test) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, data;
		Hashtable<String, String[]> hash = new Hashtable<String, String[]>();
		int snpCol, orCol, pCol, statCol, testCol, a1Col, a2Col, genoCol;
		String caseP;

		try {
			reader = new BufferedReader(new FileReader(dir + results));
			line = reader.readLine().trim().split("[\\s]+");
			snpCol = ext.indexFactors(new String[] {"SNP"}, line, false, true)[0];
			orCol = ext.indexFactors(new String[] {"OR"}, line, true, false)[0];
			statCol = ext.indexFactors(new String[] {"STAT"}, line, true, true)[0];
			pCol = ext.indexFactors(new String[] {"P"}, line, true, true)[0];
			testCol = ext.indexFactors(new String[] {"TEST"}, line, true, true)[0];
			System.out.println(ext.getTime() + "\tReading in p-values for " + test + " in " + results);
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (line[testCol].equalsIgnoreCase(test)) {
					hash.put(line[snpCol], new String[] {line[orCol], line[statCol], line[pCol]});
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + dir + results + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dir + results + "\"");
			System.exit(2);
		}

		try {
			reader = new BufferedReader(new FileReader(dir + HWE_FILE));
			writer = new PrintWriter(new FileWriter(dir + "summary.xln"));
			writer.println(ArrayUtils.toStr(FINAL_HEADER));
			line = reader.readLine().trim().split("[\\s]+");
			snpCol = ext.indexFactors(new String[] {"SNP"}, line, false, true)[0];
			a1Col = ext.indexFactors(new String[] {"A1"}, line, true, false)[0];
			a2Col = ext.indexFactors(new String[] {"A2"}, line, true, true)[0];
			genoCol = ext.indexFactors(new String[] {"GENO"}, line, true, true)[0];
			pCol = ext.indexFactors(new String[] {"P"}, line, true, true)[0];
			testCol = ext.indexFactors(new String[] {"TEST"}, line, true, true)[0];
			System.out.println(ext.getTime() + "\tReading in counts and p-values from " + HWE_FILE);
			while (reader.ready()) {
				reader.readLine();

				line = reader.readLine().trim().split("[\\s]+");
				if (!line[testCol].equalsIgnoreCase("aff")) {
					System.err.println("Error - out of sync at " + line[snpCol]);
					System.exit(1);
				}
				caseP = line[pCol];
				writer.print(CHIP_BATCH_NAME	+ "\t" + line[snpCol] + "\t" + line[a1Col] + "\t" + line[a2Col]
											+ "\t" + ArrayUtils.toStr(line[genoCol].split("/")));

				line = reader.readLine().trim().split("[\\s]+");
				if (!line[testCol].equalsIgnoreCase("unaff")) {
					System.err.println("Error - out of sync at " + line[snpCol]);
					System.exit(1);
				}
				writer.print("\t"	+ ArrayUtils.toStr(line[genoCol].split("/")) + "\t" + line[pCol] + "\t"
											+ caseP);

				data = hash.get(line[snpCol]);
				if (data == null) {
					System.err.println("Error - no p-value for " + line[snpCol]);
					data = new String[] {"NA", "NA", "NA"};
				}
				writer.println("\t" + ArrayUtils.toStr(data));

			}
			writer.close();
			reader.close();
			System.out.println(ext.getTime() + "\tDone!");
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + dir + results + "\" not found in current directory");

			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dir + results + "\"");
			System.exit(2);
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		// String dir = "";
		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\dbGaP\\";
		// String results = "plink.assoc.logistic";
		String results = "plink.assoc.recessive.logistic";
		String test = "REC";

		String usage = "\n"	+ "gwas.SummaryFor_dbGaP requires 0-1 arguments\n"
										+ "   (1) directory (i.e. dir=results/ (not the default))\n"
										+ "   (2) results file (i.e. dir=" + results + " (default))\n"
										+ "   (3) test (i.e. test=" + test + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("dir=")) {
				dir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("results=")) {
				results = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("test=")) {
				test = arg.split("=")[1];
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			summarize(dir, results, test);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
