// expecting a HapMart/BioMart export that includes position, marker_id, and reference_allele
// pedinfo2sample_CEU.txt is available at http://www.hapmap.org/downloads/samples_individuals/
package org.genvisis.dead;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Files;
import org.genvisis.common.ext;

public class HapMapParserOld {
	public static boolean USE_NUCLEOTIDES = true;
	public static boolean INCLUDE_MONOMORPHIC = true;

	public static void parse(String dir, String filename, String famstruct) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, indIDs = null, trans;
		String temp, trav;
		Hashtable<String, String[]> hash = new Hashtable<String, String[]>();
		Vector<String> v = new Vector<String>();
		int index, ones, twos;
		String root = ext.rootOf(filename);

		try {
			reader = Files.getReader(filename, dir);
			temp = reader.readLine();
			indIDs = (temp.substring(temp.indexOf("[") + 1, temp.indexOf("]"))).split("[\\s]+");

			writer = new PrintWriter(new FileWriter((new File(dir).exists() ? dir : "")	+ root
																							+ ".info"));
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				trans = new String[line.length - 3];
				ones = twos = 0;
				for (int i = 3; i < line.length; i++) {
					trans[i - 3] = "";
					for (int j = 0; j < 2; j++) {
						if (line[i].charAt(j) == 'N') {
							trans[i - 3] += "0";
						} else if (line[i].charAt(j) == line[2].charAt(0)) {
							if (USE_NUCLEOTIDES) {
								trans[i - 3] += line[i].charAt(j);
							} else {
								trans[i - 3] += "2";
							}
							twos++;
						} else {
							if (USE_NUCLEOTIDES) {
								trans[i - 3] += line[i].charAt(j);
							} else {
								trans[i - 3] += "1";
							}
							ones++;
						}
						trans[i - 3] += j == 0 ? "\t" : "";
					}
				}
				if (INCLUDE_MONOMORPHIC || (ones > 0 && twos > 0)) {
					writer.println(line[1] + "\t" + line[0]);
					v.add(line[1]);
					hash.put(line[1], trans);
				}
			}
			reader.close();
			writer.close();
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			ioe.printStackTrace();
			System.exit(2);
		}

		try {
			reader = Files.getReader(	famstruct,
																new String[] {dir, "/home/npankrat/NCBI/",
																							"C:\\Documents and Settings\\npankrat\\My Documents\\jProjects\\park\\runtime\\"});
			writer = new PrintWriter(new FileWriter((new File(dir).exists() ? dir : "") + root + ".pre"));
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				for (int i = 0; i < 5; i++) {
					writer.print(line[i] + "\t");
				}
				writer.print("1");
				index = -2;
				trav = line[6];
				if (famstruct.contains("pedinfo2sample_")) {
					line = line[6].split(":");
					if (line.length != 6) {
						System.err.println("Error - different format than expected for pedinfo2sample_***.txt file (do not alter from what's posted)");
						System.exit(1);
					}
					trav = line[4];
				}
				for (int i = 0; i < indIDs.length; i++) {
					if (indIDs[i].equals(trav)) {
						index = i;
					}
				}
				if (index == -2) {
					System.err.println("Error - Could not find sample "	+ trav
															+ " from the pedigree file in the genotype file");
				}
				for (int i = 0; i < v.size(); i++) {
					line = hash.get(v.elementAt(i));
					writer.print("\t" + line[index]);
				}
				writer.println();
			}
			reader.close();
			writer.close();
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + famstruct + "\"");
			ioe.printStackTrace();
			System.exit(2);
		}

	}

	public static void plinkMapToHaploviewInfo(String from, String to) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;

		try {
			reader = new BufferedReader(new FileReader(from));
			writer = new PrintWriter(new FileWriter(to));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				writer.println(line[1] + "\t" + line[3]);
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + from + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + from + "\"");
			System.exit(2);
		}
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		// String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\SNCA\\";
		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\followup\\APOJ\\";
		// String filename = "2qCEU-GENOS.tsv";
		// String filename = "39rGV1lKEW5.tsv";
		// String filename = "kept39rGV1lKEW2.tsv";
		// String filename = "SNCA_50kb.tsv";
		// String filename = "SNCA_50kb_43tags.tsv";
		// String filename = "SNCA_50kb_21tags.tsv";
		// String filename = "mapt.tsv";
		// String filename = "CTNNA2.tsv";
		// String filename = "deleteriousTau.tsv";
		// String filename = "GAKimputed.tsv";
		// String filename = "maptUTRsigs.tsv";
		// String filename = "gakExonic.tsv";
		// String filename = "maptUTR2.tsv";
		// String filename = "LAMP1.tsv";
		String filename = "APOJ.tsv";


		// String famstruct = "famstruct_CEU.dat";
		String famstruct = "/home/npankrat/NCBI/pedinfo2sample_CEU.txt";

		String usage = "\n"	+ "CommonTools.bioinformatics.hapmapParser requires 0-1 arguments\n"
										+ "   (1) filename (i.e. file=" + filename + " (default)\n"
										+ "   (2) famstruct (i.e. struct=" + famstruct + " (default)\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("struct=")) {
				famstruct = arg.split("=")[1];
				numArgs--;
			}
		}
		if (filename.startsWith("C:")) {
			dir = "";
		}
		for (int i = 0; i < args.length; i++) {
			System.out.println((i + 1) + ")" + args[i]);
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			parse(dir, filename, famstruct);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
