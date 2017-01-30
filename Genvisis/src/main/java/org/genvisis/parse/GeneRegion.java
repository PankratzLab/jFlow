// used to see how many genes/snps are within each of the 235 regions on the Metabochip
package org.genvisis.parse;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CountVector;
import org.genvisis.common.ext;

public class GeneRegion {
	public static final String[] HEADER = {"SNP", "Chr", "Position", "region", "Gene(s)"};

	public static void parse(String filename) {
		BufferedReader reader;
		PrintWriter writer, writer2;
		String[] line;
		Hashtable<String, CountVector> countVectors;
		CountVector cv;
		int[][] boundaries;
		String[] genes;
		int region;
		int[] counts;

		boundaries = new int[235][3];
		countVectors = new Hashtable<String, CountVector>();

		try {
			reader = new BufferedReader(new FileReader(filename));
			line = reader.readLine().trim().split("[\\s]+");
			ext.checkHeader(line, HEADER, true);
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				genes = line[4].split("\\|");
				if (genes[0].equals(".")) {
					genes[0] = "intergenic";
				}
				region = line[3].equals("NA") ? 0 : Integer.parseInt(line[3]);
				if (countVectors.containsKey(region + "")) {
					cv = countVectors.get(region + "");
				} else {
					countVectors.put(region + "", cv = new CountVector());
					cv.add("intergenic");
				}

				for (String gene : genes) {
					cv.add(gene);
				}

				if (boundaries[region][0] == 0) {
					boundaries[region][0] = Integer.parseInt(line[1]);
					boundaries[region][1] = Integer.parseInt(line[2]);
					boundaries[region][2] = Integer.parseInt(line[2]);
				} else if (boundaries[region][1] > Integer.parseInt(line[2])) {
					boundaries[region][1] = Integer.parseInt(line[2]);
				} else if (boundaries[region][2] < Integer.parseInt(line[2])) {
					boundaries[region][2] = Integer.parseInt(line[2]);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}

		try {
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false) + "_genes.xln"));
			writer.println("Region\tGene\tNumSNPsWithin");
			writer2 = new PrintWriter(new FileWriter(ext.rootOf(filename, false) + "_regions.xln"));
			writer2.println("Region\tChr\tStartPosition\tStopPosition\tNumber of Genes");
			for (int i = 0; i < boundaries.length; i++) {
				cv = countVectors.get(i + "");
				genes = cv.getValues();
				counts = cv.getCounts();
				writer2.println((i == 0 ? "NA" : i)	+ "\t" + boundaries[i][0] + "\t" + boundaries[i][1]
												+ "\t" + boundaries[i][2] + "\t" + cv.getSize() + "\t"
												+ (i == 0 ? "too many to list" : ArrayUtils.toStr(genes)));
				for (int j = 0; j < genes.length; j++) {
					writer.println(i + "\t" + genes[j] + "\t" + (j == 0 ? counts[j] - 1 : counts[j]));
				}
			}
			writer.close();
			writer2.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.rootOf(filename, false) + "_genes.xln");
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "GeneRegion.dat";

		String usage = "\n"	+ "parse.GeneRegion requires 0-1 arguments\n"
										+ "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		filename = "D:\\Myron\\CALICO\\Metabochip\\00src\\snp_position_region_gene.txt";
		try {
			parse(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
