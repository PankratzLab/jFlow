package seq;

import java.io.*;
import java.util.*;
import common.*;

public class RankRecessive {
	// you can change the name of these variables, but don't change the order without editing how the Hashtable is generated  
//	public static final String[] REQS = {"Position", "Gene", "Ref_Allele", "GT of Siblings", "rs", "merged_MAF"};
	public static final String[] REQS = {"Position", "Gene", "Ref_Allele", "GT_of_Siblings", "rs", "merged_MAF"};
	
	public static void rank(String filename, double mafForMissingValues, double mafFloor) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, header, genes;
		Hashtable<String, Vector<String[]>> hash;
		Vector<String[]> variants;
		int count;
		int[] indices, order;
		double[] mafs;
		String[] alleles;
		
		hash = new Hashtable<String, Vector<String[]>>();
		try {
			reader = new BufferedReader(new FileReader(filename));
			header = reader.readLine().trim().split("\\t", -1);
			indices = ext.indexFactors(REQS, header, false, true);
			while (reader.ready()) {
				line = reader.readLine().trim().split("\\t", -1);
				HashVec.addToHashArrayVec(hash, line[1], new String[] {line[indices[0]], line[indices[4]], line[indices[5]], line[indices[2]], line[indices[3]]});
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
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_ranked.xln"));
			genes = HashVec.getKeys(hash);
			writer.println("Gene\tProduct\tPosition1\trs1\tVariant1\tRef1\tMAF1\tPosition2\trs2\tVariant2\tRef2\tMAF2\tNumberOfOtherAllelesInGene");
			for (int i = 0; i < genes.length; i++) {
				variants = hash.get(genes[i]);
				count = variants.size();
				if (count >= 2) {
					for (int j = 0; j < count; j++) {
						line = variants.elementAt(j);
						alleles = line[4].split("/");
						if (!line[3].equals(alleles[0]) && !line[3].equals(alleles[0])) {
							variants.add(line);
						}
					}
					mafs = new double[variants.size()];
					for (int j = 0; j < variants.size(); j++) {
						line = variants.elementAt(j);
						if (line[2].equals(".")) {
							mafs[j] = mafForMissingValues; 
						} else {
							try {
								mafs[j] = Double.parseDouble(line[2]);
							} catch (NumberFormatException nfe) {
								System.err.println("Error parsing MAF ('"+line[2]+"') for "+line[0]+" in gene "+genes[i]+"; setting to missing value");
								mafs[j] = mafForMissingValues;
							}
							if (mafs[j] < mafFloor) {
								mafs[j] = mafFloor;
							}
						}
					}
					order = Sort.quicksort(mafs);
					writer.print(genes[i]+"\t"+ext.formDeci(mafs[order[0]]*mafs[order[1]], 5));
					for (int j = 0; j < 2; j++) {
						line = variants.elementAt(order[j]);
						writer.print("\t"+line[0]+"\t"+line[1]+"\t"+line[3]+"\t"+line[4].substring(line[4].indexOf("/")+1)+"\t"+line[2]);
					}
					writer.println("\t"+(variants.size()==2?".":(variants.size()-2)));
				}
								
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.rootOf(filename, false)+"_ranked.xln");
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		double mafForMissingValues = 0.01;
		double mafFloor = 0.01;

		String usage = "\n" + 
		"seq.RankRecessive requires 1-3 arguments\n" + 
		"   (1) filename (i.e. file=" + filename + " (required; not the default))\n" + 
		"   (2) MAF for any position with a missing value (i.e. mafForMissingValues=" + mafForMissingValues + " (default))\n" + 
		"   (3) minimum MAF, all values below this will be set to this value (i.e. mafFloor=" + mafFloor + " (default))\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("mafForMissingValues=")) {
				mafForMissingValues = Double.parseDouble(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("mafFloor=")) {
				mafFloor = Double.parseDouble(args[i].split("=")[1]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		
//		filename = "D:\\tWork\\FIA\\family1303_5_recessive_gene_centric_MAF_2.txt.xln";
//		filename = "D:\\tWork\\FIA\\family1003_4_recessive_gene_centric_MAF_2.txt.xln";
		
		try {
			if (filename != null) {
				rank(filename, mafForMissingValues, mafFloor);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
