package gaw17;

import java.io.*;
import java.util.*;
import common.*;

public class Parser {
	public static final String[] ALLELES = {"A", "C", "G", "T"};

	public static void parse(String geno_dir, String snpFile, String sexFile, String suffix) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String trav;
		Hashtable<String, String> snpInfo, sexInfo;
		
		snpInfo = HashVec.loadFileToHashString(snpFile, new int[] {0}, new int[] {1,2}, true, "\t", true, false, false);
		sexInfo = HashVec.loadFileToHashString(sexFile, new int[] {0}, new int[] {1}, true, "\t", true, false, false);
		
		for (int chr = 1; chr <= 22; chr++) {
			try {
				reader = new BufferedReader(new FileReader(geno_dir+"c"+chr+"_snps"+suffix));
				line = reader.readLine().trim().split(",");
				try {
					writer = new PrintWriter(new FileWriter(geno_dir+"chr"+chr+".map"));
					for (int i = 1; i < line.length; i++) {
						trav = snpInfo.get(line[i]);
						if (trav == null) {
							System.err.println("Error - no map info for marker: "+line[i]);
						}
						if (!trav.startsWith(chr+"\t")) {
							System.err.println("Error - marker is on the wrong chromosome: "+line[i]);
						}
						writer.println(chr+"\t"+line[i]+"\t0\t"+trav.split("[\\s]+")[1]);
					}
					writer.close();
				} catch (Exception e) {
					System.err.println("Error writing to " + geno_dir+"chr"+chr+".map");
					e.printStackTrace();
				}
				writer = new PrintWriter(new FileWriter(geno_dir+"chr"+chr+".ped"));
				while (reader.ready()) {
					line = reader.readLine().trim().split(",");
					trav = sexInfo.get(line[0]);
					if (trav == null) {
						System.err.println("Error - no sex info for individual: "+line[0]);
					}
					writer.print(line[0]+"\t"+line[0]+"\t0\t0\t"+trav+"\t1");
					for (int i = 1; i < line.length; i++) {
						if (line[i].length() != 3 || ext.indexOfStr(line[i].substring(0, 1), ALLELES) == -1 || ext.indexOfStr(line[i].substring(2, 3), ALLELES) == -1 || !line[i].substring(1, 2).equals("/")) {
							System.err.println("Error - invalid genotype: "+line[i]);
						}
						writer.print("\t"+line[i].charAt(0)+"\t"+line[i].charAt(2));
					}
					writer.println();
				}
				writer.close();
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + geno_dir+"" + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + geno_dir+"" + "\"");
				System.exit(2);
			}
		}
		
		try {
			writer = new PrintWriter(new FileWriter(geno_dir+"merge.txt"));
			for (int chr = 2; chr <= 22; chr++) {
				writer.println("chr"+chr+".ped chr"+chr+".map");
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + geno_dir+"merge.txt");
			e.printStackTrace();
		}

		try {
			writer = new PrintWriter(new FileWriter(geno_dir+"merge.bat"));
			writer.println("plink --file chr1 --merge-list merge.txt --make-bed");
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + geno_dir+"merge.txt");
			e.printStackTrace();
		}
	}
	
	public static void parseMarkerSets(String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Hashtable<String, Vector<String>> hash, nonsynon;
		Vector<String> genes;
		
		genes = new Vector<String>();
		hash = new Hashtable<String, Vector<String>>();
		nonsynon = new Hashtable<String, Vector<String>>();
		try {
			reader = new BufferedReader(new FileReader(filename));
			line = reader.readLine().trim().split(",");
			ext.checkHeader(line, new String[] {"snp_name", "chromosome", "position", "gene", "snp_type", "variant", "maf"}, true);
			while (reader.ready()) {
				line = reader.readLine().trim().split(",");
				HashVec.addToHashVec(hash, line[3], line[0], false);
				if (line[4].equals("Nonsynonymous")) {
					HashVec.addToHashVec(nonsynon, line[3], line[0], false);
				}
				HashVec.addIfAbsent(line[3], genes);
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
			writer = new PrintWriter(new FileWriter(ext.parseDirectoryOfFile(filename)+"markerSets.dat"));
			for (int i = 0; i < genes.size(); i++) {
				writer.println(genes.elementAt(i)+"\t"+Array.toStr(Array.toStringArray(hash.get(genes.elementAt(i)))));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.parseDirectoryOfFile(filename));
			e.printStackTrace();
		}
		
		try {
			writer = new PrintWriter(new FileWriter(ext.parseDirectoryOfFile(filename)+"nonynonMarkerSets.dat"));
			for (int i = 0; i < genes.size(); i++) {
				if (nonsynon.containsKey(genes.elementAt(i))) {
					writer.println(genes.elementAt(i)+"\t"+Array.toStr(Array.toStringArray(nonsynon.get(genes.elementAt(i)))));
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.parseDirectoryOfFile(filename));
			e.printStackTrace();
		}
		
	}
	
	public static void main(String[] args) {
		String snpFile = "D:\\GAW17\\source\\snp_info.csv";
//		String dir = "D:\\GAW17\\source\\genos_unr\\";
//		String sexFile = "D:\\GAW17\\source\\unrelateds.ped";
		String dir = "D:\\GAW17\\source\\genos_fam\\";
		String sexFile = "D:\\GAW17\\source\\families_sex.ped";
		String suffix = ".fam";
		
		try {
			parse(dir, snpFile, sexFile, suffix);
//			parseMarkerSets(snpFile);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
