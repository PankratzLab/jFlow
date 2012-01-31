package bioinformatics;

import java.io.*;
import java.util.*;

import common.*;

public class procTagsToGeneIDs {
	public static final String[] FILES_WITH_GENE = {"d1_gene_info.prn", "nbt1239-S4.prn"};
	public static final String[] FILES_WITH_GENE_ID = {"gnf1b-anntable_geneID_tag.prn", "lookup.xls.prn"};
	public static final String TAGS = "probset.txt";
	public static final String OVERRIDE = "override.prn";
	public static final String DEFAULT_Q = "C:\\Download\\seq_gene.q";

	public static void runProcTagsToGeneIDs() {
		BufferedReader reader;
		PrintWriter writer, geneids;
		String[] line, files, ids;
		String temp, trav, geneid;
		// Hashtable override = new Hashtable();
		Hashtable<String,Hashtable<String,String>> hash = new Hashtable<String,Hashtable<String,String>>();
		Hashtable<String,String> h2;
		Hashtable<String,Vector<String>> aliases = new Hashtable<String,Vector<String>>();
		Vector<String> v;

		// if (new File(OVERRIDE).exists()) {
		// override = HashVec.loadFileToHashString(OVERRIDE, false);
		// }

		System.out.println("Loading aliases...");
		try {
			reader = new BufferedReader(new FileReader(DEFAULT_Q));
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				if (line[0].startsWith("GeneID:")) {
					temp = line[0].substring(7);
					HashVec.addToHashVec(aliases, line[2], temp, true);
					line = line[5].split("\\|");
					for (int i = 0; i<line.length; i++) {
						HashVec.addToHashVec(aliases, line[i], temp, true);
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+DEFAULT_Q+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+DEFAULT_Q+"\"");
			System.exit(2);
		}

		files = new String[FILES_WITH_GENE.length+FILES_WITH_GENE_ID.length];
		for (int i = 0; i<FILES_WITH_GENE.length; i++) {
			files[i] = FILES_WITH_GENE[i];
			try {
				System.out.println("Loading "+FILES_WITH_GENE[i]+"...");
				reader = new BufferedReader(new FileReader(FILES_WITH_GENE[i]));
				writer = new PrintWriter(new FileWriter("noDirectMatch.xls"));
				reader.readLine();
				while (reader.ready()) {
					line = reader.readLine().split("[\\s]+");
					if (!line[1].equals("---")) {
						if (aliases.containsKey(line[1])) {
							v = aliases.get(line[1]);
							if (v.size()>1) {
								writer.println(line[0]+"\t"+FILES_WITH_GENE[i]+"\t"+Array.toStr(Array.toStringArray(v)));
							} else {
								HashVec.addToHashHash(hash, line[0], FILES_WITH_GENE[i], v.elementAt(0));
							}
						} else {
							writer.println(line[0]+"\t"+FILES_WITH_GENE[i]+"\tmissing\t"+line[1]);
						}
					}
				}
				reader.close();
				writer.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \""+FILES_WITH_GENE[i]+"\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \""+FILES_WITH_GENE[i]+"\"");
				System.exit(2);
			}
		}

		for (int i = 0; i<FILES_WITH_GENE_ID.length; i++) {
			files[FILES_WITH_GENE.length+i] = FILES_WITH_GENE_ID[i];
			try {
				System.out.println("Loading "+FILES_WITH_GENE_ID[i]+"...");
				reader = new BufferedReader(new FileReader(FILES_WITH_GENE_ID[i]));
				reader.readLine();
				while (reader.ready()) {
					line = reader.readLine().split("[\\s]+");
					if (!line[1].equals("---")) {
						HashVec.addToHashHash(hash, line[0], FILES_WITH_GENE_ID[i], line[1]);
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \""+FILES_WITH_GENE_ID[i]+"\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \""+FILES_WITH_GENE_ID[i]+"\"");
				System.exit(2);
			}
		}

		System.out.println("Merging with "+TAGS+"...");
		try {
			reader = new BufferedReader(new FileReader(TAGS));
			writer = new PrintWriter(new FileWriter("tagsToGeneIDs.xls"));
			geneids = new PrintWriter(new FileWriter("tagsToGeneIDs.prn"));
			writer.println("Tag\t"+Array.toStr(files));
			while (reader.ready()) {
				trav = reader.readLine();
				writer.print(trav);
				if (hash.containsKey(trav)) {
					h2 = hash.get(trav);
					ids = HashVec.getKeys(h2);
					geneid = "";
					temp = "";
					for (int i = 0; i<files.length; i++) {
						temp += "\t";
						for (int j = 0; j<ids.length; j++) {
							if (files[i].equals(ids[j])) {
								temp += h2.get(ids[j]);
								if (geneid.equals("")) {
									geneid = h2.get(ids[j]);
								} else if (!geneid.equals(h2.get(ids[j]))) {
									geneid = "-999";
								}
							}
						}
					}
					if (geneid.equals("-999")) {
						writer.print(temp+"\t\t1");
					} else {
						writer.print(temp+"\t"+geneid);
						geneids.println(trav+"\t"+geneid);
					}
				}
				writer.println();
			}
			writer.close();
			geneids.close();
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+TAGS+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+TAGS+"\"");
			System.exit(2);
		}

	}

	public static void main(String[] args) throws IOException {
		try {
			runProcTagsToGeneIDs();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
