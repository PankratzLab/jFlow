package org.genvisis.gwas;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;

public class CountProxies {
	public static final String[][] HEADER_EXPECTATIONS = {{"Region"}, {"MarkerName", "SNP", "RSID"}, {"Chr"}, {"Position", "BP"}, {"Rsq"}, {"Pval", "P-value"}, {"Replication"}};
	public static final double INDEX_THRESHOLD = 0.00000005;
	public static final double[] RSQ_THRESHOLDS = {0.80, 0.50, 0.30};
	
	public static void parse(String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp;
		Hashtable<String, Vector<String>> hash;
		Vector<String> v;
		int index, best;
		int[] indices, order;
		String[] keys;
		String region;
		String[] markerNames;
		double[] rsqs, pvalues, replications;
		
		hash = new Hashtable<String, Vector<String>>();
		line = Files.getHeaderOfFile(filename, "[\\s]+", new Logger());
		indices = ext.indexFactors(HEADER_EXPECTATIONS, line, false, true, true, true);
		
		try {
			reader = new BufferedReader(new FileReader(filename));
			reader.readLine();
			while (reader.ready()) {
				temp = reader.readLine();
				line = temp.trim().split("[\\s]+");
				region = line[indices[0]];
				if (!region.equals("NA")) {
					HashVec.addToHashVec(hash, region, temp, false);
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
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename)+"_proxies.xls"));
			writer.print("Region\tIndexName\tIndexPval\tIndexSignificant\tIndexAvailableOnReplicationArray\tIndexRsq\tIndexReplicationPval\tIndexNominallySigReplication\tIndexSigReplication");
			for (int i = 0; i < RSQ_THRESHOLDS.length; i++) {
				writer.print("\tRsq>"+RSQ_THRESHOLDS[i]+"_Available\tRsq>"+RSQ_THRESHOLDS[i]+"_IsIndex\tRsq>"+RSQ_THRESHOLDS[i]+"_Name\tRsq>"+RSQ_THRESHOLDS[i]+"_Pval\tRsq>"+RSQ_THRESHOLDS[i]+"_Rsq\tRsq>"+RSQ_THRESHOLDS[i]+"_ReplicationPval\tRsq>"+RSQ_THRESHOLDS[i]+"NominallySigReplication\tRsq>"+RSQ_THRESHOLDS[i]+"SigReplication");
			}
			writer.println();
			keys = HashVec.getKeys(hash, true, true);
			for (int keyIndex = 0; keyIndex < keys.length; keyIndex++) {
				v = hash.get(keys[keyIndex]);
//			for (int regi = 1; regi <= Integer.parseInt(keys[keys.length-1]); regi++) {
//				v = hash.get(regi+"");
				markerNames = new String[v.size()];
				rsqs = new double[v.size()];
				pvalues = new double[v.size()];
				replications = new double[v.size()];
				
				for (int i = 0; i < v.size(); i++) {
					line = v.elementAt(i).trim().split("[\\s]+");
					markerNames[i] = line[indices[1]];
					rsqs[i] = line[indices[4]].equals("NA")?-1:Double.parseDouble(line[indices[4]]);
					pvalues[i] = line[indices[5]].equals("NA")?999:Double.parseDouble(line[indices[5]]);
					replications[i] = line[indices[6]].equals("NA")?999:Double.parseDouble(line[indices[6]]);
				}
				
				order = Sort.quicksort(pvalues);
//				writer.print(regi);
				writer.print(keys[keyIndex]);
				
				index = order[0];
				if (pvalues[index] <= INDEX_THRESHOLD) {
					writer.print("\t"+markerNames[index]+"\t"+pvalues[index]+"\t1\t"+(rsqs[index]>RSQ_THRESHOLDS[RSQ_THRESHOLDS.length-1]?1:0)+"\t"+rsqs[index]+"\t"+replications[index]+"\t"+(replications[index]<0.05?1:0)+"\t"+(replications[index]<0.000373?1:0));
				} else {
					writer.print("\t.\t.\t0\t0\t.\t.\t0\t0");
				}
				for (int i = 0; i < RSQ_THRESHOLDS.length; i++) {
					best = -1;
					for (int j = 0; j < order.length; j++) {
						if (rsqs[order[j]] > RSQ_THRESHOLDS[i] && best == -1) {
							best = order[j];
						}
					}
//					writer.print("\tRsq>"+RSQ_THRESHOLDS[i]+"_Available\tRsq>"+RSQ_THRESHOLDS[i]+"_Name\tRsq>"+RSQ_THRESHOLDS[i]+"_Pval\tRsq>"+RSQ_THRESHOLDS[i]+"_Rsq\tRsq>"+RSQ_THRESHOLDS[i]+"_Replication");
					if (best == -1) {
						writer.print("\t0\t0\t.\t.\t.\t.\t0\t0");
					} else {
						writer.print("\t1\t"+(best==order[0]?1:0)+"\t"+markerNames[best]+"\t"+pvalues[best]+"\t"+rsqs[best]+"\t"+replications[best]+"\t"+(replications[best]<0.05?1:0)+"\t"+(replications[best]<0.000373?1:0));
					}
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.rootOf(filename)+"_proxies.xln");
			e.printStackTrace();
		}
		
		
		
		
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "input.txt";

		String usage = "\n" + "gwas.SelectProxies requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
//			parse(filename);
			parse("input_FE.txt");
			parse("input_RE2.txt");
			parse("input_Either.txt");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
