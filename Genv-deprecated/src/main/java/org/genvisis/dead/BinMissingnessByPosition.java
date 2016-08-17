package org.genvisis.dead;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;

public class BinMissingnessByPosition {
	public static final String MAP = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\documentation\\HumanCNV370v1 Marker Info files\\HumanCNV370_MarkerList\\gwas.map";
	public static final double TARGET_BIN_SIZE = 5000000.0;
	public static final int NUM_CHRS = 25;
	public static final double[] THRESHOLDS = {1E-2, 1E-3, 1E-4, 1E-5, 1E-10, 1E-15, 1E-20};

	public static void prepFile(String filename) {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;
		Hashtable<String,String> hash = new Hashtable<String,String>();

		try {
			reader = new BufferedReader(new FileReader(MAP));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				hash.put(line[1], line[3]);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+MAP+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+MAP+"\"");
			System.exit(2);
		}
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(filename+".out"));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				writer.println(Array.toStr(line)+"\t"+hash.get(line[1]));

			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}
	}

	public static void binByPosition(String filename) {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;
		int[] maxes = new int[NUM_CHRS];
		int chr, pos, index;
		int[] numBins = new int[NUM_CHRS];
		double[] binSizes = new double[NUM_CHRS];
		int[][][] counts = new int[NUM_CHRS][][];
		double pval, prop;

		System.out.println(ext.getTime());
		try {
			reader = new BufferedReader(new FileReader(filename));
			reader.readLine();
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				chr = Integer.parseInt(line[0]);
				pos = Integer.parseInt(line[5]);
				if (maxes[chr-1]<pos) {
					maxes[chr-1] = pos;
				}
			}
			for (int i = 0; i<maxes.length; i++) {
				numBins[i] = (int)Math.ceil((maxes[i])/TARGET_BIN_SIZE);
				binSizes[i] = (maxes[i]+1)/Math.ceil((maxes[i])/TARGET_BIN_SIZE);
				counts[i] = new int[numBins[i]][THRESHOLDS.length+1];
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}

		try {
			reader = new BufferedReader(new FileReader(filename));
			reader.readLine();
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				chr = Integer.parseInt(line[0]);
				pos = Integer.parseInt(line[5]);
				index = (int)Math.floor(pos/binSizes[chr-1]);
				if (index>=numBins[chr-1]) {
					System.err.println("Error - marker "+line[1]+" was set to bin "+index+" (max: "+numBins[chr-1]+" for chr "+chr+")");
				}
				pval = Double.parseDouble(line[4]);
				for (int i = 0; i<THRESHOLDS.length; i++) {
					if (pval<THRESHOLDS[i]) {
						counts[chr-1][index][i]++;
					}
				}
				counts[chr-1][index][THRESHOLDS.length]++;
			}
			reader.close();

			writer = new PrintWriter(new FileWriter(filename+"_binned.xln"));
			writer.println("Chr\tbin start pos\tbin stop pos\t"+Array.toStr(THRESHOLDS));
			for (int i = 0; i<counts.length; i++) {
				for (int j = 0; j<counts[i].length; j++) {
					writer.print((i+1)+"\t"+ext.formDeci((j+0)*binSizes[i], 1)+"\t"+ext.formDeci((j+1)*binSizes[i], 1));
					for (int k = 0; k<THRESHOLDS.length; k++) {
						prop = (double)counts[i][j][k]/(double)counts[i][j][THRESHOLDS.length];
						writer.print("\t"+(Double.isNaN(prop)?"-0.1":ext.formDeci(prop, 4)));
					}
					writer.println();
				}

			}
			writer.close();
			System.out.println(ext.getTime());
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		// String filename = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\blood v WGA\\blood_v_WGA.missing";
		// String filename = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\blood v WGA\\blood_v_WGA2.missing";
		String filename = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\blood v WGA\\blood_v_LCL.missing";
		// String filename = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\blood v WGA\\testmissi.missing";

		String usage = "\\n"+"park.gwa.BinMissingnessByPosition requires 0-1 arguments\n"+"   (1) filename (i.e. file="+filename+" (default))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (!new File(filename+".out").exists()) {
				System.out.println("Adding map positions");
				prepFile(filename);
			} else {
				System.out.println("Using existing converted file (map positions already added)");
			}
			binByPosition(filename+".out");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
