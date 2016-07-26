package org.genvisis.gaw17;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;
import org.genvisis.seq.WeightedSumStatistic;
import org.genvisis.stats.ProbDist;

public class RareVariants {
	public static final double[] MAFS = {0.01, 0.02, 0.05};
	
	public static void test(String root, String filename) {
		WeightedSumStatistic wss;
		String[][] subsets;
		long time;

		double mean, stdev, z;
		double[] perms;
		
		wss = new WeightedSumStatistic(root);
		subsets = HashVec.loadFileToStringMatrix(filename, false, null, false);

		time = new Date().getTime();
		for (int i = 0; i < subsets.length; i++) {
			System.out.println(subsets[i][0]);
			wss.setSubset(Array.subArray(subsets[i], 1));
			wss.setThresholdMAF(0.02);
			perms = wss.getPermutations();
	        mean = Array.mean(perms);
	        stdev = Array.stdev(perms);
	        z = (wss.getStat() - mean)/stdev;
	        
	        System.out.println(wss.getMarkerNames().length+"\t"+ext.formDeci(z, 3)+"\t"+ext.prettyP(ProbDist.NormDist(Math.abs(z)*-1))+"\t"+ext.prettyP(wss.getEmpiricalSig()));
		}
		System.out.println("Finished in " + ext.getTimeElapsed(time));
		
	}
	
	public static void run(String analysis_dir, String files, String pop, double maf) {
		PrintWriter writer;
		WeightedSumStatistic wss;
		String[][][] subsets;
		String[] markerSets;
		Hashtable<String,String> hash;
		String[][] ids;
		double[] trait;
		String trav;
		String filename;
		int[] order;
		byte[] controlDesignations;
		
		markerSets = files.split(";");

		subsets = new String[markerSets.length][][];
		for (int k = 0; k < markerSets.length; k++) {
			subsets[k] = HashVec.loadFileToStringMatrix(analysis_dir+markerSets[k], false, null, false);
		}
		System.out.println(ext.getTime()+"\t"+pop);
		wss = new WeightedSumStatistic(analysis_dir+pop);
		wss.setThresholdMAF(maf);
		for (int rep = 1; rep <= 1; rep++) {
//			for (int j = 0; j < 1; j++) {
			for (int j = 0; j < Traits.PHENOS.length; j++) {
				System.out.println(ext.getTime()+"\t"+Traits.PHENOS[j]);
//				filename = "phenos/unr_phen."+rep+".tab";
				filename = "phenos/pheno1_adj.xln";
				hash = HashVec.loadFileToHashString(analysis_dir+filename, new int[] {0,1},
						ext.indexFactors(new String[] {Traits.PHENOS[j]}, Files.getHeaderOfFile(analysis_dir+filename, "\t", new Logger()), false, true),
						false, "\t", true, false, false);
				ids = wss.getIDs();
				trait = new double[ids.length];
				for (int i = 0; i < ids.length; i++) {
					trav = hash.get(ids[i][0]+"\t"+ids[i][1]);
					if (trav == null) {
						System.err.println("Error - no pheno '"+Traits.PHENOS[j]+"' for indiviudal "+ids[i][0]+","+ids[i][1]);
						System.exit(1);
					}
					trait[i] = Double.parseDouble(trav);
				}
				wss.setTrait(trait);
				order = Sort.quicksort(trait);
				controlDesignations = new byte[ids.length];
				for (int i = 0; i < ids.length; i++) {
					controlDesignations[i] = (order[i] >= ids.length/3 && order[i] <= ids.length*2/3)?(byte)1:(byte)0;
				}
				wss.setAffectionStatus(controlDesignations);
				for (int k = 0; k < markerSets.length; k++) {
					try {
						writer = new PrintWriter(new FileWriter(analysis_dir+"mid_"+pop+"."+rep+"."+Traits.PHENOS[j]+".maf"+maf+"."+markerSets[k]));
						writer.println("Gene\tn_SNPs\tpval\tEMP1");
						for (int m = 0; m < subsets[k].length; m++) {
							if (m%30==0) {
								System.out.print(".");
								writer.flush();
							}
							wss.setSubset(Array.subArray(subsets[k][m], 1));
							writer.println(subsets[k][m][0]+"\t"+wss.getMarkerNames().length+"\t"+ext.prettyP(wss.getSig())+"\t"+ext.prettyP(wss.getEmpiricalSig()));
						}
						System.out.println();
						writer.close();
					} catch (Exception e) {
						System.err.println("Error writing to " + analysis_dir+pop+"."+rep+"."+Traits.PHENOS[j]+".maf"+maf+"."+markerSets[k]);
						e.printStackTrace();
					}
				}
			}
		}
		System.out.println(ext.getTime()+"\tDone!");
	}
	
	public static void mergeAll(String analysis_dir, String files) {
		PrintWriter writer;
		String[] markerSets;
		
		markerSets = files.split(";");

		try {
			writer = new PrintWriter(new FileWriter(analysis_dir+"mergeAll.crf"));
			writer.println("lookup");
			writer.println("answerSets.dat 0 out=allRep1Results.xln");
			for (int maf = 0; maf < MAFS.length; maf++) {
				for (int j = 0; j < Traits.PHENOS.length; j++) {
					for (int k = 0; k < markerSets.length; k++) {
						for (int i = 0; i < Traits.POPS.length; i++) {
							writer.println("mid_"+Traits.POPS[i]+"."+1+"."+Traits.PHENOS[j]+".maf"+MAFS[maf]+"."+markerSets[k]+" 0 3="+Traits.POPS[i]+"."+1+"."+Traits.PHENOS[j]+".maf"+MAFS[maf]+"."+markerSets[k]);
						}
					}
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + analysis_dir+"mergeAll.crf");
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
//		String analysis_dir = "D:\\GAW17\\analysis\\";
//		String filenames = "D:\\GAW17\\source\\answerSets.dat;D:\\GAW17\\source\\markerSets.dat;D:\\GAW17\\source\\nonynonMarkerSets.dat";
		String analysis_dir = "";
//		String filenames = "answerSets.dat";
		String filenames = "answerSetsNonsynon.dat";
		String pop = "Caucasians";
		double maf = 0.01;

		String usage = "\n" +
		"gaw17.RareVariants requires 0-1 arguments\n" +
		"   (1) directory (i.e. dir=" + analysis_dir + " (default))\n" + 
		"   (2) markerSet filenames (i.e. files=" + filenames + " (default))\n" + 
		"   (3) population (i.e. pop=" + pop + " (default))\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				analysis_dir = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("files=")) {
				filenames = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("pop=")) {
				pop = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("maf=")) {
				maf = ext.parseDoubleArg(args[i]);
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
//			test("plink", "sets.dat");
			for (int i = 0; i < Traits.POPS.length; i++) {
				run(analysis_dir, filenames, pop, maf);
				run(analysis_dir, filenames, Traits.POPS[i], maf);
			}
//			mergeAll("D:\\GAW17\\analysis\\", filenames);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
 	