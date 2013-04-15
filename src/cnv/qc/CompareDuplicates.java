package cnv.qc;

import java.io.*;

import stats.Correlation;

import cnv.filesys.*;
import common.*;

public class CompareDuplicates {

	public static String doIt(Project proj, String[] pair) {
		Sample[] fsamps;
		byte[] chrs;
		float[][] xs, ys;
//		float[][] lrrs, bafs;
		int max, count, sameGeno;
		String summary;
		double[][] correl_icin;
		boolean use;
		byte[][] genos;
		byte geno;
		
		System.out.println("Comparing "+ext.listWithCommas(pair, true));

		summary = "";
		try {
//			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+Array.toStr(pair, "_")+".xln"));
//			writer.println("X_1\tX2\tY1\tY2");
			
			
			fsamps = new Sample[pair.length];
			xs = new float[pair.length][];
			ys = new float[pair.length][];
			genos = new byte[pair.length][];
			for (int i = 0; i < pair.length; i++) {
				fsamps[i] = proj.getFullSampleFromRandomAccessFile(pair[i]);
//				xs[i]= fsamps[i].getX_Raws();
//				ys[i] = fsamps[i].getY_Raws();
				xs[i]= fsamps[i].getXs();
				ys[i] = fsamps[i].getYs();
//				genos[i] = fsamps[i].getAB_Genotypes();
				genos[i] = fsamps[i].getForwardGenotypes();
			}

			chrs = proj.getMarkerSet().getChrs();
//			max = Array.indexOfByte(chrs, (byte)23);
//			max = 1000000;
			max = chrs.length;

			geno = 0;
			count = 0;
			sameGeno = 0;
			correl_icin = new double[pair.length][max];
			for (int i = 0; i < max; i++) {
				use = true;
				for (int j = 0; j < pair.length; j++) {
//					writer.print((i==0?"":"\t")+xs[j][i]);
					correl_icin[j][count] = xs[j][i];
//					if ((xs[j][i]+"").equalsIgnoreCase("NaN")) {
//						use = false;
//					}
					if (genos[j][i] == 0) {
						use = false;
					}
					if (j==0) {
						geno = genos[j][i];
					} else if (genos[j][i] != 0 && geno != 0 && genos[j][i] == geno) {
//						sameGeno++;
						sameGeno += 2-Math.abs(genos[j][i] - geno);
					}
				}
				if (use) {
					count++;
				}
//				writer.println();
			}
			sameGeno /= 2;
			
			for (int i = 0; i < pair.length; i++) {
				correl_icin[i] = Array.subArray(correl_icin[i], 0, count);
			}
			summary += ((double)sameGeno/((double)pair.length-1)/(double)count)+"\t"+Array.toStr(Correlation.Pearson(correl_icin));
			
//			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return summary;
	}
	
	public static void allPairs(Project proj, String pairFile) {
		PrintWriter writer;
		String[][] pairs;
		
		pairs = HashVec.loadFileToStringMatrix(proj.getProjectDir()+pairFile, false, new int[] {0,1}, proj.getJarStatus());
		
		try {
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+"DuplicateQC.xln"));
			for (int i = 0; i < pairs.length; i++) {
				writer.println(doIt(proj, pairs[i]));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + proj.getProjectDir()+"DuplicateQC.xln");
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = Project.DEFAULT_PROJECT;
		String pairs = "duplicatePairs.txt";

		String usage = "\n"+
		"cnv.qc.CompareDuplicates requires 0-1 arguments\n"+
		"   (1) project file (i.e. proj="+filename+" (default))\n"+
		"   (2) file with list of sample pairs (i.e. pairs="+pairs+" (default; to be found in project directory))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("pairs=")) {
				pairs = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			filename = "D:/home/npankrat/projects/TsaiPilot.properties";
//			filename = "D:/home/npankrat/projects/SDRG.properties";
			allPairs(new Project(filename, false), pairs);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
