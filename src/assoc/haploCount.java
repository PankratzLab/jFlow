package assoc;

import java.io.*;
import java.util.*;

import stats.LogisticRegression;
import stats.ProbDist;

import common.*;

public class haploCount {
	public static boolean USING_DIRS = true;

	public haploCount(String prefix) throws IOException {
		new haploCount(prefix, null, null);
	}

	public haploCount(String prefix, String deleterious, String protective) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line = {"0", "-1"};
		String temp, global = "missed global statistic";
		Hashtable<String,String> hash = new Hashtable<String,String>();
		int count, offset, numHaps;
		double[][] counts, output;
		double[][][] genotypes;
		boolean done;
		Vector<String> haps = new Vector<String>();
		String[][] scores;
		String dir, filename, bins, bintest, phenofile;
		int[] model = {-1, -1};
		int[][] logModel;

		if (USING_DIRS) {
			dir = "";
			line = new File(".").list();
			for (int i = 0; i<line.length; i++) {
				if (line[i].startsWith(prefix)&&new File(line[i]).isDirectory()&&new File(line[i]+"/"+prefix+"_haplos.out").exists()) {
					if (dir.equals("")) {
						dir = line[i]+"/";
					} else {
						System.err.println("Error - there are 2+ directories starting with the prefix "+prefix+"; using "+dir);
					}
				}
			}

		} else {
			dir = "";
		}
		filename = dir+prefix+"_haplos.out";
		bins = dir+prefix+"_em.out";
		bintest = dir+prefix+"_bintest.out";
		phenofile = dir+prefix+"-haplostats.dat";

		if (!new File(bins).exists()) {
			System.err.println("Error - could not find "+bins+" in current directory");
			System.exit(2);
		}
		reader = new BufferedReader(new FileReader(bins));
		for (int i = 0; i<5; i++) {
			// for (int i = 0; i<6; i++) {
			reader.readLine();
		}
		count = 1;
		done = false;
		while (!done) {
			line = reader.readLine().split("[\\s]+");
			offset = (line.length>1&&line[0].equals(""))?1:0;
			if (line[0+offset].equals(count+"")) {
				temp = line[1+offset];
				for (int i = 2+offset; i<line.length-1; i++) {
					temp += "-"+line[i];
				}
				haps.add(temp);
				if (temp.equals(deleterious)) {
					model[0] = count-1;
				}
				if (temp.equals(protective)) {
					model[1] = count-1;
				}
				count++;
			} else {
				done = true;
			}
		}
		reader.close();
		numHaps = haps.size();
		if (deleterious!=null&&protective!=null) {
			if (model[0]==-1) {
				System.err.println("Deleterious haplotype not found or not specified for "+prefix+"; skipping summary for deleterious haplotype");
				deleterious = null;
			}
			if (model[1]==-1) {
				System.err.println("Protective haplotype not found or not specified for "+prefix+"; skipping summary for protective haplotype");
				protective = null;
			}
			if (model[0]==model[1]&&model[0]>=0) {
				System.err.println("Error - deleterious haplotype specified for "+prefix+" is the same as the protective haplotype; summarizing only one of them");
				protective = null;
			}
		}

		if (!new File(bintest).exists()) {
			System.err.println("Error - could not find "+bintest+" in current directory");
			System.exit(2);
		}
		reader = new BufferedReader(new FileReader(bintest));
		for (int i = 0; i<13; i++) {
			// for (int i = 0; i<15; i++) {
			if (i==4) {
				global = reader.readLine();
			} else {
				reader.readLine();
			}
		}
		scores = new String[numHaps][2];
		for (int i = 0; i<numHaps; i++) {
			scores[i][0] = ".";
			scores[i][1] = "-";
		}
		done = false;
		while (!done) {
			line = reader.readLine().split("[\\s]+");
			offset = (line.length>1&&line[0].equals(""))?1:0;
			if (line.length<2) {
				done = true;
			} else {
				temp = line[1+offset];
				for (int i = 2+offset; i<line.length-3; i++) {
					temp += "-"+line[i];
				}
				if (haps.contains(temp)) {
					scores[haps.indexOf(temp)][0] = line[line.length-2];
					scores[haps.indexOf(temp)][1] = line[line.length-1];
				} else {
					System.err.println("Error- "+temp+" was not found in "+bins);
				}
			}
		}
		reader.close();

		if (!new File(phenofile).exists()) {
			System.err.println("Error - could not find "+phenofile+" in current directory");
			System.exit(2);
		}
		reader = new BufferedReader(new FileReader(phenofile));
		count = 0;
		while (reader.ready()) {
			line = reader.readLine().split("[\\s]+");
			hash.put(line[0], line[1]);
		}
		reader.close();

		if (!new File(filename).exists()) {
			System.err.println("Error - could not find "+filename+" in current directory");
			System.exit(2);
		}
		reader = new BufferedReader(new FileReader(filename));

		count = 1;
		counts = new double[2][numHaps+1];
		genotypes = new double[2][numHaps+1][numHaps+1];
		done = false;
		while (!done) {
			do {
				line = reader.readLine().split("[\\s]+");
				offset = (line.length>1&&line[0].equals(""))?1:0;
				if (line.length>1&&line[1].equals("Number")) {
					done = true;
				}
			} while (!done&&!line[0+offset].equals(count+""));
			if (!done) {
				counts[Integer.valueOf(hash.get(line[1+offset])).intValue()][Integer.valueOf(line[2+offset]).intValue()-1] += Double.valueOf(line[4+offset]).doubleValue();
				counts[Integer.valueOf(hash.get(line[1+offset])).intValue()][Integer.valueOf(line[3+offset]).intValue()-1] += Double.valueOf(line[4+offset]).doubleValue();
				genotypes[Integer.valueOf(hash.get(line[1+offset])).intValue()][Integer.valueOf(line[2+offset]).intValue()-1][Integer.valueOf(line[3+offset]).intValue()-1] += Double.valueOf(line[4+offset]).doubleValue();
				count++;
			}
		}
		reader.close();

		writer = new PrintWriter(new FileWriter(dir+prefix+"_summary.out"));
		for (int i = 0; i<numHaps; i++) {
			counts[0][numHaps] += counts[0][i];
			counts[1][numHaps] += counts[1][i];
			for (int j = i; j<numHaps; j++) {
				genotypes[0][numHaps][numHaps] += genotypes[0][i][j];
				genotypes[1][numHaps][numHaps] += genotypes[1][i][j];
			}
		}
		writer.println("Haplostats summary");
		writer.println("Haplotype\tCases\tfreq\tControls\tfreq\tHapScore\tp-value");
		for (int i = 0; i<numHaps+1; i++) {
			writer.println((i<numHaps?"'"+haps.elementAt(i):"Sum total")+"\t"+ext.formDeci(counts[1][i], 2, true)+"\t"+ext.formDeci(counts[1][i]/counts[1][numHaps], 3, true)+"\t"+ext.formDeci(counts[0][i], 2, true)+"\t"+ext.formDeci(counts[0][i]/counts[0][numHaps], 3, true)+"\t"+(i==numHaps||scores[i][0].equals(".")?"":ext.formDeci(Double.valueOf(scores[i][0]).doubleValue(), 3, true)+"")+"\t"+(i==numHaps||scores[i][1].equals("-")?"":ext.formDeci(Double.valueOf(scores[i][1]).doubleValue(), 5, true)+""));
		}
		writer.println();
		writer.println(global);

		for (int mode = 0; mode<2; mode++) {
			if (model[mode]>=0) {
				temp = (mode==0?"del":"pro");
				counts = new double[2][3+1];
				writer.println();
				writer.println();
				writer.println("Summary for "+(mode==0?"deleterious":"protective")+" haplotype");
				writer.println("\t"+temp+"/"+temp+"\t"+temp+"/*\t*/*\t\tModel\tOR\tsig");
				for (int i = 0; i<numHaps; i++) {
					for (int j = i; j<numHaps; j++) {
						count = 1+(i==model[mode]?0:1)+(j==model[mode]?0:1);

						counts[0][count] += genotypes[0][i][j];
						counts[1][count] += genotypes[1][i][j];
						counts[0][0] += genotypes[0][i][j];
						counts[1][0] += genotypes[1][i][j];
					}
				}
				for (int i = 1; i>=-1; i--) {
					writer.print(i==1?"cases":(i==0?"controls":"diff"));
					for (int j = 1; j<=3; j++) {
						writer.print("\t"+ext.formDeci(i==-1?counts[1][j]/counts[1][0]-counts[0][j]/counts[0][0]:counts[i][j]/counts[i][0], 3, true));
					}
					logModel = new int[][] { {(i==-1?2:1)}, {(i==0?0:1)}, {0}};
					output = fakeLogistic(counts, logModel, 10);
					writer.println("\t\t"+(i==1?"dominant":(i==0?"recessive":"additive"))+"\t"+ext.formDeci(output[0][0], 3, true)+"\t"+ext.formDeci(output[1][0], 3, true));
				}
			}
		}

		if (model[0]>=0&&model[1]>=0) {
			counts = new double[2][6+1];
			writer.println();
			writer.println();
			writer.println("Summary for deleterious and protective haplotypes");
			writer.println("\tdel/del\tdel/*\tdel/pro\t*/*\tpro/*\tpro/pro\t\tModel\tdel-OR\tdel-sig\tpro-OR\tpro-sig");
			for (int i = 0; i<numHaps; i++) {
				for (int j = i; j<numHaps; j++) {
					if (i==model[0]&&j==model[0]) {
						count = 1;
					} else if (i==model[1]&&j==model[1]) {
						count = 6;
					} else if ((i==model[0]||j==model[0])&&(i!=model[1]&&j!=model[1])) {
						count = 2;
					} else if ((i==model[0]||j==model[0])&&(i==model[1]||j==model[1])) {
						count = 3;
					} else if ((i!=model[0]&&j!=model[0])&&(i!=model[1]&&j!=model[1])) {
						count = 4;
					} else if ((i==model[1]||j==model[1])&&(i!=model[0]&&j!=model[0])) {
						count = 5;
					} else {
						System.err.println("Unforseen is this");
					}

					counts[0][count] += genotypes[0][i][j];
					counts[1][count] += genotypes[1][i][j];
					counts[0][0] += genotypes[0][i][j];
					counts[1][0] += genotypes[1][i][j];
				}
			}
			for (int i = 1; i>=-1; i--) {
				writer.print(i==1?"cases":(i==0?"controls":"diff"));
				for (int j = 1; j<=6; j++) {
					writer.print("\t"+ext.formDeci(i==-1?counts[1][j]/counts[1][0]-counts[0][j]/counts[0][0]:counts[i][j]/counts[i][0], 3, true));
				}
				logModel = new int[][] { {(i==-1?2:1), 0}, {(i==0?0:1), 0}, {(i==0?0:1), (i==0?0:1)}, {0, 0}, {0, (i==0?0:1)}, {0, (i==-1?2:1)}};
				output = fakeLogistic(counts, logModel, 10);
				writer.println("\t\t"+(i==1?"dominant":(i==0?"recessive":"additive"))+"\t"+ext.formDeci(output[0][0], 3, true)+"\t"+ext.formDeci(output[1][0], 3, true)+"\t"+ext.formDeci(output[0][1], 3, true)+"\t"+ext.formDeci(output[1][1], 3, true));
			}
		}

		writer.println();
		writer.println();
		writer.println("Full distribution of haplotypes");
		writer.println("Haplotype1\tHaplotype2\tCases\tfreq\tControls\tfreq");
		for (int i = 0; i<numHaps; i++) {
			for (int j = i; j<numHaps; j++) {
				writer.println("'"+haps.elementAt(i)+"\t'"+haps.elementAt(j)+"\t"+ext.formDeci(genotypes[1][i][j], 2, true)+"\t"+ext.formDeci(genotypes[1][i][j]/genotypes[1][numHaps][numHaps], 3, true)+"\t"+ext.formDeci(genotypes[0][i][j], 2, true)+"\t"+ext.formDeci(genotypes[0][i][j]/genotypes[0][numHaps][numHaps], 3, true));
			}
		}
		writer.println("Sum total\t\t"+ext.formDeci(genotypes[1][numHaps][numHaps], 2, true)+"\t1.0\t"+ext.formDeci(genotypes[0][numHaps][numHaps], 2, true)+"\t1.0");
		writer.close();

	}

	public static double[][] fakeLogistic(double[][] counts, int[][] translations, int multiplier) throws IOException {
		LogisticRegression lr;
		Vector<String> deps = new Vector<String>();
		Vector<int[]> indeps = new Vector<int[]>();
		double[] walds;
		double[][] output, ORs;

		for (int i = 0; i<counts.length; i++) {
			for (int j = 1; j<counts[i].length; j++) {
				for (int k = 0; k<(int)Math.round(counts[i][j]*multiplier); k++) {
					deps.add(i+"");
					indeps.add(translations[j-1]);
				}
			}
		}
		lr = new LogisticRegression(deps, indeps);
		walds = lr.getWalds();
		ORs = lr.getOddsRatios();
		output = new double[2][walds.length];
		for (int i = 0; i<ORs.length; i++) {
			output[0][i] = ORs[i][0];
		}
		for (int i = 0; i<walds.length; i++) {
			output[1][i] = ProbDist.ChiDist(walds[i]/multiplier, 1);
		}

		return output;
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String prefix = "4by4", phenofile = "aft.dat";

		String usage = "\n"+"park.haploCount requires 2 arguments:\n"+"   (1) the prefix for the output from haplostats (i.e. prefix="+prefix+" (default))\n"+"       (for the individual haplotypes, the em and the binary test output)\n"+"   (2) the pheno file used in haplostats (i.e. pheno="+phenofile+" (default))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("prefix=")) {
				prefix = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("pheno=")) {
				phenofile = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		if (args.length==0) {
			System.err.println("Using defaults (prefix="+prefix+"), corresponding to "+prefix+"_haplos.out, "+prefix+"_em.out and "+prefix+"_bintest.out)");
		}

		try {
			// new haploCount("all", phenofile);
			// new haploCount("f3", phenofile);
			// new haploCount("l3", phenofile);
			// new haploCount("f2", phenofile);
			// new haploCount("m2", phenofile);
			// new haploCount("l2", phenofile);

			new haploCount("markers391,395", "2-1", "1-2");
			new haploCount("markers391,393,395", "2-1-1", "1-2-2");
			new haploCount("markers609,611", "2-2", "1-1");
			new haploCount("markers646,649,651", "2-1-2", "1-1-1");
			new haploCount("markers646,651", "2-2", "1-1");

		} catch (Exception e) {
			e.printStackTrace();
		}

	}
}
