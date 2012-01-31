package link;

import java.io.*;
import java.util.*;

import common.*;

public class alleleDistribution {
	public alleleDistribution(String filename, boolean sum) throws IOException {
		new alleleDistribution(filename, sum, 1, 23);
	}

	public alleleDistribution(String filename, boolean sum, int start, int stop) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer, suggest, oddevens, shiftoes;
		Hashtable<String,String> hash, notthese;
		String[] line;
		String temp, chrome, trav;
		Vector<String> species;
		double[][][] counts;
		String[][] names;
		int numMarkers, plate;
		double[] total, maxis, diffs;
		double score, a, b, avg;
		int numIncr, offset, off_mod;
		int[] keys;
		double[][] oe_sums;

		notthese = new Hashtable<String,String>();
		if (new File("notthese.dat").exists()) {
			reader = new BufferedReader(new FileReader("notthese.dat"));
			while (reader.ready()) {
				notthese.put(reader.readLine(), "");
			}
			reader.close();
		}

		try {
			reader = new BufferedReader(new FileReader(filename));
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error - Cannot find "+filename+" in current directory");
			System.exit(1);
		}

		hash = new Hashtable<String,String>();
		species = new Vector<String>();
		while (reader.ready()) {
			line = reader.readLine().split("[\\s]+");
			if (!species.contains(line[2])) {
				species.add(line[2]);
			}
			if (!hash.containsKey(line[0]+"\t"+line[1])) {
				hash.put(line[0]+"\t"+line[1], line[2]);
			} else if (!(hash.get(line[0]+"\t"+line[1])).equals(line[2])) {
				// System.err.println(line[0]+"\t"+line[1]+" is on both
				// "+hash.get(line[0]+"\t"+line[1])+" and "+line[2]);
			}
		}
		reader.close();

		shiftoes = new PrintWriter(new FileWriter("shiftAllelesForOddEvens.txt"));
		oddevens = new PrintWriter(new FileWriter("odd-evens.xls"));
		oddevens.print("Marker");
		for (int i = 0; i<species.size(); i++) {
			oddevens.print("\t%off "+species.elementAt(i));
		}
		oddevens.println("\t%off ALL\tmax %offr");

		suggest = new PrintWriter(new FileWriter("suggestions.out"));
		for (int i = start; i<=stop; i++) {
			chrome = (i<10)?"0"+i:""+i;
			try {
				reader = new BufferedReader(new FileReader("map"+chrome+".dat"));
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error - Cannot find "+"map"+chrome+".dat"+" in current directory");
				System.exit(1);
			}
			line = reader.readLine().split("[\\s]+");
			offset = (line[0].equals("")?1:0);
			numMarkers = Integer.valueOf(line[0+offset]).intValue()-1;
			for (int j = 0; j<6; j++) {
				reader.readLine();
			}
			if (i==23) {
				reader.readLine();
			}
			names = new String[numMarkers][];
			counts = new double[numMarkers][][];
			for (int j = 0; j<numMarkers; j++) {
				line = reader.readLine().split("[\\s]+");
				offset = (line[0].equals("")?1:0);
				names[j] = new String[Integer.valueOf(line[1+offset]).intValue()+1];
				names[j][0] = line[3+offset];
				for (int k = 1; k<names[j].length; k++) {
					if (line.length>5) {
						names[j][k] = line[k+5+offset];
					} else {
						names[j][k] = k+"";
					}
				}
				counts[j] = new double[names[j].length][species.size()+1];
				line = reader.readLine().split("[\\s]+");
				offset = (line[0].equals("")?1:0);
				for (int k = 1; k<names[j].length; k++) {
					if (sum) {
						counts[j][k][0] = 0;
					} else {
						counts[j][k][0] = Double.valueOf(line[k-1+offset]).doubleValue();
					}
				}
			}
			reader.close();

			if (!new File("re_chrom"+chrome+".pre").exists()) {
				System.err.println("Error - could not find "+"re_chrom"+chrome+".pre"+" in current directory");
				System.exit(2);
			}
			reader = new BufferedReader(new FileReader("re_chrom"+chrome+".pre"));
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				plate = species.indexOf(hash.get(line[0]+"\t"+line[1]))+1;
				if (plate==0&&Integer.valueOf(line[1]).intValue()<100) {
					System.err.println("Error - don't know which plate "+line[0]+"\t"+line[1]+" is on");
				}
				for (int j = 0; j<numMarkers; j++) {
					try {
						counts[j][Integer.valueOf(line[j*2+6]).intValue()][plate]++;
						counts[j][Integer.valueOf(line[j*2+6+1]).intValue()][plate]++;
						if (sum) {
							counts[j][Integer.valueOf(line[j*2+6]).intValue()][0]++;
							counts[j][Integer.valueOf(line[j*2+6+1]).intValue()][0]++;
						}
					} catch (Exception e) {
						System.err.println("Error - marker allele outside expected range for marker "+(j+1)+" of individual "+line[0]+"-"+line[1]);
						System.err.println("        alleles "+(Integer.valueOf(line[j*2+6]).intValue()+1)+"/"+(Integer.valueOf(line[j*2+6+1]).intValue()+1)+" outside maximum of "+counts[j].length);
						System.exit(4);
					}
				}
			}
			reader.close();

			writer = new PrintWriter(new FileWriter("alleles"+chrome+".xls"));
			for (int j = 0; j<numMarkers; j++) {
				total = new double[species.size()+1];
				writer.println(names[j][0]);
				oddevens.print(names[j][0]);
				numIncr = 0;
				for (int k = 1; k<names[j].length; k++) {
					writer.print("\t"+names[j][k]);
					if (k!=1&&Integer.valueOf(names[j][k]).intValue()-Integer.valueOf(names[j][k-1]).intValue()==1) {
						numIncr++;
					}
					for (int l = 0; l<species.size()+1; l++) {
						total[l] += counts[j][k][l];
					}
				}
				writer.println();
				diffs = new double[species.size()];
				maxis = new double[species.size()*(names[j].length-1)];
				oe_sums = new double[2][species.size()+2];
				for (int l = 0; l<species.size()+1; l++) {
					writer.print((l==0?"source":species.elementAt(l-1)));
					for (int k = 1; k<names[j].length; k++) {
						counts[j][k][l] = (counts[j][k][l]==0?0:counts[j][k][l]/(l==0&&!sum?1.0:total[l]));
						writer.print("\t"+ext.formDeci(counts[j][k][l], 6, true));
						if (l!=species.size()+1) {
							oe_sums[Integer.valueOf(names[j][k]).intValue()%2][l] += counts[j][k][l];
							oe_sums[Integer.valueOf(names[j][k]).intValue()%2][species.size()+1] += counts[j][k][l];
						}
					}
					writer.println();
				}
				score = 0;
				off_mod = (oe_sums[1][species.size()+1]>oe_sums[0][species.size()+1]?0:1);
				trav = names[j][0];
				temp = "";
				for (int k = 1; k<names[j].length; k++) {
					trav += "\t"+names[j][k];
					if (Integer.valueOf(names[j][k]).intValue()%2!=off_mod) {
						temp += "\t"+names[j][k];
					} else {
						temp += "\t0";
					}
				}

				// if (names[j][0].equals("D2S396")) {
				// boolean crap = true;
				// }
				for (int l = 1; l<species.size()+2; l++) {
					oddevens.print("\t"+oe_sums[off_mod][l]);
					if (oe_sums[off_mod][l]>score&&l<=species.size()) {
						score = oe_sums[off_mod][l];
					}
					if (oe_sums[off_mod][l]>0&&l<=species.size()&&!notthese.containsKey(names[j][0])) {
						shiftoes.println(trav);
						shiftoes.println(species.elementAt(l-1)+temp);
					}
				}
				oddevens.println("\t"+score);

				score = 0;
				for (int l = 1; l<=species.size(); l++) {
					for (int k = 1; k<names[j].length; k++) {
						avg = 0;
						for (int m = 1; m<=species.size(); m++) {
							a = (counts[j][k][l]>=counts[j][k][m]?counts[j][k][l]:counts[j][k][m]);
							b = (counts[j][k][l]>=counts[j][k][m]?counts[j][k][m]:counts[j][k][l]);

							diffs[l-1] += (a-b)*(a==0||b==0?1:a/b);
							score += (a-b)*(a==0||b==0?1:a/b);

							if (m!=l) {
								avg += counts[j][k][m];
							}
						}
						maxis[(k-1)*species.size()+(l-1)] = Math.abs(counts[j][k][l]-avg/(species.size()-1));
					}
				}

				keys = Sort.quicksort(diffs);
				suggest.print(names[j][0]+"\t"+ext.formDeci(score/(names[j].length-1), 2)+"\t"+(names[j].length-1)+"\t"+numIncr+"\t"+ext.formDeci(diffs[keys[diffs.length-1]], 3));
				keys = Sort.quicksort(maxis);
				suggest.println("\t"+ext.formDeci(maxis[keys[maxis.length-1]], 3)+"\t"+ext.formDeci(maxis[keys[maxis.length-2]], 3)+"\t"+ext.formDeci(maxis[keys[maxis.length-3]], 3));
				writer.println();
			}
			writer.close();

		}
		suggest.close();
		oddevens.close();
		shiftoes.close();

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String filename = "plateList.dat";
		boolean sum = true;

		String usage = "\n"+"park.alleleDistribution requires 1-2 arguments:\n"+"   (1) a list of DNAs and their plates (2 cols)\n"+"       (i.e. file="+filename+" (default))\n"+"   (2) calculate total (i.e. sum=true (default))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("sum=")) {
				sum = args[i].split("=")[1].equals("true");
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		if (args.length==0) {
			System.err.println("Warning: using defaults (file="+filename+")");
		}
		try {
			new alleleDistribution(filename, sum);
			// new alleleDistribution(filename, sum, 5, 5);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
