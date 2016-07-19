// should create a chi square for binomials
package stats;

import java.io.*;
import java.util.*;

import common.Array;
import common.DoubleVector;
import common.HashVec;
import common.Sort;
import common.ext;

public class PermuteOnePer {
	public static final int SIGFIGS = 2;
	public static final int NUM_REPS = 1000;
	public static final int UPPER_LIMIT = 15;
	private int numTraits;
	private double[][] nQdata;
	private String[][] variableNames;
	private double[][] aggregateCounts;
	private double[][] means;
	private double[][] stdevs;
	private double[] stats;
	private double[] pVals;
	private double[][][] results;
	private int[] counts;
	private int[] offsets;
	private boolean binary;

	public PermuteOnePer(String[] famIDs, double[] trait, int[][] variables) {
		this(famIDs, trait, variables, new String[][] {new String[] {"Trait"}, Array.stringArraySequence(variables.length, "Var ")});
	}

	public PermuteOnePer(String[] famIDs, double[] trait, int[][] variables, String[][] varNames) {
		Hashtable<String,DoubleVector[][]> hash = new Hashtable<String,DoubleVector[][]>();
		Hashtable<String,double[][]> cHash;
		DoubleVector[][] sorted;
		Ttest tt;
		Vector<String> v;
		int count;
		int[] keys, order;
		double[][] contTable;

		numTraits = variables[0].length;
		variableNames = varNames;
		counts = new int[numTraits];
		offsets = new int[numTraits];
		for (int i = 0; i<=numTraits; i++) {
			count = 0;
			v = new Vector<String>();
			while (v.size()<=UPPER_LIMIT&&count<famIDs.length) {
				if (i==numTraits&&!Double.isNaN(trait[count])) {
					HashVec.addIfAbsent(trait[count]+"", v);
				} else if (i<numTraits&&variables[count][i]!=Integer.MIN_VALUE) {
					HashVec.addIfAbsent(variables[count][i]+"", v);
				}
				count++;
			}
			if (i==numTraits) {
				binary = v.size()==2;
			} else {
				if (v.size()>UPPER_LIMIT) {
					throw new RuntimeException("Kind of hard to interpret more than "+UPPER_LIMIT+" levels, don't you think? (for '"+variableNames[1][i]+"')");
				}
				order = Array.toIntArray(Array.toStringArray(v));
				keys = Sort.quicksort(order);
				counts[i] = order[keys[keys.length-1]]-order[keys[0]]+1;
				offsets[i] = order[keys[0]];

				if (counts[i]>UPPER_LIMIT) {
					throw new RuntimeException("The spread for '"+variableNames[1][i]+"' is a bit large, don't you think ("+offsets[i]+"-"+(offsets[i]+counts[i]-1)+")?");
				}
			}
		}

		for (int i = 0; i<famIDs.length; i++) {
			if (variableNames[1].length!=variables[i].length) {
				throw new RuntimeException("Error - mismatched number of columns for row starting '"+famIDs[i]+"'");
			}
			for (int trt = 0; trt<numTraits; trt++) {
				if (variables[i][trt]!=Integer.MIN_VALUE) {
					if (hash.containsKey(famIDs[i])) {
						sorted = hash.get(famIDs[i]);
					} else {
						sorted = new DoubleVector[numTraits][];
						for (int j = 0; j<numTraits; j++) {
							sorted[j] = new DoubleVector[counts[j]];
							for (int k = 0; k<counts[j]; k++) {
								sorted[j][k] = new DoubleVector();
							}
						}
						hash.put(famIDs[i], sorted);
					}
					sorted[trt][(int)variables[i][trt]-offsets[trt]].add(trait[i]);
				}
			}
		}

		results = new double[numTraits][][];
		means = new double[numTraits][];
		stdevs = new double[numTraits][];
		aggregateCounts = new double[numTraits][];
		nQdata = new double[numTraits][3];
		stats = new double[numTraits];
		pVals = new double[numTraits];
		for (int trt = 0; trt<numTraits; trt++) {
			cHash = convertHashtable(hash, trt);
			results[trt] = permuteMeans(cHash, counts[trt]);
			means[trt] = results[trt][0];
			aggregateCounts[trt] = results[trt][1];
			stdevs[trt] = permuteStdev(cHash, counts[trt]);
			if (binary) {
				contTable = new double[2][aggregateCounts[trt].length];
				for (int i = 0; i<aggregateCounts[trt].length; i++) {
					contTable[1][i] = aggregateCounts[trt][i]*means[trt][i];
					contTable[0][i] = aggregateCounts[trt][i]-contTable[1][i];
				}
				stats[trt] = ContingencyTable.ChiSquare(contTable, false, true);
				pVals[trt] = ContingencyTable.ChiSquareOptimizedSignificance(contTable, true, true);
				nQdata[trt][0] = -999;
				nQdata[trt][1] = aggregateCounts[trt][1];
				nQdata[trt][2] = aggregateCounts[trt][0];
			} else if (counts[trt]==2) {
				tt = new Ttest(means[trt][1], stdevs[trt][1], aggregateCounts[trt][1], means[trt][0], stdevs[trt][0], aggregateCounts[trt][0]);
				stats[trt] = tt.getT();
				pVals[trt] = tt.getPvalue();
				nQdata[trt][0] = tt.getStdev();
				nQdata[trt][1] = aggregateCounts[trt][1];
				nQdata[trt][2] = aggregateCounts[trt][0];
			} else {
				System.err.println("Error - ANOVA is not yet implemented");
				// anova
			}
		}
	}

	public static PermuteOnePer procStrings(String[] data, String[][] varNames) {
		String[] line;
		int numTraits = varNames[1].length;
		String[] famIDs = new String[data.length];
		double[] trait = new double[data.length];
		int[][] variables = new int[data.length][numTraits];

		for (int i = 0; i<data.length; i++) {
			line = data[i].split("[\\s]+");
			if (line.length!=numTraits+2) {
				System.err.println("Error - line does not have the proper number of columns");
			}
			famIDs[i] = line[0];
			trait[i] = Double.parseDouble(line[1]);
			for (int j = 0; j<numTraits; j++) {
				variables[i][j] = procInt(line[2+j]);
			}
		}

		return new PermuteOnePer(famIDs, trait, variables, varNames);
	}

	public static int procInt(String str) {
		if (str.equals(".")) {
			return Integer.MIN_VALUE;
		}
		try {
			return Integer.parseInt(str);
		} catch (Exception e) {
			System.err.println("Error - '"+str+"' is an invalid number");
			return Integer.MIN_VALUE;
		}
	}

	public static void loadFile(String filename) {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String temp;
		String[] line;
		String[][] varNames = null;
		Vector<String> v = new Vector<String>();
		PermuteOnePer pop;

		try {
			reader = new BufferedReader(new FileReader(filename));
			line = reader.readLine().split("\t");
			varNames = new String[][] {new String[] {line[1]}, new String[line.length-2]};
			for (int i = 0; i<line.length-2; i++) {
				varNames[1][i] = line[2+i];
			}
			if (!line[0].equals("FamID")) {
				System.err.println("Warning - the first column was not called 'FamID'; still assuming that this column contains the grouping variable and that the second contains the trait to be permuted");
			}
			while (reader.ready()) {
				temp = reader.readLine();
				line = temp.split("\t");
				if (line[0].equals(".")) {
					System.err.println("Warning - "+"FamID"+" is missing for an individual");
				} else if (line[1].equals(".")) {
					System.err.println("Warning - "+varNames[0][0]+" is missing for an individual in "+line[0]);
				} else {
					v.add(temp);
				}
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
			pop = procStrings(Array.toStringArray(v), varNames);
			writer = new PrintWriter(new FileWriter(filename+"-summary.out"));
			writer.println(pop.getSummary());
			writer.println();
			writer.println();
			writer.println(pop.getNQueryInput());
			writer.close();
			double[][] results = pop.getResults()[0];
			for (int i = 0; i<results.length; i++) {
				System.out.println(Array.toStr(results[i]));
			}
		} catch (IOException ioe) {
			System.err.println("Error writing file \""+filename+"-summary.out"+"\"");
			System.exit(2);
		}
	}

	public String getSummary() {
		String str = "";

		str += "Variable in question is '"+variableNames[0][0]+"'\n";
		str += "\n";
		for (int trt = 0; trt<means.length; trt++) {
			str += "Trait";
			for (int i = offsets[trt]; i<offsets[trt]+counts[trt]; i++) {
				str += "\tscore for "+i+"s (num)";
			}
			str += "\t[Diff]\tp-value"+"\n";
			str += variableNames[1][trt];
			for (int i = 0; i<means[trt].length; i++) {
				str += "\t"+ext.formDeci(means[trt][i], SIGFIGS, true)+" +- "+ext.formDeci(stdevs[trt][i], SIGFIGS, true)+" ("+ext.formDeci(aggregateCounts[trt][i], 3, true)+")";
			}
			str += "\t"+(counts[trt]==2?ext.formDeci(means[trt][1]-means[trt][0], SIGFIGS, true):"")+"\t"+ext.prettyP(pVals[trt])+"\n";
			str += "\n";
		}

		return str;
	}

	public double[][][] getResults() {
		double[][][] results = new double[numTraits][][];
		double[] allCounts = new double[numTraits];

		for (int trt = 0; trt<means.length; trt++) {
			results[trt] = new double[counts[trt]+2][3];

			for (int i = 0; i<means[trt].length; i++) {
				results[trt][i][0] = means[trt][i];
				results[trt][i][1] = stdevs[trt][i];
				results[trt][i][2] = aggregateCounts[trt][i];
				allCounts[trt] += aggregateCounts[trt][i];
			}
			results[trt][counts[trt]] = new double[] {pVals[trt], stats[trt], allCounts[trt]};
		}

		return results;
	}

	public String getNQueryInput() {
		String str = "";

		for (int i = 0; i<13; i++) {
			for (int trt = 0; trt<means.length; trt++) {
				str += trt==0?"":"\t";
				switch (i) {
				case 0:
					str += variableNames[1][trt];
					break;
				case 1:
					str += "0.010";
					break;
				case 2:
					str += "2";
					break;
				case 6:
					str += nQdata[trt][0];
					break;
				case 8:
					str += "80";
					break;
				case 9:
					str += Math.round(nQdata[trt][1]);
					break;
				case 10:
					str += Math.round(nQdata[trt][2]);
					break;
				default:
					str += "";
					break;
				}
			}
			str += "\n";
		}

		return str;
	}

	public static double[][] permuteMeans(Hashtable<String,double[][]> hash, int spreadSize) {
		double[] means, aggregateSums, aggregateCounts;
		String[] fams = HashVec.getKeys(hash);
		double[][] data;
		double sum;

		aggregateSums = new double[spreadSize];
		aggregateCounts = new double[spreadSize];
		for (int i = 0; i<fams.length; i++) {
			data = hash.get(fams[i]);
			sum = 0;
			for (int j = 0; j<data.length; j++) {
				sum += data[j].length;
			}
			for (int j = 0; j<data.length; j++) {
				if (data[j].length>0) {
					aggregateSums[j] += Array.mean(data[j])*(double)data[j].length/sum;
					aggregateCounts[j] += (double)data[j].length/sum;
				}
			}
		}
		means = new double[spreadSize];
		for (int i = 0; i<means.length; i++) {
			means[i] = aggregateSums[i]/aggregateCounts[i];
		}

		return new double[][] {means, aggregateCounts};
	}

	public static double[] permuteStdev(Hashtable<String,double[][]> hash, int spreadSize) {
		double[] means, aggregateSums, aggregateCounts;
		String[] fams = HashVec.getKeys(hash);
		double[][] data;
		double sum;

		data = permuteMeans(hash, spreadSize);
		means = data[0];
		aggregateCounts = data[1];

		aggregateSums = new double[spreadSize];
		for (int i = 0; i<fams.length; i++) {
			data = hash.get(fams[i]);
			sum = 0;
			for (int j = 0; j<data.length; j++) {
				sum += data[j].length;
			}
			for (int j = 0; j<data.length; j++) {
				for (int k = 0; k<data[j].length; k++) {
					aggregateSums[j] += Math.pow(data[j][k]-means[j], 2)/sum;
				}
			}
		}
		for (int i = 0; i<aggregateSums.length; i++) {
			aggregateSums[i] = Math.sqrt(aggregateSums[i]/(aggregateCounts[i]-1));
		}

		return aggregateSums;
	}

	public static int countEm(DoubleVector[] dvs) {
		int sum = 0;

		for (int i = 0; i<dvs.length; i++) {
			sum += dvs[i].size();
		}

		return sum;
	}

	public Hashtable<String,double[][]> convertHashtable(Hashtable<String,DoubleVector[][]> hash, int col) {
		Hashtable<String,double[][]> newHash = new Hashtable<String,double[][]>();
		String[] fams = HashVec.getKeys(hash);
		DoubleVector[][] sorted;
		double[][] matrix;

		for (int i = 0; i<fams.length; i++) {
			sorted = hash.get(fams[i]);
			matrix = new double[sorted[col].length][];
			for (int j = 0; j<sorted[col].length; j++) {
				matrix[j] = sorted[col][j].toArray();
			}
			newHash.put(fams[i], matrix);
		}

		return newHash;
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		// String filename = "oneperttest.txt";
		// String filename = "oneperttest_VPD_nottreated.txt";
		// String filename = "GBA_Affected.txt";
		// String filename = "GBA_2mut_AOO.txt";
		// String filename = "GBA 2mut Affected unbiased (no NOEV sibs).txt";
		// String filename = "GBA 2+2mut Affected unbiased (no NOEV sibs).txt";
		// String filename = "Rep1_AOO_263carrier.txt";
		// String filename = "AOO-dumpAll.xls";
		// String filename = "Affected_cases.txt";
		// String filename = "Affected_controls.txt";
		// String filename = "permutedGenderCases.txt";
		// String filename = "permutedGenderControls.txt";

		// String filename = "Male.dat";
		// String filename = "AOO.dat";
		String filename = "EarlyOnset.dat";
		// String filename = "DurationFromAOO.dat";
		// String filename = "VPD.dat";
		// String filename = "AffParent.dat";
		// String filename = "Depression.dat";
		// String filename = "MMSE.dat";
		// String filename = "UPDRSmotor.dat";
		// String filename = "UPDRSliving.dat";
		// String filename = "BlessedFunctionality.dat";
		// String filename = "Education.dat";
		// String filename = "SchwabExaminer.dat";
		// String filename = "Hoehn&Yahr.dat";
		// String filename = "PDopinion.dat";
		// String filename = "PD_GTE90.dat";

		String usage = "\n"+"park.onePerTTest requires 0-1 arguments\n"+"   (1) filename (i.e. file="+filename+" (default)\n"+"       the file should be tab delimited and should contain 3+ columns (the grouping variable (FamID), the trait to be permuted, and one or more secondary groups (i.e. something like carrier status)\n"+"";

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
			loadFile(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
