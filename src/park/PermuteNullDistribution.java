package park;

import java.io.*;
import java.util.*;

import stats.LeastSquares;
import stats.LogisticRegression;
import stats.RegressionModel;

import common.*;

public class PermuteNullDistribution {
	public static final int NUM_NULL_REPS = 1000;
	public static final int NUM_ONEPER_REPS = 1000;
	public static final int NUM_BOOT_REPS = 100;
	// public static final int NUM_NULL_REPS = 10;
	// public static final int NUM_ONEPER_REPS = 10;
	// public static final int NUM_BOOT_REPS = 10;

	public static void permuteNullDistribution(String filename, int numReps, int perFamReps, int bootReps) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;
		String temp, size;
		Hashtable<String,Vector<String[]>> hash = new Hashtable<String,Vector<String[]>>();
		Hashtable<String,Vector<String>> countHash;
		Vector<String[]> vStringArray;
		Vector<String> vString;
		Vector<String> pheno;
		String[] fams, sizes, famIDs;
		double[][][] depMatrix; // fam size, fams, outcome phenotype of fam
		// members
		double[][][][] indepMatrix; // fam size, fams, fam members, member's
		// independent factors
		int numCols = -1;
		double[] deps;
		double[][] indeps;
		int count, n;
		int[] randomFams, randomPeople;
		RegressionModel model;
		boolean logistic;
		double[] stats;

		try {
			reader = new BufferedReader(new FileReader(filename));
			line = reader.readLine().split("\t");
			numCols = line.length;
			if (!line[0].equals("FamID")) {
				System.err.println("FYI, assuming the first column is FamID (the header reads '"+line[0]+"')");
			}
			while (reader.ready()) {
				temp = reader.readLine();
				line = temp.split("[\\s]+");
				if (line.length!=numCols) {
					System.err.println("Error - different number of columns starting at line:\n"+temp);
					System.exit(1);
				}
				if (ext.indexOfStr(".", line)==-1) {
					if (hash.containsKey(line[0])) {
						vStringArray = hash.get(line[0]);
					} else {
						hash.put(line[0], vStringArray = new Vector<String[]>());
					}
					vStringArray.add(line);
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

		countHash = new Hashtable<String,Vector<String>>();
		fams = HashVec.getKeys(hash);
		for (int i = 0; i<fams.length; i++) {
			size = hash.get(fams[i]).size()+"";
			if (countHash.containsKey(size)) {
				vString = countHash.get(size);
			} else {
				countHash.put(size, vString = new Vector<String>());
			}
			vString.add(fams[i]);
		}

		sizes = Array.toStr(Sort.putInOrder(Array.toIntArray(HashVec.getKeys(countHash))), "\t").split("[\\s]+");
		depMatrix = new double[sizes.length][][];
		indepMatrix = new double[sizes.length][][][];
		pheno = new Vector<String>();
		count = 0;
		for (int i = 0; i<sizes.length; i++) {
			fams = Array.toStringArray(countHash.get(sizes[i]));
			if (fams.length<5) {
				System.err.println("Warning - There are only "+fams.length+" families with a family size of "+sizes[i]+".");
			}
			depMatrix[i] = new double[fams.length][];
			indepMatrix[i] = new double[fams.length][][];
			for (int j = 0; j<fams.length; j++) {
				vStringArray = hash.get(fams[j]);
				depMatrix[i][j] = new double[vStringArray.size()];
				indepMatrix[i][j] = new double[vStringArray.size()][];
				for (int k = 0; k<vStringArray.size(); k++) {
					line = vStringArray.elementAt(k);
					try {
						depMatrix[i][j][k] = Double.parseDouble(line[1]);
						HashVec.addIfAbsent(depMatrix[i][j][k]+"", pheno);
						indepMatrix[i][j][k] = new double[line.length-2];
						for (int ins = 0; ins<line.length-2; ins++) {
							indepMatrix[i][j][k][ins] = Double.parseDouble(line[2+ins]);
						}
						count++;
					} catch (NumberFormatException nfe) {
						System.err.println("Error - could not parse one of the following numbers:\n"+Array.toStr(line));
					}
				}
			}
		}

		n = count;
		deps = new double[n];
		indeps = new double[n][];
		famIDs = new String[n];
		if (pheno.size()<2) {
			System.err.println("Error - no variation in the dependent variable; exiting");
			System.exit(1);
		}
		logistic = pheno.size()==2;

		System.out.println(ext.getTime());
		writer = new PrintWriter(new FileWriter(filename+".out"));
		writer.println("Replicate\tWald\tp-value");
		for (int rep = -1; rep<numReps; rep++) {
			if (rep%1==0) {
				System.out.println("Rep: "+(rep+1));
			}
			count = 0;
			for (int i = 0; i<depMatrix.length; i++) {
				randomFams = Array.random(depMatrix[i].length);
				for (int j = 0; j<depMatrix[i].length; j++) {
					randomPeople = Array.random(depMatrix[i][j].length);
					for (int k = 0; k<depMatrix[i][j].length; k++) {
						deps[count] = depMatrix[i][j][k];
						if (rep==-1) {
							indeps[count] = indepMatrix[i][j][k];
						} else {
							indeps[count] = indepMatrix[i][randomFams[j]][randomPeople[k]];
						}
						famIDs[count] = i+""+j;
						count++;
					}

				}
			}
			model = logistic?new LogisticRegression(deps, indeps, false, false):new LeastSquares(deps, indeps, false, false);
			model.onePerFamily(famIDs, perFamReps*(rep==-1?10:1), bootReps*(rep==-1?100:1));
			stats = model.getStats();
			for (int i = 1; i<stats.length; i++) {
				writer.print((i==1?"":"\t")+stats[i]);
			}
			writer.println();
			writer.flush();
		}
		writer.close();
		System.out.println(ext.getTime());

	}

	public static void batchAll(String filename, int numFiles, int numReps, int perFamReps, int bootReps) {
		BufferedReader reader;
		PrintWriter[] writers;
		String[] line;
		int increment;

		try {
			reader = new BufferedReader(new FileReader(filename));
			line = reader.readLine().split("[\\s]+");
			increment = (int)Math.floor((double)(line.length-2)/(double)numFiles);
			writers = new PrintWriter[numFiles];
			for (int i = 0; i<numFiles; i++) {
				writers[i] = new PrintWriter(new FileWriter("batch."+(i+1)));
				writers[i].println("mkdir "+(i+1));
				for (int j = 0; j<line.length-2; j++) {
					writers[i].println("cp "+line[2+((j+i*increment)%(line.length-2))]+".dat "+(i+1));
					writers[i].println("cd "+(i+1));
					writers[i].println("java -classpath /home/npankrat/park.jar park.PermuteNullDistribution file="+line[2+((j+i*increment)%(line.length-2))]+".dat reps="+numReps+" perFamReps="+perFamReps+" bootReps="+bootReps);
					writers[i].println("cd ..");
				}
				writers[i].close();
			}

			writers = new PrintWriter[line.length-2];
			for (int i = 0; i<writers.length; i++) {
				writers[i] = new PrintWriter(new FileWriter(line[2+i]+".dat"));
				writers[i].println(line[0]+"\t"+line[1]+"\t"+line[2+i]);
			}
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				for (int i = 0; i<writers.length; i++) {
					writers[i].println(line[0]+"\t"+line[1]+"\t"+line[2+i]);
				}
			}
			reader.close();
			for (int i = 0; i<writers.length; i++) {
				writers[i].close();
			}
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}

	}

	public static void gatherbatchedData(String filename, int numFiles) {
		BufferedReader reader;
		PrintWriter writer, grand;
		String[] traits;
		double observed;
		int countAbove, count, grandCountAbove, grandCount;

		try {
			reader = new BufferedReader(new FileReader(filename));
			traits = Array.subArray(reader.readLine().split("[\\s]+"), 2);
			reader.close();

			writer = Files.getWriter(filename+"-permutedSummary.out");
			grand = Files.getWriter(filename+"_p-values.out");
			writer.println("Trait\t"+Array.toStr(Array.stringArraySequence(numFiles, "Rep "), "\t\t\t"));
			for (int i = 0; i<traits.length; i++) {
				writer.print(traits[i]);
				grandCountAbove = grandCount = 0;
				for (int j = 0; j<numFiles; j++) {
					try {
						reader = new BufferedReader(new FileReader((j+1)+"/"+traits[i]+".dat.out"));
						count = countAbove = 0;
						if (reader.ready()) {
							reader.readLine();
							observed = Math.abs(Double.parseDouble(reader.readLine().split("[\\s]+")[0]));
						} else {
							observed = -1;
						}
						while (reader.ready()) {
							if (Math.abs(Double.parseDouble(reader.readLine().split("[\\s]+")[0]))>=observed) {
								countAbove++;
								grandCountAbove++;
							}
							count++;
							grandCount++;
						}
						reader.close();
						writer.print("\t"+countAbove+"\t"+count+"\t"+(count==0?".":ext.prettyP((double)countAbove/(double)count)));
					} catch (FileNotFoundException fnfe) {
						System.err.println("File \""+(j+1)+"/"+traits[i]+".dat.out"+"\" not found in current directory");
						writer.print("\t.\t.\t.");
					} catch (IOException ioe) {
						System.err.println("Error reading file \""+(j+1)+"/"+traits[i]+".dat.out"+"\"");
						writer.print("\t.\t.\t.");
					}
				}
				writer.println();

				grand.println(traits[i]+"\t"+grandCountAbove+"\t"+grandCount+"\t"+(grandCount==0?".":ext.prettyP((double)grandCountAbove/(double)grandCount)));
			}
			writer.close();
			grand.close();
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
		// String filename = "permuteQuantNull.dat";
		// String filename = "permuteQuantNull(both).dat";
		// String filename = "UPDRSlivingCount.dat";
		// String filename = "UPDRSmotorCount.dat";
		// String filename = "UPDRSmotorCarrier.dat";
		// String filename = "UPDRSlivingCarrier.dat";
		// String filename = "UPDRSmotorOnly.dat";
		String filename = "all.dat";
		int batch = 0;
		boolean gather = false;

		int numReps = NUM_NULL_REPS;
		int oneperReps = NUM_ONEPER_REPS;
		int bootReps = NUM_BOOT_REPS;

		String usage = "\n"+"park.permuteQuantNull requires 0-1 arguments\n"+"   (1) filename (i.e. file="+filename+" (default)\n"+"   (2) num replicates to create null distribution (i.e. reps="+numReps+" (default)\n"+"   (3) num replicates for oneperfamily calculations (i.e. perFamReps="+oneperReps+" (default)\n"+"   (4) num replicates to bootstrap (i.e. bootReps="+bootReps+" (default)\n"+"   (5) batch separately N times (as opposed to treat as covariates) (batch>0  i.e. batch="+batch+" (default)\n"+"   (6) gather batched results (requires filename for reference i.e. gather="+gather+" (default)\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("reps=")) {
				numReps = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("perFamReps=")) {
				oneperReps = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("bootReps=")) {
				bootReps = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("batch=")) {
				batch = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("gather=")) {
				gather = args[i].split("=")[1].equalsIgnoreCase("true");
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (gather) {
				gatherbatchedData(filename, batch);
			} else if (batch>0) {
				batchAll(filename, batch, numReps, oneperReps, bootReps);
			} else {
				permuteNullDistribution(filename, numReps, oneperReps, bootReps);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
