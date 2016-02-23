// mandating uniqueid, famid, indid as the first three is limiting
// -Xms1024M -Xmx1024M
package db;

import java.io.*;
import java.util.*;

import mining.Transformations;
import common.AlleleFreq;
import common.Array;
import common.DoubleVector;
import common.HashVec;
import common.IntVector;
import common.Logger;
import common.Matrix;
import common.Sort;
import common.ext;
import stats.ContingencyTable;
import stats.LeastSquares;
import stats.LogisticRegression;
import stats.PermuteOnePer;
import stats.ProbDist;
import stats.RegressionModel;
import stats.Stepwise;
import stats.Ttest;

public class comp {
	public static String DEFAULT_TRAIT = "AOO";
	public static String[][] DEFAULT_ID_NAMES = {{"UniqueID", "UID"}, {"FID", "FamID"}, {"IID", "IndID"}};
	public static String DEFAULT_INPUT = "linear.ctl";
	public static String DEFAULT_DB = "\"C:\\Documents and Settings\\npankrat\\My Documents\\1_CRFdb\\crf_db.dat\"";
	public static String[] DEFAULT_USE = {"BirthDate", "AgeAtOnset", "AgeAtExam", "Duration", "AOO", "DurationFromAOO", "EarlyOnset45", "EarlyOnset50", "EarlyOnset60", "VPD", "CONF_PD", "Male", "Caucasian", "Hispanic", "AffFather", "AffMother", "AffParent", "parkin", "polymorph", "G2019S", "Depression", "Depressed", "DSMIV_Dx_Depression", "MajorDepression", "MinorDepression", "MMSE", "Demented", "BlessedFunctionality", "Education", "MilitaryYears", "smoked", "smoker", "alcohol", "pesticides", "OnsetWithTremor", "DominantSideFirst", "LeftSideFirst", "RightSideFirst", "BothSidesFirst", "HallucinationWithDrugs", "HallucinationsWithoutDrugs", "depressionBefore", "depressionSince", "HeadInjury", "Infection", "SchwabExaminer", "SchwabSubject", "SchwabDiff", "Hoehn&Yahr", "Bradykinesia", "Rigidity", "Instability", "PersistentAsymmetry", "RestTremor", "ProgressiveDisorder", "levodopaChorea", "levodopa5PlusYears", "Course10+", "UnexplainedMotor", "Strokes", "Encephalitis", "OculogyricCrisis", "Alzheimers", "SensoryDeficits_Apraxia", "PDfromDopaDrugs", "Remission", "unilateral3PlusYears", "SupranuclearGaze", "CerebellarSigns", "Hypotension", "NoResponseLDopa", "lesionMRI", "ldopaResponse", "PDopinion", "PD>90%", "logisticE4", "logisticE2", "APOE4count", "UPDRSliving", "UPDRSmotor", "BradykinesiaScore", "BradykinesiaSubScore", "RigiditySubScore", "PIGD_score", "Tremor_score", "PIGD_scale", "PIGD_dominant", "PIGD_intermediate", "Tremor_dominant", "SpeechScore", "RestTremorScore", "ActionTremorScore", "CombinedTremorScore"};
	public static double SIG_CUTOFF = 0.05;
	public static int NUM_THREADS_DEFAULT = 1;
	public static int FAM_REPS_DEFAULT = 10000;
	public static int BOOT_REPS_DEFAULT = 10000;
	public static final String[] OPTIONS = {"dump", "dumpAll", "sw", "allsw", "predicteds", "residuals", "normalized", "inverseNormalized", "exactRegressionValues", "table", "sdtable", "trend", "oneperfamily", "verbose", "force", "noserialperm", "chis", "audit", "hwe"};
	public static final int MAX_CLASSES = 15;
	public static final int DEFAULT_SIG_FIGS = 3;
	public static final int SIG_FIGS_PERCENTAGES = 1;

	private boolean[] flags = new boolean[OPTIONS.length];

	public comp() throws IOException {
		PrintWriter writer = new PrintWriter(new FileWriter(DEFAULT_INPUT));
		writer.println(DEFAULT_DB);
		writer.println(DEFAULT_TRAIT);
		writer.println(Array.toStr(DEFAULT_USE));
		writer.println();
		writer.println("user notes: ");
		writer.close();

		new comp(DEFAULT_INPUT);
	}

	@SuppressWarnings("resource")
	public comp(String filename) throws IOException {
		BufferedReader reader, r2;
		PrintWriter writer = null;
		String[] line, traits, split;
		String temp;
		String db_file, notes;
		Vector<String> vString, deps, depCount, included = new Vector<String>();
		Vector<String[]> idV;
		String[][] ids;
		Vector<double[]> vDoubleArray;
		Vector<double[]> indeps;
		int count, N, M, numSig, n, k;
		int[] indices, limits, idIndices, sigSpots, factorNs, counts, failures;
		String[] factorNames, sigNames, validNames, idline;
		double[] dataline;
		double[][] data, sigs, means;
		RegressionModel model;
		Stepwise sw;
		boolean logistic;
		double[] predicteds = null;
		double[][] residuals = null;
		Vector<String> limitKeys = new Vector<String>(), limitTargets = new Vector<String>();
		boolean[] factorDirections;
		double[][] effectsAndConfidenceIntervals;
		DoubleVector[] trends;
		Hashtable<String,String> trendyKey;
//		Vector<String> famIDs;
		int famReps = FAM_REPS_DEFAULT;
		int bootReps = BOOT_REPS_DEFAULT;
		// int numThreads = NUM_THREADS_DEFAULT;
		boolean noMissing;
		double[][] results;
		boolean percentMe;
		int sigfigs = DEFAULT_SIG_FIGS;
		double[] finalDeps;
		double[][] finalIndeps;
		int[] order;
		String suffix = "";
		Logger log;
		String delimiter;
		
		log = new Logger(ext.rootOf(filename, false)+".log");
		if (!new File(filename).exists()) {
			log.reportError("Error - "+filename+" does not exist");
			return;
		}
		
		residuals = new double[3][]; // raw, normalized, inverseNormalized


		new File(filename+".bak").delete();
		new File(filename).renameTo(new File(filename+".bak"));
		try {
			reader = new BufferedReader(new FileReader(filename+".bak"));
		} catch (FileNotFoundException ex) {
			throw new RuntimeException(ex.getMessage());
		}
		
		temp = reader.readLine();
		line = temp.split("[\\s]+");
		db_file = line[0];
		delimiter = db_file.endsWith(".csv")?",":"\t";
		
		if (db_file.startsWith("\"")) {
			db_file = temp.substring(1, temp.substring(1).indexOf("\"")+1); // why is this +1 at the end? seems like it should be +0
			line = temp.substring(temp.substring(1).indexOf("\"")+1).split("[\\s]+");
		}
		order = null;
		for (int i = 1; i<line.length; i++) {
			if (line[i].indexOf("=")>0) {
				if (line[i].split("=")[0].equals("famreps")) {
					famReps = Integer.parseInt(line[i].split("=")[1]);
				} else if (line[i].split("=")[0].equals("sigfigs")) {
					sigfigs = Integer.parseInt(line[i].split("=")[1]);
				} else if (line[i].split("=")[0].equals("bootreps")) {
					bootReps = Integer.parseInt(line[i].split("=")[1]);
					// } else if (line[i].split("=")[0].equals("threads")) {
					// numThreads = Integer.parseInt(line[i].split("=")[1]);
				} else if (line[i].split("=")[0].equals("order")) {
					split = line[i].split("=")[1].split(",");
					order = new int[split.length];
					for (int j = 0; j < split.length; j++) {
						try {
							order[j] = Integer.parseInt(split[j]);
						} catch (Exception e) {
							System.err.println("Error - invalid order= paramter: "+line[i]);
							System.exit(1);
						}
					}
				} else {
					limitKeys.add(line[i].split("=")[0]);
					limitTargets.add(line[i].split("=")[1]);
				}
			} else {
				flagOption(line[i].toLowerCase());
			}
		}
		if (!optionFlagged("sw")) {
			flagOption("allsw", false);
		}
		if (reader.ready()) {
			traits = reader.readLine().split("\t");
			for (int i = 0; i < traits.length; i++) {
				if (traits[i].startsWith("suffix=")) {
					suffix = ext.parseStringArg(traits[i], "");
					traits = Array.removeFromArray(traits, i);
					i--;
				}
			}
			line = reader.readLine().split("\t");
			notes = reader.ready()?reader.readLine():"user notes:";
			for (int i = 0; i<line.length; i++) {
				if (!line[i].equals("")) {
					included.add(line[i]);
				}
			}
		} else {
			writer = new PrintWriter(new FileWriter(filename));
			writer.println(db_file);
			writer.println("Affected");
			try {
				writer.println(new BufferedReader(new FileReader(db_file)).readLine());
			} catch (FileNotFoundException ex) {
				reader.close();
				writer.close();
				throw new RuntimeException(ex.getMessage());
			}
			writer.println();
			writer.println();
			writer.print("Options not flagged: ");
			for (int i = 0; i<OPTIONS.length; i++) {
				if (!optionFlagged(OPTIONS[i])) {
					writer.print("\t"+OPTIONS[i]);
				}
			}
			writer.println();
			writer.close();
			traits = new String[0];
			notes = "";
		}
		reader.close();
		noMissing = optionFlagged("chis")||optionFlagged("hwe"); // not sure quite what this did when I added hwe
		for (int trt = 0; trt<traits.length; trt++) {
			try {
				reader = new BufferedReader(new FileReader(db_file));
			} catch (FileNotFoundException ex) {
				reader.close();
				throw new RuntimeException(ex.getMessage());
			}
			factorNames = reader.readLine().split(delimiter);

			indices = ext.indexFactors(Array.addStrToArray(traits[trt], Array.toStringArray(included), 0), factorNames, true, log, true, true);
			limits = ext.indexFactors(Array.addStrToArray(traits[trt], Array.toStringArray(limitKeys), 0), factorNames, true, log, true, true);
			idIndices = ext.indexFactors(DEFAULT_ID_NAMES, factorNames, false, true, false, false);
			M = indices.length-1;

			if (optionFlagged("oneperfamily") && idIndices[1] == -1) {
				System.err.println("Error - OnePerFamily option was enabled, but the default FamID name ('"+DEFAULT_ID_NAMES[1]+"') was not found in the database; disabling option");
				flagOption("oneperfamily", false);
			}

			if ((optionFlagged("chis")||optionFlagged("hwe")) && Array.min(idIndices) < 0) {
				System.err.println("Error - chis option was enabled, but family and indiviudal ids were not found; disabling option");
				flagOption("chis", false);
				flagOption("hwe", false);
			}

			System.out.println("Analyzing "+traits[trt]);
			// System.out.println("Loading data into memory");
			vDoubleArray = new Vector<double[]>();
			idV = new Vector<String[]>();
			while (reader.ready()) {
				line = reader.readLine().split(delimiter);
				if (factorNames.length!=line.length) {
					reader.close();
					throw new RuntimeException("Error - different number of values ("+line.length+" versus "+factorNames.length+" factors) for "+line[0]);
				}
				idline = new String[] {idIndices[0]==-1?null:line[idIndices[0]], idIndices[1]==-1?null:line[idIndices[1]], idIndices[2]==-1?null:line[idIndices[2]]};
				dataline = new double[M+1];
				limits[0] = 1;
				for (int i = 1; i<limits.length; i++) {
					if (!line[limits[i]].equals(limitTargets.elementAt(i-1))) {
						limits[0] = -1;
					}
				}
				for (int i = 0; i<M+1; i++) {
//					if (Double.isNaN(dataline[i])) { // important, "NaN" doesn't parse as Double.NaN for some reason
					if (ext.isMissingValue(line[indices[i]])) {
						dataline[i] = Double.MIN_VALUE;
					} else {
						try {
							dataline[i] = Double.parseDouble(line[indices[i]]);
						} catch (NumberFormatException nfe) {
							// System.out.println(line[indices[i]]);
							dataline[i] = Double.MIN_VALUE;
							if (noMissing) {
								limits[0] = -1;
							}
						}
					}
				}
				if (dataline[0]!=Double.MIN_VALUE&&limits[0]!=-1) {
					vDoubleArray.add(dataline);
					idV.add(idline);
				}
			}
			N = vDoubleArray.size();
			data = new double[N][];
			ids = Matrix.toStringArrays(idV);
			idV.removeAllElements();
			for (int i = 0; i<N; i++) {
				data[i] = vDoubleArray.elementAt(i);
			}
			if (N==1) {
				reader.close();
				throw new RuntimeException("Error - with only one subject, dependent variable has zero variation");
			}
			vDoubleArray = null;

			depCount = new Vector<String>();
			count = 0;
			while (count<N&&depCount.size()<=MAX_CLASSES) {
				HashVec.addIfAbsent(data[count][0]+"", depCount);
				count++;
			}
			k = depCount.size();
			logistic = (k==2);

			if (optionFlagged("chis")) {
				if (M!=2) {
					System.err.println("Error - \"chis\" is flagged, but there are "+M+" variables instead of 2 (the two alleles); aborting analyses");
				} else if (!logistic) {
					System.err.println("Error - \"chis\" is flagged, but the dependent variable is not binary");
				} else {
					doChis(ids, data, factorNames, indices, optionFlagged("audit"));
				}
			}

			if (optionFlagged("hwe")) {
				if (M!=2) {
					System.err.println("Error - \"hwe\" is flagged, but there are "+M+" variables instead of 2 (the two alleles); aborting analyses");
				} else if (!logistic) {
					System.err.println("Error - \"hwe\" is flagged, but the dependent variable is not binary");
				} else {
					doHWE(ids, data, factorNames, indices, optionFlagged("audit"));
				}
			}

			System.out.println("Performing serial "+(logistic?"logistic":"linear")+" regressions for variable '"+traits[trt]+"'");
			sigs = new double[M+1][3];
			means = new double[M+1][2];
			factorNs = new int[M+1];
			factorDirections = new boolean[M+1];
			failures = new int[M+1];
			effectsAndConfidenceIntervals = new double[M+1][];
			for (int factor = 1; factor<=M; factor++) {
				try {
					deps = new Vector<String>(400000);
					indeps = new Vector<double[]>(400000);
					idV = new Vector<String[]>(400000);
					for (int i = 0; i<N; i++) {
						if (data[i][factor]!=Double.MIN_VALUE) {
							idV.add(ids[i]);
							deps.add(data[i][0]+"");
							indeps.add(new double[] {(data[i][factor])});
						}
					}
					model = logistic?(RegressionModel)new LogisticRegression(deps, indeps, false, optionFlagged("verbose")):(RegressionModel)new LeastSquares(deps, indeps, false, optionFlagged("verbose"));
					if (optionFlagged("oneperfamily")&&!optionFlagged("noserialperm")) {
						Date date = new Date();
						System.out.print("Running "+factorNames[indices[factor]]+"... ");
						model.onePerFamily(Matrix.extractColumn(Matrix.toStringArrays(idV), 1), famReps, bootReps);
						System.out.println(ext.getTimeElapsed(date.getTime()));
						failures[factor] = model.getNumFailures();
					}
					if (model.analysisFailed()) {
						System.err.println("Error - analysis failed for "+included.elementAt(factor-1));
						sigs[factor][0] = 99;
						sigs[factor][1] = 0;
						effectsAndConfidenceIntervals[factor] = new double[] {0,0,0};
					} else {
						sigs[factor][0] = model.getSigs()[1];
						sigs[factor][1] = model.getRsquare();
						sigs[factor][2] = model.getStats()[1];
						factorNs[factor] = deps.size();
						factorDirections[factor] = model.getBetas()[1]>0;
						effectsAndConfidenceIntervals[factor] = model.getEffectsAndConfidenceIntervals()[1];
					}
					if (logistic) {
						finalDeps = model.getFinalDependentVariables();
						finalIndeps = model.getFinalIndependentVariables();
						counts = new int[2];
						for (int i = 0; i < finalDeps.length; i++) {
							means[factor][(int)finalDeps[i]] += finalIndeps[i][0];
							counts[(int)finalDeps[i]]++;
						}
						for (int i = 0; i < 2; i++) {
							means[factor][i] /= (double)counts[i];
						}
					}
				} catch (Exception e) {
					e.printStackTrace();
					sigs[factor][0] = 99;
					sigs[factor][1] = 0;
				}
			}

			if (optionFlagged("dumpAll")) {
				writer = new PrintWriter(new FileWriter(traits[trt]+suffix+"-dumpAll.xln"));
				writer.print("UniqueID\tFamID\tIndID");
				for (int j = 0; j<M+1; j++) {
					writer.print("\t"+factorNames[indices[j]]);
				}
				writer.println();
				for (int i = 0; i<N; i++) {
					writer.print(Array.toStr(ids[i]));
					for (int j = 0; j<M+1; j++) {
						writer.print("\t"+(data[i][j]==Double.MIN_VALUE?".":data[i][j]+""));
					}
					writer.println();
				}
				writer.close();
			}

			if (optionFlagged("trend")) {
				if (k>MAX_CLASSES) {
					System.err.println("Error - too many dependent states; skipping trend analysis");
				} else {
					try {
						writer = new PrintWriter(new FileWriter(traits[trt]+suffix+"-trend.xln"));
						for (int i = 0; i<k; i++) {
							writer.print("\t"+depCount.elementAt(i));
						}
						writer.println();
						for (int factor = 0; factor<=M; factor++) {
							trends = new DoubleVector[k];
							for (int i = 0; i<trends.length; i++) {
								trends[i] = new DoubleVector();
							}
							for (int i = 0; i<N; i++) {
								if (data[i][factor]!=Double.MIN_VALUE) {
									trends[depCount.indexOf(data[i][0]+"")].add(data[i][factor]);
								}
							}
							for (int i = 0; i<trends.length; i++) {
								if (factor==0) {
									writer.print((i==0?"N":"")+"\t"+trends[i].size());
								} else {
									writer.print((i==0?factorNames[indices[factor]]:"")+"\t"+Array.mean(trends[i].toArray()));
								}
							}
							writer.println();
						}
						writer.close();

						writer = new PrintWriter(new FileWriter(traits[trt]+suffix+"-trend.sas"));
						trends = new DoubleVector[M+1];
						validNames = new String[M+1];
						for (int factor = 0; factor<=M; factor++) {
							trends[factor] = new DoubleVector();
							validNames[factor] = ext.replaceWithLinuxSafeCharacters(factorNames[indices[factor]], true);
						}
						for (int factor = 1; factor<=M; factor++) {
							writer.println("data Rocks;");
							writer.println("   input "+validNames[0]+" "+validNames[factor]+";");
							writer.println("   datalines;");
							for (int i = 0; i<N; i++) {
								writer.println(ext.formDeci(data[i][0], 5)+" "+(data[i][factor]==Double.MIN_VALUE?".":ext.formDeci(data[i][factor], 5)));
								if (data[i][factor]!=Double.MIN_VALUE&&trends[factor].size()<3&&!trends[factor].contains(data[i][factor])) {
									trends[factor].add(data[i][factor]);
								}
							}
							writer.println(";");
							writer.println();

							if (trends[factor].size()<3) {
								writer.println("proc freq data=Rocks;");
								writer.println("   tables "+validNames[0]+"*"+validNames[factor]+" / trend;");
								writer.println("   exact trend / maxtime=60;");
								writer.println("   title1 'Cochran-Armitage Trend Test for "+factorNames[indices[factor]]+"';");
								writer.println("run;");
								writer.println("");
							} else {
								writer.println("proc glm data=Rocks;");
								writer.println("   class "+validNames[0]+";");
								writer.println("   model "+validNames[factor]+" = "+validNames[0]+" / ss3;");
								writer.println("   lsmeans "+validNames[0]+";");
								writer.println("   contrast \"Blessed linear\" "+validNames[0]+" -1 0 1;");
								writer.println("   estimate \"Blessed linear\" "+validNames[0]+" -1 0 1;");
								writer.println("   contrast \"Blessed quadratic\" "+validNames[0]+" 1 -2 1;");
								writer.println("   estimate \"Blessed quadratic\" "+validNames[0]+" 1 -2 1;");
								writer.println("   title1 'GLM Trend Test for "+factorNames[indices[factor]]+"';");
								writer.println("run;");
								writer.println("");
							}
						}
						writer.close();
					} catch (Exception e) {
						System.err.println("Sorry dude, error performing trend analysis");
						e.printStackTrace();
					}
					if (new File(traits[trt]+"-trend.lst").exists()) {
						trendyKey = new Hashtable<String,String>();
						try {
							try {
								if (!new File("trendy_key.xln").exists()) {
									throw new IOException();
								}
								r2 = new BufferedReader(new FileReader("trendy_key.xln"));
								while (r2.ready()) {
									line = reader.readLine().split("\t");
									trendyKey.put(line[0], line[1]);
								}
								r2.close();
							} catch (IOException ioe) {
								System.err.println("Error reading file \""+"trendy_key.xln"+"\"");
							}

							reader = new BufferedReader(new FileReader(traits[trt]+"-trend.lst"));
							writer = new PrintWriter(new FileWriter(traits[trt]+"-trend-results.xln"));
							writer.print(traits[trt]);
							for (int i = 0; i<k; i++) {
								writer.print("\t"+depCount.elementAt(i));
							}
							writer.println("\tn\tp-value\tcalcP\tQuad p-value\tcalc Quad p");
							counts = new int[k];
							for (int i = 0; i<N; i++) {
								counts[(int)data[i][0]]++;
							}
							writer.print("N");
							for (int i = 0; i<k; i++) {
								writer.print("\t"+counts[i]);
							}
							writer.println();
							count = 0;
							temp = "";
							while (reader.ready()) {
								temp = reader.readLine();
								if (temp.indexOf("Trend Test for")>0) {
									line = temp.split("[\\s]+");
									writer.print(trendyKey.containsKey(line[5])?trendyKey.get(line[5]):line[5]);
									count++;
									if (line[1].equals("GLM")) {
										for (int i = 0; i<12; i++)
											temp = reader.readLine();
										if (temp.indexOf("Number of Observations Used")==-1) {
											System.err.println("Error - GLM parsing needs to be updated:");
											System.err.println("'"+temp+"'");
										}
										n = Integer.parseInt(temp.split("[\\s]+")[5]);
										while (temp.indexOf("LSMEAN")==-1)
											temp = reader.readLine();
										reader.readLine();
										for (int i = 0; i<k; i++) {
											line = reader.readLine().split("[\\s]+");
											writer.print("\t"+ext.formDeci(Double.parseDouble(line[2]), 1, true));
										}
										writer.print("\t"+n);
										for (int i = 0; i<16; i++)
											temp = reader.readLine();
										if (temp.indexOf("linear")==-1) {
											System.err.println("Error - GLM parsing needs to be updated:");
											System.err.println("'"+temp+"'");
										}
										line = temp.split("[\\s]+");
										writer.print("\t"+line[6]+"\t"+"=TDIST("+Math.abs(Double.parseDouble(line[3])/Double.parseDouble(line[4]))+", "+(n-1)+", 2)");
										line = reader.readLine().split("[\\s]+");
										writer.print("\t"+line[6]+"\t"+"=TDIST("+Math.abs(Double.parseDouble(line[3])/Double.parseDouble(line[4]))+", "+(n-1)+", 2)");
										writer.println("\t"+"=IF(G"+(count+1)+">I"+(count+1)+", \"*\",\"\")");
									} else if (line[1].equals("Cochran-Armitage")) {
										while (temp.indexOf("---------+--------+--------+")==-1)
											temp = reader.readLine();
										for (int i = 0; i<k; i++) {
											line = reader.readLine().split("[\\s]+");
											writer.print("\t"+ext.formDeci((Double.parseDouble(line[5])/Double.parseDouble(line[7]))*100, 0)+"%");
											for (int j = 0; j<4; j++)
												temp = reader.readLine();
										}
										n = Integer.parseInt(reader.readLine().split("[\\s]+")[4]);
										writer.print("\t"+n);
										while (temp.indexOf("Exact Test")==-1)
											temp = reader.readLine();
										reader.readLine();
										line = reader.readLine().split("[\\s]+");
										writer.println("\t\t"+line[5]);
									} else {
										System.err.println("Error - unknown trend test: "+line[1]);
									}
								}
							}
							writer.close();
							reader.close();
						} catch (Exception e) {
							System.err.println("Sorry yo, error parsing results of trend analysis found in '"+traits[trt]+"-trend.lst"+"'");
							writer.close();
							reader.close();
						}
					}
				}
			}

			if (optionFlagged("table")) {

				if (k>MAX_CLASSES) {
					try {
						DoubleVector dv1s, dv0s;
						boolean failed;

						writer = new PrintWriter(new FileWriter(traits[trt]+suffix+"-table.xln"));
						writer.println("Trait\t1s\tNum1s\t0s\tNum0s\tp-value");
						for (int factor = 1; factor<=M; factor++) {
							dv1s = new DoubleVector();
							dv0s = new DoubleVector();
							counts = new int[2];
							failed = false;
							for (int i = 0; !failed && i<N; i++) {
								if (data[i][factor]!=Double.MIN_VALUE) {
									if (data[i][factor] == 1) {
										dv1s.add(data[i][0]);
										counts[1]++;
									} else if (data[i][factor] == 0) {
										dv0s.add(data[i][0]);
										counts[0]++;
									} else {
										System.err.println("Error - cannot create a table for "+factorNames[indices[0]]+" and "+factorNames[indices[factor]]+" (invalid value: "+data[i][factor]+"; requires 1 or 0)");
										failed = true;
									}
								}
							}
							writer.println(factorNames[indices[factor]]+"\t"+ext.formDeci(Array.mean(dv1s.toArray()), sigfigs, true)+"\t"+counts[1]+"\t"+ext.formDeci(Array.mean(dv0s.toArray()), sigfigs, true)+"\t"+counts[0]+"\t"+ext.prettyP(new Ttest(dv1s.toArray(), dv0s.toArray()).getPvalue()));
						}
						writer.close();

					} catch (Exception e) {
						System.err.println("Error creating table");
						e.printStackTrace();
					}
				} else {
					try {
						DoubleVector dv1, dv2;

						writer = new PrintWriter(new FileWriter(traits[trt]+suffix+"-table.xln"));
						for (int i = 0; i<k; i++) {
							writer.print("\t"+depCount.elementAt(i)+"\tcounts");
						}
						writer.println("\tp-value");
						for (int factor = 0; factor<=M; factor++) {
							idV = new Vector<String[]>();
							dv1 = new DoubleVector();
							dv2 = new DoubleVector();
							vString = new Vector<String>();
							counts = new int[k];
							for (int i = 0; i<N; i++) {
								if (data[i][factor]!=Double.MIN_VALUE) {
									idV.add(ids[i]);
									dv1.add(data[i][factor]);
									HashVec.addIfAbsent(data[i][factor]+"", vString);
									dv2.add(depCount.indexOf(data[i][0]+""));
								}
								counts[depCount.indexOf(data[i][0]+"")]++;
							}
							results = new PermuteOnePer(optionFlagged("oneperfamily")?Matrix.extractColumn(Matrix.toStringArrays(idV), 1):Array.stringArraySequence(idV.size(), "IND"), dv1.toArray(), dummyIntMatrix(dv2.toArray()), new String[][] { {factorNames[indices[0]]}, {factorNames[indices[factor]]}}).getResults()[0];
							for (int i = 0; i<k; i++) {
								percentMe = vString.size()==2;
								writer.print(
										(i==0?(factor==0?"N":factorNames[indices[factor]]):"")+
										"\t=CONCATENATE(\""+(percentMe?
												ext.formDeci(results[i][0]*100, SIG_FIGS_PERCENTAGES, true)+"%"
												:ext.formDeci(results[i][0], sigfigs, true)+(optionFlagged("sdtable")?" ("+ext.formDeci(results[i][1], sigfigs, true)+")":"")
												)+"\")"+
										"\t"+ext.formDeci(results[i][2], sigfigs, true)+
										(i==k-1&&factor!=0?"\t"+ext.prettyP(results[k][0])+(percentMe?"\t=CHIDIST("+Math.abs(results[k][1])+",1)":"\t=TDIST("+Math.abs(results[k][1])+","+(Math.round(results[k][2])-2)+",2)"):""));
							}
							writer.println();
						}
						writer.close();

					} catch (Exception e) {
						System.err.println("Error creating table");
						e.printStackTrace();
					}
				}
			}

			if (optionFlagged("exactRegressionValues")) {
				try {
					writer = new PrintWriter(new FileWriter(traits[trt]+suffix+"-exactRegressionValues.xln"));
					writer.println("Factor\tT\tp-value\tcalcP");
					for (int i = 1; i<=M; i++) {
						writer.println(factorNames[indices[i]]+"\t"+sigs[i][2]+"\t"+sigs[i][0]+(logistic?"\t=CHIDIST("+Math.abs(sigs[i][2])+",1)":"\t=TDIST("+Math.abs(sigs[i][2])+","+(factorNs[i]-2)+",2)")+"\t"+(factorDirections[i]?"+":"-")+"\t"+ext.formStr(ext.formDeci(sigs[i][1]*100, SIG_FIGS_PERCENTAGES, true)+"%", 5)+"");
					}
					writer.close();
				} catch (IOException ioe) {
					System.err.println("Error creating "+traits[trt]+suffix+"-exactRegressionValues.xln");
				}
			}

			if (traits.length>1) {
				if (trt==0) {
					writer = new PrintWriter(new FileWriter(filename));
					writer.print(db_file.split("[\\s]+").length>1?"\""+db_file+"\"":db_file);
					for (int i = 0; i<limitKeys.size(); i++) {
						writer.print("\t"+limitKeys.elementAt(i)+"="+limitTargets.elementAt(i));
					}
					for (int i = 0; i<OPTIONS.length; i++) {
						writer.print(optionFlagged(OPTIONS[i])?"\t"+OPTIONS[i]:"");
					}
					// if (numThreads != NUM_THREADS_DEFAULT) {
					// writer.print("\tthreads="+numThreads);
					// }
					if (sigfigs!=DEFAULT_SIG_FIGS) {
						writer.print("\tsigfigs="+sigfigs);
					}
					if (famReps!=FAM_REPS_DEFAULT) {
						writer.print("\tfamreps="+famReps);
					}
					if (bootReps!=BOOT_REPS_DEFAULT) {
						writer.print("\tbootreps="+bootReps);
					}
					if (order!=null) {
						writer.print("\torder="+Array.toStr(order, ","));
					}
					writer.println();
					writer.println(Array.toStr(traits)+(suffix.equals("")?"":"\tsuffix="+suffix));
					for (int i = 0; i<included.size(); i++) {
						writer.print((i==0?"":"\t")+included.elementAt(i));
					}
					writer.println();
					writer.close();
				}
				writer = new PrintWriter(new FileWriter(traits[trt]+suffix+".ctl"));
			} else {
				writer = new PrintWriter(new FileWriter(filename));
			}
			writer.print(db_file.split("[\\s]+").length>1?"\""+db_file+"\"":db_file);
			for (int i = 0; i<limitKeys.size(); i++) {
				writer.print("\t"+limitKeys.elementAt(i)+"="+limitTargets.elementAt(i));
			}
			for (int i = 0; i<OPTIONS.length; i++) {
				writer.print(optionFlagged(OPTIONS[i])?"\t"+OPTIONS[i]:"");
			}
			// if (numThreads != NUM_THREADS_DEFAULT) {
			// writer.print("\tthreads="+numThreads);
			// }
			if (sigfigs!=DEFAULT_SIG_FIGS) {
				writer.print("\tsigfigs="+sigfigs);
			}
			if (famReps!=FAM_REPS_DEFAULT) {
				writer.print("\tfamreps="+famReps);
			}
			if (bootReps!=BOOT_REPS_DEFAULT) {
				writer.print("\tbootreps="+bootReps);
			}
			if (order!=null) {
				writer.print("\torder="+Array.toStr(order, ","));
			}
			
			writer.println();
			writer.println(traits[trt]);
			for (int i = 0; i<included.size(); i++) {
				writer.print((i==0?"":"\t")+included.elementAt(i));
			}
			writer.println();
			writer.println(notes);
			writer.println();
			writer.print("Not used: ");
			for (int i = 1; i<factorNames.length; i++) {
				if (!included.contains(factorNames[i])) {
					writer.print("\t"+factorNames[i]);
				}
			}
			writer.println();
			writer.print("significant: ");
			numSig = 0;
			for (int i = 1; i<sigs.length; i++) {
				if (sigs[i][0]<SIG_CUTOFF) {
					writer.print("\t"+factorNames[indices[i]]);
					numSig++;
				}
			}
			writer.println();
			writer.print("Not significant: ");
			for (int i = 1; i<sigs.length; i++) {
				if (sigs[i][0]>=SIG_CUTOFF) {
					writer.print("\t"+factorNames[indices[i]]);
				}
			}
			writer.println();
			writer.print("Options not flagged: ");
			for (int i = 0; i<OPTIONS.length; i++) {
				if (!optionFlagged(OPTIONS[i])) {
					writer.print("\t"+OPTIONS[i]);
				}
			}
			writer.println();
			writer.println();
			writer.println();

			double min = Math.min(Array.max(Matrix.extractColumn(means, 0)), Array.max(Matrix.extractColumn(means, 1)));
			double max = Math.max(Array.max(Matrix.extractColumn(means, 0)), Array.max(Matrix.extractColumn(means, 1)));
			
			int maxFigs = Math.max(10, (int)Math.floor(Math.log10(Math.max(Math.max(Math.abs(min), Math.abs(max)), 1))));
			
			line = new String[] {"R-Sqr", "  Sig  ", "N", (logistic?ext.formStr("MeanAff", maxFigs, true)+ext.formStr("MeanUnaff", maxFigs, true):"")+"dir", (logistic?"Odds ratio (95% CI)   ":"Beta (95% CI)         "), "Factor", "Mean +/- SD", "Failures", "p-value"};
			writer.println(Array.toStr(order==null?line:Sort.putInOrder(line, order)));
			
			for (int i = 1; i<=M; i++) {
//				count = 0;
//				while (sigs[i][0]*Math.pow(10, count)<1&&count<20) {
//					count++;
//				}

				min = Math.min(means[i][0], means[i][1]);
				max = Math.max(means[i][0], means[i][1]);
				
				int meanFigs = 5-(int)Math.floor(Math.log10(Math.max(Math.max(Math.abs(min), Math.abs(max)), 1)));
				
				line = new String[9];
				line[0] = ext.formStr(ext.formDeci(sigs[i][1]*100, SIG_FIGS_PERCENTAGES*2, true)+"%", 5);
				line[1] = ext.formDeci(sigs[i][0], 5, true);
				line[2] = factorNs[i]+"";
				line[3] = (logistic?ext.formStr(ext.formDeci(means[i][1], meanFigs, true), maxFigs, true)+ext.formStr(ext.formDeci(means[i][0], meanFigs, true), maxFigs, true):"")+(factorDirections[i]?"+":"-");
				line[4] = ext.formDeci(effectsAndConfidenceIntervals[i][0],3,true)+" ("+ext.formDeci(effectsAndConfidenceIntervals[i][1],3,true)+", "+ext.formDeci(effectsAndConfidenceIntervals[i][2],3,true)+")";
				line[5] = factorNames[indices[i]];
				line[6] = "("+ext.formDeci(Array.mean(filterArray(data, i, Double.MIN_VALUE)), 3)+" +/- "+ext.formDeci(Array.stdev(filterArray(data, i, Double.MIN_VALUE)), 3)+")";
				line[7] = (failures[i]>0?", "+failures[i]+" failures (potential cause of bias)":"");
//				line[8] = "("+ext.formSciNot(sigs[i][0], 1, true)+")";
				line[8] = factorNs[i]==0?".":(logistic?"=NORMSDIST("+Math.sqrt(sigs[i][2])+")":"=TDIST("+Math.abs(sigs[i][2])+","+factorNs[i]+",2)");
				
				writer.println(Array.toStr(order==null?line:Sort.putInOrder(line, order)));
			}
			writer.println();
			writer.println();

			if (optionFlagged("force")) {
				deps = new Vector<String>();
				indeps = new Vector<double[]>();
				idV = new Vector<String[]>();
				for (int i = 0; i<N; i++) {
					dataline = new double[data[i].length-1];
					count = 0;
					for (int j = 1; j<data[i].length; j++) {
						dataline[j-1] = data[i][j];
						if (Math.abs(dataline[j-1]+999)<0.0001) {
							count++;
						}
					}
					if (count==0) {
						deps.add(data[i][0]+"");
						indeps.add(dataline);
						idV.add(ids[i]);
					}
				}
				model = logistic?(RegressionModel)new LogisticRegression(deps, indeps, false, optionFlagged("verbose")):(RegressionModel)new LeastSquares(deps, indeps, false, optionFlagged("verbose"));
				model.setSigFigs(sigfigs);
				line = new String[indices.length-1];
				for (int i = 0; i<M; i++) {
					line[i] = factorNames[indices[i+1]];
				}
				model.setVarNames(line);
				if (optionFlagged("oneperfamily")) {
					Date date = new Date();
					System.out.print("Forcing all variables into the model... ");
					model.onePerFamily(Matrix.extractColumn(Matrix.toStringArrays(idV), 1), famReps, bootReps);
					System.out.println(ext.getTimeElapsed(date.getTime()));
				}
				if (optionFlagged("dump")) {
					model.dumpData("force_data.xln");
				}
				writer.println(model.getSummary());
				if (model.getNumFailures()>0) {
					writer.println("There were "+model.getNumFailures()+" failures (potential cause of bias)");
				}
				writer.println();
				writer.println();
				predicteds = model.getPredicteds();
				residuals[0] = model.getResiduals();
			}

			if (optionFlagged("sw")) {
				writer.println("Stepwise regression using those significant by themselves:");
				sigSpots = new int[numSig];
				sigNames = new String[numSig];
				count = 0;
				for (int i = 1; i<sigs.length; i++) {
					if (sigs[i][0]<SIG_CUTOFF) {
						sigNames[count] = factorNames[indices[i]];
						sigSpots[count++] = i;
					}
				}
				if (count==0) {
					System.err.println("Stepwise regression canceled - no significant variables to use");
					writer.println("Stepwise regression canceled - no significant variables to use");
				} else {
					if (count==M&&optionFlagged("allsw")) {
						System.err.println(" - stepwise regression for all variables canceled since all are significant");
						flagOption("allsw", false);
					}
					deps = new Vector<String>();
					indeps = new Vector<double[]>();
					idV = new Vector<String[]>();
					for (int i = 0; i<N; i++) {
						dataline = new double[numSig];
						count = 0;
						for (int j = 0; j<sigSpots.length; j++) {
							dataline[j] = data[i][sigSpots[j]];
							if (Math.abs(dataline[j]+999)<0.0001) {
								count++;
							}
						}
						if (count==0) {
							deps.add(data[i][0]+"");
							indeps.add(dataline);
							idV.add(ids[i]);
						}
					}
					writer.println("N = "+deps.size());
					sw = new Stepwise(deps, indeps);
					sw.setVarNames(sigNames);
					if (optionFlagged("dump")) {
						sw.dumpData("sw_with_sig_data.xln");
					}
					writer.println(sw.getSummary());
					writer.println();
					writer.println("Final Model contains:");
					writer.println(sw.getFinalNames());
					predicteds = sw.getFinalPredicteds();
					residuals[0] = sw.getFinalResiduals();
					writer.println(sw.getAccuracySummary());
				}
			}

			if (optionFlagged("allsw")) {
				writer.println();
				writer.println("If all variables were considered in the stepwise regression, the final model would contain:");
				deps = new Vector<String>();
				indeps = new Vector<double[]>();
				idV = new Vector<String[]>();
				sigNames = new String[M];
				for (int i = 1; i<=M; i++) {
					sigNames[i-1] = factorNames[indices[i]];
				}
				for (int i = 0; i<N; i++) {
					dataline = new double[M];
					count = 0;
					for (int j = 0; j<M; j++) {
						dataline[j] = data[i][j+1];
						if (Math.abs(dataline[j]+999)<0.0001) {
							count++;
						}
					}
					if (count==0) {
						deps.add(data[i][0]+"");
						indeps.add(dataline);
						idV.add(ids[i]);
					}
				}
				writer.println("N = "+deps.size());
				sw = new Stepwise(deps, indeps);
				sw.setVarNames(sigNames);
				if (optionFlagged("dump")) {
					sw.dumpData("sw_with_all_data.xln");
				}
				writer.println(sw.getFinalNames());
				// writer.println();
				// writer.println(sw.getSummary());
				writer.println(sw.getAccuracySummary());
			}
			writer.close();

			if ((optionFlagged("sw") || optionFlagged("force"))&&optionFlagged("predicteds")) {
				writer = new PrintWriter(new FileWriter(ext.rootOf(filename)+"_predicteds.xln"));
				for (int i = 0; i<predicteds.length; i++) {
					writer.println(data[i][0]+"\t"+predicteds[i]);
				}
				writer.close();
			}

			if ((optionFlagged("sw") || optionFlagged("force"))&&(optionFlagged("residuals")||optionFlagged("normalized")||optionFlagged("inverseNormalized"))) {
				if (optionFlagged("normalized")) {
					residuals[1] = Transformations.transform(residuals[0], Transformations.NORMALIZE);
				}
				if (optionFlagged("inverseNormalized")) {
					residuals[2] = Transformations.transform(residuals[0], Transformations.INVERSE_NORMALIZE);
				}
				writer = new PrintWriter(new FileWriter(ext.rootOf(filename)+"_residuals.xln"));
				temp = "";
				for (int i = 0; i<DEFAULT_ID_NAMES.length; i++) {
					if (idIndices[i] != -1) {
						temp += (temp.equals("")?"":"\t")+DEFAULT_ID_NAMES[i][0];
					}
                }
				if (optionFlagged("residuals")) {
					temp += (temp.equals("")?"":"\t")+"Residuals";
				}
				if (optionFlagged("normalized")) {
					temp += (temp.equals("")?"":"\t")+"NormalizedResiduals";
				}
				if (optionFlagged("inverseNormalized")) {
					temp += (temp.equals("")?"":"\t")+"InverseNormalizedResiduals";
				}
				writer.println(temp);
				for (int i = 0; i<residuals[0].length; i++) {
					temp = "";
					for (int j = 0; j<DEFAULT_ID_NAMES.length; j++) {
						line = idV.elementAt(i);
						if (idIndices[j] != -1) {
							temp += (temp.equals("")?"":"\t")+line[j];
						}
	                }
					if (optionFlagged("residuals")) {
						temp += (temp.equals("")?"":"\t")+residuals[0][i];
					}
					if (optionFlagged("normalized")) {
						temp += (temp.equals("")?"":"\t")+residuals[1][i];
					}
					if (optionFlagged("inverseNormalized")) {
						temp += (temp.equals("")?"":"\t")+residuals[2][i];
					}
					writer.println(temp);
				}
				writer.close();
			}
			System.out.println();
		}
		new File(filename+".bak").delete();
	}

	public static void doChis(String[][] ids, double[][] data, String[] factorNames, int[] indices, boolean audit) {
		PrintWriter writer;
		String[] fams;
		int count;
		double min;
		Hashtable<String,Vector<String>> hashVec;
		Hashtable<String,Vector<double[]>> hashVecData;
		Vector<String> members;
		Vector<double[]> vDoubleArray;
		double[] dataline, subtotals;
		double[][] avgCounts, baseCounts;
		DoubleVector dv, dv1, dv2;

		fams = Array.unique(Matrix.extractColumn(ids, 1));

		hashVec = new Hashtable<String,Vector<String>>();
		hashVecData = new Hashtable<String,Vector<double[]>>();
		dv1 = new DoubleVector();
		dv2 = new DoubleVector();
		for (int i = 0; i<ids.length; i++) {
			if (hashVecData.containsKey(ids[i][1])) {
				vDoubleArray = hashVecData.get(ids[i][1]);
				members = hashVec.get(ids[i][1]+"mems");
			} else {
				hashVecData.put(ids[i][1], vDoubleArray = new Vector<double[]>());
				hashVec.put(ids[i][1]+"mems", members = new Vector<String>());
			}
			if (members.indexOf(ids[i][2])==-1) {
				vDoubleArray.add(new double[] {data[i][0], data[i][1], data[i][2]});
				dv1.addIfAbsent(data[i][1]);
				dv2.addIfAbsent(data[i][2]);
				members.add(ids[i][2]);
			} else {
				dataline = vDoubleArray.elementAt(members.indexOf(ids[i][2]));
				if (data[i][0]!=dataline[0]||data[i][1]!=dataline[1]||data[i][2]!=dataline[2]) {
					System.err.println("Error - "+ids[i][1]+"-"+ids[i][2]+" has been genotyped twice; first as "+dataline[1]+"/"+dataline[2]+" and second as "+data[1]+"/"+data[2]+"; ignoring second");
				}
			}
		}
		dv = new DoubleVector();
		for (int i = 0; i<dv1.size(); i++) {
			dv.addIfAbsent(dv1.elementAt(i));
		}
		for (int i = 0; i<dv2.size(); i++) {
			dv.addIfAbsent(dv2.elementAt(i));
		}
		if (dv.size()==dv1.size()+dv2.size()) {
			throw new RuntimeException("Error - no overlap in alleles between "+factorNames[indices[1]]+" and "+factorNames[indices[2]]);
		}
		if (dv.size()>MAX_CLASSES) {
			throw new RuntimeException("Error - more than "+MAX_CLASSES+" alleles (n="+dv.size()+") between "+factorNames[indices[1]]+" and "+factorNames[indices[2]]+"; aborting analyses.");
		}

		if (audit) {
			try {
				writer = new PrintWriter(new FileWriter("audit.out"));
				for (int i = 0; i<fams.length; i++) {
					writer.print(fams[i]);
					vDoubleArray = hashVecData.get(fams[i]);
					for (int j = 0; j<vDoubleArray.size(); j++) {
						dataline = vDoubleArray.elementAt(j);
						writer.print("\t"+dataline[1]+"\t"+dataline[2]);
					}
					writer.println();
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing audit file for chis");

			}
		}

		min = Array.min(Matrix.extractColumn(data, 0));

		avgCounts = new double[2][dv.size()];
		for (int i = 0; i<fams.length; i++) {
			vDoubleArray = hashVecData.get(fams[i]);
			baseCounts = new double[2][dv.size()];
			for (int j = 0; j<vDoubleArray.size(); j++) {
				dataline = vDoubleArray.elementAt(j);
				baseCounts[dataline[0]==min?0:1][dv.indexOf(dataline[1])]++;
				baseCounts[dataline[0]==min?0:1][dv.indexOf(dataline[2])]++;
			}
			count = Array.sum(baseCounts[1])>0?1:0;
			for (int j = 0; j<dv.size(); j++) {
				avgCounts[count][j] += (double)baseCounts[count][j]/(double)(Array.sum(baseCounts[count])/2);
			}
		}
		subtotals = new double[] {Array.sum(avgCounts[0]), Array.sum(avgCounts[1])};

		try {
			writer = new PrintWriter(new FileWriter("chis.out"));
			for (int i = 0; i<dv.size(); i++) {
				writer.print("      \t"+(int)dv.elementAt(i)+"      ");
			}
			writer.println("      \tOverall significance");
			for (int i = avgCounts.length-1; i>=0; i--) {
				writer.print(i==1?"Affected":"Unaffected");
				for (int j = 0; j<avgCounts[0].length; j++) {
					writer.print("\t"+ext.formDeci(avgCounts[i][j], 3, true)+" ("+ext.formDeci(avgCounts[i][j]/subtotals[i], 2, true)+")");
				}
				if (i==0) {
					writer.print("\tp="+ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(avgCounts), avgCounts[0].length-1)+""));
				}
				writer.println();
			}

			// writer.print("Pairwise");
			// for (int i = 0; i<avgCounts[0].length; i++) {
			// writer.print("\tp="+ProbDist.ChiDist(ContingencyTable.ChiSquare(new
			// double[][] {{avgCounts[0][i], avgCounts[1][i]}, {subtotals[0],
			// subtotals[1]}}), 1));
			//
			// }
			// writer.println();
			writer.print("Pairwise");
			for (int i = 0; i<avgCounts[0].length; i++) {
				writer.print("\tp="+ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(new double[][] { {avgCounts[0][i], avgCounts[1][i]}, {subtotals[0], subtotals[1]}}), 1)+"")+"      ");

			}
			writer.println();
			writer.println();
			writer.println();
			writer.println("Allele\tAffected\t% of Affecteds\tUnaffected\t% of Unaffecteds\tOR\tp-value");
			for (int i = 0; i<dv.size(); i++) {
				writer.print((int)dv.elementAt(i));
				writer.print("\t"+ext.formDeci(avgCounts[1][i], 3, true)+"\t"+ext.formDeci(avgCounts[1][i]/subtotals[1]*100, 1, true));
				writer.print("\t"+ext.formDeci(avgCounts[0][i], 3, true)+"\t"+ext.formDeci(avgCounts[0][i]/subtotals[0]*100, 1, true));
				writer.print("\t"+ext.formDeci((avgCounts[1][i]/subtotals[1]*(1-avgCounts[0][i]/subtotals[0]))/(avgCounts[0][i]/subtotals[0]*(1-avgCounts[1][i]/subtotals[1])), 2, true));
				writer.print("\t"+ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(new double[][] { {avgCounts[0][i], avgCounts[1][i]}, {subtotals[0], subtotals[1]}}), 1)+""));
				writer.println();
			}
			writer.println("\tOverall significance");
			writer.print("\tp="+ext.prettyP(ProbDist.ChiDist(ContingencyTable.ChiSquare(avgCounts), avgCounts[0].length-1)+""));

			writer.close();

		} catch (Exception e) {
			System.err.println("Error writing chis.out");

		}
	}

	public static void doHWE(String[][] ids, double[][] data, String[] factorNames, int[] indices, boolean audit) {
		PrintWriter writer;
		String[] fams;
		int count;
		int min;
		Hashtable<String,Vector<String>> hashVec;
		Hashtable<String,Vector<int[]>> hashVecData;
		Vector<String> members;
		Vector<int[]> genotypes;
		int[] dataline;
		double[] subtotals;
		double[][] avgCounts, baseCounts;
		IntVector iv, iv1, iv2;
		Vector<String> types = new Vector<String>();
		int[] order, orderedAlleles;
		double[] countsInOrder;
		String[] alleleLabels;

		fams = Array.unique(Matrix.extractColumn(ids, 1));

		hashVec = new Hashtable<String,Vector<String>>();
		hashVecData = new Hashtable<String,Vector<int[]>>();
		iv1 = new IntVector();
		iv2 = new IntVector();
		for (int i = 0; i<ids.length; i++) {
			if (hashVecData.containsKey(ids[i][1])) {
				genotypes = hashVecData.get(ids[i][1]);
				members = hashVec.get(ids[i][1]+"mems");
			} else {
				hashVecData.put(ids[i][1], genotypes = new Vector<int[]>());
				hashVec.put(ids[i][1]+"mems", members = new Vector<String>());
			}
			if (members.indexOf(ids[i][2])==-1) {
				genotypes.add(new int[] {(int)data[i][0], (int)data[i][1], (int)data[i][2]});
				iv1.addIfAbsent((int)data[i][1]);
				iv2.addIfAbsent((int)data[i][2]);
				HashVec.addIfAbsent((int)data[i][1]+"/"+(int)data[i][2], types);
				members.add(ids[i][2]);
			} else {
				dataline = genotypes.elementAt(members.indexOf(ids[i][2]));
				if (data[i][0]!=dataline[0]||data[i][1]!=dataline[1]||data[i][2]!=dataline[2]) {
					System.err.println("Error - "+ids[i][1]+"-"+ids[i][2]+" has been genotyped twice; first as "+dataline[1]+"/"+dataline[2]+" and second as "+data[1]+"/"+data[2]+"; ignoring second");
				}
			}
		}
		iv = new IntVector();
		for (int i = 0; i<iv1.size(); i++) {
			iv.addIfAbsent(iv1.elementAt(i));
		}
		for (int i = 0; i<iv2.size(); i++) {
			iv.addIfAbsent(iv2.elementAt(i));
		}
		if (iv.size()==iv1.size()+iv2.size()) {
			throw new RuntimeException("Error - no overlap in alleles between "+factorNames[indices[1]]+" and "+factorNames[indices[2]]);
		}
		if (iv.size()>MAX_CLASSES) {
			throw new RuntimeException("Error - more than "+MAX_CLASSES+" alleles (n="+iv.size()+") between "+factorNames[indices[1]]+" and "+factorNames[indices[2]]+"; aborting analyses.");
		}

		if (audit) {
			try {
				writer = new PrintWriter(new FileWriter("audit.out"));
				for (int i = 0; i<fams.length; i++) {
					writer.print(fams[i]);
					genotypes = hashVecData.get(fams[i]);
					for (int j = 0; j<genotypes.size(); j++) {
						dataline = genotypes.elementAt(j);
						writer.print("\t"+dataline[1]+"\t"+dataline[2]);
					}
					writer.println();
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing audit file for hwe");
			}
		}

		min = (int)Array.min(Matrix.extractColumn(data, 0));

		avgCounts = new double[2][types.size()];
		for (int i = 0; i<fams.length; i++) {
			genotypes = hashVecData.get(fams[i]);
			baseCounts = new double[2][types.size()];
			for (int j = 0; j<genotypes.size(); j++) {
				dataline = genotypes.elementAt(j);
				baseCounts[dataline[0]==min?0:1][types.indexOf(dataline[1]+"/"+dataline[2])]++;
			}
			count = Array.sum(baseCounts[1])>0?1:0;
			for (int j = 0; j<types.size(); j++) {
				avgCounts[count][j] += (double)baseCounts[count][j]/(double)Array.sum(baseCounts[count]);
			}
		}
		subtotals = new double[] {Array.sum(avgCounts[0]), Array.sum(avgCounts[1])};

		try {
			writer = new PrintWriter(new FileWriter(factorNames[indices[1]]+"_"+factorNames[indices[2]]+"_hwe.out"));
			orderedAlleles = Sort.putInOrder(iv.toArray());
			if (iv.size()==3) {
				alleleLabels = new String[] {orderedAlleles[0]+"/"+orderedAlleles[0], orderedAlleles[0]+"/"+orderedAlleles[1], orderedAlleles[0]+"/"+orderedAlleles[2], orderedAlleles[1]+"/"+orderedAlleles[1], orderedAlleles[1]+"/"+orderedAlleles[2], orderedAlleles[2]+"/"+orderedAlleles[2]};
			} else if (iv.size()<3) {
				alleleLabels = new String[] {orderedAlleles[0]+"/"+orderedAlleles[0], orderedAlleles[0]+"/"+orderedAlleles[1], orderedAlleles[1]+"/"+orderedAlleles[1]};
			} else {
				System.err.println("Error - HWE is not currently set up to handle more than 3 alleles");
				System.err.println("   (found "+Array.toStr(iv.toArray())+")");
				System.exit(1);
				alleleLabels = null;
			}
			order = new int[alleleLabels.length];
			for (int i = 0; i<alleleLabels.length; i++) {
				order[i] = types.indexOf(alleleLabels[i]);
			}

			for (int i = 0; i<alleleLabels.length; i++) {
				writer.print("      \t"+alleleLabels[i]);
			}
			writer.println();
			for (int i = avgCounts.length-1; i>=0; i--) {
				writer.print(i==1?"Affected":"Unaffected");
				for (int j = 0; j<alleleLabels.length; j++) {
					if (order[j]==-1) {
						writer.print("\t0.000 (0.00%)");
					} else {
						writer.print("\t"+ext.formDeci(avgCounts[i][order[j]], 3, true)+" ("+ext.formDeci(avgCounts[i][order[j]]/subtotals[i], 2, true)+")");
					}
				}
				if (iv.size()<=3) {
					countsInOrder = new double[order.length];
					for (int j = 0; j<order.length; j++) {
						countsInOrder[j] = order[j]==-1?0:avgCounts[i][order[j]];
					}
					writer.print("\tp="+ext.prettyP(AlleleFreq.HWEsig(countsInOrder)));
				}
				writer.println();
			}
			writer.println();
			if (iv.size()>3) {
				writer.println("Since there were more than three alleles, the genotypes were not checked for Hardy-Weinberg equilibrium.");
			} else {
				countsInOrder = new double[order.length];
				for (int j = 0; j<order.length; j++) {
					countsInOrder[j] = order[j]==-1?0:avgCounts[0][order[j]]+avgCounts[1][order[j]];
				}
				writer.println("Hardy-Weinberg for combined sample is: p="+ext.prettyP(AlleleFreq.HWEsig(countsInOrder)));
			}
			writer.println();
			writer.println();

			for (int i = avgCounts.length-1; i>=0; i--) {
				writer.print("\t"+(i==1?"Affected":"Unaffected"));
			}
			writer.println();
			for (int j = 0; j<alleleLabels.length; j++) {
				writer.print(alleleLabels[j]);
				for (int i = avgCounts.length-1; i>=0; i--) {
					if (order[j]==-1) {
						writer.print("\t0.000\t(0.00%)");
					} else {
						writer.print("\t"+ext.formDeci(avgCounts[i][order[j]], 3, true)+"\t("+ext.formDeci(avgCounts[i][order[j]]/subtotals[i]*100, SIG_FIGS_PERCENTAGES, true)+"%)");
					}
				}
				writer.println();
			}
			writer.println();

			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing hwe.out");

		}

	}

	public static double[] filterArray(double[][] data, int column, double missing) {
		DoubleVector dv = new DoubleVector();
		for (int i = 0; i<data.length; i++) {
			if (data[i][column]!=missing) {
				dv.add(data[i][column]);
			}
		}

		return dv.toArray();
	}

	public boolean containsStr(String target, String[] list) {
		for (int i = 0; i<list.length; i++) {
			if (list[i].equals(target)) {
				return true;
			}
		}
		return false;
	}

	public static String[][] separateFactorNames(String[] factorNames, int[] indices) {
		String[][] varNames = new String[2][];

		varNames[0] = new String[] {factorNames[indices[0]]};
		varNames[1] = new String[indices.length-1];

		for (int i = 1; i<indices.length; i++) {
			varNames[1][i-1] = factorNames[indices[i]];
		}

		return varNames;
	}

	public static int[][] dummyIntMatrix(double[] array) {
		int[][] matrix = new int[array.length][1];

		for (int i = 0; i<array.length; i++) {
			if (Math.abs(array[i]-(int)array[i])>0.00001) {
				System.err.println("Error - there may be less than "+MAX_CLASSES+" classes, but their not integers (found a '"+array[i]+"')");
			}
			matrix[i][0] = (int)array[i];
		}

		return matrix;
	}

	public boolean optionFlagged(String str) {
		for (int i = 0; i<OPTIONS.length; i++) {
			if (OPTIONS[i].toLowerCase().equals(str.toLowerCase())) {
				return flags[i];
			}
		}
		System.err.println("Error - tried to access an option that doesn't exist: '"+str+"'");
		return false;
	}

	public void flagOption(String str) {
		flagOption(str, true);
	}

	public void flagOption(String str, boolean state) {
		for (int i = 0; i<OPTIONS.length; i++) {
			if (OPTIONS[i].toLowerCase().equals(str.toLowerCase())) {
				flags[i] = state;
				return;
			}
		}
		System.err.println("Error - tried to flag an option that doesn't exist: '"+str+"'");
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		boolean suppress = false;
		// String filename = DEFAULT_INPUT;
		// String filename = "Dementia.MMSE-depressed.ctl";
		// String filename = "Depressed_Genetics_OnePer.ctl";
		// String filename = "Demented.ctl";
		// String filename = "Copy of APOE.ctl";
		// String filename = "Rep1.ctl";
		// String filename = "Affected.ctl";
		// String filename = "GBA_Table.ctl";
		// String filename = "all_VPD_chis.ctl";
		// String filename = "Affected_Hap_OnePer_HWE.ctl";
//		String filename = "ARIC_AboveBelow.ctl";
//		String filename = "residuals.ctl";
		String filename = "analyzeSignals.ctl";		

		String usage = "\n"+"park.comp requires 1 argument\n"+"   (1) filename (i.e. 'file="+DEFAULT_INPUT+"')\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("-suppress")) {
				suppress = true;
				numArgs--;
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
			if (new File(filename).exists()&&new File(filename).length()==0) {
				System.err.println("Using all available info to predict "+DEFAULT_TRAIT+" (also creating '"+DEFAULT_INPUT+"' if you want to alter)");
				DEFAULT_INPUT = filename;
				new comp();
			} else if (args.length==0&&(!new File(filename).exists()||new File(filename).length()==0)) {
				System.err.println("Using all available info to predict "+DEFAULT_TRAIT+" (also creating '"+DEFAULT_INPUT+"' if you want to alter)");
				new comp();
			} else {
				new comp(filename);
			}
		} catch (IOException e) {
			System.err.println("IO Exception");
		} catch (RuntimeException e) {
			System.err.println("Runtime Error: "+e.getMessage());
			e.printStackTrace();
			if (!new File(filename).exists()&&new File(filename+".bak").exists()) {
				new File(filename+".bak").renameTo(new File(filename));
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (System.getProperty("os.name").startsWith("Windows") && !suppress) {
				System.out.println("...done");
				ext.waitForResponse();
			}
		}
	}
}
