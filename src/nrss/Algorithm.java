// -Xms1024M -Xmx1024M
package nrss;

import java.io.*;
import java.util.*;

import common.*;
import filesys.LDdatabase;
import filesys.LongLDdb;
import filesys.SnpMarkerSet;
import bioinformatics.*;

public class Algorithm {
	public static final double INDEX_THRESHOLD = 0.0001;
	public static final double INCLUSION_THRESHOLD = 0.01;
	public static final int DEFAULT_WINDOW = 50000;

	// public static final String[][] MODELS = {{"0.001", "0.01", "150000"},
	// {"0.001", "0.01", "50000"},
	// {"0.001", "0.01", "25000"},
	// {"0.0001", "0.01", "250000"},
	// {"0.0001", "0.01", "150000"},
	// {"0.0001", "0.01", "25000"},
	// };

	// public static final String[][] MODELS = {{"0.01", "0.1", "150000"}};

	public static final String[][] MODELS = {{"0.001", "0.01", "150000"}};

//	public static final int EST_NUM_MARKERS_IN_LARGEST_CHR = 25000; // GWAS
	public static final int EST_NUM_MARKERS_IN_LARGEST_CHR = 320000; // imputation

//	public static final String DEFAULT_LD_LOC = "";

	// public static final String DEFAULT_LD_LOC = "C:\\Documents and
	// Settings\\npankrat\\My Documents\\gwas\\udall\\LD\\";

	public static void procFile(String dir, String filename, int p_column, int variate_column, String ldRoot, int window, String outfile, int genomeBuild) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Hashtable<String,String> chrHash;
		Vector<String> markerVector = new Vector<String>(EST_NUM_MARKERS_IN_LARGEST_CHR);
		IntVector markerLocations = new IntVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
		DoubleVector pvalueVector = new DoubleVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
		DoubleVector varVector = variate_column==-1?null:new DoubleVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
		int count, lowIndex, highIndex, countMissing;
		IntVector indexVector = new IntVector();
		Vector<String> inclusionVector = new Vector<String>();
		String chr;
		double d;
		String[] markerNames;
		int[] positions;
		int[][] markerPositions;
		double[] pvals, vars = null;
		double[][] stats;
		int[] indexSNPs, markerCounts;
		boolean done;
		IntVector cluster;
		double maxStat;
		int maxSNP;
		double minPvalue;
		Vector<String> v = new Vector<String>();
		String[] genes;
		LDdatabase lddb;
		
		lddb = new LDdatabase(ldRoot, LDdatabase.TYPE_LONG, new Logger());
		chrHash = lddb.getChrHash();

		try {
			System.out.println(ext.getTime());
			reader = new BufferedReader(new FileReader(dir+filename));
			line = reader.readLine().trim().split("[\\s]+");
			ext.checkHeader(line, new String[] {"Chr", "Position"}, new int[] {1,2}, false, new Logger(), true);
			line = reader.readLine().trim().split("[\\s]+");
			chr = "";
			count = 0;
			countMissing = 0;
			done = false;
            new File(dir+"NOT_PRESENT_IN_LD_SOURCE.txt").delete();
			while (!done) {
				if (chr.equals("")) {
					chr = line[1];
					count = 0;
				} else {
					if (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
					} else {
						done = true;
						line[1] = "done";
					}
				}

				if (line[1].equals(chr)) {
					if (!chrHash.containsKey(line[0])) {
						countMissing++;
						try {
	                        writer = new PrintWriter(new FileWriter(dir+"NOT_PRESENT_IN_LD_SOURCE.txt", true));
	                        writer.println(line[0]);
	                        writer.close();
                        } catch (Exception e) {
	                        System.err.println("Error writing to "+dir+"NOT_PRESENT_IN_LD_SOURCE.txt");
	                        e.printStackTrace();
                        }
					} else if (!line[p_column].equals(".") && !line[p_column].equals("NA")) {
						d = Double.parseDouble(line[p_column]);
						pvalueVector.add(d);
						if (varVector!=null) {
							varVector.add(Double.parseDouble(line[variate_column]));
						}
						if (d<INDEX_THRESHOLD) {
							indexVector.add(count);
						}
						if (d < INCLUSION_THRESHOLD) { // added to optimize lddb efficiency
							inclusionVector.add(line[0]);
						}
						markerVector.add(line[0]);
						markerLocations.add(Integer.parseInt(line[2]));
						count++;
					}
				} else {
					markerNames = Array.toStringArray(markerVector);
					positions = markerLocations.toArray();
					pvals = pvalueVector.toArray();
					indexSNPs = indexVector.toArray();
					markerLocations.clear();
					indexVector.clear();
					markerVector.clear();
					pvalueVector.clear();
					if (varVector!=null) {
						vars = varVector.toArray();
						varVector.clear();
					}
					
			        lddb.updateWithTheseMarkers(Array.toStringArray(inclusionVector), ext.replaceDirectoryCharsWithUnderscore(dir+filename, 2)+"."+chr);
					inclusionVector.clear();
					if (countMissing>0) {
						System.err.println("Error - number of markers on chr"+chr+" added to not present file: "+countMissing);
						countMissing = 0;
					}

					stats = new double[indexSNPs.length][];
					markerCounts = new int[indexSNPs.length];
					System.out.println("Analyzing "+indexSNPs.length+" index SNPs on chromosome "+chr);
//					if (indexSNPs.length>0) {
//						// hash = new Hashtable<String, String>();
//						hash = loadLD(DEFAULT_LD_LOC+"chr"+chr+".recode.ped.LD");
//					} else {
//						hash = null;
//					}
					
					cluster = new IntVector();

					for (int i = 0; i<indexSNPs.length; i++) {
						lowIndex = indexSNPs[i];
						minPvalue = 1;
						while (lowIndex>=0&&positions[lowIndex]>positions[indexSNPs[i]]-window) {
							if (pvals[lowIndex]<minPvalue) {
								minPvalue = pvals[lowIndex];
							}
							lowIndex--;
						}
						lowIndex++;
						highIndex = indexSNPs[i];
						while (highIndex<positions.length&&positions[highIndex]<positions[indexSNPs[i]]+window) {
							if (pvals[highIndex]<minPvalue) {
								minPvalue = pvals[highIndex];
							}
							highIndex++;
						}
						highIndex--;
						stats[i] = computeStatistic(markerNames, pvals, lowIndex, highIndex, lddb.getLongChrLDdb(Positions.chromosomeNumber(chr)), INCLUSION_THRESHOLD, vars);
						markerCounts[i] = highIndex-lowIndex+1;
						// System.out.println((i+1)+")
						// "+markerNames[indexSNPs[i]]+"\t"+positions[indexSNPs[i]]+"\t"+ext.formDeci(stats[i],
						// 2)+" (using "+(highIndex-lowIndex+1)+" markers)");

						if (i==0) {
							cluster.add(i);
						} else if (positions[indexSNPs[i]]-positions[indexSNPs[i-1]]<window) {
							cluster.add(i);
						} else {
							maxStat = 0;
							maxSNP = -1;
							for (int j = 0; j<cluster.size(); j++) {
								if (stats[cluster.elementAt(j)][0]>maxStat) {
									maxSNP = cluster.elementAt(j);
									maxStat = stats[maxSNP][0];
								}
							}
							v.add(markerNames[indexSNPs[maxSNP]]+"\t"+chr+"\t"+positions[indexSNPs[maxSNP]]+"\t"+ext.formDeci(stats[maxSNP][0], 2)+"\t"+markerCounts[maxSNP]+"\t"+ext.formDeci(stats[maxSNP][0]*Math.sqrt(markerCounts[maxSNP]), 2)+"\tchr"+chr+":"+(positions[indexSNPs[maxSNP]]-window)+"-"+(positions[indexSNPs[maxSNP]]+window)+"\t"+minPvalue+(variate_column==-1?"":"\t"+Array.toStr(Array.subArray(stats[maxSNP], 1, stats[maxSNP].length))));
							cluster.clear();
							cluster.add(i);
						}

						if (i==indexSNPs.length-1) {
							maxStat = 0;
							maxSNP = -1;
							for (int j = 0; j<cluster.size(); j++) {
								if (stats[cluster.elementAt(j)][0]>maxStat) {
									maxSNP = cluster.elementAt(j);
									maxStat = stats[maxSNP][0];
								}
							}
							v.add(markerNames[indexSNPs[maxSNP]]+"\t"+chr+"\t"+positions[indexSNPs[maxSNP]]+"\t"+ext.formDeci(stats[maxSNP][0], 2)+"\t"+markerCounts[maxSNP]+"\t"+ext.formDeci(stats[maxSNP][0]*Math.sqrt(markerCounts[maxSNP]), 2)+"\tchr"+chr+":"+(positions[indexSNPs[maxSNP]]-window)+"-"+(positions[indexSNPs[maxSNP]]+window)+"\t"+minPvalue+(variate_column==-1?"":"\t"+Array.toStr(Array.subArray(stats[maxSNP], 1, stats[maxSNP].length))));
						}
					}

					markerNames = null;
					positions = null;
					pvals = null;
					indexSNPs = null;
					markerVector = new Vector<String>(EST_NUM_MARKERS_IN_LARGEST_CHR);
					markerLocations = new IntVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
					pvalueVector = new DoubleVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
					indexVector = new IntVector();
					chr = "";
				}
			}
			reader.close();

			markerPositions = new int[v.size()][2];
			for (int i = 0; i<v.size(); i++) {
				line = v.elementAt(i).trim().split("[\\s]+");
				markerPositions[i][0] = Integer.parseInt(line[1]);
				markerPositions[i][1] = Integer.parseInt(line[2]);
			}
			genes = MapSNPsAndGenes.mapSNPsToGenesLoosely(markerPositions, 0, genomeBuild, new Logger());
			
			writer = new PrintWriter(new FileWriter(dir+outfile));
			writer.println("Marker\tChr\tPostition\tGenes(s)\tWeightedStatistic\tNumMarkers\tUnweightedStatistic\tUCSC coordinates\tMin p-value"+(variate_column==-1?"":"\tvars"));
			for (int i = 0; i<v.size(); i++) {
				line = v.elementAt(i).trim().split("[\\s]+");
				writer.println(Array.toStr(Array.subArray(line, 0, 3))+"\t"+genes[i]+"\t"+Array.toStr(Array.subArray(line, 3)));
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
		
		countMissing = Files.countLines(dir+"NOT_PRESENT_IN_LD_SOURCE.txt", false);
		if (countMissing > 0) {
			System.out.println("\nFYI - there were "+countMissing+" marker(s) exclused from inclusion in NRS scores due to lack of LD information in source files");
		}
	}

	public static Hashtable<String,String> loadLD(String filename) {
		BufferedReader reader;
		String[] line;
		Hashtable<String,String> hash = new Hashtable<String,String>();

		try {
			reader = new BufferedReader(new FileReader(filename));
			reader.readLine();
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				hash.put(line[0]+":"+line[1], line[4]);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}

		return hash;
	}

	public static void procRegion(String results, int p_column, LongLDdb chrLDdb) {
		BufferedReader reader = null;
		String[] line;
//		Hashtable<String,String> hash;
		Vector<String> markerVector = new Vector<String>();
		DoubleVector pvalueVector = new DoubleVector();

//		hash = loadLD(ld);
		
		try {
			reader = new BufferedReader(new FileReader(results));
			reader.readLine();
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				markerVector.add(line[0]);
				pvalueVector.add(Double.parseDouble(line[p_column]));
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+results+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+results+"\"");
			System.exit(2);
		}

		System.out.println("Stat is "+ext.formDeci(computeStatistic(Array.toStringArray(markerVector), pvalueVector.toArray(), chrLDdb, INCLUSION_THRESHOLD)[0], 2));
	}

	public static double[] computeStatistic(String[] markerNames, double[] pvalues, LongLDdb chrLDdb, double threshold) {
		return computeStatistic(markerNames, pvalues, 0, pvalues.length-1, chrLDdb, threshold, null);
	}

	public static double[] computeStatistic(String[] markerNames, double[] pvalues, int lowIndex, int highIndex, LongLDdb chrLDdb, double threshold, double[] vars) {
		double[] stats = new double[5];
		int[] keys;
		double maxR2, r2;
		double[] pval_window, var_window = null;
		int count;

		if (vars!=null) {
			var_window = new double[highIndex-lowIndex+1];
			for (int i = 0; i<var_window.length; i++) {
				var_window[i] = vars[lowIndex+i];
			}
		}

		pval_window = new double[highIndex-lowIndex+1];
		for (int i = 0; i<pval_window.length; i++) {
			pval_window[i] = pvalues[lowIndex+i];
		}
		keys = Sort.quicksort(pval_window);

		count = 0;
		stats[2] = 1;
		for (int i = 0; i<pval_window.length; i++) {
			if (pval_window[keys[i]]<threshold) {
				maxR2 = 0;
				for (int j = 0; j<i; j++) {
					r2 = chrLDdb.get(markerNames[lowIndex+keys[i]], markerNames[lowIndex+keys[j]]);
            		if (r2 == LDdatabase.MISSING_INFO) {
//            			if (!allMissingMarkers.containsKey(subset[index]) || !allMissingMarkers.containsKey(subset[index+k*offset])) {
            			System.err.println("Error - missing LD info for "+markerNames[lowIndex+keys[i]]+"/"+markerNames[lowIndex+keys[j]]+" pair");
//            			}
            		} else if (r2>maxR2) {
						maxR2 = r2;
					}
				}
				stats[0] += (1-maxR2)*minusLog(pval_window[keys[i]]);
				if (vars!=null) {
					stats[1] += var_window[keys[i]];
					stats[2] = Math.min(var_window[keys[i]], stats[2]);
					if (var_window[keys[i]]<0.05) {
						stats[3]++;
					}
					if (var_window[keys[i]]<0.10) {
						stats[4]++;
					}
					count++;
				}
			}
		}

		stats[0] /= Math.sqrt(pval_window.length);
		if (vars!=null) {
			stats[1] /= (double)count;
			stats[3] /= (double)count;
			stats[4] /= (double)count;
		}

		return stats;
	}

	public static double minusLog(double p) {
		return -1*Math.log10(p);
	}

	public static void batchHaploview(int numBatches) {
		String commands = "";

		// commands += "plink --bfile pd_gwas --chr [%0] --out chr[%0] --recode
		// --remove WGAs.txt\n";
		commands += "plink --bfile plink --chr [%0] --out chr[%0].recode --recode\n\n";
		commands += "jcp nrss.Nrss procMap=[%0]\n";
		commands += "java -Dsun.java2d.noddraw=true -Xmx1024m -classpath /home/npankrat/Haploview.jar -Djava.awt.headless=true edu.mit.wi.haploview.HaploView -nogui -log -pedfile chr[%0].recode.ped -info chr[%0].map -skipcheck -hwcutoff 0 -maxDistance 500 -dprime\n\n";

		Files.batchIt("haplo", null, numBatches, commands, Array.stringArraySequence(23, ""));
	}

	public static void procMapForHaploview(String in, String out) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;

		try {
			reader = new BufferedReader(new FileReader(in));
			writer = new PrintWriter(new FileWriter(out));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				writer.println(line[1]+"\t"+line[3]);
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+in+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+in+"\"");
			System.exit(2);
		}
	}

	public static void simulateNull(String pedfile, int replicates) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Vector<String> v = new Vector<String>();
		String[][] data;
		int[] order;
		int count;

		try {
			reader = new BufferedReader(new FileReader(pedfile));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				v.add(line[0]+"\t"+line[1]+"\t"+line[5]);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error - could not find "+pedfile+" in current directory");
			System.exit(2);
		} catch (IOException ioe) {
			System.err.println("Error parsing "+pedfile+"");
			System.exit(3);
		}

		data = new String[v.size()][];
		for (int i = 0; i<v.size(); i++) {
			data[i] = v.elementAt(i).split("[\\s]+");
		}

		for (int i = 1; i<=replicates; i++) {
			try {
				writer = new PrintWriter(new FileWriter("sim."+i+".dat"));
				writer.println("FID\tIID\tAff");
				order = Array.random(data.length);
				for (int j = 0; j<data.length; j++) {
					writer.println(data[j][0]+"\t"+data[j][1]+"\t"+data[order[j]][2]);
				}
				writer.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		try {
			writer = null;
			count = 0;
			for (int i = 0; i<replicates; i++) {
				if (i%250==0) {
					if (writer!=null) {
						writer.close();
						Files.chmod("batchSims."+count);
					}
					count++;
					writer = new PrintWriter(new FileWriter("batchSims."+count));
				}

				// writer.println("plink --bfile pd_gwas --pheno
				// sim."+(i+1)+".dat --logistic --out sim."+(i+1)+".additive");
				writer.println("plink --bfile plink --pheno sim."+(i+1)+".dat --logistic --out sim."+(i+1)+".additive");
			}
			writer.close();
			Files.chmod("batchSims."+count);
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static void parseSimulations(final String suffix, int p_column, int chr_column, int pos_column, int markerName_column) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String filename, outfile;

		File[] files = new File(".").listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(suffix);
			}
		});
		System.out.println(ext.getTime()+"\tFound "+files.length+" files with suffix '"+suffix+"' to parse");

		for (int f = 0; f<files.length; f++) {
			filename = files[f].getName();
			outfile = filename+".parsed";

			try {
				reader = new BufferedReader(new FileReader(filename));
				writer = new PrintWriter(new FileWriter(outfile));
				while (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
					writer.println(line[markerName_column]+"\t"+line[chr_column]+"\t"+line[pos_column]+"\t"+(line[p_column].equals("NA")?".":line[p_column]));
				}
				reader.close();
				writer.close();
			} catch (Exception e) {
				System.err.println("Error parsing '"+filename+"'");
			}
		}
		System.out.println(ext.getTime()+"\tDone");

	}

	public static void procSimulations(final String suffix, int p_column, int chr_column, int pos_column, int markerName_column, String ldRoot) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
//		Hashtable<String,String> hash;
		Vector<String> markerVector = new Vector<String>(EST_NUM_MARKERS_IN_LARGEST_CHR);
		IntVector markerLocations = new IntVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
		DoubleVector pvalueVector = new DoubleVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
		int count, lowIndex, highIndex;
		IntVector indexVector = new IntVector();
		// String chr;
		double d;
		String[] markerNames;
		int[] positions;
		double[] pvals, stats;
		int[] indexSNPs, markerCounts;
		IntVector cluster;
		double maxStat;
		int maxSNP;
		String filename, outfile;
		String[] models;
		double[] indexThresholds, inclusionThresholds;
		double maxThreshold;
		int[] windowSizes;
		boolean firstOne;
		LDdatabase lddb;

		File[] files = new File(".").listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(suffix);
			}
		});
		System.out.println("Found "+files.length+" files with suffix '"+suffix+"' to compute from");

		maxThreshold = 0;
		models = new String[MODELS.length];
		indexThresholds = new double[MODELS.length];
		inclusionThresholds = new double[MODELS.length];
		windowSizes = new int[MODELS.length];
		for (int t = 0; t<MODELS.length; t++) {
			models[t] = "index"+MODELS[t][0]+"_incl"+MODELS[t][1]+"_"+(Integer.parseInt(MODELS[t][2])/1000)+"kb/";
			new File(models[t]).mkdirs();
			indexThresholds[t] = Double.parseDouble(MODELS[t][0]);
			inclusionThresholds[t] = Double.parseDouble(MODELS[t][1]);
			windowSizes[t] = Integer.parseInt(MODELS[t][2]);
			if (indexThresholds[t]>maxThreshold) {
				maxThreshold = indexThresholds[t];
			}
		}
		System.out.println("Found "+MODELS.length+" thresholds to consider (min index pvalue of "+maxThreshold+")");

		lddb = new LDdatabase(ldRoot, LDdatabase.TYPE_LONG, new Logger());

//        lddb.updateWithTheseMarkers(superset, ext.replaceDirectoryCharsWithUnderscore(dir+filename, 2));
		
		for (int chr_target = 1; chr_target<=23; chr_target++) {
			System.out.println(ext.getTime()+"\tLoading LD information for chromosome "+chr_target);
//			hash = loadLD(DEFAULT_LD_LOC+"chr"+chr_target+".recode.ped.LD");
			// hash = new Hashtable<String, String>();

			System.out.println(ext.getTime()+"\tparsing chromosome "+chr_target);

			for (int f = 0; f<files.length; f++) {
				filename = files[f].getName();
				outfile = filename+"_cluster.xln";
				// System.out.println(filename);

				try {
					reader = new BufferedReader(new FileReader(filename));
					reader.readLine();
					count = 0;
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						if (line[chr_column].equals(chr_target+"")) {
							if (!line[p_column].equals(".")&&!line[p_column].equals("NA")) {
								d = Double.parseDouble(line[p_column]);
								pvalueVector.add(d);
								if (d<maxThreshold) {
									indexVector.add(count);
								}
								markerVector.add(line[markerName_column]);
								markerLocations.add(Integer.parseInt(line[pos_column]));
								count++;
							}
						}
					}
					reader.close();
					markerNames = Array.toStringArray(markerVector);
					positions = markerLocations.toArray();
					pvals = pvalueVector.toArray();
					indexSNPs = indexVector.toArray();
					markerLocations = null;
					indexVector = null;
					markerVector = null;
					pvalueVector = null;

					for (int t = 0; t<MODELS.length; t++) {
						writer = new PrintWriter(new FileWriter(models[t]+outfile, chr_target>1));
						if (chr_target==1) {
							writer.println("Marker\tChr\tPostition\tWeightedStatistic\tNumMarkers\tUnweightedStatistic\tUCSC coordinates");
						}

						stats = new double[indexSNPs.length];
						markerCounts = new int[indexSNPs.length];
						//
						// System.out.println("Analyzing "+indexSNPs.length+"
						// index SNPs on chromosome "+chr_target);

						cluster = new IntVector();
						firstOne = true;
						for (int i = 0; i<indexSNPs.length; i++) {
							if (indexSNPs[i]==-1) {
								System.err.println("Error - model "+t+" indexSNP "+i);
							}
							if (pvals[indexSNPs[i]]<indexThresholds[t]) {
								lowIndex = indexSNPs[i];
								while (lowIndex>=0&&positions[lowIndex]>positions[indexSNPs[i]]-windowSizes[t]) {
									lowIndex--;
								}
								lowIndex++;
								highIndex = indexSNPs[i];
								while (highIndex<positions.length&&positions[highIndex]<positions[indexSNPs[i]]+windowSizes[t]) {
									highIndex++;
								}
								highIndex--;
								stats[i] = computeStatistic(markerNames, pvals, lowIndex, highIndex, lddb.getLongChrLDdb(chr_target), inclusionThresholds[t], null)[0];
								markerCounts[i] = highIndex-lowIndex+1;
								//
								// System.out.println((i+1)+")
								// "+markerNames[indexSNPs[i]]+"\t"+positions[indexSNPs[i]]+"\t"+ext.formDeci(stats[i],
								// 2)+" (using "+(highIndex-lowIndex+1)+"
								// markers)");

								if (firstOne) {
									cluster.add(i);
									firstOne = false;
								} else if (positions[indexSNPs[i]]-positions[indexSNPs[i-1]]<windowSizes[t]) {
									cluster.add(i);
								} else {
									maxStat = 0;
									maxSNP = -1;
									for (int j = 0; j<cluster.size(); j++) {
										if (stats[cluster.elementAt(j)]>maxStat) {
											maxSNP = cluster.elementAt(j);
											maxStat = stats[maxSNP];
										}
									}
									writer.println(markerNames[indexSNPs[maxSNP]]+"\t"+chr_target+"\t"+positions[indexSNPs[maxSNP]]+"\t"+ext.formDeci(stats[maxSNP], 2)+"\t"+markerCounts[maxSNP]+"\t"+ext.formDeci(stats[maxSNP]*Math.sqrt(markerCounts[maxSNP]), 2)+"\tchr"+chr_target+":"+(positions[indexSNPs[maxSNP]]-windowSizes[t])+"-"+(positions[indexSNPs[maxSNP]]+windowSizes[t]));
									cluster.clear();
									cluster.add(i);
								}

								if (i==indexSNPs.length-1) {
									maxStat = 0;
									maxSNP = -1;
									for (int j = 0; j<cluster.size(); j++) {
										if (stats[cluster.elementAt(j)]>maxStat) {
											maxSNP = cluster.elementAt(j);
											maxStat = stats[maxSNP];
										}
									}
									writer.println(markerNames[indexSNPs[maxSNP]]+"\t"+chr_target+"\t"+positions[indexSNPs[maxSNP]]+"\t"+ext.formDeci(stats[maxSNP], 2)+"\t"+markerCounts[maxSNP]+"\t"+ext.formDeci(stats[maxSNP]*Math.sqrt(markerCounts[maxSNP]), 2)+"\tchr"+chr_target+":"+(positions[indexSNPs[maxSNP]]-windowSizes[t])+"-"+(positions[indexSNPs[maxSNP]]+windowSizes[t]));
								}

							}
						}

						writer.close();
					}

					markerNames = null;
					positions = null;
					pvals = null;
					indexSNPs = null;
					markerVector = new Vector<String>(EST_NUM_MARKERS_IN_LARGEST_CHR);
					markerLocations = new IntVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
					pvalueVector = new DoubleVector(EST_NUM_MARKERS_IN_LARGEST_CHR);
					indexVector = new IntVector();

				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""+filename+"\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \""+filename+"\"");
					System.exit(2);
				}
			}
		}
	}

	public static void formDistributions(int column, int levelsDeep) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp;
		int count;
		File[] files, dirs;
		DoubleVector[][] dvs;
		DoubleVector trav;
		int[] keys;
		String col;

		dirs = new File(".").listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return file.isDirectory()&&filename.startsWith("index");
			}
		});

		dvs = new DoubleVector[dirs.length][];
		for (int i = 0; i<dirs.length; i++) {
			files = dirs[i].listFiles(new FilenameFilter() {
				public boolean accept(File file, String filename) {
					return filename.endsWith("parsed_cluster.xln");
				}
			});
			dvs[i] = DoubleVector.newDoubleVectors(levelsDeep);
			for (int j = 0; j<files.length; j++) {
				trav = new DoubleVector();
				try {
					reader = new BufferedReader(new FileReader(files[j]));
					reader.readLine();
					while (reader.ready()) {
						line = reader.readLine().split("[\\s]+");
						trav.add(Double.parseDouble(line[column]));
					}
					reader.close();
					keys = Sort.quicksort(trav, Sort.DESCENDING);
					for (int k = 0; k<levelsDeep&&trav.size()>k; k++) {
						dvs[i][k].add(trav.elementAt(keys[k]));
					}
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""+files[j].getName()+"\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \""+files[j].getName()+"\"");
					System.exit(2);
				}
			}
		}

		try {
			writer = new PrintWriter(new FileWriter("NRSS_distributions.xln"));
			for (int i = 0; i<dirs.length; i++) {
				writer.print((i==0?"":"\t")+dirs[i].getName()+Array.toStr(Array.stringArray(levelsDeep), "\t"));
			}
			writer.println();
			for (int i = 0; i<dirs.length; i++) {
				writer.print((i==0?"":"\t")+Array.toStr(Array.stringArraySequence(levelsDeep, ""), "\t"));
			}
			writer.println();
			writer.println(Array.toStr(Array.stringArray(dirs.length*levelsDeep, "0.5")));
			temp = "";
			for (int i = 0; i<dirs.length; i++) {
				for (int j = 0; j<levelsDeep; j++) {
					col = ext.getExcelColumn(i*levelsDeep+j);
					temp += (i==0&&j==0?"":"\t")+"=COUNTIF("+col+"5:"+col+(Math.max(4+dvs[i][j].size(), 5))+", \">\"&"+col+"3)/COUNT("+col+"5:"+col+(Math.max(4+dvs[i][j].size(), 5))+")";
				}
			}

			count = 0;
			while (!temp.equals(Array.toStr(Array.stringArray(dirs.length*levelsDeep)))) {
				writer.println(temp);
				temp = "";
				for (int i = 0; i<dirs.length; i++) {
					for (int j = 0; j<levelsDeep; j++) {
						if (dvs[i][j].size()>count) {
							temp += (i==0&&j==0?"":"\t")+dvs[i][j].elementAt(count);
						} else {
							temp += (i==0&&j==0?"":"\t");
						}
					}
				}
				count++;
			}

			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void prepFile(String filename) {
		BufferedReader reader;
        PrintWriter writer;
        String[] line;
        SnpMarkerSet markerSet;
        
        System.out.println("Loading marker set...");
        markerSet = new SnpMarkerSet(filename, SnpMarkerSet.GENERIC_FORMAT_ANNOTATED_IGNORE_FIRST_LINE, true, new Logger());
        System.out.println("Sorting markers...");
        markerSet.sortMarkers();
        System.out.println("Writing sorted list...");
        markerSet.writeToFile(ext.rootOf(filename, false)+"_sorted.xln", SnpMarkerSet.GENERIC_FORMAT_ANNOTATED);
        
        System.out.println("Writing filtered list...");
        try {
	        reader = new BufferedReader(new FileReader(ext.rootOf(filename, false)+"_sorted.xln"));
	        writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_input.xln"));
	        writer.println("SNP\tChr\tPositon\tpvalue");
	        while (reader.ready()) {
	        	line = reader.readLine().trim().split("[\\s]+");
	        	if (Integer.parseInt(line[2]) > 0) {
	        		writer.println(Array.toStr(line));
	        	}
	        }
	        reader.close();
            writer.close();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+ext.rootOf(filename, false)+"_sorted.xln"+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+ext.rootOf(filename, false)+"_sorted.xln"+"\"");
	        System.exit(2);
        }
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		int batch = 0;
		int procMap = -1;
		int simulate = -1;
		boolean parseSims = false;
		boolean procSims = false;
		boolean formDist = false;
		boolean prep = false;
		String ldRoot = LDdatabase.MASTER_HAPMAP_ROOT;
		int genomeBuild = 37;

		// String filename = "plink.assoc.logistic.flipped";
		// String dir = "C:\\Documents and Settings\\npankrat\\My
		// Documents\\_GAW16\\results\\goodMarkers\\replicate\\";
		// String filename = "additive.included";
		// int col = 8;

		// String dir = "";
		// String filename = "additive.mafs.prn";
		// int col = 3;
		// int var = 4;

		// String dir = "C:\\Documents and Settings\\npankrat\\My
		// Documents\\gwas\\udall\\";
		// String filename = "MetalResults.xln";
		// int col = 7;
		// int var = -1;

//		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\Renewal\\progress report\\2q NRSS\\";
//		String filename = "SFH2q_final_perm.assoc.logistic.xln";
//		int col = 4;
//		int var = -1;

//		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\Consortium\\analysisOfImputation\\Aff\\filtered2\\";
//		String filename = "parseNRSSinput.xln";
//		int col = 3;
//		int var = -1;

		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\Consortium\\analysisOfImputation\\Aff_AAO_combo\\finalList\\";
		String filename = "v2_nrss_input.txt";
		int col = 3;
		int var = -1;

		
		
		String usage = "\\n"+
		"park.gwa.Nrss requires 0-1 arguments\n"+
		"   (1) filename (i.e. file="+filename+" (default))\n"+
		"   (1) column of p-values (i.e. col="+col+" (default))\n"+
		"   (1) column of variate (i.e. var="+var+" (default))\n"+
		"   (2) proc map for Haploview (i.e. procMap=4 for chromosome 4 (not the default))\n"+
		"   (3) batch haploview LD computation, using N files (i.e. batch=N (not the default))\n"+
		"   (4) simulate null distribution (i.e. sim=100 to do 100 replicates (not the default))\n"+
		"   (5) the genome build to determine genes (i.e. build="+genomeBuild+" (default))\n"+
		"   (3) parse simulations (i.e. -parseSims (not the default))\n"+
		"   (3) process simulation results (i.e. -procSims (not the default))\n"+
		"   (3) form distributions (i.e. -formDist (not the default))\n"+
	    "   (3) root of PLINK file to test for LD (i.e. ldRoot="+ldRoot+" (default))\n"+
	    "   (3) prep file for input (sort markers and remove unknown positions (i.e. -prep (not the default))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("col=")) {
				col = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("var=")) {
				var = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
		    } else if (args[i].startsWith("ldRoot=")) {
		    	ldRoot = ext.parseStringArg(args[i], null);
			    numArgs--;
			} else if (args[i].startsWith("procMap=")) {
				procMap = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("sim=")) {
				simulate = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("batch=")) {
				batch = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("-parseSims")) {
				parseSims = true;
				numArgs--;
			} else if (args[i].startsWith("-procSims")) {
				procSims = true;
				numArgs--;
			} else if (args[i].startsWith("-formDist")) {
				formDist = true;
				numArgs--;
			} else if (args[i].startsWith("-prep")) {
				prep = true;
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (prep) {
				prepFile(dir+filename);
			} else if (batch>0) {
				batchHaploview(batch);
			} else if (procMap>0) {
				procMapForHaploview("chr"+procMap+".recode.map", "chr"+procMap+".map");
			} else if (simulate>0) {
				// simulateNull("pd_gwas.fam", simulate);
				simulateNull("plink.fam", simulate);
			} else if (parseSims) {
				parseSimulations("logistic", 8, 0, 2, 1);
			} else if (procSims) {
				// procSimulations("logistic.parsed", 8, 0, 2, 1);
				// procSimulations("logistic.parsed", 3, 1, 2, 0);
				// procSimulations(".final", 3, 1, 2, 0);
				procSimulations(".included", 3, 1, 2, 0, ldRoot);
			} else if (formDist) {
				formDistributions(3, 5);
			} else {
				// procRegion("C:\\Documents and Settings\\npankrat\\My
				// Documents\\gwas\\postHocs\\nonredundency\\SNCA.results.txt",
				// 3, "C:\\Documents and Settings\\npankrat\\My
				// Documents\\gwas\\postHocs\\nonredundency\\SNCA.recode.ped.LD");
				// procFile("chr4.xls", 6);
				// procFile("SNCA.results.txt", 3);
				// procFile("logistic.xls", 6);
				// procFile("3plus.prn", 6);
				// procFile("sim.1.additive.assoc.logistic.parsed", 3);
				// procFile("sim.1.additive.assoc.logistic.parsed", 3);

				// procFile(dir, filename, 7, var, DEFAULT_WINDOW,
				// "IU_Add_Nrss.xln");
				// procFile(dir, filename, 10, var, DEFAULT_WINDOW,
				// "Udall_Add_Nrss.xln");
				// procFile(dir, filename, 11, var, DEFAULT_WINDOW,
				// "Meta_Add_Nrss.xln");
				// procFile(dir, filename, 13, var, DEFAULT_WINDOW,
				// "Meta_Dom_Nrss.xln");
				// procFile(dir, filename, 15, var, DEFAULT_WINDOW,
				// "Meta_Rec_Nrss.xln");

				procFile(dir, filename, col, var, ldRoot, DEFAULT_WINDOW, "2q_Add.xln", genomeBuild);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
