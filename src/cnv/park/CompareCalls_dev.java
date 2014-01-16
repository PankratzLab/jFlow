package cnv.park;

import java.io.*;
import java.util.*;

import cnv.filesys.MarkerSet;
import cnv.var.CNVariant;
import common.*;
import stats.Maths;

public class CompareCalls_dev {
	public static final String DEFAULT_ROOT = "C:/data/pennComp/";
	public static final String DEFAULT_CNV_FILE = "all_gw6.cnv";
	public static final String[] DEFAULT_COMP_LISTS = { "roots1.txt", "dups1.txt" };
	public static final String DEFAULT_LRR_LOOKUP = "LRR_full_samps.txt";
	public static final Double DEFAULT_LRR_FILTER = .30;
	public static final Double DEFAULT_CONF_FILTER = 10.0;
	public static final Integer DEFAULT_NUM_MARKERS_FILTER = 20;
	private int numberComputed;
	private double averagePercent;
	private int totalCallsAnalyzed;
	private String[] cnvFiles;
	private int[] exactMatches;
	private int[] sigOlapMatches;
	private int[][] indCNVCalls;
	private String[] inds;

	// public static final String[] DEFAULT_FILES = {"conf.cnv", "allMarkers.cnv"};
	// public static final String[] DEFAULT_FILES = {"conf_100kb_5SNP_10.0.cnv", "allMarkers_100kb_5SNP_10.0.cnv"};
	public static final String[] DEFAULT_FILES = { "conf_100kb_5SNP_10.0.cnv", "conf_100kb_20SNP_10.0.cnv" };

	public CompareCalls_dev(int numberComputed, double averagePercent, int totalCallsAnalyzed, String[] cnvFiles, int[] exactMatches, int[] sigOlapMatches, String[] inds) {
		this.numberComputed = numberComputed;
		this.averagePercent = averagePercent;
		this.totalCallsAnalyzed = totalCallsAnalyzed;
		this.cnvFiles = cnvFiles;
		this.exactMatches = exactMatches;
		this.sigOlapMatches = sigOlapMatches;
		this.inds = inds;
	}

	public CompareCalls_dev() {
		this.numberComputed = 0;
		this.averagePercent = 0;

	}

	public int getNumberComputed() {
		return numberComputed;
	}

	public double getAveragePercent() {
		return averagePercent;
	}

	public static CompareCalls_dev compare(String rootDir, String[] files) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, inds;
		Hashtable<String, Hashtable<String, Vector<CNVariant>>> hash = new Hashtable<String, Hashtable<String, Vector<CNVariant>>>();
		Hashtable<String, Vector<CNVariant>> source = new Hashtable<String, Vector<CNVariant>>();
		CNVariant[][] cnvs;
		Vector<CNVariant> v = new Vector<CNVariant>();
		int match;
		int numPassingAndPresent = 0;
		int totalCallsAnalyzed = 0;
		int[] counts;
		int[][] allPossibleCombinations = Maths.getIndicesForAllCombinations(files.length, 2);
		int[] exactMatches;
		int[] sigOlapMatches;
		String[] comparedIDs;
		double goodCallPPercent = 0;
		double averageGoodCallPercent = 0;
		CompareCalls_dev comparedCalls;
		for (int i = 0; i < files.length; i++) {
			try {
				reader = new BufferedReader(new FileReader(files[i]));
				if (!ext.checkHeader(reader.readLine().trim().split("[\\s]+"), CNVariant.PLINK_CNV_HEADER, false)) {
					reader.close();
					System.err.println("quitting comparison");

				}
				while (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
					if (hash.containsKey(line[0] + "\t" + line[1])) {
						source = hash.get(line[0] + "\t" + line[1]);
					} else {
						hash.put(line[0] + "\t" + line[1], source = new Hashtable<String, Vector<CNVariant>>());
					}
					if (source.containsKey(i + "")) {
						v = source.get(i + "");
					} else {
						source.put(i + "", v = new Vector<CNVariant>());
					}
					v.add(new CNVariant(line, i));
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + files[i] + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + files[i] + "\"");
				System.exit(2);
			}
		}
		for (int i = 0; i < files.length; i++) {

		}
		inds = HashVec.getKeys(hash);
		exactMatches = new int[inds.length];
		sigOlapMatches = new int[inds.length];
		comparedIDs = new String[inds.length];
		for (int i = 0; i < allPossibleCombinations.length; i++) {
			goodCallPPercent = 0;
			numPassingAndPresent = 0;
			System.out.println(i);

			try {

				// writer = new PrintWriter(new FileWriter(rootDir + "Compare" + ext.rootOf(files[allPossibleCombinations[i][0]]) + "And" + ext.rootOf(files[allPossibleCombinations[i][1]]) + ".xln"));
				// writer.println("FID\tIID\tTotal" + ext.rootOf(files[allPossibleCombinations[i][0]]) + "\tTotal" + ext.rootOf(files[allPossibleCombinations[i][1]]) + "\tUnique" + ext.rootOf(files[allPossibleCombinations[i][0]]) + "\tUnique" + ext.rootOf(files[allPossibleCombinations[i][1]]) + "\tOverlaps\tExactMatches\tsigOverlap\tOverlapsPlusSigOverLapPlusExact\tSigOverLapPlusExact");
				for (int j = 0; j < inds.length; j++) {
					// System.out.println(j);
					cnvs = new CNVariant[][] { CNVariant.toCNVariantArray(hash.get(inds[j]).get(allPossibleCombinations[i][0] + "")), CNVariant.toCNVariantArray(hash.get(inds[j]).get(allPossibleCombinations[i][1] + "")) };
					counts = new int[5];
					if (cnvs[0].length == 0) {
						// j++;
						// System.err.println("Warning - " + inds[j] + " not found in " + files[allPossibleCombinations[i][0]] + ", comparing what is left");
						// continue;
					}
					if (cnvs[1].length == 0) {
						// j++;
						// System.err.println("Warning - " + inds[j] + " not found in " + files[allPossibleCombinations[i][1]] + ", comparing what is left");
						// continue;
					}

					for (int a = 0; a < cnvs[0].length; a++) {
						match = 0;
						for (int b = 0; b < cnvs[1].length; b++) {
							if (cnvs[0][a].equals(cnvs[1][b])) {

								match = 3;
								cnvs[1][b].setSource(99);
							} else if (match < 2 && cnvs[0][a].significantOverlap(cnvs[1][b])) {
								match = 4;
								cnvs[1][b].setSource(99);
							} else if (match < 2 && cnvs[0][a].overlaps(cnvs[1][b])) {
								match = 2;
								cnvs[1][b].setSource(99);
							}
						}
						counts[match]++;
					}
					for (int b = 0; b < cnvs[1].length; b++) {
						match = 1;
						for (int a = 0; a < cnvs[0].length; a++) {
							if (cnvs[1][b].getSource() != 99 && cnvs[1][b].equals(cnvs[0][a])) {
								match = 3;
							} else if (match < 2 && cnvs[1][b].getSource() != 99 && cnvs[1][b].significantOverlap(cnvs[0][a])) {
								match = 4;
							} else if (match < 2 && cnvs[1][b].getSource() != 99 && cnvs[1][b].overlaps(cnvs[0][a])) {
								match = 2;
							}
						}
						if (cnvs[1][b].getSource() != 99) {
							counts[match]++;
						}
					}
					exactMatches[j] = counts[3];
					sigOlapMatches[j] = counts[4];
					comparedIDs[j] = inds[j].replaceAll("\t.*", "");

					int oLapSigoLapExactTotal = counts[2] + counts[3] + counts[4];
					int SigoLapExactTotal = counts[3] + counts[4];
					int ComparisonTotalCalls = cnvs[0].length + cnvs[1].length;
					// System.out.println(j + "\t" + comparedIDs[j] + "\t" + exactMatches[j] + "\t" + sigOlapMatches[j] + "\t" + SigoLapExactTotal);
					totalCallsAnalyzed += ComparisonTotalCalls;

					if (cnvs[0].length > 0 && cnvs[1].length > 0) {
						numPassingAndPresent++;
						goodCallPPercent += ((double) (2 * SigoLapExactTotal) / ComparisonTotalCalls);
					}
					// goodCallPPercent /= numTotalCalls;
					// writer.println(inds[j] + "\t" + cnvs[0].length + "\t" + cnvs[1].length + "\t" + Array.toStr(counts) + "\t" + oLapSigoLapExactTotal + "\t" + SigoLapExactTotal);
				}
				// writer.close();
				averageGoodCallPercent = goodCallPPercent / numPassingAndPresent;
				// System.out.println(goodCallPPercent / numPassingAndPresent);
			} catch (Exception e) {
				System.err.println("Error comparing " + files[allPossibleCombinations[i][0]] + " and " + files[allPossibleCombinations[i][1]]);
				e.printStackTrace();
			}
		}
		// System.out.println(0 + "Hellooooooo\t" + comparedIDs[1] + "\t" + exactMatches[1] + "\t" + sigOlapMatches[1]);

		comparedCalls = new CompareCalls_dev(numPassingAndPresent, averageGoodCallPercent, totalCallsAnalyzed, files, exactMatches, sigOlapMatches, comparedIDs);
		return comparedCalls;
	}

	public static String[] filterCNVs(String rootDirectory, String cnvFile, String LRR_lookup, String[] compFiles, double lrrFilter, double confFilter, int numMarkers) {
		Hashtable<String, Double> lrrs = new Hashtable<String, Double>();
		Hashtable<String, Integer> defineCompHash = new Hashtable<String, Integer>();
		BufferedReader reader;

		PrintWriter[] writers = new PrintWriter[compFiles.length];
		String[] line;
		String[] files = new String[compFiles.length];

		int[] filteredCounts = new int[compFiles.length];
		int[] includeCounts = new int[compFiles.length];
		int totalCounts = 0;

		try {
			lrrs = getLRRs(rootDirectory, LRR_lookup);
			defineCompHash = defineCompLists(rootDirectory, compFiles);
			reader = new BufferedReader(new FileReader(rootDirectory + cnvFile));
			if (!ext.checkHeader(reader.readLine().trim().split("[\\s]+"), CNVariant.PLINK_CNV_HEADER, false)) {
				reader.close();
				System.err.println("quitting comparison");
			} else {
				for (int i = 0; i < compFiles.length; i++) {
					files[i] = rootDirectory + "intermediateFiles/Filtered_LRR_" + lrrFilter + "_conf_" + confFilter + "_numMarkers_" + numMarkers + "_" + compFiles[i] + ".cnv";
					writers[i] = new PrintWriter(new FileWriter(files[i]));
					writers[i].println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
					filteredCounts[i] = 0;
				}
				while (reader.ready()) {
					line = reader.readLine().trim().split("\t");
					totalCounts++;
					String id = line[0].replaceAll("_1", "");
					if (defineCompHash.containsKey(id) && defineCompHash.containsKey(id + "_1")) {
						int fileNumber = defineCompHash.get(line[0]);

						if (lrrs.get(line[0]) <= lrrFilter && Double.parseDouble(line[6]) >= confFilter && Integer.parseInt(line[7]) >= numMarkers) {
							// TODO
							// this is very specific and only when the ids are not identical...in this case _1 is a duplicate;
							// System.out.println(line[0]);
							line[0] = line[0].replaceAll("_1", "");
							writers[fileNumber].println(Array.toStr(line));
							includeCounts[fileNumber]++;
						} else {
							filteredCounts[fileNumber]++;
						}
					}
				}
				reader.close();
				for (int i = 0; i < compFiles.length; i++) {
					writers[i].close();
					int totalInCompFile = includeCounts[i] + filteredCounts[i];
					System.out.println("Info - Out of a total of " + totalCounts + " cnv calls in " + rootDirectory + cnvFile + " , a total of " + totalInCompFile + " matched ids listed in " + compFiles[i] + ", " + filteredCounts[i] + " were filtered out, " + includeCounts[i] + " were included for comparison");
				}
			}
		} catch (Exception e) {
			System.err.println("Error comparing ");
			e.printStackTrace();
		}
		return files;

	}

	public static Hashtable<String, Double> getLRRs(String rootDirectory, String LRR_lookup) {
		Hashtable<String, Double> lrrs = new Hashtable<String, Double>();
		BufferedReader reader;
		String[] line;
		try {
			reader = new BufferedReader(new FileReader(rootDirectory + LRR_lookup));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (lrrs.containsKey(line[0])) {
					System.err.println("Warning - duplicate samples found in " + rootDirectory + LRR_lookup + " only using one of them");
				} else {
					lrrs.put(line[0], Double.parseDouble(line[1]));
					// System.out.println(line[0] + "\t" + line[1]);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + rootDirectory + LRR_lookup + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + rootDirectory + LRR_lookup + "\"");
			System.exit(2);
		}
		return lrrs;
	}

	public static Hashtable<String, Integer> defineCompLists(String rootDirectory, String[] compFiles) {
		Hashtable<String, Integer> defineCompHash = new Hashtable<String, Integer>();
		BufferedReader reader;
		String[] line;
		if (compFiles.length <= 1) {
			System.err.println("Error - two files of IDs are required for comparison...filtering on " + compFiles.length + " files");
		}

		for (int i = 0; i < compFiles.length; i++) {
			try {
				reader = new BufferedReader(new FileReader(rootDirectory + compFiles[i]));
				while (reader.ready()) {
					line = reader.readLine().trim().split("/t");
					if (defineCompHash.containsKey(line[0])) {
						System.err.println("Warning - duplicate samples found in " + rootDirectory + compFiles[i] + " only using one of them");
					} else {
						defineCompHash.put(line[0], i);
						// System.out.println(line[0]);
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + rootDirectory + compFiles[i] + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + rootDirectory + compFiles[i] + "\"");
				System.exit(2);
			}
		}

		return defineCompHash;
	}

	private static double[] binIt(double startVal, double stopVal, int numBins) {
		double[] values = new double[numBins + 1];
		double inc = (stopVal - startVal) / numBins;
		for (int i = 0; i < numBins + 1; i++) {
			// System.out.println(inc * i + startVal);
			values[i] = (inc * i) + startVal;
		}
		return values;
	}

	public static void iterate(double lrrstart, double lrrstop, int numLrrBins, double confstart, double confstop, int numConfBins, int probestart, int probeStop, String rootDirectory, String cnvFile, String LRR_lookup, String[] compFiles) {
		PrintWriter probeWriter;
		PrintWriter confWriter;
		PrintWriter lrrWriter;
		PrintWriter lrrProbeWriter;
		CompareCalls_dev comparedCalls;
		double[] lrrValues = new double[numLrrBins];
		double[] confValues = new double[numConfBins];
		lrrValues = binIt(lrrstart, lrrstop, numLrrBins);
		confValues = binIt(confstart, confstop, numConfBins);
		long time;
		time = new Date().getTime();
		String commonHeader = "_Value\tConcordance_Percent(((2X+(Exact+sigolap)/(totalcallsFile1+totalcallsFIle2))/numindividuals\tnumindividual\ttotalCallsPassingFilter\tFile1\tFile2";
		try {
			probeWriter = new PrintWriter(new FileWriter(rootDirectory + "probeConcordance_" + probestart + "_" + probeStop + "LRR" + DEFAULT_LRR_FILTER + "Conf" + DEFAULT_CONF_FILTER + ".concord"));
			confWriter = new PrintWriter(new FileWriter(rootDirectory + "confConcordance_" + confstart + "_" + confstop + "LRR" + DEFAULT_LRR_FILTER + "numMarkers" + DEFAULT_NUM_MARKERS_FILTER + ".concord"));
			lrrWriter = new PrintWriter(new FileWriter(rootDirectory + "lrrConcordance_" + lrrstart + "_" + lrrstop + "numMarkers" + DEFAULT_NUM_MARKERS_FILTER + "Conf" + DEFAULT_CONF_FILTER + ".concord"));
			// probeWriter.println("numMarker" + commonHeader);
			confWriter.println("conf" + commonHeader);
			lrrWriter.println("lrr_SD" + commonHeader);
			lrrProbeWriter = new PrintWriter(new FileWriter(rootDirectory + "lrrProbetestConcordance.concord"));
			lrrProbeWriter.println("LRR_Value\tnumMarker" + commonHeader);
			for (int i = probestart; i < probeStop + 1; i++) {
				comparedCalls = compare(rootDirectory, filterCNVs(rootDirectory, cnvFile, LRR_lookup, compFiles, DEFAULT_LRR_FILTER, DEFAULT_CONF_FILTER, i));
				if (i == probestart) {
					probeWriter.println("numMarker" + commonHeader + "\t" + Array.toStr(comparedCalls.getInds()) + "\t" + Array.toStr(comparedCalls.getInds()));
				}
				probeWriter.println(i + "\t" + comparedCalls.getAveragePercent() + "\t" + comparedCalls.getNumberComputed() + "\t" + comparedCalls.getTotalCallsAnalyzed() + "\t" + Array.toStr(comparedCalls.getCnvFiles()) + "\t" + Array.toStr(comparedCalls.getExactMatches()) + "\t" + Array.toStr(comparedCalls.getSigOlapMatches()));
			}
			// for (int i = 0; i < lrrValues.length; i++) {
			// comparedCalls = compare(rootDirectory, filterCNVs(rootDirectory, cnvFile, LRR_lookup, compFiles, lrrValues[i], DEFAULT_CONF_FILTER, DEFAULT_NUM_MARKERS_FILTER));
			// lrrWriter.println(lrrValues[i] + "\t" + comparedCalls.getAveragePercent() + "\t" + comparedCalls.getNumberComputed() + "\t" + comparedCalls.getTotalCallsAnalyzed() + "\t" + Array.toStr(comparedCalls.getCnvFiles()));
			// }
			// for (int i = 0; i < confValues.length; i++) {
			// comparedCalls = compare(rootDirectory, filterCNVs(rootDirectory, cnvFile, LRR_lookup, compFiles, DEFAULT_LRR_FILTER, confValues[i], DEFAULT_NUM_MARKERS_FILTER));
			// confWriter.println(confValues[i] + "\t" + comparedCalls.getAveragePercent() + "\t" + comparedCalls.getNumberComputed() + "\t" + comparedCalls.getTotalCallsAnalyzed() + "\t" + Array.toStr(comparedCalls.getCnvFiles()));
			// }
			// for (int j = 0; j < lrrValues.length; j++) {
			// for (int i = probestart; i < probeStop + 1; i++) {
			// comparedCalls = compare(rootDirectory, filterCNVs(rootDirectory, cnvFile, LRR_lookup, compFiles, lrrValues[j], DEFAULT_CONF_FILTER, i));
			// lrrProbeWriter.println(lrrValues[j] + "\t" + i + "\t" + comparedCalls.getAveragePercent() + "\t" + comparedCalls.getNumberComputed() + "\t" + comparedCalls.getTotalCallsAnalyzed() + "\t" + Array.toStr(comparedCalls.getCnvFiles()));
			// }
			// }

			probeWriter.close();
			confWriter.close();
			lrrWriter.close();
			lrrProbeWriter.close();
			System.out.println("Iterations took " + ext.getTimeElapsed(time));
		} catch (Exception e) {
			System.err.println("Error comparing ");
			e.printStackTrace();
		}

	}

	public int[] getExactMatches() {
		return exactMatches;
	}

	public void setExactMatches(int[] exactMatches) {
		this.exactMatches = exactMatches;
	}

	public int[] getSigOlapMatches() {
		return sigOlapMatches;
	}

	public void setSigOlapMatches(int[] sigOlapMatches) {
		this.sigOlapMatches = sigOlapMatches;
	}

	public int[][] getIndCNVCalls() {
		return indCNVCalls;
	}

	public void setIndCNVCalls(int[][] indCNVCalls) {
		this.indCNVCalls = indCNVCalls;
	}

	public String[] getInds() {
		return inds;
	}

	public void setInds(String[] inds) {
		this.inds = inds;
	}

	public void setNumberComputed(int numberComputed) {
		this.numberComputed = numberComputed;
	}

	public void setAveragePercent(double averagePercent) {
		this.averagePercent = averagePercent;
	}

	public int getTotalCallsAnalyzed() {
		return totalCallsAnalyzed;
	}

	public void setTotalCallsAnalyzed(int totalCallsAnalyzed) {
		this.totalCallsAnalyzed = totalCallsAnalyzed;
	}

	public String[] getCnvFiles() {
		return cnvFiles;
	}

	public void setCnvFiles(String[] cnvFiles) {
		this.cnvFiles = cnvFiles;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String rootDirectory = DEFAULT_ROOT;
		String cnvFile = DEFAULT_CNV_FILE;
		String[] files = DEFAULT_FILES;
		// right now we take two lists of ids to compare...since these are duplicates with unique ids
		String[] compFiles = DEFAULT_COMP_LISTS;
		String LRR_lookup = DEFAULT_LRR_LOOKUP;
		double lrrFilter = DEFAULT_LRR_FILTER;
		double confFilter = DEFAULT_CONF_FILTER;
		int numMarkers = DEFAULT_NUM_MARKERS_FILTER;
		String usage = "\\n" + "park.cnv.ComparePlinkResults requires 0-1 arguments\n" + "   (1) directory (i.e. dir=" + rootDirectory + " (default))\n" + "   (2) files to be compared (i.e. files=" + Array.toStr(files, ",") + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				rootDirectory = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("files=")) {
				files = args[i].split("=")[1].split(",");
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			// compare(rootDirectory, files);
			// defineCompLists(rootDirectory, compLists);
			// getLRRs(rootDirectory, LRR_lookup);
			iterate(0.2, 0.6, 100, 0.0, 100, 100, 0, 100, rootDirectory, cnvFile, LRR_lookup, compFiles);
			// compare(rootDirectory, filterCNVs(rootDirectory, cnvFile, LRR_lookup, compFiles, lrrFilter, confFilter, numMarkers));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
