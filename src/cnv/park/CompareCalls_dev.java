package cnv.park;

import java.io.*;
import java.util.*;

import cnv.var.CNVariant;
import common.*;
import stats.Maths;

public class CompareCalls_dev {
	public static final String DEFAULT_ROOT = "C:/data/pennComp/";
	public static final String DEFAULT_MARKER_FILE = "affygw6.hg18.pfb";
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
	private double[] indPercentMatches;
	private int[][] indCNVCalls;
	private String[] inds;

	// public static final String[] DEFAULT_FILES = {"conf.cnv", "allMarkers.cnv"};
	// public static final String[] DEFAULT_FILES = {"conf_100kb_5SNP_10.0.cnv", "allMarkers_100kb_5SNP_10.0.cnv"};
	public static final String[] DEFAULT_FILES = { "conf_100kb_5SNP_10.0.cnv", "conf_100kb_20SNP_10.0.cnv" };

	public CompareCalls_dev(int numberComputed, double averagePercent, int totalCallsAnalyzed, String[] cnvFiles, int[] exactMatches, int[] sigOlapMatches, String[] inds, int[][] indCNVCalls, double[] indPercentMatches) {
		this.numberComputed = numberComputed;
		this.averagePercent = averagePercent;
		this.totalCallsAnalyzed = totalCallsAnalyzed;
		this.cnvFiles = cnvFiles;
		this.exactMatches = exactMatches;
		this.sigOlapMatches = sigOlapMatches;
		this.inds = inds;
		this.indCNVCalls = indCNVCalls;
		this.indPercentMatches = indPercentMatches;
	}

	public double[] getindPercentMatches() {
		return indPercentMatches;
	}

	public void setindPercentMatches(double[] indPercentMatches) {
		this.indPercentMatches = indPercentMatches;
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
		String[] line, inds, comparedIDs;
		Hashtable<String, Hashtable<String, Vector<CNVariant>>> hash = new Hashtable<String, Hashtable<String, Vector<CNVariant>>>();
		Hashtable<String, Vector<CNVariant>> source = new Hashtable<String, Vector<CNVariant>>();
		CNVariant[][] cnvs;
		Vector<CNVariant> v = new Vector<CNVariant>();
		int match;
		int numPassingAndPresent = 0;
		int totalCallsAnalyzed = 0;
		int[] counts;
		int[] exactMatches;
		int[] sigOlapMatches;
		int[][] indCNVCalls;
		int[][] allPossibleCombinations = Maths.getIndicesForAllCombinations(files.length, 2);
		double[] indPercentMatches;
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
		inds = HashVec.getKeys(hash);
		exactMatches = new int[inds.length];
		sigOlapMatches = new int[inds.length];
		comparedIDs = new String[inds.length];
		indPercentMatches = new double[inds.length];
		indCNVCalls = new int[2][inds.length];

		for (int i = 0; i < allPossibleCombinations.length; i++) {
			goodCallPPercent = 0;
			numPassingAndPresent = 0;
			try {
				for (int j = 0; j < inds.length; j++) {
					cnvs = new CNVariant[][] { CNVariant.toCNVariantArray(hash.get(inds[j]).get(allPossibleCombinations[i][0] + "")), CNVariant.toCNVariantArray(hash.get(inds[j]).get(allPossibleCombinations[i][1] + "")) };
					counts = new int[5];
					// we handle samples with missing calls later
					for (int a = 0; a < cnvs[0].length; a++) {
						match = 0;
						for (int b = 0; b < cnvs[1].length; b++) {
							if (cnvs[0][a].equals(cnvs[1][b])) {
								match = 3;
								cnvs[1][b].setSource(99);
							} else if (match < 2 && cnvs[0][a].significantOverlap(cnvs[1][b])) {
								match = 4;
								// an overlap is not assumed in both directions
								// cnvs[1][b].setSource(99);
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
					indCNVCalls[0][j] = cnvs[0].length;
					indCNVCalls[1][j] = cnvs[1].length;
					totalCallsAnalyzed += cnvs[0].length + cnvs[1].length;

					if (cnvs[0].length > 0 && cnvs[1].length > 0) {
						numPassingAndPresent++;
						double callPercent = ((double) ((2 * counts[3]) + counts[4]) / (cnvs[0].length + cnvs[1].length));
						goodCallPPercent += callPercent;
						indPercentMatches[j] = callPercent;
					} else {
						indPercentMatches[j] = Double.NaN;
					}
				}
				averageGoodCallPercent = goodCallPPercent / numPassingAndPresent;
			} catch (Exception e) {
				System.err.println("Error comparing " + files[allPossibleCombinations[i][0]] + " and " + files[allPossibleCombinations[i][1]]);
				e.printStackTrace();
			}
		}
		comparedCalls = new CompareCalls_dev(numPassingAndPresent, averageGoodCallPercent, totalCallsAnalyzed, files, exactMatches, sigOlapMatches, comparedIDs, indCNVCalls, indPercentMatches);
		return comparedCalls;
	}

	public static String[] filterCNVs(String rootDirectory, String cnvFile, String LRR_lookup, String[] compFiles, double lrrFilter, double confFilter, int numMarkers) {
		Hashtable<String, Double> lrrs = new Hashtable<String, Double>();
		Hashtable<String, Integer> defineCompHash = new Hashtable<String, Integer>();
		BufferedReader reader;
		String[] line;
		PrintWriter[] writers = new PrintWriter[compFiles.length];
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
				for (int i = 0; i < compFiles.length; i++) {
					writers[i].close();
					int totalInCompFile = includeCounts[i] + filteredCounts[i];
					System.out.println("Info - Out of a total of " + totalCounts + " cnv calls in " + rootDirectory + cnvFile + " , a total of " + totalInCompFile + " matched ids listed in " + compFiles[i] + ", " + filteredCounts[i] + " were filtered out, " + includeCounts[i] + " were included for comparison");
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + rootDirectory + cnvFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + rootDirectory + cnvFile + "\"");
			System.exit(2);
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
		// TODO
		// redo for when ids are the same
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
		double inc = getIncrement(startVal, stopVal, numBins);
		for (int i = 0; i < numBins + 1; i++) {
			// System.out.println(inc * i + startVal);
			values[i] = (inc * i) + startVal;
		}
		return values;
	}

	private static double getIncrement(double startVal, double stopVal, int numBins) {
		return (stopVal - startVal) / numBins;
	}

	public static void iterate(double lrrstart, double lrrstop, int numLrrBins, double confstart, double confstop, int numConfBins, int probestart, int probeStop, String rootDirectory, String cnvFile, String LRR_lookup, String[] compFiles) {
		PrintWriter probeWriter, confWriter, lrrWriter, lrrProbeWriter;
		CompareCalls_dev comparedCalls;
		double[] lrrValues = new double[numLrrBins];
		double[] confValues = new double[numConfBins];
		long time;

		time = new Date().getTime();
		lrrValues = binIt(lrrstart, lrrstop, numLrrBins);
		confValues = binIt(confstart, confstop, numConfBins);
		String commonHeader = "_Value\tCallsFile1\tCallsFile2\tExactMatch\tSigOverlap\tindividualConcordance_Percent\tAverageConcordance_Percent\tnumindividual\ttotalCallsPassingFilter\tFile1\tFile2";

		try {
			probeWriter = new PrintWriter(new FileWriter(rootDirectory + "probeConcordance_" + probestart + "_" + probeStop + "LRR" + DEFAULT_LRR_FILTER + "Conf" + DEFAULT_CONF_FILTER + ".concord"));
			confWriter = new PrintWriter(new FileWriter(rootDirectory + "confConcordance_" + confstart + "_" + confstop + "LRR" + DEFAULT_LRR_FILTER + "numMarkers" + DEFAULT_NUM_MARKERS_FILTER + ".concord"));
			lrrWriter = new PrintWriter(new FileWriter(rootDirectory + "lrrConcordance_" + lrrstart + "_" + lrrstop + "numMarkers" + DEFAULT_NUM_MARKERS_FILTER + "Conf" + DEFAULT_CONF_FILTER + ".concord"));
			probeWriter.println("ID\tnumMarker" + commonHeader);
			confWriter.println("ID\tconf" + commonHeader);
			lrrWriter.println("ID\tlrr_SD" + commonHeader);
			// lrrProbeWriter = new PrintWriter(new FileWriter(rootDirectory + "lrrProbetestConcordance.concord"));
			// lrrProbeWriter.println("LRR_Value\tnumMarker" + commonHeader);
			for (int i = probestart; i < probeStop + 1; i++) {
				comparedCalls = compare(rootDirectory, filterCNVs(rootDirectory, cnvFile, LRR_lookup, compFiles, DEFAULT_LRR_FILTER, DEFAULT_CONF_FILTER, i));
				for (int j = 0; j < comparedCalls.getInds().length; j++) {
					if (comparedCalls.getIndCNVCalls()[0][j] > 0 && comparedCalls.getIndCNVCalls()[1][j] > 0) {
						probeWriter.println(comparedCalls.getInds()[j] + "\t" + i + "\t" + comparedCalls.getIndCNVCalls()[0][j] + "\t" + comparedCalls.getIndCNVCalls()[1][j] + "\t" + comparedCalls.getExactMatches()[j] + "\t" + comparedCalls.getSigOlapMatches()[j] + "\t" + comparedCalls.getindPercentMatches()[j] + "\t" + comparedCalls.getAveragePercent() + "\t" + comparedCalls.getNumberComputed() + "\t" + comparedCalls.getTotalCallsAnalyzed() + "\t" + Array.toStr(comparedCalls.getCnvFiles()) + "\t");
					}
				}
				// probeWriter.println(i + "\t" + comparedCalls.getAveragePercent() + "\t" + comparedCalls.getNumberComputed() + "\t" + comparedCalls.getTotalCallsAnalyzed() + "\t" + Array.toStr(comparedCalls.getCnvFiles()) + "\t" + Array.toStr(comparedCalls.getExactMatches()) + "\t" + Array.toStr(comparedCalls.getSigOlapMatches()));
			}
			for (int i = 0; i < lrrValues.length; i++) {
				comparedCalls = compare(rootDirectory, filterCNVs(rootDirectory, cnvFile, LRR_lookup, compFiles, lrrValues[i], DEFAULT_CONF_FILTER, DEFAULT_NUM_MARKERS_FILTER));
				for (int j = 0; j < comparedCalls.getInds().length; j++) {
					if (comparedCalls.getIndCNVCalls()[0][j] > 0 && comparedCalls.getIndCNVCalls()[1][j] > 0) {
						lrrWriter.println(comparedCalls.getInds()[j] + "\t" + lrrValues[i] + "\t" + comparedCalls.getIndCNVCalls()[0][j] + "\t" + comparedCalls.getIndCNVCalls()[1][j] + "\t" + comparedCalls.getExactMatches()[j] + "\t" + comparedCalls.getSigOlapMatches()[j] + "\t" + comparedCalls.getindPercentMatches()[j] + "\t" + comparedCalls.getAveragePercent() + "\t" + comparedCalls.getNumberComputed() + "\t" + comparedCalls.getTotalCallsAnalyzed() + "\t" + Array.toStr(comparedCalls.getCnvFiles()) + "\t");
					}
				}
				// lrrWriter.println(lrrValues[i] + "\t" + comparedCalls.getAveragePercent() + "\t" + comparedCalls.getNumberComputed() + "\t" + comparedCalls.getTotalCallsAnalyzed() + "\t" + Array.toStr(comparedCalls.getCnvFiles()));
			}
			for (int i = 0; i < confValues.length; i++) {
				comparedCalls = compare(rootDirectory, filterCNVs(rootDirectory, cnvFile, LRR_lookup, compFiles, DEFAULT_LRR_FILTER, confValues[i], DEFAULT_NUM_MARKERS_FILTER));
				for (int j = 0; j < comparedCalls.getInds().length; j++) {
					if (comparedCalls.getIndCNVCalls()[0][j] > 0 && comparedCalls.getIndCNVCalls()[1][j] > 0) {
						confWriter.println(comparedCalls.getInds()[j] + "\t" + confValues[i] + "\t" + comparedCalls.getIndCNVCalls()[0][j] + "\t" + comparedCalls.getIndCNVCalls()[1][j] + "\t" + comparedCalls.getExactMatches()[j] + "\t" + comparedCalls.getSigOlapMatches()[j] + "\t" + comparedCalls.getindPercentMatches()[j] + "\t" + comparedCalls.getAveragePercent() + "\t" + comparedCalls.getNumberComputed() + "\t" + comparedCalls.getTotalCallsAnalyzed() + "\t" + Array.toStr(comparedCalls.getCnvFiles()) + "\t");
					}
				}

				// confWriter.println(confValues[i] + "\t" + comparedCalls.getAveragePercent() + "\t" + comparedCalls.getNumberComputed() + "\t" + comparedCalls.getTotalCallsAnalyzed() + "\t" + Array.toStr(comparedCalls.getCnvFiles()));
			}
			// for (int j = 0; j < lrrValues.length; j++) {
			// for (int i = probestart; i < probeStop + 1; i++) {
			// comparedCalls = compare(rootDirectory, filterCNVs(rootDirectory, cnvFile, LRR_lookup, compFiles, lrrValues[j], DEFAULT_CONF_FILTER, i));
			// lrrProbeWriter.println(lrrValues[j] + "\t" + i + "\t" + comparedCalls.getAveragePercent() + "\t" + comparedCalls.getNumberComputed() + "\t" + comparedCalls.getTotalCallsAnalyzed() + "\t" + Array.toStr(comparedCalls.getCnvFiles()));
			// }
			// }

			probeWriter.close();
			confWriter.close();
			lrrWriter.close();
			// lrrProbeWriter.close();
			System.out.println("Iterations took " + ext.getTimeElapsed(time));
		} catch (Exception e) {
			System.err.println("Error comparing ");
			e.printStackTrace();
		}

	}

	public static Hashtable<Byte, List<Integer>> getMarkerLookup(String rootDirectory, String markerPositionsFile) {
		Hashtable<String, Integer> tracker = new Hashtable<String, Integer>();
		Hashtable<Byte, List<Integer>> markers = new Hashtable<Byte, List<Integer>>();
		BufferedReader reader;
		byte chr;
		String[] line;
		try {
			reader = Files.getAppropriateReader(rootDirectory + markerPositionsFile);
			do {
				line = reader.readLine().trim().split("\t", -1);
			} while (reader.ready() && !line[0].startsWith("Name"));
			// System.out.println(line[0]);
			if (!reader.ready()) {
				System.err.println("DID not FIND THE INDEX FACTORS THAT NEED TO BE ADDED");
			}
			while (reader.ready()) {
				line = reader.readLine().trim().split("\t", -1);
				if (line[1].startsWith("X")) {
					chr = (byte) 23;
				} else if (line[1].equals("Y")) {
					chr = (byte) 24;
				} else if (line[1].equals("MT") || line[1].equals("0") || line[1].equals("0")) {
					continue;
				} else {
					chr = Byte.valueOf(line[1]);
					if (!markers.containsKey(chr)) {
						List<Integer> bpList = new ArrayList<Integer>();
						markers.put(chr, bpList);
					}
					markers.get(chr).add(Integer.parseInt(line[2]));
				}
			}
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + markerPositionsFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + rootDirectory + markerPositionsFile + "\"");
			System.exit(2);
		}
		sortBPbyChromosome(markers);
		return markers;
	}

	private static void sortBPbyChromosome(Hashtable<Byte, List<Integer>> markers) {
		for (int i = 1; i <= markers.keySet().size(); i++) {
			if (markers.containsKey((byte) i)) {
				Collections.sort(markers.get((byte) i));
			} else {
				System.err.println("chromosome " + i + " is undefined");
			}
		}
	}

	public static void cleanCNVs(String rootDir, String cnvfile, String markerFile, double bpFraction, double markerFraction) {
		PrintWriter writer;
		String[] inds;
		CNVariant[] cnvs;
		Hashtable<String, Hashtable<String, Vector<CNVariant>>> hashCNVs = new Hashtable<String, Hashtable<String, Vector<CNVariant>>>();
		hashCNVs = getCNVs(cnvfile);
		inds = HashVec.getKeys(hashCNVs);
		writer = Files.getAppropriateWriter(cnvfile + ".clean.cnv");
		writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
		Hashtable<Byte, List<Integer>> markers = getMarkerLookup(rootDir, markerFile);
		int changecount = 0;
		int nochange = 0;
		int totalCount = 0;

		for (int i = 0; i < inds.length; i++) {
			cnvs = new CNVariant[hashCNVs.get(inds[i]).get(0 + "").size()];
			cnvs = CNVariant.sort(CNVariant.toCNVariantArray(hashCNVs.get(inds[i]).get(0 + "")));
			CNVariant[] cleanedCNVs = new CNVariant[cnvs.length];

			for (int j = 0; j < cnvs.length; j++) {
				totalCount++;
				if (j + 1 < cnvs.length) {
					CNVariant thisCNV = cnvs[j];
					CNVariant nextCNV = cnvs[j + 1];
					if (checkChr(thisCNV, nextCNV) && bpGapCheck(thisCNV, nextCNV, bpFraction) && checkCN(thisCNV, nextCNV) && markerGapCheck(thisCNV, nextCNV, markers, markerFraction)) {
						System.out.println(cnvs[j + 1].getIndividualID() + "\t" + cnvs[j + 1].getChr() + "\t" + cnvs[j + 1].getStart());
						changecount++;
					} else {
						nochange++;
						writer.println(thisCNV.toString());
					}
				}
			}
		}
		System.out.println(changecount + "\t" + nochange + "\t" + totalCount);
		writer.close();
	}

	private static boolean markerGapCheck(CNVariant thisCNV, CNVariant nextCNV, Hashtable<Byte, List<Integer>> markers, double markerFraction) {
		int gapLengthMarkers = getStartIndex(nextCNV, markers) - getStopIndex(thisCNV, markers) - 1;
		int proposedNewMarkerLength = getStopIndex(nextCNV, markers) - getStartIndex(thisCNV, markers) + 1;
		return testFrac(markerFraction, ((double) gapLengthMarkers / proposedNewMarkerLength));
	}

	private static int getStartIndex(CNVariant CNV, Hashtable<Byte, List<Integer>> markers) {
		return markers.get(CNV.getChr()).indexOf(CNV.getStart());
	}

	private static int getStopIndex(CNVariant CNV, Hashtable<Byte, List<Integer>> markers) {
		return markers.get(CNV.getChr()).indexOf(CNV.getStop());
	}

	private static boolean bpGapCheck(CNVariant thisCNV, CNVariant nextCNV, double bpFraction) {
		int gapLengthBP = nextCNV.getStart() - thisCNV.getStop() - 1;
		int proposedNewBPLength = nextCNV.getStop() - thisCNV.getStart() + 1;
		return testFrac(bpFraction, ((double) (gapLengthBP / proposedNewBPLength)));
	}

	private static boolean testFrac(double Fraction, double FractionToTest) {
		if (FractionToTest <= Fraction) {
			return true;
		} else {
			return false;
		}
	}

	private static boolean checkChr(CNVariant thisCNV, CNVariant nextCNV) {
		return thisCNV.getChr() == nextCNV.getChr();
	}

	private static boolean checkCN(CNVariant thisCNV, CNVariant nextCNV) {
		return thisCNV.getCN() == nextCNV.getCN();
	}

	private static Hashtable<String, Hashtable<String, Vector<CNVariant>>> getCNVs(String file) {
		Vector<CNVariant> v = new Vector<CNVariant>();
		BufferedReader reader;
		String[] line;
		Hashtable<String, Hashtable<String, Vector<CNVariant>>> hash = new Hashtable<String, Hashtable<String, Vector<CNVariant>>>();
		Hashtable<String, Vector<CNVariant>> source = new Hashtable<String, Vector<CNVariant>>();
		try {
			reader = new BufferedReader(new FileReader(file));
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
				if (source.containsKey(0 + "")) {
					v = source.get(0 + "");
				} else {
					source.put(0 + "", v = new Vector<CNVariant>());
				}
				v.add(new CNVariant(line, 0));
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + file + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + file + "\"");
			System.exit(2);
		}
		return hash;
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
		String markerFile = DEFAULT_MARKER_FILE;
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

			// getMarkerLookup(rootDirectory, markerFile);
			cleanCNVs(rootDirectory, filterCNVs(rootDirectory, cnvFile, LRR_lookup, compFiles, lrrFilter, confFilter, numMarkers)[0], markerFile, 0.2, 0.2);
			// compare(rootDirectory, files);
			// defineCompLists(rootDirectory, compLists);
			// getLRRs(rootDirectory, LRR_lookup);
			// iterate(0.2, 0.6, 100, 0.0, 100, 100, 0, 100, rootDirectory, cnvFile, LRR_lookup, compFiles);
			// compare(rootDirectory, filterCNVs(rootDirectory, cnvFile, LRR_lookup, compFiles, lrrFilter, confFilter, numMarkers));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

// if (testFrac(markerFraction, ((double) gapLengthMarkers / proposedNewMarkerLength))) {
// System.out.println(nextCNV.getIndividualID() + "\t" + thisCNV.getIndividualID() + "\t" + markers.get(nextCNV.getChr()).indexOf(nextCNV.getStart()) + "\t" + (markers.get(thisCNV.getChr()).indexOf(thisCNV.getStop()) + 1));
// System.out.println(gapLengthMarkers + "\t" + proposedNewMarkerLength + "\t" + markerFractionToTest + "\t" + testFrac(markerFraction, markerFractionToTest));
//
// }

// System.out.println(nextCNV.getIndividualID() + "\t" + thisCNV.getIndividualID() + "\t" + markers.get(nextCNV.getChr()).indexOf(nextCNV.getStart()) + "\t" + (markers.get(thisCNV.getChr()).indexOf(thisCNV.getStop()) + 1));
// System.out.println(nextCNV.getChr() + "\t" + nextCNV.getStart() + "\t" + thisCNV.getChr() + "\t" + thisCNV.getStart());

