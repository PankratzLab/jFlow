package org.genvisis.cnv.qc;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.genvisis.cnv.analysis.ProjectCNVFiltering;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CNVFilter;
import org.genvisis.common.CNVFilter.CNVFilterPass;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.CNVariantHash;

/**
 * Measure statistics for agreement (concordance) of CNV calls. There are two methods for
 * measurement:
 * <p>
 * <ul>
 * <li>Blind duplicates - the <code>pairingFile</code> should map separate samples that are expected
 * to be genetically identical (e.g. duplicate samples in the dataset), and only a single
 * <code>cnvFile</code> is input. The output concordance indicates how well the given CNV algorithm
 * does when provided similar data.</li>
 * <li>Gold standard - the <code>pairingFile</code> should map the same sample which has been called
 * in two separate CNV methods. The base <code>cnvFile</code> is taken as the test case and compared
 * against a <code>cnvControl<code> file. Output concordance indicates how well a modified CNV
 * algorithm does against a known standard.</li>
 * </ul>
 * </p>
 */
public class CNVConcordance {
	private static class ComparisionGlobalResults {
		private final CNVConcordance.ComparisionIndividualResults[] indResults;
		private int totalCNVs;
		private int totalOverlap;
		private int totalSigOverlap;
		private int totalPerfectOverlap;
		private int totalUnpaired;
		private double avgIndOverlap;
		private double avgIndSigOverlap;
		private double avgIndPerfectOverlap;
		private double avgIndOverlapScore;
		private double avgIndUnpaired;
		private double avgIndRelSigOverlap;
		private double avgIndRelPerfectOverlap;
		private double globalOverlap;
		private double globalSigOverlap;
		private double globalPerfectOverlap;
		private double avgGlobalOverlapScore;
		private final ArrayList<Double> globalOverLapScores;

		public ComparisionGlobalResults(CNVConcordance.ComparisionIndividualResults[] comparisionIndividualResults) {
			indResults = comparisionIndividualResults;
			globalOverLapScores = new ArrayList<Double>();
		}

		public String getReport() {
			StringBuilder sb = new StringBuilder();
			sb.append(totalCNVs).append("\t");
			sb.append(totalOverlap).append("\t");
			sb.append(totalSigOverlap).append("\t");
			sb.append(totalPerfectOverlap).append("\t");
			sb.append(totalUnpaired).append("\t");
			sb.append(avgIndOverlap).append("\t");
			sb.append(avgIndSigOverlap).append("\t");
			sb.append(avgIndPerfectOverlap).append("\t");
			sb.append(avgIndOverlapScore).append("\t");
			sb.append(avgIndUnpaired).append("\t");
			sb.append(avgIndRelSigOverlap).append("\t");
			sb.append(avgIndRelPerfectOverlap).append("\t");
			sb.append(globalOverlap).append("\t");
			sb.append(globalSigOverlap).append("\t");
			sb.append(globalPerfectOverlap).append("\t");
			sb.append(avgGlobalOverlapScore).append("\t");

			return sb.toString();
		}

		public void populateMetrics() {
			for (ComparisionIndividualResults indResult : indResults) {
				totalCNVs += indResult.getTotalCNVCount();
				totalOverlap += indResult.getNumOverlap();
				totalSigOverlap += indResult.getNumSigOverlap();
				totalPerfectOverlap += indResult.getNumPerfectOverlap();
				avgIndOverlap += indResult.getPercentOverlap();
				avgIndSigOverlap += indResult.getPercentSigOverlap();
				avgIndPerfectOverlap += indResult.getPercentPerfectOverlap();
				avgIndOverlapScore += indResult.getAvgOverlapScore();
				totalUnpaired += indResult.getUnpaired();
				avgIndRelSigOverlap += indResult.getRelativePercentSigOverlap();
				avgIndRelPerfectOverlap += indResult.getRelativePercentPerfectOverlap();
				indResult.addAllScores(globalOverLapScores);
			}
			avgIndOverlap /= indResults.length;
			avgIndSigOverlap /= indResults.length;
			avgIndPerfectOverlap /= indResults.length;
			avgIndOverlapScore /= indResults.length;
			avgIndUnpaired = (double) totalUnpaired / indResults.length;
			avgIndRelSigOverlap /= indResults.length;
			avgIndRelPerfectOverlap /= indResults.length;
			globalOverlap = (double) totalOverlap / totalCNVs;
			globalSigOverlap = (double) totalSigOverlap / totalCNVs;
			globalPerfectOverlap = (double) totalPerfectOverlap / totalCNVs;
			avgGlobalOverlapScore = ArrayUtils.mean(globalOverLapScores);
		}
	}
	private static class ComparisionIndividualResults {
		private String ind1;
		private String ind2;
		private final int numCNVind1;
		private final int numCNVind2;
		private int numOverlap;
		private int numSigOverlap;
		private int numPerfectOverlap;
		private final List<Double> overlapScores;

		public ComparisionIndividualResults(String ind1, String ind2, int numCNVind1, int numCNVind2) {
			this.ind1 = ind1;
			this.ind2 = ind2;
			this.numCNVind1 = numCNVind1;
			this.numCNVind2 = numCNVind2;
			overlapScores = new ArrayList<Double>(getTotalCNVCount());
		}

		public void addAllScores(List<Double> globalOverLapScores) {
			globalOverLapScores.addAll(overlapScores);
		}

		public void addOverlap() {
			numOverlap += 1;
		}

		public void addOverlapScore(double overlapScore) {
			if (!Double.isNaN(overlapScore)) {
				overlapScores.add(overlapScore);
			}
		}

		public void addPerfect() {
			numPerfectOverlap += 1;
		}

		public void addSigOverlap() {
			numSigOverlap += 1;
		}

		public double getAvgOverlapScore() {
			if (!overlapScores.isEmpty()) {
				return ArrayUtils.mean(overlapScores);
			}
			return 0.0;
		}

		public int getNumOverlap() {
			return numOverlap;
		}

		public int getNumPerfectOverlap() {
			return numPerfectOverlap;
		}

		public int getNumSigOverlap() {
			return numSigOverlap;
		}

		public double getPercentOverlap() {
			if (getTotalCNVCount() == 0) {
				return 0.0;
			}
			return (double) numOverlap / getTotalCNVCount();
		}

		public double getPercentPerfectOverlap() {
			if (getTotalCNVCount() == 0) {
				return 0.0;
			}
			return (double) numPerfectOverlap / getTotalCNVCount();
		}

		public double getPercentSigOverlap() {
			if (getTotalCNVCount() == 0) {
				return 0.0;
			}
			return (double) numSigOverlap / getTotalCNVCount();
		}

		public int getTotalCNVCount() {
			return numCNVind1 + numCNVind2;
		}

		public int getUnpaired() {
			return getTotalCNVCount() - getNumOverlap();
		}

		public double getRelativePercentPerfectOverlap() {
			return getRelative(numPerfectOverlap);
		}

		public double getRelativePercentSigOverlap() {
			return getRelative(numSigOverlap);
		}

		private double getRelative(int target) {
			if (getTotalCNVCount() == 0 || getNumOverlap() == 0) {
				return 0.0;
			}
			return (double) target / getNumOverlap();

		}

		public String getInfo() {
			// TODO for printing individual information
			return "";
		}
	}

	// TODO would prefer to integrate this with the global stats in a way that avoids parallel arrays
	public static final String[] REPORT_HEADER = {"Total CNVs Compared", "Total Overlapping CNVs",
																								"Total Significantly Overlapping CNVs",
																								"Total Perfectly overlapping cnvs",
																								"Total Unpaired CNVs", "Average Individual Overlap",
																								"Average Individual Significant Overlap",
																								"Average Individual Perfect Overlap",
																								"Average Individual Overlap Score",
																								"Average Individual Unpaired",
																								"Average Relative Significant Overlap",
																								"Average Relative Perfect Overlap",
																								"Global Overlap", "Global Significant Overlap",
																								"Global Perfect Overlap",
																								"Average Global Overlap Score"};
	public static final String COMMAND_PROJECT = "proj=";
	public static final String COMMAND_CNV = "cnvFile=";
	public static final String COMMAND_PAIR = "pairFile=";
	public static final String COMMAND_LRRSD_THRESHOLD = "lrrSD=";
	public static final String COMMAND_OUTPUT = "output=";
	public static final String COMMAND_DEFUALTS = "-defaults";
	public static final String COMMAND_DIR = "dir=";
	public static final String COMMAND_NUMCNVS = "numCNVS=";
	public static final String COMMAND_CNV_CONCORDANCE = "cnvConcordance";
	public static final String COMMAND_CNV_CONCORDANCE_DESCRIPTION = "- compute the concordance between replicates in a cnv file";
	private final Project proj;
	private final String[][] pairs;
	private final Map<String, double[]> qcMap;
	private String report;
	private final CNVariantHash testHash;
	private final CNVariantHash controlHash;
	private ComparisionIndividualResults[] comparisionResults;
	private boolean fail;
	private final CNVFilter filter;
	private final int numCNVs;
	private final double maxLrr;
	private final double minCallRate;
	private SampleData sampleData;

	public CNVConcordance(Project proj, String[][] pairs, Map<String, double[]> qcMap, double maxLrr,
												double minCallRate, CNVariantHash cNVariantHash, CNVFilter filter,
												int numCNVs) {
		this(proj, pairs, qcMap, maxLrr, minCallRate, cNVariantHash, cNVariantHash, filter, numCNVs);
	}

	// TODO need to put some thought in how markers affect scoring
	public CNVConcordance(Project proj, String[][] pairs, Map<String, double[]> qcMap, double maxLrr,
												double minCallRate, CNVariantHash testHash, CNVariantHash controlHash,
												CNVFilter filter, int numCNVs) {
		this.testHash = testHash;
		this.controlHash = controlHash;
		this.proj = proj;
		sampleData = proj.getSampleData(SampleData.MINIMAL_SAMPLE_DATA_HEADER.length, false);
		this.qcMap = qcMap;
		this.maxLrr = maxLrr;
		this.minCallRate = minCallRate;
		this.pairs = pairs;
		this.filter = filter;
		this.numCNVs = numCNVs;
	}

	public void determineConcordance() {
		comparisionResults = compareAll();
		if (!fail) {
			ComparisionGlobalResults cgr = new ComparisionGlobalResults(comparisionResults);
			cgr.populateMetrics();
			report = cgr.getReport();
		} else {
			proj.getLog().reportError("Error - comparison has failed");
		}
	}

	public String getReport() {
		return report;
	}

	private ComparisionIndividualResults[] compareAll() {
		if (!fail) {
			int numComp = 0;
			ArrayList<ComparisionIndividualResults> allResults = new ArrayList<ComparisionIndividualResults>();
			for (String[] pair : pairs) {
				String ind1 = pair[0];
				String ind2 = pair[1];
				String pairString = ind1 + ", " + ind2;
				if (checkQc(ind1, ind2)) {
					List<CNVariant> ind1CNVs = findVariants(ind1, testHash);
					List<CNVariant> ind2CNVs = findVariants(ind2, controlHash);
					ComparisionIndividualResults results = compareInds(ind1, ind1CNVs, ind2, ind2CNVs);

					if (results.getTotalCNVCount() > 0) {
						allResults.add(results);
					} else {
						System.out.println("No common cnvs found for pair: " + pairString);
					}
				}
				numComp++;
			}
			System.out.println("Total of " + numComp + " comparisons");
			int excluded = numComp - allResults.size();
			System.out.println("Total comparisons excluded= " + excluded);
			return allResults.toArray(new ComparisionIndividualResults[allResults.size()]);
		}
		return null;
	}// This command is what the kids call "janky." currently won't work if the FIDs and IIDs don't

	private boolean checkQc(String... inds) {
		int passed = 1;

		// At least one sample in this comparison must pass qc. If the rest are not found that is OK
		// (assume gold standard cnv list). If any fail, this comparison is excluded.
		for (String ind : inds) {
			passed = Math.min(passed, checkInd(ind));
		}

		return passed == 0;
	}

	/**
	 * @return 1 if qc not found, 0 if qc passed, -1 if qc failed
	 */
	private int checkInd(String ind) {
		String id = ind;
		double[] qc = qcMap.get(id);
		if (qc == null) {
			String[] lookup = sampleData.lookup(ind);
			if (lookup != null) {
				id = lookup[1];
				qc = qcMap.get(id);
			}
		}
		if (qc == null) {
			// qc not found
			return 1;
		}

		// qc failed
		if (qc[0] > maxLrr) {
			System.out.println("Warning - " + ind + " has LRR > " + maxLrr + " (" + qc[0]
												 + "). Excluded from final analysis");
			return -1;
		} else if (qc[1] < minCallRate) {
			System.out.println("Warning - " + ind + " has callrate < " + minCallRate + " (" + qc[0]
												 + "). Excluded from final analysis");
			return -1;
		}

		// qc passed
		return 0;
	}

	private List<CNVariant> findVariants(String id, CNVariantHash hash) {
		String fidiid = id;
		Hashtable<String, CNVariant[]> data = hash.getDataFor(fidiid);
		if (data.isEmpty()) {
			fidiid = id + "\t" + id;
			data = hash.getDataFor(fidiid);
		}
		if (data.isEmpty()) {
			fidiid = sampleData.lookup(id)[1];
			data = hash.getDataFor(fidiid);
		}
		return extractVariants(data);
	}

	private ComparisionIndividualResults compareInds(String ind1, List<CNVariant> ind1CNVs,
																									 String ind2, List<CNVariant> ind2CNVs) {
		if (ind1CNVs.isEmpty() || ind2CNVs.isEmpty() || ind1CNVs.size() > numCNVs
				|| ind2CNVs.size() > numCNVs) {
			return new ComparisionIndividualResults(ind1, ind2, ind1CNVs.size(), ind2CNVs.size());
		}
		ind1CNVs = filterCNVs(filter, ind1CNVs);
		ind2CNVs = filterCNVs(filter, ind2CNVs);
		// changed to 1200 from 1000 to shorten print output
		int warn = 1200;
		if (ind1CNVs.size() > warn) {
			proj.getLog().report("Warning - " + ind1 + " has more than " + warn + " CNVs ("
													 + ind1CNVs.size() + "), this could lead to skewed concordance rates");
		}
		if (ind2CNVs.size() > warn) {
			proj.getLog().report("Warning - " + ind2 + " has more than " + warn + " CNVs ("
													 + ind2CNVs.size() + "), this could lead to skewed concordance rates");
		}
		if (ind1CNVs.isEmpty() || ind2CNVs.isEmpty()) {
			return new ComparisionIndividualResults(ind1, ind2, ind1CNVs.size(), ind2CNVs.size());
		}
		return compareIndCNVs(ind1CNVs, ind2CNVs);
	}

	public static void determineConcordance(Project proj, String cnvFile, String dir, String pairFile,
																					String qcFile, String cnvControl, CNVFilter filter,
																					int numCNVs, int CN, double lrrMax, double callMin,
																					String output) {
		if ((pairFile == null) || (pairFile.equals(""))) {
			proj.getLog()
					.reportError("Error - a file of sample pairings must be provided to determine concordance");
			return;
		}
		if ((cnvFile == null) || (cnvFile.equals("")) && dir == null) {
			proj.getLog().reportError("Error - a file of cnvs must be provided to determine concordance");
			return;
		}
		String[][] pairs = loadPairs(proj.PROJECT_DIRECTORY.getValue() + pairFile);
		Map<String, double[]> qcMap = loadQC(proj.PROJECT_DIRECTORY.getValue() + qcFile);
		filter.setCN(CN);
		if (dir != null) {
			String[] cnvFiles = Files.list(proj.PROJECT_DIRECTORY.getValue() + dir, ".cnv", false);
			cnvFiles = Files.toFullPaths(cnvFiles, proj.PROJECT_DIRECTORY.getValue() + dir);
			proj.getLog().report(ArrayUtils.toStr(cnvFiles));
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + dir
																														+ output));
				int start = filter.getMinNumMarkers();
				CNVariantHash[] cNVariantHash = new CNVariantHash[cnvFiles.length];
				for (int i = 0; i < cnvFiles.length; i++) {
					cNVariantHash[i] = CNVariantHash.load(cnvFiles[i], 1, false, proj.getLog());
					for (int j = 0; j < REPORT_HEADER.length; j++) {
						writer.print(((i == 0) && (j == 0) ? "" : "\t") + REPORT_HEADER[j] + "."
												 + ext.rootOf(cnvFiles[i]));
					}
					writer.print("\tnumberOFProbes." + ext.rootOf(cnvFiles[i]));
				}
				writer.println();
				for (int j = start; j <= 100; j++) {
					filter.setMinNumMarkers(j);
					for (int i = 0; i < cnvFiles.length; i++) {
						long time = System.currentTimeMillis();
						proj.getLog()
								.report(ext.getTime() + " Info - beginning comparision for " + cnvFiles[i]);

						CNVConcordance cnvConcordance = new CNVConcordance(proj, pairs, qcMap, lrrMax, callMin,
																															 cNVariantHash[i], filter, numCNVs);
						cnvConcordance.determineConcordance();
						writer.print((i == 0 ? "" : "\t") + cnvConcordance.getReport() + "\t" + j);
						proj.getLog().report(ext.getTime() + " Info - finished comparision for " + cnvFiles[i]
																 + " and took " + ext.getTimeElapsed(time));
					}
					writer.println();
				}
				writer.close();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + proj.PROJECT_DIRECTORY.getValue() + output);
				proj.getLog().reportException(e);
			}
		} else {
			CNVariantHash cNVariantHash = CNVariantHash.load(proj.PROJECT_DIRECTORY.getValue() + cnvFile,
																											 1, false, proj.getLog());
			CNVConcordance cnvConcordance;

			if (cnvControl == null) {
				cnvConcordance = new CNVConcordance(proj, pairs, qcMap, lrrMax, callMin, cNVariantHash,
																						filter, numCNVs);
			} else {
				CNVariantHash controlHash = CNVariantHash.load(proj.PROJECT_DIRECTORY.getValue()
																											 + cnvControl, 1, false, proj.getLog());
				cnvConcordance = new CNVConcordance(proj, pairs, qcMap, lrrMax, callMin, cNVariantHash,
																						controlHash, filter, numCNVs);
			}

			cnvConcordance.determineConcordance();
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()
																														+ output));
				writer.println(ArrayUtils.toStr(REPORT_HEADER));
				writer.println(cnvConcordance.getReport());

				writer.close();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + proj.PROJECT_DIRECTORY.getValue() + output);
				proj.getLog().reportException(e);
			}
			proj.getLog().report(ArrayUtils.toStr(REPORT_HEADER));
			proj.getLog().report(cnvConcordance.getReport());
		}
	}

	public static List<CNVariant> filterCNVs(CNVFilter cnvFilter, List<CNVariant> cnvs) {
		List<CNVariant> filtered = new ArrayList<CNVariant>(cnvs.size());
		for (CNVariant cnv : cnvs) {
			CNVFilterPass filterPass = cnvFilter.getCNVFilterPass(cnv);
			if (filterPass.passedFilter()) {
				if ((filterPass.isCentromeric()) && (cnvFilter.isBreakupCentromeres())) {
					filtered.add(cnvFilter.breakUpCentromere(filterPass, cnv)[0]);
					filtered.add(cnvFilter.breakUpCentromere(filterPass, cnv)[1]);
				} else {
					filtered.add(cnv);
				}
			} else if (filterPass.isIndIsExcluded()) {
				return Collections.<CNVariant>emptyList();
			}
		}
		return filtered;
	}

	public static void fromParameters(String filename, Logger log) {
		Vector<String> params = Files.parseControlFile(filename, "cnvConcordance", getParserParams(),
																									 log);
		if (params != null) {
			main(ArrayUtils.toStringArray(params));
		}
	}

	public static String[] getParserParams() {
		String[] params = new String[15];
		params[0] = "#To intialize the cnv concorder, provide the following arguments";
		params[1] = "#the full path to a project properties file";
		params[2] = "proj=";
		params[3] = "#a path (relative to the project directory) to a cnv file ";
		params[4] = "cnvFile=";
		params[5] = "#a path (relative to the project directory) to a file of sample pairs ";
		params[6] = "pairFile=";
		params[7] = "#a threshold for log R ratio Standard deviation, but this is not implemented yet";
		params[8] = "#lrrSD=";
		params[9] = "#a path (relative to the project directory) to a file for ouput";
		params[10] = "#output=";
		params[11] = "#a path (relative to the project directory) containing multiple \".cnv\" files, all will be analyzed and the cnvFile= command will be overridden";
		params[12] = "#dir=";

		params[13] = "#maximum number of cnvs per individual";
		params[14] = "#numCNVS=";

		System.out.println(ArrayUtils.toStr(CNVFilter.getDefaultCNVParams()));
		params = ArrayUtils.concatAll(params, new String[][] {CNVFilter.getDefaultCNVParams()});

		return params;
	}

	private static ComparisionIndividualResults compareIndCNVs(List<CNVariant> ind1CNVs,
																														 List<CNVariant> ind2CNVs) {
		ComparisionIndividualResults currentComparison = new ComparisionIndividualResults(ind1CNVs.get(0)
																																															.getFamilyID()
																																											+ "\t"
																																											+ ind1CNVs.get(0)
																																																.getIndividualID(),
																																											ind2CNVs.get(0)
																																															.getFamilyID() + "\t" + ind2CNVs.get(0)
																																																															.getIndividualID(),
																																											ind1CNVs.size(),
																																											ind2CNVs.size());
		compareIndCNVHelper(currentComparison, ind1CNVs, ind2CNVs);
		compareIndCNVHelper(currentComparison, ind2CNVs, ind1CNVs);

		return currentComparison;
	}

	private static void compareIndCNVHelper(ComparisionIndividualResults currentComparison,
																					List<CNVariant> cnvs1, List<CNVariant> cnvs2) {
		for (CNVariant ind1cnv : cnvs1) {
			for (CNVariant ind2cnv : cnvs2) {
				if (ind1cnv.getCN() == ind2cnv.getCN()) {
					if (ind1cnv.overlaps(ind2cnv)) {
						currentComparison.addOverlapScore(ind1cnv.overlapScore(ind2cnv));
						currentComparison.addOverlap();
						if (ind1cnv.significantOneWayOverlap(ind2cnv)) {
							currentComparison.addSigOverlap();
						}
						if (ind1cnv.equals(ind2cnv)) {
							currentComparison.addPerfect();
						}
					}
				}
			}
		}
	}

	private static List<CNVariant> extractVariants(Hashtable<String, CNVariant[]> indCNVS) {
		List<CNVariant> cnvs = new ArrayList<CNVariant>();
		List<String> arr = Collections.list(indCNVS.keys());
		for (int i = 0; i < arr.size(); i++) {
			CNVariant[] chrCNVs = indCNVS.get(arr.get(i));
			for (CNVariant chrCNV : chrCNVs) {
				cnvs.add(chrCNV);
			}
		}
		return cnvs;
	}

	private static String[][] loadPairs(String pairFile) {
		String[] load = HashVec.loadFileToStringArray(pairFile, false, new int[] {0, 1}, false);
		String[][] pairs = new String[load.length][];
		for (int i = 0; i < load.length; i++) {
			String[] tmp = load[i].split("\t");
			ArrayList<String> compareDef = new ArrayList<String>();
			for (int j = 0; j < tmp.length; j++) {
				if (!tmp[j].equals("")) {
					compareDef.add(tmp[j]);
				}
			}
			pairs[i] = compareDef.toArray(new String[compareDef.size()]);
		}
		return pairs;
	}

	private static Map<String, double[]> loadQC(String qcFile) {
		// Known Issue - genvisis-generated LRR file might have headers sprinkled throughout
		String[] load = HashVec.loadFileToStringArray(qcFile, false, null, false);
		Map<String, double[]> qcMap = new HashMap<String, double[]>();
		for (int i = 0; i < load.length; i++) {
			String[] line = load[i].split("\t");
			// skip header lines
			if (line[2].contains("LRR")) {
				continue;
			}
			double[] qc = new double[2];
			qcMap.put(line[0], qc);
			qc[0] = Double.valueOf(line[2]);
			qc[1] = Double.valueOf(line[7]);
		}
		return qcMap;
	}

	public static void main(String[] args) {
		long startTime = System.currentTimeMillis();
		int numArgs = args.length;
		String filename = null;
		String logfile = null;
		String cnvFile = null;
		String cnvControl = null;
		String pairFile = null;
		String dir = null;
		String qcFile = null;
		int numCNVs = 2147483647;
		boolean defaults = false;
		String output = "cnv.concordance.txt";
		int CN = -1;
		double lrr = 0.32;
		double call = 0.97;

		String usage = "\njlDev.CNVConcordance requires 0-1 arguments\n";
		usage = usage + "   (1) project filename  (i.e. proj=" + filename + " (no default))\n";
		usage = usage + "   (2) cnvFile  (i.e.cnvFile=" + filename + " (no default))\n";
		usage = usage + "   (3) sample pair file (i.e. pairFile=" + filename + " (no default))\n";
		usage = usage + "   (4) quality control file (i.e. qcFile=" + filename + " (no default))\n";
		usage = usage + "   OPTIONAL:";
		usage = usage + "   (5) control cnv file (i.e. cnvControl=" + filename + " (no default))\n";
		usage = usage + "   (6) output file name  (i.e.output=" + output + " (default))\n";
		usage = usage + "   (7) log file  (i.e. log=" + filename + " (no default))\n";
		usage = usage
						+ "\t (8) For cnv filtering, use the default values (i.e. -default ( not the default))\n";
		usage = usage + "\t (9) a directory containing multiple cnv files (i.e. dir= ( no default))\n";
		usage = usage + "\t (10) maximum number of cnvs (i.e. numCNVS=" + numCNVs + " (default))\n";
		usage = usage + "\t (11) lrr cutoff (i.e. lrr=" + lrr + " (default))\n";
		usage = usage + "\t (12) minimum callrate (i.e. call=" + call + " (default))\n";

		usage = usage + "\t (13) further usage:\n" + ArrayUtils.toStr(CNVFilter.getDefaultCNVParams());
		Project proj;
		if (ext.indexOfStr("proj=", args, true, false) >= 0) {
			proj = new Project(ext.parseStringArg(args[ext.indexOfStr("proj=", args, true, false)], ""),
												 logfile, false);
		} else {
			proj = new Project(filename, logfile, false);
		}
		CNVFilter filter = ProjectCNVFiltering.setupCNVFilterFromArgs(proj, args, null, defaults,
																																	proj.getLog());
		if (ext.indexOfStr("-defaults", args) >= 0) {
			ProjectCNVFiltering.setCNVDefaults(filter, proj);
		}
		for (String arg : args) {
			if ((arg.equals("-h")) || (arg.equals("-help")) || (arg.equals("/h"))
					|| (arg.equals("/help"))) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("proj=")) {
				filename = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("cnvFile=")) {
				cnvFile = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("cnvControl=")) {
				cnvControl = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("pairFile=")) {
				pairFile = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("output=")) {
				output = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("dir=")) {
				dir = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("numCNVS=")) {
				numCNVs = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("CN=")) {
				CN = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("log=")) {
				logfile = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("qcFile=")) {
				qcFile = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("-default")) {
				defaults = true;
				numArgs--;
			} else if (arg.startsWith("lrr")) {
				lrr = ext.parseDoubleArg(arg);
				numArgs--;
			} else if (arg.startsWith("call")) {
				call = ext.parseDoubleArg(arg);
				numArgs--;
			} else if (filter.isCommandLineFilterInEffect(arg)) {
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			proj.setLog(new Logger(proj.PROJECT_DIRECTORY.getValue() + (dir == null ? "" : dir)
														 + "concordLog"));

			determineConcordance(proj, cnvFile, dir, pairFile, qcFile, cnvControl, filter, numCNVs, CN,
													 lrr, call, output);
			long endTime = System.currentTimeMillis();
			long finalTime = endTime - startTime;
			System.out.println("total time: " + finalTime);
		} catch (Exception e) {
			e.printStackTrace();

		}

	}
}

