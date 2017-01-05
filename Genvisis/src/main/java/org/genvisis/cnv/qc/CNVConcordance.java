package org.genvisis.cnv.qc;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import org.genvisis.cnv.analysis.ProjectCNVFiltering;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.CNVFilter;
import org.genvisis.common.CNVFilter.CNVFilterPass;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.CNVariantHash;

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
		private double avgIndRelOverlap;
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
			sb.append(avgIndRelOverlap).append("\t");
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
				avgIndRelOverlap += indResult.getRelativePercentOverlap();
				avgIndRelSigOverlap += indResult.getRelativePercentSigOverlap();
				avgIndRelPerfectOverlap += indResult.getRelativePercentPerfectOverlap();
				indResult.addAllScores(globalOverLapScores);
			}
			avgIndOverlap /= indResults.length;
			avgIndSigOverlap /= indResults.length;
			avgIndPerfectOverlap /= indResults.length;
			avgIndOverlapScore /= indResults.length;
			avgIndUnpaired = (double) totalUnpaired / indResults.length;
			avgIndRelOverlap /= indResults.length;
			avgIndRelSigOverlap /= indResults.length;
			avgIndRelPerfectOverlap /= indResults.length;
			globalOverlap = (double) totalOverlap / totalCNVs;
			globalSigOverlap = (double) totalSigOverlap / totalCNVs;
			globalPerfectOverlap = (double) totalPerfectOverlap / totalCNVs;
			avgGlobalOverlapScore = Array.mean(globalOverLapScores);
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
				return Array.mean(overlapScores);
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

		public double getRelativePercentOverlap() {
			return getRelative(numOverlap);
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
			//TODO for printing individual information
			return "";
		}
	}

	private static final String[] REPORT = {"Total CNVs Compared", "Total Overlapping CNVs",
	                                        "Total Significantly Overlapping CNVs",
	                                        "Total Perfectly overlapping cnvs",
	                                        "Total Unpaired CNVs",
	                                        "Average Individual Overlap",
	                                        "Average Individual Significant Overlap",
	                                        "Average Individual Perfect Overlap","Average Individual Overlap Score",
	                                        "Average Individual Unpaired",
	                                        "Average Relative Overlap",
	                                        "Average Relative Significant Overlap",
	                                        "Average Relative Perfect Overlap",
	                                        "Global Overlap",
	                                        "Global Significant Overlap", "Global Perfect Overlap",
	                                        "Average Global Overlap Score"};
	public static final String COMMAND_PROJECT = "proj=";
	public static final String COMMAND_CNV = "cnvFile=";
	public static final String COMMAND_DUPLICATE = "duplicateFile=";
	public static final String COMMAND_LRRSD_THRESHOLD = "lrrSD=";
	public static final String COMMAND_OUTPUT = "output=";
	public static final String COMMAND_DEFUALTS = "-defaults";
	public static final String COMMAND_DIR = "dir=";
	public static final String COMMAND_NUMCNVS = "numCNVS=";
	public static final String COMMAND_CNV_CONCORDANCE = "cnvConcordance";
	public static final String COMMAND_CNV_CONCORDANCE_DESCRIPTION =
	                                                               "- compute the concordance between replicates in a cnv file";
	private final Project proj;
	private final String[][] duplicates;
	private final double[] lrr;
	private final double[] callRate;
	private String report;
	private final CNVariantHash cNVariantHash;
	private ComparisionIndividualResults[] comparisionResults;
	private boolean fail;
	private final CNVFilter filter;
	private final int numCNVs;
	private final double maxLrr;
	private final double minCallRate;
	private SampleData sampleData;

	public CNVConcordance(Project proj, String[][] duplicates, double[] lrr, double maxLrr,
	                      double[] callRate, double minCallRate, CNVariantHash cNVariantHash,
	                      CNVFilter filter, int numCNVs) {
		this.proj = proj;
		sampleData = proj.getSampleData(SampleData.MINIMAL_SAMPLE_DATA_HEADER.length, false);
		this.lrr = lrr;
		this.callRate = callRate;
		this.maxLrr = maxLrr;
		this.minCallRate = minCallRate;
		this.duplicates = duplicates;
		this.cNVariantHash = cNVariantHash;
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
			ArrayList<ComparisionIndividualResults> allResults =
			                                                   new ArrayList<ComparisionIndividualResults>();
			for (String[] duplicate : duplicates) {
				if (lrr[numComp] > maxLrr) {
					System.out.println("Warning - " + duplicate[0] + " has LRR > " + maxLrr + " ("
					                   + lrr[numComp] + "). Excluded from final analysis");
				} else if (callRate[numComp] < minCallRate) {
					System.out.println("Warning - " + duplicate[0] + " has callrate < " + minCallRate + " ("
					                   + callRate[numComp] + "). Excluded from final analysis");
				} else {
					// Try dna\tdna first
					String ind1 = duplicate[0] + "\t" + duplicate[0];
					String ind2 = duplicate[1] + "\t" + duplicate[1];

					ComparisionIndividualResults results = compareInds(ind1, ind2);

					if (results.getTotalCNVCount() <= 0) {
						// Try looking up actual fid\tiid next
						ind1 = sampleData.lookup(duplicate[0])[1];
						ind2 = sampleData.lookup(duplicate[1])[1];
						results = compareInds(ind1, ind2);
					}

					if (results.getTotalCNVCount() > 0) {
						allResults.add(results);
					} else {
						System.out.println("No common cnvs found for pair: " + duplicate[0] + ", "
						                   + duplicate[1]);
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
	 // match

	private ComparisionIndividualResults compareInds(String ind1, String ind2) {
		CNVariant[] ind1CNVs = extractVariants(cNVariantHash.getDataFor(ind1));
		CNVariant[] ind2CNVs = extractVariants(cNVariantHash.getDataFor(ind2));
		if ((ind1CNVs == null) || (ind2CNVs == null)
		    || ((ind1CNVs.length == 0) && (ind2CNVs.length == 0))) {
			return new ComparisionIndividualResults(ind1, ind2, 0, 0);
		}
		if ((ind1CNVs.length > numCNVs) || (ind2CNVs.length > numCNVs)) {
			return new ComparisionIndividualResults(ind1, ind2, 0, 0);
		}
		ind1CNVs = filterCNVs(filter, ind1CNVs, proj);
		ind2CNVs = filterCNVs(filter, ind2CNVs, proj);
		if ((ind1CNVs == null) || (ind2CNVs == null)
		    || ((ind1CNVs.length == 0) && (ind2CNVs.length == 0))) {
			return new ComparisionIndividualResults(ind1, ind2, 0, 0);
		}
		if (ind1CNVs.length > 1200) {
			proj.getLog().report("Warning - " + ind1 + " has more than " + 1200 + " CNVs ("
			                     + ind1CNVs.length + "), this could lead to skewed concordance rates");
		}
		if (ind2CNVs.length > 1200) { // changed to 1200 from 1000 to shorten print output
			proj.getLog().report("Warning - " + ind2 + " has more than " + 1200 + " CNVs ("
			                     + ind2CNVs.length + "), this could lead to skewed concordance rates");
		}
		if (ind1CNVs.length == 0) {
			return new ComparisionIndividualResults(ind1, ind2, 0, ind2CNVs.length);
		}
		if (ind2CNVs.length == 0) {
			return new ComparisionIndividualResults(ind1, ind2, ind1CNVs.length, 0);
		}
		return compareIndCNVs(ind1CNVs, ind2CNVs);
	}

	/*
	 * private boolean validDuplicateFormat() { boolean valid = true; if (duplicates == null) {
	 * proj.getLog().report("Error - replicates were not defined"); valid = false; } else { int length
	 * = 0; for (int i = 0; i < duplicates.length; i++) { if (duplicates[i] == null) {
	 * proj.getLog().report("Error - replicates were not defined for replicate " + (i + 1)); valid =
	 * false; } else { length = duplicates[i].length; } } if (valid) { for (String[] duplicate :
	 * duplicates) { if (duplicate.length != length) {
	 * proj.getLog().report("Error - mismatched array sizes were found in replicate array"); valid =
	 * false; } else { for (String element : duplicate) { if (element != null) { String[] ind =
	 * sampleData.lookup(element); if (ind == null) { proj.getLog() .report("Error - did not find " +
	 * element + " in the sample data file " + proj.SAMPLE_DATA_FILENAME.getValue()); valid = false; }
	 * } } } } } } return !valid; }
	 */

	public static void determineConcordance(Project proj, String cnvFile, String dir,
	                                        String duplicateFile, String qcFile, CNVFilter filter,
	                                        int numCNVs, int CN, double lrrMax, double callMin,
	                                        String output) {
		if ((duplicateFile == null) || (duplicateFile.equals(""))) {
			proj.getLog()
			    .reportError("Error - a file of duplicates must be provided to determine concordance");
			return;
		}
		if (((cnvFile == null) || (cnvFile.equals(""))) && (dir == null)) {
			proj.getLog().reportError("Error - a file of cnvs must be provided to determine concordance");
			return;
		}
		double[][] qc = loadQC(proj.PROJECT_DIRECTORY.getValue() + qcFile);
		double[] lrr = qc[0];
		double[] callRate = qc[1];
		filter.setCN(CN);
		String[][] duplicates = loadDuplicates(proj.PROJECT_DIRECTORY.getValue() + duplicateFile);
		if (dir != null) {
			String[] cnvFiles = Files.list(proj.PROJECT_DIRECTORY.getValue() + dir, ".cnv", false);
			cnvFiles = Files.toFullPaths(cnvFiles, proj.PROJECT_DIRECTORY.getValue() + dir);
			proj.getLog().report(Array.toStr(cnvFiles));
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + dir

				                                                    + output));
				int start = filter.getMinNumMarkers();
				CNVariantHash[] cNVariantHash = new CNVariantHash[cnvFiles.length];
				for (int i = 0; i < cnvFiles.length; i++) {
					cNVariantHash[i] = CNVariantHash.load(cnvFiles[i], 1, false, proj.getLog());
					for (int j = 0; j < REPORT.length; j++) {
						writer.print(((i == 0) && (j == 0) ? "" : "\t") + REPORT[j] + "."
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

						CNVConcordance cnvConcordance = new CNVConcordance(proj, duplicates, lrr, lrrMax,
						                                                   callRate, callMin, cNVariantHash[i],
						                                                   filter, numCNVs);
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
			CNVConcordance cnvConcordance = new CNVConcordance(proj, duplicates, lrr, lrrMax, callRate,
			                                                   callMin, cNVariantHash, filter, numCNVs);
			cnvConcordance.determineConcordance();
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()
				                                                    + output));
				writer.println(Array.toStr(REPORT));
				writer.println(cnvConcordance.getReport());

				writer.close();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + proj.PROJECT_DIRECTORY.getValue() + output);
				proj.getLog().reportException(e);
			}
			proj.getLog().report(Array.toStr(REPORT));
			proj.getLog().report(cnvConcordance.getReport());
		}
	}

	public static CNVariant[] filterCNVs(CNVFilter cnvFilter, CNVariant[] cnvs, Project proj) {
		ArrayList<CNVariant> filtered = new ArrayList<CNVariant>(cnvs.length);
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
				return null;
			}
		}
		return filtered.toArray(new CNVariant[filtered.size()]);
	}

	public static void fromParameters(String filename, Logger log) {
		Vector<String> params = Files.parseControlFile(filename, "cnvConcordance", getParserParams(),
		                                               log);
		if (params != null) {
			main(Array.toStringArray(params));
		}
	}

	public static String[] getParserParams() {
		String[] params = new String[15];
		params[0] = "#To intialize the cnv concorder, provide the following arguments";
		params[1] = "#the full path to a project properties file";
		params[2] = "proj=";
		params[3] = "#a path (relative to the project directory) to a cnv file ";
		params[4] = "cnvFile=";
		params[5] = "#a path (relative to the project directory) to a file of duplicates ";
		params[6] = "duplicateFile=";
		params[7] = "#a threshold for log R ratio Standard deviation, but this is not implemented yet";
		params[8] = "#lrrSD=";
		params[9] = "#a path (relative to the project directory) to a file for ouput";
		params[10] = "#output=";
		params[11] =
		           "#a path (relative to the project directory) containing multiple \".cnv\" files, all will be analyzed and the cnvFile= command will be overridden";
		params[12] = "#dir=";

		params[13] = "#maximum number of cnvs per individual";
		params[14] = "#numCNVS=";

		System.out.println(Array.toStr(CNVFilter.getDefaultCNVParams()));
		params = Array.concatAll(params, new String[][] {CNVFilter.getDefaultCNVParams()});

		return params;
	}

	private static ComparisionIndividualResults compareIndCNVs(CNVariant[] ind1CNVs,
	                                                           CNVariant[] ind2CNVs) {
		ComparisionIndividualResults currentComparison =
		                                               new ComparisionIndividualResults(ind1CNVs[0].getFamilyID()
		                                                                                + "\t"
		                                                                                + ind1CNVs[0].getIndividualID(),
		                                                                                ind2CNVs[0].getFamilyID() + "\t" + ind2CNVs[0].getIndividualID(),
		                                                                                ind1CNVs.length,
		                                                                                ind2CNVs.length);
		compareIndCNVHelper(currentComparison, ind1CNVs, ind2CNVs);
		compareIndCNVHelper(currentComparison, ind2CNVs, ind1CNVs);

		return currentComparison;
	}

	private static void compareIndCNVHelper(ComparisionIndividualResults currentComparison,
	                                        CNVariant[] cnvs1, CNVariant[] cnvs2) {
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

	private static CNVariant[] extractVariants(Hashtable<String, CNVariant[]> indCNVS) {
		ArrayList<CNVariant> cnvs = new ArrayList<CNVariant>();
		ArrayList<String> arr = Collections.list(indCNVS.keys());
		for (int i = 0; i < arr.size(); i++) {
			CNVariant[] chrCNVs = indCNVS.get(arr.get(i));
			for (CNVariant chrCNV : chrCNVs) {
				cnvs.add(chrCNV);
			}
		}
		return cnvs.toArray(new CNVariant[cnvs.size()]);
	}

	private static String[][] loadDuplicates(String duplicateFile) {
		String[] load = HashVec.loadFileToStringArray(duplicateFile, false, new int[] {0, 1}, false);
		String[][] duplicates = new String[load.length][];
		for (int i = 0; i < load.length; i++) {
			String[] tmp = load[i].split("\t");
			ArrayList<String> compareDef = new ArrayList<String>();
			for (int j = 0; j < tmp.length; j++) {
				if (!tmp[j].equals("")) {
					compareDef.add(tmp[j]);
				}
			}
			duplicates[i] = (compareDef.toArray(new String[compareDef.size()]));
		}
		return duplicates;
	}

	private static double[][] loadQC(String qcFile) { // Known Issue - genvisis-generated LRR file
	                                                  // might have headers sprinkled throughout
		String[] load = HashVec.loadFileToStringArray(qcFile, false, null, false);
		double[][] qc = new double[2][load.length];
		for (int i = 0; i < load.length; i++) {
			String[] line = load[i].split("\t");
			// skip header lines
			if (line[2].contains("LRR")) {
				continue;
			}
			qc[0][i] = Double.valueOf(line[2]);
			qc[1][i] = Double.valueOf(line[7]);
		}
		return qc;
	}

	/*
	 * private static double[] loadCallRate(String qcFile) { String[][] load =
	 * HashVec.loadFileToStringMatrix(qcFile, false, null, false); double[] callRate = new
	 * double[load.length]; for (int i=0; i< load.length; i++) { String[] parts =
	 * load[i][7].split("/t"); callRate[i]=Double.valueOf(parts[0]); } return callRate; }
	 */

	public static void main(String[] args) {
		long startTime = System.currentTimeMillis();
		int numArgs = args.length;
		String filename = null;
		String logfile = null;
		String cnvFile = null;
		String duplicateFile = null;
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
		usage = usage + "   (3) duplicate file  (i.e. duplicateFile=" + filename + " (no default))\n";
		usage = usage + "   (4) quality control file (i.e. qcFile=" + filename + " (no default))\n";
		usage = usage + "   OPTIONAL:";
		usage = usage + "   (5) output file name  (i.e.output=" + output + " (default))\n";
		usage = usage + "   (6) log file  (i.e. log=" + filename + " (no default))\n";
		usage = usage
		        + "\t (7) For cnv filtering, use the default values (i.e. -default ( not the default))\n";
		usage = usage + "\t (8) a directory containing multiple cnv files (i.e. dir= ( no default))\n";
		usage = usage + "\t (9) maximum number of cnvs (i.e. numCNVS=" + numCNVs + " (default))\n";
		usage = usage + "\t (10) lrr cutoff (i.e. lrr=" + lrr + " (default))\n";
		usage = usage + "\t (11) minimum callrate (i.e. call=" + call + " (default))\n";

		usage = usage + "\t (12) further usage:\n" + Array.toStr(CNVFilter.getDefaultCNVParams());
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
			} else if (arg.startsWith("duplicateFile=")) {
				duplicateFile = ext.parseStringArg(arg, "");
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

			determineConcordance(proj, cnvFile, dir, duplicateFile, qcFile, filter, numCNVs, CN, lrr,
			                     call, output);
			long endTime = System.currentTimeMillis();
			long finalTime = endTime - startTime;
			System.out.println("total time: " + finalTime);
		} catch (Exception e) {
			e.printStackTrace();

		}

	}
}

