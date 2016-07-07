package cnv.qc;

import cnv.analysis.ProjectCNVFiltering;
import cnv.filesys.Project;
import cnv.var.SampleData;
import common.Array;
import common.CNVFilter;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;
import common.CNVFilter.CNVFilterPass;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Vector;

import filesys.CNVariant;
import filesys.CNVariantHash;

public class CNVConcordance {
	private static final String[] REPORT = { "Total CNVs Compared", "Total Overlapping CNVs", "Total Significantly Overlapping CNVs", "Total Perfectly overlapping cnvs", "Average Individual Overlap", "Average Individual Significant Overlap", "Average Individual Perfect Overlap", "Average Individual Overlap Score", "Global Overlap", "Global Significant Overlap", "Global Perfect Overlap", "Average Global Overlap Score" };
	// private static final int WARN_NUM_CNVs = 1000;
	public static final String COMMAND_PROJECT = "proj=";
	public static final String COMMAND_CNV = "cnvFile=";
	public static final String COMMAND_DUPLICATE = "duplicateFile=";
	public static final String COMMAND_LRRSD_THRESHOLD = "lrrSD=";
	public static final String COMMAND_OUTPUT = "output=";
	public static final String COMMAND_DEFUALTS = "-defaults";
	public static final String COMMAND_DIR = "dir=";
	public static final String COMMAND_NUMCNVS = "numCNVS=";
	public static final String COMMAND_CNV_CONCORDANCE = "cnvConcordance";
	public static final String COMMAND_CNV_CONCORDANCE_DESCRIPTION = "- compute the concordance between replicates in a cnv file";
	private Project proj;
	private String[][] duplicates;
	private String report;
	private CNVariantHash cNVariantHash;
	private SampleData sampleData;
	private ComparisionIndividualResults[] comparisionResults;
	private boolean fail;
	private CNVFilter filter;
	private int numCNVs;

	public CNVConcordance(Project proj, String[][] duplicates, CNVariantHash cNVariantHash, CNVFilter filter, int numCNVs) {
		this.proj = proj;
		this.sampleData = proj.getSampleData(0, false);
		this.duplicates = duplicates;
		this.cNVariantHash = cNVariantHash;
		this.filter = filter;
		this.fail = validDuplicateFormat();
		this.numCNVs = numCNVs;
	}

	public void determineConcordance() {
		this.comparisionResults = compareAll();
		if (!this.fail) {
			ComparisionGlobalResults cgr = new ComparisionGlobalResults(this.comparisionResults);
			cgr.populateMetrics();
			this.report = cgr.getReport();
		} else {
			this.proj.getLog().reportError("Error - comparison has failed");
		}
	}

	public String getReport() {
		return this.report;
	}

	private ComparisionIndividualResults[] compareAll() {
		if (!this.fail) {
			int numComp = 0;
			ArrayList<ComparisionIndividualResults> allResults = new ArrayList<ComparisionIndividualResults>();
			for (int i = 0; i < this.duplicates.length; i++) {
				Hashtable<String, String> track = new Hashtable<String, String>();
				for (int j = 0; j < this.duplicates[i].length; j++) {
					for (int k = 0; k < this.duplicates[i].length; k++) {
						String currentComp = j + "V" + k;
						if ((k != j) && (!track.containsKey(currentComp))) {
							numComp++;

							track.put(currentComp, currentComp);
							track.put(k + "V" + j, currentComp);
							String ind1 = this.sampleData.lookup(this.duplicates[i][j])[1];
							String ind2 = this.sampleData.lookup(this.duplicates[i][k])[1];

							ComparisionIndividualResults results = compareInds(ind1, ind2);
							if (results.getTotalCNVCount() > 0) {
								allResults.add(results);
							}
						}
					}
				}
			}
			System.out.println(numComp + " This many comparisons");
			return (ComparisionIndividualResults[]) allResults.toArray(new ComparisionIndividualResults[allResults.size()]);
		}
		return null;
	}

	private ComparisionIndividualResults compareInds(String ind1, String ind2) {
		CNVariant[] ind1CNVs = extractVariants(this.cNVariantHash.getDataFor(ind1));
		CNVariant[] ind2CNVs = extractVariants(this.cNVariantHash.getDataFor(ind2));
		if ((ind1CNVs == null) || (ind2CNVs == null) || ((ind1CNVs.length == 0) && (ind2CNVs.length == 0))) {
			return new ComparisionIndividualResults(ind1, ind2, 0, 0);
		}
		if ((ind1CNVs.length > this.numCNVs) || (ind2CNVs.length > this.numCNVs)) {
			return new ComparisionIndividualResults(ind1, ind2, 0, 0);
		}
		ind1CNVs = filterCNVs(this.filter, ind1CNVs, this.proj);
		ind2CNVs = filterCNVs(this.filter, ind2CNVs, this.proj);
		if ((ind1CNVs == null) || (ind2CNVs == null) || ((ind1CNVs.length == 0) && (ind2CNVs.length == 0))) {
			return new ComparisionIndividualResults(ind1, ind2, 0, 0);
		}
		if (ind1CNVs.length > 1000) {
			this.proj.getLog().report("Warning - " + ind1 + " has more than " + 1000 + " CNVs (" + ind1CNVs.length + "), this could lead to skewed concordance rates");
		}
		if (ind2CNVs.length > 1000) {
			this.proj.getLog().report("Warning - " + ind2 + " has more than " + 1000 + " CNVs (" + ind1CNVs.length + "), this could lead to skewed concordance rates");
		}
		if (ind1CNVs.length == 0) {
			return new ComparisionIndividualResults(ind1, ind2, 0, ind2CNVs.length);
		}
		if (ind2CNVs.length == 0) {
			return new ComparisionIndividualResults(ind1, ind2, ind1CNVs.length, 0);
		}
		return compareIndCNVs(ind1CNVs, ind2CNVs);
	}

	private boolean validDuplicateFormat() {
		boolean valid = true;
		if (this.duplicates == null) {
			this.proj.getLog().report("Error - replicates were not defined");
			valid = false;
		} else {
			int length = 0;
			for (int i = 0; i < this.duplicates.length; i++) {
				if (this.duplicates[i] == null) {
					this.proj.getLog().report("Error - replicates were not defined for replicate " + (i + 1));
					valid = false;
				} else {
					length = this.duplicates[i].length;
				}
			}
			if (valid) {
				for (int i = 0; i < this.duplicates.length; i++) {
					if (this.duplicates[i].length != length) {
						this.proj.getLog().report("Error - mismatched array sizes were found in replicate array");
						valid = false;
					} else {
						for (int j = 0; j < this.duplicates[i].length; j++) {
							if (this.duplicates[i][j] != null) {
								String[] ind = this.sampleData.lookup(this.duplicates[i][j]);
								if (ind == null) {
									this.proj.getLog().report("Error - did not find " + this.duplicates[i][j] + " in the sample data file " + proj.SAMPLE_DATA_FILENAME.getValue());
									valid = false;
								}
							}
						}
					}
				}
			}
		}
		return !valid;
	}

	private static CNVariant[] extractVariants(Hashtable<String, CNVariant[]> indCNVS) {
		ArrayList<CNVariant> cnvs = new ArrayList<CNVariant>();
		ArrayList<String> arr = Collections.list(indCNVS.keys());
		for (int i = 0; i < arr.size(); i++) {
			CNVariant[] chrCNVs = (CNVariant[]) indCNVS.get(arr.get(i));
			for (int j = 0; j < chrCNVs.length; j++) {
				cnvs.add(chrCNVs[j]);
			}
		}
		return (CNVariant[]) cnvs.toArray(new CNVariant[cnvs.size()]);
	}

	private static ComparisionIndividualResults compareIndCNVs(CNVariant[] ind1CNVs, CNVariant[] ind2CNVs) {
		ComparisionIndividualResults currentComparison = new ComparisionIndividualResults(ind1CNVs[0].getFamilyID() + "\t" + ind1CNVs[0].getIndividualID(), ind2CNVs[0].getFamilyID() + "\t" + ind2CNVs[0].getIndividualID(), ind1CNVs.length, ind2CNVs.length);
		boolean[] foundInd2 = new boolean[ind2CNVs.length];
		Arrays.fill(foundInd2, false);
		for (int i = 0; i < ind1CNVs.length; i++) {
			boolean foundSignificant = false;
			for (int j = 0; j < ind2CNVs.length; j++) {
				if (ind1CNVs[i].getCN() == ind2CNVs[j].getCN()) {
					if (ind1CNVs[i].overlaps(ind2CNVs[j])) {
						double minOverLapScore = Math.min(ind1CNVs[i].overlapScore(ind2CNVs[j]), ind2CNVs[j].overlapScore(ind1CNVs[i]));
						currentComparison.addOverlapScore(minOverLapScore);
						if ((!foundSignificant) && (foundInd2[j])) {
							currentComparison.addOverlap();
						}
					}
					if ((ind1CNVs[i].significantOverlap(ind2CNVs[j])) && (ind2CNVs[j].significantOverlap(ind1CNVs[i]))) {
						if ((!foundSignificant) && (foundInd2[j])) {
							currentComparison.addSigOverlap();
						}
						foundSignificant = true;
						foundInd2[j] = true;
					}
					if (ind1CNVs[i].equals(ind2CNVs[j])) {
						currentComparison.addPerfect();
					}
				}
			}
		}
		return currentComparison;
	}

	private static class ComparisionIndividualResults {
//		private String ind1;
//		private String ind2;
		private int numCNVind1;
		private int numCNVind2;
		private int numOverlap;
		private int numSigOverlap;
		private int numPerfectOverlap;
		private ArrayList<Double> overlapScores;

		public ComparisionIndividualResults(String ind1, String ind2, int numCNVind1, int numCNVind2) {
//			this.ind1 = ind1;
//			this.ind2 = ind2;
			this.numCNVind1 = numCNVind1;
			this.numCNVind2 = numCNVind2;
			this.numOverlap = 0;
			this.numSigOverlap = 0;
			this.numPerfectOverlap = 0;
			this.overlapScores = new ArrayList<Double>(getTotalCNVCount());
		}

		public void addOverlap() {
			this.numOverlap += 1;
		}

		public void addPerfect() {
			this.numPerfectOverlap += 1;
		}

		public int getNumPerfectOverlap() {
			return this.numPerfectOverlap;
		}

		public void addOverlapScore(double overlapScore) {
			if (!Double.isNaN(overlapScore)) {
				this.overlapScores.add(Double.valueOf(overlapScore));
			}
		}

		public double getAvgOverlapScore() {
			if (this.overlapScores.size() > 0) {
				return Array.mean(Array.toDoubleArray(this.overlapScores));
			}
			return 0.0D;
		}

		public void addAllScores(ArrayList<Double> globalOverLapScores) {
			for (int i = 0; i < this.overlapScores.size(); i++) {
				globalOverLapScores.add((Double) this.overlapScores.get(i));
			}
		}

		public void addSigOverlap() {
			this.numSigOverlap += 1;
		}

		public int getNumOverlap() {
			return this.numOverlap;
		}

		public int getNumSigOverlap() {
			return this.numSigOverlap;
		}

		public int getTotalCNVCount() {
			return this.numCNVind1 + this.numCNVind2;
		}

		public double getPercentOverlap() {
			if (getTotalCNVCount() == 0) {
				return 0.0D;
			}
			return 2.0D * this.numOverlap / getTotalCNVCount();
		}

		public double getPercentSigOverlap() {
			if (getTotalCNVCount() == 0) {
				return 0.0D;
			}
			return 2.0D * this.numSigOverlap / getTotalCNVCount();
		}

		public double getPercentPerfectOverlap() {
			if (getTotalCNVCount() == 0) {
				return 0.0D;
			}
			return 2.0D * this.numPerfectOverlap / getTotalCNVCount();
		}
	}

	private static class ComparisionGlobalResults {
		private CNVConcordance.ComparisionIndividualResults[] indResults;
		private int totalCNVs;
		private int totalOverlap;
		private int totalSigOverlap;
		private int totalPerfectOverlap;
		private double avgIndOverlap;
		private double avgIndSigOverlap;
		private double avgIndPerfectOverlap;
		private double avgIndOverlapScore;
		private double globalOverlap;
		private double globalSigOverlap;
		private double globalPerfectOverLap;
		private double avgGlobalOverlapScore;
		private ArrayList<Double> globalOverLapScores;

		public ComparisionGlobalResults(CNVConcordance.ComparisionIndividualResults[] indResults, int totalCNVs, int totalOverlap, int totalSigOverlap, int totalPerfectOverlap, double avgIndOverlap, double avgIndSigOverlap, double avgIndPerfectOverlap, double avgIndOverlapScore, double globalOverlap, double globalSigOverlap, double avgGlobalOverlapScore) {
			this.indResults = indResults;
			this.totalCNVs = totalCNVs;
			this.totalOverlap = totalOverlap;
			this.totalSigOverlap = totalSigOverlap;
			this.totalPerfectOverlap = totalPerfectOverlap;
			this.avgIndOverlap = avgIndOverlap;
			this.avgIndSigOverlap = avgIndSigOverlap;
			this.avgIndPerfectOverlap = avgIndPerfectOverlap;
			this.avgIndOverlapScore = avgIndOverlapScore;
			this.globalOverlap = globalOverlap;
			this.globalSigOverlap = globalSigOverlap;
			this.avgGlobalOverlapScore = avgGlobalOverlapScore;
			this.globalOverLapScores = new ArrayList<Double>();
		}

		public ComparisionGlobalResults(CNVConcordance.ComparisionIndividualResults[] comparisionIndividualResults) {
			this(comparisionIndividualResults, 0, 0, 0, 0, 0.0D, 0.0D, 0.0D, 0.0D, 0.0D, 0.0D, 0.0D);
		}

		public void populateMetrics() {
			for (int i = 0; i < this.indResults.length; i++) {
				this.totalCNVs += this.indResults[i].getTotalCNVCount();
				this.totalOverlap += this.indResults[i].getNumOverlap();
				this.totalSigOverlap += this.indResults[i].getNumSigOverlap();
				this.totalPerfectOverlap += this.indResults[i].getNumPerfectOverlap();
				this.avgIndOverlap += this.indResults[i].getPercentOverlap();
				this.avgIndSigOverlap += this.indResults[i].getPercentSigOverlap();
				this.avgIndPerfectOverlap += this.indResults[i].getPercentPerfectOverlap();
				this.avgIndOverlapScore += this.indResults[i].getAvgOverlapScore();
				this.indResults[i].addAllScores(this.globalOverLapScores);
			}
			this.avgIndOverlap /= this.indResults.length;
			this.avgIndSigOverlap /= this.indResults.length;
			this.avgIndPerfectOverlap /= this.indResults.length;
			this.avgIndOverlapScore /= this.indResults.length;
			this.globalOverlap = getGlobalPercentOverlap();
			this.globalSigOverlap = getGlobalPercentSigOverlap();
			this.globalPerfectOverLap = getGlobalPercentPerfectOverlap();
			this.avgGlobalOverlapScore = getAvgGlobalOverlapScore();
		}

		private double getAvgGlobalOverlapScore() {
			return Array.mean(Array.toDoubleArray(this.globalOverLapScores));
		}

		private double getGlobalPercentOverlap() {
			return 2.0D * this.totalOverlap / this.totalCNVs;
		}

		private double getGlobalPercentSigOverlap() {
			return 2.0D * this.totalSigOverlap / this.totalCNVs;
		}

		private double getGlobalPercentPerfectOverlap() {
			return 2.0D * this.totalPerfectOverlap / this.totalCNVs;
		}

		public String getReport() {
			return this.totalCNVs + "\t" + this.totalOverlap + "\t" + this.totalSigOverlap + "\t" + this.totalPerfectOverlap + "\t" + this.avgIndOverlap + "\t" + this.avgIndSigOverlap + "\t" + this.avgIndPerfectOverlap + "\t" + this.avgIndOverlapScore + "\t" + this.globalOverlap + "\t" + this.globalSigOverlap + "\t" + this.globalPerfectOverLap + "\t" + this.avgGlobalOverlapScore;
		}
	}

	private static String[][] loadDuplicates(String duplicateFile) {
		String[] load = HashVec.loadFileToStringArray(duplicateFile, false, null, false);
		String[][] duplicates = new String[load.length][];
		for (int i = 0; i < load.length; i++) {
			String[] tmp = load[i].split("\t");
			ArrayList<String> compareDef = new ArrayList<String>();
			for (int j = 0; j < tmp.length; j++) {
				if (!tmp[j].equals("")) {
					compareDef.add(tmp[j]);
				}
			}
			duplicates[i] = ((String[]) compareDef.toArray(new String[compareDef.size()]));
		}
		return duplicates;
	}

	public static CNVariant[] filterCNVs(CNVFilter cnvFilter, CNVariant[] cnvs, Project proj) {
		ArrayList<CNVariant> filtered = new ArrayList<CNVariant>(cnvs.length);
		for (int i = 0; i < cnvs.length; i++) {
			CNVFilterPass filterPass = cnvFilter.getCNVFilterPass(cnvs[i]);
			if (filterPass.passedFilter()) {
				if ((filterPass.isCentromeric()) && (cnvFilter.isBreakupCentromeres())) {
					filtered.add(cnvFilter.breakUpCentromere(filterPass, cnvs[i])[0]);
					filtered.add(cnvFilter.breakUpCentromere(filterPass, cnvs[i])[1]);
				} else {
					filtered.add(cnvs[i]);
				}
			} else if (filterPass.isIndIsExcluded()) {
				return null;
			}
		}
		return (CNVariant[]) filtered.toArray(new CNVariant[filtered.size()]);
	}

	public static void determineConcordance(Project proj, String cnvFile, String dir, String duplicateFile, CNVFilter filter, int numCNVs, int CN, String output) {
		if ((duplicateFile == null) || (duplicateFile.equals(""))) {
			proj.getLog().reportError("Error - a file of duplicates must be provided to determine concordance");
			return;
		}
		if (((cnvFile == null) || (cnvFile.equals(""))) && (dir == null)) {
			proj.getLog().reportError("Error - a file of cnvs must be provided to determine concordance");
			return;
		}
		filter.setCN(CN);
		String[][] duplicates = loadDuplicates(proj.PROJECT_DIRECTORY.getValue() + duplicateFile);
		if (dir != null) {
			String[] cnvFiles = Files.list(proj.PROJECT_DIRECTORY.getValue() + dir, ".cnv", false);
			cnvFiles = Files.toFullPaths(cnvFiles, proj.PROJECT_DIRECTORY.getValue() + dir);
			proj.getLog().report(Array.toStr(cnvFiles));
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + dir + output));
				int start = filter.getMinNumMarkers();
				//int stop = filter.getMaxNumMarkers();
				CNVariantHash[] cNVariantHash = new CNVariantHash[cnvFiles.length];
				for (int i = 0; i < cnvFiles.length; i++) {
					cNVariantHash[i] = CNVariantHash.load(cnvFiles[i], 1, false, proj.getLog());
					for (int j = 0; j < REPORT.length; j++) {
						writer.print(((i == 0) && (j == 0) ? "" : "\t") + REPORT[j] + "." + ext.rootOf(cnvFiles[i]));
					}
					writer.print("\tnumberOFProbes." + ext.rootOf(cnvFiles[i]));
				}
				writer.println();
				for (int j = start; j <= 100; j++) {
					filter.setMinNumMarkers(j);
					for (int i = 0; i < cnvFiles.length; i++) {
						long time = System.currentTimeMillis();
						proj.getLog().report(ext.getTime() + " Info - beginning comparision for " + cnvFiles[i]);

						CNVConcordance cnvConcordance = new CNVConcordance(proj, duplicates, cNVariantHash[i], filter, numCNVs);
						cnvConcordance.determineConcordance();
						writer.print((i == 0 ? "" : "\t") + cnvConcordance.getReport() + "\t" + j);
						proj.getLog().report(ext.getTime() + " Info - finished comparision for " + cnvFiles[i] + " and took " + ext.getTimeElapsed(time));
					}
					writer.println();
				}
				writer.close();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + proj.PROJECT_DIRECTORY.getValue() + output);
				proj.getLog().reportException(e);
			}
		} else {
			CNVariantHash cNVariantHash = CNVariantHash.load(proj.PROJECT_DIRECTORY.getValue() + cnvFile, 1, false, proj.getLog());
			CNVConcordance cnvConcordance = new CNVConcordance(proj, duplicates, cNVariantHash, filter, numCNVs);
			cnvConcordance.determineConcordance();
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + output));
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
		params[11] = "#a path (relative to the project directory) containing multiple \".cnv\" files, all will be analyzed and the cnvFile= command will be overridden";
		params[12] = "#dir=";

		params[13] = "#maximum number of cnvs per individual";
		params[14] = "#numCNVS=";

		System.out.println(Array.toStr(CNVFilter.getDefaultCNVParams()));
		params = (String[]) Array.concatAll(params, new String[][] { CNVFilter.getDefaultCNVParams() });

		return params;
	}

	public static void fromParameters(String filename, Logger log) {
		Vector<String> params = Files.parseControlFile(filename, "cnvConcordance", getParserParams(), log);
		if (params != null) {
			main(Array.toStringArray(params));
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String logfile = null;
		String cnvFile = null;
		String duplicateFile = null;
		String dir = null;
		int numCNVs = 2147483647;
		boolean defaults = false;
		String output = "cnv.concordance.txt";
		int CN = -1;

		String usage = "\njlDev.CNVConcordance requires 0-1 arguments\n";
		usage = usage + "   (1) project filename  (i.e. proj=" + filename + " (no default))\n";
		usage = usage + "   (2) cnvFile  (i.e.cnvFile=" + filename + " (no default))\n";
		usage = usage + "   (3) duplicate file  (i.e. duplicateFile=" + filename + " (no default))\n";
		usage = usage + "   OPTIONAL:";
		usage = usage + "   (4) output file name  (i.e.output=" + output + " (default))\n";
		usage = usage + "   (5) log file  (i.e. log=" + filename + " (no default))\n";
		usage = usage + "\t (6) For cnv filtering, use the default values (i.e. -default ( not the default))\n";
		usage = usage + "\t (7) a directory containing multiple cnv files (i.e. dir= ( no default))\n";
		usage = usage + "\t (8) maximum number of cnvs (i.e. numCNVS=" + numCNVs + " (default))\n";

		usage = usage + "\t (9) further usage:\n" + Array.toStr(CNVFilter.getDefaultCNVParams());
		Project proj;
		if (ext.indexOfStr("proj=", args, true, false) >= 0) {
			proj = new Project(ext.parseStringArg(args[ext.indexOfStr("proj=", args, true, false)], ""), logfile, false);
		} else {
			proj = new Project(filename, logfile, false);
		}
		CNVFilter filter = ProjectCNVFiltering.setupCNVFilterFromArgs(proj, args, null, defaults, proj.getLog());
		if (ext.indexOfStr("-defaults", args) >= 0) {
		    ProjectCNVFiltering.setCNVDefaults(filter, proj);
		}
		for (int i = 0; i < args.length; i++) {
			if ((args[i].equals("-h")) || (args[i].equals("-help")) || (args[i].equals("/h")) || (args[i].equals("/help"))) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("cnvFile=")) {
				cnvFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("duplicateFile=")) {
				duplicateFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("output=")) {
				output = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("dir=")) {
				dir = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("numCNVS=")) {
				numCNVs = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("CN=")) {
				CN = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("-default")) {
				defaults = true;
				numArgs--;
			} else if (filter.isCommandLineFilterInEffect(args[i])) {
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			proj.setLog(new Logger(proj.PROJECT_DIRECTORY.getValue() + (dir == null ? "" : dir) + "concordLog"));

			determineConcordance(proj, cnvFile, dir, duplicateFile, filter, numCNVs, CN, output);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
