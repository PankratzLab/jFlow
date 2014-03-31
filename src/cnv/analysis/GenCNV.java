package cnv.analysis;

import gwas.PhenoPrep;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;

import stats.LeastSquares;
import stats.RegressionModel;
import common.Array;
import common.Files;
import common.Logger;
import common.ext;
import cnv.filesys.MarkerFreqs;
import cnv.manage.ExportCNVsToPedFormat;

public class GenCNV implements Runnable {

	public static final String[] ANALYSIS_MODELS = { "ANY_CNV_COLLAPSED", "ANY_CNV_ORDERED_NOT_COLLAPSED", "HOMOZYGOUS_DELETIONS_ONLY" };
	public static final String[] ANALYSIS_TYPES = { "POSITION", "WINDOW" };
	private static final int[] DEFAULT_WINDOWS = { 0, 200000 };
	private static final int DEFAULT_LOCI_TO_RETURN = 100;
	public static final String[] ANALYSIS_SUMMARY_HEADER = { "MODEL", "TYPE", "PhenoType", "MinPvalue", "MinPvalueLocus", "NumPassing pval cutoff", "lociTested/lociTotal", "Top SignificantLoci..." };
	public static final boolean[][] ANALYSIS_MODEL_PARAMS = { { true, true, false, true, false }, { true, true, true, true, false }, { true, false, false, false, true } };
	public static final String[] GPHENO_HEADERS = { "FID", "IID" };
	public static final String[] GPED_HEADERS = { "markerName" };
	private static final int PHENO_START = 2;
	private static final String DEFAULT_LOG_NAME = "GenCNV.log";
	private static final String DEFAULT_BACKUP = "BackUp/";
	private static final String DEFAULT_RESULTS_BACKUP = "resultsBackUp/";
	private static final String[] DEFUALT_MISSING_VALUES = ext.MISSING_VALUES;
	private static final String ID_DELIMITER = "-";
	private Analysis[] analyses;
	private int threadID;

	public void run() {
		for (int i = 0; i < analyses.length; i++) {
			String directory = analyses[i].getDirectory();
			String logFile = directory + DEFAULT_LOG_NAME;
			if (Files.exists(logFile)) {
				System.out.println(logFile + "\n " + directory + DEFAULT_BACKUP);
				Files.backup(DEFAULT_LOG_NAME, directory, directory + DEFAULT_BACKUP);
			}
			Logger log = new Logger(logFile);
			String[] analysisFiles = analyses[i].getFiles();
			for (int k = 0; k < analysisFiles.length; k++) {
				log.report(ext.getTime() + " Thread " + threadID + ": analysis " + i + "/" + analyses.length + ": file " + k + "/" + analysisFiles.length);
				String file = directory + analysisFiles[k];
				GpedCNVGenos gpedCNVGenos = getGpedCNVGenos(file, log);
				log.report("Info - Found " + gpedCNVGenos.getGpedCNVGenos().length + " CNV genotypes in file " + file + " using thread " + threadID);
				runAnalysis(analyses[i], gpedCNVGenos, log, file);
			}
			computeBurden(analyses[i]);
		}
	}

	private static void computeBurden(Analysis analysis) {
		Pheno[] phenos = analysis.getPhenos();
		for (int i = 0; i < phenos.length; i++) {
			if (hasVariance(phenos[i])) {
				RegressionModel model = (RegressionModel) new LeastSquares(phenos[i].getArrayPheno(), Array.toDoubleArray(analysis.getBurdens()[i].getCounts()));
				if (!model.analysisFailed()) {
					analysis.getBurdens()[i].setBurdenPvalue(model.getOverallSig());
				}
			}
		}
	}

	// assign genotypes for inds with this phenotype
	// double[] indeps = assignedGenos(headers, gpedGenos[k], phenos[i], analysis.getBurdens()[i]);

	private static void runAnalysis(Analysis analysis, GpedCNVGenos gpedCNVGenos, Logger log, String file) {
		String[] headers = gpedCNVGenos.getHeaders();
		Pheno[] phenos = analysis.getPhenos();
		String[][] gpedGenos = gpedCNVGenos.getGpedCNVGenos();
		// iterate phenos
		for (int i = 0; i < phenos.length; i++) {
			if (hasVariance(phenos[i])) {
				int count = 0;
				int[] indices = ext.indexFactors(phenos[i].getArrayInds(), headers, true, true);
				for (int k = 0; k < gpedGenos.length; k++) {
					double[] indeps = Array.toDoubleArray(Array.subArray(gpedGenos[k], indices));
					analysis.getSignificance()[i].addNumTotal();
					assignBurdens(indeps, analysis.getBurdens()[i], log);
					// check if frequency of cnvs passes cutoff
					if (checkIndepsFreq(analysis, indeps)) {
						RegressionModel model = (RegressionModel) new LeastSquares(phenos[i].getArrayPheno(), indeps);
						if (!model.analysisFailed()) {
							count++;
							double overAllSig = model.getOverallSig();
							double minSig = analysis.getSignificance()[i].getMinPvalue();
							analysis.getSignificance()[i].addLoci(gpedGenos[k][0]);
							analysis.getSignificance()[i].addPval(overAllSig);
							analysis.getSignificance()[i].addNumTest();
							if (overAllSig < minSig) {
								analysis.getSignificance()[i].setMinPvalue(overAllSig);
								analysis.getSignificance()[i].setMinPvalLocus(gpedGenos[k][0]);
							}
							if (overAllSig < analysis.getPvalCutoff()) {
								analysis.getSignificance()[i].addNumPassingThreshold();
							}
						} else {
							log.reportError("Analysis failed for locus " + gpedGenos[k][0] + " , " + analysis.getSignificance()[0].getAnalysisModel() + "\t" + analysis.getSignificance()[0].getAnalysisType());
							log.reportError("Failed analysis for data in file" + file + ", position " + k + ", phenotype " + phenos[i].getPhenoName());
						}
					}
				}
				log.report(ext.getTime() + " Info - tested " + count + " individuals in file " + file + " with phenotype " + phenos[i].getPhenoName());
			} else {
				System.out.println("NO VAR " + phenos[i].getPhenoName());
			}

		}
	}

	private static boolean checkIndepsFreq(Analysis analysis, double[] indeps) {
		double excludeFreqBelow = analysis.getExcludeFreqBelow();
		double numNeeded = excludeFreqBelow * indeps.length;
		double numHave = 0;
		for (int i = 0; i < indeps.length; i++) {
			if (indeps[i] != 0) {
				numHave++;
			}
			if (numHave >= numNeeded) {
				return true;
			}
		}
		return false;
	}

	private static void assignBurdens(double[] genos, Burden burden, Logger log) {
		if (genos.length != burden.getCounts().length) {
			log.reportError("Error - the number of genotypes do not match the number of individuals in the burden test");
			System.exit(1);
		} else {
			for (int i = 0; i < genos.length; i++) {
				if (genos[i] != (double) 0) {
					burden.countIt(i);
					// System.out.println(burden.getCounts()[i]);
				}
			}
		}
	}

	public GenCNV(Analysis[] analyses, int threadID) {
		this.analyses = analyses;
		this.threadID = threadID;

	}

	public int getThreadID() {
		return threadID;
	}

	public void setAnalyses(Analysis[] analyses) {
		this.analyses = analyses;
	}

	public Analysis[] getAnalyses() {
		return analyses;
	}

	public static GpedCNVGenos getGpedCNVGenos(String gpedFile, Logger log) {
		int numLines = getNumLines(gpedFile, log);
		if (numLines < 2) {
			log.reportError("Error -data not found in " + gpedFile);
			System.exit(1);

		}
		String[] headers = null;
		String[][] gpedCNVGenos = new String[numLines - 1][];
		try {
			int index = 0;
			BufferedReader reader = new BufferedReader(new FileReader(gpedFile));
			if (!reader.ready()) {
				log.reportError("Error reading file \"" + gpedFile + "\"");
				System.exit(1);
			}
			headers = reader.readLine().trim().split("\t");
			if (!headers[0].equals(GPED_HEADERS[0])) {
				log.reportError("Error - Need a column header ending with the following suffixes; missing at least one");
				log.reportError("        " + Array.toStr(GPED_HEADERS, "  "));
				System.exit(1);
			}
			while (reader.ready()) {
				gpedCNVGenos[index] = reader.readLine().trim().split("\t");
				index++;
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + gpedFile + "\" not found in current directory");
			fnfe.printStackTrace();
			System.exit(1);

		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + gpedFile + "\"");
			ioe.printStackTrace();
			System.exit(1);
		}
		if (headers == null || gpedCNVGenos[0] == null) {
			log.reportError("Error retrieving data from file \"" + gpedFile + "\"");
			System.exit(1);
		}
		return new GpedCNVGenos(headers, gpedCNVGenos);
	}

	private static int getNumLines(String file, Logger log) {
		int numLines = 0;
		try {
			BufferedReader reader = new BufferedReader(new FileReader(file));
			if (!reader.ready()) {
				log.reportError("Error reading file \"" + file + "\"");
				System.exit(1);
			}
			while (reader.ready()) {
				reader.readLine();
				numLines++;
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + file + "\" not found in current directory");
			fnfe.printStackTrace();
			System.exit(1);

		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + file + "\"");
			ioe.printStackTrace();
			System.exit(1);

		}
		return numLines;
	}

	public static Pheno[] loadGPHENO(String gPhenoFIle, Logger log) {
		Pheno[] phenos = null;
		try {
			BufferedReader reader = Files.getAppropriateReader(gPhenoFIle);
			if (!reader.ready()) {
				log.reportError("Error reading file \"" + gPhenoFIle + "\"");
				System.exit(1);
			}
			String[] header = reader.readLine().trim().split("\t");
			checkHeader(header, log);
			checkPheno(header, log);
			phenos = getPhenos(header);
			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("\t");
				String sample = line[0] + "-" + line[1];
				for (int i = PHENO_START; i < header.length; i++) {
					// skip misses right away;
					if (!ext.isMissingValue(line[i])) {
						try {
							Double.parseDouble(line[i]);
							phenos[i - PHENO_START].setIndPheno(sample, Double.parseDouble(line[i]));

						} catch (NumberFormatException nfe) {
							// set to NAN for compatability with phenoprep
							phenos[i - PHENO_START].setIndPheno(sample, Double.NaN);

						}
					}
				}
			}
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + gPhenoFIle + "\" not found in current directory");
			fnfe.printStackTrace();
			System.exit(1);
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + gPhenoFIle + "\"");
			ioe.printStackTrace();
			System.exit(1);

		}
		return phenos;
	}

	private static boolean isMissing(String pheno) {
		boolean miss = false;
		for (int i = 0; i < DEFUALT_MISSING_VALUES.length; i++) {
			if (pheno.equals(DEFUALT_MISSING_VALUES[i])) {
				miss = true;
			}
		}
		return miss;
	}

	private static Pheno[] getPhenos(String[] header) {
		Pheno[] phenos = new Pheno[header.length - PHENO_START];
		for (int i = PHENO_START; i < header.length; i++) {
			phenos[i - PHENO_START] = new Pheno(header[i]);
		}
		return phenos;
	}

	private static boolean checkPheno(String[] header, Logger log) {
		boolean apheno = false;
		if (header.length < 3) {
			log.reportError("Error - Need at least one phenotype...");
			System.exit(1);
		} else {
			apheno = true;
		}
		return apheno;
	}

	private static boolean checkHeader(String[] header, Logger log) {
		boolean head = false;
		if (header[0].equals(GPHENO_HEADERS[0]) && header[1].equals(GPHENO_HEADERS[1])) {
			head = true;
		} else {
			log.reportError("Error - Need a column header ending with the following suffixes; missing at least one");
			log.reportError("        " + Array.toStr(GPHENO_HEADERS, "  "));
			System.exit(1);
		}
		return head;
	}

	public static class GpedCNVGenos {
		private String[] headers;
		private String[][] gpedCNVGenos;

		public GpedCNVGenos(String[] headers, String[][] gpedCNVGenos) {
			this.headers = headers;
			this.gpedCNVGenos = gpedCNVGenos;
		}

		public String[] getHeaders() {
			return headers;
		}

		public String[][] getGpedCNVGenos() {
			return gpedCNVGenos;
		}
	}

	public static class Pheno {
		private String phenoName;
		private ArrayList<String> inds;
		private ArrayList<Double> phenos;
		private Hashtable<String, Integer> hasPheno;

		public Pheno(String phenoName) {
			this.phenoName = phenoName;
			this.inds = new ArrayList<String>();
			this.phenos = new ArrayList<Double>();
			this.hasPheno = new Hashtable<String, Integer>();
		}

		public Hashtable<String, Integer> getHasPheno() {
			return hasPheno;
		}

		public String getPhenoName() {
			return phenoName;
		}

		public double[] getArrayPheno() {
			return toDoubleArray(phenos);
		}

		public String[] getArrayInds() {
			return toStringArray(inds);
		}

		public void setIndPheno(String ind, double value) {
			addPhenoValue(value);
			addInd(ind);
			definedPheno(ind, inds.size() - 1);
		}

		public void definedPheno(String ind, int index) {
			hasPheno.put(ind, index);
		}

		public void addPhenoValue(double value) {
			phenos.add(value);
		}

		public void addInd(String ind) {
			inds.add(ind);
		}

		public ArrayList<String> getInds() {
			return inds;
		}
	}

	public static class Burden {
		private String analysisModel;
		private String analysisType;
		private String pheno;
		private int[] counts;
		private double burdenPvalue;

		public Burden(String analysisModel, String analysisType, String pheno, int numInds) {
			super();
			this.analysisModel = analysisModel;
			this.analysisType = analysisType;
			this.pheno = pheno;
			this.counts = new int[numInds];
			this.burdenPvalue = 1;
		}

		public void countIt(int sampleIndex) {
			counts[sampleIndex]++;
		}

		public int[] getCounts() {
			return counts;
		}

		public String getPheno() {
			return pheno;
		}

		public double getBurdenPvalue() {
			return burdenPvalue;
		}

		public void setBurdenPvalue(double burdenPvalue) {
			this.burdenPvalue = burdenPvalue;
		}

		public String getAnalysisModel() {
			return analysisModel;
		}

		public String getAnalysisType() {
			return analysisType;
		}

	}

	public static class Analysis {
		private String directory;
		private String[] files;
		private Significance[] significance;
		private Pheno[] phenos;
		private Burden[] burdens;
		private double excludeFreqBelow;
		private double pvalCutoff;

		public Analysis(String directory, String[] files, Significance[] significance, Pheno[] phenos, Burden[] burdens, double excludeFreqBelow, double pvalCutoff) {
			this.directory = directory;
			this.files = files;
			this.significance = significance;
			this.phenos = phenos;
			this.burdens = burdens;
			this.excludeFreqBelow = excludeFreqBelow;
			this.pvalCutoff = pvalCutoff;
		}

		public Burden[] getBurdens() {
			return burdens;
		}

		public void setSignificance(Significance[] significance) {
			this.significance = significance;
		}

		public double getPvalCutoff() {
			return pvalCutoff;
		}

		public double getExcludeFreqBelow() {
			return excludeFreqBelow;
		}

		public String getDirectory() {
			return directory;
		}

		public String[] getFiles() {
			return files;
		}

		public Significance[] getSignificance() {
			return significance;
		}

		public Pheno[] getPhenos() {
			return phenos;
		}
	}

	public static class AllSigs implements Serializable {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		private Significance[] sigs;

		public AllSigs(Significance[] sigs) {
			this.sigs = sigs;
		}

		public void serialize(String filename) {
			Files.writeSerial(this, filename);
		}

		public static AllSigs load(String filename, boolean jar) {
			return (AllSigs) Files.readSerial(filename, jar, true);
		}

	}

	// tracks significant results for each phenotype
	public static class Significance implements Serializable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		private String analysisModel;
		private String analysisType;
		private String phenoType;
		private String minPvalLocus;
		private double minPvalue;
		private int numPassingThreshold;
		private int numTests;
		private int numTotal;
		private ArrayList<String> lociTested;
		private ArrayList<Double> lociTestedPvalue;

		public Significance(String analysisModel, String analysisType, String phenoType) {
			this.analysisModel = analysisModel;
			this.analysisType = analysisType;
			this.phenoType = phenoType;
			this.minPvalLocus = "NA";
			this.minPvalue = 1;
			this.numPassingThreshold = 0;
			this.numTests = 0;
			this.numTotal = 0;
			this.lociTested = new ArrayList<String>();
			this.lociTestedPvalue = new ArrayList<Double>();
		}

		public String getMinPvalLocus() {
			return minPvalLocus;
		}

		public void setMinPvalLocus(String minPvalLocus) {
			this.minPvalLocus = minPvalLocus;
		}

		public String getFullAnalysis() {
			return this.analysisModel + "\t" + this.analysisType + "\t" + this.phenoType + "\t" + this.minPvalue + "\t" + this.minPvalLocus + "\t" + this.numPassingThreshold + "\t" + this.numTests + "/" + this.numTotal;
		}

		public String getAllloci() {
			return Array.toStr(this.lociTested.toArray(new String[lociTested.size()]));
		}

		public void addLoci(String loci) {
			lociTested.add(loci);
		}

		public void addPval(double pval) {
			lociTestedPvalue.add(pval);
		}

		public String getAnalysisModel() {
			return analysisModel;
		}

		public String getAnalysisType() {
			return analysisType;
		}

		public String getPhenoType() {
			return phenoType;
		}

		public double getMinPvalue() {
			return minPvalue;
		}

		public void setMinPvalue(double minPvalue) {
			this.minPvalue = minPvalue;
		}

		public int getNumPassingThreshold() {
			return numPassingThreshold;
		}

		public void setNumPassingThreshold(int numPassingThreshold) {
			this.numPassingThreshold = numPassingThreshold;
		}

		public void addNumPassingThreshold() {
			this.numPassingThreshold = numPassingThreshold + 1;
		}

		public void addNumTest() {
			this.numTests = numTests + 1;
		}

		public void addNumTotal() {
			this.numTotal = numTotal + 1;
		}

		public int getNumTotal() {
			return numTotal;
		}

		// pvalue ties count toward total
		public String getTopLociByPvalue(int numberToReturn, double pvalCutoff) {
			String[] sortedLoci = new String[numberToReturn];
			Hashtable<Double, ArrayList<String>> map = new Hashtable<Double, ArrayList<String>>();
			for (int i = 0; i < lociTestedPvalue.size(); i++) {
				if (map.containsKey(lociTestedPvalue.get(i))) {
					map.get(lociTestedPvalue.get(i)).add(lociTested.get(i));
				} else {
					map.put(lociTestedPvalue.get(i), new ArrayList<String>());
					map.get(lociTestedPvalue.get(i)).add(lociTested.get(i));
				}
			}
			double[] sortedPvalys = toDoubleArray(lociTestedPvalue);
			Arrays.sort(sortedPvalys);
			int numReturned = 0;
			for (int i = 0; i < numberToReturn; i++) {
				ArrayList<String> tops = map.get(sortedPvalys[i]);
				for (int j = 0; j < tops.size(); j++) {
					if (numReturned >= numberToReturn || sortedPvalys[j] > pvalCutoff) {
						if (sortedPvalys[j] > pvalCutoff) {
							sortedLoci[numReturned] = "";
						}
						continue;
					} else {
						sortedLoci[numReturned] = tops.get(j);
						numReturned++;
					}
				}
			}
			return Array.toStr(sortedLoci);
		}
	}

	public static class PrepResults {
		private String phenotype;
		private double[] phenoValues;
		Hashtable<String, Integer> tracker;
		boolean error = false;

		public Hashtable<String, Integer> getTracker() {
			return tracker;
		}

		public String getPhenotype() {
			return phenotype;
		}

		public double[] getPhenoValues() {
			return phenoValues;
		}

		public boolean isError() {
			return error;
		}

		public PrepResults(String phenotype, String[] finalIDs, double[][] database) {
			this.phenotype = phenotype;
			this.phenoValues = new double[finalIDs.length];
			this.tracker = new Hashtable<String, Integer>();
			assignPhenos(finalIDs, database);
		}

		public void assignPhenos(String[] finalIDs, double[][] database) {

			for (int i = 0; i < finalIDs.length; i++) {
				if (tracker.containsKey(finalIDs[i])) {
					error = true;
					System.err.println("Error - duplicate ids were provided in the phenotype file being prepared");
					System.exit(1);
				} else {
					// only retreive one phenotype;
					tracker.put(finalIDs[i], i);
					phenoValues[i] = database[i][0];

				}
			}
		}
	}

	private static double[] toDoubleArray(ArrayList<Double> al) {
		double[] d = new double[al.size()];
		for (int i = 0; i < al.size(); i++) {
			d[i] = al.get(i);
		}
		return d;
	}

	private static String[] toStringArray(ArrayList<String> al) {
		String[] s = new String[al.size()];
		for (int i = 0; i < al.size(); i++) {
			s[i] = al.get(i);
		}
		return s;
	}

	private static GenCNV[] runAnalysis(ArrayList<ArrayList<Analysis>> cabinet, int numThreads, Thread[] threads, Logger log) {
		GenCNV[] genCNVs = new GenCNV[numThreads];
		for (int i = 0; i < numThreads; i++) {
			genCNVs[i] = new GenCNV(cabinet.get(i).toArray(new Analysis[cabinet.get(i).size()]), i);
			threads[i] = new Thread(genCNVs[i]);
			threads[i].start();
		}
		checkThreadStatus(numThreads, threads);
		return genCNVs;
	}

	private static ArrayList<ArrayList<Analysis>> getcabinet(ArrayList<Analysis> analyses, int numThreads) {
		ArrayList<ArrayList<Analysis>> cabinet = new ArrayList<ArrayList<Analysis>>();
		for (int i = 0; i < numThreads; i++) {
			cabinet.add(new ArrayList<Analysis>());
		}
		for (int i = 0; i < analyses.size(); i++) {
			cabinet.get(i % numThreads).add(analyses.get(i));
		}
		return cabinet;
	}

	private static void checkThreadStatus(int numThreads, Thread[] threads) {
		boolean complete;
		complete = false;
		while (!complete) {
			complete = true;
			for (int i = 0; i < numThreads; i++) {
				if (threads[i].isAlive()) {
					complete = false;
				}
			}
			if (!complete) {
				try {
					Thread.sleep(1000L);
				} catch (InterruptedException ex) {
				}
			}
		}
	}

	private static String getFileBase(int modelNum, int analyisNum) {
		return ANALYSIS_MODELS[modelNum] + "_" + ANALYSIS_TYPES[analyisNum];
	}

	private static String getOutputDir(String dir, int modelNum, int analyisNum) {
		return dir + ANALYSIS_MODELS[modelNum] + "/" + ANALYSIS_TYPES[analyisNum];
	}

	public static void analyzeALL(String dir, String gPhenoFile, int numThreads, double excludeFreqBelow, double pvalCutoff, Logger log) {
		Thread[] threads = new Thread[numThreads];
		Pheno[] phenos = loadGPHENO(dir + gPhenoFile, log);
		ArrayList<Analysis> analyses = new ArrayList<Analysis>();
		log.report(ext.getTime() + " Info - positions with frequency less than " + excludeFreqBelow + " will be removed");
		log.report(ext.getTime() + " Info - p-value cutoff set to " + pvalCutoff);
		for (int i = 0; i < ANALYSIS_MODELS.length; i++) {
			for (int j = 0; j < ANALYSIS_TYPES.length; j++) {
				String prefix = getFileBase(i, j);
				String resultDir = getOutputDir(dir, i, j) + "/";
				String[] files = Files.list(resultDir, prefix, null, true, false);
				Significance[] significance = new Significance[phenos.length];
				Burden[] burdens = new Burden[phenos.length];
				for (int k = 0; k < phenos.length; k++) {
					significance[k] = new Significance(ANALYSIS_MODELS[i], ANALYSIS_TYPES[j], phenos[k].getPhenoName());
					burdens[k] = new Burden(ANALYSIS_MODELS[i], ANALYSIS_TYPES[j], phenos[k].getPhenoName(), phenos[k].getArrayInds().length);
				}
				analyses.add(new Analysis(resultDir, files, significance, phenos, burdens, excludeFreqBelow, pvalCutoff));
			}
		}
		ArrayList<ArrayList<Analysis>> cabinet = getcabinet(analyses, numThreads);
		GenCNV[] genCNVs = runAnalysis(cabinet, numThreads, threads, log);
		summerizeAll(genCNVs, dir, gPhenoFile, log);
	}

	private static void summerizeAll(GenCNV[] genCNVs, String dir, String gPhenoFile, Logger log) {
		log.report(Array.toStr(ANALYSIS_SUMMARY_HEADER));
		ArrayList<Significance> allSigs = new ArrayList<Significance>();
		for (int i = 0; i < genCNVs.length; i++) {
			Analysis[] analyses = genCNVs[i].getAnalyses();
			for (int k = 0; k < analyses.length; k++) {
				Significance[] significances = analyses[k].getSignificance();
				Burden[] burdens = analyses[k].getBurdens();
				for (int j = 0; j < significances.length; j++) {
					allSigs.add(significances[j]);
					log.report(significances[j].getFullAnalysis() + "\t" + significances[j].getTopLociByPvalue(DEFAULT_LOCI_TO_RETURN, analyses[k].getPvalCutoff()));
					log.report("BURDEN_" + burdens[j].getAnalysisModel() + "\t" + burdens[j].getAnalysisType() + "\t" + burdens[j].getPheno() + "\t" + burdens[j].getBurdenPvalue());
				}
			}
		}
		String output = dir + ext.rootOf(gPhenoFile) + "_results.ser";
		if (Files.exists(output)) {
			Files.backup(ext.rootOf(gPhenoFile) + "_results.ser", dir, dir + DEFAULT_BACKUP);
		}
		new AllSigs(allSigs.toArray(new Significance[allSigs.size()])).serialize(output);

	}

	public static void runALL(String dir, String cnvFilename, String pedFilename, int numMarkersPerFile, Logger log) {
		for (int i = 0; i < ANALYSIS_MODELS.length; i++) {
			for (int j = 0; j < ANALYSIS_TYPES.length; j++) {
				String prefix = getFileBase(i, j);
				String outputRoot = getOutputDir(dir, i, j) + "/" + prefix;
				ExportCNVsToPedFormat.export(dir + cnvFilename, dir + pedFilename, outputRoot, "\n", false, ANALYSIS_MODEL_PARAMS[i][0], ANALYSIS_MODEL_PARAMS[i][1], ANALYSIS_MODEL_PARAMS[i][2], ANALYSIS_MODEL_PARAMS[i][3], ANALYSIS_MODEL_PARAMS[i][4], true, numMarkersPerFile, DEFAULT_WINDOWS[j]);
			}
		}
	}

	// TODO will make this more modular, but now just setting a specific analysis Type;
	// basically run the prep with each phenotype and print to one file;
	public static String prepPhenos(String dir, String gPhenoFIle, String idFile, String[] covars, Logger log) {
		String newGPhenoFile = ext.rootOf(gPhenoFIle) + ".gprep";
		log.report(dir + gPhenoFIle);
		Pheno[] phenos = loadGPHENO(dir + gPhenoFIle, log);

		String[] uniqInds = getUniqInds(phenos);
		Hashtable<String, String> hashcovars = defineCovars(covars, log);
		ArrayList<PrepResults> prepResults = new ArrayList<PrepResults>();
		int phenoNumber = 0;
		for (int i = 0; i < phenos.length; i++) {
			// log.report("" + phenos[i].getArrayInds().length + "\t" + phenos[i].getPhenoName());
			if (hashcovars.containsKey(phenos[i].getPhenoName()) || !hasVariance(phenos[i])) {
				if (!hasVariance(phenos[i])) {
					log.report("Warning - no variance detected in phenotype " + phenos[i].getPhenoName() + ", removing from analysis");
				}
				continue;
			} else {
				PhenoPrep prep = new PhenoPrep(dir + gPhenoFIle, idFile == null ? null : dir + idFile, GPHENO_HEADERS[0], phenos[i].getPhenoName(), covars, log);
				prep.computeResiduals();
				prep.inverseNormalize();
				prepResults.add(new PrepResults(phenos[i].getPhenoName(), prep.getFinalIDs(), prep.getDatabase()));
				// log.report("" + prep.getFinalIDs().length + "\t" + phenos[i].getPhenoName());
				phenoNumber++;
			}
		}
		printNewGPheno(dir, newGPhenoFile, prepResults.toArray(new PrepResults[prepResults.size()]), uniqInds, log);
		return newGPhenoFile;
	}

	private static boolean hasVariance(Pheno pheno) {
		Hashtable<Double, Boolean> tracker = new Hashtable<Double, Boolean>();
		double[] phenos = pheno.getArrayPheno();
		for (int i = 0; i < phenos.length; i++) {
			tracker.put(phenos[i], true);
			if (tracker.size() > 1) {
				return true;
			}
		}
		return tracker.size() > 1;
	}

	private static void printNewGPheno(String dir, String newGPhenoFile, PrepResults[] prepResults, String[] uniqInds, Logger log) {
		String output = dir + newGPhenoFile;
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			writer.print("FID\tIID");
			for (int i = 0; i < prepResults.length; i++) {
				writer.print("\t" + prepResults[i].getPhenotype());
			}
			writer.println();
			for (int i = 0; i < uniqInds.length; i++) {
				String indKey = uniqInds[i].split(ID_DELIMITER)[0];
				writer.print(Array.toStr(uniqInds[i].split(ID_DELIMITER)));
				for (int j = 0; j < prepResults.length; j++) {
					if (prepResults[j].getTracker().containsKey(indKey)) {
						writer.print("\t" + prepResults[j].getPhenoValues()[prepResults[j].getTracker().get(indKey)]);
					} else {
						writer.print("\t.");
					}
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + output);
			e.printStackTrace();
		}
		log.report(ext.getTime() + " Finished prepping " + prepResults.length + " phenotypes for " + uniqInds.length + " individuals");

	}

	private static String[] getUniqInds(Pheno[] phenos) {
		ArrayList<String> inds = new ArrayList<String>();
		Hashtable<String, Boolean> tracker = new Hashtable<String, Boolean>();
		for (int i = 0; i < phenos.length; i++) {
			String[] phenoinds = phenos[i].getArrayInds();
			for (int j = 0; j < phenoinds.length; j++) {
				if (tracker.containsKey(phenoinds[j])) {
					continue;
				} else {
					inds.add(phenoinds[j]);
					tracker.put(phenoinds[j], true);
				}
			}
		}
		return inds.toArray(new String[inds.size()]);
	}

	// for (int j = 0; j < finalIDs.length; j++) {
	// if (hash.containsKey(finalIDs[j])) {
	// writer.println(hash.get(finalIDs[j]) + "\t" + finalIDs[j] + "\t" + database[j][0]);
	// }
	// for (int m = 0; m < finalIDs.length; m++) {
	// writer.println(finalIDs[m] + delimiter + Array.toStr(database[m], -1, -1, delimiter));
	// }

	public static Hashtable<String, String> defineCovars(String[] covars, Logger log) {
		Hashtable<String, String> hash = new Hashtable<String, String>();
		for (int i = 0; i < covars.length; i++) {
			if (hash.containsKey(covars[i])) {
				log.reportError("Error - covariate " + covars[i] + " was defined twice");
			} else {
				hash.put(covars[i], covars[i]);
			}
		}
		return hash;
	}

	// public static void export(String cnvFilename, String pedFilename, String outputRoot, String endOfLine, boolean rfglsOutput, boolean includeDele, boolean includeDupl, boolean ordered, boolean collapsed, boolean homozygous, boolean excludeMonomorphicLoci, int markersPerFile, int win) {
	// Pheno[] phenos = loadGPHENO(dir + gPhenoFIle, log);
	// ExportCNVsToPedFormat.export(dir + cnvFile, dir + pedFilename, dir + ANALYSIS_MODELS[0] + "/" + ANALYSIS_MODELS[0], "\n", false, true, true, false, false, false, false, 5000, 0);

	public static void main(String[] args) {
		String dir = "C:/data/ARIC/ARIC_Pheno/PheWas/PheWASW/";
		// String cnvFile = "low.cnv";
		String pedFilename = "pedW.gped";
		String cnvFile = "low_probRemoved.cnv";

		// String gPhenoFIle = "phenoB.gfam";
		String gSummary = "summary.txt";
		// String gPhenoFIle = "phenoB_test.gfam";
		String gPhenoFIle = "phenoW.gpheno";
		String covars = "Male,Age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20";
		int numThreads = 8;
		int numMarkersPerFile = 5000;
		double excludeFreqBelow = 0.01;
		double pvalCutoff = 0.0001;
		boolean prepPhenos = true;
		// int window =

		if (Files.exists(dir + DEFAULT_LOG_NAME)) {
			Files.backup(DEFAULT_LOG_NAME, dir, dir + DEFAULT_BACKUP);
		}

		Logger log = new Logger(dir + DEFAULT_LOG_NAME);

		gPhenoFIle = prepPhenos(dir, gPhenoFIle, null, covars.split(","), log);
		log.report(ext.getTime() + " Generating CNV ped files ");
		// runALL(dir, cnvFile, pedFilename, numMarkersPerFile, log);
		log.report(ext.getTime() + " Finished generating CNV ped files ");
		log.report(ext.getTime() + " Analyzing files CNV ped files");
		analyzeALL(dir, gPhenoFIle, numThreads, excludeFreqBelow, pvalCutoff, log);
		log.report(ext.getTime() + " Finished analyzing  CNV ped files");

		// ExportCNVsToPedFormat.export(dir + cnvFile, dir + pedFilename, dir + ANALYSIS_MODELS[0] + "/" + ANALYSIS_MODELS[0], "\n", false, true, true, false, false, false, false, 5000, 0);
		// parseCNVGenotypes(dir + cnvFile, dir + pedFilename, dir + ANALYSIS_MODELS[0] + "/" + ANALYSIS_MODELS[0], "\n", false, true, true, false, false, false, false, 5000, 0, phenos);
		// runLinearRegression(dir + "ANY_CNV_COLLAPSED/" + testgpedFile, phenos, log);
		// public static void export(String cnvFilename, String pedFilename, String outputRoot, String endOfLine, boolean rfglsOutput, boolean includeDele, boolean includeDupl, boolean ordered, boolean collapsed, boolean homozygous, boolean excludeMonomorphicLoci, int markersPerFile, int win) {
		// Pheno[] phenos = loadGPHENO(dir + gPhenoFIle, log);
		// ExportCNVsToPedFormat.export(dir + cnvFile, dir + pedFilename, dir + ANALYSIS_MODELS[0] + "/" + ANALYSIS_MODELS[0], "\n", false, true, true, false, false, false, false, 5000, 0);
	}
}
