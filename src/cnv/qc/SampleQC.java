package cnv.qc;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.Hashtable;

import stats.Quantiles;
import common.Array;
import common.Files;
import common.HashVec;
import common.Sort;
import common.ext;
import cnv.analysis.pca.PrincipalComponentsResiduals;
import cnv.filesys.Project;
import cnv.var.SampleData;

/**
 * Class that automates adding sample QC metrics to sample data, and also parses the qc metrics to the desired quantile class
 * <p>
 * Currently requires that all samples are present
 *
 */
public class SampleQC {
	private Project proj;
	private String[] samples;
	private boolean[] excludes;
	private String[] excludeNotes;
	private double[][] qcMatrix;
	private String[] qctitles;
	private int numAdded;

	public SampleQC(Project proj, String[] qctitles) {
		super();
		this.proj = proj;
		this.samples = proj.getSamples();
		this.excludes = new boolean[samples.length];
		Arrays.fill(this.excludes, false);
		this.excludeNotes = new String[samples.length];
		Arrays.fill(this.excludeNotes, ".");
		this.qcMatrix = new double[qctitles.length][samples.length];
		this.qctitles = qctitles;
		this.numAdded = 0;
	}
	
	/**
	 * @param base existing SampleQC to start from
	 * @param useCols boolean array of cols to keep from base
	 */
	public SampleQC(SampleQC base, boolean[] useCols){
		super();
		this.proj = base.proj;
		this.samples = base.samples;
		if (!base.verify()) {
			proj.getLog().reportTimeError("Could not verify that all data was added");
		}
		this.excludes = base.excludes;
		this.excludeNotes = base.excludeNotes;
		this.qcMatrix = Array.subArray(base.getQcMatrix(), useCols);
		this.qctitles = Array.subArray(base.getQctitles(), useCols);
		this.numAdded = qctitles.length * samples.length;
	}

	public String[] getSamples() {
		return samples;
	}

	public String[] getQctitles() {
		return qctitles;
	}

	public double[][] getQcMatrix() {
		return qcMatrix;
	}

	public boolean verify() {
		return numAdded == qctitles.length * samples.length;
	}

	public double[] getDataFor(String qctitle) {
		int indexToExtract = ext.indexOfStr(qctitle, qctitles);
		if (indexToExtract < 0) {
			proj.getLog().reportTimeError("Invalid title " + qctitle + " returning null");
			return null;
		} else {
			return qcMatrix[indexToExtract];
		}

	}

	private void addToMatrix(String sample, int qcTitleIndex, double data) {
		int sampleIndex = ext.indexOfStr(sample, samples);
		addToMatrix(sampleIndex, qcTitleIndex, data);
	}

	private void addToMatrix(int sampleIndex, int qcTitleIndex, double data) {
		qcMatrix[qcTitleIndex][sampleIndex] = data;
		numAdded++;
		if (numAdded > qctitles.length * samples.length) {
			proj.getLog().reportTimeError("Internal Error; too many QC metrics have been added");
		}
	}
	
	private boolean addToExclude(int sampleIndex, String reason){
		if (!excludes[sampleIndex]){
			excludes[sampleIndex] = true;
			excludeNotes[sampleIndex] = reason;
			return true;
		}
		excludeNotes[sampleIndex] += ("; " + reason);
		return false;
	}
	
	/**
	 * Prepares a hash to be used in the addition to sample data
	 */
	private Hashtable<String, String> developQCHash(Quantiles[] quantiles, boolean justClasses) {
		Hashtable<String, String> hashtable = new Hashtable<String, String>();
		for (int i = 0; i < samples.length; i++) {
			String qcInfo = "";
			for (int j = 0; j < quantiles.length; j++) {
				qcInfo += (j == 0 ? "" : "\t");
				if (!justClasses) {
					qcInfo += qcMatrix[j][i] + "\t";
				}
				qcInfo += quantiles[j].getQuantileMembershipAsRoundedInt()[i];
			}
			qcInfo += ("\t" + (excludes[i] ? "1" : "0") + "\t" + excludeNotes[i]);
			hashtable.put(samples[i], qcInfo);
		}
		return hashtable;
	}
	
	private String[] developQCHeader(Quantiles[] quantiles, int numQ, boolean justClasses) {
		return Array.combine(developHeader(quantiles, qctitles, numQ, justClasses), new String[] {"CLASS=Exclude","ExcludeNote"});
	}

	public void addQCsToSampleData(int numQ, boolean justClasses) {
		Quantiles[] quantiles = Quantiles.qetQuantilesFor(numQ, qcMatrix, qctitles, proj.getLog());
		Hashtable<String, String> hashtable = developQCHash(quantiles, justClasses);
		String[] header = developQCHeader(quantiles, numQ, justClasses);
		addToSampleData(proj, hashtable, header, numQ, justClasses);
	}

	public void addPCsToSampleData(int numQ, int numPCs, boolean justClasses) {
		proj.getLog().reportTimeInfo("Adding " + numPCs + " to sample data");
		PrincipalComponentsResiduals pcResiduals = proj.loadPcResids();
		double[][] pcBasisSubset = new double[numPCs][];
		String[] pcTitles = new String[numPCs];
		if (pcResiduals.isSortedByProject()) {
			for (int i = 0; i < pcTitles.length; i++) {
				pcTitles[i] = "PC" + (i + 1);
				pcBasisSubset[i] = pcResiduals.getBasisAt((i + 1));
			}
			Quantiles[] quantiles = Quantiles.qetQuantilesFor(numQ, pcBasisSubset, pcTitles, proj.getLog());
			Hashtable<String, String> hashtable = developHash(quantiles, proj.getSamples(), pcBasisSubset, justClasses);
			String[] header = developHeader(quantiles, pcTitles, numQ, justClasses);
			addToSampleData(proj, hashtable, header, numQ, justClasses);
		} else {
			proj.getLog().reportTimeError("PCs are not sorted by project, currently this is not supported");
		}
	}

	private static void addToSampleData(Project proj, Hashtable<String, String> hashtable, String[] header, int numQ, boolean justClasses) {
		SampleData sampledata = proj.getSampleData(0, false);
		proj.getLog().reportTimeInfo("Adding " + header.length + " columns to sample data based on sample QC");
		sampledata.addData(hashtable, "DNA", header, "NaN", "\t", proj.getLog());
	}

	/**
	 * Prepares a hash to be used in the addition to sample data
	 */
	private static Hashtable<String, String> developHash(Quantiles[] quantiles, String[] samples, double[][] variableDomMatrix, boolean justClasses) {
		Hashtable<String, String> hashtable = new Hashtable<String, String>();
		for (int i = 0; i < samples.length; i++) {
			String qcInfo = "";
			for (int j = 0; j < quantiles.length; j++) {
				qcInfo += (j == 0 ? "" : "\t");
				if (!justClasses) {
					qcInfo += variableDomMatrix[j][i] + "\t";
				}
				qcInfo += quantiles[j].getQuantileMembershipAsRoundedInt()[i];
			}
			hashtable.put(samples[i], qcInfo);
		}
		return hashtable;
	}

	private static String[] developHeader(Quantiles[] quantiles, String[] titles, int numQ, boolean justClasses) {
		String[] qcHeader = new String[justClasses ? quantiles.length : quantiles.length * 2];
		int curIndex = 0;
		for (int i = 0; i < quantiles.length; i++) {
			if (!justClasses) {
				qcHeader[curIndex] = titles[i];
				curIndex++;
			}
			qcHeader[curIndex] = getClassForQuantile(quantiles[i], titles[i], numQ);
			curIndex++;
		}
		return qcHeader;
	}

	/**
	 * Creates sample data friendly class titles for column headers
	 */
	private static String getClassForQuantile(Quantiles quantiles, String qcTitle, int numQ) {
		String thisClass = "CLASS=QUANTILE_" + numQ + "_" + qcTitle;
		int[] uniqLabels = Array.toIntArray(Array.unique(Array.toStringArray(quantiles.getQuantileMembershipAsRoundedInt())));
		double[] uniqQs = Array.toDoubleArray(Array.unique(Array.toStringArray(quantiles.getQuantileMembership())));
		int[] orderLabels = Sort.quicksort(uniqLabels);
		int[] orderQs = Sort.quicksort(uniqQs);
		for (int i = 0; i < orderLabels.length; i++) {
			thisClass += ";" + uniqLabels[orderLabels[i]] + "=q_" + ext.formDeci(uniqQs[orderQs[i]], 3);
		}
		return thisClass;
	}

	private static boolean verifyAllProjectSamples(Project proj, String lrrSdToLoad, int sampleColumn) {
		String[] projSamples = proj.getSamples();
		String[] fileSamples = HashVec.loadFileToStringArray(lrrSdToLoad, false, new int[] { sampleColumn }, false);
		int[] indices = ext.indexLargeFactors(fileSamples, projSamples, true, proj.getLog(), false, false);
		return fileSamples.length - Array.countIf(indices, -1) == projSamples.length;
	}

	public static SampleQC loadSampleQC(Project proj) {
		return loadSampleQC(proj, LrrSd.SAMPLE_COLUMN, LrrSd.NUMERIC_COLUMNS, false);
	}

	public static SampleQC loadSampleQC(Project proj, boolean generate) {
		return loadSampleQC(proj, LrrSd.SAMPLE_COLUMN, LrrSd.NUMERIC_COLUMNS, generate);
	}

	/**
	 * @param proj
	 * @param sampleColumnName
	 *            header of the column containing sample names
	 * @param qcTitlesToLoad
	 *            qc titles to load from the sample QC file
	 * @param generate
	 *            generate sampleQC if missing
	 * @return
	 */
	public static SampleQC loadSampleQC(Project proj, String sampleColumnName, String[] qcTitlesToLoad, boolean generate) {
		// String lrrSdToLoad = proj.getFilename(proj.SAMPLE_QC_FILENAME);
		String lrrSdToLoad = proj.SAMPLE_QC_FILENAME.getValue();
		SampleQC sampleQC = null;
		if (!Files.exists(lrrSdToLoad) && !generate) {
			proj.getLog().reportTimeError("Could not find sample QC file " + lrrSdToLoad);

		}

		else {
			if (!Files.exists(lrrSdToLoad) && !generate) {
				proj.getLog().reportTimeInfo("Attempting to generate QC file " + lrrSdToLoad);
				LrrSd.init(proj, null, null, proj.NUM_THREADS.getValue());
				if (!Files.exists(lrrSdToLoad)) {
					proj.getLog().reportTimeError("Could not find generate QC file " + lrrSdToLoad);
				}
			}
			proj.getLog().reportTimeInfo("Loading qc data from " + lrrSdToLoad);
			try {
				BufferedReader reader = Files.getAppropriateReader(lrrSdToLoad);
				String[] header = reader.readLine().trim().split("[\\s]+");
				int[] indicesToLoad = ext.indexFactors(qcTitlesToLoad, header, true, proj.getLog(), true, false);
				int sampleColumn = ext.indexOfStr(sampleColumnName, header);
				if (Array.countIf(indicesToLoad, -1) > 0 || sampleColumn < 0) {
					proj.getLog().reportTimeError("Could not find all desired columns in qc file " + lrrSdToLoad);
					proj.getLog().reportTimeError("Consider re-creating " + lrrSdToLoad + " if sample qc has been updated");

				} else if (!verifyAllProjectSamples(proj, lrrSdToLoad, sampleColumn)) {
					proj.getLog().reportTimeError("Could not find all of the projects samples in qc file " + lrrSdToLoad);
				} else {
					sampleQC = new SampleQC(proj, qcTitlesToLoad);

					while (reader.ready()) {
						String[] line = reader.readLine().trim().split("[\\s]+");
						String sample = line[sampleColumn];
						if (ext.indexOfStr(sampleColumnName, line) < 0) {
							for (int i = 0; i < indicesToLoad.length; i++) {
								double data = Double.NaN;
								try {
									data = Double.parseDouble(line[indicesToLoad[i]]);
								} catch (NumberFormatException e) {
									proj.getLog().reportTimeWarning("line " + Array.toStr(line) + " contained an invalid number for qc column " + qcTitlesToLoad[i]);
								}
								sampleQC.addToMatrix(sample, i, data);
							}
						}
					}
					if (!sampleQC.verify()) {
						proj.getLog().reportTimeError("Could not verify that all data has been added");
						sampleQC = null;
					}
				}
				reader.close();
				proj.getLog().reportTimeInfo("Finished loading qc data from " + lrrSdToLoad);

			} catch (FileNotFoundException fnfe) {
				proj.getLog().reportTimeError("file \"" + lrrSdToLoad + "\" not found in current directory");
				return null;
			} catch (IOException ioe) {
				proj.getLog().reportTimeError("Error reading file \"" + lrrSdToLoad + "\"");
				return null;
			}

			if (sampleQC != null) {
				proj.getLog().reportTimeInfo("Filtering empty columns from " + lrrSdToLoad);

				boolean[] useCols = new boolean[sampleQC.getQcMatrix().length];
				Arrays.fill(useCols, false);
				int numFiltered = sampleQC.getQcMatrix().length;
				for (int i = 0; i < sampleQC.getQcMatrix().length; i++) {
					for (int j = 0; j < sampleQC.getQcMatrix()[i].length; j++){
						if (!Double.isNaN(sampleQC.getQcMatrix()[i][j])) {
							useCols[i] = true;
							numFiltered--;
							break;
						}
					}
				}
				sampleQC = new SampleQC(sampleQC, useCols);

				proj.getLog().reportTimeInfo("Filtered " + numFiltered + " empty columns from " + lrrSdToLoad);
				
				proj.getLog().reportTimeInfo("Finding samples to exclude");
				
				int numExcluded = 0;
				int lrr_sdIndex = ext.indexOfStr("LRR_SD_Post_Correction", sampleQC.getQctitles());
				if (lrr_sdIndex == -1) lrr_sdIndex = ext.indexOfStr("LRR_SD", sampleQC.getQctitles());
				int callrateIndex = ext.indexOfStr("Genotype_callrate", sampleQC.getQctitles());
				
				if (lrr_sdIndex != -1){
					for (int i = 0; i < sampleQC.getQcMatrix()[lrr_sdIndex].length; i++){
						if (sampleQC.getQcMatrix()[lrr_sdIndex][i] > proj.LRRSD_CUTOFF.getValue()){
							if (sampleQC.addToExclude(i,"LRR_SD > " + proj.LRRSD_CUTOFF.getValue())) numExcluded++;
						}
					}
				}
				if (callrateIndex != -1){
					for (int i = 0; i < sampleQC.getQcMatrix()[callrateIndex].length; i++){
						if (sampleQC.getQcMatrix()[callrateIndex][i] < proj.SAMPLE_CALLRATE_THRESHOLD.getValue()){
							if (sampleQC.addToExclude(i,"Callrate < " + proj.SAMPLE_CALLRATE_THRESHOLD.getValue())) numExcluded++;
						}
					}
				}
				
				proj.getLog().reportTimeInfo("Found " + numExcluded + " samples to exclude");
			}
			
		}
		return sampleQC;
	}

	public static void parseAndAddToSampleData(Project proj, int numQ, int numPCs, boolean justClasses) {
		SampleQC sampleQC = loadSampleQC(proj, false);
		sampleQC.addQCsToSampleData(numQ, justClasses);
		if (numPCs > 0) {
			sampleQC.addPCsToSampleData(numQ, numPCs, justClasses);
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		int numQ = 5;
		int numPCs = 0;
		boolean justClasses = true;

		String usage = "\n" + "cnv.qc.SampleQC requires 0-1 arguments\n";
		usage += "   (1) filename (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) number of quantiles to divide the sample QC metrics into (5 = quintiles, 100 = percentiles) (i.e. numQ=" + numQ + " (default))\n" + "";
		usage += "   (3) add only class (quantiled) qc data to sample data  (i.e. justClasses=" + justClasses + " (default))\n" + "";
		usage += "   (4) if a pc file is available, add this many pcs to the sample data file , must be set to >=1 to be added  (i.e. numPCs=" + numPCs + " (default,no addition))\n" + "";

		usage += "   NOTE: the projects sample qc file must be present, for the qc metrics to be added to sample data";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("numQ=")) {
				numQ = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("numPCs=")) {
				numPCs = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("justClasses=")) {
				justClasses = ext.parseBooleanArg(args[i]);
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
			Project proj = new Project(filename, false);
			parseAndAddToSampleData(proj, numQ, numPCs, justClasses);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
