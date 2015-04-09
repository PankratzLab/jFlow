package cnv.qc;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
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
	private double[][] qcMatrix;
	private String[] qctitles;
	private int numAdded;

	public SampleQC(Project proj, String[] qctitles) {
		super();
		this.proj = proj;
		this.samples = proj.getSamples();
		this.qcMatrix = new double[qctitles.length][samples.length];
		this.qctitles = qctitles;
		this.numAdded = 0;
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

	public void addQCsToSampleData(int numQ, boolean justClasses) {
		addToSampleData(proj, qcMatrix, qctitles, numQ, justClasses);
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
			addToSampleData(proj, pcBasisSubset, pcTitles, numQ, justClasses);
		} else {
			proj.getLog().reportTimeError("PCs are not sorted by project, currently this is not supported");
		}
	}

	private static void addToSampleData(Project proj, double[][] variableDomMatrix, String[] titles, int numQ, boolean justClasses) {
		SampleData sampledata = proj.getSampleData(0, false);
		Quantiles[] quantiles = Quantiles.qetQuantilesFor(numQ, variableDomMatrix, titles, proj.getLog());
		Hashtable<String, String> hashtable = developQCHash(quantiles, proj.getSamples(), variableDomMatrix, justClasses);
		String[] header = developQCHeader(quantiles, titles, numQ, justClasses);
		proj.getLog().reportTimeInfo("Adding " + titles.length + " sample qc metric(s) to sample data");
		sampledata.addData(hashtable, "DNA", header, "NaN", "\t", proj.getLog());
	}

	/**
	 * Prepares a hash to be used in the addition to sample data
	 */
	private static Hashtable<String, String> developQCHash(Quantiles[] quantiles, String[] samples, double[][] variableDomMatrix, boolean justClasses) {
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

	private static String[] developQCHeader(Quantiles[] quantiles, String[] titles, int numQ, boolean justClasses) {
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
		return loadSampleQC(proj, LrrSd.SAMPLE_COLUMN, LrrSd.NUMERIC_COLUMNS);
	}

	public static SampleQC loadSampleQC(Project proj, String sampleColumnName, String[] qcTitlesToLoad) {
		String lrrSdToLoad = proj.getFilename(proj.SAMPLE_QC_FILENAME);
		SampleQC sampleQC = null;
		if (!Files.exists(lrrSdToLoad)) {
			proj.getLog().reportTimeError("Could not find sample QC file " + lrrSdToLoad);

		} else {
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
									proj.getLog().reportTimeWarning("line " + Array.toStr(line) + " contained and invalid number for qc column " + qcTitlesToLoad[i]);
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
		}
		return sampleQC;
	}

	public static void parseAndAddToSampleData(Project proj, int numQ, int numPCs, boolean justClasses) {
		SampleQC sampleQC = loadSampleQC(proj, LrrSd.SAMPLE_COLUMN, LrrSd.NUMERIC_COLUMNS);
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
