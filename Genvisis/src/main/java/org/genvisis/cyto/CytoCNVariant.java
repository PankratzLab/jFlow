package org.genvisis.cyto;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;

/**
 * Class to store cyto variants and related info. Currently variants are loaded from files with header containing the following:
 * <p>
 * "AberrationNo	CytoBand	ChrName	ProbeName	Start	Stop	Description	Genes		Logratio	Amplification	Deletion		Logratio	Amplification	Deletion		Logratio	Amplification	Deletion		Logratio	Amplification	Deletion" but only CYTO_ABERATION_HEADER headers are needed
 * <p>
 * above each Logratio column should be a sample name
 * <p>
 * Data can look like:
 * <p>
 * 1.1 p36.11 chr1 A_18_P10061443 25599273 25599326 ref|Homo sapiens Rh blood group, D antigen (RHD), transcript variant 2, mRNA. RHD 0.31111443 0.45292675 -0.46087968 0 -0.096791215 0 0.2777727 0
 */
public class CytoCNVariant extends CNVariant {

	public static final String[][] CYTO_ABERATION_HEADER = { { "AberrationNo", "CytoBand", "ChrName", "ProbeName", "Start", "Stop", "Genes" } };
	public static final String[] CYTO_SAMPLE_HEADER = { "Logratio", "Amplification", "Deletion" };
	public static final String[] SPLITS = { "\t" };
	public static final String[] NOTHINGS = { "0", "0.0", "" };
	public static final String[] INTERPS = { "loss in ", "gain in " };
	public static final String GENES_DELIM = ",";
	public static final String GENES_BLANK = "---";
	public static final String CNV_EXT = ".cnv";
	public static final int DEFAULT_SCORE = 127;

	public static final int[] CN_NUMBER = { 1, 3 };

	private static final long serialVersionUID = 1L;
	private String cytoBand;
	private String[] probeNames, genes;
	private ArrayList<String> probeNamestmp, genestmp, cytoBandtmp;
	private double avgLogRatio;

	/**
	 * Constructor when all the pieces are known
	 */
	public CytoCNVariant(String familyID, String individualID, byte chr, int start, int stop, int cn, double score, int numMarkers, int source, String cytoBand, String[] probeNames, String[] genes, double avgLogratio) {
		super(familyID, individualID, chr, start, stop, cn, score, numMarkers, source);
		this.cytoBand = cytoBand;
		this.probeNames = probeNames;
		this.genes = genes;
		this.avgLogRatio = avgLogratio;
	}

	/**
	 * Constructor when things may need to updated, i.e scanning and collecting from an input file with multiple samples per line, call finalize to erase the ArrayLists and convert to Array[]
	 */
	public CytoCNVariant(String familyID, String individualID, byte chr, int start, int stop, int cn, double score, int numMarkers, int source, String cytoBand, double avgLogratio, String firstProbe, String firstGene) {
		super(familyID, individualID, chr, start, stop, cn, score, numMarkers, source);
		this.cytoBand = cytoBand;
		this.probeNamestmp = new ArrayList<String>();
		this.genestmp = new ArrayList<String>();
		this.cytoBandtmp = new ArrayList<String>();
		this.avgLogRatio = avgLogratio;
		probeNamestmp.add(firstProbe);
		addGene(firstGene);
	}

	/**
	 * Returns ISCN format such as "p22.33(1,770,348-1,783,772)x1" or "p11.2(57,189,194-57,596,150)x3"
	 * <p>
	 * x1 = deletion
	 * <p>
	 * x3 = amplification
	 */
	public String getISCN() {
		String ISCN = "";
		ISCN += getCytoBand() + "(" + ext.addCommas(getStart()) + "-" + ext.addCommas(getStop()) + ")x" + getCN();
		return ISCN;
	}

	public String[] getGenes() {
		return genes;
	}

	/**
	 * The interpretation such as "loss in p22.33 (13.4 kb)" or "gain in q11.2 (842.0 kb)"
	 * <p>
	 * Warning - returns "" if invalid copy number, must be either 1 or 3
	 */
	public String getInterpretation() {
		String interp = "";
		if (getCN() == 1) {
			interp += INTERPS[0];
		} else if (getCN() == 3) {
			interp += INTERPS[1];
		} else {
			System.err.println("Error - invalid copy number " + getCN() + " for CytoCNVariants, copy number must be " + Array.toStr(CN_NUMBER, ", or"));
			return "";
		}
		interp += getCytoBand() + " (" + ext.prettyUpDistance(getSize(), 1) + ")";
		return interp;
	}

	public double getAvgLogRatio() {
		return avgLogRatio;
	}

	public String getCytoBand() {
		return cytoBand;
	}

	/**
	 * @return a String representing the genes of the aberration, GENES_DELIM delimited
	 */
	public String getStringGenes() {
		return (genes.length > 0 ? Array.toStr(genes, GENES_DELIM) : GENES_BLANK);
	}

	/**
	 * Compares a new Average Log Ratio to the existing one for the variant
	 * <p>
	 * We assume that a different Log Ratio indicates a different variant
	 */
	public boolean isSame(double newAvgLogratio) {
		return newAvgLogratio == avgLogRatio;
	}

	/**
	 * @return a summary of the region (as requested by the cytogenics lab)
	 */
	public String getMyReport() {
		return getUCSClocation() + "\t" + getUCSCLink("hg18") + "\t" + getISCN() + "\t" + getInterpretation() + "\t" + getStringGenes() + "\t" + getAvgLogRatio() + "\t" + getNumMarkers();
	}

	/**
	 * A method that will update a CytoCNVariants as we scan through a file
	 */
	public boolean update(String newCytoBand, byte newChr, int newstop, String newProbeName, String newGene, double newAvgLogratio, Logger log) {
		boolean updated = true;
		if (probeNamestmp == null) {
			log.reportError("Error - this method can only be used with the temporary constructor");
			updated = false;
		}
		// TODO, how to properly format this when there are multiple cytobands, need to know for proper ISCN formatting?
		if (!cytoBandtmp.contains(newCytoBand)) {
			cytoBandtmp.add(newCytoBand);
		}
		if (!isSameChr(newChr)) {
			log.reportError("Error - this aberation spans multiple chromosomes, this should not happen");
			updated = false;
		}
		if (newstop <= stop) {
			log.reportError("Error - the current stop being added is less than the previous stop position, this should not happen");
			updated = false;
		}
		if (newstop <= start) {
			log.reportError("Error - the current stop being added is less than the start position, this should not happen");
			updated = false;
		}
		if (!isSame(newAvgLogratio)) {
			log.reportError("Error - the current logRatio does not equal the previous, this should not happen");
			updated = false;
		}
		if (updated) {
			probeNamestmp.add(newProbeName);
			addGene(newGene);
			this.stop = newstop;
		}
		return updated;
	}

	public boolean isSameChr(byte newChr) {
		return chr == newChr;
	}

	/**
	 * We only add the gene if it is non-missing, and not already added.
	 * <p>
	 * TODO the test for existence could be sped up, but likely not worth it
	 */
	public void addGene(String gene) {
		if (ext.indexOfStr(gene, NOTHINGS) == -1 && !genestmp.contains(gene)) {
			genestmp.add(gene);
		}
	}

	/**
	 * Format the cytoband as chromosome(UCSC,not int if appropriate) and band1band2....
	 */
	private void formatCytoBand() {
		this.cytoBand = Positions.getChromosomeUCSC(chr,false) + Array.toStr(cytoBandtmp.toArray(new String[cytoBandtmp.size()]), "");
	}

	/**
	 * Call this method to clear ArrayLists and convert to arrays. Assigns number of probes found. Used for loading from file when unknown number of probes are being added
	 */
	public void finalizeVariant() {
		this.probeNames = probeNamestmp.toArray(new String[probeNamestmp.size()]);
		this.probeNamestmp.clear();
		this.genes = genestmp.toArray(new String[genestmp.size()]);
		genestmp.clear();
		this.numMarkers = probeNames.length;
		formatCytoBand();
		cytoBandtmp.clear();
	}

	/**
	 * Parse a CytoCNVariant[] to a CytoCNVariant[][] of individuals
	 * <p>
	 * Warning - this method assumes that FID\tIID ids are unique to a sample, this is not always the case
	 * <p>
	 * Warning - samples are ordered on a "first-seen basis" when iterating over CytoCNVariant[] cytoCNVariants
	 * 
	 * @param cytoCNVariants
	 *            the CytoCNVariant[] to parse
	 * @param log
	 * @return CytoCNVariant[][] organized as CytoCNVariant[sample0][CytoCNVariantsForSample0]
	 */
	public static CytoCNVariant[][] toIndividuals(CytoCNVariant[] cytoCNVariants, Logger log) {
		Hashtable<String, ArrayList<CytoCNVariant>> track = new Hashtable<String, ArrayList<CytoCNVariant>>();
		ArrayList<String> inds = new ArrayList<String>();

		for (int i = 0; i < cytoCNVariants.length; i++) {
			String key = cytoCNVariants[i].getFamilyID() + "\t" + cytoCNVariants[i].getIndividualID();
			if (!track.containsKey(key)) {
				track.put(key, new ArrayList<CytoCNVariant>());
				inds.add(key);
			}
			track.get(key).add(cytoCNVariants[i]);
		}

		CytoCNVariant[][] cytoCNVariantsInds = new CytoCNVariant[inds.size()][];
		for (int i = 0; i < inds.size(); i++) {
			cytoCNVariantsInds[i] = track.get(inds.get(i)).toArray(new CytoCNVariant[track.get(inds.get(i)).size()]);
		}
		return cytoCNVariantsInds;
	}

	/**
	 * Load a cyto aberration file to a CytoCNVariant[]
	 */
	public static CytoCNVariant[] loadCytoCNVariant(String cytoCNVariantFile, Logger log) {
		String[] line;
		String[] header, sampleNames;
		// the common indices for all samples on each line in the file
		int[] indicesCommon;
		// the unique indices for the samples in the file
		int[][] indicesSamples;
		// temporary storage for the parsed (finalized) aberrations as we move through the file
		ArrayList<CytoCNVariant> cytoCNVariant = new ArrayList<CytoCNVariant>();
		int scan = 0;
		int count = 0;
		// temporary storage for aberrations that are currently being parsed (one per sample)
		CytoCNVariant[] tmps;
		try {
			BufferedReader reader = Files.getAppropriateReader(cytoCNVariantFile);
			do {
				line = reader.readLine().trim().split(SPLITS[0], -1);
				count++;
			} while (reader.ready() && (ext.indexFactors(CYTO_ABERATION_HEADER, line, false, true, false, false)[0] == -1));
			if (!reader.ready()) {
				log.reportError("Error - did not find the neccesary column headers in file " + cytoCNVariantFile);
				return null;
			}
			header = line;
			indicesCommon = ext.indexFactors(CYTO_ABERATION_HEADER[0], header, false, false);
			indicesSamples = determineSampleIndices(header, log);
			reader.close();

			reader = Files.getAppropriateReader(cytoCNVariantFile);
			// the sample names are located right above the data header, so we re-scan the reader to that location
			do {
				line = reader.readLine().trim().split(SPLITS[0], -1);
				scan++;
			} while (scan < count - 1);

			sampleNames = getSamples(line);

			if (indicesSamples.length == sampleNames.length) {
				log.report(ext.getTime() + " Info - found the following samples to parse: " + Array.toStr(sampleNames) + " in file " + cytoCNVariantFile);
			} else {
				log.reportError("Error - data columns were found for " + indicesSamples.length + "samples, but only " + Array.toStr(sampleNames) + " names were found");
				return null;
			}
			// skip the already parsed data header
			reader.readLine();
			tmps = new CytoCNVariant[sampleNames.length];

			while (reader.ready()) {
				line = reader.readLine().trim().split(SPLITS[0], -1);
				parseInputLine(sampleNames, line, indicesCommon, indicesSamples, tmps, cytoCNVariant, log);
			}
			reader.close();
			// finalize any remaining CytoCNVariants in tmps and add to cytoCNVariant
			finalAction(cytoCNVariant, tmps);
			log.report(ext.getTime() + " Info - found a total of " + cytoCNVariant.size() + " consolidated aberrations across " + sampleNames.length + " samples in file " + cytoCNVariantFile);

		} catch (FileNotFoundException e) {
			log.reportError("Error - could not find file " + cytoCNVariantFile);
			log.reportException(e);
		} catch (IOException e) {
			log.reportError("Error - reading file" + cytoCNVariantFile);
			log.reportException(e);
		}
		return cytoCNVariant.toArray(new CytoCNVariant[cytoCNVariant.size()]);
	}

	/**
	 * Get just the sample names present in the file
	 */
	public static String[] getSampleNames(String cytoCNVariantFile, Logger log) {
		String[] line;
		String[] sampleNames = new String[0];
		int scan = 0;
		int count = 0;
		try {
			BufferedReader reader = Files.getAppropriateReader(cytoCNVariantFile);
			do {
				line = reader.readLine().trim().split(SPLITS[0], -1);
				count++;
			} while (reader.ready() && (ext.indexFactors(CYTO_ABERATION_HEADER, line, false, true, false, false)[0] == -1));
			if (!reader.ready()) {
				log.reportError("Error - did not find the neccesary column headers in file " + cytoCNVariantFile);
				return sampleNames;
			}
			reader.close();

			reader = Files.getAppropriateReader(cytoCNVariantFile);
			// the sample names are located right above the data header, so we re-scan the reader to that location
			do {
				line = reader.readLine().trim().split(SPLITS[0], -1);
				scan++;
			} while (scan < count - 1);

			sampleNames = getSamples(line);
		} catch (FileNotFoundException e) {
			log.reportError("Error - could not find file " + cytoCNVariantFile);
			log.reportException(e);
		} catch (IOException e) {
			log.reportError("Error - reading file" + cytoCNVariantFile);
			log.reportException(e);
		}
		return sampleNames;
	}

	/**
	 * @param cytoCNVariantFile
	 * @param log
	 * @return the cytoCNVariantFile loaded directly to individual format
	 */
	public static CytoCNVariant[][] directToInds(String cytoCNVariantFile, Logger log) {
		return toIndividuals(loadCytoCNVariant(cytoCNVariantFile, log), log);
	}

	/**
	 * The base file name for each individual is FID_IID.cnv, we assume that FID_IID are unique
	 * 
	 * @param cytoCNVariantInds
	 *            CytoCNVariant[][] organized CytoCNVariant[sample0][CytoCNVariantsForSample0]
	 * @param dir
	 *            added to base file name of FID_IID,if not null
	 * @param log
	 * 
	 * @return String[] of the output files (full paths) written to: Final output will be dir (if not null)+FID_IID.cnv
	 * 
	 */
	public static String[] writeIndCNVariantFiles(CytoCNVariant[][] cytoCNVariantInds, String dir, Logger log) {
		String[] outputs = new String[cytoCNVariantInds.length];
		for (int i = 0; i < cytoCNVariantInds.length; i++) {
			if (cytoCNVariantInds[i] != null && cytoCNVariantInds[i].length > 0) {
				String FID_IID = cytoCNVariantInds[i][0].getFamilyID() + "_" + cytoCNVariantInds[i][0].getIndividualID();
				String fileName = ext.replaceWithLinuxSafeCharacters(FID_IID, true);
				if (dir != null) {
					fileName = dir + fileName;
				}
				fileName = fileName + CNV_EXT;
				outputs[i] = fileName;
				writeToPlink(cytoCNVariantInds[i], fileName, log);
			} else {
				log.reportError("Warning - index " + i + " did not contain any variants");
				outputs[i] = "";
			}
		}
		return outputs;
	}

	/**
	 * @param CytoCNVariant
	 *            CytoCNVariant[] to write
	 * @param fileName
	 *            full path
	 * @param log
	 */
	public static void writeToPlink(CytoCNVariant[] CytoCNVariant, String fileName, Logger log) {
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(fileName));
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
			for (int i = 0; i < CytoCNVariant.length; i++) {
				writer.println(CytoCNVariant[i].toPlinkFormat());
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + fileName);
			log.reportException(e);
		}
	}

	/**
	 * Here we check the corresponding amplification/deletions for each sample.
	 * <p>
	 * If it is adding to a previous aberation (same logRatio score, we update the cytoCNVariant to reflect the new info (new stop position, add any new markers/genes)
	 * <p>
	 * If it is not the same, we create a new cytoCNVariant to track for the corresponding sample. If the aberration has ended and a new one is not immediately following, we finalize the variant's arrays and add it to tmps
	 * 
	 * @param sampleNames
	 * @param line
	 * @param indicesCommon
	 * @param indicesSamples
	 * @param tmps
	 *            temporary variants corresponding to one per sample as we scan the rows. Assuming no amplification/duplication combinations
	 * @param cytoCNVariant
	 *            the final parsed variants
	 * @param log
	 */
	private static void parseInputLine(String[] sampleNames, String[] line, int[] indicesCommon, int[][] indicesSamples, CytoCNVariant[] tmps, ArrayList<CytoCNVariant> cytoCNVariant, Logger log) {

		try {
			// parse the common info for all samples on this line
			String cytoBand = line[indicesCommon[1]];
			byte chr = Positions.chromosomeNumber(line[indicesCommon[2]], log);
			String probe = line[indicesCommon[3]];
			int start = Integer.parseInt(line[indicesCommon[4]]);
			int stop = Integer.parseInt(line[indicesCommon[5]]);
			String gene = line[indicesCommon[6]];
			for (int i = 0; i < indicesSamples.length; i++) {
				String sampAmp = getSafeIndex(line, indicesSamples[i][1]);
				String sampDel = getSafeIndex(line, indicesSamples[i][2]);

				double logRatio = Double.NaN;
				int CN = 0;

				// make sure there is not both a deletion and amplification for a particular sample
				if (ext.indexOfStr(sampAmp, NOTHINGS, true, true) == -1 && ext.indexOfStr(sampDel, NOTHINGS, true, true) == -1) {
					log.reportError("Error - found an amplification and a deletion on line " + Array.toStr(line));
					return;
				}
				if (ext.indexOfStr(sampAmp, NOTHINGS, true, true) == -1) {
					logRatio = Double.parseDouble(sampAmp);
					CN = CN_NUMBER[1];
				} else if (ext.indexOfStr(sampDel, NOTHINGS, true, true) == -1) {
					logRatio = Double.parseDouble(sampDel);
					CN = CN_NUMBER[0];
				}
				if (!Double.isNaN(logRatio) && CN != 0) {// if there is an amplification/deletion for the current sample on this line,
					if (tmps[i] != null) { // if one existed on the previous line
						if (tmps[i].isSame(logRatio) && tmps[i].isSameChr((byte) chr)) {// if it is the same aberration as the previous line we update it
							if (!tmps[i].update(cytoBand, chr, stop, probe, gene, logRatio, log)) {
								log.reportError("Error - could not update the variant for sample " + sampleNames[i] + " on line " + Array.toStr(line));
							}
						} else {// else we finalize the previous and create a new one
							tmps[i].finalizeVariant();
							cytoCNVariant.add(tmps[i]);
							tmps[i] = new CytoCNVariant("" + (i + 1), ext.replaceWithLinuxSafeCharacters(sampleNames[i],true), chr, start, stop, CN, DEFAULT_SCORE, 1, 0, cytoBand, logRatio, probe, gene);
						}
					} else {// else we create a new one
						tmps[i] = new CytoCNVariant("" + (i + 1), ext.replaceWithLinuxSafeCharacters(sampleNames[i],true), chr, start, stop, CN, DEFAULT_SCORE, 1, 0, cytoBand, logRatio, probe, gene);
					}
				} else if (tmps[i] != null) {// else we check if there was one on the previous line and finalize it
					tmps[i].finalizeVariant();
					cytoCNVariant.add(tmps[i]);
					tmps[i] = null;
				}
			}
		} catch (NumberFormatException nfe) {
			log.reportError("Error - found an invalid number on that could not be parsed on line " + Array.toStr(line));
			return;
		}
	}

	/**
	 * To finalize the last aberration(s) and add to the ArrayList of parsed aberrations
	 * 
	 * @param cytoCNVariant
	 * @param tmps
	 */
	private static void finalAction(ArrayList<CytoCNVariant> cytoCNVariant, CytoCNVariant[] tmps) {
		for (int i = 0; i < tmps.length; i++) {
			if (tmps[i] != null) {
				tmps[i].finalizeVariant();
				cytoCNVariant.add(tmps[i]);
				tmps[i] = null;
			}
		}
	}

	/**
	 * Need this due to strange blank formatting where the line length can vary, causing index out of bounds exceptions
	 * 
	 * @param line
	 * @param index
	 * @return "" if out of index
	 */
	private static String getSafeIndex(String[] line, int index) {
		String st;
		if (index > line.length - 1) {
			st = NOTHINGS[2];
		} else {
			st = line[index];
		}
		return st;
	}

	/**
	 * The assumption here is that the sample names are the only non-blank entries on the line right above the data header. The non-blank entries will be used for the FID/IID
	 * 
	 * @param line
	 * @return
	 */
	private static String[] getSamples(String[] line) {
		ArrayList<String> tmpSamps = new ArrayList<String>();
		for (int i = 0; i < line.length; i++) {
			if (!line[i].equals("")) {
				tmpSamps.add(line[i]);
			}
		}
		return tmpSamps.toArray(new String[tmpSamps.size()]);
	}

	/**
	 * 
	 * We use this function to get the indices of each sample's (CYTO_SAMPLE_HEADER) info for the file. Also, each of the samples must have all three columns present
	 * 
	 * @param header
	 * @param log
	 * @return
	 */
	private static int[][] determineSampleIndices(String[] header, Logger log) {
		int[][] sampleIndices;
		int count = 0;
		for (int i = 0; i < header.length; i++) {
			if (header[i].equals(CYTO_SAMPLE_HEADER[0])) {
				count++;
			}
		}
		sampleIndices = new int[count][3];
		int sampleIndex = 0;
		int numFields = 0;
		for (int i = 0; i < header.length; i++) {
			int index = ext.indexOfStr(header[i], CYTO_SAMPLE_HEADER);
			if (index >= 0) {
				sampleIndices[sampleIndex][index] = i;
				numFields++;
			}
			if (numFields > 2) {
				sampleIndex++;
				numFields = 0;
			}
		}
		if (sampleIndex != sampleIndices.length) {
			log.reportError("Error - did not find the necessary column headers \"" + Array.toStr(CYTO_SAMPLE_HEADER) + "\" for all samples");
		}
		return sampleIndices;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String cytoCNVariantFile = "D:/data/Hirch_CYTO/Hirch_Cytolab/5 6 7 8 workbench.xls";
		String logFile = "D:/data/Hirch_CYTO/Hirch_Cytolab/testlog.log";
		String usage = "\n" + "jlDev.CytoCNVariant requires 0-9 arguments\n";
		usage += "   (1) cyto call file (i.e. filename=" + cytoCNVariantFile + " (default))\n";
		usage += "   (2) name of the log file (i.e. logFile=" + cytoCNVariantFile + " (default))\n";
		usage += "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("filename=")) {
				cytoCNVariantFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("logFile=")) {
				logFile = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		Logger log = new Logger(logFile);
		test(cytoCNVariantFile, log);
	}

	public static void test(String cytoCNVariantFile, Logger log) {
		CytoCNVariant[] tests = loadCytoCNVariant(cytoCNVariantFile, log);
		for (int i = 0; i < tests.length; i++) {
			log.report(tests[i].getMyReport());
		}
	}
}
