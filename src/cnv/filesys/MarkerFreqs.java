package cnv.filesys;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.Hashtable;

import common.Array;
import common.Files;
import common.Logger;
import common.ext;

//class to store markerQC metrics across all (or a subset of) samples , currently only Minor Allele Frequency
public class MarkerFreqs implements Serializable {
	public static final long serialVersionUID = 1L;
	public static final String[] MARKER_FREQ_FIELDS = { "Name", "MAF" };

	private double[] mafs;
	private long fingerprint;

	public MarkerFreqs(double[] mafs, long fingerprint) {
		this.mafs = mafs;
		this.fingerprint = fingerprint;
	}

	public double[] getMafs() {
		return mafs;
	}

	public long getFingerprint() {
		return fingerprint;
	}

	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static MarkerFreqs load(String filename, boolean jar) {
		if (Files.exists(filename)) {
			return (MarkerFreqs) Files.readSerial(filename, jar, true);
		} else {
			return null;
		}
	}

	public static void exportToText(String filename, String mafFilename, String exportFilename) {
		Project proj = new Project(filename, false);
		MarkerFreqs markerMAF = load(proj.getProjectDir() + mafFilename, false);
		MarkerSet markerSet = proj.getMarkerSet();
		String[] markerNames = markerSet.getMarkerNames();
		double[] mafs = markerMAF.getMafs();
		Logger log;
		
		log = proj.getLog();
		if (markerNames.length != mafs.length) {
			log.reportError("Error - mismatched number of markers in the project's marker set and the imported AlleleFrequency file (" + mafFilename + "); aborting");
			System.exit(1);
		}
		if (markerSet.getFingerprint() != markerMAF.getFingerprint()) {
			log.reportError("Error - mismatched marker fingerprints in the project's marker set and the imported AlleleFrequency file (" + mafFilename + "); aborting");
			System.exit(1);
		} else {
			dump(proj.getProjectDir() + exportFilename, mafs, markerNames, markerSet.getFingerprint(), log);
		}
	}

	private static void dump(String exportFilename, double[] mafs, String[] markerNames, long fingerprint, Logger log) {
		PrintWriter writer;
		try {
			writer = new PrintWriter(new FileWriter(exportFilename));
			writer.println("marker_fingerprint=" + fingerprint);
			writer.println(Array.toStr(MARKER_FREQ_FIELDS));
			for (int i = 0; i < markerNames.length; i++) {
				writer.println(markerNames[i] + "\t" + mafs[i]);
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + exportFilename);
			e.printStackTrace();
			System.exit(1);
		}
	}

	public static void convertMarkerFreqsFromTxt(Project proj, String Freqfilename, String outputFileNameSer) {
		BufferedReader reader;
		String[] line, header;
		Hashtable<String, String> hash;
		int[] indices;
		String[] markerNames;
		int index;
		MarkerSet markerSet;
		Logger log;
		
		log = proj.getLog();
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		double[] mafs = new double[markerNames.length];
		hash = new Hashtable<String, String>();
		for (int i = 0; i < markerNames.length; i++) {
			hash.put(markerNames[i], i + "");
		}
		try {
			int numMafs = 0;
			reader = new BufferedReader(new FileReader(proj.getProjectDir() + Freqfilename));
			header = reader.readLine().trim().split("[\\s]+");
			indices = Array.intArray(MARKER_FREQ_FIELDS.length, -1);
			for (int i = 0; i < header.length; i++) {
				index = ext.indexOfEndsWith(header[i], MARKER_FREQ_FIELDS, true);
				if (index >= 0) {
					indices[index] = i;
				}
			}
			if (Array.min(indices) == -1) {
				log.reportError("Error - Need a column header ending with the following suffixes; missing at least one");
				log.reportError("        " + Array.toStr(MARKER_FREQ_FIELDS, "  "));
				System.exit(1);
			}
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (!hash.containsKey(line[indices[0]])) {
					log.reportError("Error - marker '" + line[indices[0]] + "' was not found in MarkerSet");
					System.exit(1);
				}
				numMafs++;
				index = Integer.parseInt(hash.get(line[indices[0]]));
				mafs[index] = Double.parseDouble(line[1]);
			}
			if (numMafs != markerNames.length) {
				log.reportError("Error - " + (markerNames.length - numMafs) + "markers were not found in MarkerSet");
				System.exit(1);
			}

			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + Freqfilename + "\" not found in " + proj.getProjectDir());
			System.exit(1);
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + Freqfilename + "\"");
			System.exit(2);
		}
		new MarkerFreqs(mafs, markerSet.getFingerprint()).serialize(proj.getProjectDir() + outputFileNameSer);
	}

}
