package org.genvisis.filesys;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;

import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;

public class ResultSet implements Serializable {
	public static final long serialVersionUID = 1L;
	public static final int PVALUES_ONLY = 0;
	public static final int GENERIC_FORMAT = 1;
	public static final int METAL_TBL_FORMAT = 2;
	public static final int PROBABEL_FORMAT = 3;

	/** Marker name, Chr, Position, centiMorgans, A1, A2, annotation, header */
	public static final int[][] SNP_INDICES = {	{0, -1, -1, -1, -1, -1, -1, 1},
																							{0, 1, 2, -1, -1, -1, -1, 1},
																							{0, -1, -1, -1, 1, 2, -1, 1},
																							{0, -1, -1, -1, 1, 2, -1, 1},};

	/** Effect, Stderr, p-value, annotation, header */
	public static final int[][] RESULT_INDICES = {{-1, -1, 1, -1, 1}, {-1, -1, 1, -1, 1},
																								{3, 4, 5, 6, 1}, {10, 11, 12, -1, 1},};

	private final SnpMarkerSet markerSet;
	private float[] effects;
	private float[] stderrs;
	private float[] pvals;
	private String[] annotation;


	public ResultSet(String filename, int format) {
		this(filename, SNP_INDICES[format], SNP_INDICES[format][7] == 1, RESULT_INDICES[format]);
	}

	public ResultSet(String filename, int[] snp_indices, boolean header, int[] result_indices) {
		BufferedReader reader;
		String[] line;
		int numMarkers;

		markerSet = new SnpMarkerSet(filename, snp_indices, header, true, new Logger());
		numMarkers = markerSet.getMarkerNames().length;
		try {
			reader = new BufferedReader(new FileReader(filename));
			if (result_indices[4] == 1) {
				reader.readLine();
			}
			if (result_indices[0] == -1) {
				effects = null;
			} else {
				effects = new float[numMarkers];
			}
			if (result_indices[1] == -1) {
				stderrs = null;
			} else {
				stderrs = new float[numMarkers];
			}
			if (result_indices[2] == -1) {
				pvals = null;
			} else {
				pvals = new float[numMarkers];
			}
			if (result_indices[3] == -1) {
				annotation = null;
			} else {
				annotation = new String[numMarkers];
			}
			for (int i = 0; i < numMarkers; i++) {
				line = reader.readLine().trim().split("[\\s]+");
				if (Array.max(result_indices) >= line.length) {
					System.err.println("Error - required column does not exist for the marker in row "
															+ (i + 1));
					System.exit(1);
				}

				if (result_indices[0] != -1) {
					effects[i] = Float.parseFloat(line[result_indices[0]]);
				}
				if (result_indices[1] != -1) {
					stderrs[i] = Float.parseFloat(line[result_indices[1]]);
				}
				if (result_indices[2] != -1) {
					pvals[i] = Float.parseFloat(line[result_indices[2]]);
				}
				if (result_indices[3] != -1) {
					annotation[i] = line[result_indices[3]];
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
	}

	public String[] getMarkerNames() {
		return markerSet.getMarkerNames();
	}

	public byte[] getChrs() {
		return markerSet.getChrs();
	}

	public int[] getPositions() {
		return markerSet.getPositions();
	}


	public float[] getEffects() {
		return effects;
	}

	public float[] getStderrs() {
		return stderrs;
	}

	public float[] getPvals() {
		return pvals;
	}

	public String[] getAnnotation() {
		return annotation;
	}

	public void writeToFile(String filename, int format) {
		writeToFile(filename, SNP_INDICES[format], RESULT_INDICES[format]);
	}

	public void writeToFile(String filename, int[] snp_indices, int[] result_indices) {
		PrintWriter writer;
		String[] line, markerNames;
		byte[] chrs;
		int[] positions;
		char[][] alleles;

		markerNames = markerSet.getMarkerNames();
		chrs = markerSet.getChrs();
		positions = markerSet.getPositions();
		alleles = markerSet.getAlleles();

		line = new String[Math.max(Array.max(snp_indices), Array.max(result_indices)) + 1];
		try {
			writer = new PrintWriter(new FileWriter(filename));
			for (int i = 0; i < markerNames.length; i++) {
				for (int j = 0; j < line.length; j++) {
					line[j] = "0";
				}
				if (snp_indices[0] != -1 && markerNames != null) {
					line[snp_indices[0]] = markerNames[i];
				}
				if (snp_indices[1] != -1 && chrs != null) {
					line[snp_indices[1]] = chrs[i] + "";
				}
				if (snp_indices[2] != -1 && positions != null) {
					line[snp_indices[2]] = positions[i] + "";
				}
				if (snp_indices[3] != -1 && alleles != null) {
					line[snp_indices[3]] = alleles[i][0] + "";
				}
				if (snp_indices[4] != -1 && alleles != null) {
					line[snp_indices[4]] = alleles[i][1] + "";
				}
				if (result_indices[0] != -1 && effects != null) {
					line[result_indices[0]] = effects[i] + "";
				}
				if (result_indices[1] != -1 && stderrs != null) {
					line[result_indices[1]] = stderrs[i] + "";
				}
				if (result_indices[2] != -1 && pvals != null) {
					line[result_indices[2]] = pvals[i] + "";
				}
				if (result_indices[3] != -1 && annotation != null) {
					line[result_indices[3]] = annotation[i];
				}
				writer.println(Array.toStr(line));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing " + filename);
			e.printStackTrace();
		}
	}

	public void serialize(String filename) {
		SerializedFiles.writeSerial(this, filename);
	}

	public static ResultSet load(String filename, boolean jar, boolean kill) {
		return (ResultSet) SerializedFiles.readSerial(filename, jar, kill);
	}

}
