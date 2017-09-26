package org.genvisis.filesys;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.List;
import java.util.Map.Entry;

import org.genvisis.common.Aliases;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

public class Hits {

	private final Hashtable<String, Double> hash;

	public Hits() {
		hash = new Hashtable<String, Double>();
	}

	public void incorporateFromFile(String filename, double threshold, Logger log) {
		int[] indices;

		indices = ext.indexFactors(new String[][] {Aliases.MARKER_NAMES, Aliases.PVALUES},
															 Files.getHeaderOfFile(filename, log), false, true, true, log);
		incorporateFromFile(filename, indices, threshold, log);
	}

	public void incorporateFromFile(String filename, int[] indices, double threshold, Logger log) {
		BufferedReader reader;
		String[] line;
		String delimiter;
		double value;
		int testColumn;

		try {
			reader = new BufferedReader(new FileReader(filename));
			delimiter = Files.determineDelimiter(filename, log);
			line = ext.splitLine(reader.readLine(), delimiter, log);
			testColumn = ext.indexOfStr("TEST", line);
			while (reader.ready()) {
				line = ext.splitLine(reader.readLine(), delimiter, log);
				if (!ext.isMissingValue(line[indices[1]])
						&& (testColumn == -1 || line[testColumn].equals("ADD"))) {
					try {
						value = Double.parseDouble(line[indices[1]]);
					} catch (NumberFormatException nfe) {
						System.err.println("Error - invalid p-value (" + line[indices[1]] + ") for marker "
															 + line[indices[0]]);
						value = 1;
					}
					if (value <= threshold) {
						if (hash.containsKey(line[indices[0]])) {
							value = Math.min(value, hash.get(line[indices[0]]));
						}
						hash.put(line[indices[0]], value);
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			return;
		}
	}

	public void writeHits(String filename) {
		PrintWriter writer;

		List<Entry<String, Double>> entries = Sort.entriesSortedByValues(hash);

		try {
			writer = Files.openAppropriateWriter(filename);
			for (Entry<String, Double> e : entries) {
				writer.println(e.getKey() + "\t" + e.getValue());
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + filename);
			e.printStackTrace();
		}
	}

	public String[] getKeys() {
		return HashVec.getKeys(hash);
	}

}
