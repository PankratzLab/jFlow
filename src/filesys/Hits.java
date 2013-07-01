package filesys;

import java.io.*;
import java.util.*;
import common.*;
import gwas.Metal;

public class Hits {

	private Hashtable<String,Double> hash;
	
	public Hits() {
		hash = new Hashtable<String, Double>();
	}
	
	public void incorporateFromFile(String filename, double threshold, Logger log) {
		int[] indices;

		indices = ext.indexFactors(new String[][] {Metal.MARKER_NAMES, Metal.PVALUES}, Files.getHeaderOfFile(filename, log), false, true, true, log, true);
		incorporateFromFile(filename, indices, threshold, log);		
	}
	
	public void incorporateFromFile(String filename, int[] indices, double threshold, Logger log) {
		BufferedReader reader;
		String[] line;
		String delimiter;
		double value;
		
		try {
			reader = new BufferedReader(new FileReader(filename));
			delimiter = Files.determineDelimiter(filename, log);
			reader.readLine();
			while (reader.ready()) {
				if (delimiter.equals("\t")) {
					line = reader.readLine().split(delimiter, -1);
				} else if (delimiter.equals(",")) {
					line = ext.splitCommasIntelligently(reader.readLine(), true, log);
				} else {
					line = reader.readLine().trim().split(delimiter);
				}
				if (!ext.isMissingValue(line[indices[1]])) {
					try {
						value = Double.parseDouble(line[indices[1]]);
					} catch (NumberFormatException nfe) {
						System.err.println("Error - invalid p-value ("+line[indices[1]]+") for marker "+line[indices[0]]);
						value = 1;
					}
					if (value < threshold) {
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
		String[] keys;
		int[] order;
		double[] values;
		
		keys = HashVec.getKeys(hash, false, false);
		values = new double[keys.length];
		for (int i = 0; i < keys.length; i++) {
			values[i] = hash.get(keys[i]);
		}
		order = Sort.quicksort(values);
		
		try {
			writer = new PrintWriter(new FileWriter(filename));
			for (int i = 0; i < order.length; i++) {
				writer.println(keys[order[i]] +"\t"+ values[order[i]]);
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
