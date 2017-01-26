package org.genvisis.common;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Vector;

public class HashVec {
	public static boolean addIfAbsent(String str, Vector<String> v) {
		return addIfAbsent(str, v, false);
	}

	/**
	 * @param inAlphabeticalOrder - if the vector is already sorted and the new element should be at
	 *        the appropriate sorted index.
	 * @return true if added
	 */
	public static boolean addIfAbsent(String str, Vector<String> v, boolean inAlphabeticalOrder) {

		if (inAlphabeticalOrder) {
			int idx = Collections.binarySearch(v, str);
			// negative value indicates value was not present.
			// adjust index and insert
			if (idx < 0) {
				idx = -(idx + 1);
				v.insertElementAt(str, idx);
				return true;
			}
		} else if (!v.contains(str)) {
			v.add(str);
			return true;
		}

		return false;
	}

	public static void addAllInArrayToVector(String[] array, Vector<String> vector) {
		for (String element : array) {
			vector.add(element);
		}
	}

	public static Object addToHashIfAbsent(Hashtable<String, Object> hash, String key, Object value) {
		Object o;

		if (hash.containsKey(key)) {
			o = hash.get(key);
		} else {
			hash.put(key, o = value);
		}

		return o;
	}

	public static void addToHashIfAbsent(Hashtable<String, String> hash, String key, String value) {
		if (!hash.containsKey(key)) {
			hash.put(key, value);
		}
	}

	public static void addToHashIfAbsent(	Hashtable<String, String> hash, String[] keys,
																				String[] values) {
		for (int i = 0; i < keys.length; i++) {
			if (!hash.containsKey(keys[i])) {
				hash.put(keys[i], values == null ? "" : values[i]);
			}
		}
	}

	public static void addToHashVec(Hashtable<String, Vector<String>> hash, String key, String value,
																	boolean onlyifabsent) {
		Vector<String> v;

		if (hash.containsKey(key)) {
			v = hash.get(key);
		} else {
			hash.put(key, v = new Vector<String>());
		}

		if (onlyifabsent) {
			if (!v.contains(value)) {
				v.add(value);
			}
		} else {
			v.add(value);
		}
	}

	public static void addToHashArrayVec(	Hashtable<String, Vector<String[]>> hash, String key,
																				String[] array) {
		Vector<String[]> v;

		if (hash.containsKey(key)) {
			v = hash.get(key);
		} else {
			hash.put(key, v = new Vector<String[]>());
		}

		v.add(array);
	}

	public static void addToHashHash(	Hashtable<String, Hashtable<String, Object>> hash1, String key1,
																		String key2, Object value) {
		Hashtable<String, Object> hash2;

		if (hash1.containsKey(key1)) {
			hash2 = hash1.get(key1);
		} else {
			hash1.put(key1, hash2 = new Hashtable<String, Object>());
		}
		hash2.put(key2, value);
	}

	public static void addToHashHash(	Hashtable<String, Hashtable<String, String>> hash1, String key1,
																		String key2, String value) {
		Hashtable<String, String> hash2;

		if (hash1.containsKey(key1)) {
			hash2 = hash1.get(key1);
		} else {
			hash1.put(key1, hash2 = new Hashtable<String, String>());
		}
		hash2.put(key2, value);
	}

	public static void addToHashHashVec(Hashtable<String, Hashtable<String, Vector<String>>> hash1,
																			String key1, String key2, String value,
																			boolean onlyifabsent) {
		Hashtable<String, Vector<String>> hash2;
		Vector<String> v;

		if (hash1.containsKey(key1)) {
			hash2 = hash1.get(key1);
		} else {
			hash1.put(key1, hash2 = new Hashtable<String, Vector<String>>());
		}
		if (hash2.containsKey(key2)) {
			v = hash2.get(key2);
		} else {
			hash2.put(key2, v = new Vector<String>());
		}
		if (onlyifabsent) {
			addIfAbsent(value, v);
		} else {
			v.add(value);
		}
	}

	public static Hashtable<String, String> loadFileToHashString(	String filename,
																																boolean ignoreFirstLine) {
		return loadFileToHashString(filename, 0, new int[] {1}, null, ignoreFirstLine);
	}

	public static HashSet<String> loadFileToHashSet(String filename, boolean ignoreFirstLine) {
		return convertHashNullToHashSet(loadFileToHashString(filename, 0, null, null, ignoreFirstLine));
	}

	public static HashSet<String> loadFileToHashSet(String filename, int[] keyIndices,
																									String delimiterWithinHash,
																									boolean ignoreFirstLine) {
		return convertHashNullToHashSet(loadFileToHashString(	filename, keyIndices, null, false,
																													delimiterWithinHash, ignoreFirstLine,
																													false, false));
	}

	public static HashSet<String> loadToHashSet(String[] list) {
		if (list == null)
			return new HashSet<String>();
		HashSet<String> hash = new HashSet<String>(list.length);

		for (String s : list) {
			hash.add(s);
		}

		return hash;
	}

	@SuppressWarnings("rawtypes")
	public static HashSet<String> convertHashNullToHashSet(Hashtable hash) {
		HashSet<String> set = new HashSet<String>();
		for (Object k : hash.keySet()) {
			set.add(k.toString());
		}
		return set;
	}

	public static Hashtable<String, Integer> loadToHashIndices(String[] list, Logger log) {
		Hashtable<String, Integer> hash = new Hashtable<String, Integer>((list == null	? 10
																																										: list.length));

		for (int i = 0; list != null && i < list.length; i++) {
			if (list[i] == null) {
				log.report("Cannot add null as a key in a hashIndices object");
			} else {
				hash.put(list[i], i);
			}
		}

		return hash;
	}

	public static Vector<String> loadFileToVec(	String filename, boolean ignoreFirstLine,
																							boolean onlyFirstColumn, boolean onlyIfAbsent) {
		return loadFileToVec(	filename, ignoreFirstLine, onlyFirstColumn ? new int[] {0} : null,
													onlyIfAbsent, false);
	}

	public static Vector<String> loadFileToVec(	String filename, boolean ignoreFirstLine,
																							boolean onlyFirstColumn, boolean onlyIfAbsent,
																							boolean jar) {
		return loadFileToVec(	filename, ignoreFirstLine, onlyFirstColumn ? new int[] {0} : null,
													onlyIfAbsent, jar);
	}

	public static String[] loadFileToStringArray(	String filename, boolean ignoreFirstLine, int[] cols,
																								boolean onlyIfAbsent) {
		Vector<String> v = loadFileToVec(filename, ignoreFirstLine, cols, onlyIfAbsent, false);
		return v == null ? null : ArrayUtils.toStringArray(v);
	}

	public static String[] loadFileToStringArray(	String filename, boolean jar,
																								boolean ignoreFirstLine, int[] cols,
																								boolean onlyIfAbsent) {
		Vector<String> v = loadFileToVec(filename, ignoreFirstLine, cols, onlyIfAbsent, jar);
		return v == null ? null : ArrayUtils.toStringArray(v);
	}

	public static String[] loadFileToStringArray(	String filename, boolean jar,
																								boolean ignoreFirstLine, int[] cols,
																								boolean trimFirst, boolean onlyIfAbsent,
																								String delimiter) {
		Vector<String> v = loadFileToVec(	filename, ignoreFirstLine, cols, trimFirst, onlyIfAbsent, jar,
																			delimiter);
		return v == null ? null : ArrayUtils.toStringArray(v);
	}

	public static Vector<String> loadFileToVec(	String filename, boolean ignoreFirstLine, int[] cols,
																							boolean onlyIfAbsent, boolean jar) {
		return loadFileToVec(	filename, ignoreFirstLine, cols, true, onlyIfAbsent, jar,
													Files.determineDelimiter(filename, new Logger()));
	}

	public static Vector<String> loadFileToVec(	String filename, boolean ignoreFirstLine, int[] cols,
																							boolean trimFirst, boolean onlyIfAbsent, boolean jar,
																							String delimiter) {
		BufferedReader reader = null;
		Vector<String> v = new Vector<String>();
		HashSet<String> onlyIfAbsentHash = new HashSet<String>();
		String trav;
		String[] line;
		int count;

		try {
			reader = Files.getReader(filename, jar, true, false);
			if (reader == null) {
				return null; // Should return empty? No - empty could be valid, we need to show something
											// invalid, so null or exception
			}
			if (ignoreFirstLine) {
				reader.readLine();
			}
			count = 1;
			while (reader.ready()) {
				trav = reader.readLine();
				if (cols != null) {
					if (trimFirst) {
						trav = trav.trim(); // trim() needed for all PLINK files
					}
					if (delimiter.equals(",")) {
						line = ext.splitCommasIntelligently(trav, true, new Logger());
					} else {
						line = trav.split(delimiter, -1);
					}
					trav = "";
					for (int i = 0; i < cols.length; i++) {
						if (line.length <= cols[i]) {
							System.err.println("Error - not enough columns at line "+ count + " of file "
																	+ filename + ": " + ArrayUtils.toStr(line));
						}
						trav += (i == 0 ? "" : "\t") + line[cols[i]];
					}
				}
				if (!onlyIfAbsent || !onlyIfAbsentHash.contains(trav)) {
					v.add(trav);
					if (onlyIfAbsent) {
						onlyIfAbsentHash.add(trav);
					}
				}
				count++;
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}

		return v;
	}

	public static String[][] loadFileToStringMatrix(String filename, boolean ignoreFirstLine,
																									int[] cols, boolean jar) {
		return loadFileToStringMatrix(filename, ignoreFirstLine, cols,
																	Files.determineDelimiter(filename, new Logger()), jar, 1000,
																	false);
	}

	public static String[][] loadFileToStringMatrix(String filename, boolean ignoreFirstLine,
																									int[] cols, String delimiter, boolean jar,
																									int initialCapacity, boolean allowMissingData) {
		BufferedReader reader = null;
		Vector<String[]> v = new Vector<String[]>(initialCapacity);
		String temp;
		String[] line, data;

		try {
			reader = Files.getReader(filename, jar, true, false);
			if (reader == null) {
				return null;
			}
			if (ignoreFirstLine) {
				reader.readLine();
			}
			while ((temp = reader.readLine()) != null) {
				line = ext.splitLine(temp, delimiter, null);
				if (cols == null) {
					// if (!temp.equals("") || allowMissingData) {
					v.add(line);
					// }
				} else {
					data = new String[cols.length];
					for (int i = 0; i < cols.length; i++) {
						if (allowMissingData && line.length <= cols[i]) {
							data[i] = null;
						} else {
							data[i] = line[cols[i]];
						}
					}
					v.add(data);
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

		return Matrix.toStringArrays(v);
	}

	public static Hashtable<String, String> loadFileToHashString(	String filename, String keyHeader,
																																String[] valueHeaders,
																																String delimiterWithinHash) {
		BufferedReader reader = null;
		String[] line;
		int keyIndex;
		int[] valueIndices;

		try {
			reader = new BufferedReader(new FileReader(filename));
			line = reader.readLine().trim().split("\t", -1);
			keyIndex = ext.indexOfStr(keyHeader, line);
			if (keyIndex == -1) {
				System.err.println("Error - '" + keyHeader + "' not found in " + filename);
				System.exit(1);
			}
			valueIndices = ext.indexFactors(valueHeaders, line, false, true);

			reader.close();
			return loadFileToHashString(filename, keyIndex, valueIndices, delimiterWithinHash, true);
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}

		return null;
	}

	public static Hashtable<String, String> loadFileToHashString(	String filename, int keyIndex,
																																int[] valueIndices,
																																String delimiterWithinHash,
																																boolean ignoreFirstLine) {
		return loadFileToHashString(filename, keyIndex, valueIndices, delimiterWithinHash,
																ignoreFirstLine, false);
	}

	public static Hashtable<String, String> loadFileToHashString(	String filename, int keyIndex,
																																int[] valueIndices,
																																String delimiterWithinHash,
																																boolean ignoreFirstLine,
																																boolean jar) {
		return loadFileToHashString(filename, new int[] {keyIndex}, valueIndices, false,
																delimiterWithinHash, ignoreFirstLine, jar, false);
	}

	public static Hashtable<String, String> loadFileToHashString(	String filename, int[] keyIndices,
																																int[] valueIndices,
																																boolean commaDelimitedFile,
																																String delimiterWithinHash,
																																boolean ignoreFirstLine,
																																boolean jar,
																																boolean allowMissingData) {
		BufferedReader reader = null;
		String[] line;
		Hashtable<String, String> hash = new Hashtable<String, String>();
		String key, temp;

		if (valueIndices == null) {
			valueIndices = new int[] {};
		}

		try {
			reader = Files.getReader(filename, jar, true, true);
			if (ignoreFirstLine) {
				reader.readLine();
			}
			while (reader.ready()) {
				temp = reader.readLine();
				if (filename.endsWith(".csv") || commaDelimitedFile) {
					// line = temp.split(",", -1);
					line = ext.splitCommasIntelligently(temp, true, new Logger());
				} else if (temp.indexOf("\t") == -1) {
					line = temp.trim().split("[\\s]+");
				} else {
					line = temp.split("\t", -1);
				}
				key = "";
				for (int i = 0; i < keyIndices.length; i++) {
					key += ((i == 0) ? "" : delimiterWithinHash) + line[keyIndices[i]];
				}
				temp = "";
				for (int i = 0; i < valueIndices.length; i++) {
					if (valueIndices[i] == -7) {
						temp += ((i == 0) ? "" : delimiterWithinHash) + hash.size();
					} else if (valueIndices[i] < line.length) {
						temp += ((i == 0) ? "" : delimiterWithinHash) + line[valueIndices[i]];
					} else if (allowMissingData) {
						temp += ((i == 0) ? "" : delimiterWithinHash) + ".";
					} else {
						System.err.println("Error - not enough columns for key '"+ key + "' in file '"
																+ filename + "'; and allowMissingData was not flagged");
					}
				}
				hash.put(key, temp);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
		}

		return hash;
	}

	/**
	 * As {@link #loadFileToHashVec(String, int, int[], String, boolean, boolean)} but all columns of
	 * interest can be enumerated, avoiding the need to look up their indices. This method assumes
	 * {@code ignoreFirstLine} to be true (that is, this file does have a header) and detects delimiters automatically.
	 */
	public static Hashtable<String, Vector<String>> loadFileToHashVec(String filename, boolean onlyIfAbsent, String key, String... values) {
		String firstLine = Files.getFirstLineOfFile(filename, null);
		String delimiter = ext.determineDelimiter(firstLine);
		String[] header = firstLine.split(delimiter);

		int keyCol = ext.indexOfStr(key, header);
		int[] valueCols = ext.indexFactors(values, header, true, true);

		return loadFileToHashVec(	filename, keyCol, valueCols, delimiter,
		                         	true, onlyIfAbsent);
	}

	/**
	 * As {@link #loadFileToHashVec(String, int[], int[], String, boolean, boolean)} with a single key
	 * index.
	 */
	public static Hashtable<String, Vector<String>> loadFileToHashVec(String filename, int keyIndex,
																																		int[] valueIndices,
																																		String delimiter,
																																		boolean ignoreFirstLine,
																																		boolean onlyIfAbsent) {
		return loadFileToHashVec(	filename, new int[] {keyIndex}, valueIndices, delimiter,
															ignoreFirstLine, onlyIfAbsent);
	}

	/**
	 * Loads a file such that each line is a set of columns, separated by the specified delimiter. An
	 * arbitrary number of key and valu columns may be specified.
	 */
	public static Hashtable<String, Vector<String>> loadFileToHashVec(String filename,
																																		int[] keyIndices,
																																		int[] valueIndices,
																																		String delimiter,
																																		boolean ignoreFirstLine,
																																		boolean onlyIfAbsent) {
		BufferedReader reader = null;
		String[] line;
		Hashtable<String, Vector<String>> hash = new Hashtable<String, Vector<String>>();
		String key, temp;

		try {
			reader = new BufferedReader(new FileReader(filename));
			if (ignoreFirstLine) {
				reader.readLine();
			}
			while (reader.ready()) {
				temp = reader.readLine();
				if (temp.indexOf("\t") == -1) {
					line = temp.trim().split("[\\s]+");
				} else {
					line = temp.split("\t", -1);
				}
				key = "";
				for (int i = 0; i < keyIndices.length; i++) {
					key += ((i == 0) ? "" : delimiter) + line[keyIndices[i]];
				}
				temp = "";
				for (int i = 0; i < valueIndices.length; i++) {
					temp += ((i == 0) ? "" : delimiter) + line[valueIndices[i]];
				}
				addToHashVec(hash, key, temp, onlyIfAbsent);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}

		return hash;
	}

	/**
	 * As {@link #loadFileToHashHashVec(String, int, int, int[], boolean, boolean)} but only the
	 * targetIndex column is stored.
	 */
	public static Hashtable<String, Hashtable<String, String>> loadFileToHashHash(String filename,
																																								int key1Index,
																																								int key2Index,
																																								int targetIndex,
																																								boolean ignoreFirstLine) {
		BufferedReader reader = null;
		String[] line;
		Hashtable<String, Hashtable<String, String>> hashes =
																												new Hashtable<String, Hashtable<String, String>>();
		String temp;

		try {
			reader = new BufferedReader(new FileReader(filename));
			if (ignoreFirstLine) {
				reader.readLine();
			}
			while (reader.ready()) {
				temp = reader.readLine();
				if (temp.indexOf("\t") == -1) {
					line = temp.trim().split("[\\s]+");
				} else {
					line = temp.split("\t", -1);
				}
				addToHashHash(hashes, line[key1Index], line[key2Index], line[targetIndex]);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}

		return hashes;
	}

	/**
	 * A file is read such that each line is an entry with the first two columns as keys.
	 */
	public static Hashtable<String, Hashtable<String, Vector<String>>> loadFileToHashHashVec(	String filename,
																																														int key1Index,
																																														int key2Index,
																																														int[] targetIndices,
																																														boolean ignoreFirstLine,
																																														boolean onlyifabsent) {
		BufferedReader reader = null;
		String[] line;
		Hashtable<String, Hashtable<String, Vector<String>>> hashes =
																																new Hashtable<String, Hashtable<String, Vector<String>>>();
		String temp, str;

		try {
			reader = new BufferedReader(new FileReader(filename));
			if (ignoreFirstLine) {
				reader.readLine();
			}
			while (reader.ready()) {
				temp = reader.readLine();
				if (temp.indexOf("\t") == -1) {
					line = temp.trim().split("[\\s]+");
				} else {
					line = temp.split("\t", -1);
				}
				str = "";
				for (int i = 0; i < targetIndices.length; i++) {
					str += (i == 0 ? "" : "\t") + line[targetIndices[i]];
				}
				addToHashHashVec(hashes, line[key1Index], line[key2Index], str, onlyifabsent);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}

		return hashes;
	}

	/**
	 * Create a vector of String vectors, pre-allocated to the specified size
	 */
	public static Vector<Vector<String>> newVecVecString(int size) {
		Vector<Vector<String>> v = new Vector<Vector<String>>();

		for (int i = 0; i < size; i++) {
			v.add(new Vector<String>());
		}

		return v;
	}

	/**
	 * Create a vector of hashtables, pre-allocated to the specified size.
	 */
	public static Vector<Hashtable<String, String>> newVecHashStringString(int size) {
		Vector<Hashtable<String, String>> hashes = new Vector<Hashtable<String, String>>();

		for (int i = 0; i < size; i++) {
			hashes.add(new Hashtable<String, String>());
		}

		return hashes;
	}

	/**
	 * Copy the specified vector
	 *
	 * @return the new copy
	 */
	public static Vector<String> cloneVectorString(Vector<String> v) {
		Vector<String> vec = new Vector<String>();

		for (String s : v) {
			vec.add(s);
		}

		return vec;
	}

	/**
	 * Lookup the key in the given table. If not found, use missingValue instead.
	 */
	public static String get(Hashtable<String, String> hash, String key, String missingValue) {
		// TODO this would be nicer in a subclass
		if (hash.containsKey(key)) {
			return hash.get(key);
		}
		System.err.println("Error - no valid value in hash for " + key);
		return missingValue;
	}

	/**
	 * @return A sorted array of keys for the given map
	 */
	public static String[] getKeys(Map<String, ?> map) {
		return getKeys(map, true);
	}

	/**
	 * @return An array, optionally sorted, of the keys for the given map
	 */
	public static String[] getKeys(Map<String, ?> map, boolean sort) {
		String[] keys = map.keySet().toArray(new String[map.size()]);
		if (sort) {
			Arrays.sort(keys);
		}
		return keys;
	}

	/**
	 * @return A sorted list of keys for the given map
	 */
	public static List<String> getKeyList(Map<String, ?> map) {
		List<String> keys = new ArrayList<String>(map.keySet());
		Collections.sort(keys);
		return keys;
	}

	/**
	 * Sorting method when keys may contain numeric values. See <a href=
	 * "http://stackoverflow.com/questions/104599/sort-on-a-string-that-may-contain-a-number">this SO
	 * post</a> for information on how these keys are sorted.
	 *
	 * @return A sorted list of keys for the given hashtable, sorted alphanumerically
	 */
	public static List<String> getNumericKeyList(Map<String, ?> map) {
		List<String> keys = new ArrayList<String>(map.keySet());
		Collections.sort(keys, new SciStringComparator());
		return keys;
	}

	/**
	 * As {@link #getNumericKeys(Map)} but returns an array instead
	 */
	public static String[] getNumericKeys(Map<String, ?> map) {
		String[] keys = map.keySet().toArray(new String[map.size()]);
		Arrays.sort(keys, new SciStringComparator());
		return keys;
	}

	public static void main(String...strings) {
		Map<String, String> map = new HashMap<String, String>();
		for (int i=0; i<1000000; i++) {
			map.put(String.valueOf(i), "");
		}

		long t = System.currentTimeMillis();
		getNumericKeys(map);
		t = System.currentTimeMillis() - t;
		System.out.println(t);
	}
}
