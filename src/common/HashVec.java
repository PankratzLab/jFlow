package common;

import java.io.*;
import java.text.*;
import java.util.*;

public class HashVec {
	public static boolean addIfAbsent(String str, Vector<String> v) {
		return addIfAbsent(str, v, false);
	}

	public static boolean addIfAbsent(String str, Vector<String> v, boolean inAlphabeticalOrder) {
		RuleBasedCollator collator;

		if (v.contains(str)) {
			return false;
		} else {
			if (inAlphabeticalOrder) {
				try {
					collator = new RuleBasedCollator(Sort.RULES);
					for (int i = 0; i<=v.size(); i++) {
						if (i==v.size()||collator.compare(v.elementAt(i), str)>0) {
							v.insertElementAt(str, i);
							return true;
						}
					}
				} catch (ParseException pe) {}
			} else {
				v.add(str);
			}
			return true;
		}
	}

	@SuppressWarnings({ "rawtypes" })
	public static String[] getKeys(Hashtable hash) {
		return getKeys(hash, true, false);
	}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	public static String[] getKeys(Hashtable hash, boolean sort, boolean treatAsNumbers) {
		String[] array = new String[hash.size()];
		Enumeration<String> enumer = hash.keys();
		int count = 0;

		try {
			while (enumer.hasMoreElements()) {
				array[count++] = enumer.nextElement();
			}
		} catch (Exception e) {
			System.err.println("Error - hash keys were not Strings");
			e.printStackTrace();
		}

		if (sort) {
			return Sort.putInOrder(array, treatAsNumbers);
		} else {
			return array;
		}
	}

	public static Object addToHashIfAbsent(Hashtable<String,Object> hash, String key, Object value) {
		Object o;

		if (hash.containsKey(key)) {
			o = hash.get(key);
		} else {
			hash.put(key, o = value);
		}

		return o;
	}

	public static void addToHashIfAbsent(Hashtable<String,String> hash, String key, String value) {
		if (!hash.containsKey(key)) {
			hash.put(key, value);
		}
	}

	public static void addToHashIfAbsent(Hashtable<String,String> hash, String[] keys, String[] values) {
		for (int i = 0; i < keys.length; i++) {
			if (!hash.containsKey(keys[i])) {
				hash.put(keys[i], values==null?"":values[i]);
			}
		}
	}

	public static void addToHashVec(Hashtable<String,Vector<String>> hash, String key, String value, boolean onlyifabsent) {
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

	public static void addToHashArrayVec(Hashtable<String,Vector<String[]>> hash, String key, String[] array) {
		Vector<String[]> v;

		if (hash.containsKey(key)) {
			v = hash.get(key);
		} else {
			hash.put(key, v = new Vector<String[]>());
		}

		v.add(array);
	}

	public static void addToHashHash(Hashtable<String,Hashtable<String,Object>> hash1, String key1, String key2, Object value) {
		Hashtable<String,Object> hash2;

		if (hash1.containsKey(key1)) {
			hash2 = hash1.get(key1);
		} else {
			hash1.put(key1, hash2 = new Hashtable<String,Object>());
		}
		hash2.put(key2, value);
	}

	public static void addToHashHash(Hashtable<String,Hashtable<String,String>> hash1, String key1, String key2, String value) {
		Hashtable<String,String> hash2;

		if (hash1.containsKey(key1)) {
			hash2 = hash1.get(key1);
		} else {
			hash1.put(key1, hash2 = new Hashtable<String,String>());
		}
		hash2.put(key2, value);
	}

	public static void addToHashHashVec(Hashtable<String,Hashtable<String,Vector<String>>> hash1, String key1, String key2, String value, boolean addIfAbsent) {
		Hashtable<String,Vector<String>> hash2;
		Vector<String> v;
		
		if (hash1.containsKey(key1)) {
			hash2 = hash1.get(key1);
		} else {
			hash1.put(key1, hash2 = new Hashtable<String,Vector<String>>());
		}
		if (hash2.containsKey(key2)) {
			v = hash2.get(key2);
		} else {
			hash2.put(key2, v = new Vector<String>());
		}
		if (addIfAbsent) {
			addIfAbsent(value, v);
		} else {
			v.add(value);
		}
	}

	public static Hashtable<String,String> loadFileToHashString(String filename, boolean ignoreFirstLine) {
		return loadFileToHashString(filename, 0, new int[] {1}, null, ignoreFirstLine);
	}

	public static Hashtable<String,String> loadFileToHashNull(String filename, boolean ignoreFirstLine) {
		return loadFileToHashString(filename, 0, null, null, ignoreFirstLine);
	}

	public static Hashtable<String,String> loadFileToHashNull(String[] list) {
		Hashtable<String,String> hash = new Hashtable<String,String>((list == null?10:list.length));
		
		for (int i = 0; list != null && i<list.length; i++) {
			hash.put(list[i], "");
		}

		return hash;
	}

	public static Hashtable<String,Integer> loadFileToHashIndex(String[] list) {
		Hashtable<String,Integer> hash = new Hashtable<String,Integer>((list == null?10:list.length));
		
		for (int i = 0; list != null && i<list.length; i++) {
			hash.put(list[i], i);
		}

		return hash;
	}

	public static Vector<String> loadFileToVec(String filename, boolean ignoreFirstLine, boolean onlyFirstColumn, boolean onlyIfAbsent) {
		return loadFileToVec(filename, ignoreFirstLine, onlyFirstColumn?new int[] {0}:null, onlyIfAbsent, false);
	}

	public static Vector<String> loadFileToVec(String filename, boolean ignoreFirstLine, boolean onlyFirstColumn, boolean onlyIfAbsent, boolean jar) {
		return loadFileToVec(filename, ignoreFirstLine, onlyFirstColumn?new int[] {0}:null, onlyIfAbsent, jar);
	}

	public static String[] loadFileToStringArray(String filename, boolean ignoreFirstLine, int[] cols, boolean onlyIfAbsent) {
		return Array.toStringArray(loadFileToVec(filename, ignoreFirstLine, cols, onlyIfAbsent, false));
	}

	public static String[] loadFileToStringArray(String filename, boolean jar, boolean ignoreFirstLine, int[] cols, boolean onlyIfAbsent) {
		return Array.toStringArray(loadFileToVec(filename, ignoreFirstLine, cols, onlyIfAbsent, jar));
	}

	public static String[] loadFileToStringArray(String filename, boolean jar, boolean ignoreFirstLine, int[] cols, boolean onlyIfAbsent, String delimiter) {
		return Array.toStringArray(loadFileToVec(filename, ignoreFirstLine, cols, onlyIfAbsent, jar, delimiter));
	}

	public static Vector<String> loadFileToVec(String filename, boolean ignoreFirstLine, int[] cols, boolean onlyIfAbsent, boolean jar) {
		return loadFileToVec(filename, ignoreFirstLine, cols, onlyIfAbsent, jar, "[\\s]+");
	}

	public static Vector<String> loadFileToVec(String filename, boolean ignoreFirstLine, int[] cols, boolean onlyIfAbsent, boolean jar, String delimiter) {
		BufferedReader reader = null;
		Vector<String> v = new Vector<String>();
		String trav;
		String[] line;
		int count;

		try {
			reader = Files.getReader(filename, jar, true, false);
			if (ignoreFirstLine) {
				reader.readLine();
			}
			count = 1;
			while (reader.ready()) {
				trav = reader.readLine();
				if (cols!=null) {
					line = trav.trim().split(delimiter, -1);
					trav = "";
					for (int i = 0; i<cols.length; i++) {
						if (line.length <= cols[i]) {
							System.err.println("Error - not enough columns at line "+count+": "+Array.toStr(line));
						}
						trav += (i==0?"":"\t")+line[cols[i]];
					}
				}
				if (onlyIfAbsent) {
					addIfAbsent(trav, v);
				} else {
					v.add(trav);
				}
				count++;
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}

		return v;
	}

	public static String[][] loadFileToStringMatrix(String filename, boolean ignoreFirstLine, int[] cols, boolean jar) {
		return loadFileToStringMatrix(filename, ignoreFirstLine, cols, "[\\s]+", jar, 1000, false);
	}
	
	public static String[][] loadFileToStringMatrix(String filename, boolean ignoreFirstLine, int[] cols, String delimiter, boolean jar, int initialCapacity, boolean allowMissingData) {
		BufferedReader reader = null;
		Vector<String[]> v = new Vector<String[]>(initialCapacity);
		String[] line, data;

		try {
			reader = Files.getReader(filename, jar, true, false);
			if (ignoreFirstLine) {
				reader.readLine();
			}
			while (reader.ready()) {
				line = reader.readLine().trim().split(delimiter);
				if (cols == null) {
					v.add(line);
				} else {
					data = new String[cols.length];
					for (int i = 0; i<cols.length; i++) {
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
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}

		return Matrix.toStringArrays(v);
	}

	public static Hashtable<String,String> loadFileToHashString(String filename, String keyHeader, String[] valueHeaders, String delimiterWithinHash) {
		BufferedReader reader = null;
		String[] line;
		int keyIndex;
		int[] valueIndices;
		
		try {
	        reader = new BufferedReader(new FileReader(filename));
	        line = reader.readLine().trim().split("\t", -1);
	        keyIndex = ext.indexOfStr(keyHeader, line);
	        if (keyIndex == -1) {
	        	System.err.println("Error - '"+keyHeader+"' not found in "+filename);
	        	System.exit(1);
	        }
	        valueIndices = ext.indexFactors(valueHeaders, line, false, true);
	        
	        reader.close();
	        return loadFileToHashString(filename, keyIndex, valueIndices, delimiterWithinHash, true);
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+filename+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+filename+"\"");
	        System.exit(2);
        }
        
		return null;
	}
	
	public static Hashtable<String,String> loadFileToHashString(String filename, int keyIndex, int[] valueIndices, String delimiterWithinHash, boolean ignoreFirstLine) {
		return loadFileToHashString(filename, keyIndex, valueIndices, delimiterWithinHash, ignoreFirstLine, false);
	}
	
	public static Hashtable<String,String> loadFileToHashString(String filename, int keyIndex, int[] valueIndices, String delimiterWithinHash, boolean ignoreFirstLine, boolean jar) {
		return loadFileToHashString(filename, new int[] {keyIndex}, valueIndices, false, delimiterWithinHash, ignoreFirstLine, jar, false);
	}
	
	public static Hashtable<String,String> loadFileToHashString(String filename, int[] keyIndices, int[] valueIndices, boolean commaDelimitedFile, String delimiterWithinHash, boolean ignoreFirstLine, boolean jar, boolean allowMissingData) {
		BufferedReader reader = null;
		String[] line;
		Hashtable<String,String> hash = new Hashtable<String,String>();
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
				if (commaDelimitedFile) {
					line = temp.split(",");
				} else if (temp.indexOf("\t")==-1) {
					line = temp.trim().split("[\\s]+");
				} else {
					line = temp.split("\t", -1);
				}
				key = "";
				for (int i = 0; i<keyIndices.length; i++) {
					key += ((i==0)?"":delimiterWithinHash)+line[keyIndices[i]];
				}
				temp = "";
				for (int i = 0; i<valueIndices.length; i++) {
					if (valueIndices[i] == -7) {
						temp += ((i==0)?"":delimiterWithinHash)+hash.size();
					} else if (valueIndices[i]<line.length) {
						temp += ((i==0)?"":delimiterWithinHash)+line[valueIndices[i]];
					} else if (allowMissingData) {
						temp += ((i==0)?"":delimiterWithinHash)+".";
					} else {
						System.err.println("Error - not enough columns for key '"+key+"'; and allMissingData was not flagged");
					}
				}
				hash.put(key, temp);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}

		return hash;
	}

	public static Hashtable<String,Vector<String>> loadFileToHashVec(String filename, int keyIndex, int[] valueIndices, String delimiter, boolean ignoreFirstLine, boolean onlyIfAbsent) {
		return loadFileToHashVec(filename, new int[] {keyIndex}, valueIndices, delimiter, ignoreFirstLine, onlyIfAbsent);
	}
	
	public static Hashtable<String,Vector<String>> loadFileToHashVec(String filename, int[] keyIndices, int[] valueIndices, String delimiter, boolean ignoreFirstLine, boolean onlyIfAbsent) {
		BufferedReader reader = null;
		String[] line;
		Hashtable<String,Vector<String>> hash = new Hashtable<String,Vector<String>>();
		String key, temp;

		try {
			reader = new BufferedReader(new FileReader(filename));
			if (ignoreFirstLine) {
				reader.readLine();
			}
			while (reader.ready()) {
				temp = reader.readLine();
				if (temp.indexOf("\t")==-1) {
					line = temp.trim().split("[\\s]+");
				} else {
					line = temp.split("\t", -1);
				}
				key = "";
				for (int i = 0; i<keyIndices.length; i++) {
					key += ((i==0)?"":delimiter)+line[keyIndices[i]];
				}
				temp = "";
				for (int i = 0; i<valueIndices.length; i++) {
					temp += ((i==0)?"":delimiter)+line[valueIndices[i]];
				}
				addToHashVec(hash, key, temp, onlyIfAbsent);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}

		return hash;
	}

	public static Hashtable<String,Hashtable<String,String>> loadFileToHashHash(String filename, int key1Index, int key2Index, int targetIndex, boolean ignoreFirstLine) {
		BufferedReader reader = null;
		String[] line;
		Hashtable<String,Hashtable<String,String>> hashes = new Hashtable<String,Hashtable<String,String>>();
		String temp;

		try {
			reader = new BufferedReader(new FileReader(filename));
			if (ignoreFirstLine) {
				reader.readLine();
			}
			while (reader.ready()) {
				temp = reader.readLine();
				if (temp.indexOf("\t")==-1) {
					line = temp.trim().split("[\\s]+");
				} else {
					line = temp.split("\t", -1);
				}
				addToHashHash(hashes, line[key1Index], line[key2Index], line[targetIndex]);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}

		return hashes;
	}

	public static Hashtable<String,Hashtable<String,Vector<String>>> loadFileToHashHashVec(String filename, int key1Index, int key2Index, int[] targetIndices, String delimiter, boolean ignoreFirstLine, boolean addIfAbsent) {
		BufferedReader reader = null;
		String[] line;
		Hashtable<String,Hashtable<String,Vector<String>>> hashes = new Hashtable<String,Hashtable<String,Vector<String>>>();
		String temp, str;

		try {
			reader = new BufferedReader(new FileReader(filename));
			if (ignoreFirstLine) {
				reader.readLine();
			}
			while (reader.ready()) {
				temp = reader.readLine();
				if (temp.indexOf("\t")==-1) {
					line = temp.trim().split("[\\s]+");
				} else {
					line = temp.split("\t", -1);
				}
				str = "";
				for (int i = 0; i<targetIndices.length; i++) {
					str += (i==0?"":"\t")+line[targetIndices[i]];
                }
				addToHashHashVec(hashes, line[key1Index], line[key2Index], str, addIfAbsent);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}

		return hashes;
	}

	public static Vector<Vector<String>> newVecVecString(int size) {
		Vector<Vector<String>> v = new Vector<Vector<String>>();

		for (int i = 0; i<size; i++) {
			v.add(new Vector<String>());
		}

		return v;
	}

	public static Vector<Hashtable<String,String>> newVecHashStringString(int size) {
		Vector<Hashtable<String,String>> hashes = new Vector<Hashtable<String,String>>();

		for (int i = 0; i<size; i++) {
			hashes.add(new Hashtable<String,String>());
		}

		return hashes;
	}

	public static Vector<String> cloneVectorString(Vector<String> v) {
		Vector<String> vec = new Vector<String>();

		for (int i = 0; i<v.size(); i++) {
			vec.add(v.elementAt(i));
		}

		return vec;
	}

	public static void addToVector(Vector<String> container, Vector<String> additions) {
		for (int i = 0; i<additions.size(); i++) {
			container.add(additions.elementAt(i));
		}
	}

	public static String get(Hashtable<String,String> hash, String key, String missingValue) {
		if (hash.containsKey(key)) {
			return hash.get(key);
		} else {
			System.err.println("Error - no valid value in hash for "+key);
			return missingValue;
		}
	}
	
	public static void mergeHash2IntoHash1(Hashtable<String,String> hash1, Hashtable<String,String> hash2) {
		String[] keys;
		
		keys = getKeys(hash2, true, false);
		for (int i = 0; i < keys.length; i++) {
			hash1.put(keys[i], hash2.get(keys[i]));
		}
	}
	
}
