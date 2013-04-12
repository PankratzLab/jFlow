package common;

/**
 * Contains code from the following classes:
 * 			common.Array
 * 			common.ext
 * 			common.Files
 * 			common.HashVec
 * 			common.Sort
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.URI;
import java.net.URL;
import java.text.ParseException;
import java.text.RuleBasedCollator;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;

public class Collapsed {
	public static final int GREATER = 0;
	public static final int LESSTHAN = 1;
	public static final int GREATEROREQUAL = 2;
	public static final int LESSOREQUAL = 3;

	public static final String RULES = "<0<1<2<3<4<5<6<7<8<9"+"<a,A<b,B<c,C<d,D<e,E<f,F<g,G<h,H<i,I<j,J"+"<k,K<l,L<m,M<n,N<o,O<p,P<q,Q<r,R<s,S<t,T"+"<u,U<v,V<w,W<x,X<y,Y<z,Z";

	public static String[] list(String directory, final String suffix, boolean jar) {
		if (jar) {
			try {
//				System.err.println("I haven't been able to get listFiles() to work inside a jar file yet");

				URL repositoryURL = ClassLoader.getSystemResource("common/Files.class");
				String repositoryPath = repositoryURL.getPath();
				URI jarURI = new URI(repositoryPath.substring(0, repositoryPath.indexOf('!')));
				JarFile jarFile = new JarFile(new File(jarURI));
				Vector<String> v = new Vector<String>();

				Enumeration<JarEntry> entries = jarFile.entries();
				while (entries.hasMoreElements()) {
					String entryName = entries.nextElement().getName();

					if (entryName.startsWith(directory)&&entryName.endsWith(suffix)) {
						String trav = entryName.substring(directory.length());
						if (trav.startsWith("/")) {
							trav = trav.substring(1);
						}
						if (!trav.contains("/")) {
							v.add(trav);
						}
					}
				}

				return toStringArray(v);
			} catch (Exception e) {
				System.err.println("Error reading files in jar file");
				e.printStackTrace();
				return null;
			}
		} else {
			String[] files;
			
			files = new File(directory).list(new FilenameFilter() {
				public boolean accept(File file, String filename) {
					return filename.endsWith(suffix);
				}
			});
			
			if (files == null) {
				return new String[0];
			} else {
				return files;
			}
		}
	}

	public static Hashtable<String,String> loadFileToHashString(String filename, int keyIndex, int[] valueIndices, String delimiterWithinHash, boolean ignoreFirstLine) {
		return loadFileToHashString(filename, keyIndex, valueIndices, delimiterWithinHash, ignoreFirstLine, false);
	}
	
	public static Hashtable<String,String> loadFileToHashString(String filename, int keyIndex, int[] valueIndices, String delimiterWithinHash, boolean ignoreFirstLine, boolean jar) {
		return loadFileToHashString(filename, new int[] {keyIndex}, valueIndices, delimiterWithinHash, ignoreFirstLine, jar);
	}
	
	public static Hashtable<String,String> loadFileToHashString(String filename, int[] keyIndices, int[] valueIndices, String delimiterWithinHash, boolean ignoreFirstLine, boolean jar) {
		BufferedReader reader = null;
		String[] line;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		String key, temp;
		
		if (valueIndices == null) {
			valueIndices = new int[] {};
		}

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
					key += ((i==0)?"":delimiterWithinHash)+line[keyIndices[i]];
				}
				temp = "";
				for (int i = 0; i<valueIndices.length; i++) {
					temp += ((i==0)?"":delimiterWithinHash)+line[valueIndices[i]];
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

	public static void writeList(String[] list, String filename) {
        PrintWriter writer;
        
		try {
	        writer = new PrintWriter(new FileWriter(filename));
	        for (int i = 0; i<list.length; i++) {
	        	writer.println(list[i]);
            }
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+filename);
	        e.printStackTrace();
        }
	}
	
	public static boolean exists(String filename, boolean jar) {
		if (jar) {
			try {
				ClassLoader.getSystemResourceAsStream(filename).close();
				return true;
			} catch (Exception e) {
				return false;
			}
		} else {
			return new File(filename).exists();
		}
	}

	public static int[] indexFactors(String[] subset, String[] superset, boolean casesensitive, boolean kill) {
		return indexFactors(subset, superset, casesensitive, new Logger(null), kill);
	}
	
	public static int[] indexFactors(String[] subset, String[] superset, boolean casesensitive, Logger log, boolean kill) {
		int[] indices = new int[subset.length];
		boolean err = false;

		for (int i = 0; i<subset.length; i++) {
			indices[i] = -1;
			for (int j = 0; j<superset.length; j++) {
				if (casesensitive?subset[i].equals(superset[j]):subset[i].equalsIgnoreCase(superset[j])) {
					if (indices[i]==-1) {
						indices[i] = j;
					} else {
						log.reportError("Error - more than one factor was named '"+subset[i]+"'");
						err = true;
					}
				}
			}
			if (indices[i]==-1) {
				log.reportError("Error - no factor was named '"+subset[i]+"'");
				err = true;
			}
		}

		if (kill&&err) {
			System.exit(1);
		}

		return indices;
	}

	public static int indexOfStr(String target, String[] array) {
		return indexOfStr(target, array, true, true);
	}

	public static int indexOfStr(String target, String[] array, boolean caseSensitive, boolean exactMatch) {
		for (int i = 0; i<array.length; i++) {
			if (exactMatch) {
				if (caseSensitive?array[i].equals(target):array[i].toLowerCase().equals(target.toLowerCase())) {
					return i;
				}
			} else {
				if (caseSensitive?array[i].contains(target)||target.contains(array[i]):array[i].toLowerCase().contains(target.toLowerCase())||target.toLowerCase().contains(array[i].toLowerCase())) {
					return i;
				}
			}
		}

		// System.err.println("Error - '"+target+"' was not found in array");
		//		
		return -1;
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

	@SuppressWarnings({ "unchecked" })
	public static String[] getKeys(Hashtable hash) {
		return getKeys(hash, true, false);
	}

	@SuppressWarnings({ "unchecked" })
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

	/**
	 * Prints an array of objects separated by a tab
	 * 
	 * @param array
	 *            an array of objects
	 * @return String of printed objects
	 */
	public static String toStr(String[] array) {
		return toStr(array, "\t");
	}

	/**
	 * Prints an array of objects separated by the specified delimiter
	 * 
	 * @param array
	 *            an array of objects
	 * @param delimiter
	 *            String delimiter
	 * @return String of printed objects
	 */
	public static String toStr(String[] array, String delimiter) {
		String str = "";

		for (int i = 0; i<array.length; i++) {
			str += (i==0?"":delimiter)+array[i];
		}

		return str;
	}

	/**
	 * Creates an array of Strings and copies the contents of a Vector into it
	 * 
	 * @param v
	 *            vector of Strings
	 * @return an array of Strings from the Vector
	 */
	public static String[] toStringArray(Vector<String> v) {
		return v.toArray(new String[v.size()]);
	}
	
	/**
	 * Sorts an array of strings and returns the sorted array
	 * 
	 * @param array
	 *            the array to be sorted
	 * @return the sorted array
	 */
	public static String[] putInOrder(String[] array, boolean treatAsNumbers) {
		int[] keys = quicksort(array);
		String[] newArray = new String[array.length];

		for (int i = 0; i<array.length; i++) {
			newArray[i] = array[keys[i]];
		}

		return newArray;
	}

	/**
	 * Sort an array of Strings into ASCENDING alphabetical order
	 * 
	 * @param arr
	 *            the array of Strings to index
	 * @return an integer array of the sorted indices
	 */
	public static int[] quicksort(String[] arr) {
		int n = arr.length;
		int[] indx = new int[arr.length];
		int i, indxt, ir = n-1, j, k, l = 0;
		int jstack = 0, istack[];
		String a;

		istack = new int[arr.length+1];

		RuleBasedCollator theCollation;
		try {
			theCollation = new RuleBasedCollator(RULES);
			for (j = 0; j<n; j++)
				indx[j] = j;
			for (;;) {
				if (ir-l<7) {
					for (j = l+1; j<=ir; j++) {
						indxt = indx[j];
						a = arr[indxt];
						for (i = j-1; i>=0; i--) {
							if (compare(theCollation, arr[indx[i]], a, LESSOREQUAL))
								break;
							indx[i+1] = indx[i];
						}
						indx[i+1] = indxt;
					}
					if (jstack==0)
						break;
					ir = istack[jstack--];
					l = istack[jstack--];
				} else {
					k = (l+ir)>>1;
					SWAP(indx, k, l+1);
					if (compare(theCollation, arr[indx[l+1]], arr[indx[ir]], GREATER)) {
						SWAP(indx, l+1, ir);
					}
					if (compare(theCollation, arr[indx[l]], arr[indx[ir]], GREATER)) {
						SWAP(indx, l, ir);
					}
					if (compare(theCollation, arr[indx[l+1]], arr[indx[l]], GREATER)) {
						SWAP(indx, l+1, l);
					}
					i = l+1;
					j = ir;
					indxt = indx[l];
					a = arr[indxt];
					for (;;) {
						do
							i++;
						while (compare(theCollation, arr[indx[i]], a, LESSTHAN));
						do
							j--;
						while (compare(theCollation, arr[indx[j]], a, GREATER));
						if (j<i)
							break;
						SWAP(indx, i, j);
					}
					indx[l] = indx[j];
					indx[j] = indxt;
					jstack += 2;
					if (ir-i+1>=j-l) {
						istack[jstack] = ir;
						istack[jstack-1] = i;
						ir = j-1;
					} else {
						istack[jstack] = j-1;
						istack[jstack-1] = l;
						l = i;
					}
				}
			}

		} catch (ParseException pe) {
			System.err.println("Error initializing RuleBasedCollator for quicksort(String[] arr)");
			pe.printStackTrace();
			System.exit(-1);
		}

		return indx;
	}

	/**
	 * Swaps the positions of data in an array
	 * 
	 * @param array
	 *            the array containg the data to swap
	 * @param index1
	 *            position 1 to swap
	 * @param index2
	 *            position 2 to swap
	 */
	private static void SWAP(int[] array, int index1, int index2) {
		int temp = array[index1];
		array[index1] = array[index2];
		array[index2] = temp;
	}

	/**
	 * Rule comparitor to compare 2 strings and return if it matches the
	 * comparison type
	 * 
	 * @param collator
	 *            the rule collator for the comparison
	 * @param arg1
	 *            string 1
	 * @param arg2
	 *            string 2
	 * @param type
	 *            the type of comparison
	 *            GREATER/LESSTHAN/GREATEROREQUAL/LESSOREQUAL
	 * @return the result of the comparison
	 */
	private static boolean compare(RuleBasedCollator collator, String arg1, String arg2, int type) {
		int compareResult = collator.compare(arg1, arg2);
		switch (type) {
		case (GREATER):
			if (compareResult>0) {
				return true;
			}
			break;
		case (LESSTHAN):
			if (compareResult<0) {
				return true;
			}
			break;
		case (GREATEROREQUAL):
			if (compareResult>=0) {
				return true;
			}
			break;
		case (LESSOREQUAL):
			if (compareResult<=0) {
				return true;
			}
			break;
		}
		return false;
	}
}
