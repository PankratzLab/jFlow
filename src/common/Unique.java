package common;

import java.io.*;
import java.util.*;

public class Unique {
	public static void proc(String[] filenames, int[] skips, String[] delimiters, String uniquesFile, String countsFile, boolean noInput) {
		BufferedReader reader;
		PrintWriter writer = null;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		String[] line, keys;
		String temp;
		int count;
		int[] cols;

		if (skips == null) {
			skips = Array.intArray(filenames.length, 0);
		}
		if (delimiters == null) {
			delimiters = Array.stringArray(filenames.length, "[\\s]+");
		}
		if (uniquesFile != null && uniquesFile.equalsIgnoreCase("null")) {
			uniquesFile = null;
		}
		if (countsFile != null && countsFile.equalsIgnoreCase("null")) {
			countsFile = null;
		}

		// loop through to check that all parameters are valid before wasting time processing data
		cols = new int[filenames.length];
		for (int i = 0; i < filenames.length; i++) {
			if (ext.removeDirectoryInfo(filenames[i]).contains(":")) {
				try {
					cols[i] = Integer.parseInt(filenames[i].substring(filenames[i].lastIndexOf(":")+1));
				} catch (NumberFormatException nfe) {
					System.err.println("Error - cannot parse column from '"+filenames[i]+"'");
					return;
				}
				filenames[i] = filenames[i].substring(0, filenames[i].lastIndexOf(":"));
			} else {
				cols[i] = 0;							
				filenames[i] = filenames[i];
			}
			if (!new File(filenames[i]).exists()) {
				System.err.println("Error - file '"+filenames[i]+"' does not exist");
				return;
			}
		}
		
		try {
			if (uniquesFile != null) {
				writer = new PrintWriter(new FileWriter(uniquesFile));
			}
			count = 0;

			for (int i = 0; i < filenames.length; i++) {
				System.out.println("Parsing column index "+cols[i]+" of "+ext.removeDirectoryInfo(filenames[i]));
				
				reader = new BufferedReader(new FileReader(filenames[i]));
				if (skips != null) {
					for (int j = 0; j < skips[i]; j++) {
						reader.readLine();
					}
				}
				while (reader.ready()) {
					temp = reader.readLine();
					if (count == 0 && (double)(new File(filenames[i]).length())/(double)temp.length() > 1000) {
						hash = new Hashtable<String, String>((int)((double)(new File(filenames[i]).length())/(double)temp.length()*1.2));
					}						
					line = temp.trim().split(delimiters[i]);

					if (line.length<=cols[i]) {
						System.err.println("Error - Not enough columns for line:\n"+temp);
						reader.close();
						return;
					}

					if (!hash.containsKey(line[cols[i]])) {
						if (uniquesFile != null) {
							writer.println(line[cols[i]]);
						}
						hash.put(line[cols[i]], "1");
					} else {
						hash.put(line[cols[i]], (Integer.parseInt(hash.get(line[cols[i]]))+1)+"");
					}
					count++;
				};

				reader.close();

			}
			if (uniquesFile != null) {
				writer.close();
			}
			System.out.println("Found "+hash.size()+" unique records among "+count+" total records");			
		} catch (Exception e) {
			System.err.println("Error writing to \""+uniquesFile+"\"");
			e.printStackTrace();
			return;
		}
		
		if (countsFile != null) {
			try {
				writer = new PrintWriter(new FileWriter(countsFile));
				keys = HashVec.getKeys(hash, false, false);
				writer.println(keys.length);
				writer.println();
				for (int i = 0; i<keys.length; i++) {
					writer.println(hash.get(keys[i])+"\t"+keys[i]+"\t"+hash.get(keys[i]));
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to " + countsFile);
				e.printStackTrace();
			}
		}		

		System.out.println("...done");

		if (!noInput) {
			ext.waitForResponse();
		}
	}

	public static void fromParamters(String filename, Logger log) {
		Vector<String> params;
		String[] line, files;
		int[] skips;
		String[] delimiters;
		String out, outCounts;
		
		params = Files.parseControlFile(filename, "unique", new String[] {"unique.out", "uniqueCounts.out", "file1.txt", "file2.txt:3 TAB", "file3.txt:1 skip=9 ,"}, log);
		if (params != null) {
    		out = params.remove(0);
    		outCounts = params.remove(0);
    		files = new String[params.size()];
    		skips = Array.intArray(params.size(), 0);
    		delimiters = Array.stringArray(params.size(), "[\\s]+");
    		for (int i = 0; i < files.length; i++) {
    			line = params.elementAt(i).trim().split("[\\s]+");
    			files[i] = line[0];
    			for (int j = 1; j < line.length; j++) {
        			if (line[j].startsWith("skip=") || line[j].startsWith("skips=")) {
            			skips[i] = Integer.parseInt(line[j].split("=")[1]);
        			} else if (line[j].equals(",")) {
            			delimiters[i] = ",";
        			} else if (line[j].equalsIgnoreCase("tab")) {
            			delimiters[i] = "\t";
        			} else {
        				System.err.println("Error - invalid argument (\""+line[j]+"\") ");
        			}
				}
			}
			proc(files, skips, delimiters, out, outCounts, true);
		}
	}
	
	public static String proc(String[] array) {
		Hashtable<String,String> hash = new Hashtable<String,String>();
	    StringBuilder stringBuilder;
		String[] keys;
	    String ls, temp;
		int count;

		temp = null;
		count = 0;
		for (int i = 0; i<array.length; i++) {
			if (!hash.containsKey(array[i])) {
				hash.put(array[i], "1");
			} else {
				hash.put(array[i], (Integer.parseInt(hash.get(array[i]))+1)+"");
			}
			count++;
		};

		System.out.println("Found "+hash.size()+" unique records among "+count+" total records");
		
		keys = HashVec.getKeys(hash, false, false);
		
	    stringBuilder = new StringBuilder();
	    ls = System.getProperty("line.separator");
        stringBuilder.append(keys.length);
        stringBuilder.append(ls);
        stringBuilder.append(ls);
		for (int i = 0; i<keys.length; i++) {
	        stringBuilder.append(keys[i]+"\t"+hash.get(keys[i]));
	        stringBuilder.append(ls);
		}
		temp = stringBuilder.toString();

		System.out.println("...done");
		return temp;
	}

	public static String[][] proc(String[][] arrays, boolean verbose) {
		Hashtable<String,int[]> hash = new Hashtable<String,int[]>();
		String[] keys;
		int[] counts;
		String[][] allCounts;
		int count;

		count = 0;
		for (int i = 0; i<arrays.length; i++) {
			for (int j = 0; j < arrays[i].length; j++) {
				if (hash.containsKey(arrays[i][j])) {
					counts = hash.get(arrays[i][j]);
				} else {
					counts = new int[arrays.length+1];
				}
				count++;
				counts[0]++;
				counts[i+1]++;
				hash.put(arrays[i][j], counts);
			}
		};

		if (verbose) {
			System.out.println("Found "+hash.size()+" unique records among "+count+" total records found within "+arrays.length+" array"+(arrays.length==1?"":"s"));
		}

		keys = HashVec.getKeys(hash, false, false);
		allCounts = new String[keys.length][arrays.length+2];
		for (int i = 0; i<keys.length; i++) {
			counts = hash.get(keys[i]);
			allCounts[i][0] = keys[i];
			for (int j = 0; j < arrays.length+1; j++) {
				allCounts[i][j+1] = counts[j]+"";
			}
		}

		return allCounts;
	}

	public static void main(String[] args) throws IOException {
		String uniquesFile = null;
		String countsFile = null;
		Vector<String> filenames;

		String usage = "\n"+
		"park.findUnique requires 0-3 arguments\n"+
		"   (1) filenames delimited by a space and including a :# suffix if the column index is not 0 (i.e. file1.txt:4 file2.txt file2.txt:1 (not the default)\n"+
		"   (2) output filename for uniques (i.e. outUniques=unique.out (default for multiple files; default for single files is [filename]-unique.out)\n"+
		"   (3) output filename for counts (i.e. outCounts=uniqueCounts.out (default for multiple files; default for single files is [filename]-uniqueCounts.out)\n"+
		"\n"+
		"   Unique is also available in the Launch.crf variety and includes options for comma/tab-delimited and skipping headers\n"+
		"";

		filenames = new Vector<String>();
		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("outUniques=")) {
				uniquesFile = args[i].split("=")[1];
			} else if (args[i].startsWith("outCounts=")) {
				countsFile = args[i].split("=")[1];
			} else {
				filenames.add(args[i]);
			}
		}
		if (filenames.size() == 0) {
			System.err.println(usage);
			return;
		}
		if (uniquesFile == null) {
			if (filenames.size() == 1) {
				uniquesFile = filenames.elementAt(0)+"-unique.out";
			} else {
				uniquesFile = "unique.out";
			}
		}
		if (countsFile == null) {
			if (filenames.size() == 1) {
				countsFile = filenames.elementAt(0)+"-uniqueCounts.out";
			} else {
				countsFile = "uniqueCounts.out";
			}
		}
		try {
			proc(Array.toStringArray(filenames), null, null, uniquesFile, countsFile, true);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
