package common;

import java.io.*;
import java.util.*;

public class Unique {
	public static void proc(String[] filenames, String uniquesFile, String countsFile, boolean noInput) { // , boolean trashTheRestOfTheLine
		BufferedReader reader = null;
		PrintWriter writer = null;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		String[] line, keys;
		String temp;
		int count;
		String filename;
		int col;

		for (int i = 0; i < filenames.length; i++) {
			if (ext.removeDirectoryInfo(filenames[i]).contains(":")) {
				filename = filenames[i].substring(0, filenames[i].lastIndexOf(":"));
				if (!new File(filename).exists()) {
					System.err.println("Error - file '"+filename+"' does not exist");
					return;
				}
				try {
					col = Integer.parseInt(filenames[i].substring(filenames[i].lastIndexOf(":")+1));
				} catch (NumberFormatException nfe) {
					System.err.println("Error - cannot parse column from '"+filenames[i]+"'");
					return;
				}
			} else {
				filename = filenames[i];
				col = 0;							
			}
		}
		
		try {
			writer = new PrintWriter(new FileWriter(uniquesFile));
			count = 0;

			for (int i = 0; i < filenames.length; i++) {
				if (ext.removeDirectoryInfo(filenames[i]).contains(":")) {
					filename = filenames[i].substring(0, filenames[i].lastIndexOf(":"));
					try {
						col = Integer.parseInt(filenames[i].substring(filenames[i].lastIndexOf(":")+1));
					} catch (NumberFormatException nfe) {
						System.err.println("Error - cannot parse column from '"+filenames[i]+"'");
						return;
					}
				} else {
					filename = filenames[i];
					col = 0;							
				}
				System.out.println("Parsing column index "+col+" of "+ext.removeDirectoryInfo(filename));
				
				reader = new BufferedReader(new FileReader(filename));
				while (reader.ready()) {
					temp = reader.readLine();
					if (count == 0 && (double)(new File(filename).length())/(double)temp.length() > 1000) {
						hash = new Hashtable<String, String>((int)((double)(new File(filename).length())/(double)temp.length()*1.2));
					}						
					line = temp.trim().split("[\\s]+");

					if (line.length<=col) {
						System.err.println("Error - Not enough columns for line:\n"+temp);
						return;
					}

					if (!hash.containsKey(line[col])) {
//						if (trashTheRestOfTheLine) {
							writer.println(line[col]);
////							writer.flush();
//						} else {
//							writer.println(temp);
////							writer.flush();
//						}
						hash.put(line[col], "1");
					} else {
						hash.put(line[col], (Integer.parseInt(hash.get(line[col]))+1)+"");
					}
					count++;
				};

				reader.close();

			}
			writer.close();
			System.out.println("Found "+hash.size()+" unique records among "+count+" total records");			
		} catch (Exception e) {
			System.err.println("Error writing to \""+uniquesFile+"\"");
			e.printStackTrace();
			return;
		}
		
		if (filenames.length == 1) {
			filename = filenames[0]+"-uniqueCount.out";
		} else {
			filename = "unique.out";
		}
		try {
			writer = new PrintWriter(new FileWriter(filename));
			keys = HashVec.getKeys(hash, false, false);
			writer.println(keys.length);
			writer.println();
			for (int i = 0; i<keys.length; i++) {
				writer.println(hash.get(keys[i])+"\t"+keys[i]+"\t"+hash.get(keys[i]));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + filename);
			e.printStackTrace();
		}

		System.out.println("...done");

		if (!noInput) {
			CmdLine.stdin();
		}
	}

	public static String proc(String[] array) {
		Hashtable<String,String> hash = new Hashtable<String,String>();
		String[] keys;
		String temp;
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
		temp = keys.length +"\n\n";
		for (int i = 0; i<keys.length; i++) {
			temp += keys[i]+"\t"+hash.get(keys[i])+"\n";
		}
		System.out.println("...done");

		return temp;
	}

	public static void main(String[] args) throws IOException {
//		int numArgs = args.length;
		String uniquesFile = null;
		String countsFile = null;
		Vector<String> filenames;
//		boolean trash = false;
		boolean noInput = false;

		String usage = "\n"+
		"park.findUnique requires 0-1 arguments\n"+
		"   (1) filenames delimited by a space and including a :# suffix if the column index is not 0 (i.e. file1.txt:4 file2.txt file2.txt:1 (not the default)\n"+
		"   (2) output filename for uniques (i.e. outUniques=unique.out (default for multiple files; default for single files is [filename]-unique.out)\n"+
		"   (3) output filename for counts (i.e. outCounts=uniqueCounts.out (default for multiple files; default for single files is [filename]-uniqueCounts.out)\n"+
//		"   (optional) trash the rest of the line (i.e. -trash (not the default)\n"+
		"   (optional) don't wait for a keyboard response before closing (i.e. -noInput (not the default)\n"+
		"";

		filenames = new Vector<String>();
		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("outUniques=")) {
				uniquesFile = args[i].split("=")[1];
//				numArgs--;
			} else if (args[i].startsWith("outCounts=")) {
				countsFile = args[i].split("=")[1];
//				numArgs--;
//			} else if (args[i].equalsIgnoreCase("-trash")) {
//				trash = true;
//				numArgs--;
			} else if (args[i].equalsIgnoreCase("-noInput")) {
				noInput = true;
//				numArgs--;
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
			proc(Array.toStringArray(filenames), uniquesFile, countsFile, noInput);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
