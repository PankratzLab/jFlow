package common;

import java.io.*;
import java.util.*;

public class FilterByLists {
	public static void process(String filename, String keeps, String deletes, int col, String outfile, boolean keepFirstLine, Logger log) {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;
		String temp;
		Hashtable<String,String> keepsHash, deletesHash;
		Vector<String> vKeeps, vDeletes;
		boolean isFirstLine, problem;

		keepsHash = new Hashtable<String,String>();
		vKeeps = new Vector<String>();
		if (keeps != null) {
			if (!new File(keeps).exists()) {
				System.err.println("Since '"+keeps+"' is not a filename, assuming this is the key to keep");
				vKeeps.add(keeps);
			} else {
				vKeeps = HashVec.loadFileToVec(keeps, false, true, true);
			}	
			keepsHash = HashVec.loadFileToHashNull(Array.toStringArray(vKeeps));
		}

		deletesHash = new Hashtable<String,String>();
		vDeletes = new Vector<String>();
		if (deletes != null) {
			if (!new File(deletes).exists()) {
				System.err.println("Since '"+deletes+"' is not a filename, assuming this is the key to delete");
				vDeletes.add(deletes);
			} else {
				vDeletes = HashVec.loadFileToVec(deletes, false, true, true);
			}
			deletesHash = HashVec.loadFileToHashNull(Array.toStringArray(vDeletes));
		}
		
		problem = false;
		for (int i = 0; i < vDeletes.size(); i++) {
			if (keepsHash.containsKey(vDeletes.elementAt(i))) {
				System.err.println("Error - token '"+vDeletes.elementAt(i)+"' was found in both the keeps file and the deletes file");
				problem = true;
			}
		}
		if (problem) {
			return;
		}

		isFirstLine = true;
		if (outfile == null) {
			new File(filename).renameTo(new File(filename+".bak"));
			outfile = filename;
			filename = filename+".bak";
		}
		
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(outfile));
			while (reader.ready()) {
				temp = reader.readLine();
				line = temp.trim().split("[\\s]+");
				if (isFirstLine&&keepFirstLine) {
					writer.println(temp);
					isFirstLine = false;
				} else if (deletes != null && !deletesHash.containsKey(line[col])) {
					if (vDeletes.contains(line[col])) {
						vDeletes.remove(line[col]);
					}
				} else if (keeps == null || keepsHash.containsKey(line[col])) {
					writer.println(temp);
					if (vKeeps.contains(line[col])) {
						vKeeps.remove(line[col]);
					}
				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
		
		if (log.getLevel() > 0) {
			if (vKeeps.size()>0) {
				System.err.println("Warning - the following were found in the keeps list but not in the data file:");
				for (int i = 0; i<vKeeps.size(); i++) {
					System.err.println(vKeeps.elementAt(i));
				}
			}
	
			if (vDeletes.size()>0) {
				System.err.println("Warning - the following were found in the deletes list but not in the data file:");
				for (int i = 0; i<vDeletes.size(); i++) {
					System.err.println(vDeletes.elementAt(i));
				}
			}
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String keepsFile = null;
		String deletesFile = null;
		String outfile = null;
		boolean keepFirstLine = false;
		int col = 0;
		
		String usage = "\n" + 
		"common.FilterByLists requires 0-1 arguments\n" + 
		"   (1) input filename (i.e. file=" + filename + " (default))\n" + 
		"   (2) filename containing list of tokens to keep (i.e. keeps=keeps.txt (not the default; if token is not a filename, assuming it is the key to keep))\n" + 
		"   (3) filename containing list of tokens to remove (i.e. deletes=deltetes.txt (not the default; if token is not a filename, assuming it is the key to delete))\n" + 
		"   (4) index of the column of the token to match on (i.e. col="+col+" (default))\n" + 
		"   (5) output filename (i.e. out=output.out (not the default; default is to backup and replace the input file))\n" + 
		"   (6) always keep first line (i.e. -keepFirst (not the default))\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("keeps=")) {
				keepsFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("deletes=")) {
				deletesFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("col=")) {
				col = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				outfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-keepFirst")) {
				keepFirstLine = true;
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			process(filename, keepsFile, deletesFile, col, outfile, keepFirstLine, new Logger());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
