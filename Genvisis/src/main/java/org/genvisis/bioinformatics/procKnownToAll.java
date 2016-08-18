package org.genvisis.bioinformatics;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;

public class procKnownToAll {
	public static final String[] DEFAULT_FILES = {"knownToGNfAtlas2.prn", "knownToGnf1h.prn", "knownToU133.prn", "knownToU133Plus2.prn", "knownToU95.prn"};

	public procKnownToAll(String[] filenames) {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line, genes;
		Hashtable<String,Vector<String>> hash = new Hashtable<String,Vector<String>>();
		Vector<String> v = new Vector<String>();
		for (int i = 0; i<filenames.length; i++) {
			try {
				reader = new BufferedReader(new FileReader(filenames[i]));
				reader.readLine();
				while (reader.ready()) {
					line = reader.readLine().split("[\\s]+");
					HashVec.addToHashVec(hash, line[0], line[1], true);
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \""+filenames[i]+"\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \""+filenames[i]+"\"");
				System.exit(2);
			}
		}
		try {
			genes = HashVec.getKeys(hash);
			writer = new PrintWriter(new FileWriter("knownToAll.prn"));
			writer.println("#name\tvalue");
			for (int i = 0; i<genes.length; i++) {
				v = hash.get(genes[i]);
				for (int j = 0; j<v.size(); j++) {
					writer.println(genes[i]+"\t"+v.elementAt(j));
				}
			}
			writer.close();
		} catch (IOException ioe) {
			System.err.println("Error writing "+"knownToAll.prn");
			ioe.printStackTrace();
		}

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		Vector<String> v = null;

		String usage = "\n"+"park.procKnownToAll requires 0+ arguments\n"+"   filenames (i.e. "+Array.toStr(DEFAULT_FILES, " ")+" (default))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else {
				if (v==null) {
					v = new Vector<String>();
				}
				v.add(args[i]);
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			new procKnownToAll(v==null?DEFAULT_FILES:Array.toStringArray(v));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
