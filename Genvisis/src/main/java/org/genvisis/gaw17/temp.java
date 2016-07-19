package org.genvisis.gaw17;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;

public class temp {
	public static void run(String filename) {
		PrintWriter writer;
		String[] keys;
		Hashtable<String, Vector<String>> hash;
		
		hash = HashVec.loadFileToHashVec(filename, 0, new int[] {1}, "", false, false);
		keys = HashVec.getKeys(hash);
		try {
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename)+"_genes.dat"));
			for (int i = 0; i < keys.length; i++) {
				writer.println(keys[i]+"\t"+Array.toStr(Array.toStringArray(hash.get(keys[i]))));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.rootOf(filename)+"_genes.dat");
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "Aff.txt";

		String usage = "\n" + "gaw17.temp requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			run(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
