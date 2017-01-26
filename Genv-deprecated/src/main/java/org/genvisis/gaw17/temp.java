package org.genvisis.gaw17;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;

public class temp {
	public static void run(String filename) {
		PrintWriter writer;
		String[] keys;
		Hashtable<String, Vector<String>> hash;

		hash = HashVec.loadFileToHashVec(filename, 0, new int[] {1}, "", false, false);
		keys = HashVec.getKeys(hash);
		try {
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename) + "_genes.dat"));
			for (String key : keys) {
				writer.println(key + "\t" + ArrayUtils.toStr(ArrayUtils.toStringArray(hash.get(key))));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.rootOf(filename) + "_genes.dat");
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "Aff.txt";

		String usage = "\n"	+ "gaw17.temp requires 0-1 arguments\n" + "   (1) filename (i.e. file="
										+ filename + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
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
