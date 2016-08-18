package org.genvisis.widgets;

import java.io.*;

import org.genvisis.common.*;

public class TabVersion {
	public static void make(String filename) {
		BufferedReader reader;
		PrintWriter writer;
		
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(filename+".xln"));
			while (reader.ready()) {
				writer.println(Array.toStr(reader.readLine().trim().split("[\\s]+")));
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
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "MakeTabVersion.dat";

		String usage = "\n"+
		"widgets.TabVersion requires 0-1 arguments\n"+
		"   (1) filename (i.e. file=" + filename + " (default))\n"+
		"";

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
			make(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
