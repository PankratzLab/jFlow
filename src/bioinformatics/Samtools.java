package bioinformatics;

import java.io.*;
import java.util.*;
import java.util.Map.Entry;

import common.*;

public class Samtools {

	public static void extractRegions(String filename, int windowInBp, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp, trav;
		Hashtable<String, Vector<String>> hash = new Hashtable<String, Vector<String>>();
		Set<String> keySet;
		int count;
		long time;
		String key;

		hash = loadFromFile(filename, log);
		keySet = hash.keySet();
		for (String entry : keySet) {
			
		}
//		writer = new PrintWriter(new FileWriter(filename + ".IGV_Script"));
//		writer.println("samtools view " + bamFilenames[i] + "-b chr" + chr + ":" + (pos - window) + "-" + (pos + window));
	}
	
	public static Hashtable<String, Vector<String>> loadFromFile(String filename, Logger log) {
		BufferedReader reader;
		String[] line;
		Hashtable<String, Vector<String>> hash;
		Vector<String> v;

		hash = new Hashtable<String, Vector<String>>();
		try {
			reader = new BufferedReader(new FileReader(filename));
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				if(hash.containsKey(line[2])) {
					v = hash.get(line[2]);
				} else {
					v = new Vector<String>();
				}
				v.add(line[4]);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return hash;
	}
	
	public static Hashtable<String, String> writerToFile(String filename, Logger log) {
		BufferedReader reader;
		String[] line;
		String temp, trav;
		Hashtable<String, String> hash = new Hashtable<String, String>();
		Vector<String> v = new Vector<String>();
		int count;
		long time;

		hash = new Hashtable<String, String>();
		try {
			reader = new BufferedReader(new FileReader(filename));
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				hash.put(line[0], line[1]);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return hash;
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "Samtools.dat";
		int windowInBp = 5000;
		String logfile = null;
		Logger log;

		String usage = "\n" + "bioinformatics.Samtools requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("window=")) {
				windowInBp = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = ext.parseStringArg(args[i], null);
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
			log = new Logger(logfile);
			extractRegions(filename, windowInBp, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
