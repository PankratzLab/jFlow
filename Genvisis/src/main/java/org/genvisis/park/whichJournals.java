package org.genvisis.park;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;

public class whichJournals {
	public static final String IMPACT_DB = "impact_factors.dat";
	public static final String JOURNAL_KEY = "journal_key.dat";
	public static final String DB_FILE = "articles.db.dat";

	public whichJournals(String filename) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line, js;
		String temp, trav, journal;
		Hashtable<String,String> jk, hash;
		Hashtable<String,String[]> impacts;
		Vector<String> v;

		jk = new Hashtable<String,String>();
		try {
			for (reader = new BufferedReader(new FileReader(JOURNAL_KEY)); reader.ready(); jk.put(line[1].toLowerCase(), line[0].toLowerCase())) {
				line = reader.readLine().split("\\|");
			}
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+JOURNAL_KEY+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+JOURNAL_KEY+"\"");
			System.exit(2);
		}

		impacts = new Hashtable<String,String[]>();
		try {
			reader = new BufferedReader(new FileReader(IMPACT_DB));
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				impacts.put(line[0].toLowerCase(), line);
				if (jk.containsKey(line[0].toLowerCase())) {
					impacts.put(jk.get(line[0].toLowerCase()), line);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+IMPACT_DB+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+IMPACT_DB+"\"");
			System.exit(2);
		}

		v = new Vector<String>();
		try {
			for (reader = new BufferedReader(new FileReader(filename)); reader.ready();) {
				temp = reader.readLine();
				if (temp.indexOf("|")>0) {
					line = temp.split("\\|");
					for (int i = 0; i<line.length; i++) {
						v.add(line[i]);
					}
				} else {
					v.add(temp);
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

		hash = new Hashtable<String,String>();
		try {
			reader = new BufferedReader(new FileReader(DB_FILE));
			while (reader.ready()) {
				trav = reader.readLine();
				reader.readLine();
				reader.readLine();
				journal = reader.readLine().split("\\|")[0].toLowerCase();
				reader.readLine();
				reader.readLine();
				reader.readLine();
				reader.readLine();
				reader.readLine();
				if (v.contains(trav)) {
					if (hash.containsKey(journal)) {
						hash.put(journal, hash.get(journal)+"\t"+trav);
					} else {
						hash.put(journal, trav);
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+DB_FILE+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+DB_FILE+"\"");
			System.exit(2);
		}

		writer = new PrintWriter(new FileWriter(filename+"-out.xls"));
		writer.println("Journal\timpact factor\ttimes used");
		js = HashVec.getKeys(hash);
		for (int i = 0; i<js.length; i++) {
			line = (impacts.containsKey(js[i])?impacts.get(js[i]):new String[] {"?"});
			writer.println(js[i]+"\t"+line[line.length-1]+"\t"+(hash.get(js[i])).split("[\\s]+").length+"\t"+hash.get(js[i]));
		}
		writer.close();

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		// String filename = "depression_papers.txt";
		String filename = "REP1_papers.txt";

		String usage = "\n"+"park.whichJournals requires 0-1 arguments\n"+"   (1) filename (i.e. file="+filename+" (default)\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			new whichJournals(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
