package org.genvisis.kaput;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.StringTokenizer;
import java.util.Vector;

public class keepFams {
	public static boolean KEEPFIRSTLINE;

	public keepFams(String struct, String keeps) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		StringTokenizer st;
		String temp, trav;
		Hashtable<String, String> hash = new Hashtable<String, String>();
		Vector<String> duplicateList = new Vector<String>();
		boolean isFirstLine; // workaround for Leah

		if (keeps.startsWith("key=")) {
			hash.put(keeps.substring(4), "null");
		} else {
			reader = new BufferedReader(new FileReader(keeps));
			while (reader.ready()) {
				st = new StringTokenizer(reader.readLine(), " \t\n\r\f,");
				if (st.hasMoreTokens()) {
					trav = st.nextToken();
					if (!hash.containsKey(trav)) {
						hash.put(trav, "null");
						duplicateList.add(trav);
					}
				}
			}
			reader.close();
		}

		isFirstLine = true;

		(new File(struct)).renameTo(new File(struct + ".bak"));
		reader = new BufferedReader(new FileReader(struct + ".bak"));
		writer = new PrintWriter(new FileWriter(struct));
		while (reader.ready()) {
			temp = reader.readLine();
			st = new StringTokenizer(temp, " \t\n\r\f,");
			trav = "";
			if (st.hasMoreTokens()) {
				trav = st.nextToken();
			}
			if (hash.containsKey(trav)) {
				writer.println(temp);
				if (duplicateList.contains(trav)) {
					duplicateList.remove(trav);
				}
			} else if (isFirstLine && KEEPFIRSTLINE) {
				writer.println(temp);
				isFirstLine = false;
			}
		}
		reader.close();
		writer.close();
		if (duplicateList.size() > 0) {
			System.err.println("Warning - the following were found in the list but not in the source:");
			for (int i = 0; i < duplicateList.size(); i++) {
				System.err.println(duplicateList.elementAt(i));
			}
		}

	}

	public static void main(String[] args) throws IOException {
		if (args.length < 2 || args.length > 3) {
			System.out.println("Expecting 2-3 arguments: source_filename and keys_to_keep filename (or \"key=key_you_want_to_keep\").");
			System.out.println("Third optional argument: -keepfirstline");
		} else {
			try {
				if (args.length == 3 && args[2].equals("-keepfirstline")) {
					KEEPFIRSTLINE = true;
				} else {
					KEEPFIRSTLINE = false;
				}
				new keepFams(args[0], args[1]);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
}
