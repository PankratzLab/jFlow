package org.genvisis.park;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;

public class removeFams {

	public removeFams(String struct, String deletes) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		StringTokenizer st;
		String temp;
		Hashtable<String,String> hash = new Hashtable<String,String>();

		reader = new BufferedReader(new FileReader(deletes));
		while (reader.ready()) {
			st = new StringTokenizer(reader.readLine(), " \t\n\r\f,");
			hash.put(st.nextToken(), "null");
		}
		reader.close();

		String bakFilename = Files.getBakFilename(struct, super.getClass().getName());
		(new File(struct)).renameTo(new File(bakFilename));
		reader = new BufferedReader(new FileReader(bakFilename));
		writer = new PrintWriter(new FileWriter(struct));
		while (reader.ready()) {
			temp = reader.readLine();
			st = new StringTokenizer(temp, " \t\n\r\f,");
			if (!hash.containsKey(st.nextToken())) {
				writer.println(temp);
			}
		}
		reader.close();
		writer.close();
	}

	public static void main(String[] args) throws IOException {
		if (args.length!=2) {
			System.out.println("Expecting 2 arguments: struct filename and fams_to_delete filename.");
		} else {
			try {
				new removeFams(args[0], args[1]);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
}
