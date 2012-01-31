package park;

import java.io.*;
import java.util.*;

import common.*;

public class keepInds {
	public static boolean KEEPFIRSTLINE;

	public keepInds(String target, String keepsfile) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		StringTokenizer st;
		String temp, id;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		boolean isFirstLine; // workaround for Leah
		String bakFilename;

		reader = new BufferedReader(new FileReader(keepsfile));
		while (reader.ready()) {
			st = new StringTokenizer(reader.readLine());
			id = st.nextToken()+":"+st.nextToken();
			hash.put(id, "");
		}
		reader.close();

		isFirstLine = true;

		bakFilename = Files.getBakFilename(target, super.getClass().getName());
		(new File(target)).renameTo(new File(bakFilename));
		reader = new BufferedReader(new FileReader(bakFilename));
		writer = new PrintWriter(new FileWriter(target));
		while (reader.ready()) {
			temp = reader.readLine();
			st = new StringTokenizer(temp);
			id = st.nextToken()+":"+st.nextToken();
			if (hash.containsKey(id)) {
				writer.println(temp);
			} else if (isFirstLine&&KEEPFIRSTLINE) {
				writer.println(temp);
				isFirstLine = false;
			}
		}
		reader.close();
		writer.close();
	}

	public static void main(String[] args) throws IOException {
		if (args.length<2||args.length>3) {
			System.out.println("Expecting 2-3 arguments: struct filename and fams_to_keep filename (or \"key=key_you_want_to_keep\").");
			System.out.println("Third optional argument: -keepfirstline");
		} else {
			try {
				if (args.length==3&&args[2].equals("-keepfirstline")) {
					KEEPFIRSTLINE = true;
				} else {
					KEEPFIRSTLINE = false;
				}
				new keepInds(args[0], args[1]);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
}
