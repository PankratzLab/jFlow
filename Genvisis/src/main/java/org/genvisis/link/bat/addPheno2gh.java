package org.genvisis.link.bat;

import java.io.*;
import java.util.*;

public class addPheno2gh {

	public static String MISSING_VALUE = "-";

	public addPheno2gh(int chromosome, int columnOfStruct) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer;
		StringTokenizer st;
		String temp, chrome, id;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		double aao;

		try {
			reader = new BufferedReader(new FileReader("struct.dat"));
			while (reader.ready()) {
				st = new StringTokenizer(reader.readLine());
				id = st.nextToken()+":"+st.nextToken();
				for (int i = 0; i<columnOfStruct-3; i++) {
					st.nextToken();
				}
				temp = st.nextToken();
				if (!temp.equals(".")) {
					aao = Double.parseDouble(temp);
					// aao = Integer.valueOf(temp).intValue();
					if (aao>0&&aao<99) {
						hash.put(id, temp);
					}
				}
			}
			reader.close();
		} catch (Exception e) {
			System.err.println("Error trying to read column "+columnOfStruct+" of struct.dat");
			System.exit(1);
		}

		chrome = (chromosome<10)?"0"+chromosome:""+chromosome;
		try {
			reader = new BufferedReader(new FileReader("re_chrom"+chrome+".pre"));
		} catch (Exception e) {
			System.err.println("Could not open the pedigree file: re_chrom"+chrome+".pre");
			System.err.println("Please rectify");
			System.exit(1);
		}
		writer = new PrintWriter(new FileWriter("gh_chrom"+chrome+".pre"));

		while (reader.ready()) {
			temp = reader.readLine();
			st = new StringTokenizer(temp);
			id = st.nextToken()+":"+st.nextToken();
			if (hash.containsKey(id)) {
				writer.println(temp+"\t"+hash.get(id));
			} else {
				writer.println(temp+"\t"+MISSING_VALUE);
			}
		}
		reader.close();
		writer.close();

		try {
			reader = new BufferedReader(new FileReader("map"+chrome+".dat"));
		} catch (Exception e) {
			System.err.println("Could not open the map file: map"+chrome+".dat");
			System.err.println("Please rectify");
			System.exit(1);
		}
		writer = new PrintWriter(new FileWriter("gh_map"+chrome+".dat"));

		temp = reader.readLine();
		st = new StringTokenizer(temp);
		int numMarkers = Integer.valueOf(st.nextToken()).intValue();
		writer.println((numMarkers+1)+temp.substring((numMarkers<10)?1:2));
		writer.println(reader.readLine());
		for (int j = 1; j<=numMarkers+1; j++) {
			writer.print(j+" ");
		}
		writer.println();
		reader.readLine();

		writer.println(reader.readLine());
		writer.println(reader.readLine());
		writer.println(reader.readLine());
		writer.println(reader.readLine());
		if (chromosome==23) {
			reader.readLine();
		}

		temp = reader.readLine();
		while (!temp.startsWith("0 0 ")) {
			writer.println(temp);
			temp = reader.readLine();
		}
		writer.println("0 2");
		writer.println("");
		writer.println("");
		writer.println("");
		writer.println("");
		writer.println("");
		writer.println(temp);
		st = new StringTokenizer(reader.readLine());
		for (int i = 0; i<numMarkers-1; i++) {
			writer.print(st.nextToken()+" ");
		}
		writer.print("10.0 ");
		while (st.hasMoreTokens()) {
			writer.print(st.nextToken()+" ");
		}
		writer.println();
		writer.println(reader.readLine());

		reader.close();
		writer.close();
	}

	public static void main(String[] args) throws IOException {
		// new addPheno2gh(5,8);
		// BufferedReader stdin = new BufferedReader(new
		// InputStreamReader(System.in));
		if (args.length!=2) {
			System.out.println("Expecting 2 arguments: chromosome number and column of struct to add.");
		} else {
			try {
				new addPheno2gh(Integer.valueOf(args[0]).intValue(), Integer.valueOf(args[1]).intValue());
			} catch (Exception e) {
				e.printStackTrace();
				System.err.println("Error in processing chromosome "+args[0]);
			}
		}
	}
}
