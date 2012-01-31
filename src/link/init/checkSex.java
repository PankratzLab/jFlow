// Make PD independent: change so that this used re_chrom23.pre instead of chromosome23.dat
package link.init;

import java.io.*;
import java.util.*;
import park.tools;

public class checkSex {
	public checkSex() throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String id, first, second, dna;
		StringTokenizer st;
		int ctSame, ctDiff;
		Hashtable<String,String> hash = new Hashtable<String,String>();

		reader = tools.getNinfoReader(1, false);
		reader.readLine();
		while (reader.ready()) {
			st = new StringTokenizer(reader.readLine());
			id = st.nextToken()+"-"+st.nextToken();
			st.nextToken();
			st.nextToken();
			st.nextToken();
			hash.put(id, st.nextToken());
		}
		reader.close();

		try {
			reader = new BufferedReader(new FileReader("chromosome23.dat"));
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: could not find the file \"chromosome23.dat\"");
			System.exit(1);
		}

		writer = new PrintWriter(new FileWriter("checkSex.out", true));
		System.out.println("Verifying sex of males and females using X chromosome marker data.");
		reader.readLine();
		reader.readLine();
		while (reader.ready()) {
			st = new StringTokenizer(reader.readLine());
			dna = st.nextToken();
			id = st.nextToken()+"-"+st.nextToken();
			ctSame = ctDiff = 0;
			while (st.hasMoreTokens()) {
				first = st.nextToken();
				second = st.nextToken();
				if (!first.equals("0")||!first.equals("0")) {
					if (first.equals(second)) {
						ctSame++;
					} else {
						ctDiff++;
					}
				}
			}

			writer.println(id+"\t"+hash.get(id)+"\t"+dna+"\t"+ctSame+"\t"+ctDiff+"\t"+((ctDiff>3)?"female":"male"));
			if (ctDiff+ctSame>3&&((ctDiff>3&&!(hash.get(id)).equals("F"))||(ctDiff<=3&&(hash.get(id)).equals("F")))) {
				System.err.println("Check: "+id+"\t"+hash.get(id)+"\t"+dna+"\t"+ctSame+"\t"+ctDiff+"\t"+((ctDiff>3)?"female":"male"));
			}
		}
		reader.close();
		writer.close();
	}

	public static void main(String[] args) throws IOException {
		if (args.length>0) {
			System.out.println("Expecting no arguments.");
		}
		try {
			new checkSex();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
