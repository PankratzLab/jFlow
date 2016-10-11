package org.genvisis.assoc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

public class haplostats {
	public class haploguy {
		public String id;

		public String sex;

		public int[] data;

		public haploguy(String new_id, String new_sex, int[] new_data) {
			id = new_id;
			sex = new_sex;
			data = new_data;
		}
	}

	public haplostats(String filename) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null, pheno = null;
		String[] line;
		Hashtable<String, Vector<haploguy>> hash = new Hashtable<String, Vector<haploguy>>();
		Vector<String> fams = new Vector<String>();
		Vector<haploguy> v;
		String temp;
		int[] data;
		int score;
		haploguy dude;
		int count;

		if (!new File(filename).exists()) {
			System.err.println("Error - could not find " + filename + " in current directory");
			System.exit(2);
		}
		reader = new BufferedReader(new FileReader(filename));
		reader.readLine();
		while (reader.ready()) {
			line = reader.readLine().split("[\\s]+");
			data = new int[line.length - 6];
			if (data.length % 2 != 0) {
				System.err.println("Error - Expecting "	+ (line.length - 6)
														+ " columns (including affection status) followed by an even number of alleles");
				System.exit(4);
			}
			score = 0;
			for (int i = 0; i < line.length - 6; i++) {
				data[i] = Integer.valueOf(line[i + 6]).intValue();
				if (data[i] != 0) {
					score++;
				}
			}
			// if (score > 1 && line[5].equals("2")) {
			if (score == 8 && line[5].equals("2")) {
				if (hash.containsKey(line[0])) {
					v = hash.get(line[0]);
				} else {
					v = new Vector<haploguy>();
					fams.add(line[0]);
				}
				dude = new haploguy(""	+ (Integer.valueOf(line[0]).intValue() * 1000
																	+ Integer.valueOf(line[1]).intValue()),
														line[4], data);
				v.add(dude);
				hash.put(line[0], v);
			}
		}
		reader.close();

		writer = new PrintWriter(new FileWriter("haplostats.pre.prn"));
		pheno = new PrintWriter(new FileWriter("haplostats.dat"));

		for (int i = 0; i < fams.size(); i++) {
			v = hash.get(fams.elementAt(i));
			dude = v.elementAt(0);
			writer.print(fams.elementAt(i) + "\t" + dude.id + "\t0\t0\t" + dude.sex);
			pheno.println(dude.id + "\t1");
			data = dude.data;
			for (int element : data) {
				writer.print("\t" + element);
			}
			writer.println();
		}

		if (!new File("negneuros.ped").exists()) {
			System.err.println("Error - could not find " + "negneuros.ped" + " in current directory");
			System.exit(2);
		}
		reader = new BufferedReader(new FileReader("negneuros.ped"));
		count = 0;
		while (reader.ready()) {
			count++;
			line = reader.readLine().split("[\\s]+");
			temp = (99000 + count) + "";
			writer.print(temp + "\t" + temp + "\t0\t0\t1");
			pheno.println(temp + "\t0");
			for (int j = 6; j < line.length; j++) {
				writer.print("\t" + line[j]);
			}
			writer.println();
		}
		reader.close();
		writer.close();
		pheno.close();
	}

	public String translate(int[] data, String[][] lookup) throws IOException {
		String str = "";
		for (int i = 0; i < data.length / 2; i++) {
			str += (i == 0 ? "" : ",") + lookup[i][data[i * 2]] + lookup[i][data[i * 2 + 1]];
		}
		return str;
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		// String filename = "progeni-haploview.ped";
		String filename = "markers390-395.pre";

		String usage = "\n"	+ "park.haplostats requires 1 arguments:\n"
										+ "   (1) a pre file (i.e. file=" + filename + " (default))\n" + "";

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
		if (args.length == 0) {
			System.err.println("Using defaults (file=" + filename + ")");
		}

		try {
			new haplostats(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}

	}
}
