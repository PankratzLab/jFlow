package org.genvisis.link.bat;

import java.io.*;
import java.util.*;

public class solarGatherByFam {
	public solarGatherByFam(String[] args) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		StringTokenizer st;
		String temp, fam;
		String[] chrome;
		int[] chromosome, position;
		Vector<String> fams = new Vector<String>();
		Hashtable<String,String> hash = new Hashtable<String,String>();

		reader = new BufferedReader(new FileReader("solar_ped"));
		reader.readLine();
		while (reader.ready()) {
			st = new StringTokenizer(reader.readLine(), ",");
			temp = st.nextToken();
			if (!hash.containsKey(temp)) {
				hash.put(temp, ";");
				fams.add(temp);
			}
		}

		st = new StringTokenizer((new File(".")).getAbsolutePath(), "/");
		while (st.countTokens()>3) {
			st.nextElement();
		}
		writer = new PrintWriter(new FileWriter(st.nextToken()+"-"+st.nextToken()+" linkage summary.xls"));

		chromosome = new int[args.length];
		position = new int[args.length];
		chrome = new String[args.length];

		writer.print("chr");
		for (int i = 0; i<args.length; i++) {
			chromosome[i] = Integer.valueOf(args[i].substring(0, args[i].indexOf("@"))).intValue();
			position[i] = Integer.valueOf(args[i].substring(args[i].indexOf("@")+1)).intValue();
			chrome[i] = (chromosome[i]<10)?"0"+chromosome[i]:""+chromosome[i];
			writer.print("\t"+args[i]);
		}
		writer.println();

		for (int huh = 0; huh<fams.size(); huh++) {
			fam = fams.elementAt(huh);
			writer.print(fam);
			for (int i = 0; i<args.length; i++) {
				try {
					reader = new BufferedReader(new FileReader("fam"+fam+"/chrom"+chrome[i]+"/AAO/multipoint1.out"));

					for (int j = 0; j<position[i]+2; j++) {
						reader.readLine();
					}
					st = new StringTokenizer(reader.readLine());
					for (int j = 0; j<4; j++) {
						st.nextToken();
					}
					writer.print("\t"+st.nextToken());
					reader.close();
				} catch (Exception e) {}
			}
			writer.println();
		}
		writer.close();
	}

	public static void main(String[] args) throws IOException {
		if (args.length==0) {
			System.err.println("Correct usage : java gatherByFam 2@170 12@1 21@42 ... etc.");
		}
		try {
			new solarGatherByFam(args);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
