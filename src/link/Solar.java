// to do - missing first individual
// to do - don't put 99+ in ptypes file
// to do - any reason why h2 is so far off?

package link;

import java.io.*;
import java.util.*;

import common.*;

import park.tools;

public class Solar {
	public static void createFiles(int chr, String trait) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null, famtastic = null;
		String temp, famID, indID;
		StringTokenizer st;
		String[] line;
		Vector<String> markerV = new Vector<String>();
		int numMarkers, numAlleles;
		double total;
		Hashtable<String,String> hash;

		reader = new BufferedReader(new FileReader("map"+ext.chrome(chr)+".dat"));
		writer = new PrintWriter(new FileWriter("solar.freqs."+chr));

		st = new StringTokenizer(reader.readLine());
		numMarkers = Integer.valueOf(st.nextToken()).intValue()-1;
		for (int i = 0; i<6; i++) {
			reader.readLine();
		}
		if (chr==23) {
			reader.readLine();
		}
		for (int i = 0; i<numMarkers; i++) {
			st = new StringTokenizer(reader.readLine());
			st.nextToken();
			numAlleles = Integer.valueOf(st.nextToken()).intValue();
			st.nextToken();
			temp = st.nextToken();
			markerV.add(temp);
			writer.print(ext.formStr(temp, 9, true));
			st = new StringTokenizer(reader.readLine());
			for (int j = 1; j<=numAlleles; j++) {
				writer.print(" "+j+" "+st.nextToken());
			}
			writer.println();
		}
		writer.close();

		reader.readLine();
		st = new StringTokenizer(reader.readLine());
		st.nextToken();

		writer = new PrintWriter(new FileWriter("solar.map."+chr));
		writer.println(chr);
		total = 0;
		writer.println(ext.formStr(markerV.elementAt(0), 9, true)+"   0.00");
		for (int i = 1; i<numMarkers; i++) {
			total += Double.valueOf(st.nextToken()).doubleValue();
			writer.println(ext.formStr(markerV.elementAt(i), 9, true)+" "+ext.formStr(ext.formDeci(total, 2, true), 6));
		}
		reader.close();
		writer.close();

		reader = new BufferedReader(new FileReader("re_chrom"+ext.chrome(chr)+".pre"));
		writer = new PrintWriter(new FileWriter("solar.gtypes."+chr));
		famtastic = new PrintWriter(new FileWriter("solar.fam"));
		writer.print("FAMID,ID");
		famtastic.println("FAMID,ID,FA,MO,SEX");
		for (int i = 0; i<markerV.size(); i++) {
			writer.print(","+markerV.elementAt(i));
		}
		writer.println();
		while (reader.ready()) {
			st = new StringTokenizer(reader.readLine());
			famID = st.nextToken();
			indID = st.nextToken();
			writer.print(famID+",");
			writer.print(indID);
			famtastic.println(famID+","+indID+","+st.nextToken()+","+st.nextToken()+","+st.nextToken());
			st.nextToken();
			do {
				writer.print(","+st.nextToken()+"/"+st.nextToken());
			} while (st.hasMoreElements());
			writer.println();
		}
		reader.close();
		writer.close();
		famtastic.close();

		hash = tools.pullTraitFromDB(trait);

		reader = new BufferedReader(new FileReader("re_chrom"+ext.chrome(chr)+".pre"));
		writer = new PrintWriter(new FileWriter("solar.ptypes"));
		writer.println("FAMID,ID,"+trait);

		while (reader.ready()) {
			line = reader.readLine().split("[\\s]+");
			temp = hash.containsKey(line[0]+"\t"+line[1])?hash.get(line[0]+"\t"+line[1]):"";
			writer.println(line[0]+","+line[1]+","+(temp.equals(".")?"":temp));
		}
		reader.close();
		writer.close();
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		int chr = -1;
		String trait = "AOO";

		String usage = "\n"+"park.createSolar requires 1-2 arguments\n"+"   (1) chromosome number (i.e. chr=2)\n"+"   (2) trait (i.e. trait="+trait+" (default)\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("chr=")) {
				chr = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("trait=")) {
				trait = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		if (chr==-1) {
			System.err.println("Error - createSolar requires you to specify which chromosome to do");
			System.exit(2);
		}
		try {
			createFiles(chr, trait);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
