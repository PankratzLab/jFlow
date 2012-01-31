package kaput;

import java.io.*;
import java.util.*;

import common.*;

public class borel {

	public borel(int chromosome) throws IOException {
		BufferedReader ped = null;
		BufferedReader map = null;
		PrintWriter writer;
		StringTokenizer st;
		String temp, chrome;
		int numMarkers;
		double num;

		chrome = (chromosome<10)?"0"+chromosome:""+chromosome;

		try {
			ped = new BufferedReader(new FileReader("chrom"+chrome+".ped"));
		} catch (Exception e) {
			System.err.println("Could not open the pedigree file: chrom"+chrome+".ped");
			System.err.println("Please rectify");
			System.exit(1);
		}
		try {
			// map = new BufferedReader(new
			// FileReader("/home/npankrat/park/00masters/map"+chrome+".dat"));
			map = new BufferedReader(new FileReader("map"+chrome+".dat"));
		} catch (Exception e) {
			System.err.println("Could not open the map file for chromosome "+chromosome+": map"+chrome+".dat");
			System.err.println("Please rectify");
			System.exit(1);
		}

		writer = new PrintWriter(new FileWriter("r_chrom"+chrome+".ped"));

		while (ped.ready()) {
			temp = ped.readLine();
			writer.println(temp.substring(0, 22)+temp.substring(24));
		}
		ped.close();
		writer.close();

		writer = new PrintWriter(new FileWriter("r_map"+chrome+".dat"));

		temp = map.readLine();
		st = new StringTokenizer(temp);
		numMarkers = Integer.valueOf(st.nextToken()).intValue();
		writer.println((numMarkers-1)+temp.substring((numMarkers<10)?3:2));
		writer.println(map.readLine());
		for (int j = 1; j<numMarkers; j++) {
			writer.print(j+" ");
		}
		writer.println();
		map.readLine();
		map.readLine();
		map.readLine();
		map.readLine();
		map.readLine();

		temp = map.readLine();
		while (temp.startsWith("3 ")) {
			writer.println(temp);
			writer.println(map.readLine());
			temp = map.readLine();
		}
		writer.println(temp);
		st = new StringTokenizer(map.readLine());
		for (int j = 0; j<numMarkers-1; j++) {
			num = Double.valueOf(st.nextToken()).doubleValue();
			writer.print(ext.formDeci(num/100, 5)+" ");
		}
		writer.println(" << RECOMB VALUES");
		writer.println(map.readLine());

		// System.out.println(numMarkers);
		map.close();
		writer.close();
	}

	public static void main(String[] args) throws IOException {
		if (args.length!=1) {
			System.out.println("Expecting 1 argument: chromosome number.");
		} else {
			try {
				new borel(Integer.valueOf(args[0]).intValue());
			} catch (Exception e) {
				e.printStackTrace();
				System.err.println("Error in processing chromosome "+args[0]);
			}
		}
	}
}
