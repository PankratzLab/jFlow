package one;

import java.io.*;
//import java.util.*;

import common.*;

public class GEDI_Analyses {

	public static void dummyMarkerPositions() {
		BufferedReader reader;
		PrintWriter writer;
		int count;
		
		String filename = "D:/GEDI/SDRG/markerList.txt";
		try {
			writer = new PrintWriter(new FileWriter("D:/GEDI/SDRG/markerPositions.txt"));
			writer.println("Name\tChr\tMapInfo");
			try {
				count = 0;
				reader = new BufferedReader(new FileReader(filename));
				while (reader.ready()) {
					count++;
					writer.println(reader.readLine()+"\t"+(count/105000+1)+"\t"+count);
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + filename
						+ "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + filename + "\"");
				System.exit(2);
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + "D:/GEDI/SDRG/00src/markerPositions.txt");
			e.printStackTrace();
		}
	}
	
	private static void platelist(String pattern) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp;
		String filename;
		int[] indices;
		
		String[] columns_needed = new String[] {"Institute Sample Label", "Row", "Institute Plate Label", "Well"};
		
		try {
			writer = new PrintWriter(new FileWriter(pattern));
			writer.println("DNA\tPlate#\tRow\tPlateLabel\tWell");
			for (int i = 0; i < 120; i++) {
				filename = ext.insertNumbers(pattern, i);
				if (Files.exists(filename, false)) {
					try {
						reader = new BufferedReader(new FileReader(filename));
						do {
							temp = reader.readLine();
						} while (!temp.contains("Well"));
						indices = ext.indexFactors(columns_needed, temp.trim().split(","), false, true);
						while (reader.ready()) {
							line = reader.readLine().trim().split(",");
							writer.println(line[indices[0]]+"\t"+i+"\t"+line[indices[1]]+"\t"+line[indices[2]]+"\t"+line[indices[3]]);							
						}
						reader.close();
					} catch (FileNotFoundException fnfe) {
						System.err.println("Error: file \"" + filename + "\" not found in current directory");
						System.exit(1);
					} catch (IOException ioe) {
						System.err.println("Error reading file \"" + filename + "\"");
						System.exit(2);
					}
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + pattern);
			e.printStackTrace();
		}
		
		
		
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "GEDI_Analyses.dat";

		String usage = "\n" + "one.GEDI_Analyses requires 0-1 arguments\n"
				+ "   (1) filename (i.e. file=" + filename + " (default))\n"
				+ "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help")
					|| args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
//			dummyMarkerPositions();
			platelist("D:/data/GEDI_exome/00src/McGue11001_Exome_LIMS_manifest/McGue11001_Exome_p##.csv");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
