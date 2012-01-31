package widgets;

import java.io.*;
import java.util.*;

public class Load {
	public static void createHeader(String root, String output, boolean split) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, batch_dna;
		Vector<String> v = new Vector<String>();
		int count;
		String[] markerNames;
		
		try {
			reader = new BufferedReader(new FileReader(root+".map"));
			while (reader.ready()) {
				v.add(reader.readLine().trim().split("[\\s]+")[1]);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + root+".map" + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + root+".map" + "\"");
			try {
				writer = new PrintWriter(new FileWriter("ERROR_COULD_NOT_OPEN_"+root+".map_FILE.TXT"));
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to " + "ERROR_COULD_NOT_OPEN_"+root+".map_FILE.TXT");
				e.printStackTrace();
			}
			return;
		}
		
		markerNames = new String[v.size()];
		for (int i = 0; i < markerNames.length; i++) {
			markerNames[i] = v.elementAt(i);
		}
		
		try {
			reader = new BufferedReader(new FileReader(root+".ped"));
			writer = new PrintWriter(new FileWriter(output));
			writer.print((split?"BATCH\t":"")+"FID\tIID\tFATHER\tMOTHER\tSEX\tAFFECTED");
			for (int i = 0; i < markerNames.length; i++) {
				writer.print("\t"+markerNames[i]+"A"+"\t"+markerNames[i]+"B");
			}
			writer.println();
			count = 0;
			while (reader.ready()) {
				count++;
				line = reader.readLine().trim().split("[\\s]+");
				if (line.length != markerNames.length*2+6) {
					System.err.println("Error - line "+count+" has "+line.length+" columns instead of the expected "+markerNames.length*2+6+" (#markers*2+6)");
					try {
						writer = new PrintWriter(new FileWriter("ERROR_NUMBER_OF_COLUMNS_DON'T_MATCHUP_FOR_THESE_LINES.TXT", true));
						writer.println(count);
						writer.close();
					} catch (Exception e) {
						System.err.println("Error writing to " + "ERROR_COULD_NOT_OPEN_"+root+".map_FILE.TXT");
						e.printStackTrace();
					}
				}
				for (int i = 0; i < line.length; i++) {
					if (i==0) {
						if (split) {
							batch_dna = line[i].split("_");
							if (batch_dna.length != 2) {
								try {
									writer = new PrintWriter(new FileWriter("ERROR_THESE_FIDS_AREN'T_CONCATENATED_BY_AN_UNDERSCORE.TXT", true));
									writer.println("line"+count+"\t"+line[i]);
									writer.close();
								} catch (Exception e) {
									System.err.println("Error writing to " + "ERROR_COULD_NOT_OPEN_"+root+".map_FILE.TXT");
									e.printStackTrace();
								}
								writer.print(line[i]+"\t"+line[i]);
							} else {
								writer.print(batch_dna[0]+"\t"+batch_dna[1]);
							}
						} else {
							writer.print(line[i]);
						}
					} else {
						writer.print("\t"+line[i]);
					}
				}
				writer.println();
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + root+".map" + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + root+".map" + "\"");
			try {
				writer = new PrintWriter(new FileWriter("ERROR_COULD_NOT_OPEN_"+root+".map_FILE.TXT"));
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to " + "ERROR_COULD_NOT_OPEN_"+root+".map_FILE.TXT");
				e.printStackTrace();
			}
			return;
		}
		
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String root = "plink";
		String output = "plink.txt";
		boolean split = false;

		String usage = "\n" + 
		"widgets.Load requires 0-1 arguments\n" + 
		"   (1) root of plink file (i.e. root=" + root + " (default))\n" + 
		"   (2) name of output file (i.e. output=" + output + " (default))\n" + 
		"   (3) split the FID (i.e. -split (not the default))\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("root=")) {
				root = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("output=")) {
				output = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-split")) {
				split = true;
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			createHeader(root, output, split);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
