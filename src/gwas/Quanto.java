package gwas;

import java.io.*;
//import java.util.*;
import common.*;

public class Quanto {
	public static final String[] EXPECTED_TRIADS = {"Frequency", "RG", "Gene", "kP"};
	
	private static void parseResults(String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String freq;
		double prev, trav;
		
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+".out"));
			line = reader.readLine().trim().split("[\\s]+");
			ext.checkHeader(line, EXPECTED_TRIADS, true);
			freq = null;
			prev = 999;
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (line.length > 3) {
					freq = line[0];
					prev = 999;
					line = Array.subArray(line, 1);
				}
				trav = Double.parseDouble(line[1]);
				if (prev < 0.8 && trav >= 0.8) {
					writer.println(freq+"\t"+Array.toStr(line));
				}			
				
				prev = trav;
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "Quanto.dat";

		String usage = "\n" + "gwas.Quanto requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
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
			filename = "D:/grants/2012.06 Hepatoblastoma/QuantoPower/input.txt";
			filename = "D:/grants/2012.06 Hepatoblastoma/QuantoPower/470trios.txt";
			filename = "D:/grants/2012.06 Hepatoblastoma/QuantoPower/290_4cc.txt";
			filename = "D:/grants/2012.06 Hepatoblastoma/QuantoPower/520_4cc.txt";
			filename = "D:/grants/2012.06 Hepatoblastoma/QuantoPower/600_4cc.txt";
			filename = "D:/grants/2012.06 Hepatoblastoma/QuantoPower/100_4cc_to5.txt";
			filename = "D:/grants/2012.06 Hepatoblastoma/QuantoPower/520_4cc_5E6.txt";
			
			parseResults(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
