package one;

import gwas.Metal;

import java.io.*;
import java.util.*;

import javax.swing.plaf.metal.MetalScrollBarUI;

import bioinformatics.Sequence;

import common.*;

public class IndianDiabetes {
	private static void parseGenotypes(String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp, trav;
		Hashtable<String, String> hash = new Hashtable<String, String>();
		Vector<String> v = new Vector<String>();
		int count;
		long time;
		Logger log;
		int numMarkers;
		String[] markerNames;
		
		log = new Logger(filename+".log");
		
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+".ped"));
			line = reader.readLine().trim().split(",", -1);
			ext.checkHeader(line, new String[] {"NUMBER ", "Sex", "Affected"}, new int[] {0,1,2}, false, log, true);
			numMarkers = line.length - 3;
			markerNames = Array.subArray(line, 3);
			while (reader.ready()) {
				line = reader.readLine().trim().split(",", -1);
				writer.print(line[0]+"\t"+line[0]+"\t0\t0"+"\t"+line[1]+"\t"+line[2]);
				if (line.length != numMarkers + 3) {
					log.reportError("Error - mismatched number of alleles for "+line[0]);
					System.exit(1);
				}
				for (int i = 0; i < numMarkers; i++) {
					if (line[i+3].equals("")) {
						writer.print("\t0\t0");
					} else if (line[i+3].length() != 2) {
						log.reportError("Error - mishapen alleles for indiviudal "+line[0]+" marker "+markerNames[i]+": '"+line[i+3]+"'");
					} else {
						for (int j = 0; j < 2; j++) {
							if (Array.indexOfChar(Sequence.ALLELES, line[i+3].charAt(j)) == -1) {
								log.reportError("Error - invalid allele for indiviudal "+line[0]+" marker "+markerNames[i]+": '"+line[i+3]+"'");
							}
							writer.print("\t"+line[i+3].charAt(j));
						}
					}
				}
				writer.println();
			}
			reader.close();
			writer.close();
			
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+".map"));
			for (int i = 0; i < markerNames.length; i++) {
				if (markerNames[i].toLowerCase().startsWith("chr")) {
					line = markerNames[i].substring(3).split("_");
					writer.println(line[0]+"\t"+markerNames[i]+"\t0\t"+line[1]);
				} else {
					System.err.println("Error - don't know how to parse the position for marker: "+markerNames[i]);
					writer.println("0\t"+markerNames[i]+"\t0\t0");
				}
			}
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
		String dir = "D:/Myron/Indian_Diabetes/SequencingPilot/Replication_Sequenome_Indian/";
		String genotypes = "IndianReplication_genotypes.csv";
		
		parseGenotypes(dir+genotypes);
	}
}
