package link;

import java.io.*;
import java.util.*;

import common.*;

public class shiftAlleles {

	public shiftAlleles(String keyfile, String plateList) throws IOException {
		BufferedReader reader = null;
		String[] line;
		String temp;
		int numMarkers, chr;
		Vector<String> oldAlleles, newAlleles;
		Hashtable<String,String> hash, markers, integrityCheck = new Hashtable<String,String>();
		Vector<Hashtable<String,String>> plates = new Vector<Hashtable<String,String>>();
		Vector<String> plateNames = new Vector<String>();
		boolean[] touched = new boolean[23];
		String markerName, plateAffected;

		if (!new File(plateList).exists()) {
			System.err.println("Error - could not find "+plateList+" in current directory");
			System.exit(2);
		}
		reader = new BufferedReader(new FileReader(plateList));
		while (reader.ready()) {
			line = reader.readLine().split("[\\s]+");
			if (line.length<3&&!line[0].equals("")) {
				System.err.println("Error - file with plate listing requires three columns (FamID IndID plate)");
				System.err.println("      - could not process line beginning with '"+line[0]+"'");
				System.exit(1);
			}
			if (plateNames.contains(line[2])) {
				hash = plates.elementAt(plateNames.indexOf(line[2]));
			} else {
				plateNames.add(line[2]);
				plates.add(hash = new Hashtable<String,String>());
			}
			hash.put(line[0]+"\t"+line[1], line[2]);
		}
		reader.close();

		try {
			reader = new BufferedReader(new FileReader("marker.database"));
		} catch (FileNotFoundException fnfe) {
			try {
				reader = new BufferedReader(new FileReader("/home/npankrat/marker.database"));
			} catch (FileNotFoundException fnfe2) {
				System.err.println("Error - could not find marker.database in current or root directory");
				System.exit(1);
			}
		}
		markers = new Hashtable<String,String>();
		chr = 0;
		while (reader.ready()) {
			temp = reader.readLine();
			if (temp.equals("")||temp.startsWith("Marker")) {
				continue;
			}
			if (temp.startsWith("chromosome")) {
				chr = Integer.valueOf(temp.substring(10, 12)).intValue();
				continue;
			}
			line = temp.split("[ ]+");
			markers.put(line[0].toUpperCase(), chr+"");
		}
		reader.close();

		if (!new File(keyfile).exists()) {
			System.err.println("Error - could not find "+keyfile+" in current directory");
			System.exit(2);
		}
		reader = new BufferedReader(new FileReader(keyfile));
		numMarkers = 0;
		while (reader.ready()) {
			reader.readLine();
			numMarkers++;
		}
		reader.close();
		if (numMarkers%2!=0) {
			System.err.println("Error - Odd number of lines in "+keyfile);
			System.err.println("      - shiftAlleles requires 2 lines for each marker, the old and new alleles");
			System.err.println("      - start the old alleles with the marker name and the new alleles with the plates affected");
		}
		numMarkers /= 2;

		reader = new BufferedReader(new FileReader(keyfile));
		for (int i = 0; i<numMarkers; i++) {
			line = reader.readLine().split("[\\s]+");
			markerName = line[0].toUpperCase();
			line[0] = "0";
			oldAlleles = new Vector<String>();
			for (int j = 0; j<line.length; j++) {
				oldAlleles.add(line[j]);
			}
			line = reader.readLine().split("[\\s]+");
			plateAffected = line[0];
			if (integrityCheck.containsKey(markerName+":"+plateAffected)) {
				// if you wanted to be really fancy you could find out which
				// alleles were changed and then report only those markers where
				// this warning is pertinent
				// System.err.println("Warning - Plate(s) "+plateAffected+" have
				// more than one key for marker "+markerName+". Possible repeat
				// of shifts occuring.");
			} else {
				integrityCheck.put(markerName+":"+plateAffected, "");
			}
			line[0] = "0";
			newAlleles = new Vector<String>();
			for (int j = 0; j<line.length; j++) {
				newAlleles.add(line[j]);
			}
			if (!markers.containsKey(markerName)) {
				System.err.println("Error - '"+markerName+"' (number "+(i+1)+" in the keyfile) is not listed as a marker in the Marshfield database");
				System.exit(4);
			}
			chr = Integer.valueOf(markers.get(markerName)).intValue();
			if (!new File("chromosome"+chr+".dat").exists()) {
				System.err.println("Error - could not find "+"chromosome"+chr+".dat"+" in current directory");
				System.err.println("      - necessary to shift alleles for marker "+markerName);
				System.exit(5);
			}
			if (!touched[chr-1]) {
				temp = Files.getBakFilename("chromosome"+chr+".dat", super.getClass().getName());
				new File("chromosome"+chr+".dat").renameTo(new File(temp));
				if (!Files.copyFile(temp, "chromosome"+chr+".dat")) {
					System.err.println("Error - could not backup "+"chromosome"+chr+".dat");
					System.exit(6);
				}
				touched[chr-1] = true;
			}
			if (!plateNames.contains(plateAffected)) {
				System.err.println("Error - '"+plateAffected+"' (from one of the second lines in the key file) was not the name of a plate seen in file "+plateList);
				System.exit(7);
			}
			shift(markerName, chr, oldAlleles, newAlleles, plates.elementAt(plateNames.indexOf(plateAffected)));

		}
		reader.close();
	}

	public void shift(String markerName, int chr, Vector<String> oldAlleles, Vector<String> newAlleles, Hashtable<String,String> impacted) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer;
		String[] line;
		String temp;
		int index = -1, numMarkers;

		new File("chromosome"+chr+".dat").renameTo(new File("temp"));
		reader = new BufferedReader(new FileReader("temp"));
		writer = new PrintWriter(new FileWriter("chromosome"+chr+".dat"));

		writer.println(reader.readLine());
		temp = reader.readLine();
		writer.println(temp);
		line = temp.split("[\\s]+");
		numMarkers = line.length-3;
		for (int i = 3; i<line.length; i++) {
			if (line[i].equals(markerName)) {
				index = i-3;
			}
		}
		if (index==-1) {
			System.err.println("Error - could not find '"+markerName+"' in "+"chromosome"+chr+".dat");
			System.err.println("        skipping marker and continuing with file");
			return;
		}

		while (reader.ready()) {
			line = reader.readLine().split("[\\s]+");
			if (line[0].equals("")) {
				continue;
			}
			if (line.length!=3+2*numMarkers) {
				System.err.println("Error - invalid number of columns for "+line[0]+"-"+line[1]);
				System.err.println("        expecting 3 + 2*numMarkers ("+(3+2*numMarkers)+" columns) and found "+line.length);
				System.exit(3);
			}
			for (int i = 0; i<3; i++) {
				writer.print((i==0?"":"\t")+line[i]);
			}
			for (int i = 0; i<numMarkers; i++) {
				if (i==index&&impacted.containsKey(line[1]+"\t"+line[2])) {
					for (int j = 0; j<2; j++) {
						if (oldAlleles.indexOf(line[3+i*2+j])==-1) {
							System.err.println("Error - Allele "+line[3+i*2+j]+" was not found in the key for marker "+markerName);
							System.err.println("      - and therefore was not shifted");
							writer.print("\t"+line[3+i*2+j]);
						} else {
							writer.print("\t"+newAlleles.elementAt(oldAlleles.indexOf(line[3+i*2+j])));
						}
					}
				} else {
					writer.print("\t"+line[3+i*2]);
					writer.print("\t"+line[3+i*2+1]);
				}
			}
			writer.println();
		}
		reader.close();
		writer.close();
		new File("temp").delete();

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String keyfile = "shiftAllelesForPlates.txt";
		String plateList = "plateList.dat";

		String usage = "\n"+"park.shiftAlleles requires 2 arguments:\n"+"   (1) list of markers, plates affected, and the key to the new alleles (i.e. key="+keyfile+" (default))\n"+"       start the old alleles with the marker name and the new alleles with the plates affected"+"   (2) list of individuals and their plates (i.e. pl="+plateList+" (default))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("key=")) {
				keyfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("pl=")) {
				plateList = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		System.out.println("Using "+keyfile+" to recode appropriate chromosomeX.dat files using plateList file "+plateList);
		try {
			new shiftAlleles(keyfile, plateList);
		} catch (Exception e) {
			e.printStackTrace();
		}

	}
}
