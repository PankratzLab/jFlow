package org.genvisis.link;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;

public class addMarkers {

	public addMarkers(String filename) throws IOException {
		BufferedReader reader = null;
		Hashtable<String,String> markers = new Hashtable<String,String>(), hash;
		String[] line;
		String temp, newMarker;
		int chr, numberOfNewMarkers;

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
			markers.put(line[0].toUpperCase(), chr+"\t"+line[2]);
		}
		reader.close();

		if (!new File(filename).exists()) {
			System.err.println("Error - could not find "+filename+" in current directory");
			System.exit(2);
		}
		reader = new BufferedReader(new FileReader(filename));
		reader.readLine();
		numberOfNewMarkers = reader.readLine().split("[\\s]+").length-3;
		for (int i = 0; i<numberOfNewMarkers; i++) {
			reader = new BufferedReader(new FileReader(filename));
			reader.readLine();
			newMarker = reader.readLine().split("[\\s]+")[i+3];
			hash = new Hashtable<String,String>();
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				hash.put(line[0], line[1]+"\t"+line[2]+"\t"+line[i*2+3]+"\t"+line[i*2+4]);
			}
			reader.close();

			markers.put("currentMarker", i+"");
			new addMarkers(markers, newMarker.toUpperCase(), hash);
		}
	}

	public addMarkers(Hashtable<String,String> markers, String newMarker, Hashtable<String,String> hash) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer;
		double trav, prev, mark;
		String[] line, data;
		int offset;
		boolean replacement = false;
		String bakFilename;
		int chr;

		if (!markers.containsKey(newMarker.toUpperCase())) {
			System.err.println("Error - "+newMarker.toUpperCase()+" was not found in the Marshfield database");
			System.err.println("      - make sure you have your placeholder line holding its place");
			System.exit(3);
		}
		line = (markers.get(newMarker.toUpperCase())).split("[\\s]+");
		chr = Integer.valueOf(line[0]).intValue();
		mark = Double.valueOf(line[1]).doubleValue();
		bakFilename = Files.getBakFilename("chromosome"+chr+".dat", super.getClass().getName());
		if (new File("chromosome"+chr+".dat").exists()) {
			new File("chromosome"+chr+".dat").renameTo(new File(bakFilename));
		} else {
			System.err.println("Error - "+"chromosome"+chr+".dat"+" was not found in current directory");
			System.exit(5);
		}

		reader = new BufferedReader(new FileReader(bakFilename));
		writer = new PrintWriter(new FileWriter("chromosome"+chr+".dat"));

		writer.println(reader.readLine());
		line = reader.readLine().split("[\t]+");
		writer.print(line[0]+"\t"+line[1]+"\t"+line[2]);
		offset = -1;
		trav = prev = 0.0;
		for (int i = 3; i<line.length; i++) {
			if (offset==-1) {
				data = (markers.get(line[i].toUpperCase())).split("[\\s]+");
				if (Integer.valueOf(data[0]).intValue()!=chr) {
					System.err.println("Error - Marker "+line[i]+" in chromosome"+chr+".dat is not on chromosome "+chr);
					quitProperly(reader, writer, bakFilename, "chromosome"+chr+".dat");
					System.exit(2);
				}
				trav = Double.valueOf(data[1]).doubleValue();
				if (prev>trav) {
					System.err.println("Error - the markers in chromosome"+chr+".dat are not in order");
					quitProperly(reader, writer, bakFilename, "chromosome"+chr+".dat");
					System.exit(3);
				}
				if (mark==trav&&line[i].toUpperCase().equals(newMarker.toUpperCase())) {
					offset = i-3;
					replacement = true;
				} else if (mark<=trav&&mark>=prev) {
					writer.print("\t"+newMarker+"\t");
					offset = i-3;
				}
				prev = trav;
			}
			writer.print("\t"+line[i]+"\t");
		}
		writer.println();
		while (reader.ready()) {
			line = reader.readLine().split("[\\s]+");
			writer.print(line[0]+"\t"+line[1]+"\t"+line[2]);
			for (int i = 3; i<line.length; i++) {
				if (i==3+offset*2) {
					if (hash.containsKey(line[0])) {
						data = (hash.get(line[0])).split("[\\s]+");
						if (!(line[1]+"\t"+line[2]).equals(data[0]+"\t"+data[1])) {
							System.err.println("Error - DNA mismatch with regards to "+line[0]+" ("+line[1]+"\t"+line[2]+" in chromosome"+chr+".dat and "+data[0]+"\t"+data[1]+" in file with new marker(s)");
							System.err.println("      - but letting it slide... THIS time, matching up using DNA numbers");
							// quitProperly(reader, writer, bakFilename,
							// "chromosome"+chr+".dat");
							// System.exit(4);
						}
						writer.print("\t"+data[2]+"\t"+data[3]);
					} else {
						if ((markers.get("currentMarker")).equals("0")) {
							System.err.println("Warning - "+line[1]+"\t"+line[2]+" does not have data for the new marker.");
						}
						writer.print("\t0\t0");
					}
					if (!replacement) {
						writer.print("\t"+line[i]);
						i++;
						writer.print("\t"+line[i]);
					} else {
						i++;
					}
				} else {
					writer.print("\t"+line[i]);
				}
			}
			writer.println();
		}
		reader.close();
		writer.close();
		if (!(markers.get("currentMarker")).equals("0")) {
			new File(bakFilename).delete();
		}
	}

	public void quitProperly(BufferedReader reader, PrintWriter writer, String bakFilename, String filename) throws IOException {
		reader.close();
		writer.close();
		new File(filename).delete();
		new File(bakFilename).renameTo(new File(filename));
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		// String filename = "corrected_d6s292.dat";
		String filename = "add_D6S305.dat";
		// String filename = "finemarkers.dat";
		// String filename = "chr5fine.dat";

		String usage = "\n"+"park.addMarker requires 1 argument:\n"+"   (1) a chromosom-like file - DNA FAMID INDID genotypes...\n"+"       (i.e. file="+filename+" (default))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		if (args.length==0) {
			System.err.println("Warning: using defaults (file="+filename+")");
		}
		try {
			new addMarkers(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}

	}
}
