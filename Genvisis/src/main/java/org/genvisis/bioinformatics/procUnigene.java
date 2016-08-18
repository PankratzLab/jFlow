package org.genvisis.bioinformatics;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;

public class procUnigene {
	public static final String[] VARS = {"ID", "GENE ", "TITLE", "GENE_ID", "CHROMOSOME", "CYTOBAND", "EXPRESS"};

	public procUnigene(String filename) {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;
		String temp = "", trav;
		Vector<String> v = new Vector<String>();
		IntVector iv = new IntVector();
		int index, count = 0;
		String[] data = Array.stringArray(VARS.length+1, ".");
		int[] keys;

		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(filename+".xls"));
			writer.println("NM_ID\t"+Array.toStr(VARS));
			while (reader.ready()) {
				temp = reader.readLine();
				if (temp.equals("//")) {
					writer.println(Array.toStr(data));
					line = data[ext.indexOfStr("EXPRESS", VARS)+1].split("\\|");
					if (!line[0].equals(".")&&!data[3+1].equals(".")) {
						count++;
						for (int i = 0; i<line.length; i++) {
							trav = line[i].trim();
							if (v.contains(trav)) {
								index = v.indexOf(trav);
								iv.setElementAt(iv.elementAt(index)+1, index);
							} else {
								v.add(trav);
								iv.add(1);
							}
						}
					}
					data = Array.stringArray(VARS.length+1, ".");
				}
				for (int i = 0; i<VARS.length; i++) {
					if (temp.startsWith(VARS[i])) {
						data[i+1] = temp.substring(VARS[i].length()).trim();
					}
				}
				if (temp.indexOf("ACC=NM_")>0) {
					temp = temp.indexOf(";")>0?temp.substring(temp.indexOf("ACC=NM_")+4, temp.indexOf(";")):temp.substring(temp.indexOf("ACC=NM_")+4).split("[\\s]+")[0];
					if (data[0].equals(".")) {
						data[0] = temp;
					} else {
						data[0] += ";"+temp;
					}
				}

			}
			writer.close();
			reader.close();

			writer = new PrintWriter(new FileWriter("possibleTissues.out"));
			keys = Sort.quicksort(Array.toStringArray(v));
			for (int i = 0; i<v.size(); i++) {
				writer.println(v.elementAt(keys[i])+"\t"+iv.elementAt(keys[i])+"\t"+((double)iv.elementAt(keys[i])/(double)count));
			}
			writer.println();
			writer.println("Total\t"+count);
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String filename = "C:\\Download\\Hs.data";

		String usage = "\n"+"park.procUnigene requires 0-1 arguments\n"+"   (1) filename (i.e. file="+filename+" (default))\n"+"";

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
		try {
			new procUnigene(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
