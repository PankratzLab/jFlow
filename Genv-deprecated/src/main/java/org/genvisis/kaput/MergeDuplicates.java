// -Xms1024M -Xmx1024M
package org.genvisis.kaput;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;

public class MergeDuplicates {
	public static final String WINDOWS_DIRECTORY = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\duplicates\\";
	public static final String LINUX_DIRECTORY = "";
	public static final String PEDS = "peds/";
	public static final String MERGED = "merged/";

	public static void mergeDuplicates(String filename, String markerfile) {
		BufferedReader reader;
		PrintWriter writer, log;
		String[] line;
		String temp;
		Vector<String> v = new Vector<String>(), markers;
		int count = 0;
		String dir = new File(WINDOWS_DIRECTORY).exists()?WINDOWS_DIRECTORY:LINUX_DIRECTORY;
		String[][] samples;
		int numSamples, numMarkers, source;

		markers = HashVec.loadFileToVec(dir+markerfile, false, false, false);
		numMarkers = markers.size();

		try {
			reader = new BufferedReader(new FileReader(dir+filename));
			while (reader.ready()) {
				temp = reader.readLine();
				v.add(temp);
				line = temp.split("[\\s]+");
				if (line.length>count) {
					count = line.length;
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+dir+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+dir+filename+"\"");
			System.exit(2);
		}

		if (!new File(dir+MERGED).exists()) {
			new File(dir+MERGED).mkdirs();
		}

		log = Files.getWriter(dir+"merge.log");
		log.println("Sample\tAlternate(s)\tMarker\tFinal Allele1\tFinal Allele2\tCall Source\tmethod"+Array.toStr(Array.stringArraySequence(count, "\tSample #")));
		for (int i = 0; i<v.size(); i++) {
			System.out.println("Processing: "+v.elementAt(i));
			line = v.elementAt(i).split("[\\s]+");
			numSamples = line.length;
			samples = new String[numSamples+1][];
			for (int j = 0; j<numSamples; j++) {
				if (!new File(dir+PEDS+line[j]+".ped").exists()) {
					System.err.println("Error - file '"+line[j]+".ped' not found in directory '"+dir+PEDS+"'");
					System.exit(1);
				}
				try {
					reader = new BufferedReader(new FileReader(dir+PEDS+line[j]+".ped"));
					samples[j] = Array.subArray(reader.readLine().split("[\\s]+"), 1);
					if (j==0) {
						count = samples.length;
					} else if (samples.length!=count) {
						System.err.println("Error - mismatched number of columns for DNAs: "+v.elementAt(i));
						System.exit(1);
					}
					if (samples[j].length!=numMarkers*2) {
						System.err.println("Error - number of columns ("+(samples[j].length+1)+") does not equal 2*#markers for DNA ("+numMarkers*2+"): "+line[j]);
						System.exit(1);
					}
					reader.close();
				} catch (IOException ioe) {
					System.err.println("Error reading file \""+dir+PEDS+line[j]+".ped"+"\"");
					System.exit(2);
				}
			}
			samples[numSamples] = new String[numMarkers*2];
			for (int j = 0; j<numMarkers; j++) {
				temp = "";
				source = -1;
				for (int k = 0; k<numSamples; k++) {
					if (k==0) {
						samples[numSamples][j*2+0] = samples[k][j*2+0];
						samples[numSamples][j*2+1] = samples[k][j*2+1];
						source = k;
					} else if ((samples[numSamples][j*2+0].equals("0")||samples[numSamples][j*2+1].equals("0"))&&!samples[k][j*2+0].equals("0")&&!samples[k][j*2+1].equals("0")) {
						samples[numSamples][j*2+0] = samples[k][j*2+0];
						samples[numSamples][j*2+1] = samples[k][j*2+1];
						temp = line[0]+"\t"+Array.toStr(Array.subArray(line, 1), "|")+"\t"+markers.elementAt(j)+"\t"+samples[numSamples][j*2+0]+"\t"+samples[numSamples][j*2+1]+"\t"+line[k]+"\tmerged";
						for (int m = 0; m<numSamples; m++) {
							temp += "\t"+samples[m][j*2+0]+"\t"+samples[m][j*2+1];
						}
						source = k;
					} else if ((samples[k][j*2+0].equals("0")||samples[k][j*2+1].equals("0"))&&!samples[numSamples][j*2+0].equals("0")&&!samples[numSamples][j*2+1].equals("0")) {
						temp = line[0]+"\t"+Array.toStr(Array.subArray(line, 1), "|")+"\t"+markers.elementAt(j)+"\t"+samples[numSamples][j*2+0]+"\t"+samples[numSamples][j*2+1]+"\t"+line[source]+"\tmerged";
						for (int m = 0; m<numSamples; m++) {
							temp += "\t"+samples[m][j*2+0]+"\t"+samples[m][j*2+1];
						}
					} else if (!samples[numSamples][j*2+0].equals(samples[k][j*2+0])||!samples[numSamples][j*2+1].equals(samples[k][j*2+1])) {
						samples[numSamples][j*2+0] = "X";
						samples[numSamples][j*2+1] = "X";
					}
				}
				if (samples[numSamples][j*2+0].equals("X")) {
					samples[numSamples][j*2+0] = "0";
					samples[numSamples][j*2+1] = "0";
					temp = line[0]+"\t"+Array.toStr(Array.subArray(line, 1), "|")+"\t"+markers.elementAt(j)+"\t0\t0\tnull\tmismatch";
					for (int m = 0; m<numSamples; m++) {
						temp += "\t"+samples[m][j*2+0]+"\t"+samples[m][j*2+1];
					}
				}
				if (!temp.equals("")) {
					log.println(temp);
					log.flush();
				}
			}
			writer = Files.getWriter(dir+MERGED+line[0]+".ped");
			writer.print(line[0]);
			for (int j = 0; j<numMarkers; j++) {
				writer.print("\t"+samples[numSamples][j*2+0]+"\t"+samples[numSamples][j*2+1]);
			}
			writer.println();
			writer.close();

		}
		log.close();

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "MergeDuplicates.dat";
		String markers = "markers.dat";

		String usage = "\\n"+"park.gwa.MergeDuplicates requires 0-1 arguments\n"+"   (1) filename (i.e. file="+filename+" (default))\n"+"   (2) filename of markers (i.e. markers="+markers+" (default))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("markers=")) {
				markers = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			mergeDuplicates(filename, markers);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
