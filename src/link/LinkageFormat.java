package link;

import java.io.*;
import java.util.*;

import common.*;

public class LinkageFormat {
	public static String STRUCT_EXAMPLE = ""+"70001	101	0	0	1	0	NoDNAprsnt	.	.\n"+"70001	102	0	0	2	0	NoDNAprsnt	.	.\n"+"70001	1	101	102	1	2	1999PD0004	55	WHTE\n"+"70001	11	101	102	1	2	1999PD0005	59	WHTE"+"70002	101	0	0	1	0	NoDNAprsnt	.	."+"70002	102	0	0	2	0	NoDNAprsnt	.	."+"70002	1	101	102	1	2	1999PD0014	75	WHTE"+"70002	13	101	102	1	2	1999PD0015	60	WHTE";

	public static void create(int chr, boolean makePheno) throws IOException {
		String chrome = ext.chrome(chr);

		if (!new File("struct.dat").exists()) {
			System.err.println("Error - could not find "+"struct.dat"+" in current directory");
			System.err.println("Example:\n"+STRUCT_EXAMPLE);
			System.exit(2);
		}

		if (new File("mrkr"+chrome+".dat").exists()) {
			create("mrkr"+chrome+".dat", "chrom"+chrome+".pre", makePheno);
		} else if (new File("/home/npankrat/park/00masters/mrkr"+chrome+".dat").exists()) {
			create("/home/npankrat/park/00masters/mrkr"+chrome+".dat", "chrom"+chrome+".pre", makePheno);
		} else {
			System.err.println("Could not open the genotype file for chromosome "+chr+": mrkr"+chrome+".dat\n"+"in present directory or /home/npankrat/park/00masters/");
		}
	}

	public static void create(String filein, String fileout, boolean makePheno) {
		BufferedReader reader;
		PrintWriter writer, phenophile = null;
		String[] line = null, data;
		Hashtable<String,String[]> hash = new Hashtable<String,String[]>();

		try {
			reader = new BufferedReader(new FileReader(filein));
			if (!reader.ready()) {
				System.err.println("Error - no data in "+filein);
				System.exit(4);
			}
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				hash.put(line[0]+"\t"+line[1], line);
			}
			line = new String[line.length];
			line[0] = line[1] = "blank";
			for (int i = 2; i<line.length; i++) {
				line[i] = "0";
			}
			hash.put(line[0], line);
	
			writer = new PrintWriter(new FileWriter(fileout));
			if (makePheno) {
				phenophile = new PrintWriter(new FileWriter("pheno.sibs"));
			}
			reader.close();
	
			reader = new BufferedReader(new FileReader("struct.dat"));
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				for (int i = 0; i<6; i++) {
					writer.print((i==0?"":"\t")+line[i]);
				}
				data = hash.get(line[0]+"\t"+line[1]);
				if (data==null&&line[0].contains("_")) {
					data = hash.get(line[0].substring(0, line[0].indexOf("_"))+"\t"+line[1]);
				}
				if (data==null) {
					if (!line[6].equals("NoDNAprsnt")&&!line[6].equals("0")) {
						System.err.println("Error - expecting genotypes for "+line[0]+"\t"+line[1]+" (i.e. contains DNA# '"+line[6]+"'), but was not found in genotype file");
						System.exit(6);
					}
					data = hash.get("blank");
				}
				if (makePheno&&line.length>7) {
					phenophile.println(line[0]+"\t"+line[1]+"\t"+line[7]);
				}
				for (int i = 2; i<data.length; i++) {
					writer.print("\t"+data[i]);
				}
				writer.println();
			}
	
			reader.close();
			writer.close();
			if (makePheno) {
				phenophile.close();
			}
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filein + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filein + "\"");
			System.exit(2);
		}

	}

	public static int filterMarkers(String pedin, String pedout, String mapin, String mapout, String markersToKeepFile, String markersToDropFile) {
		return filterMarkers(pedin, pedout, mapin, mapout, markersToKeepFile == null?null:HashVec.loadFileToStringArray(markersToKeepFile, false, new int[] {0}, true), markersToDropFile == null?null:HashVec.loadFileToStringArray(markersToDropFile, false, new int[] {0}, true));
	}
	
	public static int filterMarkers(String pedin, String pedout, String mapin, String mapout, String[] markersToKeep, String[] markersToDrop) {
		boolean[] keeps;
		int numDropped;
		
		String[] markerNames;
		LinkageMap map;

		if ((markersToKeep == null || markersToKeep.length == 0) && (markersToDrop == null || markersToDrop.length == 0)) {
			System.err.println("Error - nothing to filter if both the keeps list and the drops list are null or empty");
			return -1;
		}
		
		map = new LinkageMap(mapin);
		markerNames = map.getMarkerNames();
		keeps = Array.booleanArray(markerNames.length, false);		
		for (int i = 0; i<markerNames.length; i++) {
			if (markersToKeep != null && ext.indexOfStr(markerNames[i], markersToKeep) >= 0 && markersToDrop != null && ext.indexOfStr(markerNames[i], markersToDrop) >= 0) {
				System.err.println("Error - marker '"+markerNames[i]+"' is listed in the keeps list and the drops list");
			} else if (markersToKeep == null) {
				keeps[i] = ext.indexOfStr(markerNames[i], markersToDrop) == -1; 
			} else if (markersToDrop == null) {
				keeps[i] = ext.indexOfStr(markerNames[i], markersToKeep) >= 0;
			} else {
				keeps[i] = ext.indexOfStr(markerNames[i], markersToDrop) == -1 && ext.indexOfStr(markerNames[i], markersToKeep) >= 0;
			}
        }

		numDropped = keeps.length - Array.booleanArraySum(keeps);
		
		if (numDropped > 0) {
			if (pedout.equals(pedin)) {
				pedin = Files.backup(pedin, "", "");
			}

			filterChromosome(pedin, pedout, keeps);
			map.filter(keeps);
			map.createFile(mapout);
		} else {
			if (!pedout.equals(pedin)) {
				Files.copyFile(pedin, pedout);
			}
			if (!mapout.equals(mapin)) {
				Files.copyFile(mapin, mapout);
			}
		}
		
		return numDropped;
	}
	
	public static void filterOutMarkersFromGenome(String dropFile, String prefix) {
		CountVector[] counts = CountVector.initArray(3);
		int numDropped;

		for (int chr = 1; chr<=23; chr++) {
			if (new File("map"+ext.chrome(chr)+".dat").exists()) {
				numDropped = filterMarkers(prefix+ext.chrome(chr)+".pre", prefix+ext.chrome(chr)+".pre", "map"+ext.chrome(chr)+".dat", "map"+ext.chrome(chr)+".dat", null, dropFile); 
				if (numDropped>0) {
					for (int i = 0; i<numDropped; i++) {
						counts[0].add(chr+"");
					}
				} else {
					counts[1].add(chr+"");
				}
			} else {
				counts[2].add(chr+"");
			}
		}

		for (int i = 0; i<counts.length; i++) {
			System.out.println("Chromosomes "+(i==0?"altered":(i==1?"unaltered":"missing"))+": "+(counts[i].getSize()==0?"none":Array.toStr(i==0?counts[i].list():counts[i].getValues(), ", ")));
		}
	}

	public static void filterChromosome(String filein, String fileout, boolean[] keeps) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;

		try {
			reader = new BufferedReader(new FileReader(filein));
			writer = new PrintWriter(new FileWriter(fileout));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (line.length!=keeps.length*2+6) {
					System.err.println("I do not think you are keeping the markers you think you are keeping...");
					System.exit(1);
				}
				writer.print(Array.toStr(Array.subArray(line, 0, 6)));
				for (int i = 0; i<keeps.length; i++) {
					if (keeps[i]) {
						writer.print("\t"+line[6+2*i+0]+"\t"+line[6+2*i+1]);
					}
				}
				writer.println();
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filein+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filein+"\"");
			System.exit(2);
		}
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		int chr = -1;
		String filename = "mrkr02.dat";
		String fileout = (filename.lastIndexOf(".")!=-1?filename.substring(0, filename.lastIndexOf("."))+".pre":filename+".pre");
		boolean pheno = false;
//		String drops = "drops.dat";
		String drops = "";
		String prefix = "re_chrom";

		String usage = "\n"+"link.LinkageFormat requires 0-2 arguments\n"+"   (1) chromosome number (i.e. chr=2 (no default))\n"+"   (2) add '-pheno' if you want to create \"pheno.sibs\" (not done by default)\n"+"   (3) file with list of markers you want to drop (i.e. drop=drops.dat (no default))\n"+"alternative [default] usage:\n"+"   (1) take from file other than a mrkr##.dat file (i.e. file="+filename+" (default))\n"+"   (2) specify the name of the .pre file (i.e. out="+fileout+" (default))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("chr=")) {
				chr = Positions.chromosomeNumber(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].equals("-pheno")) {
				pheno = true;
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("drop=")) {
				drops = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				fileout = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (drops!=null&&!drops.equals("")) {
				filterOutMarkersFromGenome(drops, prefix);
			} else if (chr==-1) {
				System.out.println("Attempting to load genotypes from '"+filename+"'");
				create(filename, fileout, pheno);
			} else {
				create(chr, pheno);
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("Error in processing chromosome "+chr);
		}
	}
}
