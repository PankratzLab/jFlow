package cnv.manage;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import cnv.filesys.*;
import common.*;

public class Markers {

	public static int[] orderMarkers(String[] markerNames, String markerDatabase, String output) {
		Hashtable<String,String> snpPositions;
		byte[] chrs;
		int[] positions, keys;
		String[] line;
		Vector<String> v;

		snpPositions = HashVec.loadFileToHashString(markerDatabase, 0, new int[] {1, 2}, "\t", false);
		if (markerNames == null) {
			// TODO
			// replace with Metal.MARKER_NAMES after merge
			String[] MARKER_NAMES = new String[] {"MarkerName", "Marker", "Name", "SNP"};
			for (int i = 0; i < MARKER_NAMES.length; i++) {
				if (snpPositions.containsKey(MARKER_NAMES[i])) {
					snpPositions.remove(MARKER_NAMES[i]);
				}
			}
			markerNames = HashVec.getKeys(snpPositions);
		}

		v = new Vector<String>();
		System.out.println(ext.getTime()+"\tOrdering markers");
		chrs = new byte[markerNames.length];
		positions = new int[markerNames.length];
		for (int i = 0; i<markerNames.length; i++) {
			if (snpPositions.containsKey(markerNames[i])) {
				line = snpPositions.get(markerNames[i]).split("[\\s]+");
				chrs[i] = Positions.chromosomeNumber(line[0]);
				positions[i] = Integer.parseInt(line[1]);
			} else {
				v.add(markerNames[i]);
			}
		}
		if (v.size() > 0) {
			System.err.println("Error - The there "+(v.size()==1?"was one":"were "+v.size())+" markers found in the FinalReport file that were not listed in the file of SNP positions; halting parse operation");
			Files.writeList(Array.toStringArray(v), ext.parseDirectoryOfFile(markerDatabase)+"markersNotInPositionsFile.txt");
			return null;			
		}

		keys = Sort.orderTwoLayers(chrs, positions);

		new MarkerSet(markerNames, chrs, positions, keys).serialize(output);

		System.out.println(ext.getTime()+"\tFinished ordering");

		return keys;
	}

	public static void generateMarkerPositions(Project proj, String snpTable) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		int[] indices;
		String delimiter;
		long time;
	
		time = new Date().getTime();
		delimiter = proj.getSourceFileDelimiter();
		try {
			reader = new BufferedReader(new FileReader(proj.getProjectDir()+snpTable));
			writer = new PrintWriter(new FileWriter(proj.getFilename(Project.MARKER_POSITION_FILENAME, true, true)));
			indices = ext.indexFactors(ParseIllumina.SNP_TABLE_FIELDS, reader.readLine().trim().split(delimiter), false, true, true, true);
			writer.println("Marker\tChr\tPosition");
			while (reader.ready()) {
				line = reader.readLine().trim().split(delimiter);
				writer.println(line[indices[0]]+"\t"+line[indices[1]]+"\t"+line[indices[2]]);
			}
			writer.close();
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+snpTable+"\" not found in "+proj.getProjectDir());
			return;
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+snpTable+"\"");
			return;
		}
		System.out.println("Finished parsing "+proj.getFilename(Project.MARKER_POSITION_FILENAME, false, false)+" in " + ext.getTimeElapsed(time));
	}

	public static void useAlleleLookup(String filename, int alleleCol, String lookupFile, int setFrom, int setTo) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Hashtable<String, String> hash;
		String[] header, alleles;
		
		header = Files.getHeaderOfFile(lookupFile, "\t", new Logger(null));
		System.out.println("Converting from columns "+header[1+setFrom*2+0]+"/"+header[1+setFrom*2+1]+" to columns "+header[1+setTo*2+0]+"/"+header[1+setTo*2+1]);
		hash = HashVec.loadFileToHashString(lookupFile, 0, Array.intArray(header.length), "\t", true);
		
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_"+header[1+setTo*2+0]+"_"+header[1+setTo*2+1]+".xln"));
			line = reader.readLine().trim().split("[\\s]+");
			line = Array.insertStringAt(header[1+setTo*2+0], line, alleleCol+2);
			line = Array.insertStringAt(header[1+setTo*2+1], line, alleleCol+3);
			writer.println(Array.toStr(line));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (hash.containsKey(line[0])) {
					alleles = Array.subArray(hash.get(line[0]).split("[\\s]+"), 1);
					if (line[alleleCol].equals(alleles[setFrom*2+0]) && line[alleleCol+1].equals(alleles[setFrom*2+1])) {
						line = Array.insertStringAt(alleles[setTo*2+0], line, alleleCol+2);
						line = Array.insertStringAt(alleles[setTo*2+1], line, alleleCol+3);
					} else if (line[alleleCol].equals(alleles[setFrom*2+1]) && line[alleleCol+1].equals(alleles[setFrom*2+0])) {
						line = Array.insertStringAt(alleles[setTo*2+1], line, alleleCol+2);
						line = Array.insertStringAt(alleles[setTo*2+0], line, alleleCol+3);
					} else {
						System.err.println("Error - snp '"+line[0]+"' has alleles "+line[alleleCol]+"/"+line[alleleCol+1]+" in the file and "+alleles[setFrom*2+0]+"/"+alleles[setFrom*2+1]+" in allele lookup table");
						line = Array.insertStringAt("XXXXX", line, alleleCol+2);
						line = Array.insertStringAt("XXXXX", line, alleleCol+3);
					}
				} else {
					System.err.println("Error - snp '"+line[0]+"' not found in allele lookup table");
					line = Array.insertStringAt("XXXXX", line, alleleCol+2);
					line = Array.insertStringAt("XXXXX", line, alleleCol+3);
				}				
				writer.println(Array.toStr(line));
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
		Project proj;
		String filename = Project.DEFAULT_PROJECT;
		String snpTable = "";
		String fileToConvert = "";
		String lookupFile = "alleleLookup.txt";
		int alleleCol = 6;
		int setFrom = 1;
		int setTo = 2;

		String usage = "\n" +
		"cnv.manage.Markers requires 0-1 arguments\n" +
		"   (1) project file (i.e. proj="+filename+" (default))\n"+
		"   (2) filename of SNP Table (i.e. snps=Table.csv (not the default))\n"+
		" OR:\n"+
		"   (2) use allele lookup to convert a file form forward to TOP strand (i.e. convert=file.txt (not the default))\n"+
		"   (3) column of A1 in file (i.e. col="+alleleCol+" (default))\n"+
		"   (4) allele set to lookup from (i.e. from="+setFrom+" (default; if 1+8 columns fo 4 pairs, then use indices 0-3))\n"+
		"   (5) allele set to lookup to (i.e. to="+setTo+" (default))\n"+
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("snps=")) {
				snpTable = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

		proj = new Project(filename, false);

		try {
			if (!snpTable.equals("")) {
				Markers.generateMarkerPositions(proj, snpTable);
			} else if (!fileToConvert.equals("")) {
				Markers.useAlleleLookup(fileToConvert, alleleCol, lookupFile, setFrom, setTo);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
