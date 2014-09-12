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
	public static final int MAX_ERRORS_TO_REPORT = 30;

	public static int[] orderMarkers(String[] markerNames, String markerDatabase, String output, Logger log) {
		Hashtable<String,String> snpPositions;
		byte[] chrs;
		int[] positions, keys;
		String[] line;
		Vector<String> v;
		long time;
		int logLevel;

		logLevel = log.getLevel();
		log.setLevel(9);
		time = new Date().getTime();
		snpPositions = loadFileToHashString(markerDatabase, log);
		if (snpPositions == null) {
			return null;
		}
		if (markerNames == null) {
			for (int i = 0; i < Aliases.MARKER_NAMES.length; i++) {
				if (snpPositions.containsKey(Aliases.MARKER_NAMES[i])) {
					snpPositions.remove(Aliases.MARKER_NAMES[i]);
				}
			}
			markerNames = HashVec.getKeys(snpPositions);
		}

		v = new Vector<String>();
		log.report(ext.getTime()+"\tSorting markers by chromosome and position");
		chrs = new byte[markerNames.length];
		positions = new int[markerNames.length];
		for (int i = 0; i<markerNames.length; i++) {
			if (snpPositions.containsKey(markerNames[i])) {
				line = snpPositions.get(markerNames[i]).split("[\\s]+");
				chrs[i] = Positions.chromosomeNumber(line[0], log);
				positions[i] = Integer.parseInt(line[1]);
			} else {
				v.add(markerNames[i]);
			}
		}
		if (v.size() > 0) {
			log.reportError("Error - There "+(v.size()==1?"was one":"were "+v.size())+" markers found in the FinalReport file that were not listed in the file of marker positions; halting parse operation");
			log.reportError("\nThe best source of complete marker positions is the SNP manifest (e.g., SNP_Map.csv from Illumina's GenomeStudio that should be exported along with the FinalReport files)");
			Files.writeList(Array.toStringArray(v), ext.parseDirectoryOfFile(markerDatabase)+"markersNotInPositionsFile.txt");
			return null;			
		}

		keys = Sort.orderTwoLayers(chrs, positions, log);

		new MarkerSet(markerNames, chrs, positions, keys).serialize(output);

		log.report(ext.getTime()+"\tFinished sorting in " + ext.getTimeElapsed(time));
		log.setLevel(logLevel);
		
		return keys;
	}
	
	public static Hashtable<String,String> loadFileToHashString(String filename, Logger log) {
		BufferedReader reader = null;
		String[] line;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		String markerName, chr, position, delimiter, temp;
		byte chrValue;
		int count, countBad, numBlankNames, numBlankChrs, numBlankPositions, numRepeatedNames, numInvalidChrs, numInvalidPositions, numIncompleteLines;
		
		delimiter = Files.determineDelimiter(filename, log);

		count = countBad = 0;
		numBlankNames = numBlankChrs = numBlankPositions = numRepeatedNames = numInvalidChrs = numInvalidPositions = numIncompleteLines = 0;
		try {
			reader = Files.getAppropriateReader(filename);
			while (reader.ready()) {
				temp = reader.readLine();
				if (delimiter.equals(",")) {
					line = ext.splitCommasIntelligently(temp, true, new Logger());
				} else if (temp.indexOf("\t") == -1) {
					line = temp.trim().split("[\\s]+");
				} else {
					line = temp.split("\t", -1);
				}
				if (count == 0 && ext.indexOfStr(line[0], Aliases.MARKER_NAMES) >= 0) {
					
				} else if (line.length < 3) {
					if (countBad < MAX_ERRORS_TO_REPORT) {
						log.report("Error - incomplete line at row "+count+" for marker \""+line[0]+"\"; line will not be added");
					}
					numIncompleteLines++;
				} else {
					markerName = line[0];
					chr = line[1];
					position = line[2];
					
					if (markerName.equals("")) {
						if (countBad < MAX_ERRORS_TO_REPORT) {
							log.reportError("Error - blank marker name at line "+count+" of "+filename);
						}
						numBlankNames++;
						countBad++;
					} else if (chr.equals("")) {
						if (countBad < MAX_ERRORS_TO_REPORT) {
							log.reportError("Error - blank chr for marker '"+markerName+"' at line "+count+" of "+filename);
						}
						numBlankChrs++;
						countBad++;
					} else if (position.equals("")) {
						if (countBad < MAX_ERRORS_TO_REPORT) {
							log.reportError("Error - blank position for marker '"+markerName+"' at line "+count+" of "+filename);
						}
						numBlankPositions++;
						countBad++;
					} else {
						if (hash.containsKey(markerName)) {
							log.reportError("Error - marker '"+markerName+"' is already listed in the markerPositions file and is seen again at line "+count+" of "+filename);
							numRepeatedNames++;
							countBad++;
						}
						chrValue = Positions.chromosomeNumber(chr, log);
						if (chrValue < 0 || chrValue > 26) {
							numInvalidChrs++;
							countBad++;
						}
						try {
							Integer.parseInt(position);
						} catch (NumberFormatException nfe) {
							if (countBad < MAX_ERRORS_TO_REPORT) {
								log.reportError("Error - invalid position ("+position+") for marker '"+markerName+"' at line "+count+" of "+filename);
							}
							numInvalidPositions++;
							countBad++;
						}
					}
					
					hash.put(markerName, chr+"\t"+position);
				}
				if (countBad == MAX_ERRORS_TO_REPORT) {
					log.reportError("...");
					countBad++;
				}
				count++;
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}

		log.report("\nRead in "+count+" markers from the markerPositions file");
		if (countBad > MAX_ERRORS_TO_REPORT) {
			countBad--;
		}

		if (countBad > 0) {
			log.report("...with a total of "+countBad+" problem"+(countBad==1?"":"s"));
		}
		if (numIncompleteLines > 0) {
			log.report("...including "+numIncompleteLines+" incomplete line"+(countBad==1?"":"s"));
		}
		if (numBlankNames > 0) {
			log.report("Number of blank marker names: "+numBlankNames);
		}
		if (numBlankChrs > 0) {
			log.report("Number of blank chromosomes: "+numBlankChrs);
		}
		if (numBlankPositions > 0) {
			log.report("Number of blank marker positions: "+numBlankPositions);
		}
		if (numRepeatedNames > 0) {
			log.report("Number of repeated marker names: "+numRepeatedNames);
		}
		if (numInvalidChrs > 0) {
			log.report("Number of invalid chromosomes: "+numInvalidChrs);
		}
		if (numInvalidPositions > 0) {
			log.report("Number of invalid positions: "+numInvalidPositions);
		}
		log.report("");
		
		return hash;
	}	

	public static void generateMarkerPositions(Project proj, String snpTable) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		int[] indices;
		String delimiter;
		long time;
		Logger log;
		
		log = proj.getLog();
		time = new Date().getTime();
		delimiter = Files.determineDelimiter(snpTable, log);
		try {
			if (!Files.exists(snpTable) && Files.exists(proj.getProjectDir()+snpTable)) {
				snpTable = proj.getProjectDir()+snpTable;
			}
			reader = Files.getAppropriateReader(snpTable);
			writer = new PrintWriter(new FileWriter(proj.getFilename(Project.MARKER_POSITION_FILENAME, false, false)));
			indices = ext.indexFactors(ParseIllumina.SNP_TABLE_FIELDS, reader.readLine().trim().split(delimiter), false, true, true, true);
			writer.println("Marker\tChr\tPosition");
			while (reader.ready()) {
				line = reader.readLine().trim().split(delimiter);
				writer.println(line[indices[0]]+"\t"+line[indices[1]]+"\t"+line[indices[2]]);
			}
			writer.close();
			reader.close();
		} catch (FileNotFoundException fnfe) {
			proj.message("Error: file \""+snpTable+"\" not found in "+proj.getProjectDir());
			return;
		} catch (IOException ioe) {
			proj.message("Error reading file \""+snpTable+"\"");
			return;
		}
		log.report("Finished parsing "+proj.getFilename(Project.MARKER_POSITION_FILENAME, false, false)+" in " + ext.getTimeElapsed(time));
	}

	public static void useAlleleLookup(String filename, int alleleCol, String lookupFile, int setFrom, int setTo) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Hashtable<String, String> hash;
		String[] header, alleles;
		
		header = Files.getHeaderOfFile(lookupFile, "\t", new Logger());
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
		String filename = null;
		String snpTable = "";
		String fileToConvert = "";
		String lookupFile = "alleleLookup.txt";
		int alleleCol = 6;
		int setFrom = 1;
		int setTo = 2;

		String usage = "\n" +
		"cnv.manage.Markers requires 0-1 arguments\n" +
		"   (1) project properties filename (i.e. proj="+cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n"+
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
