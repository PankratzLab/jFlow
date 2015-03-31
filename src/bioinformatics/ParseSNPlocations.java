// -Xmx4g
// -Xmx1024M
// found b132 in ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/database/organism_data/
package bioinformatics;

import java.io.*;
import java.util.Hashtable;

import common.*;
import filesys.SerialHash;
import filesys.SnpMarkerSet;

public class ParseSNPlocations {
	public static final String DEFAULT_NCBI_DIRECTORY = "C:/bin/NCBI/";
//	public static final String DEFAULT_NCBI_DIRECTORY = "N:/statgen/NCBI/";

	public static final String DEFAULT_B36_SOURCE = DEFAULT_NCBI_DIRECTORY+"b130_SNPChrPosOnRef_36_3.bcp.gz";
	public static final String DEFAULT_B36_DB = DEFAULT_NCBI_DIRECTORY+"b130_36_3.ser";

	public static final String DEFAULT_B37_SOURCE = DEFAULT_NCBI_DIRECTORY+"b138_SNPChrPosOnRef.bcp.gz";
	public static final String DEFAULT_B37_DB = DEFAULT_NCBI_DIRECTORY+"b138_37_3.ser";

	public static final String DEFAULT_MERGE_SOURCE = DEFAULT_NCBI_DIRECTORY+"RsMergeArch.bcp.gz";
	public static final String DEFAULT_MERGE = DEFAULT_NCBI_DIRECTORY+"RsMerge.ser";
	
	
	public static final int OFFSET = 1;
	public static final int MULTIPLE_POSITIONS = -2;
	public static final int UNMAPPED = -9;

	public static void lowMemParse(String snpListFile, boolean useExistingPositions, Logger log) {
		lowMemParse(snpListFile, ParseSNPlocations.DEFAULT_B37_DB, ParseSNPlocations.DEFAULT_MERGE, useExistingPositions, log);
	}
	
	public static void lowMemParse(String snpListFile, String db, String mergeDB, boolean useExistingPositions, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		SnpMarkerSet dbMarkerSet;
		int[] dbRSnumbers, dbPositions;
		byte[] dbChrs;
		int index;
		Hashtable<Integer,Integer> mergeHash;
		Integer trav, next;
		boolean first;
		byte chr;
		int position;

		log.report(ext.getTime());

		dbMarkerSet = null;
		dbRSnumbers = null;
		dbPositions = null;
		dbChrs = null;
		mergeHash = null;
		try {
	        reader = new BufferedReader(new FileReader(snpListFile));
	        writer = new PrintWriter(new FileWriter(ext.rootOf(snpListFile, false)+"_positions.xln"));
			writer.println("SNP\tChr\tPosition");
			first = true;
	        while (reader.ready()) {
	        	line = reader.readLine().trim().split("[\\s]+");
	        	chr = -2;
	        	position = -2;
	        	if (first) {
	        		if (line[0].equalsIgnoreCase("snp") || line[0].toLowerCase().startsWith("marker")) {
	    	        	line = reader.readLine().trim().split("[\\s]+");
	        		}
	        		first = false;
	        	}
	        	if (line.length > 2 && !line[1].equals(".") && !line[2].equals(".")) {
					writer.println(line[0]+"\t"+line[1]+"\t"+line[2]);
	        	} else {
	        		if (dbMarkerSet == null) {
	        			log.report("Loading database...");
	        			log.report("(if file is not found and it's searching for the wrong file, then force a recompile of MapSNPsAndGenes so that it gets the new constant from ParseSNPlocations)");

	        			dbMarkerSet = SnpMarkerSet.load(db, false, log);
	        			dbRSnumbers = dbMarkerSet.getRSnumbers();
	        			dbPositions = dbMarkerSet.getPositions();
	        			dbChrs = dbMarkerSet.getChrs();

	        			log.report("Searching database...");
	        			
	        		}
	        		if (line[0].startsWith("rs")) {
						index = Array.binarySearch(dbRSnumbers, Integer.parseInt(line[0].substring(2)), true);
						if (index == -1) {
							if (mergeDB != null) {
								if (mergeHash == null) {
									log.report("Could not find a marker, loading merge database...", false, true);
									if (new File(mergeDB).exists()) {
										mergeHash = SerialHash.loadSerializedIntHash(mergeDB);
										log.report("done");
									} else {
										log.report("failed; "+mergeDB+" not found");
										mergeHash = new Hashtable<Integer,Integer>();
									}
								}
								trav = Integer.parseInt(line[0].substring(2));
								next = trav;
								while (next != null) {
									trav = next;
									next = mergeHash.get(trav);
								}
								if (line[0].equals("rs"+trav)) {
									log.reportError("\n\n****ERROR**** failed to find "+line[0]+" in any NCBI database ****ERROR****\n\n");
								} else if (trav != null) {
									log.reportError("FYI - "+line[0]+" has merged with rs"+trav);
									index = Array.binarySearch(dbRSnumbers, trav.intValue(), true);
									chr = dbChrs[index];
									position = dbPositions[index]+ParseSNPlocations.OFFSET;
									index = Array.binarySearch(dbRSnumbers, trav.intValue(), true);
								} else {
									log.reportError("Error - tried and failed to find a merge record for "+line[0]);
								}
							}
							if (index == -1) {
								log.reportError("Error - could not find "+line[0]+" in "+db);
								chr = (byte)0;
								position = 0;
							}
						} else {
							if (dbChrs[index] == ParseSNPlocations.UNMAPPED) {
								log.reportError("Error - marker "+line[0]+" is not mapped to a chromosome");
								chr = 0;
								position = -1;
							} else if (dbPositions[index] == ParseSNPlocations.MULTIPLE_POSITIONS) {
								log.reportError("Warning - marker "+line[0]+" likely has multiple positions on chromosome "+dbChrs[index]);
								chr = dbChrs[index];
								position = -1;
							} else {
								chr = dbChrs[index];
								position = dbPositions[index]+ParseSNPlocations.OFFSET;
							}
						}
						writer.println(line[0]+"\t"+chr+"\t"+position);
					} else {
						log.reportError("Error - can't look up a SNP without an rs number ("+line[0]+")");
						writer.println(line[0]+"\t0\t0");
					}
	        	}
	        }
	        reader.close();
	        writer.close();
        } catch (FileNotFoundException fnfe) {
        	log.reportError("Error: file \""+snpListFile+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
        	log.reportError("Error reading file \""+snpListFile+"\"");
	        System.exit(2);
        }

		log.report("Done with the positions! "+ext.getTime());
	}

	public static void hashParse(String snpListFile, String db, String mergeDB, boolean useExistingPositions, Logger log) {
//		boolean reminderToCode = true;
		// TODO Auto-generated catch block
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		SnpMarkerSet dbMarkerSet;
		int[] dbRSnumbers, dbPositions;
		byte[] dbChrs;
		int index;
		Hashtable<Integer,Integer> mergeHash;
		Integer trav, next;
		boolean first;
		byte chr;
		int position;

		log.report(ext.getTime());

		dbMarkerSet = null;
		dbRSnumbers = null;
		dbPositions = null;
		dbChrs = null;
		mergeHash = null;
		try {
	        reader = new BufferedReader(new FileReader(snpListFile));
	        writer = new PrintWriter(new FileWriter(ext.rootOf(snpListFile, false)+"_positions.xln"));
			writer.println("SNP\tChr\tPosition");
			first = true;
	        while (reader.ready()) {
	        	line = reader.readLine().trim().split("[\\s]+");
	        	chr = -2;
	        	position = -2;
	        	if (first) {
	        		if (line[0].equalsIgnoreCase("snp") || line[0].toLowerCase().startsWith("marker")) {
	    	        	line = reader.readLine().trim().split("[\\s]+");
	        		}
	        		first = false;
	        	}
	        	if (line.length > 2 && !line[1].equals(".") && !line[2].equals(".")) {
					writer.println(line[0]+"\t"+line[1]+"\t"+line[2]);
	        	} else {
	        		if (dbMarkerSet == null) {
	        			log.report("Loading database...");

	        			dbMarkerSet = SnpMarkerSet.load(db, false);
	        			dbRSnumbers = dbMarkerSet.getRSnumbers();
	        			dbPositions = dbMarkerSet.getPositions();
	        			dbChrs = dbMarkerSet.getChrs();

	        			log.report("Searching database...");
	        			
	        		}
	        		if (line[0].startsWith("rs")) {
						index = Array.binarySearch(dbRSnumbers, Integer.parseInt(line[0].substring(2)), true);
						if (index == -1) {
							if (mergeDB != null) {
								if (mergeHash == null) {
									log.report("Could not find a marker, loading merge database...", false, true);
									if (new File(mergeDB).exists()) {
										mergeHash = SerialHash.loadSerializedIntHash(mergeDB);
										log.report("done");
									} else {
										log.report("failed; "+mergeDB+" not found");
										mergeHash = new Hashtable<Integer,Integer>();
									}
								}
								trav = Integer.parseInt(line[0].substring(2));
								next = trav;
								while (next != null) {
									trav = next;
									next = mergeHash.get(trav);
								}
								if (line[0].equals("rs"+trav)) {
									log.reportError("\n\n****ERROR**** failed to find "+line[0]+" in any NCBI database ****ERROR****\n\n");
								} else if (trav != null) {
									log.reportError("FYI - "+line[0]+" has merged with rs"+trav);
									index = Array.binarySearch(dbRSnumbers, trav.intValue(), true);
									chr = dbChrs[index];
									position = dbPositions[index]+ParseSNPlocations.OFFSET;
									index = Array.binarySearch(dbRSnumbers, trav.intValue(), true);
								} else {
									log.reportError("Error - tried and failed to find a merge record for "+line[0]);
								}
							}
							if (index == -1) {
								log.reportError("Error - could not find "+line[0]+" in "+db);
								chr = (byte)0;
								position = 0;
							}
						} else {
							if (dbChrs[index] == ParseSNPlocations.UNMAPPED) {
								log.reportError("Error - marker "+line[0]+" is not mapped to a chromosome");
								chr = 0;
								position = -1;
							} else if (dbPositions[index] == ParseSNPlocations.MULTIPLE_POSITIONS) {
								log.reportError("Warning - marker "+line[0]+" likely has multiple positions on chromosome "+dbChrs[index]);
								chr = dbChrs[index];
								position = -1;
							} else {
								chr = dbChrs[index];
								position = dbPositions[index]+ParseSNPlocations.OFFSET;
							}
						}
						writer.println(line[0]+"\t"+chr+"\t"+position);
					} else {
						log.reportError("Error - can't look up a SNP without an rs number ("+line[0]+")");
						writer.println(line[0]+"\t0\t0");
					}
	        	}
	        }
	        reader.close();
	        writer.close();
        } catch (FileNotFoundException fnfe) {
        	log.reportError("Error: file \""+snpListFile+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
        	log.reportError("Error reading file \""+snpListFile+"\"");
	        System.exit(2);
        }

		log.report("Done with the positions! "+ext.getTime());
	}
	
	public static void createFromSource(String source, String db) {
		BufferedReader reader;
		String[] line;
		int count, countPARys;
		String[] rsNumbers;
		int[] positions;
		byte[] chrs;
		String temp;

		try {
			System.out.println("Reading in "+source);
//			reader = new BufferedReader(new FileReader(source));
			reader = Files.getAppropriateReader(source);
			count = 0;
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
//				if (line.length>2) {
				if (line.length>1) {
					count++;
				}
			}
			reader.close();
			System.out.println("Found "+count+" markers");

//			reader = new BufferedReader(new FileReader(source));
			reader = Files.getAppropriateReader(source);
			rsNumbers = new String[count];
			positions = new int[count];
			chrs = new byte[count];
			count = 0;
			countPARys = 0;
			while (reader.ready()) {
				temp = reader.readLine();
				try {
					line = temp.trim().split("[\\s]+");
//					if (line.length>2) {
					if (line.length>1) {
						rsNumbers[count] = "rs"+line[0];
						chrs[count] = Positions.chromosomeNumber(line[1]);
						if (line[1].equals("PAR") && line[2].equals("y")) {
							countPARys++;
						} else if (line.length > 2) { // new
							try {
								positions[count] = Integer.parseInt(line[2]);
							} catch (NumberFormatException nfe) {
								System.err.println("Error - parsing line #"+(count+1)+": "+Array.toStr(line));
							}
						} else {
							positions[count] = MULTIPLE_POSITIONS;
						}
						count++;
					} else { // new
						chrs[count] = 0;
						positions[count] = Integer.parseInt(line[2]);
						count++;
					}
				} catch (Exception e) {
					System.err.println("Error parsing "+source+" at the following line:");
					System.err.println(temp);
					e.printStackTrace();
					reader.close();
					return;
				}
				
			}
			reader.close();
			if (countPARys > 0) {
				System.err.println("Warning - there were "+countPARys+" instances where the chromosome was listed as PAR and the position was nothing but 'y'");
			}
			// already in order, d'uh!
//			System.out.println("Sorting...");
//			keys = Sort.quicksort(rsNumbers);
			System.out.println("Serializing...");
			new SnpMarkerSet(rsNumbers, chrs, positions, null, null, false, true).serialize(db);
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+source+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+source+"\"");
			System.exit(2);
		}
	}

	public static void createMergeDBfromSource(String source, String db) {
		BufferedReader reader;
		String[] line;
		int count;
		Hashtable<Integer,Integer> hash;

		count = 0;
		try {
			System.out.println(ext.getTime()+"\tReading in "+source);
			reader = new BufferedReader(new FileReader(source));
			while (reader.ready()) {
				reader.readLine();
				count++;
			}
			reader.close();
			System.out.println(ext.getTime()+"\tFound "+count+" merged markers");
			
			hash = new Hashtable<Integer,Integer>(count);

			count = 0;
			reader = Files.getAppropriateReader(source);
			while (reader.ready()) {
				line = reader.readLine().trim().split("\\t", -1);
//				hash.put(line[0], line[6]); // not complete
				hash.put(Integer.parseInt(line[0]), Integer.parseInt(line[1]));
				count++;
			}
			reader.close();
			System.out.println(ext.getTime()+"\tSerializing...");
			SerialHash.createSerializedIntHash(db, hash);
			System.out.println(ext.getTime()+"\tDone!");
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+source+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+source+"\"");
			System.exit(2);
		} catch (OutOfMemoryError oome) {
			System.err.println("Error - Ran out of memory at line "+(count+1));
			oome.printStackTrace();
		}
	}	

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String source = "";
		String mergeSource = "";
		String db = DEFAULT_B37_DB;
		String merge = DEFAULT_MERGE;

		// uncomment one of these to compile
//		source = DEFAULT_B37_SOURCE;
//		db = DEFAULT_B37_DB;

//		source = DEFAULT_B36_SOURCE;
//		db = DEFAULT_B36_DB;
		
//		mergeSource = DEFAULT_MERGE_SOURCE;

		// String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\PD-singleton\\";
		// String filename = "pd.recode.map";
		// String filename = "whathappened.dat";
		// boolean plinkFormat = true;

		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\";
//		String filename = "top_snplist_neck.txt";
//		String filename = "top_aa_snps.txt";
//		String filename = "aff.txt";
//		String filename = "CARE_hits.txt";
//		String filename = "CARE_X_hits.txt";
//		String filename = "X_diffs.txt";
		String filename = "list.txt";
		
		boolean plinkFormat = false;

		String usage = "\n"+
		"bioinformatics.ParseSNPlocations requires 0-1 arguments\n"+
		"   (0) file from which to create db (i.e. source="+DEFAULT_B37_SOURCE+" (not the default))\n"+
		"   (1) snpDB (i.e. db="+db+" (default))\n"+
		"   (2) merge DB (i.e. merge="+merge+" (default))\n"+
		"   (3) dir (i.e. dir="+dir+" (default))\n"+
		"   (4) filename (i.e. file="+filename+" (default))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("dir=")) {
				dir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("db=")) {
				db = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("merge=")) {
				merge = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (!source.equals("")) {
				createFromSource(source, db);
			} else if (!mergeSource.equals("")) {
				createMergeDBfromSource(mergeSource, merge);
			} else {
				SnpMarkerSet map = new SnpMarkerSet(dir+filename, plinkFormat?SnpMarkerSet.PLINK_MAP_FORMAT:SnpMarkerSet.NAMES_ONLY, true, new Logger());
				map.parseSNPlocations(db, merge, new Logger());
				map.writeToFile(dir+ext.rootOf(filename)+"_newPositions.out", SnpMarkerSet.GENERIC_FORMAT);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
