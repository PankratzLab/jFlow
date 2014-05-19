//-Xms2048M -Xmx2048M
package bioinformatics;

import java.io.*;

import common.*;

// TODO was not a clean merge, could be undetected issues
public class MapSNPsAndGenes {
//	public static final boolean AUTOMATICALLY_ADD_ONE = true;
//	public static final boolean AUTOMATICALLY_ADD_ONE = false;
	public static final String GENES36 = "N:/statgen/NCBI/genes36.xln";
	public static final String GENES37 = "N:/statgen/NCBI/genes37.xln";
//	public static final int DEFAULT_WIGGLE_ROOM = 0;
	public static final int DEFAULT_WIGGLE_ROOM = 15000;
//	public static final int DEFAULT_WIGGLE_ROOM = 100000;
	public static final int UCSC_WINDOW = DEFAULT_WIGGLE_ROOM;

//	public static final String DEFAULT_DIRECTORY = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\udall\\";
//	public static final String DEFAULT_FILENAME = "combinedMarkerPositions.xln";

	public static final String DEFAULT_DIRECTORY = "C:\\Documents and Settings\\npankrat\\My Documents\\";

	// public static final String DEFAULT_FILENAME = "COGA.txt";
	// public static final String DEFAULT_FILENAME = "FIA_snp_list.txt";
	// public static final String DEFAULT_FILENAME = "FIA_snps.xln";
	// public static final String DEFAULT_FILENAME = "top5000nathan.txt";
	// public static final String DEFAULT_FILENAME = "bone snps.txt";
	// public static final String DEFAULT_FILENAME = "top_snps_for_nathan_2009_04_16.txt";
	// public static final String DEFAULT_FILENAME = "top_aa_snps_positions.txt";
//	public static final String DEFAULT_FILENAME = "cpru.txt";
//	public static final String DEFAULT_FILENAME = "aff.txt";
//	public static final String DEFAULT_FILENAME = "onset.txt";
//	public static final String DEFAULT_FILENAME = "CARE_hits.txt";
//	public static final String DEFAULT_FILENAME = "CARE_X_hits.txt";
	public static final String DEFAULT_FILENAME = "SNPs.txt";
//	public static final String DEFAULT_FILENAME = "E3E3_LTE0.05.txt";
	
	public static String[] mapSNPsToGenes(int[][] markerPositions, int wiggleRoom, Logger log) {
		BufferedReader reader;
		String[] line, genes;
		int[] chr_start_stop;

		genes = Array.stringArray(markerPositions.length, "");
		try {
			reader = new BufferedReader(new FileReader(GENES36));
			reader.readLine();
			chr_start_stop = new int[3];
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				if (!line[3].equals(".")) {
					chr_start_stop[0] = line[2].equals("X")?23:(line[2].equals("Y")?24:(line[2].equals("XY")?25:(line[2].equals("MT")?26:(line[2].equals("Un")?27:Integer.parseInt(line[2])))));
					chr_start_stop[1] = Integer.parseInt(line[3]);
					chr_start_stop[2] = Integer.parseInt(line[4]);
	
					for (int j = 0; j<markerPositions.length; j++) {
						if (markerPositions[j][0]==chr_start_stop[0]&&markerPositions[j][1]>=chr_start_stop[1]-wiggleRoom&&markerPositions[j][1]<=chr_start_stop[2]+wiggleRoom) {
							genes[j] += (genes[j].equals("")?"":"|")+line[1];
						}
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \""+GENES36+"\" not found");
			System.exit(1);
		} catch (IOException ioe) {
			log.reportError("Error reading file \""+GENES36+"\"");
			System.exit(2);
		}

		return genes;
	}

	public static String getGeneDB(byte build, Logger log) {
		String geneDB;
		
		geneDB = GENES37;
		switch (build) {
		case 36:
			geneDB = GENES36;
			break;
		case 37:
			geneDB = GENES37;
			break;
		default:
			log.reportError("Error - unknown build '"+build+"'; using the default instead (b37/hg19)");
			break;
		}

		return geneDB;
	}
	
	public static String getSNPDB(byte build, Logger log) {
		String snpDB;
		
		snpDB = ParseSNPlocations.DEFAULT_B37_DB;
		switch (build) {
		case 36:
			snpDB = ParseSNPlocations.DEFAULT_B36_DB;
			break;
		case 37:
			snpDB = ParseSNPlocations.DEFAULT_B37_DB;
			break;
		default:
			log.reportError("Error - unknown build '"+build+"'; using the default instead (b37/hg19)");
			break;
		}

		return snpDB;
	}
	
	public static String[][] mapSNPsToGenes(int[][] markerPositions, byte build, int wiggleRoom, Logger log) {
		return mapSNPsToGenes(markerPositions, getGeneDB(build, log), wiggleRoom, log);
	}
	
	/**
	 * Map chromosome positions to genes.
	 * @param markerPositions a two dimensional array of marker positions. An example of the format is
	 * 			0	1
	 * 		0	chr	pos
	 * 		1	chr	pos
	 * 		...
	 * @param geneDB
	 * @param wiggleRoom
	 * @param log
	 * @return
	 */
	public static String[][] mapSNPsToGenes(int[][] markerPositions, String geneDB, int wiggleRoom, Logger log) {
		BufferedReader reader;
		String[] line, genes, distances;
		int[] chr_start_stop;
		int dist;
		String[][] finalGenes;
		String[] geneNames;
		int[] dists;
		int[] order;
	
		genes = Array.stringArray(markerPositions.length, "");
		distances = Array.stringArray(markerPositions.length, "");
		try {
			reader = new BufferedReader(new FileReader(geneDB));
			reader.readLine();
			chr_start_stop = new int[3];
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				if (!line[3].equals(".")) {
					chr_start_stop[0] = line[2].equals("X")?23:(line[2].equals("Y")?24:(line[2].equals("XY")?25:(line[2].equals("MT")?26:(line[2].equals("Un")?27:Integer.parseInt(line[2])))));
					chr_start_stop[1] = Integer.parseInt(line[3]);
					chr_start_stop[2] = Integer.parseInt(line[4]);
	
					for (int j = 0; j<markerPositions.length; j++) {
						if (markerPositions[j][0]==chr_start_stop[0]&&markerPositions[j][1]>=chr_start_stop[1]-wiggleRoom&&markerPositions[j][1]<=chr_start_stop[2]+wiggleRoom) {
							genes[j] += (genes[j].equals("")?"":"|")+line[1];
							if (markerPositions[j][1]<chr_start_stop[1]) {
								dist = (line[5].equals("+")?-1:1)*(chr_start_stop[1] - markerPositions[j][1]);
							} else if (markerPositions[j][1]>chr_start_stop[2]) {
								dist = (line[5].equals("+")?1:-1)*(markerPositions[j][1] - chr_start_stop[2]);
							} else {
								dist = 0;
							}
	
//							distances[j] += (distances[j].equals("")?"":"|")+dist+","+line[5];
							distances[j] += (distances[j].equals("")?"":"|")+dist;
						}
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \""+geneDB+"\" not found");
			System.exit(1);
		} catch (IOException ioe) {
			log.reportError("Error reading file \""+geneDB+"\"");
			System.exit(2);
		}
		
		finalGenes = new String[markerPositions.length][3];
		for (int i = 0; i < markerPositions.length; i++) {
			if (genes[i].equals("")) {
				finalGenes[i][0] = ".";
				finalGenes[i][1] = ".";
				finalGenes[i][2] = ".";
			} else {
				geneNames = genes[i].split("\\|");
				dists = Array.toIntArray(distances[i].split("\\|"));
				for (int j = 0; j < geneNames.length; j++) {
					if (dists[j] < 0) {
						geneNames[j] += " "+ext.prettyUpDistance(-1*dists[j], 1)+" upstream";
					}
					if (dists[j] > 0) {
						geneNames[j] += " "+ext.prettyUpDistance(dists[j], 1)+" downstream";
					}
				}
				for (int j = 0; j < dists.length; j++) {
					dists[j] = Math.abs(dists[j]);
				}
				order = Sort.quicksort(dists);
				for (int j = 0; j < geneNames.length; j++) {
					if (dists[j] == 0) {
						finalGenes[i][0] = finalGenes[i][0] == null?geneNames[j]:finalGenes[i][0]+"|"+geneNames[j];
						finalGenes[i][1] = finalGenes[i][1] == null?geneNames[j]:finalGenes[i][1]+"|"+geneNames[j];
					}
				}
				if (finalGenes[i][1] == null) {
					finalGenes[i][1] = geneNames[order[0]];
				}
				finalGenes[i][2] = Array.toStr(Sort.putInOrder(geneNames, order), "|");
			}
		}
		
		return finalGenes;
	}


	public static void procSNPsToGenes(String dir, String snps, int wiggleRoom, byte build, Logger log) {
		PrintWriter writer;
		String[] line;
		String[] data, markers;
		String[][] genes;
		int[][] markerPositions;
		String snpDB, geneDB;

		snpDB = getSNPDB(build, log);
		geneDB = getGeneDB(build, log);
		ParseSNPlocations.lowMemParse(dir+snps, snpDB, ParseSNPlocations.DEFAULT_MERGE, true, log);
//		ParseSNPlocations.lowMemParse(dir+snps, true, log);
		
		data = Array.toStringArray(HashVec.loadFileToVec(ext.rootOf(dir+snps, false)+"_positions.xln", false, false, false));

		markers = new String[data.length-1];
		markerPositions = new int[markers.length][2];
		for (int i = 0; i<markers.length; i++) {
			line = data[i+1].trim().split("[\\s]+");
			markers[i] = line[0];
			markerPositions[i][0] = Positions.chromosomeNumber(line[1]);
			markerPositions[i][1] = Integer.parseInt(line[2]);
		}

//		genes = mapSNPsToGenes1_zx(markerPositions, wiggleRoom, log);
		genes = mapSNPsToGenes(markerPositions, geneDB, wiggleRoom, log);

		try {
			System.out.println(dir+ext.rootOf(snps)+"_genes.xln");
			writer = new PrintWriter(new FileWriter(dir+ext.rootOf(snps)+"_genes.xln"));
			writer.println("SNP\tChr\tPosition\tGene(s)\t" + UCSC_WINDOW + "\tClosest" + "\t<- dynamically linked basepair buffer in UCSC hyperlink");
			for (int i = 0; i<markers.length; i++) {
//				writer.println(markers[i]+"\t"+markerPositions[i][0]+"\t"+(AUTOMATICALLY_ADD_ONE?markerPositions[i][1]+1:markerPositions[i][1])+"\t"+(genes[i].equals("")?".":genes[i])+"\t=HYPERLINK(CONCATENATE(\"http://genome.ucsc.edu/cgi-bin/hgTracks?position=chr"+markerPositions[i][0]+":\","+markerPositions[i][1]+"-$E$1+1,\"-\","+markerPositions[i][1]+"+$E$1),CONCATENATE(\"chr"+markerPositions[i][0]+":\","+markerPositions[i][1]+"-$E$1+1,\"-\","+markerPositions[i][1]+"+$E$1))");
				writer.println(markers[i]+"\t"+markerPositions[i][0]+"\t"+markerPositions[i][1]+"\t"+(genes[i][0]==null?".":genes[i][0])+"\t"+(genes[i][1]==null?".":genes[i][1])+"\t"+(genes[i][2]==null?".":genes[i][2])+"\t=HYPERLINK(CONCATENATE(\"http://genome.ucsc.edu/cgi-bin/hgTracks?position=chr"+markerPositions[i][0]+":\","+markerPositions[i][1]+"-$E$1+1,\"-\","+markerPositions[i][1]+"+$E$1),CONCATENATE(\"chr"+markerPositions[i][0]+":\","+markerPositions[i][1]+"-$E$1+1,\"-\","+markerPositions[i][1]+"+$E$1))");
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing genes file");
			log.reportException(e);
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = DEFAULT_DIRECTORY;
		String filename = DEFAULT_FILENAME;
		int wiggleRoom = DEFAULT_WIGGLE_ROOM;
		String temp;
		byte build = 37;
		Logger log;

		String usage = "\\n"+
		"bioinformatics.MapSNPsAndGenes requires 0-1 arguments\n"+
		"   (1) directory (i.e. dir="+dir+" (default))\n"+
		"   (2) filename (i.e. file="+filename+" (default))\n"+
		"   (3) # bp up and down stream to count as an associated gene (i.e. win="+wiggleRoom+" (default))\n"+
		"   (4) build # of the NCBI gene map file (i.e. build="+build+" (default))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				dir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("win=")) {
				wiggleRoom = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("build=")) {
				build = Byte.parseByte(args[i].split("=")[1]);
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}

		try {
			if (filename.endsWith(".snps")) {
				log = new Logger(filename+".log");
				temp = filename.substring(0, filename.length()-(".snps").length());
				if (filename.indexOf(".b36.")>0) {
					build = 36;
				} else if (filename.indexOf(".b37.")>0) {
					build = 37;
				} else {
					log.reportError("Warning - using the default build (b37/hg19) since the file '"+filename+"' does not explicitly specify one using the convention \"filename.b37.20000.snps\"");
					build = 37;
				}
				if (temp.indexOf(".")>0) {
					try {
						wiggleRoom = Integer.parseInt(temp.substring(temp.lastIndexOf(".")+1));
					} catch (NumberFormatException nfe) {
						wiggleRoom = DEFAULT_WIGGLE_ROOM;
					}
				}				
//				procSNPsToGenes("", filename, wiggleRoom, build, new Logger(filename+".log"));
				procSNPsToGenes("", filename, wiggleRoom, build, log);
			} else {
//				procSNPsToGenes(dir, filename, wiggleRoom, build, new Logger());
				procSNPsToGenes(dir, filename, wiggleRoom, build, new Logger());
			}
		} catch (Exception e) {
			e.printStackTrace();
			if (filename.endsWith(".snps")) {
				try {
					new BufferedReader(new InputStreamReader(System.in)).readLine();
				} catch (Exception e2) {}					
			}			
		}
	}

}
