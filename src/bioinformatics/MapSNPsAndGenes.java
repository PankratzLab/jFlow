//-Xms2048M -Xmx2048M
package bioinformatics;

import java.io.*;
import common.*;

public class MapSNPsAndGenes {
//	public static final boolean AUTOMATICALLY_ADD_ONE = true;
//	public static final boolean AUTOMATICALLY_ADD_ONE = false;
	public static final String GENES = "/home/npankrat/NCBI/genes.xls";
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
//	public static final String DEFAULT_FILENAME = "ea_final_list.txt";
//	public static final String DEFAULT_FILENAME = "supp_t1_snps.txt";
//	public static final String DEFAULT_FILENAME = "table1_snps.txt";
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
			reader = new BufferedReader(new FileReader(GENES));
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
			log.reportError("Error: file \""+GENES+"\" not found");
			System.exit(1);
		} catch (IOException ioe) {
			log.reportError("Error reading file \""+GENES+"\"");
			System.exit(2);
		}

		return genes;
	}

	public static void procSNPsToGenes(String dir, String snps, int wiggleRoom, Logger log) {
		PrintWriter writer;
		String[] line;
		String[] data, markers, genes;
		int[][] markerPositions;

		ParseSNPlocations.lowMemParse(dir+snps, true, log);
		
		data = Array.toStringArray(HashVec.loadFileToVec(ext.rootOf(dir+snps, false)+"_positions.xln", false, false, false));

		markers = new String[data.length-1];
		markerPositions = new int[markers.length][2];
		for (int i = 0; i<markers.length; i++) {
			line = data[i+1].trim().split("[\\s]+");
			markers[i] = line[0];
			markerPositions[i][0] = Positions.chromosomeNumber(line[1]);
			markerPositions[i][1] = Integer.parseInt(line[2]);
		}

		genes = mapSNPsToGenes(markerPositions, wiggleRoom, log);

		try {
			writer = new PrintWriter(new FileWriter(dir+ext.rootOf(snps)+"_genes.xln"));
			writer.println("SNP\tChr\tPosition\tGene(s)\t"+UCSC_WINDOW+"\t<- dynamically linked basepair buffer in UCSC hyperlink");
			for (int i = 0; i<markers.length; i++) {
//				writer.println(markers[i]+"\t"+markerPositions[i][0]+"\t"+(AUTOMATICALLY_ADD_ONE?markerPositions[i][1]+1:markerPositions[i][1])+"\t"+(genes[i].equals("")?".":genes[i])+"\t=HYPERLINK(CONCATENATE(\"http://genome.ucsc.edu/cgi-bin/hgTracks?position=chr"+markerPositions[i][0]+":\","+markerPositions[i][1]+"-$E$1+1,\"-\","+markerPositions[i][1]+"+$E$1),CONCATENATE(\"chr"+markerPositions[i][0]+":\","+markerPositions[i][1]+"-$E$1+1,\"-\","+markerPositions[i][1]+"+$E$1))");
				writer.println(markers[i]+"\t"+markerPositions[i][0]+"\t"+markerPositions[i][1]+"\t"+(genes[i].equals("")?".":genes[i])+"\t=HYPERLINK(CONCATENATE(\"http://genome.ucsc.edu/cgi-bin/hgTracks?position=chr"+markerPositions[i][0]+":\","+markerPositions[i][1]+"-$E$1+1,\"-\","+markerPositions[i][1]+"+$E$1),CONCATENATE(\"chr"+markerPositions[i][0]+":\","+markerPositions[i][1]+"-$E$1+1,\"-\","+markerPositions[i][1]+"+$E$1))");
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

		String usage = "\\n"+
		"bioinformatics.MapSNPsAndGenes requires 0-1 arguments\n"+
		"   (1) directory (i.e. dir="+dir+" (default))\n"+
		"   (2) filename (i.e. file="+filename+" (default))\n"+
		"   (3) # bp up and down stream to count as an associated gene (i.e. win="+wiggleRoom+" (default))\n"+
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
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}

		try {
			if (filename.endsWith(".snps")) {
				temp = filename.substring(0, filename.length()-(".snps").length());
				if (temp.indexOf(".")>0) {
					try {
						wiggleRoom = Integer.parseInt(temp.substring(temp.lastIndexOf(".")+1));
					} catch (NumberFormatException nfe) {
						wiggleRoom = DEFAULT_WIGGLE_ROOM;
					}
				}				
				procSNPsToGenes("", filename, wiggleRoom, new Logger(filename+".log"));
			} else {
				procSNPsToGenes(dir, filename, wiggleRoom, new Logger());
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
