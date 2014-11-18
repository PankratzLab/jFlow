// expecting a HapMart/BioMart export that includes position, marker_id, reference_allele and of course CEU genotype 
// pedinfo2sample_CEU.txt is available at http://www.hapmap.org/downloads/samples_individuals/
package bioinformatics;

import java.io.*;
import java.util.*;

import common.*;
import filesys.*;

public class HapMapParser {
	public static final boolean USE_NUCLEOTIDES = true;
	public static final boolean INCLUDE_MONOMORPHIC = true;
	
	public static final String[][] TARGETS_WITH_ALTS = {{"genotype"}, {"position"}, {"reference allele"}, {"marker id"}, {"chromosome"}};
	
	public static void parse(String dir, String filename, String famstruct, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, trans;
		String temp, trav;
		Hashtable<String,String[]> hash = new Hashtable<String,String[]>();
		Vector<String> v = new Vector<String>();
		int index, ones, twos;
		String root = ext.rootOf(filename);
		int[] indices;
		String refAllele, position, chr, markerName;
		String[] indIDs = null;

		try {
			reader = Files.getReader(filename, dir);
			temp = reader.readLine();
			indices = ext.indexFactors(TARGETS_WITH_ALTS, temp.trim().split("\\t", -1), false, false, true, log, true);

			writer = new PrintWriter(new FileWriter((new File(dir).exists()?dir:"")+root+".map"));
			trans = null;
			while (reader.ready()) {
				line = reader.readLine().split("\\t", -1);
				refAllele = line[indices[2]];
				position = line[indices[1]];
				chr = line[indices[4]].substring(3);
				markerName = line[indices[3]];
				
				line = line[indices[0]].split("[\\s]+");
				trans = new String[line.length];
				ones = twos = 0;
				for (int i = 0; i<line.length; i++) {
					trans[i] = "";
					for (int j = 0; j<2; j++) {
						if (line[i].charAt(j)=='N') {
							trans[i] += "0";
						} else if (line[i].charAt(j)==refAllele.charAt(0)) {
							if (USE_NUCLEOTIDES) {
								trans[i] += line[i].charAt(j);
							} else {
								trans[i] += "2";
							}
							twos++;
						} else {
							if (USE_NUCLEOTIDES) {
								trans[i] += line[i].charAt(j);
							} else {
								trans[i] += "1";
							}
							ones++;
						}
						trans[i] += j==0?"\t":"";
					}
				}
				if ((ones>0&&twos>0) || INCLUDE_MONOMORPHIC) {
					writer.println(chr+"\t"+markerName+"\t0\t"+position);
					v.add(markerName);
					hash.put(markerName, trans);
				}
			}
			reader.close();
			writer.close();
			indIDs = Array.stringArraySequence(trans.length, "");
		} catch (IOException ioe) {
			log.reportError("Error reading file \""+filename+"\"");
			log.reportException(ioe);
			return;
		}

		try {
			reader = Files.getReader(famstruct, new String[] {dir, "/home/npankrat/NCBI/", "C:\\Documents and Settings\\npankrat\\My Documents\\jProjects\\park\\runtime\\"});
			writer = new PrintWriter(new FileWriter((new File(dir).exists()?dir:"")+root+".pre"));
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				for (int i = 0; i<5; i++) {
					writer.print(line[i]+"\t");
				}
				writer.print("1");
				index = -2;
				trav = line[6];
				if (famstruct.contains("pedinfo2sample_")) {
					line = line[6].split(":");
					if (line.length!=6) {
						log.reportError("Error - different format than expected for pedinfo2sample_***.txt file (do not alter from what's posted)");
						reader.close();
						writer.close();
						return;
					}
					trav = line[4];
				}
				for (int i = 0; i<indIDs.length; i++) {
					if (indIDs[i].equals(trav)) {
						index = i;
					}
				}
				if (index==-2) {
					log.reportError("Error - Could not find sample "+trav+" from the pedigree file in the genotype file");
				}
				for (int i = 0; i<v.size(); i++) {
					line = hash.get(v.elementAt(i));
					writer.print("\t"+line[index]);
				}
				writer.println();
			}
			reader.close();
			writer.close();
		} catch (IOException ioe) {
			log.reportError("Error reading file \""+famstruct+"\"");
			log.reportException(ioe);
			return;
		}
		
		generateHaploviewBatch(dir, root, true, log);		
	}

	public static void generateHaploviewBatch(String dir, String root, boolean preNotPed, Logger log) {
		PrintWriter writer;

		new SnpMarkerSet(dir+root+".map").writeToFile(dir+root+".info", SnpMarkerSet.HAPLOVIEW_INFO_FORMAT);
		try {
			writer = new PrintWriter(new FileWriter((new File(dir).exists()?dir:"")+root+".bat"));
			writer.println("java -jar /home/npankrat/Haploview.jar -pedfile "+root+"."+(preNotPed?"pre":"ped")+" -info "+root+".info");
			writer.close();
        } catch (Exception e) {
        	log.reportError("Error writing batch file");
        	log.reportException(e);
        }
	}

	public static void plinkMapToHaploviewInfo(String from, String to) {
		new SnpMarkerSet(from).writeToFile(to, SnpMarkerSet.HAPLOVIEW_INFO_FORMAT);
	}
	
	public static void splitBedByChromosome(String root) {
		for (int chr = 1; chr<=23; chr++) {
			System.out.print(".");
            CmdLine.run("plink --noweb --bfile "+root+" --chr "+chr+" --make-bed --out "+root+".chr"+chr, "./");
        }
		System.out.println();
	}

	public static void splitPedByChromosome(String root) {
		for (int chr = 1; chr<=23; chr++) {
			System.out.print(".");
            CmdLine.run("plink --noweb --file "+root+" --chr "+chr+" --recode --out "+root+".chr"+chr, "./");
        }
		System.out.println();
	}

	public static void fixAffectionStatusInBed(String root) {
        FamilyStructure struct;
        long timestamp;
        
		for (int chr = 1; chr<=23; chr++) {
			timestamp = new File(root+".chr"+chr+".fam").lastModified();
			struct = new FamilyStructure(root+".chr"+chr+".fam");
			struct.setAffections(Array.byteArray(struct.getAffections().length, (byte)1));
			struct.writeToFile(root+".chr"+chr+".fam", false);
			new File(root+".chr"+chr+".fam").setLastModified(timestamp);
			System.out.print(".");
        }
		System.out.println();
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String dir = "";
//		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\SNCA\\";
//		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\followup\\APOJ\\";
		// String filename = "2qCEU-GENOS.tsv";
		// String filename = "39rGV1lKEW5.tsv";
		// String filename = "kept39rGV1lKEW2.tsv";
		// String filename = "SNCA_50kb.tsv";
		// String filename = "SNCA_50kb_43tags.tsv";
		// String filename = "SNCA_50kb_21tags.tsv";
		// String filename = "mapt.tsv";
		// String filename = "CTNNA2.tsv";
		// String filename = "deleteriousTau.tsv";
		// String filename = "GAKimputed.tsv";
		// String filename = "maptUTRsigs.tsv";
		// String filename = "gakExonic.tsv";
		// String filename = "maptUTR2.tsv";
//		String filename = "LAMP1.tsv";
//		String filename = "APOJ.tsv";
		String filename = null;
		String map = "";		

//		String famstruct = "famstruct_CEU.dat";
//		String famstruct = "/home/npankrat/NCBI/pedinfo2sample_CEU.txt";
//		String famstruct = "/home/npankrat/NCBI/fakeIDs_for_CEU.txt";
		String famstruct = "/home/npankrat/NCBI/CEU_structUnrelated.dat";
		
		String bed = "";
		String ped = "";
		String fix = "";
		
		Logger log;

		String usage = "\n"+
		"bioinformatics.HapMapParser requires 0-1 arguments\n"+
		"   (1) filename (i.e. file=MAPT.tsv (not the default)\n"+
		"   (2) famstruct (i.e. struct="+famstruct+" (default)\n"+
		"   (3) split bed by chromosome (i.e. bed=plink (not the default)\n"+
		"   (4) split ped by chromosome (i.e. ped=plink (not the default)\n"+
		"   (5) fix affection status in fam files (i.e. fix=plink (not the default)\n"+
		"   (6) generate Haploview batch file from PLINK files (i.e. map=plink.map (not the default)\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("struct=")) {
				famstruct = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("bed=")) {
				bed = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("ped=")) {
				ped = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("fix=")) {
				fix = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("map=")) {
				map = args[i].split("=")[1];
				numArgs--;
			}
		}
		
		if (filename == null && map.equals("")) {
			System.err.println("Error - need to pass a filename as an argument (i.e. file=MAPT.tsv)");
			ext.waitForResponse();
			return;
		}
		
		log = new Logger(ext.rootOf(filename==null?map:filename, false)+"_hapmap_parser.log");
				
//		if (filename.startsWith("C:")) {
//			dir = "";
//		}
		
		for (int i = 0; i<args.length; i++) {
			log.report((i+1)+") "+args[i]);
		}
		if (numArgs!=0) {
			log.reportError(usage);
			System.exit(1);
		}
		try {
			if (!map.equals("")) {
				generateHaploviewBatch("", ext.removeDirectoryInfo(map.substring(0, map.lastIndexOf("."))), false, log);
			} else if (!bed.equals("")) {
				splitBedByChromosome(bed);
			} else if (!ped.equals("")) {
				splitPedByChromosome(ped);
			} else if (!fix.equals("")) {
				fixAffectionStatusInBed(fix);
			} else {
				parse(dir, filename, famstruct, log);
			}
		} catch (Exception e) {
			e.printStackTrace();
			ext.waitForResponse();
		}
	}
}
