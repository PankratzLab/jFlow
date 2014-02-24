package bioinformatics;

import java.io.*;
import java.util.*;

import common.*;

public class SeattleSeq {
	public static final String[][] NEEDS = {{"Chr"}, {"MapInfo", "Position"}, {"MarkerName", "SNP"}, {"RefStrand"}};
	public static final String[][] RELEVANTS = {{"chromosome"}, {"position"}, {"sampleAlleles"}, {"accession"}, {"functionGVS"}, {"aminoAcids"}, {"geneList"}, {"inDBSNPOrNot"}, {"microRNAs"}};
//	"# inDBSNPOrNot", "", "position", "referenceBase", "sampleGenotype", "sampleAlleles", "allelesDBSNP", "accession", "functionGVS", "functionDBSNP", "rsID", "aminoAcids", "proteinPosition", "cDNAPosition", "polyPhen", "granthamScore", "scorePhastCons", "consScoreGERP", "chimpAllele", "CNV", "geneList", "AfricanHapMapFreq", "EuropeanHapMapFreq", "AsianHapMapFreq", "hasGenotypes", "dbSNPValidation", "repeatMasker", "tandemRepeat", "clinicalAssociation", "distanceToSplice", "microRNAs", "proteinSequence"
	public static final String[] ORDER = {"nonsense", "missense", "splice-5", "splice-3", "coding-synonymous", "coding-notMod3", "utr-5", "utr-3", "intron", "near-gene-5", "near-gene-3", "intergenic"};
	public static final String[] BAD = {"missense", "stop-gained", "stop-lost", "missense-near-splice", "splice-donor", "splice-acceptor"};
	public static final String[] NEUTRAL = {"intron-near-splice", "5-prime-UTR", "downstream-gene", "upstream-gene", "synonymous", "coding-synonymous", "intergenic", "non-coding-exon", "3-prime-UTR", "intron", "coding-notMod3"};
	
	public static void proc(String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		char temp;
		int[] indices;
		char[] alleles;
		
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+".input"));
			line = reader.readLine().trim().split("[\\s]+");
			indices = ext.indexFactors(NEEDS, line, false, true, true, true);
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (!line[indices[0]].startsWith("chr")) {
					line[indices[0]] = "chr"+line[indices[0]];
				}
				try {
					Integer.parseInt(line[indices[1]]);
				} catch (Exception e) {
					System.err.println("Error - invalid position ('"+line[indices[1]]+"')");
				}
				alleles = determineAlleles(line[indices[2]]);
				if (line[indices[3]].equals("-")) {
					temp = alleles[0];
					alleles[0] = alleles[1];
					alleles[1] = temp;
				} else if (!line[indices[3]].equals("+")) {
					System.err.println("Error - invalid strand ('"+line[indices[3]]+"')");
				}
				for (int i = 0; i < alleles.length; i++) {
					if (alleles[i] == 'I' || alleles[i] == 'D') {
						alleles[i] = 'A';
					}
				}
				writer.println(line[indices[0]]+"\t"+line[indices[1]]+"\t0\t"+alleles[0]+"\t"+alleles[1]);
			}
			writer.close();
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
		
	}
	
	public static char[] determineAlleles(String str) {
		String[] line;
		
		if (str.startsWith("[")) {
			str = str.substring(1);
		}
		
		if (str.startsWith("]")) {
			str = str.substring(0, str.length()-1);
		}
		
		line = str.split("/");
		
		return new char[] {line[0].charAt(0), line[1].charAt(0)};		
	}
	
	// the indels are not captured by SeattleSeq
	// there are triallelic markers on the exome chip, hence then Hashtable, but that doesn't necessarily caputre all, still missing a little over a hundred
	public static void summarize(String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, trav;
		String temp;
		Hashtable<String, Vector<String[]>> hash = new Hashtable<String, Vector<String[]>>();
		Vector<String[]> v;
		int[] indices;
		Logger log;
		String prev;
		boolean done;
		int type, worstType, worst;
		String[] keys;
		Hashtable<String,String> hashFreq;
		double freq;
		
		hashFreq = HashVec.loadFileToHashString("D:/home/npankrat/NCBI/ESP_exome_chip/EVS/freqInfo_proc.dat", new int[] {0}, new int[] {3}, false, "", true, false, false);
		
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_summary.out"));
			writer.println(Array.toStr(Matrix.extractColumn(RELEVANTS, 0)));
			log = new Logger(ext.rootOf(filename, false)+".log");
			temp = reader.readLine().trim();
			if (temp.startsWith("#"));
			temp = temp.substring(1).trim();
			line = temp.split("[\\s]+");
			indices = ext.indexFactors(RELEVANTS, line, false, true, true, log, true);
//			public static final String[][] RELEVANTS = {{"chromosome"}, {"position"}, {"accession"}, {"functionGVS"}, {"aminoAcids"}, {"geneList"}, {"inDBSNPOrNot"}, {"microRNAs"}};
			
			prev = "";
			
			done = false;
			while (!done) {
				if (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
				} else {
					done = true;
				}
				if (hash.size() > 0 && !(line[indices[0]]+"\t"+line[indices[1]]).equals(prev) || done) {
					keys = HashVec.getKeys(hash);
					for (int key = 0; key < keys.length; key++) {
						v = hash.get(keys[key]);
						worst = -1;
						worstType = ORDER.length;
						type = -1;
						for (int i = 0; i < v.size(); i++) {
							trav = v.elementAt(i);
							type = ext.indexOfStr(trav[indices[4]], ORDER);
							if (type == -1) {
								System.out.println(Array.toStr(trav));
								System.exit(1);
							}
							if (type < worstType) {
								worstType = type;
								worst = i;
							}						
						}
						trav = v.elementAt(worst);
						if (worstType == 0) {
							if (trav[indices[5]].startsWith("stop")) {
								trav[indices[4]] = "stoploss";
							} else {
								trav[indices[4]] = "stopgain";
							}
						}
						for (int i = 0; i < indices.length; i++) {
							writer.print((i==0?"":"\t")+trav[indices[i]]);
						}
						if (hashFreq.containsKey(line[indices[0]].toLowerCase()+"^"+line[indices[1]])) {
							freq = Double.parseDouble(hashFreq.get(line[indices[0]].toLowerCase()+"^"+line[indices[1]]));
							writer.print("\t"+freq+"\t"+(freq<0.01?1:0)+"\t"+(freq<0.05?1:0));
						} else {
							writer.print("\t.\t.\t.");
						}
						writer.println();
					}
					hash = new Hashtable<String, Vector<String[]>>();
				}
				HashVec.addToHashArrayVec(hash, line[indices[2]], line);
				prev = line[indices[0]]+"\t"+line[indices[1]];
			}
			writer.close();
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
	}
	
	public static void parseFreq(String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		double a1, a2, temp;
		
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_proc.dat"));
			writer.println(reader.readLine());
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				for (int i = 1; i < line.length; i++) {
					a1 = ext.parseDoubleArg(line[i].split("/")[0]);
					a2 = ext.parseDoubleArg(line[i].split("/")[1]);
					if (a2<a1) {
						temp = a1;
						a1 = a2;
						a2 = temp;
					}
					line[i] = a1/(a1+a2)+"";
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

	public static void geneCounts(String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, trav;
		String temp;
		Hashtable<String, Vector<String[]>> hash = new Hashtable<String, Vector<String[]>>();
		Vector<String[]> v;
		int[] indices;
		Logger log;
		String prev;
		boolean done;
		int type, worstType, worst;
		String[] keys;
		Hashtable<String,String> hashFreq;
		double freq;
		
		hashFreq = HashVec.loadFileToHashString("D:/home/npankrat/NCBI/ESP_exome_chip/EVS/freqInfo_proc.dat", new int[] {0}, new int[] {3}, false, "", true, false, false);
		
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_summary.out"));
			writer.println(Array.toStr(Matrix.extractColumn(RELEVANTS, 0)));
			log = new Logger(ext.rootOf(filename, false)+".log");
			temp = reader.readLine().trim();
			if (temp.startsWith("#"));
			temp = temp.substring(1).trim();
			line = temp.split("[\\s]+");
			indices = ext.indexFactors(RELEVANTS, line, false, true, true, log, true);
//			public static final String[][] RELEVANTS = {{"chromosome"}, {"position"}, {"accession"}, {"functionGVS"}, {"aminoAcids"}, {"geneList"}, {"inDBSNPOrNot"}, {"microRNAs"}};
			
			prev = "";
			
			done = false;
			while (!done) {
				if (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
				} else {
					done = true;
				}
				if (hash.size() > 0 && !(line[indices[0]]+"\t"+line[indices[1]]).equals(prev) || done) {
					keys = HashVec.getKeys(hash);
					for (int key = 0; key < keys.length; key++) {
						v = hash.get(keys[key]);
						worst = -1;
						worstType = ORDER.length;
						type = -1;
						for (int i = 0; i < v.size(); i++) {
							trav = v.elementAt(i);
							type = ext.indexOfStr(trav[indices[4]], ORDER);
							if (type == -1) {
								System.out.println(Array.toStr(trav));
								System.exit(1);
							}
							if (type < worstType) {
								worstType = type;
								worst = i;
							}						
						}
						trav = v.elementAt(worst);
						if (worstType == 0) {
							if (trav[indices[5]].startsWith("stop")) {
								trav[indices[4]] = "stoploss";
							} else {
								trav[indices[4]] = "stopgain";
							}
						}
						for (int i = 0; i < indices.length; i++) {
							writer.print((i==0?"":"\t")+trav[indices[i]]);
						}
						if (hashFreq.containsKey(line[indices[0]].toLowerCase()+"^"+line[indices[1]])) {
							freq = Double.parseDouble(hashFreq.get(line[indices[0]].toLowerCase()+"^"+line[indices[1]]));
							writer.print("\t"+freq+"\t"+(freq<0.01?1:0)+"\t"+(freq<0.05?1:0));
						} else {
							writer.print("\t.\t.\t.");
						}
						writer.println();
					}
					hash = new Hashtable<String, Vector<String[]>>();
				}
				HashVec.addToHashArrayVec(hash, line[indices[2]], line);
				prev = line[indices[0]]+"\t"+line[indices[1]];
			}
			writer.close();
			reader.close();
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
		String filename = "SeattleSeq.dat";

		String usage = "\n" + 
		"bioinformatics.SeattleSeq requires 0-1 arguments\n" + 
				"   (1) filename (i.e. file=" + filename + " (default))\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
//		filename = "D:/home/npankrat/NCBI/ESP_exome_chip/seattleInput.txt";
		filename = "D:/home/npankrat/NCBI/ESP_exome_chip/SeattleSeqAnnotation131.seattleInput.input.240042724786.txt";
//		filename = "D:/home/npankrat/NCBI/ESP_exome_chip/FGG.txt";
//		filename = "D:/home/npankrat/NCBI/ESP_exome_chip/F11.txt";
		try {
//			proc(filename);
//			summarize(filename);
			geneCounts(filename);
//			parseFreq("D:/home/npankrat/NCBI/ESP_exome_chip/EVS/freqInfo.out");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static Hashtable<String, String[]> loadAllAnnotationInDir(String directory, Logger log) {
		BufferedReader reader;
		String[] files, line;
		Hashtable<String, String[]> hash;
		String markerName;
		String function;

		hash = new Hashtable<String, String[]>();
		
		if (directory == null) {
			log.reportError("The SeattleSeq annotation directory was null; returning an empty hashtable");
		} else if (!Files.exists(directory) || !Files.isDirectory(directory)) {
			log.reportError("Error - SeattleSeq annotation directory directory not found: "+directory);
			log.reportError("        returning an empty hashtable");
		} else {
			files = Files.list(directory, "SeattleSeqAnnotation", ".txt.gz", false, false);
			log.report("Found "+files.length+" file(s) with a .SeattleSeq extension to include");
			for (int i = 0; i < files.length; i++) {
				try {
					reader = Files.getAppropriateReader(directory+files[i]);
					while (reader.ready()) {
						line = reader.readLine().trim().split("\t", -1);
						if (line.length > 1) {
							markerName = "chr"+line[1]+":"+line[2]+"_"+line[3]+"_"+line[4];
							if (!hash.containsKey(markerName) || hash.get(markerName) == null) {
								if (ext.indexOfStr(line[8], BAD) >= 0) {
									function = line[8];
									if (!line[11].equals("none")) {
										function += " "+ext.replaceAllWith(line[11], ",", line[12].substring(0, line[12].indexOf("/")));
									}
									function += "\t"+line[0];
									hash.put(markerName, new String[] {function});
								} else {
									hash.put(markerName, new String[0]);
								}
							}
						}
						
						// TODO
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \"" + directory+files[i]
							+ "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + directory+files[i] + "\"");
					System.exit(2);
				}
			}
		}		

		return hash;
	}
}
