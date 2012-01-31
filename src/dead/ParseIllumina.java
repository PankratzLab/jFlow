package dead;

import java.io.*;
import java.util.*;

import common.*;

public class ParseIllumina {
	// public static final String WINDOWS_DIRECTORY = "C:\\Documents and Settings\\npankrat\\My Documents\\jProjects\\park\\runtime";
	// public static final String WINDOWS_DIRECTORY = "C:\\Documents and Settings\\npankrat\\Desktop";
	// public static final String WINDOWS_DIRECTORY = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\00src\\";
	public static final String WINDOWS_DIRECTORY = "C:\\Documents and Settings\\npankrat\\My Documents\\jProjects\\park\\runtime\\adni\\";
	// public static final String LINUX_DIRECTORY = "/archive/parkinsons/gwas/Genotype and intensity data files";
	// public static final String LINUX_DIRECTORY = "/work/parkinsons/gwas/Genotype and intensity data files/Final_Reports";
	// public static final String LINUX_DIRECTORY = "/home/npankrat/cnv/oldClusters/raw";
	// public static final String LINUX_DIRECTORY = "/home/npankrat/gwas/reloadedIndividuals_n4";
	// public static final String LINUX_DIRECTORY = "/home/npankrat/gwas/peds/newnewcalls/source";
	// public static final String LINUX_DIRECTORY = "/archive/parkinsons/gwas/Genotype and intensity data files/mosaics/deleted_X";
	// public static final String LINUX_DIRECTORY = "/work/adnigwas/ADNI 676 Subjects";
	// public static final String LINUX_DIRECTORY = "/work/osteo/GWAS_data/Genotype_and_Intensity_Data_Files/Econs_610Quad_genotype_files_010609";
	public static final String LINUX_DIRECTORY = "/net/collabs/genetics_anals/LONI_subjs";
	public static final String PED_DIRECTORY = "peds/";
	public static final String CNV_DIRECTORY = "cnvs/";
	public static final String MARKERS_IN_INDEX_ORDER = "markersInIndexOrder.dat";
	public static final String PED_FILE = "pedfile.dat";
	public static final String CNV_MARKER_POSITIONS = "cnvMarkerPositions.xln"; // only include those markers you want in the CNV files
	public static final String TARGET_MARKERS = "markerListWithIndices.dat";
	// public static final double GC_SCORE_CUTOFF = 0.25;
	public static final double GC_SCORE_CUTOFF = 0.15;
	public static final double[] POSSIBLE_CUTOFFS = {0.15, 0.25};
	public static final boolean MAKE_CNV_FILES = true;
	public static final String[] FIELDS = {"SNP Name", "Sample ID", "Sample Name", "GC Score", "Allele1 - Forward", "Allele2 - Forward", "Allele1 - AB", "Allele2 - AB", "B Allele Freq", "Log R Ratio"};
	public static final String[] INDEXING = {"SNP Name", "Chr", "Position", "SNP Index"};

	public static void parseToIndividualPedFiles() {
		BufferedReader reader;
		PrintWriter writer, missWriter, cnvWriter = null;
		String[] line, snpNames = null, header;
		int count, version, index, indexIndex;
		String id = "", trav, temp;
		boolean pdgwas;
		int[] indices;
		int[][] countMissAcrossMarkers;
		int[] countMissAcrossIndividuals;
		double score;
		char[][] data;
		Hashtable<String,String> cnvPoslar = new Hashtable<String,String>();

		if (new File(WINDOWS_DIRECTORY).exists()) {
			trav = WINDOWS_DIRECTORY;
		} else if (new File(LINUX_DIRECTORY).exists()) {
			trav = LINUX_DIRECTORY;
		} else {
			trav = "";
			System.err.println("Error - could not resolve directory to parse (none of the following worked)");
			System.err.println(WINDOWS_DIRECTORY);
			System.err.println(LINUX_DIRECTORY);
			System.exit(1);
		}

		File[] files = new File(trav).listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				// return filename.startsWith("Myers") &&
				// filename.endsWith(".csv");
				return filename.endsWith(".csv");
			}
		});

		System.out.println(ext.getTime());
		System.out.println("Found "+files.length+" files to parse");

		snpNames = Array.toStringArray(HashVec.loadFileToVec(MARKERS_IN_INDEX_ORDER, false, true, false));
		System.out.println("There are "+snpNames.length+" markers being read in");

		countMissAcrossMarkers = new int[snpNames.length][POSSIBLE_CUTOFFS.length+2];

		new File(PED_DIRECTORY).mkdirs();
		if (MAKE_CNV_FILES) {
			new File(CNV_DIRECTORY).mkdirs();
		}
		cnvPoslar = HashVec.loadFileToHashString(CNV_MARKER_POSITIONS, 0, new int[] {1, 2}, "\t", false);

		System.out.println(ext.getTime());
		try {
			missWriter = new PrintWriter(new FileWriter("missingByIndividual.xln"));
			missWriter.println("Marker\t"+Array.toStr(POSSIBLE_CUTOFFS)+"\tNotCalled\tNotPresent");
			for (int i = 0; i<files.length; i++) {
				pdgwas = files[i].getName().contains("Myers");
				try {
					reader = new BufferedReader(new FileReader(files[i]));
					do {
						temp = reader.readLine();
					} while (reader.ready()&&!temp.contains("SNP Name")&&!temp.contains("Sample ID"));

					if (!reader.ready()) {
						System.err.println("Error - went through enitre file without finding a line containing both 'SNP Name' and 'Sample ID'");
						System.exit(1);
					}
					header = temp.trim().split(",");
					indices = ext.indexFactors(FIELDS, header, false, true);
					indexIndex = ext.indexOfStr("SNP Index", header);

					count = 0;
					countMissAcrossIndividuals = new int[POSSIBLE_CUTOFFS.length+2];
					data = new char[snpNames.length][2];
					version = 0;
					while (reader.ready()) {
						line = reader.readLine().split(",");
						if (pdgwas) {
							trav = line[indices[ext.indexOfStr("Sample Name", FIELDS)]].substring(0, line[indices[ext.indexOfStr("Sample Name", FIELDS)]].indexOf("@"));
						} else {
							trav = line[indices[ext.indexOfStr("Sample ID", FIELDS)]];
						}
						if (count==0) {
							id = trav;
							while (new File(PED_DIRECTORY+id+(version==0?"":"__"+version)+".ped").exists()) {
								version++;
							}
							if (MAKE_CNV_FILES) {
								cnvWriter = new PrintWriter(new FileWriter(CNV_DIRECTORY+id+(version==0?"":"."+version)));
								cnvWriter.println("Name\tChr\tPosition\t"+id+".GType\t"+id+".Log R Ratio\t"+id+".B Allele Freq");
							}
						} else if (!trav.equals(id)) {
							System.err.println("Found "+trav+" -- expecting "+id+" in file "+files[i].getName());
							System.exit(1);
						}

						if (indexIndex==-1) {
							index = count;
							if (!snpNames[count].equals(line[indices[ext.indexOfStr("SNP Name", FIELDS)]])) {
								System.err.println("Found "+line[indices[ext.indexOfStr("SNP Name", FIELDS)]]+" -- expecting "+snpNames[count]+" in file "+files[i].getName());
								System.exit(1);
							}
						} else {
							index = Integer.parseInt(line[indexIndex])-1;
							if (!line[indices[ext.indexOfStr("SNP Name", FIELDS)]].equals(snpNames[index])) {
								System.err.println("Error - two SNPs with the same index ("+index+") in file "+files[i].getName()+": "+snpNames[index]+" in "+MARKERS_IN_INDEX_ORDER+" and "+line[indices[ext.indexOfStr("SNP Name", FIELDS)]]+" here");
								System.exit(1);
							}
						}

						score = Double.parseDouble(line[indices[ext.indexOfStr("GC Score", FIELDS)]]);
						for (int j = 0; j<POSSIBLE_CUTOFFS.length; j++) {
							if (score<POSSIBLE_CUTOFFS[j]) {
								countMissAcrossMarkers[index][j]++;
								countMissAcrossIndividuals[j]++;
							}
						}
						if (score<GC_SCORE_CUTOFF) {
							line[indices[ext.indexOfStr("Allele1 - Forward", FIELDS)]] = "0";
							line[indices[ext.indexOfStr("Allele2 - Forward", FIELDS)]] = "0";
							countMissAcrossMarkers[index][POSSIBLE_CUTOFFS.length]++;
							countMissAcrossIndividuals[POSSIBLE_CUTOFFS.length]++;
						}

						data[index][0] = line[indices[ext.indexOfStr("Allele1 - Forward", FIELDS)]].charAt(0);
						data[index][1] = line[indices[ext.indexOfStr("Allele2 - Forward", FIELDS)]].charAt(0);

						if (MAKE_CNV_FILES&&cnvPoslar.containsKey(snpNames[index])&&!line[indices[ext.indexOfStr("Log R Ratio", FIELDS)]].equals("NaN")) {
							cnvWriter.print(snpNames[index]+"\t"+cnvPoslar.get(snpNames[index]));
							temp = line[indices[ext.indexOfStr("Allele1 - AB", FIELDS)]]+line[indices[ext.indexOfStr("Allele2 - AB", FIELDS)]];
							if (temp.contains("-")) {
								temp = "NC";
							}
							cnvWriter.print("\t"+temp);
							cnvWriter.print("\t"+line[indices[ext.indexOfStr("Log R Ratio", FIELDS)]]);
							cnvWriter.print("\t"+line[indices[ext.indexOfStr("B Allele Freq", FIELDS)]]);
							cnvWriter.println();
						}
						count++;
					}
					reader.close();
					if (MAKE_CNV_FILES) {
						cnvWriter.close();
					}

					writer = new PrintWriter(new FileWriter(PED_DIRECTORY+id+(version==0?"":"__"+version)+".ped"));
					writer.print(id);
					for (int j = 0; j<data.length; j++) {
						if (data[j][0]=='-'||data[j][1]=='-') {
							writer.print("\t0\t0");
							countMissAcrossIndividuals[POSSIBLE_CUTOFFS.length]++;
							countMissAcrossMarkers[j][POSSIBLE_CUTOFFS.length]++;
						} else if ((int)data[j][0]==0) {
							writer.print("\t0\t0");
							countMissAcrossIndividuals[POSSIBLE_CUTOFFS.length+1]++;
							countMissAcrossMarkers[j][POSSIBLE_CUTOFFS.length+1]++;
						} else {
							writer.print("\t"+data[j][0]+"\t"+data[j][1]);
						}
					}
					writer.println();
					writer.close();

					missWriter.println(trav+"\t"+Array.toStr(countMissAcrossIndividuals));
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""+files[i].getName()+"\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error processing file \""+files[i].getName()+"\"");
					writer = new PrintWriter(new FileWriter("ERRORS FOUND WHEN PARSING ILLUMINA DATA.txt", true));
					writer.println("Error processing file \""+files[i].getName()+"\"");
					writer.close();
				}
			}
			missWriter.close();

			missWriter = new PrintWriter(new FileWriter("missingByMarker.xln"));
			missWriter.println("SNP\t"+Array.toStr(POSSIBLE_CUTOFFS)+"\tNotCalled\tNotPresent");
			for (int i = 0; i<snpNames.length; i++) {
				missWriter.println(snpNames[i]+"\t"+Array.toStr(countMissAcrossMarkers[i]));
			}
			missWriter.close();

			System.out.println(ext.getTime());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void createPlink(String pedigreeFile, String targetMarkers) {
		BufferedReader reader, pedreader;
		PrintWriter writer;
		Vector<String> targets = new Vector<String>();
		IntVector indexVec = new IntVector();
		int count;
		String[] line;
		int[] keys, indices;
		String trav;

		System.out.println(ext.getTime());
		try {
			reader = new BufferedReader(new FileReader(targetMarkers));
			writer = new PrintWriter(new FileWriter("gwas.map"));
			ext.checkHeader(reader.readLine().split("[\\s]+"), new String[] {"Name", "GenomeBuild", "Chr", "MapInfo", "Index"}, false);
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				targets.add(line[0]);
				indexVec.add(Integer.parseInt(line[4]));
				writer.println((line[2].equals("X")?"23":(line[2].equals("Y")?"24":(line[2].equals("XY")?"25":(line[2].equals("MT")?"26":line[2]))))+"\t"+line[0]+"\t"+(line.length>5?line[5]:"0")+"\t"+line[3]);
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+targetMarkers+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+targetMarkers+"\"");
			System.exit(2);
		}
		indices = indexVec.toArray();
		keys = Sort.quicksort(indices);
		try {
			reader = new BufferedReader(new FileReader(MARKERS_IN_INDEX_ORDER));
			count = 0;
			for (int i = 0; i<indices.length; i++) {
				while (count<indices[keys[i]]) {
					reader.readLine();
					count++;
				}
				trav = reader.readLine().split("[\\s]+")[0];
				if (!trav.equals(targets.elementAt(keys[i]))) {
					System.err.println("Error - marker "+targets.elementAt(keys[i])+" was not at the index you said it would be (found "+trav+"). For shame...");
				}
				count++;
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+MARKERS_IN_INDEX_ORDER+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+MARKERS_IN_INDEX_ORDER+"\"");
			System.exit(2);
		}
		try {
			reader = new BufferedReader(new FileReader(pedigreeFile));
			writer = new PrintWriter(new FileWriter("gwas.ped"));
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				writer.print(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]);
				if (line[6].equals(".")) {
					for (int i = 0; i<indices.length; i++) {
						writer.print("\t0\t0");
					}
				} else if (!new File(PED_DIRECTORY+line[6]+".ped").exists()) {
					System.err.println("Error - the DNA# "+line[6]+" was listed in the pedigree file but the following file was nowhere in sight: "+PED_DIRECTORY+line[6]+".ped");
					for (int i = 0; i<indices.length; i++) {
						writer.print("\t0\t0");
					}
				} else {
					try {
						pedreader = new BufferedReader(new FileReader(PED_DIRECTORY+line[6]+".ped"));
						line = pedreader.readLine().split("[\\s]+");
						for (int i = 0; i<indices.length; i++) {
							writer.print("\t"+line[1+indices[i]*2+0]+"\t"+line[1+indices[i]*2+1]);
						}
						pedreader.close();
					} catch (FileNotFoundException fnfe) {
						System.err.println("Error: file \""+line[6]+".ped"+"\" not found in current directory");
						System.exit(1);
					} catch (IOException ioe) {
						System.err.println("Error reading file \""+line[6]+".ped"+"\"");
						System.exit(2);
					}
				}
				writer.println();

			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+pedigreeFile+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+pedigreeFile+"\"");
			System.exit(2);
		}

		System.out.println(ext.getTime());

	}

	public static void getListOfSamples() {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String trav, temp;
		int[] indices;

		if (new File(WINDOWS_DIRECTORY).exists()) {
			trav = WINDOWS_DIRECTORY;
		} else if (new File(LINUX_DIRECTORY).exists()) {
			trav = LINUX_DIRECTORY;
		} else {
			trav = "";
			System.err.println("Error - could not resolve directory to parse (none of the following worked)");
			System.err.println(WINDOWS_DIRECTORY);
			System.err.println(LINUX_DIRECTORY);
			System.exit(1);
		}

		File[] files = new File(trav).listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(".csv");
			}
		});

		System.out.println(ext.getTime());
		System.out.println("Found "+files.length+" files");
		try {
			writer = new PrintWriter(new FileWriter("individuals.xln"));
			for (int i = 0; i<files.length; i++) {
				try {
					reader = new BufferedReader(new FileReader(files[i]));
					do {
						temp = reader.readLine();
					} while (reader.ready()&&!temp.contains("SNP Name")&&!temp.contains("Sample ID"));

					if (!reader.ready()) {
						System.err.println("Error - went through enitre file without finding a line containing both 'SNP Name' and 'Sample ID'");
						System.exit(1);
					}
					indices = ext.indexFactors(FIELDS, temp.trim().split(","), false, true);

					line = reader.readLine().split(",");
					writer.println(files[i].getName()+"\t"+line[indices[ext.indexOfStr("Sample ID", FIELDS)]]+"\t"+line[indices[ext.indexOfStr("Sample Name", FIELDS)]]);
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""+files[i].getName()+"\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \""+files[i].getName()+"\"");
					System.exit(2);
				}
			}
			System.out.println(ext.getTime());
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void indexMarkers(String filename) {
		// originally written with the following file in mind
		// String filename = "C:\\Documents and Settings\\npankrat\\My
		// Documents\\tWork\\COGA_CIDR_GWAS\\SNP_Table_042808.csv";
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String trav, temp;
		byte[] chrs;
		int[] indices, poslar, keys;
		Vector<String[]> v = new Vector<String[]>();
		File[] files = null;
		int maxIndex;

		if (filename==null) {
			if (new File(WINDOWS_DIRECTORY).exists()) {
				trav = WINDOWS_DIRECTORY;
			} else if (new File(LINUX_DIRECTORY).exists()) {
				trav = LINUX_DIRECTORY;
			} else {
				trav = "";
				System.err.println("Error - could not resolve directory to parse (none of the following worked)");
				System.err.println(WINDOWS_DIRECTORY);
				System.err.println(LINUX_DIRECTORY);
				System.exit(1);
			}

			files = new File(trav).listFiles(new FilenameFilter() {
				public boolean accept(File file, String filename) {
					return filename.endsWith(".csv");
				}
			});
			if (files.length==0) {
				System.err.println("Error - no .csv files found in directory: "+trav);
				System.exit(1);
			}
			System.out.println(ext.getTime()+"\tparsing marker indices from '"+files[0].getName()+"'");
		} else {
			System.out.println(ext.getTime()+"\tparsing marker indices from '"+filename+"'");
			System.out.println("Though I must warn you that I haven't yet coded to filter out intensity only markers from "+TARGET_MARKERS);
		}
		try {
			reader = new BufferedReader(filename==null?new FileReader(files[0]):new FileReader(filename));

			do {
				temp = reader.readLine();
			} while (reader.ready()&&!temp.contains("SNP Name")&&!temp.contains("Sample ID"));

			if (!reader.ready()) {
				System.err.println("Error - went through enitre file without finding a line containing both 'SNP Name' and 'Sample ID'");
				System.exit(1);
			}
			indices = ext.indexFactors(INDEXING, temp.trim().split(","), false, true);
			maxIndex = 0;
			while (reader.ready()) {
				line = reader.readLine().split(",");
				v.add(new String[] {line[indices[ext.indexOfStr("SNP Name", INDEXING)]], Positions.chromosomeNumber(line[indices[ext.indexOfStr("Chr", INDEXING)]])+"", line[indices[ext.indexOfStr("Position", INDEXING)]], (Integer.parseInt(line[indices[ext.indexOfStr("SNP Index", INDEXING)]])-1)+""});
				if (Integer.parseInt(line[indices[ext.indexOfStr("SNP Index", INDEXING)]])>maxIndex) {
					maxIndex = Integer.parseInt(line[indices[ext.indexOfStr("SNP Index", INDEXING)]]);
				}
			}
			reader.close();

			if (maxIndex!=v.size()) {
				System.err.println("Error - your maximum Index is "+maxIndex+" but there are only "+v.size()+" records in this file. You might want to pick a file with complete data or start form another source...");
			} else {
				System.out.println("Maximum Index was "+maxIndex+", found "+v.size()+" records.");
			}

			chrs = new byte[v.size()];
			poslar = new int[v.size()];
			indices = new int[v.size()];
			for (int i = 0; i<v.size(); i++) {
				chrs[i] = Byte.parseByte(v.elementAt(i)[1]);
				poslar[i] = Integer.parseInt(v.elementAt(i)[2]);
				indices[i] = Integer.parseInt(v.elementAt(i)[3]);
			}

			keys = Sort.quicksort(indices);

			writer = new PrintWriter(new FileWriter(MARKERS_IN_INDEX_ORDER));
			for (int i = 0; i<v.size(); i++) {
				line = v.elementAt(keys[i]);
				writer.println(line[0]);
			}
			writer.close();

			keys = Sort.orderTwoLayers(chrs, poslar);

			writer = new PrintWriter(new FileWriter("potential_"+TARGET_MARKERS));
			writer.println("Name\tGenomeBuild\tChr\tMapInfo\tIndex");
			for (int i = 0; i<v.size(); i++) {
				line = v.elementAt(keys[i]);
				writer.println(line[0]+"\t36\t"+line[1]+"\t"+line[2]+"\t"+line[3]);
			}
			writer.close();

			writer = new PrintWriter(new FileWriter("potential_"+CNV_MARKER_POSITIONS));
			writer.println("Name\tChr\tPosition");
			for (int i = 0; i<v.size(); i++) {
				line = v.elementAt(keys[i]);
				writer.println(line[0]+"\t"+line[1]+"\t"+line[2]);
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+(filename==null?files[0]:filename)+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+(filename==null?files[0]:filename)+"\"");
			System.exit(2);
		}
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String pedfile = PED_FILE;
		String markers = TARGET_MARKERS;
		boolean parse = false;
		boolean listIDs = false;
		boolean indexMarkers = false;
		String indexfile = null;

		String usage = "\n"+"park.gwa.ParseIllumina requires 0-4 arguments\n"+"   (1) filename (i.e. ped="+pedfile+" (default))\n"+"   (2) marker list with indices (i.e. markers="+markers+" (default))\n"+
		// " (3) directory of csv files (no spaces allowed) (i.e.
		// dir="+dir+" (default))\n" +
		"   (4) -parse (not the default))\n"+"   (5) -listIDs (not the default))\n"+"   (6) -indexMarkers (not the default; if no file is specified, then it takes the first in the source directory))\n"+"   (7) file to index (i.e. index=Table.csv (not the default))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("ped=")) {
				pedfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("markers=")) {
				markers = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("index=")) {
				indexfile = args[i].split("=")[1];
				indexMarkers = true;
				numArgs--;
				// } else if (args[i].startsWith("dir=")) {
				// dir = args[i].split("=")[1];
				// numArgs--;
			} else if (args[i].equals("-parse")) {
				parse = true;
				numArgs--;
			} else if (args[i].equals("-listIDs")) {
				listIDs = true;
				numArgs--;
			} else if (args[i].equals("-indexMarkers")) {
				indexMarkers = true;
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (listIDs) {
				getListOfSamples();
			} else if (indexMarkers) {
				indexMarkers(indexfile);
			} else if (parse) {
				parseToIndividualPedFiles();
			} else {
				createPlink(pedfile, markers);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
