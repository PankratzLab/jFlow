package dead;

import java.io.*;
import java.util.*;

import common.*;

public class ParseADNI2Plink {
	public static final String WINDOWS_DIRECTORY = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\00src\\";

	public static final String LINUX_DIRECTORY = "/work/adnigwas/ADNI 676 Subjects";

	public static final String PED_DIRECTORY = "peds/";

	public static final String CNV_DIRECTORY = "cnvs/";

	public static final String MARKER_FILE = "markers.dat";

	public static final String PED_FILE = "pedfile.dat";

	public static final String ALL_MARKERS = "snpPositions.xln";

	public static final String TARGET_MARKERS = "markerListWithIndices.dat";

	public static final double GC_SCORE_CUTOFF = 0.25;

	public static final double[] POSSIBLE_CUTOFFS = {0.15, 0.25};

	public static final int NUM_MARKERS = 620901;

	public static final boolean MAKE_CNV_FILES = true;

	public static final String[] FIELDS = {"SNP Name", "Sample ID", "Sample Name", "GC Score", "SNP Index", "Chr", "Position", "Allele1 - Forward", "Allele2 - Forward", "Allele1 - AB", "Allele2 - AB", "B Allele Freq", "Log R Ratio"};

	public static void parseCIDRtoIndividualPedFiles() {
		BufferedReader reader;
		PrintWriter writer = null, missWriter, cnvWriter = null;
		String[] line;
		int count, index;
		String id = "", trav, temp;
		int[] indices;
		int[][] countMissAcrossMarkers;
		int[] countMissAcrossIndividuals;
		double score;
		String[][] snps;
		char[][] data;

		if (new File(WINDOWS_DIRECTORY).exists()) {
			trav = WINDOWS_DIRECTORY;
		} else if (new File(LINUX_DIRECTORY).exists()) {
			trav = LINUX_DIRECTORY;
		} else {
			trav = "";
			System.err.println("Error - could not resolve directory to parse (none of the following worked):");
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
		System.out.println("Found "+files.length+" files to parse");
		System.out.println("Assuming there are "+NUM_MARKERS+" markers being read in");

		snps = new String[NUM_MARKERS][3];
		countMissAcrossMarkers = new int[snps.length][POSSIBLE_CUTOFFS.length+2];

		new File(PED_DIRECTORY).mkdirs();
		if (MAKE_CNV_FILES) {
			new File(CNV_DIRECTORY).mkdirs();
		}

		System.out.println(ext.getTime());
		try {
			missWriter = new PrintWriter(new FileWriter("missingByIndividual.xln"));
			missWriter.println("Marker\t"+Array.toStr(POSSIBLE_CUTOFFS)+"\tNotCalled\tNotPresent");
			for (int i = 0; i<files.length; i++) {
				try {
					reader = new BufferedReader(new FileReader(files[i]));
					do {
						temp = reader.readLine();
					} while (reader.ready()&&!temp.contains("SNP Name")&&!temp.contains("Sample ID"));

					if (!reader.ready()) {
						System.err.println("Error - went through entire file without finding a line containing both 'SNP Name' and 'Sample ID'");
						System.exit(1);
					}
					indices = ext.indexFactors(FIELDS, temp.trim().split(","), false, true);

					count = 0;
					countMissAcrossIndividuals = new int[POSSIBLE_CUTOFFS.length+2];
					data = new char[snps.length][2];
					while (reader.ready()) {
						line = reader.readLine().split(",");
						trav = line[indices[ext.indexOfStr("Sample ID", FIELDS)]];

						if (count==0) {
							id = trav;
							if (MAKE_CNV_FILES) {
								cnvWriter = new PrintWriter(new FileWriter(CNV_DIRECTORY+id));
								cnvWriter.println("Name\tChr\tPosition\t"+id+".GType\t"+id+".Log R Ratio\t"+id+".B Allele Freq");
							}
						} else if (!trav.equals(id)) {
							System.err.println("Found "+trav+" -- expecting "+id+" in file "+files[i].getName());
							System.exit(1);
						}
						index = Integer.parseInt(line[indices[ext.indexOfStr("SNP Index", FIELDS)]])-1;
						if (snps[index][0]==null) {
							snps[index][0] = line[indices[ext.indexOfStr("SNP Name", FIELDS)]];
							snps[index][1] = line[indices[ext.indexOfStr("Chr", FIELDS)]];
							snps[index][2] = line[indices[ext.indexOfStr("Position", FIELDS)]];
						} else if (!line[indices[ext.indexOfStr("SNP Name", FIELDS)]].equals(snps[index][0])) {
							System.err.println("Error - two SNPs with the same index ("+index+") in file "+files[i].getName()+": "+snps[index][0]+" and "+line[indices[ext.indexOfStr("SNP name", FIELDS)]]);
							System.exit(1);
						}

						score = Double.parseDouble(line[indices[ext.indexOfStr("GC Score", FIELDS)]]);
						for (int j = 0; j<POSSIBLE_CUTOFFS.length; j++) {
							if (score<POSSIBLE_CUTOFFS[j]) {
								countMissAcrossMarkers[count][j]++;
								countMissAcrossIndividuals[j]++;
							}
						}
						if (score<GC_SCORE_CUTOFF) {
							line[indices[ext.indexOfStr("Allele1 - Forward", FIELDS)]] = "0";
							line[indices[ext.indexOfStr("Allele2 - Forward", FIELDS)]] = "0";
							countMissAcrossMarkers[count][POSSIBLE_CUTOFFS.length]++;
							countMissAcrossIndividuals[POSSIBLE_CUTOFFS.length]++;
						}

						data[index][0] = line[indices[ext.indexOfStr("Allele1 - Forward", FIELDS)]].charAt(0);
						data[index][1] = line[indices[ext.indexOfStr("Allele2 - Forward", FIELDS)]].charAt(0);

						if (MAKE_CNV_FILES) {
							cnvWriter.print(snps[index][0]+"\t"+snps[index][1]+"\t"+snps[index][2]);
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

					writer = new PrintWriter(new FileWriter(PED_DIRECTORY+id+".ped"));
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
					System.err.println("Error reading file \""+files[i].getName()+"\"");
					System.exit(2);
				}
			}
			missWriter.close();

			writer = new PrintWriter(new FileWriter(MARKER_FILE));
			for (int i = 0; i<snps.length; i++) {
				writer.println(snps[i][0]);
			}
			writer.close();

			missWriter = new PrintWriter(new FileWriter("missingByMarker.xln"));
			missWriter.println("SNP Index\tSNP Name\tChr\tPosition\t"+Array.toStr(POSSIBLE_CUTOFFS)+"\tNotCalled\tNotPresent");
			for (int i = 0; i<snps.length; i++) {
				missWriter.println((i+1)+"\t"+snps[i][0]+"\t"+snps[i][1]+"\t"+snps[i][2]+"\t"+Array.toStr(countMissAcrossMarkers[i]));
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
			writer = new PrintWriter(new FileWriter("plink.map"));
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
			reader = new BufferedReader(new FileReader(MARKER_FILE));
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
			System.err.println("Error: file \""+MARKER_FILE+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+MARKER_FILE+"\"");
			System.exit(2);
		}
		try {
			reader = new BufferedReader(new FileReader(pedigreeFile));
			writer = new PrintWriter(new FileWriter("plink.ped"));
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				writer.print(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]);
				if (line[6].equals(".")) {
					for (int i = 0; i<indices.length; i++) {
						writer.print("\t0\t0");
					}
				} else if (!new File(PED_DIRECTORY+line[6]+".ped").exists()) {
					System.err.println("Error - the DNA# "+line[6]+" was listed in the pedigree file but the following file was nowhere to be found: "+PED_DIRECTORY+line[6]+".ped");
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
			System.err.println("Error - could not resolve directory to parse (none of the following worked):");
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
						System.err.println("Error - went through entire file without finding a line containing both 'SNP Name' and 'Sample ID'");
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

	public static void indexMarkers() {
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
			System.err.println("Error - could not resolve directory to parse (none of the following worked):");
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
						System.err.println("Error - went through entire file without finding a line containing both 'SNP Name' and 'Sample ID'");
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

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String pedfile = PED_FILE;
		String markers = TARGET_MARKERS;
		boolean parseCIDR = false;
		boolean listIDs = false;
		boolean indexMarkers = false;

		String usage = "\n"+"park.gwa.ParseADNI2Plink requires 0-4 arguments\n"+"   (1) filename (i.e. ped="+pedfile+" (default))\n"+"   (2) marker list with indices (i.e. markers="+markers+" (default))\n"+"   (3) -parseCIDR (not the default))\n"+"   (4) -listIDs (not the default))\n"+"   (5) -indexMarkers (not the default))\n"+"";

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
			} else if (args[i].equals("-parseCIDR")) {
				parseCIDR = true;
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
				indexMarkers();
			} else if (parseCIDR) {
				parseCIDRtoIndividualPedFiles();
			} else {
				createPlink(pedfile, markers);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
