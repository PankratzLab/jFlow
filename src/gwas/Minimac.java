package gwas;

import java.io.*;
import java.util.*;

import stats.ContingencyTable;
import bioinformatics.Sequence;
import common.*;
import filesys.SnpMarkerSet;

public class Minimac {
	public static final String REF_HAP_ROOT = "/share/archive/1000G_phased/hap/EUR/EUR.chr#.hap.gz";
	public static final String REF_SNPS_ROOT = "/share/archive/1000G_phased/snps/EUR.chr#.snps";
	public static final String REF_FREQ_ROOT = "/share/archive/1000G_phased/freq/EUR.chr#_freq.xln";
//	public static final String REF_HAP_ROOT = "/share/archive/1000G_phased/hg18/0908_CEU/hap/0908_CEU_NoSingleton_chr#.hap";
//	public static final String REF_SNPS_ROOT = "/share/archive/1000G_phased/hg18/0908_CEU/snps/chr#.snps";
//	public static final String REF_FREQ_ROOT = "/share/archive/1000G_phased/hg18/0908_CEU/freq/0908_CEU_NoSingleton_chr#_freq.xln";
	public static final String MINIMAC = "/home/npankrat/bin/minimac";
	public static final String MACH2DAT = "/share/apps/bin/mach2dat"; // needed to be recompiled, per the recommendation of the website
	public static final String BGL_TO_PED = "/home/npankrat/bin/bgl_to_ped";
	public static final String[] FREQ_HEADER = {"Marker", "numMissing", "numA", "numC", "numG", "numT", "%A", "%C", "%G", "%T"};
	public static final double CHISQ_THRESHOLD = 15.05;

	public static void splitBglToPedOutputToMachHaplotypes(String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter((filename.lastIndexOf(".")>0?filename.substring(0, filename.lastIndexOf(".")):filename)+".haps"));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				for (int i = 0; i < 2; i++) {
					writer.print(line[0]+"->"+line[1]+" HAPLO"+(i+1)+" ");
					for (int j = 6+i; j < line.length; j=j+2) {
						writer.print(line[j]);
					}
					writer.println();
				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename+"\"");
			System.exit(2);
		}
	}
	
	public static void parseLogfile(String logfile) {
		BufferedReader reader;
		String[] line;
		String temp;
		Vector<String> flips, mismatches, others;
		boolean done;
		
		flips = new Vector<String>();
		mismatches = new Vector<String>();
		others = new Vector<String>();
		try {
			reader = new BufferedReader(new FileReader(logfile));
			do {
				temp = reader.readLine();
			} while (reader.ready() && !temp.startsWith("Loading target haplotypes"));
			
			done = false;
			while (!done) {
				temp = reader.readLine();
				line = temp.trim().split("[\\s]+");
				
				try {
				
					if (temp.indexOf("Target Haplotypes Loaded") > 0 || temp.indexOf("markers") > 0) {
						break;
					} else if (temp.substring(temp.indexOf("[")+1,temp.indexOf("]")).split(",").length == 3) {
						if (temp.indexOf("Possible strand flip") > 0) {
							others.add(line[4].substring(1, line[4].indexOf("':"))+"\t"+line[11]+"\t"+line[5]+"\t"+line[7]+"\t"+line[9].substring(0,line[9].length()-1));
						} else if (temp.indexOf("Mismatched frequencies") > 0) {
							others.add(line[3].substring(1, line[3].indexOf("':"))+"\t"+line[10]+"\t"+line[4]+"\t"+line[6]+"\t"+line[8].substring(0,line[8].length()-1));
						} else {
							System.err.println("Error - don't know what to do with: "+temp);
							others.add(temp);
						}
					} else if (temp.indexOf("Possible strand flip") > 0) {
						flips.add(line[4].substring(1, line[4].indexOf("':")));
					} else if (temp.indexOf("Mismatched frequencies") > 0) {
						mismatches.add(line[3].substring(1, line[3].indexOf("':"))+"\t"+line[10]+"\t"+line[4]+"\t"+line[6]+"\t"+line[8].substring(0,line[8].length()-1));
					} else {
						System.err.println("Error - don't know what to do with: "+temp);
						others.add(temp);
					}
				
				} catch (Exception e) {
					System.err.println("Error - problem parsing: "+temp);
					e.printStackTrace();
				}
			}
			
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + logfile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + logfile + "\"");
			System.exit(2);
		}

		Files.writeList(Array.toStringArray(flips), "flips.txt");
		Files.writeList(Array.toStringArray(mismatches), "mismatches.txt");
		Files.writeList(Array.toStringArray(others), "problems.txt");
	}
	
	public static void filterHaplotypes(String hapFile, String mapFile, String flips, String drops, String newHapFilename, String newMapFilename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Hashtable<String, String> flipHash, dropHash;
		int count;
		String[] markerNames, finalMarkerNames;
		int mapType;
		boolean[] flip, drop;
		char[] alleles;
		
		mapType = SnpMarkerSet.determineType(mapFile);
		if (mapType == -1) {
			System.err.println("Assuming map file ('"+mapFile+"') is just a list of marker names");
			mapType = SnpMarkerSet.NAMES_ONLY;
		}
		markerNames = new SnpMarkerSet(mapFile, mapType, false, new Logger()).getMarkerNames();
		
		flipHash = flips==null?new Hashtable<String, String>():HashVec.loadFileToHashNull(flips, false);
		dropHash = drops==null?new Hashtable<String, String>():HashVec.loadFileToHashNull(drops, false);
		flip = new boolean[markerNames.length];
		drop = new boolean[markerNames.length];
		for (int i = 0; i < markerNames.length; i++) {
			flip[i] = flipHash.containsKey(markerNames[i]);
			drop[i] = dropHash.containsKey(markerNames[i]);
		}
		
		try {
			reader = new BufferedReader(new FileReader(hapFile));
			writer = new PrintWriter(new FileWriter(newHapFilename));
			count = 0;
			while (reader.ready()) {
				count++;
				line = reader.readLine().trim().split("[\\s]+");
				if (line[2].length() != markerNames.length) {
					if (line[2].length() == 1) {
						System.err.println("Error - remove all whitespace between alleles in file '"+hapFile+"'");
					} else {
						System.err.println("Error - mismatched number of alleles at line "+count+" (expecting "+markerNames.length+"; found "+line[2].length()+")");
					}
					return;
				}
				writer.print(line[0]+" "+line[1]+" ");
				alleles = line[2].toCharArray();
				for (int i = 0; i < markerNames.length; i++) {
					if (flip[i]) {
						writer.print(Sequence.flip(alleles[i]));
					} else if (!drop[i]) {
						writer.print(alleles[i]);
					}
				}
				writer.println();
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + hapFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + hapFile + "\"");
			System.exit(2);
		}
			
		finalMarkerNames = new String[markerNames.length-Array.booleanArraySum(drop)];
		count = 0;
		for (int i = 0; i < markerNames.length; i++) {
			if (!drop[i]) {
				finalMarkerNames[count] = markerNames[i];
				count++;
			}
		}
		Files.writeList(finalMarkerNames, newMapFilename);
	}
	
	public static void extractMarkerGenotypesFromMultipleHaplotypes(String hapFileFormat, String mapFileFormat, String markerFile) {
		BufferedReader reader;
		String[] line, markers, list;
		int[] indices;
		int chr;
		Vector<String> merges;
		String command;
		
		markers = new String[24];
		try {
			reader = new BufferedReader(new FileReader(markerFile));
			reader.mark(5000);
			line = reader.readLine().trim().split("[\\s]+");
			indices = ext.indexFactors(new String[][] {{"MarkerName", "SNP", "RSID"}, {"Chr"}}, line, false, true, true, false);
			if (Array.min(indices) == -1) {
				System.err.println("Error - assuming there is no header; found '"+Array.toStr(line)+"' as first line");
				reader.reset();
			}
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				chr = Integer.parseInt(line[indices[1]]);
				markers[chr] = markers[chr]==null?line[0]:markers[chr]+","+line[0];
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + markerFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + markerFile + "\"");
			System.exit(2);
		}
		
		merges = new Vector<String>();
		for (int i = 0; i < markers.length; i++) {
			if (markers[i] != null) {
				list = markers[i].split(",");
				System.out.println("Extracting "+list.length+" marker"+(list.length>1?"s":"")+" from chr"+i);
//				extractMarkerGenotypesFromHaplotypes(ext.replaceAllWith(hapFileFormat, "#", i+""), ext.replaceAllWith(mapFileFormat, "#", i+""), list, ext.rootOf(markerFile)+"_chr"+i);
				merges.add(ext.rootOf(markerFile)+"_chr"+i+".ped "+ext.rootOf(markerFile)+"_chr"+i+".map");
			}
		}
		if (merges.size() > 1) {
			command = "plink --file "+ext.rootOf(merges.remove(0).split("[\\s]+")[0])+" --merge-list merges.txt --make-bed --out "+ext.rootOf(markerFile);
			Files.writeList(Array.toStringArray(merges), "merges.txt");
			CmdLine.run(command, "./");
			
		}
	}
	
	public static void extractMarkerGenotypesFromHaplotypes(String hapFile, String mapFile, String[] markers, String rootForOutput) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String trav;
		Hashtable<String, String> hash;
		int[] indices;
		int max, count;
		String[][][] alleles;
		String[] ids;
		String[] map;
		
		indices = ext.indexFactors(markers, HashVec.loadFileToStringArray(mapFile, false, new int[] {1}, false), false, false);
		hash = new Hashtable<String, String>();
		for (int i = 0; i < indices.length; i++) {
			hash.put(indices[i]+"", i+"");
		}
		max = Array.max(indices);
		try {
			count = Files.countLines(hapFile, false);
			System.out.println("Found "+count+" indiviudals to parse");
			reader = new BufferedReader(new FileReader(hapFile));
			alleles = new String[count][markers.length][2];
			ids = new String[count/2];
			count = 0;
			while (reader.ready()) {
				for (int i = 0; i < 2; i++) {
					line = reader.readLine().trim().split("[\\s]+");
					if (line.length == 1) {
						ids[count] = (count+1)+"";
						trav = line[0];
					} else {
						ids[count] = line[0];
						trav = line[1];
						if (trav.startsWith("HAPLO")) {
							trav = line[2];
						}
					}
					
					for (int j = 0; j <= max; j++) {
						if (hash.containsKey(j+"")) {
							alleles[count][Integer.parseInt(hash.get(j+""))][i] = trav.charAt(j)+"";
						}
					}
				}
				count++;
			}
			reader.close();
			
			writer = new PrintWriter(new FileWriter(rootForOutput+".ped"));
			for (int i = 0; i < ids.length; i++) {
				writer.print(ids[i]+"\t"+ids[i]+"\t0\t0\t1\t1");
				for (int j = 0; j < markers.length; j++) {
					if (alleles[i][j][0] == null) {
						writer.print("\t0\t0");
					} else {
						writer.print("\t"+alleles[i][j][0]+"\t"+alleles[i][j][1]);
					}
				}
				writer.println();
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + hapFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + hapFile + "\"");
			System.exit(2);
		}

		try {
			reader = new BufferedReader(new FileReader(mapFile));
			count = 0;
			map = new String[markers.length];
			while (reader.ready()) {
				trav = reader.readLine();
				if (hash.containsKey(count+"")) {
					map[Integer.parseInt(hash.get(count+""))] = trav;
				}
				count++;
			}
			reader.close();
			
			writer = new PrintWriter(new FileWriter(rootForOutput+".map"));
			for (int i = 0; i < markers.length; i++) {
				if (map[i] == null) {
					writer.println("0\t"+markers[i]+"\t0\t"+i);
				} else {
					line = Array.insertStringAt("0", map[i].trim().split("[\\s]+"), 2);
					writer.println(Array.toStr(line));
				}
			}
			writer.close();

			writer = new PrintWriter(new FileWriter(rootForOutput+".info"));
			for (int i = 0; i < markers.length; i++) {
				if (map[i] == null) {
					writer.println(markers[i]+"\t"+i);
				} else {
					line = Array.removeFromArray(map[i].trim().split("[\\s]+"), 0);
					writer.println(Array.toStr(line));
				}
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + mapFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + mapFile + "\"");
			System.exit(2);
		}
		
		Files.writeList(new String[] {"java -jar /home/npankrat/Haploview.jar -pedfile "+rootForOutput+".ped -info "+rootForOutput+".info"}, rootForOutput+".bat");
	}
	
	public static void batch(boolean beagle, boolean update, int memRequiredInGb) {
		String commands;
		
		// might consider using a local temp directory, especially for the BGL_TO_PED conversion, to free up packets back and forth, and then moving the final .dose file over to wherever it needs to be
		// other options would be to process the BGL_TO_PED conversion serially before the major jobs are started or tie it to the end of the Beagle run
		
		commands = "cd chr#\n";
		if (beagle) {
			commands += "gunzip phased.pre_phase.bgl.phased.gz\n"+
			   			(update?"cat phased.pre_phase.bgl.phased | fgrep -v id | cat ../new_header - > updated.phased\n"+
			   			BGL_TO_PED+" updated.phased ../new_plink.fam 0 > target.ped\n"+ // transform updated
			   			"rm updated.phased\n"
			   			:BGL_TO_PED+" phased.pre_phase.bgl.phased plink.fam 0 > target.ped\n")+ // transform original
			   			"gzip phased.pre_phase.bgl.phased\n"+
			   			Files.JCP+"gwas.Minimac split=target.ped\n"+
			   			"awk '{print $2}' plink.map > target.snps\n";		// split

		} else {
			commands += "cp MACH_step2_chr# target.haps\n"+
						"awk '{print $2}' all.chr#.map > target.snps\n";
		}
		
		commands += Files.JCP+"gwas.Minimac -freq hapFile=target.haps mapFile=target.snps\n"+	// compute allele frequencies
		       		Files.JCP+"gwas.Minimac compStrand=target_freq.xln compRef="+REF_FREQ_ROOT+"/\n"+	// process strand issues
		       		Files.JCP+"gwas.Minimac -filter hapFile=target.haps mapFile=target.snps flips=flips.txt drops=problems.txt newHapFile=final.haps newMapFile=final.snps\n"+	// filter
		       		"rm target.ped target.haps pre_phase.bgl plink.ped\n"+		// delete large files
		       		"gzip final.haps\n"+	// compress
		       		//			       "mkdir \n"+	// create temp directory
		       		MINIMAC+" --refHaps "+REF_HAP_ROOT+" --refSnps "+REF_SNPS_ROOT+" --haps final.haps.gz --snps final.snps --rounds 5 --states 200 --prefix chr# > final_mini.log \n"+ // final run
		       		"\n";	// post process?

		Files.qsub("mini", 1, 22, commands, memRequiredInGb);
	}
	
	public static void compareHaps(String filenames, String map, int numToCompare) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Hashtable<String, String> hash;
		String[] markerNames, files;
		int[][] markerIndices;
		char[][][] allelesCubed;
		char[] alleles;
		String[][] ids;
		String consensus;
		
		markerNames = new SnpMarkerSet(map).getMarkerNames();
		files = filenames.split(",");
		markerIndices = new int[files.length][markerNames.length];
		for (int i = 0; i < files.length; i++) {
			hash = HashVec.loadFileToHashString(ext.rootOf(files[i], false)+".snps", new int[] {0}, new int[] {-7}, false, "\t", false, false, false);
			for (int j = 0; j < markerNames.length; j++) {
				markerIndices[i][j] = hash.containsKey(markerNames[j])?Integer.parseInt(hash.get(markerNames[j])):-1;
			}
		}
		ids = new String[numToCompare][2];
		allelesCubed = new char[numToCompare][files.length][markerNames.length];
		for (int i = 0; i < files.length; i++) {
			try {
				reader = new BufferedReader(new FileReader(files[i]));
				for (int j = 0; j < numToCompare && reader.ready(); j++) {
					line = reader.readLine().trim().split("[\\s]+");
					if (i==0) {
						ids[j] = new String[] {line[0], line[1]};
					} else if (!ids[j][0].equals(line[0]) || !ids[j][1].equals(line[1])) {
						System.err.println("Error - mismatched IDs ("+Array.toStr(ids[j], "/")+") in "+files[0]+" and "+files[i]);
					}
					alleles = line[2].toCharArray();
					for (int k = 0; k < markerNames.length; k++) {
						allelesCubed[j][i][k] = markerIndices[i][k]==-1?' ':alleles[markerIndices[i][k]];
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + files[i] + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + files[i] + "\"");
				System.exit(2);
			}
		}
		
		try {
			writer = new PrintWriter(new FileWriter(ext.parseDirectoryOfFile(files[0])+"haplotype_comparison.xln"));
			writer.print("Marker");
			for (int i = 0; i < ids.length; i++) {
				for (int j = 0; j < files.length+1; j++) {
					writer.print("\t"+ids[i][0]);
				}
			}
			writer.println();
			for (int i = 0; i < ids.length; i++) {
				for (int j = 0; j < files.length+1; j++) {
					writer.print("\t"+ids[i][1]);
				}
			}
			writer.println();
			for (int i = 0; i < ids.length; i++) {
				for (int j = 0; j < files.length; j++) {
					writer.print("\t"+ext.rootOf(files[j]));
				}
				writer.print("\tconsensus");
			}
			writer.println();
			for (int i = 0; i < markerNames.length; i++) {
				writer.print(markerNames[i]);
				for (int j = 0; j < ids.length; j++) {
					consensus = "";
					for (int k = 0; k < files.length; k++) {
						writer.print("\t"+allelesCubed[j][k][i]);
					}
					writer.print("\t"+consensus);
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.parseDirectoryOfFile(files[0])+"haplotype_comparison.xln");
			e.printStackTrace();
		}
	}
	
	public static void freq(String hapFile, String mapFile) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String trav;
		int count;
		String[] markerNames;
		char[] alleles;
		int alleleIndex;
		int[][] alleleCounts;
		double sum;
		boolean makeUpperCase;
		
		makeUpperCase = false;
		
		markerNames = HashVec.loadFileToStringArray(mapFile, false, new int[] {0}, false);
		alleleCounts = new int[markerNames.length][5];
		try {
			reader = new BufferedReader(new FileReader(hapFile));
			count = 0;
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				count++;
				if (line.length == 1) {
					trav = line[0];
				} else {
					trav = line[1];
					if (trav.startsWith("HAPLO")) {
						trav = line[2];
					}
				}
				alleles = (makeUpperCase?trav.toUpperCase():trav).toCharArray();
				if (alleles.length != markerNames.length) {
					System.err.println("Error - number of markers in line "+count+" of "+ext.removeDirectoryInfo(hapFile)+" ("+alleles.length+") does not match the number of markers in "+ext.rootOf(hapFile, true)+".snps ("+markerNames.length+")");
					return;
				}
				for (int i = 0; i < markerNames.length; i++) {
					alleleIndex = ext.indexOfChar(alleles[i], Sequence.ALLELES);
					if (alleleIndex == -1) {
						if (ext.indexOfChar((alleles[i]+"").toUpperCase().charAt(0), Sequence.ALLELES) == -1) {
							System.err.println("Error - unknown allele code '"+alleles[i]+"'; add to null chars if this is a valid missing value");
						} else {
							makeUpperCase = true;
							alleles = trav.toUpperCase().toCharArray();
							alleleIndex = ext.indexOfChar(alleles[i], Sequence.ALLELES);
						}
					}
					alleleCounts[i][alleleIndex+1]++;
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + hapFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + hapFile + "\"");
			System.exit(2);
		}
		
		try {
			writer = new PrintWriter(new FileWriter(ext.rootOf(hapFile, true)+"_freq.xln"));
			writer.println(Array.toStr(FREQ_HEADER));
			for (int i = 0; i < markerNames.length; i++) {
				writer.print(markerNames[i]+"\t"+alleleCounts[i][0]);
				sum = 0;
				for (int j = 0; j < 4; j++) {
					writer.print("\t"+alleleCounts[i][j+1]);
					sum += alleleCounts[i][j+1];
				}
				for (int j = 0; j < 4; j++) {
					writer.print("\t"+ext.formDeci((double)alleleCounts[i][j+1]/sum, 4));
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.rootOf(hapFile, true)+"_freq.xln");
			e.printStackTrace();
		}
	}
	
	public static void allFreq() {
		PrintWriter writer;
		
		try {
			writer = new PrintWriter(new FileWriter("batchFreq"));
			writer.println("mkdir /share/archive/1000G_phased/freq/");
			writer.println("cd /share/archive/1000G_phased/freq/");
			writer.println("mkdir /share/archive/1000G_phased/hap/EUR/unzipped/");
			for (int chr = 1; chr <= 22; chr++) {
				writer.println("gunzip -c /share/archive/1000G_phased/hap/EUR/EUR.chr"+chr+".hap.gz > /share/archive/1000G_phased/hap/EUR/unzipped/EUR.chr"+chr+".hap");
				writer.println(Files.JCP+"gwas.Minimac -freq hapFile=/share/archive/1000G_phased/hap/EUR/unzipped/EUR.chr"+chr+".hap mapFile=/share/archive/1000G_phased/snps/EUR.chr"+chr+".snps");
				writer.println("rm /share/archive/1000G_phased/hap/EUR/unzipped/EUR.chr"+chr+".hap");
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + "batchFreq");
			e.printStackTrace();
		}
	}
	
	public static void compareStrands(String compStrand, String compRef) {
		BufferedReader reader;
		String[] line, target, ref;
		Hashtable<String, String> hash;
		boolean[][] matchup;
		Vector<String> flips, mismatches, others;
		double chisq;
		int numObserved;
		
		flips = new Vector<String>();
		mismatches = new Vector<String>();
		others = new Vector<String>();
		
		hash = HashVec.loadFileToHashString(compStrand, new int[] {0}, new int[] {2,3,4,5}, false, "\t", true, false, false);
		try {
			reader = new BufferedReader(new FileReader(compRef));
			line = reader.readLine().trim().split("[\\s]+");
			ext.checkHeader(line, FREQ_HEADER, false);
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (hash.containsKey(line[0])) {
					target = hash.get(line[0]).split("[\\s]+");
					ref = Array.subArray(line, 2, 6);
					matchup = new boolean[2][4];
					for (int i = 0; i < 2; i++) {
						for (int j = 0; j < 4; j++) {
							matchup[i][j] = !(i==0?target:ref)[j].equals("0");
						}
					}
					chisq = ContingencyTable.ChiSquare(Matrix.prune(new int[][] {Array.toIntArray(target), Array.toIntArray(ref)}), false);
					numObserved = Array.booleanArraySum(allelesUsed(matchup));
					if (Array.booleanArraySum(matchup[0]) > 2 || Array.booleanArraySum(matchup[1]) > 2) {
						others.add(line[0]+"\t"+ext.formDeci(chisq, 1, true)+"\t"+displayFreqs(matchup, target, ref));
					} else if (numObserved > 2) {
						if (Array.booleanArraySum(allelesUsed(new boolean[][] {flip(matchup[0]), matchup[1]})) > 2) {
							others.add(line[0]+"\t"+ext.formDeci(chisq, 1, true)+"\t"+displayFreqs(matchup, target, ref));
						} else {
							flips.add(line[0]+"\t"+displayFreqs(matchup, target, ref));
						}
					} else if (chisq >= CHISQ_THRESHOLD) {
						mismatches.add(line[0]+"\t"+ext.formDeci(chisq, 1, true)+"\t"+displayFreqs(matchup, target, ref));
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + compRef + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + compRef + "\"");
			System.exit(2);
		}
		
		Files.writeList(Array.toStringArray(flips), "flips.txt");
		Files.writeList(Array.toStringArray(mismatches), "mismatches.txt");
		Files.writeList(Array.toStringArray(others), "problems.txt");

	}

	public static String[] flip(String[] array) {
		String[] flipped;
		
		flipped = new String[array.length];
		flipped[0] = array[3];
		flipped[1] = array[2];
		flipped[2] = array[1];
		flipped[3] = array[0];
		
		return flipped;
	}
	
	public static boolean[] flip(boolean[] array) {
		boolean[] flipped;
		
		flipped = new boolean[array.length];
		flipped[0] = array[3];
		flipped[1] = array[2];
		flipped[2] = array[1];
		flipped[3] = array[0];
		
		return flipped;
	}
	
	public static boolean[] allelesUsed(boolean[][] matchup) {
		boolean[] allelesUsed;

		allelesUsed = Array.booleanArray(4, false);
		for (int i = 0; i < matchup.length; i++) {
			for (int j = 0; j < matchup[i].length; j++) {
				if (matchup[i][j]) {
					allelesUsed[j] = true;
				}
			}
		}

		return allelesUsed;
	}
	
	public static String displayFreqs(boolean[][] matchup, String[] target, String[] ref) {
		boolean[] allelesUsed;
		String[] strs;
		boolean first;
		double[] denoms;

		allelesUsed = allelesUsed(matchup);

		strs = new String[3];
		strs[0] = "f[";
		strs[1] = "[";
		strs[2] = "[";

		first = true;
		denoms = new double[] {Array.sum(Array.toDoubleArray(target)), Array.sum(Array.toDoubleArray(ref))};
		for (int i = 0; i < allelesUsed.length; i++) {
			if (allelesUsed[i]) {
				if (first) {
					first = false;
				} else {
					for (int j = 0; j < 3; j++) {
						strs[j] += ",";
					}
				}
				strs[0] += Sequence.ALLELES[i];
				strs[1] += ext.formDeci(Double.parseDouble(target[i])/denoms[0], 2, true)+"";
				strs[2] += ext.formDeci(Double.parseDouble(ref[i])/denoms[1], 2, true)+"";
			}
		}
		for (int j = 0; j < 3; j++) {
			strs[j] += "]";
		}
		
		return Array.toStr(strs);
	}
	
	//  the mach2dat executable needed to be recompiled, per the recommendation of the website
	public static void qsubMach2dat(String pheno) {
		PrintWriter writer;
		String[] line;
		String commands;

		try {
			writer = new PrintWriter(new FileWriter(ext.addToRoot(pheno, "_desc")));
			line = Files.getHeaderOfFile(pheno, "[\\s]+", new Logger());
			writer.println("A affected");
			for (int i = 6; i < line.length; i++) {
				writer.println("C covariate"+(i-5));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.addToRoot(pheno, "_desc"));
			e.printStackTrace();
		}		

		commands = "echo \"start at: \" `date`\n";
		commands += "cd chr#\n";
//		commands += "awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7}' chr#.info > chr#.pinfo\n";
		commands += MACH2DAT+" -d ../"+ext.addToRoot(pheno, "_desc")+" -p ../"+pheno+" -i chr#.info -g chr#.dose > "+ext.rootOf(pheno)+"_chr#.out\n";
		commands += "cd ..\n";
		commands += "echo \"end at: \" `date`\n";
		Files.qsub("chr#_"+ext.rootOf(pheno), 1, 22, commands, -1);
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String extract = null;
		String extractSet = null;
		String[] markerNames;
		String hapFile = "chr6.haps";
		String mapFile = "chr6.snps";
		String hapFileFormat = "hap/EUR/unzipped/EUR.chr#.hap";
		String mapFileFormat = "map/EUR.chr#.map";
		String logfile = null;
		boolean filter = false;
		String flips = "flips.txt";
		String drops = "mismatches.txt";
		String split = null;
		String newHapFile = null;
		String newMapFile = null;
		boolean batch = false;
		boolean update = false;
		int memRequiredInGb = -1;
		String comp = null;
		String map = "plink.bim";
		int numToCompare = 4;
		boolean freq = false;
		boolean freqAll = false;
		String compStrand = null;
		String compRef = null;
		String rootForOutput;
		String mach2dat = null;
		boolean beagle = false;

		String usage = "\n" +
		"gwas.Minimac requires 0-1 arguments\n" +
		"   (1) extract specific markers from haplotype files (i.e. extract=filename.txt OR extract=rs3129882,rs2395163,chr6:32588205 (not the default))\n" + 
		"   (2) haplotype filename (i.e. hapFile=" + hapFile + " (default))\n" + 
		"   (3) map filename (i.e.  mapFile=" + mapFile + " (default))\n" + 
		" OR\n" + 
		"   (1) extract markers from a set of haplotype files using the names and chrs in this file (i.e. extractSet=filename.txt (not the default))\n" + 
		"   (2) haplotype filename format (i.e. hapFileFormat=" + hapFileFormat + " (default))\n" + 
		"   (3) map filename (i.e.  mapFileFormat=" + mapFileFormat + " (default))\n" + 
		" OR\n" + 
		"   (1) compute allele frequencies from haplotype files (i.e. -freq (not the default))\n" + 
		"   (2) haplotype filename (i.e. hapFile=" + hapFile + " (default)))\n" + 
		"   (3) map filename (i.e.  mapFile=" + mapFile + " (default))\n" + 
		" OR\n" +
		"   (1) generate frequency reports for full reference panel (i.e. -freqAll (not the default))\n" + 
		" OR\n" + 
		"   (1) minimac logfile to be parsed for flipped/mismatched markers (i.e. logfile=chr6_mini.log (not the default))\n" + 
		" OR\n" + 
		"   (1) compare alleles/strand directly to a reference freq report (i.e. compStrand=chr6_freq.xln (not the default))\n" + 
		"   (2) reference freq report to which to compare (i.e. compRef=ref/EUR.chr6_freq.xln (not the default))\n" + 
		" OR\n" + 
		"   (1) split and compress the results of a bgl_to_ped conversion (i.e. split=chr6.ped (not the default))\n" + 
		" OR\n" + 
		"   (1) filter haplotypes (i.e. -filter (not the default))\n" + 
		"   (2) haplotype filename (i.e. hapFile=" + hapFile + " (default))\n" + 
		"   (3) map filename (i.e.  mapFile=" + mapFile + " (default))\n" + 
		"   (4) list of markers to flip (i.e. flips="+flips+" (default))\n" + 
		"   (5) list of markers to drop (i.e. drop="+drops+" (default))\n" + 
		"   (6) name of resulting haplotype file (i.e. newHapFile=[oldfileRoot]_filtered.[oldFileExtension] (default))\n" + 
		"   (7) name of resulting map filename (i.e. newMapFile=[oldfileRoot]_filtered.[oldFileExtension] (default))\n" + 
		" OR\n" +
		"   (1) compare haplotypes (i.e. comp=target.haps,final.haps (not the default; need a .snps file for each root))\n" + 
		"   (2) master map file with positions (i.e. map="+map+" (default))\n" + 
		"   (3) number of lines/haplotypes to compare (i.e. numToCompare="+numToCompare+" (default))\n" + 
		" OR\n" +
		"   (1) batch convert and run (i.e. -batch (not the default))\n" + 
		"   (2) used phased beagle data instead of phased mach data (i.e. -beagle (not the default))\n" + 
		"   (3) (optional) update IDs with ../new_header and ../new_plink.fam when running (i.e. -update (not the default))\n" +
		"   (4) (optional) amount of memory in Gb to reserve on each node (i.e. mem=4 (not the default))\n" +
		" OR\n"+
		"   (1) set up qsubs for mach2dat using pheno file (i.e. mach2dat=pheno.dat (not the default))\n"+
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("extract=")) {
				extract = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("hapFile=")) {
				hapFile = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("mapFile=")) {
				mapFile = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("extractSet=")) {
				extractSet = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("hapFileFormat=")) {
				hapFileFormat = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("mapFileFormat=")) {
				mapFileFormat = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("-freq")) {
				freq = true;
				numArgs--;
			} else if (args[i].startsWith("logfile=")) {
				logfile = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("compStrand=")) {
				compStrand = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("compRef=")) {
				compRef = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("split=")) {
				split = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("-filter")) {
				filter = true;
				numArgs--;
			} else if (args[i].startsWith("flips=")) {
				flips = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("drops=")) {
				drops = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("newHapFile=")) {
				newHapFile = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("newMapFile=")) {
				newMapFile = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("-batch")) {
				batch = true;
				numArgs--;
			} else if (args[i].startsWith("-beagle")) {
				beagle = true;
				numArgs--;
			} else if (args[i].startsWith("-update")) {
				update = true;
				numArgs--;
			} else if (args[i].startsWith("mem=")) {
				memRequiredInGb = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("comp=")) {
				comp = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("map=")) {
				map = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("numToCompare=")) {
				numToCompare = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("-freqAll")) {
				freqAll = true;
				numArgs--;
			} else if (args[i].startsWith("mach2dat=")) {
				mach2dat = ext.parseStringArg(args[i], null);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
//		hapFile = "C:\\Users\\npankrat\\Downloads\\hap\\chr6.hap";
//		mapFile = "C:\\Users\\npankrat\\Downloads\\map\\chr6.map";
//		extract = "rs3129882,rs2395163,chr6:32588205";
		
//		split = "D:/target.ped";

//		logfile = "D:/logfile.out";

//		filter = true;
//		hapFile = "D:/targets.haps";
//		mapFile = "D:/targets.snps";
		
//		batch = true;
		
//		logfile = "first_mini.log";
		
//		comp = "D:\\home\\npankrat\\jProjects\\master\\chr16\\target.haps,D:\\home\\npankrat\\jProjects\\master\\chr16\\final.haps";
//		map = "D:\\home\\npankrat\\jProjects\\master\\chr16\\plink.bim";

//		freq = true;
//		hapFile = "D:\\home\\npankrat\\jProjects\\master\\chr16\\target.haps";
//		hapFile = "D:\\home\\npankrat\\jProjects\\master\\chr16\\final.haps";
//		hapFile = "D:\\home\\npankrat\\jProjects\\master\\chr16\\EUR.chr16.haps";

//		compStrand = "D:\\home\\npankrat\\jProjects\\master\\chr16\\target_freq.xln";
//		compRef = "D:\\home\\npankrat\\jProjects\\master\\chr16\\all_EUR.chr16_freq.xln";

//		extractSet = "map.txt";
		
		try {
			if (split != null) {
				splitBglToPedOutputToMachHaplotypes(split);
			} else if (extractSet != null) { 
				extractMarkerGenotypesFromMultipleHaplotypes(hapFileFormat, mapFileFormat, extractSet);
			} else if (extract != null) { 
				if (Files.exists(extract, false)) {
					markerNames = HashVec.loadFileToStringArray(extract, false, new int[] {0}, true);
					rootForOutput = ext.rootOf(extract);
				} else {
					markerNames = extract.split(",");
					rootForOutput = ext.rootOf(hapFile)+"_subset";
				}
				extractMarkerGenotypesFromHaplotypes(hapFile, mapFile, markerNames, rootForOutput);
			} else if (freq) {
				freq(hapFile, mapFile);
			} else if (freqAll) {
				allFreq();
			} else if (logfile != null) {
				parseLogfile(logfile);
			} else if (filter) {
				filterHaplotypes(hapFile, mapFile, flips, drops, newHapFile==null?ext.addToRoot(hapFile, "_filtered"):newHapFile, newMapFile==null?ext.addToRoot(mapFile, "_filtered"):newMapFile);
			} else if (batch) {
				batch(beagle, update, memRequiredInGb);
			} else if (comp != null) {
				compareHaps(comp, map, numToCompare);
			} else if (compStrand != null) {
				if (compRef == null) {
					System.err.println("Error - reference freq file not delineated");
				} else {
					compareStrands(compStrand, compRef);
				}
			} else if (mach2dat != null) {
				qsubMach2dat(mach2dat);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
