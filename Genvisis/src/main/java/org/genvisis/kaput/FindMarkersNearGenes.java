// -Xms1024M -Xmx1024M
package kaput;

import java.io.*;
import java.util.*;

import common.*;

public class FindMarkersNearGenes {
	// public static final String DEFAULT_DIR = "C:\\Documents and
	// Settings\\npankrat\\My Documents\\gwas\\merged\\results\\";
	public static final String DEFAULT_DIR = "";

	public static final String DEFAULT_ROOT = "hapmap-ceu";

	public static final String GENES_ALT = "C:\\Documents and Settings\\npankrat\\My Documents\\jProjects\\CommonTools\\runtime";

	public static final String GENES = "genes.xls";

	// public static final String DEFAULT_HIT_LIST = "hitGenes.txt";
	// public static final String DEFAULT_HIT_LIST = "Top29.txt";
	// public static final String DEFAULT_HIT_LIST = "Top86.txt";
	// public static final String DEFAULT_HIT_LIST = "chr21.txt";
	public static final String DEFAULT_HIT_LIST = "chr2G.txt";

	public static final double DEFAULT_MAF_CUTOFF = 0.05;

	public static final int DEFAULT_WINDOW_SIZE = 2000;

	// public static final String[] TAG_SUFFIXES = {"", "2", "23", "-0.9",
	// "2-0.9", "23-0.9"};
	// public static final String[] TAG_SUFFIXES = {"", "-0.9"};
	public static final String[] TAG_SUFFIXES = {"", "23", "-0.9", "23-0.9"};

	public static final int SPLIT_SIZE = 200;

	public static final int GAP_MIN = 100000;

	public static final int NUM_BATCHES = 1;

	public static void findMarkers(String dir, String filename, String root, double mafCutoff, int windowSize) {
		BufferedReader reader;
		PrintWriter writer;
		String trav;
		String[] line;
		String[] geneNames;
		int[][] genePositions;
		Hashtable<String,String> geneLookup = new Hashtable<String,String>();
		Hashtable<String,String> markerPositions = new Hashtable<String,String>();
		Hashtable<String,String> mafLookup;
		int index, chr, pos;
		Vector<String> markers;
		String[] markerNames;

		System.out.println(ext.getTime()+"  Loading genes...");
		geneNames = Array.toStringArray(HashVec.loadFileToVec(dir+filename, false, true, true));
		System.out.println("            Found "+geneNames.length+" genes...");

		System.out.println(ext.getTime()+"  Determining gene boundaries...");
		genePositions = new int[geneNames.length][];
		try {
			reader = Files.getReader(GENES, GENES_ALT);
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				index = ext.indexOfStr(line[1], geneNames);
				if (index>=0) {
					genePositions[index] = new int[3];
					genePositions[index][0] = Positions.chromosomeNumber(line[2]);
					genePositions[index][1] = Integer.parseInt(line[3])-windowSize;
					genePositions[index][2] = Integer.parseInt(line[4])+windowSize;
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+GENES+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+GENES+"\"");
			System.exit(2);
		}
		trav = "";
		for (int i = 0; i<genePositions.length; i++) {
			if (genePositions[i]==null) {
				trav += (trav.equals("")?"":", ")+geneNames[i];
			}
		}
		if (!trav.equals("")) {
			System.err.println("Error - the following genes were not found in the genes database: "+trav);
			System.exit(1);
		}

		System.out.println(ext.getTime()+"  Finding qualified SNPs...");
		markers = new Vector<String>();
		try {
			reader = new BufferedReader(new FileReader(dir+root+".bim"));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				chr = Integer.parseInt(line[0]);
				pos = Integer.parseInt(line[3]);
				for (int i = 0; i<genePositions.length; i++) {
					if (chr==genePositions[i][0]&&pos>=genePositions[i][1]&&pos<=genePositions[i][2]) {
						HashVec.addIfAbsent(line[1], markers);
						trav = geneLookup.containsKey(line[1])?geneLookup.get(line[1]):"";
						trav += (trav.equals("")?"":"; ")+geneNames[i];
						geneLookup.put(line[1], trav);
						markerPositions.put(line[1], line[0]+"\t"+line[3]);
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+dir+root+".bim"+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+dir+root+".bim"+"\"");
			System.exit(2);
		}

		System.out.println(ext.getTime()+"  Loading MAFs...");
		mafLookup = HashVec.loadFileToHashString(dir+root+".frq", 1, new int[] {4}, null, false);
		// mafLookup = new Hashtable<String, String>();

		markerNames = Array.toStringArray(markers);
		System.out.println(ext.getTime()+"  Filtering "+markers.size()+" SNPs; discarding MAF <"+mafCutoff);
		try {
			writer = new PrintWriter(new FileWriter(dir+ext.rootOf(filename)+"_markers.xln"));
			for (int i = 0; i<markerNames.length; i++) {
				trav = mafLookup.get(markerNames[i]);
				if (trav==null) {
					System.err.println("Error - "+markerNames[i]+" was listed in the map file but not the frq file");
				} else if (Double.parseDouble(trav)<mafCutoff) {
					markers.remove(markerNames[i]);
				}
				writer.println(markerNames[i]+"\t"+markerPositions.get(markerNames[i])+"\t"+(trav==null?"??":trav)+"\t"+geneLookup.get(markerNames[i]));
			}
			writer.close();
		} catch (IOException ioe) {
			System.err.println("Error writing to file \""+dir+ext.rootOf(filename)+"_markers.xln"+"\"");
			System.exit(2);
		}
		System.out.println(ext.getTime()+"  "+markers.size()+" SNPs survived filter");

		try {
			writer = new PrintWriter(new FileWriter(dir+ext.rootOf(filename)+"_qualified.txt"));
			for (int j = 0; j<markers.size(); j++) {
				writer.println(markers.elementAt(j));
			}
			writer.close();

			if (markers.size()==0) {
				writer = new PrintWriter(new FileWriter(dir+ext.rootOf(filename)+"_tags.txt"));
				writer.close();
				writer = new PrintWriter(new FileWriter(dir+ext.rootOf(filename)+"_tags-0.9.txt"));
				writer.close();
			}
		} catch (IOException ioe) {
			System.err.println("Error writing to file \""+dir+ext.rootOf(filename)+"_qualified.txt"+"\"");
			System.exit(2);
		}
	}

	public static void createPlink(String dir, String filename, String root) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;

		if (!new File(dir+ext.rootOf(filename)+".ped").exists()||!new File(dir+ext.rootOf(filename)+".info").exists()) {
			System.out.println(ext.getTime()+"  Creating plink files for root "+ext.rootOf(filename));
			CmdLine.run("plink --noweb --bfile "+root+" --extract "+ext.rootOf(filename)+"_qualified.txt --out "+ext.rootOf(filename)+" --recode", "./");
			try {
				new File(dir+ext.rootOf(filename)+".ped").renameTo(new File(dir+ext.rootOf(filename)+".ped.bak"));
				reader = new BufferedReader(new FileReader(dir+ext.rootOf(filename)+".ped.bak"));
				writer = new PrintWriter(new FileWriter(dir+ext.rootOf(filename)+".ped"));
				while (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
					line[5] = "0";
					writer.println(Array.toStr(line, " "));
				}
				reader.close();
				writer.close();
			} catch (IOException e) {
				System.err.println("Error fixing affection status "+ext.rootOf(filename)+".ped");
				System.exit(2);
			}
			new File(dir+ext.rootOf(filename)+".ped.bak").delete();
			try {
				reader = new BufferedReader(new FileReader(dir+ext.rootOf(filename)+".map"));
				writer = new PrintWriter(new FileWriter(dir+ext.rootOf(filename)+".info"));
				while (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
					writer.println(line[1]+"\t"+line[3]);
				}
				reader.close();
				writer.close();
			} catch (IOException e) {
				System.err.println("Error converting "+ext.rootOf(filename)+".map to "+ext.rootOf(filename)+".info");
				System.exit(2);
			}
			System.out.println(ext.getTime()+" ...all set!");
		}

	}

	public static String prettyUpDistance(int dist) {
		String str = dist+"";
		if (str.length()>6) {
			return str.substring(0, str.length()-6)+"Mb";
		}
		if (str.length()>3) {
			return str.substring(0, str.length()-3)+"kb";
		}
		return str;
	}

	public static int[] parseChrInfo(String[] line) {
		int[] chr_start_stop = new int[3];

		chr_start_stop[0] = line[2].equals("X")?23:(line[2].equals("Y")?24:(line[2].equals("XY")?25:(line[2].equals("MT")?26:(line[2].equals("Un")?27:Integer.parseInt(line[2])))));
		chr_start_stop[1] = Integer.parseInt(line[3]);
		chr_start_stop[2] = Integer.parseInt(line[4]);

		return chr_start_stop;
	}

	public static void summarizeMarkers(String dir, String filename, double mafCutoff, int windowSize) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String trav;
		Hashtable<String,String> seen = new Hashtable<String,String>();
		String[] genes;
		String[][] tags;
		int[][] genePositions, geneCounts;
		int index, chr, pos;

		System.out.println(ext.getTime()+"  Loading genes...");
		genes = Array.toStringArray(HashVec.loadFileToVec(dir+filename, false, true, true));
		System.out.println("            Found "+genes.length+" genes...");

		System.out.println(ext.getTime()+"  Loading tags...");
		tags = new String[TAG_SUFFIXES.length][];
		for (int i = 0; i<TAG_SUFFIXES.length; i++) {
			if (new File(dir+ext.rootOf(filename)+"_tags"+TAG_SUFFIXES[i]+".txt").exists()) {
				tags[i] = Array.toStringArray(HashVec.loadFileToVec(dir+ext.rootOf(filename)+"_tags"+TAG_SUFFIXES[i]+".txt", false, true, true));
			} else {
				System.err.println("Error - "+ext.rootOf(filename)+"_tags"+TAG_SUFFIXES[i]+".txt failed");
			}
		}
		System.out.println("            Found "+tags.length+" tags...");

		System.out.println(ext.getTime()+"  Determining gene boundaries...");
		genePositions = new int[genes.length][];
		try {
			reader = Files.getReader(GENES, GENES_ALT);
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				index = ext.indexOfStr(line[1], genes);
				if (index>=0) {
					genePositions[index] = new int[3];
					genePositions[index][0] = Positions.chromosomeNumber(line[2]);
					genePositions[index][1] = Integer.parseInt(line[3])-windowSize;
					genePositions[index][2] = Integer.parseInt(line[4])+windowSize;
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+GENES+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+GENES+"\"");
			System.exit(2);
		}
		trav = "";
		for (int i = 0; i<genePositions.length; i++) {
			if (genePositions[i]==null) {
				trav += (trav.equals("")?"":", ")+genes[i];
			}
		}
		if (!trav.equals("")) {
			System.err.println("Error - the following genes were not found in the genes database: "+trav);
			System.exit(1);
		}

		System.out.println(ext.getTime()+"  Sifting through qualified genes");
		geneCounts = Matrix.intMatrix(genes.length, 3+TAG_SUFFIXES.length, 0);
		try {
			reader = new BufferedReader(new FileReader(dir+ext.rootOf(filename)+"_markers.xln"));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				// System.err.println(line[0]);
				chr = Integer.parseInt(line[1]);
				pos = Integer.parseInt(line[2]);
				for (int i = 0; i<genePositions.length; i++) {
					if (chr==genePositions[i][0]&&pos>=genePositions[i][1]&&pos<=genePositions[i][2]) {
						geneCounts[i][0]++;
						if (seen.containsKey(line[0])) {
							geneCounts[i][1]++;
						}
						if (!line[3].equals("??")&&Double.parseDouble(line[3])>=mafCutoff) {
							geneCounts[i][2]++;
						}
						for (int j = 0; j<TAG_SUFFIXES.length; j++) {
							if (tags[j]!=null&&ext.containsAny(line[0], tags[j])) {
								geneCounts[i][3+j]++;
							}
						}
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("\""+dir+ext.rootOf(filename)+"_markers.xln"+"\" was not found; please re-run findMarkers algorithm");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+dir+ext.rootOf(filename)+"_markers.xln"+"\"");
			System.exit(2);
		}

		try {
			writer = new PrintWriter(new FileWriter(dir+ext.rootOf(filename)+"_summary.xln"));
			writer.print("Marker\tchr\tstart\tstop\tSize\tNumMarkers\tNumOverlap\tNumMAF<"+mafCutoff);
			for (int i = 0; i<TAG_SUFFIXES.length; i++) {
				writer.print("\tNumTags"+TAG_SUFFIXES[i]);
			}
			writer.println();
			for (int i = 0; i<genes.length; i++) {
				writer.print(genes[i]+"\t"+genePositions[i][0]+"\t"+genePositions[i][1]+"\t"+genePositions[i][2]+"\t"+(genePositions[i][2]-genePositions[i][1]));
				for (int j = 0; j<3; j++) {
					writer.print("\t"+geneCounts[i][j]);
				}
				for (int j = 0; j<tags.length; j++) {
					writer.print("\t"+(tags[j]==null?"??":geneCounts[i][3+j]));
				}
				writer.println();
			}
			writer.close();
		} catch (IOException ioe) {
			System.err.println("Error writing to file \""+dir+ext.rootOf(filename)+"_summary.xln"+"\"");
			System.exit(2);
		}
	}

	public static void splitGenesAndBatch() {
		BufferedReader reader;
		PrintWriter writer = null, lookup;
		String[] line = null;
		int chr = 0;
		int prev = 0, count, suffix;
		boolean done = false;
		Vector<String[]> v = new Vector<String[]>();
		String commands;

		try {
			reader = Files.getReader(GENES, GENES_ALT);
			lookup = new PrintWriter(new FileWriter("lookup.xln"));

			chr = count = suffix = 0;
			reader.readLine();

			while (!done) {
				if (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
					line[2] = Positions.chromosomeNumber(line[2])+"";
				} else {
					done = true;
				}
				if (done||Integer.parseInt(line[2])>chr||(count>=SPLIT_SIZE&&Integer.parseInt(line[3])-prev>GAP_MIN)) {
					if (count>0) {
						writer.close();
						System.out.println("Filled chr"+chr+ext.getExcelColumn(suffix)+".txt (n="+count+")");
						v.add(new String[] {"chr"+chr+ext.getExcelColumn(suffix)});
					}
					if (Integer.parseInt(line[2])>chr) {
						chr = Integer.parseInt(line[2]);
						suffix = 0;
					} else {
						suffix++;
					}
					if (chr>24) {
						done = true;
					}
					if (!done) {
						writer = new PrintWriter(new FileWriter("chr"+chr+ext.getExcelColumn(suffix)+".txt"));
						count = 0;
					}
				}
				if (!done) {
					writer.println(line[1]);
					lookup.println(line[1]+"\t"+"chr"+chr+ext.getExcelColumn(suffix));
					prev = Integer.parseInt(line[4]);
					count++;
				}
			}
			reader.close();
			lookup.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+GENES+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+GENES+"\"");
			System.exit(2);
		}

		commands = "jcpm gwa.FindMarkersNearGenes dir=./ file=[%0].txt\n"+"jcp gwa.FindMarkersNearGenes dir=./ file=[%0].txt -createPlink\n"+"java -Dsun.java2d.noddraw=true -Xms2048M -Xmx2048M -classpath Haploview.jar -Djava.awt.headless=true edu.mit.wi.haploview.HaploView -nogui -pedfile [%0].ped -info [%0].info -pairwiseTagging\n"+"mv [%0].ped.TESTS [%0]_tags.txt\n"+"java -Dsun.java2d.noddraw=true -Xms2048M -Xmx2048M -classpath Haploview.jar -Djava.awt.headless=true edu.mit.wi.haploview.HaploView -nogui -pedfile [%0].ped -info [%0].info -aggressiveTagging -aggressiveNumMarkers 3\n"+"mv [%0].ped.TESTS [%0]_tags23.txt\n"+"java -Dsun.java2d.noddraw=true -Xms2048M -Xmx2048M -classpath Haploview.jar -Djava.awt.headless=true edu.mit.wi.haploview.HaploView -nogui -pedfile [%0].ped -info [%0].info -tagrsqcutoff 0.9 -pairwiseTagging\n"+"mv [%0].ped.TESTS [%0]_tags-0.9.txt\n"+"java -Dsun.java2d.noddraw=true -Xms2048M -Xmx2048M -classpath Haploview.jar -Djava.awt.headless=true edu.mit.wi.haploview.HaploView -nogui -pedfile [%0].ped -info [%0].info -tagrsqcutoff 0.9 -aggressiveTagging -aggressiveNumMarkers 3\n"+"mv [%0].ped.TESTS [%0]_tags23-0.9.txt\n"+"rm [%0].ped.TAGS\n"+"jcpm gwa.FindMarkersNearGenes dir=./ file=[%0].txt\n"+"";
		commands = "jcpm gwa.FindMarkersNearGenes dir=./ file=[%0].txt\n";
		Files.batchIt("batch", null, NUM_BATCHES, commands, Matrix.toStringArrays(v));
	}

	public static void fillIn() {
		PrintWriter writer;
		String trav;
		int count;
		String commands;

		File[] files = new File(".").listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.startsWith("chr")&&filename.endsWith(".txt")&&!filename.contains("_");
			}
		});

		try {
			writer = new PrintWriter(new FileWriter("fillIn"));

			for (int i = 0; i<files.length; i++) {
				trav = ext.rootOf(files[i].getName());
				commands = "";
				count = 0;
				if (!new File(trav+"_tags.txt").exists()) {
					commands += "java -Dsun.java2d.noddraw=true -Xms2048M -Xmx2048M -classpath Haploview.jar -Djava.awt.headless=true edu.mit.wi.haploview.HaploView -nogui -pedfile "+trav+".ped -info "+trav+".info -pairwiseTagging\n"+"mv "+trav+".ped.TESTS "+trav+"_tags.txt\n";
					count++;
				}
				if (!new File(trav+"_tags23.txt").exists()) {
					commands += "java -Dsun.java2d.noddraw=true -Xms2048M -Xmx2048M -classpath Haploview.jar -Djava.awt.headless=true edu.mit.wi.haploview.HaploView -nogui -pedfile "+trav+".ped -info "+trav+".info -aggressiveTagging -aggressiveNumMarkers 3\n"+"mv "+trav+".ped.TESTS "+trav+"_tags23.txt\n";
					count++;
				}
				if (!new File(trav+"_tags-0.9.txt").exists()) {
					commands += "java -Dsun.java2d.noddraw=true -Xms2048M -Xmx2048M -classpath Haploview.jar -Djava.awt.headless=true edu.mit.wi.haploview.HaploView -nogui -pedfile "+trav+".ped -info "+trav+".info -tagrsqcutoff 0.9 -pairwiseTagging\n"+"mv "+trav+".ped.TESTS "+trav+"_tags-0.9.txt\n";
					count++;
				}
				if (!new File(trav+"_tags23-0.9.txt").exists()) {
					commands += "java -Dsun.java2d.noddraw=true -Xms2048M -Xmx2048M -classpath Haploview.jar -Djava.awt.headless=true edu.mit.wi.haploview.HaploView -nogui -pedfile "+trav+".ped -info "+trav+".info -tagrsqcutoff 0.9 -aggressiveTagging -aggressiveNumMarkers 3\n"+"mv "+trav+".ped.TESTS "+trav+"_tags23-0.9.txt\n";
					count++;
				}

				if (count>0) {
					commands += "rm "+trav+".ped.TAGS\n"+"jcpm gwa.FindMarkersNearGenes dir=./ file="+trav+".txt\n\n";
				}
				writer.print(commands);
			}
			writer.println();
			writer.close();
		} catch (IOException ioe) {
			System.err.println("Error writing fillIn");

		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = DEFAULT_DIR;
		String filename = DEFAULT_HIT_LIST;
		String root = DEFAULT_ROOT;
		double mafCutoff = DEFAULT_MAF_CUTOFF;
		int windowSize = DEFAULT_WINDOW_SIZE;
		boolean split = false;
		boolean createPlink = false;
		boolean fillIn = false;

		String usage = "\n"+"park.gwa.FindNearestGenes requires 0-1 arguments\n"+"   (1) directory (i.e. dir="+dir+" (default))\n"+"   (2) list of genes (i.e. file="+filename+" (default))\n"+"   (3) plink root for hapmap data (i.e. root="+root+" (default))\n"+"   (4) MAF cutoff (i.e. maf="+mafCutoff+" (default))\n"+"   (5) window size (i.e. win="+windowSize+" (default))\n"+"   (6) (optional) split genes (i.e. -splitGenes (not the default))\n"+"   (7) (optional) creates plink file if needed (i.e. -createPlink (not the default))\n"+"   (8) (optional) fill in missing tags (i.e. -fillIn (not the default))\n"+"";

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
			} else if (args[i].startsWith("root=")) {
				root = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("maf=")) {
				mafCutoff = Double.parseDouble(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("win=")) {
				windowSize = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("-splitGenes")) {
				split = true;
				numArgs--;
			} else if (args[i].startsWith("-createPlink")) {
				createPlink = true;
				numArgs--;
			} else if (args[i].startsWith("-fillIn")) {
				fillIn = true;
				numArgs--;
			} else {
				System.err.println("Error - '"+args[i]+"' is not a valid flag");
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (split) {
				splitGenesAndBatch();
			} else if (fillIn) {
				fillIn();
			} else if (createPlink) {
				createPlink(dir, filename, root);
			} else if (!new File(dir+ext.rootOf(filename)+"_qualified.txt").exists()||!new File(dir+ext.rootOf(filename)+"_tags.txt").exists()) {
				findMarkers(dir, filename, root, mafCutoff, windowSize);
			} else {
				System.out.println("Found qualified.txt and tags.txt; processing results");
				summarizeMarkers(dir, filename, mafCutoff, windowSize);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
