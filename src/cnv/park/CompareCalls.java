package cnv.park;

import java.io.*;
import java.util.*;

import common.*;
import filesys.CNVariant;
import stats.Maths;

public class CompareCalls {
	public static final String DEFAULT_ROOT = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV\\allCalls\\";

	// public static final String[] DEFAULT_FILES = {"conf.cnv", "allMarkers.cnv"};
	// public static final String[] DEFAULT_FILES = {"conf_100kb_5SNP_10.0.cnv", "allMarkers_100kb_5SNP_10.0.cnv"};
	public static final String[] DEFAULT_FILES = {"conf_100kb_5SNP_10.0.cnv", "conf_100kb_20SNP_10.0.cnv"};

	public static void compare(String rootDir, String[] files) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, inds;
		Hashtable<String,Hashtable<String,Vector<CNVariant>>> hash = new Hashtable<String,Hashtable<String,Vector<CNVariant>>>();
		Hashtable<String,Vector<CNVariant>> source = new Hashtable<String,Vector<CNVariant>>();
		CNVariant[][] cnvs;
		Vector<CNVariant> v = new Vector<CNVariant>();
		int match;
		int[] counts;
		int[][] allPossibleCombinations = Maths.getIndicesForAllCombinations(files.length, 2);

		for (int i = 0; i<files.length; i++) {
			try {
				reader = new BufferedReader(new FileReader(rootDir+files[i]));
				if (!ext.checkHeader(reader.readLine().trim().split("[\\s]+"), CNVariant.PLINK_CNV_HEADER, false)) {
					reader.close();
					return;
				}
				while (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
					if (hash.containsKey(line[0]+"\t"+line[1])) {
						source = hash.get(line[0]+"\t"+line[1]);
					} else {
						hash.put(line[0]+"\t"+line[1], source = new Hashtable<String,Vector<CNVariant>>());
					}
					if (source.containsKey(i+"")) {
						v = source.get(i+"");
					} else {
						source.put(i+"", v = new Vector<CNVariant>());
					}
					v.add(new CNVariant(line, i));
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \""+rootDir+files[i]+"\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \""+rootDir+files[i]+"\"");
				System.exit(2);
			}
		}

		inds = HashVec.getKeys(hash);
		for (int i = 0; i<allPossibleCombinations.length; i++) {
			try {
				writer = new PrintWriter(new FileWriter(rootDir+"Compare "+ext.rootOf(files[allPossibleCombinations[i][0]])+" and "+ext.rootOf(files[allPossibleCombinations[i][1]])+".xln"));
				writer.println("FID\tIID\tTotal"+ext.rootOf(files[allPossibleCombinations[i][0]])+"\tTotal"+ext.rootOf(files[allPossibleCombinations[i][1]])+"\tUnique"+ext.rootOf(files[allPossibleCombinations[i][0]])+"\tUnique"+ext.rootOf(files[allPossibleCombinations[i][1]])+"\tOverlaps\tExactMatches");
				for (int j = 0; j<inds.length; j++) {
					cnvs = new CNVariant[][] {CNVariant.toCNVariantArray(hash.get(inds[j]).get(allPossibleCombinations[i][0]+"")), CNVariant.toCNVariantArray(hash.get(inds[j]).get(allPossibleCombinations[i][1]+""))};
					counts = new int[4];
					if (cnvs[0].length==0) {
						System.err.println("Error - "+inds[j]+" not found in "+files[allPossibleCombinations[i][0]]);
					}
					if (cnvs[1].length==0) {
						System.err.println("Error - "+inds[j]+" not found in "+files[allPossibleCombinations[i][1]]);
					}
					for (int a = 0; a<cnvs[0].length; a++) {
						match = 0;
						for (int b = 0; b<cnvs[1].length; b++) {
							if (cnvs[0][a].equals(cnvs[1][b])) {
								match = 3;
								cnvs[1][b].setSource(99);
							} else if (match<2&&cnvs[0][a].overlaps(cnvs[1][b])) {
								match = 2;
								cnvs[1][b].setSource(99);
							}
						}
						counts[match]++;
					}
					for (int b = 0; b<cnvs[1].length; b++) {
						match = 1;
						for (int a = 0; a<cnvs[0].length; a++) {
							if (cnvs[1][b].getSource()!=99&&cnvs[1][b].equals(cnvs[0][a])) {
								match = 3;
							} else if (match<2&&cnvs[1][b].getSource()!=99&&cnvs[1][b].overlaps(cnvs[0][a])) {
								match = 2;
							}
						}
						if (cnvs[1][b].getSource()!=99) {
							counts[match]++;
						}
					}
					writer.println(inds[j]+"\t"+cnvs[0].length+"\t"+cnvs[1].length+"\t"+Array.toStr(counts));
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error comparing "+files[allPossibleCombinations[i][0]]+" and "+files[allPossibleCombinations[i][1]]);
				e.printStackTrace();
			}
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String rootDirectory = DEFAULT_ROOT;
		String[] files = DEFAULT_FILES;

		String usage = "\\n"+"park.cnv.ComparePlinkResults requires 0-1 arguments\n"+"   (1) directory (i.e. dir="+rootDirectory+" (default))\n"+"   (2) files to be compared (i.e. files="+Array.toStr(files, ",")+" (default))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				rootDirectory = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("files=")) {
				files = args[i].split("=")[1].split(",");
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			compare(rootDirectory, files);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
