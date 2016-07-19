package org.genvisis.link;

import java.io.*;
import java.util.*;

import org.genvisis.common.ext;

public class gatherQuant {
	// public static String[] GH_FILE_TYPES = {"em-he", "ml_ndv", "ml",
	// "trad-he", "vc"};
	public static String[] GH_FILE_TYPES = {"vc"};

	// public static int[] GH_COLUMN = {3, 2, 2, 3, 2};
	public static int[] GH_COLUMN = {2};

	public static double LOD_THRESHOLD = 2.0;

	public gatherQuant(String trait, String[] dirs) throws IOException {
		BufferedReader reader = null;
		PrintWriter writer = null;
		PrintWriter peaks = null;
		PrintWriter lot = null;
		String[] line, vLine;
		Vector<String[]> vLines = new Vector<String[]>();
		String chrome;
		int numDirs, pos, max_pos;
		double lod, max_lod;
		boolean[] solar;
		String[] missing;
		File f;

		numDirs = dirs.length;
		solar = new boolean[numDirs];
		missing = new String[numDirs];
		for (int i = 0; i<numDirs; i++) {
			if (!new File(dirs[i]).exists()) {
				System.err.println("Error - "+dirs[i]+" does not exist.");
				System.exit(1);
			}
			if (!dirs[i].endsWith("/")) {
				dirs[i] += "/";
			}
			for (int j = 1; j<23; j++) {
				if (new File(dirs[i]+"chrom"+(j<10?"0"+j:""+j)).isDirectory()) {
					solar[i] = true;
				}
			}
			missing[i] = "Missing chrs ";
		}

		new File(trait+"_summary").mkdir();
		peaks = new PrintWriter(new FileWriter(trait+"_summary/linkage peaks.xls"));
		lot = new PrintWriter(new FileWriter(trait+"_summary/all chromosomes.xls"));
		for (int i = 0; i<numDirs; i++) {
			peaks.print("\t"+dirs[i]+"\tcM");
			lot.print("\t"+dirs[i]);
		}
		peaks.println();
		lot.println();
		for (int chr = 1; chr<=23; chr++) {
			chrome = chr<10?"0"+chr:""+chr;
			peaks.print(chr);
			writer = new PrintWriter(new FileWriter(trait+"_summary/summary"+chrome+".xls"));
			vLines.removeAllElements();
			for (int i = 0; i<numDirs; i++) {
				if (solar[i]) {
					max_lod = -999;
					max_pos = -1;
					if (chr<23) { // solar don't play X
						f = new File(dirs[i]+"chrom"+chrome+"/"+trait+"/multipoint1.out");
						if (!f.exists()||f.length()<1000) {
							missing[i] += " "+chrome;
						} else {
							missing[i] += "   ";
							reader = new BufferedReader(new FileReader(f));
							reader.readLine();
							reader.readLine();
							while (reader.ready()) {
								line = reader.readLine().split("[\\s]+");
								pos = Integer.valueOf(line[3]).intValue();
								if (ext.isMissingValue(line[4])) {
									lod = -998;
								} else {
									lod = Double.valueOf(line[4]).doubleValue();
								}
								if (lod>max_lod) {
									max_lod = lod;
									max_pos = pos;
								}

								if (vLines.size()<pos+1) {
									vLine = new String[numDirs+1];
									vLine[0] = pos+"";
									vLine[runIndex(i, 0, solar)] = lod+"";
									vLines.add(vLine);
								} else {
									vLine = vLines.elementAt(pos);
									vLine[runIndex(i, 0, solar)] = lod+"";
								}
							}
							reader.close();
						}
					}

					if (max_lod>LOD_THRESHOLD) {
						peaks.print("\t"+max_lod+"\t"+max_pos);
					} else {
						peaks.print("\t"+max_lod+"\t");
					}
				} else {
					try {
						for (int j = 0; j<GH_FILE_TYPES.length; j++) {
							f = new File(dirs[i]+GH_FILE_TYPES[j]+chrome+".out");
							if (!f.exists()||f.length()==0) {
								missing[i] += " "+chrome;
							} else {
								missing[i] += "   ";
								reader = new BufferedReader(new FileReader(f));
								reader.readLine();

								max_lod = -999;
								max_pos = -1;
								while (reader.ready()) {
									line = reader.readLine().split("[\\s]+");
									if (line.length<2) {
										line = reader.readLine().split("[\\s]+");
									}
									if (line[0].startsWith("==")||line[0].startsWith("Parameter")) {
										break;
									}
									pos = (int)Double.valueOf(line[0]).doubleValue();
									lod = Double.valueOf(line[GH_COLUMN[j]-1]).doubleValue();
									if (lod>max_lod) {
										max_lod = lod;
										max_pos = pos;
									}
									if (vLines.size()<pos+1) {
										vLine = new String[numDirs+1];
										vLine[0] = pos+"";
										vLine[runIndex(i, j, solar)] = lod+"";
										vLines.add(vLine);
									} else {
										vLine = vLines.elementAt(pos);
										vLine[runIndex(i, j, solar)] = lod+"";
									}
								}

								if (max_lod>LOD_THRESHOLD) {
									peaks.print("\t"+max_lod+"\t"+max_pos);
								} else {
									peaks.print("\t"+max_lod+"\t");
								}
								reader.close();
							}
						}
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
			}
			peaks.println();
			for (int i = 0; i<numDirs; i++) {
				if (solar[i]) {
					writer.print("\t"+dirs[i]);
				} else {
					for (int j = 0; j<GH_FILE_TYPES.length; j++) {
						writer.print("\t"+GH_FILE_TYPES[j]+"_"+dirs[i]);
					}
				}
			}
			writer.println();
			for (int i = 0; i<vLines.size(); i++) {
				vLine = vLines.elementAt(i);
				lot.print(i==0?chr:"");
				for (int j = 0; j<vLine.length; j++) {
					writer.print((j==0?"":"\t")+vLine[j]);
					lot.print((j==0?"":"\t")+vLine[j]);
				}
				writer.println();
				lot.println();
			}
			writer.close();
		}
		for (int i = 0; i<numDirs; i++) {
			System.err.println(missing[i]+" for "+dirs[i]);
		}
		peaks.close();
		lot.close();
	}

	public int runIndex(int dir, int type, boolean[] solar) {
		int count = 1;

		for (int i = 0; i<solar.length; i++) {
			for (int j = 0; j<(solar[i]?1:GH_FILE_TYPES.length); j++) {
				if (dir==i&&type==j) {
					return count;
				}
				count++;
			}
		}

		return count;
	}

	public static void main(String[] args) throws IOException {
		String trait = "AOO";
		Vector<String> directories = new Vector<String>();
		String[] dirs;

		String usage = "\n"+"park.gatherQuant requires at least 1 argument\n"+"   names of directories to be parsed\n"+"   (optional) trait (i.e. 'trait="+trait+"' (default))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("trait=")) {
				trait = args[i].split("=")[1];
			} else {
				directories.add(args[i]);
				if (!new File(args[i]).isDirectory()) {
					System.err.println("Error - directory '"+args[i]+"' does not exist");
					System.exit(2);
				}
			}
		}
		dirs = new String[directories.size()];
		for (int i = 0; i<dirs.length; i++) {
			dirs[i] = directories.elementAt(i);
		}
		if (dirs.length==0) {
			System.err.println(usage);
			System.exit(1);
		}

		try {
			new gatherQuant(trait, dirs);
			// new gatherQuant("AOO", new String[] {"111+/", "111-/", "100+/",
			// "100-/"});
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
