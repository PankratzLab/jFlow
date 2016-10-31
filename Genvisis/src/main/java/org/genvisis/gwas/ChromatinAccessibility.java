package org.genvisis.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Matrix;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.filesys.SegmentLists;
import org.genvisis.stats.ProbDist;

public class ChromatinAccessibility {
	// public static final String SRC_DIR = "D:/grants/2012.06 Hepatoblastoma/Chromatin
	// accessibility/rawData/";
	// public static final String SRC_DIR = "D:/grants/2012.06 Hepatoblastoma/Chromatin
	// accessibility/rawData/heme/";
	// public static final String SRC_DIR = "D:/grants/2012.06 Hepatoblastoma/Chromatin
	// accessibility/rawData/cancers/";
	// public static final String SRC_DIR = "D:/grants/2012.06 Hepatoblastoma/Chromatin
	// accessibility/rawData/other/";
	public static final String SRC_DIR = "D:/grants/2012.06 Hepatoblastoma/Chromatin accessibility/rawData/all/";
	public static final String[][] PVAL_HEADER = {{"Marker", "SNP"}, {"Chr"}, {"Pos"},
																								{"pval", "p-val"}};
	public static final String[] CLASSES = {"All", "Valid", "DNaseOnly_OpenChrom",
																					"FAIREOnly_OpenChrom", "ChIPOnly"};

	public static String[] coverage(Segment[][] segs) {
		double sumIn, sumOut;
		int prevPos;

		sumOut = sumIn = 0;
		for (int chr = 0; chr <= 25; chr++) {
			prevPos = 0;
			for (int i = 0; i < segs[chr].length; i++) {
				sumOut += segs[chr][i].getStart() - prevPos;
				sumIn += segs[chr][i].getSize();
				prevPos = segs[chr][i].getStop();
			}
		}

		// return new String[] {ext.formPercent(sumIn/(sumIn+sumOut),2), sumIn+"", sumOut+""};
		return new String[] {ext.formPercent(sumIn / (sumIn + sumOut), 2)};
	}

	public static void stats(String dir) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, keys;
		String trav;
		Hashtable<String, Vector<double[]>> hash = new Hashtable<String, Vector<double[]>>();
		Vector<double[]> v;
		String[] files;
		double[][] matrix;

		files = Files.list(dir, ".bed", false);
		try {
			writer = new PrintWriter(new FileWriter("odd.xln"));
			for (String file : files) {
				try {
					reader = new BufferedReader(new FileReader(dir + file));
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						trav = line[3].substring(0, line[3].indexOf("_"));
						if (hash.containsKey(trav)) {
							v = hash.get(trav);
						} else {
							hash.put(trav, v = new Vector<double[]>());
						}
						v.add(Array.toDoubleArray(Array.subArray(line, 9)));
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \"" + dir + file + "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + dir + file + "\"");
					ioe.printStackTrace();
					System.exit(2);
				}

				System.out.println("Parsed: " + ext.rootOf(file));
				writer.println(ext.rootOf(file));
				keys = HashVec.getKeys(hash);
				for (int j = 0; j < keys.length; j++) {
					matrix = Matrix.toDoubleArrays(hash.get(keys[j]));
					writer.print("mean");
					for (int k = 0; k < matrix[j].length; k++) {
						writer.print("\t" + Array.mean(Matrix.extractColumn(matrix, k)));
					}
					writer.println("\t" + keys[j]);
					writer.print("min");
					for (int k = 0; k < matrix[j].length; k++) {
						writer.print("\t" + Array.min(Matrix.extractColumn(matrix, k)));
					}
					writer.println("\t" + keys[j]);
					writer.print("max");
					for (int k = 0; k < matrix[j].length; k++) {
						writer.print("\t" + Array.max(Matrix.extractColumn(matrix, k)));
					}
					writer.println("\t" + keys[j]);
				}
				writer.println();
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + "odd.xln");
			e.printStackTrace();
		}
	}

	public static Segment[][][][] loadSegments(String dir) {
		BufferedReader reader;
		Vector<Vector<Vector<Segment>>> allSeg;
		String[] line, files;
		Segment[][][][] segs;
		long time;
		Vector<Vector<Segment>> mSeg;
		Vector<Segment> vSeg;
		byte chr;

		files = Files.list(dir, ".bed", false);
		segs = new Segment[files.length + 1][CLASSES.length][26][];
		try {
			System.out.println("Tissue\tSerialized\ttimeToLoad\tClass\tProportion affected\tnumIn\tnumOut\tClass\tProportion affected\tnumIn\tnumOut");

			for (int i = 0; i < files.length; i++) {
				System.out.print(ext.rootOf(files[i]));
				if (Files.exists(dir + files[i] + "_" + CLASSES[0] + ".ser", false)) {
					time = new Date().getTime();
					for (int j = 0; j < CLASSES.length; j++) {
						segs[i][j] = SegmentLists	.load(dir + files[i] + "_" + CLASSES[j] + ".ser", false)
																			.getLists();
					}
					// System.out.println("Re-loaded serialized version of '"+ext.rootOf(files[i])+"' in " +
					// ext.getTimeElapsed(time));
					System.out.print("\t1\t" + ext.getTimeElapsed(time));
				} else {
					time = new Date().getTime();
					allSeg = new Vector<Vector<Vector<Segment>>>();
					for (String element : CLASSES) {
						allSeg.add(mSeg = new Vector<Vector<Segment>>());
						for (int k = 0; k <= 25; k++) {
							mSeg.add(new Vector<Segment>());
						}
					}

					try {
						// Segment.loadRegions(dir+files[i], 0, 1, 2);
						reader = new BufferedReader(new FileReader(dir + files[i]));
						System.out.println("Loading '" + files[i] + "'");
						while (reader.ready()) {
							line = reader.readLine().trim().split("[\\s]+");
							chr = Positions.chromosomeNumber(line[0]);
							allSeg.elementAt(0).elementAt(chr)
										.add(new Segment(chr, Integer.parseInt(line[1]), Integer.parseInt(line[2])));
							allSeg.elementAt(Integer.parseInt(line[20])).elementAt(chr)
										.add(new Segment(chr, Integer.parseInt(line[1]), Integer.parseInt(line[2])));
							for (int j = 2; line[20].equals("1") && j < CLASSES.length; j++) {
								allSeg.elementAt(j).elementAt(chr)
											.add(new Segment(chr, Integer.parseInt(line[1]), Integer.parseInt(line[2])));
							}
						}
						reader.close();
					} catch (FileNotFoundException fnfe) {
						System.err.println("Error: file \""	+ dir + files[i]
																+ "\" not found in current directory");
						System.exit(1);
					} catch (Exception ioe) {
						System.err.println("Error reading file \"" + dir + files[i] + "\"");
						ioe.printStackTrace();
						System.exit(2);
					}

					for (int j = 0; j < CLASSES.length; j++) {
						for (int k = 0; k <= 25; k++) {
							segs[i][j][k] = Segment.toArray(allSeg.elementAt(j).elementAt(k));
							Arrays.sort(segs[i][j][k]);
						}
						new SegmentLists(segs[i][j]).serialize(dir + files[i] + "_" + CLASSES[j] + ".ser");
					}
					// System.out.println("Parsed and serialized '"+ext.rootOf(files[i])+"' in " +
					// ext.getTimeElapsed(time));
					System.out.print("\t0\t" + ext.getTimeElapsed(time));
				}
				// System.out.println("\tClass\tProportion affected\tnumIn\tnumOut");
				for (int j = 0; j < 2; j++) {
					// for (int j = 0; j < CLASSES.length; j++) {
					System.out.print("\t" + CLASSES[j] + "\t" + Array.toStr(coverage(segs[i][j])));
				}
				System.out.println();
			}
		} catch (Exception e) {
			System.err.println("Error writing to " + "odd.xln");
			e.printStackTrace();
		}

		time = new Date().getTime();
		System.out.print("AllMerged");
		if (Files.exists(dir + "allFilesMerged_" + CLASSES[0] + ".ser", false)) {
			for (int j = 0; j < CLASSES.length; j++) {
				segs[files.length][j] = SegmentLists.load(dir	+ "allFilesMerged_" + CLASSES[j] + ".ser",
																									false)
																						.getLists();
			}
			// System.out.println("Re-loaded serialized version of the merged summaries in " +
			// ext.getTimeElapsed(time));
			System.out.print("\t1\t" + ext.getTimeElapsed(time));
		} else {
			// boolean hi = true;
			// if (hi) {
			// return segs;
			// }
			System.out.println("Merging segments across " + files.length + " files");
			for (int j = 0; j < CLASSES.length; j++) {
				for (int k = 0; k <= 25; k++) {
					System.out.print(".");
					vSeg = new Vector<Segment>();
					for (int i = 0; i < files.length; i++) {
						for (int n = 0; n < segs[i][j][k].length; n++) {
							vSeg.add(segs[i][j][k][n]);
						}
					}
					Segment.mergeOverlapsAndSort(vSeg);
					segs[files.length][j][k] = Segment.toArray(vSeg);
					Arrays.sort(segs[files.length][j][k]);
				}
				System.out.println();
				new SegmentLists(segs[files.length][j]).serialize(dir	+ "allFilesMerged_" + CLASSES[j]
																													+ ".ser");
			}
			// System.out.println("Merged all files serialized the summaries in " +
			// ext.getTimeElapsed(time));
			System.out.print("\t0\t" + ext.getTimeElapsed(time));
		}
		// System.out.println("\tClass\tProportion affected\tnumIn\tnumOut");
		for (int j = 0; j < 2; j++) {
			// for (int j = 0; j < CLASSES.length; j++) {
			System.out.print("\t" + CLASSES[j] + "\t" + Array.toStr(coverage(segs[files.length][j])));
		}
		System.out.println();

		return segs;
	}


	public static void filter(String bed_dir, String filename, String levels) {
		BufferedReader reader;
		PrintWriter[][][] writers;
		String[] line;
		String temp;
		int count;
		int[] indices;
		Segment[][][][] segs;
		Segment seg;
		String file_dir;
		String[] files;
		byte chr;
		boolean in;

		files = Files.list(bed_dir, ".bed", false);
		segs = loadSegments(bed_dir);

		file_dir = ext.parseDirectoryOfFile(filename);
		filename = ext.removeDirectoryInfo(filename);
		try {
			writers = new PrintWriter[files.length + 1][CLASSES.length][2];
			try {
				for (int i = 0; i < files.length + 1; i++) {
					for (int j = 0; j < CLASSES.length; j++) {
						new File(file_dir + ext.rootOf(filename) + "/").mkdirs();
						if (i == files.length) {
							writers[i][j][0] = new PrintWriter(new FileWriter(file_dir	+ ext.rootOf(filename)
																																+ "/"
																																+ ext.addToRoot(filename,
																																								"_allMerged_" + CLASSES[j] + "_in")));
							writers[i][j][1] = new PrintWriter(new FileWriter(file_dir	+ ext.rootOf(filename)
																																+ "/"
																																+ ext.addToRoot(filename,
																																								"_allMerged_" + CLASSES[j] + "_out")));
						} else {
							writers[i][j][0] = new PrintWriter(new FileWriter(file_dir	+ ext.rootOf(filename)
																																+ "/"
																																+ ext.addToRoot(filename,
																																								"_"				+ ext.rootOf(files[i])
																																													+ "_"
																																													+ CLASSES[j]
																																													+ "_in")));
							writers[i][j][1] = new PrintWriter(new FileWriter(file_dir	+ ext.rootOf(filename)
																																+ "/"
																																+ ext.addToRoot(filename,
																																								"_"				+ ext.rootOf(files[i])
																																													+ "_"
																																													+ CLASSES[j]
																																													+ "_out")));
						}
					}
				}
			} catch (Exception e) {
				System.err.println("Error - setting up writers");
			}
			try {
				reader = new BufferedReader(new FileReader(file_dir + filename));
				indices = ext.indexFactors(	PVAL_HEADER, reader.readLine().trim().split("[\\s]+"), false,
																		false, true, true);
				count = 0;
				while (reader.ready()) {
					count++;
					if (count % 10000 == 0) {
						System.out.print(".");
						// System.out.println(count);
					}
					temp = reader.readLine();
					line = temp.trim().split("[\\s]+");
					chr = Positions.chromosomeNumber(line[indices[1]]);
					seg = new Segment(chr, Integer.parseInt(line[indices[2]]),
														Integer.parseInt(line[indices[2]]));
					for (int i = 0; i < files.length + 1; i++) {
						for (int j = 0; j < CLASSES.length; j++) {
							if (i < files.length) {
								// in = Segment.binarySearch(seg, segs[i][j][chr]) != -1;
								in = Segment.binarySearchForOverlap(seg, segs[i][j][chr]) != -1; // this code was
																																									// altered when
																																									// merging with
																																									// other code
							} else {
								in = false;
								for (int i2 = 0; i2 < files.length; i2++) {
									if (Segment.binarySearchForOverlap(seg, segs[i2][j][chr]) != -1) { // this code
																																											// was altered
																																											// when
																																											// merging
																																											// with other
																																											// code
										in = true;
									}
								}
							}

							if (in) {
								writers[i][j][0].println(temp);
							} else {
								writers[i][j][1].println(temp);
							}
						}
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + filename + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + filename + "\"");
				ioe.printStackTrace();
				System.exit(2);
			}
			System.out.println();
			System.out.println();

			System.out.println("Tissue\tnumInAll\tnumOutAll\t%InAll\tnumInValid\tnumOutValid\t%InValid");
			for (int i = 0; i < files.length + 1; i++) {
				System.out.print(i == files.length ? "AllMerged" : ext.rootOf(files[i]));
				for (int j = 0; j < CLASSES.length; j++) {
					writers[i][j][0].close();
					writers[i][j][1].close();

					if (j < 2) {
						double[] pvals;

						double numIn = Files.countLines(file_dir	+ ext.rootOf(filename) + "/"
																						+ ext.addToRoot(filename,
																														(i == files.length	? "_allMerged_"
																																								: "_"
																																										+ ext.rootOf(files[i])
																																									+ "_")
																																			+ CLASSES[j] + "_in"),
																						0);
						// System.out.print("\t"+numIn);
						double numOut = Files.countLines(file_dir	+ ext.rootOf(filename) + "/"
																							+ ext.addToRoot(filename,
																															(i == files.length	? "_allMerged_"
																																									: "_"
																																											+ ext.rootOf(files[i])
																																										+ "_")
																																				+ CLASSES[j] + "_out"),
																							0);
						// System.out.print("\t"+numOut+"\t"+ext.formPercent(numIn/(numIn+numOut), 2));
						System.out.print("\t" + ext.formPercent(numIn / (numIn + numOut), 2));

						pvals = Array.toDoubleArray(HashVec.loadFileToStringArray(file_dir
																																				+ ext.rootOf(filename) + "/"
																																			+ ext.addToRoot(filename,
																																											(i == files.length	? "_allMerged_"
																																																					: "_"
																																																							+ ext.rootOf(files[i])
																																																						+ "_")
																																																	+ CLASSES[j]
																																																+ "_in"),
																																			false, new int[] {3}, false));
						System.out.print("\t" + ext.formDeci(
																									ProbDist.ChiDistReverse(Array.median(Array.removeNaN(pvals)),
																																					1)
																									/ ProbDist.ChiDistReverse(0.50, 1), 4));
						System.out.print("\t" + ext.formDeci(
																									ProbDist.ChiDistReverse(Array.quantExclusive(	Array.removeNaN(pvals),
																																																0.10),
																																					1)
																									/ ProbDist.ChiDistReverse(0.10, 1), 4));
						pvals = Array.toDoubleArray(HashVec.loadFileToStringArray(file_dir
																																				+ ext.rootOf(filename) + "/"
																																			+ ext.addToRoot(filename,
																																											(i == files.length	? "_allMerged_"
																																																					: "_"
																																																							+ ext.rootOf(files[i])
																																																						+ "_")
																																																	+ CLASSES[j]
																																																+ "_out"),
																																			false, new int[] {3}, false));
						System.out.print("\t" + ext.formDeci(
																									ProbDist.ChiDistReverse(Array.median(Array.removeNaN(pvals)),
																																					1)
																									/ ProbDist.ChiDistReverse(0.50, 1), 4));
						System.out.print("\t" + ext.formDeci(
																									ProbDist.ChiDistReverse(Array.quantExclusive(	Array.removeNaN(pvals),
																																																0.10),
																																					1)
																									/ ProbDist.ChiDistReverse(0.10, 1), 4));
					}
				}
				System.out.println();
			}
		} catch (Exception e) {
			System.err.println("Error filtering " + filename);
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String bed_dir = SRC_DIR;
		String filename = "ChromatinAccessibility.dat";
		String levels = "1,2,3,4";

		String usage = "\n"	+ "gwas.ChromatinAccessibility requires 0-1 arguments\n"
										+ "   (1) directory with bed files (i.e. dir=" + bed_dir + " (default))\n"
										+ "   (2) levels to use (i.e. levels=" + levels + " (default))\n"
										+ "   (3) filename with chr/pos/pvals (i.e. file=" + levels + " (default))\n"
										+ "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("dir=")) {
				bed_dir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("file=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("levels=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			// stats(SRC_DIR);
			filename = "D:/CARe/conditionalMeta/normalized/slimResultsSorted.dat";
			filter(bed_dir, filename, levels);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
