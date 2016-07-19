package cnv.park;

import java.io.*;
import java.util.*;

import cnv.var.*;
import common.*;
import filesys.CNVariant;

public class ComparePlinkResults {
//	public static final String DEFAULT_ROOT = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV\\allresults\\";
//	public static final String DEFAULT_ROOT = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV_PD\\results\\homozygosity\\";
	public static final String DEFAULT_ROOT = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV_PD\\results\\paperComp\\";
//	public static final String DEFAULT_ROOT = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV_PD\\results\\HLA\\";
//	public static final String DEFAULT_ROOT = "C:\\Documents and Settings\\npankrat\\My Documents\\Singleton\\results\\";
	
	
	// public static final String DEFAULT_ROOT = "C:\\Documents and Settings\\npankrat\\My Documents\\osteo\\comparisons\\";
	public static final String[] DEFAULT_DIRS = {"position\\", "window\\"};
	public static final String DEFAULT_CALLS_DIR = "calls\\";
	public static final String DEFAULT_GLOBAL_DIR = "global\\";
	public static final String[] INDIV_HEADER = {"FID", "IID", "PHE", "NSEG", "KB", "KBAVG"};
	public static final String[] SUMMARY_MPERM_HEADER = {"CHR", "SNP", "EMP1"};
	public static final String[] TESTS = {"Total number of CNVs", "Proportion with at least 1 CNV", "Summed length of CNVs", "Average length of CNVs"};
	public static final String[] TEST_SUFFIX = {" CNVs", "%", " kb", " kb"};

	public static final String[] MPERM_REQUIRED = {"SNP", "EMP2"};
//	public static final String[] MPERM_HEADER = {"CHR", "SNP", "STAT", "EMP1", "EMP2"};
	// public static final String[] MPERM_HEADER = {"CHR", "SNP", "BP", "EMP1", "EMP2"};
//	public static final int BLUR = 1000000;
	public static final int BLUR = 50000;

//	public static final double SIGNIFICANCE_THRESHOLD = 0.10;
	public static final double SIGNIFICANCE_THRESHOLD = 0.99999;
//	public static final double SIGNIFICANCE_THRESHOLD = 0.20;

	public static void compare(String rootDirectory, String[] dirs, String globalDirectory) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String begin, end;
		File[] files;
		Vector<int[]> hits, composite;
		int chr, start, stop, pos;
		int[] region, keys, indices;
		double[][] minPs;
		double[][][] sums;
		double[][] pvals;
		int aff;
		String match;

		for (int i = 0; i<dirs.length; i++) {
			files = new File(rootDirectory+dirs[i]).listFiles(new FilenameFilter() {
				public boolean accept(File file, String filename) {
					return filename.endsWith(".mperm");
				}
			});
			System.out.println("found "+files.length+" files in "+dirs[i]);

			hits = new Vector<int[]>();
			for (int j = 0; j<files.length; j++) {
				try {
					reader = new BufferedReader(new FileReader(files[j]));
					System.out.println(files[j].getName());
					indices = ext.indexFactors(MPERM_REQUIRED, reader.readLine().trim().split("[\\s]+"), false, true);
					begin = end = "";
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						if (Double.parseDouble(line[indices[1]])<SIGNIFICANCE_THRESHOLD) {
							if (begin.equals("")) {
								begin = line[indices[0]];
							}
							end = line[indices[0]];
						} else if (!begin.equals("")) {
							hits.add(new int[] {Integer.parseInt(begin.substring(1, begin.indexOf("-"))), Integer.parseInt(begin.substring(begin.indexOf("-")+1)), Integer.parseInt(end.substring(end.indexOf("-")+1))});
							begin = end = "";
						}
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""+files[j]+"\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \""+files[j]+"\"");
					System.exit(2);
				}
			}

			keys = Sort.orderTwoLayers(Matrix.toMatrix(hits), new Logger());

			composite = new Vector<int[]>();
			for (int j = 0; j<hits.size(); j++) {
				chr = hits.elementAt(keys[j])[0];
				start = hits.elementAt(keys[j])[1];
				stop = hits.elementAt(keys[j])[2];
				for (int l = 0; l<composite.size(); l++) {
					region = composite.elementAt(l);
					if (chr==region[0]&&(isBetween(start, region[1]-BLUR, region[2]+BLUR)||isBetween(stop, region[1]-BLUR, region[2]+BLUR))) {
						region[1] = Math.min(start, region[1]);
						region[2] = Math.max(stop, region[2]);
						hits.elementAt(keys[j])[0] = -1;
					}
				}
				if (hits.elementAt(keys[j])[0]!=-1) {
					composite.add(new int[] {chr, start, stop});
				}
			}

			minPs = Matrix.doubleMatrix(composite.size(), files.length, 1);
			for (int j = 0; j<files.length; j++) {
				try {
					reader = new BufferedReader(new FileReader(files[j]));
					indices = ext.indexFactors(MPERM_REQUIRED, reader.readLine().trim().split("[\\s]+"), false, true);
					begin = end = "";
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						chr = Integer.parseInt(line[indices[0]].substring(1, line[indices[0]].indexOf("-")));
						pos = Integer.parseInt(line[indices[0]].substring(line[indices[0]].indexOf("-")+1));
						for (int k = 0; k<composite.size(); k++) {
							region = composite.elementAt(k);
							if (chr==region[0]&&pos>=region[1]&&pos<=region[2]) {
								if (Double.parseDouble(line[indices[1]])<minPs[k][j]) {
									minPs[k][j] = Double.parseDouble(line[indices[1]]);
								}
							}
						}
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""+files[j]+"\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \""+files[j]+"\"");
					System.exit(2);
				}
			}

			try {
				writer = new PrintWriter(new FileWriter(rootDirectory+dirs[i].substring(0, dirs[i].length()-1)+"_comparison.xln"));
				writer.print("Chromosome\tStart\tStop\tUCSClink");
				for (int j = 0; j<files.length; j++) {
					writer.print("\t"+ext.rootOf(files[j].getName()));
				}
				writer.println();
				for (int j = 0; j<composite.size(); j++) {
					writer.print(Array.toStr(composite.elementAt(j))+"\t=HYPERLINK(\""+Positions.getUCSClink(composite.elementAt(j))+"\", \"link\")");
					for (int k = 0; k<files.length; k++) {
						writer.print("\t"+ext.prettyP(minPs[j][k], 2, 5, 1, false));
					}
					writer.println();
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to file");
				e.printStackTrace();
			}
			System.out.println();
		}

		try {
			files = new File(rootDirectory+globalDirectory).listFiles(new FilenameFilter() {
				public boolean accept(File file, String filename) {
					return filename.endsWith(".cnv.indiv");
				}
			});
			System.out.println("found "+files.length+" files in "+globalDirectory);

            sums = new double[files.length][2][5];
            pvals = new double[files.length][4];
			for (int i = 0; i<files.length; i++) {
				System.out.println(files[i].getName());
				try {
	                reader = new BufferedReader(new FileReader(files[i]));
	                ext.checkHeader(reader.readLine().trim().split("[\\s]+"), INDIV_HEADER, true);
	                while (reader.ready()) {
	                	line = reader.readLine().trim().split("[\\s]+");
	                	aff = Integer.parseInt(line[2])-1;
	                	sums[i][aff][0] += 1;
	                	sums[i][aff][1] += Integer.parseInt(line[3]);
	                	sums[i][aff][2] += Integer.parseInt(line[3])>0?100:0;
	                	sums[i][aff][3] += Double.parseDouble(line[4]);
	                	sums[i][aff][4] += Double.parseDouble(line[5]);
	                }
	                reader.close();
                } catch (FileNotFoundException fnfe) {
	                System.err.println("Error: file \""+files[i].getName()+"\" not found in current directory");
	                System.exit(1);
                } catch (IOException ioe) {
	                System.err.println("Error reading file \""+files[i].getName()+"\"");
	                System.exit(2);
                }
            	match = ext.replaceAllWith(files[i].getAbsolutePath(), ".cnv.indiv", ".cnv.summary.mperm");
            	try {
	                reader = new BufferedReader(new FileReader(match));
	                ext.checkHeader(reader.readLine().trim().split("[\\s]+"), SUMMARY_MPERM_HEADER, true);
	                for (int j = 0; j<4; j++) {
	                	line = reader.readLine().trim().split("[\\s]+");
	                	pvals[i][j] = Double.parseDouble(line[2]);
                    }
	                reader.close();
                } catch (FileNotFoundException fnfe) {
	                System.err.println("Error: file \""+match+"\" not found in current directory");
	                System.exit(1);
                } catch (IOException ioe) {
	                System.err.println("Error reading file \""+match+"\"");
	                System.exit(2);
                }
            }
			writer = new PrintWriter(new FileWriter(rootDirectory+globalDirectory.substring(0, globalDirectory.length()-1)+"_comparison.xln"));
			writer.print("Test\t");
			for (int i = 0; i<files.length; i++) {
				writer.print("\t"+ext.rootRootOf(ext.rootOf(files[i].getName(), true)));
			}
			writer.println();
			for (int cc = 1; cc>=0; cc--) {
				writer.print("\t"+(cc==1?"#cases":"#controls"));
				for (int i = 0; i<files.length; i++) {
					writer.print("\t"+(int)sums[i][cc][0]);
	            }
				writer.println();
            }
			for (int test = 0; test<TESTS.length; test++) {
				for (int cc = 1; cc>=0; cc--) {
					writer.print(cc==1?TESTS[test]+"\tcases":"\tcontrols");
					for (int i = 0; i<files.length; i++) {
						writer.print("\t"+ext.formDeci(sums[i][cc][test+1]/sums[i][cc][0], 2)+TEST_SUFFIX[test]);
		            }
					writer.println();
				}
				writer.print("\tp-value");
				for (int i = 0; i<files.length; i++) {
					writer.print("\t"+ext.prettyP(pvals[i][test]));
	            }
				writer.println();
            }
			
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+rootDirectory+globalDirectory);
	        e.printStackTrace();
        }
	}

	public static boolean isBetween(int value, int beginBoundary, int endBoundary) {
		return value>=beginBoundary&&value<=endBoundary;
	}
	
	public static void determineOverlapBetweenFilters(String rootDirectory, String callDirectory) {
		PrintWriter writer;
		Hashtable<String,boolean[]> hash;
		File[] files;
		String trav;
		CNVariant[] cnvs;
		boolean[] rollcall;
		String[] list;

		try {
			files = new File(rootDirectory+callDirectory).listFiles(new FilenameFilter() {
				public boolean accept(File file, String filename) {
					return filename.endsWith(".cnv");
				}
			});
			System.out.println("found "+files.length+" files in "+callDirectory);
		
			hash = new Hashtable<String,boolean[]>();
			for (int i = 0; i<files.length; i++) {
				System.out.println(files[i].getName());
				cnvs = CNVariant.loadPlinkFile(files[i].getAbsolutePath(), false);
				for (int j = 0; j<cnvs.length; j++) {
					trav = cnvs[j].getFingerprint();
					if (hash.containsKey(trav)) {
						rollcall = hash.get(trav);
					} else {
						hash.put(trav, rollcall = new boolean[files.length]);
					}
					rollcall[i] = true;
                }
            }
			
			list = HashVec.getKeys(hash);
			writer = new PrintWriter(new FileWriter(rootDirectory+callDirectory.substring(0, callDirectory.length()-1)+"_comparison.xln"));
			writer.print("CNV");
			for (int i = 0; i<files.length; i++) {
				writer.print("\t"+ext.rootOf(files[i].getName(), true));
            }
			writer.println();
			for (int i = 0; i<list.length; i++) {
				writer.println(list[i]+"\t"+Array.toStr(Array.booleanArrayToStringArray(hash.get(list[i]))));
            }
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+rootDirectory+callDirectory.substring(0, callDirectory.length()-1)+"_comparison.xln");
	        e.printStackTrace();
        }
		
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "ComparePlinkResults.dat";
		String rootDirectory = DEFAULT_ROOT;
		String[] dirs = DEFAULT_DIRS;
		String calls = DEFAULT_CALLS_DIR;
		String global = DEFAULT_GLOBAL_DIR;

		String usage = "\\n"+"park.cnv.ComparePlinkResults requires 0-1 arguments\n"+"   (1) filename (i.e. file="+filename+" (default))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (true) {
				compare(rootDirectory, dirs, global);
			}
			determineOverlapBetweenFilters(rootDirectory, calls);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
