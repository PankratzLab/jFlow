package cnv.park;

import java.io.*;
import java.util.*;

import cnv.var.CNVariant;

import common.*;

public class DNAcopy {
//	public static final String ROOT_DIR = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV\\Software\\DNAcopy\\";
	public static final String[] ROOT_DIRS = {"C:\\Documents and Settings\\npankrat\\My Documents\\CNV\\dnacopy\\", "/work/parkinsons/cnvs/dnacopy/", "/home/npankrat/penncnv/again/dnacopy/"};
//	public static final String CNV_DIRECTORY = "../quantisnp/source/";
	public static final String CNV_DIRECTORY = "../cnvs/";
	public static final String DEFAULT_OUTPUT = "results/";
	public static final String DEFAULT_BATCH = "batch/";
	public static final String LOOKUP_PLUS_GENDER_FILE = "LookupPlusGender.xln";
	public static final double ABSOLUT_THRESHOLD = 1.5;
	public static final String[] LOOKUP_HEADER = {"Sample_ID", "FID", "IID", "Gender"};
	public static final String[] RESULTS_HEADER = {"ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean"};

	public static void createBatch(String inputDirectory, String batchDirectory, String outputDirectory, int numBatches) {
		PrintWriter writer;
		Vector<String[]> v = new Vector<String[]>();
		String[] inputs, outputs;
		String trav, commands;

		new File(batchDirectory).mkdirs();
		new File(outputDirectory).mkdirs();

		inputs = new File(inputDirectory).list(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				// return filename.endsWith(".qs");
				return !filename.contains(".");
			}
		});

		outputs = new File(outputDirectory).list(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(".segments");
			}
		});

		System.out.println("Found "+inputs.length+" samples, as well as results for "+outputs.length+" that have been done (not necessarily the same ones)");

		for (int i = 0; i<inputs.length; i++) {
			// trav = ext.rootOf(inputs[i]);
			trav = inputs[i];
			if (ext.indexOfStr(trav+"_output.out", outputs)==-1) {
				try {
					writer = new PrintWriter(new FileWriter(batchDirectory+trav+".batch"));
					writer.println("library(DNAcopy)");
					// writer.println("quant<-read.table(\""+inputDirectory+trav+".qs\",
					// header=TRUE, sep=\"\\t\")");
					writer.println("quant<-read.table(\""+inputDirectory+trav+"\", header=TRUE, sep=\"\\t\")");
					writer.println("cna <- CNA(cbind(quant$"+rName(trav)+".Log.R.Ratio), quant$Chr, quant$Position, data.type = \"logratio\", sampleid = \""+trav+"\")");
					writer.println("smoothed <- smooth.CNA(cna)");
					writer.println("write.table(smoothed, \""+outputDirectory+trav+".smooth\")");
					writer.println("segmented <- segment(smoothed, verbose = 1)");
					writer.println("sink(\""+outputDirectory+trav+".segments\")");
					writer.println("segmented");
					writer.println("sink()");
					writer.println("quit(\"no\")");
					writer.close();
				} catch (Exception e) {
					System.err.println("Error writing to file "+batchDirectory+trav+".batch");
					e.printStackTrace();
				}
				v.add(new String[] {trav});
			}
		}

		System.out.println("Made "+numBatches+" batch files that will take care of the "+v.size()+" files yet to analyze");

		commands = "R CMD BATCH "+batchDirectory+"[%0].batch\n";
		Files.batchIt("batch", null, numBatches, commands, Matrix.toStringArrays(v));
	}

	public static String rName(String dna) {
		return (dna.charAt(0)>='0'&&dna.charAt(0)<='9'?"X":"")+dna;
	}

	public static void parseResults(String rootDirectory, String outputDirectory, String included, double threshold) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Hashtable<String,String> lookupDNAtoSubject = new Hashtable<String,String>();
		Hashtable<String,String> lookupSubjectToDNA = new Hashtable<String,String>();
		String[] outputs, inds;
		int[][] counts;
		int[] array;

		try {
			reader = new BufferedReader(new FileReader(rootDirectory+LOOKUP_PLUS_GENDER_FILE));
			ext.checkHeader(reader.readLine().trim().split("[\\s]+"), LOOKUP_HEADER, true);
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				lookupDNAtoSubject.put(line[0], line[1]+"\t"+line[2]);
				lookupSubjectToDNA.put(line[1]+"\t"+line[2], line[0]);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+rootDirectory+LOOKUP_PLUS_GENDER_FILE+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+rootDirectory+LOOKUP_PLUS_GENDER_FILE+"\"");
			System.exit(2);
		}

		outputs = new File(rootDirectory+outputDirectory).list(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(".segments");
			}
		});

		if (included==null) {
			System.out.println("Using all indivduals with result files (n="+outputs.length+")");
			inds = new String[outputs.length];
			for (int i = 0; i<outputs.length; i++) {
				inds[i] = ext.rootOf(outputs[i]);
			}
		} else {
			inds = Array.toStringArray(HashVec.loadFileToVec(rootDirectory+included, false, false, true, true));
			System.out.println("Using all indivduals listed in '"+included+"' (n="+inds.length+")");
			for (int i = 0; i<inds.length; i++) {
				if (ext.indexOfStr(inds[i]+".segments", outputs)==-1) {
					System.err.println("Error - '"+inds[i]+"' was found in the list of indiviudals to parse, but no results were found for this sample");
				}
			}
		}

		try {
			writer = new PrintWriter(new FileWriter(rootDirectory+"dnacopy.cnv"));
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
			counts = new int[inds.length][2];
			for (int i = 0; i<inds.length; i++) {
				try {
					reader = new BufferedReader(new FileReader(rootDirectory+outputDirectory+inds[i]+".segments"));
					reader.readLine();
					reader.readLine();
					reader.readLine();
					ext.checkHeader(reader.readLine().trim().split("[\\s]+"), RESULTS_HEADER, true);

					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						if (!line[1].equals(rName(inds[i]))) {
							System.err.println("Error - sample mismatch for "+inds[i]);
							break;
						}
						counts[i][0]++;
						if (Math.abs(Double.parseDouble(line[6]))>threshold) {
							writer.println(lookupDNAtoSubject.get(inds[i])+"\t"+Positions.chromosomeNumber(line[2])+"\t"+line[3]+"\t"+line[4]+"\t"+(Double.parseDouble(line[6])<0?1:3)+"\t"+line[6]+"\t0");
							counts[i][1]++;
						}
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""+rootDirectory+outputDirectory+inds[i]+".segments"+"\" not found in directory... won't be in final plink.cnv file");
				} catch (IOException ioe) {
					System.err.println("Error reading file \""+rootDirectory+outputDirectory+inds[i]+".segments"+"\"");
					System.exit(2);
				}
				// writer.println(inds[i]);
			}
			array = Matrix.extractColumn(counts, 0);
			System.out.println("Mean/Total number of CNVs per person is: "+Array.mean(array)+" / "+Array.sum(array));
			array = Matrix.extractColumn(counts, 1);
			System.out.println("Mean/Total number of CNVs per person exceeding the threshold is: "+Array.mean(array)+" / "+Array.sum(array));
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("Error writing to file");

		}
	}

	public static void analyzeDistribution(String rootDirectory, String outputDirectory, String included, double min, double max, int sigfigs) {
		BufferedReader reader;
		PrintWriter writer;
		String[] outputs, inds, line;
		double d, obsMin, obsMax;
		boolean faulty = false;
		int[] counts;
		String output;

		outputs = new File(rootDirectory+outputDirectory).list(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(".smooth");
			}
		});

		if (included==null||!new File(included).exists()) {
			if (included!=null) {
				System.err.println("Error - file '"+included+"' was not found in directory "+rootDirectory);
			}
			System.out.println("Using all indivduals with result files (n="+outputs.length+")");
			inds = new String[outputs.length];
			for (int i = 0; i<outputs.length; i++) {
				inds[i] = ext.rootOf(outputs[i]);
			}
			output = "Distribution_for_ALL."+sigfigs+"sigfig"+".xln";
		} else {
			inds = Array.toStringArray(HashVec.loadFileToVec(included, false, false, true, true));
			System.out.println("Using all dnas listed in '"+included+"' (n="+inds.length+")");
			for (int i = 0; i<inds.length; i++) {
				if (ext.indexOfStr(inds[i]+".smooth", outputs)==-1) {
					System.err.println("Error - '"+inds[i]+"' was found in the list of indiviudals to parse, but no results were found for this sample");
				}
			}
			output = "Distribution_for_"+ext.rootOf(included)+"."+sigfigs+"sigfig"+".xln";
		}

		obsMin = Double.POSITIVE_INFINITY;
		obsMax = Double.NEGATIVE_INFINITY;
		counts = new int[(int)((max-min)*Math.pow(10, sigfigs))];
		for (int i = 0; i<inds.length; i++) {
			System.out.println((i+1)+" of "+inds.length+": "+inds[i]);
			try {
				reader = new BufferedReader(new FileReader(rootDirectory+outputDirectory+inds[i]+".smooth"));
				reader.readLine();
				while (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
					// if (!line[1].contains("X") && !line[1].contains("Y")) {
					if (!line[1].contains("X")&&!line[1].contains("Y")&&!line[1].equals("\"6\"")) {
						d = Double.parseDouble(line[3]);
						if (d<obsMin) {
							obsMin = d;
							System.out.println("New observed minimum: "+d);
						}
						if (d>obsMax) {
							obsMax = d;
							System.out.println("New observed maximum: "+d);
						}
						if (!faulty) {
							if (d<min) {
								System.err.println("Error - '"+d+"' exceeds minimum ("+min+"); continuing on to determine absolute maximum and minimum");
								faulty = true;
							} else if (d>max) {
								System.err.println("Error - '"+d+"' exceeds minimum ("+min+"); continuing on to determine absolute maximum and minimum");
								faulty = true;
							} else {
								counts[(int)(d*Math.pow(10, sigfigs)-min*Math.pow(10, sigfigs))]++;
							}
						}
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \""+rootDirectory+outputDirectory+inds[i]+".smooth"+"\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \""+rootDirectory+outputDirectory+inds[i]+".smooth"+"\"");
				System.exit(2);
			}
		}
		try {
			writer = new PrintWriter(new FileWriter(output));
			for (int i = 0; i<counts.length; i++) {
				writer.println(ext.formDeci(min+Math.pow(0.1, sigfigs)*i, sigfigs)+"\t"+ext.formDeci(min+Math.pow(0.1, sigfigs)*(i+1), sigfigs)+"\t"+counts[i]);
			}
			writer.close();

		} catch (Exception e) {
			e.printStackTrace();
		}

		System.out.println("Observed minimum: "+obsMin);
		System.out.println("Observed maximum: "+obsMax);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String root = Files.firstDirectoryThatExists(ROOT_DIRS, true, true);
		String inputs = CNV_DIRECTORY;
		String batches = DEFAULT_BATCH;
		String outputs = DEFAULT_OUTPUT;
		String include = null;
		int numBatches = 6;
		boolean createFiles = false;
		boolean createBatch = false;
		boolean parseResults = false;
		boolean analyzeDis = false;
		double min = -15;
		double max = 5;
		int sigfigs = 2;
		double threshold = ABSOLUT_THRESHOLD;

		String usage = "\\n"+"park.cnv.QuantiSNP requires 0-1 arguments\n"+"   (1) root directory of CNV files (i.e. root="+root+" (default))\n"+"   (2) -createFiles (not the default)\n"+"   (3) number of batches to create (i.e. batches="+numBatches+" (not the default))\n"+"   (4) -parseResults (not the default)\n"+"   (5) file with list of inidiuvdals to parse (i.e. include=keeps.txt (default is to include all with results))\n"+"   (6) -dist (analyze the distribution of smoothed LRR values) (not the default)\n"+"   (7) number of significant digits to histogram (i.e. sig="+sigfigs+" (default))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("root=")) {
				root = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("batches=")) {
				numBatches = Integer.parseInt(args[i].split("=")[1]);
				createBatch = true;
				numArgs--;
			} else if (args[i].startsWith("-createFiles")) {
				createFiles = true;
				numArgs--;
			} else if (args[i].startsWith("-parseResults")) {
				parseResults = true;
				numArgs--;
			} else if (args[i].startsWith("-dist")) {
				analyzeDis = true;
				numArgs--;
			} else if (args[i].startsWith("include=")) {
				include = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("sig=")) {
				sigfigs = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("thresh=")) {
				threshold = Double.parseDouble(args[i].split("=")[1]);
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (createFiles) {
				QuantiSNP.createFiles();
			}
			if (createBatch) {
				createBatch(inputs, batches, outputs, numBatches);
			}
			if (parseResults) {
				parseResults(root, outputs, include, threshold);
			}
			if (analyzeDis) {
				analyzeDistribution(root, outputs, include, min, max, sigfigs);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
