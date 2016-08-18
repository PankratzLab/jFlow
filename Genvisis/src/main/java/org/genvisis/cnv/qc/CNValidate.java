package org.genvisis.cnv.qc;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.manage.Transforms;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

public class CNValidate implements Runnable {
	private static final double SCALE_FACTOR_MAD = 0.134894516;
	private Project proj;
	private String[] inds;
	private Hashtable<String, CNVariantQC[]> allIndcnVariantQCs;
	private SampleData sampleData;
	private MarkerSet markerSet;
	private CNVariantQC[][] indcnVariantQCs;
	private Logger log;

	public CNValidate(Project proj, String[] inds, Hashtable<String, CNVariantQC[]> allIndcnVariantQCs, MarkerSet markerSet) {
		this.proj = proj;
		this.inds = inds;
		this.allIndcnVariantQCs = allIndcnVariantQCs;
		this.sampleData = proj.getSampleData(2, false);
		this.markerSet = markerSet;
		this.indcnVariantQCs = new CNVariantQC[inds.length][];
		this.log = proj.getLog();
	}

	public void run() {

		for (int i = 0; i < inds.length; i++) {
			log.report(ext.getTime() + "\t" + (i + 1) + " of " + inds.length);
			CNVariantQC[] cnVariantQCs = allIndcnVariantQCs.get(inds[i]);
			String[] ids = sampleData.lookup(cnVariantQCs[0].getCnVariant().getFamilyID() + "\t" + cnVariantQCs[0].getCnVariant().getIndividualID());
			if (ids != null) {
				Sample samp = proj.getFullSampleFromRandomAccessFile(ids[0]);
				if (samp != null) {
					log.report(ext.getTime() + " Loaded Sample..." + samp.getSampleName());
					if (markerSet.getFingerprint() != samp.getFingerprint()) {
						log.reportError("Error - sample and marker fingerprints do not match");
						System.exit(1);
					} else {
						log.report(ext.getTime() + " Computing validations for sample..." + samp.getSampleName());
						indcnVariantQCs[i] = computeValidations(samp, markerSet, cnVariantQCs, log);
						log.report(ext.getTime() + " Finshed computing validations for sample..." + samp.getSampleName());
					}
				}
			} else {
				log.reportError("Error - Sample data for " + inds[i] + " not found in samples directory");
				System.exit(1);
			}
		}
	}

	public static CNVariantQC[][] computeMultiThreadedValidations(Project proj, String[] inds, Hashtable<String, CNVariantQC[]> allIndcnVariantQCs, MarkerSet markerSet, int processors) {
		if (processors == 0) {
			processors = Runtime.getRuntime().availableProcessors();
		}
		Thread[] threads = new Thread[processors];
		Vector<Vector<String>> cabinet = getcabinet(inds, processors);
		CNValidate[] cnvals = processValidations(proj, processors, threads, cabinet, allIndcnVariantQCs, markerSet);
		return collectAllValidations(processors, cnvals, inds, proj.getLog());

	}

	public Hashtable<String, CNVariantQC[]> getAllIndcnVariantQCs() {
		return allIndcnVariantQCs;
	}

	public CNVariantQC[][] getIndcnVariantQCs() {
		return indcnVariantQCs;
	}

	public String[] getInds() {
		return inds;
	}

	private static CNVariantQC[] computeValidations(Sample samp, MarkerSet markerSet, CNVariantQC[] cnVariantQCs, Logger log) {
		int[][] indices = markerSet.getIndicesByChr();
		float[] LRRsInvTransformedByChr = Transforms.transform(samp.getLRRs(), 2, true, markerSet);
		double[] chrLRRMediansMADScaled = getChrLRRMediansMADScaled(indices, LRRsInvTransformedByChr, SCALE_FACTOR_MAD, log);
		double[] LRRsInvTransformedByChrMADScaled = scaleMAD(indices, LRRsInvTransformedByChr, chrLRRMediansMADScaled);
		return evaluateCNVariantQCs(samp, markerSet, cnVariantQCs, LRRsInvTransformedByChrMADScaled, indices, log);
	}

	// the main event
	private static CNVariantQC[] evaluateCNVariantQCs(Sample samp, MarkerSet markerSet, CNVariantQC[] cnVariantQCs, double[] LRRsInvTransformedByChrMADScaled, int[][] indices, Logger log) {
		String[] markerNames = markerSet.getMarkerNames();
		float[] bafs = samp.getBAFs();
		byte[] abGenotypes = samp.getAB_Genotypes();
		double callRate =getCallRate(markerNames , abGenotypes);
		for (int i = 0; i < cnVariantQCs.length; i++) {
			Hashtable<String, Integer> markersIncnVariant = cnVariantQCs[i].getMarkersIncnVariant();
			String[] markerList = cnVariantQCs[i].getMarkerList();
			double[] variantLRRs = new double[markerList.length];
			double[] variantbafs = new double[markerList.length];
			byte[] variantGenotypes = new byte[markerList.length];
			for (int k = 0; k < markerList.length; k++) {
				int markerindex = markersIncnVariant.get(markerList[k]);
				if (markerNames[markerindex].equals(markerList[k])) {
					if (Double.isNaN(LRRsInvTransformedByChrMADScaled[markerindex])) {
						log.reportError("Error - the cnv " + cnVariantQCs[i].getCnVariant().toPlinkFormat() + " contained a NaN LRR value at marker " + markerList[k] + " the height will be set to zero. QC metrics may be inaccurate");
						variantLRRs[k] = 0;
					} else {
						variantLRRs[k] = LRRsInvTransformedByChrMADScaled[markerindex];
						variantbafs[k] = bafs[markerindex];
						variantGenotypes[k] = abGenotypes[markerindex];
					}
				} else {
					log.reportError("Error - Received unmatched indices for marker " + markerList[k] + ", got " + markerNames[markerindex] + "this should not happen");
					System.exit(1);
				}
			}
			if (markerList.length == cnVariantQCs[i].getCnVariant().getNumMarkers() && markerList.length == variantLRRs.length) {
				cnVariantQCs[i].setHeight(Array.median(variantLRRs));
				cnVariantQCs[i].setBafs(variantbafs);
				cnVariantQCs[i].setGenotypes(variantGenotypes);
				cnVariantQCs[i].setLrrs(variantLRRs);
				cnVariantQCs[i].setSampleCallRate(callRate);
				cnVariantQCs[i].setSourceFile(samp.getSampleName());
			} else {
				log.reportError("Error - there were less markers contained in the cnv region " + cnVariantQCs[i].getCnVariant().toPlinkFormat() + "  than markers in the position file ");
				System.exit(1);
			}
		}
		return cnVariantQCs;
	}

	// TODO
	// Affy AND Illumina specific, need to determine SNP markers
	private static double getCallRate(String[] markerNames, byte[] abGenotypes) {
		double calls = 0;
		double snpMarkers = 0;
		for (int i = 0; i < markerNames.length; i++) {
			if (!markerNames[i].startsWith("CN_") || markerNames[i].startsWith("cnv")) {
				snpMarkers++;
				if (abGenotypes[i] != -1) {
					calls++;
				}
			}
		}
		return (calls / snpMarkers);
	}

	private static double[] scaleMAD(int[][] indices, float[] LRRsInvTransformedByChr, double[] chrLRRMediansMADScaled) {
		double[] LRRsInvTransformedByChrMADScaled = new double[LRRsInvTransformedByChr.length];
		for (int i = 0; i < indices.length; i++) {
			if (indices[i].length > 0) {
				for (int j = 0; j < indices[i].length; j++) {
					LRRsInvTransformedByChrMADScaled[indices[i][j]] = (double) (LRRsInvTransformedByChr[indices[i][j]] / (chrLRRMediansMADScaled[i]));
				}
			}
		}
		return LRRsInvTransformedByChrMADScaled;
	}

	private static double[] getChrLRRMediansMADScaled(int[][] indices, float[] LRRsInvTransformedByChr, double SCALE_FACTOR_MAD, Logger log) {
		double[] chrLRRMediansMADScaled = new double[indices.length];
		for (int i = 0; i < indices.length; i++) {
			if (indices[i].length > 0) {
				ArrayList<Float> chrLRRal = new ArrayList<Float>();
				for (int j = 0; j < indices[i].length; j++) {
					// check for NaN
					if (LRRsInvTransformedByChr[j] == LRRsInvTransformedByChr[j]) {
						chrLRRal.add(Math.abs(LRRsInvTransformedByChr[j]));
					}
				}
				chrLRRMediansMADScaled[i] = (double) (median(getFloatArray(chrLRRal)) / SCALE_FACTOR_MAD);
			} else {
				log.report("Warning - not analyzing chromomosome " + i + " , did not find any markers");
				continue;
			}
		}
		return chrLRRMediansMADScaled;
	}

	private static float[] getFloatArray(ArrayList<Float> al) {
		float[] array = new float[al.size()];
		for (int i = 0; i < al.size(); i++) {
			array[i] = al.get(i);
		}
		return array;
	}

	private static CNValidate[] processValidations(Project proj, int processors, Thread[] threads, Vector<Vector<String>> cabinet, Hashtable<String, CNVariantQC[]> allIndcnVariantQCs, MarkerSet markerSet) {
		CNValidate[] cnvals = new CNValidate[processors];
		for (int i = 0; i < processors; i++) {
			cnvals[i] = new CNValidate(proj, cabinet.elementAt(i).toArray(new String[cabinet.elementAt(i).size()]), allIndcnVariantQCs, markerSet);
			threads[i] = new Thread(cnvals[i]);
			threads[i].start();
		}
		checkThreadStatus(processors, threads);
		return cnvals;
	}

	private static CNVariantQC[][] collectAllValidations(int processors, CNValidate[] cnvals, String[] inds, Logger log) {
		CNVariantQC[][] cnVariantQCs = new CNVariantQC[inds.length][];
		int indIndex = 0;
		int counter = 0;
		for (int i = 0; i < inds.length; i++) {
			counter++;
			if (counter > processors) {
				indIndex += 1;
				counter = 1;
			}
			if (cnvals[i % processors].getInds()[indIndex].equals(inds[i])) {
				if (cnvals[i % processors].getIndcnVariantQCs()[indIndex].length == cnvals[i % processors].getAllIndcnVariantQCs().get(inds[i]).length) {
					cnVariantQCs[i] = cnvals[i % processors].getIndcnVariantQCs()[indIndex];
				} else {
					log.reportError("Error - recieved unmatched cnv numbers while collecting results for  " + inds[i] + ": " + cnvals[i % processors].getIndcnVariantQCs()[indIndex].length + " and " + cnvals[i % processors].getAllIndcnVariantQCs().get(inds[i]).length);
					System.exit(1);
				}
			} else {
				log.reportError("Error - recieved unmatched ids while collecting results for " + cnvals[i % processors].getInds()[indIndex] + "\t" + inds[i]);
				System.exit(1);
			}
		}
		return cnVariantQCs;
	}

	private static Vector<Vector<String>> getcabinet(String[] inds, int processors) {
		Vector<Vector<String>> cabinet = new Vector<Vector<String>>();
		for (int i = 0; i < processors; i++) {
			cabinet.add(new Vector<String>());
		}
		for (int i = 0; i < inds.length; i++) {
			cabinet.elementAt(i % processors).add(inds[i]);
		}
		return cabinet;
	}

	private static float median(float[] array) {
		return (quant(array, (float) 0.50));
	}

	private static float quant(float[] array, float q) {
		int keys[] = Sort.quicksort(array);
		try {
			if (q > 1 || q < 0) {
				return (0);
			} else {
				double index = (array.length + 1) * q;
				if (index - (int) index == 0) {
					return array[keys[(int) index - 1]];
				} else {
					return q * array[keys[(int) Math.floor(index) - 1]] + (1 - q) * array[keys[(int) Math.ceil(index) - 1]];
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			return -1234567890;
		}
	}

	private static void checkThreadStatus(int processors, Thread[] threads) {
		boolean complete;
		complete = false;
		while (!complete) {
			complete = true;
			for (int i = 0; i < processors; i++) {
				if (threads[i].isAlive()) {
					complete = false;
				}
			}
			if (!complete) {
				try {
					Thread.sleep(1000L);
				} catch (InterruptedException ex) {
				}
			}
		}
	}
}

// TODO for testing purposes , delete later
// if (samp.getSampleName().matches("BURDS_p_ARIC_batch4_011_affy_GenomeWideSNP_6_G05_213086.CEL.IND.txt")) {
// try {
// int nanacount = 0;
// PrintWriter cnvOutput = new PrintWriter(new FileWriter("D:/data/compare/BURDS_p_ARIC_batch4_011_affy_GenomeWideSNP_6_G05_213086.CEL.IND.LRRE.txt", false));
// System.err.println("HIDDFMGenomeWideSNP_6_F04_213068");
// cnvOutput.println("Name\tBURDS_p_ARIC_batch4_011_affy_GenomeWideSNP_6_G05_213086.CEL.Log R Ratio");
// for (int i = 0; i < markerSet.getMarkerNames().length; i++) {
// if ((samp.getLRRs()[i] + "").equals("NaN")) {
// nanacount++;
// }
// cnvOutput.println(markerSet.getMarkerNames()[i] + "\t" + samp.getLRRs()[i]);
//
// }
// System.err.println(nanacount + " NANs");
// cnvOutput.close();
// } catch (IOException e) {
// System.err.println("Could Not Open");
// e.printStackTrace();
// System.exit(1);
// }
//
// }

// 6 min_gap
// 30 max_gap
// 0.25 ratio_min
// 50000 block_size
// 3 number of delimiter characters before location column
// 8 number of delimiter characters before logratio column
// datafile.dat
// 1 indicator for whether file is sorted, 1 if sorted 0 otherwise
// 0.4 exponent for score
// 46 ascii value of the delimiter character
// 1 flag for whether to quantile normalize data
// 0 output file cnv.pred.res with predicted values at each probe if 1"
// summary.res name of summary output file
// pred.res name of predicted values file

// public static void validate(String plinkCnvs, Project proj) {
// int processors = Runtime.getRuntime().availableProcessors();
// Thread[] threads = new Thread[processors];
// CNVariant[] fileCNVs = CNVariant.loadPlinkFile(plinkCnvs, false);
// String[] inds = CompareCalls_dev.toStringArray(CompareCalls_dev.getIDList(fileCNVs));
// Hashtable<String, ArrayList<CNVariant>> allIndCNVs = cnv.park.CompareCalls_dev.getindCNVsfromFile(fileCNVs);
// Vector<Vector<String>> cabinet = getcabinet(inds, processors);
//
// startJobs(proj, processors, threads, cabinet, allIndCNVs);
// checkThreadStatus(processors, threads);
//
// }

// public static void main(String[] args) {
// int numArgs = args.length;
// // String filename = cnv.Launch.getDefaultDebugProjectFile();
// String filename = "C:/workspace/Genvisis/projects/inSilicoValidate.properties";
// String plinkCnvs = "C:/data/inSilicoValidate/Filtered_LRR_35.35_conf_0.0_numMarkers_0_RootsANDDups.cnv";
//
// String loggerFilename = "C:/data/inSilicoValidate/SureIllLog.txt";
// Logger log;
// String usage = "\n" + "cnv.qc.CNValidate requires 0-1 arguments\n" + "   (0) project properties filename (i.e. proj=" + filename + " (default))\n" + "   (1) plink CNV format .cnv file to validate(i.e. cnvs=all_cnvs.cnv\n" + "   (2) clean the .cnv file before validation (i.e -clean\n" + "";
// // String usage = "\n"+
// // "cnv.park.PennCNV requires 0-1 arguments\n"+
// // "   (0) project properties filename (i.e. proj="+filename+" (default))\n"+
// // "   (1) plink CNV format .cnv file to validate(i.e. cnvs=all_cnvs.cnv\n"+
// // "   (2) clean the .cnv file before validation (i.e -clean)\n"+
// // "   (3) filter the .cnv file before validation (i.e -filter)\n"+
// // "";
//
// for (int i = 0; i < args.length; i++) {
// if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
// System.err.println(usage);
// System.exit(1);
// } else if (args[i].startsWith("proj=")) {
// filename = args[i].split("=")[1];
// numArgs--;
// } else if (args[i].startsWith("cnvs=")) {
// plinkCnvs = args[i].split("=")[1];
// numArgs--;
// } else if (args[i].startsWith("-clean")) {
// clean = true;
// numArgs--;
// } else if (args[i].startsWith("-filter")) {
// filter = true;
// numArgs--;
// }
// }
// if (numArgs != 0) {
// System.err.println(usage);
// System.exit(1);
// }
// proj = new Project(filename, false);
// log = new Logger(loggerFilename);
//
// try {
//
// proj = new Project(filename, false);
//
// log = new Logger(loggerFilename);
//
// if (plinkCnvs != null) {
// // validate(plinkCnvs, proj);
// }
//
// } catch (Exception e) {
// e.printStackTrace();
// }
// }
// }

// try {
// // PrintWriter writercnv = Files.getAppropriateWriter("BMI_NEGR1_BEASTYSummary");
// BufferedReader reader = Files.getAppropriateReader(cnvFile + ".summary.res");
// PrintWriter writer = new PrintWriter(new FileWriter("D:/data/inSilicoValidate/NEGR1_all_BEAST.txt", true));
// try {
// reader.readLine();
// Vector<Double> thierTransforms = new Vector<Double>();
// if (!reader.ready()) {
// continue;
// }
// while (reader.ready()) {
// String[] line = reader.readLine().trim().split("\\s+");
// thierTransforms.add(Double.parseDouble(line[0]));
//
// }
// double[] thierTransformsdubs = new double[thierTransforms.size()];
// for (int j = 1; j < thierTransforms.size(); j++) {
// thierTransformsdubs[j] = Math.abs(thierTransforms.get(j));
// }
// System.out.println(cnvFile + ".pred.res\tmedian\t" + Array.median(thierTransformsdubs));
//
// } catch (IOException e) {
// // TODO Auto-generated catch block
// e.printStackTrace();
// }
//
// } catch (FileNotFoundException e) {
// // TODO Auto-generated catch block
// e.printStackTrace();
// }
// // System.out.println("BEASTIFIED " + cnvFile);
// catch (IOException e1) {
// // TODO Auto-generated catch block
// e1.printStackTrace();
// }
//
// }

// private static void unleashTheBeast(Sample samp, MarkerSet markerSet) {
// int[][] indices = markerSet.getIndicesByChr();
// int[] positions = markerSet.getPositions();
// byte[] chr = markerSet.getChrs();
// String[] markerNames = markerSet.getMarkerNames();
// float[] LRRs = samp.getLRRs();
//
// for (int i = 1; i < indices.length; i++) {
// if (i > 10) {
// // System.out.println(i);
// continue;
// } else {
// String cnvFile = "D:/data/inSilicoValidate/UnleashTheBEAST/" + samp.getSampleName().replaceAll(".*affy", "") + "_chr" + i + ".dat";
// String configFile = "D:/data/inSilicoValidate/UnleashTheBEAST/" + samp.getSampleName().replaceAll(".*affy", "") + "_chr" + i + ".config";
// PrintWriter writercnv = Files.getAppropriateWriter(cnvFile);
// writercnv.println("probeset_id\tchr\tposition\tLRR");
// PrintWriter writerconfig = Files.getAppropriateWriter(configFile);
// writerconfig.println("6\n30\n0.25\n50000\n2\n3\n" + cnvFile + "\n1\n0.4\n9\n1\n1\n" + cnvFile + ".summary.res\n" + cnvFile + ".pred.res\n");
// for (int k = 0; k < indices[i].length; k++) {
// writercnv.println(markerNames[indices[i][k]] + "\t" + i + "\t" + positions[indices[i][k]] + "\t" + LRRs[indices[i][k]]);
// // System.out.println(markerNames[indices[i][k]] + "\t" + i + "\t" + positions[indices[i][k]] + "\t" + LRRs[indices[i][k]]);
//
// }
// writercnv.close();
// writerconfig.close();
//
// CmdLine.run("D:/data/inSilicoValidate/cnv.beast.exe " + configFile, "");
//
// // // PrintWriter writercnv = Files.getAppropriateWriter("BMI_NEGR1_BEASTYSummary");
// //
//
// try {
// // PrintWriter writercnv = Files.getAppropriateWriter("BMI_NEGR1_BEASTYSummary");
// BufferedReader reader = Files.getAppropriateReader(cnvFile + ".pred.res");
// PrintWriter writer = new PrintWriter(new FileWriter("D:/data/inSilicoValidate/NEGR1_all_BEAST.txt", true));
// try {
// reader.readLine();
// Vector<Double> thierTransforms = new Vector<Double>();
// if (!reader.ready()) {
// continue;
// }
// while (reader.ready()) {
// String[] line = reader.readLine().trim().split("\\s+");
// thierTransforms.add(Double.parseDouble(line[0]));
//
// }
// double[] thierTransformsdubs = new double[thierTransforms.size()];
// for (int j = 1; j < thierTransforms.size(); j++) {
// thierTransformsdubs[j] = Math.abs(thierTransforms.get(j));
// }
// System.out.println(cnvFile + ".pred.res\tmedian\t" + Array.median(thierTransformsdubs));
//
// } catch (IOException e) {
// // TODO Auto-generated catch block
// e.printStackTrace();
// }
//
// } catch (FileNotFoundException e) {
// // TODO Auto-generated catch block
// e.printStackTrace();
// }
// // System.out.println("BEASTIFIED " + cnvFile);
// catch (IOException e1) {
// // TODO Auto-generated catch block
// e1.printStackTrace();
// }
// PrintWriter writer = new PrintWriter(new FileWriter("D:/data/inSilicoValidate/NEGR1_all_BEAST.txt", true));
// BufferedReader reader = Files.getAppropriateReader(cnvFile + ".summary.res");
// String testline = reader.readLine();
// if (testline == null) {
// System.err.println("there is a blank line on   " + cnvFile);
// break;
// }
// while (reader.ready()) {
// String[] line = reader.readLine().trim().split("\\s+");
// if (Integer.parseInt(line[0]) < 68000000 || Integer.parseInt(line[1]) > 77000000) {
// continue;
// }
// System.out.println(Array.toStr(line));
// writer.println(samp.getSampleName() + "\t" + Array.toStr(line));
//
// }
// writer.close();
//
// } catch (IOException e) {
// // TODO Auto-generated catch block
// e.printStackTrace();
// }
// }
// }
//
// }

// private static void evaluateCNVs(Sample samp, MarkerSet markerSet, CNVariant[] cnvs, double[] LRRsInvTransformedByChrMADScaled, int[][] indices, double alpha, float[] LRRsInvTransformedByChr, double[] chrLRRMediansMADScaled) {
// int[] positions = markerSet.getPositions();
// CNVariant[] beastConfcnvs = new CNVariant[cnvs.length];
// String[] markerNames = markerSet.getMarkerNames();
// System.out.println(markerNames.length);
// String cnvFile = "D:/data/inSilicoValidate/Roots_AND_Dups_BEAST_ALL.cnv";
//
// try {
// PrintWriter writer = new PrintWriter(new FileWriter(cnvFile, true));
// for (int i = 0; i < cnvs.length; i++) {
//
// double score = 0;
// double[] cnvLRRs = new double[cnvs[i].getNumMarkers()];
// int index = 0;
// for (int k = 0; k < indices[cnvs[i].getChr()].length; k++) {
// // System.out.println(cnvs[i].getChr());
// int position = positions[indices[cnvs[i].getChr()][k]];
// if (position >= cnvs[i].getStart() && position <= cnvs[i].getStop()) {
// cnvLRRs[index] = LRRsInvTransformedByChrMADScaled[indices[cnvs[i].getChr()][k]];
// // writer.println((index) + "\t" + cnvs[i].toPlinkFormat() + "\t" + cnvLRRs.length + "\t" + markerNames[indices[cnvs[i].getChr()][k]] + "\t" + position + "\t" + LRRsInvTransformedByChrMADScaled[indices[cnvs[i].getChr()][k]] + "\t" + LRRsInvTransformedByChr[indices[cnvs[i].getChr()][k]] + "\t" + chrLRRMediansMADScaled[cnvs[i].getChr()]);
// index++;
// }
// }
//
// // the main event
//
// if (index == cnvs[i].getNumMarkers() && index == cnvLRRs.length) {
// score = score(cnvs[i], alpha, cnvLRRs);
// writer.println(cnvs[i].toPlinkFormat() + "\t" + cnvLRRs.length + "\t" + Array.median(cnvLRRs) + "\t" + score);
// } else {
// System.err.println("Error - mismatched number of markers contained in cnv " + cnvs[i].toPlinkFormat() + " and markers in the position file ");
// }
// }
//
// writer.close();
// } catch (FileNotFoundException fnfe) {
// System.err.println("Error: file \"" + cnvFile + "\" could not be written to (it's probably open)");
// System.exit(1);
// } catch (IOException ioe) {
// System.err.println("Error reading file \"" + cnvFile + "\"");
// System.exit(2);
// }
//
// }

// private static void startJobs(Project proj, int processors, Thread[] threads, Vector<Vector<String>> cabinet, Hashtable<String, CNVariant[]> allIndCNVsArrays) {
// CNValidate cnval;
// // cnval = new CNValidate(proj, cabinet.elementAt(i).toArray(new String[cabinet.elementAt(i).size()]), allIndCNVs);
// for (int i = 0; i < processors; i++) {
//
// cnval = new CNValidate(proj, cabinet.elementAt(i).toArray(new String[cabinet.elementAt(i).size()]), allIndCNVsArrays);
//
// // threads[i] = new Thread();
// threads[i] = new Thread(cnval);
// threads[i].start();
// }
// }

// private static void unleashTheBeast(Sample samp, MarkerSet markerSet, SampleData sampleData) {
// int[][] indices = markerSet.getIndicesByChr();
// int[] positions = markerSet.getPositions();
//
// String[] markerNames = markerSet.getMarkerNames();
// float[] LRRs = samp.getLRRs();
//
// for (int i = 1; i < indices.length; i++) {
// String cnvFile = "D:/data/inSilicoValidate/UnleashTheBEAST/" + samp.getSampleName().replaceAll(".*affy", "") + "_chr" + i + ".dat";
// String configFile = "D:/data/inSilicoValidate/UnleashTheBEAST/" + samp.getSampleName().replaceAll(".*affy", "") + "_chr" + i + ".config";
// try {
//
// PrintWriter writercnv = new PrintWriter(new FileWriter(cnvFile, false));
// PrintWriter writerconfig = new PrintWriter(new FileWriter(configFile, false));
// writercnv.println("probeset_id\tchr\tposition\tLRR");
// writerconfig.println("6\n30\n0.25\n5000000\n2\n3\n" + cnvFile + "\n1\n0.4\n9\n1\n1\n" + cnvFile + ".summary.res\n" + cnvFile + ".pred.res\n");
// Hashtable<Integer, Float> tracker = new Hashtable<Integer, Float>();
//
// for (int k = 0; k < indices[i].length; k++) {
// if (tracker.containsKey(positions[indices[i][k]]) && (!markerNames[indices[i][k]].startsWith("AFFX") || !markerNames[indices[i][k]].startsWith("CN_"))) {
// tracker.put(positions[indices[i][k]], LRRs[indices[i][k]]);
// } else {
// tracker.put(positions[indices[i][k]], LRRs[indices[i][k]]);
// }
// }
//
// for (int k = 0; k < indices[i].length; k++) {
// if (tracker.containsKey(positions[indices[i][k]]) && tracker.get(positions[indices[i][k]]) == LRRs[indices[i][k]]) {
// writercnv.println(markerNames[indices[i][k]] + "\t" + i + "\t" + positions[indices[i][k]] + "\t" + LRRs[indices[i][k]]);
// tracker.remove(positions[indices[i][k]]);
// } else {
// // System.err.println("Skipping this marker since it was on a duplicate location " + markerNames[indices[i][k]] + "\t" + i + "\t" + positions[indices[i][k]] + "\t" + LRRs[indices[i][k]]);
// }
// }
//
// writercnv.close();
// writerconfig.close();
//
// CmdLine.run("D:/data/inSilicoValidate/cnv.beast.exe " + configFile, "");
// Thread.sleep(10);
//
// } catch (IOException e1) {
// System.err.println("Error - could not create file " + cnvFile);
// e1.printStackTrace();
// } catch (InterruptedException e2) {
// }
//
// // System.out.println("BEASTIFIED " + cnvFile);
// // // PrintWriter writercnv = Files.getAppropriateWriter("BMI_NEGR1_BEASTYSummary");
// //
// tameTheBeast(samp, sampleData, i, cnvFile);
// }
//
// }
// private static void tameTheBeast(Sample samp, SampleData sampleData, int chr, String cnvFile) {
// try {
// // PrintWriter writercnv = Files.getAppropriateWriter("BMI_NEGR1_BEASTYSummary");
// BufferedReader reader = Files.getAppropriateReader(cnvFile + ".summary.res");
// int numTests = 20;
// PrintWriter[] writer = new PrintWriter[numTests];
// double constant = 10;
// for (int k = 0; k < numTests; k++) {
// double alpha = (double) ((k / constant));
// // System.out.println(alpha);
// String file = "D:/data/inSilicoValidate/BEAST_alpha_" + alpha + ".cnv";
// if (Files.exists(file)) {
// writer[k] = new PrintWriter(new FileWriter(file, true));
// } else {
// writer[k] = new PrintWriter(new FileWriter(file, false));
// writer[k].println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
// }
// }
//
// // skip beast header
// reader.readLine();
// if (!reader.ready()) {
// System.err.println("Warning - no Beast calls were found in " + cnvFile + ".summary.res");
// return;
// }
// while (reader.ready()) {
// String[] line = reader.readLine().trim().split("\\s+");
// int type = 2;
// // convert beast calls to 1=deletion and 3 =duplication, removing any strange calls
// // with mixed types (Beast type!= 1)
// if (line[4].equals("1")) {
// if (Double.parseDouble(line[2]) > 0) {
// type = 3;
// } else {
// type = 1;
// }
// // System.out.println(line[4] + "\t" + type + "\t" + line[2]);
// for (int k = 0; k < numTests; k++) {
// double alpha = (double) ((k / constant));
// double score = Math.abs(Math.pow(Double.parseDouble(line[3]), alpha) * Double.parseDouble(line[2]));
// writer[k].println(sampleData.lookup(samp.getSampleName())[1] + "\t" + chr + "\t" + line[0] + "\t" + line[1] + "\t" + type + "\t" + score + "\t" + line[3]);
// }
// }
// }
// for (int k = 0; k < numTests; k++) {
// writer[k].close();
// }
//
// } catch (FileNotFoundException e) {
//
// e.printStackTrace();
// }
//
// catch (IOException e1) {
//
// e1.printStackTrace();
// }
// }