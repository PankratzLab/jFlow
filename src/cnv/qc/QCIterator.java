package cnv.qc;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;

import stats.Maths;
import common.Array;
import common.Files;
import common.Logger;
import common.ext;
import cnv.filesys.CNVQC;
import cnv.filesys.MarkerFreqs;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.var.CNVariant;

public class QCIterator implements Runnable {
	private Hashtable<String, CNVSampleQC> cnvSampleQCHash;
	private CNVariantQC[][][] cnvQCsAssigned;
	private double[] targetPercentages;
	private OptimizedQCThresholds[][] optqcs;
	private int optimizationType;
	private Logger log;

	public QCIterator(CNVariantQC[][][] cnvQCsAssigned, Hashtable<String, CNVSampleQC> cnvSampleQCHash, double[] targetPercentages, int optimizationType, Logger log) {
		this.cnvQCsAssigned = cnvQCsAssigned;
		this.cnvSampleQCHash = cnvSampleQCHash;
		this.targetPercentages = targetPercentages;
		this.optqcs = new OptimizedQCThresholds[targetPercentages.length][];
		this.optimizationType = optimizationType;
		this.log = log;
	}

	public void run() {
		for (int i = 0; i < targetPercentages.length; i++) {
			log.report(ext.getTime() + " Beginning iterations for target percentage " + targetPercentages[i]);
			optqcs[i] = iterate(cnvSampleQCHash, cnvQCsAssigned, targetPercentages[i], optimizationType, log);
			log.report(ext.getTime() + " Finished iterations for target percentage " + targetPercentages[i]);
		}
	}

	// preps CNV file for comparision, loads CNV heights, and sends jobs out for parameter iteration at each percent target
	public static void optimizeQCThresholds(Project proj, String plinkCnvQCs, String duplicatesFile, String SampleQCFile, String output, Logger log, int optimizationType) {
		int processors = Runtime.getRuntime().availableProcessors();
		double[] allPercents = binIt(0, 1, 15);
		getcabinet(allPercents, processors);
		Thread[] threads = new Thread[processors];
		log.report(ext.getTime() + " Prepping cnvQCs in " + plinkCnvQCs + " for comparisons");
		Hashtable<String, Hashtable<String, Integer>> defineCompHash = defineCompLists(proj.getProjectDir(), duplicatesFile, log);
		CNVariantQC[][][] cnvQCsAssigned = CNVariantQC.prepCNVQCsForComparison(proj, plinkCnvQCs, defineCompHash, log);
		log.report(ext.getTime() + " Finished prepping cnvQCs in " + plinkCnvQCs + " for " + cnvQCsAssigned[0].length + " comparisons");
		log.report(ext.getTime() + " Beginning iterations for " + allPercents.length + " target percentages");
		Hashtable<String, CNVSampleQC> cnvSampleQCHash = CNVSampleQC.getSampleQCs(proj, SampleQCFile, log);
		QCIterator[] qcits = iteratePercents(proj, processors, threads, getcabinet(allPercents, processors), cnvQCsAssigned, cnvSampleQCHash, optimizationType, log);
		log.report(ext.getTime() + " Finished iterations for " + allPercents.length + " target percentages");
		OptimizedQCThresholds[][] optqcs = collectOptqcs(processors, qcits, allPercents, log);
		summarizeOptqcs(proj, optqcs, output, log);
	}

	public static void filterByComparison(Project proj, String plinkCnvQCs, String duplicatesFile ,Logger log) {
		Hashtable<String, Hashtable<String, Integer>> defineCompHash = defineCompLists(proj.getProjectDir(), duplicatesFile, log);
		CNVariantQC.filterCNVQCsByComparison(proj, plinkCnvQCs, defineCompHash, log);

	}
	
	public static void filterCNVsByQCThresholds(Project proj, String plinkCnvQCs, String SampleQCFile, String qcThresholdFileName, String output, String QCsubset, int optimizationType, Logger log) {
		OptimizedQCThresholds qcThresholds = OptimizedQCThresholds.loadThresholdsFromTxt(proj.getProjectDir() + qcThresholdFileName, log);
		log.report(ext.getTime() + " Loaded thresholds from " + proj.getProjectDir() + qcThresholdFileName);
		log.report(ext.getTime() + " Retrieving sample qc data from " + proj.getProjectDir() + SampleQCFile);
		Hashtable<String, CNVSampleQC> cnvSampleQCHash = CNVSampleQC.getSampleQCs(proj, SampleQCFile, log);
		log.report(ext.getTime() + " Loading cnvQCs from " + proj.getProjectDir() + plinkCnvQCs);
		CNVariantQC[] cnVariantQCs = CNVQC.load(proj.getProjectDir() + plinkCnvQCs, false).getCnVariantQCs();
		log.report(ext.getTime() + " Finished loading cnvQCs from " + proj.getProjectDir() + plinkCnvQCs);
		String[] inds = CNVariantQC.getIDList(cnVariantQCs, null);
		CNVariantQC[][] unfilteredcnvsQCs = new CNVariantQC[inds.length][];
		Hashtable<String, CNVariantQC[]> allIndcnVariantQCs = CNVariantQC.getIndCNVQCs(inds, cnVariantQCs);
		for (int i = 0; i < inds.length; i++) {
			unfilteredcnvsQCs[i] = allIndcnVariantQCs.get(inds[i]);
		}
		log.report(ext.getTime() + " Beginning to filter " + plinkCnvQCs + " using filter type " + CNVComparison.QC_PARAMETERs[optimizationType]);
		CNVComparison filteredOnly = new CNVComparison(unfilteredcnvsQCs, cnvSampleQCHash, qcThresholds, optimizationType, log);
		log.report(ext.getTime() + " Finished Filtering " + plinkCnvQCs + " using filter type " + CNVComparison.QC_PARAMETERs[optimizationType]);
		CNVariantQC[][] filteredcnvsQCs = filteredOnly.getFilteredcnvQCs1();
		log.report(ext.getTime() + " Creating output files");
		summarizeFiltering(proj, inds, filteredcnvsQCs, cnvSampleQCHash, output, QCsubset, log);
		log.report(ext.getTime() + " Finished");

	}



	public static void convertToQCFormat(Project proj, String plinkCnvs, String markerMAFser, String output, String QCsubset, int threads, Logger log) {
		String[] inds;
		MarkerSet markerSet = proj.getMarkerSet();
		int[][] indices = markerSet.getIndicesByChr();
		int[] positions = markerSet.getPositions();
		String[] markerNames = markerSet.getMarkerNames();
		MarkerFreqs markerFreqs = MarkerFreqs.load(proj.getProjectDir() + markerMAFser, false);
		double[] mafs = markerFreqs.getMafs();
		if (markerFreqs.getFingerprint() != markerSet.getFingerprint()) {
			log.reportError("Error - mismatched marker fingerprints in the project's marker set and the imported AlleleFrequency file (" + markerMAFser + "); aborting");
			System.exit(1);
		}
		log.report(ext.getTime() + " Retrieving cnvs from " + proj.getProjectDir() + plinkCnvs);
		CNVariantQC[] cnVariantQCs = CNVariantQC.getCNVariantQCFromPlinkFile(proj, plinkCnvs);
		log.report(ext.getTime() + " Finished retrieving cnvs from " + proj.getProjectDir() + plinkCnvs);
		if (QCsubset != null) {
			log.report(ext.getTime() + " Filtering cnvs by the QC subset file");
			inds = CNVariantQC.getIDList(cnVariantQCs, defineCompLists(proj.getProjectDir(), QCsubset, log));
			log.report(ext.getTime() + " Finished filtering cnvs by the QC subset file , computing cnv QC metrics for " + inds.length + " individuals");
		} else {
			inds = CNVariantQC.getIDList(cnVariantQCs, null);
			log.report(ext.getTime() + " Using " + inds.length + " individuals to compute cnv QC metrics ");
		}
		Hashtable<String, CNVariantQC[]> allIndcnVariantQCs = CNVariantQC.getIndCNVQCs(inds, cnVariantQCs);
		Hashtable<String, Double> markerMAFhash = hashMAFs(markerNames, mafs);
		log.report(ext.getTime() + " Retrieving markers in cnvs and assigning MAFs");
		for (int i = 0; i < inds.length; i++) {
			for (int j = 0; j < allIndcnVariantQCs.get(inds[i]).length; j++) {
				allIndcnVariantQCs.get(inds[i])[j].findMarkerNamesinCNV(proj, indices, positions, markerNames, log);
				allIndcnVariantQCs.get(inds[i])[j].assignMAFs(markerMAFhash, log);
			}
		}
		CNVariantQC[][] allIndcnVariantQCsArrays = CNValidate.computeMultiThreadedValidations(proj, inds, allIndcnVariantQCs, markerSet, threads, log);
		printNewCNVariantQCFile(proj, output, inds, allIndcnVariantQCsArrays);
		log.report(ext.getTime() + " Completed QC computations");
	}

	// TODO to compare two different plink format files;
	public static void compareFiles(Project proj, String[] plinkCnvQCs, Logger log) {
		Hashtable<Integer, Hashtable<String, CNVariantQC[]>> fileIndcnVariantQCs = new Hashtable<Integer, Hashtable<String, CNVariantQC[]>>();
		String[][] inds = new String[plinkCnvQCs.length][];
		int[][] allPossibleCombinations = Maths.getIndicesForAllCombinations(plinkCnvQCs.length, 2);
		//CNVComparison cnvComparison;
		for (int i = 0; i < plinkCnvQCs.length; i++) {
			CNVariantQC[] cnVariantQCs = CNVariantQC.getCNVariantQCFromPlinkFile(proj, plinkCnvQCs[i]);
			inds[i] = CNVariantQC.getIDList(cnVariantQCs, null);
			fileIndcnVariantQCs.put(i, CNVariantQC.getIndCNVQCs(inds[i], cnVariantQCs));
		}
		for (int i = 0; i < allPossibleCombinations.length; i++) {
			//ArrayList<CNVariantQC[]> cnvQCX = new ArrayList<CNVariantQC[]>();
			//ArrayList<CNVariantQC[]> cnvQCY = new ArrayList<CNVariantQC[]>();
			int fileX = allPossibleCombinations[i][0];
			int fileY = allPossibleCombinations[i][1];
			String[] indsx = inds[fileX];
			String[] indsy = inds[fileY];
			for (int j = 0; j < indsx.length; j++) {
				for (int k = 0; k < indsy.length; k++) {
					if (fileIndcnVariantQCs.get(fileX).containsKey(indsx[j]) && fileIndcnVariantQCs.get(fileY).containsKey(indsy[k])) {

					}
				}
			}
		}
	}

	public OptimizedQCThresholds[][] getOptqcs() {
		return optqcs;
	}

	public double[] getTargetPercentages() {
		return targetPercentages;
	}

	private static Hashtable<String, Double> hashMAFs(String[] markerNames, double[] mafs) {
		Hashtable<String, Double> markerMAFs = new Hashtable<String, Double>();
		for (int i = 0; i < markerNames.length; i++) {
			markerMAFs.put(markerNames[i], mafs[i]);
		}
		return markerMAFs;
	}

	private static void summarizeFiltering(Project proj, String[] inds, CNVariantQC[][] filteredcnvsQCs, Hashtable<String, CNVSampleQC> cnvSampleQCHash, String output, String QCsubset, Logger log) {
		PrintWriter sampleWriter, cnvWriter, famWriter;
		Hashtable<String, Hashtable<String, Integer>> defineCompHash = null;
		if (QCsubset != null) {
			defineCompHash = defineCompLists(proj.getProjectDir(), QCsubset, log);
		}
		try {
			sampleWriter = new PrintWriter(new FileWriter(proj.getProjectDir() + output + ".txt"));
			famWriter = new PrintWriter(new FileWriter(proj.getProjectDir() + output + ".fam"));
			cnvWriter = new PrintWriter(new FileWriter(proj.getProjectDir() + output + ".cnv"));
			sampleWriter.println("Sample\tPassQC?0:1");
			cnvWriter.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
			for (int i = 0; i < inds.length; i++) {
				if (defineCompHash.containsKey(inds[i]) || QCsubset == null) {
					if (filteredcnvsQCs[i] != null) {
						String lookup = filteredcnvsQCs[i][0].getSourceFile();
						int sex = proj.getSampleData(0, false).getSexForIndividual(lookup);
						famWriter.println(filteredcnvsQCs[i][0].getCnVariant().getFamilyID() + "\t" + filteredcnvsQCs[i][0].getCnVariant().getIndividualID() + "\t0\t0\t" + sex + "\t1");
						sampleWriter.println(inds[i] + "\t0");
						for (int j = 0; j < filteredcnvsQCs[i].length; j++) {
							cnvWriter.println(filteredcnvsQCs[i][j].getCnVariant().toPlinkFormat());
						}
					} else {
						sampleWriter.println(inds[i] + "\t1");
					}
				}
			}
			sampleWriter.close();
			famWriter.close();
			cnvWriter.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + proj.getProjectDir() + output);
			e.printStackTrace();
			System.exit(1);
		}
	}
	private OptimizedQCThresholds[] iterate(Hashtable<String, CNVSampleQC> cnvSampleQCHash, CNVariantQC[][][] cnvQCsAssigned, double targetConcordancePercentage, int optimizationType, Logger log) {
		OptimizedQCThresholds[] bestOptqcs = OptimizedQCThresholds.getNewOptqcs(targetConcordancePercentage, 6);
		OptimizedQCThresholds[] qcThresholds = getQCIterations();
		log.report(ext.getTime() + " Iterating " + qcThresholds.length + " parameter combinations");
		for (int i = 0; i < qcThresholds.length; i++) {
			CNVComparison cnvComp = new CNVComparison(cnvQCsAssigned[0], cnvQCsAssigned[1], cnvSampleQCHash, qcThresholds[i], optimizationType, log);
			bestOptqcs = checkThresholds(cnvComp, bestOptqcs, targetConcordancePercentage, qcThresholds[i]);
		}
		return bestOptqcs;
	}

	private static OptimizedQCThresholds[] getQCIterations(){
		// TODO these will be input parameters, with defaults if not provided
		// CNV -specific
		int numMarkerStart = 20;
		int numMarkerStop = 20;
		double[] confCutoffs = binIt(.75, .75, 300);
		double[] alphas = binIt(0.4, .4, 50);
		// less than
		double[] BAFQCcutoffs = binIt(200000, 200000, 1000);
		// greater than
		double[] twopqCutoffs = binIt(0, 0, 10);
		// less than
		double[] hetCutoffs = binIt(1, 1, 10);
		// less than
		double[] cnvLRRSTDevs = binIt(100, 100, 100);
		double[] pennconf = binIt(0, 0, 10);
		double[] bafDrifts = binIt(0.03, .03, 10);
		double[] kbSize = binIt(0, 0, 10);
		double[] kbDensity = binIt(0, 0, 10);
		// SampleSpecific
		int numSampleCNVsStart = 100;
		int numSampleCNVsStop = 100;
		double[] lrrCutoffs = binIt(0.35, .35, 20);
		double[] GCWFCutoffs = binIt(0.02, 0.02, 20);
		double[] sampleCallRates = binIt(0.96, .96, 20);
		
		ArrayList<OptimizedQCThresholds> qcThresholds = new ArrayList<OptimizedQCThresholds>();
		// qcThresholds.add(new OptimizedQCThresholds(alphas[i], confCutoffs[j], lrrCutoffs[k], l, BAFQCcutoffs[m], twopqCutoffs[n], hetCutoffs[o], GCWFCutoffs[p], q));
		// OptimizedQCThresholds noFilter = new OptimizedQCThresholds(0.5, 1.5, 20.0, 200000, 1, 0.0, 1.0, 10.0, 80000, 10.0);
		// qcThresholds.add(noFilter);
		// public OptimizedQCThresholds(double alpha, double confCutoff, double lrrCutoff, int numMarkers, double BAFQCcutoff, double twopqCutoff, double hetCutoff, double GCWF, int numSampleCNVs, double cnvLRRSTDev) {
		for (int i = 0; i < alphas.length; i++) {
			for (int j = 0; j < confCutoffs.length; j++) {
				for (int k = 0; k < lrrCutoffs.length; k++) {
					for (int l = numMarkerStart; l <= numMarkerStop; l++) {
						for (int m = 0; m < BAFQCcutoffs.length; m++) {
							for (int n = 0; n < twopqCutoffs.length; n++) {
								for (int o = 0; o < hetCutoffs.length; o++) {
									for (int p = 0; p < GCWFCutoffs.length; p++) {
										for (int q = numSampleCNVsStart; q <= numSampleCNVsStop; q++) {
											for (int r = 0; r < cnvLRRSTDevs.length; r++) {
												for (int s = 0; s < sampleCallRates.length; s++) {
													for (int t= 0; t < pennconf.length; t++) {	
														for (int u = 0; u < bafDrifts.length; u++) {
															for (int v = 0; v < kbSize.length; v++) {
																for (int w = 0; w < kbDensity.length; w++) {
																	qcThresholds.add(new OptimizedQCThresholds(alphas[i], confCutoffs[j], lrrCutoffs[k], l, BAFQCcutoffs[m], twopqCutoffs[n], hetCutoffs[o], GCWFCutoffs[p], q, cnvLRRSTDevs[r], sampleCallRates[s], pennconf[t], bafDrifts[u], kbSize[v], kbDensity[w]));
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		return qcThresholds.toArray(new OptimizedQCThresholds[qcThresholds.size()]);
	}

	// This is the main qc tracker, if there are more calls passing the thresholds for a given percent target, the qc parameters will be updated
	// copy number specific
	private static OptimizedQCThresholds[] checkThresholds(CNVComparison cnvComp, OptimizedQCThresholds[] bestOptqcs, double targetConcordancePercentage, OptimizedQCThresholds qcIteration) {
		double[] averageCNPercent = cnvComp.getAverageCNPercent();
		int[] numCallsPassing = cnvComp.getFilteredCallsAvailable();
		for (int i = 0; i < averageCNPercent.length; i++) {
			if (averageCNPercent[i] >= targetConcordancePercentage && bestOptqcs[i].getCallsPassingFilter() < numCallsPassing[i]) {
				bestOptqcs[i] = new OptimizedQCThresholds(qcIteration, targetConcordancePercentage, averageCNPercent[i], numCallsPassing[i], cnvComp.getGoodCalls()[i], cnvComp.getTotalCallsAvailable()[i], cnvComp.getNumberIndsCompared(), i);
			}
		}
		return bestOptqcs;
	}

	private static void printNewCNVariantQCFile(Project proj, String output, String[] inds, CNVariantQC[][] allIndcnVariantQCsArrays) {
			ArrayList<CNVariantQC> toSerialize = new ArrayList<CNVariantQC>();
			for (int i = 0; i < inds.length; i++) {
				CNVariantQC[] indcnVariantQCs = allIndcnVariantQCsArrays[i];
				for (int k = 0; k < indcnVariantQCs.length; k++) {
					toSerialize.add(indcnVariantQCs[k]);
				}
			}

		// CNVariantQC[] toTest = toSerialize.toArray(new CNVariantQC[toSerialize.size()]);

		new CNVQC(toSerialize.toArray(new CNVariantQC[toSerialize.size()])).serialize(proj.getProjectDir() + output + ".ser");
	}

	private static QCIterator[] iteratePercents(Project proj, int processors, Thread[] threads, ArrayList<ArrayList<Double>> cabinet, CNVariantQC[][][] cnvQCsAssigned, Hashtable<String, CNVSampleQC> cnvSampleQCHash, int optimizationType, Logger log) {
		QCIterator[] qcIts = new QCIterator[processors];
		for (int i = 0; i < processors; i++) {
			qcIts[i] = new QCIterator(cnvQCsAssigned, cnvSampleQCHash, getDoubleThreadPercents(cabinet.get(i)), optimizationType, log);
			threads[i] = new Thread(qcIts[i]);
			threads[i].start();
		}
		checkThreadStatus(processors, threads);
		return qcIts;
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



	private static double[] getDoubleThreadPercents(ArrayList<Double> doubles) {
		double[] threadPercents = new double[doubles.size()];
		for (int i = 0; i < doubles.size(); i++) {
			threadPercents[i] = doubles.get(i);
		}
		return threadPercents;
	}

	private static ArrayList<ArrayList<Double>> getcabinet(double[] percents, int processors) {
		ArrayList<ArrayList<Double>> cabinet = new ArrayList<ArrayList<Double>>();

		for (int i = 0; i < processors; i++) {
			cabinet.add(new ArrayList<Double>());
		}
		for (int i = 0; i < percents.length; i++) {
			cabinet.get(i % processors).add(percents[i]);
		}
		return cabinet;
	}

	// Defines the comparisions of interest
	private static Hashtable<String, Hashtable<String, Integer>> defineCompLists(String rootDir, String compFile, Logger log) {
		Hashtable<String, Hashtable<String, Integer>> defineCompHash = new Hashtable<String, Hashtable<String, Integer>>();
		BufferedReader reader;
		String[] line;
		int maxNumComparisions = 0;
		try {
			reader = new BufferedReader(new FileReader(rootDir + compFile));
			while (reader.ready()) {
				line = reader.readLine().trim().split("\t");
				if (line.length > maxNumComparisions) {
					maxNumComparisions = line.length;
				}
				// same id in multiple lines not allowed
				if (alreadyDefinedID(defineCompHash, line, log)) {
					continue;
				}
				defineComparisons(defineCompHash, line, log);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + rootDir + compFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + rootDir + compFile + "\"");
			System.exit(2);
		}
		return defineCompHash;
	}

	private static boolean alreadyDefinedID(Hashtable<String, Hashtable<String, Integer>> defineCompHash, String[] line, Logger log) {
		boolean alreadyDefined = false;
		for (int i = 0; i < line.length; i++) {
			// same id not allowed in same line
			if (defineCompHash.containsKey(line[i])) {
				log.reportError("Warning - duplicate IDs were detectected in the comparision file: " + line[i] + " was seen twice, only the first comparison will be used ");
				alreadyDefined = true;
			}
		}
		return alreadyDefined;
	}

	private static void defineComparisons(Hashtable<String, Hashtable<String, Integer>> defineCompHash, String[] line, Logger log) {
		Hashtable<String, Integer> defined = new Hashtable<String, Integer>();
		for (int i = 0; i < line.length; i++) {
			for (int j = 0; j < line.length; j++) {
				if (i != j && line[i] == line[j]) {
					log.reportError("Warning - duplicate IDs were detectected in the comparision file: " + line[i] + " in column " + i + "and " + line[j] + " in column " + j);
				}
				defined.put(line[j], j);
				defineCompHash.put(line[i], defined);
			}
		}
	}

	private static double[] binIt(double startVal, double stopVal, int numBins) {
		double[] values;
		if(startVal==stopVal){
			values = new double[1];
			values[0] = startVal;
		}else{
			values = new double[numBins + 1];
			double inc = getIncrement(startVal, stopVal, numBins);
			for (int i = 0; i < numBins + 1; i++) {
				values[i] = (inc * i) + startVal;
			}
		}
		return values;
	}

	private static double getIncrement(double startVal, double stopVal, int numBins) {
		return (stopVal - startVal) / numBins;
	}

	// collects results from all threads
	private static OptimizedQCThresholds[][] collectOptqcs(int processors, QCIterator[] qcIts, double[] allPercents, Logger log) {
		OptimizedQCThresholds[][] allOptQcs = new OptimizedQCThresholds[allPercents.length][];
		log.report(ext.getTime() + " Collecting QC thresholds for " + allPercents.length + " target percents from available threads...");
		int indIndex = 0;
		int counter = 0;
		for (int i = 0; i < allPercents.length; i++) {
			counter++;
			if (counter > processors) {
				indIndex += 1;
				counter = 1;
			}

			if (qcIts[i % processors].getTargetPercentages()[indIndex] != allPercents[i]) {
				log.reportError("Error - recieved unmatched results while collecting results for " + qcIts[i % processors].getOptqcs()[indIndex] + "\t" + allPercents[i]);
				System.exit(1);
			} else if (qcIts[i % processors].getTargetPercentages()[indIndex] == allPercents[i]) {
				allOptQcs[i] = qcIts[i % processors].getOptqcs()[indIndex];
			}
		}
		log.report(ext.getTime() + " Sucessfully collected QC thresholds for " + allPercents.length + " target percents from available threads...");
		return allOptQcs;
	}

	private static void summarizeOptqcs(Project proj, OptimizedQCThresholds[][] optqcs, String output, Logger log) {
		try {
			PrintWriter qcoutput = new PrintWriter(new FileWriter(proj.getProjectDir() + output, false));
			qcoutput.println(Array.toStr(OptimizedQCThresholds.OPT_QC_HEADS));
			for (int i = 0; i < optqcs.length; i++) {
				qcoutput.println(Array.toStr(OptimizedQCThresholds.OPT_QC_HEADS));
				for (int k = 0; k < optqcs[i].length; k++) {
					qcoutput.print(optqcs[i][k].getDisplayString());
					qcoutput.print("\n");
				}
			}
			qcoutput.close();
		} catch (IOException e) {
			log.reportError("Could Not Open " + proj.getProjectDir() + output);
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		// String filename = Project.DEFAULT_PROJECT;
		int numArgs = args.length;
		String usage = "TODO";
		// String usage = "\n"+
		// "cnv.qc.QCIterator requires 0-4 arguments\n"+
		// "   (0) project properties filename (i.e. proj="+filename+" (default))\n"+
		// "   (1) plink CNV format .cnv file to test qc thresholds (i.e. cnvs=all_cnvs.cnv\n"+
		// "   (2) a tab delimited file defining the comparisons (i.e. comp=replicates.txt (default)) "+
		// "   (3) name of the output file (i.e. out=qcThresholds.txt (default)"+
		// "";
		String filename = "C:/workspace/Genvisis/projects/ARICGenvisis_CEL_11908.properties";
		// C:\workspace\Genvisis\projects\ARICGenvisis_CEL_11908.properties
		String duplicatesFile = "rootsANDdubs.comp.txt";
		String plinkCnvs = "Not_all_gw6.cnv";
		String MarkerFreqs = "MarkerFreq.ser";
		// String plinkCnvs = "all_gw6.cnv ";
		// String beastHeights = plinkCnvs.replaceAll(".cnv", "heights");
		String plinkCnvQCs = "qcThresholdsFull.txt.beast.ser";
		String logfile = null;
		String QCsubset = null;
		// String output;
		String qcThresholdFileName = "qcParmaters.txt";
		String output = "qcThresholds.txt";
		String SampleQCFile = "Sample_QC.xln";
		int optimizationType = 1;
		boolean convert = false;
		boolean filter = false;
		int threads = 0;
		// String plinkCnvs = null;
		Logger log;
		Project proj;
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("cnvs=")) {
				plinkCnvs = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("comp=")) {
				duplicatesFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				output = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("cnvQCs=")) {
				plinkCnvQCs = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("qcsubset=")) {
				QCsubset = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			}else if (args[i].startsWith("mafs=")) {
				MarkerFreqs = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("params=")) {
				qcThresholdFileName = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("type=")) {
				optimizationType = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("threads=")) {
				threads = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("-convert")) {
				convert = true;
				numArgs--;
			} else if (args[i].startsWith("-filter")) {
				filter = true;
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			proj = new Project(filename, false);
			if(logfile==null){
				logfile = "QClog.txt";
			}
			log = new Logger(proj.getProjectDir() + logfile);
			System.out.println("Logging progress to "+proj.getProjectDir() + logfile);
			if (convert) {
				output += ".GENQC";
				convertToQCFormat(proj, plinkCnvs, MarkerFreqs, output, QCsubset, threads, log);
			} else if (filter) {
				filterCNVsByQCThresholds(proj, plinkCnvQCs, SampleQCFile, qcThresholdFileName, output, QCsubset, optimizationType, log);
			} else {
				output += ".summary";
				optimizeQCThresholds(proj, plinkCnvQCs, duplicatesFile, SampleQCFile, output, log, optimizationType);
			}
			Files.backup(logfile, proj.getProjectDir(), proj.getProjectDir());

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
