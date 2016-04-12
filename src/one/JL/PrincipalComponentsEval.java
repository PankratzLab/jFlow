package one.JL;

import cnv.analysis.pca.PrincipalComponentsIntensity;
import cnv.analysis.pca.PrincipalComponentsResiduals;
import cnv.filesys.MarkerData;
import cnv.filesys.Project;
import cnv.manage.MarkerDataLoader;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import stats.Correlation;
import stats.ICC;
import stats.IrrTable;
import stats.Quantiles;
import stats.LeastSquares.LS_TYPE;

public class PrincipalComponentsEval {
	public static void evaluatePCs(Project proj, String markers, String evalFile, String output, String samplesFile, int jumpPC, int numQs,LS_TYPE lType) {
		PrincipalComponentsResiduals pcResiduals = proj.loadPcResids();
		double[] eval = loadDataFile(proj, evalFile, proj.getLog());
		boolean[] samplesToUse = proj.getSamplesToInclude( proj.PROJECT_DIRECTORY.getValue()+ samplesFile);
		if (samplesFile != null) {
			boolean[] tmpsamplesToUse = proj.getSamplesToInclude(null);
			for (int i = 0; i < tmpsamplesToUse.length; i++) {
				if (!tmpsamplesToUse[i]) {
					samplesToUse[i] = false;
				}
			}
		}
		proj.getLog().report(ext.getTime() + " Info - models will be constructed with " + Array.booleanArraySum(samplesToUse) + " samples");

		pcResiduals.setMarkersToAssessFile(markers);
		pcResiduals.setRecomputeLRR(true);
		pcResiduals.computeAssessmentDataMedians();
		pcResiduals.setHomozygousOnly(true);

		Files.writeMatrix(new String[][] { proj.getSamples(), Array.toStringArray(pcResiduals.getMedians()) }, ext.parseDirectoryOfFile(evalFile) + ext.rootOf(output) + ".medianValues", "\t");
		try {
			output = ext.parseDirectoryOfFile(evalFile) + ext.rootOf(output) + ".evaluated";
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			for (int i = 0; i <= pcResiduals.getNumComponents(); i += jumpPC) {
				double[] data;
				if (i == 0) {
					data = pcResiduals.getMedians();
				} else {
					data = pcResiduals.getCorrectedDataAt(pcResiduals.getMedians(), samplesToUse, i, lType, "PC" + i, true).getResiduals();
				}
				evalAndPrint(proj, eval, writer, i, data, numQs, i == 0, evalFile);
			}
			writer.close();
		} catch (Exception e) {
			proj.getLog().reportError("Error writing to " + output);
			proj.getLog().reportException(e);
		}
	}

	private static void evalAndPrint(Project proj, double[] eval, PrintWriter writer, int PC, double[] data, int numQs, boolean header, String evalFile) {
		if (header) {
			writer.print("PC\t" + ext.rootOf(evalFile) + ".PearsonCorrel" + "\t" + ext.rootOf(evalFile) + ".PearsonP\t" + ext.rootOf(evalFile) + ".SpearmanCorrel" + "\t" + ext.rootOf(evalFile) + ".SpearmanP" + "\t" + ext.rootOf(evalFile) + ".ICCMETHODS" + "\t" + ext.rootOf(evalFile) + ".ICCMATCHED");
		}
		ArrayList<Double> evalHave = new ArrayList<Double>();
		ArrayList<Double> residHave = new ArrayList<Double>();
		for (int j = 0; j < data.length; j++) {
			if ((!Double.isNaN(data[j])) && (!Double.isNaN(eval[j]))) {
				evalHave.add(Double.valueOf(eval[j]));
				residHave.add(Double.valueOf(data[j]));
			}
		}
		double[] evals = Array.toDoubleArray(evalHave);
		double[] resids = Array.toDoubleArray(residHave);
		double[] iccData = new double[evals.length * 2];
		String[] ICCDef_METHODS = new String[evals.length * 2];
		String[] ICCDef_MATCHED = new String[evals.length * 2];
		double[] normResids = Array.normalize(resids);
		double[] normevals = Array.normalize(evals);
		int index = 0;
		int match = 1;
		for (int j = 0; j < evals.length; j++) {
			iccData[index] = normevals[j];
			ICCDef_METHODS[index] = "QPCR";
			ICCDef_MATCHED[index] = match + "";
			index++;
			iccData[index] = normResids[j];
			ICCDef_METHODS[index] = "GEN_ESTIMATE";
			ICCDef_MATCHED[index] = match + "";
			match++;
			index++;
		}
		ICC icc = new ICC(iccData, ICCDef_METHODS, null, null, true, proj.getLog());
		icc.computeICC();
		double iccM = icc.getICC();
		icc = new ICC(iccData, ICCDef_MATCHED, null, null, true, proj.getLog());
		icc.computeICC();
		double iccMa = icc.getICC();
		double[] correlP = Correlation.Pearson(evals, resids);
		double[] correlS = Correlation.Spearman(new double[][] { evals, resids });
		double[][] pearsonCorrelQs = new double[numQs][];
		double[][] quantileAgreement = new double[numQs][];
		int[] qs = new int[numQs];
		if (numQs >= 2) {
			for (int i = 0; i < numQs; i++) {
				qs[i] = (i + 2);
				String q = "Q" + qs[i];
				if (header) {
					writer.print("\t" + q + ".PearsonCorrel\t" + q + ".PearsonPvalue\t" + q + ".PercentAgreementHigh\t" + q + ".PercentAgreementLow\t" + q + ".KAPPA");
				}
			}
			if (header) {
				writer.println();
			}
			Quantiles[] evalQs = Quantiles.qetQuantilesFor(qs, evals, "QPCR", proj.getLog());
			Quantiles[] residQs = Quantiles.qetQuantilesFor(qs, resids, "US", proj.getLog());
			for (int i = 0; i < qs.length; i++) {
				pearsonCorrelQs[i] = Correlation.Pearson(Array.toDoubleArray(evalQs[i].getQuantileMembershipAsRoundedInt()), Array.toDoubleArray(residQs[i].getQuantileMembershipAsRoundedInt()));
				IrrTable irrTable = new IrrTable(2, evals.length, true, proj.getLog());
				irrTable.addRatings(0, evalQs[i].getQuantileMembershipAsRoundedInt());
				irrTable.addRatings(1, residQs[i].getQuantileMembershipAsRoundedInt());
				irrTable.parseAgreement();
				quantileAgreement[i] = new double[] { irrTable.getPercentAgreementFor(Array.min(irrTable.getUniqRatings())), irrTable.getPercentAgreementFor(Array.max(irrTable.getUniqRatings())), irrTable.getCohensKappa() };
			}
		}
		proj.getLog().report(ext.getTime() + " Info - finshed PC " + PC);
		proj.getLog().report(ext.getTime() + " Pearson " + Array.toStr(correlP) + "\tNumSamples: " + evalHave.size());
		proj.getLog().report(ext.getTime() + " Spearman " + Array.toStr(correlS) + "\tNumSamples: " + evalHave.size());
		proj.getLog().report(ext.getTime() + " ICCMethods " + iccM + "\tNumSamples: " + evalHave.size());
		proj.getLog().report(ext.getTime() + " ICCMatched " + iccMa + "\tNumSamples: " + evalHave.size());

		writer.print(PC + "\t" + Array.toStr(correlP) + "\t" + Array.toStr(correlS) + "\t" + iccM + "\t" + iccMa);
		if (numQs >= 2) {
			for (int i = 0; i < pearsonCorrelQs.length; i++) {
				writer.print("\t" + Array.toStr(pearsonCorrelQs[i]) + "\t" + quantileAgreement[i][0] + "\t" + quantileAgreement[i][1] + "\t" + quantileAgreement[i][2]);
			}
		}
		writer.println();
		writer.flush();
	}

	public static void evaluatePCsIntensityCorrection(Project proj, String markers, String evalFile, String output, String samplesFile, LS_TYPE lType, int numMarkerThreads, int numCorrectionThreads, int jump, int numQs) {
		PrincipalComponentsResiduals pcResiduals = proj.loadPcResids();
		pcResiduals.setHomozygousOnly(true);
		double[] eval = loadDataFile(proj, evalFile, proj.getLog());
		boolean[] samplesToUse = proj.getSamplesToInclude( proj.PROJECT_DIRECTORY.getValue()+ samplesFile);
		String[] markersToUse = HashVec.loadFileToStringArray( proj.PROJECT_DIRECTORY.getValue() + markers, false, new int[1], true);
		MarkerData[] markerDatas = new MarkerData[markersToUse.length];
		MarkerDataLoader markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markersToUse);
		byte[][] abGentypes = new byte[markersToUse.length][proj.getSamples().length];
		boolean[][] homozygousMarkersToUse = new boolean[markerDatas.length][proj.getSamples().length];
		for (int i = 0; i < markerDatas.length; i++) {
			markerDatas[i] = markerDataLoader.requestMarkerData(i);
			abGentypes[i] = markerDatas[i].getAbGenotypes();
			for (int j = 0; j < abGentypes[i].length; j++) {
				if ((samplesToUse[j]) && ((abGentypes[i][j] == 0) || (abGentypes[i][j] == 2))) {
					homozygousMarkersToUse[i][j] = true;
				} else {
					homozygousMarkersToUse[i][j] = false;
				}
			}
			markerDataLoader.releaseIndex(i);
		}
		output = ext.parseDirectoryOfFile(evalFile) + ext.rootOf(output) + ".evaluatedByMarker";
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			writer.println("PC\t" + ext.rootOf(evalFile) + ".PearsonCorrel" + "\t" + ext.rootOf(evalFile) + ".PearsonP\t" + ext.rootOf(evalFile) + ".SpearmanCorrel" + "\t" + ext.rootOf(evalFile) + ".SpearmanP" + "\t" + ext.rootOf(evalFile) + ".ICCMETHODS" + "\t" + ext.rootOf(evalFile) + ".ICCMATCHED");
			for (int i = 0; i <= pcResiduals.getNumComponents(); i += jump) {
				double[][] curData = new double[markersToUse.length][proj.getSamples().length];
				double[][] tmpData = new double[proj.getSamples().length][markersToUse.length];
				double[] cnEstimate = new double[proj.getSamples().length];
				ExecutorService executor = Executors.newFixedThreadPool(numMarkerThreads);
				Hashtable<String, Future<double[]>> tmpResults = new Hashtable<String, Future<double[]>>();
				if (i > 0) {
					for (int j = 0; j < markersToUse.length; j++) {
						PrincipalComponentsIntensity pcComponentsIntensity = new PrincipalComponentsIntensity(pcResiduals, markerDatas[j], true, null, homozygousMarkersToUse[j], 1.0D, 0.0D, null, false, lType, 2, 5, 0.0D, 0.1D, numCorrectionThreads, false, null);
						tmpResults.put(j + "", executor.submit(new WorkerCorrection(pcComponentsIntensity, i, markersToUse[j], proj.getLog())));
					}
					for (int j = 0; j < markersToUse.length; j++) {
						try {
							curData[j] = tmpResults.get(j).get();
							for (int k = 0; k < proj.getSamples().length; k++) {
								tmpData[k][j] = curData[j][k];
							}
						} catch (InterruptedException e) {
							proj.getLog().reportError("Error - could running Correction on internal index " + j);
							proj.getLog().reportException(e);
						} catch (ExecutionException e) {
							proj.getLog().reportError("Error - could running Correction on internal index " + j);
							proj.getLog().reportException(e);
						}
					}
					executor.shutdown();
					try {
						executor.awaitTermination(10L, TimeUnit.DAYS);
					} catch (InterruptedException e) {
						proj.getLog().reportException(e);
					}
				} else {
					for (int j = 0; j < markersToUse.length; j++) {
						curData[j] = Array.toDoubleArray(markerDatas[j].getLRRs());
					}
				}
				pcResiduals.setFullData(curData);
				pcResiduals.setAbGenotypesAfterFilters(abGentypes);
				cnEstimate = pcResiduals.getLRRMedian();
				for (int j = 0; j < proj.getSamples().length; j++) {
				}
				evalAndPrint(proj, eval, writer, i, cnEstimate, numQs, i == 0, evalFile);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static class WorkerCorrection implements Callable<double[]> {
		private PrincipalComponentsIntensity pcComponentsIntensity;
		private int atComponent;
		private String currentMarker;
		private Logger log;

		public WorkerCorrection(PrincipalComponentsIntensity pcComponentsIntensity, int atComponent, String currentMarker, Logger log) {
			this.atComponent = atComponent;
			this.pcComponentsIntensity = pcComponentsIntensity;
			this.currentMarker = currentMarker;
			this.log = log;
		}

		public double[] call() {
			this.log.reportTimeInfo("Starting corrections for PC" + this.atComponent + " for marker " + this.currentMarker + " on thread" + Thread.currentThread().getName());
			this.pcComponentsIntensity.correctXYAt(this.atComponent);
			this.log.reportTimeInfo("Finished corrections for PC" + this.atComponent + " for marker " + this.currentMarker + " on thread" + Thread.currentThread().getName());
			return Array.toDoubleArray(this.pcComponentsIntensity.getCorrectedIntensity("BAF_LRR", true)[1]);
		}
	}

	private static double[] loadDataFile(Project proj, String dataFile, Logger log) {
		double[] data = new double[proj.getSamples().length];
		Arrays.fill(data, Double.NaN);
		try {
			BufferedReader reader = Files.getAppropriateReader(dataFile);
			reader.readLine();
			ArrayList<String> samples = new ArrayList<String>();
			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("[\\s]+");
				samples.add(line[0]);
			}
			reader.close();
			int[] indices = ext.indexLargeFactors((String[]) samples.toArray(new String[samples.size()]), proj.getSamples(), true, log, true, true);
			reader = Files.getAppropriateReader(dataFile);
			reader.readLine();
			int index = 0;
			while (reader.ready()) {
				if (indices[index] >= 0) {
					try {
						data[indices[index]] = Double.parseDouble(reader.readLine().trim().split("[\\s]+")[1]);
					} catch (NumberFormatException numberFormatException) {
						data[indices[index]] = Double.NaN;
					}
				} else {
					data[indices[index]] = Double.NaN;
				}
				index++;
			}
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + dataFile + "\" not found in current directory");
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + dataFile + "\"");
		}
		return data;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String markers = null;
		String evalFile = null;
		String sampleFile = null;
		String output = "Mitos.Evaluted";

		boolean svdRegression = false;
		boolean separate = false;
		int numCorrectionThreads = 1;
		int numMarkerThreads = 4;
		int jumpPC = 5;
		int numQs = 0;
		String usage = "\nseq.BWA_Analysis requires 2 argument\n";
		usage = usage + "   (1) project to use (i.e. proj=" + filename + " (no default))\n";
		usage = usage + "   (2) data file of markers to compute median from (i.e. markers=" + markers + " (no default))\n";
		usage = usage + "   (3) data file to evaluate (i.e. eval=" + markers + " (no default))\n";

		usage = usage + "   (4) use svd regression (i.e. -svd (not the default))\n";
		usage = usage + "   (5) number of threads for markers (i.e. numMarkerThreads=" + numMarkerThreads + " ( default))\n";
		usage = usage + "   (5) number of threads for each marker (i.e. numCorrectionThreads=" + numCorrectionThreads + " ( default))\n";

		usage = usage + "   (6) the jump for each pc tested(i.e. jump=" + jumpPC + " ( default))\n";
		usage = usage + "   (7) file of samples to use for model building (i.e. samples=" + markers + " (no default))\n";
		usage = usage + "   (8) correct each marker separately (takes a long time) (not the default))\n";
		usage = usage + "   (9) test the quantiles from 2 to this number for correlation, must be >=2 (i.e. numQs=" + numQs + " (no default))\n";
		for (int i = 0; i < args.length; i++) {
			if ((args[i].equals("-h")) || (args[i].equals("-help")) || (args[i].equals("/h")) || (args[i].equals("/help"))) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("markers=")) {
				markers = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("eval=")) {
				evalFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("samples=")) {
				sampleFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("-svd")) {
				svdRegression = true;
				numArgs--;
			} else if (args[i].startsWith("numMarkerThreads=")) {
				numMarkerThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("numCorrectionThreads=")) {
				numCorrectionThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("-separate")) {
				separate = true;
				numArgs--;
			} else if (args[i].startsWith("jump=")) {
				jumpPC = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("numQs=")) {
				numQs = ext.parseIntArg(args[i]);
				numArgs--;
			} else {
				System.out.println("Invalid argument " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.out.println(usage);
			System.exit(1);
		}
		Project proj = new Project(filename, null, false);
		if (separate) {
			evaluatePCsIntensityCorrection(proj, markers, evalFile, output, sampleFile, svdRegression ? LS_TYPE.SVD : LS_TYPE.REGULAR, numMarkerThreads, numCorrectionThreads, jumpPC, numQs);
		} else {
			evaluatePCs(proj, markers, evalFile, output, sampleFile, jumpPC, numQs, svdRegression ? LS_TYPE.SVD : LS_TYPE.REGULAR);
		}
	}
}
