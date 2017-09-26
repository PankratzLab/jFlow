package org.genvisis.one.JL;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity.CHROMOSOME_X_STRATEGY;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsResiduals;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.MarkerDataLoader;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;
import org.genvisis.stats.CrossValidation;
import org.genvisis.stats.ICC;
import org.genvisis.stats.LeastSquares.LS_TYPE;

import com.google.common.primitives.Doubles;

public class IntensityCorrectionQC {
	public static final String OUTPUT_EXT = ".icc";
	public static final String[] MASK = {".", "0", "e"};
	public static final String INCLUDED_IN_PC = "INCLUDE_IN_MODEL";
	public static final int INCLUDED_IN_PC_INT = 1;

	public static void ICCtheClasses(Project proj, double[] data, String output, String dir,
																	 int startPC, int stopPC, int jumpPC, LS_TYPE lType,
																	 int numThreads) {
		new File(dir).mkdirs();
		PrincipalComponentsResiduals pcComponentsResiduals = proj.loadPcResids();
		ClassDefinition[] classDefinitions = ClassDefinition.getClassDefinitionsFromSampleData(proj);
		// int[] sampleSex = getSexDef(classDefinitions);
		boolean[] modelDefMask = getModelDefMask(proj, classDefinitions);

		// String[] allClasses = getAllClasses(classDefinitions);
		output = proj.PROJECT_DIRECTORY.getValue() + dir + output + ".icc";
		new File(proj.PROJECT_DIRECTORY.getValue() + dir).mkdirs();
		int numTests = 0;
		for (int j = startPC; j < stopPC; j += jumpPC) {
			numTests++;
		}
		int[] pcsTested = new int[numTests];
		int index = 0;
		for (int j = startPC; j < stopPC; j += jumpPC) {
			pcsTested[index] = j;
			index++;
		}
		ExecutorService executor = Executors.newFixedThreadPool(numThreads);
		Hashtable<String, Future<double[]>> tmpResults = new Hashtable<String, Future<double[]>>();
		double[][] allResults = new double[numTests][];
		for (int i = startPC; i < stopPC; i += jumpPC) {
			tmpResults.put(i + "",
										 executor.submit(new WorkerResidual(data, lType, pcComponentsResiduals,
																												classDefinitions, modelDefMask, i,
																												proj.getLog())));
		}
		index = 0;
		for (int i = startPC; i < stopPC; i += jumpPC) {
			try {
				allResults[index] = tmpResults.get(i + "").get();
				index++;
			} catch (InterruptedException e) {
				proj.getLog()
						.reportError("Error - when running GATK Base recalibraion on internal index " + i);
				proj.getLog().reportException(e);
			} catch (ExecutionException e) {
				proj.getLog()
						.reportError("Error - when running GATK Base recalibraion on internal index " + i);
				proj.getLog().reportException(e);
			}
		}
		executor.shutdown();
		try {
			executor.awaitTermination(10L, TimeUnit.DAYS);
		} catch (InterruptedException e) {
			proj.getLog().reportException(e);
		}
		try {
			PrintWriter writer = Files.openAppropriateWriter(output);
			writer.print("PC");
			for (ClassDefinition classDefinition : classDefinitions) {
				writer.print("\tICC." + classDefinition.getClassTitle());
			}
			writer.println();
			int report = 0;
			for (int i = startPC; i < stopPC; i += jumpPC) {
				writer.println(i + "\t" + ArrayUtils.toStr(allResults[report]));
				report++;
			}
			writer.close();
		} catch (Exception e) {
			proj.getLog().reportError("Error writing to " + output);
			proj.getLog().reportException(e);
		}
	}

	private static class WorkerResidual implements Callable<double[]> {
		private final double[] data;
		private final LS_TYPE lType;
		private final PrincipalComponentsResiduals pcComponentsResiduals;
		private final IntensityCorrectionQC.ClassDefinition[] classDefinitions;
		private final boolean[] modelDefMask;
		private final int pc;
		private final Logger log;

		public WorkerResidual(double[] data, LS_TYPE lType,
													PrincipalComponentsResiduals pcComponentsResiduals,
													IntensityCorrectionQC.ClassDefinition[] classDefinitions,
													boolean[] modelDefMask, int pc, Logger log) {
			this.data = data;
			this.lType = lType;
			this.pcComponentsResiduals = pcComponentsResiduals;
			this.classDefinitions = classDefinitions;
			this.modelDefMask = modelDefMask;
			this.pc = pc;
			this.log = log;
		}

		@Override
		public double[] call() {
			return IntensityCorrectionQC.computeAt(data, lType, pcComponentsResiduals, classDefinitions,
																						 modelDefMask, pc, log);
		}
	}

	public static double[] computeAt(double[] data, LS_TYPE lType,
																	 PrincipalComponentsResiduals pcComponentsResiduals,
																	 ClassDefinition[] classDefinitions, boolean[] modelDefMask,
																	 int pc, Logger log) {
		double[] currentData = null;
		if (pc == 0) {
			currentData = data;
		} else {
			CrossValidation crossValidation = pcComponentsResiduals.getCorrectedDataAt(data, modelDefMask,
																																								 pc, lType,
																																								 "PC" + pc, false);
			currentData = crossValidation.getResiduals();
			if (pc == 100) {
				Files.writeArray(ArrayUtils.toStringArray(currentData),
												 pcComponentsResiduals.getProj().PROJECT_DIRECTORY.getValue() + "DFSD.txt");
				System.exit(1);
			}
		}
		log.report("Info - test on PC " + pc);
		double[] iccs = new double[classDefinitions.length];
		String summary = "";
		summary = summary + pc;
		String title = "";
		for (int j = 0; j < classDefinitions.length; j++) {
			double icc = Double.NaN;
			if (classDefinitions[j].isValidForICC()) {
				String[] classDefs = classDefinitions[j].getClassDefs();
				ICC iccComp = new ICC(currentData, classDefs, MASK, null, false, log);
				iccComp.computeICC();
				icc = iccComp.getICC();
			}
			summary = summary + "\t" + icc;
			iccs[j] = icc;
			title += "\t" + classDefinitions[j].getClassTitle();
		}
		log.report(title);

		log.report(summary);
		return iccs;
	}

	public static void ICCtheClasses(Project proj, String[] markers, int numCorrectionThreads,
																	 int numMarkerThreads, String output, String dir, int startPC,
																	 int stopPC, int jumpPC, boolean mitoMode) {
		new File(dir).mkdirs();
		PrincipalComponentsResiduals pcComponentsResiduals = proj.loadPcResids();
		ClassDefinition[] classDefinitions = ClassDefinition.getClassDefinitionsFromSampleData(proj);
		MarkerDataLoader markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj,
																																																markers);
		int[] sampleSex = getSexDef(classDefinitions);
		boolean[] samplesToUseCluster = proj.getSamplesToInclude(null);
		output = proj.PROJECT_DIRECTORY.getValue() + dir + output + ".icc";
		new File(proj.PROJECT_DIRECTORY.getValue() + dir).mkdirs();
		int numTests = 0;
		for (int j = startPC; j <= stopPC; j += jumpPC) {
			numTests++;
		}
		int[] pcsTested = new int[numTests];
		int index = 0;
		for (int j = startPC; j <= stopPC; j += jumpPC) {
			pcsTested[index] = j;
			index++;
		}
		String[] allClasses = getAllClasses(classDefinitions);

		ICCMarkerResultsBatch icMarkerResultsBatch = new ICCMarkerResultsBatch(markers.length);
		for (int i = 0; i < markers.length; i++) {
			icMarkerResultsBatch.addICCMarkerResults(i, new ICCMarkerResults(markers[i], allClasses,
																																			 pcsTested));
			MarkerData markerData = markerDataLoader.requestMarkerData(i);
			proj.getLog().report("Info - marker " + i + " of " + markers.length);
			float[] lrrs = markerData.getLRRs();
			int pcIndex = 0;
			byte[] abGenotypes = markerData.getAbGenotypes();
			boolean[] tmpSampleFilter = new boolean[samplesToUseCluster.length];
			Arrays.fill(tmpSampleFilter, true);
			System.out.println(ArrayUtils.booleanArraySum(samplesToUseCluster));
			for (int k = 0; k < abGenotypes.length; k++) {
				if ((mitoMode) && (abGenotypes[k] != 2) && (abGenotypes[k] != 0)) {
					tmpSampleFilter[k] = false;
				}
				if (!samplesToUseCluster[k]) {
					tmpSampleFilter[k] = false;
				}
			}
			for (int j = startPC; j <= stopPC; j += jumpPC) {
				proj.getLog().report("Info - marker " + i + " of " + markers.length + " on PC " + j);
				float[] lrrICC;
				if (j == 0) {
					lrrICC = lrrs;
				} else {
					PrincipalComponentsIntensity principalComponentsIntensity = new PrincipalComponentsIntensity(pcComponentsResiduals,
																																																			 markerData,
																																																			 true,
																																																			 sampleSex,
																																																			 tmpSampleFilter,
																																																			 1.0D,
																																																			 0.0D,
																																																			 null,
																																																			 true,
																																																			 LS_TYPE.REGULAR,
																																																			 2,
																																																			 5,
																																																			 0.0D,
																																																			 0.1D,
																																																			 numCorrectionThreads,
																																																			 false,
																																																			 null,
																																																			 CHROMOSOME_X_STRATEGY.BIOLOGICAL);
					principalComponentsIntensity.correctXYAt(j);
					if (!principalComponentsIntensity.isFail()) {
						lrrICC = principalComponentsIntensity.getCorrectedIntensity("BAF_LRR", true)[1];
					} else {
						lrrICC = lrrs;
					}
				}
				for (int k = 0; k < classDefinitions.length; k++) {
					double icc = (0.0D / 0.0D);
					if (classDefinitions[k].isValidForICC()) {
						String[] classDefs = new String[classDefinitions[k].getClassDefs().length];
						if (mitoMode) {
							for (int l = 0; l < tmpSampleFilter.length; l++) {
								if (tmpSampleFilter[l]) {
									classDefs[l] = MASK[1];
								} else {
									classDefs[l] = classDefinitions[k].getClassDefs()[l];
								}
							}
						}
						ICC iccComp = new ICC(ArrayUtils.toDoubleArray(lrrICC), classDefs, MASK, null, false,
																	proj.getLog());
						iccComp.computeICC();
						icc = iccComp.getICC();
					}
					icMarkerResultsBatch.addICC(i, k, pcIndex, icc);
				}
				pcIndex++;
			}
		}
		icMarkerResultsBatch.serialize(output);
	}
	//
	// private static void getICC(Project proj, int numCorrectionThreads, PrincipalComponentsResiduals
	// pcComponentsResiduals, ClassDefinition[] classDefinitions, int[] sampleSex, boolean[]
	// samplesToUseCluster, ICCMarkerResultsBatch icMarkerResultsBatch, int i, MarkerData markerData,
	// float[] lrrs, int pcIndex, int j) {
	// float[] lrrICC;
	// if (j == 0) {
	// lrrICC = lrrs;
	// } else {
	// PrincipalComponentsIntensity principalComponentsIntensity = new
	// PrincipalComponentsIntensity(pcComponentsResiduals, markerData, true, sampleSex,
	// samplesToUseCluster, 1.0D, 0.0D, null, true, false, 2, 5, 0.0D, 0.1D, numCorrectionThreads,
	// false, null);
	// principalComponentsIntensity.correctXYAt(j);
	// if (!principalComponentsIntensity.isFail()) {
	// lrrICC = principalComponentsIntensity.getCorrectedIntensity("BAF_LRR", true)[1];
	// } else {
	// lrrICC = lrrs;
	// }
	// }
	// for (int k = 0; k < classDefinitions.length; k++) {
	// double icc = (0.0D / 0.0D);
	// if (classDefinitions[k].isValidForICC()) {
	// ICC iccComp = new ICC(Array.toDoubleArray(lrrICC), classDefinitions[k].getClassDefs(), MASK,
	// null, false, proj.getLog());
	// iccComp.computeICC();
	// icc = iccComp.getICC();
	// }
	// icMarkerResultsBatch.addICC(i, k, pcIndex, icc);
	// }
	// }

	// private static class WorkerICC implements Callable<CrossValidation> {
	// private PrincipalComponentsResiduals principalComponentsResiduals;
	// private double[] data;
	// private boolean[] samplesTobuildModel;
	// private int clusterComponent;
	// private boolean svdRegression;
	// private String title;
	//
	// public WorkerICC(PrincipalComponentsResiduals principalComponentsResiduals, double[] data,
	// boolean[] samplesTobuildModel, int clusterComponent, boolean svdRegression, String title,
	// Logger log) {
	// this.principalComponentsResiduals = principalComponentsResiduals;
	// this.data = data;
	// this.samplesTobuildModel = samplesTobuildModel;
	// this.clusterComponent = clusterComponent;
	// this.svdRegression = svdRegression;
	// this.title = title;
	// }
	//
	// public CrossValidation call() {
	// return this.principalComponentsResiduals.getCorrectedDataAt(this.data,
	// this.samplesTobuildModel, this.clusterComponent, this.svdRegression, this.title, true);
	// }
	// }

	public static void dumpToText(Project proj, String dir) {
		String[] batches = Files.list(proj.PROJECT_DIRECTORY.getValue() + dir, ".icc");
		ICCMarkerResultsBatch[] icMarkerResultsBatchs = new ICCMarkerResultsBatch[batches.length];
		String[] classDefs = null;
		int[] pcsTested = null;
		// double[][] classPCAverages = null;
		// int[][] classPCCounts = null;
		if (icMarkerResultsBatchs.length < 1) {
			proj.getLog().reportError("Error - did not find any result files");
			return;
		}
		icMarkerResultsBatchs[0] = ICCMarkerResultsBatch.loadSerial(proj.PROJECT_DIRECTORY.getValue()
																																+ dir + batches[0]);

		classDefs = icMarkerResultsBatchs[0].getClasses();
		pcsTested = icMarkerResultsBatchs[0].getPcsTested();
		// classPCAverages = new double[classDefs.length][pcsTested.length];
		// classPCCounts = new int[classDefs.length][pcsTested.length];

		ICCClassResults[] icClassResults = new ICCClassResults[classDefs.length];
		for (int i = 0; i < icClassResults.length; i++) {
			icClassResults[i] = new ICCClassResults(pcsTested, classDefs[i]);
		}
		for (int i = 0; i < icMarkerResultsBatchs.length; i++) {
			icMarkerResultsBatchs[i] = ICCMarkerResultsBatch.loadSerial(proj.PROJECT_DIRECTORY.getValue()
																																	+ dir + batches[i]);
			System.out.println(batches[i]);
			if (icMarkerResultsBatchs[i].verify(classDefs, pcsTested)) {
				String[] classes = icMarkerResultsBatchs[i].getClasses();
				for (int j = 0; j < classes.length; j++) {
					String currentOutput = proj.PROJECT_DIRECTORY.getValue() + dir + classes[j]
																 + "_fullResults.txt";
					try {
						double[][] iccs = icMarkerResultsBatchs[i].getAllICCsForClass(j);
						for (int k = 0; k < iccs.length; k++) {
							for (int k2 = 0; k2 < iccs[k].length; k2++) {
								if (!Double.isNaN(iccs[k][k2])) {
									icClassResults[j].addPCICC(k2, iccs[k][k2]);
								}
							}
						}
					} catch (Exception e) {
						proj.getLog().reportError("Error writing to " + currentOutput);
						proj.getLog().reportException(e);
					}
				}
			} else {
				proj.getLog()
						.reportError("Error - could not verify that class definitions and pcs were in the same order");
				return;
			}
		}
		for (int i = 0; i < icClassResults.length; i++) {
			icClassResults[i].finalizeArrays();
			String currentOutput = proj.PROJECT_DIRECTORY.getValue() + dir + classDefs[i]
														 + "_summary.txt";
			try {
				PrintWriter writer = Files.openAppropriateWriter(currentOutput);
				writer.println("NumMarkers\tPC\tAvgICC\tMedianICC\tStdevICC");
				for (int j = 0; j < pcsTested.length; j++) {
					writer.println(icClassResults[i].getSizeAt(j) + "\t" + pcsTested[j] + "\t"
												 + icClassResults[i].getAvgAt(j) + "\t" + icClassResults[i].getMedianAt(j)
												 + "\t" + icClassResults[i].getStdevAt(j));
				}
				writer.close();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + currentOutput);
				proj.getLog().reportException(e);
			}
		}
	}

	private static String[] getAllClasses(ClassDefinition[] classDefinitions) {
		String[] allClasses = new String[classDefinitions.length];
		for (int i = 0; i < classDefinitions.length; i++) {
			allClasses[i] = classDefinitions[i].getClassTitle();
		}
		return allClasses;
	}

	public static class ICCClassResults {
		private IntensityCorrectionQC.ArraySpecialLists iccs;
		private final double[][] finalIccs;

		public ICCClassResults(int[] pcsTested, String classDef) {
			iccs = new IntensityCorrectionQC.ArraySpecialLists(pcsTested.length, 100000);
			finalIccs = new double[pcsTested.length][];
		}

		public void addPCICC(int pcIndex, double icc) {
			iccs.addTo(pcIndex, icc);
		}

		public double getAvgAt(int pcIndex) {
			return ArrayUtils.mean(finalIccs[pcIndex], true);
		}

		public int getSizeAt(int pcIndex) {
			return finalIccs[pcIndex].length;
		}

		public double getMedianAt(int pcIndex) {
			if (finalIccs[pcIndex].length < 2) {
				return (0.0D / 0.0D);
			}
			return ArrayUtils.median(finalIccs[pcIndex]);
		}

		public double getStdevAt(int pcIndex) {
			return ArrayUtils.stdev(finalIccs[pcIndex], true);
		}

		public void finalizeArrays() {
			for (int i = 0; i < finalIccs.length; i++) {
				finalIccs[i] = iccs.getAt(i);
			}
			iccs = null;
		}
	}

	private static class ArraySpecialLists {
		private final IntensityCorrectionQC.ArraySpecialList[] arraySpecialLists;

		public ArraySpecialLists(int num, int initCapacity) {
			arraySpecialLists = new IntensityCorrectionQC.ArraySpecialList[num];
			init(initCapacity);
		}

		private void init(int initCapacity) {
			for (int i = 0; i < arraySpecialLists.length; i++) {
				arraySpecialLists[i] = new IntensityCorrectionQC.ArraySpecialList(initCapacity);
			}
		}

		public void addTo(int index, double d) {
			arraySpecialLists[index].add(Double.valueOf(d));
		}

		public double[] getAt(int index) {
			return Doubles.toArray(arraySpecialLists[index]);
		}

		// public double getMedianAt(int index) {
		// return Array.median(getAt(index));
		// }
	}

	private static class ArraySpecialList extends ArrayList<Double> {
		private static final long serialVersionUID = 1L;

		public ArraySpecialList(int initCapacity) {
			super();
		}
	}

	public static class ICCMarkerResultsBatch implements Serializable {
		private static final long serialVersionUID = 1L;
		private final IntensityCorrectionQC.ICCMarkerResults[] icMarkerResults;

		public ICCMarkerResultsBatch(int init) {
			icMarkerResults = new IntensityCorrectionQC.ICCMarkerResults[init];
		}

		public void addICC(int markerIndex, int classDefIndex, int pcIndex, double icc) {
			icMarkerResults[markerIndex].addICC(classDefIndex, pcIndex, icc);
		}

		public String getMarkerAt(int markerIndex) {
			return icMarkerResults[markerIndex].getMarker();
		}

		public void addICCMarkerResults(int index,
																		IntensityCorrectionQC.ICCMarkerResults icMarkerResult) {
			icMarkerResults[index] = icMarkerResult;
		}

		public String[] getClasses() {
			return icMarkerResults[0].getClassDefs();
		}

		public int[] getPcsTested() {
			return icMarkerResults[0].getPcsTested();
		}

		public boolean verify(String[] classDefs, int[] pcsTested) {
			for (int i = 0; i < icMarkerResults.length; i++) {
				if (!ArrayUtils.equals(icMarkerResults[i].getClassDefs(), classDefs, true)) {
					System.out.println(i + " classDefs");

					return false;
				}
				if (!Arrays.equals(pcsTested, icMarkerResults[i].getPcsTested())) {
					System.out.println(i + " pcs");
					System.out.println(ArrayUtils.toStr(pcsTested));
					System.out.println(ArrayUtils.toStr(icMarkerResults[i].getPcsTested()));
					return false;
				}
			}
			return true;
		}

		public double[][] getAllICCsForClass(int classDefIndex) {
			double[][] iccs = new double[icMarkerResults.length][];
			for (int i = 0; i < iccs.length; i++) {
				iccs[i] = icMarkerResults[i].getIccsForClass(classDefIndex);
			}
			return iccs;
		}

		public void serialize(String fullPathToOutput) {
			SerializedFiles.writeSerial(this, fullPathToOutput);
		}

		public static ICCMarkerResultsBatch loadSerial(String fullPathToSerial) {
			return (ICCMarkerResultsBatch) SerializedFiles.readSerial(fullPathToSerial);
		}
	}

	public static class ICCMarkerResults implements Serializable {
		private static final long serialVersionUID = 1L;
		private final String marker;
		private final String[] classDefs;
		private final double[][] ICCPcResults;
		private final int[] pcsTested;

		public ICCMarkerResults(String marker, String[] classDefs, int[] pcsTested) {
			this.marker = marker;
			this.classDefs = classDefs;
			this.pcsTested = pcsTested;
			ICCPcResults = new double[classDefs.length][pcsTested.length];
		}

		public void addICC(int classDefIndex, int PCIndex, double icc) {
			ICCPcResults[classDefIndex][PCIndex] = icc;
		}

		public String getMarker() {
			return marker;
		}

		public int[] getPcsTested() {
			return pcsTested;
		}

		public double[] getIccsForClass(int classDefIndex) {
			return ICCPcResults[classDefIndex];
		}

		public String[] getClassDefs() {
			return classDefs;
		}
	}

	public static class ClassDefinition {
		private final String[] classDefs;
		private final String classTitle;
		private final boolean sexClass;
		private boolean validForICC;
		private final boolean includedInPCDef;

		public ClassDefinition(Project proj, String classTitle, int numSamples) {
			this.classTitle = classTitle;
			classDefs = new String[numSamples];
			sexClass = (ext.indexOfStr(classTitle, SampleData.EUPHEMISMS, false, true) >= 0);
			includedInPCDef = (ext.indexOfStr(classTitle, new String[] {"INCLUDE_IN_MODEL"}, false,
																				true) >= 0);
		}

		public boolean isValidForICC() {
			return validForICC;
		}

		public boolean isIncludedInPCDef() {
			return includedInPCDef;
		}

		public void addClassDef(String classDef, int index) {
			classDefs[index] = classDef;
		}

		public void determineValidForICC() {
			Hashtable<String, String> track = new Hashtable<String, String>();
			validForICC = false;
			for (String classDef : classDefs) {
				if (ext.indexOfStr(classDef, IntensityCorrectionQC.MASK, true, true) < 0) {
					track.put(classDef, classDef);
					if (track.size() >= 2) {
						validForICC = true;
						break;
					}
				}
			}
		}

		public String getClassTitle() {
			return classTitle;
		}

		public String[] getClassDefs() {
			return classDefs;
		}

		public boolean isSexClass() {
			return sexClass;
		}

		public static ClassDefinition[] getClassDefinitionsFromSampleData(Project proj) {
			SampleData sampleData = proj.getSampleData(false);
			String[] samples = proj.getSamples();
			int numClasses = sampleData.getNumActualClasses();
			boolean[] samplesToUse = proj.getSamplesToInclude(null);
			ClassDefinition[] classDefinitions = new ClassDefinition[numClasses];
			for (int i = 0; i < numClasses; i++) {
				classDefinitions[i] = new ClassDefinition(proj, sampleData.getClassName(i),
																									proj.getSamples().length);
				for (int j = 0; j < samples.length; j++) {
					if (samplesToUse[j]) {
						classDefinitions[i].addClassDef(sampleData.getClassForInd(samples[j], i) + "", j);
					} else {
						classDefinitions[i].addClassDef(IntensityCorrectionQC.MASK[1], j);
					}
				}
			}
			for (ClassDefinition classDefinition : classDefinitions) {
				classDefinition.determineValidForICC();
			}
			return classDefinitions;
		}
	}

	private static int[] getSexDef(ClassDefinition[] classDefinitions) {
		for (ClassDefinition classDefinition : classDefinitions) {
			if (classDefinition.isSexClass()) {
				return ArrayUtils.toIntArray(classDefinition.getClassDefs());
			}
		}
		return null;
	}

	private static boolean[] getModelDefMask(Project proj, ClassDefinition[] classDefinitions) {
		int[] modelDef = new int[proj.getSamples().length];
		Arrays.fill(modelDef, 1);
		for (ClassDefinition classDefinition : classDefinitions) {
			if (classDefinition.isIncludedInPCDef()) {
				modelDef = ArrayUtils.toIntArray(classDefinition.getClassDefs());
			}
		}
		boolean[] modelDefMask = new boolean[modelDef.length];
		for (int i = 0; i < modelDefMask.length; i++) {
			modelDefMask[i] = (modelDef[i] == 1 ? true : false);
		}
		return modelDefMask;
	}

	private static double[] loadDataFile(Project proj, String dataFile, Logger log) {
		double[] data = new double[proj.getSamples().length];
		Arrays.fill(data, (0.0D / 0.0D));
		try {
			BufferedReader reader = Files.getAppropriateReader(dataFile);
			reader.readLine();
			ArrayList<String> samples = new ArrayList<String>();
			while (reader.ready()) {
				String[] line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
				samples.add(line[0]);
			}
			reader.close();
			int[] indices = ext.indexLargeFactors(samples.toArray(new String[samples.size()]),
																						proj.getSamples(), true, log, true);

			reader = Files.getAppropriateReader(dataFile);
			reader.readLine();
			int index = 0;
			while (reader.ready()) {
				if (indices[index] >= 0) {
					try {
						data[indices[index]] = Double.parseDouble(reader.readLine().trim()
																														.split(PSF.Regex.GREEDY_WHITESPACE)[1]);
					} catch (NumberFormatException numberFormatException) {
						data[indices[index]] = (0.0D / 0.0D);
					}
				} else {
					data[indices[index]] = (0.0D / 0.0D);
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

	public static void test2(Project proj, String dataFile, LS_TYPE lType, int numThreads,
													 int jumpPC) {
		double[] data = loadDataFile(proj, dataFile, proj.getLog());

		ICCtheClasses(proj, data, "Mito", "mitos/", 0, proj.INTENSITY_PC_NUM_COMPONENTS.getValue(),
									jumpPC, lType, numThreads);
	}

	public static void test(Project proj) {
		MarkerSetInfo markerSet = proj.getMarkerSet();
		int[][] chrInd = markerSet.getIndicesByChr();
		// String[][] chunkMarkers = Array.splitUpStringArray(Array.subArray(proj.getMarkerNames(),
		// chrInd[3]), 300, proj.getLog());

		ICCtheClasses(proj, ArrayUtils.subArray(proj.getMarkerNames(), chrInd[26]), 6, 1, "Mito",
									"mitos/", 0, 1500, 5, true);
		dumpToText(proj, "mitos/");
		for (int i = 0; i < 25; i++) {
		}
		dumpToText(proj, "test/");
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String dataFile = "D:/data/gedi_gwas/ICCTest/MitoMedian.txt";
		boolean svdRegression = false;
		int numThreads = 4;
		int jumpPC = 5;
		String usage = "\nseq.BWA_Analysis requires 2 argument\n";
		usage = usage + "   (1) project to use (i.e. proj=" + filename + " (no default))\n";
		usage = usage + "   (2) data file to evaluate (i.e. data=" + dataFile + " (no default))\n";
		usage = usage + "   (3) use svd regression (i.e. -svd (not the default))\n";
		usage = usage + "   (4) number of threads (i.e. numThreads=" + numThreads + " ( default))\n";
		usage = usage + "   (5) the jump for each pc tested(i.e. jump=" + jumpPC + " ( default))\n";
		for (String arg : args) {
			if ((arg.equals("-h")) || (arg.equals("-help")) || (arg.equals("/h"))
					|| (arg.equals("/help"))) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("proj=")) {
				filename = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("data=")) {
				dataFile = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("-svd")) {
				svdRegression = true;
				numArgs--;
			} else if (arg.startsWith("numThreads=")) {
				numThreads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("jump=")) {
				jumpPC = ext.parseIntArg(arg);
				numArgs--;
			} else {
				System.out.println("Invalid argument " + arg);
			}
		}
		if (numArgs != 0) {
			System.out.println(usage);
			System.exit(1);
		}
		Project proj = new Project(filename, null);
		test2(proj, dataFile, svdRegression ? LS_TYPE.SVD : LS_TYPE.REGULAR, numThreads, jumpPC);
	}
}
