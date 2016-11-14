package org.genvisis.cnv.analysis.pca;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Random;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.stats.CrossValidation;
import org.genvisis.stats.LeastSquares.LS_TYPE;

/**
 * Class to cross validate a series of pc files (and to assist in batching for computing them)
 * -currently focusing on the number of PCs to include for a given sample size
 *
 */
public class PrincipalComponentsValidate {
	public static final String OUTPUT = ".crossValidate";
	public static final String IND = ".ind.txt";
	public static final String REP = ".rep.txt";
	public static final String OUTSAMP = ".outsample.txt";
	private static final String BATCH = "b";
	private static final String JAVA = "/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/bin/java ";
	private static final String JAR = "-jar /home/pankrat2/lanej/park2.jar ";
	private static final String MT_PIPE = " cnv.manage.MitoPipeline ";
	private static final String VAL_PIPE = " cnv.analysis.pca.PrincipalComponentsValidate ";

	private static final String XMX = "-Xmx23g ";
	private static final int MEMORY_MB = 23000;
	private static final double WALLTIME = 48;

	private static final String[] PC_TYPES = {PrincipalComponentsCompute.OUTPUT_EXT[0],
																						PCA.FILE_EXTs[0]};
	private static final String PC = "PC number from pc type";

	/**
	 * @param proj
	 * @param dir directory under the project directory containing pc files
	 * @param startAtComponent start the validations at this component
	 * @param stopAtComponent stop the validations at this component
	 * @param numPcSamplings how many chunks to validate tha pcs at
	 * @param kfolds how many folds for the cross validation
	 * @param pcType the type of PC ({@link PrincipalComponentsValidate#PC_TYPES}
	 * @param mtMarkers the file listing the set of markers to use as dependent variables
	 * @param numThreads
	 * @param svdRegression use an svd-based regression, much faster for large number of pcs
	 * @param compareAcrossFiles do an all versus all comparison using one file as the training set
	 *        and another as the validation
	 * @param recomputeLRR recompute Log R Ratios for median markers
	 * @param output
	 * @param log
	 * @return
	 */
	public static PrincipalComponentsResiduals[] analyzeRegressions(Project proj, String dir,
																																	int startAtComponent,
																																	int stopAtComponent,
																																	int numPcSamplings, int kfolds,
																																	int pcType, String mtMarkers,
																																	int numThreads, LS_TYPE lType,
																																	boolean compareAcrossFiles,
																																	boolean recomputeLRR,
																																	String output, Logger log) {
		if (numPcSamplings > stopAtComponent) {
			log.reportError("Error - cannot break "	+ stopAtComponent + " components into "
											+ numPcSamplings + " chunks");
			return new PrincipalComponentsResiduals[0];
		} else {
			log.report(ext.getTime() + " Info - initializing regressions");
			PrincipalComponentsResiduals[] principalComponentsResiduals =
																																	initRegressions(proj, dir,
																																									stopAtComponent,
																																									pcType, mtMarkers,
																																									recomputeLRR,
																																									log);
			log.report(ext.getTime() + " Info - finished initializing regressions");
			if (principalComponentsResiduals.length > 0) {
				if (compareAcrossFiles) {
					if (principalComponentsResiduals.length > 1) {
						compareAcrossFiles(	proj, dir, principalComponentsResiduals, startAtComponent,
																stopAtComponent, numPcSamplings, kfolds, pcType, mtMarkers,
																numThreads, lType, output, log);
					} else {
						log.reportError("Error - must have at least two PC files of type "	+ PC_TYPES[pcType]
														+ " in directory " + proj.PROJECT_DIRECTORY.getValue() + dir
														+ " to compare across files");
					}
				} else {
					compareWithinFile(proj, dir, startAtComponent, stopAtComponent, numPcSamplings, kfolds,
														pcType, numThreads, lType, output, log, principalComponentsResiduals);
				}
			} else {
				log.reportError("Error - did not find any PC files of type "	+ PC_TYPES[pcType]
												+ " in directory " + proj.PROJECT_DIRECTORY.getValue() + dir);
			}
			return principalComponentsResiduals;
		}
	}

	/**
	 * Batch the regression validations by pc number. Subdirectories will be created for the analysis
	 * to run in
	 */
	public static PrincipalComponentsResiduals[] analyzeRegressionsBatched(	Project proj, String dir,
																																					int startAtComponent,
																																					int stopAtComponent,
																																					int numPcSamplings,
																																					int kfolds, int pcType,
																																					String mtMarkers,
																																					int numThreads,
																																					boolean svdRegression,
																																					boolean compareAcrossFiles,
																																					boolean recomputeLRR,
																																					String output,
																																					Logger log) {
		PrincipalComponentsResiduals[] principalComponentsResiduals =
																																initRegressions(proj, dir,
																																								stopAtComponent,
																																								pcType, mtMarkers,
																																								recomputeLRR, log);
		ArrayList<String> componentBatches = new ArrayList<String>();
		for (int i = 0; i < principalComponentsResiduals.length; i++) {
			int numTotalComponents = principalComponentsResiduals[i].getTotalNumComponents();
			componentBatches.add(numTotalComponents + "");
			String subDir = dir + BATCH + numTotalComponents + "/";
			if (!new File(proj.PROJECT_DIRECTORY.getValue() + subDir).exists()) {
				new File(proj.PROJECT_DIRECTORY.getValue() + subDir).mkdirs();
			}
			if (!Files.exists(proj.PROJECT_DIRECTORY.getValue()	+ subDir
												+ ext.removeDirectoryInfo(principalComponentsResiduals[i].getPcFile()))) {
				Files.copyFile(proj.PROJECT_DIRECTORY.getValue()
													+ principalComponentsResiduals[i].getPcFile(),
												proj.PROJECT_DIRECTORY.getValue() + subDir + ext.removeDirectoryInfo(principalComponentsResiduals[i].getPcFile()));
			}
		}

		int[] uniqBatches =
											Array.toIntArray(Array.unique(componentBatches.toArray(new String[componentBatches.size()])));
		Arrays.sort(uniqBatches);
		String[][] batches = new String[uniqBatches.length][1];
		for (int i = 0; i < uniqBatches.length; i++) {
			batches[i][0] = uniqBatches[i] + "";
		}
		// String command = JAVA + CP + XMX + VAL_PIPE + "proj=" +
		// proj.getFilename(proj.PROJECT_PROPERTIES_FILENAME) + " dir=" + dir + BATCH + "[%0]/" + "
		// numThreads=" + numThreads;
		String command = JAVA	+ XMX + JAR + VAL_PIPE + "proj="
											+ proj.PROJECT_PROPERTIES_FILENAME.getValue() + " dir=" + dir + BATCH
											+ "[%0]/" + "  numThreads=" + numThreads;
		command += " startAtComponent="	+ startAtComponent + " stopAtComponent=" + stopAtComponent
								+ " numPcIterations=" + numPcSamplings + " pcType=" + pcType;
		command += " kfolds=" + kfolds + " mtMarkers=" + mtMarkers + " output=regress[%0]";
		Files.qsub("Regression", command, batches, MEMORY_MB, WALLTIME, numThreads);
		return principalComponentsResiduals;
	}

	private static void compareWithinFile(Project proj, String dir, int startAtComponent,
																				int stopAtComponent, int numPcSamplings, int kfolds,
																				int pcType, int numThreads, LS_TYPE lType, String output,
																				Logger log,
																				PrincipalComponentsResiduals[] principalComponentsResiduals) {
		assignMedians(principalComponentsResiduals, log);// compute the medians for samples in each PC
		ValidationResults[] validationResults = computeRegressions(	principalComponentsResiduals, null,
																																startAtComponent, stopAtComponent,
																																numPcSamplings, kfolds,
																																numThreads,
																																proj.PROJECT_DIRECTORY.getValue() + dir + "tmp_" + output + OUTPUT,
																																lType, log);
		summarize(validationResults, principalComponentsResiduals, pcType,
							proj.PROJECT_DIRECTORY.getValue() + dir + output + OUTPUT, log);
	}

	private static void summarize(ValidationResults[] validationResults,
																PrincipalComponentsResiduals[] principalComponentsResiduals,
																int pcType, String output, Logger log) {
		PrintWriter writer = Files.getAppropriateWriter(output + IND);

		writer.print(PC + "\"" + PC_TYPES[pcType] + "\"");
		for (PrincipalComponentsResiduals principalComponentsResidual : principalComponentsResiduals) {
			writer.print("\t" + principalComponentsResidual.getOutput());
		}
		writer.println();
		for (ValidationResults validationResult : validationResults) {
			writer.println(validationResult.getSummary());
		}
		writer.close();
		if (validationResults[0].getSortedUniqNumComponents().length < principalComponentsResiduals.length) {// detect
																																																					// equal
																																																					// number
																																																					// of
																																																					// samples
			log.report("Detected that multiple files contained the same number of principal components, these will be treated as replicates");
			log.report("Warning - it is assumed that the number of PCs +1 is equal to the number of samples use to compute components (i.e a file with 99 PCs was computed from 100 samples)");
			summarizeReplicates(validationResults, principalComponentsResiduals, pcType, output, log);
		}
	}

	/**
	 * This runs a cross-validation where the validation data comes from another files PCs and sample
	 * medians.
	 * <p>
	 * We use the same {@link PrincipalComponentsValidate#computeRegressions()} as
	 * {@link PrincipalComponentsValidate#compareWithinFile()} but provide a validation set that comes
	 * from another pc file
	 *
	 */
	private static ValidationResults[][][] compareAcrossFiles(Project proj, String dir,
																														PrincipalComponentsResiduals[] principalComponentsResiduals,
																														int startAtComponent,
																														int stopAtComponent, int numPcSamplings,
																														int kfolds, int pcType,
																														String mtMarkers, int numThreads,
																														LS_TYPE lType, String output,
																														Logger log) {
		assignMedians(principalComponentsResiduals, log);
		ValidationResults[] inSamplevalidationResults =
																									computeRegressions(	principalComponentsResiduals,
																																			null, startAtComponent,
																																			stopAtComponent,
																																			numPcSamplings, kfolds,
																																			numThreads,
																																			proj.PROJECT_DIRECTORY.getValue()
																																										+ dir + "tmp_"
																																									+ output + OUTPUT,
																																			lType, log);
		summarize(inSamplevalidationResults, principalComponentsResiduals, pcType,
							proj.PROJECT_DIRECTORY.getValue() + dir + output + OUTPUT, log);

		// train #,validation #, PC #
		ValidationResults[][][] outOfSampleValidationResults = new ValidationResults[principalComponentsResiduals.length][principalComponentsResiduals.length
																																																											- 1][];// all
																																																															// v
																																																															// all
		String[][] comparisons = new String[outOfSampleValidationResults.length][outOfSampleValidationResults[0].length];

		for (int i = 0; i < principalComponentsResiduals.length; i++) {
			PrincipalComponentsResiduals[] train_pc = {principalComponentsResiduals[i]};
			int index = 0;
			for (int j = 0; j < principalComponentsResiduals.length; j++) {
				if (j != i) {
					// PrincipalComponentsResiduals val_pc = filterExclude(train_pc[0],
					// principalComponentsResiduals[j], log);// remove training individuals from the
					// validation
					PrincipalComponentsResiduals val_pc = principalComponentsResiduals[j];
					String file = train_pc[0].getOutput()	+ "_val_on_" + ext.rootOf(val_pc.getOutput()) + "_"
												+ val_pc.getPcBasis()[0].length + "_individuals";
					if (val_pc.getPcBasis()[0].length > 0) {
						outOfSampleValidationResults[i][index] = computeRegressions(train_pc, val_pc,
																																				startAtComponent,
																																				stopAtComponent,
																																				numPcSamplings, kfolds,
																																				numThreads, null, lType,
																																				log);// this
																																							// can
																																							// create
																																							// alot
																																							// of
																																							// tempory
																																							// files
																																							// we
																																							// skip
																																							// the
																																							// tmp
																																							// reporting
					} else {
						log.reportError("Warning - all individuals from the validation pc file "
															+ principalComponentsResiduals[j].getOutput()
														+ " were present in the training file " + train_pc[0].getOutput());
					}
					comparisons[i][index] = file;
					index++;
				}
			}
		}
		summarizePCGenerators(inSamplevalidationResults, outOfSampleValidationResults,
													principalComponentsResiduals, comparisons, pcType,
													proj.PROJECT_DIRECTORY.getValue() + dir + output, log);
		return outOfSampleValidationResults;
	}

	/**
	 * This summarizes both the within file and out of file validation results
	 */
	private static void summarizePCGenerators(ValidationResults[] inSamplevalidationResults,
																						ValidationResults[][][] outOfSampleValidationResults,
																						PrincipalComponentsResiduals[] principalComponentsResiduals,
																						String[][] comparisons, int pcType, String output,
																						Logger log) {
		PrintWriter writer = Files.getAppropriateWriter(output + OUTSAMP);

		writer.print(PC + "\"" + PC_TYPES[pcType] + "\"");
		for (int i = 0; i < principalComponentsResiduals.length; i++) {
			writer.print("\tinSample_"	+ principalComponentsResiduals[i].getOutput() + "\t"
										+ Array.toStr(comparisons[i]));
		}
		writer.println();
		for (int i = 0; i < inSamplevalidationResults.length; i++) {// for all components tested
			writer.print(inSamplevalidationResults[i].getRegressAtComponent());
			for (int j = 0; j < outOfSampleValidationResults.length; j++) {// for all training file used
				writer.print("\t" + inSamplevalidationResults[i].getAverageSSerrAt(j));
				for (int j2 = 0; j2 < outOfSampleValidationResults[j].length; j2++) {// for all validation
																																							// files used
					if (outOfSampleValidationResults[j][j2] == null) {
						writer.print("\t" + Double.NaN);
					} else {
						writer.print("\t" + outOfSampleValidationResults[j][j2][i].getAverageSSerrAt(0));// only
																																															// one
																																															// train
																																															// file
					}
				}
			}
			writer.println();
		}
		writer.close();
	}

	// /**
	// * We remove individuals contained in toExclude from toFilter by changing the basis and median
	// arrays.
	// */
	// private static PrincipalComponentsResiduals filterExclude(PrincipalComponentsResiduals
	// toExclude, PrincipalComponentsResiduals toFilter, Logger log) {
	// String[] toFilterInds = Collections.list(toFilter.getSamplesInPc().keys()).toArray(new
	// String[toFilter.getSamplesInPc().size()]);
	//
	// Hashtable<String, Integer> toExcludeSamps = toExclude.getSamplesInPc();
	// double[] tmpMedians = toFilter.getMedians();
	// double[][] tmpBasis = toFilter.getPcBasis();
	// Hashtable<Integer, Integer> indicesToKeep = new Hashtable<Integer, Integer>();
	//
	// for (int i = 0; i < toFilterInds.length; i++) {
	// if (!toExcludeSamps.containsKey(toFilterInds[i])) {// we keep an individual if they are not in
	// the toExclude PC
	// int indexToKeep = toFilter.getSamplesInPc().get(toFilterInds[i]);
	// indicesToKeep.put(indexToKeep, indexToKeep);
	// }
	// }
	// double[] filteredMedians = new double[indicesToKeep.size()];
	// double[][] filteredBasis = new double[tmpBasis.length][indicesToKeep.size()];
	// int index = 0;
	// for (int i = 0; i < tmpMedians.length; i++) {
	// if (indicesToKeep.containsKey(i)) {// keep it
	// filteredMedians[index] = tmpMedians[i];
	// index++;
	// }
	// }
	// for (int i = 0; i < tmpBasis.length; i++) {// tmpBasis[i][sample1...] is in same order as
	// medians
	// index = 0;
	// for (int j = 0; j < tmpBasis[0].length; j++) {
	// if (indicesToKeep.containsKey(j)) {
	// filteredBasis[i][index] = tmpBasis[i][j];
	// index++;
	// }
	// }
	// }
	// log.report(ext.getTime() + " Kept " + indicesToKeep.size() + " individuals from " +
	// toFilter.getOutput() + " that were not in the " + toExcludeSamps.size() + " samples from" +
	// toExclude.getOutput() + " for the validation set");
	// PrincipalComponentsResiduals filtered = toFilter.clone();// need to de-reference and create a
	// new one
	// if (indicesToKeep.size() > 0) { // all individuals present in toExclude
	// filtered.setMedians(filteredMedians);
	// filtered.setPcBasis(filteredBasis);
	// } else {
	// filtered.setMedians(new double[0]);
	// filtered.setPcBasis(new double[1][0]);// 0 individuals
	// }
	// return filtered;
	// }

	/**
	 * Warning - pc files with the same number of components are summarize as replicates
	 */
	private static void summarizeReplicates(ValidationResults[] validationResults,
																					PrincipalComponentsResiduals[] principalComponentsResiduals,
																					int pcType, String output, Logger log) {
		PrintWriter writer = Files.getAppropriateWriter(output + REP);
		int[] uniqNumComponents = validationResults[0].getSortedUniqNumComponents();
		writer.print(PC + "\"" + PC_TYPES[pcType] + "\"");
		for (int uniqNumComponent : uniqNumComponents) {
			String numSamples = (uniqNumComponent + 1) + "";// the assumption
			String repString = "("	+ validationResults[0].getNumReplicatesForPCNumber(uniqNumComponent)
													+ ")X ";
			writer.print("\t"	+ repString + numSamples + ".samples.percentRegressed" + "\t" + repString
										+ numSamples + ".samples.averageSSerr" + "\t" + repString + numSamples
										+ ".samples.averageR2");
		}
		writer.println();
		for (ValidationResults validationResult : validationResults) {
			writer.println(validationResult.getReplicateSummary());
		}
		writer.close();
	}

	private static ValidationResults[] computeRegressions(PrincipalComponentsResiduals[] principalComponentsResiduals,
																												PrincipalComponentsResiduals val_pcs,
																												int startAtComponent, int stopAtComponent,
																												int numPcSamplings, int kFolds,
																												int numThreads, String tmpOutput,
																												LS_TYPE lType, Logger log) {
		int[] pcChunks = Array.splitUp((stopAtComponent + 1 - startAtComponent), numPcSamplings);
		int[] componentsToTest = getComponentsToTest(pcChunks, startAtComponent);
		ValidationResults[] validationResults = new ValidationResults[componentsToTest.length];
		initValidations(validationResults, principalComponentsResiduals.length, kFolds,
										componentsToTest);
		for (int i = 0; i < principalComponentsResiduals.length; i++) {
			int totalNumComponents = principalComponentsResiduals[i].getTotalNumComponents();// total
																																												// number of
																																												// PCs in
																																												// the
																																												// file,
																																												// assumed
																																												// to
																																												// be number
																																												// of
																																												// samples
																																												// minus 1
			int numSamples = principalComponentsResiduals[i].getNumSamples();
			addToValidations(	principalComponentsResiduals[i].crossValidate(kFolds, componentsToTest,
																																			numThreads, lType,
																																			(tmpOutput == null	? null
																																													: tmpOutput
																																														+ ext.rootOf(principalComponentsResiduals[i].getOutput())),
																																			val_pcs),
												validationResults, i, totalNumComponents, numSamples);
		}
		return validationResults;
	}

	/**
	 * We add all the cross validations for a particular file to each of the validation results,
	 * results are PC oriented
	 */
	/**
	 * @param crossValidations fully computed validations
	 * @param validationResults the storage container
	 * @param fileIndex
	 * @param totalNumComponents total number of components in the file
	 * @param numSamples total number of samples in the file
	 */
	private static void addToValidations(	CrossValidation[][] crossValidations,
																				ValidationResults[] validationResults, int fileIndex,
																				int totalNumComponents, int numSamples) {
		for (int i = 0; i < crossValidations.length; i++) {
			validationResults[i].addCrossValidation(fileIndex, totalNumComponents, numSamples,
																							crossValidations[i]);
		}
	}

	/**
	 * set up the validations for PC based multi-threading
	 */
	private static void initValidations(ValidationResults[] validationResults, int numFiles,
																			int kFolds, int[] componentsToTest) {
		for (int i = 0; i < validationResults.length; i++) {
			validationResults[i] = new ValidationResults(numFiles, kFolds, componentsToTest[i]);
		}
	}

	/**
	 * This method picks the PCs to test in the cross-validation
	 */
	private static int[] getComponentsToTest(int[] pcChunks, int startAt) {
		int[] componentsToTest = new int[pcChunks.length + 1];
		int currentComponent = 1;
		while (currentComponent < startAt) {
			currentComponent++;
		}
		componentsToTest[0] = currentComponent;
		for (int i = 0; i < pcChunks.length; i++) {
			currentComponent += pcChunks[i];
			componentsToTest[i + 1] = currentComponent - 1;
		}
		return componentsToTest;
	}

	/**
	 * Since files may contain different individuals, we assign the medians separately
	 */
	private static void assignMedians(PrincipalComponentsResiduals[] principalComponentsResiduals,
																		Logger log) {
		for (PrincipalComponentsResiduals principalComponentsResidual : principalComponentsResiduals) {
			principalComponentsResidual.computeAssessmentDataMedians();
		}
	}

	/**
	 * Class that stores results for a particular PC (which may contain multiple files), and
	 * facilitates summarization
	 *
	 */
	private static class ValidationResults {
		private final int regressAtComponent;
		private final int[] numSamples, totalNumComponents;
		private final CrossValidation[][] crossValidations;

		public ValidationResults(int numFiles, int kfolds, int regressAtComponent) {
			super();
			crossValidations = new CrossValidation[numFiles][kfolds];
			numSamples = new int[numFiles];
			totalNumComponents = new int[numFiles];
			this.regressAtComponent = regressAtComponent;
		}

		public void addCrossValidation(	int fileIndex, int totalNumComponent, int numSample,
																		CrossValidation[] crossValidation) {
			crossValidations[fileIndex] = crossValidation;
			totalNumComponents[fileIndex] = totalNumComponent;// how we determine the replicates
			numSamples[fileIndex] = numSample;// in case we ever want to group by number of samples
		}

		public double getAverageSSerrAt(int fileIndex) {
			return CrossValidation.getEstimateError(crossValidations[fileIndex]);
		}

		public double getAverageR2At(int fileIndex) {
			return CrossValidation.getAverageR2(crossValidations[fileIndex]);
		}

		/**
		 * @return number of unique TOTAL number of components between files
		 */
		private Hashtable<Integer, Integer> getNumReps() {// We assume files with the same number of
																											// components were computed from the same
																											// number of samples
			Hashtable<Integer, Integer> numReps = new Hashtable<Integer, Integer>();
			for (int i = 0; i < totalNumComponents.length; i++) {
				if (!numReps.containsKey(totalNumComponents[i])) {
					numReps.put(totalNumComponents[i], 1);
				} else {
					numReps.put(totalNumComponents[i], numReps.get(totalNumComponents[i]) + 1);
				}
			}
			return numReps;
		}

		public int getNumReplicatesForPCNumber(int pcNumber) {
			return getNumReps().get(pcNumber);// a bit inefficient
		}

		/**
		 * across cross-validations and replicates
		 */
		private Hashtable<Integer, Double> getReplicateAvgSSerr() {
			Hashtable<Integer, Double> average = new Hashtable<Integer, Double>();
			Hashtable<Integer, Integer> numReps = getNumReps();
			for (int i = 0; i < totalNumComponents.length; i++) {
				if (!average.containsKey(totalNumComponents[i])) {
					average.put(totalNumComponents[i], getAverageSSerrAt(i));
				} else {
					average.put(totalNumComponents[i],
											average.get(totalNumComponents[i]) + getAverageSSerrAt(i));
				}
			}
			int[] uniqTotalComponents = getUniqNumComponents();
			for (int uniqTotalComponent : uniqTotalComponents) {
				average.put(uniqTotalComponent,
										average.get(uniqTotalComponent) / (double) (numReps.get(uniqTotalComponent)));
			}
			return average;
		}

		/**
		 * across cross-validations and replicates
		 */
		private Hashtable<Integer, Double> getReplicateAvgR2() {
			Hashtable<Integer, Double> average = new Hashtable<Integer, Double>();
			Hashtable<Integer, Integer> numReps = getNumReps();
			for (int i = 0; i < totalNumComponents.length; i++) {
				if (!average.containsKey(totalNumComponents[i])) {
					average.put(totalNumComponents[i], getAverageR2At(i));
				} else {
					average.put(totalNumComponents[i],
											average.get(totalNumComponents[i]) + getAverageR2At(i));
				}
			}
			int[] uniqTotalNumComponents = getUniqNumComponents();
			for (int uniqTotalNumComponent : uniqTotalNumComponents) {
				average.put(uniqTotalNumComponent, average.get(uniqTotalNumComponent)
																						/ (double) (numReps.get(uniqTotalNumComponent)));
			}
			return average;
		}

		public int[] getUniqNumComponents() {
			String[] array = Array.toStringArray(totalNumComponents);
			String[] unique = Array.unique(array);
			return Array.toIntArray(unique);
		}

		public int[] getSortedUniqNumComponents() {
			int[] sort = getUniqNumComponents();
			return sort;
		}

		public int getRegressAtComponent() {
			return regressAtComponent;
		}

		public String getSummary() {
			String summary = "" + getRegressAtComponent();
			for (int i = 0; i < crossValidations.length; i++) {
				summary += "\t" + getAverageSSerrAt(i);
			}
			return summary;
		}

		public String getReplicateSummary() {
			int[] uniqSortedComponents = getSortedUniqNumComponents();
			Hashtable<Integer, Double> averageSSerr = getReplicateAvgSSerr();
			Hashtable<Integer, Double> averageR2 = getReplicateAvgR2();
			String summary = "" + getRegressAtComponent();
			for (int uniqSortedComponent : uniqSortedComponents) {
				summary += "\t"	+ ((double) (regressAtComponent) / (double) (uniqSortedComponent)) + "\t"
										+ averageSSerr.get(uniqSortedComponent) + "\t"
										+ averageR2.get(uniqSortedComponent);
			}
			return summary;
		}
	}

	/**
	 * @param proj
	 * @param dir under the project directory containing pc files
	 * @param numComponents max number of components to load from each pc file
	 * @param pcType
	 * @param mtMarkers (markers for medians and residuals)
	 * @param log
	 * @return
	 */

	private static PrincipalComponentsResiduals[] initRegressions(Project proj, String dir,
																																int numComponents, int pcType,
																																String mtMarkers,
																																boolean recomputeLRR, Logger log) {
		String[] files = Files.toFullPaths(	Files.list(proj.PROJECT_DIRECTORY.getValue()	+ dir,
																									PC_TYPES[pcType], false),
																				dir);
		if (files == null || files.length < 1) {
			log.reportError("Error - could not find any files of type " + PC_TYPES[pcType]);
			return new PrincipalComponentsResiduals[0];
		}
		PrincipalComponentsResiduals[] principalComponentsResiduals =
																																new PrincipalComponentsResiduals[files.length];
		log.report(ext.getTime()	+ " Info - Using the following files for validation\n"
								+ Array.toStr(files, "\n"));
		for (int i = 0; i < files.length; i++) {
			principalComponentsResiduals[i] = new PrincipalComponentsResiduals(	proj, files[i], mtMarkers,
																																					numComponents, false, 0,
																																					true, recomputeLRR,
																																					ext.rootOf(files[i]));
		}
		return principalComponentsResiduals;
	}

	/**
	 * Note: use this to set up a bunch of pca runs Warning: your project should be fully parsed prior
	 * to using this Warning: you should have lrr_sd and markerMetrics computed prior to using this,
	 * otherwise the multiple runs will get confused when filtering files etc..
	 *
	 * @param proj current project
	 * @param dir directory within the project directory
	 * @param indsToValidateFile list of inds to chunk at project directory/dir/
	 * @param numberOfChunks number to have a subset removed (total jobs is numberOfChunks +1)
	 * @param mtMarkers mitochondrial markers to use
	 * @param pcMarkers pc markers to use
	 * @param log
	 */
	public static void batchJobs(	Project proj, String dir, String indsToValidateFile,
																int numberOfChunks, int numComponents, String mtMarkers,
																String pcMarkers, int pcReplicates, Logger log) {
		String curDir = proj.PROJECT_DIRECTORY.getValue() + dir;
		String[] inds = HashVec.loadFileToStringArray(curDir	+ indsToValidateFile, false, new int[] {0},
																									false);
		int[] chunks = Array.splitUp(inds.length, numberOfChunks);
		String[][] batches = getBatches(inds, chunks, pcReplicates, log);
		// Files.writeList(inds, curDir + BATCH + 0);
		for (int i = 0; i < batches.length; i++) {
			String batchName = BATCH + "_" + i + "_" + batches[i].length;
			String batchFile = curDir + batchName;
			Files.writeArray(batches[i], batchFile);
		}
		// String command = JAVA + CP + XMX + MT_PIPE + "proj=" +
		// proj.getFilename(proj.PROJECT_PROPERTIES_FILENAME) + " PCmarkers=" + pcMarkers + "
		// numComponents=[%1] medianMarkers=" + proj.getProjectDir() + mtMarkers + " useFile=" + curDir
		// + BATCH + "_[%0] output=" + dir + BATCH + "_[%0]";
		String command = JAVA	+ XMX + JAR + MT_PIPE + "proj="
											+ proj.PROJECT_PROPERTIES_FILENAME.getValue() + " PCmarkers=" + pcMarkers
											+ " numComponents=[%1]  medianMarkers=" + proj.PROJECT_DIRECTORY.getValue()
											+ mtMarkers + " useFile=" + curDir + BATCH + "_[%0] output=" + dir + BATCH
											+ "_[%0]";
		Files.qsub("PCA", command, getIters(batches, log));
	}

	/**
	 * Just creating a String[][] with increasing index (String[0][0]=0,String[1][1]=1 etc)
	 */
	private static String[][] getIters(String[][] batches, Logger log) {
		String[][] iters = new String[batches.length][2];
		for (int i = 0; i < batches.length; i++) {
			iters[i][0] = i + "_" + batches[i].length;
			iters[i][1] = "" + (batches[i].length - 1);
		}
		return iters;
	}

	/**
	 * We exclude one chunk of inds from each batch, so that each individual appears chunks.length -1
	 * times in batches
	 */

	private static String[][] getBatches(String[] inds, int[] chunks, int replicates, Logger log) {
		int num = 0;
		int batchIndex = 0;
		String[][] batches = new String[chunks.length * replicates][];
		for (int chunk : chunks) {
			num += (chunk + 50) / 100 * 100;
			for (int j = 0; j < replicates; j++) {
				batches[batchIndex] = getRandomSubset(inds, num, log);
				batchIndex++;
			}
		}
		return batches;
	}

	private static String[] getRandomSubset(String[] a, int num, Logger log) {
		Random randomGenerator = new Random();
		String[] randomSubset = new String[num];
		if (a.length == num) {
			randomSubset = a;
		} else if (a.length < num) {
			randomSubset = a;
			log.reportError("Error - cannot choose "	+ num + " random strings from " + a.length
											+ " strings, using all available instead");
		} else {
			Hashtable<Integer, Boolean> tracker = new Hashtable<Integer, Boolean>();
			for (int i = 0; i < num; i++) {
				int index = randomGenerator.nextInt(a.length);
				if (tracker.containsKey(index)) {
					while (tracker.containsKey(index)) {
						index = randomGenerator.nextInt(a.length);
					}
				}
				randomSubset[i] = a[index];
				tracker.put(index, true);
			}
		}
		return randomSubset;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String fileName = "/home/pankrat2/lanej/projects/gedi_exomechip_original.properties";
		String indsToValidateFile = "3050W_unrelateds.txt";
		String dir = "cross_validate/";
		String output = "test_validate";
		String pcMarkers = "/home/pankrat2/shared/gedi_exomechip/targetMarkers.txt";
		String mtMarkers = "annotation_u_use.xln";
		String logfile = "Validation.log";
		int numberOfBatches = 10;// how many groups to break the validation individuals into if running
															// PC files
		int pcReplicates = 3;
		int stopAtComponent = 3049;// number to compute and load
		int startAtComponent = 1;// go from start to stop
		int numPcIterations = 30;// number of samplings across the number of components
		int kfolds = 10;
		int pcType = 1;
		int numThreads = 1;
		boolean batch = false;
		boolean regressionBatch = false;
		boolean svdRegression = true;
		boolean compareAcrossFiles = false;
		boolean recomputeLRR = false;
		Logger log;

		String usage = "\n" + "cnv.analysis.pca.PrincipalComponentsValidate requires 0-1 arguments\n";
		usage += "   For creating jobs to compute PCs:";
		usage += "   (1) Directory (i.e. dir=" + dir + " (default))\n" + "";
		usage +=
					"   (2) Filename of individuals to validate by computing batches of PCs (i.e. indsToValidateFile="
							+ indsToValidateFile + " (default))\n" + "";
		usage += "   (3) Project filename  (i.e. proj=" + fileName + " (default))\n" + "";
		usage += "   (4) PCMarkers to validate on (i.e. pcMarkers=" + pcMarkers + " (default))\n" + "";
		usage += "   (5) MT markers to validate on (i.e. mtMarkers=" + mtMarkers + " (default))\n" + "";
		usage += "   (6) Number of subset batches to compute separate PCs with (i.e. numberOfBatches="
							+ numberOfBatches + " (default))\n" + "";
		usage += "   (7) Number of PC replicates for the batches(i.e. pcReplicates="	+ pcReplicates
							+ " (default))\n" + "";
		usage +=
					"   Note: for each batch the number of PCs computed will be one less than the number of samples";
		usage += "   For cross validating pre-computed PC files:";
		usage += "   (1) project filename  (i.e. proj=" + fileName + " (default))\n" + "";
		usage += "   (2) directory under the project directory with pc files(i.e. dir="	+ dir
							+ " (default))\n" + "";
		usage += "   (3) Max number of PCs to iterate with (i.e. stopAtComponent="	+ stopAtComponent
							+ " (default))\n" + "";
		usage += "   (4) Minimum number of PCs to begin iterations with (i.e. startAtComponent="
							+ startAtComponent + " (default))\n" + "";
		usage += "   (5) Create the batches to run (i.e. -batch (not the default))\n" + "";
		usage += "   (6) Type of PC (i.e. pcType=" + pcType + " (default))\n" + "";
		usage +=
					"   (7) Number of evenly spaced samplings of PCs to compute residuals with (e.g 100 PCs and 10 samplings = PC1,PC10,PC20...PC100) (i.e. numPcIterations="
							+ numPcIterations + " (default))\n" + "";
		usage += "   (8) Number of k-folds to cross validate with  (i.e. kfolds="	+ kfolds
							+ " (default))\n" + "";
		usage += "   (9) Number of threads to use in the regressions (i.e. numThreads="	+ numThreads
							+ " (default))\n" + "";
		usage += "   (10) do not use SVD to compute betas (i.e. -noSVD (not the default))\n" + "";
		usage += "   (11) output base name (i.e. output=" + output + " (default))\n";
		usage +=
					"   (12) along with the default within-file validation, train on each file, and predict the others(i.e -compareAcrossFiles)\n";
		usage += "   (13) batch the cross-validations(i.e -regressionBatch (not the default))\n";
		usage +=
					"   (14) Recompute Log R Ratios for each marker from genotypes/intensities when computing median values(i.e. -recomputeLRR (not the default))\n";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("indsToValidateFile=")) {
				indsToValidateFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("dir=")) {
				dir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("output=")) {
				output = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("pcMarkers=")) {
				pcMarkers = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("mtMarkers=")) {
				mtMarkers = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("proj=")) {
				fileName = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("numberOfBatches=")) {
				numberOfBatches = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("pcReplicates=")) {
				pcReplicates = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("stopAtComponent=")) {
				stopAtComponent = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("startAtComponent=")) {
				startAtComponent = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("pcType=")) {
				pcType = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("numPcIterations=")) {
				numPcIterations = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("kfolds=")) {
				kfolds = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("numThreads=")) {
				numThreads = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("-batch")) {
				batch = true;
				numArgs--;
			} else if (arg.startsWith("-recomputeLRR")) {
				recomputeLRR = true;
				numArgs--;
			} else if (arg.startsWith("-regressionBatch")) {
				regressionBatch = true;
				numArgs--;
			} else if (arg.startsWith("-noSVD")) {
				svdRegression = false;
				numArgs--;
			} else if (arg.startsWith("-compareAcrossFiles")) {
				compareAcrossFiles = true;
				numArgs--;
			} else if (arg.startsWith("log=")) {
				logfile = arg.split("=")[1];
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
			Project proj = new Project(fileName, false);
			log = new Logger(logfile);
			if (batch) {
				batchJobs(proj, dir, indsToValidateFile, numberOfBatches, stopAtComponent, mtMarkers,
									pcMarkers, pcReplicates, log);
			} else if (regressionBatch) {
				analyzeRegressionsBatched(proj, dir, startAtComponent, stopAtComponent, numPcIterations,
																	kfolds, pcType, mtMarkers, numThreads, svdRegression,
																	compareAcrossFiles, recomputeLRR, output, log);
			} else {
				analyzeRegressions(	proj, dir, startAtComponent, stopAtComponent, numPcIterations, kfolds,
														pcType, mtMarkers, numThreads,
														svdRegression ? LS_TYPE.SVD : LS_TYPE.REGULAR, compareAcrossFiles,
														recomputeLRR, output, log);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

