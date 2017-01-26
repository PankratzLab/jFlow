package org.genvisis.cnv.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity.CHROMOSOME_X_STRATEGY;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity.CORRECTION_TYPE;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity.PcCorrectionProducer;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsResiduals;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.ext;
import org.genvisis.stats.LeastSquares.LS_TYPE;

/**
 * This is a temporary fix to export corrected intensities, geared toward PennCNV output and shadow
 * samples. It currently seems to work, but be weary
 *
 */
public class PennCNVPrep {
	private static final String[] PENN_STRINGS = {"Name", ".GType", ".Log R Ratio", ".B Allele Freq"};
	private static final String STORAGE_BASE = "firstMarkerIndex_";
	private static final String STORAGE_EXT = ".ser";

	private static volatile Map<String, MarkerDataStorage> fileToMarkerMap = new HashMap<String, MarkerDataStorage>();

	private final Project proj;
	private final PrincipalComponentsResiduals principalComponentsResiduals;
	private final boolean[] samplesToExport, samplesToUseCluster;
	private final int numCorrectionThreads, numMarkerThreads, numComponents;
	private final int[] sampleSex;
	private final String[] markers;
	private final String dir;
	private final LS_TYPE lType;

	public PennCNVPrep(	Project proj, PrincipalComponentsResiduals principalComponentsResiduals,
											boolean[] samplesToExport, boolean[] samplesToUseCluster, int[] sampleSex,
											String[] markers, int numComponents, String dir, LS_TYPE lType,
											int numThreads, int numMarkerThreads) {
		super();
		this.proj = proj;
		this.principalComponentsResiduals = principalComponentsResiduals;
		this.samplesToExport = samplesToExport;
		this.samplesToUseCluster = samplesToUseCluster;
		this.sampleSex = sampleSex;
		this.numComponents = numComponents;
		numCorrectionThreads = numThreads;
		this.markers = markers;
		this.dir = dir;
		this.numMarkerThreads = numMarkerThreads;
		this.lType = lType;
	}

	/**
	 * This creates the temporary serialized {@link MarkerData} lists stored in
	 * {@link MarkerDataStorage} objects, and contain only LRR/BAF and genotypes (currently original
	 * genotypes)
	 *
	 */
	public void exportSpecialMarkerDataMoreThreads(	String tmpDir, boolean preserveBafs,
																									CORRECTION_TYPE correctionType,
																									CHROMOSOME_X_STRATEGY sexStrategy) {
		String output = (tmpDir == null ? proj.PROJECT_DIRECTORY.getValue() : tmpDir)+ dir
										+ STORAGE_BASE + ext.indexLargeFactors(	markers, proj.getMarkerNames(), true,
																														proj.getLog(), true, true)[0]
										+ STORAGE_EXT;
		if (!Files.exists(output)) {
			new File(ext.parseDirectoryOfFile(output)).mkdirs();
			System.out.println("writing to " + output);
			PcCorrectionProducer producer = new PcCorrectionProducer(	principalComponentsResiduals,
																																numComponents, sampleSex,
																																samplesToUseCluster, lType,
																																numCorrectionThreads, 1, markers,
																																correctionType, sexStrategy);
			proj.getLog().reportTimeInfo("Using correction type " + correctionType);
			WorkerTrain<PrincipalComponentsIntensity> train =
																											new WorkerTrain<PrincipalComponentsIntensity>(producer,
																																																		numMarkerThreads,
																																																		10,
																																																		proj.getLog());
			ArrayList<String> notCorrected = new ArrayList<String>();
			MarkerDataStorage markerDataStorage = new MarkerDataStorage(markers.length);
			int index = 0;
			while (train.hasNext()) {
				MarkerData markerDataToStore = null;
				PrincipalComponentsIntensity principalComponentsIntensity = train.next();
				MarkerData markerData = principalComponentsIntensity.getCentroidCompute().getMarkerData();
				if (principalComponentsIntensity.isFail()) {
					notCorrected.add(markers[index]);
					markerDataToStore = markerData;
				} else {
					byte[] abGenotypes = principalComponentsIntensity.getGenotypesUsed();
					float[][] correctedXY =
																principalComponentsIntensity.getCorrectedIntensity(	PrincipalComponentsIntensity.XY_RETURN,
																																										true);
					float[][] correctedLRRBAF =
																		principalComponentsIntensity.getCorrectedIntensity(	PrincipalComponentsIntensity.BAF_LRR_RETURN,
																																												true);// for
																																															// now
					markerDataToStore = new MarkerData(	markerData.getMarkerName(), markerData.getChr(),
																							markerData.getPosition(), markerData.getFingerprint(),
																							markerData.getGCs(), null, null, correctedXY[0],
																							correctedXY[1], null, null,
																							preserveBafs	? markerData.getBAFs()
																														: correctedLRRBAF[0],
																							correctionType == CORRECTION_TYPE.XY	? correctedLRRBAF[1]
																																										: principalComponentsIntensity.getCorrectedLRR(),
																							abGenotypes, abGenotypes);
				}
				markerDataStorage.addToNextIndex(markerDataToStore);
				index++;
			}
			markerDataStorage.serialize(output);
			if (notCorrected.size() > 0) {
				Files.writeArray(	notCorrected.toArray(new String[notCorrected.size()]),
													output.replaceAll("\\.ser",
																						"_") + notCorrected.size() + "_markersThatFailedCorrection.txt");
			}
		} else {
			proj.getLog().reportFileExists(output);
		}
	}



	/**
	 * This function exports to PennCNV format starting from the temporary files made in
	 * {@link PennCNVPrep#exportSpecialMarkerData()} Currently, all individuals will be exported
	 *
	 * @param fileNamesOfMarkerDataInOrder files of serialized {@link MarkerDataStorage};
	 */
	public void exportSpecialPennCNVData(String[] fileNamesOfMarkerDataInOrder) {
		int[] sampleIndicesInProject = ext.indexLargeFactors(	ArrayUtils.subArray(proj.getSamples(),
																																				samplesToExport),
																													proj.getSamples(), true, proj.getLog(),
																													true, true);
		// int numMarkersPerWrite =
		// Integer.parseInt(proj.getProperty(Project.MAX_MARKERS_LOADED_PER_CYCLE));
		int numMarkersPerWrite = proj.getProperty(proj.MAX_MARKERS_LOADED_PER_CYCLE);
		int numMarkersThisRound = 0;
		String[] subSamples = ArrayUtils.subArray(proj.getSamples(), samplesToExport);
		PennCNVIndividual[] pennCNVIndividuals = initSamples(	proj, dir, subSamples, false,
																													numMarkersPerWrite, proj.getLog());
		for (int i = 0; i < fileNamesOfMarkerDataInOrder.length; i++) {
			MarkerDataStorage markerDataStorage = MarkerDataStorage.load(	fileNamesOfMarkerDataInOrder[i],
																																		false);
			MarkerData[] markerDatas = markerDataStorage.getMarkerDatas();
			for (int j = 0; j < markerDatas.length; j++) {
				addData(numMarkersThisRound, pennCNVIndividuals, markerDatas[j].getBAFs(),
								markerDatas[j].getLRRs(), markerDatas[j].getAbGenotypes(),
								markerDatas[j].getMarkerName(), sampleIndicesInProject);
				numMarkersThisRound++;
				if (numMarkersThisRound == numMarkersPerWrite
						|| (j == markerDatas.length - 1 && i == fileNamesOfMarkerDataInOrder.length - 1)) {
					numMarkersThisRound = 0;
					for (PennCNVIndividual pennCNVIndividual : pennCNVIndividuals) {
						pennCNVIndividual.dump(proj);
					}
					pennCNVIndividuals = initSamples(	proj, dir, subSamples, true, numMarkersPerWrite,
																						proj.getLog());
				}
			}
		}
	}

	/**
	 * Helper method to control access to a volatile map of marker data file names to marker datas.
	 * This method is thread-safe and will ensure only one thread actually loads a given marker data
	 * file. This allows marker data to be shared across threads to prevent unnecessary re-loading.
	 */
	private static MarkerDataStorage loadMarkersIfNeeded(	String markerDataFile, boolean jar,
																												Logger log) {
		MarkerDataStorage markerDataStorage = fileToMarkerMap.get(markerDataFile);
		if (markerDataStorage == null) {
			// Load the MarkerDataStorage for this file in a thread-safe way
			markerDataStorage = loadMarkers(markerDataFile, jar, log);
		}
		return markerDataStorage;
	}

	/**
	 * Synchronized to ensure the actual loading only happens once
	 */
	private static synchronized MarkerDataStorage loadMarkers(String markerDataFile, boolean jar,
																														Logger log) {
		// Double-check to ensure another thread didn't already load this file
		// while we were waiting for a lock
		MarkerDataStorage markerDataStorage;
		if (!fileToMarkerMap.containsKey(markerDataFile)) {
			log.reportTimeInfo("Loading " + markerDataFile);
			markerDataStorage = MarkerDataStorage.load(markerDataFile, jar);
			log.reportTimeInfo("Finished loading " + markerDataFile);
			fileToMarkerMap.put(markerDataFile, markerDataStorage);
		} else {
			// If it was loaded, just look up the result again
			markerDataStorage = fileToMarkerMap.get(markerDataFile);
		}
		return markerDataStorage;
	}

	/**
	 * This function exports to sampraf format starting from the temporary files made in
	 * {@link PennCNVPrep#exportSpecialMarkerData()}
	 *
	 * @param fileNamesOfMarkerDataInOrder files of serialized {@link MarkerDataStorage};
	 * @param forceLoadFromFiles If a large project is being exported, force individual .ser files to
	 *        be loaded each time.
	 */
	public Hashtable<String, Float> exportSpecialSamples(	String[] fileNamesOfMarkerDataInOrder,
																												boolean[] samplesToExport,
																												boolean forceLoadFromFiles) {
		Hashtable<String, Float> allOutliers = new Hashtable<String, Float>();
		int[] subSampleIndicesInProject = ext.indexLargeFactors(
																														ArrayUtils.subArray(	proj.getSamples(),
																																						samplesToExport),
																														proj.getSamples(), true, proj.getLog(),
																														true, true);
		String[] subSamples = ArrayUtils.subArray(proj.getSamples(), samplesToExport);
		String dir = sampleDir(proj);
		proj.getLog().report("Info - checking for existing files in " + dir + "...");
		boolean allExist = true;
		for (int i = 0; i < subSamples.length; i++) {
			if (!Files.exists(dir + subSamples[i] + Sample.SAMPLE_FILE_EXTENSION)) {
				allExist = false;
			}
		}

		if (!allExist) {// if all files exist we skip this export
			proj.getLog().report("Info - detected that not all files exist");

			ShadowSample[] shadowSamples = new ShadowSample[ArrayUtils.booleanArraySum(samplesToExport)];
			for (int i = 0; i < shadowSamples.length; i++) {
				shadowSamples[i] = new ShadowSample(subSamples[i], proj.getMarkerNames());
				if (i % 200 == 0) {
					proj.getLog()
							.report(ext.getTime()+ "\tData loaded = "
											+ Math.round(((double) i / shadowSamples.length) * 100.0)
											+ "%\tFree memory: "
											+ Math.round(((double) MemUtils.availableMem() / Runtime.getRuntime().maxMemory()) * 100.0) + "%");

				}
			}
			int currentIndex = 0;
			for (int i = 0; i < fileNamesOfMarkerDataInOrder.length; i++) {
				MarkerDataStorage markerDataStorage;
				if (forceLoadFromFiles) {
					markerDataStorage = MarkerDataStorage.load(fileNamesOfMarkerDataInOrder[i], false);
				} else {
					markerDataStorage = loadMarkersIfNeeded(fileNamesOfMarkerDataInOrder[i], false,
																									proj.getLog());
				}

				MarkerData[] markerDatas = markerDataStorage.getMarkerDatas();
				for (int j = 0; j < markerDatas.length; j++) {
					if ((j + 1) % 100 == 0) {
						proj.getLog().report("Info - exporting marker "+ (j + 1) + " of " + markerDatas.length
																	+ " from file " + fileNamesOfMarkerDataInOrder[i]);
						if ((j + 1) % 100000 == 0) {

							float usedMemory = Runtime.getRuntime().totalMemory()
																	- Runtime.getRuntime().freeMemory();
							float freeMemory = Runtime.getRuntime().maxMemory() - usedMemory;
							float maxMemory = Runtime.getRuntime().maxMemory();
							proj.getLog().report(ext.getTime()+ "\tData loaded = "
																		+ Math.round(((double) i/ (double) proj.getMarkerNames().length
																									* 100.0))
																		+ "%\tFree memory: "
																		+ Math.round(((double) freeMemory / (double) maxMemory * 100.0))
																		+ "%");
						}
					}
					MarkerData markerData = markerDatas[j];
					for (int j2 = 0; j2 < subSampleIndicesInProject.length; j2++) {
						int sampIndex = subSampleIndicesInProject[j2];
						shadowSamples[j2].addData(proj.getSamples()[sampIndex], currentIndex,
																			markerData.getMarkerName(), markerData.getXs()[sampIndex],
																			markerData.getYs()[sampIndex], markerData.getGCs()[sampIndex],
																			markerData.getBAFs()[sampIndex],
																			markerData.getLRRs()[sampIndex],
																			markerData.getAbGenotypes()[sampIndex], proj.getLog());
					}
					currentIndex++;
				}

			}
			for (ShadowSample shadowSample : shadowSamples) {
				shadowSample.writeShadow(proj, dir, proj.getMarkerSet().getFingerprint(), allOutliers);
			}
		} else {
			proj.getLog()
					.report("Info - detected that all "+ subSamples.length + " shadow samples exist in "
									+ dir + " for the current batch, skipping export...");
		}
		return allOutliers;
	}

	private static void addData(int index, PennCNVIndividual[] pennCNVIndividuals, float[] baf,
															float[] lrr, byte[] genotypes, String marker,
															int[] sampleIndicesInProject) {
		for (int i = 0; i < sampleIndicesInProject.length; i++) {
			pennCNVIndividuals[i].addLine(index, marker, genotypes[sampleIndicesInProject[i]],
																		lrr[sampleIndicesInProject[i]], baf[sampleIndicesInProject[i]]);
		}
	}

	/**
	 * Currently only stores new LRR/BAF values and original genotypes
	 *
	 */
	private static class MarkerDataStorage implements Serializable {
		private static final long serialVersionUID = 1L;
		private final MarkerData[] markerDatas;
		private int currentIndex;

		public MarkerDataStorage(int numMarkers) {
			markerDatas = new MarkerData[numMarkers];
			currentIndex = 0;
		}

		public void addToNextIndex(MarkerData markerData) {
			markerDatas[currentIndex] = markerData;
			currentIndex++;
		}

		public void serialize(String filename) {
			SerializedFiles.writeSerial(this, filename);
		}

		public static MarkerDataStorage load(String filename, boolean jar) {
			return (MarkerDataStorage) SerializedFiles.readSerial(filename, jar, true);
		}

		public MarkerData[] getMarkerDatas() {
			return markerDatas;
		}

	}

	// public Sample(String sampleName, long fingerprint, float[][] data, byte[][] genotypes, boolean
	// canXYBeNegative) {
	// this.sampleName = sampleName;
	// this.fingerprint = fingerprint;
	// this.gcs = data[0];
	// this.xs = data[3];
	// this.ys = data[4];
	// this.thetas = null;
	// this.rs = null;
	// this.bafs = data[7];
	// this.lrrs = data[8];
	// this.forwardGenotypes = genotypes[0];
	// this.abGenotypes = genotypes[1];
	// this.canXYBeNegative = canXYBeNegative;
	// updateNullStatus();
	// }
	public static class ShadowSample {
		private final String sampleName;
		private final String[] shadowMarkers;
		private final float[] shadowGCs;
		private final float[] shadowXs;
		private final float[] shadowYs;
		private final float[] shadowLrrs;
		private final float[] shadowBafs;
		private final byte[] shadowABGenotypes;
		private int numAdded;
		private boolean valid;

		public ShadowSample(String sampleName, String[] markers) {
			this.sampleName = sampleName;
			shadowMarkers = markers;
			shadowGCs = new float[shadowMarkers.length];
			shadowXs = new float[shadowMarkers.length];
			shadowYs = new float[shadowMarkers.length];
			shadowLrrs = new float[shadowMarkers.length];
			shadowBafs = new float[shadowMarkers.length];
			shadowABGenotypes = new byte[shadowMarkers.length];
			numAdded = 0;
			valid = true;
		}

		public void addData(String sample, int index, String marker, float X, float Y, float gc,
												float baf, float lrr, byte abGenotype, Logger log) {
			if (!sample.equals(sampleName)) {
				log.reportError("Error - incorrect sample is being added, trying to add "+ sample
												+ " and should be adding " + sampleName);
				valid = false;
			}
			if (!shadowMarkers[index].equals(marker)) {
				log.reportError("Error - data is not in correct order, got "+ marker
												+ " and should have been " + shadowMarkers[index]);
				valid = false;
			} else if (valid) {
				numAdded++;
				shadowGCs[index] = gc;
				shadowXs[index] = X;
				shadowYs[index] = Y;
				shadowBafs[index] = baf;
				shadowLrrs[index] = lrr;
				shadowABGenotypes[index] = abGenotype;
			}
		}

		public Hashtable<String, Float> writeShadow(Project proj, String dir, long fingerprint,
																								Hashtable<String, Float> allOutliers) {
			if (numAdded != shadowMarkers.length) {
				proj.getLog().reportError("Error not all data was added");
				return null;
			}
			if (!valid) {
				proj.getLog().reportError("Error data was not valid");
				return null;
			}
			new File(dir).mkdirs();
			Sample samp = new Sample(	sampleName, fingerprint, shadowGCs, shadowXs, shadowYs, shadowBafs,
																shadowLrrs, shadowABGenotypes, shadowABGenotypes, false);
			samp.saveToRandomAccessFile(dir+ sampleName + Sample.SAMPLE_FILE_EXTENSION, allOutliers,
																	sampleName);
			return allOutliers;
		}
	}

	/**
	 * Aids in processing to PennCNV output when iterating over markers
	 *
	 */
	private static class PennCNVIndividual {
		private final String ouputFile, header;
		private final String[] markers;
		private final boolean append;
		private final byte[] genotypes;
		private final float[] lrrs;
		private final float[] bafs;
		private int numAdded;

		public PennCNVIndividual(String outputFile, String header, int numMarkers, boolean append) {
			ouputFile = outputFile;
			this.header = header;
			markers = new String[numMarkers];
			genotypes = new byte[numMarkers];
			lrrs = new float[numMarkers];
			bafs = new float[numMarkers];
			this.append = append;
			numAdded = 0;
		}

		public void addLine(int projectIndex, String marker, byte genotype, float lrr, float baf) {
			markers[projectIndex] = marker;
			genotypes[projectIndex] = genotype;
			lrrs[projectIndex] = lrr;
			bafs[projectIndex] = baf;
			numAdded++;
		}

		private String getStringAt(int index) {
			return markers[index]+ "\t"
							+ (genotypes[index] == -1 ? "NC" : Sample.AB_PAIRS[genotypes[index]]) + "\t"
							+ lrrs[index] + "\t" + bafs[index];
		}

		public void dump(Project proj) {
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(ouputFile, append));
				if (!append) {
					writer.println(header);
				}
				for (int i = 0; i < numAdded; i++) {
					writer.println(getStringAt(i));
				}
				writer.close();

			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + ouputFile);
				proj.getLog().reportException(e);
			}
		}
	}

	/**
	 * Grab the {@link PrincipalComponentsResiduals} from {@link Project#INTENSITY_PC_FILENAME}
	 */
	public static PrincipalComponentsResiduals loadPcResids(Project proj, int numComponents) {
		// String pcFile = proj.getFilename(proj.INTENSITY_PC_FILENAME);
		String pcFile = proj.INTENSITY_PC_FILENAME.getValue();
		if (numComponents <= 0) {
			numComponents = proj.getProperty(proj.INTENSITY_PC_NUM_COMPONENTS);
		}
		PrincipalComponentsResiduals pcResids;
		if (Files.exists(proj.PROJECT_DIRECTORY.getValue() + ext.removeDirectoryInfo(pcFile))) {
			proj.getLog().report("Info - loading " + ext.removeDirectoryInfo(pcFile));
			// pcResids = new PrincipalComponentsResiduals(proj, ext.removeDirectoryInfo(pcFile), null,
			// Integer.parseInt(proj.getProperty(Project.INTENSITY_PC_NUM_COMPONENTS)), false, 0, false,
			// false, null);
			pcResids = new PrincipalComponentsResiduals(proj, ext.removeDirectoryInfo(pcFile), null,
																									numComponents, false, 0, false, false, null);
		} else {
			proj.getLog()
					.reportError("Error - did not find Intensity PC File "+ proj.PROJECT_DIRECTORY.getValue()
												+ ext.removeDirectoryInfo(pcFile) + " as defined by"
												+ proj.INTENSITY_PC_FILENAME);
			pcResids = null;
		}
		return pcResids;
	}

	private static PennCNVIndividual[] initSamples(	Project proj, String dir, String[] samplesToExport,
																									boolean append, int numMarkers, Logger log) {
		PennCNVIndividual[] samples = new PennCNVIndividual[samplesToExport.length];
		for (int i = 0; i < samplesToExport.length; i++) {
			String header = PENN_STRINGS[0]+ "\t" + samplesToExport[i] + PENN_STRINGS[1] + "\t"
											+ samplesToExport[i] + PENN_STRINGS[2] + "\t" + samplesToExport[i]
											+ PENN_STRINGS[3];
			samples[i] = new PennCNVIndividual(samplesToExport[i] + ".txt", header, numMarkers, append);
		}
		return samples;
	}

	public static int[] getSampleSex(Project proj) {
		String[] samples = proj.getSamples();
		SampleData sampleData = proj.getSampleData(0, false);
		int[] sex = new int[samples.length];
		for (int i = 0; i < samples.length; i++) {
			sex[i] = sampleData.getSexForIndividual(samples[i]);
		}
		int numMales = ArrayUtils.countIf(sex, 1);
		int numFemales = ArrayUtils.countIf(sex, 2);
		double percentDefined = (double) (numFemales + numMales) / samples.length;
		if (percentDefined > .90) {
			proj.getLog().report("Info - detected "+ numMales + " males and " + numFemales
														+ " females to use for sex-specific reclustering");
		}
		if (numFemales == 0) {
			proj.getLog()
					.report("Error - no females were detected in sample data, chr X will not be reclustered properly");
			return null;
		}
		if (numMales == 0) {
			proj.getLog()
					.report("Error - no males were detected in sample data, chr Y will not be reclustered properly");
			return null;
		}
		return sex;
	}

	private static String sampleDir(Project proj) {
		return ensureExists(proj.PROJECT_DIRECTORY.getValue() + "shadowSamples/");
	}

	private static String transposedDir(Project proj) {
		return ensureExists(proj.PROJECT_DIRECTORY.getValue() + "shadowTransposed/");
	}

	private static String ensureExists(String dir) {
		new File(dir).mkdirs();
		return dir;
	}

	private static String[] getSortedFileNames(Project proj, String dir, String tmpDir) {
		String[] markers = proj.getMarkerNames();
		ArrayList<String> files = new ArrayList<String>();
		int diff = -1;
		for (int i = 0; i < markers.length; i++) {
			String possibleExist = (tmpDir == null ? proj.PROJECT_DIRECTORY.getValue() : tmpDir)+ dir
															+ STORAGE_BASE + i + STORAGE_EXT;
			if (Files.exists(possibleExist)) {

				files.add(possibleExist);
				proj.getLog().reportTimeInfo("Found file "+ possibleExist + " Diff " + (diff - i) + " Num "
																			+ files.size());
				diff = i;
			}
		}
		proj.getLog()
				.report("Info - detected " + files.size() + " files to use for the PennCNV export");
		return files.toArray(new String[files.size()]);
	}

	/**
	 * @param proj
	 * @param dir directory under the project directory
	 * @param numComponents number of components to correct for
	 * @param markerFile path (full path) to a list of markers (single column) to correct for
	 * @param numThreads number of threads for each correction, maximum of 6 will be used
	 * @param exportToPennCNV flag if the directory supplied has serialized files already
	 */
	public static void exportSpecialPennCNV(Project proj, String dir, String tmpDir,
																					int numComponents, String markerFile, int numThreads,
																					int numMarkerThreads,
																					/* boolean exportToPennCNV, */ boolean shadowSamples,
																					LS_TYPE lType, int numSampleChunks, boolean preserveBafs,
																					boolean forceLoadFromFiles,
																					CORRECTION_TYPE correctionType,
																					CHROMOSOME_X_STRATEGY sexStrategy) {
		new File(proj.PROJECT_DIRECTORY.getValue() + dir).mkdirs();
		// if (exportToPennCNV) {
		// boolean[] exportThese = new boolean[proj.getSamples().length];
		// Arrays.fill(exportThese, true);
		// // TODO, samples for Clustering!
		//
		// PennCNVPrep specialPennCNVFormat = new PennCNVPrep(proj, null, exportThese, null, null, null,
		// numComponents, dir, lType, numThreads, numMarkerThreads);
		// String[] sortedFileNames = getSortedFileNames(proj, dir, tmpDir);
		// if (sortedFileNames == null || sortedFileNames.length == 0) {
		// proj.getLog().reportError("Error - did not find any files to export from");
		//
		// } else {
		// specialPennCNVFormat.exportSpecialPennCNVData(sortedFileNames);
		// }
		// }
		if (shadowSamples) {
			if (numSampleChunks == 0) {
				// Conservatively use 75% of available memory for chunking
				double memToUse = 0.75 * MemUtils.availableMem();
				// Estimate worst-case memory use per thread
				long worstCaseSize = Files.worstCaseDirSize(proj.SAMPLE_DATA_FILENAME.getValue(), "sampRAF", false);
				numSampleChunks = (int) Math.round((worstCaseSize / memToUse) / numThreads);
			}


			boolean[][] batches = ArrayUtils.splitUpStringArrayToBoolean(proj.getSamples(), numSampleChunks,
																															proj.getLog());

			ExecutorService executor = Executors.newFixedThreadPool(numThreads);
			Hashtable<String, Future<Hashtable<String, Float>>> tmpResults = new Hashtable<String, Future<Hashtable<String, Float>>>();
			String[] sortedFileNames = getSortedFileNames(proj, dir, tmpDir);
			Hashtable<String, Float> outliers = new Hashtable<String, Float>();
			String outlierFile = sampleDir(proj) + "outliers.ser";

			for (int i = 0; i < batches.length; i++) {

				PennCNVPrep specialPennCNVFormat = new PennCNVPrep(	proj, null, null, null, null, null,
																														numComponents, dir, lType, 1,
																														numMarkerThreads);

				if (sortedFileNames == null || sortedFileNames.length == 0) {
					proj.getLog().reportError("Error - did not find any files to export from");
				} else {
					proj.getLog().reportTimeInfo("Found " + sortedFileNames.length + " special files");
					tmpResults.put(i+ "",
													executor.submit(new WorkerShadow(	specialPennCNVFormat, sortedFileNames,
																														batches[i], i, forceLoadFromFiles,
																														proj.getLog())));
				}
			}

			for (int i = 0; i < batches.length; i++) {
				proj.getLog().reportTimeInfo("Found " + sortedFileNames.length + " special files");
				try {
					if (sortedFileNames == null || sortedFileNames.length == 0) {
						proj.getLog().reportError("Error - did not find any files to export from");
					} else {
						outliers.putAll(tmpResults.get(i + "").get());
					}
				} catch (InterruptedException e) {
					proj.getLog().reportError("Error - interrupted when processing file " + i);
					proj.getLog().reportException(e);
				} catch (ExecutionException e) {
					proj.getLog().reportError("Error - in file on internal index " + i);
					proj.getLog().reportException(e);
				}
			}
			executor.shutdown();
			try {
				executor.awaitTermination(10, TimeUnit.DAYS);
			} catch (InterruptedException e) {
				proj.getLog().reportException(e);
			}
			proj.SAMPLE_DIRECTORY.setValue(sampleDir(proj));
			proj.MARKER_DATA_DIRECTORY.setValue(transposedDir(proj));
			if (outliers.size() == 0) {// usually caused by skipping sample export, so will generate it
				proj.NUM_THREADS.setValue(numMarkerThreads * numThreads);
				proj.verifyAndGenerateOutliers(true);
			} else {
				SerializedFiles.writeSerial(outliers, outlierFile);
			}

			proj.getLog().report("Saving shadow project properties to: "
														+ proj.PROJECT_DIRECTORY.getValue() + "shadow.properties");
			proj.saveProperties(proj.PROJECT_DIRECTORY.getValue() + "shadow.properties");
		} else {
			prepExport(	proj, dir, tmpDir, numComponents, markerFile, numThreads, numMarkerThreads, lType,
									preserveBafs, correctionType, sexStrategy);
		}
	}

	public static void prepExport(Project proj, String dir, String tmpDir, int numComponents,
																String markerFile, int numThreads, int numMarkerThreads,
																LS_TYPE lType, boolean preserveBafs, CORRECTION_TYPE correctionType,
																CHROMOSOME_X_STRATEGY sexStrategy) {
		String[] markers;
		PrincipalComponentsResiduals principalComponentsResiduals = loadPcResids(proj, numComponents);
		if (principalComponentsResiduals == null) {
			proj.getLog().reportError("Error: no principal component residuals for project");
			return;
		}
		int[] sex = getSampleSex(proj);
		if (sex == null && proj.getAutosomalMarkers().length != proj.getMarkerNames().length) {
			proj.getLog().reportTimeWarning("missing sex codes");
			// return;
		}
		if (markerFile == null) {
			markers = proj.getMarkerNames();
			proj.getLog()
					.report("Info - a file of markers was not provided, exporting all in this batch");
		} else {
			markers = HashVec.loadFileToStringArray(markerFile, false, new int[] {0}, false);
			proj.getLog()
					.report("Info - loaded " + markers.length + " markers from " + markerFile + " to export");

		}
		PennCNVPrep specialPennCNVFormat = new PennCNVPrep(	proj, principalComponentsResiduals, null,
																												proj.getSamplesToInclude(null), sex,
																												markers, numComponents, dir, lType,
																												numThreads, numMarkerThreads);
		specialPennCNVFormat.exportSpecialMarkerDataMoreThreads(tmpDir, preserveBafs, correctionType,
																														sexStrategy);
	}

	public static void batchCorrections(Project proj, String java, String classPath, int memoryInMB,
																			int wallTimeInHours, String dir, String tmpDir,
																			int numBatches, int numThreads, int numMarkerThreads,
																			int numComponents) {
		String[] allMarkers = proj.getMarkerNames();
		int[] chunks = ArrayUtils.splitUp(allMarkers.length, numBatches);
		int index = 0;
		String[][] batches = new String[numBatches][1];
		String thisDir = (tmpDir == null ? proj.PROJECT_DIRECTORY.getValue() : tmpDir);
		new File(thisDir + dir).mkdirs();
		for (int i = 0; i < chunks.length; i++) {
			ArrayList<String> chunk = new ArrayList<String>(chunks[i]);
			for (int j = 0; j < chunks[i]; j++) {
				chunk.add(allMarkers[index]);
				index++;
			}
			batches[i][0] = "batch_" + i + "_" + chunks[i] + "_markers";
			Files.writeArray(	chunk.toArray(new String[chunk.size()]),
												thisDir + dir + batches[i][0] + ".txt");
		}
		StringBuilder cmd = new StringBuilder("module load java\n");
		cmd	.append("java").append(" -Xmx").append(memoryInMB).append("M -jar ").append(classPath)
				.append(" cnv.analysis.PennCNVPrep proj=").append(proj.getPropertyFilename())
				.append(" dir=").append(dir);
		Files.qsub(	"PennCNVPrepFormatExport", cmd.toString() + " -create", new String[][] {{""}},
								memoryInMB, 3 * wallTimeInHours, 1);
		cmd.append(" tmpDir=").append(thisDir);
		Files.qsub(	"ShadowCNVPrepFormatExport",
								cmd.toString() + " -shadow sampleChunks=NeedToFillThisIn numThreads=1 -forceLoadFromFiles",
								new String[][] {{""}}, memoryInMB, 3 * wallTimeInHours, 1);
		cmd	.append(" numMarkerThreads=").append(numMarkerThreads).append(" numThreads=")
				.append(numThreads).append(" numComponents=").append(numComponents).append(" markers=")
				.append(thisDir).append(dir).append("[%0].txt");
		Files.qsub(	"PennCNVPrepFormatTmpFiles", cmd.toString(), batches, memoryInMB, wallTimeInHours,
								numThreads * numMarkerThreads);
		if (!Files.exists(proj.INTENSITY_PC_FILENAME.getValue())) {
			proj.getLog()
					.report("Warning - all jobs will fail if the property "
									+ proj.INTENSITY_PC_FILENAME.getName() + " in " + proj.getPropertyFilename()
									+ " is not set to an existing file");
			proj.getLog().report("		  - did not find " + proj.INTENSITY_PC_FILENAME.getValue());
		}
		if (getSampleSex(proj) == null) {
			proj.getLog().report("Warning - all jobs will fail if sample sex is not provided in "
														+ proj.SAMPLE_DATA_FILENAME.getValue());
			proj.getLog().report("		  - please specify sex for as many individuals as possible");
		}
	}

	private static class WorkerShadow implements Callable<Hashtable<String, Float>> {

		private final PennCNVPrep specialPennCNVFormat;
		private final String[] sortedFileNames;
		private final boolean[] batch;
		private final boolean forceLoadFromFiles;
		private final int batchIndex;
		private final Logger log;

		public WorkerShadow(PennCNVPrep specialPennCNVFormat, String[] sortedFileNames, boolean[] batch,
												int batchIndex, boolean forceLoadFromFiles, Logger log) {
			super();
			this.specialPennCNVFormat = specialPennCNVFormat;
			this.sortedFileNames = sortedFileNames;
			this.batch = batch;
			this.forceLoadFromFiles = forceLoadFromFiles;
			this.batchIndex = batchIndex;
			this.log = log;
		}

		@Override
		public Hashtable<String, Float> call() {
			log.report("Info - exporting batch "+ batchIndex + " with thread "
									+ Thread.currentThread().getName());

			Hashtable<String, Float> outliers =
																				specialPennCNVFormat.exportSpecialSamples(sortedFileNames,
																																									batch,
																																									forceLoadFromFiles);
			return outliers;
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String classPath = "/home/pankrat2/coleb/" + org.genvisis.common.PSF.Java.GENVISIS;
		String java = "/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/bin/java";// for lab
		// String java ="/soft/java/jdk1.7.0_45/bin/java"; //for itasca
		int memoryInMB = 22000;
		int wallTimeInHours = 15;
		String filename = null;
		String logfile = null;
		String dir = "PennCNVPrep/";
		String tmpDir = null;
		int numThreads = 6;// can only utilize 6 (3 genotype clusters by X/Y)
		int numMarkerThreads = 2;
		int numComponents = 40;
		String markers = null;
		boolean shadowSamples = false;
		boolean svdRegression = false;
		int batch = 0;
		int sampleChunks = 0;
		boolean forceLoadFromFiles = false;
		CHROMOSOME_X_STRATEGY strategy = CHROMOSOME_X_STRATEGY.BIOLOGICAL;
		// Ex - Recommend modifying this to run the corrections

		// java -jar " + common.PSF.Java.GENVISIS + " cnv.analysis.PennCNVPrep batch=100
		// proj=/home/usr/projects/x.properties classPath=/yourPathTo/" + common.PSF.Java.GENVISIS + "
		// dir=PennCNVPrep/ numComponents=40

		// then run recommended
		// ./master.PennCNVPrepFormatTmpFiles

		// and after completion of all batches
		// ./master.PennCNVPrepFormatExport

		// also, if you submit a bunch of jobs, this will kill them all
		// showq -u Your_Usr_Name_Here |cut -d ' ' -f 1|xargs qdel

		String usage = "\n" + "cnv.analysis.PennCNVPrep requires 1 argument\n";
		usage += "   (1) Project (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) logfile (i.e. log=" + logfile + " ( no default))\n" + "";
		usage += "   (3) directory (relative to the project directory) for output (i.e. dir="+ dir
							+ " ( no default))\n" + "";
		usage += "   (6) number of principal components for correction, (i.e. numComponents="
							+ numComponents + " (default))\n" + "";
		usage +=
					"   (7) a full path to a file listing markers to export in the current batch, (i.e. markers="
							+ numComponents + " (default))\n" + "";
		usage +=
					"   (8) create PennCNV files from the tempory markerData files, (i.e. -create ( not the default))\n"
							+ "";
		usage +=
					"   (9) set this up for a batch run, which is recommended. Set to 0 if batch is not wanted (i.e. batch="
							+ batch + " (default))\n" + "";
		// usage += " (9) java location for batch run (i.e. java=" + java + " (default))\n" + "";
		usage += "   (10) classPath for batch run (i.e. classPath=" + classPath + " (default))\n" + "";
		usage +=
					"   (12) export shadow samples for quickly comparing the current sample data to corrected data(i.e. -shadow ( not the default))\n"
							+ "";
		usage +=
					"   (12) if exporting shadow samples, the number of chunks to export at once (this many samples*the number of threads will be held in memory). A value <= 0 will make chunking determined automatically. (i.e. sampleChunks="
							+ sampleChunks + " (default))\n" + "";
		usage += "   (13) walltime in hours for batched run (i.e. walltime="+ wallTimeInHours
							+ " (default))\n" + "";
		usage +=
					"   (14) memory in mb for batched run (i.e. memory=" + memoryInMB + " (default))\n" + "";
		usage +=
					"   (15) if using a large number of PCs (>150) use a svd regression method (i.e. -svd (not the default))\n"
							+ "";
		usage +=
					"   (16) number of threads for a single marker (correction within a marker) (i.e. numThreads="
							+ numThreads + " (default))\n" + "";
		usage +=
					"   (17) number of threads for between a marker  (correction between a marker)(i.e. numMarkerThreads="
							+ numThreads + " (default))\n" + "";
		usage += "   (18) full path to a temporary directory (i.e. tmpDir= (no default))\n" + "";
		usage += "   (7) Chromosome X correction strategy.  Options include: "
							+ ArrayUtils.toStr(CHROMOSOME_X_STRATEGY.values(), ", ") + " (i.e. sexStrategy=" + strategy
							+ " (default))\n";
		usage += "   NOTE: the total number of threads is numThreads*numMarkerThreads";
		usage += "   NOTE: aprox 50 *(numSamples/5000) batches per 500,000 markers" + "";
		usage +=
					"   NOTE: If using ChrX/ChrY markers, it is important that the projects sample data file has sex defined for all samples that are used for clustering"
							+ "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("proj=")) {
				filename = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("log=")) {
				logfile = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("dir=")) {
				dir = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("tmpDir=")) {
				tmpDir = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("numThreads=")) {
				numThreads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("numMarkerThreads=")) {
				numMarkerThreads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("numComponents=")) {
				numComponents = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("sampleChunks=")) {
				sampleChunks = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("markers=")) {
				markers = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("java=")) {
				java = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("classPath=")) {
				classPath = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("-create")) {
				numArgs--;
			} else if (arg.startsWith("-shadow")) {
				shadowSamples = true;
				numArgs--;
			} else if (arg.startsWith("-svd")) {
				svdRegression = true;
				numArgs--;
			} else if (arg.startsWith("-forceLoadFromFiles")) {
				forceLoadFromFiles = true;
				numArgs--;
			} else if (arg.startsWith("batch=")) {
				batch = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("walltime=")) {
				wallTimeInHours = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("sexStrategy=")) {
				strategy = CHROMOSOME_X_STRATEGY.valueOf(ext.parseStringArg(arg, strategy.toString()));
				numArgs--;
			} else if (arg.startsWith("memory=")) {
				memoryInMB = ext.parseIntArg(arg);
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
			Project proj = new Project(filename, logfile, false);
			if (batch > 0) {
				batchCorrections(	proj, java, classPath, memoryInMB, wallTimeInHours, dir, tmpDir, batch,
													numThreads, numMarkerThreads, numComponents);
			} else {
				exportSpecialPennCNV(	proj, dir, tmpDir, numComponents, markers, numThreads,
															numMarkerThreads, shadowSamples,
															svdRegression ? LS_TYPE.SVD : LS_TYPE.REGULAR, sampleChunks, false,
															forceLoadFromFiles, CORRECTION_TYPE.XY,
															strategy);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
