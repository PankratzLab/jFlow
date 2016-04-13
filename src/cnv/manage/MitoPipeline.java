package cnv.manage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;

import cnv.LaunchProperties;
import cnv.analysis.PennCNV;
import cnv.analysis.pca.PCA;
import cnv.analysis.pca.PrincipalComponentsApply;
import cnv.analysis.pca.PrincipalComponentsCompute;
import cnv.analysis.pca.PrincipalComponentsResiduals;
import cnv.filesys.ABLookup;
import cnv.filesys.MarkerLookup;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.filesys.SampleList;
import cnv.manage.Resources.GENOME_BUILD;
import cnv.manage.Resources.Resource;
import cnv.qc.GcAdjustor.GCAdjustorBuilder;
import cnv.qc.GcAdjustorParameter.GcAdjustorParameters;
import cnv.qc.GcAdjustor;
import cnv.qc.GcAdjustorParameter;
import cnv.qc.MarkerMetrics;
import cnv.var.SampleData;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;

/**
 * A class that serves as the outer wrapper for PCA related happenings...does import, marker QC, sample QC, PCA, etc
 * 
 */
public class MitoPipeline {
	public static final String FILE_BASE = "PCA_GENVISIS";
	public static final String[] PED_INPUT = { "DNA", "FID", "IID", "FA", "MO", "SEX", "AFF" };
	public static final String[] SAMPLEMAP_INPUT = { "Index", "Name", "ID", "Gender", "Plate", "Well", "Group", "Parent1", "Parent2", "Replicate", "SentrixPosition" };
	public static final String[] QC_COLUMNS = { "Sample", "LRR_SD", "Genotype_callrate" };
	public static final String[] SAMPLE_QC_SUMMARY = { "DNA", "LRR_SD", "Genotype_callrate", "Included in PC?" };
	public static final String[] SEX = { "female", "male" };
	public static final String[] SAMPLE_DATA_ADDITION_HEADERS = { "LRR_SD", "Genotype_callrate", "CLASS=Exclude" };

	static final String MARKERS_TO_QC_FILE = "markers_to_QC.txt";
	static final String MARKERS_FOR_ABCALLRATE = "markers_ABCallRate.txt";
	public static final String DNA_LINKER = "DNA";
	public static final String PCA_SAMPLES = ".samples.USED_PC.txt";
	public static final String PCA_SAMPLES_SUMMARY = ".samples.QC_Summary.txt";
	public static final String PROJECT_EXT = ".properties";
	private static final String PCA_FINAL_REPORT = ".finalReport.txt";
	private static final String PC_MARKER_COMMAND = "PCmarkers=";
	private static final String MITO_MARKER_COMMAND = "mitochondrialMarkers=";
	private static final String USE_FILE_COMMAND = "useFile=";
	private static final String IMPORT_EXTENSION = "dirExt=";
	private static final long RECOMMENDED_MEMORY = 1000000000;

	private String projectName;
	private String filename;
	private Project proj;
	private Logger log;
	private String projectDirectory, sourceDirectory, dataExtension, defaultLRRSdFilter, defaultCallRateFilter, targetMarkers, medianMarkers, markerPositions, idHeader, abLookup;

	/**
	 * 
	 * 
	 * @param projectName
	 *            future name of the project
	 * @param projectDirectory
	 *            future project directory
	 * @param sourceDirectory
	 *            where the data is
	 * @param dataExtension
	 * @param idHeader
	 *            header of the sample name in the data files
	 * @param abLookup
	 *            if necessary
	 * @param defaultLRRSdFilter
	 *            to perform sample qc
	 * @param defaultCallRateFilter
	 *            to perform sample qc
	 * @param targetMarkers
	 *            markers that will be used in PC computation
	 * @param markerPositions
	 *            positions of all markers in the project
	 * @param medianMarkers
	 *            markers that will have their median values summarized. PCs will be regressed out of the median values
	 * @param log
	 */
	public MitoPipeline(Project existingProj, String projectName, String projectDirectory, String sourceDirectory, String dataExtension, String idHeader, String abLookup, String defaultLRRSdFilter, String defaultCallRateFilter, String targetMarkers, String markerPositions, String medianMarkers, String logfile) {
		String path = initGenvisisProject();
		this.projectName = projectName;
		this.filename = path + projectName + PROJECT_EXT;
		this.projectDirectory = projectDirectory;
		this.sourceDirectory = sourceDirectory;
		this.dataExtension = dataExtension;
		this.idHeader = idHeader;
		this.abLookup = abLookup;
		this.defaultLRRSdFilter = defaultLRRSdFilter;
		this.defaultCallRateFilter = defaultCallRateFilter;
		this.targetMarkers = targetMarkers;
		this.markerPositions = markerPositions;
		this.medianMarkers = medianMarkers;
		this.log = new Logger(logfile);
		initProjectDir();
		System.out.println(path);
		if (existingProj == null) {
			initProject(path);
		} else {
			this.proj = existingProj;
		}
		if (logfile != null) {
			proj.setLog(log);
		} else {
			log = proj.getLog();
		}
		initMainProperties();
		copyAuxFiles();
		initDataDir();
	}

	/**
	 * Sets up the location for projects
	 */
	public static String initGenvisisProject() {
		String launchPropertiesFile = LaunchProperties.DEFAULT_PROPERTIES_FILE;
		String path = LaunchProperties.directoryOfLaunchProperties(launchPropertiesFile);
		path = (new LaunchProperties(path + launchPropertiesFile)).getDirectory();
		// if (!new File(path + "projects/").exists()) {
		// new File(path + "projects/").mkdirs();
		// }
		if (!new File(path).exists()) {
			new File(path).mkdirs();
		}
		if (!new File(launchPropertiesFile).exists()) {
			new File(path + "example/").mkdirs();
			Files.writeList(new String[] { "LAST_PROJECT_OPENED=example.properties", "PROJECTS_DIR=" + path }, launchPropertiesFile);
			if (!new File(path + "example.properties").exists()) {
				Files.writeList(new String[] { "PROJECT_NAME=Example", "PROJECT_DIRECTORY=example/", "SOURCE_DIRECTORY=sourceFiles/" }, path + "example.properties");
			}
		}
		return path;
		// if (!new File(launchPropertiesFile).exists()) {
		// new File(path + "example/").mkdirs();
		// Files.writeList(new String[] { "LAST_PROJECT_OPENED=example.properties", "PROJECTS_DIR=" + path + "projects/" }, launchPropertiesFile);
		// if (!new File(path + "projects/example.properties").exists()) {
		// Files.writeList(new String[] { "PROJECT_NAME=Example", "PROJECT_DIRECTORY=example/", "SOURCE_DIRECTORY=sourceFiles/" }, path + "projects/example.properties");
		// }
		// }
		// return path + "projects/";
	}

	/**
	 * Copies the default project to the project directory if the desired fileName does not already exist
	 */
	public void initProject(String path) {
		if (Files.exists(filename)) {
			Files.backup(ext.removeDirectoryInfo(filename), path, path + "backup/", false);
			log.report("Using project file " + filename + ", you may also specify project filename using the command line argument \"proj=\"");
		} else {
			// if (proj != null) {
			//
			// }
			log.report("Project properties file can be found at " + filename);
			Files.write("PROJECT_NAME=" + projectName, filename);
		}
		this.proj = new Project(filename, null, false, false);
	}

	public void initProjectDir() {
		mkdir(projectDirectory, log);
	}

	public void initDataDir() {
		mkdir(proj.getProperty(proj.DATA_DIRECTORY), log);
	}

	public Project getProj() {
		return proj;
	}

	/**
	 * Copy all auxiliary files to project directory if not already present
	 */
	public void copyAuxFiles() {
		boolean missingFile;

		missingFile = false;
		log.report("Preparing auxiliary files in " + projectDirectory);
		if (!Files.exists(projectDirectory + ext.removeDirectoryInfo(targetMarkers))) {
			if (!Files.copyFile(targetMarkers, projectDirectory + ext.removeDirectoryInfo(targetMarkers))) {
				log.reportError("Error - the filename specified for targetMarkers (\"" + targetMarkers + "\") was not found");
				missingFile = true;
			}
		}
		if (markerPositions == null) {
			log.report("Info - markerPositions was set to null, so attempting to auto-create from SNP_Map.csv");
			String snpMapFilename = proj.getLocationOfSNP_Map(true);
			if (snpMapFilename == null) {
				log.reportError("Error - unfortunately a SNP_Map.csv file could not be found; either create a markerPositions.txt file with MarkerName/Chr/Position, or supply a SNP_Map.csv file");
				log.reportError("aborting MitoPipeline");
				System.exit(1);
			} else {
				log.report("\nParsing " + proj.MARKER_POSITION_FILENAME.getValue(false, false) + " using " + snpMapFilename);
				cnv.manage.Markers.generateMarkerPositions(proj, snpMapFilename);
			}
		} else if (!Files.exists(projectDirectory + ext.removeDirectoryInfo(markerPositions))) {
			if (!Files.copyFile(markerPositions, projectDirectory + ext.removeDirectoryInfo(markerPositions))) {
				log.reportError("Error - the filename specified for markerPositions (\"" + markerPositions + "\") was not found");
				missingFile = true;
			}
		}
		if (!Files.exists(projectDirectory + ext.removeDirectoryInfo(medianMarkers))) {
			if (!Files.copyFile(medianMarkers, projectDirectory + ext.removeDirectoryInfo(medianMarkers))) {
				log.reportError("Error - the filename specified for medianMarkers (\"" + medianMarkers + "\") was not found");
				missingFile = true;
			}
		}
		if (abLookup != null && Files.exists(abLookup) && !Files.exists(projectDirectory + abLookup)) {
			if (!Files.copyFile(abLookup, projectDirectory + ext.removeDirectoryInfo(abLookup))) {
				log.reportError("Error - the filename specified for abLookup (\"" + abLookup + "\") was not found");
				missingFile = true;
			}
		}

		if (missingFile) {
			log.reportError("Error - critical file missing; aborting MitoPipeline");
			System.exit(1);
		}
	}

	/**
	 * Sets the project properties as defined at the command line
	 */
	public void initMainProperties() {
		proj.setProperty(proj.PROJECT_NAME, projectName);
		proj.setProperty(proj.PROJECT_DIRECTORY, projectDirectory);
		proj.setProperty(proj.SOURCE_DIRECTORY, sourceDirectory);
		proj.setProperty(proj.SOURCE_FILENAME_EXTENSION, dataExtension);
		proj.setProperty(proj.ID_HEADER, idHeader);
		// proj.setProperty(proj.LRRSD_CUTOFF, defaultLRRSdFilter);
		proj.setProperty(proj.LRRSD_CUTOFF, Double.valueOf(defaultLRRSdFilter));
		proj.setProperty(proj.SAMPLE_CALLRATE_THRESHOLD, Double.valueOf(defaultCallRateFilter));
		proj.setProperty(proj.INTENSITY_PC_MARKERS_FILENAME, ext.removeDirectoryInfo(targetMarkers));
		if (markerPositions != null) {
			proj.setProperty(proj.MARKER_POSITION_FILENAME, ext.removeDirectoryInfo(markerPositions));
		}
		if (abLookup != null && Files.exists(projectDirectory + abLookup)) {
			proj.setProperty(proj.AB_LOOKUP_FILENAME, ext.removeDirectoryInfo(abLookup));
		}
		proj.saveProperties();
	}

	public static Project initializeProject(Project proj, String projectName, String projectDirectory, String sourceDirectory, String dataExtension, String idHeader, String abLookup, String targetMarkers, String medianMarkers, String markerPositions, String defaultLRRSdFilter, String defaultCallRateFilter, String logfile) {
		MitoPipeline projG = new MitoPipeline(proj, projectName, projectDirectory, sourceDirectory, dataExtension, idHeader, abLookup, defaultLRRSdFilter, defaultCallRateFilter, targetMarkers, markerPositions, medianMarkers, logfile);
		return projG.getProj();
	}

	/**
	 * The main event. Takes the samples from raw data through import and PCA
	 */

	public static int catAndCaboodle(Project proj, int numThreads, String medianMarkers, int numComponents, String outputBase, boolean homosygousOnly, boolean markerQC, double markerCallRateFilter, String useFile, String pedFile, String sampleMapCsv, boolean recomputeLRR_PCs, boolean recomputeLRR_Median, boolean sampLrr, boolean doAbLookup, boolean imputeMeanForNaN, boolean gcCorrect, String refGenomeFasta, int bpGcModel, int regressionDistance, GENOME_BUILD build) {
		String sampleDirectory;
		SampleList sampleList;
		int[] counts;
		Logger log;
		long memoryAvailable;
		int result;
		log = proj.getLog();

		memoryAvailable = Runtime.getRuntime().maxMemory();
		log.report("Memory available = " + memoryAvailable + "  (" + ext.prettyUpSize(memoryAvailable, 1) + ")");
		if (memoryAvailable < RECOMMENDED_MEMORY) {
			log.reportError("\nWarning - " + ext.prettyUpSize(memoryAvailable, 1) + " may not be enough RAM to get the job done; add the following -Xmx argument at the command line to increase the amount of memory to the Java virtual machine");
			log.reportError("java -Xmx10g -cp /path/to/genvisis.jar ... ");
			log.reportError("which will allocate 10 Gb of RAM (likewise, you can set it to -Xmx2g for 2 GB or -Xmx250 for 250Gb)\n");
		}

		sampleDirectory = proj.SAMPLE_DIRECTORY.getValue(false, false);
		if (Files.exists(sampleDirectory) && Files.list(sampleDirectory, Sample.SAMPLE_DATA_FILE_EXTENSION, false).length > 0 && proj.getSampleList() != null && proj.getSampleList().getSamples().length > 0) {
			sampleList = proj.getSampleList();
			log.report("Detected that " + (sampleList.getSamples().length > 1 ? sampleList.getSamples().length + " samples have" : sampleList.getSamples().length + " sample has") + " already been parsed");
			// log.report("Skipping sample import step for the analysis. If this is an incorrect number of samples, please remove (or change the name of) " + proj.getFilename(proj.SAMPLELIST_FILENAME) + " and " + proj.getDir(proj.SAMPLE_DIRECTORY));
			log.report("Skipping sample import step for the analysis. If this is an incorrect number of samples, please remove (or change the name of) " + proj.SAMPLELIST_FILENAME.getValue() + " and " + proj.SAMPLE_DIRECTORY.getValue(false, true));
		} else {
			result = cnv.manage.SourceFileParser.createFiles(proj, numThreads);
			if (result == 0) {
				return 0;
			} else if (result == 7) {
				log.reportError("\nAlternatively, mitoPipeline can do this for you, if you set the markerPositions argument to be null or blank (i.e., \"markerPositions=\")");
				return 0;
			} else if (result == 6) {
				doAbLookup = true;
			}
			if (Files.exists(sampleDirectory) && (Files.list(sampleDirectory, null, false).length == 0)) {
				log.reportError("\nMake sure your " + IMPORT_EXTENSION + " argument is set to the right file extension");
			}
		}
		try {
			Thread.sleep(5000); // Got hit with the error below for no reason twice now
		} catch (InterruptedException ie) {
		}
		sampleList = proj.getSampleList();
		if (sampleList == null || sampleList.getSamples().length == 0) {
			log.report("\n" + ext.getTime() + "\tError - samples were not imported properly, halting MitoPipeline");
			if (doAbLookup) {
				return 40;// we return 40 so that the next attempt will remember to create an ab Lookup
			} else {
				return 41;
			}
		}

		SampleData.createSampleData(pedFile, sampleMapCsv, proj);
		// we require that every sample that has been parsed has an entry in sampleData
		if (verifyAllSamples(proj, sampleList.getSamples())) {
			if (doAbLookup) {
				log.report("Info - determined that an AB lookup is required and was not provided, attempting to generate one now");
				if (generateABLookup(proj, log) == 0) {
					return 0;
				}
			}
			// if a useFile is given, all samples must be available
			if (verifyUseFile(proj, sampleList.getSamples(), useFile)) {
				if (new File(proj.MARKER_DATA_DIRECTORY.getValue(false, false) + "markers.0.mdRAF").exists()) {
					log.report("Marker data (at least the first file 'markers.0.mdRAF') have already been parsed");
					log.report("Skipping transpose step for the analysis. If you would like to re-transpose the data, please remove (or change the name of) " + proj.MARKER_DATA_DIRECTORY.getValue(false, false));
				} else {
					TransposeData.transposeData(proj, 2000000000, false);
				}

				// we make sure each marker has an entry in the projects Markerlookup. I am doing this in case previous steps have already failed, and this should catch it
				if (verifyAllProjectMarkersAreAvailable(proj)) {
					String markersForABCallRate = null;
					String markersForEverythingElse = null;
					// check that all target markers are available
					// if (verifyAuxMarkers(proj, proj.getFilename(proj.TARGET_MARKERS_FILENAME), PC_MARKER_COMMAND)) {
					if (verifyAuxMarkers(proj, proj.INTENSITY_PC_MARKERS_FILENAME.getValue(), PC_MARKER_COMMAND)) {
						// if marker QC is not flagged, sample qc is based on all target markers by default
						String markersToQC = proj.PROJECT_DIRECTORY.getValue() + outputBase + "_" + MARKERS_TO_QC_FILE;
						String markersABCallrate = proj.PROJECT_DIRECTORY.getValue() + outputBase + "_" + MARKERS_FOR_ABCALLRATE;
						String baseLineMarkers = proj.PROJECT_DIRECTORY.getValue() + outputBase + "_baselineMarkers.txt";
						String[] auto = proj.getAutosomalMarkers();
						ArrayList<String> tmp = new ArrayList<String>();
						for (int i = 0; i < auto.length; i++) {
							if (!proj.getArrayType().isCNOnly(auto[i])) {
								tmp.add(auto[i]);
							}
						}
						Files.writeArrayList(tmp, baseLineMarkers);
						if (markerQC) {
							String markerQCFile = outputBase + "_markerQC.txt";
							proj.MARKER_METRICS_FILENAME.setValue(markerQCFile);
							qcMarkers(proj, baseLineMarkers, markersToQC, markersABCallrate, markerCallRateFilter, numThreads);
							markersForABCallRate = markersABCallrate;
							if (!Files.exists(markersForABCallRate)) {
								log.reportError("Error - markerQC was flagged but the file " + markersABCallrate + " could not be found");
								return 1;
							}
						} else {
							markersForABCallRate = baseLineMarkers;
							writeMarkersToQC(proj, baseLineMarkers, markersToQC);
						}
						markersForEverythingElse = markersToQC;
						String qcFile = outputBase + "_lrr_sd.txt";
						proj.SAMPLE_QC_FILENAME.setValue(qcFile);

						counts = cnv.qc.LrrSd.filterSamples(proj, outputBase, markersForABCallRate, markersForEverythingElse, numThreads, useFile,false);
						if (counts == null || counts[1] != sampleList.getSamples().length) {
							if (counts == null || counts[1] == 0 && Files.exists(proj.SAMPLE_QC_FILENAME.getValue())) {
								log.reportError("Error - was unable to parse QC file " + proj.SAMPLE_QC_FILENAME.getValue() + ", backing up this file to " + proj.BACKUP_DIRECTORY.getValue(false, false) + " and re-starting sample qc");
								Files.backup(ext.removeDirectoryInfo(proj.SAMPLE_QC_FILENAME.getValue()), proj.PROJECT_DIRECTORY.getValue(), proj.BACKUP_DIRECTORY.getValue(true, false), true);
							}
							counts = cnv.qc.LrrSd.filterSamples(proj, outputBase, markersForABCallRate, markersForEverythingElse, numThreads, useFile, false);
							if (counts == null || counts[1] != sampleList.getSamples().length) {
								if (counts == null) {
									log.reportError("Error - could not parse QC file (" + proj.SAMPLE_QC_FILENAME.getValue() + ")");
								} else {
									log.reportError("Error - different number of samples (n=" + counts[1] + ") listed in the QC file (" + proj.SAMPLE_QC_FILENAME.getValue() + ") compared to the number of samples in the project (n=" + sampleList.getSamples().length + ")");
									log.reportError("      - delete the QC file (" + proj.SAMPLE_QC_FILENAME.getValue() + ") to regenerate it with the correct number of samples");
								}
								log.reportError("aborting...");
								return 2;
							}
						}
						if (counts == null || counts[0] == 0) {// no samples passed threshold, null case shouldn't happen but we will test anyway
							return 2;// message handled already
						}
						// check that all median markers are available
						if (verifyAuxMarkers(proj, medianMarkers, MITO_MARKER_COMMAND)) {
							// compute PCs with samples passing QC

							GcAdjustorParameters params = null;
							if (gcCorrect) {// TODO, non gc sampleSpecific recomputing
								String samps = proj.PROJECT_DIRECTORY.getValue() + outputBase + PCA_SAMPLES;
								boolean[] sampsToUseRecompute = null;

								if (sampLrr) {
									if (!recomputeLRR_Median || !recomputeLRR_PCs) {
										proj.getLog().reportTimeWarning("Sample specific LRR flagged, overiding other recompute LRR options");
										recomputeLRR_Median = true;
										recomputeLRR_PCs = true;
									}
									sampsToUseRecompute = Array.booleanArray(proj.getSamples().length, false);
									int[] indices = ext.indexFactors(HashVec.loadFileToStringArray(samps, false, new int[] { 0 }, false), proj.getSamples(), true, false);
									for (int i = 0; i < indices.length; i++) {
										sampsToUseRecompute[indices[i]] = true;
									}
									proj.getLog().reportTimeInfo("LRR will be recomputed with " + Array.booleanArraySum(sampsToUseRecompute) + " samples from " + samps);
								}
								if ((refGenomeFasta != null && !Files.exists(refGenomeFasta)) && Files.exists(proj.REFERENCE_GENOME_FASTA_FILENAME.getValue())) {
									proj.getLog().reportTimeWarning("Command line reference genome did not exist or was not provided, using default " + proj.REFERENCE_GENOME_FASTA_FILENAME.getValue());
									refGenomeFasta = proj.REFERENCE_GENOME_FASTA_FILENAME.getValue();
								}
								Resource gmodelBase = Resources.getGCBaseResource(build);
								if (!Files.exists(proj.GC_MODEL_FILENAME.getValue()) && (refGenomeFasta == null || !Files.exists(refGenomeFasta)) && gmodelBase.isAvailable()) {
									log.reportTimeWarning("Generating gcModel for " + build.getBuild() + " at " + proj.GC_MODEL_FILENAME.getValue() + " from " + gmodelBase.getResource());
									proj.getLog().setLevel(3);
									PennCNV.gcModel(proj, gmodelBase.getResource(), proj.GC_MODEL_FILENAME.getValue(), 100);
									refGenomeFasta = null;
								}
								if (Files.exists(refGenomeFasta) || Files.exists(proj.GC_MODEL_FILENAME.getValue())) {// TODO, after evaluating reference genome based gc model files, will demand a refGenome
									if (refGenomeFasta != null && Files.exists(refGenomeFasta)) {
										proj.REFERENCE_GENOME_FASTA_FILENAME.setValue(refGenomeFasta);
									}
									// try {
									GCAdjustorBuilder gAdjustorBuilder = new GCAdjustorBuilder();
									gAdjustorBuilder.regressionDistance(regressionDistance);

									params = GcAdjustorParameter.generate(proj, outputBase + "_GC_ADJUSTMENT/", refGenomeFasta, gAdjustorBuilder, sampsToUseRecompute, recomputeLRR_Median || recomputeLRR_PCs, bpGcModel, numThreads);
									if ((recomputeLRR_Median || recomputeLRR_PCs) && params.getCentroids() == null) {
										throw new IllegalStateException("Internal error, did not recieve centroids");
									} else if ((!recomputeLRR_Median && !recomputeLRR_PCs) && params.getCentroids() != null) {
										throw new IllegalStateException("Internal error, should not have recieved centroids");
									}
									recomputeLRR_Median = false;// recomputed if params has centroid
									recomputeLRR_PCs = false;
									// } catch (IllegalStateException e) {
									//
									// proj.getLog().reportTimeError("GC adjustment was flagged, but could not generate neccesary files");
									// }
								} else {
									proj.getLog().reportTimeError("Can not gc correct values without a valid reference genome");
									proj.getLog().reportTimeError("please supply a valid reference genome (full path) with the \"ref=\" argument");
								}
							}
							log.report("\nReady to perform the principal components analysis (PCA)\n");
							PrincipalComponentsCompute pcs = PCA.computePrincipalComponents(proj, false, numComponents, false, false, true, true, imputeMeanForNaN, recomputeLRR_PCs, proj.PROJECT_DIRECTORY.getValue() + outputBase + PCA_SAMPLES, outputBase, params);
							if (pcs == null) {
								return 3;
							}
							// apply PCs to everyone, we set useFile to null and excludeSamples to false to get all samples in the current project.
							// TODO, if we ever want to apply to only a subset of the project, we can do that here.....
							log.report("\nApplying the loadings from the principal components analysis to all samples\n");
							PrincipalComponentsApply pcApply = PCA.applyLoadings(proj, numComponents, pcs.getSingularValuesFile(), pcs.getMarkerLoadingFile(), null, false, imputeMeanForNaN, recomputeLRR_PCs, outputBase, params);
							// Compute Medians for (MT) markers and compute residuals from PCs for everyone
							log.report("\nComputing residuals after regressing out " + numComponents + " principal component" + (numComponents == 1 ? "" : "s") + "\n");
							PrincipalComponentsResiduals pcResids = PCA.computeResiduals(proj, pcApply.getExtrapolatedPCsFile(), ext.removeDirectoryInfo(medianMarkers), numComponents, true, 0f, homosygousOnly, recomputeLRR_Median, outputBase, params);
							generateFinalReport(proj, outputBase, pcResids.getResidOutput());
							proj.setProperty(proj.INTENSITY_PC_FILENAME, pcApply.getExtrapolatedPCsFile());
							proj.setProperty(proj.INTENSITY_PC_NUM_COMPONENTS, numComponents);
						}
					}
				}
			}
		}
		return 42;
	}

	/**
	 * if a subset of individuals is provided, we try to verify that they are all present in the project The String[] samples should be retrieved from the sample list so that it reflects parsed samples
	 */
	public static boolean verifyUseFile(Project proj, String[] samples, String useFile) {
		boolean allParsed = true;
		Logger log = proj.getLog();

		if (useFile != null) {
			String[] samplesToVerify = HashVec.loadFileToStringArray(useFile, false, new int[] { 0 }, false);
			Hashtable<String, String> track = new Hashtable<String, String>();
			ArrayList<String> notAvailable = new ArrayList<String>();
			ArrayList<String> available = new ArrayList<String>();
			for (int i = 0; i < samples.length; i++) {
				track.put(samples[i], samples[i]);
			}
			for (int i = 0; i < samplesToVerify.length; i++) {
				if (track.containsKey(samplesToVerify[i])) {
					available.add(samplesToVerify[i]);
				} else {
					notAvailable.add(samplesToVerify[i]);
					allParsed = false;
				}
			}
			if (notAvailable.size() > 0) {
				String missingFile = useFile + ".missing";
				String haveFile = useFile + ".have";
				Files.writeList(notAvailable.toArray(new String[notAvailable.size()]), missingFile);
				Files.writeList(available.toArray(new String[available.size()]), haveFile);
				log.reportError("Error - detected that not all samples (missing " + notAvailable.size() + ") from " + useFile + " are availble in the current project");
				log.reportError("	   - Please review the missing samples in " + missingFile + " and the samples available in " + haveFile + ". If you wish to continue after review,  change the argument \"" + USE_FILE_COMMAND + useFile + "\" to \"" + USE_FILE_COMMAND + haveFile + "\"");
				log.reportError("	   - Please note that sample names should correspond to the \"DNA\" column in the sample data file " + proj.PROJECT_DIRECTORY.getValue() + proj.getProperty(proj.SAMPLE_DATA_FILENAME) + ", and the name of the file (directory and extension removed) in " + proj.SAMPLE_DIRECTORY.getValue(false, false));
			}
		}
		return allParsed;
	}

	/**
	 * We try to verify that every sample has an entry in sample data, since we will need to look up ids from the PC file later.
	 * <p>
	 * The String[] samples should be retrieved from the sample list so that it reflects parsed samples
	 */
	private static boolean verifyAllSamples(Project proj, String[] samples) {
		boolean allParsed = true;
		SampleData sampleData = proj.getSampleData(0, false);
		ArrayList<String> notInSampleData = new ArrayList<String>();
		Logger log = proj.getLog();

		for (int i = 0; i < samples.length; i++) {
			if (sampleData.lookup(samples[i]) == null) {
				notInSampleData.add(samples[i]);
				allParsed = false;
			}
		}
		if (notInSampleData.size() > 0) {
			// log.reportError("Error - detected that some samples (missing " + notInSampleData.size() + ") do not have an entry in the sample data file " + proj.getFilename(proj.SAMPLE_DATA_FILENAME) + ", halting");
			log.reportError("Error - detected that some samples (missing " + notInSampleData.size() + ") do not have an entry in the sample data file " + proj.SAMPLE_DATA_FILENAME.getValue() + ", halting");
			log.reportError("	   - Please make sure the following samples have entries: " + Array.toStr(notInSampleData.toArray(new String[notInSampleData.size()]), "\n"));
		}
		return allParsed;
	}

	/**
	 * We try to ensure that all the markers contained in each aux file (PC markers/median markers) are available in the current project.
	 * <p>
	 * If they are not, we create lists of the markers available and missing TODO this could be consolidated with verify methods, but currently I wanted some specific reporting
	 */
	private static boolean verifyAuxMarkers(Project proj, String fileOfMarkers, String command) {
		boolean allAvailable = true;
		MarkerLookup markerLookup = proj.getMarkerLookup();
		String[] markersToVerify = HashVec.loadFileToStringArray(fileOfMarkers, false, new int[] { 0 }, false);
		ArrayList<String> notAvailable = new ArrayList<String>();
		ArrayList<String> available = new ArrayList<String>();
		Logger log = proj.getLog();

		if (markerLookup == null) {
			allAvailable = false;
		} else {
			for (int i = 0; i < markersToVerify.length; i++) {
				if (!markerLookup.contains(markersToVerify[i])) {
					notAvailable.add(markersToVerify[i]);
					allAvailable = false;
				} else {
					available.add(markersToVerify[i]);
				}
			}
		}
		if (notAvailable.size() > 0) {
			String missingFile = fileOfMarkers + ".missing";
			String haveFile = fileOfMarkers + ".have";
			Files.writeList(notAvailable.toArray(new String[notAvailable.size()]), missingFile);
			Files.writeList(available.toArray(new String[available.size()]), haveFile);
			log.reportError("Error - detected that not all markers (missing " + notAvailable.size() + ") from " + fileOfMarkers + " are availble in the current project");
			log.reportError("	   - Please review the missing markers in " + missingFile + " and the markers available in " + haveFile + ". If you wish to continue after review,  change the argument \"" + command + fileOfMarkers + "\" to \"" + command + haveFile + "\"");
		}
		return allAvailable;
	}

	/**
	 * We try to verify that every marker has a place in a transposed file
	 */
	private static boolean verifyAllProjectMarkersAreAvailable(Project proj) {
		boolean allParsed = true;
		ArrayList<String> notParsed = new ArrayList<String>();
		String[] markers = proj.getMarkerNames();
		MarkerLookup markerLookup = proj.getMarkerLookup();
		if (markerLookup == null) {
			allParsed = false;
		} else {
			for (int i = 0; i < markers.length; i++) {
				if (!markerLookup.contains(markers[i])) {
					notParsed.add(markers[i]);
					allParsed = false;
				}
			}
		}
		if (notParsed.size() > 0) {
			proj.getLog().reportError("Error - detected that not all markers (missing " + notParsed.size() + ") were properly parsed, halting: This should not happen");
		}
		return allParsed;
	}

	/**
	 * @param proj
	 * @param targetMarkersFile
	 * @param markersToQCFile
	 * @param markersABCallrate
	 * @param markerCallRateFilter
	 * @param numthreads
	 */
	public static void qcMarkers(Project proj, String targetMarkersFile, String markersToQCFile, String markersABCallrate, double markerCallRateFilter, int numthreads) {
		Logger log;
		String markerMetricsFilename;

		log = proj.getLog();
		markerMetricsFilename = proj.MARKER_METRICS_FILENAME.getValue(true, false);
		// skip if marker qc file exists
		if (Files.exists(markerMetricsFilename) && new File(markerMetricsFilename).length() > 0 && Files.exists(markersToQCFile) && Files.countLines(markerMetricsFilename, 1) >= Files.countLines(markersToQCFile, 0)) {
			log.report("Marker QC file " + proj.MARKER_METRICS_FILENAME.getValue(true, false) + " exists");
			log.report("Skipping Marker QC computation for the analysis, filtering on existing file");
		} else {
			log.report("Computing marker QC for " + targetMarkersFile == null ? "all markers in project." : "markers in " + targetMarkersFile);
			writeMarkersToQC(proj, targetMarkersFile, markersToQCFile);
			boolean[] samplesToExclude = new boolean[proj.getSamples().length];
			Arrays.fill(samplesToExclude, false);
			if (Thread.currentThread().isInterrupted()) {
				throw new RuntimeException(new InterruptedException());
			}
			MarkerMetrics.fullQC(proj, samplesToExclude, ext.removeDirectoryInfo(markersToQCFile), false, numthreads);

			// MarkerMetrics.fullQC(proj, samplesToExclude, null, false, numthreads);
		}
		if (Thread.currentThread().isInterrupted()) {
			throw new RuntimeException(new InterruptedException());
		}
		filterMarkerMetricsFile(proj, markerCallRateFilter, markersABCallrate);
	}

	/**
	 * Currently un-neccesary, but it is set up in case we want to QC the median markers at the same time
	 */
	private static void writeMarkersToQC(Project proj, String targetMarkersFile, String markersToQCFile) {
		String[] markers = null;
		if (targetMarkersFile == null) {
			markers = proj.getMarkerNames();
		} else {
			String[] markersToQC = { targetMarkersFile };
			markers = setMarkersToQC(proj, markersToQC);
		}
		// TODO remove CNVi probe markers (here?)
		Files.writeList(markers, markersToQCFile);
	}

	/**
	 * Similar to above, currently un-neccesary, but it is set up in case we want to QC the median markers at the same time
	 */
	private static String[] setMarkersToQC(Project proj, String[] files) {
		ArrayList<String> markersToUse = new ArrayList<String>();
		for (int i = 0; i < files.length; i++) {
			String[] markers = HashVec.loadFileToStringArray(files[i], false, new int[] { 0 }, false);
			for (int j = 0; j < markers.length; j++) {
				markersToUse.add(markers[j]);
			}
		}
		return PrincipalComponentsCompute.sortByProjectMarkers(proj, markersToUse.toArray(new String[markersToUse.size()]));
	}

	/**
	 * Write a filtered list of markers to use for ab callRate in sample QC
	 */
	private static boolean filterMarkerMetricsFile(Project proj, double markerCallRateFilter, String markersABCallrate) {
		ArrayList<String> abMarkersToUse = new ArrayList<String>();
		BufferedReader reader;
		Logger log = proj.getLog();

		try {
			// reader = Files.getReader(proj.getFilename(proj.MARKER_METRICS_FILENAME), false, true, false);
			reader = Files.getReader(proj.MARKER_METRICS_FILENAME.getValue(), false, true, false);
			String[] header = reader.readLine().trim().split("\t");
			int abIndex = ext.indexOfStr(MarkerMetrics.FULL_QC_HEADER[2], header);
			if (abIndex == -1) {
				// log.reportError("Error - the necessary marker metrics header " + MarkerMetrics.FULL_QC_HEADER[2] + " was not found in the marker metrics file" + proj.getFilename(proj.MARKER_METRICS_FILENAME));
				log.reportError("Error - the necessary marker metrics header " + MarkerMetrics.FULL_QC_HEADER[2] + " was not found in the marker metrics file" + proj.MARKER_METRICS_FILENAME.getValue());
				return false;
			} else {
				String[] metrics;
				while (reader.ready()) {
					// proj.getArrayType().isCNOnly(markerName)
					metrics = reader.readLine().trim().split("\t");
					try {
						double callRate = Double.parseDouble(metrics[abIndex]);
						if (callRate >= markerCallRateFilter) {
							abMarkersToUse.add(metrics[0]);
						}
					} catch (NumberFormatException nfe) {
						log.report("Warning - found an invalid number " + metrics[abIndex] + " for marker" + metrics[0] + " skipping this marker");
					}
				}
			}
			if (abMarkersToUse.size() == 0) {
				log.reportError("Error - no markers passed the callRate threshold. Please consider lowering threshold, or ensure that markers can have call rates (not cnv only probes)");
				return false;
			} else {
				log.report("Sample call rate will be computed with " + abMarkersToUse.size() + " markers");
				Files.writeList(abMarkersToUse.toArray(new String[abMarkersToUse.size()]), markersABCallrate);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			// log.reportError("Error: file \"" + proj.getFilename(proj.MARKER_METRICS_FILENAME) + "\" not found in current directory");
			log.reportError("Error: file \"" + proj.MARKER_METRICS_FILENAME.getValue() + "\" not found in current directory");
		} catch (IOException ioe) {
			// log.reportError("Error reading file \"" + proj.getFilename(proj.MARKER_METRICS_FILENAME) + "\"");
			log.reportError("Error reading file \"" + proj.MARKER_METRICS_FILENAME.getValue() + "\"");
		}

		return true;
	}

	/**
	 * The last thing to do, gives us a report of summarized median markers, qc stats, and pcs..
	 */
	public static void generateFinalReport(Project proj, String outputBase, String residualFile) {
		BufferedReader reader;
		PrintWriter writer;
		int DNAIndex = 0;
		Hashtable<String, String> qcLookup = loadQC(proj, outputBase);
		Logger log = proj.getLog();

		try {
			String finalReport = ext.rootOf(residualFile) + PCA_FINAL_REPORT;
			if (Files.exists(proj.PROJECT_DIRECTORY.getValue() + finalReport)) {
				Files.backup(finalReport, proj.PROJECT_DIRECTORY.getValue(), proj.PROJECT_DIRECTORY.getValue() + proj.getProperty(proj.BACKUP_DIRECTORY));
			}

			DNAIndex = getDNAIndex(proj, proj.PROJECT_DIRECTORY.getValue() + residualFile);
			reader = Files.getReader(proj.PROJECT_DIRECTORY.getValue() + residualFile, false, true, false);
			writer = Files.getAppropriateWriter(proj.PROJECT_DIRECTORY.getValue() + finalReport);

			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("\t");
				String key = line[DNAIndex];
				if (qcLookup.containsKey(key)) {
					writer.println(Array.toStr(line) + "\t" + qcLookup.get(key));
				} else {
					log.reportError("Error - could not match ids " + key + " in the qc file to produce final report");
				}

			}
			reader.close();
			writer.close();
			log.report(ext.getTime() + " Finished summarizing to final report at " + finalReport);

		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + residualFile + "\" not found in current directory");
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + residualFile + "\"");
		}
	}

	/**
	 * Loads the PCA_SAMPLES_SUMMARY File
	 */
	private static Hashtable<String, String> loadQC(Project proj, String outputBase) {
		BufferedReader reader;
		Hashtable<String, String> qcLookup = new Hashtable<String, String>();
		int DNAIndex = 0;
		Logger log = proj.getLog();
		String sampleQcFile = proj.PROJECT_DIRECTORY.getValue() + outputBase + PCA_SAMPLES_SUMMARY;
		try {
			DNAIndex = getDNAIndex(proj, sampleQcFile);
			reader = Files.getReader(sampleQcFile, false, true, false);
			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("\t");
				qcLookup.put(line[DNAIndex], Array.toStr(Array.subArray(line, 1)));
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + sampleQcFile + "\" not found.");
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + sampleQcFile + "\"");
		}
		return qcLookup;
	}

	/**
	 * Tries to get the DNA index so we always match up
	 */
	private static int getDNAIndex(Project proj, String fileName) {
		Logger log = proj.getLog();
		int DNAIndex;
		DNAIndex = ext.indexOfEndsWith(DNA_LINKER, Files.getHeaderOfFile(fileName, log), false);
		if (DNAIndex == -1) {
			log.report("Warning - could not find the linker " + DNA_LINKER + " in file " + fileName + ", assuming it is the first index");
			DNAIndex = 0;
		}
		return DNAIndex;
	}

	/**
	 * This function will attempt to generate an AB lookup for the project. It should only be called if it is required
	 */
	static int generateABLookup(Project proj, Logger log) {
		ABLookup abLookup;
		String snpMapFile;
		abLookup = new ABLookup();
		abLookup.parseFromGenotypeClusterCenters(proj);
		abLookup.writeToFile(proj.AB_LOOKUP_FILENAME.getValue(false, false), proj.getLog());
		if (Files.exists(proj.AB_LOOKUP_FILENAME.getValue(false, false))) {
			snpMapFile = proj.getLocationOfSNP_Map(true);
			if (snpMapFile != null) {
				log.report("Info - attempting to fill in missing alleles from " + snpMapFile);
				ABLookup.fillInMissingAlleles(proj, proj.AB_LOOKUP_FILENAME.getValue(false, false), snpMapFile, false);
				if (Files.exists(ext.addToRoot(proj.AB_LOOKUP_FILENAME.getValue(false, false), "_filledIn"))) {
					proj.setProperty(proj.AB_LOOKUP_FILENAME, ext.addToRoot(proj.AB_LOOKUP_FILENAME.getValue(false, false), "_filledIn"));
					proj.saveProperties();
				} else {
					log.reportError("Error - detected " + snpMapFile + ", but could not fill in missing AB Lookup values. Reverting to " + proj.AB_LOOKUP_FILENAME.getValue(false, false));
				}
			} else {
				log.report("Warning - will not be able to fill in missing alleles from " + proj.AB_LOOKUP_FILENAME.getValue(false, false) + ", a SNP_MAP file could not be found");
			}
			ABLookup.applyABLookupToFullSampleFiles(proj);
			return 1;
		} else {
			log.reportError("Error - detected that an AB Lookup file is required, but failed to create AB Lookup file " + proj.AB_LOOKUP_FILENAME.getValue(false, false));
			return 0;
		}
	}

	private static void mkdir(String outputDirectory, Logger log) {
		File file = new File(outputDirectory);
		if (!file.exists()) {
			if (file.mkdirs()) {
				log.report("\n" + ext.getTime() + " Created directory " + outputDirectory);
			} else {
				log.reportError("Error - failed to create  " + outputDirectory + ", please manually create it unless it already exists");
			}
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String projectName = "GENVISIS_project1";

		boolean[] requiredArray = new boolean[5];
		Arrays.fill(requiredArray, false);
		String[] requiredArgs = { "dirProj=", "dirSrc=", PC_MARKER_COMMAND, MITO_MARKER_COMMAND, "markerPositions=" };

		String projectDirectory = null;
		String sourceDirectory = null;
		String targetMarkers = null;
		String markerPositions = null;
		String medianMarkers = null;

		String output = FILE_BASE;
		String idHeader = "Sample ID";
		String abLookup = null;
		// String pedFile = "D:/data/TestAuto/testped.txt";
		String pedFile = null;
		// String sampleMapCsv = null;
		String sampleMapCsv = null;
		String dataExtension = ".csv";
		String sampleLRRSdFilter = "0.5";
		String sampleCallRateFilter = "0.95";
		double markerCallRateFilter = 0.98;
		boolean markerQC = true;
		boolean recomputeLRR_PCs = true;
		boolean recomputeLRR_Median = true;
		String useFile = null;
		String logfile = null;
		Project proj;

		int numThreads = 1;
		int numComponents = 100;
		boolean imputeMeanForNaN = true;
		boolean homosygousOnly = true;
		boolean doAbLookup = false;

		String referenceGenomeFasta = null;
		String gcmodel = null;
		int regressionDistance = GcAdjustor.DEFAULT_REGRESSION_DISTANCE[0];
		int bpGcModel = GcAdjustor.GcModel.DEFAULT_GC_MODEL_BIN_FASTA;
		boolean gcCorrect = true;
		int attempts, result;
		GENOME_BUILD build = GENOME_BUILD.HG19;
		boolean recompSampleSpecific = true;

		String usage = "\n";
		usage += "The MitoPipeline currently requires 5 arguments and allows for many more optional arguments:\n";
		usage += "  \n";
		usage += "   (1) The full path for the project directory (where results will be stored) (i.e. dirProj=/home/usr/projects/)\n";
		usage += "   (2) The full path for the source data directory  (where final report files are located) (i.e. dirSrc=/home/usr/data/project1/)\n";
		usage += "   (3) The full path for a file with a list of markers (one per line) to use for computing PCs (i.e. " + PC_MARKER_COMMAND + "/home/usr/auxFiles/exomeChip.PC_Markers.txt)\n";
		usage += "   (4) The full path for a file with a list of markers (one per line,mitochondrial markers) to use for computing computing median Log R Ratios (i.e. " + MITO_MARKER_COMMAND + "/home/usr/auxFiles/exomeChip.MT_Markers.txt)\n";
		usage += "   (5) The full path for a tab-delimited file with marker positions (with columns \"Marker\", \"Chr\", and \"Position\")  (i.e. markerPositions=/home/usr/auxFiles/exomeChip.Positions.txt)\n";

		usage += "   OPTIONAL:\n";
		usage += "	 (6) A file listing a subset of samples (DNA ID) to use for QC and PC computation portions of the analysis, often a list of unrelated individuals. If a list is not provided, all samples in the source directory will be analyzed (i.e. " + USE_FILE_COMMAND + useFile + " (no default))\n";
		usage += "   (7) The full path for a tab-delimited .PED format file with header \"" + Array.toStr(PED_INPUT) + "\" (i.e. pedFile=" + pedFile + "(no default))\n";
		usage += "   OR:\n";
		usage += "   (8) The full path for a Sample_Map.csv file, with at least two columns having headers \"" + SAMPLEMAP_INPUT[1] + "\" and \"" + SAMPLEMAP_INPUT[2] + "\"(i.e. mapFile=" + sampleMapCsv + " (default))\n\n";
		usage += "   NOTE:\n";
		usage += "   All samples to be analyzed must be contained in the sample manifest (.PED format file, or Sample_Map.csv file)\n";
		usage += "   (9) The desired name of the project (i.e. projName=" + projectName + " (default))\n";
		usage += "   (10) Data extension for files contained in the source data directory (i.e. dirExt=" + dataExtension + " (default))\n";
		usage += "   (11) Log R Ratio standard deviation filter to exclude samples from PCs (i.e. LRRSD=" + sampleLRRSdFilter + " (default))\n";
		usage += "   (12) Call rate filter to exclude samples from PCs (i.e. sampleCallRate=" + sampleCallRateFilter + " (default))\n";
		usage += "   (13) Number of principal components to compute (must be less than the number of samples AND the number of markers) (i.e. numComponents=" + numComponents + " (default))\n";
		usage += "   (14) Number of threads to use for multi-threaded portions of the analysis (i.e. numThreads=" + numThreads + " (default))\n";
		usage += "   (15) Output file full path and baseName (i.e. output=" + output + " (default))\n";
		usage += "   (16) Project filename (if you manually created a project properties file, or edited an existing project). Note that default arguments available here can overide existing project properties (i.e. proj=" + filename + " (no default))\n";
		usage += "   (17) The header of the column containing sample ids in the final report files (for command-line interpretability, space characters must be replaced with \"_\". Common options are \"Sample_ID\" and \"Sample_Name\", corresponding to \"Sample ID\" and \"Sample Name\")  (i.e. idHeader=" + idHeader + " (default))\n";
		// usage += "   (18) A file specifying the AB allele lookup for markers, often times required  (i.e. abLookup=" + idHeader + " (default))\n";
		usage += "   (18) Do not perform a marker qc step to select higher quality markers (or remove cnv-only markers) to use for computing the sample call rate (i.e. -nomarkerQC (not the default))\n";
		usage += "   (19) If marker qc is performed, the call rate cutoff for markers to be passed on to the sample QC step (i.e. markerCallRate=" + markerCallRateFilter + " (default))\n";
		usage += "   (20) Name of the log file (i.e. log=[project_directory]/logs/Genvisis_[date].log (default))\n";
		usage += "   (21) Recompute Log R Ratios for each marker from genotypes/intensities when computing AND extrapolating PCs(i.e. recomputeLRR_PCs=" + recomputeLRR_PCs + " (default))\n";
		usage += "   (22) Recompute Log R Ratios for each marker from genotypes/intensities when computing median values(i.e. recomputeLRR_Median=" + recomputeLRR_Median + " (default))\n";
		usage += "   (23) Impute mean for markers with NaN data, if false markers with NaN values for any sample will be skipped (i.e. imputeMeanForNaN=" + imputeMeanForNaN + " (default))\n";
		usage += "   (24) gc correct Log R Ratios, cannot be used with recomputeLRR options (i.e. gcCorrect=" + gcCorrect + " (default))\n";
		usage += "   (25) A reference genome file used to assign gc content to each marker (i.e. ref= (no default))\n";
		usage += "   (26) base-pair bins for the gc model generated from the reference (i.e. bpGcModel=" + bpGcModel + " (default))\n";
		usage += "   (27) regression distance for the gc adjustment (i.e. regressionDistance=" + regressionDistance + " (default))\n";
		usage += "   (28) full path to a .gcmodel file, this model will take precedence over base-pair bins, and the reference genome will not be used (i.e. gcmodel=" + gcmodel + " (default))\n";
		usage += "   (29) if relying on a gc5Base.txt (default), the genomic build to use (i.e. build=" + build + " (default))\n";
		usage += "   (30) recompute LRR using only those samples that pass QC, and are in the use file (i.e. sampLRR=" + recompSampleSpecific + " (default))\n";

		usage += "   NOTE:\n";
		usage += "   Project properties can be manually edited in the .properties file for the project. If you would like to use an existing project properties file, please specify the filename using the \"proj=\" argument\n";
		usage += "   Editing the project properties file can be useful when a command line option is not available\n";

		usage += "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("projName=")) {
				projectName = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("sampLRR=")) {
				recompSampleSpecific = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("bpGcModel=")) {
				bpGcModel = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("regressionDistance=")) {
				regressionDistance = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("proj=")) {
				filename = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("dirProj=")) {
				projectDirectory = ext.parseStringArg(args[i], null);
				numArgs--;
				requiredArray[0] = true;
			} else if (args[i].startsWith("dirSrc=")) {
				sourceDirectory = ext.parseStringArg(args[i], null);
				numArgs--;
				requiredArray[1] = true;
			} else if (args[i].startsWith("pedFile=")) {
				pedFile = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("ref=")) {
				referenceGenomeFasta = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("mapFile=")) {
				sampleMapCsv = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith(IMPORT_EXTENSION)) {
				dataExtension = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("LRRSD=")) {
				sampleLRRSdFilter = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("sampleCallRate=")) {
				sampleCallRateFilter = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("gcmodel=")) {
				gcmodel = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith(PC_MARKER_COMMAND)) {
				targetMarkers = ext.parseStringArg(args[i], null);
				numArgs--;
				requiredArray[2] = true;
			} else if (args[i].startsWith("output=")) {
				output = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith(MITO_MARKER_COMMAND) || args[i].startsWith("medianMarkers=")) {
				medianMarkers = ext.parseStringArg(args[i], null);
				numArgs--;
				requiredArray[3] = true;
			} else if (args[i].startsWith("markerPositions=")) {
				markerPositions = ext.parseStringArg(args[i], null);
				numArgs--;
				requiredArray[4] = true;
			} else if (args[i].startsWith("numThreads=")) {
				numThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("numComponents=")) {
				numComponents = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("-allCalls")) {
				homosygousOnly = false;
				numArgs--;
			} else if (args[i].startsWith("-nomarkerQC")) {
				markerQC = true;
				numArgs--;
			} else if (args[i].startsWith("recomputeLRR_PCs=")) {
				recomputeLRR_PCs = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("recomputeLRR_Median=")) {
				recomputeLRR_Median = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("gcCorrect=")) {
				gcCorrect = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("imputeMeanForNaN=")) {
				imputeMeanForNaN = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("markerCallRate=")) {
				markerCallRateFilter = ext.parseDoubleArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith(USE_FILE_COMMAND)) {
				useFile = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("abLookup=")) {
				abLookup = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("idHeader=")) {
				System.out.println(args[i]);
				idHeader = ext.parseStringArg(args[i].replaceAll("_", " "), null);
				numArgs--;
			} else if (args[i].startsWith("build=")) {
				try {
					build = GENOME_BUILD.valueOf(ext.parseStringArg(args[i], ""));
					numArgs--;
				} catch (IllegalArgumentException ile) {
					System.err.println("Invalid build " + ext.parseStringArg(args[i], ""));
					System.err.println("Options Are: ");
					for (int j = 0; j < GENOME_BUILD.values().length; j++) {
						System.err.println(GENOME_BUILD.values()[j]);
					}
				}
			} else if (args[i].startsWith("log=")) {
				logfile = ext.parseStringArg(args[i], null);
				numArgs--;
			} else {
				System.err.println("\n\nError - invalid argument " + args[i] + "\n\n");
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

		if (Array.booleanArraySum(requiredArray) != requiredArray.length) {
			System.err.println(usage + "\n\n");
			System.err.println("The MitoPipeline currently requires " + requiredArray.length + " arguments and we only detected " + Array.booleanArraySum(requiredArray) + " of the " + requiredArray.length);
			System.err.println("Here is a list of missing arguments...");
			for (int i = 0; i < requiredArgs.length; i++) {
				if (!requiredArray[i]) {
					System.err.println("\"" + requiredArgs[i] + "\"");
				}
			}
			System.exit(1);
		}

		initGenvisisProject();
		proj = null;
		if (filename == null) {
			proj = initializeProject(proj, projectName, projectDirectory, sourceDirectory, dataExtension, idHeader, abLookup, targetMarkers, medianMarkers, markerPositions, sampleLRRSdFilter, sampleCallRateFilter, logfile);
		} else {
			proj = new Project(filename, logfile, false);
			proj = initializeProject(proj, projectName, projectDirectory, sourceDirectory, dataExtension, idHeader, abLookup, targetMarkers, medianMarkers, markerPositions, sampleLRRSdFilter, sampleCallRateFilter, logfile);
		}
		attempts = 0;
		while (attempts < 2) {
			if (gcmodel != null) {
				if (!proj.GC_MODEL_FILENAME.getValue().equals(gcmodel)) {
					Files.copyFileUsingFileChannels(gcmodel, proj.GC_MODEL_FILENAME.getValue(), proj.getLog());
				}
				// proj.GC_MODEL_FILENAME.setValue(gcmodel);
				if (referenceGenomeFasta != null) {
					proj.getLog().reportTimeWarning("Ignoring reference genome " + referenceGenomeFasta);
					proj.REFERENCE_GENOME_FASTA_FILENAME.setValue("Ignore");
				}
			}
			result = catAndCaboodle(proj, numThreads, medianMarkers, numComponents, output, homosygousOnly, markerQC, markerCallRateFilter, useFile, pedFile, sampleMapCsv, recomputeLRR_PCs, recomputeLRR_Median, recompSampleSpecific, doAbLookup, imputeMeanForNaN, gcCorrect, referenceGenomeFasta, bpGcModel, regressionDistance, build);
			attempts++;
			if (result == 41 || result == 40) {
				proj.getLog().report("Attempting to restart pipeline once to fix SampleList problem");
				if (result == 40) {// ParseIllumina determined an AB Lookup was necessary if result is 40
					doAbLookup = true;
				}
			} else {
				attempts++;
			}
		}
	}
}
