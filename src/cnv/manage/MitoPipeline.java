package cnv.manage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;

import cnv.LaunchProperties;
import cnv.analysis.pca.PCA;
import cnv.analysis.pca.PrincipalComponentsApply;
import cnv.analysis.pca.PrincipalComponentsCompute;
import cnv.analysis.pca.PrincipalComponentsResiduals;
import cnv.filesys.ABLookup;
import cnv.filesys.MarkerLookup;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.filesys.SampleList;
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
 * */
public class MitoPipeline {
	public static final String[] PED_INPUT = { "DNA", "FID", "IID", "FA", "MO", "SEX", "AFF" };
	public static final String[] SAMPLEMAP_INPUT = { "Index", "Name", "ID", "Gender", "Plate", "Well", "Group", "Parent1", "Parent2", "Replicate", "SentrixPosition" };
	public static final String[] QC_COLUMNS = { "Sample", "LRR_SD", "Genotype_callrate" };
	public static final String[] SAMPLE_QC_SUMMARY = { "DNA", "LRR_SD", "Genotype_callrate", "Included in PC?" };
	public static final String[] SEX = { "female", "male" };
	public static final String[] SAMPLE_DATA_ADDITION_HEADERS = { "LRR_SD", "Genotype_callrate", "CLASS=Exclude" };

	private static final String DNA_LINKER = "DNA";
	private static final String MARKERS_TO_QC_FILE = "markers_to_QC.txt";
	static final String MARKERS_FOR_ABCALLRATE = "markers_ABCallRate.txt";
	protected static final String PCA_SAMPLES = ".samples.USED_PC.txt";
	private static final String PCA_SAMPLES_SUMMARY = ".samples.QC_Summary.txt";
	private static final String PCA_FINAL_REPORT = ".finalReport.txt";
	private static final String[] SPLITS = { "\t", "," };
	public static final String PROJECT_EXT = ".properties";
	private static final String DEFAULT_QC_FILE = "lrr_sd.xln";
	private static final String PC_MARKER_COMMAND = "PCmarkers=";
	private static final String MEDIAN_MARKER_COMMAND = "medianMarkers=";
	private static final String USE_FILE_COMMAND = "useFile=";
	private static final String IMPORT_EXTENSION = "dirExt=";
	private static final long RECOMMENDED_MEMORY = 1000000000;

	private String projectName;
	private String filename;
	private Project proj;
	private Logger log;
	private String projectDirectory, sourceDirectory, dataExtension, defaultLRRSdFilter, targetMarkers, medianMarkers, markerPositions, idHeader, abLookup;

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
		this.targetMarkers = targetMarkers;
		this.markerPositions = markerPositions;
		this.medianMarkers = medianMarkers;
		this.log = new Logger(logfile);
		initProjectDir();
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
//		if (!new File(path + "projects/").exists()) {
//			new File(path + "projects/").mkdirs();
//		}
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
//		if (!new File(launchPropertiesFile).exists()) {
//		    new File(path + "example/").mkdirs();
//		    Files.writeList(new String[] { "LAST_PROJECT_OPENED=example.properties", "PROJECTS_DIR=" + path + "projects/" }, launchPropertiesFile);
//		    if (!new File(path + "projects/example.properties").exists()) {
//		        Files.writeList(new String[] { "PROJECT_NAME=Example", "PROJECT_DIRECTORY=example/", "SOURCE_DIRECTORY=sourceFiles/" }, path + "projects/example.properties");
//		    }
//		}
//		return path + "projects/";
	}

	/**
	 * Copies the default project to the project directory if the desired fileName does not already exist
	 */
	public void initProject(String path) {
		if (Files.exists(filename)) {
			Files.backup(ext.removeDirectoryInfo(filename), path, path + "backup/", false);
			log.report("Using project file " + filename + ", you may also specify project filename using the command line argument \"proj=\"");
		} else {
			log.report("Project properties file can be found at " + filename);
			if (proj != null) {
				Files.write(proj.PROJECT_NAME.getName() + "=" + projectName, filename);
			}
		}
		this.proj = new Project(filename, false);

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
//		proj.setProperty(proj.LRRSD_CUTOFF, defaultLRRSdFilter);
		proj.setProperty(proj.LRRSD_CUTOFF, Double.valueOf(defaultLRRSdFilter));
		proj.setProperty(proj.TARGET_MARKERS_FILENAME, ext.removeDirectoryInfo(targetMarkers));
		if (markerPositions != null) {
			proj.setProperty(proj.MARKER_POSITION_FILENAME, ext.removeDirectoryInfo(markerPositions));
		}
		if (abLookup != null && Files.exists(projectDirectory + abLookup)) {
			proj.setProperty(proj.AB_LOOKUP_FILENAME, ext.removeDirectoryInfo(abLookup));
		}
		proj.saveProperties();
	}

	public static void generateSampleDataMap(Project proj, String sampleMapCsv) {
		generateSampleData(proj, loadSampleMapFile(sampleMapCsv, proj.getLog()));
	}

	public static void generateSampleDataPed(Project proj, String pedFile) {
		generateSampleData(proj, loadPedFile(pedFile, proj.getLog()));
	}

	public static Project initializeProject(Project proj, String projectName, String projectDirectory, String sourceDirectory, String dataExtension, String idHeader, String abLookup, String targetMarkers, String medianMarkers, String markerPositions, String defaultLRRSdFilter, String defaultCallRateFilter, String logfile) {
		MitoPipeline projG = new MitoPipeline(proj, projectName, projectDirectory, sourceDirectory, dataExtension, idHeader, abLookup, defaultLRRSdFilter, defaultCallRateFilter, targetMarkers, markerPositions, medianMarkers, logfile);
		return projG.getProj();
	}

	/**
	 * @param pedFile
	 *            ped format file to create sample data,
	 * @param sampleMapCsv
	 *            sample_Map.csv format file to create sample data
	 * @param proj
	 *            an existing, or newly created project
	 * @param log
	 *            Note: if the pedFile and sampleMapCsv file are both null, we create a minimal sample data instead Note: if sample data already exists, we leave it alone
	 */
	public static int createSampleData(String pedFile, String sampleMapCsv, Project proj) {
		String sampleDataFilename = proj.SAMPLE_DATA_FILENAME.getValue(false, false);
		if ((sampleMapCsv == null || "".equals(sampleMapCsv) || !Files.exists(sampleMapCsv)) && (pedFile == null || "".equals(pedFile) || !Files.exists(pedFile)) && !Files.exists(sampleDataFilename)) {
			proj.getLog().report("Neither a sample manifest nor a sample map file was provided; generating sample data file at: " + sampleDataFilename);
			SampleData.createMinimalSampleData(proj);
			return 0;
		} else if (!Files.exists(sampleDataFilename)) {
			if (pedFile != null) {
				generateSampleDataPed(proj, pedFile);
				return 1;
			} else {
				generateSampleDataMap(proj, sampleMapCsv);
				return 2;
			}
		} else {
			proj.getLog().report("Detected that a sampleData file already exists at " + sampleDataFilename + ", skipping sampleData creation");
			return -1;
		}
	}

	/**
	 * The main event. Takes the samples from raw data through import and PCA
	 */

	public static int catAndCaboodle(Project proj, int numThreads, String sampleCallRateFilter, String medianMarkers, int numComponents, String outputBase, boolean homosygousOnly, boolean markerQC, double markerCallRateFilter, String useFile, String pedFile, String sampleMapCsv, boolean recomputeLRR_PCs, boolean recomputeLRR_Median, boolean doAbLookup, boolean imputeMeanForNaN) {
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
//			log.report("Skipping sample import step for the analysis. If this is an incorrect number of samples, please remove (or change the name of) " + proj.getFilename(proj.SAMPLELIST_FILENAME) + " and " + proj.getDir(proj.SAMPLE_DIRECTORY));
			log.report("Skipping sample import step for the analysis. If this is an incorrect number of samples, please remove (or change the name of) " + proj.SAMPLELIST_FILENAME.getValue() + " and " + proj.SAMPLE_DIRECTORY.getValue(false, true));
		} else {
			result = cnv.manage.ParseIllumina.createFiles(proj, numThreads);
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

		createSampleData(pedFile, sampleMapCsv, proj);
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
//					if (verifyAuxMarkers(proj, proj.getFilename(proj.TARGET_MARKERS_FILENAME), PC_MARKER_COMMAND)) {
					if (verifyAuxMarkers(proj, proj.TARGET_MARKERS_FILENAME.getValue(), PC_MARKER_COMMAND)) {
						// if marker QC is not flagged, sample qc is based on all target markers by default
						if (markerQC) {
							qcMarkers(proj, markerCallRateFilter, numThreads);
							markersForABCallRate = proj.PROJECT_DIRECTORY.getValue() + MARKERS_FOR_ABCALLRATE;
							if (!Files.exists(markersForABCallRate)) {
								log.reportError("Error - markerQC was flagged but the file " + proj.PROJECT_DIRECTORY.getValue() + MARKERS_FOR_ABCALLRATE + " could not be found");
								return 1;
							}
						} else {
							markersForABCallRate = proj.PROJECT_DIRECTORY.getValue() + MARKERS_TO_QC_FILE;
							writeMarkersToQC(proj);
						}

						markersForEverythingElse = proj.PROJECT_DIRECTORY.getValue() + MARKERS_TO_QC_FILE;

						counts = filterSamples(proj, outputBase, markersForABCallRate, markersForEverythingElse, numThreads, sampleCallRateFilter, useFile);
						if (counts == null || counts[1] != sampleList.getSamples().length) {
							if (counts == null || counts[1] == 0 && Files.exists(proj.PROJECT_DIRECTORY.getValue() + DEFAULT_QC_FILE)) {
								log.reportError("Error - was unable to parse QC file " + proj.PROJECT_DIRECTORY.getValue() + DEFAULT_QC_FILE + ", backing up this file to " + proj.BACKUP_DIRECTORY.getValue(false, false) + " and re-starting sample qc");
								Files.backup(DEFAULT_QC_FILE, proj.PROJECT_DIRECTORY.getValue(), proj.BACKUP_DIRECTORY.getValue(true, false), true);
							}
							counts = filterSamples(proj, outputBase, markersForABCallRate, markersForEverythingElse, numThreads, sampleCallRateFilter, useFile);
							if (counts == null || counts[1] != sampleList.getSamples().length) {
								if (counts == null) {
									log.reportError("Error - could not parse QC file (" + DEFAULT_QC_FILE + ")");
								} else {
									log.reportError("Error - different number of samples (n=" + counts[1] + ") listed in the QC file (" + DEFAULT_QC_FILE + ") compared to the number of samples in the project (n=" + sampleList.getSamples().length + ")");
									log.reportError("      - delete the QC file (" + DEFAULT_QC_FILE + ") to regenerate it with the correct number of samples");
								}
								log.reportError("aborting...");
								return 2;
							}
						}
						if (counts == null || counts[0] == 0) {// no samples passed threshold, null case shouldn't happen but we will test anyway
							return 2;// message handled already
						}
						// check that all median markers are available
						if (verifyAuxMarkers(proj, medianMarkers, MEDIAN_MARKER_COMMAND)) {
							// compute PCs with samples passing QC
							log.report("\nReady to perform the principal components analysis (PCA)\n");
							PrincipalComponentsCompute pcs = PCA.computePrincipalComponents(proj, false, numComponents, false, false, true, true, imputeMeanForNaN, recomputeLRR_PCs, outputBase + PCA_SAMPLES, outputBase);
							if (pcs == null) {
								return 3;
							}
							// apply PCs to everyone, we set useFile to null and excludeSamples to false to get all samples in the current project.
							// TODO, if we ever want to apply to only a subset of the project, we can do that here.....
							log.report("\nApplying the loadings from the principal components analysis to all samples\n");
							PrincipalComponentsApply pcApply = PCA.applyLoadings(proj, numComponents, pcs.getSingularValuesFile(), pcs.getMarkerLoadingFile(), null, false, imputeMeanForNaN, recomputeLRR_PCs, outputBase);
							// Compute Medians for (MT) markers and compute residuals from PCs for everyone
							log.report("\nComputing residuals after regressing out " + numComponents + " principal component" + (numComponents == 1 ? "" : "s") + "\n");
							PrincipalComponentsResiduals pcResids = PCA.computeResiduals(proj, pcApply.getExtrapolatedPCsFile(), ext.removeDirectoryInfo(medianMarkers), numComponents, true, 0f, homosygousOnly, recomputeLRR_Median, outputBase);
							generateFinalReport(proj, outputBase, pcResids.getResidOutput());
							proj.setProperty(proj.INTENSITY_PC_FILENAME, pcApply.getExtrapolatedPCsFile());
//							proj.setProperty(proj.INTENSITY_PC_NUM_COMPONENTS, numComponents + "");
							proj.setProperty(proj.INTENSITY_PC_NUM_COMPONENTS, numComponents);
							proj.saveProperties();
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
//			log.reportError("Error - detected that some samples (missing " + notInSampleData.size() + ") do not have an entry in the sample data file " + proj.getFilename(proj.SAMPLE_DATA_FILENAME) + ", halting");
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

	public static void qcMarkers(Project proj, double markerCallRateFilter, int numthreads) {
		Logger log;
		String markerMetricsFilename;

		log = proj.getLog();
		markerMetricsFilename = proj.MARKER_METRICS_FILENAME.getValue(true, false);
		// skip if marker qc file exists
		if (Files.exists(markerMetricsFilename) && new File(markerMetricsFilename).length() > 0 && Files.exists(proj.PROJECT_DIRECTORY.getValue() + MARKERS_TO_QC_FILE) && Files.countLines(markerMetricsFilename, true) >= Files.countLines(proj.PROJECT_DIRECTORY.getValue() + MARKERS_TO_QC_FILE, false)) {
			log.report("Marker QC file " + proj.MARKER_METRICS_FILENAME.getValue(true, false) + " exists");
			log.report("Skipping Marker QC computation for the analysis, filtering on existing file");
		} else {
//			log.report("Computing marker QC for markers in " + proj.getFilename(proj.TARGET_MARKERS_FILENAME));
			log.report("Computing marker QC for markers in " + proj.TARGET_MARKERS_FILENAME.getValue());
			writeMarkersToQC(proj);
			boolean[] samplesToExclude = new boolean[proj.getSamples().length];
			Arrays.fill(samplesToExclude, false);
			MarkerMetrics.fullQC(proj, samplesToExclude, MARKERS_TO_QC_FILE, numthreads);
		}
		filterMarkerMetricsFile(proj, markerCallRateFilter);
	}

	/**
	 * Currently un-neccesary, but it is set up in case we want to QC the median markers at the same time
	 */
	private static void writeMarkersToQC(Project proj) {
//		String[] markersToQC = { proj.getFilename(proj.TARGET_MARKERS_FILENAME, true, false) };
		String[] markersToQC = { proj.TARGET_MARKERS_FILENAME.getValue(true, false) };
		Files.writeList(setMarkersToQC(proj, markersToQC), proj.PROJECT_DIRECTORY.getValue() + MARKERS_TO_QC_FILE);
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
	private static boolean filterMarkerMetricsFile(Project proj, double markerCallRateFilter) {
		ArrayList<String> abMarkersToUse = new ArrayList<String>();
		BufferedReader reader;
		Logger log = proj.getLog();

		try {
//			reader = Files.getReader(proj.getFilename(proj.MARKER_METRICS_FILENAME), false, true, false);
			reader = Files.getReader(proj.MARKER_METRICS_FILENAME.getValue(), false, true, false);
			String[] header = reader.readLine().trim().split("\t");
			int abIndex = ext.indexOfStr(MarkerMetrics.FULL_QC_HEADER[2], header);
			if (abIndex == -1) {
//				log.reportError("Error - the necessary marker metrics header " + MarkerMetrics.FULL_QC_HEADER[2] + " was not found in the marker metrics file" + proj.getFilename(proj.MARKER_METRICS_FILENAME));
				log.reportError("Error - the necessary marker metrics header " + MarkerMetrics.FULL_QC_HEADER[2] + " was not found in the marker metrics file" + proj.MARKER_METRICS_FILENAME.getValue());
				return false;
			} else {
				String[] metrics;
				while (reader.ready()) {
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
				Files.writeList(abMarkersToUse.toArray(new String[abMarkersToUse.size()]), proj.PROJECT_DIRECTORY.getValue() + MARKERS_FOR_ABCALLRATE);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
//			log.reportError("Error: file \"" + proj.getFilename(proj.MARKER_METRICS_FILENAME) + "\" not found in current directory");
			log.reportError("Error: file \"" + proj.MARKER_METRICS_FILENAME.getValue() + "\" not found in current directory");
		} catch (IOException ioe) {
//			log.reportError("Error reading file \"" + proj.getFilename(proj.MARKER_METRICS_FILENAME) + "\"");
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

		try {
			DNAIndex = getDNAIndex(proj, proj.PROJECT_DIRECTORY.getValue() + outputBase + PCA_SAMPLES_SUMMARY);
			reader = Files.getReader(proj.PROJECT_DIRECTORY.getValue() + outputBase + PCA_SAMPLES_SUMMARY, false, true, false);
			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("\t");
				qcLookup.put(line[DNAIndex], Array.toStr(Array.subArray(line, 1)));
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + proj.PROJECT_DIRECTORY.getValue() + PCA_SAMPLES_SUMMARY + "\" not found in current directory");
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + proj.PROJECT_DIRECTORY.getValue() + PCA_SAMPLES_SUMMARY + "\"");
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
	 * @param proj
	 *            current project
	 * @param outputBase
	 * @param markersForABCallRate
	 * @param markersForEverythingElse
	 * @param numThreads
	 *            threads for LRR_SD
	 * @param sampleCallRateFilter
	 *            filter samples by this call rate (LRR_SD filter is set in project)
	 * @param computeQC
	 * @param useFile
	 *            a further filter of samples that will be used
	 * @param log
	 */
	public static int[] filterSamples(Project proj, String outputBase, String markersForABCallRate, String markersForEverythingElse, int numThreads, String sampleCallRateFilter, String useFile) {
		Hashtable<String, String> sampDataQC = new Hashtable<String, String>();
		int[] indices;
		String[] line;
		int numPassing, count;
		Logger log = proj.getLog();

		SampleData sampleData = proj.getSampleData(0, false);
//		double lrrSdFilter = Double.parseDouble(proj.getProperty(proj.LRRSD_CUTOFF));
		double lrrSdFilter = proj.getProperty(proj.LRRSD_CUTOFF);
		double callRateFilter = Double.parseDouble(sampleCallRateFilter);
		boolean addToSampleData = checkSampleData(proj, sampleData);
		Hashtable<String, String> subset = checkSubset(useFile, log);

		if (Files.exists(proj.PROJECT_DIRECTORY.getValue() + DEFAULT_QC_FILE)) {
			log.report("The sample qc file " + proj.PROJECT_DIRECTORY.getValue() + DEFAULT_QC_FILE + " already exists");
			log.report("Skipping qc computation, filtering on existing qc file " + proj.PROJECT_DIRECTORY.getValue() + DEFAULT_QC_FILE);
		} else {
			log.report("Computing sample QC for all samples...");
			cnv.qc.LrrSd.init(proj, null, markersForABCallRate, markersForEverythingElse, null, numThreads);
		}

		count = 0;
		numPassing = 0;
		try {
			BufferedReader reader = Files.getReader(proj.PROJECT_DIRECTORY.getValue() + DEFAULT_QC_FILE, false, true, false);
			PrintWriter writerUse = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + outputBase + PCA_SAMPLES));
			PrintWriter writerSummary = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + outputBase + PCA_SAMPLES_SUMMARY));

			writerSummary.println(Array.toStr(SAMPLE_QC_SUMMARY));
			if (!reader.ready()) {
				writerUse.close();
				writerSummary.close();
				reader.close();
				log.reportError("Error - QC file (" + DEFAULT_QC_FILE + ") was empty");
				return new int[] { numPassing, count };
			}
			line = reader.readLine().trim().split(SPLITS[0]);
			indices = ext.indexFactors(QC_COLUMNS, line, true, log, true, false);

			if (!checkIndices(proj, indices)) {
				writerUse.close();
				writerSummary.close();
				reader.close();
				log.reportError("Error - could not detect proper header in QC file (" + DEFAULT_QC_FILE + ")");
				return null;
			}

			while (reader.ready()) {
				line = reader.readLine().trim().split(SPLITS[0]);
				// skip any headers as a result of concatenating the qc results
				if (!line[indices[0]].equals(QC_COLUMNS[0])) {
					// if passes qc
					if (Double.parseDouble(line[indices[1]]) < lrrSdFilter && Double.parseDouble(line[indices[2]]) > callRateFilter) {
						sampDataQC.put(line[indices[0]], line[indices[1]] + "\t" + line[indices[2]] + "\t" + "0");
						// check the subset
						if (subset.size() == 0 || subset.containsKey(line[indices[0]])) {
							writerUse.println(line[indices[0]]);
							writerSummary.println(line[indices[0]] + "\t" + line[indices[1]] + "\t" + line[indices[2]] + "\t" + "TRUE");
							numPassing++;
						} else {
							// sampDataQC.put(line[indices[0]], line[indices[1]] + "\t" + line[indices[2]] + "\t" + "1");
							writerSummary.println(line[indices[0]] + "\t" + line[indices[1]] + "\t" + line[indices[2]] + "\t" + "FALSE");
						}
					} else {
						sampDataQC.put(line[indices[0]], line[indices[1]] + "\t" + line[indices[2]] + "\t" + "1");
						writerSummary.println(line[indices[0]] + "\t" + line[indices[1]] + "\t" + line[indices[2]] + "\t" + "FALSE");
					}
					count++;
				}
			}

			reader.close();
			writerUse.close();
			writerSummary.close();

			if (numPassing == 0) {
				log.reportError("Error - all Samples were filtered out by the QC step");
				log.reportError("If there are a large number of cnv-only probes on the array, try lowering the call rate threshold for samples or use the \"-markerQC\" option to only select markers with high quality call rates");
				return new int[] { numPassing, count };
			}
			log.report("Info - " + numPassing + " " + (numPassing == 1 ? "sample" : "samples") + " passed the QC threshold" + (subset.size() > 0 ? " and were present in the subset file " + useFile : ""));
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + DEFAULT_QC_FILE + "\" not found in current directory");
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + DEFAULT_QC_FILE + "\"");
		}

		if (addToSampleData) {
			sampleData.addData(sampDataQC, DNA_LINKER, SAMPLE_DATA_ADDITION_HEADERS, ext.MISSING_VALUES[1], SPLITS[0], log);
		}
		return new int[] { numPassing, count };
	}

	/**
	 * If the useFile is not null, we return a hash with the subset of individuals. Can return null if useFile does not exist or does not contain any individuals
	 * 
	 * @param useFile
	 * @param log
	 * @return
	 */
	private static Hashtable<String, String> checkSubset(String useFile, Logger log) {
		Hashtable<String, String> subset = new Hashtable<String, String>();
		if (useFile != null) {
			if (Files.exists(useFile)) {
				subset = HashVec.loadFileToHashString(useFile, 0, new int[] { 0 }, null, false, false);
				if (subset.size() == 0) {
					log.reportError("Error - did not find any samples in the subset file " + useFile);
					return null;
				} else {
					log.report("Analysis will be performed starting with the subset of " + subset.size() + " samples found in " + useFile);
				}
			} else {
				log.reportError("Error - a file list of samples to use was provided, but the file " + useFile + " does not exist");
				return null;
			}
		} else {
			log.report("Info - A subset of samples was not provided with the \"useFile=\" argument, using all parsed samples as input to the SVD");
		}
		return subset;
	}

	/**
	 * Check to make sure that sample data has DNA header, and that the QC has not already been added
	 */
	private static boolean checkSampleData(Project proj, SampleData sampleData) {
		// This should not happen, but if it is the case we will not attempt to add qc metrics
		boolean addToSampleData = true;
		Logger log = proj.getLog();

		if (!sampleData.containsDNA()) {
			addToSampleData = false;
			log.reportError("Error - sample data did not contain column with header \"DNA\", not adding sample qc summaries to sample data");
		}
		if (qcAdded(proj)) {
			addToSampleData = false;
			log.reportError("Detected that sample data QC metrics have been added already, will not add these again");
//			log.reportError("If new thresholds were used, please remove the columns [" + ext.listWithCommas(SAMPLE_DATA_ADDITION_HEADERS, true) + "] in " + proj.getFilename(proj.SAMPLE_DATA_FILENAME));
			log.reportError("If new thresholds were used, please remove the columns [" + ext.listWithCommas(SAMPLE_DATA_ADDITION_HEADERS, true) + "] in " + proj.SAMPLE_DATA_FILENAME.getValue());
		}
		return addToSampleData;
	}

	/**
	 * Check all indices for -1 status
	 */
	private static boolean checkIndices(Project proj, int[] indices) {
		boolean allGood = true;
		for (int i = 0; i < indices.length; i++) {
			if (indices[i] == -1) {
				allGood = false;
				proj.getLog().reportError("Error - The sample QC file " + proj.PROJECT_DIRECTORY.getValue() + DEFAULT_QC_FILE + " did not contain the proper headings, this should not happen");
			}
		}
		return allGood;
	}

	/**
	 * Check the header of the sample data file to see if the sample data qc headers are present
	 */
	private static boolean qcAdded(Project proj) {
		boolean added = true;
//		String[] header = Files.getHeaderOfFile(proj.getFilename(proj.SAMPLE_DATA_FILENAME), proj.getLog());
		String[] header = Files.getHeaderOfFile(proj.SAMPLE_DATA_FILENAME.getValue(), proj.getLog());
		int[] indices = ext.indexFactors(SAMPLE_DATA_ADDITION_HEADERS, header, true, proj.getLog(), false, false);
		for (int i = 0; i < indices.length; i++) {
			if (indices[i] < 0) {
				added = false;
			}
		}
		return added;
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

	/**
	 * We use the Individual class as input so that we only need one method to generate the sample data
	 */
	public static void generateSampleData(Project proj, Individual[] inds) {
//		String sampleDataFile = proj.PROJECT_DIRECTORY.getValue() + proj.getProperty(proj.SAMPLE_DATA_FILENAME);
		String sampleDataFile = proj.SAMPLE_DATA_FILENAME.getValue(false, true);
		Logger log = proj.getLog();

		try {
			PrintWriter writer = new PrintWriter(new FileWriter(sampleDataFile));
			String[] classed = PED_INPUT;
			classed[5] = "Class=Sex";
			writer.println(Array.toStr(classed) + (inds[0].getSampleMapHeader() == null ? "" : "\t" + Array.toStr(inds[0].getSampleMapHeader())));
			for (int i = 0; i < inds.length; i++) {
				writer.println(inds[i].getSampDataFormat());
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + sampleDataFile + "\" could not be written to (it's probably open)");
			log.reportException(fnfe);
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + sampleDataFile + "\"");
			log.reportException(ioe);
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

	public static Individual[] loadPedFile(String pedFile, Logger log) {
		String[] line;
		ArrayList<Individual> al = new ArrayList<Individual>();
		try {
			BufferedReader reader = Files.getReader(pedFile, false, true, false);
			String temp = reader.readLine().trim();
			String delim = ext.determineDelimiter(temp);
			line = temp.split(delim);
			int[] indices = ext.indexFactors(PED_INPUT, line, true, false);
			for (int i = 0; i < indices.length; i++) {
				if (indices[i] < 0) {
					log.reportError("Error - Improper formatting of the pedigree file, can not generate sampleData");
					log.reportError("Warning - Parsing can proceed, but a sample data file is needed to generate principal components");
				}
			}
			while (reader.ready()) {
				line = reader.readLine().trim().split(delim);
				al.add(new Individual(line[indices[0]], line[indices[1]], line[indices[2]], line[indices[3]], line[indices[4]], line[indices[5]], line[indices[6]]));
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + pedFile + "\" not found in current directory");
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + pedFile + "\"");
		}
		return al.toArray(new Individual[al.size()]);
	}

	public static Individual[] loadSampleMapFile(String sampleMapCsv, Logger log) {
		String[] line;
		ArrayList<Individual> al = new ArrayList<Individual>();
		log.report("Using Sample Map file " + sampleMapCsv);
		try {
			BufferedReader reader = Files.getReader(sampleMapCsv, false, true, false);
			line = reader.readLine().trim().split(SPLITS[1]);
			String[] header = line;
			int[] indices = ext.indexFactors(SAMPLEMAP_INPUT, line, true, false);
			if (indices[1] == -1 || indices[2] == -1) {
				log.reportError("Error - Columns \"" + SAMPLEMAP_INPUT[1] + "\" and \"" + SAMPLEMAP_INPUT[2] + "\" must be provided in .csv format " + sampleMapCsv);
			}
			while (reader.ready()) {
				line = reader.readLine().trim().split(SPLITS[1]);
				al.add(new Individual(indices, line, header));
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + sampleMapCsv + "\" not found in current directory");
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + sampleMapCsv + "\"");
		}
		return al.toArray(new Individual[al.size()]);
	}

	/**
	 * A helper class to facilitate creating sample data from either .ped or Sample_Map.csv formats
	 * 
	 */
	public static class Individual {
		private String fid;
		private String iid;
		private String fa;

		private String mo;
		private String sex;
		private String aff;
		private String dna;
		private String[] sampleMapLine;
		private String[] sampleMapHeader;

		public Individual(String dna, String fid, String iid, String fa, String mo, String sex, String aff) {
			this.fid = fid;
			this.iid = iid;
			this.fa = fa;
			this.mo = mo;
			this.sex = sex;
			this.aff = aff;
			this.dna = dna;
		}

		/**
		 * @param indices
		 *            of the required columns
		 * @param sampleMapLine
		 *            a line from the sampleMap file
		 * @param header
		 */
		public Individual(int[] indices, String[] sampleMapLine, String[] header) {
			this.fid = sampleMapLine[indices[1]];
			this.iid = sampleMapLine[indices[2]];
			this.fa = "NA";
			this.mo = "NA";
			this.sex = indices[3] == -1 ? "NA" : parseSex(sampleMapLine[indices[3]]);
			this.dna = sampleMapLine[indices[2]];
			this.aff = "NA";
			this.sampleMapHeader = header;
			this.sampleMapLine = sampleMapLine;
		}

		public String getDna() {
			return dna;
		}

		public void setDna(String dna) {
			this.dna = dna;
		}

		public String parseSex(String sex) {
			String s = "-1";
			if (sex.toLowerCase().equals(SEX[0])) {
				s = "2";
			} else if (sex.toLowerCase().equals(SEX[1])) {
				s = "1";
			}
			return s;
		}

		public String[] getSampleMapHeader() {
			return sampleMapHeader;
		}

		public String getSampDataFormat() {
			if (sampleMapHeader == null) {
				return this.dna + "\t" + this.fid + "\t" + this.iid + "\t" + this.fa + "\t" + this.mo + "\t" + this.sex + "\t" + this.aff;
			} else {
				return this.dna + "\t" + this.fid + "\t" + this.iid + "\t" + this.fa + "\t" + this.mo + "\t" + this.sex + "\t" + this.aff + "\t" + Array.toStr(sampleMapLine);
			}
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String projectName = "GENVISIS_project1";

		boolean[] requiredArray = new boolean[5];
		Arrays.fill(requiredArray, false);
		String[] requiredArgs = { "dirProj=", "dirSrc=", PC_MARKER_COMMAND, MEDIAN_MARKER_COMMAND, "markerPositions=" };

		String projectDirectory = null;
		String sourceDirectory = null;
		String targetMarkers = null;
		String markerPositions = null;
		String medianMarkers = null;

		String output = "PCA_GENVISIS";
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

		int attempts, result;

		String usage = "\n";
		usage += "The MitoPipeline currently requires 5 arguments and allows for many more optional arguments:\n";
		usage += "  \n";
		usage += "   (1) The full path for the project directory (where results will be stored) (i.e. dirProj=/home/usr/projects/)\n";
		usage += "   (2) The full path for the source data directory  (where final report files are located) (i.e. dirSrc=/home/usr/data/project1/)\n";
		usage += "   (3) The full path for a file with a list of markers (one per line) to use for computing PCs (i.e. " + PC_MARKER_COMMAND + "/home/usr/auxFiles/exomeChip.PC_Markers.txt)\n";
		usage += "   (4) The full path for a file with a list of markers (one per line, typically mitochondrial markers) to use for computing computing median Log R Ratios (i.e. " + MEDIAN_MARKER_COMMAND + "/home/usr/auxFiles/exomeChip.MT_Markers.txt)\n";
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
		usage += "   (15) Output file baseName (i.e. output=" + output + " (default))\n";
		usage += "   (16) Project filename (if you manually created a project properties file, or edited an existing project). Note that default arguments available here can overide existing project properties (i.e. proj=" + filename + " (no default))\n";
		usage += "   (17) The header of the column containing sample ids in the final report files (for command-line interpretability, space characters must be replaced with \"_\". Common options are \"Sample_ID\" and \"Sample_Name\", corresponding to \"Sample ID\" and \"Sample Name\")  (i.e. idHeader=" + idHeader + " (default))\n";
		// usage += "   (18) A file specifying the AB allele lookup for markers, often times required  (i.e. abLookup=" + idHeader + " (default))\n";
		usage += "   (18) Do not perform a marker qc step to select higher quality markers (or remove cnv-only markers) to use for computing the sample call rate (i.e. -nomarkerQC (not the default))\n";
		usage += "   (19) If marker qc is performed, the call rate cutoff for markers to be passed on to the sample QC step (i.e. markerCallRate=" + markerCallRateFilter + " (default))\n";
		usage += "   (20) Name of the log file (i.e. log=[project_directory]/logs/Genvisis_[date].log (default))\n";
		usage += "   (21) Recompute Log R Ratios for each marker from genotypes/intensities when computing AND extrapolating PCs(i.e. recomputeLRR_PCs=" + recomputeLRR_PCs + " (default))\n";
		usage += "   (22) Recompute Log R Ratios for each marker from genotypes/intensities when computing median values(i.e. recomputeLRR_Median=" + recomputeLRR_Median + " (default))\n";
		usage += "   (23) Impute mean for markers with NaN data, if false markers with NaN values for any sample will be skipped (i.e. imputeMeanForNaN=" + imputeMeanForNaN + " (default))\n";

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
			} else if (args[i].startsWith(PC_MARKER_COMMAND)) {
				targetMarkers = ext.parseStringArg(args[i], null);
				numArgs--;
				requiredArray[2] = true;
			} else if (args[i].startsWith("output=")) {
				output = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith(MEDIAN_MARKER_COMMAND)) {
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

		proj = null;
		if (filename == null) {
			proj = initializeProject(proj, projectName, projectDirectory, sourceDirectory, dataExtension, idHeader, abLookup, targetMarkers, medianMarkers, markerPositions, sampleLRRSdFilter, sampleCallRateFilter, logfile);
		} else {
			proj = new Project(filename, logfile, false);
			proj = initializeProject(proj, projectName, projectDirectory, sourceDirectory, dataExtension, idHeader, abLookup, targetMarkers, medianMarkers, markerPositions, sampleLRRSdFilter, sampleCallRateFilter, logfile);
		}
		attempts = 0;
		while (attempts < 2) {
			result = catAndCaboodle(proj, numThreads, sampleCallRateFilter, medianMarkers, numComponents, output, homosygousOnly, markerQC, markerCallRateFilter, useFile, pedFile, sampleMapCsv, recomputeLRR_PCs, recomputeLRR_Median, doAbLookup, imputeMeanForNaN);
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
