package cnv.manage;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;

import common.Array;
import common.Files;
import common.ext;
import cnv.filesys.Project;
import cnv.filesys.SampleList;
import cnv.manage.ExtProjectDataParser.Builder;

public class DemoProject extends Project {
	private static final long serialVersionUID = 1L;
	private Project proj;

	public static enum DEMO_TYPE {
		MARKER_FOCUS, SAMPLE_FOCUS
	}

	/**
	 * where the new project directory will be placed
	 */
	private String demoDirectory;
	private boolean overwriteExisting, fail;
	private DEMO_TYPE dType;

	public DemoProject(Project proj, String demoDirectory, boolean overwriteExisting, DEMO_TYPE dType) {
		super();
		this.proj = new Project();//(Project) proj.clone();
		this.demoDirectory = demoDirectory;
		this.overwriteExisting = overwriteExisting;
		this.dType = dType;
		this.fail = false;
		init();
	}

	public boolean isFail() {
		return fail;
	}

	public DEMO_TYPE getdType() {
		return dType;
	}

	private void init() {
		String demoProjectDirectory = demoDirectory + proj.getNameOfProject() + "_" + dType + "/";
		setProperty(proj.PROJECT_DIRECTORY, demoProjectDirectory);
		if (!Files.exists(getProjectDir()) || overwriteExisting) {
			new File(getProjectDir()).mkdirs();
			getDir(proj.DATA_DIRECTORY, true, false);
			getDir(proj.SAMPLE_DIRECTORY, true, false);
			getDir(proj.RESULTS_DIRECTORY, true, false);

			if (dType == DEMO_TYPE.MARKER_FOCUS) {
				getDir(proj.MARKER_DATA_DIRECTORY, true, false);
			}

			// single file Copies
//			copyFileIfExists(proj.MARKER_POSITION_FILENAME);
//			copyFileIfExists(proj.SAMPLE_DATA_FILENAME);
//			copyFileIfExists(proj.SAMPLE_QC_FILENAME);
//			copyFileIfExists(proj.MOSAIC_RESULTS_FILENAME);
//			copyFileIfExists(proj.SEXCHECK_RESULTS_FILENAME);
//			copyFileIfExists(proj.INTENSITY_PC_FILENAME);
//			copyFileIfExists(proj.CLUSTER_FILTER_COLLECTION_FILENAME);
//			copyFileIfExists(proj.ANNOTATION_FILENAME);
			copyFileIfExists(proj.MARKER_POSITION_FILENAME.getName());
			copyFileIfExists(proj.SAMPLE_DATA_FILENAME.getName());
			copyFileIfExists(proj.SAMPLE_QC_FILENAME.getName());
			copyFileIfExists(proj.MOSAIC_RESULTS_FILENAME.getName());
			copyFileIfExists(proj.SEXCHECK_RESULTS_FILENAME.getName());
			copyFileIfExists(proj.INTENSITY_PC_FILENAME.getName());
			copyFileIfExists(proj.CLUSTER_FILTER_COLLECTION_FILENAME.getName());
			copyFileIfExists(proj.ANNOTATION_FILENAME.getName());
			// multi-file Copies
//			copyIfFilesExists(proj.TWOD_LOADED_FILENAMES);
			copyIfFilesExists(proj.TWOD_LOADED_FILENAMES.getName());
			
//			copyIfFilesExists(proj.INDIVIDUAL_CNV_LIST_FILENAMES);//copying to both marker and sample focus in case of marker focused cnvs
			copyIfFilesExists(proj.INDIVIDUAL_CNV_LIST_FILENAMES.getName());//copying to both marker and sample focus in case of marker focused cnvs

			copyStratResults(proj, this);
			copyGeneTrack(proj, this);

		} else {
			proj.getLog().reportTimeError(demoProjectDirectory + " exists and the overwrite option was not flagged, halting");
			fail = true;
		}
	}

	private static void copyGeneTrack(Project original, Project demo) {
		String geneTrack = original.getGeneTrackFileName(true);
		if (geneTrack != null && Files.exists(geneTrack)) {
			demo.setProperty(demo.GENETRACK_FILENAME, ext.removeDirectoryInfo(geneTrack));
			original.getLog().reportTimeInfo("Detected " + geneTrack + ", copying to " + demo.getFilename(demo.GENETRACK_FILENAME, false, false) + "\n\t (this takes a while due to byte by byte copying)");
			Files.copyFileExactlyByteByByte(geneTrack, demo.getFilename(demo.GENETRACK_FILENAME, false, false));
		}
	}

	private static void copyStratResults(Project original, Project demo) {
		String[] strats = Array.toStringArray(original.getStratResults());
		if (strats != null && strats.length > 0) {
			for (int i = 0; i < strats.length; i++) {
				Files.copyFile(strats[i], demo.getProjectDir() + ext.removeDirectoryInfo(strats[i]));
			}
		} else {
			original.getLog().reportTimeWarning("Did not find any stratification results");
		}

	}

	private void copyIfFilesExists(String propertyWithMultipleFiles) {
//		String propertyValue = proj.getProperty(propertyWithMultipleFiles);
		String propertyValue = proj.getProperty(propertyWithMultipleFiles).getValueString();
		String propertyValueNew = propertyValue;
		propertyValueNew = propertyValueNew.replace(proj.getProjectDir(), getProjectDir());
		setProperty(propertyWithMultipleFiles, propertyValueNew);

		String[] filesOriginal = proj.getFilenames(propertyWithMultipleFiles);
		String[] filesNew = getFilenames(propertyWithMultipleFiles, true);
		copyFiles(filesOriginal, filesNew);
	}

	private void copyFiles(String[] filesOriginal, String[] filesNew) {
		if (filesOriginal != null && filesOriginal.length > 0) {
			for (int i = 0; i < filesOriginal.length; i++) {
				if (filesOriginal[i] != null && Files.exists(filesOriginal[i])) {
					Files.copyFile(filesOriginal[i], filesNew[i]);
				} else {
					proj.getLog().reportTimeWarning("Did not find file " + filesOriginal[i] + " cannot copy to demo");
				}
			}
		}
	}

	private void copyFileIfExists(String property) {
		String file = proj.getFilename(property, false, false);
		if (file != null && Files.exists(file)) {
			Files.copyFile(file, getFilename(property, false, false));
		} else {
			proj.getLog().reportTimeWarning("Did not find file " + file + " cannot copy to demo");
		}
	}

	private static boolean[] loadSamples(Project proj, String samplesFile, DEMO_TYPE dType) {
		boolean[] samplesToUse = new boolean[proj.getSamples().length];
		Arrays.fill(samplesToUse, true);// defualt to all samples;
		String sampleSubsetFile = proj.getFilename(proj.SAMPLE_SUBSET_FILENAME);
		if (samplesFile != null) {
			samplesToUse = getParserFor(proj, samplesFile, true).getDataPresent();
		} else if (Files.exists(sampleSubsetFile) && dType == DEMO_TYPE.SAMPLE_FOCUS) {
			samplesToUse = getParserFor(proj, sampleSubsetFile, true).getDataPresent();
		} else {
			if (dType == DEMO_TYPE.SAMPLE_FOCUS) {
				proj.getLog().reportTimeWarning("A sample subset file was not provided and " + sampleSubsetFile + " did not exist, exporting all samples");
			} else {
				proj.getLog().reportTimeInfo("A sample subset file was not provided, exporting all samples for demotype " + DEMO_TYPE.MARKER_FOCUS);
			}
		}
		return samplesToUse;
	}

	private static String[] loadMarkers(Project proj, String markersFile, DEMO_TYPE dType) {
		String[] markersToUse = proj.getMarkerNames();
		String targetMarkerFile = proj.getFilename(proj.TARGET_MARKERS_FILENAME);
		if (markersFile != null) {
			markersToUse = getParserFor(proj, markersFile, false).getStringDataAt(0, true);
		} else if (Files.exists(targetMarkerFile) && dType == DEMO_TYPE.MARKER_FOCUS) {
			markersToUse = getParserFor(proj, targetMarkerFile, false).getStringDataAt(0, true);
		} else {
			if (dType == DEMO_TYPE.MARKER_FOCUS) {
				proj.getLog().reportTimeWarning("A marker subset file was not provided and " + markersFile + " did not exist, exporting all markers");
			} else {
				proj.getLog().reportTimeInfo("A marker subset file was not provided, exporting all markers for demotype " + DEMO_TYPE.SAMPLE_FOCUS);

			}
		}
		return markersToUse;
	}

	private static ExtProjectDataParser getParserFor(Project proj, String file, boolean sampleBased) {
		ExtProjectDataParser.Builder builder = new Builder();
		builder.dataKeyColumnIndex(0);
		builder.stringDataIndices(new int[][] { { 0 } });
		builder.treatAllNumeric(false);
		builder.hasHeader(false);
		builder.skipUnFound(false);
		builder.sampleBased(sampleBased);
		ExtProjectDataParser sampleParser = null;
		try {
			sampleParser = builder.build(proj, file);
		} catch (FileNotFoundException e) {
			proj.getLog().reportException(e);
			e.printStackTrace();
		}
		sampleParser.loadData();
		return sampleParser;
	}

	public boolean createProjectDemo(String markersFile, String samplesFile, int numThreads) {
		if (!fail) {
			return createProjectDemo(loadMarkers(proj, markersFile, dType), loadSamples(proj, samplesFile, dType), numThreads);
		} else {
			return false;
		}
	}

	public boolean createProjectDemo(String[] markersToExport, boolean[] samplesToUse, int numThreads) {

		boolean created = true;
		if (!fail) {
			if (markersToExport == null) {
				markersToExport = proj.getMarkerNames();
				proj.getLog().reportTimeWarning("Since the markers to export were not provided, all " + markersToExport.length + " markers will be exported");
			}
			if (samplesToUse == null) {
				samplesToUse = new boolean[proj.getSamples().length];
				Arrays.fill(samplesToUse, true);
				proj.getLog().reportTimeWarning("Since the samples to export were not provided, all " + samplesToUse.length + " samples will be exported");
			}
			if (samplesToUse.length != proj.getSamples().length) {
				proj.getLog().reportTimeError("The array length provided does (" + samplesToUse.length + ") does not contain boolean values for all samples");
				created = false;
				return created;
			} else {
				Markers.orderMarkers(markersToExport, getFilename(proj.MARKER_POSITION_FILENAME), getFilename(proj.MARKERSET_FILENAME, true, true), proj.getLog());
				long fingerPrint = getMarkerSet().getFingerprint();
				proj.getLog().reportTimeInfo("Attempting to subset the samples...");
				new File(getFilename(SAMPLELIST_FILENAME, false, false)).delete();
				created = FocusedSample.focusAProject(proj, this, markersToExport, samplesToUse, fingerPrint, numThreads, overwriteExisting, getLog());
				if (created) {
					if (dType == DEMO_TYPE.MARKER_FOCUS) {
						proj.getLog().reportTimeInfo("Finished subsetting the samples...Attempting to transpose the data");
						TransposeData.transposeData(this, 2000000000, false);
						Files.writeList(markersToExport, getFilename(proj.DISPLAY_MARKERS_FILENAME));
					}
					SampleList.generateSampleList(this).writeToTextFile(getProjectDir() + "ListOfSamples.txt");
					getSamples();
				} else {
					proj.getLog().reportTimeError("Could not create focused samples, halting...");
				}
			}
		} else {
			created = false;
		}
		return created;
	}

	public static void test() {
		Project proj = new Project(null, false);
		DemoProject demoProject = new DemoProject(proj, proj.getProjectDir() + "Demo/", true, DEMO_TYPE.MARKER_FOCUS);
		String[] adfa = null;
		demoProject.createProjectDemo(adfa, null, 8);
	}

	public static void main(String[] args) {
		test();
	}

}
