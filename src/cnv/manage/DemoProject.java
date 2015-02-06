package cnv.manage;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;

import common.Array;
import common.Files;
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
		this.proj = (Project) proj.clone();
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
		setProperty(Project.PROJECT_DIRECTORY, demoProjectDirectory);
		if (!Files.exists(getProjectDir()) || overwriteExisting) {
			new File(getProjectDir()).mkdirs();
			getDir(Project.DATA_DIRECTORY, true, false);
			getDir(Project.SAMPLE_DIRECTORY, true, false);
			if (dType == DEMO_TYPE.MARKER_FOCUS) {
				getDir(Project.MARKER_DATA_DIRECTORY, true, false);
			}
			copyFileIfExists(Project.MARKER_POSITION_FILENAME);
			copyFileIfExists(Project.SAMPLE_DATA_FILENAME);
			copyFileIfExists(Project.SAMPLE_QC_FILENAME);
			// copyFileIfExists(Project.DISPLAY_MARKERS_FILENAME);
			copyFileIfExists(Project.INTENSITY_PC_FILENAME);

		} else {
			proj.getLog().reportTimeError(demoProjectDirectory + " exists and the overwrite option was not flagged, halting");
			fail = true;
		}
	}

	private void copyFileIfExists(String property) {
		if (Files.exists(proj.getFilename(property))) {
			Files.copyFile(proj.getFilename(property, false, false), getFilename(property, false, false));
		} else {
			proj.getLog().reportTimeWarning("Did not find file " + proj.getFilename(property, false, false) + " cannot copy to demo");
		}
	}

	private static boolean[] loadSamples(Project proj, String samplesFile, DEMO_TYPE dType) {
		boolean[] samplesToUse = new boolean[proj.getSamples().length];
		Arrays.fill(samplesToUse, true);// defualt to all samples;
		String sampleSubsetFile = proj.getFilename(Project.SAMPLE_SUBSET_FILENAME);
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
		String targetMarkerFile = proj.getFilename(Project.TARGET_MARKERS_FILENAME);
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
				Markers.orderMarkers(markersToExport, getFilename(Project.MARKER_POSITION_FILENAME), getFilename(Project.MARKERSET_FILENAME, true, true), proj.getLog());
				long fingerPrint = getMarkerSet().getFingerprint();
				proj.getLog().reportTimeInfo("Attempting to subset the samples...");
				new File(getFilename(SAMPLELIST_FILENAME, false, false)).delete();
				created = FocusedSample.focusAProject(proj, this, markersToExport, samplesToUse, fingerPrint, numThreads, overwriteExisting, getLog());
				if (created) {
					if (dType == DEMO_TYPE.MARKER_FOCUS) {
						proj.getLog().reportTimeInfo("Finished subsetting the samples...Attempting to transpose the data");
						TransposeData.transposeData(this, 2000000000, false);
						Files.writeList(markersToExport, getFilename(Project.DISPLAY_MARKERS_FILENAME));
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
