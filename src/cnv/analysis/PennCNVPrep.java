package cnv.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

import cnv.analysis.pca.PrincipalComponentsIntensity;
import cnv.analysis.pca.PrincipalComponentsResiduals;
import cnv.filesys.MarkerData;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.manage.MarkerDataLoader;
import cnv.var.SampleData;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;

/**
 * This is a temporary fix to export corrected intensities, geared toward PennCNV output. It currently seems to work, but be weary
 *
 */
public class PennCNVPrep {
	private static final String[] PENN_STRINGS = { "Name", ".GType", ".Log R Ratio", ".B Allele Freq" };
	private static final String STORAGE_BASE = "firstMarkerIndex_";
	private static final String STORAGE_EXT = ".ser";

	private Project proj;
	private PrincipalComponentsResiduals principalComponentsResiduals;
	private boolean[] samplesToExport, samplesToUseCluster;
	private int numCorrectionThreads, numComponents;
	private int[] sampleSex;
	private String[] markers;
	private String dir;

	public PennCNVPrep(Project proj, PrincipalComponentsResiduals principalComponentsResiduals, boolean[] samplesToExport, boolean[] samplesToUseCluster, int[] sampleSex, String[] markers, int numComponents, String dir, int numThreads) {
		super();
		this.proj = proj;
		this.principalComponentsResiduals = principalComponentsResiduals;
		this.samplesToExport = samplesToExport;
		this.samplesToUseCluster = samplesToUseCluster;
		this.sampleSex = sampleSex;
		this.numComponents = numComponents;
		this.numCorrectionThreads = numThreads;
		this.markers = markers;
		this.dir = dir;
	}

	/**
	 * This creates the temporary serialized {@link MarkerData} lists stored in {@link MarkerDataStorage} objects, and contain only LRR/BAF and genotypes (currently original genotypes)
	 * 
	 */
	public void exportSpecialMarkerData() {
		MarkerDataLoader markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markers);
		MarkerDataStorage markerDataStorage = new MarkerDataStorage(markers.length);
		ArrayList<String> notCorrected = new ArrayList<String>();
		for (int i = 0; i < markers.length; i++) {
			if ((i + 1) % 10 == 0) {
				proj.getLog().report("Info - correcting marker " + (i + 1) + " of " + markers.length);
			}
			// TODO, if we have more than 6 threads we could probably do that here
			MarkerData markerData = markerDataLoader.requestMarkerData(i);
			PrincipalComponentsIntensity principalComponentsIntensity = new PrincipalComponentsIntensity(principalComponentsResiduals, markerData, true, sampleSex, samplesToUseCluster, 1, 0, null, true, false, 2, 5, 0, numCorrectionThreads, false, null);
			byte[] genotypes = markerData.getAbGenotypes();
			markerDataLoader.releaseIndex(i);
			principalComponentsIntensity.correctXYAt(numComponents);
			if (principalComponentsIntensity.isFail()) {
				notCorrected.add(markers[i]);
			}
			float[][] correctedLRRBAF = principalComponentsIntensity.getCorrectedLRRBAF(true);// for now
			MarkerData markerDataToStore = new MarkerData(markerData.getMarkerName(), markerData.getChr(), markerData.getPosition(), markerData.getFingerprint(), null, null, null, null, null, null, null, correctedLRRBAF[0], correctedLRRBAF[1], genotypes, genotypes);
			markerDataStorage.addToNextIndex(markerDataToStore);
		}
		String output = proj.getProjectDir() + dir + STORAGE_BASE + ext.indexLargeFactors(markers, proj.getMarkerNames(), true, proj.getLog(), true, true)[0] + STORAGE_EXT;
		markerDataStorage.serialize(output);
		if (notCorrected.size() > 0) {
			Files.writeList(notCorrected.toArray(new String[notCorrected.size()]), output.replaceAll("\\.ser", "_") + notCorrected.size() + "_markersThatFailedCorrection.txt");
		}
	}

	/**
	 * This function exports to PennCNV format starting from the temporary files made in {@link PennCNVPrep#exportSpecialMarkerData()} Currently, all individuals will be exported
	 * 
	 * @param fileNamesOfMarkerDataInOrder
	 *            files of serialized {@link MarkerDataStorage};
	 */
	public void exportSpecialPennCNVData(String[] fileNamesOfMarkerDataInOrder) {
		int[] sampleIndicesInProject = ext.indexLargeFactors(Array.subArray(proj.getSamples(), samplesToExport), proj.getSamples(), true, proj.getLog(), true, true);
		int numMarkersPerWrite = Integer.parseInt(proj.getProperty(Project.MAX_MARKERS_LOADED_PER_CYCLE));
		int numMarkersThisRound = 0;
		String[] subSamples = Array.subArray(proj.getSamples(), samplesToExport);
		PennCNVIndividual[] pennCNVIndividuals = initSamples(proj, dir, subSamples, false, numMarkersPerWrite, proj.getLog());
		for (int i = 0; i < fileNamesOfMarkerDataInOrder.length; i++) {
			MarkerDataStorage markerDataStorage = MarkerDataStorage.load(fileNamesOfMarkerDataInOrder[i], false);
			MarkerData[] markerDatas = markerDataStorage.getMarkerDatas();
			for (int j = 0; j < markerDatas.length; j++) {
				if ((j + 1) % 100 == 0) {
					proj.getLog().report("Info - exporting marker " + (j + 1) + " of " + markerDatas.length + " from file " + fileNamesOfMarkerDataInOrder[i]);
				}
				addData(numMarkersThisRound, pennCNVIndividuals, markerDatas[j].getBAFs(), markerDatas[j].getLRRs(), markerDatas[j].getAbGenotypes(), markerDatas[j].getMarkerName(), sampleIndicesInProject);
				numMarkersThisRound++;
				if (numMarkersThisRound == numMarkersPerWrite || (j == markerDatas.length - 1 && i == fileNamesOfMarkerDataInOrder.length - 1)) {
					numMarkersThisRound = 0;
					for (int k = 0; k < pennCNVIndividuals.length; k++) {
						pennCNVIndividuals[k].dump(proj);
					}
					pennCNVIndividuals = initSamples(proj, dir, subSamples, true, numMarkersPerWrite, proj.getLog());
				}
			}
		}
	}

	private static void addData(int index, PennCNVIndividual[] pennCNVIndividuals, float[] baf, float[] lrr, byte[] genotypes, String marker, int[] sampleIndicesInProject) {
		for (int i = 0; i < sampleIndicesInProject.length; i++) {
			pennCNVIndividuals[i].addLine(index, marker, genotypes[sampleIndicesInProject[i]], lrr[sampleIndicesInProject[i]], baf[sampleIndicesInProject[i]]);
		}
	}

	/**
	 * Currently only stores new LRR/BAF values and original genotypes
	 *
	 */
	private static class MarkerDataStorage implements Serializable {
		private static final long serialVersionUID = 1L;
		private MarkerData[] markerDatas;
		private int currentIndex;

		public MarkerDataStorage(int numMarkers) {
			this.markerDatas = new MarkerData[numMarkers];
			this.currentIndex = 0;
		}

		public void addToNextIndex(MarkerData markerData) {
			markerDatas[currentIndex] = markerData;
			currentIndex++;
		}

		public void serialize(String filename) {
			Files.writeSerial(this, filename);
		}

		public static MarkerDataStorage load(String filename, boolean jar) {
			return (MarkerDataStorage) Files.readSerial(filename, jar, true);
		}

		public MarkerData[] getMarkerDatas() {
			return markerDatas;
		}

	}

	/**
	 * Aids in processing to PennCNV output when iterating over markers
	 *
	 */
	private static class PennCNVIndividual {
		private String ouputFile, header;
		private String[] markers;
		private boolean append;
		private byte[] genotypes;
		private float[] lrrs;
		private float[] bafs;
		private int numAdded;

		public PennCNVIndividual(String outputFile, String header, int numMarkers, boolean append) {
			this.ouputFile = outputFile;
			this.header = header;
			this.markers = new String[numMarkers];
			this.genotypes = new byte[numMarkers];
			this.lrrs = new float[numMarkers];
			this.bafs = new float[numMarkers];
			this.append = append;
			this.numAdded = 0;
		}

		public void addLine(int projectIndex, String marker, byte genotype, float lrr, float baf) {
			markers[projectIndex] = marker;
			genotypes[projectIndex] = genotype;
			lrrs[projectIndex] = lrr;
			bafs[projectIndex] = baf;
			numAdded++;
		}

		private String getStringAt(int index) {
			return markers[index] + "\t" + (genotypes[index] == -1 ? "NC" : Sample.AB_PAIRS[genotypes[index]]) + "\t" + lrrs[index] + "\t" + bafs[index];
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
	private static PrincipalComponentsResiduals loadPcResids(Project proj, int numComponents) {
		String pcFile = proj.getFilename(Project.INTENSITY_PC_FILENAME);
		PrincipalComponentsResiduals pcResids;
		if (Files.exists(proj.getProjectDir() + ext.removeDirectoryInfo(pcFile))) {
			proj.getLog().report("Info - loading " + ext.removeDirectoryInfo(pcFile));
			pcResids = new PrincipalComponentsResiduals(proj, ext.removeDirectoryInfo(pcFile), null, Integer.parseInt(proj.getProperty(Project.INTENSITY_PC_NUM_COMPONENTS)), false, 0, false, false, null);
		} else {
			proj.getLog().reportError("Error - did not find Intensity PC File " + proj.getProjectDir() + ext.removeDirectoryInfo(pcFile) + " as defined by" + Project.INTENSITY_PC_FILENAME);
			pcResids = null;
		}
		return pcResids;
	}

	private static PennCNVIndividual[] initSamples(Project proj, String dir, String[] samplesToExport, boolean append, int numMarkers, Logger log) {
		PennCNVIndividual[] samples = new PennCNVIndividual[samplesToExport.length];
		for (int i = 0; i < samplesToExport.length; i++) {
			String header = PENN_STRINGS[0] + "\t" + samplesToExport[i] + PENN_STRINGS[1] + "\t" + samplesToExport[i] + PENN_STRINGS[2] + "\t" + samplesToExport[i] + PENN_STRINGS[3];
			samples[i] = new PennCNVIndividual(samplesToExport[i] + ".txt", header, numMarkers, append);
		}
		return samples;
	}

	private static int[] getSampleSex(Project proj) {
		String[] samples = proj.getSamples();
		SampleData sampleData = proj.getSampleData(0, false);
		int[] sex = new int[samples.length];
		for (int i = 0; i < samples.length; i++) {
			sex[i] = sampleData.getSexForIndividual(samples[i]);
		}
		int numMales = Array.countIf(sex, 1);
		int numFemales = Array.countIf(sex, 2);
		double percentDefined = (double) (numFemales + numMales) / samples.length;
		if (percentDefined > .90) {
			proj.getLog().report("Info - detected " + numMales + " males and " + numFemales + " females to use for sex-specific reclustering");
		}
		if (numFemales == 0) {
			proj.getLog().report("Error - no females were detected in sample data, chr X will not be reclustered properly");
			return null;
		}
		if (numMales == 0) {
			proj.getLog().report("Error - no males were detected in sample data, chr Y will not be reclustered properly");
			return null;
		}
		return sex;
	}

	private static String[] getSortedFileNames(Project proj, String dir) {
		String[] markers = proj.getMarkerNames();
		ArrayList<String> files = new ArrayList<String>();
		for (int i = 0; i < markers.length; i++) {
			String possibleExist = proj.getProjectDir() + dir + STORAGE_BASE + i + STORAGE_EXT;
			if (Files.exists(possibleExist)) {
				files.add(possibleExist);
			}
		}
		proj.getLog().report("Info - detected " + files.size() + " files to use for the PennCNV export");
		return files.toArray(new String[files.size()]);
	}

	/**
	 * @param proj
	 * @param dir
	 *            directory under the project directory
	 * @param numComponents
	 *            number of components to correct for
	 * @param markerFile
	 *            path (full path) to a list of markers (single column) to correct for
	 * @param numThreads
	 *            number of threads for each correction, maximum of 6 will be used
	 * @param exportFromSpecialMarkerData
	 *            flag if the directory supplied has serialized files already
	 */
	public static void exportSpecialPennCNV(Project proj, String dir, int numComponents, String markerFile, int numThreads, boolean exportFromSpecialMarkerData) {
		String[] markers;
		new File(proj.getProjectDir() + dir).mkdirs();
		if (!exportFromSpecialMarkerData) {
			PrincipalComponentsResiduals principalComponentsResiduals = loadPcResids(proj, numComponents);
			if (principalComponentsResiduals == null) {
				return;
			}
			int[] sex = getSampleSex(proj);
			if (sex == null) {
				proj.getLog().reportError("Error - missing sex codes");
				return;
			}
			if (markerFile == null) {
				markers = proj.getMarkerNames();
				proj.getLog().report("Info - a file of markers was not provided, exporting all in this batch");
			} else {
				markers = HashVec.loadFileToStringArray(markerFile, false, new int[] { 0 }, false);
				proj.getLog().report("Info - loaded " + markers.length + " markers from " + markerFile + " to export");

			}
			PennCNVPrep specialPennCNVFormat = new PennCNVPrep(proj, principalComponentsResiduals, null, proj.getSamplesToInclude(null), sex, markers, numComponents, dir, numThreads);
			specialPennCNVFormat.exportSpecialMarkerData();
		} else {
			boolean[] exportThese = new boolean[proj.getSamples().length];
			Arrays.fill(exportThese, true);
			PennCNVPrep specialPennCNVFormat = new PennCNVPrep(proj, null, exportThese, null, null, null, numComponents, dir, numThreads);
			String[] sortedFileNames = getSortedFileNames(proj, dir);
			if (sortedFileNames == null || sortedFileNames.length == 0) {
				proj.getLog().reportError("Error - did not find any files to export");

			} else {
				specialPennCNVFormat.exportSpecialPennCNVData(sortedFileNames);
			}
		}
	}

	public static void batchCorrections(Project proj, String java, String classPath, int memoryInMB, int wallTimeInHours, String dir, int numBatches, int numThreads, int numComponents) {
		String[] allMarkers = proj.getMarkerNames();
		int[] chunks = Array.splitUp(allMarkers.length, numBatches);
		int index = 0;
		String[][] batches = new String[numBatches][1];
		new File(proj.getProjectDir() + dir).mkdirs();
		for (int i = 0; i < chunks.length; i++) {
			ArrayList<String> chunk = new ArrayList<String>(chunks[i]);
			for (int j = 0; j < chunks[i]; j++) {
				chunk.add(allMarkers[index]);
				index++;
			}
			batches[i][0] = "batch_" + i + "_" + chunks[i] + "_markers";
			Files.writeList(chunk.toArray(new String[chunk.size()]), proj.getProjectDir() + dir + batches[i][0] + ".txt");
		}
		String command = "module load java\njava -cp  " + classPath + " -Xmx" + memoryInMB + "m cnv.analysis.PennCNVPrep proj=" + proj.getFilename(Project.PROJECT_PROPERTIES_FILENAME) + " dir=" + dir;
		Files.qsub("PennCNVPrepFormatExport", command + " -create", new String[][] { { "" } }, memoryInMB, 3 * wallTimeInHours, 1);
		command += " numThreads=" + numThreads + " numComponents=" + numComponents + " markers=" + proj.getProjectDir() + dir + "[%0].txt";
		Files.qsub("PennCNVPrepFormatTmpFiles", command, batches, memoryInMB, wallTimeInHours, numThreads);
		if (!Files.exists(proj.getFilename(Project.INTENSITY_PC_FILENAME))) {
			proj.getLog().report("Warning - all jobs will fail if the property " + Project.INTENSITY_PC_FILENAME + " in " + proj.getPropertyFilename() + " is not set to an existing file");
			proj.getLog().report("		  - did not find " + proj.getFilename(Project.INTENSITY_PC_FILENAME));
		}
		if (getSampleSex(proj) == null) {
			proj.getLog().report("Warning - all jobs will fail if sample sex is not provided in " + proj.getFilename(Project.SAMPLE_DATA_FILENAME));
			proj.getLog().report("		  - please specify sex for as many individuals as possible");
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String classPath = "/home/pankrat2/lanej/park.jar";
		String java = "/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/bin/java";// for lab
		// String java ="/soft/java/jdk1.7.0_45/bin/java"; //for itasca
		int memoryInMB = 22000;
		int wallTimeInHours = 8;
		String filename = null;
		String logfile = null;
		String dir = "PennCNVPrep/";
		int numThreads = 6;// can only utilize 6 (3 genotype clusters by X/Y)
		int numComponents = 40;
		String markers = null;
		boolean exportToPennCNV = false;
		int batch = 0;

		// Ex - Recommend modifying this to run the corrections

		// java -cp park.jar cnv.analysis.PennCNVPrep batch=100 proj=/home/usr/projects/x.properties classPath=/yourPathTo/park.jar dir=PennCNVPrep/ numComponents=40

		// then run recommended
		// ./master.PennCNVPrepFormatTmpFiles

		// and after completion of all batches
		// ./master.PennCNVPrepFormatExport

		// also, if you submit a bunch of jobs, this will kill them all
		// showq -u Your_Usr_Name_Here |cut -d ' ' -f 1|xargs qdel

		String usage = "\n" + "cnv.analysis.PennCNVPrep requires 1 argument\n";
		usage += "   (1) Project (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) logfile (i.e. log=" + logfile + " ( no default))\n" + "";
		usage += "   (3) directory (relative to the project directory) for output (i.e. dir=" + dir + " ( no default))\n" + "";
		usage += "   (4) number of threads, no more than 6 will used (i.e. numThreads=" + numThreads + " (default))\n" + "";
		usage += "   (5) number of principal components for correction, (i.e. numComponents=" + numComponents + " (default))\n" + "";
		usage += "   (6) a full path to a file listing markers to export in the current batch, (i.e. markers=" + numComponents + " (default))\n" + "";
		usage += "   (7) create PennCNV files from the tempory markerData files, (i.e. -create ( not the default))\n" + "";
		usage += "   (8) set this up for a batch run, which is recommended. Set to 0 if batch is not wanted (i.e. batch=" + batch + " (default))\n" + "";
		// usage += "   (9) java location for batch run (i.e. java=" + java + " (default))\n" + "";
		usage += "   (9) classPath for batch run (i.e. classPath=" + classPath + " (default))\n" + "";

		usage += "   NOTE: aprox 50 *(numSamples/5000) batches per 500,000 markers" + "";
		usage += "   NOTE: If using ChrX/ChrY markers, it is important that the projects sample data file has sex defined for all samples that are used for clustering" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("dir=")) {
				dir = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("numThreads=")) {
				numThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("numComponents=")) {
				numComponents = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("markers=")) {
				markers = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("java=")) {
				java = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("classPath=")) {
				classPath = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("-create")) {
				exportToPennCNV = true;
				numArgs--;
			} else if (args[i].startsWith("batch=")) {
				batch = ext.parseIntArg(args[i]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			Project proj = new Project(filename, logfile, false);
			if (batch > 0) {
				batchCorrections(proj, java, classPath, memoryInMB, wallTimeInHours, dir, batch, numThreads, numComponents);
			} else {
				exportSpecialPennCNV(proj, dir, numComponents, markers, numThreads, exportToPennCNV);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
