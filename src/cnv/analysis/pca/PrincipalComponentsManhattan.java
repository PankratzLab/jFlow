package cnv.analysis.pca;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import stats.CrossValidation;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.Matrix;
import common.PSF;
import common.ext;
import cnv.filesys.MarkerData;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.manage.MarkerDataLoader;

/**
 * Generates manhattan plots for PCs and other optional data...can be loaded directly into IGV or haploview
 *
 */
public class PrincipalComponentsManhattan extends PrincipalComponentsResiduals {
	/**
	 * 
	 */
	private static final String[] HEADER = { "Index", "SNP", "CHR", "BP", "P", "T_STAT" };
	private static final String EXT = ".hat.linear";

	private String[] markersToTest;
	private double[][][] results;// type,marker,result
	private ManhattanTest[] manhattanTests;
	private int numPCs;

	/**
	 * @param proj
	 * @param markersToTest
	 *            these markers will be tested against the pcs, and any other data provided
	 * @param numPcs
	 */
	public PrincipalComponentsManhattan(Project proj, String[] markersToTest, int numPcs) {
		super(proj.loadPcResids());
		this.markersToTest = PrincipalComponentsCompute.sortByProjectMarkers(getProj(), markersToTest);
		this.manhattanTests = null;
		this.results = null;
		this.numPCs = numPcs > 0 ? numPcs : getNumComponents();

	}

	/**
	 * @param mTests
	 *            can add more tests here,just call before populating results
	 */
	public void addManhattanTest(ManhattanTest[] mTests) {
		manhattanTests = Array.concatAll(manhattanTests, mTests);
		updateResultSize();
	}

	/**
	 * computes for all markers
	 */
	public void populateResults(String fullPathToaltDataFile, int numThreads, boolean verbose, boolean svdRegression) {
		initTests(fullPathToaltDataFile);
		getProj().getLog().reportTimeInfo("Generating Manhattan plot(s) from " + markersToTest.length + " markers");
		MarkerDataLoader markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(getProj(), markersToTest);
		for (int i = 0; i < markersToTest.length; i++) {
			if (i != 0 && i % 50 == 0) {
				getProj().getLog().reportTimeInfo("marker " + i + " of " + markersToTest.length);
			}
			MarkerData markerData = markerDataLoader.requestMarkerData(i);
			markerDataLoader.releaseIndex(i);
			populateDataForMarker(numThreads, verbose, svdRegression, i, markerData);
		}
	}

	/**
	 * So that this can be done in conjuntion with other marker based things
	 */
	public void populateDataForMarker(int numThreads, boolean verbose, boolean svdRegression, int i, MarkerData markerData) {
		double[] lrrs = Array.toDoubleArray(markerData.getLRRs());
		ExecutorService executor = Executors.newFixedThreadPool(numThreads);
		Hashtable<String, Future<double[]>> tmpResults = new Hashtable<String, Future<double[]>>();
		for (int j = 0; j < manhattanTests.length; j++) {
			tmpResults.put(j + "", executor.submit(new ManhattanTestWorker(manhattanTests[j], lrrs, verbose, svdRegression, getProj().getLog())));
		}
		for (int j = 0; j < manhattanTests.length; j++) {
			try {
				results[j][i] = tmpResults.get(j + "").get();
			} catch (InterruptedException e) {
				getProj().getLog().reportError("Error - could running Picard on internal index " + j);
				getProj().getLog().reportException(e);
			} catch (ExecutionException e) {
				getProj().getLog().reportError("Error - could running Picard on internal index " + j);
				getProj().getLog().reportException(e);
			}
		}
		executor.shutdown();
		try {
			executor.awaitTermination(10, TimeUnit.DAYS);
		} catch (InterruptedException e) {
			getProj().getLog().reportException(e);
		}
	}

	/**
	 * dumps all tests to separate files
	 */
	public void dumpResults(String fullPathToOutputBase) {
		MarkerSet markerSet = proj.getMarkerSet();
		int[] markerIndicesInProject = ext.indexLargeFactors(markersToTest, markerSet.getMarkerNames(), true, getProj().getLog(), true, false);
		byte[] chrs = markerSet.getChrs();
		int[] pos = markerSet.getPositions();
		for (int i = 0; i < manhattanTests.length; i++) {
			String output = fullPathToOutputBase + ext.replaceWithLinuxSafeCharacters(manhattanTests[i].getTitle(), true) + EXT;

			try {
				PrintWriter writer = new PrintWriter(new FileWriter(output));
				writer.println(Array.toStr(HEADER));
				for (int j = 0; j < markersToTest.length; j++) {
					writer.println((j + 1) + "\t" + markersToTest[j] + "\t" + chrs[markerIndicesInProject[j]] + "\t" + pos[markerIndicesInProject[j]] + "\t" + (results[i][j] == null ? "NaN\tNaN" : results[i][j][0] + "\t" + Math.abs(results[i][j][1])));
				}
				writer.close();
			} catch (Exception e) {
				getProj().getLog().reportError("Error writing to " + output);
				getProj().getLog().reportException(e);
			}
		}

	}

	/**
	 * adds tests from file if not null
	 */
	private void initTests(String fullPathToaltDataFile) {
		manhattanTests = getManhattanTestsForPCs(this, numPCs);
		if (fullPathToaltDataFile != null) {
			proj.getLog().reportTimeInfo("Attempting to load " + fullPathToaltDataFile);
			ManhattanTest[] tmp = loadManhattanTestFromFile(getProj(), fullPathToaltDataFile);
			manhattanTests = Array.concatAll(tmp, manhattanTests);
		}
		updateResultSize();
	}

	private void updateResultSize() {
		results = new double[manhattanTests.length][markersToTest.length][];
	}

	private static ManhattanTest[] getManhattanTestsForPCs(PrincipalComponentsManhattan pcManhattan, int numPCs) {
		ManhattanTest[] manhattanTests = new ManhattanTest[numPCs];
		for (int i = 0; i < numPCs; i++) {
			manhattanTests[i] = new ManhattanTest("PC" + (i + 1), pcManhattan.getBasisAt((i + 1)), null);
		}
		return manhattanTests;
	}

	private static class ManhattanTestWorker implements Callable<double[]> {
		private ManhattanTest mTest;
		private double[] dataToTest;
		private boolean verbose;
		private boolean svdRegression;
		private Logger log;

		public ManhattanTestWorker(ManhattanTest mTest, double[] dataToTest, boolean verbose, boolean svdRegression, Logger log) {
			super();
			this.mTest = mTest;
			this.dataToTest = dataToTest;
			this.verbose = verbose;
			this.svdRegression = svdRegression;
			this.log = log;
		}

		@Override
		public double[] call() throws Exception {
			if (verbose) {
				log.reportTimeInfo("Executing on thread" + Thread.currentThread().getName());
			}
			return mTest.testData(dataToTest, verbose, svdRegression, log);
		}
	}

	public static class ManhattanTest {
		private static final String SAMPLE = "DNA";
		private String title;
		private double[] dataTest;
		private boolean[] dataMask;

		public ManhattanTest(String title, double[] dataTest, boolean[] dataMask) {
			super();
			this.title = title;
			this.dataTest = dataTest;
			this.dataMask = dataMask;
		}

		public double[] testData(double[] dataToTest, boolean verbose, boolean svdRegression, Logger log) {
			if (dataToTest.length != dataTest.length) {
				log.reportTimeError("Mismatched array sizes for regression input");
			}
			double[] stats = null;
			double[] tmpdep = dataTest;
			double[] tmpInd = dataToTest;
			if (dataMask != null) {
				tmpdep = Array.subArray(tmpdep, dataMask);
				tmpInd = Array.subArray(tmpInd, dataMask);
			}
			CrossValidation crossValidation = new CrossValidation(tmpdep, Matrix.toMatrix(tmpInd), tmpdep, Matrix.toMatrix(tmpInd), verbose, svdRegression, log);
			crossValidation.train();
			crossValidation.computePredictedValues();
			crossValidation.computeResiduals();
			if (!crossValidation.analysisFailed()) {
				stats = crossValidation.getStats();
			}
			return stats;
		}

		public String getTitle() {
			return title;
		}

		public double[] getDataTest() {
			return dataTest;
		}

	}

	/**
	 * @param proj
	 * @param fullPathToFile
	 *            , must have a DNA column corresponding to samples in the project (can be missing samples) , all other columns will be tested
	 * @return
	 */
	private static ManhattanTest[] loadManhattanTestFromFile(Project proj, String fullPathToFile) {
		ManhattanTest[] maTests = null;
		String[] header = Files.getHeaderOfFile(fullPathToFile, proj.getLog());
		int sampIndex = ext.indexOfStr(ManhattanTest.SAMPLE, header);
		if (sampIndex < 0) {
			proj.getLog().reportTimeError(fullPathToFile + " must have a header with " + ManhattanTest.SAMPLE);
			return null;
		} else {
			String[] titles = new String[header.length - 1];
			int titleIndex = 0;
			for (int i = 0; i < header.length; i++) {
				if (i != sampIndex) {
					titles[titleIndex] = header[i];
				}
			}
			int[] titleIndices = ext.indexFactors(titles, header, true, true);
			double[][] data = new double[titles.length][proj.getSamples().length];
			boolean[][] masks = new boolean[titles.length][proj.getSamples().length];
			for (int i = 0; i < masks.length; i++) {
				Arrays.fill(data[i], Double.NaN);
				Arrays.fill(masks[i], false);
			}
			try {
				BufferedReader reader = Files.getAppropriateReader(fullPathToFile);
				reader.readLine();
				while (reader.ready()) {
					String[] line = reader.readLine().trim().split("[\\s]+");
					int curSamp = ext.indexOfStr(line[sampIndex], proj.getSamples());
					if (curSamp < 0) {
						proj.getLog().reportTimeError("did not find sample " + line[sampIndex] + " in this project");
						return null;
					} else {
						for (int i = 0; i < titleIndices.length; i++) {
							data[i][curSamp] = Double.parseDouble(line[titleIndices[i]]);
							masks[i][curSamp] = true;
						}
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				proj.getLog().reportError("Error: file \"" + fullPathToFile + "\" not found in current directory");
			} catch (IOException ioe) {
				proj.getLog().reportError("Error reading file \"" + fullPathToFile + "\"");
			}
			maTests = new ManhattanTest[titleIndices.length];
			for (int i = 0; i < maTests.length; i++) {
				proj.getLog().reportTimeInfo("Found " + Array.booleanArraySum(masks[i]) + " samples for column " + titles[i]);
				maTests[i] = new ManhattanTest(titles[i], data[i], masks[i]);
			}
		}

		return maTests;
	}

	public static void createManhattans(Project proj, String outputBase, String altDataFile, int numPCs, int numThreads, boolean verbose, boolean svdRegression) {
		proj.getLog().reportTimeInfo("Loading markers from " + proj.getFilename(Project.TARGET_MARKERS_FILENAME));
		String[] markers = HashVec.loadFileToStringArray(proj.getFilename(Project.TARGET_MARKERS_FILENAME), false, new int[] { 0 }, true);
		PrincipalComponentsManhattan principalComponentsManhattan = new PrincipalComponentsManhattan(proj, markers, numPCs);
		principalComponentsManhattan.populateResults(altDataFile == null ? null : proj.getProjectDir() + altDataFile, numThreads, verbose, svdRegression);
		outputBase = proj.getProjectDir() + outputBase;
		new File(ext.parseDirectoryOfFile(outputBase)).mkdirs();
		principalComponentsManhattan.dumpResults(outputBase);

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String logfile = null;
		int numPCs = 5;
		String altDataFile = "PC1Test.txt";
		int numThreads = 8;
		String outputBase = "Manhattan/manhattan";

		String usage = "\n" + "cnv.analysis.pca.PrincipalComponentsManhattan requires 0-1 arguments\n";
		usage += "   (1) project filename (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) number of PCs to test (i.e. numPCs=" + numPCs + " (default))\n" + "";
		usage += "   (3) another data file (under the project's directory) to use in the tests, must have a header with at least this column " + ManhattanTest.SAMPLE + " (i.e. altDataFile=" + altDataFile + " (no default))\n" + "";
		usage += "   (4) number of threads (i.e. " + PSF.Ext.NUM_THREADS_COMMAND + numThreads + " (default))\n" + "";
		usage += "   (5) output base name (i.e.  outputBase=" + outputBase + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("altDataFile=")) {
				altDataFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("outputBase=")) {
				outputBase = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("numPCs=")) {
				numPCs = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(args[i]);
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
			createManhattans(proj, outputBase, altDataFile, numPCs, numThreads, false, false);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
