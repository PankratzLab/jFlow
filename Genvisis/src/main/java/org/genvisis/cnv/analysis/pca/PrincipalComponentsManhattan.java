package org.genvisis.cnv.analysis.pca;

import java.awt.Component;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
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

import javax.swing.JOptionPane;

import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.MarkerDataLoader;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.stats.CrossValidation;
import org.genvisis.stats.LeastSquares.LS_TYPE;

/**
 * Generates manhattan plots for PCs and other optional data...can be loaded directly into IGV or
 * haploview
 *
 */
public class PrincipalComponentsManhattan extends PrincipalComponentsResiduals {
	/**
	 *
	 */
	private static final long serialVersionUID = 1L;
	/**
	 *
	 */
	public static final String PRINCIPAL_MANHATTAN_MI = "Generate manhattan plots";
	private static final String[] HEADER = {"Index", "SNP", "CHR", "BP", "P", "T_STAT_ABS"};
	private static final String EXT = ".hat.linear";

	private final String[] markersToTest;
	private double[][][] results;// type,marker,result
	private ManhattanTest[] manhattanTests;
	private final int numPCs;

	/**
	 * @param proj
	 * @param markersToTest these markers will be tested against the pcs, and any other data provided
	 * @param fullPathToaltDataFile besides manhattans for the pcs, generate for this data file as
	 *        well
	 * @param numPcs
	 */
	public PrincipalComponentsManhattan(Project proj, String[] markersToTest,
																			String fullPathToaltDataFile, int numPcs) {
		super(proj.loadPcResids());
		this.markersToTest = PrincipalComponentsCompute.sortByProjectMarkers(getProj(), markersToTest);
		numPCs = numPcs > 0 ? numPcs : getNumComponents();
		initTests(fullPathToaltDataFile);

	}

	/**
	 * @param mTests can add more tests here,just call before populating results
	 */
	public void addManhattanTest(ManhattanTest[] mTests) {
		manhattanTests = ArrayUtils.concatAll(manhattanTests, mTests);
		updateResultSize();
	}

	/**
	 * computes for all markers
	 */
	public void populateResults(int numThreads, boolean verbose, LS_TYPE lType) {
		getProj().getLog().reportTimeInfo("Generating Manhattan plot(s) from " + markersToTest.length
																			+ " markers");
		MarkerDataLoader markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(getProj(),
																																																markersToTest);
		for (int i = 0; i < markersToTest.length; i++) {
			if (i != 0 && i % 50 == 0) {
				getProj().getLog().reportTimeInfo("marker " + i + " of " + markersToTest.length);
			}
			MarkerData markerData = markerDataLoader.requestMarkerData(i);
			markerDataLoader.releaseIndex(i);
			populateDataForMarker(numThreads, verbose, lType, i, markerData);
		}
	}

	/**
	 * So that this can be done in conjunction with other marker based things
	 */
	public void populateDataForMarker(int numThreads, boolean verbose, LS_TYPE lType, int i,
																		MarkerData markerData) {
		double[] lrrs = ArrayUtils.toDoubleArray(markerData.getLRRs());
		ExecutorService executor = Executors.newFixedThreadPool(numThreads);
		Hashtable<String, Future<double[]>> tmpResults = new Hashtable<String, Future<double[]>>();
		for (int j = 0; j < manhattanTests.length; j++) {
			tmpResults.put(j + "",
										 executor.submit(new ManhattanTestWorker(manhattanTests[j], lrrs, verbose,
																														 lType, getProj().getLog())));
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
		MarkerSetInfo markerSet = proj.getMarkerSet();
		int[] markerIndicesInProject = ext.indexLargeFactors(markersToTest, markerSet.getMarkerNames(),
																												 true, getProj().getLog(), true);
		byte[] chrs = markerSet.getChrs();
		int[] pos = markerSet.getPositions();
		for (int i = 0; i < manhattanTests.length; i++) {
			String output = fullPathToOutputBase
											+ ext.replaceWithLinuxSafeCharacters(manhattanTests[i].getTitle(), true)
											+ EXT;

			try {
				PrintWriter writer = Files.openAppropriateWriter(output);
				writer.println(ArrayUtils.toStr(HEADER));
				for (int j = 0; j < markersToTest.length; j++) {
					writer.println((j + 1) + "\t" + markersToTest[j] + "\t" + chrs[markerIndicesInProject[j]]
												 + "\t" + pos[markerIndicesInProject[j]] + "\t"
												 + (results[i][j] == null ? "NaN\tNaN" : results[i][j][0] + "\t"
																																 + Math.abs(results[i][j][1])));
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
			manhattanTests = ArrayUtils.concatAll(tmp, manhattanTests);
		}
		updateResultSize();
	}

	private void updateResultSize() {
		results = new double[manhattanTests.length][markersToTest.length][];
	}

	private static ManhattanTest[] getManhattanTestsForPCs(PrincipalComponentsManhattan pcManhattan,
																												 int numPCs) {
		ManhattanTest[] manhattanTests = new ManhattanTest[numPCs];
		for (int i = 0; i < numPCs; i++) {
			manhattanTests[i] = new ManhattanTest("PC" + (i + 1), pcManhattan.getBasisAt((i + 1)), null);
		}
		return manhattanTests;
	}

	private static class ManhattanTestWorker implements Callable<double[]> {
		private final ManhattanTest mTest;
		private final double[] dataToTest;
		private final boolean verbose;
		private final LS_TYPE lType;
		private final Logger log;

		public ManhattanTestWorker(ManhattanTest mTest, double[] dataToTest, boolean verbose,
															 LS_TYPE lType, Logger log) {
			super();
			this.mTest = mTest;
			this.dataToTest = dataToTest;
			this.verbose = verbose;
			this.lType = lType;
			this.log = log;
		}

		@Override
		public double[] call() throws Exception {
			if (verbose) {
				log.reportTimeInfo("Executing on thread" + Thread.currentThread().getName());
			}
			return mTest.testData(dataToTest, verbose, lType, log);
		}
	}

	public static class ManhattanTest {
		private static final String SAMPLE = "DNA";
		private final String title;
		private final double[] dataTest;
		private final boolean[] dataMask;

		public ManhattanTest(String title, double[] dataTest, boolean[] dataMask) {
			super();
			this.title = title;
			this.dataTest = dataTest;
			this.dataMask = dataMask;
		}

		public double[] testData(double[] dataToTest, boolean verbose, LS_TYPE lType, Logger log) {
			if (dataToTest.length != dataTest.length) {
				log.reportError("Mismatched array sizes for regression input");
			}
			double[] result = new double[2];
			Arrays.fill(result, Double.NaN);
			double[] tmpdep = dataTest;
			double[] tmpInd = dataToTest;
			if (dataMask != null) {
				tmpdep = ArrayUtils.subArray(tmpdep, dataMask);
				tmpInd = ArrayUtils.subArray(tmpInd, dataMask);
			}
			CrossValidation crossValidation = new CrossValidation(tmpdep, Matrix.toMatrix(tmpInd), tmpdep,
																														Matrix.toMatrix(tmpInd), verbose, lType,
																														log);
			crossValidation.train();
			crossValidation.computePredictedValues();
			crossValidation.computeResiduals();
			if (!crossValidation.analysisFailed()) {
				result[0] = crossValidation.getSigs()[1];
				result[1] = crossValidation.getStats()[1];
			}
			return result;
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
	 * @param fullPathToFile , must have a DNA column corresponding to samples in the project (can be
	 *        missing samples) , all other columns will be tested
	 * @return
	 */
	private static ManhattanTest[] loadManhattanTestFromFile(Project proj, String fullPathToFile) {
		ManhattanTest[] maTests = null;
		String[] header = Files.getHeaderOfFile(fullPathToFile, proj.getLog());
		int sampIndex = ext.indexOfStr(ManhattanTest.SAMPLE, header);
		if (sampIndex < 0) {
			proj.getLog()
					.reportError(fullPathToFile + " must have a header with " + ManhattanTest.SAMPLE);
			return null;
		} else {
			String[] titles = new String[header.length - 1];
			int titleIndex = 0;
			for (int i = 0; i < header.length; i++) {
				if (i != sampIndex) {
					titles[titleIndex] = header[i];
				}
			}
			int[] titleIndices = ext.indexFactors(titles, header, true);
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
					String[] line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
					int curSamp = ext.indexOfStr(line[sampIndex], proj.getSamples());
					if (curSamp < 0) {
						proj.getLog()
								.reportError("did not find sample " + line[sampIndex] + " in this project");
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
				proj.getLog()
						.reportError("Error: file \"" + fullPathToFile + "\" not found in current directory");
			} catch (IOException ioe) {
				proj.getLog().reportError("Error reading file \"" + fullPathToFile + "\"");
			}
			maTests = new ManhattanTest[titleIndices.length];
			for (int i = 0; i < maTests.length; i++) {
				proj.getLog().reportTimeInfo("Found " + ArrayUtils.booleanArraySum(masks[i])
																		 + " samples for column " + titles[i]);
				maTests[i] = new ManhattanTest(titles[i], data[i], masks[i]);
			}
		}
		return maTests;
	}

	public static void guiAccess(Project proj, Component parentComponent) {
		String pcFile = proj.INTENSITY_PC_FILENAME.getValue(false, false);
		String[] targetMarkers = proj.getTargetMarkers();
		if (targetMarkers == null) {
			// JOptionPane.showMessageDialog(parentComponent, "Failed to load target markers '" +
			// proj.getFilename(proj.TARGET_MARKERS_FILENAME) + "'; this is the designated marker in the
			// project properties file", "Error", JOptionPane.ERROR_MESSAGE);
			JOptionPane.showMessageDialog(parentComponent,
																		"Failed to load target markers '"
																										 + proj.TARGET_MARKERS_FILENAMES.getValue()[0]
																										 + "'; this is the designated marker in the project properties file",
																		"Error", JOptionPane.ERROR_MESSAGE);
			return;
		}
		int numMarkers = targetMarkers.length;
		if (Files.exists(pcFile)) {
			String ObjButtons[] = {"OK", "Cancel"};
			// int promptResult = JOptionPane.showOptionDialog(parentComponent, "Generate manhattan plots
			// for " + numMarkers + " marker(s) over " + proj.getInt(proj.INTENSITY_PC_NUM_COMPONENTS) + "
			// component(s)?", "Manhattan Plot", JOptionPane.DEFAULT_OPTION, JOptionPane.WARNING_MESSAGE,
			// null, ObjButtons, ObjButtons[1]);
			int promptResult = JOptionPane.showOptionDialog(parentComponent,
																											"Generate manhattan plots for " + numMarkers
																																			 + " marker(s) over "
																																			 + proj.INTENSITY_PC_NUM_COMPONENTS.getValue()
																																			 + " component(s)?",
																											"Manhattan Plot", JOptionPane.DEFAULT_OPTION,
																											JOptionPane.WARNING_MESSAGE, null, ObjButtons,
																											ObjButtons[1]);
			if (promptResult == 0) {
				PrincipalComponentsManhattan.createManhattans(proj);
			}
		} else {
			JOptionPane.showMessageDialog(parentComponent,
																		"Failed to detect " + proj.INTENSITY_PC_FILENAME + " " + pcFile
																										 + " ; this is the designated intensity pc filename in the project properties file",
																		"Error", JOptionPane.ERROR_MESSAGE);
		}
	}

	/**
	 * Mainly for gui access using all defualts
	 */
	public static void createManhattans(Project proj) {
		// int numComponents = proj.getInt(proj.INTENSITY_PC_NUM_COMPONENTS);
		int numComponents = proj.INTENSITY_PC_NUM_COMPONENTS.getValue();
		createManhattans(proj, "Manhattan/manhattan", null, numComponents, 1, false,
										 numComponents >= 250 ? LS_TYPE.SVD : LS_TYPE.REGULAR);
	}

	public static void createManhattans(Project proj, String outputBase, String altDataFile,
																			int numPCs, int numThreads, boolean verbose, LS_TYPE ltType) {
		// proj.getLog().reportTimeInfo("Loading markers from " +
		// proj.getFilename(proj.TARGET_MARKERS_FILENAME));
		// String[] markers =
		// HashVec.loadFileToStringArray(proj.getFilename(proj.TARGET_MARKERS_FILENAME), false, new
		// int[] { 0 }, true);
		proj.getLog()
				.reportTimeInfo("Loading markers from " + proj.TARGET_MARKERS_FILENAMES.getValue()[0]);
		String[] markers = HashVec.loadFileToStringArray(proj.TARGET_MARKERS_FILENAMES.getValue()[0],
																										 false, new int[] {0}, true);
		PrincipalComponentsManhattan principalComponentsManhattan = new PrincipalComponentsManhattan(proj,
																																																 markers,
																																																 altDataFile == null ? null
																																																										 : proj.PROJECT_DIRECTORY.getValue()
																																																											 + altDataFile,
																																																 numPCs);
		principalComponentsManhattan.populateResults(numThreads, verbose, ltType);
		outputBase = proj.PROJECT_DIRECTORY.getValue() + outputBase;
		new File(ext.parseDirectoryOfFile(outputBase)).mkdirs();
		principalComponentsManhattan.dumpResults(outputBase);

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String logfile = null;
		int numPCs = 5;
		String altDataFile = "PC1Test.txt";
		int numThreads = 1;
		String outputBase = "Manhattan/manhattan";

		String usage = "\n" + "cnv.analysis.pca.PrincipalComponentsManhattan requires 0-1 arguments\n";
		usage += "   (1) project filename (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) number of PCs to test (i.e. numPCs=" + numPCs + " (default))\n" + "";
		usage += "   (3) another data file (under the project's directory) to use in the tests, must have a header with at least this column "
						 + ManhattanTest.SAMPLE + " (i.e. altDataFile=" + altDataFile + " (no default))\n" + "";
		usage += "   (4) number of threads (i.e. " + PSF.Ext.NUM_THREADS_COMMAND + numThreads
						 + " (default))\n" + "";
		usage += "   (5) output base name (i.e.  outputBase=" + outputBase + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("proj=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("altDataFile=")) {
				altDataFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("outputBase=")) {
				outputBase = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("numPCs=")) {
				numPCs = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(arg);
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
			Project proj = new Project(filename, logfile);
			createManhattans(proj, outputBase, altDataFile, numPCs, numThreads, false, LS_TYPE.REGULAR);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
