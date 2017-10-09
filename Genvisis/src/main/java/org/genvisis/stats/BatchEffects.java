package org.genvisis.stats;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import javax.annotation.Nonnull;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.stat.inference.TTest;
import org.genvisis.CLI;
import org.genvisis.cnv.plots.TwoDPlot;
import org.genvisis.cnv.plots.TwoDPlot.ScreenToCapture;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Images;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

/**
 * This class is designed to look for batch effects within sample factor data. It accomplishes this goal by running t-tests on sample factor data organized by batch, then visualizing the results.
 * 
 * @author Travis Rogers
 */
public class BatchEffects {
	private final List<String> batchLabels;
	private final List<String> factorLabels;
	private final List<String> batchValuesToInclude;
	private final List<List<Double>> factorValuesToInclude;
	private double maxNegLog10PValue;
	private final Logger logger;
	private double pValueTruncation = 1.0E-300;

	private BatchEffects(Logger logger) {
		batchLabels = new ArrayList<>();
		factorLabels = new ArrayList<>();
		batchValuesToInclude = new ArrayList<>();
		factorValuesToInclude = new ArrayList<>();
		this.logger = logger;
	}

	/**
	 * Returns BatchEffects instance containing parsed batch and factor data.
	 * see {@link #parseFiles(String,String)} for more information on how to format input files.
	 * 
	 * @param batchFilePath String path to file containing labels in first row, sample identifiers in first column, and batch values in second column.
	 * @param factorFilePath String path to file containing labels in first row, sample identifiers in first column, and factor values in subsequent columns.
	 * @param logger Logger used by instance.
	 * @return BatchEffects instance containing parsed batch and factor data.
	 * @throws FileNotFoundException if either String path parameter does not lead to a valid file.
	 * @throws ParseException if file with factor data contains fewer than 2 column labels.
	 * @throws IOException
	 */
	public static BatchEffects getInstance(@Nonnull String batchFilePath, @Nonnull String factorFilePath, @Nonnull Logger logger) throws FileNotFoundException, ParseException, IOException {
		BatchEffects instance = new BatchEffects(logger);
		instance.parseFiles(batchFilePath, factorFilePath);
		return instance;
	}

	/**
	 * Parses factor and batch data from files, then creates batch effect screenshots.
	 * see {@link #parseFiles(String,String)} for more information on how to format input files.
	 * 
	 * @param args <ul><li>first argument (required) is full path to file with batch data. <li>second argument (required) is full path to file with factor data. <li>third argument (optional) is lower truncation boundary for p-values. </ul>
	 * @throws FileNotFoundException if either String path parameter does not lead to a valid file.
	 * @throws ParseException instance containing parsed batch and factor data.
	 * @throws IOException
	 * @throws org.apache.commons.cli.ParseException 
	 */
	public static void main(String[] args) throws FileNotFoundException, ParseException, IOException, org.apache.commons.cli.ParseException {
		CLI cli = new CLI("BatchEffect-Plots");
		cli.addArgWithDefault("pValueTruncation", "Batch effect p-values less than this double value will be set to this value.", 1.0E-300);
		cli.addArg("batchFilePath", "String path to file containing labels in first row, sample identifiers in first column, and batch values in second column.", true, CLI.Arg.FILE);
		cli.addArg("factorFilePath", "String path to file containing labels in first row, sample identifiers in first column, and factor values in subsequent columns.", true, CLI.Arg.FILE);
		cli.parse(args);
		if (cli.has("batchFilePath") && cli.has("factorFilePath")) {
			BatchEffects instance = BatchEffects.getInstance(cli.get("batchFilePath"), cli.get("factorFilePath"), new Logger());
			instance.pValueTruncation = cli.getD("pValueTruncation");
			instance.createNegLog10PValueScatterPlotScreenshots();
		} 
	}

	/**
	 * Creates scatter plot images of negative Log10 p-values obtained by comparing factor values between batches.
	 * Saves the following to current working directory:
	 * <ul>
	 * <li>text file containing data used for images
	 * <li>individual scatter plot images
	 * <li>aggregate image containing individual images stitched together
	 * </ul>
	 */
	public void createNegLog10PValueScatterPlotScreenshots() {
		// create negative log10 p-value matrix from class-level batch and factor data
		String[][] negLog10PValueMatrix;
		negLog10PValueMatrix = this.getNegLog10PValueMatrix(pValueTruncation);
		String outputFileName = System.currentTimeMillis() + ".txt";

		// write matrix to file because twoDPlot.createScreenshots requires data to be in a file 
		logger.report("writing negative log10 p-value matrix to file...");
		Files.writeMatrix(negLog10PValueMatrix, outputFileName, "\t");
		logger.report("new text file: " + outputFileName);
		logger.report("writing to file complete");

		// set up and configure TwoDPlot and ScreenToCapture objects needed for image creation
		TwoDPlot twoDPlot = TwoDPlot.createGUI(null, false, false);
		twoDPlot.getPanel().setSize(800, 600); // this line determines the size of each individual image
		twoDPlot.getPanel().overrideAxisLabels("Factors", "-Log10(P-Value)");
		List<ScreenToCapture> screens = new ArrayList<>();
		String[] plotFileNames = new String[negLog10PValueMatrix[0].length - 1];
		for (int i = 1; i < negLog10PValueMatrix[0].length; i++) {
			String[] files = {outputFileName, outputFileName, null};
			int[] dataIndices = {0, i, -1};
			int[] idIndices = {0, 0, -1};
			float[] displayWindow = {0, negLog10PValueMatrix.length, 0, (float)(Math.ceil(this.maxNegLog10PValue))};
			String plotFileName = System.currentTimeMillis() + "_" + i;
			plotFileNames[i - 1] = plotFileName;
			String title = (i == negLog10PValueMatrix[0].length - 1)?"Min.P-Value":this.batchLabels.get(i - 1);
			ScreenToCapture screen = new ScreenToCapture(files, dataIndices, idIndices, displayWindow, false, false, false, false, plotFileName, title);
			screens.add(screen);
		}

		// create individual images
		String targetDirectory = System.getProperty("user.dir") + "/";
		logger.report("creating screenshots...");
		twoDPlot.createScreenshots(targetDirectory, screens);
		for (int i = 0; i < plotFileNames.length; i++) {
			plotFileNames[i] += ".png";
			logger.report("new image: " + plotFileNames[i]);
		}

		// create aggregate image
		String summaryImageFileName = System.currentTimeMillis() + "_Aggregate.png";
		Images.stitchImages(targetDirectory, plotFileNames, summaryImageFileName, null, false, false);
		logger.report("new aggregate image: " + summaryImageFileName);
		logger.report("batchEffects has completed running");
	}

	/**
	 * Returns a two-dimensional array of negative log10 p-values. 
	 * Each row represents a factor. Each column represents a batch. 
	 * The first row contains batch labels. The first column contains an arbitrary numerical index. 
	 * For a given factor, each p-value represents the result of a t-test of the values for one batch compared with all other batches combined.
	 * 
	 * @param pValueTruncation p-values less than this double value will be set to this value.
	 * @return two-dimensional array of negative log10 p-values.
	 */
	public String[][] getNegLog10PValueMatrix(double pValueTruncation) {
		// create boolean matrix representing sample membership in each of the batches detected. each sub-array index of the same position represents the same sample.
		logger.report("generating p-value matrix...");
		boolean[][] batchMembershipMatrix = ArrayUtils.classifyStringsToBoolean(batchValuesToInclude.toArray(new String[] {}), new String[0]).getClassified();

		// instantiate negative log10 p-value matrix and add column headers
		String[][] negLog10PValueMatrix = new String[factorLabels.size() + 1][];
		negLog10PValueMatrix[0] = new String[batchLabels.size() + 2];
		negLog10PValueMatrix[0][0] = "Factors";
		negLog10PValueMatrix[0][negLog10PValueMatrix[0].length - 1] = "Min.P-Value";
		for (int i = 1; i < negLog10PValueMatrix[0].length - 1; i++) {
			negLog10PValueMatrix[0][i] = batchLabels.get(i - 1);
		}

		maxNegLog10PValue = Double.MIN_VALUE; // this is updated and used to determine the maximum y-axis value on scatter plot images

		// run t-tests and populate matrix that will be returned
		// loop for each factor
		for (int i = 0; i < factorValuesToInclude.size(); i++) {
			negLog10PValueMatrix[i + 1] = new String[batchLabels.size() + 2]; // 2 is added to the number of batches because index 0 will contain an arbitrary row index and the last index will contain the minimum batch effects p-value detected for that factor
			negLog10PValueMatrix[i + 1][0] = String.valueOf(i + 1);  // this line sets arbitrary indices as row identifiers. these values are used for x-axis values on scatter plot images
			double minPValue = Double.MAX_VALUE;
			// loop for each batch
			for (int j = 0; j < batchMembershipMatrix.length; j++) {
				double[][] tTestGroups = createTTestGroups(batchMembershipMatrix[j], factorValuesToInclude.get(i));
				TTest tTest = new TTest();
				// run t-test for current factor values of current batch compared with current factor values of all other batches combined
				double pValue = tTest.tTest(tTestGroups[0], tTestGroups[1]);
				if (pValue < minPValue) {
					minPValue = pValue;
				}
				if (pValue < pValueTruncation) {
					pValue = pValueTruncation;
				}
				double negLog10PValue = -1.0 * Math.log10(pValue);
				negLog10PValueMatrix[i + 1][j + 1] = String.valueOf(negLog10PValue); // add negative log10 p-value to matrix
				if (negLog10PValue > maxNegLog10PValue) {
					maxNegLog10PValue = negLog10PValue;
				}
			}
			negLog10PValueMatrix[i + 1][negLog10PValueMatrix[i].length - 1] = String.valueOf(-1.0 * Math.log10(minPValue)); // add minimum batch effect p-value for current factor to last index of array
		}
		logger.report("p-value matrix complete");

		return negLog10PValueMatrix;
	}

	/**
	 * Parses sample batch and factor information from files to class-level data structures.
	 * 
	 * @param batchFilePath String path to file containing labels in first row, sample identifiers in first column, and batch values in second column.
	 * @param factorFilePath String path to file containing labels in first row, sample identifiers in first column, and factor values in subsequent columns.
	 * @throws FileNotFoundException if either String path parameter does not lead to a valid file.
	 * @throws ParseException if file with factor data contains fewer than 2 column labels.
	 * @throws IOException
	 */
	private void parseFiles(@Nonnull String batchFilePath, @Nonnull String factorFilePath) throws FileNotFoundException, ParseException, IOException {
		logger.report("parsing data from files...");

		// parse batch file
		Map<String, String> sampleBatchMap = HashVec.loadFileToHashString(batchFilePath, new int[] {0}, new int[] {1}, false, null, true, false);

		// parse factor file
		BufferedReader reader = Files.getAppropriateReader(factorFilePath);
		String[] rowData;
		int expectedNumColumns = -1;
		// read factor column headers
		String line = reader.readLine();
		if (line != null) {
			rowData = line.split("\t");
			if (rowData.length >= 2) {
				expectedNumColumns = rowData.length;
				for (int i = 1; i < rowData.length; i++) { // this loop starts at 1 to exclude the sample id column
					factorLabels.add(rowData[i]);
					factorValuesToInclude.add(new ArrayList<>()); // each sub-list represents a different factor
				}
			} else {
				throw new ParseException("Error: expected at least 2 column headers in factor data. Found " + rowData.length + ".", 0);
			}
		}
		// read factor sample data
		long numSampleLinesFromFactorData = 0;
		long numSamplesIncluded = 0;
		while ((line = reader.readLine()) != null) {
			numSampleLinesFromFactorData++;
			rowData = line.split("\t");
			boolean sampleIncluded = addSampleIfValid(rowData, sampleBatchMap, expectedNumColumns); // add current sample's batch and factor data to class level data structures if passes validation
			if (sampleIncluded) { 
				numSamplesIncluded++; 
			} else {
				logger.report("line excluded from factor file: " + line);
			}
		}
		logger.report("parsing complete");
		logger.report(numSamplesIncluded + " samples will be included in analysis. " + (numSampleLinesFromFactorData - numSamplesIncluded) + " lines of potential sample data were excluded from the factor file. " 
				+ (sampleBatchMap.size() - numSamplesIncluded) + " lines of potential sample data were excluded from the batch file. ");
	}

	/**
	 * If sample passes validation, adds sample values to class-level data structures that will be used for analysis.
	 * Sample will fail validation if it is missing any data or if any of its factor values cannot be parsed as a double.
	 * 
	 * @param factorDataForSample array of factor data. index 0 should be the sample identifier. subsequent indices should be factor values for that sample.
	 * @param sampleBatchMap map of sample identifier to batch value.
	 * @param expectedLength expected length of factorDataForSample parameter.
	 * @return a boolean value indicating whether the sample was added to class-level data structures that will be used for analysis.
	 */
	private boolean addSampleIfValid(String[] factorDataForSample, Map<String, String> sampleBatchMap, int expectedLength) {
		boolean includeSample = true;

		if (factorDataForSample.length != expectedLength) { 
			includeSample = false; 
			}

		// exclude sample if the sample id from the factor file does not match any from the batch file OR if the batch value associated with that sample id is missing 
		if (includeSample) {
			String batchValue = sampleBatchMap.get(factorDataForSample[0]);
			if (batchValue == null || ext.isMissingValue(batchValue)) { includeSample = false; }
		}

		// exclude sample if any of its factor values are missing or do not parse into doubles
		List<Double> parsedValues = new ArrayList<>();
		if (includeSample) {
			for (int i = 0; includeSample && i < factorDataForSample.length; i++) {
				if (ext.isMissingValue(factorDataForSample[i])) { includeSample = false; }
				if (includeSample && i != 0) {
					try {
						double d = Double.parseDouble(factorDataForSample[i]);
						parsedValues.add(d);
					} catch (NumberFormatException e) {
						includeSample = false;
					}
				}
			}
		}

		// if sample passes validation, add to class-level data structures that will be used for analysis
		if (includeSample) {
			String batchValue = sampleBatchMap.get(factorDataForSample[0]); // inputValues[0] is the sample id
			if (!batchLabels.contains(batchValue)) {
				batchLabels.add(batchValue);
			}
			batchValuesToInclude.add(batchValue); 
			for (int i = 0; i < parsedValues.size(); i++) {
				double currentValue = parsedValues.get(i);
				factorValuesToInclude.get(i).add(currentValue);
			}
		}

		// return boolean indication of whether sample was included
		return includeSample;
	}

	/**
	 * Splits the values from sampleValues into two arrays, based on the true/false indices found in sampleMembership. 
	 * 
	 * sampleMembership and sampleValues should have parallel indices, each index representing a specific sample. 
	 * if sampleMembership and sampleValues arguments are both size 0, then a two-dimensional array will be returned that contains two arrays of length 0.
	 * 
	 * @param sampleMembership boolean array that conveys sample membership in one of two groups (true or false). parallel sample indices with sampleValues. used to assign values from sampleValues to one of two arrays. not null.
	 * @param sampleValues list of values that will be split into two groups for T-Test use. not null.
	 * @return a two dimensional array that contains only two arrays, representing the two T-Test groups.
	 * @throws DimensionMismatchException if sampleMembership length and sampleValues length are not equal.
	 */
	private double[][] createTTestGroups(@Nonnull boolean[] sampleMembership, @Nonnull List<Double> sampleValues) throws DimensionMismatchException {
		//incoming data structures should be same size because each parallel index represents the same sample
		if (sampleMembership.length != sampleValues.size()) {
			throw new DimensionMismatchException(sampleMembership.length, sampleValues.size());
		}

		// assign sample values to two groups, based on boolean sample membership
		List<Double> group1 = new LinkedList<>();
		List<Double> group2 = new LinkedList<>();
		for (int i = 0; i < sampleValues.size(); i++) {
			if (sampleMembership[i] == true) {
				group1.add(sampleValues.get(i));
			} else {
				group2.add(sampleValues.get(i));
			}
		}

		// convert two lists into a single two-dimensional array
		double[][] tTestGroups = new double[2][];
		tTestGroups[0] = toPrimitive(group1);
		tTestGroups[1] = toPrimitive(group2);
		return tTestGroups;
	}

	/**
	 * Returns an array of primitive doubles in the same order as the list received.
	 * 
	 * @param list the list to convert to an array. not null. empty list returns empty array.
	 * @return an array containing the elements from the list argument.
	 */
	private double[] toPrimitive(@Nonnull List<Double> list) {
		double[] array = new double[list.size()];
		for (int i = 0; i < list.size(); i++) {
			array[i] = list.get(i);
		}
		return array;
	}

	/**
	 * Returns a view on batchLabels list.
	 */
	public List<String> getBatchLabels() {
		return Collections.unmodifiableList(batchLabels);
	}

	/**
	 * Returns a view on factorLabels list.
	 */
	public List<String> getFactorLabels() {
		return Collections.unmodifiableList(factorLabels);
	}

	/**
	 * Returns a view on batchValuesToInclude.
	 */
	public List<String> getBatchValuesToInclude() {
		return Collections.unmodifiableList(batchValuesToInclude);
	}

	/**
	 * Returns a view on factorValuesToInclude.
	 */
	public List<List<Double>> getFactorValuesToInclude() {
		return Collections.unmodifiableList(factorValuesToInclude);
	}

}
