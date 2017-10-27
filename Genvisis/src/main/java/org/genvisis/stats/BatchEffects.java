package org.genvisis.stats;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Paths;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import javax.annotation.Nonnull;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.stat.inference.TTest;
import org.genvisis.CLI;
import org.genvisis.cnv.plots.TwoDPlot;
import org.genvisis.cnv.plots.TwoDPlot.ScreenToCapture;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Images;
import org.genvisis.common.Logger;

import com.google.common.collect.ImmutableList;

/**
 * Looks for batch effects within sample factor data. Accomplished by running t-tests on sample
 * factor data organized by batch, then visualizing the results.
 * 
 * @author Travis Rogers
 */
public class BatchEffects {
	private final List<String> batchLabels;
	private final List<String> factorLabels;
	/*
	 * the batchValuesToInclude list and factorValuesToInclude sub-lists have parallel indices, each
	 * index representing a value (batch value or factor value, respectively) for a particular sample
	 */
	private final List<String> batchValuesToInclude;
	// each list in factorValuesToInclude represents the sample values for a particular factor
	private final List<List<Double>> factorValuesToInclude;
	private double maxNegLog10PValue;
	private final Logger logger;
	private static final double DEFAULT_PVALUE_TRUNCTATION = 1.0E-300;



	BatchEffects(BatchEffectsBuilder b) {
		this.batchLabels = ImmutableList.copyOf(b.getBatchLabels());
		this.factorLabels = ImmutableList.copyOf(b.getFactorLabels());
		this.batchValuesToInclude = ImmutableList.copyOf(b.getBatchValuesToInclude());
		this.factorValuesToInclude = ImmutableList.copyOf(b.getFactorValuesToInclude());
		this.logger = b.getLogger();
	}

	/**
	 * Parses factor and batch data from files, then creates batch effect screenshots. See
	 * {@link #parseFiles(String,String)} for more information on how to format input files. Command
	 * line arguments use unordered name=value format.
	 * 
	 * @param args unordered arguments that follow name=value format
	 *        <ul>
	 *        <li>"batchFilePath": (required) String path to file containing labels in first row,
	 *        sample identifiers in first column, and batch values in second column
	 *        <li>"factorFilePath": (required) String path to file containing labels in first row,
	 *        sample identifiers in first column, and factor values in subsequent columns
	 *        <li>"pValueTruncation": (optional) batch effect p-values less than this double value
	 *        will be set to this value. default value is 1.0E-300
	 *        </ul>
	 * @throws FileNotFoundException if either String path parameter does not lead to a valid file.
	 * @throws ParseException instance containing parsed batch and factor data.
	 * @throws IOException
	 * @throws org.apache.commons.cli.ParseException
	 */
	public static void main(String[] args) throws FileNotFoundException, ParseException, IOException,
																				 org.apache.commons.cli.ParseException {
		CLI cli = new CLI("BatchEffect-Plots");
		cli.addArgWithDefault("pValueTruncation",
													"Batch effect p-values less than this double value will be set to this value.",
													DEFAULT_PVALUE_TRUNCTATION);
		cli.addArg("batchFilePath",
							 "Path to file containing labels in first row, sample identifiers in first column, and batch values in second column.",
							 true);
		cli.addArg("factorFilePath",
							 "Path to file containing labels in first row, sample identifiers in first column, and factor values in subsequent columns.",
							 true);
		cli.parse(args);
		if (cli.has("batchFilePath") && cli.has("factorFilePath")) {
			BatchEffectsBuilder builder = new BatchEffectsBuilder(new Logger());
			BatchEffects instance = builder.build(cli.get("batchFilePath"),
																						cli.get("factorFilePath"));
			instance.createNegLog10PValueScatterPlotScreenshots(cli.getD("pValueTruncation"));
		}
	}

	/**
	 * Calls {@link BatchEffects#createNegLog10PValueScatterPlotScreenshots(double)} with default
	 * p-value truncation.
	 */
	public void createNegLog10PValueScatterPlotScreenshots() {
		createNegLog10PValueScatterPlotScreenshots(DEFAULT_PVALUE_TRUNCTATION);
	}

	/**
	 * Creates scatter plot images of negative Log10 p-values obtained by comparing factor values
	 * between batches. Saves the following to current working directory:
	 * <ul>
	 * <li>text file containing data used for images
	 * <li>individual scatter plot images
	 * <li>aggregate image containing individual images stitched together
	 * </ul>
	 * 
	 * @param pValueTruncation p-values less than this double value will be set to this value
	 */
	public void createNegLog10PValueScatterPlotScreenshots(double pValueTruncation) {
		// create negative log10 p-value matrix from class-level batch and factor data
		String[][] negLog10PValueMatrix = getNegLog10PValueMatrix(pValueTruncation);
		String outputFileName = System.currentTimeMillis() + "_" + pValueTruncation + ".txt";

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
			float[] displayWindow = {0, negLog10PValueMatrix.length, 0,
															 (float) (Math.ceil(this.maxNegLog10PValue))};
			String plotFileName = System.currentTimeMillis() + "_" + i;
			plotFileNames[i - 1] = plotFileName;
			String title = (i == negLog10PValueMatrix[0].length - 1) ? "Min.P-Value"
																															 : this.batchLabels.get(i - 1);
			ScreenToCapture screen = new ScreenToCapture(files, dataIndices, idIndices, displayWindow,
																									 false, false, false, false, plotFileName, title);
			screens.add(screen);
		}

		// create individual images
		String targetDirectory = System.getProperty("user.dir") + "/";
		logger.report("creating screenshots...");
		twoDPlot.createScreenshots(targetDirectory, screens);
		for (int i = 0; i < plotFileNames.length; i++) {
			plotFileNames[i] += ".png";
		}

		// create aggregate image
		String outFile = System.currentTimeMillis() + "_" + pValueTruncation + ".png";
		Images.stitchImages(targetDirectory, plotFileNames, outFile, null, false, false);
		logger.report("new aggregate image: " + outFile);

		// delete individual image files
		for (String fileName : plotFileNames) {
			try {
				java.nio.file.Files.deleteIfExists(Paths.get(targetDirectory, fileName));
			} catch (IOException e) {
				logger.reportError("Error during clean up - " + targetDirectory + fileName
													 + " could not be deleted.");
			}
		}

		logger.report("batchEffects has completed running");
	}

	/**
	 * Calls {@link BatchEffects#getNegLog10PValueMatrix(double)} with default p-value truncation.
	 * 
	 * @return two-dimensional array of negative log10 p-values
	 */
	public String[][] getNegLog10PValueMatrix() {
		return getNegLog10PValueMatrix(DEFAULT_PVALUE_TRUNCTATION);
	}

	/**
	 * Returns a two-dimensional array of negative log10 p-values. Each row represents a factor. Each
	 * column represents a batch. The first row contains batch labels. The first column contains an
	 * arbitrary numerical index. For a given factor, each p-value represents the result of a t-test
	 * of the values for one batch compared with all other batches combined.
	 * 
	 * @param pValueTruncation p-values less than this double value will be set to this value
	 * @return two-dimensional array of negative log10 p-values
	 */
	public String[][] getNegLog10PValueMatrix(double pValueTruncation) {
		// create boolean matrix representing sample membership in each of the batches detected. each
		// sub-array index of the same position represents the same sample.
		logger.report("generating p-value matrix...");
		boolean[][] batchMembershipMatrix = ArrayUtils.classifyStringsToBoolean(batchValuesToInclude.toArray(new String[] {}),
																																						new String[0])
																									.getClassified();

		// instantiate negative log10 p-value matrix and add column headers
		String[][] negLog10PValueMatrix = new String[factorLabels.size() + 1][];
		negLog10PValueMatrix[0] = new String[batchLabels.size() + 2];
		negLog10PValueMatrix[0][0] = "Factors";
		negLog10PValueMatrix[0][negLog10PValueMatrix[0].length - 1] = "Min.P-Value";
		for (int i = 1; i < negLog10PValueMatrix[0].length - 1; i++) {
			negLog10PValueMatrix[0][i] = batchLabels.get(i - 1);
		}

		// this variable will be updated and used for the maximum y-axis value on scatter plot
		maxNegLog10PValue = Double.MIN_VALUE;

		// run t-tests and populate matrix that will be returned
		// loop for each factor
		for (int i = 0; i < factorValuesToInclude.size(); i++) {
			/*
			 * on the following line, 2 is added to the number of batches because two extra columns are
			 * desired in the matrix. the first column will hold factor identifiers and the last will hold
			 * the min p-value for each factor
			 */
			negLog10PValueMatrix[i + 1] = new String[batchLabels.size() + 2];
			/*
			 * the following line sets arbitrary indices as factor identifiers. if creating a scatter plot
			 * image, these values will be used on the x-axis
			 */
			negLog10PValueMatrix[i + 1][0] = String.valueOf(i + 1);
			double minPValue = Double.MAX_VALUE;
			// loop for each batch
			TTest tTest = new TTest();
			for (int j = 0; j < batchMembershipMatrix.length; j++) {
				double[][] tTestGroups = createTTestGroups(batchMembershipMatrix[j],
																									 factorValuesToInclude.get(i));
				// run t-test for current factor values of current batch compared with current factor values
				// of all other batches combined
				double pValue = tTest.tTest(tTestGroups[0], tTestGroups[1]);
				if (pValue < minPValue) {
					minPValue = pValue;
				}
				if (pValue < pValueTruncation) {
					pValue = pValueTruncation;
				}
				double negLog10PValue = -1.0 * Math.log10(pValue);
				negLog10PValueMatrix[i + 1][j + 1] = String.valueOf(negLog10PValue);
				if (negLog10PValue > maxNegLog10PValue) {
					maxNegLog10PValue = negLog10PValue;
				}
			}
			// include minimum p-values to last column
			negLog10PValueMatrix[i + 1][negLog10PValueMatrix[i].length
																	- 1] = String.valueOf(-1.0 * Math.log10(minPValue));
		}
		logger.report("p-value matrix complete");

		return negLog10PValueMatrix;
	}

	/**
	 * Splits the values from sampleValues into two arrays, based on the true/false indices found in
	 * sampleMembership.
	 * 
	 * sampleMembership and sampleValues should have parallel indices, each index representing a
	 * specific sample. if sampleMembership and sampleValues arguments are both size 0, then a
	 * two-dimensional array will be returned that contains two arrays of length 0.
	 * 
	 * @param sampleMembership boolean array that conveys sample membership in one of two groups (true
	 *        or false). parallel sample indices with sampleValues. used to assign values from
	 *        sampleValues to one of two arrays. not null.
	 * @param sampleValues list of values that will be split into two groups for T-Test use. not null.
	 * @return a two dimensional array that contains only two arrays, representing the two T-Test
	 *         groups.
	 * @throws DimensionMismatchException if sampleMembership length and sampleValues length are not
	 *         equal.
	 */
	private double[][] createTTestGroups(@Nonnull boolean[] sampleMembership,
																			 @Nonnull List<Double> sampleValues) throws DimensionMismatchException {
		// incoming data structures should be same size because each parallel index represents the same
		// sample
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
		return batchLabels;
	}

	/**
	 * Returns a view on factorLabels list.
	 */
	public List<String> getFactorLabels() {
		return factorLabels;
	}

	/**
	 * Returns a view on batchValuesToInclude.
	 */
	public List<String> getBatchValuesToInclude() {
		return batchValuesToInclude;
	}

	/**
	 * Returns a view on factorValuesToInclude.
	 */
	public List<List<Double>> getFactorValuesToInclude() {
		return factorValuesToInclude;
	}

}
