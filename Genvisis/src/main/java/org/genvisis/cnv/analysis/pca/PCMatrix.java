package org.genvisis.cnv.analysis.pca;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import org.apache.commons.math3.stat.inference.TestUtils;
import org.genvisis.CLI;
import org.genvisis.CLI.Arg;
import org.genvisis.cnv.Launch;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.stats.Ttest;

/**
 * Generates matrices plotting principal components against {@link SampleData} columns. This can be
 * used to identify which components sufficiently capture biases such as batch or plate effects, by
 * plotting against the respective column.
 */
public class PCMatrix {

	public static final String ARGS_PC_FILE = "pcFile";
	public static final String ARGS_COLUMN = "columnTitle";
	public static final String COLUMN_DEFAULT = "PLATE";
	public static final String OUT_DEFAULT = "results/";
	public static final String MENU_ENTRY = "Generate data column matrix";

	/**
	 * Create PC vs column matrices with default PC intensity file
	 *
	 * @see #createMatrix(Project, String, String, String)
	 *
	 * @param proj Base project
	 * @param column Column of interest in the sample data.
	 * @param pcFile Principal component file path
	 * @return {@code true} iff the matrices are output without error
	 */
	public static boolean createMatrix(Project proj, String column) {
		return createMatrix(proj, column, proj.INTENSITY_PC_FILENAME.getValue());
	}

	/**
	 * Create PC vs column matrices with default output directory
	 *
	 * @see #createMatrix(Project, String, String, String)
	 *
	 * @param proj Base project
	 * @param column Column of interest in the sample data.
	 * @param pcFile Principal component file path
	 * @return {@code true} iff the matrices are output without error
	 */
	public static boolean createMatrix(Project proj, String column, String pcFile) {
		return createMatrix(proj, column, pcFile, OUT_DEFAULT);
	}

	/**
	 * Create a matrix of principal components vs. sample data columns. The output data has 1 row per
	 * PC in the specified file and one column per unique value in the specified data column, with one
	 * extra column for ANOVA data. To compute the value of each position we split the sample data for
	 * the current row's PC based on the particular value of the current column. We compute a t-test
	 * on these data groups and record the -log10(p_val) of the test. This data can then be used to
	 * select a number of components that fully captures any bias effects of the given column.
	 *
	 * @param proj Base project
	 * @param column Column of interest in the sample data.
	 * @param pcFile Principal component file path
	 * @param outDir Directory to write generated files
	 * @return {@code true} iff the matrices are output without error
	 */
	public static boolean createMatrix(Project proj, String column, String pcFile, String outDir) {
		Logger log = proj.getLog();
		if (!outDir.endsWith("/")) {
			outDir = outDir + "/";
		}

		String matrixOutPath = proj.PROJECT_DIRECTORY.getValue() + outDir + column + "_PC_matrix.csv";
		String valuesOutPath = proj.PROJECT_DIRECTORY.getValue() + outDir + column + "_values.txt";

		log.report("Creating matrix of Principal Component vs " + column + "...");
		PrincipalComponentsResiduals pcResids = new PrincipalComponentsResiduals(	proj, pcFile, null,
																																							proj.INTENSITY_PC_NUM_COMPONENTS.getValue(),
																																							false, 0, false,
																																							false, null);

		// PC values by pc index, sample
		double[][] pcsForSample = pcResids.getPcBasis();
		String[] pcSamples = pcResids.getAllProjSamples();

		final String file = proj.SAMPLE_DATA_FILENAME.getValue();
		// Find the index of the column of interest
		String[] headerOfFile = Files.getHeaderOfFile(file, log);
		int[] valCols = ext.indexFactors(new String[] {column}, headerOfFile, false, false);
		int[] keyCols = ext.indexFactors(new String[] {"FID", "IID"}, headerOfFile, false, false);

		if (valCols == null || valCols.length == 0) {
			log.reportError(PCMatrix.class	+ ": requested column (" + column
											+ ") not found in current project.");
			return false;
		}

		// Get column values for all samples
		Hashtable<String, Vector<String>> sampleVec = HashVec.loadFileToHashVec(file, keyCols, valCols,
																																						"\t", true, false);

		// Align the PC values with sample data column values
		String[][] filtered = filterSamples(proj.getSampleData(0, false), pcSamples, sampleVec);
		String[] colValsBySample = filtered[0]; // All column values for samples with PC data
		String[] valList = filtered[1]; // Unique column values


		// output matrix for the resulting t-test values
		// Extra column for anova stat
		String[][] matrix = new String[pcsForSample.length][valList.length + 1];

		// Compute t-tests
		for (int i = 0; i < matrix.length; i++) {
			List<double[]> colData = new ArrayList<double[]>();
			for (int j = 0; j < valList.length; j++) {
				// we will split data for t-test based on this column
				String colOfInterest = valList[j];
				double[][] splitData = Ttest.splitOut(colValsBySample, pcsForSample[i], colOfInterest);
				// Store the column-specific data for this PC for later anova
				colData.add(splitData[0]);
				Ttest t = new Ttest(splitData);
				matrix[i][j] = getValue(t.getPvalue());
			}

			// Compute anova
			matrix[i][matrix[i].length - 1] = getValue(TestUtils.oneWayAnovaPValue(colData));
		}

		Files.writeArray(colValsBySample, valuesOutPath);
		Files.writeMatrix(matrix, matrixOutPath, ",");

		log.report("PC Matrix written to: " + matrixOutPath);

		return true;
	}

	private static String[][] filterSamples(SampleData sampleData, String[] pcSamples,
																		Hashtable<String, Vector<String>> sampleVec) {
		String[][] values = new String[2][];
		values[0] = new String[pcSamples.length];
		Set<String> unique = new HashSet<String>();

		for (int i=0; i<pcSamples.length; i++) {
			// Convert the sample ID to FID\tIID key
			String sample = sampleData.lookup(pcSamples[i])[1];
			Vector<String> vector = sampleVec.get(sample);
			values[0][i] = vector.get(0);
			unique.add(vector.get(0));
		}

		values[1] = unique.toArray(new String[unique.size()]);

		return values;
	}

	/**
	 * Convert doubles to printable values
	 *
	 * @return String representation of -log10(value)
	 */
	private static String getValue(double v) {
		return Double.toString(-Math.log10(v));
	}

	/**
	 * Command-line entry point
	 */
	public static void main(String... args) {
		CLI c = new CLI(PCMatrix.class);
		c.addArgWithDefault(CLI.ARG_PROJ, CLI.DESC_PROJ,
												Launch.getDefaultDebugProjectFile(false));
		c.addArgWithDefault(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, OUT_DEFAULT);
		c.addArgWithDefault(ARGS_COLUMN, "Sample data column to compare with PCs", COLUMN_DEFAULT);
		c.addArg(	ARGS_PC_FILE,
							"Input principal component file. If not specified, will use project default.", false,
							Arg.FILE);

		c.parseWithExit(args);

		Project proj = new Project(c.get(CLI.ARG_PROJ), false);
		String column = c.get(ARGS_COLUMN);
		String pcFile =
									c.has(ARGS_PC_FILE) ? c.get(ARGS_PC_FILE) : proj.INTENSITY_PC_FILENAME.getValue();
		String outDir = c.get(CLI.ARG_OUTDIR);

		createMatrix(proj, column, pcFile, outDir);
	}
}
