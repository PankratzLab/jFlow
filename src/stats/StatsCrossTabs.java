package stats;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;

import common.Array;
import common.Logger;
import common.ext;

/**
 * Generates a matrix for comparing multiple columns of data,
 *
 */
public class StatsCrossTabs {
	public static enum STAT_TYPE {
		SPEARMAN_CORREL, PEARSON_CORREL, LIN_REGRESSION

	}

	public static enum VALUE_TYPE {
		STAT, PVALUE
	}

	private double[][] data, indeps;
	private double[][] statisticTable, sigTable;
	private String[] dataTitles;
	private boolean verify;
	private boolean verbose;
	private Logger log;
	private STAT_TYPE sType;

	/**
	 * @param data
	 *            organized as data[variable][dataForVariable]
	 * @param data
	 *            organized as data[dataForVariable][variable], must be the same length as data
	 * 
	 * @param dataTitles
	 *            must be the same length as data
	 * 
	 * @param sType
	 *            type of test to run, see {@link STAT_TYPE}
	 * 
	 * @param log
	 */
	public StatsCrossTabs(double[][] data, double[][] indeps, String[] dataTitles, STAT_TYPE sType, boolean verbose, Logger log) {
		super();
		this.verbose = verbose;
		this.data = data;
		this.dataTitles = dataTitles;
		this.log = log;
		this.statisticTable = new double[data.length][data.length];
		this.sigTable = new double[data.length][data.length];
		this.sType = sType;
		this.indeps = indeps;
		this.verify = verify();
	}

	public void computeTable() {
		if (verify) {
			Hashtable<Integer, Integer> complete = new Hashtable<Integer, Integer>();
			for (int i = 0; i < statisticTable.length; i++) {
				for (int j = 0; j < statisticTable.length; j++) {
					double stat = Double.NaN;
					double sig = Double.NaN;

					if (!complete.containsKey(j)) {
						double[][] dataToCorrel = cleanNaNs(new double[][] { data[i], data[j] }, new String[] { dataTitles[i], dataTitles[j] }, verbose, log);
						double[] result;
						switch (sType) {
						case PEARSON_CORREL:
							result = stats.Correlation.Pearson(dataToCorrel);

							stat = result[0];
							sig = result[1];
							break;
						case SPEARMAN_CORREL:
							result = stats.Correlation.Spearman(dataToCorrel);
							stat = result[0];
							sig = result[1];
							break;
						case LIN_REGRESSION:
							double[][] newIndeps = new double[data[i].length][indeps == null ? 1 : 1 + indeps[0].length];
							for (int k = 0; k < newIndeps.length; k++) {
								for (int k2 = 0; k2 < newIndeps[0].length; k2++) {
									newIndeps[k][0] = data[1][k];
									newIndeps[k][k2 + 1] = indeps[k][k2];
								}
							}
							log.reportTimeInfo("Using " + newIndeps[0].length + " independant predictors for regression");
							RegressionModel model = new LeastSquares(data[i], newIndeps, null, false, verbose, false);
							stat = model.getRsquare();
							sig = model.getOverallSig();
							break;
						default:
							break;
						}
					}
					statisticTable[j][i] = stat;
					sigTable[j][i] = sig;
				}
				complete.put(i, i);
			}
		}
	}

//	public int[] getInOrder(VALUE_TYPE type) {
//		switch (type) {
//		case PVALUE:
//			return Array
//			break;
//		case STAT:
//			break;
//		default:
//			log.reportTimeError("Invalid type " + type);
//			break;
//
//		}
//	}

	public void dumpTables(String fullPathToOutputBase) {
		dump(ext.addToRoot(fullPathToOutputBase, "." + sType + "_" + VALUE_TYPE.STAT), dataTitles, statisticTable, sType, VALUE_TYPE.STAT, log);
		dump(ext.addToRoot(fullPathToOutputBase, "." + sType + "_" + VALUE_TYPE.PVALUE), dataTitles, sigTable, sType, VALUE_TYPE.PVALUE, log);

	}

	private static void dump(String fullPathToFile, String[] dataTitles, double[][] dataTable, STAT_TYPE cType, VALUE_TYPE vType, Logger log) {
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(fullPathToFile));
			writer.println(cType + "_" + vType + "\t" + Array.toStr(dataTitles));
			for (int i = 0; i < dataTable.length; i++) {
				writer.println(dataTitles[i] + "\t" + Array.toStr(dataTable[i]));
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + fullPathToFile);
			log.reportException(e);
		}
	}

	private boolean verify() {
		boolean verify = true;
		if (dataTitles.length != data.length) {
			log.reportTimeError("Data titles and data matrix must be the same size");
			verify = false;
		}
		if (indeps != null && data[0].length != indeps.length) {
			log.reportTimeError("Independant predictors and data must be the same size");
			verify = false;
		}
		if (sType == STAT_TYPE.LIN_REGRESSION && indeps == null) {
			log.reportTimeWarning("Independent predictors were not provided for stat type " + sType);
		} else if (sType != STAT_TYPE.LIN_REGRESSION && indeps != null) {
			log.reportTimeWarning("Independent predictors were provided for stat type " + sType + " and will be skipped");

		}
		return verify;
	}

	private static double[][] cleanNaNs(double[][] dataToClean, String[] titles, boolean verbose, Logger log) {
		int numNansC1 = 0;
		int numNansC2 = 0;

		ArrayList<Double> c1 = new ArrayList<Double>(dataToClean[0].length);
		ArrayList<Double> c2 = new ArrayList<Double>(dataToClean[0].length);

		for (int i = 0; i < dataToClean[0].length; i++) {
			boolean keep = true;
			if (Double.isNaN(dataToClean[0][i])) {
				numNansC1++;
				keep = false;
			}
			if (Double.isNaN(dataToClean[1][i])) {
				numNansC2++;
				keep = false;
			}
			if (keep) {
				c1.add(dataToClean[0][i]);
				c2.add(dataToClean[1][i]);
			}
		}
		if ((numNansC1 > 0 || numNansC2 > 0) && verbose) {
			log.reportTimeWarning(titles[0] + " had " + numNansC1 + " NaN value(s)");
			log.reportTimeWarning(titles[1] + " had " + numNansC2 + " NaN value(s)");
			log.reportTimeWarning("Retaining " + c1.size() + " matched values for correlations");
		}
		return new double[][] { Array.toDoubleArray(c1), Array.toDoubleArray(c2) };
	}

	public static void test() {
		StatsCrossTabs cTable = new StatsCrossTabs(new double[][] { { 1, 2, 3, 4 }, { 1, 2, 3, 4 }, { 1, 6, 3, 5 }, { 1, 6, 3, Double.NaN } }, null, new String[] { "1", "2", "3", "hasNaN" }, STAT_TYPE.SPEARMAN_CORREL, true, new Logger());
		cTable.computeTable();
		cTable.dumpTables("D:/data/singapore_PCs/Manhattan/test.correl");
	}

	public static void main(String[] args) {
		test();
	}

}
