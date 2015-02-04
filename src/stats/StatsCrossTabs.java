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
	public static enum CORREL_TYPE {
		SPEARMAN, PEARSON
	}

	private static enum VALUE_TYPE {
		STAT, PVALUE
	}

	private double[][] data;
	private double[][] statisticTable, sigTable;
	private String[] dataTitles;
	private boolean verify;
	private boolean verbose;
	private Logger log;
	private CORREL_TYPE cType;

	/**
	 * @param data
	 *            organized as data[variable][dataForVariable]
	 * @param dataTitles
	 *            must be the same length as data
	 * @param cType
	 *            type of test to run, see {@link CORREL_TYPE}
	 * 
	 * @param log
	 */
	public StatsCrossTabs(double[][] data, String[] dataTitles, CORREL_TYPE cType, boolean verbose, Logger log) {
		super();
		this.verbose = verbose;
		this.data = data;
		this.dataTitles = dataTitles;
		this.log = log;
		this.statisticTable = new double[data.length][data.length];
		this.sigTable = new double[data.length][data.length];
		this.cType = cType;
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
						switch (cType) {
						case PEARSON:
							result = stats.Correlation.Pearson(dataToCorrel);
							stat = result[0];
							sig = result[1];
							break;
						case SPEARMAN:
							result = stats.Correlation.Spearman(dataToCorrel);
							stat = result[0];
							sig = result[1];
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

	public void dumpTables(String fullPathToOutputBase) {
		dump(ext.addToRoot(fullPathToOutputBase, "." + cType + "_" + VALUE_TYPE.STAT), dataTitles, statisticTable, cType, VALUE_TYPE.STAT, log);
		dump(ext.addToRoot(fullPathToOutputBase, "." + cType + "_" + VALUE_TYPE.PVALUE), dataTitles, sigTable, cType, VALUE_TYPE.PVALUE, log);

	}

	private static void dump(String fullPathToFile, String[] dataTitles, double[][] dataTable, CORREL_TYPE cType, VALUE_TYPE vType, Logger log) {
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
		StatsCrossTabs cTable = new StatsCrossTabs(new double[][] { { 1, 2, 3, 4 }, { 1, 2, 3, 4 }, { 1, 6, 3, 5 }, { 1, 6, 3, Double.NaN } }, new String[] { "1", "2", "3", "hasNaN" }, CORREL_TYPE.SPEARMAN, true, new Logger());
		cTable.computeTable();
		cTable.dumpTables("D:/data/singapore_PCs/Manhattan/test.correl");
	}

	public static void main(String[] args) {
		test();
	}

}
