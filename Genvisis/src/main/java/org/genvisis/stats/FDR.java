/**
 * 
 */
package org.genvisis.stats;

import java.util.Arrays;

import org.genvisis.CLI;
import org.genvisis.CLI.Arg;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;



/**
 * methods to estimate false discovery rate, from Stolen from
 * https://forums.roguewave.com/showthread.php?248-False-Discovery-Rate-methods-in-JMSL
 *
 */
public class FDR {
	double t;
	double n;

	/**
	 * Create a new FDR object.
	 *
	 * @param p a double array, the vector of p-values
	 * @param q a double, the False Discovery Rate level
	 */
	public FDR(double[] p, double q) {
		Arrays.sort(p);
		int v = p.length;
		double cvn = 0;
		for (int i = 0; i < v; i++) {
			cvn += 1.0 / (i + 1);
		}
		double cvid = 1;

		double[] test1 = new double[v];
		double[] test2 = new double[v];
		for (int i = 0; i < v; i++) {
			test1[i] = (i + 1) / (double) v * q / cvid;
			test2[i] = (i + 1) / (double) v * q / cvn;
		}

		int ixt = findlemax(p, test1);
		if (ixt != -1)
			t = p[ixt];
		int ixn = findlemax(p, test2);
		if (ixn != -1)
			n = p[ixn];
	}

	/**
	 * @return p-value threshold based on independence or positive dependence
	 */
	public double getThreshold() {
		return t;
	}

	/**
	 * @return Nonparametric p-value threshold
	 */
	public double getNThreshold() {
		return n;
	}

	int findlemax(double[] a, double[] b) {
		final int len = a.length;
		int ix = -1;
		for (int i = 0; i < len; i++) {
			if (a[i] <= b[i])
				ix = i;
		}
		return ix;
	}


	/**
	 * @param p array of pvalues
	 * @param q false discovery rate
	 * @return
	 */
	public static FDR compute(double[] p, double q) {
		return new FDR(p, q);

	}


	private static void computeFromFile(String pvalFile, double q, String output, Logger log) {
		double[] pvals = Array.toDoubleArray(HashVec.loadFileToStringArray(	pvalFile, false,
																																				new int[] {0}, false));
		log.reportTimeInfo("False discovery rate set to " + q);
		FDR f = compute(pvals, q);
		StringBuilder builder = new StringBuilder();
		builder.append("\nThreshold =\t" + f.getThreshold() + "\n");
		builder.append("N Threshold =\t" + f.getNThreshold());

		log.reportTimeInfo(builder.toString());
		Files.write(builder.toString(), output == null ? pvalFile + ".fdr" : output);

	}



	public static void main(String[] arg) {
		CLI c = new CLI(FDR.class);

		c.addArg("pvals", "full path to a file of p-values in the first column (no header)", true);
		c.addArgWithDefault("qvalue", "the False Discovery Rate level", "0.05", Arg.NUMBER);
		c.addArgWithDefault("output", "output file (defaults to file of p-values + \".fdr\"", null);

		c.parseWithExit(arg);

		computeFromFile(c.get("pvals"), c.getD("qvalue"), c.get("output"), new Logger());


	}


}
