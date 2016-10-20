/**
 * 
 */
package org.genvisis.stats;

import java.text.DecimalFormat;

import org.apache.commons.math3.random.EmpiricalDistribution;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

/**
 * @author Kitty
 * 
 *         Apache commons based histogram
 *
 */
public class BasicHistogram {

	private long[] counts;
	private double[] binMax;
	private double[] binMin;

	private BasicHistogram(long[] counts, double[] binMax, double[] binMin) {
		super();
		this.counts = counts;
		this.binMax = binMax;
		this.binMin = binMin;
	}

	public long[] getCounts() {
		return counts;
	}

	public double[] getBinMax() {
		return binMax;
	}

	public double[] getBinMin() {
		return binMin;
	}

	public String[] getBinLabels() {
		return getBinLabels("#.#");
	}

	private String[] getBinLabels(String format) {
		DecimalFormat d = new DecimalFormat(format);
		String[] labels = new String[getCounts().length];
		for (int i = 0; i < labels.length; i++) {
			String key = d.format(binMin[i]) + "-" + d.format(binMax[i]);
			labels[i] = key;
		}

		return labels;
	}

	/**
	 * @param binCount number of bins to use
	 * @param data data to create a histogram for
	 * @return
	 */
	public static BasicHistogram getHistogram(int binCount, double[] data) {
		long[] counts = new long[binCount];
		double[] binMax = new double[binCount];
		double[] binMin = new double[binCount];

		EmpiricalDistribution distribution = new EmpiricalDistribution(binCount);
		distribution.load(data);
		int k = 0;
		for (SummaryStatistics stats : distribution.getBinStats()) {
			counts[k] = stats.getN();
			binMax[k] = stats.getMax();
			binMin[k] = stats.getMin();
			k++;
		}
		return new BasicHistogram(counts, binMax, binMin);
	}
}


