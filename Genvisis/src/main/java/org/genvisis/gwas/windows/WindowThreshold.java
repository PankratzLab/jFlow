package org.genvisis.gwas.windows;

/**
 * Container structure for hit window detection parameters: index p-value, suggestive p-value, and
 * window size (in base pairs, one-sided)
 */
public class WindowThreshold {
	public static final double DEFAULT_INDEX_PVAL = 0.0000001;
	public static final double DEFAULT_SUGGESTIVE_PVAL = 0.00001;
	public static final int DEFAULT_WINDOW = 250000;

	private double iPval;
	private double sPval;
	private int window;

	/**
	 * Constructs a WindowThreshold with {@link DEFAULT_INDEX_PVAL}, {@link #DEFAULT_SUGGESTIVE_PVAL},
	 * and {@link #DEFAULT_WINDOW}
	 */
	public WindowThreshold() {
		iPval = DEFAULT_INDEX_PVAL;
		sPval = DEFAULT_SUGGESTIVE_PVAL;
		window = DEFAULT_WINDOW;
	}

	/**
	 * @param indexPval New index p-value
	 * @return This WindowThreshold, for method chaining
	 */
	public WindowThreshold index(double indexPval) {
		iPval = indexPval;
		return this;
	}

	/**
	 * @param suggestivePval New suggestive p-value
	 * @return This WindowThreshold, for method chaining
	 */
	public WindowThreshold suggestive(double suggestivePval) {
		sPval = suggestivePval;
		return this;
	}

	/**
	 * @param window New window size
	 * @return This WindowThreshold, for method chaining
	 */
	public WindowThreshold window(int window) {
		this.window = window;
		return this;
	}

	public double getIndexPval() {
		return iPval;
	}

	public double getSugPval() {
		return sPval;
	}

	public int getWindow() {
		return window;
	}
}
