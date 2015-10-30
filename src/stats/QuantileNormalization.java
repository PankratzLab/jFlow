package stats;

import java.util.Arrays;

import common.Array;
import common.Logger;
import common.Sort;
import common.ext;

public class QuantileNormalization {
	private double thresholdFactor;
	private double[][] dataToNorm;
	private double[][] normData;
	private Logger log;

	/**
	 * @param dataToNorm
	 * @param thresholdDifference
	 *            the difference between the
	 * @param log
	 */
	public QuantileNormalization(double[][] dataToNorm, Logger log) {
		super();
		this.dataToNorm = dataToNorm;
		this.thresholdFactor = Double.NaN;
		this.log = log;
	}

	public void setThresholdFactor(double thresholdFactor) {
		this.thresholdFactor = thresholdFactor;
	}

	public double[][] getNormData() {
		return normData;
	}

	public void verify() {
		if (dataToNorm == null || dataToNorm.length < 2) {
			String error = "Data must be present and length greater than 1";
			log.reportTimeError(error);
			throw new IllegalArgumentException(error);
		} else {
			int base = 0;
			for (int i = 0; i < dataToNorm.length; i++) {
				if (dataToNorm[i] == null) {
					String error = "Data arrays invalid";
					log.reportTimeError(error);
					throw new IllegalArgumentException(error);
				}
				if (i == 0) {
					base = dataToNorm[i].length;
				} else if (base != dataToNorm[i].length) {
					String error = "Data arrays have inconsistant length";
					log.reportTimeError(error);
					throw new IllegalArgumentException(error);
				}
			}
			for (int i = 0; i < dataToNorm[0].length; i++) {
				double[] tmp = new double[dataToNorm.length];
				for (int j = 0; j < tmp.length; j++) {
					tmp[j] = dataToNorm[j][i];
				}
				int numNan = tmp.length - Array.removeNaN(tmp).length;
				if (numNan != 0 && numNan != tmp.length) {
					String error = i + "\tAll indices must either have complete (non-NaN) data, or all NaN data " + numNan + " Nans out of " + tmp.length;
					log.reportTimeError(error);
					throw new IllegalArgumentException(error);
				}
			}
		}
	}

	public void normalize() {
		int[][] ranks = new int[dataToNorm.length][];
		for (int i = 0; i < ranks.length; i++) {
			ranks[i] = Sort.trickSort(dataToNorm[i]);
		}
		double[] rankMean = new double[dataToNorm[0].length];
		for (int i = 0; i < dataToNorm[0].length; i++) {
			double[] tmp = new double[dataToNorm.length];
			for (int j = 0; j < tmp.length; j++) {
				int nextRank = ranks[j][i];
				tmp[j] = dataToNorm[j][nextRank];
			}
			rankMean[i] = Array.mean(tmp, true);
		}

		this.normData = new double[dataToNorm.length][dataToNorm[0].length];
		for (int i = 0; i < dataToNorm.length; i++) {
			for (int j = 0; j < dataToNorm[0].length; j++) {

				normData[i][j] = rankMean[ranks[i][j]];
				if (!Double.isNaN(thresholdFactor)) {
					if (Double.isNaN(dataToNorm[i][j])) {
						normData[i][j] = Double.NaN;
					} else if (dataToNorm[i][j] * thresholdFactor < normData[i][j]) {
						normData[i][j] = thresholdFactor * dataToNorm[i][j];
					}
				}
			}
		}
	}
}
