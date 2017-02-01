package org.genvisis.stats;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.genvisis.common.ArrayUtils;

/**
 * Class for computing a moving average of binned data (represented by lists). For example, if
 * computing a moving average income of a population, this class would allow you to bin people by
 * age and have a window that covers N years. This saves on computation, as the mean value for a
 * given year only needs to be computed once.
 * <p>
 * Use: you must decide how your data is to be binned. For each sample, in bin-sorted order, compute
 * the bin number for that sample. Check {@link #lastBin(int)} to see if that sample will cause a
 * new bin to be created. If so, check {@link #getValue()} for the average value of all current bins
 * (which will be {@code -1} if we do not have sufficient bins yet). If you are using a two-sided
 * moving average, {@link #mid()} will give you the current median bin value. Finally, use
 * {@link #add(Number, int)} to add the value to the specified bin. If you need to check the moving
 * average as though the "in process" bin has been filled, e.g. if your last sample did not trigger addition of a new bin, use
 * {@link #getValWithCurrentBin()}
 * </p>
 * <p>
 * Inspired by
 * http://commons.apache.org/proper/commons-math/apidocs/org/apache/commons/math4/stat/descriptive/DescriptiveStatistics.html.
 * </p>
 */
public class BinnedMovingAverage<T extends Number> {

	private final int window;
	private double movingSum;
	private double movingCount;

	private List<BinStat<T>> bins;
	private List<T> currentBin;
	private List<Integer> binVals;

	/**
	 * Create a moving average which retains the specified number of lists in scope at any given time.
	 */
	public BinnedMovingAverage(int window) {
		this.window = window;
		bins = new LinkedList<BinStat<T>>();
		binVals = new LinkedList<Integer>();
		currentBin = new ArrayList<T>();
	}

	/**
	 * Add a bin to this statistics collection, and return the current average value.
	 */
	public double add(T val, int bin) {
		if (binVals.isEmpty()) {
			// First bin value, so start tracking
			binVals.add(bin);
		} else if (bin != binVals.get(binVals.size() - 1)) {
			// Current bin is done. Add it to our running stats
			BinStat<T> stat = new BinStat<T>(currentBin);
			movingSum += stat.sum();
			movingCount += stat.count();
			bins.add(stat);
			// If we've passed the window size, pop off the oldest bin
			if (bins.size() > window) {
				BinStat<T> removed = bins.remove(0);
				movingSum -= removed.sum();
				movingCount -= removed.count();
				binVals.remove(0);
			}
			// Start building the next bin
			binVals.add(bin);
			currentBin = new ArrayList<T>();
		}
		currentBin.add(val);

		return getValue();
	}

	/**
	 * Returns the current running average not including the bin currently under construction. To
	 * include those values, use {@link #getValWithCurrentBin()}
	 *
	 * @return The mean of all values in all lists in scope, or -1 if insufficient values have been
	 *         added to this moving average.
	 */
	public double getValue() {
		if (movingCount < window) {
			return -1;
		}
		return movingSum / movingCount;
	}

	/**
	 * @return The average of all bins, temporarily treating the bin under construction as complete.
	 */
	public double getValWithCurrentBin() {
		if (currentBin.isEmpty()) {
			return getValue();
		}
		BinStat<T> stat = new BinStat<T>(currentBin);
		// Replace the oldest bin's values with the current.
		return (movingSum + stat.sum() - bins.get(0).sum())
		       / (movingCount + stat.count() - bins.get(0).count());
	}

	/**
	 * Remove all values from this moving average
	 */
	public void clear() {
		bins.clear();
		binVals.clear();
		movingSum = 0;
		movingCount = 0;
	}

	/**
	 * @return True if a) the given bin value matches the most recent bin, or b) there are currently
	 *         no bins. In either case, there is nothing interesting about the current moving average.
	 */
	public boolean lastBin(int bin) {
		return binVals.isEmpty() || binVals.get(binVals.size() - 1).compareTo(bin) == 0;
	}

	/**
	 * @return Median bin value
	 */
	public int mid() {
		return binVals.get(binVals.size() / 2);
	}

	/**
	 * Helper class to avoid recomputation of statistics for an individual list.
	 */
	private static class BinStat<T extends Number> {

		private final double sum;
		private final int count;

		public BinStat(List<T> bin) {
			count = bin.size();
			sum = ArrayUtils.sum(bin);
		}

		/**
		 * @return Sum of all values in this list
		 */
		public double sum() {
			return sum;
		}

		/**
		 * @return Number of values in this list
		 */
		public int count() {
			return count;
		}
	}
}
