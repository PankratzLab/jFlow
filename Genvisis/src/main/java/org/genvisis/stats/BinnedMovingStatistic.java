package org.genvisis.stats;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.collections4.list.TreeList;
import org.genvisis.common.ArrayUtils;

/**
 * Class for computing a moving statistic of binned data (represented by lists). For example, if
 * computing a moving mean income of a population, this class would allow you to bin people by age
 * and have a window that covers N years. This saves on computation, as the mean value for a given
 * year only needs to be computed once. Which statistic is computed is determined by the
 * {@link MovingStat} provided at construction.
 * <p>
 * Use: you must decide how your data is to be binned. For each sample, in bin-sorted order, compute
 * the bin number for that sample. Check {@link #lastBin(int)} to see if that sample will cause a
 * new bin to be created. If so, check {@link #getValue()} for the average value of all current bins
 * (which will be {@code -1} if we do not have sufficient bins yet). If you are using a two-sided
 * moving average, {@link #mid()} will give you the current median bin value. Finally, use
 * {@link #add(Number, int)} to add the value to the specified bin.
 * </p>
 * <p>
 * Inspired by
 * http://commons.apache.org/proper/commons-math/apidocs/org/apache/commons/math4/stat/descriptive/DescriptiveStatistics.html.
 * </p>
 */
public class BinnedMovingStatistic<T extends Number> {

	/**
	 * Available moving statistics to compute
	 */
	public enum MovingStat {
														MEAN, MAD
	}

	private BinManager<T> binManager;
	private List<Integer> binVals;
	private List<T> currentBin;

	/**
	 * Create a moving statistic which retains the specified number of bins in scope at any given
	 * time.
	 */
	public BinnedMovingStatistic(int window, MovingStat movingType) {
		switch (movingType) {
			case MEAN:
				binManager = new MeanBinManager<T>(window);
				break;
			case MAD:
				binManager = new MADBinManager<T>(window);
				break;
		}
		binVals = new LinkedList<Integer>();
		currentBin = new ArrayList<T>();
	}

	/**
	 * Add a bin to this statistics collection, and return the current stat value.
	 */
	public double add(T val, int bin) {
		double rVal = -1;
		if (binVals.isEmpty()) {
			// First bin value, so start tracking
			binVals.add(bin);
		} else if (bin != binVals.get(binVals.size() - 1)) {
			// Current bin is done. Add it to our running stats

			// If we've passed the window size, pop off the oldest bin
			if (binManager.add(currentBin)) {
				binVals.remove(0);
				rVal = getValue();
			}
			// Start building the next bin
			binVals.add(bin);
			currentBin = new ArrayList<T>();
		}
		currentBin.add(val);

		return rVal;
	}

	/**
	 * Returns the current moving statistic, not including the bin currently under construction. To
	 * get that value, first explicitly call {@link #binBreak()}.
	 *
	 * @return The mean of all values in all lists in scope, or -1 if insufficient values have been
	 *         added to this moving average.
	 */
	public double getValue() {
		return binManager.hasStat() ? binManager.getStat() : -1.0;
	}

	/**
	 * Forcibly stops construction of the current bin, adding its contents to the running total. This
	 * method does not need to be called under normal circumstances, but may be useful at the end of
	 * iteration.
	 */
	public void binBreak() {

	}

	/**
	 * Remove all values from this moving statistic collection
	 */
	public void clear() {
		binVals.clear();
		currentBin.clear();
		binManager.clear();
	}

	/**
	 * @return True if a) the given bin value matches the most recent bin, or b) there are currently
	 *         no bins. In either case, there is nothing interesting about the current statistics.
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
	 * Marker interface for moving statistic types
	 */
	private static interface BinManager<T extends Number> {

		/**
		 * Add the given bin to the current statistics
		 *
		 * @return true if the number of bins has exceeded the window
		 */
		boolean add(List<T> bin);

		/**
		 * Remove the earliest bin
		 */
		void evict();


		/**
		 * @return true if we have as many bins as our window size. False otherwise.
		 */
		boolean hasStat();

		/**
		 * @return The appropriate statistic covering all values in all bins
		 */
		double getStat();

		/**
		 * Drop all values from this statistic
		 */
		void clear();
	}

	/**
	 * Abstract superclass for {@link BinManager} implementations. Covers some boilerplate. Subclasses
	 * should calling {@code super.add} so they do not need to manually determine if {@link #evict()}
	 * is necessary. All subclasses should call {@code super.clear}.
	 */
	private abstract static class AbstractBinManager<T extends Number> implements BinManager<T> {

		private final int window;
		private int bins;

		public AbstractBinManager(int window) {
			this.window = window;
			bins = 0;
		}

		@Override
		public boolean add(List<T> bin) {
			bins++;
			if (bins > window) {
				evict();
				bins--;
				return true;
			}
			return false;
		}

		@Override
		public void clear() {
			bins = 0;
		}

		@Override
		public boolean hasStat() {
			return bins == window;
		}
	}

	/**
	 * {@link BinManager} implementations that returns the mean value in its {@link #getStat()}
	 * method.
	 */
	private static class MeanBinManager<T extends Number> extends AbstractBinManager<T> {

		private List<Double> sums;
		private List<Integer> counts;
		private int movingCount;
		private double movingSum;

		public MeanBinManager(int window) {
			super(window);
			sums = new LinkedList<Double>();
			counts = new LinkedList<Integer>();
		}

		@Override
		public boolean add(List<T> bin) {
			double sum = ArrayUtils.sum(bin);
			sums.add(sum);
			counts.add(bin.size());
			movingSum += sum;
			movingCount += bin.size();
			return super.add(bin);
		}

		@Override
		public double getStat() {
			return movingSum / movingCount;
		}

		@Override
		public void evict() {
			Double removedSum = sums.remove(0);
			Integer removedCount = counts.remove(0);
			movingSum -= removedSum;
			movingCount -= removedCount;
		}

		@Override
		public void clear() {
			super.clear();
			sums.clear();
			counts.clear();
			movingCount = 0;
			movingSum = 0;
		}
	}

	/**
	 * {@link BinManager} implementations that returns the median absolute difference value in its
	 * {@link #getStat()} method.
	 */
	private static class MADBinManager<T extends Number> extends AbstractBinManager<T> {

		private List<T> values;
		// A list of all the bins covered in this window
		private List<List<T>> bins;

		public MADBinManager(int window) {
			super(window);
			// Need to store the raw bin unfortunately, as it is necessary on eviction
			bins = new LinkedList<List<T>>();
			values = new TreeList<T>();
		}

		@Override
		public boolean add(List<T> bin) {
			// Convert the bin to a map of value > counts. Using a map reduces how frequently we have to
			// sort data, as normally the MAD would require two sorts.
			values.addAll(bin);
			bins.add(bin);
			return super.add(bin);
		}

		@Override
		public double getStat() {
			double median = ArrayUtils.medianSorted(values);

			TreeList<Double> absDiffs = new TreeList<Double>();
			for (T v : values) {
				absDiffs.add(Math.abs(v.doubleValue() - median));
			}

			return ArrayUtils.medianSorted(absDiffs);
		}

		@Override
		public void evict() {
			List<T> removed = bins.remove(0);
			// Have to update counts for each value in the removed bin
			for (T val : removed) {
				values.remove(val);
			}
		}

		@Override
		public void clear() {
			super.clear();
			bins.clear();
			values.clear();
		}
	}
}
