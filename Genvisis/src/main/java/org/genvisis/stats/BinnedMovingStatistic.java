package org.genvisis.stats;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Queue;
import java.util.TreeMap;

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
 * the bin number for that sample. Check {@link #inBin(int)} to see if that sample will cause a new
 * bin to be created. If so, check {@link #getValue()} for the average value of all current bins
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
	private int currentBinVal;
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
		binVals = new TreeList<Integer>();
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
			currentBinVal = bin;
		} else if (bin != currentBinVal) {
			// Current bin is done. Add it to our running stats

			// If we've passed the window size, pop off the oldest bin
			if (binManager.add(currentBin)) {
				binVals.remove(0);
				rVal = getValue();
			}
			// Start building the next bin
			binVals.add(bin);
			currentBinVal = bin;
			currentBin = new ArrayList<T>();
		}
		currentBin.add(val);

		return rVal;
	}

	/**
	 * Returns the current moving statistic, not including the bin currently under construction. To
	 * get that value, first explicitly call {@link #forceBinBreak()}.
	 *
	 * @return The mean of all values in all lists in scope, or -1 if insufficient values have been
	 *         added to this moving average.
	 */
	public double getValue() {
		return binManager.hasStat() ? binManager.getStat() : -1.0;
	}

	/**
	 * Forcibly stops construction of the current bin, adding its contents to the running total. This
	 * method does not need to be called under normal circumstances, but may be useful to end current
	 * iteration.
	 *
	 * @return true the current statistic has a valid value
	 */
	public boolean forceBinBreak() {
		if (!currentBin.isEmpty()) {
			binVals.add(currentBinVal);
			binManager.add(currentBin);
			currentBin = new ArrayList<T>();
		}
		return binManager.hasStat();
	}

	/**
	 * Forcibly removes the oldest bin. This method does not need to be called under normal
	 * circumstances, but may be useful to end current iteration.
	 *
	 * @return true the current statistic has a valid value
	 */
	public boolean forceBinPop() {
		if (!binVals.isEmpty()) {
			binManager.evict();
			binVals.remove(0);
		}
		return binManager.hasStat();
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
	public boolean inBin(int bin) {
		return binVals.isEmpty() || currentBinVal == bin;
	}

	/**
	 * @return bin value corresponding to current stat
	 */
	public int mid() {
		return binVals.get(binManager.fillOffset());
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

		/**
		 * @return Offset from the current bin to the bin position that corresponds to the value of
		 *         {@link #getStat()}. Reasonable expected values: 0 - use only preceding bins,
		 *         window_size - use only following bins, window_size/2 - use half the window on both
		 *         side..
		 */
		int fillOffset();
	}

	/**
	 * Abstract superclass for {@link BinManager} implementations. Covers some boilerplate. Subclasses
	 * should calling {@code super.add} so they do not need to manually determine if {@link #evict()}
	 * is necessary. All subclasses should call {@code super.clear} and {@code super.evict}.
	 * <p>
	 * NB: this class implements {@link #hasStat()} and {@link #fillOffset()} in a way to facilitate
	 * two-sided moving statistics that are "valid" when either window margin is satisfied.
	 * </p>
	 * <p>
	 * This superclass assumes moving bins over a continuous data range. If is designed to handle a
	 * single turning point in window size (e.g. effective window size grows at the beginning as you
	 * are adding more data, and shrinks at the tail end). {@link #clear()} will reset the statistic.
	 * </p>
	 */
	private abstract static class AbstractBinManager<T extends Number> implements BinManager<T> {

		private final int window;
		private final int minFill;
		private int bins;
		// Flag for whether or not we have ever reached the full window size
		private boolean filledWindow = false;

		public AbstractBinManager(int window) {
			this.window = window;
			minFill = (window / 2) + 1;
			bins = 0;
		}

		@Override
		public boolean add(List<T> bin) {
			bins++;
			if (bins > window) {
				evict();
				filledWindow = true;
				return true;
			}
			return false;
		}

		@Override
		public void evict() {
			bins--;
		}

		@Override
		public void clear() {
			filledWindow = false;
			bins = 0;
		}

		@Override
		public boolean hasStat() {
			return bins >= minFill;
		}

		@Override
		public int fillOffset() {
			// This behavior has to depend on whether or not the moving stat is
			// growing or shrinking. As you start to build the stat, you must use points ahead of you. As
			// your stat moves into terminal boundaries, it must use points behind it.
			return filledWindow ? minFill : (bins - minFill);
		}
	}

	/**
	 * {@link BinManager} implementations that returns the mean value in its {@link #getStat()}
	 * method.
	 */
	private static class MeanBinManager<T extends Number> extends AbstractBinManager<T> {

		private Queue<Double> sums;
		private Queue<Integer> counts;
		private int movingCount;
		private double movingSum;

		public MeanBinManager(int window) {
			super(window);
			sums = new ArrayDeque<Double>();
			counts = new ArrayDeque<Integer>();
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
			Double removedSum = sums.remove();
			Integer removedCount = counts.remove();
			movingSum -= removedSum;
			movingCount -= removedCount;
			super.evict();
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

		private Map<T, int[]> values;
		// A list of all the bins covered in this window
		private List<List<T>> bins;
		private int valueCount = 0;

		public MADBinManager(int window) {
			super(window);
			// Need to store the raw bin unfortunately, as it is necessary on eviction
			bins = new TreeList<List<T>>();
			values = new TreeMap<T, int[]>();
		}

		@Override
		public boolean add(List<T> bin) {
			// Convert the bin to a map of value > counts. Using a map reduces how frequently we have to
			// sort data, as normally the MAD would require two sorts.
			for (T v : bin) {
				add(v, values);
				valueCount++;
			}
			bins.add(bin);
			return super.add(bin);
		}

		@Override
		public double getStat() {
			double median = median(values);

			Map<Double, int[]> absDiffs = new TreeMap<Double, int[]>();
			for (Entry<T, int[]> entry : values.entrySet()) {
				Double diffKey = Math.abs(entry.getKey().doubleValue() - median);
				if (absDiffs.containsKey(diffKey)) {
					absDiffs.get(diffKey)[0] += entry.getValue()[0];
				} else {
					absDiffs.put(diffKey, new int[] {entry.getValue()[0]});
				}
			}

			return median(absDiffs);
		}

		@Override
		public void evict() {
			List<T> removed = bins.remove(0);
			// Have to update counts for each value in the removed bin
			for (T val : removed) {
				remove(val, values);
			}
			super.evict();
		}

		@Override
		public void clear() {
			super.clear();
			bins.clear();
			values.clear();
			valueCount = 0;
		}

		private void add(T value, Map<T, int[]> valueMap) {
			int[] counts = valueMap.get(value);
			if (counts == null) {
				counts = new int[1];
				valueMap.put(value, counts);
			}
			counts[0]++;
		}

		private void remove(T value, Map<T, int[]> valueMap) {
			int[] counts = valueMap.get(value);
			counts[0]--;
			if (counts[0] == 0) {
				valueMap.remove(value);
			}
			valueCount--;
		}

		private <N extends Number> double median(Map<N, int[]> valueMap) {
			N p1 = null;
			N p2 = null;
			final int idx1 = valueCount / 2;
			final int idx2 = valueCount % 2 == 0 ? (valueCount / 2) - 1 : idx1;
			int pos = 0;
			for (Entry<N, int[]> entry : valueMap.entrySet()) {
				pos += entry.getValue()[0];
				if (p1 == null && idx1 < pos) {
					p1 = entry.getKey();
				}
				if (p2 == null && idx2 < pos) {
					p2 = entry.getKey();
				}
				if (p1 != null && p2 != null) {
					break;
				}
			}
			return (p1.doubleValue() + p2.doubleValue()) / 2;
		}
	}
}
