/**
 * 
 */
package org.genvisis.cnv.wbs;

import java.util.Comparator;

class ChangePoint implements Comparable<ChangePoint> {



	private int start;
	private int end;
	private int changePoint;
	private double cusum;
	private double minth;
	private int scale;

	ChangePoint(int start, int end, int changePoint, double cusum, double minth,
							int scale) {
		super();
		this.start = start;
		this.end = end;
		this.changePoint = changePoint;
		this.cusum = cusum;
		this.minth = minth;
		this.scale = scale;
	}

	static double[] extract() {
		return null;
	}

	int getStart() {
		return start;
	}

	int getEnd() {
		return end;
	}

	int getChangePoint() {
		return changePoint;
	}

	double getCusum() {
		return cusum;
	}

	double getMinth() {
		return minth;
	}

	int getScale() {
		return scale;
	}

	@Override
	public String toString() {
		return "WBSChangePoint [start=" + start + ", end=" + end + ", changePoint=" + changePoint
					 + ", cusum=" + cusum + ", minth=" + minth + ", scale=" + scale + "]";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(ChangePoint o) {
		return -1 * Double.compare(Math.abs(getCusum()), Math.abs(o.getCusum()));
	}

	static class ChangePointComparable implements Comparator<ChangePoint> {

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
		 */
		@Override
		public int compare(ChangePoint o1, ChangePoint o2) {
			// increasing
			return Double.compare(o1.getChangePoint(), o2.getChangePoint());
		}

	}

	static class MinthComparable implements Comparator<ChangePoint> {

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
		 */
		@Override
		public int compare(ChangePoint o1, ChangePoint o2) {
			return -1 * Double.compare(o1.getMinth(), o2.getMinth());
		}

	}

}
