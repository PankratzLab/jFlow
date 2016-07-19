package org.genvisis.common;

import htsjdk.samtools.Cigar;

import java.util.ArrayList;

import org.genvisis.cnv.annotation.BlastAnnotationTypes.BlastAnnotation;

public class ArraySpecialList {
	public static class ArrayStringList extends ArrayList<String> {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		public ArrayStringList(int capacity) {
			super(capacity);
		}

	}
	
	
	
	public static class ArrayIntList extends ArrayList<Integer> {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		public ArrayIntList(int capacity) {
			super(capacity);
		}

	}

	public static class ArrayCigarList extends ArrayList<Cigar> {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		public ArrayCigarList(int capacity) {
			super(capacity);
		}

	}

	public static class ArrayBlastAnnotationList extends ArrayList<BlastAnnotation> {
		private double maxEval;
		private double minEval;
		private int maxEvalIndex;
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		public ArrayBlastAnnotationList(int capacity) {
			super(capacity);
			this.maxEval = -1;
			this.minEval = Double.MAX_VALUE;
			this.maxEvalIndex = -1;
		}

		public void update() {
			this.maxEvalIndex = -1;
			this.minEval = Double.MAX_VALUE;
			this.maxEval = -1;
			for (int i = 0; i < size(); i++) {
				BlastAnnotation tmpBa = get(i);
				double tmp = tmpBa.geteValue();
				if (tmp > maxEval) {
					maxEval = tmp;
					maxEvalIndex = i;
				}
				if (tmp < minEval) {
					minEval = tmp;
				}
			}
		}

		public boolean add(BlastAnnotation blastAnnotation) {
			boolean add = super.add(blastAnnotation);
			if (add) {
				double tmp = blastAnnotation.geteValue();
				if (tmp > maxEval) {
					maxEval = tmp;
					maxEvalIndex = size() - 1;

				}
				if (tmp < minEval) {
					minEval = tmp;
				}
			}
			return true;
		}

		public double getMaxEval() {
			return maxEval;
		}

		public double getMinEval() {
			return minEval;
		}

		public int getMaxEvalIndex() {
			return maxEvalIndex;
		}

	}

}
