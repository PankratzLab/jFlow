package common;

import htsjdk.samtools.Cigar;

import java.util.ArrayList;

import cnv.annotation.BlastAnnotationTypes.BlastAnnotation;

public class ArraySpecialList {

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

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		public ArrayBlastAnnotationList(int capacity) {
			super(capacity);
		}

	}

}
