package common;

import htsjdk.samtools.Cigar;

import java.util.ArrayList;

public class ArraySpecialList {

	public static class ArrayIntList extends ArrayList<Integer> {
		private int capacity;
		public ArrayIntList(int capacity) {
			super(capacity);
		}

	}

	public static class ArrayCigarList extends ArrayList<Cigar> {
		private int capacity;

		public ArrayCigarList(int capacity) {
			super(capacity);
		}

	}

}
