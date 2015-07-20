package common;

import java.util.ArrayList;

public class ArraySpecialList {

	public static class ArrayIntList extends ArrayList<Integer> {
		private int capacity;
		public ArrayIntList(int capacity) {
			super(capacity);
		}

	}

}
