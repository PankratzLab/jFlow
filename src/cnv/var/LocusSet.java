package cnv.var;

import java.io.Serializable;
import java.util.Arrays;

import common.Array;
import common.Logger;
import filesys.Segment;

public abstract class LocusSet<T extends Segment> implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private T[] loci;
	private boolean sorted;
	private Logger log;

	public LocusSet(final T[] loci, boolean sort, final Logger log) {
		super();
		this.loci = loci;
		this.sorted = false;
		this.log = log;
		if (sort) {
			sort();
		}
	}

	public void sort() {
		loci = putInOrder(loci, Segment.quicksort(loci));
		sorted = true;
	}

	public T[] getLoci() {
		return loci;
	}

	public T[] getOverLappingLoci(final Segment seg) {
		if (!sorted) {
			log.reportTimeError("Internal error: must sort internal segment array prior to overlap search");
			return null;
		} else {
			int[] indices = Segment.binarySearchForAllOverLappingIndices(seg, loci);
			if (indices == null) {
				return null;
			} else {
				return Array.subArray(loci, indices);
			}
		}
	}

	private T[] putInOrder(final T[] array, int[] order) {
		T[] newArray;

		newArray = Arrays.copyOf(array, array.length);
		for (int i = 0; i < order.length; i++) {
			newArray[i] = array[order[i]];
		}

		return newArray;
	}

}