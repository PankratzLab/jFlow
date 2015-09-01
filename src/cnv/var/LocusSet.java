package cnv.var;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Vector;

import common.Array;
import common.Logger;
import filesys.Segment;

public abstract class LocusSet<T extends Segment> implements Serializable {
	/**
	 * 
	 */
	protected static final long serialVersionUID = 1L;
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
		loci = putInOrder(loci, Segment.quicksort(loci, false));
		sorted = true;
	}

	/**
	 * @return simply the sum of the segments length, should be a merged set if you want the unique coverage
	 */
	public int getBpCovered() {
		int sum = 0;
		for (int i = 0; i < loci.length; i++) {
			sum += loci[i].getSize();
		}
		return sum;
	}

	public T[] getLoci() {
		return loci;
	}

	public LocusSet<Segment> mergeOverlapping() {

		if (!sorted) {
			log.reportTimeError("Internal error: must sort internal segment array prior to merge");
			return null;
		} else {
			byte currentChr = -1;
			ArrayList<Segment> merged = new ArrayList<Segment>();
			Vector<Segment> tmp = new Vector<Segment>();
			int originalSize = loci.length;
			for (int i = 0; i < loci.length; i++) {
				if (loci[i].getChr() != currentChr) {
					if (tmp.size() > 0) {
						Segment.mergeOverlapsAndSort(tmp);
						for (int j = 0; j < tmp.size(); j++) {
							merged.add(tmp.get(j));
						}
						tmp.clear();
					}
				}
				currentChr = loci[i].getChr();
				tmp.add(loci[i]);

			}
			if (tmp.size() > 0) {
				Segment.mergeOverlapsAndSort(tmp);
				for (int j = 0; j < tmp.size(); j++) {
					merged.add(tmp.get(j));
				}
				tmp.clear();
			}
			log.reportTimeInfo("Merged "+originalSize +" segments to "+merged.size());
			LocusSet<Segment> mergedSet = new LocusSet<Segment>(merged.toArray(new Segment[merged.size()]), true, log) {

				/**
				 * 
				 */
				private static final long serialVersionUID = 1L;

			};
			return mergedSet;
		}
	}

	public int[] getOverlappingIndices(final Segment seg) {
		if (!sorted) {
			log.reportTimeError("Internal error: must sort internal segment array prior to overlap search");
			return null;
		} else {
			int[] indices = Segment.binarySearchForAllOverLappingIndices(seg, loci);
			return indices;

		}
	}

	public T[] getOverLappingLoci(final Segment seg) {
		int[] indices = getOverlappingIndices(seg);
		if (indices == null) {
			return null;
		} else {
			return Array.subArray(loci, indices);
		}
	}

	public enum TO_STRING_TYPE {
		/**
		 * calls the {@link Segment#getUCSClocation()}, or any overide
		 */
		UCSC, /**
		 * calls the {@link Segment#toString()}, or any override
		 */
		REGULAR;
	}

	/**
	 * @param array
	 *            writes according to the to string method of the class
	 * @param filename
	 * @param log
	 * @return
	 */
	public boolean writeRegions(String filename, TO_STRING_TYPE type, boolean header, Logger log) {
		log.reportTimeInfo("Writing " + loci.length + " loci to " + filename);
		boolean written = true;
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(filename));
			for (int i = 0; i < loci.length; i++) {
				if (i == 0 && header) {
					writer.println(Array.toStr(loci[i].getHeader()));
				}
				switch (type) {
				case REGULAR:
					writer.println(loci[i].toAnalysisString());
					break;
				case UCSC:
					writer.println(loci[i].getUCSClocation());
					break;
				default:
					log.reportTimeError("Invalid type " + type);
					written = false;
					break;
				}
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + filename);
			log.reportException(e);
			written = false;
		}
		return written;
	}

	private T[] putInOrder(final T[] array, int[] order) {
		T[] newArray;

		newArray = Arrays.copyOf(array, array.length);
		for (int i = 0; i < order.length; i++) {
			newArray[i] = array[order[i]];
		}

		return newArray;
	}

	public static LocusSet<Segment> loadSegmentSetFromFile(String file, int chrCol, int startCol, int stopCol, int skipNumLines, boolean inclusiveStart, boolean inclusiveStop, int bpBuffer, Logger log) {
		Segment[] segs = Segment.loadRegions(file, chrCol, startCol, stopCol, skipNumLines, inclusiveStart, inclusiveStop, bpBuffer);
		LocusSet<Segment> lSet = new LocusSet<Segment>(segs, true, log) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};
		return lSet;
	}

}