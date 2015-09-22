package cnv.var;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Vector;

import common.Array;
import common.Files;
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

	public LocusSet(final List<T> lociAl, boolean sort, final Logger log) {
		super();
		if (lociAl == null || lociAl.size() < 1) {
			String error = "BUG: constructor call with 0 size or null in set contruction";
			log.reportTimeError(error);
			throw new IllegalArgumentException(error);
		}

		@SuppressWarnings("unchecked")
		T[] loci = (T[]) java.lang.reflect.Array.newInstance(lociAl.get(0).getClass(), lociAl.size());
		for (int i = 0; i < loci.length; i++) {
			loci[i] = lociAl.get(i);
		}
		this.loci = loci;
		this.sorted = false;
		this.log = log;
		if (sort) {
			sort();
		}
	}

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
	public long getBpCovered() {
		long sum = 0;
		for (int i = 0; i < loci.length; i++) {
			sum += loci[i].getSize();
		}
		return sum;
	}

	public T[] getLoci() {
		return loci;
	}

	/**
	 * @param setToRemove
	 *            remove the bases represented by this set from the current loci
	 * @param bpBuffer
	 *            the buffer to extend on either side of the removing set's loci
	 * @return
	 * 
	 * 
	 */
	public <E extends Segment> LocusSet<Segment> removeThese(final LocusSet<E> setToRemove, int bpBuffer) {
		ArrayList<Segment> newLoci = new ArrayList<Segment>();
		LocusSet<Segment> operateSet = setToRemove.getStrictSegmentSet();
		if (bpBuffer > 0) {
			operateSet = setToRemove.getBufferedSegmentSet(bpBuffer);
		}
		for (int i = 0; i < loci.length; i++) {
			Segment[] overlaps = operateSet.getOverLappingLoci(loci[i]);// exons in the bin
			if (overlaps == null || overlaps.length == 0) {
				newLoci.add(loci[i]);
			} else {
				LocusSet<Segment> overlapsRemoved = loci[i].removeAll(overlaps, log);
				for (int j = 0; j < overlapsRemoved.getLoci().length; j++) {
					newLoci.add(overlapsRemoved.getLoci()[j]);
				}
			}
		}
		LocusSet<Segment> toReturn = new LocusSet<Segment>(newLoci.toArray(new Segment[newLoci.size()]), true, log) {
			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};
		for (int i = 0; i < operateSet.getLoci().length; i++) {
			if (toReturn.getOverLappingLoci(operateSet.getLoci()[i]) != null) {
				String error = "BUG: found overlapping loci from the removed set in the set to be returned";
				log.reportTimeError(error);
				throw new IllegalStateException(error);
			}
		}
		return toReturn.mergeOverlapping();
	}

	public LocusSet<Segment> getBufferedSegmentSet(int bpBuffer) {
		ArrayList<Segment> buffered = new ArrayList<Segment>();
		for (int i = 0; i < loci.length; i++) {
			buffered.add(loci[i].getBufferedSegment(bpBuffer));
		}

		LocusSet<Segment> bufSet = new LocusSet<Segment>(buffered.toArray(new Segment[buffered.size()]), true, log) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};
		return bufSet;
	}

	public boolean hasNoOverlap() {
		boolean hasOverlap = false;
		out: for (int i = 0; i < loci.length; i++) {
			T[] overlaps = getOverLappingLoci(loci[i]);
			if (overlaps != null && overlaps.length > 0) {
				for (int j = 0; j < overlaps.length; j++) {
					if (!overlaps[j].equals(loci[i])) {
						hasOverlap = true;
						break out;
					}
				}
			}
		}
		return hasOverlap;
	}

	public LocusSet<Segment> mergeOverlapping() {
		return mergeOverlapping(false);
	}

	public LocusSet<Segment> mergeOverlapping(boolean verbose) {// TODO, use clone and set positions instead to get same type returned

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
			if (verbose) {
				log.reportTimeInfo("Merged " + originalSize + " segments to " + merged.size());
			}
			LocusSet<Segment> mergedSet = new LocusSet<Segment>(merged.toArray(new Segment[merged.size()]), true, log) {

				/**
				 * 
				 */
				private static final long serialVersionUID = 1L;

			};
			return mergedSet;
		}
	}

	public LocusSet<Segment> getStrictSegmentSet() {
		LocusSet<Segment> segSet = new LocusSet<Segment>(getStrictSegments(), true, log) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};
		return segSet;
	}

	public Segment[] getStrictSegments() {
		Segment[] segs = new Segment[loci.length];
		for (int i = 0; i < loci.length; i++) {
			segs[i] = new Segment(loci[i].getChr(), loci[i].getStart(), loci[i].getStop());
		}
		return segs;
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

	@SuppressWarnings("unchecked")
	public static <T extends Segment> LocusSet<T> combine(LocusSet<T> one, LocusSet<T> two, boolean sort, Logger log) {
		T[] combinedLoci = Array.concatAll(one.getLoci(), two.getLoci());
		LocusSet<T> combined = new LocusSet<T>(combinedLoci, sort, log) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};
		return combined;
	}

	/**
	 * @param array
	 *            writes according to the to string method of the class
	 * @param filename
	 * @param type
	 *            see {@link TO_STRING_TYPE}
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

	public void writeSerial(String fileName) {
		Files.writeSerial(this, fileName, true);
	}

	@SuppressWarnings("unchecked")
	public static LocusSet<CNVariant> readSerialCnvSet(String filename, Logger log) {
		return ((LocusSet<CNVariant>) Files.readSerial(filename, false, log, false, true));
	}

}