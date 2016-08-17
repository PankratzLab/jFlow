package org.genvisis.filesys;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import org.genvisis.cnv.var.MosaicRegion;
import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.SerializedFiles;

import com.google.common.primitives.Bytes;
import com.google.common.primitives.Ints;

public abstract class LocusSet<T extends Segment> implements Serializable {
  public enum TO_STRING_TYPE {
    /**
     * calls the {@link Segment#getUCSClocation()}, or any overide
     */
    UCSC,
    /**
     * calls the {@link Segment#toString()}, or any override
     */
    REGULAR;
  }

  /**
   * 
   */
  protected static final long serialVersionUID = 1L;

  @SuppressWarnings("unchecked")
  public static <T extends Segment> LocusSet<T> combine(LocusSet<T> one, LocusSet<T> two,
      boolean sort, Logger log) {
    T[] combinedLoci = Array.concatAll(one.getLoci(), two.getLoci());
    LocusSet<T> combined = new LocusSet<T>(combinedLoci, sort, log) {

      /**
       * 
       */
      private static final long serialVersionUID = 1L;

    };
    return combined;
  }

  public static LocusSet<Segment> loadSegmentSetFromFile(String file, int chrCol, int startCol,
      int stopCol, int skipNumLines, boolean inclusiveStart, boolean inclusiveStop, int bpBuffer,
      Logger log) {
    Segment[] segs = Segment.loadRegions(file, chrCol, startCol, stopCol, skipNumLines,
        inclusiveStart, inclusiveStop, bpBuffer);
    LocusSet<Segment> lSet = new LocusSet<Segment>(segs, true, log) {

      /**
       * 
       */
      private static final long serialVersionUID = 1L;

    };
    return lSet;
  }

  @SuppressWarnings("unchecked")
  public static LocusSet<CNVariant> readSerialCnvSet(String filename, Logger log) {
    return ((LocusSet<CNVariant>) SerializedFiles.readSerial(filename, false, log, false, true));
  }

  @SuppressWarnings("unchecked")
  public static LocusSet<MosaicRegion> readSerialMRSet(String filename, Logger log) {
    return ((LocusSet<MosaicRegion>) SerializedFiles.readSerial(filename, false, log, false, true));
  }

  private T[] loci;

  private boolean sorted;

  private final Logger log;

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

  // /**
  // * Warning, uses equals to remove
  // */
  // public LocusSet<T> getUnique() {
  // throw new IllegalStateException("Not ready yet");
  // if (loci == null || loci.length < 1) {
  // return this;
  // } else {
  // ArrayList<T> tmp = new ArrayList<T>();
  // T[] tmpSearch = loci;
  // Hashtable<String, String> used = new Hashtable<String, String>();
  // for (int i = 0; i < tmpSearch.length; i++) {
  // int[] indices = getOverlappingIndices(tmpSearch[i]);
  // if (indices == null || indices.length == 0) {
  // tmp.add(tmpSearch[i]);
  // }
  // }
  //
  // LocusSet<T> returnSet = new LocusSet<T>(tmp, true, log) {
  //
  // /**
  // *
  // */
  // private static final long serialVersionUID = 1L;
  //
  // };
  // return returnSet;
  // }
  // };

  public LocusSet(final T[] loci, boolean sort, final Logger log) {
    super();
    this.loci = loci;
    this.sorted = false;
    this.log = log;
    if (sort) {
      sort();
    }
  }

  public void addAll(ArrayList<T> toAdd) {
    for (T element : loci) {
      toAdd.add(element);
    }
  }

  public LocusSet<T> autosomal(boolean sort, Logger log) {
    ArrayList<T> auto = new ArrayList<T>();
    for (int i = 0; i < getLoci().length; i++) {
      if (getLoci()[i].getChr() < 23) {
        auto.add(getLoci()[i]);
      }
    }
    if (auto.size() < 1) {
      throw new IllegalArgumentException("no autosomals T found");
    }
    LocusSet<T> aut = new LocusSet<T>(auto, sort, log) {

      /**
       * 
       */
      private static final long serialVersionUID = 1L;

    };
    return aut;
  }

  /**
   * @return simply the sum of the segments length, should be a merged set if you want the unique
   *         coverage
   */
  public long getBpCovered() {
    long sum = 0;
    for (T element : loci) {
      sum += element.getSize();
    }
    return sum;
  }

  public LocusSet<Segment> getBufferedSegmentSet(int bpBuffer) {
    ArrayList<Segment> buffered = new ArrayList<Segment>();
    for (T element : loci) {
      buffered.add(element.getBufferedSegment(bpBuffer));
    }

    LocusSet<Segment> bufSet =
        new LocusSet<Segment>(buffered.toArray(new Segment[buffered.size()]), true, log) {

          /**
           * 
           */
          private static final long serialVersionUID = 1L;

        };
    return bufSet;
  }

  public int[] getExactMatch(final Segment seg) {
    int[] overlaps = getOverlappingIndices(seg);
    ArrayList<Integer> exacts = new ArrayList<Integer>();
    if (overlaps == null || overlaps.length == 0) {
      return null;
    } else {
      for (int overlap : overlaps) {
        if (loci[overlap].equals(seg)) {
          exacts.add(overlap);
        }
      }
      return Ints.toArray(exacts);
    }
  }

  public T[] getLoci() {
    return loci;
  }

  public int[] getOverlappingIndices(final Segment seg) {
    if (!sorted) {
      log.reportTimeError(
          "Internal error: must sort internal segment array prior to overlap search");
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

  /**
   * @return a sorted int[][] of all the starts and stops of these loci
   */
  public int[][] getStartsAndStopsByChromosome() {
    Hashtable<Byte, ArrayList<Integer>> tracks = new Hashtable<Byte, ArrayList<Integer>>();
    ArrayList<Byte> uniqueChrs = new ArrayList<Byte>();
    for (int i = 0; i < loci.length; i++) {
      if (!tracks.containsKey(loci[i].getChr())) {
        tracks.put(loci[i].getChr(), new ArrayList<Integer>());
        uniqueChrs.add(loci[i].getChr());
      }
      tracks.get(loci[i].getChr()).add(loci[i].getStart());
      tracks.get(loci[i].getChr()).add(loci[i].getStop());
    }
    byte[] sortedChr = Bytes.toArray(uniqueChrs);
    Arrays.sort(sortedChr);
    int[][] startsStopsByChr = new int[27][0];

    for (int i = 0; i < sortedChr.length; i++) {
      int[] tmp = Ints.toArray(tracks.get(sortedChr[i]));
      Arrays.sort(tmp);
      startsStopsByChr[sortedChr[i]] = tmp;
    }
    return startsStopsByChr;
  }

  public Segment[] getStrictSegments() {
    Segment[] segs = new Segment[loci.length];
    for (int i = 0; i < loci.length; i++) {
      segs[i] = new Segment(loci[i].getChr(), loci[i].getStart(), loci[i].getStop());
    }
    return segs;
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

  public LocusSet<Segment> mergeOverlapping(boolean verbose) {// TODO, use clone and set positions
                                                              // instead to get same type returned

    if (!sorted) {
      log.reportTimeError("Internal error: must sort internal segment array prior to merge");
      return null;
    } else {
      byte currentChr = -1;
      ArrayList<Segment> merged = new ArrayList<Segment>();
      Vector<Segment> tmp = new Vector<Segment>();
      int originalSize = loci.length;
      for (T element : loci) {
        if (element.getChr() != currentChr) {
          if (tmp.size() > 0) {
            Segment.mergeOverlapsAndSort(tmp);
            for (int j = 0; j < tmp.size(); j++) {
              merged.add(tmp.get(j));
            }
            tmp.clear();
          }
        }
        currentChr = element.getChr();
        tmp.add(element);

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
      LocusSet<Segment> mergedSet =
          new LocusSet<Segment>(merged.toArray(new Segment[merged.size()]), true, log) {

            /**
             * 
             */
            private static final long serialVersionUID = 1L;

          };
      return mergedSet;
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

  /**
   * @param setToRemove remove the bases represented by this set from the current loci
   * @param bpBuffer the buffer to extend on either side of the removing set's loci
   * @return
   * 
   * 
   */
  public <E extends Segment> LocusSet<Segment> removeThese(final LocusSet<E> setToRemove,
      int bpBuffer) {
    ArrayList<Segment> newLoci = new ArrayList<Segment>();
    LocusSet<Segment> operateSet = setToRemove.getStrictSegmentSet();
    if (bpBuffer > 0) {
      operateSet = setToRemove.getBufferedSegmentSet(bpBuffer);
    }
    for (T element : loci) {
      Segment[] overlaps = operateSet.getOverLappingLoci(element);// exons in the bin
      if (overlaps == null || overlaps.length == 0) {
        newLoci.add(element);
      } else {
        LocusSet<Segment> overlapsRemoved = element.removeAll(overlaps, log);
        for (int j = 0; j < overlapsRemoved.getLoci().length; j++) {
          newLoci.add(overlapsRemoved.getLoci()[j]);
        }
      }
    }
    LocusSet<Segment> toReturn =
        new LocusSet<Segment>(newLoci.toArray(new Segment[newLoci.size()]), true, log) {
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

  public void sort() {
    loci = putInOrder(loci, Segment.quicksort(loci, false));
    sorted = true;
  }

  public boolean verifyPositiveLength() {
    for (T element : loci) {
      if (element.getSize() < 1) {
        return false;
      }
    }
    return true;
  }

  /**
   * @param array writes according to the to string method of the class
   * @param filename
   * @param type see {@link TO_STRING_TYPE}
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

  /**
   * Useful for converting bed files between b37 and hg19, etc
   */
  public boolean writeSegmentRegions(String filename, boolean numericChrs, Logger log) {
    log.reportTimeInfo("Writing " + loci.length + " loci to " + filename);
    boolean written = true;
    try {
      PrintWriter writer = new PrintWriter(new FileWriter(filename));
      for (T seg : loci) {
        writer.println(Positions.getChromosomeUCSC(seg.getChr(), !numericChrs) + "\t"
            + seg.getStart() + "\t" + seg.getStop());
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + filename);
      log.reportException(e);
      written = false;
    }
    return written;
  }

  public void writeSerial(String fileName) {
    SerializedFiles.writeSerial(this, fileName, true);
  }

}
