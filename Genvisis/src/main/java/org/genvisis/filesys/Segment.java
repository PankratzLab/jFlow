package org.genvisis.filesys;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.Vector;
import org.genvisis.common.Files;
import org.genvisis.common.GenomicPosition;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import com.google.common.primitives.Ints;

public class Segment implements Serializable, Comparable<Segment> {

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + getChr();
    result = prime * result + getStart();
    result = prime * result + getStop();
    return result;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    if (getClass() != obj.getClass()) {
      return false;
    }
    Segment other = (Segment) obj;
    if (getChr() != other.getChr()) {
      return false;
    }
    if (getStart() != other.getStart()) {
      return false;
    }
    if (getStop() != other.getStop()) {
      return false;
    }
    return true;
  }

  public static final long serialVersionUID = 1L;

  private final byte chr;
  private final int start;
  private final int stop;

  public Segment(byte chr, int start, int stop) {
    this.chr = chr;
    this.start = start;
    this.stop = stop;
  }

  /**
   * Construct a {@link Segment} for a single base
   * 
   * @param chr
   * @param position
   */
  public Segment(byte chr, int position) {
    this(chr, position, position);
  }

  /**
   * Construct a {@link Segment} for a single {@link GenomicPosition}
   * 
   * @param genomicPosition
   */
  public Segment(GenomicPosition genomicPosition) {
    this(genomicPosition.getChr(), genomicPosition.getPosition());
  }

  public Segment(int start, int stop) {
    this((byte) 0, start, stop);
  }

  public Segment(String ucscLocation) {
    this(Positions.parseUCSClocation(ucscLocation));
  }

  private Segment(int[] ucscLocation) {
    this((byte) ucscLocation[0], ucscLocation[1], ucscLocation[2]);
  }

  public Segment(String chr, String start, String stop) {
    this(Positions.chromosomeNumber(chr), Integer.parseInt(start), Integer.parseInt(stop));
  }

  public Segment(String chr, int start, int stop) {
    this(Positions.chromosomeNumber(chr), start, stop);
  }

  public byte getChr() {
    return chr;
  }

  public int getStart() {
    return start;
  }

  public int getStop() {
    return stop;
  }

  /**
   * @return difference between stop and start positions
   */
  public int getLength() {
    return getStop() - getStart();
  }

  public String getUCSClocation() {
    return Positions.getUCSCformat(new int[] {getChr(), getStart(), getStop()});
  }

  public String getChromosomeUCSC() {
    return Positions.getChromosomeUCSC(getChr(), true);
  }

  public String getUCSCLink(String hgBuild) {
    return Positions.getUCSClinkInExcel(new int[] {getChr(), getStart(), getStop()}, hgBuild);
  }

  public int getSize() {
    return getStop() - getStart() + 1;
  }

  public boolean equals(Segment seg) {
    return chr == seg.chr && start == seg.start && stop == seg.stop;
  }

  public int amountOfOverlapInBasepairs(Segment seg) {
    if (getChr() == seg.getChr() || getChr() == -1 || seg.getChr() == -1) {
      if (getStart() >= seg.getStart() && getStop() <= seg.getStop()) {
        return getSize();
      }
      if (seg.getStart() >= getStart() && seg.getStop() <= getStop()) {
        return seg.getSize();
      }
      if (getStart() >= seg.getStart() && getStart() <= seg.getStop()) {
        return seg.getStop() - getStart() + 1;
      }
      if (getStop() >= seg.getStart() && getStop() <= seg.getStop()) {
        return getStop() - seg.getStart() + 1;
      }
    }

    return -1;
  }

  public boolean overlaps(Segment seg) {
    return amountOfOverlapInBasepairs(seg) > 0;
  }

  public Segment getBufferedSegment(int buffer) {
    return new Segment(getChr(), getStart() - buffer, getStop() + buffer);
  }

  /**
   * @return true if this and target Segment contain more overlapping positions than half the
   *         smaller Segment's length (e.g., given a segment of length 10 and length 3, if 2
   *         positions overlap this method will return true).
   */
  public boolean significantOverlap(Segment seg) {
    return significantOverlap(seg, false);
  }

  /**
   * @param checkLarger If true, the significance threshold is stricter - using the smaller value of
   *          [larger segment's length / 2] and [smaller segment's length]. If false, significance
   *          is tested against 1/2 the smaller segment's length.
   * @return true if the number of overlapping positions in this and the target segment exceeds the
   *         requested threshold
   */
  public boolean significantOverlap(Segment seg, boolean checkLarger) {
    int smallSize = Math.min(getSize(), seg.getSize());
    int bigSize = Math.max(getSize(), seg.getSize()) / 2;

    return significantOverlap(seg, checkLarger ? Math.min(bigSize, smallSize) : smallSize / 2);
  }

  /**
   * @return true if the number of overlapping positions between this and the target segment exceed
   *         half this segment's length (that is, we don't care about the target segment's length -
   *         if it's too short it can't be called as overlapping, and if it's very long we can still
   *         call significance).
   */
  public boolean significantOneWayOverlap(Segment seg) {
    return significantOverlap(seg, getSize() / 2);
  }

  /**
   * @return true if the {@link #amountOfOverlapInBasepairs(Segment)} exceeds the specified
   *         threshold
   */
  public boolean significantOverlap(Segment seg, int threshold) {
    return amountOfOverlapInBasepairs(seg) > threshold;
  }

  public double overlapScore(Segment seg) {
    int overLap = amountOfOverlapInBasepairs(seg);
    return ((double) overLap / getSize());
  }

  // public boolean significantOverlapOld(Segment seg) {
  // return
  // chr==seg.chr&&((start>=seg.start&&start<=seg.stop&&seg.stop-start>=getSize()/2)||(stop>=seg.start&&stop<=seg.stop&&stop-seg.start>=getSize()/2)||(seg.start>=start&&seg.start<=stop&&seg.stop-seg.start>=getSize()/2)||(seg.stop>=start&&seg.stop<=stop&&seg.stop-seg.start>=getSize()/2));
  // }

  public Segment merge(Segment seg) {
    if (getChr() != seg.getChr()) {
      System.err.println("Error - merging segments on different chromosomes");
    }
    return new Segment(getChr(), Math.min(getStart(), seg.getStart()),
                       Math.max(getStop(), seg.getStop()));
  }

  /**
   * @param seg another segment that must overlap
   * @param log
   * @return the intersection of the two segments
   */
  public Segment getIntersection(Segment seg, Logger log) {
    if (getChr() != seg.getChr()) {
      String error = "merging segments on different chromosomes";
      log.reportError(error);
      throw new IllegalArgumentException(error);
    }
    if (!overlaps(seg)) {
      String error = "segments do not overlap";
      log.reportError(error);
      throw new IllegalArgumentException(error);
    }
    return new Segment(getChr(), Math.max(getStart(), seg.getStart()),
                       Math.min(getStop(), seg.getStop()));
  }

  /**
   * @param segsToRemove remove all of these segments from the current
   * @param log
   * @return the list of remaining dust NOTE: In the interest of speed, only pass segments that
   *         actually overlap
   */
  public LocusSet<Segment> removeAll(Segment[] segsToRemove, Logger log) {
    if (segsToRemove == null || segsToRemove.length == 0) {
      LocusSet<Segment> original = new LocusSet<Segment>(new Segment[] {this}, true, log) {

        private static final long serialVersionUID = 1L;
      };
      return original;
    } else {

      LocusSet<Segment> removers = new LocusSet<Segment>(segsToRemove, true, log) {

        private static final long serialVersionUID = 1L;
      };

      Segment[] finalRemovers = removers.mergeOverlapping().getLoci();// deal with overlapping
                                                                      // removal segments

      ArrayList<Segment> currentSegs = new ArrayList<>();

      Segment[] removed = remove(finalRemovers[0], log);// seed removal
      if (removed != null) {
        for (Segment element : removed) {
          currentSegs.add(element);
        }
      }

      int totalBasePairsToRemove = 0;
      for (Segment finalRemover : finalRemovers) {
        totalBasePairsToRemove += Math.max(amountOfOverlapInBasepairs(finalRemover), 0);
      }

      int currentIndex = 1;
      while (currentIndex < finalRemovers.length) {// branch removal
        ArrayList<Segment> tmp = new ArrayList<>();

        for (int i = 0; i < currentSegs.size(); i++) {
          Segment[] removedMore = currentSegs.get(i).remove(finalRemovers[currentIndex], log);
          if (removedMore != null) {
            for (Segment element : removedMore) {
              tmp.add(element);
            }
          }
        }
        currentSegs = new ArrayList<>();
        currentSegs.addAll(tmp);
        currentIndex++;
      }

      LocusSet<Segment> finalSet = new LocusSet<Segment>(currentSegs.toArray(new Segment[currentSegs.size()]),
                                                         true, log) {

        /**
        * 
        */
        private static final long serialVersionUID = 1L;

      };
      int totalBpRemaining = 0;
      LocusSet<Segment> finalMergedSet = finalSet.mergeOverlapping();

      for (int i = 0; i < finalMergedSet.getLoci().length; i++) {
        totalBpRemaining += finalMergedSet.getLoci()[i].getSize();
      }
      if (getSize() - totalBasePairsToRemove != totalBpRemaining) {

        String error = "BUG: did not remove segments properly";
        log.reportError(error);
        throw new IllegalStateException(error);
      }

      for (Segment element : segsToRemove) {
        if (finalMergedSet.getOverLappingLoci(element) != null) {
          String error = "BUG: not all segments were properly removed";
          log.reportError(error);
          throw new IllegalStateException(error);
        }
      }

      return finalMergedSet;
    }
  }

  /**
   * **Warning, not really tested
   *
   * @param seg
   * @param log
   * @return
   */
  public Segment[] remove(Segment seg, Logger log) {
    Segment[] cleaned = null;
    if (!overlaps(seg)) {
      cleaned = new Segment[] {this};
    } else {
      Segment intersection = getIntersection(seg, log);
      if (equals(intersection)) {
        cleaned = null;// removed all
      } else if (intersection.getStart() > getStart() && intersection.getStop() < getStop()) {// split
        Segment first = new Segment(getChr(), getStart(), intersection.getStart() - 1);
        Segment second = new Segment(getChr(), intersection.getStop() + 1, getStop());
        cleaned = new Segment[] {first, second};
      } else if (intersection.getStart() > getStart() && intersection.getStop() >= getStop()) {
        Segment head = new Segment(getChr(), getStart(), intersection.getStart() - 1);
        cleaned = new Segment[] {head};
      } else if (intersection.getStart() <= getStart() && intersection.getStop() < getStop()) {
        Segment tail = new Segment(getChr(), intersection.getStop() + 1, getStop());
        cleaned = new Segment[] {tail};
      } else {
        String error = "Un accounted for remove" + getUCSClocation() + " trying to remove "
                       + seg.getUCSClocation();
        log.reportError(error);
        throw new IllegalStateException(error);
      }
    }

    int numBpRemaining = 0;
    int bpShouldHaveBeenRemoved = Math.max(amountOfOverlapInBasepairs(seg), 0);
    if (cleaned != null) {
      for (Segment element : cleaned) {
        numBpRemaining += element.getSize();
      }
    }
    int numBpRemoved = getSize() - numBpRemaining;

    if (numBpRemoved != bpShouldHaveBeenRemoved) {
      String error = "BUG: " + numBpRemoved + " base pairs were removed, but "
                     + bpShouldHaveBeenRemoved + " should have been removed";
      error += "\nOriginal: " + getUCSClocation() + " Removed: " + seg.getUCSClocation();
      if (cleaned != null) {
        for (Segment element : cleaned) {
          error += "\n New: " + element.getUCSClocation();
        }
      }
      log.reportError(error);
      throw new IllegalStateException(error);
    }
    return cleaned;
  }

  public String toAnalysisString() {
    return Positions.getChromosomeUCSC(getChr(), true) + "\t" + getStart() + "\t" + getStop()
           + "\t";
  }

  public String[] getHeader() {
    return new String[] {"chr", "start", "stop"};
  }

  @Override
  public int compareTo(Segment o) {
    int c = getChr() - o.getChr();
    if (c == 0) {
      c = Integer.compare(this.getStart(), o.getStart());
    }
    return c;
  }

  public static boolean addIfAbsent(Segment seg, List<Segment> exons) {
    for (int i = 0; i < exons.size(); i++) {
      if (seg.equals(exons.get(i))) {
        return false;
      }
    }
    exons.add(seg);
    return true;
  }

  public static Segment[] mergeOverlapsAndSortAllChromosomes(Segment[] segments, int buffer) {
    Hashtable<String, Vector<Segment>> splits = new Hashtable<>();
    for (int i = 0; i < segments.length; i++) {
      if (!splits.containsKey(segments[i].getChr() + "")) {
        splits.put(segments[i].getChr() + "", new Vector<Segment>());
      }
      splits.get(segments[i].getChr() + "")
            .add(buffer > 0 ? segments[i].getBufferedSegment(buffer) : segments[i]);
    }
    ArrayList<Segment> merged = new ArrayList<>();
    for (String chr : splits.keySet()) {
      Vector<Segment> tmp = splits.get(chr);
      mergeOverlapsAndSort(tmp);
      merged.addAll(tmp);
    }

    Segment[] sorted = merged.toArray(new Segment[merged.size()]);
    Arrays.sort(sorted);
    return sorted;
  }

  // this method must be run separately for each chromosome
  public static void mergeOverlapsAndSort(List<Segment> segments) {
    byte chr;
    int[][] segBoundaries;
    int count, start, stop;

    if (segments.size() == 0) {
      return;
    }

    chr = segments.get(0).getChr();
    for (int i = 0; i < segments.size(); i++) {
      if (segments.get(i).getChr() != chr) {
        System.err.println("Mismatched chromosmes for merging...");
        segments = null;
        return;
      }
    }
    segBoundaries = convertListToSortedBoundaries(segments);

    segments.clear();
    for (int i = 0; i < segBoundaries.length; i++) {
      if (segBoundaries[i][0] != -1) {
        count = 0;
        start = segBoundaries[i][0];
        stop = segBoundaries[i][1];
        while (i + count < segBoundaries.length
               && (segBoundaries[i + count][0] <= stop || segBoundaries[i + count][0] == -1)) {
          stop = Math.max(stop, segBoundaries[i + count][1]);
          segBoundaries[i + count][0] = -1;
          count++;
        }
        segments.add(new Segment(chr, start, stop));
      }
    }
  }

  public static void mergeOverlapsOld(List<Segment> segments) {
    boolean newlyAdded = true;
    Segment seg1, seg2;

    while (newlyAdded) {
      newlyAdded = false;
      for (int j = 0; j < segments.size(); j++) {
        for (int k = j + 1; k < segments.size(); k++) {
          if (segments.get(j).overlaps(segments.get(k))) {
            seg2 = segments.remove(k);
            seg1 = segments.remove(j);
            segments.add(j, seg1.merge(seg2));
            j = segments.size();
            k = segments.size();
            newlyAdded = true;
          }
        }
      }
    }
  }

  public static int[][] convertListToSortedBoundaries(List<Segment> segs) {
    int[][] segBoundaries = new int[segs.size()][2];

    Segment[] arr = segs.toArray(new Segment[segs.size()]);
    Arrays.sort(arr);

    for (int i = 0; i < segs.size(); i++) {
      segBoundaries[i][0] = arr[i].getStart();
      segBoundaries[i][1] = arr[i].getStop();
    }

    return segBoundaries;
  }

  public static Segment[] toArray(List<Segment> setVec) {
    Segment[] list = new Segment[setVec.size()];

    for (int i = 0; i < setVec.size(); i++) {
      list[i] = setVec.get(i);
    }

    return list;
  }

  public static int binarySearchForOverlap(Segment seg, Segment[] orderedList) {
    int low, high, mid;

    low = 0;
    high = orderedList.length - 1;
    while (low <= high) {
      mid = low + ((high - low) / 2);
      Segment inspect = orderedList[mid];
      if (inspect.overlaps(seg)) {
        return mid;
      } else if (seg.getStop() < inspect.getStart()) {
        high = mid - 1;
      } else {
        low = mid + 1;
      }
    }

    return -1;
  }

  /**
   * This function searches for a seed index using
   * {@link Segment#binarySearchForOverlapChromosomeAware(Segment, Segment[])} and then scans up and
   * down to
   * <p>
   * retrieve any other matches
   *
   * @param seg segment to search for
   * @param orderedList orderList of segments, in order by chromosome and then position
   * @return index/indices of the overlapping segment, or null if not found
   */
  public static int[] binarySearchForAllOverLappingIndices(Segment seg, Segment[] orderedList) {
    int index = binarySearchForOverlapChromosomeAware(seg, orderedList);
    if (index < 0) {
      return null;
    } else {
      ArrayList<Integer> overlaps = new ArrayList<>();
      overlaps.add(index);
      for (int i = index + 1; i < orderedList.length; i++) {
        if (seg.overlaps(orderedList[i])) {
          overlaps.add(i);
        } else {
          break;
        }
      }
      for (int i = index - 1; i >= 0; i--) {
        if (seg.overlaps(orderedList[i])) {
          overlaps.add(i);
        } else {
          break;
        }
      }
      return Ints.toArray(overlaps);
    }
  }

  /**
   * Note that this function is chromosome aware and
   * {@link Segment#binarySearchForOverlap(Segment, Segment[])} is not
   *
   * @param seg segment to search for
   * @param orderedList orderList of segments, in order by chromosome and then position
   * @return index of the overlapping segment, or -1 if not found
   */
  public static int binarySearchForOverlapChromosomeAware(Segment seg, Segment[] orderedList) {
    int low, high, mid;

    low = 0;
    high = orderedList.length - 1;
    while (low <= high) {
      mid = low + (high - low) / 2;
      if (orderedList[mid].overlaps(seg)) {
        return mid;
      } else if (seg.getChr() < orderedList[mid].getChr()
                 || (seg.getChr() == orderedList[mid].getChr()
                     && seg.getStart() < orderedList[mid].getStart())) {
        high = mid - 1;
      } else {
        low = mid + 1;
      }
    }
    return -1;
  }

  public static int binarySearchForStartPositions(Segment seg, Segment[] orderedList) {
    int low, high, mid;

    low = 0;
    high = orderedList.length - 1;
    while (low <= high) {
      mid = low + (high - low) / 2;
      if (orderedList[mid].getChr() == seg.getChr()
          && orderedList[mid].getStart() == seg.getStart()) {
        return mid;
      } else if (seg.getChr() < orderedList[mid].getChr()
                 || (seg.getChr() == orderedList[mid].getChr()
                     && seg.getStart() < orderedList[mid].getStart())) {
        high = mid - 1;
      } else {
        low = mid + 1;
      }
    }

    return -1;
  }

  public static int binarySearchForOverlapMinIndex(Segment seg, Segment[] orderedList) {
    int low, high, mid;
    boolean overlapping = false;
    low = 0;
    high = orderedList.length - 1;
    while (low <= high) {
      mid = low + (high - low) / 2;
      if (orderedList[mid].getChr() == seg.getChr() && orderedList[mid].overlaps(seg)) {
        overlapping = true;
        // scan for minimum index
        while (overlapping) {
          if ((mid - 1) >= 0) {
            if (orderedList[mid - 1].getChr() == seg.getChr()
                && orderedList[mid - 1].overlaps(seg)) {
              mid = mid - 1;
            } else {
              overlapping = false;
            }
          } else {
            overlapping = false;
          }
        }
        return mid;
      } else if (seg.getChr() < orderedList[mid].getChr()
                 || (seg.getChr() == orderedList[mid].getChr()
                     && seg.getStart() < orderedList[mid].getStart())) {
        high = mid - 1;
      } else {
        low = mid + 1;
      }
    }
    return -1;
  }

  public static boolean overlapsAny(Segment seg, Segment[] orderedList) {
    return binarySearchForOverlap(seg, orderedList) >= 0;
  }

  public static boolean contains(Segment seg, Segment[] unorderedList) {
    for (Segment element : unorderedList) {
      if (seg.equals(element)) {
        return true;
      }
    }
    return false;
  }

  public static Segment[] getSegments(String[] ucscLocations) {
    Segment[] segs;

    segs = new Segment[ucscLocations.length];
    for (int i = 0; i < segs.length; i++) {
      segs[i] = new Segment(ucscLocations[i]);
    }

    return segs;
  }

  public static Segment[] loadUCSCregions(String filename, boolean ignoreFirstLine) {
    return loadUCSCregions(filename, 0, ignoreFirstLine, new Logger());
  }

  public static Segment[] loadUCSCregions(String filename, int column, boolean ignoreFirstLine,
                                          Logger log) {
    BufferedReader reader;
    Vector<Segment> v = new Vector<>();

    try {
      reader = new BufferedReader(new FileReader(filename));
      if (ignoreFirstLine) {
        reader.readLine();
      }
      while (reader.ready()) {
        v.add(new Segment(reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE)[column]));
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + filename + "\" not found");
      return null;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + filename + "\"");
      return null;
    }

    return Segment.toArray(v);
  }

  /**
   * Removes duplicates based on {@link Segment#getUCSClocation()};
   */
  public static Segment[] unique(Segment[] segments) {
    Hashtable<String, String> track = new Hashtable<>();
    ArrayList<Segment> unique = new ArrayList<>();
    for (int i = 0; i < segments.length; i++) {
      if (!track.containsKey(segments[i].getUCSClocation())) {
        track.put(segments[i].getUCSClocation(), segments[i].getUCSClocation());
        unique.add(segments[i]);
      }
    }
    return unique.toArray(new Segment[unique.size()]);
  }

  public static Segment[] loadRegions(String filename, int chrCol, int startCol, int stopCol,
                                      boolean ignoreFirstLine) {
    return loadRegions(filename, chrCol, startCol, stopCol, ignoreFirstLine ? 1 : 0, true, true, 0);
  }

  public static Segment[] loadRegions(String filename, int chrCol, int startCol, int stopCol,
                                      int skipNumLines, boolean sorted, boolean inclusiveStart,
                                      boolean inclusiveStop, int bpBuffer) {
    Segment[] regions = loadRegions(filename, chrCol, startCol, stopCol, skipNumLines,
                                    inclusiveStart, inclusiveStop, bpBuffer);
    if (sorted) {
      Arrays.sort(regions);
    }

    return regions;
  }

  public static Segment[] loadRegions(String filename, int chrCol, int startCol, int stopCol,
                                      int skipNumLines, boolean inclusiveStart,
                                      boolean inclusiveStop, int bpBuffer) {
    BufferedReader reader;
    Vector<Segment> v = new Vector<>();
    String[] line;

    try {
      reader = new BufferedReader(new FileReader(filename));
      for (int i = 0; i < skipNumLines; i++) {
        reader.readLine();
      }
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        v.add(new Segment(Positions.chromosomeNumber(line[chrCol]),
                          (inclusiveStart ? Integer.parseInt(line[startCol])
                                          : Integer.parseInt(line[startCol]) + 1) - bpBuffer,
                          (inclusiveStop ? Integer.parseInt(line[stopCol])
                                         : Integer.parseInt(line[stopCol]) - 1) + bpBuffer));
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }

    return Segment.toArray(v);
  }

  public static Segment[][] parseToChromosomeArrays(Segment[] segs, Logger log) {
    Hashtable<Integer, Integer> track = new Hashtable<>();
    int index = 0;
    for (int i = 0; i < segs.length; i++) {
      if (!track.containsKey((int) segs[i].getChr())) {
        track.put((int) segs[i].getChr(), index);
        index++;
      }
    }

    ArrayList<ArrayList<Segment>> tmp = new ArrayList<>();
    Set<Integer> indices = track.keySet();
    for (int i = 0; i < indices.size(); i++) {
      tmp.add(new ArrayList<Segment>());

    }
    for (Segment seg : segs) {
      tmp.get(track.get((int) seg.getChr())).add(seg);
    }

    Segment[][] parsed = new Segment[tmp.size()][];
    for (int i = 0; i < parsed.length; i++) {
      Segment[] chrSegs = tmp.get(i).toArray(new Segment[tmp.get(i).size()]);
      parsed[i] = chrSegs;
    }
    return parsed;
  }

  public static Vector<Segment> toVector(Segment[] segs) {
    Vector<Segment> v = new Vector<>(segs.length);
    for (Segment seg : segs) {
      v.add(seg);
    }
    return v;
  }

  public static void parseFirstInSecond(String firstFile, String secondFile, Logger log) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String temp;
    Segment[][] secondList;
    Segment seg;

    secondList = SegmentLists.parseSegmentList(secondFile, 0, 1, 2, false).getLists();

    try {
      reader = Files.getAppropriateReader(firstFile);
      writer = Files.openAppropriateWriter(firstFile + "_filteredOn_"
                                           + ext.removeDirectoryInfo(secondFile) + ".out");
      while (reader.ready()) {
        temp = reader.readLine();
        line = temp.trim().split(PSF.Regex.GREEDY_WHITESPACE);
        seg = new Segment(line[0], line[1], line[2]);
        if (secondList[seg.getChr()] != null && overlapsAny(seg, secondList[seg.getChr()])) {
          writer.println(temp);
        }
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + firstFile + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + firstFile + "\"");
      System.exit(2);
    }
  }

  /**
   * Nested class to compare a segment to a list of other segments using
   * {@link Segment#overlapScore}
   */
  public class SegmentCompare {

    private final Segment[] compareSegments;
    private double maximumOverlapScore;
    private final double scoreThreshold;
    private Segment maxScoreSegment;
    private int numOverlapingPastThreshold;
    private double avgOverlapScore;

    /**
     * Usage: Segment.SegmentCompare segmentCompare = segment.new SegmentCompare(compareSegments,
     * scoreThreshold, log);
     * <p>
     *
     * @param compareSegments Segment[] to compare this Segment to
     * @param scoreThreshold score threshold for the {@link Segment#overlapScore}. The number of
     *          segments exceeding this threshold will be stored.
     * @param log place holder, not currently used
     */
    public SegmentCompare(Segment[] compareSegments, double scoreThreshold, Logger log) {
      this.compareSegments = compareSegments;
      this.scoreThreshold = scoreThreshold;
      maximumOverlapScore = 0;
      numOverlapingPastThreshold = 0;
      maxScoreSegment = new Segment((byte) 0, 0, 0);
      avgOverlapScore = 0;
      // this.log = log;
    }

    public void compare() {
      for (Segment compareSegment : compareSegments) {
        double score = overlapScore(compareSegment);
        if (score > maximumOverlapScore) {
          maximumOverlapScore = score;
          maxScoreSegment = compareSegment;
        }
        if (score > scoreThreshold) {
          numOverlapingPastThreshold++;
          avgOverlapScore += score;
        }
      }
      avgOverlapScore = avgOverlapScore / numOverlapingPastThreshold;
    }

    public double getAvgOverlapScore() {
      return avgOverlapScore;
    }

    public Segment getMaxScoreSegment() {
      return maxScoreSegment;
    }

    public double getMaximumOverlapScore() {
      return maximumOverlapScore;
    }

    public int getNumOverlapingPastThreshold() {
      return numOverlapingPastThreshold;
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String logfile = null;
    Logger log;
    String firstFile = null;
    String secondFile = null;
    boolean firstInSecond = false;

    firstFile = "D:/Logan/Cosmic/S04380219_Regions_noHeader.bed";
    secondFile = "D:/Logan/Cosmic/cosmic_gene_positions.txt";
    firstInSecond = true;

    String usage = "\n" + "filesys.SegmentLists requires 0-1 arguments\n"
                   + "   (1) first .bed filename (i.e. firstFile=onTarget.bed (default))\n"
                   + "   (2) second .bed filename (i.e. secondFile=genesOfInterest.bed (default))\n"
                   + "   (3) find segments in first that overlap any segment in second (i.e. -firstInSecond (not the default))\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("firstFile=")) {
        firstFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("secondFile=")) {
        secondFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-firstInSecond")) {
        firstInSecond = true;
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      log = new Logger(logfile);
      if (firstInSecond) {
        parseFirstInSecond(firstFile, secondFile, log);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
