package org.genvisis.filesys;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Set;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

public class Segment implements Serializable {
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + chr;
		result = prime * result + start;
		result = prime * result + stop;
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Segment other = (Segment) obj;
		if (chr != other.chr)
			return false;
		if (start != other.start)
			return false;
		if (stop != other.stop)
			return false;
		return true;
	}

	public static final long serialVersionUID = 1L;

	protected byte chr;
	protected int start;
	protected int stop;
	
	public Segment() {}
	
	public Segment(byte chr, int start, int stop) {
		this.chr = chr;
		this.start = start;
		this.stop = stop;
	}

	public Segment(int start, int stop) {
		this.chr = 0;
		this.start = start;
		this.stop = stop;
	}

	public Segment(String ucscLocation) {
		int[] location = Positions.parseUCSClocation(ucscLocation);
		
		this.chr = (byte)location[0];
		this.start = location[1];
		this.stop = location[2];
	}
	
	public Segment(String chr, String start, String stop) {
		this.chr = Positions.chromosomeNumber(chr);
		this.start = Integer.parseInt(start);
		this.stop = Integer.parseInt(stop);
	}

	public Segment(String chr, int start, int stop) {
		this.chr = Positions.chromosomeNumber(chr);
		this.start = start;
		this.stop = stop;
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

	public String getUCSClocation() {
		return Positions.getUCSCformat(new int[] { getChr(), getStart(), getStop() });
	}
	
	public String getChromosomeUCSC() {
		return Positions.getChromosomeUCSC(chr, true);
	}

	public String getUCSCLink(String hgBuild) {
		return Positions.getUCSClinkInExcel(new int[] { getChr(), getStart(), getStop() }, hgBuild);
	}

	public int getSize() {
		return stop - start + 1;
	}

	public boolean equals(Segment seg) {
		return chr==seg.chr&&start==seg.start && stop==seg.stop;
	}

	public int amountOfOverlapInBasepairs(Segment seg) {
		if (chr == seg.chr || chr == -1 || seg.chr == -1) {
			if (start>=seg.start&&stop<=seg.stop) {
				return getSize();
			}
			if (seg.start>=start&&seg.stop<=stop) {
				return seg.getSize();
			}
			if (start>=seg.start&&start<=seg.stop) {
				return seg.stop-start+1;
			}
			if (stop>=seg.start&&stop<=seg.stop) {
				return stop-seg.start+1;
			}
		}
		
		return -1;
	}

	public boolean overlaps(Segment seg) {
		return amountOfOverlapInBasepairs(seg)>0;
	}

	public Segment getBufferedSegment(int buffer) {
		return new Segment(chr, start - buffer, stop + buffer);
	}

	public boolean significantOverlap(Segment seg) {
		return amountOfOverlapInBasepairs(seg) > Math.min(getSize(), seg.getSize()) / 2;
	}
	
	public boolean significantOverlap(Segment seg, boolean checkLarger) {
		int overlap = amountOfOverlapInBasepairs(seg);
		
		int mySize = getSize();
		int segSize = seg.getSize();

		/* 
		 * Get the threshold that overlap needs to pass to be significant:
		 * 		First check which is larger, then
		 * 		check argument flag,
		 * 		then set the threshold to the minimum of either 
		 * 			half the size of the larger CNV or
		 * 			the entirety of the smaller CNV 
		 * 				(if the smaller CNV is smaller than half the size of the larger CNV - otherwise we'd never call significance) 
		 */
		int threshold = mySize > segSize ? 
							(checkLarger ? Math.min(mySize / 2, segSize) : Math.min(segSize / 2, mySize))
						: 
							(checkLarger ? Math.min(segSize / 2, mySize ) : Math.min(mySize  / 2, segSize));
		
		return overlap > threshold;
	}
	
	public double overlapScore(Segment seg) {
		int overLap = amountOfOverlapInBasepairs(seg);
		return ((double) overLap / getSize());
	}
	
	

//	public boolean significantOverlapOld(Segment seg) {
//		return chr==seg.chr&&((start>=seg.start&&start<=seg.stop&&seg.stop-start>=getSize()/2)||(stop>=seg.start&&stop<=seg.stop&&stop-seg.start>=getSize()/2)||(seg.start>=start&&seg.start<=stop&&seg.stop-seg.start>=getSize()/2)||(seg.stop>=start&&seg.stop<=stop&&seg.stop-seg.start>=getSize()/2));
//	}

	public Segment merge(Segment seg) {
		if (chr != seg.chr) {
			System.err.println("Error - merging segments on different chromosomes");
		}
		return new Segment(chr, Math.min(start, seg.start), Math.max(stop, seg.stop));
	}
	
	/**
	 * @param seg
	 *            another segment that must overlap
	 * @param log
	 * @return the intersection of the two segments
	 */
	public Segment getIntersection(Segment seg, Logger log) {
		if (chr != seg.chr) {
			String error = "merging segments on different chromosomes";
			log.reportTimeError(error);
			throw new IllegalArgumentException(error);
		}
		if (!overlaps(seg)) {
			String error = "segments do not overlap";
			log.reportTimeError(error);
			throw new IllegalArgumentException(error);
		}
		return new Segment(chr, Math.max(start, seg.start), Math.min(stop, seg.stop));
	}

	/**
	 * @param segsToRemove
	 *            remove all of these segments from the current
	 * @param log
	 * @return the list of remaining dust
	 * 
	 *         NOTE: In the interest of speed, only pass segments that actually overlap
	 */
	public LocusSet<Segment> removeAll(Segment[] segsToRemove, Logger log) {
		if (segsToRemove == null || segsToRemove.length == 0) {
			LocusSet<Segment> original = new LocusSet<Segment>(new Segment[] { this }, true, log) {
				private static final long serialVersionUID = 1L;
			};
			return original;
		} else {

			LocusSet<Segment> removers = new LocusSet<Segment>(segsToRemove, true, log) {
				private static final long serialVersionUID = 1L;
			};

			Segment[] finalRemovers = removers.mergeOverlapping().getLoci();// deal with overlapping removal segments

			ArrayList<Segment> currentSegs = new ArrayList<Segment>();

			Segment[] removed = remove(finalRemovers[0], log);//seed removal
			if (removed != null) {
				for (int i = 0; i < removed.length; i++) {
					currentSegs.add(removed[i]);
				}
			}

			int totalBasePairsToRemove = 0;
			for (int i = 0; i < finalRemovers.length; i++) {
				totalBasePairsToRemove += Math.max(amountOfOverlapInBasepairs(finalRemovers[i]), 0);
			}

			int currentIndex = 1;
			while (currentIndex < finalRemovers.length) {// branch removal
				ArrayList<Segment> tmp = new ArrayList<Segment>();

				for (int i = 0; i < currentSegs.size(); i++) {
					Segment[] removedMore = currentSegs.get(i).remove(finalRemovers[currentIndex], log);
					if (removedMore != null) {
						for (int j = 0; j < removedMore.length; j++) {
							tmp.add(removedMore[j]);
						}
					}
				}
				currentSegs = new ArrayList<Segment>();
				currentSegs.addAll(tmp);
				currentIndex++;
			}

			LocusSet<Segment> finalSet = new LocusSet<Segment>(currentSegs.toArray(new Segment[currentSegs.size()]), true, log) {
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
				log.reportTimeError(error);
				throw new IllegalStateException(error);
			}

			for (int i = 0; i < segsToRemove.length; i++) {
				if (finalMergedSet.getOverLappingLoci(segsToRemove[i]) != null) {
					String error = "BUG: not all segments were properly removed";
					log.reportTimeError(error);
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
			cleaned = new Segment[] { this };
		} else {
			Segment intersection = getIntersection(seg, log);
			if (equals(intersection)) {
				cleaned = null;// removed all
			} else if (intersection.getStart() > getStart() && intersection.getStop() < getStop()) {// split
				Segment first = new Segment(getChr(), getStart(), intersection.getStart() - 1);
				Segment second = new Segment(getChr(), intersection.getStop() + 1, getStop());
				cleaned = new Segment[] { first, second };
			} else if (intersection.getStart() > getStart() && intersection.getStop() >= getStop()) {
				Segment head = new Segment(getChr(), getStart(), intersection.getStart() - 1);
				cleaned = new Segment[] { head };
			} else if (intersection.getStart() <= getStart() && intersection.getStop() < getStop()) {
				Segment tail = new Segment(getChr(), intersection.getStop() + 1, getStop());
				cleaned = new Segment[] { tail };
			} else {
				String error = "Un accounted for remove" + getUCSClocation() + " trying to remove " + seg.getUCSClocation();
				log.reportTimeError(error);
				throw new IllegalStateException(error);
			}
		}

		int numBpRemaining = 0;
		int bpShouldHaveBeenRemoved = Math.max(amountOfOverlapInBasepairs(seg), 0);
		if (cleaned != null) {
			for (int i = 0; i < cleaned.length; i++) {
				numBpRemaining += cleaned[i].getSize();
			}
		}
		int numBpRemoved = getSize() - numBpRemaining;

		if (numBpRemoved != bpShouldHaveBeenRemoved) {
			String error = "BUG: " + numBpRemoved + " base pairs were removed, but " + bpShouldHaveBeenRemoved + " should have been removed";
			error += "\nOriginal: " + getUCSClocation() + " Removed: " + seg.getUCSClocation();
			if (cleaned != null) {
				for (int i = 0; i < cleaned.length; i++) {
					error += "\n New: " + cleaned[i].getUCSClocation();
				}
			}
			log.reportTimeError(error);
			throw new IllegalStateException(error);
		}
		return cleaned;
	}
	public String toAnalysisString() {
		return Positions.getChromosomeUCSC(chr, true) + "\t" + start + "\t" + stop + "\t";
	}
	
	public String[] getHeader() {
		return new String[] { "chr", "start", "stop" };
	}

	public static boolean addIfAbsent(Segment seg, Vector<Segment> exons) {
		for (int i = 0; i<exons.size(); i++) {
			if (seg.equals(exons.elementAt(i))) {
				return false;
			}
        }
		exons.add(seg);
		return true;
	}
	
	
	public static Segment[] mergeOverlapsAndSortAllChromosomes(Segment[] segments,int buffer) {
		Hashtable<String, Vector<Segment>> splits = new Hashtable<String, Vector<Segment>>();
		for (int i = 0; i < segments.length; i++) {
			if (!splits.containsKey(segments[i].getChr() + "")) {
				splits.put(segments[i].getChr() + "", new Vector<Segment>());
			}
			splits.get(segments[i].getChr() + "").add(buffer > 0 ? segments[i].getBufferedSegment(buffer) : segments[i]);
		}
		ArrayList<Segment> merged = new ArrayList<Segment>();
		for (String chr : splits.keySet()) {
			Vector<Segment> tmp = splits.get(chr);
			mergeOverlapsAndSort(tmp);
			merged.addAll(tmp);
		}
		
		return sortSegments(merged.toArray(new Segment[merged.size()]));
	}

	// this method must be run separately for each chromosome
	public static void mergeOverlapsAndSort(Vector<Segment> segments) {
		byte chr;
		int[][] segBoundaries;
		int count, start, stop;
		
		if (segments.size() == 0) {
			return;
		}

		chr = segments.elementAt(0).getChr();
		for (int i = 0; i < segments.size(); i++) {
			if (segments.get(i).getChr() != chr) {
				System.err.println("Mismatched chromosmes for merging...");
				segments = null;
				return;
			}
		}
		segBoundaries = convertListToSortedBoundaries(segments);

		segments.clear();
		for (int i = 0; i<segBoundaries.length; i++) {
			if (segBoundaries[i][0] != -1) {
				count = 0;
				start = segBoundaries[i][0];
				stop = segBoundaries[i][1];
				while (i+count < segBoundaries.length && (segBoundaries[i+count][0] <= stop || segBoundaries[i+count][0] == -1)) {
					stop = Math.max(stop, segBoundaries[i+count][1]);
					segBoundaries[i+count][0] = -1;
					count++;
				}
				segments.add(new Segment(chr, start, stop));
			}
        }
	}

	public static void mergeOverlapsOld(Vector<Segment> segments) {
		boolean newlyAdded = true;
		Segment seg1, seg2;
		
		while (newlyAdded) {
			newlyAdded = false;
			for (int j = 0; j<segments.size(); j++) {
				for (int k = j+1; k<segments.size(); k++) {
					if (segments.elementAt(j).overlaps(segments.elementAt(k))) {
						seg2 = segments.remove(k);
						seg1 = segments.remove(j);
						segments.insertElementAt(seg1.merge(seg2), j);
						j = segments.size();
						k = segments.size();
		    			newlyAdded = true;
					}
                }
            }
		}
	}

	public static int[][] convertListToSortedBoundaries(Vector<Segment> segs) {
		int[][] segBoundaries = new int[segs.size()][2];
		int[] starts, keys;
		
		starts = new int[segs.size()];
		for (int i = 0; i<segs.size(); i++) {
			starts[i] = segs.elementAt(i).getStart();
		}
		keys = Sort.quicksort(starts);
		
		for (int i = 0; i<segs.size(); i++) {
			segBoundaries[i][0] = segs.elementAt(keys[i]).getStart();
			segBoundaries[i][1] = segs.elementAt(keys[i]).getStop();
        }

		return segBoundaries;
	}
	
	public static Segment[] toArray(Vector<Segment> setVec) {
		Segment[] list = new Segment[setVec.size()];
		
		for (int i = 0; i<setVec.size(); i++) {
			list[i] = setVec.elementAt(i);
        }
		
		return list;
	}
	
	public static int[] quicksort(Segment[] array) {
		return quicksort(array, true);
	}

	public static int[] quicksort(Segment[] array, boolean verbose) {
		byte[] chrs;
		int[] positions;

		chrs = new byte[array.length];
		positions = new int[array.length];
		for (int i = 0; i<array.length; i++) {
			chrs[i] = (byte)array[i].chr;
			positions[i] = array[i].start;
		}

		return Sort.orderTwoLayers(chrs, positions, verbose, new Logger());
	}

	public static Segment[] putInOrder(Segment[] array, int[] order) {
		Segment[] newArray;

		newArray = new Segment[array.length];
		for (int i = 0; i<order.length; i++) {
			newArray[i] = array[order[i]];
		}

		return newArray;
	}
	
	public static Segment[] sortSegments(Segment[] array) {
		return putInOrder(array, quicksort(array));
	}
	
	public static int binarySearchForOverlap(Segment seg, Segment[] orderedList) {
		int low, high, mid;
		
		low = 0;
		high = orderedList.length-1;
		while (low <= high) {
			mid = low + ((high - low) / 2);
			Segment inspect = orderedList[mid];
			if (inspect.overlaps(seg)) {
				return mid;
			} else if (seg.stop < inspect.start) {
				high = mid - 1;
			} else {
				low = mid + 1;
			}
		}

		return -1;
	}

	/**
	 * This function searches for a seed index using {@link Segment#binarySearchForOverlapChromosomeAware(Segment, Segment[])} and then scans up and down to
	 * <p>
	 * retrieve any other matches
	 * 
	 * @param seg
	 *            segment to search for
	 * @param orderedList
	 *            orderList of segments, in order by chromosome and then position
	 * @return index/indices of the overlapping segment, or null if not found
	 */
	public static int[] binarySearchForAllOverLappingIndices(Segment seg, Segment[] orderedList) {
		int index = binarySearchForOverlapChromosomeAware(seg, orderedList);
		if (index < 0) {
			return null;
		} else {
			ArrayList<Integer> overlaps = new ArrayList<Integer>();
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
			return Array.toIntegerArray(overlaps);
		}
	}
	
	/**
	 * Note that this function is chromosome aware and {@link Segment#binarySearchForOverlap(Segment, Segment[])} is not
	 * 
	 * @param seg
	 *            segment to search for
	 * @param orderedList
	 *            orderList of segments, in order by chromosome and then position
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
			} else if (seg.chr < orderedList[mid].chr || (seg.chr == orderedList[mid].chr && seg.start < orderedList[mid].start)) {
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
		high = orderedList.length-1;
		while (low<=high) {
			mid = low+(high-low)/2;
			if (orderedList[mid].getChr() == seg.getChr() && orderedList[mid].getStart() == seg.getStart()) {
				return mid;
			} else if (seg.chr < orderedList[mid].chr || (seg.chr == orderedList[mid].chr && seg.start < orderedList[mid].start)) {
				high = mid-1;
			} else {
				low = mid+1;
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
						if (orderedList[mid - 1].getChr() == seg.getChr() && orderedList[mid - 1].overlaps(seg)) {
							mid = mid - 1;
						} else {
							overlapping = false;
						}
					} else {
						overlapping = false;
					}
				}
				return mid;
			} else if (seg.chr < orderedList[mid].chr || (seg.chr == orderedList[mid].chr && seg.start < orderedList[mid].start)) {
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
		for (int i = 0; i < unorderedList.length; i++) {
			if (seg.equals(unorderedList[i])) {
				return true;
			}
		}
		return false;
	}
	
	public static Segment[] getSegments(String[] ucscLocations) {
		Segment[] segs;
		
		segs = new Segment[ucscLocations.length];
        for (int i = 0; i<segs.length; i++) {
        	segs[i] = new Segment(ucscLocations[i]);
        }
        
        return segs;
	}
	
	public static Segment[] loadUCSCregions(String filename, boolean ignoreFirstLine) {
		return loadUCSCregions(filename, 0, ignoreFirstLine, new Logger());
	}
	
	public static Segment[] loadUCSCregions(String filename, int column, boolean ignoreFirstLine, Logger log) {
		BufferedReader reader;
		Vector<Segment> v = new Vector<Segment>();

		try {
	        reader = new BufferedReader(new FileReader(filename));
	        if (ignoreFirstLine) {
	        	reader.readLine();
	        }
	        while (reader.ready()) {
				v.add(new Segment(reader.readLine().trim().split("[\\s]+")[column]));
	        }
	        reader.close();
        } catch (FileNotFoundException fnfe) {
        	log.reportError("Error: file \""+filename+"\" not found");
	        return null;
        } catch (IOException ioe) {
        	log.reportError("Error reading file \""+filename+"\"");
	        return null;
        }

		return Segment.toArray(v);
	}

	/**
	 * Removes duplicates based on {@link Segment#getUCSClocation()};
	 */
	public static Segment[] unique(Segment[] segments) {
		Hashtable<String, String> track = new Hashtable<String, String>();
		ArrayList<Segment> unique = new ArrayList<Segment>();
		for (int i = 0; i < segments.length; i++) {
			if (!track.containsKey(segments[i].getUCSClocation())) {
				track.put(segments[i].getUCSClocation(), segments[i].getUCSClocation());
				unique.add(segments[i]);
			}
		}
		return unique.toArray(new Segment[unique.size()]);
	}
	
	public static Segment[] loadRegions(String filename, int chrCol, int startCol, int stopCol, boolean ignoreFirstLine) {
		return loadRegions(filename, chrCol, startCol, stopCol, ignoreFirstLine ? 1 : 0, true, true, 0);
	}

	public static Segment[] loadRegions(String filename, int chrCol, int startCol, int stopCol, int skipNumLines, boolean sorted, boolean inclusiveStart, boolean inclusiveStop, int bpBuffer) {
		if (sorted) {
			return sortSegments(loadRegions(filename, chrCol, startCol, stopCol, skipNumLines, inclusiveStart, inclusiveStop, bpBuffer));
		} else {
			return loadRegions(filename, chrCol, startCol, stopCol, skipNumLines, inclusiveStart, inclusiveStop, bpBuffer);
		}
	}

	public static Segment[] loadRegions(String filename, int chrCol, int startCol, int stopCol, int skipNumLines, boolean inclusiveStart, boolean inclusiveStop, int bpBuffer) {
		BufferedReader reader;
		Vector<Segment> v = new Vector<Segment>();
		String[] line;

		try {
			reader = new BufferedReader(new FileReader(filename));
			for (int i = 0; i < skipNumLines; i++) {
				reader.readLine();
			}
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				v.add(new Segment(Positions.chromosomeNumber(line[chrCol]), (inclusiveStart ? Integer.parseInt(line[startCol]) : Integer.parseInt(line[startCol]) + 1) - bpBuffer, (inclusiveStop ? Integer.parseInt(line[stopCol]) : Integer.parseInt(line[stopCol]) - 1) + bpBuffer));
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
		Hashtable<Integer, Integer> track = new Hashtable<Integer, Integer>();
		int index = 0;
		for (int i = 0; i < segs.length; i++) {
			if (!track.containsKey((int) segs[i].getChr())) {
				track.put((int) segs[i].getChr(), index);
				index++;
			}
		}

		ArrayList<ArrayList<Segment>> tmp = new ArrayList<ArrayList<Segment>>();
		Set<Integer> indices = track.keySet();
		for (int i = 0; i < indices.size(); i++) {
			tmp.add( new ArrayList<Segment>());

		}
		for (int i = 0; i < segs.length; i++) {
			tmp.get(track.get((int) segs[i].getChr())).add(segs[i]);
		}

		Segment[][] parsed = new Segment[tmp.size()][];
		for (int i = 0; i < parsed.length; i++) {
			Segment[] chrSegs = tmp.get(i).toArray(new Segment[tmp.get(i).size()]);
			parsed[i] = chrSegs;
		}
		return parsed;
	}
	
	public static Vector<Segment> toVector(Segment[] segs) {
		Vector<Segment> v = new Vector<Segment>(segs.length);
		for (int i = 0; i < segs.length; i++) {
			v.add(segs[i]);
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
			writer = new PrintWriter(new FileWriter(firstFile+"_filteredOn_"+ext.removeDirectoryInfo(secondFile)+".out"));
			while (reader.ready()) {
				temp = reader.readLine();
				line = temp.trim().split("[\\s]+");
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
	 * Nested class to compare a segment to a list of other segments using {@link Segment#overlapScore}
	 * 
	 */
	public class SegmentCompare {

		private Segment[] compareSegments;
		private double maximumOverlapScore, scoreThreshold;
		private Segment maxScoreSegment;
		private int numOverlapingPastThreshold;
		private double avgOverlapScore;

		/**
		 * 
		 * Usage: Segment.SegmentCompare segmentCompare = segment.new SegmentCompare(compareSegments, scoreThreshold, log);
		 * <p>
		 * 
		 * @param compareSegments
		 *            Segment[] to compare this Segment to
		 * @param scoreThreshold
		 *            score threshold for the {@link Segment#overlapScore}. The number of segments exceeding this threshold will be stored.
		 * @param log
		 *            place holder, not currently used
		 */
		public SegmentCompare(Segment[] compareSegments, double scoreThreshold, Logger log) {
			this.compareSegments = compareSegments;
			this.scoreThreshold = scoreThreshold;
			this.maximumOverlapScore = 0;
			this.numOverlapingPastThreshold = 0;
			this.maxScoreSegment = new Segment((byte) 0, 0, 0);
			this.avgOverlapScore=0;
			// this.log = log;
		}

		public void compare() {
			for (int i = 0; i < compareSegments.length; i++) {
				double score = Segment.this.overlapScore(compareSegments[i]);
				if (score > maximumOverlapScore) {
					maximumOverlapScore = score;
					maxScoreSegment = compareSegments[i];
				}
				if (score > scoreThreshold) {
					numOverlapingPastThreshold++;
					avgOverlapScore += score;
				}
			}
			avgOverlapScore = (double) avgOverlapScore / numOverlapingPastThreshold;
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

		String usage = "\n" +
		"filesys.SegmentLists requires 0-1 arguments\n" +
		"   (1) first .bed filename (i.e. firstFile=onTarget.bed (default))\n" + 
		"   (2) second .bed filename (i.e. secondFile=genesOfInterest.bed (default))\n" + 
		"   (3) find segments in first that overlap any segment in second (i.e. -firstInSecond (not the default))\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("firstFile=")) {
				firstFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("secondFile=")) {
				secondFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-firstInSecond")) {
				firstInSecond = true;
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
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