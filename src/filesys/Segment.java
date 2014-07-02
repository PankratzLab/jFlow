package filesys;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.Vector;

import common.Files;
import common.Logger;
import common.Positions;
import common.Sort;
import common.ext;

public class Segment implements Serializable {
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
		return "chr"+chr+":"+start+"-"+stop;
	}
	
	public String getUCSCLink() {
		return Positions.getUCSClinkInExcel(new int[] { getChr(), getStart(), getStop() });
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

	public boolean significantOverlap(Segment seg) {
		return amountOfOverlapInBasepairs(seg)>Math.min(getSize(), seg.getSize())/2;
	}

	public double overlapScore(Segment seg) {
		int overLap = amountOfOverlapInBasepairs(seg);
		return (((double) overLap / getSize()) + ((double) overLap / seg.getSize())) / 2;
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

	public static boolean addIfAbsent(Segment seg, Vector<Segment> exons) {
		for (int i = 0; i<exons.size(); i++) {
			if (seg.equals(exons.elementAt(i))) {
				return false;
			}
        }
		exons.add(seg);
		return true;
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
		byte[] chrs;
		int[] positions;

		chrs = new byte[array.length];
		positions = new int[array.length];
		for (int i = 0; i<array.length; i++) {
			chrs[i] = (byte)array[i].chr;
			positions[i] = array[i].start;
		}

		return Sort.orderTwoLayers(chrs, positions);
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
		while (low<=high) {
			mid = low+(high-low)/2;
			if (orderedList[mid].overlaps(seg)) {
				return mid;
			} else if (seg.stop < orderedList[mid].start) {
				high = mid-1;
			} else {
				low = mid+1;
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

	public static Segment[] loadRegions(String filename, int chrCol, int startCol, int stopCol, boolean ignoreFirstLine) {
		BufferedReader reader;
		Vector<Segment> v = new Vector<Segment>();
		String[] line;

		try {
	        reader = new BufferedReader(new FileReader(filename));
	        if (ignoreFirstLine) {
	        	reader.readLine();
	        }
	        while (reader.ready()) {
	        	line = reader.readLine().trim().split("[\\s]+");
				v.add(new Segment(Positions.chromosomeNumber(line[chrCol]), Integer.parseInt(line[startCol]), Integer.parseInt(line[stopCol])));
	        }
	        reader.close();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+filename+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+filename+"\"");
	        System.exit(2);
        }

		return Segment.toArray(v);
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
