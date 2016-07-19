package org.genvisis.filesys;

import java.io.PrintWriter;
import java.io.Serializable;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.filesys.PlainTextExport;

public class SegmentLists implements Serializable, PlainTextExport {
	public static final long serialVersionUID = 1L;
	
	private Segment[][] lists;
	
	public SegmentLists(Segment[][] lists) {
		this.lists = lists;
	}

	public Segment[][] getLists() {
		return lists;
	}

	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}
	
	@Override
	public void exportToText(String outputFile, Logger log) {
	    PrintWriter writer;
	    
	    writer = Files.getAppropriateWriter(outputFile);
	    writer.println("Chr\tStart\tStop");
	    for (Segment[] segList : lists) {
	        for (Segment seg : segList) {
	            writer.println(seg.toAnalysisString());
	        }
	    }
	    writer.flush();
	    writer.close();
	}
	
	public static SegmentLists parseUCSCSegmentList(String filename, boolean ignoreFirstLine) {
		Hashtable<String, Vector<Segment>> hash = new Hashtable<String,Vector<Segment>>();
		Vector<Segment> vSegs;
		Segment[][] lists;
		int[] chrs;
		Segment[] segs = null;
		
		segs = Segment.loadUCSCregions(filename, ignoreFirstLine);		

		for (int i = 0; i<segs.length; i++) {
			if (hash.containsKey(segs[i].getChr()+"")) {
				vSegs = hash.get(segs[i].getChr()+"");
			} else {
				hash.put(segs[i].getChr()+"", vSegs = new Vector<Segment>());
			}
			vSegs.add(new Segment(segs[i].getChr(), segs[i].getStart(), segs[i].getStop()));
        }
		chrs = Array.toIntArray(HashVec.getKeys(hash));
//		lists = new Segment[Array.max(chrs)+1][];
		lists = new Segment[27][];
		for (int i = 0; i<chrs.length; i++) {
			vSegs = hash.get(chrs[i]+"");
        	Segment.mergeOverlapsAndSort(vSegs);
        	lists[chrs[i]] = Segment.toArray(vSegs);
        }
		
		return new SegmentLists(lists);
	}
	
    /**
     * Overloaded function to make a {@link SegmentList} of {@link Segment} sorted by arranged by chr
     * Use getLists() function to get the list of segments as Segment[][]
     *
     * @param filename:        name of the file containing segments
     * @param chrCol:          column number of chr
     * @param startCol:        column number of starting position
     * @param stopCol:         column number of stopping position
     * @param ignoreFirstLine: boolean to decide whether to ignore first line or not
     * @return a {@link SegmentList} of {@link Segment} sorted by arranged by chr
     */
    public static SegmentLists parseSegmentList(String filename, int chrCol, int startCol, int stopCol,
                                                boolean ignoreFirstLine) {
        Hashtable<String, Vector<Segment>> hash = new Hashtable<String, Vector<Segment>>();
        Vector<Segment> vSegs;
        Segment[][] lists;
        int[] chrs;
        Segment[] segs = null;

        segs = Segment.loadRegions(filename, chrCol, startCol, stopCol, ignoreFirstLine);

        for (int i = 0; i < segs.length; i++) {
            if (hash.containsKey(segs[i].getChr() + "")) {
                vSegs = hash.get(segs[i].getChr() + "");
            } else {
                hash.put(segs[i].getChr() + "", vSegs = new Vector<Segment>());
            }
            vSegs.add(new Segment(segs[i].getChr(), segs[i].getStart(), segs[i].getStop()));
        }
        chrs = Array.toIntArray(HashVec.getKeys(hash));
        lists = new Segment[27][];
        for (int i = 0; i < chrs.length; i++) {
            vSegs = hash.get(chrs[i] + "");
            Segment.mergeOverlapsAndSort(vSegs);
            lists[chrs[i]] = Segment.toArray(vSegs);
        }

        return new SegmentLists(lists);
    }
	
	public static SegmentLists load(String filename, boolean jar) {
		return (SegmentLists)Files.readSerial(filename, jar, true);
	}
}
