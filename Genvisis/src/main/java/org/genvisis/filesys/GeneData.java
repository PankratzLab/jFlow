package org.genvisis.filesys;

import java.io.Serializable;
import java.util.Vector;

import org.genvisis.common.*;

public class GeneData extends Segment implements Serializable {
	public static final long serialVersionUID = 1L;
	public static final int PLUS_STRAND = 1;
	public static final int MINUS_STRAND = 0;
	public static final int BOTH_STRANDS = 2;

	private String geneName;
	private String[] ncbiAssessionNumbers;
	private boolean positionFinalized;
	private boolean collapsed;
	private byte strand;
	private int[][] exonBoundaries;
	private byte multiLoc;

	public GeneData(String geneName, String[] ncbiAssessionNumbers, byte chr, boolean positionFinalized, byte strand, int startTranscription, int stop, int[][] exonBoundaries, byte multiLoc, boolean collapsedIsoformGene) {
		this.geneName = geneName;
		this.ncbiAssessionNumbers = ncbiAssessionNumbers;
		this.chr = chr;
		this.positionFinalized = positionFinalized;
		this.strand = strand;
		this.start = startTranscription;
		this.stop = stop;
		this.exonBoundaries = exonBoundaries;
		this.multiLoc = multiLoc;
		this.collapsed = collapsedIsoformGene;
	}
	
	public GeneData(String refGeneLine) {
		String[] line, starts, stops;
		
		line = refGeneLine.split("\t", -1);
		this.geneName = line[12];
		this.ncbiAssessionNumbers = new String[] {line[1]};
		
		if (!line[2].startsWith("chr")) {
			System.err.println("Error - don't know how to parse the chromosome from '"+line[2]+"'");
		} else {
			try {
				if (line[2].contains("_")) {
					this.chr = (byte)Positions.chromosomeNumber(line[2].substring(3, line[2].indexOf("_")));
					this.positionFinalized = false;
				} else {
					this.chr = (byte)Positions.chromosomeNumber(line[2].substring(3));
					this.positionFinalized = true;
				}
			} catch (Exception e) {
				System.err.println("Error - parsing chr ("+line[2]+") for accession "+line[1]);
				this.chr = -1;
				this.positionFinalized = false;
			}
		}
		if (line[3].equals("+")) {
			this.strand = PLUS_STRAND;
		} else if (line[3].equals("-")) {
			this.strand = MINUS_STRAND;
		} else {
			System.err.println("Error - unknown strand '"+line[3]+"' for accession '"+line[1]+"'");
		}

		this.start = Integer.parseInt(line[4]);
		this.stop = Integer.parseInt(line[5]);
		
		starts = line[9].trim().split(",", -1);
		stops = line[10].trim().split(",", -1);
		if (starts.length != stops.length || starts.length-1 != Integer.parseInt(line[8])) {
			System.err.println("Error - file format error: different number of start (n="+starts.length+") and stop (n="+stops.length+") exon boundaries for "+line[1]+"/"+line[12]);
		} else {
			this.exonBoundaries = new int[starts.length-1][2];
			for (int i = 0; i<starts.length-1; i++) {
				if (starts[i].equals("")) {
					System.err.println("Error - missing start position for exon "+(i+1)+" of "+line[1]);
				} else {
					this.exonBoundaries[i][0] = Integer.parseInt(starts[i]);
				}
				if (stops[i].equals("")) {
					System.err.println("Error - missing stop position for exon "+(i+1)+" of "+line[1]);
				} else {
					this.exonBoundaries[i][1] = Integer.parseInt(stops[i]);
				}
            }
		}
		
		this.collapsed = false;
	}
	
	public int getSize() {
		return stop-start+1;
	}

	public int amountOfOverlapInBasepairs(GeneData gene) {
		if (chr==gene.chr) {
			if (start>=gene.start&&stop<=gene.stop) {
				return getSize();
			}
			if (gene.start>=start&&gene.stop<=stop) {
				return gene.getSize();
			}
			if (start>=gene.start&&start<=gene.stop) {
				return gene.stop-start+1;
			}
			if (stop>=gene.start&&stop<=gene.stop) {
				return stop-gene.start+1;
			}
		}

		return -1;
	}

	public boolean overlaps(GeneData gene) {
		return amountOfOverlapInBasepairs(gene)>0;
	}

	public boolean significantOverlap(GeneData gene) {
		return amountOfOverlapInBasepairs(gene)>Math.min(getSize(), gene.getSize())/2;
	}

	public String getGeneName() {
		return geneName;
	}

	public String[] getNcbiAssessionNumbers() {
		return ncbiAssessionNumbers;
	}

	public boolean isFinalized() {
		return positionFinalized;
	}

	public byte getStrand() {
		return strand;
	}

	public int[][] getExonBoundaries() {
		return exonBoundaries;
	}

	public int getMultiLoc() {
		return multiLoc;
	}
	
	public boolean isCollapsedIsoforms() {
	    return collapsed;
	}
	
	// used to be just toArray, which overrode Segment.toArray, but Ant thinks this is a compile error for some reason
	public static GeneData[] toGeneDataArray(Vector<GeneData> setVec) {
		GeneData[] list = new GeneData[setVec.size()];
		
		for (int i = 0; i<setVec.size(); i++) {
			list[i] = setVec.elementAt(i);
        }
		
		return list;
	}
}
