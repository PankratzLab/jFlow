package org.pankratzlab.common.filesys;

import java.io.Serializable;
import java.util.Vector;

public class GeneData extends Segment implements Serializable {

  public static final long serialVersionUID = 1L;
  public static final int PLUS_STRAND = 1;
  public static final int MINUS_STRAND = 0;
  public static final int BOTH_STRANDS = 2;

  private final String geneName;
  private final String[] ncbiAssessionNumbers;
  private boolean positionFinalized;
  private final boolean collapsed;
  private byte strand;
  private int[][] exonBoundaries;
  private byte multiLoc;

  public GeneData(String geneName, String[] ncbiAssessionNumbers, byte chr,
                  boolean positionFinalized, byte strand, int startTranscription, int stop,
                  int[][] exonBoundaries, byte multiLoc, boolean collapsedIsoformGene) {
    super(chr, startTranscription, stop);
    this.geneName = geneName;
    this.ncbiAssessionNumbers = ncbiAssessionNumbers;
    this.positionFinalized = positionFinalized;
    this.strand = strand;
    this.exonBoundaries = exonBoundaries;
    this.multiLoc = multiLoc;
    collapsed = collapsedIsoformGene;
  }

  public GeneData(String refGeneLine) {
    this(refGeneLine.split("\t", -1));

  }

  private GeneData(String[] refGeneLine) {
    super(parseRefGeneChr(refGeneLine), Integer.parseInt(refGeneLine[4]),
          Integer.parseInt(refGeneLine[5]));
    geneName = refGeneLine[12];
    ncbiAssessionNumbers = new String[] {refGeneLine[1]};

    positionFinalized = getChr() != -1 && !refGeneLine[2].contains("_");

    if (refGeneLine[3].equals("+")) {
      strand = PLUS_STRAND;
    } else if (refGeneLine[3].equals("-")) {
      strand = MINUS_STRAND;
    } else {
      System.err.println("Error - unknown strand '" + refGeneLine[3] + "' for accession '"
                         + refGeneLine[1] + "'");
    }

    String[] starts = refGeneLine[9].trim().split(",", -1);
    String[] stops = refGeneLine[10].trim().split(",", -1);
    if (starts.length != stops.length || starts.length - 1 != Integer.parseInt(refGeneLine[8])) {
      System.err.println("Error - file format error: different number of start (n=" + starts.length
                         + ") and stop (n=" + stops.length + ") exon boundaries for "
                         + refGeneLine[1] + "/" + refGeneLine[12]);
    } else {
      exonBoundaries = new int[starts.length - 1][2];
      for (int i = 0; i < starts.length - 1; i++) {
        if (starts[i].equals("")) {
          System.err.println("Error - missing start position for exon " + (i + 1) + " of "
                             + refGeneLine[1]);
        } else {
          exonBoundaries[i][0] = Integer.parseInt(starts[i]);
        }
        if (stops[i].equals("")) {
          System.err.println("Error - missing stop position for exon " + (i + 1) + " of "
                             + refGeneLine[1]);
        } else {
          exonBoundaries[i][1] = Integer.parseInt(stops[i]);
        }
      }
    }

    collapsed = false;
  }

  private static byte parseRefGeneChr(String[] line) {

    byte chr = -1;

    if (!line[2].startsWith("chr")) {
      System.err.println("Error - don't know how to parse the chromosome from '" + line[2] + "'");
    } else {
      try {
        if (line[2].contains("_")) {
          chr = Positions.chromosomeNumber(line[2].substring(3, line[2].indexOf("_")));
        } else {
          chr = Positions.chromosomeNumber(line[2].substring(3));
        }
      } catch (Exception e) {
        System.err.println("Error - parsing chr (" + line[2] + ") for accession " + line[1]);
      }
    }

    return chr;

  }

  @Override
  public int getSize() {
    return getStop() - getStart() + 1;
  }

  public int amountOfOverlapInBasepairs(GeneData gene) {
    if (getChr() == gene.getChr()) {
      if (getStart() >= gene.getStart() && getStop() <= gene.getStop()) {
        return getSize();
      }
      if (gene.getStart() >= getStart() && gene.getStop() <= getStop()) {
        return gene.getSize();
      }
      if (getStart() >= gene.getStart() && getStart() <= gene.getStop()) {
        return gene.getStop() - getStart() + 1;
      }
      if (getStop() >= gene.getStart() && getStop() <= gene.getStop()) {
        return getStop() - gene.getStart() + 1;
      }
    }

    return -1;
  }

  public boolean overlaps(GeneData gene) {
    return amountOfOverlapInBasepairs(gene) > 0;
  }

  public boolean significantOverlap(GeneData gene) {
    return amountOfOverlapInBasepairs(gene) > Math.min(getSize(), gene.getSize()) / 2;
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

  // used to be just toArray, which overrode Segment.toArray, but Ant thinks this is a compile error
  // for some reason
  public static GeneData[] toGeneDataArray(Vector<GeneData> setVec) {
    GeneData[] list = new GeneData[setVec.size()];

    for (int i = 0; i < setVec.size(); i++) {
      list[i] = setVec.elementAt(i);
    }

    return list;
  }
}
