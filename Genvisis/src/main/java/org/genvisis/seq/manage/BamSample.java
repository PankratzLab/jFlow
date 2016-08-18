package org.genvisis.seq.manage;

import java.util.Hashtable;

import org.genvisis.cnv.analysis.BeastScore;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.common.Array;
import org.genvisis.common.Numbers;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.BamImport.NGS_MARKER_TYPE;

/**
 * @author lane0212 Handles the data storage and normalization prior to conversion to {@link Sample}
 */
public class BamSample {
  private static final double MAX_MAPQ = 60;
  private static final double SCALE_FACTOR_NUM_READS = 1000000;
  private static final double MAD_FACTOR = 1.4826;
  private static final int MIN_NUM_MISMATCH = 5;// can set this later;
  private final String bamFile;
  private final String sampleName;
  private final BamPile[] bamPiles;
  private double[] rawDepth;
  private double[] normDepth;
  private double[] mapQs;
  private double[] percentWithMismatch;
  private final Project proj;

  public BamSample(Project proj, String bamFile, BamPile[] bamPiles) {
    super();
    this.proj = proj;
    this.bamFile = bamFile;
    sampleName = BamOps.getSampleName(bamFile);
    this.bamPiles = bamPiles;
    process();
  }

  public String getSampleName() {
    return sampleName;
  }

  /**
   * 
   * See http://www.partek.com/Tutorials/microarray/User_Guides/UnderstandingReads.pdf , page 3
   * 
   * @param numMappedReads
   * @param bin
   * @param numTotalMappedReads
   * @return
   */
  private static double computeRPKM(int numMappedReads, Segment bin, int numTotalMappedReads) {
    double data = (double) numMappedReads / bin.getSize();
    if (Double.isNaN(data)) {
      throw new IllegalArgumentException("Size and num mapped reads cannot be NaN, size cannot be 0"
                                         + bin.getUCSClocation());
    }
    double scale = numTotalMappedReads > 0 ? SCALE_FACTOR_NUM_READS / numTotalMappedReads : 0;
    return data * scale;
  }

  private static class BamPileParams {
    private final NGS_MARKER_TYPE type;
    private final boolean[] mask;
    private int rpkmVal;

    private BamPileParams(NGS_MARKER_TYPE type, boolean[] mask) {
      super();
      this.type = type;
      this.mask = mask;
      rpkmVal = 0;
    }

    public int getRpkmVal() {
      return rpkmVal;
    }

    public void setRpkmVal(int rpkmVal) {
      this.rpkmVal = rpkmVal;
    }

    public boolean[] getMask() {
      return mask;
    }

    public NGS_MARKER_TYPE getType() {
      return type;
    }

  }



  private void process() {
    MarkerSet markerSet = proj.getMarkerSet();
    String[] markerNames = markerSet.getMarkerNames();
    if (markerNames.length != bamPiles.length) {
      throw new IllegalArgumentException("Mismatched marker sizes, this is bad");
    }
    rawDepth = new double[bamPiles.length];
    mapQs = new double[bamPiles.length];
    percentWithMismatch = new double[bamPiles.length];

    BamPileParams[] params = new BamPileParams[NGS_MARKER_TYPE.values().length];
    for (int i = 0; i < params.length; i++) {
      params[i] = new BamPileParams(NGS_MARKER_TYPE.values()[i],
                                    Array.booleanArray(markerNames.length, false));
    }
    proj.getLog().reportTimeInfo("Mapping queried segments");
    int[] traversalOrder = Array.intArray(markerNames.length, -1);
    Hashtable<String, Integer> lookup = new Hashtable<String, Integer>();
    for (int i = 0; i < bamPiles.length; i++) {
      lookup.put(bamPiles[i].getBin().getUCSClocation(), i);// if we store marker positions by start
                                                            // instead of midpoint, this will not be
                                                            // a problem
    }

    proj.getLog().reportTimeInfo("Computing rpkm mask");
    for (int i = 0; i < bamPiles.length; i++) {

      Segment markerSeg = new Segment(markerNames[i].split("\\|")[0]);
      if (!lookup.containsKey(markerSeg.getUCSClocation())) {
        System.err.println(markerSeg.getUCSClocation());
        throw new IllegalArgumentException("A major mismatching issue");
      }
      int bamPileIndex = lookup.get(markerSeg.getUCSClocation());// if more than 1, identical so
                                                                 // does'nt matter
      traversalOrder[i] = bamPileIndex;
      BamPile currentPile = bamPiles[bamPileIndex];

      NGS_MARKER_TYPE current = NGS_MARKER_TYPE.getType(markerSet.getMarkerNames()[i]);
      for (int j = 0; j < params.length; j++) {
        if (current == params[j].getType()) {
          params[j].getMask()[i] = true;
          params[j].setRpkmVal(params[j].getRpkmVal() + currentPile.getNumOverlappingReads());
          break;
        }
      }
    }

    proj.getLog().reportTimeInfo("Computing Normalized depths");
    proj.getLog().reportTimeInfo(
                                 "Percent het will be reported at variant sites with alt depth greater than "
                                 + MIN_NUM_MISMATCH);
    if (Array.countIf(traversalOrder, -1) > 0) {
      throw new IllegalArgumentException("Not all indices accounted for");
    }

    for (int i = 0; i < traversalOrder.length; i++) {
      BamPile currentPile = bamPiles[traversalOrder[i]];

      NGS_MARKER_TYPE current = NGS_MARKER_TYPE.getType(markerSet.getMarkerNames()[i]);
      for (BamPileParams param : params) {
        if (current == param.getType()) {
          rawDepth[i] = computeRPKM(currentPile.getNumOverlappingReads(), currentPile.getBin(),
                                    param.getRpkmVal());
          break;
        }
      }

      mapQs[i] = Math.min(currentPile.getOverallAvgMapQ() / MAX_MAPQ, 1);
      if (Double.isNaN(rawDepth[i])) {
        String warning = "Found invalid scale raw depth for " + bamFile + ", bin "
                         + markerSet.getMarkerNames()[i];
        proj.getLog().reportTimeWarning(warning);
        throw new IllegalArgumentException(warning);
      }
      if (current == NGS_MARKER_TYPE.VARIANT_SITE) {
        double normBasesOverlap = currentPile.getNumBasesOverlap();
        double normBasesMiss = currentPile.getNumBasesWithMismatch();
        double percentMiss = 0;
        if (normBasesMiss > MIN_NUM_MISMATCH) {
          percentMiss = normBasesMiss / normBasesOverlap;
        }
        percentWithMismatch[i] = percentMiss;
      } else {
        percentWithMismatch[i] = 0;
      }
    }
    int[][] chrIndices = markerSet.getIndicesByChr();

    BeastScore beastScoreFirst =
        new BeastScore(Array.toFloatArray(rawDepth), chrIndices, null, proj.getLog());
    beastScoreFirst.setUse(params[0].getMask());
    float[] scaleMAD = beastScoreFirst.getScaleMadRawData(MAD_FACTOR);// http://www.genomebiology.com/2014/15/12/550

    for (int i = 1; i < params.length; i++) {
      BeastScore beastScoretmp = new BeastScore(scaleMAD, chrIndices, null, proj.getLog());
      beastScoretmp.setUse(params[i].getMask());
      scaleMAD = beastScoretmp.getScaleMadRawData(MAD_FACTOR);
    }

    for (int i = 0; i < chrIndices.length; i++) {
      boolean error = false;
      for (int j = 0; j < chrIndices[i].length; j++) {
        int index = chrIndices[i][j];

        if (!Numbers.isFinite(scaleMAD[index])) {// should only happen if the MAD is NaN
          if (!error) {
            String warning = "Found invalid scale MAD depth for " + bamFile + ", bin "
                             + markerSet.getMarkerNames()[chrIndices[i][j]];
            warning += "Setting all of chr" + i + " to 0";
            warning += "This is usually caused by having zero passing reads for chr " + i;
            proj.getLog().reportTimeWarning(warning);
            error = true;
          }
        }
      }
      if (error) {
        for (int j = 0; j < chrIndices[i].length; j++) {
          scaleMAD[chrIndices[i][j]] = 0;
        }

      }
    }
    normDepth = Array.scaleMinTo(Array.toDoubleArray(scaleMAD), 1);
    for (int j = 0; j < normDepth.length; j++) {
      if (Double.isNaN(normDepth[j])) {
        String error = "Found invalid normalized depth for " + bamFile + ", bin "
                       + markerSet.getMarkerNames()[j];
        proj.getLog().reportTimeError(error);
        if (markerSet.getMarkerNames()[j].contains(NGS_MARKER_TYPE.OFF_TARGET.getFlag())) {
          normDepth[j] = 1;
        } else {
          throw new IllegalStateException(error);
        }
      }
    }
  }

  public Hashtable<String, Float> writeSample(long fingerprint) {
    Hashtable<String, Float> outliers = new Hashtable<String, Float>();
    byte[] genos = Array.byteArray(bamPiles.length, (byte) 1);
    float[] blankLRRs = Array.floatArray(bamPiles.length, 1);
    String sampleFile =
        proj.SAMPLE_DIRECTORY.getValue() + sampleName + Sample.SAMPLE_FILE_EXTENSION;
    Sample sample =
        new Sample(sampleFile, fingerprint, Array.toFloatArray(mapQs),
                   Array.toFloatArray(normDepth), Array.toFloatArray(normDepth),
                   Array.toFloatArray(percentWithMismatch), blankLRRs, genos, genos, false);
    sample.saveToRandomAccessFile(sampleFile, outliers, sampleName);
    return outliers;
  }
}

//
// // // if(Array.min(normDepth)<0){
// // // System.err.println("less dan 0");
// // // System.exit(1);
// // // }
// double maxNorm = Array.max(normDepth);
// double minNorm = Array.min(normDepth);
// for (int i = 0; i < normDepth.length; i++) {
// normDepth[i] = (normDepth[i] - minNorm) / (maxNorm - minNorm);
// normDepth[i] += 1;
// // normDepth[i] *= 10;
// }
// System.out.println(Array.mean(normDepth));
// System.out.println(Array.min(normDepth));
// System.out.println(Array.max(normDepth))
