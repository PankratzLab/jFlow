package org.genvisis.cnv.seq.manage;

import java.util.HashMap;
import java.util.Hashtable;
import java.util.Map;
import org.genvisis.cnv.analysis.BeastScore;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.seq.manage.BamImport.NGS_MARKER_TYPE;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.BamPile;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Segment;

/**
 * @author lane0212 Handles the data storage and normalization prior to conversion to {@link Sample}
 */
public class BamSample {

  /**
   * The scope of the normalization scaling factor derivation
   */
  public enum NORMALIZATON_METHOD {
    /**
     * Typical, but if whole chromosomal arms/large events exist, chromosome-based normalization can
     * become very skewed
     */
    CHROMOSOME,
    /**
     * Not sure if this is the preffered method yet, but might be a safer bet
     */
    GENOME;
  }

  private static final double MAX_MAPQ = 60;
  private static final double SCALE_FACTOR_NUM_READS = 1000000;
  private static final double MAD_FACTOR = 1.4826;
  private static final int MIN_NUM_MISMATCH = 1;// can set this later;
  private final String bamFile;
  private final String sampleName;
  private final BamPile[] bamPiles;
  private double[] rawDepth;
  private double[] normDepth;
  private double[] mapQs;
  private double[] percentWithMismatch;
  private final Project proj;
  private final NORMALIZATON_METHOD normMethod;

  /**
   * @param proj
   * @param bamFile the full path to the bamFile
   * @param bamPiles {@link BamPile} array to turn into a {@link BamSample}
   */
  BamSample(Project proj, String bamFile, BamPile[] bamPiles, NORMALIZATON_METHOD normMethod) {
    super();
    this.proj = proj;
    this.bamFile = bamFile;
    // bam file does not exist when processing cleaned sra runs
    sampleName = Files.exists(bamFile) ? BamOps.getSampleName(bamFile, proj.getLog())
                                       : ext.rootOf(bamFile);
    this.bamPiles = bamPiles;
    this.normMethod = normMethod;
    process();
  }

  String getSampleName() {
    return sampleName;
  }

  /**
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

    int getRpkmVal() {
      return rpkmVal;
    }

    void setRpkmVal(int rpkmVal) {
      this.rpkmVal = rpkmVal;
    }

    boolean[] getMask() {
      return mask;
    }

    NGS_MARKER_TYPE getType() {
      return type;
    }

  }

  private void process() {
    MarkerSetInfo markerSet = proj.getMarkerSet();
    String[] markerNames = markerSet.getMarkerNames();
    if (markerNames.length != bamPiles.length) {
      throw new IllegalArgumentException("Mismatched marker sizes, this is bad, was expecting "
                                         + markerNames.length + " and found " + bamPiles.length);
    }
    rawDepth = new double[bamPiles.length];
    mapQs = new double[bamPiles.length];
    percentWithMismatch = new double[bamPiles.length];

    BamPileParams[] params = new BamPileParams[NGS_MARKER_TYPE.values().length];
    for (int i = 0; i < params.length; i++) {
      params[i] = new BamPileParams(NGS_MARKER_TYPE.values()[i],
                                    ArrayUtils.booleanArray(markerNames.length, false));
    }
    proj.getLog().reportTimeInfo("Mapping queried segments");
    int[] traversalOrder = ArrayUtils.intArray(markerNames.length, -1);
    Map<String, Integer> lookup = new HashMap<>();
    for (int i = 0; i < bamPiles.length; i++) {
      lookup.put(bamPiles[i].getBin().getUCSClocation(), i);// if we store marker positions by start
                                                            // instead of midpoint, this will not be
                                                            // a problem
    }

    proj.getLog().reportTimeInfo("Computing rpkm mask");
    for (int i = 0; i < bamPiles.length; i++) {

      Segment markerSeg = new Segment(markerNames[i].split("\\|")[0]);
      if (!lookup.containsKey(markerSeg.getUCSClocation())) {
        proj.getLog().reportError("A major mismatching issue with " + markerSeg.getUCSClocation());
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

    proj.getLog().reportTimeInfo("Computing Normalized depths using method: " + normMethod);
    proj.getLog()
        .reportTimeInfo("Percent het will be reported at variant sites with alt depth greater than "
                        + MIN_NUM_MISMATCH);
    if (ArrayUtils.countIf(traversalOrder, -1) > 0) {
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
    int[][] normalizationIndices;
    switch (normMethod) {
      case CHROMOSOME:
        normalizationIndices = markerSet.getIndicesByChr();
        break;
      case GENOME:
        normalizationIndices = new int[][] {proj.getAutosomalMarkerIndices()};
        break;
      default:
        throw new IllegalArgumentException("Invalid normalization method: " + normMethod);
    }

    BeastScore beastScoreFirst = new BeastScore(ArrayUtils.toFloatArray(rawDepth),
                                                normalizationIndices, null, proj.getLog());
    int[] allProjectIndices = ArrayUtils.arrayOfIndices(proj.getMarkerNames().length);
    beastScoreFirst.setUse(params[0].getMask());
    float[] scaleMAD = getAppropriateScaling(beastScoreFirst, allProjectIndices, normMethod);

    for (int i = 1; i < params.length; i++) {
      BeastScore beastScoretmp = new BeastScore(scaleMAD, normalizationIndices, null,
                                                proj.getLog());
      beastScoretmp.setUse(params[i].getMask());
      scaleMAD = getAppropriateScaling(beastScoretmp, allProjectIndices, normMethod);
    }

    for (int i = 0; i < normalizationIndices.length; i++) {
      boolean error = false;
      for (int j = 0; j < normalizationIndices[i].length; j++) {
        int index = normalizationIndices[i][j];

        if (!Float.isFinite(scaleMAD[index])) {// should only happen if the MAD is NaN
          if (!error) {
            String warning = "Found invalid scale MAD depth for " + bamFile + ", bin "
                             + markerSet.getMarkerNames()[normalizationIndices[i][j]];
            warning += "Setting all of chr" + i + " to 0";
            warning += "This is usually caused by having zero passing reads for chr " + i;
            proj.getLog().reportTimeWarning(warning);
            error = true;
          }
        }
      }
      if (error) {
        for (int j = 0; j < normalizationIndices[i].length; j++) {
          scaleMAD[normalizationIndices[i][j]] = 0;
        }

      }
    }
    normDepth = ArrayUtils.scaleMinTo(ArrayUtils.toDoubleArray(scaleMAD), 1);
    for (int j = 0; j < normDepth.length; j++) {
      if (Double.isNaN(normDepth[j])) {
        String error = "Found invalid normalized depth for " + bamFile + ", bin "
                       + markerSet.getMarkerNames()[j];
        proj.getLog().reportError(error);
        if (markerSet.getMarkerNames()[j].contains(NGS_MARKER_TYPE.OFF_TARGET.getFlag())) {
          normDepth[j] = 1;
        } else {
          throw new IllegalStateException(error);
        }
      }
    }
  }

  Hashtable<String, Float> writeSample(long fingerprint) {
    Hashtable<String, Float> outliers = new Hashtable<>();
    byte[] genos = ArrayUtils.byteArray(bamPiles.length, (byte) 1);
    float[] blankLRRs = ArrayUtils.floatArray(bamPiles.length, 1);
    String sampleFile = proj.SAMPLE_DIRECTORY.getValue() + sampleName
                        + Sample.SAMPLE_FILE_EXTENSION;
    Sample sample = new Sample(sampleFile, fingerprint, ArrayUtils.toFloatArray(mapQs),
                               ArrayUtils.toFloatArray(normDepth),
                               ArrayUtils.toFloatArray(normDepth),
                               ArrayUtils.toFloatArray(percentWithMismatch), blankLRRs, genos,
                               genos, proj.getArrayType().getCanXYBeNegative());
    sample.saveToRandomAccessFile(sampleFile, outliers, sampleName);
    return outliers;
  }

  private static float[] getAppropriateScaling(BeastScore beastScore, int[] allProjectIndices,
                                               NORMALIZATON_METHOD normMethod) {
    // http://www.genomebiology.com/2014/15/12/550
    switch (normMethod) {
      case CHROMOSOME:
        return beastScore.getScaleMadRawData(MAD_FACTOR);
      case GENOME:
        return beastScore.getPropagatedScaling(MAD_FACTOR, allProjectIndices);
      default:
        throw new IllegalArgumentException("Invalid normalization method: " + normMethod);
    }
  }
}
