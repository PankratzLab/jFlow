package org.genvisis.cnv.analysis;

import java.util.ArrayList;

import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.Transforms;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.Segment;

import com.google.common.primitives.Doubles;

/**
 * Class to compute beast scores similar to the beast algorithm
 * http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0086272. Values are always
 * inverse normalized with 5df, and scaled by the empirically derived SCALE_FACTOR_MAD. NaN data is
 * discouraged, but is ignored; Note: Since the input data is always inverse normalized, the MAD for
 * a region (such as a cnv) is the median of the (absolute value) inverse normalized data across
 * that region (scaled by the MAD of the superset i.e the chromosome) Note: one difference to the
 * score scheme is that a min/max length parameter has not been added. i.e(length<=min -> score=0,
 * length>=max -> height goes to min height. This may make long cnv calls have inflated scores.
 */
public class BeastScore {
  /**
   * The scale factor was empirically derived and scales all MAD values for consistent comparisons
   * across samples
   */
  public static final double SCALE_FACTOR_MAD = 0.134894516;
  /**
   * The alpha suggested by the beast folks
   */
  public static final float DEFAULT_ALPHA = 0.5f;
  private final float[] inputData;
  private float[] beastHeights;
  private float[] beastScores;
  private float[] invTransScaledStdev;
  private final int[][] indicesToChunk;
  private final int[][] indicesForScores;
  private int[] beastLengths;
  private boolean[] use;
  private final Logger log;

  /**
   * @param inputData Most often LRR values for a single sample
   * @param indicesToChunk Most often the indices of chromosomes, if null, all indices will be used
   *        as a group
   * @param indicesForScores Most often the indices of markers in cnvs for a particular sample, if
   *        null each indice will be scored
   * @param log You know, a log
   */
  public BeastScore(float[] inputData, int[][] indicesToChunk, int[][] indicesForScores,
                    Logger log) {
    super();
    this.inputData = inputData;
    this.indicesToChunk =
                        indicesToChunk == null ? new int[][] {Array.arrayOfIndices(inputData.length)}
                                               : indicesToChunk;
    this.indicesForScores =
                          indicesForScores == null ? new int[][] {Array.arrayOfIndices(inputData.length)}
                                                   : indicesForScores;
    beastHeights = new float[this.indicesForScores.length];
    beastScores = new float[this.indicesForScores.length];
    beastLengths = new int[this.indicesForScores.length];
    use = null;
    this.log = log;
  }

  public float[] getBeastScores() {
    return beastScores;
  }

  public float[] getBeastHeights() {
    return beastHeights;
  }

  public int[] getBeastLengths() {
    return beastLengths;
  }

  public float[] getInvTransScaledStdevChunks() {
    return invTransScaledStdev;
  }

  public void setUse(boolean[] useForMedian) {
    use = useForMedian;
  }

  public String getSummaryAt(int index) {
    if (index >= beastScores.length) {
      log.reportError("Error - requested a summary for an index that is too big");
    }
    return getBeastHeights()[index] + "\t" + getBeastLengths()[index] + "\t"
           + getBeastScores()[index];
  }

  /**
   * Computes beast score using the default alpha
   */
  public void computeBeastScores() {
    computeBeastScores(DEFAULT_ALPHA);
  }

  /**
   * Inverse transforms data, scales chunks (usually chromosomes) by MAD scale factor, uses scaled
   * chromosome values to obtain scaled MAD values for input data
   *
   * @param alpha
   * @param computeSTDevRaw compute standard deviation of each chunk corresponding to the
   *        indicesForScores from raw data
   * @param computeSTDevtransformed compute standard deviation of each chunk corresponding to the
   *        indicesForScores from inverseTransformed and scaled data
   *
   */
  public void computeBeastScores(float alpha) {
    float[] inverseTransformedDataScaleMAD = getinverseTransformedDataScaleMAD(SCALE_FACTOR_MAD);
    beastHeights = getBeastHeights(inverseTransformedDataScaleMAD, indicesForScores, use, log);
    beastLengths = getBeastLengths(inverseTransformedDataScaleMAD, indicesForScores, use, log);
    beastScores = getBeastScores(beastHeights, beastLengths, alpha, log);
  }

  /**
   * Consolidates the data transformations, just in case you want to operate on the transformed data
   */
  public float[] getinverseTransformedDataScaleMAD(double scaleFactorMAD) {
    float[] inverseTransformedData = transformData(inputData, indicesToChunk, use, log);
    float[] indicesMADScaled = getscaleMADIndices(indicesToChunk, inverseTransformedData, use,
                                                  scaleFactorMAD, log);
    float[] inverseTransformedDataScaleMAD = getscaleMADData(inverseTransformedData, indicesToChunk,
                                                             use, indicesMADScaled, log);
    return inverseTransformedDataScaleMAD;
  }// JOHN hijack this

  public float[] getScaleMadRawData(double scaleFactorMAD) {
    float[] indicesMADScaled = getscaleMADIndices(indicesToChunk, inputData, use, scaleFactorMAD,
                                                  log);
    float[] madScale = getscaleMADData(inputData, indicesToChunk, use, indicesMADScaled, log);
    return madScale;
  }

  /**
   *
   * @param inputData
   * @param indicesToTransform
   * @return inputData inverse Transformed with 5 df according to the indicesToChunk
   */
  public static float[] transformData(float[] inputData, int[][] indicesToChunk, boolean[] use,
                                      Logger log) {
    return Transforms.transform(inputData, 2, indicesToChunk,
                                use == null ? Array.booleanArray(indicesToChunk.length, true)
                                            : use);
  }

  /**
   *
   * Computes the median across indicesToScale of inverseTransformedData and scales by
   * SCALE_FACTOR_MAD
   *
   * @param indicesToScale
   * @param inverseTransformedData
   * @param useScaleFactorMAD using {@link BeastScore#SCALE_FACTOR_MAD} mimics the computation
   *        performed by beast
   * @param log
   * @return
   */
  public static float[] getscaleMADIndices(int[][] indicesToScale, float[] inverseTransformedData,
                                           boolean[] use, double scaleFactorMAD, Logger log) {
    float[] indicesMADScaled = new float[indicesToScale.length];
    for (int i = 0; i < indicesToScale.length; i++) {
      if (indicesToScale[i] != null && indicesToScale[i].length > 0) {
        ArrayList<Float> medianIndices = new ArrayList<Float>();
        for (int j = 0; j < indicesToScale[i].length; j++) {
          int index = indicesToScale[i][j];
          if (!Double.isNaN(inverseTransformedData[index]) && (use == null || use[index])) {
            medianIndices.add(Math.abs(inverseTransformedData[index]));
          }
        }
        indicesMADScaled[i] =
                            (float) (Array.median(Doubles.toArray(medianIndices)) / scaleFactorMAD);
      }
    }
    return indicesMADScaled;
  }

  /**
   *
   * scales inverseTransformedData according to indicesToScale and by the factors in scaleMAD
   *
   * @param data
   * @param indicesToScale
   * @param scaleMAD
   * @param log
   * @return the inverse transformed data scaled by the MAD according to indicesToScale
   */
  public static float[] getscaleMADData(float[] data, int[][] indicesToScale, boolean[] use,
                                        float[] scaleMAD, Logger log) {
    if (scaleMAD.length != indicesToScale.length) {
      log.reportError("Error - the indices to scale and the factors to scale by must have the same length");
      return null;
    }
    float[] inverseTransformedDataScaleMAD = new float[data.length];
    for (int i = 0; i < indicesToScale.length; i++) {
      if (indicesToScale != null && indicesToScale[i].length > 0) {
        for (int j = 0; j < indicesToScale[i].length; j++) {
          int index = indicesToScale[i][j];
          if (use == null || use[index]) {
            inverseTransformedDataScaleMAD[index] = data[index]
                                                    / (scaleMAD[i] == 0 ? 1 : scaleMAD[i]);
          } else {
            inverseTransformedDataScaleMAD[index] = data[index];
          }
        }
      } else if (i < 24) {
        // log.reportError("Warning - the index " + i + " was missing data");
      }
    }
    return inverseTransformedDataScaleMAD;
  }

  /**
   * @param inverseTransformedDataScaleMAD
   * @param indicesForHeights
   * @return the median of inverseTransformedDataScaleMAD according to the indices in
   *         indicesForHeights
   */
  public static float[] getBeastHeights(float[] inverseTransformedDataScaleMAD,
                                        int[][] indicesForHeights, boolean[] use, Logger log) {
    float[] beastHeights = new float[indicesForHeights.length];
    for (int i = 0; i < indicesForHeights.length; i++) {
      ArrayList<Float> medianHeightIndices = new ArrayList<Float>();
      for (int j = 0; j < indicesForHeights[i].length; j++) {
        if ((use == null || use[indicesForHeights[i][j]])
            && !Float.isNaN(inverseTransformedDataScaleMAD[indicesForHeights[i][j]])) {
          medianHeightIndices.add(inverseTransformedDataScaleMAD[indicesForHeights[i][j]]);
        }
      }
      if (medianHeightIndices.size() > 0) {
        beastHeights[i] = (float) (Array.median(Doubles.toArray(medianHeightIndices)));
      } else {
        beastHeights[i] = Float.NaN;
      }
    }
    return beastHeights;
  }

  /**
   * @param inverseTransformedDataScaleMAD
   * @param indicesForLengths
   * @return the number of non - NaN inputs in inverseTransformedDataScaleMAD according to the
   *         indices in indicesForLengths
   */
  public static int[] getBeastLengths(float[] inverseTransformedDataScaleMAD,
                                      int[][] indicesForLengths, boolean[] use, Logger log) {
    int[] beastLengths = new int[indicesForLengths.length];
    for (int i = 0; i < indicesForLengths.length; i++) {
      beastLengths[i] = 0;
      for (int j = 0; j < indicesForLengths[i].length; j++) {
        if ((use == null || use[indicesForLengths[i][j]])
            && !Float.isNaN(inverseTransformedDataScaleMAD[indicesForLengths[i][j]])) {
          beastLengths[i]++;
        }
      }
    }
    return beastLengths;
  }

  /**
   * @param beastHeights
   * @param beastLengths
   * @param alpha
   * @param log
   * @return the beast scores for each index of beastHeights and beastLengths
   */
  public static float[] getBeastScores(float[] beastHeights, int[] beastLengths, float alpha,
                                       Logger log) {
    float[] beastScores = new float[beastHeights.length];
    if (beastHeights.length != beastLengths.length) {
      log.reportError("Error - heights and lengths must contain the same number of inputs to compute a Beast Score");
      return null;
    } else {
      for (int i = 0; i < beastHeights.length; i++) {
        beastScores[i] = scoreBeast(beastLengths[i], alpha, beastHeights[i], log);
      }
    }
    return beastScores;
  }

  /**
   * @param length
   * @param alpha
   * @param height
   * @param log
   * @return the beast score using using length^alpha * height We do not implement a max/min height
   */
  public static float scoreBeast(int length, float alpha, float height, Logger log) {
    return (float) Math.abs(Math.pow(length, alpha) * height);
  }

  /**
   * Helper function to extract indices of markers contained in a CNVariant
   *
   * @param chr
   * @param positions
   * @param cnVariant
   * @param log
   * @return
   */
  public static int[] getCNVMarkerIndices(byte[] chr, int[] positions, CNVariant cnVariant,
                                          Logger log) {
    int[] indices = new int[cnVariant.getNumMarkers()];
    int count = 0;
    if (chr.length != positions.length) {
      log.reportError("Error - the chromosome and position arrays must be the same length");
      return null;
    }
    for (int i = 0; i < chr.length; i++) {
      if (cnVariant.overlaps(new Segment(chr[i], positions[i], positions[i]))) {
        if (count < indices.length) {
          indices[count] = i;
        }
        count++;
      }
    }
    if (count > indices.length) {
      log.reportError("Warning - found too many markers (" + count + ") for cnVariant "
                      + cnVariant.toPlinkFormat()
                      + "\n Perhaps some markers were filtered out prior to calling?");
      log.reportError("Warning - the beast score will be computed for the first "
                      + cnVariant.getNumMarkers() + " found in this region, and may be inaccurate");

    }
    if (count < indices.length) {
      log.reportError("Error - not enough markers were found for cnVariant"
                      + cnVariant.toPlinkFormat() + " using the current chr/position data");
      return null;
    }
    return indices;
  }

  /**
   * Helper Function to compute the beast scores for cNVariantInds[][], where cNVariantInds.lenght =
   * number of individuals and cNVariantInds[i].length is the number of cnvs per indivdual
   *
   * @return an array of beastScores, 1 per with scores computed across the individuals cnvs
   */
  public static BeastScore[] beastInds(Project proj, CNVariant[][] cNVariantInds) {
    BeastScore[] beastScores = new BeastScore[cNVariantInds.length];
    MarkerSet markerSet = proj.getMarkerSet();
    byte[] chr = markerSet.getChrs();
    int[] positions = markerSet.getPositions();
    int[][] indicesByChr = markerSet.getIndicesByChr();
    SampleData sampleData = proj.getSampleData(0, false);

    for (int i = 0; i < cNVariantInds.length; i++) {
      beastScores[i] = beastInd(proj, sampleData, cNVariantInds[i], chr, positions, indicesByChr);
    }
    return beastScores;
  }

  public static BeastScore beastInd(Project proj, SampleData sampleData, CNVariant[] cNVariantInd,
                                    byte[] chr, int[] positions, int[][] indicesByChr) {
    return beastInd(proj, sampleData, null, cNVariantInd, chr, positions, indicesByChr);
  }


  /**
   * Helper function to compute beast scores for CNVariant[] representing a single individual
   *
   * @return BeastScore (for all the individual's cnvs)
   */
  public static BeastScore beastInd(Project proj, SampleData sampleData, float[] lrrs,
                                    CNVariant[] cNVariantInd, byte[] chr, int[] positions,
                                    int[][] indicesByChr) {
    BeastScore score = null;
    int[][] indices = new int[cNVariantInd.length][];
    String ind = null;
    Logger log = proj.getLog();

    if (cNVariantInd.length > 0) {
      String key = cNVariantInd[0].getFamilyID() + "\t" + cNVariantInd[0].getIndividualID();
      if (lrrs == null) {
        try {
          ind = sampleData.lookup(key)[0];
        } catch (NullPointerException npe) {
          log.reportError("Error - could not look up the sample " + key
                          + " in the sample data file " + proj.SAMPLE_DATA_FILENAME.getValue()
                          + ", cannot load sample to compute beast score");
          log.reportError("Error - please ensure that the sample names correspond to the varaints being processed with FID="
                          + cNVariantInd[0].getFamilyID() + " and IID="
                          + cNVariantInd[0].getIndividualID());
          return score;
        }
      }
      for (int i = 0; i < cNVariantInd.length; i++) {
        if (!key.equals(cNVariantInd[i].getFamilyID() + "\t" + cNVariantInd[i].getIndividualID())) {
          log.reportError("Error - this method can only be used with cnvs from the same individual according to FID\tIID identifiers, returning null score");
          return score;
        } else {
          indices[i] = getCNVMarkerIndices(chr, positions, cNVariantInd[i], log);
        }
      }
      if (lrrs == null) {
        try {
          lrrs = proj.getFullSampleFromRandomAccessFile(ind).getLRRs();
        } catch (NullPointerException npe) {
          log.reportError("Error - could not load data for the sample " + ind + "\t" + key
                          + ", please ensure samples have been parsed prior to computing beast score");
          log.report("Skipping beast score for sample " + ind + "\t" + key);
          return score;
        }
      }
      score = new BeastScore(lrrs, indicesByChr, indices, log);
      score.computeBeastScores();
    } else {
      log.reportError("Warning - no CNVariants were found for an individual, returning a null BeastScore");
    }
    return score;
  }
}
