package org.genvisis.cnv.manage;

import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.mining.Transformations;

public class Transforms {

  public static enum TRANSFORMATION {
    RAW("Raw Values", -3, 3),
    INV_NORM("Inverse normalized", 0, 1),
    INV_TD_5DF("Inverse T-distribution with 5 df", -8, 8),
    BEASTV("BEAST vision", -12, 12),
    MULT_5X("5X multiply", -1, -1);

    private TRANSFORMATION(String desc, int dispMin, int dispMax) {
      this.desc = desc;
      this.displayMin = dispMin;
      this.displayMax = dispMax;
    }

    private final String desc;
    private final int displayMin;
    private final int displayMax;

    public String getDescription() {
      return desc;
    }

    public int getDisplayMin() {
      return displayMin;
    }

    public int getDisplayMax() {
      return displayMax;
    }
  }

  public static final String[] TRANFORMATIONS = {"Raw Values", "Inverse normalized",
                                                 "Inverse T-distribution with 5 df", "BEAST vision",
                                                 "5X multiply"};
  public static final int[] TRANSFORMATION_TYPES = {Transformations.IDENTITY,
                                                    Transformations.INVERSE_NORMALIZE,
                                                    Transformations.INVERSE_TDIST_5DF,
                                                    Transformations.MAD_SCALED, Transformations.X5};
  public static final String[] SCOPES = {"Chromosome", "Genome"};

  public static float[] transform(float[] input, int transformation_type,
                                  boolean transformSeparatelyByChromosome,
                                  MarkerSetInfo markerSet) {
    int[][] indices;

    if (transformSeparatelyByChromosome) {
      if (markerSet.getChrs().length != input.length) {
        System.err.println("Error - cannot transform by chromosome; mismatched number of records between array and MarkerSet");
        return null;
      }
      indices = markerSet.getIndicesByChr();
    } else {
      indices = new int[][] {ArrayUtils.arrayOfIndices(input.length)};
    }

    return transform(input, transformation_type, indices,
                     ArrayUtils.booleanArray(indices.length, true));
  }

  public static float[] transform(float[] input, int transformation_type, int[][] indices,
                                  boolean[] transform) {
    indices = ArrayUtils.deepCopy(indices);
    float[] output, trav;
    int count;
    // long time;
    // String timeString;

    // time = new Date().getTime();
    output = new float[input.length];
    for (int i = 0; i < indices.length; i++) {
      if (transform[i] && indices[i].length > 0) {
        count = 0;
        trav = new float[indices[i].length];
        for (int j = 0; j < indices[i].length; j++) {
          if (Float.isNaN(input[indices[i][j]])) {
            if (indices[i][j] == 0) {
              indices[i][j] = Integer.MAX_VALUE;
            } else {
              indices[i][j] *= -1;
            }
          } else {
            trav[count] = input[indices[i][j]];
            count++;
          }
        }

        trav = Transformations.transform(ArrayUtils.subArray(trav, 0, count),
                                         TRANSFORMATION_TYPES[transformation_type]);

        count = 0;
        for (int j = 0; j < indices[i].length; j++) {
          if (indices[i][j] < 0) {
            indices[i][j] *= -1;
            output[indices[i][j]] = Float.NaN;
          } else if (indices[i][j] == Integer.MAX_VALUE) {
            indices[i][j] = 0;
            output[indices[i][j]] = Float.NaN;
          } else {
            output[indices[i][j]] = trav[count];
            count++;
          }
        }
      }

    }
    // timeString = "sorting: "+ext.getTimeElapsed(time);

    // System.out.println(transformation_type+"\t"+timeString);

    return output;
  }

}
