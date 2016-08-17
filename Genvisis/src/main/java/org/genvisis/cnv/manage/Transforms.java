package org.genvisis.cnv.manage;

import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.common.Array;
import org.genvisis.mining.Transformations;

public class Transforms {
  // public static final String[] TRANFORMATIONS = {"Raw Values", "Quantile", "Inverse normalized",
  // "Inverse T-distribution with 5 df"};
  // public static final int[] TRANSFORMATION_TYPES = {Transformations.IDENTITY,
  // Transformations.QUANTILE, Transformations.INVERSE_NORMALIZE,
  // Transformations.INVERSE_TDIST_5DF};
  public static final String[] TRANFORMATIONS = {"Raw Values", "Inverse normalized",
      "Inverse T-distribution with 5 df", "BEAST vision", "5X multiply"};
  public static final int[] TRANSFORMATION_TYPES =
      {Transformations.IDENTITY, Transformations.INVERSE_NORMALIZE,
          Transformations.INVERSE_TDIST_5DF, Transformations.MAD_SCALED, Transformations.X5};
  public static final String[] SCOPES = {"Chromosome", "Genome"};

  public static float[] transform(float[] input, int transformation_type,
      boolean transformSeparatelyByChromosome, MarkerSet markerSet) {
    int[][] indices;

    if (transformSeparatelyByChromosome) {
      if (markerSet.getChrs().length != input.length) {
        System.err.println(
            "Error - cannot transform by chromosome; mismatched number of records between array and MarkerSet");
        return null;
      }
      indices = markerSet.getIndicesByChr();
    } else {
      indices = new int[][] {Array.arrayOfIndices(input.length)};
    }

    return transform(input, transformation_type, indices, Array.booleanArray(indices.length, true));
  }

  public static float[] transform(float[] input, int transformation_type, int[][] indices,
      boolean[] transform) {
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

        trav = Transformations.transform(Array.subArray(trav, 0, count),
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
