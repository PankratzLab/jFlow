package org.genvisis.cnv.manage;

import org.genvisis.cnv.filesys.Sample;
import org.genvisis.common.ext;

public class SourceParserDep {

  public static byte[] parseGenotypes(String[] line, int[] genotypeIndices, boolean ignoreAB,
                                      char[][] abLookup, int count, String sampleName,
                                      String markerName, String filename) {
    String genotype;
    byte genoForward, genoAB;

    if (genotypeIndices[4] >= 0) {
      genotype = line[genotypeIndices[4]];
    } else {
      genotype = line[genotypeIndices[0]] + line[genotypeIndices[1]];
    }
    genoForward = (byte) ext.indexOfStr(genotype, Sample.ALLELE_PAIRS);
    if (genoForward == -1) {
      if (ext.indexOfStr(genotype, Sample.ALT_NULLS) == -1) {
        // System.err.println("Error - failed to lookup "+genotype+" for marker "+markerName+" of
        // sample "+filename);
        genoForward = -9;
      } else {
        genoForward = 0;
      }
    }
    if (ignoreAB) {
      genoAB = -1;
      // do nothing, will need to use these files to determine AB lookup table
    } else if (abLookup == null) {
      if (genotypeIndices[5] >= 0) {
        genotype = line[genotypeIndices[5]];
      } else {
        genotype = line[genotypeIndices[2]] + line[genotypeIndices[3]];
      }
      genoAB = (byte) ext.indexOfStr(genotype, Sample.AB_PAIRS);
      if (genoAB == -1) {
        if (ext.indexOfStr(genotype, Sample.ALT_NULLS) == -1) {
          System.err.println("Error - failed to lookup " + genotype + " for marker " + markerName
                             + " of sample " + filename);
        }
      }
    } else {
      if (genoForward == 0) {
        genoAB = -1;
      } else {
        genoAB = 0;
        if (genotypeIndices[4] >= 0) {
          genotype = line[genotypeIndices[4]];
        } else {
          genotype = line[genotypeIndices[0]] + line[genotypeIndices[1]];
        }
        for (int j = 0; j < 2; j++) {
          if (genotype.charAt(j) == abLookup[count][1]) {
            genoAB++;
          } else if (genotype.charAt(j) != abLookup[count][0]) {
            System.err.println("Error - alleles for individual '" + sampleName + "' (" + genotype
                               + ") do not match up with the defined AB lookup alleles ("
                               + abLookup[count][0] + "/" + abLookup[count][1] + ") for marker "
                               + markerName);
          }
        }
      }
    }

    return new byte[] {genoForward == -9 ? 1 : genoForward, genoAB};
  }

}
