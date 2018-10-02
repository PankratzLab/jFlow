package org.genvisis.cnv.util;

import java.io.PrintWriter;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.ext;
import org.pankratzlab.shared.filesys.CNVariant;
import htsjdk.tribble.annotation.Strand;

/**
 * Static utility helper methods for working with CNVs
 */
public final class CNVHelper {

  public static void generateRegionsFileFromCNVFile(String cnvFile) {
    CNVariant[] cnvs = CNVariant.loadPlinkFile(cnvFile);
    String outFile = ext.rootOf(cnvFile, false) + "_regions.txt";
    PrintWriter writer = Files.getAppropriateWriter(outFile);

    for (CNVariant cnv : cnvs) {
      String ucsc = cnv.getUCSClocation();
      if (ucsc == null) {
        ucsc = "chr" + cnv.getChr() + ":" + cnv.getStart() + "-" + cnv.getStop();
      }
      writer.println(cnv.getIndividualID() + "\t" + ucsc);
    }

    writer.flush();
    writer.close();
  }

  /**
   * @param strand
   * @return Single-character string equivalent of the given Strand enum.
   */
  public static String decode(Strand strand) {
    switch (strand) {
      case NEGATIVE:
        return "-";
      case POSITIVE:
        return "+";
      case NONE:
      default:
        return "!";
    }

  }
}
