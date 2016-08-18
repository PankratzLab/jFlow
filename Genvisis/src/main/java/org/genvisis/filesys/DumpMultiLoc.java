package org.genvisis.filesys;

import java.io.FileWriter;
import java.io.PrintWriter;

import org.genvisis.common.Logger;

public class DumpMultiLoc {

  public static void dumpMultiLoc(String geneTrackFile, String outputFile, Logger log) {
    GeneTrack geneTrack = GeneTrack.load(geneTrackFile, false);
    GeneData[][] genes = geneTrack.getGenes();

    try {
      PrintWriter writer = new PrintWriter(new FileWriter(outputFile));
      for (GeneData[] gene : genes) {
        for (int j = 0; j < gene.length; j++) {
          if (gene[j].getMultiLoc() > 0) {
            writer.println(gene[j].getChr() + "\t" + gene[j].getStart() + "\t" + gene[j].getStop());
          }
        }
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + outputFile);
      log.reportException(e);
    }

  }

}
