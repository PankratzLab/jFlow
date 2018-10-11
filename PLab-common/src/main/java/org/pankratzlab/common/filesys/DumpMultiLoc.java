package org.pankratzlab.common.filesys;

import java.io.PrintWriter;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;

public class DumpMultiLoc {

  public static void dumpMultiLoc(String geneTrackFile, String outputFile, Logger log) {
    GeneTrack geneTrack = GeneTrack.load(geneTrackFile);
    GeneData[][] genes = geneTrack.getGenes();

    try {
      PrintWriter writer = Files.openAppropriateWriter(outputFile);
      for (GeneData[] gene : genes) {
        for (GeneData element : gene) {
          if (element.getMultiLoc() > 0) {
            writer.println(element.getChr() + "\t" + element.getStart() + "\t" + element.getStop());
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
