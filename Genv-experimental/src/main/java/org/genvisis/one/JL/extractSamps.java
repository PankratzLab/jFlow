package org.genvisis.one.JL;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.LocusSet;

public class extractSamps {
  public static void main(String[] args) {
    String filename = "C:/data/ARIC/shadowCNVs/combinedMF.cnv";
    String sampFile = "C:/data/ARIC/shadowCNVs/VTECasesOnly.txt";
    String sampCNVs = ext.addToRoot(filename, "." + ext.rootOf(sampFile));
    HashSet<String> sampSet = HashVec.loadFileToHashSet(sampFile, false);
    if (!Files.exists(sampCNVs)) {
      LocusSet<CNVariant> cnvs = CNVariant.loadLocSet(filename, new Logger());

      try {
        int numTotal = 0;
        int numSamps = 0;
        PrintWriter writer = new PrintWriter(new FileWriter(sampCNVs));
        writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
        for (int i = 0; i < cnvs.getLoci().length; i++) {
          numTotal++;
          if (sampSet.contains(cnvs.getLoci()[i].getIndividualID())) {
            writer.println(cnvs.getLoci()[i].toPlinkFormat());
            numSamps++;
          }
        }
        writer.close();
        new Logger().reportTimeInfo(numTotal + " - > " + numSamps);
      } catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
    }
  }
}
