package org.genvisis.one.JL;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import org.genvisis.cnv.LocusSet;
import org.genvisis.cnv.filesys.CNVariant;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;

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
        PrintWriter writer = Files.openAppropriateWriter(sampCNVs);
        writer.println(ArrayUtils.toStr(CNVariant.PLINK_CNV_HEADER));
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
