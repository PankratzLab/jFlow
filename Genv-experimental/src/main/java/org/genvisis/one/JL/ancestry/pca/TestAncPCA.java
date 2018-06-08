package org.genvisis.one.JL.ancestry.pca;

import java.io.File;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.matrix.SVD;
import org.genvisis.pca.ancestry.AncestryPCA;
import org.genvisis.pca.ancestry.PlinkDataLoader;

public class TestAncPCA {

  public static void main(String[] args) {
    long time = System.currentTimeMillis();
    String dir = "/Volumes/Beta2/Poynter/";
    //    String dir = "/Volumes/Beta2/NGS/topmed/aims/plinkfreeze.5b.aims.pass_and_fail.gtonly.minDP10.vcf/quality_control/ancestryFull/unrelateds/";

    String plinkRoot = "plink";

    String outDir = dir + "testOut/";
    new File(outDir).mkdirs();
    Logger log = new Logger(outDir + "test.log");
    String ser = outDir + "svdBase.ser.gz";
    if (!Files.exists(ser)) {
      AncestryPCA.generatePCs(new PlinkDataLoader(dir, plinkRoot, log), 10, log).writeSerial(ser);

    }
    SVD svd = SVD.readSerial(ser);
    svd.dumpPCsToText(outDir + "testMean0_posEstAF", "SAMPLE", log);
    svd.dumpLoadingsToText(outDir + "testMean0_posEstAF", "MARKER", log);
    log.reportTimeElapsed(time);
    AncestryPCA.extrapolatePCs(svd, new PlinkDataLoader(dir, plinkRoot, log), log)
               .dumpToText(outDir + "testMean0_posEstAF.extrapolatedPcs.txt", "SAMPLE", log);
  }

}
