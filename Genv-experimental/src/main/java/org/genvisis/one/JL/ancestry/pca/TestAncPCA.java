package org.genvisis.one.JL.ancestry.pca;

import java.io.File;
import org.genvisis.common.Logger;
import org.genvisis.common.matrix.SVD;
import org.genvisis.pca.ancestry.AncestryPCA;
import org.genvisis.pca.ancestry.PlinkDataLoader;

public class TestAncPCA {

  public static void main(String[] args) {
    long time = System.currentTimeMillis();
    String dir = "/Volumes/Beta2/Poynter/";
    String plinkRoot = "plink";

    String outDir = dir + "testOut/";
    new File(outDir).mkdirs();
    Logger log = new Logger(outDir + "test.log");

    SVD svd = AncestryPCA.generatePCs(new PlinkDataLoader(dir, plinkRoot, log), log);
    svd.dumpPCsToText(dir + "testPCS/testMean0_posEstAF", log);

    log.reportTimeElapsed(time);

  }

}
