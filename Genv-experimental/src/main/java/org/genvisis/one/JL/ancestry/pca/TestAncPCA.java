package org.genvisis.one.JL.ancestry.pca;

import java.io.File;
import org.genvisis.cnv.analysis.pca.ancestry.AncestryPCA;
import org.genvisis.cnv.analysis.pca.ancestry.PlinkDataLoader;
import org.genvisis.common.Logger;

public class TestAncPCA {

  public static void main(String[] args) {
    String dir = "/Volumes/Beta2/Poynter/";
    String plinkRoot = "plink";

    String outDir = dir + "testOut/";
    new File(outDir).mkdirs();
    Logger log = new Logger(outDir + "test.log");

    AncestryPCA.computePCA(new PlinkDataLoader(dir, plinkRoot, log).getData(), log);

  }

}
