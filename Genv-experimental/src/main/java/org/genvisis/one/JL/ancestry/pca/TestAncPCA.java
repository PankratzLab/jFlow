package org.genvisis.one.JL.ancestry.pca;

import java.io.File;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.matrix.SVD;
import org.genvisis.pca.ancestry.AncestryPCA;
import org.genvisis.pca.ancestry.PlinkDataMatrixLoader;

public class TestAncPCA {

  public static void main(String[] args) {
    long time = System.currentTimeMillis();
    String dir = "/Volumes/Beta2/Poynter/";
    //    dir = "/Volumes/Beta2/NGS/topmed/aims/plinkfreeze.5b.aims.pass_and_fail.gtonly.minDP10.vcf/quality_control/ancestryFull/unrelateds/";

    String plinkRoot = "plink";

    String outDir = dir + "testOutCenterMean0NormEJML_Testing2/";
    new File(outDir).mkdirs();
    Logger log = new Logger(outDir + "test.log");
    String ser = outDir + "svdBase.ser.gz";
    if (!Files.exists(ser)) {
      AncestryPCA.generatePCs(new PlinkDataMatrixLoader(dir, plinkRoot, log), 10, log)
                 .writeSerial(ser);

    }
    AncestryPCA svd = AncestryPCA.readSerial(ser);
    svd.getSvd().dumpPCsToText(outDir + "test", "FID\tIID", log);
    svd.getSvd().dumpLoadingsToText(outDir + "test", "MARKER", log);

    AncestryPCA.extrapolatePCs(svd, new PlinkDataMatrixLoader(dir, plinkRoot, log), log)
               .dumpToText(outDir + "test.extrapolatedPcs.txt.gz", "FID\tIID", log);
    log.reportTimeElapsed(time);

  }

}
