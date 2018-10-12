package org.genvisis.one.JL;

import java.io.File;
import java.io.PrintWriter;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.shared.stats.Rscript;
import org.pankratzlab.shared.stats.Rscript.RScatter;
import org.pankratzlab.shared.stats.Rscript.SCATTER_TYPE;

public class plotContam {

  public static void main(String[] args) {
    String[] files = new String[] {"D:/data/Project_Tsai_21_25_26_spector/QC/bamContam/fullBamContam.txt",
                                   "D:/data/Tsai_Project_028/bamContam/fileOfBams.txtcontamSummary.txt"};
    String dir = "D:/data/Project_Tsai_21_25_26_spector/QC/bamContam/images/";
    new File(dir).mkdirs();
    for (String file : files) {
      String[] samps = Files.getHeaderOfFile(file, new Logger());
      String[] sampsToPlot = HashVec.loadFileToStringArray("D:/data/Project_Tsai_21_25_26_spector/QC/bamContam/CushingSamps.txt",
                                                           false, new int[] {0}, true);

      String[][] all = HashVec.loadFileToStringMatrix(file, false, null);
      int[] indices = ext.indexFactors(sampsToPlot, samps, false);
      indices = ArrayUtils.removeAllValues(indices, -1);

      String withMedian = ext.rootOf(file, false) + "withMedian.txt";
      try {
        PrintWriter writer = Files.openAppropriateWriter(withMedian);
        writer.println(all[0][0] + "\t" + ArrayUtils.toStr(ArrayUtils.subArray(all[0], indices))
                       + "\tPopulationMedian");
        for (int j = 1; j < all.length; j++) {
          double[] counts = ArrayUtils.toDoubleArray(ArrayUtils.subArray(all[j], indices));
          double median = ArrayUtils.median(counts);
          writer.println(all[j][0] + "\t" + ArrayUtils.toStr(ArrayUtils.subArray(all[j], indices))
                         + "\t" + median);
        }
        writer.close();
      } catch (Exception e) {
        // log.reportError("Error writing to " + name);
        // log.reportException(e);
      }

      String[] newheader = Files.getHeaderOfFile(withMedian, new Logger());
      for (int j = 1; j < samps.length; j++) {
        int index = ext.indexOfStr(samps[j], sampsToPlot);
        if (index >= 0) {
          String root = dir + ext.rootOf(file) + Rscript.makeRSafe(samps[j]);
          double[] data = ArrayUtils.toDoubleArray(HashVec.loadFileToStringArray(withMedian, true,
                                                                                 new int[] {ext.indexOfStr(samps[j],
                                                                                                           newheader)},
                                                                                 false));
          double kurt = ArrayUtils.kurtosis(data);
          double skew = ArrayUtils.skewness(data);
          RScatter rscScatter = new RScatter(withMedian, root + ".rscript",
                                             ext.removeDirectoryInfo(root), root + ".jpeg",
                                             "PROP_REF_BIN",
                                             new String[] {samps[j], "PopulationMedian"},
                                             SCATTER_TYPE.POINT, new Logger());
          rscScatter.setTitle(samps[j] + "\nkurtosis=" + ext.formDeci(kurt, 3) + "\nskew="
                              + ext.formDeci(skew, 3));
          rscScatter.setyLabel("Number of Reads");
          rscScatter.setxLabel("Proportion Reference");
          rscScatter.setOverWriteExisting(true);
          rscScatter.execute();
        }
      }
    }
  }
}
