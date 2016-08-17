package org.genvisis.cnv.qc;

import java.io.FileWriter;
import java.io.PrintWriter;

import org.genvisis.cnv.filesys.ClusterFilterCollection;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.common.Array;
import org.genvisis.common.HashVec;
import org.genvisis.stats.Correlation;

public class CompareDuplicates {

  public static void allPairs(Project proj, String pairFile) {
    PrintWriter writer;
    String[][] pairs;
    String[] markerNames;
    int[] discordantCounts;

    markerNames = proj.getMarkerNames();
    discordantCounts = new int[markerNames.length];

    pairs = HashVec.loadFileToStringMatrix(proj.PROJECT_DIRECTORY.getValue() + pairFile, false,
        new int[] {0, 1}, "\t", proj.JAR_STATUS.getValue(), 100, false);

    try {
      writer =
          new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + "DuplicateQC.xln"));
      writer.println("DNA1\tDNA2\tgenotypeConcordance\tx_correlation\tcorr_pval\t#diffGeno");
      for (int i = 0; i < pairs.length; i++) {
        System.out.println((i + 1) + " of " + pairs.length);
        writer.println(
            Array.toStr(pairs[i]) + "\t" + doIt(proj, pairs[i], markerNames, discordantCounts));
      }
      writer.close();
    } catch (Exception e) {
      System.err
          .println("Error writing to " + proj.PROJECT_DIRECTORY.getValue() + "DuplicateQC.xln");
      e.printStackTrace();
    }

    try {
      writer = new PrintWriter(
          new FileWriter(proj.PROJECT_DIRECTORY.getValue() + "DuplicateErrors.xln"));
      writer.println("MarkerName\t#discordant");
      for (int i = 0; i < markerNames.length; i++) {
        writer.println(markerNames[i] + "\t" + discordantCounts[i]);
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + "DuplicateErrors.xln");
      e.printStackTrace();
    }
  }

  public static String doIt(Project proj, String[] pair, String[] markerNames,
      int[] discordantCounts) {
    Sample[] fsamps;
    byte[] chrs;
    float[][] xs, ys;
    int max, count, sameGeno, diffGeno;
    String summary;
    double[][] correl_icin;
    boolean use;
    byte[][] genos;
    byte geno;
    ClusterFilterCollection clusterFilterCollection;
    float gcThreshold;

    clusterFilterCollection = proj.getClusterFilterCollection();
    // gcThreshold = Float.parseFloat(proj.getProperty(Project.GC_THRESHOLD));
    gcThreshold = proj.getProperty(proj.GC_THRESHOLD).floatValue();

    // System.out.println("Comparing "+ext.listWithCommas(pair, true));

    summary = "";
    try {
      fsamps = new Sample[pair.length];
      xs = new float[pair.length][];
      ys = new float[pair.length][];
      genos = new byte[pair.length][];
      for (int i = 0; i < pair.length; i++) {
        fsamps[i] = proj.getFullSampleFromRandomAccessFile(pair[i]);
        xs[i] = fsamps[i].getXs();
        ys[i] = fsamps[i].getYs();
        genos[i] = fsamps[i].getAB_GenotypesAfterFilters(markerNames, clusterFilterCollection,
            gcThreshold);
        // genos[i] = fsamps[i].getForwardGenotypes();
      }

      chrs = proj.getMarkerSet().getChrs();
      // max = Array.indexOfByte(chrs, (byte)23);
      // max = 1000000;
      max = chrs.length;

      geno = 0;
      count = 0;
      sameGeno = 0;
      diffGeno = 0;
      correl_icin = new double[pair.length][max];
      for (int i = 0; i < max; i++) {
        use = true;
        for (int j = 0; j < pair.length; j++) {
          // writer.print((i==0?"":"\t")+xs[j][i]);
          correl_icin[j][count] = xs[j][i];
          // if (genos[j][i] == 0) { // for ForwardGenotypes
          if (genos[j][i] < 0) { // for AB genotypes
            use = false;
          }
          if (j == 0) {
            geno = genos[j][i];
            // } else if (genos[j][i] != 0 && geno != 0) { // for ForwardGenotypes
          } else if (genos[j][i] != -1 && geno != -1) {
            if (genos[j][i] == geno) {
              // sameGeno++;
              sameGeno += 2 - Math.abs(genos[j][i] - geno);
            } else {
              diffGeno++;
              discordantCounts[i]++;
            }
          }
        }
        if (use) {
          count++;
        }
        // writer.println();
      }
      sameGeno /= 2;

      for (int i = 0; i < pair.length; i++) {
        correl_icin[i] = Array.subArray(correl_icin[i], 0, count);
      }
      summary += (sameGeno / ((double) pair.length - 1) / count) + "\t"
          + Array.toStr(Correlation.Pearson(correl_icin)) + "\t" + diffGeno;

      // writer.close();
    } catch (Exception e) {
      e.printStackTrace();
    }

    return summary;
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = null;
    String pairs = "duplicatePairs.txt";

    String usage = "\n" + "cnv.qc.CompareDuplicates requires 0-1 arguments\n"
        + "   (1) project properties filename (i.e. proj="
        + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
        + "   (2) file with list of sample pairs (i.e. pairs=" + pairs
        + " (default; to be found in project directory))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("pairs=")) {
        pairs = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      // filename = "D:/home/npankrat/projects/TsaiPilot.properties";
      // filename = "D:/home/npankrat/projects/SDRG.properties";
      filename = "D:/home/npankrat/projects/GEDI_exomeRAF.properties";
      allPairs(new Project(filename, false), pairs);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
