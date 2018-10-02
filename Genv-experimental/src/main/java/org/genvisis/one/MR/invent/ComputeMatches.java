package org.genvisis.one.MR.invent;

import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.gwas.MatchSamples;
import org.pankratzlab.gwas.MatchesVisualized;

public class ComputeMatches {

  public static void main(String[] args) {
    String d = "D:\\tWork\\SequencingProjectWithCIDR\\MatchingControls\\MatchingForMito\\";
    String anchors = "anchor_cases.dat";
    String barnaclesFile = "barnacle_controls.dat";
    String factors = "ARIC_EA_VTPEDXall.txt";
    // int[] factorIndices = new int[] {1,2,3,4};
    // int[] factorIndices = new int[] {1,2,3,4,5,6,7,8,9,10};
    // int[] factorIndices = new int[] {1,2};

    // String[] factorTargets = new String[] {"C1_norm", "C2_norm",
    // "Age_norm", "AgeAtExam_norm"};
    // double[] factorLoadings = new double[] {2, 2, 4, 1};

    String clusterfile = "cluster.genome";
    String file, pairs;

    String usage = "\\n" + "gwas.MatchSamples requires 0-1 arguments\n"
                   + "   (0) directory (i.e. dir=" + d + " (default))\n"
                   + "   (1) anchors (i.e. anchors=" + anchors + " (default))\n"
                   + "   (2) barnacles (i.e. barnacles=" + barnaclesFile + " (default))\n"
                   + "   (3) file with factors (i.e. factors=" + factors + " (default))\n" +
                   // " (4) indices of factors in clusterfile (i.e.
                   // indices="+Array.toStr(factorIndices, ",") +" (default))\n" +
                   "   (5) clusterfile (i.e. clusterfile=" + clusterfile + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        d = arg.split("=")[1];
      } else if (arg.startsWith("anchors=")) {
        anchors = arg.split("=")[1];
      } else if (arg.startsWith("barnacles=")) {
        barnaclesFile = arg.split("=")[1];
      } else if (arg.startsWith("factors=")) {
        factors = arg.split("=")[1];
      } else if (arg.startsWith("clusterfile=")) {
        clusterfile = arg.split("=")[1];
      }
    }

    System.out.println("Factors: " + factors);

    try {

      int iterations = 4;

      for (int i = 1; i <= iterations; i++) {
        String dir = d + "/matches_" + i + "/";
        Files.ensurePathExists(dir);

        String[] barnacles = HashVec.loadFileToStringArray(d + "/" + barnaclesFile, false, null,
                                                           true);

        file = MatchSamples.matchMaker(dir, "/../" + anchors, "/../" + barnaclesFile,
                                       "/../" + factors,
                                       new String[] {"PC1", "PC2", "AGE01", "SEX"},
                                       new double[] {16, 16, 4, 1}, true);

        pairs = MatchSamples.matchPairs(dir, file, true);

        System.out.println(pairs);

        MatchSamples.evalAgeSex_and_MDS_separately(dir, pairs,
                                                   "distances_PC1x16,PC2x16,AGE01x4,SEXx1.xln",
                                                   "/../" + factors, "AGE01", "SEX");
        MatchSamples.eval(dir, pairs, "/../" + factors, new String[] {"AGE01", "SEX=concordance"});

        new MatchesVisualized(dir, "/../" + anchors, "/../" + barnaclesFile, "/../" + factors,
                              new int[] {4, 5}, pairs);

        if (i < iterations) {
          String[] barns = HashVec.loadFileToStringArray(dir + pairs, true, new int[] {1}, true);
          String[] updatedBarns = ArrayUtils.removeFromArray(barnacles, barns);
          Files.writeArray(updatedBarns, d + "/" + "barnacles_" + i + ".dat");
          // update the barnacles array to the new file we just made
          barnaclesFile = "barnacles_" + i + ".dat";
        }
      }

    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
