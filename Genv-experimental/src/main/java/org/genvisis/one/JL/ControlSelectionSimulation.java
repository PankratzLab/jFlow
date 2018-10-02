package org.genvisis.one.JL;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CmdLine;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ext;
import org.pankratzlab.gwas.MatchSamples;
import com.google.common.primitives.Doubles;

public class ControlSelectionSimulation {

  private static int SAMPLE_SIZE = 550;
  // private static int SAMPLE_SIZE = 200;

  private static double defC1 = 0.006;
  private static double defC2 = 0.002;

  private static List<String> selectCases(List<String> all, Hashtable<String, String> pheno,
                                          String cent, Hashtable<String, String> have, int num,
                                          double percent) {
    ArrayList<String> finalSelect = new ArrayList<>();
    System.out.println(all.size());
    int numNotIn = 0;
    for (String ca : all) {
      String currnt = ca.split("\t")[0];
      if (pheno.containsKey(currnt) && have.containsKey(currnt)) {
        // System.out.println("FDS");
        double c1 = Double.parseDouble(have.get(currnt).split("\t")[3]);
        double c2 = Double.parseDouble(have.get(currnt).split("\t")[4]);
        if (pheno.get(currnt).split("\t")[7].equals(cent)) {
          double rand = Math.random();
          if (c1 < defC1 && c2 < defC2) {
            if (rand <= percent) {
              finalSelect.add(ca);
            } else {
              numNotIn++;
            }
          } else {
            finalSelect.add(ca);
          }
        }
      }

    }
    System.out.println(finalSelect.size() + " after excluding " + numNotIn);
    Collections.shuffle(finalSelect);
    return finalSelect.subList(0, Math.min(num - 1, finalSelect.size() - 1));
  }

  public static void main(String[] args) throws IOException {
    // plink2 --bfile ../plink --keep keeps.txt --make-bed --out ./plink
    // plink2 --bfile ../plink --extract range
    // /Volumes/Beta/data/controlSelection/plink/meregedCaptureRegions.txt
    // --keep cushingPheno.txt --make-bed --out plink
    // ./plink
    // QQPlot qqPlot1 = QQPlot.loadPvals(
    // new String[] {
    // "/Volumes/Beta/data/controlSelection/withCushing/antimatch/CUSHING_AND_OTHERS_v_ARIC/CUSHING_AND_OTHERS_v_ARICantioptimal.p.txt",
    // "/Volumes/Beta/data/controlSelection/withCushing/randomMatch/CUSHING_AND_OTHERS_v_ARIC/CUSHING_AND_OTHERS_v_ARICrandom.p.txt",
    // "/Volumes/Beta/data/controlSelection/withCushing/match/CUSHING_AND_OTHERS_v_ARIC/CUSHING_AND_OTHERS_v_ARICoptimal.p.txt"
    // },
    // "Hi", false, true, false, -1, false, Float.MAX_VALUE, new Logger());
    // qqPlot1.screenCap("/Volumes/Beta/data/controlSelection/withCushing/"
    // + "test.png");
    // System.exit(1);
    // //whites
    // String dir = "/Volumes/Beta/data/controlSelection/plink/";
    // String plink1 =
    // "/Volumes/Beta/data/controlSelection/plinkARIC/quality_control/further_analysis_QC/";
    // String plink2 =
    // "/Volumes/Beta/data/controlSelection/plinkCushing/quality_control/further_analysis_QC/";
    // String plinkNew = dir + "plinkMerge/";
    String root = "/Users/Kitty/temp/controlSelect/ARIC_Only_SpecialSimulation_maf_0_05/";

    double[] thresholds = new double[] {.4, .6};
    int reps = 100;

    for (int numS = 500; numS < SAMPLE_SIZE; numS += 52) {
      for (int k = 0; k < thresholds.length; k++) {
        ArrayList<String> lambdaReport = new ArrayList<>();
        lambdaReport.add("Comparison\tmethod\tnumSamples\tlambda\tphenoFile");
        final double threshold = thresholds[k];
        final int i = numS;
        // new Thread(new Runnable() {
        //
        // public void run() {
        String phenoFile = "/Volumes/Beta/data/controlSelection/plink/centerPheno2.txt";

        // String[] phenoFiles =
        // Files.listFullPaths("/Volumes/Beta/data/controlSelection/centerPhenos/",
        // ".txt",
        // false);
        String pheno = ext.rootOf(phenoFile);
        String rootOut = root + pheno + "_select" + i + "threshold_" + threshold + "/";
        // new File(plinkNew).mkdirs();
        new File(rootOut).mkdirs();
        String mds = "/Volumes/Beta/data/controlSelection/plinkARIC/quality_control/genome/mds20.mds";
        String keepSamples = "/Volumes/Beta/data/controlSelection/plinkARIC/quality_control/genome/plink.genome_keep.dat";

        String plinkQc = "/Volumes/Beta/data/controlSelection/plinkARIC/quality_control/further_analysis_QC/plink_QCd";

        Hashtable<String, String> have = HashVec.loadFileToHashString(mds, new int[] {0},
                                                                      new int[] {0, 1, 2, 3, 4, 5,
                                                                                 6, 7},
                                                                      false, "\t", false, true);
        Hashtable<String, String> keeps = HashVec.loadFileToHashString(keepSamples, new int[] {0},
                                                                       new int[] {0, 1, 2, 3, 4, 5,
                                                                                  6, 7},
                                                                       false, "\t", false, true);

        Hashtable<String, String> hash = HashVec.loadFileToHashString(phenoFile, new int[] {0},
                                                                      new int[] {0, 1, 2, 3, 4, 5,
                                                                                 6, 7},
                                                                      false, "\t", false, true);

        HashMap<String, ArrayList<String>> centers = new HashMap<>();
        for (String sample : hash.keySet()) {
          String center = hash.get(sample).split("\t")[7];
          if (!center.equals("-1")) {
            if (!centers.containsKey(center)) {
              centers.put(center, new ArrayList<>());
            }
            if (have.containsKey(sample) && keeps.containsKey(sample)) {
              centers.get(center).add(hash.get(sample));
            }
          }
        }

        int min = i + 1;
        for (String cent : centers.keySet()) {
          System.out.println(cent + "\t" + centers.get(cent).size());

          if (centers.get(cent).size() < min) {
            min = centers.get(cent).size();
          }
        }
        System.out.println(min);

        ArrayList<String> qqFiles = new ArrayList<>();
        HashSet<String> comps = new HashSet<>();
        for (String cent1 : centers.keySet()) {
          System.out.println(cent1);
          for (int j = 0; j < reps; j++) {

            System.out.println(centers.get(cent1).size());

            if (cent1.startsWith("W")) {
              comps.add(cent1);
              List<String> caseSelect = selectCases(centers.get(cent1), hash, cent1, have, min,
                                                    threshold);
              // System.exit(1);
              ArrayList<String> vAllBarn = new ArrayList<>();
              for (String cent2 : centers.keySet()) {
                if (!cent1.equals(cent2)) {
                  vAllBarn.addAll(centers.get(cent2));
                }
              }
              runit(lambdaReport, caseSelect, i, rootOut, mds, plinkQc, hash, centers, min, qqFiles,
                    cent1, cent1 + "_v_ALL", vAllBarn, pheno, j, threshold);
              // if (i < 1100) {
              // for (String cent2 : centers.keySet()) {
              // String centComp = cent1 + "_v_" + cent2;
              //
              // if (!cent1.equals(cent2) &&
              // !comps.contains(centComp)
              // && !cent2.startsWith("F")) {
              // comps.add(centComp);
              // comps.add(cent2 + "_v_" + cent1);
              //
              // ArrayList<String> barns = new
              // ArrayList<>();
              // barns.addAll(centers.get(cent2));
              // runit(lambdaReport, caseSelect, i,
              // rootOut, mds,
              // plinkQc, hash, centers, min,
              // qqFiles, cent1, centComp, barns,
              // pheno,j);
              // }
              //
              // }
              // }
            }
          }

          System.out.println(ArrayUtils.toStr(qqFiles));
          // QQPlot qqPlot =
          // QQPlot.loadPvals(ArrayUtils.toStringArray(qqFiles),
          // "Hi",
          // false, true, false, -1, false,
          // Float.MAX_VALUE, new Logger());
          // qqPlot.screenCap(rootOut + "test_pheno" + pheno +
          // "_" + i
          // +
          // "samples.png");
          // Files.copyFileUsingFileChannels(rootOut + "test_pheno" +
          // pheno + "_" + i + "samples.png",
          // root + "test_pheno" + pheno + "_" + i + "samples.png",
          // new Logger());
        }
        Files.writeIterable(lambdaReport, root + threshold + "lambdaReport.txt");

      }

      // }).start();

      // }
    }

  }

  private static void runit(ArrayList<String> lambdaReport, List<String> casesL, int i,
                            String rootOut, String mds, String plinkQc,
                            Hashtable<String, String> hash,
                            HashMap<String, ArrayList<String>> centers, int min,
                            ArrayList<String> qqFiles, String cent1, String centComp,
                            ArrayList<String> barns, String pheno, int rep, double threshold) {
    String matchDir = rootOut + "match/" + centComp + "/";

    String[] run = new String[] {"C1", "C2"};
    double[] w = new double[] {1, 1};
    String antiOptimalDir = rootOut + "antimatch/" + centComp + "/";
    new File(antiOptimalDir).mkdirs();
    if (!Files.exists(antiOptimalDir + "mds20.mds")) {
      Files.copyFileUsingFileChannels(mds, antiOptimalDir + "mds20.mds", new Logger());
    }

    String pairsAnti = selectControls(casesL, centComp, barns, run, w, antiOptimalDir, true);

    String[] casesAnti = HashVec.loadFileToStringArray(pairsAnti, true, new int[] {0}, true);
    String[] controlsAnti = HashVec.loadFileToStringArray(pairsAnti, true, new int[] {1}, true);
    if (!Files.exists(antiOptimalDir + "mds20.mds")) {
      Files.copyFileUsingFileChannels(mds, antiOptimalDir + "mds20.mds", new Logger());
    }
    run(plinkQc, mds, hash, qqFiles, centComp, antiOptimalDir, "antioptimal", casesAnti,
        controlsAnti, lambdaReport, i, pheno, rep, threshold);

    new File(matchDir).mkdirs();
    if (!Files.exists(matchDir + "mds20.mds")) {

      Files.copyFileUsingFileChannels(mds, matchDir + "mds20.mds", new Logger());
    }
    String pairs = selectControls(casesL, centComp, barns, run, w, matchDir, false);

    String[] cases = HashVec.loadFileToStringArray(pairs, true, new int[] {0}, true);
    String[] controls = HashVec.loadFileToStringArray(pairs, true, new int[] {1}, true);

    run(plinkQc, mds, hash, qqFiles, centComp, matchDir, "optimal", cases, controls, lambdaReport,
        i, pheno, rep, threshold);

    String matchDirRandom = rootOut + "randomMatch/" + centComp + "/";
    new File(matchDirRandom).mkdirs();
    Files.copyFileUsingFileChannels(mds, matchDirRandom + "mds20.mds", new Logger());
    Collections.shuffle(barns);
    List<String> sub = barns.subList(0, Math.min(min, cases.length) - 1);
    Collections.sort(sub);
    run(plinkQc, mds, hash, qqFiles, centComp, matchDirRandom, "random", cases,
        ArrayUtils.toStringArray(sub), lambdaReport, i, pheno, rep, threshold);
  }

  private static String selectControls(List<String> casesL, String centComp,
                                       ArrayList<String> barns, String[] run, double[] w,
                                       String dir, boolean anti) {
    String[] previous = Files.toFullPaths(Files.list(dir, "distance", null, false), dir);
    for (String file : previous) {
      new File(file).delete();
    }
    Files.writeIterable(casesL, dir + centComp + ".anchors.txt");
    Files.writeIterable(barns, dir + centComp + ".barns.txt");
    System.out.println("RUNNING match1");
    String matchFileanti = MatchSamples.matchMaker(dir, centComp + ".anchors.txt",
                                                   centComp + ".barns.txt",
                                                   ext.removeDirectoryInfo("mds20.mds"), run, w,
                                                   false);
    System.out.println("RUNNING match2");
    matchFileanti = MatchSamples.normalizeDistances(dir, matchFileanti, 0, 100);
    System.out.println("RUNNING match3");

    String pairsAnti = dir + MatchSamples.matchPairs(dir, matchFileanti, true, anti);
    System.out.println("RUNNING match4");
    return pairsAnti;
  }

  private static void run(String plinkQc, String mds, Hashtable<String, String> hash,
                          ArrayList<String> qqFiles, String cent, String matchDir, String tag,
                          String[] cases, String[] controls, ArrayList<String> lambdaReport, int i,
                          String pheno, int rep, double threshold) {

    ArrayList<String> newFam = new ArrayList<>();
    HashSet<String> cas = new HashSet<>();
    HashSet<String> cont = new HashSet<>();

    if (!Files.exists(matchDir + "plink.bed")) {
      ArrayList<String> command = new ArrayList<>();
      command.add("plink2");
      command.add("--bfile");
      command.add(plinkQc);
      // command.add("--keep");
      // command.add(matchDir + "keeps.txt");
      command.add("--make-bed");

      command.add("--out");

      command.add(matchDir + "plink");
      CmdLine.run(command, matchDir, null, null, new Logger(), false);
    }

    for (String sample : cases) {
      String[] c = hash.get(sample.split("\t")[0]).split("\t");
      c[5] = "2";
      newFam.add(ArrayUtils.toStr(ArrayUtils.subArray(c, 0, 6)));
      cas.add(sample.split("\t")[0]);
    }
    for (String sample : controls) {
      // System.out.println(sample);
      String[] c = hash.get(sample.split("\t")[0]).split("\t");
      c[5] = "1";
      newFam.add(ArrayUtils.toStr(ArrayUtils.subArray(c, 0, 6)));
      cont.add(sample.split("\t")[0]);
    }
    Files.writeIterable(newFam, matchDir + "keeps.txt");
    // double maf = 40.0 / (controls.length + cases.length);
    double maf = 0.05;
    String baseFam = matchDir + "plinkBase.fam";
    if (!Files.exists(baseFam)) {
      Files.copyFile(matchDir + "plink.fam", baseFam);
    }
    String[] ins = HashVec.loadFileToStringArray(baseFam, false, null, false);

    newFam = new ArrayList<>();
    HashSet<String> have = new HashSet<>();
    for (String in : ins) {
      have.add(in);
    }
    for (String sample : ins) {
      if (cas.contains(sample.split(PSF.Regex.GREEDY_WHITESPACE)[0])) {
        String[] c = hash.get(sample.split(PSF.Regex.GREEDY_WHITESPACE)[0]).split("\t");
        c[5] = "2";
        newFam.add(ArrayUtils.toStr(ArrayUtils.subArray(c, 0, 6)));
      } else if (cont.contains(sample.split(PSF.Regex.GREEDY_WHITESPACE)[0])) {
        String[] c = hash.get(sample.split(PSF.Regex.GREEDY_WHITESPACE)[0]).split("\t");
        c[5] = "1";
        newFam.add(ArrayUtils.toStr(ArrayUtils.subArray(c, 0, 6)));
      } else {
        String[] c = hash.get(sample.split(PSF.Regex.GREEDY_WHITESPACE)[0]).split("\t");
        c[5] = "-9";
        newFam.add(ArrayUtils.toStr(ArrayUtils.subArray(c, 0, 6)));
      }

    }

    // for (String sample : controls) {
    // // System.out.println(sample);
    // if (have.contains(sample)) {
    // String[] c = hash.get(sample.split("\t")[0]).split("\t");
    // c[5] = "1";
    // newFam.add(ArrayUtils.toStr(ArrayUtils.subArray(c, 0, 6)));
    // }
    // }
    // Collections.sort(newFam);
    Files.writeIterable(newFam, matchDir + "plink.fam");

    // if (!Files.exists(matchDir + "plink.assoc.logistic.tabs")) {

    String outCovar = matchDir + cent + tag + rep + "covarTrue.p.txt";
    runPlink(mds, matchDir, maf, outCovar, true);
    writeReport(cent, tag, cases, lambdaReport, pheno, rep, threshold, outCovar, true, maf);

    String outNoCovar = matchDir + cent + tag + rep + "covarFalse.p.txt";
    runPlink(mds, matchDir, maf, outNoCovar, false);
    writeReport(cent, tag, cases, lambdaReport, pheno, rep, threshold, outNoCovar, false, maf);
    // System.exit(1);
    Files.copyFile(matchDir + "plink.fam", matchDir + "plink.fam_" + rep);

    // qqFiles.add(out);

    // ArrayUtils.lambda(pvals[i])
  }

  private static void writeReport(String cent, String tag, String[] cases,
                                  ArrayList<String> lambdaReport, String pheno, int rep,
                                  double threshold, String out, boolean hadCovars, double maf) {
    String[] ps = HashVec.loadFileToStringArray(out, true, new int[] {0}, false);
    ArrayList<Double> p = new ArrayList<>();
    for (String a : ps) {
      try {
        p.add(Double.parseDouble(a));
      } catch (NumberFormatException nfe) {

      }
    }
    String report = cent + "\t" + tag + "\t" + cases.length + "\t" + rep + "\t" + threshold + "\t"
                    + ArrayUtils.lambda(Doubles.toArray(p)) + "\t" + pheno + "\t" + hadCovars + "\t"
                    + maf;
    lambdaReport.add(report);
    Files.write(report, out + ".lambda");
  }

  private static void runPlink(String mds, String matchDir, double maf, String out, boolean covar) {
    ArrayList<String> plink = new ArrayList<>();

    plink.add("plink2");

    plink.add("--logistic");

    plink.add("--maf");
    plink.add(maf + "");
    plink.add("--bfile");
    plink.add(matchDir + "plink");
    // plink.add(plinkQc);
    plink.add("--keep");
    plink.add(matchDir + "keeps.txt");
    plink.add("--geno");
    plink.add(".01");
    plink.add("--mind");
    plink.add(".05");
    plink.add("--out");
    plink.add(matchDir + "plink");
    plink.add("--threads");
    plink.add(4 + "");
    if (covar) {
      plink.add("--covar");
      plink.add(mds);
      plink.add("--covar-name");
      plink.add("C1,C2");
    }
    // System.out.println(maf);
    CmdLine.run(plink, matchDir, null, null, new Logger(), false);
    // System.exit(1);

    // }
    String t1 = "cat " + matchDir
                + "plink.assoc.logistic |grep \"ADD\\|SNP\"| tr -s ' ' '\t'|cut -f 2- > " + matchDir
                + "plink.assoc.logistic.tabs";
    String t2 = "cut -f9 " + matchDir + "plink.assoc.logistic.tabs > " + out;
    org.pankratzlab.common.Files.write(t1 + "\n" + t2, matchDir + "tab.sh");
    org.pankratzlab.common.Files.chmod(matchDir + "tab.sh");
    CmdLine.run(matchDir + "tab.sh", matchDir);
  }

}
