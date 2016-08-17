package org.genvisis.gaw17;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;
import org.genvisis.mining.Transformations;
import org.genvisis.stats.LeastSquares;
import org.genvisis.stats.LogisticRegression;

public class Traits {
  // public static final String[] POPS = {"Caucasians", "Asians", "Africans"};
  // public static final String[] POPS = {"Caucasians", "Asians", "Africans", "All"};
  // public static final String[] POPS = {"all697"};
  public static final String[] POPS =
      {"CEPH", "Denver", "Han", "Japanese", "Luhya", "Tuscan", "Yoruban"};
  // public static final String[] PHENOS = {"aff", "q1", "q2", "q4"};
  public static final String[] PHENOS = {"q1"};
  // public static final String[] PHENOS = {"affected", "q1", "q2", "q4", "AffResids",
  // "Aff_minusQ1Q4", "Aff_minusQ1Q2Q4", "Q1resids", "Q1residsMinusQ2", "Q2residsQ1", "Q4resids"};
  // public static final String[] PHENOS = {"AffResids", "Aff_minusQ1Q4", "Aff_minusQ1Q2Q4",
  // "Q1resids", "Q1residsMinusQ2", "Q2residsQ1", "Q4resids"};
  // public static final String[][] COVARS = {{"SEX", "AGE", "SMOKE", "C1", "C2"}, {"SEX", "AGE",
  // "SMOKE", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"}};
  public static final String[][] COVARS =
      {{"SEX", "SMOKE"}, {"SEX", "SMOKE"}, {}, {"SEX", "AGE", "SMOKE"}};

  public static void analyzeCovariates(String analysis_dir, String pheno_dir) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    int count;
    long time;
    double[][] deps, indeps;

    time = new Date().getTime();
    try {
      writer = new PrintWriter(new FileWriter(analysis_dir + "allReps.xls"));
      writer.println("Rep\tAff_age\tAff_gender\tAff_smoking\tQ1_age\tQ1_gender\tQ1_smoking\tQ2_age\tQ2_gender\tQ2_smoking\tQ4_age\tQ4_gender\tQ4_smoking");
      for (int rep = 1; rep <= 200; rep++) {
        writer.print(rep);
        try {
          reader = new BufferedReader(new FileReader(pheno_dir + "unr_phen." + rep));
          if (!reader.readLine().equals("ID,SEX,AGE,SMOKE,Q1,Q2,Q4,AFFECTED")) {
            System.err.println("Error - invalid header for replicate: " + rep);
          }
          deps = new double[4][697];
          indeps = new double[697][3];
          count = 0;
          while (reader.ready()) {
            line = reader.readLine().trim().split(",");
            deps[0][count] = Double.parseDouble(line[7]);
            deps[1][count] = Double.parseDouble(line[4]);
            deps[2][count] = Double.parseDouble(line[5]);
            deps[3][count] = Double.parseDouble(line[6]);
            indeps[count][0] = Double.parseDouble(line[2]);
            indeps[count][1] = Double.parseDouble(line[1]);
            indeps[count][2] = Double.parseDouble(line[3]);
            count++;
          }
          writer.print("\t" + Array.toStr(Array.subArray(
                                                         new LogisticRegression(deps[0], indeps)
                                                                                                .getSigs(),
                                                         1, 4)));
          writer.print("\t"
                       + Array.toStr(Array.subArray(new LeastSquares(deps[1], indeps).getSigs(), 1,
                                                    4)));
          writer.print("\t"
                       + Array.toStr(Array.subArray(new LeastSquares(deps[2], indeps).getSigs(), 1,
                                                    4)));
          writer.print("\t"
                       + Array.toStr(Array.subArray(new LeastSquares(deps[3], indeps).getSigs(), 1,
                                                    4)));
          writer.println();
          writer.flush();
          reader.close();
        } catch (FileNotFoundException fnfe) {
          System.err.println("Error: file \"" + "unr_phen." + rep
                             + "\" not found in current directory");
          System.exit(1);
        } catch (IOException ioe) {
          System.err.println("Error reading file \"" + "unr_phen." + rep + "\"");
          System.exit(2);
        }
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + analysis_dir + "allReps.xls");
      e.printStackTrace();
    }
    System.out.println("Finished in " + ext.getTimeElapsed(time));
  }

  public static void analyzeCovariates2(String analysis_dir, String pheno_dir) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    int count;
    long time;
    double[][] deps;
    double[][][] indeps;

    time = new Date().getTime();
    try {
      writer = new PrintWriter(new FileWriter(analysis_dir + "allReps2.xls"));
      writer.println("Rep\tAff_age\tAff_gender\tAff_smoking\tAff_Q1\tAff_Q2\tAff_Q4\tQ1_age\tQ1_gender\tQ1_smoking\tQ1_Q2\tQ1_Q4\tQ2_age\tQ2_gender\tQ2_smoking\tQ2_Q1\tQ2_Q4\tQ4_age\tQ4_gender\tQ4_smoking\tQ4_Q1\tQ4_Q2");
      for (int rep = 1; rep <= 200; rep++) {
        writer.print(rep);
        try {
          reader = new BufferedReader(new FileReader(pheno_dir + "unr_phen." + rep));
          if (!reader.readLine().equals("ID,SEX,AGE,SMOKE,Q1,Q2,Q4,AFFECTED")) {
            System.err.println("Error - invalid header for replicate: " + rep);
          }
          deps = new double[4][697];
          indeps = new double[4][697][];
          count = 0;
          while (reader.ready()) {
            line = reader.readLine().trim().split(",");
            deps[0][count] = Double.parseDouble(line[7]);
            deps[1][count] = Double.parseDouble(line[4]);
            deps[2][count] = Double.parseDouble(line[5]);
            deps[3][count] = Double.parseDouble(line[6]);
            indeps[0][count] = new double[6];
            indeps[1][count] = new double[5];
            indeps[2][count] = new double[5];
            indeps[3][count] = new double[5];
            for (int i = 0; i < 4; i++) {
              indeps[i][count][0] = Double.parseDouble(line[2]);
              indeps[i][count][1] = Double.parseDouble(line[1]);
              indeps[i][count][2] = Double.parseDouble(line[3]);
            }
            indeps[0][count][3] = Double.parseDouble(line[4]);
            indeps[0][count][4] = Double.parseDouble(line[5]);
            indeps[0][count][5] = Double.parseDouble(line[6]);
            indeps[1][count][3] = Double.parseDouble(line[5]);
            indeps[1][count][4] = Double.parseDouble(line[6]);
            indeps[2][count][3] = Double.parseDouble(line[4]);
            indeps[2][count][4] = Double.parseDouble(line[6]);
            indeps[3][count][3] = Double.parseDouble(line[4]);
            indeps[3][count][4] = Double.parseDouble(line[5]);
            count++;
          }
          writer.print("\t" + Array.toStr(Array.subArray(
                                                         new LogisticRegression(deps[0], indeps[0])
                                                                                                   .getSigs(),
                                                         1, 7)));
          writer.print("\t"
                       + Array.toStr(Array.subArray(new LeastSquares(deps[1], indeps[1]).getSigs(),
                                                    1, 6)));
          writer.print("\t"
                       + Array.toStr(Array.subArray(new LeastSquares(deps[2], indeps[2]).getSigs(),
                                                    1, 6)));
          writer.print("\t"
                       + Array.toStr(Array.subArray(new LeastSquares(deps[3], indeps[3]).getSigs(),
                                                    1, 6)));
          writer.println();
          writer.flush();
          reader.close();
        } catch (FileNotFoundException fnfe) {
          System.err.println("Error: file \"" + "unr_phen." + rep
                             + "\" not found in current directory");
          System.exit(1);
        } catch (IOException ioe) {
          System.err.println("Error reading file \"" + "unr_phen." + rep + "\"");
          System.exit(2);
        }
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + analysis_dir + "allReps.xls");
      e.printStackTrace();
    }
    System.out.println("Finished in " + ext.getTimeElapsed(time));
  }

  // convert to tab-delimited files
  public static void convert(String pheno_dir) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;

    for (int rep = 1; rep <= 200; rep++) {
      try {
        reader = new BufferedReader(new FileReader(pheno_dir + "unr_phen." + rep));
        writer = new PrintWriter(new FileWriter(pheno_dir + "unr_phen." + rep + ".tab"));
        writer.println("FID\tI" + Array.toStr(reader.readLine().trim().split(",")));
        while (reader.ready()) {
          line = reader.readLine().trim().split(",");
          writer.println(line[0] + "\t" + line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3]
                         + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t"
                         + (line[7].equals("1") ? "2" : "1"));
        }
        writer.close();
        reader.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + pheno_dir + "unr_phen." + rep
                           + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + pheno_dir + "unr_phen." + rep + "\"");
        System.exit(2);
      }
    }

  }

  public static void analyzeUnivariates(String analysis_dir, String pheno_dir) {
    PrintWriter writer;
    Logger log;

    try {
      writer = new PrintWriter(new FileWriter(analysis_dir + "analysisUnivariate.bat"));
      for (int rep = 1; rep <= 200; rep++) {
        for (int i = 0; i < POPS.length; i++) {
          // writer.println("plink --bfile ../"+POPS[i]+" --pheno "+pheno_dir+"unr_phen."+rep+".tab
          // --pheno-name AFFECTED --covar "+pheno_dir+"unr_phen."+rep+".tab --covar-name SEX,SMOKE
          // --logistic --ci 0.95 --out "+POPS[i]+".aff."+rep+"");
          writer.println("plink --bfile ../" + POPS[i] + " --pheno " + pheno_dir + "unr_phen." + rep
                         + ".tab --pheno-name Q1 --covar " + pheno_dir + "unr_phen." + rep
                         + ".tab --covar-name SEX,SMOKE --logistic --ci 0.95 --out " + POPS[i]
                         + ".q1." + rep + "");
          // writer.println("plink --bfile ../"+POPS[i]+" --pheno "+pheno_dir+"unr_phen."+rep+".tab
          // --pheno-name Q2 --logistic --ci 0.95 --out "+POPS[i]+".q2."+rep+"");
          // writer.println("plink --bfile ../"+POPS[i]+" --pheno "+pheno_dir+"unr_phen."+rep+".tab
          // --pheno-name Q4 --covar "+pheno_dir+"unr_phen."+rep+".tab --covar-name SEX,AGE,SMOKE
          // --logistic --ci 0.95 --out "+POPS[i]+".q4."+rep+"");
        }
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + analysis_dir + "analysisUnivariate.bat");
      e.printStackTrace();
    }

    log = new Logger(analysis_dir + "RunMetal.bat");
    for (int i = 0; i < POPS.length; i++) {
      for (int pheno = 0; pheno < PHENOS.length; pheno++) {
        try {
          writer = new PrintWriter(new FileWriter(analysis_dir + "metal_" + POPS[i] + "_"
                                                  + PHENOS[pheno] + ".txt"));
          writer.println("MARKER SNP");
          writer.println("ALLELE A1 TEST");
          writer.println("EFFECT " + (!PHENOS[pheno].startsWith("q") ? "OR" : "BETA"));
          writer.println("STDERR SE");
          writer.println("SCHEME STDERR");
          writer.println("ADDFILTER TEST IN (ADD)");
          writer.println("");
          for (int rep = 1; rep <= 200; rep++) {
            writer.println("PROCESS " + POPS[i] + ".aff." + rep + ".assoc."
                           + (!PHENOS[pheno].startsWith("q") ? "logistic" : "linear"));
          }
          writer.println("");
          writer.println("OUTFILE meta_" + POPS[i] + " .tbl");
          writer.println("ANALYZE");

          writer.close();
        } catch (Exception e) {
          System.err.println("Error writing to " + analysis_dir + "metal_" + POPS[i] + "_"
                             + PHENOS[pheno] + ".txt");
          e.printStackTrace();
        }
        log.report("metal < metal_" + POPS[i] + "_" + PHENOS[pheno] + ".txt");
      }
    }

    log = new Logger(analysis_dir + "RunMetals.bat");
    for (int pheno = 0; pheno < PHENOS.length; pheno++) {
      for (int rep = 1; rep <= 200; rep++) {
        try {
          writer = new PrintWriter(new FileWriter(analysis_dir + "metal_" + PHENOS[pheno] + "."
                                                  + rep + ".txt"));
          writer.println("MARKER SNP");
          writer.println("ALLELE A1 TEST");
          writer.println("EFFECT " + (!PHENOS[pheno].startsWith("q") ? "OR" : "BETA"));
          writer.println("STDERR SE");
          writer.println("SCHEME STDERR");
          writer.println("ADDFILTER TEST IN (ADD)");
          writer.println("");
          for (int i = 0; i < POPS.length; i++) {
            writer.println("PROCESS " + POPS[i] + "." + PHENOS[pheno] + "." + rep + ".assoc."
                           + (!PHENOS[pheno].startsWith("q") ? "logistic" : "linear"));
          }
          writer.println("");
          writer.println("OUTFILE meta_" + rep + " .tbl");
          writer.println("ANALYZE");

          writer.close();
        } catch (Exception e) {
          System.err.println("Error writing to " + analysis_dir + "metal_" + "_" + PHENOS[pheno]
                             + "." + rep + ".txt");
          e.printStackTrace();
        }
        log.report("metal < metal_" + PHENOS[pheno] + "." + rep + ".txt");
      }
    }
  }

  public static void metaAnalyze(String analysis_dir) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    int count;
    String[] snps;
    Hashtable<String, String> hash;

    for (int rep = 1; rep <= 200; rep++) {
      try {
        writer = new PrintWriter(new FileWriter(analysis_dir + rep + ".metal"));
        writer.println("MARKER SNP");
        writer.println("ALLELE A1 A2");
        // writer.println("EFFECT beta");
        // writer.println("STDERR se");
        // writer.println("SCHEME STDERR");
        // writer.println("GENOMICCONTROL ON");
        writer.println("EFFECT beta");
        writer.println("PVALUE pval");
        writer.println("WEIGHT weight");
        // writer.println("GENOMICCONTROL ON");
        writer.println("");
        for (int i = 0; i < POPS.length; i++) {
          writer.println("PROCESS " + POPS[i] + ".q1." + rep + ".in");
        }
        writer.println("");
        writer.println("");
        writer.println("OUTFILE " + rep + ".meta .tbl");
        writer.println("ANALYZE");
        writer.close();
      } catch (Exception e) {
        System.err.println("Error writing to " + rep + ".metal");
        e.printStackTrace();
      }
    }
    try {
      writer = new PrintWriter(new FileWriter(analysis_dir + "runMetals.bat"));
      for (int rep = 1; rep <= 200; rep++) {
        writer.println("metal < " + rep + ".metal > " + rep + ".log");
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + "metal.bat");
      e.printStackTrace();
    }
    System.exit(1);

    snps = HashVec.loadFileToStringArray(analysis_dir + "snps.txt", false, new int[] {0}, false);

    for (int i = 0; i < POPS.length; i++) {
      hash = HashVec.loadFileToHashString(analysis_dir + POPS[i] + ".frq", 1, new int[] {2, 3},
                                          "\t", true);
      System.out.print(POPS[i]);
      for (int rep = 1; rep <= 200; rep++) {
        if (rep % 10 == 0) {
          System.out.print(".");
        }
        for (int pheno = 0; pheno < PHENOS.length; pheno++) {
          try {
            reader =
                new BufferedReader(new FileReader(analysis_dir + POPS[i] + "." + PHENOS[pheno] + "."
                                                  + rep + ".assoc."
                                                  + (!PHENOS[pheno].startsWith("q") ? "logistic"
                                                                                    : "linear")));
            writer = new PrintWriter(new FileWriter(analysis_dir + POPS[i] + "." + PHENOS[pheno]
                                                    + "." + rep + ".in"));
            // writer.println("SNP\tA1\tA2\tbeta\tse");
            writer.println("SNP\tA1\tA2\tbeta\tpval\tweight");
            reader.readLine();
            count = 0;
            while (reader.ready()) {
              line = reader.readLine().trim().split("[\\s]+");
              if (!line[1].equals(snps[count])) {
                System.err.println("Error - out of sync in " + POPS[i] + ".aff." + rep + ".assoc."
                                   + (!PHENOS[pheno].startsWith("q") ? "logistic" : "linear")
                                   + " expecting " + snps[count] + " found " + line[1]);
              }
              if (!line[11].equals("NA")) {
                String alleles;
                alleles = hash.get(line[1]);
                if (alleles == null) {
                  System.err.println("Error - no frq for " + line[1] + " in " + POPS[i]);
                } else if (!alleles.startsWith(line[3])) {
                  System.err.println("Error - wrong A1 for " + line[1] + " found " + line[3]
                                     + " expecting " + alleles);
                }
                // writer.println(line[1]+"\t"+alleles+"\t"+line[6]+"\t"+line[7]);
                // writer.println(line[1]+"\t"+alleles+"\t"+line[6]+"\t"+line[11]+"\t"+(Double.parseDouble(line[6])<0?"-":"+")+"\t"+line[5]);
                writer.println(line[1] + "\t" + alleles + "\t" + line[6] + "\t" + line[11] + "\t"
                               + "\t" + line[5]);
              }
              for (int j = 0; j < COVARS[pheno].length; j++) {
                line = reader.readLine().trim().split("[\\s]+");
                if (!line[4].equals(COVARS[pheno][j])) {
                  System.err.println("Error - out of sync in " + POPS[i] + ".aff." + rep + ".assoc."
                                     + (!PHENOS[pheno].startsWith("q") ? "logistic" : "linear")
                                     + " expecting " + snps[count] + " found " + line[1]);
                }
              }
              count++;
            }
            reader.close();
            writer.close();
          } catch (FileNotFoundException fnfe) {
            System.err.println("Error: file \"" + analysis_dir + POPS[i] + "." + PHENOS[pheno] + "."
                               + rep + ".assoc."
                               + (!PHENOS[pheno].startsWith("q") ? "logistic" : "linear")
                               + "\" not found in current directory");
            System.exit(1);
          } catch (IOException ioe) {
            System.err.println("Error reading file \"" + analysis_dir + POPS[i] + "."
                               + PHENOS[pheno] + "." + rep + ".assoc."
                               + (!PHENOS[pheno].startsWith("q") ? "logistic" : "linear") + "\"");
            System.exit(2);
          }
        }
      }
      System.out.println();
    }
  }

  public static void addMDS(String pheno_dir, String mds_file, int numComps) {
    String comps = " 0";

    for (int i = 1; i <= numComps; i++) {
      comps += " " + (i + 2) + "=C" + i;
    }

    for (int rep = 1; rep <= 200; rep++) {
      // Files.combine(HashVec.loadFileToStringArray(pheno_dir+"unr_phen."+rep, false, true, new
      // int[] {0}, true, false, ","), new String[] {pheno_dir+"unr_phen."+rep+" 0 0=IID 4=Q1 1=SEX
      // 2=AGE 3=SMOKE ,", mds_file+" 0 3=C1 4=C2"}, "FID", pheno_dir+"pheno_C2."+rep, new Logger(),
      // false);
      // Files.combine(HashVec.loadFileToStringArray(pheno_dir+"unr_phen."+rep, false, true, new
      // int[] {0}, true, false, ","), new String[] {pheno_dir+"unr_phen."+rep+" 0 0=IID 4=Q1 1=SEX
      // 2=AGE 3=SMOKE ,", mds_file+" 0 3=C1 4=C2 5=C3 6=C4 7=C5 8=C6 9=C7 10=C8 11=C9 12=C10"},
      // "FID", pheno_dir+"pheno_C10."+rep, new Logger(), false);
      Files.combine(HashVec.loadFileToStringArray(pheno_dir + "unr_phen." + rep, false, true,
                                                  new int[] {0}, true, false, ","),
                    new String[] {pheno_dir + "unr_phen." + rep
                                  + " 0 0=IID 7=Affected 4=Q1 5=Q2 1=SEX 2=AGE 3=SMOKE ,",
                                  mds_file + comps},
                    "FID", pheno_dir + "pheno_C" + numComps + "." + rep, new Logger(), false);
    }
  }

  public static void normPheno(String pheno_dir, String output_dir, String filename_root,
                               int column, int normType) {
    output_dir = ext.verifyDirFormat(output_dir);
    new File(output_dir).mkdirs();
    for (int rep = 1; rep <= 200; rep++) {
      Transformations.transformFile(pheno_dir + filename_root + "." + rep,
                                    output_dir + filename_root + "." + rep, true, column, false,
                                    true, normType, new Logger());
    }
  }

  public static void analyzeMDS(String rootPheno, int numComps, String trait) {
    String comps = "";

    for (int i = 1; i <= numComps; i++) {
      comps += ",C" + i;
    }

    // Files.qsub("C2", 1, 200, "/home/npankrat/bin/plink --bfile all697 --pheno phenos/pheno_C10.#
    // --pheno-name Q1 --covar phenos/pheno_C10.# --covar-name SEX,AGE,SMOKE,C1,C2 --linear --ci
    // 0.95 --out C2.#");
    // Files.qsub("C10", 1, 200, "/home/npankrat/bin/plink --bfile all697 --pheno phenos/pheno_C10.#
    // --pheno-name Q1 --covar phenos/pheno_C10.# --covar-name
    // SEX,AGE,SMOKE,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10 --linear --ci 0.95 --out C10.#");
    Files.qsub("C" + numComps, 1, 200,
               "/home/npankrat/bin/plink --bfile all697 --pheno phenos/" + rootPheno
                                       + ".# --pheno-name " + trait + " --covar phenos/" + rootPheno
                                       + ".# --covar-name SEX,AGE,SMOKE" + comps
                                       + " --linear --ci 0.95 --out C" + numComps + ".#");
  }

  public static void sumSigs(String analysis_dir, double cutoff) {
    BufferedReader reader;
    String[] line;
    int count;
    String[] snps;
    int[][] counts;

    snps = HashVec.loadFileToStringArray(analysis_dir + "snps.txt", false, new int[] {0}, false);
    counts = new int[snps.length][8];

    for (int i = 0; i < POPS.length; i++) {
      System.out.print(POPS[i]);
      for (int rep = 1; rep <= 200; rep++) {
        if (rep % 10 == 0) {
          System.out.print(".");
        }
        for (int pheno = 0; pheno < PHENOS.length; pheno++) {
          try {
            reader =
                new BufferedReader(new FileReader(analysis_dir + POPS[i] + "." + PHENOS[pheno] + "."
                                                  + rep + ".assoc."
                                                  + (!PHENOS[pheno].startsWith("q") ? "logistic"
                                                                                    : "linear")));
            reader.readLine();
            count = 0;
            while (reader.ready()) {
              line = reader.readLine().trim().split("[\\s]+");
              if (!line[1].equals(snps[count])) {
                System.err.println("Error - out of sync in " + POPS[i] + ".aff." + rep + ".assoc."
                                   + (!PHENOS[pheno].startsWith("q") ? "logistic" : "linear")
                                   + " expecting " + snps[count] + " found " + line[1]);
              }
              if (!line[11].equals("NA")) {
                if (Double.parseDouble(line[11]) < cutoff) {
                  counts[count][pheno * 2 + 0]++;
                }
                counts[count][pheno * 2 + 1]++;
              }
              for (int j = 0; j < COVARS[pheno].length; j++) {
                line = reader.readLine().trim().split("[\\s]+");
                if (!line[4].equals(COVARS[pheno][j])) {
                  System.err.println("Error - out of sync in " + POPS[i] + ".aff." + rep + ".assoc."
                                     + (!PHENOS[pheno].startsWith("q") ? "logistic" : "linear")
                                     + " expecting " + snps[count] + " found " + line[1]);
                }
              }
              count++;
            }
            reader.close();
          } catch (FileNotFoundException fnfe) {
            System.err.println("Error: file \"" + analysis_dir + POPS[i] + "." + PHENOS[pheno] + "."
                               + rep + ".assoc."
                               + (!PHENOS[pheno].startsWith("q") ? "logistic" : "linear")
                               + "\" not found in current directory");
            System.exit(1);
          } catch (IOException ioe) {
            System.err.println("Error reading file \"" + analysis_dir + POPS[i] + "."
                               + PHENOS[pheno] + "." + rep + ".assoc."
                               + (!PHENOS[pheno].startsWith("q") ? "logistic" : "linear") + "\"");
            System.exit(2);
          }
        }
      }
      System.out.println();
    }
  }

  public static void sumMDS(String prefix, double cutoff, String method) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    int count;
    String[] snps;
    int[][] counts;

    snps = HashVec.loadFileToStringArray("snps.txt", false, new int[] {0}, false);
    counts = new int[snps.length][8];

    for (int rep = 1; rep <= 200; rep++) {
      if (rep % 10 == 0) {
        System.out.print(".");
      }
      try {
        reader = new BufferedReader(new FileReader(prefix + "." + rep + ".assoc." + method));
        reader.readLine();
        count = 0;
        while (reader.ready()) {
          line = reader.readLine().trim().split("[\\s]+");
          if (line[4].equals("ADD")) {
            if (!line[1].equals(snps[count])) {
              System.err.println("Error - out of sync in " + prefix + "." + rep + ".assoc." + method
                                 + " expecting " + snps[count] + " found " + line[1]);
            }
            if (!line[11].equals("NA")) {
              if (Double.parseDouble(line[11]) < cutoff) {
                counts[count][0]++;
              }
              counts[count][1]++;
            }
            // for (int j = 0; j < COVARS[prefix.length()-2].length; j++) {
            // line = reader.readLine().trim().split("[\\s]+");
            // if (!line[4].equals(COVARS[prefix.length()-2][j])) {
            // System.err.println("Error - out of sync in "+prefix+"."+rep+".assoc."+method+"
            // expecting "+snps[count]+" found "+line[1]);
            // }
            // }
            count++;
          }
        }
        if (count != snps.length) {
          System.err.println("Error - only " + count + " observations for replicate " + rep);
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + prefix + "." + rep + ".assoc." + method
                           + "\" not found in current directory");
        // System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + prefix + "." + rep + ".assoc." + method
                           + "\"");
        // System.exit(2);
      }
    }
    try {
      writer = new PrintWriter(new FileWriter(prefix + "_sigCounts_" + cutoff + ".xln"));
      writer.println("SNP\tsigs\tcount\t%");
      for (int j = 0; j < snps.length; j++) {
        writer.println(snps[j] + "\t" + counts[j][0] + "\t" + counts[j][1] + "\t"
                       + ext.formDeci((double) counts[j][0] / (double) counts[j][1] * 100, 1)
                       + "%");
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + prefix + "_sigCounts.xln");
      e.printStackTrace();
    }
  }

  public static void sumMetal(String analysis_dir, double cutoff) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    // int count;
    String[] snps;
    CountVector cv;

    snps = HashVec.loadFileToStringArray(analysis_dir + "snps.txt", false, new int[] {0}, false);
    cv = new CountVector();

    for (int rep = 1; rep <= 200; rep++) {
      if (rep % 10 == 0) {
        System.out.print(".");
      }
      try {
        reader = new BufferedReader(new FileReader(analysis_dir + rep + ".meta1.tbl"));
        reader.readLine();
        // count = 0;
        while (reader.ready()) {
          line = reader.readLine().trim().split("[\\s]+");
          if (!line[5].equals("NA")) {
            if (Double.parseDouble(line[5]) < cutoff) {
              cv.add(line[0]);
            }
          }
          // count++;
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + rep + ".meta1.tbl\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + rep + ".meta1.tbl\"");
        System.exit(2);
      }
    }
    try {
      writer = new PrintWriter(new FileWriter(analysis_dir + "meta_sigCounts_" + cutoff + ".xln"));
      String[] values = cv.getValues();
      int[] counts = cv.getCounts();
      int index;
      writer.println("SNP\tQ1_sigs\tQ1_count\tQ1_%");
      for (int j = 0; j < snps.length; j++) {
        index = ext.indexOfStr(snps[j], values);
        writer.println(snps[j] + "\t" + (index == -1 ? 0 : counts[index]) + "\t200\t"
                       + ext.formDeci((double) (index == -1 ? 0 : counts[index]) / (double) 200
                                      * 100, 1)
                       + "%");
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + "meta_sigCounts.xln");
      e.printStackTrace();
    }
  }

  public static void splitPops(String analysis_dir) {
    for (int i = 0; i < POPS.length; i++) {
      CmdLine.run("plink --bfile all697 --keep data/" + POPS[i] + ".txt --make-bed --out "
                  + POPS[i], analysis_dir);
      CmdLine.run("plink --bfile " + POPS[i] + " --freq --out " + POPS[i], analysis_dir);
    }
  }

  public static void analyzeDistribution(String pheno_dir, String pheno_root, int col,
                                         boolean commaDelimited, double[] stddevThresholds,
                                         String output_file) {
    PrintWriter writer;
    String[] ids;
    double[] values;
    int[][] counts;
    double mean, stdev;
    double[][] matrix;

    ids = HashVec.loadFileToStringArray(pheno_dir + pheno_root + ".1", false, true, new int[] {0},
                                        true, false, commaDelimited ? "," : "[\\s]+");

    counts = new int[ids.length][stddevThresholds.length];
    matrix = new double[ids.length][200];
    for (int rep = 1; rep <= 200; rep++) {
      values =
          Array.toDoubleArray(Array.toStringArray(HashVec.loadFileToVec(pheno_dir + pheno_root + "."
                                                                        + rep, true,
                                                                        new int[] {col}, true,
                                                                        false, false,
                                                                        commaDelimited ? ","
                                                                                       : "[\\s]+")));
      mean = Array.mean(values);
      stdev = Array.stdev(values);
      for (int i = 0; i < ids.length; i++) {
        matrix[i][rep - 1] = values[i];
        for (int j = 0; j < stddevThresholds.length; j++) {
          if (Math.abs(values[i] - mean) > stddevThresholds[j] * stdev) {
            counts[i][j]++;
          }
        }
      }
    }

    try {
      writer = new PrintWriter(new FileWriter(output_file));
      writer.print("UID\tmeanValue\tstdev");
      for (int i = 0; i < stddevThresholds.length; i++) {
        writer.print("\t" + stddevThresholds[i] + " SD");
      }
      writer.println();
      for (int i = 0; i < ids.length; i++) {
        writer.println(ids[i] + "\t" + Array.mean(matrix[i]) + "\t" + Array.stdev(matrix[i]) + "\t"
                       + Array.toStr(counts[i]));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + output_file);
      e.printStackTrace();
    }

    try {
      writer = new PrintWriter(new FileWriter(ext.rootOf(output_file) + "_matrix.xln"));
      writer.print("UID");
      for (int i = 1; i <= 200; i++) {
        writer.print("\tRep" + i);
      }
      writer.println();
      for (int i = 0; i < ids.length; i++) {
        writer.println(ids[i] + "\t" + Array.toStr(matrix[i]));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + output_file);
      e.printStackTrace();
    }


  }

  public static void main(String[] args) {
    // String analysis_dir = "";
    // String analysis_dir = "D:\\GAW17\\analysis\\allUnivAll\\";
    // String analysis_dir = "D:\\GAW17\\analysis\\univariate\\";
    String pheno_dir = "D:\\GAW17\\source\\phenos_unr\\";
    // String mds_file = "D:\\GAW17\\mds10.mds";
    // String analysis_dir = "D:\\GAW17\\splitUni\\";

    try {
      // analyzeCovariates(analysis_dir, pheno_dir);
      // analyzeCovariates2(analysis_dir, pheno_dir);
      // convert(pheno_dir);
      // analyzeUnivariates(analysis_dir, pheno_dir);
      // sumSigs(analysis_dir, 0.0000000001);
      // sumSigs(analysis_dir, 0.0000001);
      // sumSigs(analysis_dir, 0.000001);
      // sumSigs(analysis_dir, 0.00001);
      // sumSigs(analysis_dir, 0.0001);
      // sumSigs(analysis_dir, 0.001);
      // addMDS(pheno_dir, mds_file);
      // analyzeMDS();

      // addMDS(pheno_dir, "D:\\GAW17\\mds50.mds", 50);
      // analyzeMDS("pheno_C50", 0, "Affected");
      // analyzeMDS("pheno_C50", 2);
      // analyzeMDS("pheno_C50", 10, "Affected");
      // analyzeMDS("pheno_C50", 20, "Q1");
      // analyzeMDS("pheno_C50", 30);
      // analyzeMDS("pheno_C50", 40);
      // analyzeMDS("pheno_C50", 50);

      // normPheno(pheno_dir, "normed", "pheno_C50", 3, Transformations.INVERSE_NORMALIZE);

      // sumMDS(args[0], Double.parseDouble(args[1]), args[2]);
      // sumMDS("C10", 0.00001);
      // metaAnalyze(analysis_dir);
      // sumMetal(analysis_dir, 0.00001);
      // splitPops(analysis_dir);

      // analyzeDistribution(pheno_dir, "unr_phen", 1, true, new double[] {1, 2, 3, 4, 5, 6, 7},
      // "Sex_distribution.xln");
      // analyzeDistribution(pheno_dir, "unr_phen", 2, true, new double[] {1, 2, 3, 4, 5, 6, 7},
      // "Age_distribution.xln");
      analyzeDistribution(pheno_dir, "unr_phen", 3, true, new double[] {1, 2, 3, 4, 5, 6, 7},
                          "Smoke_distribution.xln");
      analyzeDistribution(pheno_dir, "unr_phen", 4, true, new double[] {1, 2, 3, 4, 5, 6, 7},
                          "Q1_distribution.xln");
      analyzeDistribution(pheno_dir, "unr_phen", 5, true, new double[] {1, 2, 3, 4, 5, 6, 7},
                          "Q2_distribution.xln");
      analyzeDistribution(pheno_dir, "unr_phen", 6, true, new double[] {1, 2, 3, 4, 5, 6, 7},
                          "Q5_distribution.xln");
      analyzeDistribution(pheno_dir, "unr_phen", 7, true, new double[] {1, 2, 3, 4, 5, 6, 7},
                          "Aff_distribution.xln");
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
