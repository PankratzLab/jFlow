package org.genvisis.link;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.ext;

public class Merlin {
  public static final int CHR_START = 1;
  public static final int CHR_STOP = 22;
  // public static final int CHR_START = 16;
  // public static final int CHR_STOP = 16;
  public static final int UNKNOWN_TRAIT = 0;
  public static final int BINARY_TRAIT = 1;
  public static final int QUANTITATIVE_TRAIT = 2;

  public static void batchAllBinary(String qsub, boolean blade) {
    String commands;
    // PrintWriter writer;
    // String chrome;
    //
    // try {
    // writer = new PrintWriter(new FileWriter("MerlinAll.bat"));
    // for (int i = CHR_START; i<=CHR_STOP; i++) {
    // chrome = ext.chrome(i);
    //// writer.println("java -cp \"C:\\Documents and Settings\\npankrat\\My Documents\\jProjects\\"
    // + common.PSF.Java.GENVISIS + "\" -Xms1024M -Xmx1024M link.Merlin chr="+i);
    // writer.println("jcp link.Merlin chr="+i);
    // writer.println("merlin -d chr"+chrome+".dat -p re_chrom"+chrome+".pre -m chr"+chrome+".map
    // --npl --tabulate --step 5 --markerNames --information --ibd --prefix merlin-chr"+chrome+"");
    // writer.println();
    // }
    // writer.close();
    // } catch (Exception e) {
    // System.err.println("Error writing to "+"MerlinAll.bat");
    // e.printStackTrace();
    // }
    if (qsub != null) {
      if (blade) {
        commands =
            "merlin -d chr##.dat -p re_chrom##.pre -m chr##.map --npl --tabulate --step 5 --markerNames --information --ibd --prefix merlin-chr##";
        Files.qsub(qsub, "java", CHR_START, CHR_STOP, commands, 1, 10, null);
      } else {
        commands =
            "/share/apps/bin/merlin -d chr##.dat -p re_chrom##.pre -m chr##.map --npl --tabulate --step 5 --markerNames --information --ibd --prefix merlin-chr##";
        Files.qsub(qsub, CHR_START, CHR_STOP, commands);
      }
    } else {
      commands = "echo \"Starting chr# at...\"\n" + "date\n" + "jcp link.Merlin chr=#\n"
          + "merlin -d chr##.dat -p re_chrom##.pre -m chr##.map --npl --tabulate --step 5 --markerNames --information --ibd --prefix merlin-chr##\n"
          + "echo \"Finished chr# at...\"\n" + "date\n";
      Files.batchIt("batch", 0, CHR_START, CHR_STOP, 4, commands);
    }
  }

  public static void batchAllQuant(double[] quant, String qsub, boolean blade) {
    String commands;
    boolean win = Files.isWindows();

    if (qsub != null) {
      if (blade) {
        commands = "cd " + qsub + "\n" + "java -cp /home/bc2/pankratz/"
            + org.genvisis.common.PSF.Java.GENVISIS + " link.Merlin chr=#\n"
            + "merlin-regress -d chr##.dat -p re_chrom##.pre -m chr##.map --mean " + quant[0]
            + " --var " + quant[1] + " --her " + quant[2]
            + " --tabulate --step 5 --prefix regress-chr##";
        Files.qsub(qsub + "_regress", null, CHR_START, CHR_STOP, commands, 1, 48, null);
      } else {
        commands = "cd " + qsub + "\n" + Files.getRunString() + " link.Merlin chr=#\n"
            + "/share/apps/bin/merlin-regress -d chr##.dat -p re_chrom##.pre -m chr##.map --mean "
            + quant[0] + " --var " + quant[1] + " --her " + quant[2]
            + " --tabulate --step 5 --prefix regress-chr##";
        Files.qsub(qsub + "_regress", CHR_START, CHR_STOP, commands);
      }
    } else {
      commands = "echo \"Starting chr# at...\"\n" + (win ? "date /t\ntime /t\n" : "date\n")
          + "java -cp " + (win ? "C:" : "") + "/home/npankrat/"
          + org.genvisis.common.PSF.Java.GENVISIS + " link.Merlin chr=#\n"
          + "merlin-regress -d chr##.dat -p re_chrom##.pre -m chr##.map --mean " + quant[0]
          + " --var " + quant[1] + " --her " + quant[2]
          + " --tabulate --step 5 --information --ibd --prefix regress-chr##\n"
          + "echo \"Finished chr# at...\"\n" + (win ? "date /t\ntime /t\n" : "date\n");
      Files.batchIt("regress", 0, CHR_START, CHR_STOP, 4, commands);
    }
  }

  public static void batchAllVC(String qsub, boolean blade) {
    String commands;
    boolean win = Files.isWindows();

    if (qsub != null) {
      if (blade) {
        commands = "" + "cd " + qsub + "\n" + "java -cp /home/bc2/pankratz/"
            + org.genvisis.common.PSF.Java.GENVISIS + " link.Merlin chr=#\n"
            + "merlin -d chr##.dat -p re_chrom##.pre -m chr##.map --vc --tabulate --step 5 --markerNames --information --prefix vc-chr## > vc-chr##.log";
        Files.qsub(qsub + "_vc", "java", CHR_START, CHR_STOP, commands, 1, 10, null);
      } else {
        commands = "" + "cd " + qsub + "\n" + Files.getRunString() + " link.Merlin chr=#\n"
            + "/share/apps/bin/merlin -d chr##.dat -p re_chrom##.pre -m chr##.map --vc --tabulate --step 5 --markerNames --information --prefix vc-chr## > vc-chr##.log";
        Files.qsub(qsub + "_vc", CHR_START, CHR_STOP, commands);
      }
    } else {
      commands = "echo \"Starting chr# at...\"\n" + (win ? "date /t\ntime /t\n" : "date\n")
          + "java -cp " + (win ? "C:" : "") + "/home/npankrat/"
          + org.genvisis.common.PSF.Java.GENVISIS + " link.Merlin chr=#\n"
          + "merlin -d chr##.dat -p re_chrom##.pre -m chr##.map --vc --tabulate --step 5 --markerNames --information --prefix vc-chr## > vc-chr##.log\n"
          + "echo \"Finished chr# at...\"\n" + (win ? "date /t\ntime /t\n" : "date\n");
      Files.batchIt("vc", 0, 8, CHR_STOP, 5, commands);
    }
  }

  public static boolean checkTrait(int chr) {
    BufferedReader reader;
    boolean binary;
    String[] line;
    int conf;

    binary = true;
    conf = 0;
    try {
      reader = new BufferedReader(new FileReader("re_chrom" + ext.chrome(chr) + ".pre"));
      while (binary && conf < 100 && reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        if (line[5].equals("0")) {

        } else if (line[5].equals("1") || line[5].equals("2")) {
          conf++;
        } else {
          binary = false;
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + "re_chrom" + ext.chrome(chr) + ".pre"
          + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + "re_chrom" + ext.chrome(chr) + ".pre" + "\"");
      System.exit(2);
    }

    return binary;
  }

  public static void createMerlinFiles(int chr, int trait) {
    createMerlinFiles(chr, "chr" + ext.chrome(chr), trait, null);
  }

  public static void createMerlinFiles(int chr, String root, int trait, String[] covariateNames) {
    createMerlinFiles(null, chr, root, trait, covariateNames);
  }

  public static void createMerlinFiles(String dir, int chr, String root, int trait,
      String[] covariateNames) {
    PrintWriter writer;
    LinkageMap map;
    String[] markerNames;
    double[] positions;
    double[][] alleleFreqs;

    map = dir == null ? new LinkageMap(chr) : new LinkageMap(dir, chr);
    markerNames = map.getMarkerNames();
    positions = map.getCumulativePositions(false);
    alleleFreqs = map.getAlleleFreqs();

    try {
      writer = new PrintWriter(new FileWriter(root + ".dat"));
      switch (trait) {
        case UNKNOWN_TRAIT:
          writer.println(checkTrait(chr) ? "A affection_status" : "T trait");
          break;
        case BINARY_TRAIT:
          writer.println("A affection_status");
          break;
        case QUANTITATIVE_TRAIT:
          writer.println("T trait");
          break;
      }
      for (int i = 0; covariateNames != null && i < covariateNames.length; i++) {
        writer.println("C " + covariateNames[i]);
      }
      for (String markerName : markerNames) {
        writer.println("M " + markerName);
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing " + root + ".dat");
      e.printStackTrace();
    }

    try {
      writer = new PrintWriter(new FileWriter(root + ".map"));
      writer.println("CHROMOSOME\tMARKER\tPOSITION");
      for (int i = 0; i < markerNames.length; i++) {
        writer.println(chr + "\t" + markerNames[i] + "\t" + ext.formDeci(positions[i], 10, false));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing " + root + ".map");
      e.printStackTrace();
    }

    try {
      writer = new PrintWriter(new FileWriter(root + ".freq"));
      for (int i = 0; i < markerNames.length; i++) {
        writer.println("M " + markerNames[i]);
        writer.println("F " + Array.toStr(alleleFreqs[i], 6, 6, " "));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing " + root + ".freq");
      e.printStackTrace();
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    int chr = -1;
    boolean batch = false;
    String qsub = null;
    boolean vc = false;
    boolean prep = false;
    int trait = UNKNOWN_TRAIT;
    double[] quant = null;
    boolean blade = false;

    // batchAllBinary("HL");
    // System.exit(1);

    // batchAllVC("pta", blade);
    // batchAllVC("Prin1", blade);
    // batchAllVC("Prin2", blade);
    // batchAllVC("Prin3", blade);
    // System.exit(1);

    // vc = true;
    // qsub = "PCA1";
    // chr=1;

    // PCA1
    // qsub = "PCA1";
    // quant = new double[] {0.070835749, 3.30957236, 0.4079};

    // PCA1 boxcox
    // qsub = "PCA1_boxcox";
    // quant = new double[] {3.209175862, 0.541843525, 0.4079};

    // PCA2
    // quant = new double[] {0.010311256, 1.196156432, 0.2450};

    // PCA3
    // quant = new double[] {0.016727564, 0.587978393, 0.4974};

    String usage = "\n" + "link.Merlin requires 0-1 arguments\n"
        + "   (1) generate files for specific chromosome (i.e. chr=2 (not the default))\n"
        + "   (2) binary trait (i.e. -binary)\n" + "   (3) quantitative trait (i.e. -quant)\n"
        + "  OR\n" + "   (1) generate files for all chromsomes (i.e. -prep (not the default))\n"
        + "  OR\n"
        + "   (1) batch create and analyze a binary trait (i.e. -batch (not the default))\n"
        + "  OR\n"
        + "   (1) batch create and analyze a quantitative trait using vc (i.e. -vc (not the default))\n"
        + "   (2) make .qsub files for this directory instead of batch files (i.e. qsub=dir (not the default; don't use slash))\n"
        + "   (3) use blade header for qsub instead of alc (i.e. -blade (not the default))\n"
        + "  OR\n"
        + "   (1) batch create and analyze a quantitative trait using h-e (i.e. quant=mean,variance,heritability)\n"
        + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("chr=")) {
        chr = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("-batch")) {
        batch = true;
        trait = BINARY_TRAIT;
        numArgs--;
      } else if (arg.startsWith("-prep")) {
        prep = true;
        numArgs--;
      } else if (arg.startsWith("qsub=")) {
        qsub = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-binary")) {
        trait = BINARY_TRAIT;
        numArgs--;
      } else if (arg.startsWith("-vc")) {
        vc = true;
        numArgs--;
      } else if (arg.startsWith("-quant")) {
        trait = QUANTITATIVE_TRAIT;
        numArgs--;
      } else if (arg.startsWith("quant=")) {
        quant = Array.toDoubleArray(arg.split("=")[1].split(","));
        if (quant.length == 3) {
          numArgs--;
        } else {
          System.err.println("Error - batchAllQuant requires three values");
        }
      } else if (arg.startsWith("-blade")) {
        blade = true;
        numArgs--;
      } else {
        System.err.println("Error - don't know what to do with argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    try {
      if (prep) {
        prepAll();
      } else if (batch) {
        batchAllBinary(qsub, blade);
      } else if (vc) {
        batchAllVC(qsub, blade);
      } else if (quant != null) {
        batchAllQuant(quant, qsub, blade);
      } else if (chr > 0) {
        createMerlinFiles(chr, trait);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void prepAll() {
    for (int i = CHR_START; i <= CHR_STOP; i++) {
      createMerlinFiles(i, UNKNOWN_TRAIT);
    }
  }
}
