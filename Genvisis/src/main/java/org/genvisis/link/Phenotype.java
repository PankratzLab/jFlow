package org.genvisis.link;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

import com.google.common.primitives.Ints;

public class Phenotype {
  public static final String[] OTHER_FILES =
      {"map##.dat", "chr##.dat", "chr##.map", "chr##.freq", "run#.qsub", "chr#_vc.qsub"};

  public static void main(String[] args) {
    int numArgs = args.length;
    String dir = "";
    String pheno = "pheno.dat";
    String pattern = "re_chrom##.pre";
    int index = 5;
    String missingValue = "x";
    String logfile = null;
    Logger log;

    String usage = "\n" + "link.Phenotype requires 0-1 arguments\n" + "   (1) directory (i.e. dir="
        + dir + " (default))\n" + "   (2) phenotype filename (i.e. pheno=" + pheno + " (default))\n"
        + "   (3) pattern of files to update (i.e. pattern=" + pattern + " (default))\n"
        + "   (4) index of column to swap out (i.e. index=" + index + " (default))\n"
        + "   (5) missing value (i.e. missingValue=" + missingValue + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("pheno=")) {
        pheno = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("pattern=")) {
        pattern = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("index=")) {
        index = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("missingValue=")) {
        missingValue = ext.parseStringArg(arg, "x");
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    dir = "D:/BOSS/Linkage/PCA_hits/";
    dir = "D:/BOSS/Linkage/PCA_all_files/";

    try {
      log = new Logger(logfile);
      update(dir, pheno, pattern, index, missingValue, log);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  private static void update(String dir, String pheno, String pattern, int index,
      String missingValue, Logger log) {
    BufferedReader reader;
    PrintWriter[] writers;
    String[] line, phenoNames, phenos;
    String trav;
    Hashtable<String, String> hash;
    String filename;
    IntVector skips;

    dir = ext.verifyDirFormat(dir);
    phenoNames = Files.getHeaderOfFile(dir + pheno, log);
    ext.checkHeader(phenoNames, new String[] {"FID", "IID"}, new int[] {0, 1}, false, log, false);
    System.out.println("Loading " + pheno);
    hash = HashVec.loadFileToHashString(dir + pheno, new int[] {0, 1},
        Array.subArray(Array.arrayOfIndices(phenoNames.length), 2), pheno.endsWith(".csv"), "\t",
        true, false, false);
    phenoNames = Array.subArray(phenoNames, 2);
    writers = new PrintWriter[phenoNames.length];
    for (int i = 0; i < writers.length; i++) {
      new File(dir + phenoNames[i] + "/").mkdirs();
    }

    skips = new IntVector();
    for (int chr = 1; chr <= 23; chr++) {
      filename = ext.insertNumbers(pattern, chr);
      if (Files.exists(dir + filename)) {
        try {
          reader = new BufferedReader(new FileReader(dir + filename));
          for (int i = 0; i < writers.length; i++) {
            writers[i] = new PrintWriter(new FileWriter(dir + phenoNames[i] + "/" + filename));
          }
          while (reader.ready()) {
            line = reader.readLine().trim().split("[\\s]+");
            trav = hash.get(line[0] + "\t" + line[1]);
            if (trav == null) {
              phenos = Array.stringArray(phenoNames.length, missingValue);
            } else {
              phenos = trav.split("[\\s]+");
            }
            for (int i = 0; i < writers.length; i++) {
              line[index] = phenos[i];
              writers[i].println(Array.toStr(line));
            }

          }
          reader.close();
          for (int i = 0; i < writers.length; i++) {
            writers[i].close();
            for (String element : OTHER_FILES) {
              filename = ext.insertNumbers(element, chr);
              if (Files.exists(dir + filename)) {
                Files.copyFile(dir + filename, dir + phenoNames[i] + "/" + filename);
              }
            }
          }
        } catch (FileNotFoundException fnfe) {
          System.err
              .println("Error: file \"" + dir + filename + "\" not found in current directory");
          System.exit(1);
        } catch (IOException ioe) {
          System.err.println("Error reading file \"" + dir + filename + "\"");
          System.exit(2);
        }
        log.report("Updating chr" + chr);

      } else {
        skips.add(chr);
      }
    }

    if (skips.size() > 0) {
      log.report("skipped chromosomes: " + ext.listRanges(Ints.toArray(skips)));
    }



  }
}
