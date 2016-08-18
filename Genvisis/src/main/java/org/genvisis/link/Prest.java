// commandline = prest
package org.genvisis.link;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
// import java.util.*;
import java.text.DecimalFormat;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.ext;

public class Prest {
  public static final String[] RELATIONSHIP_TRANSLATION =
      {".", "FS", "HS", "GG", "AV", "CO", "UN", "HA", "2C", "2C", "PO"};

  public static void createFiles() {
    BufferedReader reader;
    PrintWriter writer, chromfiles;
    String[] line;
    String temp, chrome;
    DecimalFormat myFormatter = new DecimalFormat("0.######");
    int chr, count;
    double max;

    chr = 1;
    while (!new File("re_chrom" + ext.chrome(chr) + ".pre").exists() && chr++ < 22) {
      ;
    }
    if (chr == 23) {
      System.err.println("There were no re_chrom##.pre files that could be used.");
      System.exit(1);
    }

    chrome = ext.chrome(chr);
    try {
      reader = new BufferedReader(new FileReader("re_chrom" + chrome + ".pre"));
      writer = new PrintWriter(new FileWriter("pedigrees"));
      while (reader.ready()) {
        writer.println(Array.toStr(Array.subArray(reader.readLine().trim().split("[\\s]+"), 0, 6)));
      }
      reader.close();
      writer.close();
    } catch (IOException ioe) {
      System.err.println("Error parsing \"" + "re_chrom" + chrome + ".pre" + "\"");
      ioe.printStackTrace();
    }

    try {
      chromfiles = new PrintWriter(new FileWriter("chromfiles"));
      for (int chromosome = 1; chromosome <= 22; chromosome++) {
        chrome = ext.chrome(chromosome);

        try {
          reader = new BufferedReader(new FileReader("map" + chrome + ".dat"));
          writer = new PrintWriter(new FileWriter("map" + chrome + ".idx"));

          writer.println((Integer.parseInt(reader.readLine().trim().split("[\\s]+")[0]) - 1) + "");
          do {
            temp = reader.readLine();
          } while (!temp.startsWith("3 "));
          count = 0;
          do {
            writer.println(temp);
            writer.println(reader.readLine());
            count++;
            temp = reader.readLine();
          } while (temp.startsWith("3 "));
          line = reader.readLine().trim().split("[\\s]+");

          max = 0;
          for (int i = 0; i < count; i++) {
            max = Math.max(Double.parseDouble(line[i]), max);
          }

          for (int i = 1; i < count; i++) {
            if (max > 1) {
              writer.print(myFormatter.format(Double.parseDouble(line[i]) / 100) + " ");
            } else {
              writer.print(line[i] + " ");
            }
          }
          writer.println();

          reader.close();
          writer.close();
          chromfiles.println("map" + chrome + ".idx re_chrom" + chrome + ".pre");
        } catch (Exception e) {
          System.err.println("skipped chromosome " + chromosome);
        }
      }
      chromfiles.close();
    } catch (Exception e) {
      System.err.println("Error parsing chromfiles");
      e.printStackTrace();
    }

    System.out.println("Try: prest pedigrees chromfiles 1");
  }

  public static void slimGenome(String dir) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String genome = "plink.genome";

    try {
      reader = new BufferedReader(new FileReader(dir + genome));
      writer = new PrintWriter(new FileWriter(dir + genome + ".slim"));
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        writer.println(line[0] + "-" + line[1] + ":" + line[2] + "-" + line[3] + "\t" + line[5]
                       + "\t" + line[6] + "\t" + line[7] + "\t" + line[8] + "\t" + line[9]);
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + genome + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + genome + "\"");
      System.exit(2);
    }
  }

  public static void slimPrest(String dir) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String prest = "prest_out2";
    double min;

    try {
      reader = new BufferedReader(new FileReader(dir + prest));
      writer = new PrintWriter(new FileWriter(dir + prest + ".slim"));
      writer.println("FID1-IID1:FID2-IID2\tPutative#\tPutativeRel\t#markers\tIBD\tp0\tp1\tp2\tEIBD-pval\tAIBS-pval\tIBS-pval\tmin-pval");
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        min = 999;
        for (int i = 9; i <= 11; i++) {
          if (!line[i].equals("NA")) {
            min = Math.min(Double.parseDouble(line[i]), min);
          }
        }
        writer.println(line[0] + "-" + line[1] + ":" + line[0] + "-" + line[2] + "\t" + line[3]
                       + "\t" + RELATIONSHIP_TRANSLATION[Integer.parseInt(line[3])] + "\t" + line[4]
                       + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + line[8] + "\t"
                       + line[9] + "\t" + line[10] + "\t" + line[11] + "\t" + min);
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + prest + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + prest + "\"");
      System.exit(2);
    }
  }

  public static void mergePrestWithPlink(String dir) {
    try {
      Files.merge(dir + "plink.genome.slim", 0, new int[] {0, 1, 2, 3, 4, 5}, false,
                  dir + "prest_out2.slim", 0, new int[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, true,
                  dir + "PlinkVersusPrest.xln");
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    // String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\Family structure
    // issues\\secondPass\\prest\\";
    String dir =
        "C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\Family structure issues\\testFinalComplete\\";
    boolean create = false;
    boolean slimGenome = true;
    boolean slimPrest = true;
    boolean merge = true;

    String usage = "\n" + "link.bat.Prest requires 0-1 arguments\n"
                   + "   (1) create files for Prest (i.e. -create (not the default))\n"
                   + "   (2) slim genome (i.e. -slimGenome (not the default))\n"
                   + "   (3) slim prest2 file (i.e. -slimPrest (not the default))\n"
                   + "   (4) merge Prest results with slimmed plink results (i.e. -merge (not the default))\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("-create")) {
        create = true;
        numArgs--;
      } else if (arg.startsWith("-slimGenome")) {
        slimGenome = true;
        numArgs--;
      } else if (arg.startsWith("-slimPrest")) {
        slimPrest = true;
        numArgs--;
      } else if (arg.startsWith("-merge")) {
        merge = true;
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    if (!new File(dir).exists()) {
      dir = "";
      System.err.println("Error - directory '" + dir + "' not found, trying current directory...");
    }

    try {
      if (create) {
        createFiles();
      }
      if (slimGenome) {
        slimGenome(dir);
      }
      if (slimPrest) {
        slimPrest(dir);
      }
      if (merge) {
        mergePrestWithPlink(dir);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
