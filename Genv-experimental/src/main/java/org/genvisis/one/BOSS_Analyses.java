package org.genvisis.one;

// import java.util.*;
// import common.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.Hits;
import org.genvisis.gwas.CreateDatabaseFromPlink;

public class BOSS_Analyses {
  public static void batch(String dir) {

  }

  public static void generate(String dir) {
    new File(dir + "split/").mkdirs();
    for (int chr = 1; chr <= 23; chr++) {
      // CmdLine.run("plink --bfile ../plink --chr "+chr+" --recode --maf 0.001 --geno 0.2 --out
      // chr"+chr, dir+"split/");
      CmdLine.run("plink --file chr" + chr + " --freq --out chr" + chr, dir + "split/");
    }

    new File(dir + "gwaf/").mkdirs();
    for (int chr = 1; chr <= 23; chr++) {
      CreateDatabaseFromPlink.toGWAF(dir + "split/chr" + chr + ".ped",
          dir + "split/chr" + chr + ".map", dir + "split/chr" + chr + ".frq", null,
          dir + "gwaf/chr" + chr + ".fhsR");
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String dir, resultDir, mafDir, geneListFile;
    Boolean createGwaf, isGenerate, isStandardize, isParse, isProcessMacs;

    dir = "D:/temp/BOSS/E1/";
    resultDir = "/home/pankrat2/shared/boss/emmax/prpidsk/results/ver5_newPCs_residual/";
    mafDir = "/home/pankrat2/shared/boss/emmax/prpidsk/phenocovars/ver5_newPCs_residual/";
    geneListFile = "/home/pankrat2/shared/boss/emmax/analyzeThis.b36_genes.xln";
    createGwaf = false;
    isGenerate = false;
    isStandardize = false;
    isParse = true;
    isProcessMacs = false;
    String usage = "\n" + "one.BOSS_Analyses requires 0-1 arguments\n"
        + "   (1) directory (i.e. dir=" + dir + " (default))\n"
        + "   (2) to parse result (i.e. parseresult=" + isParse + " (default))\n"
        + "   (3) result directory (i.e. resultdir=" + resultDir + " (default))\n"
        + "   (4) minor allele frequency file directory (i.e. mafdir=" + mafDir + " (default))\n"
        + "   (5) gene list file full path (i.e. genelist=" + geneListFile + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("resultdir=")) {
        resultDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("mafdir=")) {
        mafDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("genelist=")) {
        geneListFile = arg.split("=")[1];
        numArgs--;
      } else {
        System.out.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }


    if (createGwaf) {
      CreateDatabaseFromPlink.toGWAF(dir + "ICAM1/plink.ped", dir + "ICAM1/plink.map",
          dir + "ICAM1/plink.frq", null, dir + "ICAM1/gwaf.fhsR");
    } else if (isGenerate) {
      generate(dir);
    } else if (isStandardize) {
      standardize(dir);
    } else if (isParse) {
      parseResults(resultDir, mafDir, geneListFile, null);
    } else if (isProcessMacs) {
      processMACs(dir + "mafCheck.frq", dir + "minor_allele_counts.xln");
    } else {
      System.out.println("Please select a command to run.");
    }
  }

  public static void parseResults(String resultDir, String mafDir, String geneListFileFullPath,
      Logger log) {
    String[] filenames = null, array_of_hits, params;
    String root;
    Hits hits;

    // filenames = Files.list(dir, ".ps", false);
    filenames = new File(resultDir).list(new FilenameFilter() {
      @Override
      public boolean accept(File file, String filename) {
        if (!filename.toLowerCase().endsWith(".ps")) {
          return false;
        }
        if (filename.toLowerCase().contains("_hibs")) {
          return false;
        }
        if (new File(file, filename).isDirectory()) {
          return false;
        }
        return true;
      }
    });

    if (log == null) {
      log = new Logger();
    }
    if (filenames == null) {
      log.report("Did not find any .ps file in the directory " + resultDir);
    } else {
      hits = new Hits();
      for (byte i = 0; i < filenames.length; i++) {
        hits.incorporateFromFile(resultDir + filenames[i], new int[] {0, 2}, 0.001, log);
      }

      hits.writeHits(resultDir + "top_snps.txt");
      array_of_hits =
          HashVec.loadFileToStringArray(resultDir + "top_snps.txt", false, new int[] {0}, false);

      params = new String[filenames.length * 2 + 2];
      // params[0] = "gene_position.xln 0 1 2 3 4 5 6";
      params[0] = geneListFileFullPath + " 0 1 2 3 4 5 6 tab";
      params[1] = resultDir + "top_snps.txt 0 1=minPval skip=0";

      for (byte i = 0; i < filenames.length; i++) {
        root = ext.rootOf(filenames[i]);
        params[2 + i * 2 + 0] =
            resultDir + filenames[i] + " 0 1=beta_" + root + " 2=pval_" + root + " skip=0";
        params[2 + i * 2 + 1] = resultDir + root + ".mac 0 1=MAC_" + root;

        processMACs(mafDir + root + "/mafCheck.frq", resultDir + root + ".mac");
        // processMACs(dir + "../phenocovars/" + root + "/mafCheck.frq", dir + root + ".mac");
        // processMACs(dir + root + "/mafCheck.frq", dir + root + ".mac");
      }

      Files.combine(array_of_hits, params, null, "MarkerName", ".", resultDir + "top_hits.xln", log,
          true, true, false);
    }
  }

  private static void processMACs(String plinkFrqFile, String outputFile) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;

    try {
      writer = new PrintWriter(new FileOutputStream(outputFile));
      writer.println("SNP\tMAC");
      reader = new BufferedReader(new FileReader(plinkFrqFile));
      reader.readLine();
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        writer.println(line[1] + "\t" + (line[4].equalsIgnoreCase("NA") ? "."
            : (int) (Double.parseDouble(line[4]) * Integer.parseInt(line[5]))));
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  public static void standardize(String dir) {
    BufferedReader reader, mapReader, infoReader;
    PrintWriter writer, mapWriter, infoWriter;
    String header, root;
    int count, fileNum;

    root = "gwaf/chr";

    count = 0;
    fileNum = 0;
    new File(dir + "std/").mkdirs();
    try {
      writer = new PrintWriter(new FileWriter(dir + "std/file" + fileNum + ".gen"));
      mapWriter = new PrintWriter(new FileWriter(dir + "std/file" + fileNum + ".pmap"));
      infoWriter = new PrintWriter(new FileWriter(dir + "std/file" + fileNum + ".mlinfo"));
      for (int chr = 1; chr <= 22; chr++) {
        try {
          reader = new BufferedReader(new FileReader(dir + root + chr + ".gen"));
          mapReader = new BufferedReader(new FileReader(dir + root + chr + ".pmap"));
          infoReader = new BufferedReader(new FileReader(dir + root + chr + ".mlinfo"));
          header = infoReader.readLine();
          while (reader.ready()) {
            if (count == 0) {
              infoWriter.println(header);
            }
            writer.println(reader.readLine());
            mapWriter.println(mapReader.readLine());
            infoWriter.println(infoReader.readLine());
            count++;
            if (count == 10000) {
              writer.close();
              mapWriter.close();
              infoWriter.close();
              fileNum++;
              writer = new PrintWriter(new FileWriter(dir + "std/file" + fileNum + ".gen"));
              mapWriter = new PrintWriter(new FileWriter(dir + "std/file" + fileNum + ".pmap"));
              infoWriter = new PrintWriter(new FileWriter(dir + "std/file" + fileNum + ".mlinfo"));
              count = 0;
            }
          }
          reader.close();
          mapReader.close();
          infoReader.close();
        } catch (FileNotFoundException fnfe) {
          System.err.println("Error: file \"" + dir + "leslie_lange.FHS.IBC.CEU.chr1.gen"
              + "\" not found in current directory");
          System.exit(1);
        } catch (IOException ioe) {
          System.err
              .println("Error reading file \"" + dir + "leslie_lange.FHS.IBC.CEU.chr1.gen" + "\"");
          System.exit(2);
        }
      }
      writer.close();
      mapWriter.close();
    } catch (Exception e) {
      System.err.println("Error writing to file #" + fileNum);
      e.printStackTrace();
    }
  }

}
