// -Xms1024M -Xmx1024M
package org.genvisis.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;
import java.util.zip.GZIPInputStream;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;
import org.genvisis.filesys.FamilyStructure;
import org.genvisis.filesys.SnpMarkerSet;

import com.google.common.primitives.Ints;

public class Mach {
  // public static final String EXECUTABLE = "mach1";
  public static final String EXECUTABLE = "/share/apps/bin/mach1";
  public static final double P_THRESHOLD = 0.0001;
  public static final int BETA_SIG = 5;
  public static final int SIG_SIG = 6;
  public static final String DEFAULT_SUBSET = "subset";
  public static final String DEFAULT_ALL = "all";
  // public static final int IMPUTATION_AWARE_COLUMN = 5; // Quality (always higher)
  public static final int IMPUTATION_AWARE_COLUMN = 6; // Rsq
  // public static final int PICK_BEST_SIZE = 50;
  public static final int PICK_BEST_SIZE = 125;
  public static final String[] MLINFO_HEADER = {"SNP", "Al1", "Al2", "Freq1", "MAF", "Quality",
                                                "Rsq"};
  public static final String[] MINFO_HEADER = {"SNP", "Al1", "Al2", "Freq1", "MAF", "AvgCall",
                                               "Rsq", "Genotyped", "LooRsq", "EmpR", "EmpRsq",
                                               "Dose1", "Dose2"};
  public static final String[] MACH2DAT_HEADER = {"TRAIT", "MARKER", "ALLELES", "FREQ1", "RSQR",
                                                  "EFFECT1", "OR", "STDERR", "WALDCHISQ", "PVALUE",
                                                  "LRCHISQ", "LRPVAL", "NCASES", "NCONTROLS"};
  public static final char[] DELIMITERS = {' ', (char) 9};

  public static void batchFileCreation(String root, String keeps, String excludes,
                                       String prefixSubset, String prefixAll) {
    PrintWriter writer;

    try {
      writer = new PrintWriter(new FileWriter(root + ".batch"));
      if (keeps == null) {
        System.err.println("Warning - subset and all will be the same");
      }
      for (int i = 1; i <= 22; i++) {
        writer.println("mkdir chr" + i);
        writer.println("cd chr" + i);
        writer.println("plink --noweb --bfile ../" + root
                       + (excludes == null ? "" : " --exclude ../" + excludes) + " --chr " + i
                       + " --out " + prefixAll + ".chr" + i + " --recode");
        writer.println("plink --noweb --bfile ../" + root
                       + (excludes == null ? "" : " --exclude ../" + excludes) + " --chr " + i
                       + " --out " + prefixSubset + ".chr" + i + " --recode"
                       + (keeps == null ? "" : " --keep ../" + keeps));
        writer.println(Files.getRunString() + " gwas.Mach convert=" + prefixSubset + ".chr" + i);
        writer.println(Files.getRunString() + " gwas.Mach convert=" + prefixAll + ".chr" + i);
        writer.println("cd ..");
        writer.println();
      }

      writer.close();
      Files.chmod(root + ".batch");
    } catch (Exception e) {
      System.err.println("Error writing batchSplit");
    }
  }

  public static void convertMapToDat(String prefix) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;

    try {
      reader = new BufferedReader(new FileReader(prefix + ".map"));
      writer = new PrintWriter(new FileWriter(prefix + ".dat"));
      writer.println("A\tPD_Affection");
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        writer.println("M\t" + line[1]);
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + prefix + ".map" + "\" not found in current directory");
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + prefix + ".map" + "\"");
    }
  }

  public static void steps(String prefixSubset, String prefixAll, String[] nodesToUse) {
    String commands;
    Vector<String> v;
    // int step;
    String[] list;

    commands = "cd chr#\n";
    // commands += EXECUTABLE+" -d "+prefixSubset+".chr#.dat -p "+prefixSubset+".chr#.ped --autoFlip
    // --hapmapFormat --snps
    // "+(trimming?"truncated":"genotypes")+"_chr#_CEU_r22_nr.b36_fwd_legend.txt --haps
    // "+(trimming?"truncated":"genotypes")+"_chr#_CEU_r22_nr.b36_fwd.phase --rounds 100 --greedy
    // --prefix MACH_step1_chr#\n";
    // commands += EXECUTABLE+" -d "+prefixAll+".chr#.dat -p "+prefixAll+".chr#.ped --autoFlip
    // --hapmapFormat --snps
    // "+(trimming?"truncated":"genotypes")+"_chr#_CEU_r22_nr.b36_fwd_legend.txt --haps
    // "+(trimming?"truncated":"genotypes")+"_chr#_CEU_r22_nr.b36_fwd.phase --crossover
    // MACH_step1_chr#.rec --errormap MACH_step1_chr#.erate --greedy --mle --mldetails --prefix
    // MACH_step2_chr#\n";
    commands += EXECUTABLE + " -d " + prefixSubset + ".chr#.dat -p " + prefixSubset
                + ".chr#.ped --rounds 30 --states 400 --phase --interim 5 --sample 5 --prefix MACH_step1_chr#\n";
    commands += EXECUTABLE + " -d " + prefixAll + ".chr#.dat -p " + prefixAll
                + ".chr#.ped --crossover MACH_step1_chr#.rec --errormap MACH_step1_chr#.erate --phase --prefix MACH_step2_chr#\n";
    commands += "cd .." + "\n";
    commands += "\n";

    v = new Vector<String>();
    if (nodesToUse == null) {
      v = Array.toStringVector(Files.qsub("", "chr#", 1, 22, commands, null, 10000, 48, null));
    } else {
      v = new Vector<String>();
      // step = (int)Math.ceil((double)(22)/(double)nodesToUse.length);
      // for (int i = 0; i < nodesToUse.length; i++) {
      // list = Files.qsub("", null, i*step+0, i==nodesToUse.length-1?22:((i+1)*step-1), commands,
      // "chr", null, -1, nodesToUse[i]);
      for (int chr = 1; chr <= 22; chr++) {
        list = Files.qsub("", "chr#", chr, chr, commands, null, 10000, 48,
                          nodesToUse[chr % nodesToUse.length]);
        for (String element : list) {
          v.add(element);
        }
      }
    }
    Files.writeList(Array.toStringArray(v), "master.mach");
    Files.chmod("master.mach");
  }

  public static void trimReference(String prefixAll, int chr) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String trav;
    int count;
    int pos, start, stop, startAtIndex, stopBeforeIndex;

    start = Integer.MAX_VALUE;
    stop = -1;
    try {
      reader = new BufferedReader(new FileReader(prefixAll + ".chr" + chr + ".map"));
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        pos = Integer.parseInt(line[3]);
        start = Math.min(start, pos);
        stop = Math.max(stop, pos);
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + prefixAll + ".chr" + chr + ".map"
                         + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + prefixAll + ".chr" + chr + ".map" + "\"");
      System.exit(2);
    }

    startAtIndex = Integer.MAX_VALUE;
    stopBeforeIndex = -1;
    count = 0;
    try {
      reader = new BufferedReader(new FileReader("genotypes_chr" + chr
                                                 + "_CEU_r22_nr.b36_fwd_legend.txt"));
      writer = new PrintWriter(new FileWriter("truncated_chr" + chr
                                              + "_CEU_r22_nr.b36_fwd_legend.txt"));
      writer.println(reader.readLine());
      while (reader.ready()) {
        trav = reader.readLine();
        line = trav.trim().split("[\\s]+");
        pos = Integer.parseInt(line[1]);
        if (pos >= start && pos <= stop) {
          writer.println(trav);
          if (startAtIndex == Integer.MAX_VALUE) {
            startAtIndex = count;
          }
        } else if (startAtIndex != Integer.MAX_VALUE && stopBeforeIndex == -1) {
          stopBeforeIndex = count;
        }
        count++;
      }
      writer.close();
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + "genotypes_chr" + chr + "_CEU_r22_nr.b36_fwd_legend.txt"
                         + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + "genotypes_chr" + chr
                         + "_CEU_r22_nr.b36_fwd_legend.txt" + "\"");
      System.exit(2);
    }

    System.out.println("Truncating HapMap files to chr" + chr + ":" + start + "-" + stop
                       + " (indices: " + startAtIndex + " to " + stopBeforeIndex + ")");

    try {
      FileReader in;
      FileWriter out;
      int c;
      boolean done = false;

      in = new FileReader("genotypes_chr" + chr + "_CEU_r22_nr.b36_fwd.phase");
      out = new FileWriter("truncated_chr" + chr + "_CEU_r22_nr.b36_fwd.phase");
      while (!done) {
        for (int i = 0; i < startAtIndex && !done; i++) {
          c = in.read();
          if (c == -1) {
            done = true;
          }
          in.read();
        }
        for (int i = startAtIndex; i < stopBeforeIndex && !done; i++) {
          c = in.read();
          if (c == -1) {
            done = true;
          } else {
            out.write(c);
            out.write(in.read());
          }
        }
        for (int i = stopBeforeIndex; i < count; i++) {
          in.read();
          in.read();
        }
        c = in.read();
        if (c == 32) {
          c = in.read();
        }
        if (c == -1) {
          done = true;
        } else {
          out.write(c);
        }
      }

      in.close();
      out.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + "genotypes_chr" + chr + "_CEU_r22_nr.b36_fwd.phase"
                         + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + "genotypes_chr" + chr
                         + "_CEU_r22_nr.b36_fwd.phase" + "\"");
      System.exit(2);
    }
  }

  public static void decodePhasedHapMap(String dir, int chr) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, snpData, markerNames, positions;
    String[][] alleles;

    snpData = HashVec.loadFileToStringArray(dir + "truncated_chr" + chr
                                            + "_CEU_r22_nr.b36_fwd_legend.txt", false, true,
                                            new int[] {0, 1, 2, 3}, false);
    markerNames = new String[snpData.length];
    positions = new String[snpData.length];
    alleles = new String[snpData.length][2];
    for (int i = 0; i < snpData.length; i++) {
      line = snpData[i].trim().split("[\\s]+");
      markerNames[i] = line[0];
      positions[i] = line[1];
      alleles[i][0] = line[2];
      alleles[i][1] = line[3];
    }
    try {
      reader = new BufferedReader(new FileReader(dir + "truncated_chr" + chr
                                                 + "_CEU_r22_nr.b36_fwd.phase"));
      writer = new PrintWriter(new FileWriter(dir + "chr" + chr + "_phased.xln"));
      writer.println(Array.toStr(markerNames));
      writer.println(Array.toStr(positions));
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        for (int i = 0; i < markerNames.length; i++) {
          writer.print((i == 0 ? "" : "\t") + alleles[i][Integer.parseInt(line[i])]);
        }
        writer.println();
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + "truncated_chr" + chr
                         + "_CEU_r22_nr.b36_fwd.phase" + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + "truncated_chr" + chr
                         + "_CEU_r22_nr.b36_fwd.phase" + "\"");
      System.exit(2);
    }
  }

  public static void listIndividualsInMldoseSlowly(String filename, String outfile) {
    PrintWriter writer;
    FileReader in;
    String trav;
    long time;
    int c;

    try {
      time = new Date().getTime();
      in = new FileReader(filename);
      writer = new PrintWriter(new FileWriter(outfile));
      System.out.println(ext.getTime() + " Reading dosage file");
      c = 0;
      while (c >= 0) {
        trav = "";
        while (ext.indexOfChar((char) (c = in.read()), DELIMITERS) == -1 && c >= 0) {
          trav += (char) c;
        }
        if (c >= 0) {
          writer.println(trav.substring(0, trav.indexOf("->")) + "\t"
                         + trav.substring(trav.indexOf("->") + 2));
          writer.flush();
        }

        while ((char) (c = in.read()) != '\n' && c >= 0) {
          ;
        }
      }
      System.out.println(ext.getTime() + " Done in " + ext.getTimeElapsed(time));
      writer.close();
      in.close();
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }
  }

  public static void listIndividualsInMldose(String filename, String outfile) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    long time;

    try {
      time = new Date().getTime();
      // reader = new BufferedReader(new FileReader(filename));
      reader = Files.getAppropriateReader(filename);
      System.out.println(ext.getTime() + " Reading dosage file");
      writer = new PrintWriter(new FileWriter(outfile));
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        writer.println(line[0].substring(0, line[0].indexOf("->")) + "\t"
                       + line[0].substring(line[0].indexOf("->") + 2));
        writer.flush();
      }
      writer.close();
      reader.close();
      System.out.println(ext.getTime() + " Done in " + ext.getTimeElapsed(time));
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }
  }

  public static void extractDosageForPlink(String dir, String pedfile, String machDosage,
                                           String machMarkers, String markersToTest) {
    BufferedReader reader;
    PrintWriter writer, mapWriter;
    String[] line;
    Hashtable<String, String> hash = new Hashtable<String, String>();
    Vector<String> v;
    String[] indIDs;
    int count;
    String[] markerNames;
    boolean allMarkers;
    String[] subset = null;

    if (!dir.equals("") && !new File(dir).exists()) {
      System.err.println("Could not find directory: " + dir);
      System.err.println("  using current directory instead");
      dir = "";
    }

    markerNames = Array.toStringArray(HashVec.loadFileToVec(dir + machMarkers, true,
                                                            new int[] {0, 1, 2}, false, false));

    count = 0;
    try {
      reader = new BufferedReader(new FileReader(dir + pedfile));
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        if (Integer.parseInt(line[5]) > 0) {
          if (hash.containsKey(line[0] + "->" + line[1])) {
            System.err.println("Warning - multiple entries for individual " + line[0] + "-"
                               + line[1] + " in the .fam file (only the last one will be used)");
          }
          hash.put(line[0] + "->" + line[1],
                   line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5]);
        }
        count++;
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + pedfile + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + pedfile + "\"");
      System.exit(2);
    }

    if (new File(dir + markersToTest).exists()) {
      System.out.println("Extracting the markers listed in " + markersToTest);
      subset = Array.toStringArray(HashVec.loadFileToVec(dir + markersToTest, false, true, false));
      allMarkers = false;
    } else {
      System.out.println("Did not find a " + markersToTest + " file; extracting all markers...");
      allMarkers = true;
    }

    try {
      FileReader in;
      int c = -1;
      String[][] data;
      int inc = 0;
      String trav;
      int rowCount, step;

      // try {
      // reader = new BufferedReader(new FileReader(dir+machDosage));
      // System.out.println("There are "+reader.readLine().trim().split("[\\s]+").length+" columns
      // in the first row");
      //
      // reader.close();
      // } catch (FileNotFoundException fnfe) {
      // System.err.println("Error: file \""+dir+machDosage+"\" not found in current directory");
      // System.exit(1);
      // } catch (IOException ioe) {
      // System.err.println("Error reading file \""+dir+machDosage+"\"");
      // System.exit(2);
      // }
      in = new FileReader(dir + machDosage);
      rowCount = 0;
      count = 0;
      c = 0;
      v = new Vector<String>();
      while (c >= 0) {
        if (count == 0) {
          trav = "";
          while (ext.indexOfChar((char) (c = in.read()), DELIMITERS) == -1 && c >= 0) {
            trav += (char) c;
          }
          if (c >= 0) {
            v.add(trav);
          }
          count = 1;
        }

        c = in.read();
        if (ext.indexOfChar((char) c, DELIMITERS) >= 0) {
          count++;
        }
        if ((char) c == '\n') {
          rowCount++;
          if (count - 1 != markerNames.length) {
            System.err.println("Error - mismatched number of columns in row " + (rowCount + 1));
            System.exit(1);
          }
          count = 0;
        }
      }
      in.close();

      System.out.println("Transposing data for " + rowCount + " individuals");
      step = 0;
      if (step == 0) {
        int bytesPerMarker = 7 * rowCount; // should be 5x, but it's actually closer to 7x
        long heapMaxSize = Runtime.getRuntime().maxMemory();
        step = (int) ((heapMaxSize * 0.75) / bytesPerMarker);
        System.out.println("Heap size (memory) was set to "
                           + ext.formDeci((double) heapMaxSize / 1024 / 1024 / 1024, 2)
                           + " Gb; which allows us to do " + step
                           + " markers at a time comfortably");
        System.out.println("Suggestion: Use 90% of available memory (i.e. -Xmx15g if you have 16GB of memory available) and use a 64-bit version of java (and the -d64 option) if you plan to use more than 2GB of memory");
      }
      int numCycles = (int) Math.ceil((double) markerNames.length / (double) step);
      long freeMemory = Runtime.getRuntime().freeMemory();
      System.out.println("Max memory: " + Runtime.getRuntime().maxMemory());
      Runtime.getRuntime().gc();
      Runtime.getRuntime().runFinalization();
      Runtime.getRuntime().gc();
      if (step < markerNames.length) {
        System.out.println("Running using " + numCycles + " cycles (" + step + " markers per run)");
      } else {
        System.out.println("Running and doing all " + markerNames.length
                           + " markers in a single pass");
        step = markerNames.length + 2;
      }

      writer = new PrintWriter(new FileWriter(dir + machDosage + ".dose"));
      writer.print("SNP\tA1\tA2");
      indIDs = Array.toStringArray(v);
      for (int i = 0; i < rowCount; i++) {
        writer.print("\t" + indIDs[i].substring(0, indIDs[i].indexOf("->")) + " "
                     + indIDs[i].substring(indIDs[i].indexOf("->") + 2));
      }
      writer.println();
      mapWriter = new PrintWriter(new FileWriter(dir + machDosage + ".map"));
      data = new String[rowCount][step];
      while (inc < markerNames.length) {
        System.out.println(ext.getTime() + " " + inc);

        in = new FileReader(dir + machDosage);
        for (int i = 0; i < rowCount; i++) {
          count = -2; // one for ID and one for "MLDOSE"
          while (count < inc) {
            if (ext.indexOfChar((char) (c = in.read()), DELIMITERS) >= 0) {
              count++;
            }
          }
          for (int j = 0; j < (inc + step > markerNames.length ? markerNames.length - inc
                                                               : step); j++) {
            trav = "";
            while (ext.indexOfChar((char) (c = in.read()), DELIMITERS) == -1 && (char) c != '\n') {
              trav += (char) c;
            }
            if (trav.endsWith("\r")) {
              trav = trav.substring(0, trav.length() - 1);
            }
            data[i][j] = trav;
            count++;
          }
          while (count < markerNames.length) {
            if (ext.indexOfChar((char) (c = in.read()), DELIMITERS) >= 0 || (char) c == '\n') {
              count++;
            }
          }
          if ((char) c != '\n') {
            System.err.println("Error transposing file: did not find " + (markerNames.length + 2)
                               + " columns for row " + (i + 1));
            System.err.println("   Make sure that the delimiter is set correctly");
          }
        }

        for (int i = 0; i < (inc + step > markerNames.length ? markerNames.length - inc
                                                             : step); i++) {
          if (allMarkers || ext.indexOfStr(markerNames[inc + i].split("[\\s]+")[0], subset) >= 0) {
            writer.print(markerNames[inc + i]);
            for (String[] element : data) {
              writer.print("\t" + element[i]);
            }
            writer.println();
            writer.flush();
            mapWriter.println("0\t" + markerNames[inc + i] + "\t0\t0");
            mapWriter.flush();
          }
        }
        inc += step;
      }
      in.close();
      writer.close();

      try {
        writer = new PrintWriter(new FileWriter(dir + machDosage + ".fam"));
        for (String indID : indIDs) {
          writer.print(indID.substring(0, indID.indexOf("->")) + "\t"
                       + indID.substring(indID.indexOf("->") + 2));
          if (hash.containsKey(indID)) {
            writer.println("\t" + hash.get(indID));
          } else {
            writer.println("\t0\t0\t0\t-9");
            System.err.println("Error - missing gender and phenotype for " + indID);
          }
        }
        writer.close();
      } catch (Exception e) {
        System.err.println("Error writing to " + machDosage + ".fam");
        e.printStackTrace();
      }
      System.out.println("Looks like we used " + (Runtime.getRuntime().maxMemory() - freeMemory)
                         + " in this run");
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + machDosage + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + machDosage + "\"");
      System.exit(2);
    }
  }

  // returns true if nothing went wrong
  public static boolean extractSpecificMarkers(String dir, String markerList, String dosageFormat,
                                               String markerInfoFormat, boolean verbose,
                                               Logger log) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    Hashtable<String, Integer> indexHash;
    Vector<String> v;
    String[] indIDs;
    int count;
    String[] markerNames;
    int[] indices;
    Hashtable<String, Vector<String>> chrHash;
    Hashtable<String, int[]> indicesHash;
    Hashtable<String, String[]> infoHash;
    String[] chrKeys;
    Integer index;
    int chr;
    boolean prob;
    long time;

    // FileReader in;
    InputStreamReader in;
    int c = -1;
    String[][] data;
    String trav;

    String dosageFile, markerInfoFile;
    String pedfile;
    boolean machFormat;
    String filename, descriptor, delimiter;

    time = new Date().getTime();
    System.out.println("Beginning extraction at " + ext.getTime());

    trav = dosageFormat.endsWith(".gz") ? dosageFormat.substring(0, dosageFormat.lastIndexOf("."))
                                        : dosageFormat;
    machFormat = trav.endsWith(".mldose");
    trav = markerInfoFormat.endsWith(".gz")
                                            ? markerInfoFormat.substring(0,
                                                                         markerInfoFormat.lastIndexOf("."))
                                            : markerInfoFormat;
    if (machFormat && !trav.endsWith(".mlinfo")) {
      System.err.println("Error - mismatched format patterns, assuming these are output from minimac");
    }

    indexHash = new Hashtable<String, Integer>();
    chrHash = new Hashtable<String, Vector<String>>();
    v = new Vector<String>();
    try {
      // reader = new BufferedReader(new FileReader(markerList));
      reader = Files.getAppropriateReader(markerList);
      count = 0;
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        if (line.length < 2) {
          log.reportError("Error - the second column of the marker list must contain the chromosome number");
          return false;
        }
        if (indexHash.containsKey(line[0])) {
          log.reportError("Error - cannot have the same marker (" + line[0]
                          + ") in the file twice");
          return false;
        }
        indexHash.put(line[0], Integer.valueOf(count));
        HashVec.addToHashVec(chrHash, line[1], line[0], false);
        v.add(line[0]);
        count++;
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + markerList + "\" not found in current directory");
      return false;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + markerList + "\"");
      return false;
    }
    markerNames = Array.toStringArray(v);
    log.report("Will try to export data for " + markerNames.length + " markers");

    pedfile = dir + "list.txt";
    if (!new File(pedfile).exists()) {
      log.reportError("Could not find 'list.txt' generating...");
      int chrom = 22;
      while (!(new File(ext.insertNumbers(dosageFormat, chrom)).exists()) && chrom >= 1) {
        chrom--;
      }
      if (chrom == 0) {
        log.reportError("Error - no valid " + dosageFormat + " file found anywhere; aborting");
        return false;
      }
      listIndividualsInMldose(ext.insertNumbers(dosageFormat, chrom), dir + "list.txt");
    }
    indIDs = HashVec.loadFileToStringArray(pedfile, false, false, new int[] {0, 1}, true, false,
                                           "[\\s]+");
    log.report("Will be pulling information from " + chrHash.size() + " chromosomes");


    chrKeys = HashVec.getKeys(chrHash);
    indicesHash = new Hashtable<String, int[]>();
    infoHash = new Hashtable<String, String[]>();
    for (String chrKey : chrKeys) {
      try {
        chr = Integer.parseInt(chrKey);
      } catch (NumberFormatException nfe) {
        log.reportError("Error - invalid chromosome number listed in column 2 for marker '"
                        + chrHash.get(chrKey).elementAt(0) + "'");
        chr = -1;
        return false;
      }
      if (chr < 1 || chr > 22) {
        log.reportError("Error - invalid chromosome number ('" + chr
                        + "') listed in column 2 for marker '" + chrHash.get(chrKey).elementAt(0)
                        + "'");
      }
      dosageFile = ext.insertNumbers(dosageFormat, chr);
      markerInfoFile = ext.insertNumbers(markerInfoFormat, chr);

      if (!Files.exists(dosageFile, false)) {
        log.reportError("Error - missing file '" + dosageFile + "'");
      }
      if (!Files.exists(markerInfoFile, false)) {
        log.reportError("Error - missing file '" + markerInfoFile + "'");
      }

      log.report("Parsing indices for chr" + chr);
      indices = new int[Files.countLines(markerInfoFile, 1)];
      try {
        reader = Files.getAppropriateReader(markerInfoFile);
        line = reader.readLine().trim().split("[\\s]+");
        if (machFormat) {
          ext.checkHeader(line, MLINFO_HEADER, true);
        } else {
          ext.checkHeader(line, MINFO_HEADER, true);
        }
        count = 0;
        while (reader.ready()) {
          line = reader.readLine().trim().split("[\\s]+");
          index = indexHash.get(line[0]);
          if (index == null) {
            indices[count] = -1;
          } else {
            indices[count] = index.intValue();
            infoHash.put(line[0], line);
          }
          count++;
        }
        reader.close();
      } catch (Exception e) {
        log.reportError("Error reading file \"" + markerInfoFile + "\"");
        log.reportException(e);
        return false;
      }
      indicesHash.put(chrKey, indices);
    }

    prob = false;
    for (int i = 0; i < markerNames.length; i++) {
      if (!infoHash.containsKey(markerNames[i])) {
        log.reportError("Error - marker '" + markerNames[i] + "' was not found in any "
                        + markerInfoFormat + " file", true, verbose);
        prob = true;
      }
    }
    if (prob) {
      return false;
    }

    data = new String[indIDs.length][markerNames.length];
    log.report("Parsing " + markerNames.length + " markers for " + indIDs.length + " individuals",
               true, verbose);

    for (String chrKey : chrKeys) {
      // log.report(".", false, verbose);
      chr = Integer.parseInt(chrKey);
      dosageFile = ext.insertNumbers(dosageFormat, chr);
      indices = indicesHash.get(chrKey);
      System.out.println("Parsing " + indices.length + " markers for chromosome " + chrKey
                         + " from " + dosageFile);

      try {
        if (dosageFile.endsWith(".gz")) {
          in = new InputStreamReader(new GZIPInputStream(new FileInputStream(dosageFile)));
        } else {
          in = new FileReader(dosageFile);
        }
        count = -2;
        c = 0;
        for (int ind = 0; ind < indIDs.length; ind++) {
          trav = "";
          while (ext.indexOfChar((char) (c = in.read()), DELIMITERS) == -1) {
            trav += (char) c;
          }
          line = indIDs[ind].split("[\\s]+");
          if (!trav.equals(line[0] + "->" + line[1])) {
            log.reportError("Error - dosagee file (" + dosageFile + ") does not match pedfile ("
                            + pedfile + "); expecting '" + line[0] + "->" + line[1] + "', found '"
                            + trav + "'");
            return false;
          }

          while (ext.indexOfChar((char) (c = in.read()), DELIMITERS) == -1 && c >= 0) {
            ; // MLDOSE column
          }

          for (int j = 0; j < indices.length; j++) {
            if (indices[j] == -1) {
              while (ext.indexOfChar((char) (c = in.read()), DELIMITERS) == -1
                     && (char) c != '\n') {
                ;
              }
            } else {
              trav = "";
              while (ext.indexOfChar((char) (c = in.read()), DELIMITERS) == -1
                     && (char) c != '\n') {
                trav += (char) c;
              }
              if (trav.endsWith("\r")) {
                trav = trav.substring(0, trav.length() - 1);
              }
              data[ind][indices[j]] = trav;
            }
          }
          if ((char) c != '\n') {
            log.reportError("Error reading '" + dosageFile + "': did not find "
                            + (indices.length + 2) + " columns for row " + (ind + 1));
            log.reportError("   Make sure that the delimiter is set correctly");
          }
          if (c < 0 && ind != indIDs.length - 1) {
            log.reportError("Error - premature truncation of file '" + dosageFile + "' at line "
                            + (ind + 1));
            return false;
          }
        }
        if (in.read() != -1) {
          log.reportError("Error - stopped reading before the end of file for '" + dosageFile
                          + "' (found character '" + c + "')");
          in.close();
          return false;
        }
        in.close();
      } catch (IOException ioe) {
        log.reportError("Error reading file \"" + dosageFile + "\"");
        return false;
      }
    }
    log.report("", true, verbose);

    if (machFormat) {
      filename = ext.rootOf(markerList, false) + ".mldose";
      descriptor = "MLDOSE";
      delimiter = " ";
    } else {
      filename = ext.rootOf(markerList, false) + ".dose";
      descriptor = "DOSE";
      delimiter = "\t";
    }
    try {
      writer = new PrintWriter(new FileWriter(filename));
      for (int i = 0; i < indIDs.length; i++) {
        line = indIDs[i].split("[\\s]+");
        writer.print(line[0] + "->" + line[1] + delimiter + descriptor);
        for (int j = 0; j < markerNames.length; j++) {
          writer.print(delimiter + data[i][j]);
        }
        writer.println();
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + filename);
      log.reportException(e);
    }

    filename = ext.rootOf(markerList, false) + (machFormat ? ".mlinfo" : ".info");
    try {
      writer = new PrintWriter(new FileWriter(filename));
      if (machFormat) {
        writer.println(Array.toStr(MLINFO_HEADER));
      } else {
        writer.println(Array.toStr(MINFO_HEADER));
      }
      for (String markerName : markerNames) {
        writer.println(Array.toStr(infoHash.get(markerName)));
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + filename);
      log.reportException(e);
    }
    log.report("Finished in " + ext.getTimeElapsed(time), true, verbose);

    return true;
  }

  // // using reader.readline does not appear to be all that faster or slower
  // public static boolean extractSpecificMarkers(String dir, String markerList, String
  // dosageFormat, String markerInfoFormat, boolean verbose, Logger log) {
  // BufferedReader reader;
  // PrintWriter writer;
  // String[] line, ids;
  // Hashtable<String,Integer> indexHash;
  // Vector<String> v;
  // String[] indIDs;
  // int count;
  // String[] markerNames;
  // int[] indices;
  // Hashtable<String,Vector<String>> chrHash;
  // Hashtable<String,int[]> indicesHash;
  // Hashtable<String,String[]> infoHash;
  // String[] chrKeys;
  // Integer index;
  // int chr;
  // boolean prob;
  // long time;
  //
  //// FileReader in;
  // InputStreamReader in;
  // String[][] data;
  // String trav;
  //
  // String dosageFile, markerInfoFile;
  // String pedfile;
  // boolean machFormat;
  // String filename, descriptor, delimiter;
  //
  // time = new Date().getTime();
  // System.out.println("Beginning at "+ext.getTime());
  //
  // trav = dosageFormat.endsWith(".gz")?dosageFormat.substring(0,
  // dosageFormat.lastIndexOf(".")):dosageFormat;
  // machFormat = trav.endsWith(".mldose");
  // trav = markerInfoFormat.endsWith(".gz")?markerInfoFormat.substring(0,
  // markerInfoFormat.lastIndexOf(".")):markerInfoFormat;
  // if (machFormat && !trav.endsWith(".mlinfo")) {
  // System.err.println("Error - mismatched format patterns, assuming these are output from
  // minimac");
  // }
  //
  // indexHash = new Hashtable<String,Integer>();
  // chrHash = new Hashtable<String,Vector<String>>();
  // v = new Vector<String>();
  // try {
  //// reader = new BufferedReader(new FileReader(markerList));
  // reader = Files.getAppropriateReader(markerList);
  // count = 0;
  // while (reader.ready()) {
  // line = reader.readLine().trim().split("[\\s]+");
  // if (line.length < 2) {
  // log.reportError("Error - the second column of the marker list must contain the chromosome
  // number");
  // return false;
  // }
  // if (indexHash.containsKey(line[0])) {
  // log.reportError("Error - cannot have the same marker ("+line[0]+") in the file twice");
  // return false;
  // }
  // indexHash.put(line[0], Integer.valueOf(count));
  // HashVec.addToHashVec(chrHash, line[1], line[0], false);
  // v.add(line[0]);
  // count++;
  // }
  // reader.close();
  // } catch (FileNotFoundException fnfe) {
  // log.reportError("Error: file \""+markerList+"\" not found in current directory");
  // return false;
  // } catch (IOException ioe) {
  // log.reportError("Error reading file \""+markerList+"\"");
  // return false;
  // }
  // markerNames = Array.toStringArray(v);
  // log.report("Will try to export data for "+markerNames.length+" markers");
  //
  // pedfile = dir+"list.txt";
  // if (!new File(pedfile).exists()) {
  // log.reportError("Could not find 'list.txt' generating...");
  // int chrom = 22;
  // while (!(new File(ext.insertNumbers(dosageFormat, chrom)).exists()) && chrom >= 1) {
  // chrom--;
  // }
  // if (chrom == 0) {
  // log.reportError("Error - no valid "+dosageFormat+" file found anywhere; aborting");
  // return false;
  // }
  // listIndividualsInMldose(ext.insertNumbers(dosageFormat, chrom), dir+"list.txt");
  // }
  // indIDs = HashVec.loadFileToStringArray(pedfile, false, false, new int[] {0,1}, true, false,
  // "[\\s]+");
  // log.report("Will be pulling information from "+chrHash.size()+" chromosomes");
  //
  //
  // chrKeys = HashVec.getKeys(chrHash);
  // indicesHash = new Hashtable<String,int[]>();
  // infoHash = new Hashtable<String,String[]>();
  // for (int i = 0; i<chrKeys.length; i++) {
  // try {
  // chr = Integer.parseInt(chrKeys[i]);
  // } catch (NumberFormatException nfe) {
  // log.reportError("Error - invalid chromosome number listed in column 2 for marker
  // '"+chrHash.get(chrKeys[i]).elementAt(0)+"'");
  // chr = -1;
  // return false;
  // }
  // if (chr < 1 || chr > 22) {
  // log.reportError("Error - invalid chromosome number ('"+chr+"') listed in column 2 for marker
  // '"+chrHash.get(chrKeys[i]).elementAt(0)+"'");
  // }
  // dosageFile = ext.insertNumbers(dosageFormat, chr);
  // markerInfoFile = ext.insertNumbers(markerInfoFormat, chr);
  //
  // if (!Files.exists(dosageFile, false)) {
  // log.reportError("Error - missing file '"+dosageFile+"'");
  // }
  // if (!Files.exists(markerInfoFile, false)) {
  // log.reportError("Error - missing file '"+markerInfoFile+"'");
  // }
  //
  // log.report("Parsing indices for chr"+chr);
  // indices = new int[Files.countLines(markerInfoFile, true)];
  // try {
  // reader = Files.getAppropriateReader(markerInfoFile);
  // line = reader.readLine().trim().split("[\\s]+");
  // if (machFormat) {
  // ext.checkHeader(line, MLINFO_HEADER, true);
  // } else {
  // ext.checkHeader(line, MINFO_HEADER, true);
  // }
  // count = 0;
  // while (reader.ready()) {
  // line = reader.readLine().trim().split("[\\s]+");
  // index = indexHash.get(line[0]);
  // if (index == null) {
  // indices[count] = -1;
  // } else {
  // indices[count] = index.intValue();
  // infoHash.put(line[0], line);
  // }
  // count++;
  // }
  // reader.close();
  // } catch (Exception e) {
  // log.reportError("Error reading file \""+markerInfoFile+"\"");
  // log.reportException(e);
  // return false;
  // }
  // indicesHash.put(chrKeys[i], indices);
  // }
  //
  // prob = false;
  // for (int i = 0; i<markerNames.length; i++) {
  // if (!infoHash.containsKey(markerNames[i])) {
  // log.reportError("Error - marker '"+markerNames[i]+"' was not found in any "+markerInfoFormat+"
  // file", true, verbose);
  // prob = true;
  // }
  // }
  // if (prob) {
  // return false;
  // }
  //
  // data = new String[indIDs.length][markerNames.length];
  // log.report("Parsing "+markerNames.length+" markers for "+indIDs.length+" individuals", true,
  // verbose);
  //
  // for (int i = 0; i<chrKeys.length; i++) {
  //// log.report(".", false, verbose);
  // chr = Integer.parseInt(chrKeys[i]);
  // dosageFile = ext.insertNumbers(dosageFormat, chr);
  // indices = indicesHash.get(chrKeys[i]);
  // System.out.println("Parsing "+indices.length+" markers for chromosome "+chrKeys[i]+" from
  // "+dosageFile);
  //
  // try {
  // if (dosageFile.endsWith(".gz")) {
  // in = new InputStreamReader(new GZIPInputStream(new FileInputStream(dosageFile)));
  // } else {
  // in = new FileReader(dosageFile);
  // }
  // count = -2;
  // for (int ind = 0; ind<indIDs.length; ind++) {
  // if (!reader.ready()) {
  // log.reportError("Error - premature truncation of file '"+dosageFile+"' at line "+(ind+1));
  // in.close();
  // return false;
  // }
  // line = reader.readLine().trim().split("[\\s]+");
  // ids = indIDs[ind].split("[\\s]+");
  // if (!line[0].equals(ids[0]+"->"+ids[1])) {
  // log.reportError("Error - dosagee file ("+dosageFile+") does not match pedfile ("+pedfile+");
  // expecting '"+line[0]+"->"+line[1]+"', found '"+trav+"'");
  // in.close();
  // return false;
  // }
  //
  // if (line.length != indices.length+2) {
  // log.reportError("Error reading '"+dosageFile+"': did not find "+(indices.length+2)+" columns
  // for row "+(ind+1));
  // log.reportError(" Make sure that the delimiter is set correctly");
  // }
  //
  //// line[1] == 'MLDOSE' column
  //
  // for (int j = 0; j<indices.length; j++) {
  // if (indices[j] != -1) {
  // data[ind][indices[j]] = line[2+j];
  // }
  // }
  // }
  // in.close();
  // } catch (IOException ioe) {
  // log.reportError("Error when reading file \""+dosageFile+"\"");
  // log.reportException(ioe);
  // return false;
  // } catch (Exception e) {
  // log.reportError("Error processing file \""+dosageFile+"\"");
  // log.reportException(e);
  // return false;
  // }
  // }
  // log.report("", true, verbose);
  //
  // if (machFormat) {
  // filename = ext.rootOf(markerList, false)+".mldose";
  // descriptor = "MLDOSE";
  // delimiter = " ";
  // } else {
  // filename = ext.rootOf(markerList, false)+".dose";
  // descriptor = "DOSE";
  // delimiter = "\t";
  // }
  // try {
  // writer = new PrintWriter(new FileWriter(filename));
  // for (int i = 0; i<indIDs.length; i++) {
  // line = indIDs[i].split("[\\s]+");
  // writer.print(line[0]+"->"+line[1]+delimiter+descriptor);
  // for (int j = 0; j<markerNames.length; j++) {
  // writer.print(delimiter+data[i][j]);
  // }
  // writer.println();
  // }
  // writer.close();
  // } catch (Exception e) {
  // log.reportError("Error writing to "+filename);
  // log.reportException(e);
  // }
  //
  // filename = ext.rootOf(markerList, false)+(machFormat?".mlinfo":".info");
  // try {
  // writer = new PrintWriter(new FileWriter(filename));
  // if (machFormat) {
  // writer.println(Array.toStr(MLINFO_HEADER));
  // } else {
  // writer.println(Array.toStr(MINFO_HEADER));
  // }
  // for (int i = 0; i<markerNames.length; i++) {
  // writer.println(Array.toStr(infoHash.get(markerNames[i])));
  // }
  // writer.close();
  // } catch (Exception e) {
  // log.reportError("Error writing to "+filename);
  // log.reportException(e);
  // }
  // log.report("Finished in "+ext.getTimeElapsed(time), true, verbose);
  //
  // return true;
  // }

  public static void extractIntermediate(String markerList, int regionNameIndex,
                                         String dosageFormat, String infoFormat) {
    BufferedReader reader;
    PrintWriter writer, w2;
    String[] line;
    Hashtable<String, Integer> chrHash;
    Hashtable<String, String> markerHash;
    Vector<String> v;
    int count;
    Hashtable<String, Hashtable<String, String>> regionMarkerHashes;
    Hashtable<String, String> regionHash;
    String[] regionKeys;
    int chr;
    long time;
    String dosageFile, markersFile;
    String batchFilename, dir;
    boolean regions;

    if (regionNameIndex == -1) {
      regions = false;
      regionNameIndex = 1;
    } else {
      regions = true;
    }

    time = new Date().getTime();

    markerHash = new Hashtable<String, String>();
    chrHash = new Hashtable<String, Integer>();
    regionMarkerHashes = new Hashtable<String, Hashtable<String, String>>();
    v = new Vector<String>();
    try {
      reader = new BufferedReader(new FileReader(markerList));
      System.out.println("Reading data from " + markerList);
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        if (line.length < 3) {
          System.err.println("Error - the second column of the marker list must contain the chromosome number");
          System.exit(1);
        } else if (line.length < regionNameIndex + 1) {
          System.err.println("Error - each row must contain the region name in the specified index ('"
                             + regionNameIndex + "')");
          System.exit(1);
        }
        if (!regions && markerHash.containsKey(line[0])) {
          System.err.println("Error - cannot have the same marker (" + line[0]
                             + ") in the file twice");
          System.exit(1);
        }
        markerHash.put(line[0], "");
        try {
          chrHash.put(line[0], Integer.valueOf(line[1]));
        } catch (NumberFormatException nfe) {
          System.err.println("Error - invalid chromosome number listed in column 2 for marker '"
                             + line[0] + "'");
          System.exit(1);
        }
        line[regionNameIndex] = ext.replaceAllWith(line[regionNameIndex], ":", "_");
        HashVec.addToHashHash(regionMarkerHashes, line[regionNameIndex], line[0], "");
        HashVec.addIfAbsent(line[regionNameIndex], v);
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + markerList + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + markerList + "\"");
      System.exit(2);
    }

    if (regions) {
      batchFilename = "splitByRegion";
      // regionKeys = HashVec.getKeys(regionMarkerHashes, true, false);
      regionKeys = Array.toStringArray(v);
      Files.writeList(regionKeys, "regionList.dat");
      System.out.println("Found " + regionKeys.length + " regions to parse");
    } else {
      batchFilename = "intExtract";
      regionKeys = HashVec.getKeys(regionMarkerHashes, true, true);
      System.out.println("Markers are distributed across " + regionKeys.length + " chromosomes");
    }

    try {
      writer = new PrintWriter(new FileWriter(batchFilename));
      for (int i = 0; i < regionKeys.length; i++) {
        regionHash = regionMarkerHashes.get(regionKeys[i]);
        chr = chrHash.get(regionHash.keys().nextElement()).intValue();
        if (regions) {
          dir = regionKeys[i] + "/";
          System.out.println(regionKeys[i] + " on chr" + chr + " (region " + (i + 1) + " of "
                             + regionKeys.length + ")");
        } else {
          dir = "intermediateExtraction/";
          System.out.println("chr" + chr + " (batch " + (i + 1) + " of " + regionKeys.length + ")");
        }
        new File(dir).mkdir();
        if (chr < 1 || chr > 22) {
          System.err.println("Error - invalid chromosome number ('" + chr
                             + "') listed in column 2 for marker '"
                             + regionHash.keys().nextElement() + "'");
        }

        dosageFile = ext.insertNumbers(dosageFormat, chr);
        markersFile = ext.insertNumbers(infoFormat, chr);

        if (!Files.exists(dosageFile, false)) {
          System.err.println("Error - missing file '" + dosageFile + "'");
          System.exit(1);
        }
        if (!Files.exists(markersFile, false)) {
          System.err.println("Error - missing file '" + markersFile + "'");
          System.exit(1);
        }

        try {
          reader = Files.getAppropriateReader(markersFile);
          w2 =
             new PrintWriter(new FileWriter(dir + (regions ? "data.info" : "chr" + chr + ".info")));
          line = reader.readLine().trim().split("[\\s]+");
          w2.println(Array.toStr(line));
          if (markersFile.endsWith(".gz") || markersFile.endsWith(".zip")) {
            markersFile = markersFile.substring(0, markersFile.lastIndexOf("."));
          }
          if (markersFile.endsWith(".mlinfo")) {
            System.out.println("Checking against MLINFO header");
            ext.checkHeader(line, MLINFO_HEADER, true);
          } else if (markersFile.endsWith(".info") || markersFile.endsWith(".minfo")
                     || markersFile.endsWith(".pinfo")) {
            System.out.println("Checking against INFO header");
            ext.checkHeader(line, MINFO_HEADER, true);
          } else {
            System.err.println("Warning - unknown marker info format extension (" + markersFile
                               + "); skipping header verification");
          }
          v = new Vector<String>();
          v.add("1");
          v.add("2");
          count = 2;
          while (reader.ready()) {
            count++;
            line = reader.readLine().trim().split("[\\s]+");
            if (regionHash.containsKey(line[0])) {
              regionHash.remove(line[0]);
              v.add(count + "");
              w2.println(Array.toStr(line));
            }
          }
          reader.close();
          w2.close();

          if (regionHash.size() > 0) {
            System.err.println("Error - missed " + regionHash.size()
                               + " markers that were supposed to be "
                               + (regions ? "part of region " : "on chromosome ") + regionKeys[i]
                               + ":");
            System.err.println(Array.toStr(HashVec.getKeys(regionHash), ","));
          }
        } catch (IOException ioe) {
          System.err.println("Error reading file \"" + markersFile + "\"");
          System.exit(2);
        }
        if (regions) {
          writer.println("echo \"" + regionKeys[i] + " on chr" + chr + " (region " + (i + 1)
                         + " of " + regionKeys.length + ")...\"");
        } else {
          writer.println("echo \"parsing chr" + chr + " (batch " + (i + 1) + " of "
                         + regionKeys.length + ")...\"");
        }
        Files.write("{print $" + Array.toStr(Array.toStringArray(v), "\"\\t\"$") + "}",
                    dir + "split" + (i + 1) + ".awk");
        if (dosageFile.endsWith(".gz")) {
          writer.println("gunzip -c " + dosageFile + " | awk -f " + dir + "split" + (i + 1)
                         + ".awk > " + dir + (regions ? "data.dose" : "chr" + chr + ".dose"));
        } else {
          writer.println("awk -f " + dir + "split" + (i + 1) + ".awk " + dosageFile + " > " + dir
                         + (regions ? "data.dose" : "chr" + chr + ".dose"));
        }
      }
      if (regions) {
        writer.println("tar -zcvf intermediateExtraction/region_awks.tar.gz intermediateExtraction/*.awk");
        writer.println("rm intermediateExtraction/*.awk");
      } else {
        writer.println("tar -zcvf intermediateExtraction/intEx_awks.tar.gz intermediateExtraction/*.awk");
        writer.println("rm intermediateExtraction/*.awk");
      }
      writer.close();
      Files.chmod(batchFilename);
    } catch (Exception e) {
      System.err.println("Error writing to " + batchFilename);
      e.printStackTrace();
    }

    System.out.println("Finished setting up " + (regions ? "split-by-region" : "intermediate")
                       + " extraction in " + ext.getTimeElapsed(time)
                       + "; use the following command to have awk do the actual parsing:");
    System.out.println("./" + batchFilename);
  }

  public static void extractGenotypesForPlink(String dir, String pedfile, String machProb,
                                              String machMarkers, String markersToTest,
                                              double rsqThreshold, double probRequirement) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    Hashtable<String, String> hash, subset;
    int count;
    String[] markerNames;
    boolean allMarkers;

    SnpMarkerSet map;
    boolean[] use;
    String[] rsqs;
    char[][] alleles;

    if (!dir.equals("") && !new File(dir).exists()) {
      System.err.println("Could not find directory: " + dir);
      System.err.println("  using current cirectory instead");
      dir = "";
    }

    if (new File(dir + markersToTest).exists()) {
      System.out.println("Extracting the markers listed in " + markersToTest);
      subset = HashVec.loadFileToHashString(dir + markersToTest, 0, new int[] {0}, "", false);
      allMarkers = false;
    } else {
      System.out.println("Did not find a " + markersToTest + " file; extracting all markers...");
      allMarkers = true;
      subset = null;
    }

    map = new SnpMarkerSet(dir + machMarkers);

    markerNames = map.getMarkerNames();
    alleles = map.getAlleles();
    use = new boolean[markerNames.length];
    rsqs = Matrix.extractColumn(map.getAnnotation(), 3);
    try {
      writer = new PrintWriter(new FileWriter(dir + machProb + ".map"));
      for (int i = 0; i < use.length; i++) {
        use[i] = (allMarkers || subset.containsKey(markerNames[i]))
                 && Double.parseDouble(rsqs[i]) > rsqThreshold;
        if (use[i]) {
          writer.println("0\t" + markerNames[i] + "\t0\t0");
        }
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + machProb + ".map");
      e.printStackTrace();
    }

    count = 0;
    hash = new Hashtable<String, String>();
    try {
      reader = new BufferedReader(new FileReader(dir + pedfile));
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        if (Integer.parseInt(line[5]) > 0) {
          if (hash.containsKey(line[0] + "->" + line[1])) {
            System.err.println("Warning - multiple entries for individual " + line[0] + "-"
                               + line[1] + " in the .fam file (only the last one will be used)");
          }
          hash.put(line[0] + "->" + line[1],
                   line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5]);
        }
        count++;
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + pedfile + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + pedfile + "\"");
      System.exit(2);
    }

    try {
      FileReader in;
      int c = -1;
      String trav;
      int rowCount;
      double[] probs;

      in = new FileReader(dir + machProb);
      writer = new PrintWriter(new FileWriter(dir + machProb + ".ped"));
      rowCount = 0;
      count = -2;
      c = 0;
      probs = new double[2];
      while (c >= 0) {
        if (count == -2) {
          trav = "";
          while (ext.indexOfChar((char) (c = in.read()), DELIMITERS) == -1 && c >= 0) {
            trav += (char) c;
          }
          if (c >= 0) { // final run does pass through here
            writer.print(trav.substring(0, trav.indexOf("->")) + "\t"
                         + trav.substring(trav.indexOf("->") + 2));
            if (hash.containsKey(trav)) {
              writer.print("\t" + hash.get(trav));
            } else {
              writer.print("\t0\t0\t0\t-9");
              System.err.println("Error - missing gender and phenotype for " + trav);
            }
          }
          count++;

          while (c >= 0 && count < 0) { // and here
            if (ext.indexOfChar((char) in.read(), DELIMITERS) >= 0) {
              count++;
            }
          }
        }

        if (c >= 0) {
          for (int i = 0; i < 2; i++) {
            trav = "";
            while (ext.indexOfChar((char) (c = in.read()), DELIMITERS) == -1 && (char) c != '\n') {
              trav += (char) c;
            }
            if (trav.endsWith("\r")) {
              trav = trav.substring(0, trav.length() - 1);
            }
            probs[i] = Double.parseDouble(trav);
          }

          if (use[count]) {
            if (probs[0] >= probRequirement) {
              writer.print("\t" + alleles[count][0] + "\t" + alleles[count][0]);
            } else if (probs[1] >= probRequirement) {
              writer.print("\t" + alleles[count][0] + "\t" + alleles[count][1]);
            } else if (1 - (probs[0] + probs[1]) >= probRequirement) {
              writer.print("\t" + alleles[count][1] + "\t" + alleles[count][1]);
            } else {
              writer.print("\t0\t0");
            }
          }
          count++;

          if ((char) c == '\n' || c < 0) {
            rowCount++;
            if (count != markerNames.length) {
              System.err.println("Error - mismatched number of columns/markers in row "
                                 + (rowCount + 1));
              System.exit(1);
            }
            count = -2;
            writer.println();
          }
        }
      }
      in.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + machProb + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + machProb + "\"");
      System.exit(2);
    }
  }

  public static void parseGenotypesTogether(String machRoot) {
    PrintWriter writer;

    for (int chr = 1; chr <= 22; chr++) {
      System.out.println(Files.getRunString()
                         + " gwas.Mach -extractGenotypes rsq=0.98 prob=0.98 mach=MACH_step2_chr"
                         + chr + " ped=../../plink.fam");
      CmdLine.run(Files.getRunString()
                  + " gwas.Mach -extractGenotypes rsq=0.98 prob=0.98 mach=MACH_step2_chr" + chr
                  + " ped=../../plink.fam", "chr" + chr + "/");
    }
    try {
      writer = new PrintWriter(new FileWriter("merge.dat"));
      for (int chr = 2; chr <= 22; chr++) {
        writer.println("chr" + chr + "/" + machRoot + "_chr" + chr + ".mlprob.ped chr" + chr + "/"
                       + machRoot + "_chr" + chr + ".mlprob.map");
      }
      System.out.println("plink --noweb --file chr1/" + machRoot
                         + "_chr1.mlprob --merge-list merge.dat --make-bed --out imputedGenotypes");

      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + "merge.dat");
      e.printStackTrace();
    }
  }

  public static void pickBest(int size) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    Hashtable<String, String> hash;
    Hashtable<String, Vector<String>> hashV;
    Vector<String> v = new Vector<String>();
    FamilyStructure struct;
    String[][] ids;
    byte[] genders, affections;
    String[] keys;
    double[] values;
    int[] order;

    struct = new FamilyStructure("plink.fam");
    ids = struct.getIDs();
    genders = struct.getGenders();
    affections = struct.getAffections();

    hash = new Hashtable<String, String>();
    try {
      reader = new BufferedReader(new FileReader("plink.imiss"));
      ext.checkHeader(reader.readLine().trim().split("[\\s]+"),
                      new String[] {"FID", "IID", "MISS_PHENO", "N_MISS", "N_GENO", "F_MISS"}, true,
                      true);
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        hash.put(line[0] + "\t" + line[1], line[3]);
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + "plink.imiss" + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + "plink.imiss" + "\"");
      System.exit(2);
    }

    hashV = new Hashtable<String, Vector<String>>();
    for (int i = 0; i < ids.length; i++) {
      HashVec.addToHashVec(hashV,
                           genders[i] + "_"
                                  + affections[i],
                           ids[i][0] + "\t" + ids[i][1] + "\t"
                                                   + hash.get(ids[i][0] + "\t" + ids[i][1]),
                           false);
    }

    try {
      writer = new PrintWriter(new FileWriter("subset.txt"));
      keys = HashVec.getKeys(hashV);
      for (String key : keys) {
        v = hashV.get(key);
        values = new double[v.size()];
        for (int j = 0; j < v.size(); j++) {
          line = v.elementAt(j).split("[\\s]+");
          values[j] = Double.parseDouble(line[2]);
        }
        order = Sort.quicksort(values);
        for (int j = 0; j < Math.min(order.length, size); j++) {
          line = v.elementAt(order[j]).split("[\\s]+");
          writer.println(line[0] + "\t" + line[1]);
        }
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + "subset.txt");
      e.printStackTrace();
    }
  }

  public static void createDatabase(String machRoot) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String[] markerNames;
    long time;

    time = new Date().getTime();

    markerNames = HashVec.loadFileToStringArray(machRoot + ".mlinfo", true, new int[] {0}, false);
    try {
      reader = new BufferedReader(new FileReader(machRoot + ".mldose"));
      writer = new PrintWriter(new FileWriter(machRoot + ".xln"));
      writer.println("FID\tIID\t" + Array.toStr(markerNames));
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        line[1] = line[0].substring(line[0].indexOf("->") + 2);
        line[0] = line[0].substring(0, line[0].indexOf("->"));
        writer.println(Array.toStr(line));
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + machRoot + ".mldose"
                         + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + machRoot + ".mldose" + "\"");
      System.exit(2);
    }

    System.out.println("Done in " + ext.getTimeElapsed(time));

  }

  public static void checkResults() {
    IntVector step1Pres, step1Miss, step2Pres, step2Miss;

    step1Pres = new IntVector();
    step1Miss = new IntVector();
    step2Pres = new IntVector();
    step2Miss = new IntVector();
    for (int i = 1; i <= 22; i++) {
      if (new File("chr" + i + "/MACH_step1_chr" + i + ".erate").exists()) {
        step1Pres.add(i);
      } else {
        step1Miss.add(i);
      }
      if (new File("chr" + i + "/MACH_step2_chr" + i + ".erate").exists()) {
        step2Pres.add(i);
      } else {
        step2Miss.add(i);
      }
    }
    System.out.println("Finished the following chromosomes for step1: "
                       + ext.listRanges(Ints.toArray(step1Pres)));
    System.out.println("Still missing: " + ext.listRanges(Ints.toArray(step1Miss)));
    System.out.println();
    System.out.println("Finished the following chromosomes for step2: "
                       + ext.listRanges(Ints.toArray(step2Pres)));
    System.out.println("Still missing: " + ext.listRanges(Ints.toArray(step2Miss)));

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String root = "plink";
    String keeps = "subset.txt";
    // String excludes = "CIDR_bad_snps_final.txt.prn";
    String excludes = null;
    String convert = "";
    boolean createfiles = false;
    boolean batch = false;
    // String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\candis\\Confidence
    // Intervals\\CI95_every_which_way\\";
    String dir = "";
    // String pedfile = "pd_gwas.fam";
    // String pedfile = "plink.fam";
    String pedfile = "confidence.ped";
    String machRoot = "MACH_step2";
    String markerSubset = "testThese.txt";
    int trim = -1;
    String list = null;
    int decode = -1;
    boolean extractDosage = false;
    boolean extractGenotypes = false;
    double rsqThreshold = 0.98;
    double probRequirement = 0.98;
    boolean extractAllGenotypes = false;
    boolean pickBest = false;
    int pickSize = PICK_BEST_SIZE;
    boolean checkResults = false;
    String extract = "";
    String database = "";
    String doseFormat = "chr#/MACH_step2_chr#.mldose";
    String infoFormat = "chr#/MACH_step2_chr#.mlinfo";
    String probFormat = "MACH_step2_chr21.mlprob";
    boolean intEx = false;
    String[] nodesToUse = null;
    int splitByRegion = -1;

    // dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\tWork\\Expression\\LD_versus_GWAS\\MAPT_again\\";
    // decode = 17;

    // dir = "C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\candis\\chr11\\";
    // pedfile = "confidence.ped";
    // machRoot = "MACH_step2_chr11";
    // extractGenotypes = true;

    // extractAllGenotypes = true;

    // extract = "list11.txt";

    String usage = "\n" + "gwas.Mach requires 0-1 arguments\n"
                   + "   (1) root of plink file system (i.e. root=" + root + " (default))\n"
                   + "   (2) individuals to keep (i.e. keeps=" + keeps + " (default))\n"
                   + "   (3) markers to exclude (i.e. exclude=" + excludes + " (default))\n"
                   + "   (4) convert .bim to .dat using prefix provided (i.e. convert=pd200 (not the default))\n"
                   + " OR\n" + "   -createFiles (creates batch that creates files for MACH)\n"
                   + " OR\n" + "   -gethapmap (downloads data)\n" + " OR\n"
                   + "   (1) create batch files for MACH phasing (i.e. -batches (not the default))\n"
                   + "   (2) nodes to use (i.e. nodesToUse="
                   + Array.toStr(new String[] {"v1", "v2", "..."}, ",")
                   + " (default; qsubs only; full names, comma-delimited))\n" + " OR\n"
                   + "   chromosome to trim from source phase files (i.e. trim=1 (not the default))\n"
                   + " OR\n" + "   list IDs from .mldose file (i.e. list=chr21.mldose)\n" + " OR\n"
                   + "   (1) trim the phased HapMap data for targetted imputation (i.e. trim=17 (not the default))\n"
                   + " OR\n"
                   + "   (1) decode trimmed phased HapMap data (i.e. decode=17 (not the default))\n"
                   + " OR\n"
                   + "   (1) extract PLINK dosage files from MACH .mldose file (i.e. -extractDosage (not the default))\n"
                   + "   (2) pedigree file (i.e. ped=" + pedfile + " (default))\n"
                   + "   (3) format of dosage filenames (i.e. doseFormat=" + doseFormat
                   + " (default))\n" + "   (4) format of marker info filenames (i.e. infoFormat="
                   + infoFormat + " (default))\n"
                   + "   (5) (optional) list of subset of markers to test (i.e. subset="
                   + markerSubset + " (default))\n" + " OR\n"
                   + "   (1) recreate .mldose/.mlinfo with a subset of markers (i.e. extract=markerList.txt (not the default; list requires chr in second column))\n"
                   + "   (2) format of dosage filenames (i.e. doseFormat=" + doseFormat
                   + " (default))\n" + "   (3) format of marker info filenames (i.e. infoFormat="
                   + infoFormat + " (default))\n"
                   + "   (4) (optional) create batch to quickly parse .mldose files (i.e. -intEx (not the default))\n"
                   + "   (5) (optional) split files up by region name into separate directories (i.e. splitByRegion=[index of region name in markerList file] (not the default))\n"
                   + " OR\n"
                   + "   (1) create database .xln file from a .mldose/.mlinfo pair (i.e. database=fileRoot (not the default))\n"
                   + " OR\n"
                   + "   (1) extract PLINK genotype files from MACH .mlprob file (i.e. -extractGenotypes (not the default))\n"
                   + "   (2) pedigree file (i.e. ped=" + pedfile + " (default))\n"
                   + "   (3) MACH file name root (i.e. probFormat=" + probFormat + " (default))\n"
                   + "   (4) rsq threshold for marker inclusion (i.e. rsq=" + rsqThreshold
                   + " (default))\n" + "   (5) prob threshold for genotype calling (i.e. prob="
                   + probRequirement + " (default))\n"
                   + "   (6) (optional) list of subset of markers to test (i.e. subset="
                   + markerSubset + " (default))\n" + " OR\n"
                   + "   (1) extract PLINK genotype files for all chromosomes using defaults (i.e. -extractAllGenotypes (not the default))\n"
                   + "   (2) MACH file name root (i.e. mach=" + machRoot + " (default))\n" + " OR\n"
                   + "   (1) pick best call rate from gender/affection crosstab (i.e. -pickBest (not the default))\n"
                   + "   (2) size of each cell (i.e. pickSize=" + pickSize + " (default))\n"
                   + " OR\n"
                   + "   (1) check which chromosomes have been completed (i.e. -check (not the default))\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("root=")) {
        root = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("keep=")) {
        keeps = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("exclude=")) {
        excludes = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("convert=")) {
        convert = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-createFiles")) {
        createfiles = true;
        numArgs--;
      } else if (arg.startsWith("-batches")) {
        batch = true;
        numArgs--;
      } else if (arg.startsWith("nodesToUse=")) {
        nodesToUse = arg.split("=")[1].split(",");
        numArgs--;
      } else if (arg.startsWith("trim=")) {
        trim = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("decode=")) {
        decode = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("list=")) {
        list = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("ped=")) {
        pedfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("mach=")) {
        machRoot = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("subset=")) {
        markerSubset = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-extractDosage")) {
        extractDosage = true;
        numArgs--;
      } else if (arg.startsWith("extract=")) {
        extract = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-intEx")) {
        intEx = true;
        numArgs--;
      } else if (arg.startsWith("splitByRegion=")) {
        splitByRegion = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("doseFormat=")) {
        doseFormat = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("infoFormat=")) {
        infoFormat = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("probFormat=")) {
        probFormat = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("database=")) {
        database = arg.split("=")[1];
        if (database.endsWith(".mldose")) {
          database = database.substring(0, database.lastIndexOf("."));
        }
        numArgs--;
      } else if (arg.startsWith("-extractGenotypes")) {
        extractGenotypes = true;
        numArgs--;
      } else if (arg.startsWith("-extractAllGenotypes")) {
        extractAllGenotypes = true;
        numArgs--;
      } else if (arg.startsWith("rsq=")) {
        rsqThreshold = Double.parseDouble(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("prob=")) {
        probRequirement = Double.parseDouble(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("-pickBest")) {
        pickBest = true;
        numArgs--;
      } else if (arg.startsWith("pickSize=")) {
        pickSize = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("-check")) {
        checkResults = true;
        numArgs--;
      } else {
        System.err.println("Error - don't know what to do with: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    // database = "aric_overlap";

    try {
      if (!convert.equals("")) {
        convertMapToDat(convert);
      } else if (list != null) {
        listIndividualsInMldose(list, "list.txt");
        // long time;
        // time = new Date().getTime();
        // listIndividualsInMldoseSlowly(list, "listS1.txt");
        // System.out.println("S1 in " + ext.getTimeElapsed(time));
        // time = new Date().getTime();
        // listIndividualsInMldose(list, "listF1.txt");
        // System.out.println("F1 in " + ext.getTimeElapsed(time));
        // time = new Date().getTime();
        // listIndividualsInMldose(list, "listF2.txt");
        // System.out.println("F2 in " + ext.getTimeElapsed(time));
        // time = new Date().getTime();
        // listIndividualsInMldoseSlowly(list, "listS2.txt");
        // System.out.println("S2 in " + ext.getTimeElapsed(time));
      } else if (createfiles) {
        // batchFileCreation(root, keeps, excludes, DEFAULT_SUBSET, DEFAULT_ALL);
        batchFileCreation(root, keeps, excludes, ext.rootOf(keeps), DEFAULT_ALL);
      } else if (batch) {
        steps(ext.rootOf(keeps), DEFAULT_ALL, nodesToUse);
        // } else if (batchAllInOne) {
        // System.err.println("Using defaults:");
        // System.err.println(" root="+root);
        // System.err.println(" keep="+keeps);
        // System.err.println(" exclude="+excludes);
        // batchAllInOne(root, keeps, excludes, "subset", "all", trim > 0);
      } else if (trim > 0) {
        trimReference(DEFAULT_ALL, trim);
      } else if (decode > 0) {
        decodePhasedHapMap(dir, decode);
      } else if (pickBest) {
        pickBest(pickSize);
      } else if (extractDosage) {
        extractDosageForPlink(dir, pedfile, doseFormat, infoFormat, markerSubset);
      } else if (!extract.equals("")) {
        if (intEx || splitByRegion >= 0) {
          extractIntermediate(extract, splitByRegion, doseFormat, infoFormat);
        } else {
          extractSpecificMarkers("", extract, doseFormat, infoFormat, true, new Logger());
        }
      } else if (!database.equals("")) {
        createDatabase(database);
      } else if (extractGenotypes) {
        extractGenotypesForPlink(dir, pedfile, probFormat, infoFormat, markerSubset, rsqThreshold,
                                 probRequirement);
      } else if (extractAllGenotypes) {
        parseGenotypesTogether(machRoot);
      } else if (checkResults) {
        checkResults();
      } else {
        System.err.println("Error - no option was selected...");
        System.err.println(usage);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
