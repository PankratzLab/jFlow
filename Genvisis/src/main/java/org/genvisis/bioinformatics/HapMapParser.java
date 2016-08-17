// expecting a HapMart/BioMart export that includes position, marker_id, reference_allele and of
// course CEU genotype
// pedinfo2sample_CEU.txt is available at http://www.hapmap.org/downloads/samples_individuals/
package org.genvisis.bioinformatics;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Aliases;
import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.FamilyStructure;
import org.genvisis.filesys.SnpMarkerSet;

public class HapMapParser {
  public static final boolean USE_NUCLEOTIDES = true;
  public static final boolean INCLUDE_MONOMORPHIC = true;

  public static final String[][] TARGETS_WITH_ALTS =
      {{"genotype"}, {"position"}, {"reference allele"}, {"marker id"}, {"chromosome"}};

  public static void fixAffectionStatusInBed(String root) {
    FamilyStructure struct;
    long timestamp;

    for (int chr = 1; chr <= 23; chr++) {
      timestamp = new File(root + ".chr" + chr + ".fam").lastModified();
      struct = new FamilyStructure(root + ".chr" + chr + ".fam");
      struct.setAffections(Array.byteArray(struct.getAffections().length, (byte) 1));
      struct.writeToFile(root + ".chr" + chr + ".fam", false);
      new File(root + ".chr" + chr + ".fam").setLastModified(timestamp);
      System.out.print(".");
    }
    System.out.println();
  }

  public static void generateHaploviewBatch(String dir, String root, boolean preNotPed,
      Logger log) {
    PrintWriter writer;

    new SnpMarkerSet(dir + root + ".map", true, log).writeToFile(dir + root + ".info",
        SnpMarkerSet.HAPLOVIEW_INFO_FORMAT, log);
    try {
      writer = new PrintWriter(new FileWriter((new File(dir).exists() ? dir : "") + root + ".bat"));
      writer.println("java -jar /home/npankrat/Haploview.jar -pedfile " + root + "."
          + (preNotPed ? "pre" : "ped") + " -info " + root + ".info");
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing batch file");
      log.reportException(e);
    }
  }

  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    String dir = "";
    String filename = null;
    String map = "";
    String famstruct;

    // famstruct = "famstruct_CEU.dat";
    // famstruct = "pedinfo2sample_CEU.txt";
    // famstruct = "fakeIDs_for_CEU.txt";
    famstruct = "CEU_structUnrelated.dat";
    famstruct = Aliases.getPathToFileInReferenceDirectory(famstruct, false, new Logger());

    String bed = "";
    String ped = "";
    String fix = "";

    String usage = "\n" + "bioinformatics.HapMapParser requires 0-1 arguments\n"
        + "   (1) filename (i.e. file=MAPT.tsv (not the default)\n"
        + "   (2) famstruct (i.e. struct=" + famstruct + " (default)\n"
        + "   (3) split bed by chromosome (i.e. bed=plink (not the default)\n"
        + "   (4) split ped by chromosome (i.e. ped=plink (not the default)\n"
        + "   (5) fix affection status in fam files (i.e. fix=plink (not the default)\n" + "  OR:\n"
        + "   (1) generate Haploview batch file from PLINK files (i.e. map=plink.map (not the default)\n"
        + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("struct=")) {
        famstruct = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("bed=")) {
        bed = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("ped=")) {
        ped = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("fix=")) {
        fix = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("map=")) {
        map = arg.split("=")[1];
        numArgs--;
      }
    }

    if (filename == null && map.equals("")) {
      System.err.println("Error - need to pass a filename as an argument (i.e. file=MAPT.tsv)");
      ext.waitForResponse();
      return;
    }

    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      if (filename == null) {
        filename = map;
      }
      if (!map.equals("")) {
        generateHaploviewBatch("", ext.removeDirectoryInfo(map.substring(0, map.lastIndexOf("."))),
            false, new Logger(
                ext.rootOf(filename, false) + "_haploview_prep.log"));
      } else if (!bed.equals("")) {
        splitBedByChromosome(bed);
      } else if (!ped.equals("")) {
        splitPedByChromosome(ped);
      } else if (!fix.equals("")) {
        fixAffectionStatusInBed(fix);
      } else {
        parse(dir, filename, famstruct, new Logger(
            ext.rootOf(filename, false) + "_hapmap_parser.log"));
      }
    } catch (Exception e) {
      e.printStackTrace();
      ext.waitForResponse();
    }
  }

  public static void parse(String dir, String filename, String famstruct, Logger log) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, trans;
    String temp, trav;
    Hashtable<String, String[]> hash = new Hashtable<String, String[]>();
    Vector<String> v = new Vector<String>();
    int index, ones, twos;
    String root = ext.rootOf(filename);
    int[] indices;
    String refAllele, position, chr, markerName;
    String[] indIDs = null;

    try {
      reader = Files.getReader(filename, dir);
      temp = reader.readLine();
      indices = ext.indexFactors(TARGETS_WITH_ALTS, temp.trim().split("\\t", -1), false, false,
          true, log, true);

      writer = new PrintWriter(new FileWriter((new File(dir).exists() ? dir : "") + root + ".map"));
      trans = null;
      while (reader.ready()) {
        line = reader.readLine().split("\\t", -1);
        refAllele = line[indices[2]];
        position = line[indices[1]];
        chr = line[indices[4]].substring(3);
        markerName = line[indices[3]];

        line = line[indices[0]].split("[\\s]+");
        trans = new String[line.length];
        ones = twos = 0;
        for (int i = 0; i < line.length; i++) {
          trans[i] = "";
          for (int j = 0; j < 2; j++) {
            if (line[i].charAt(j) == 'N') {
              trans[i] += "0";
            } else if (line[i].charAt(j) == refAllele.charAt(0)) {
              if (USE_NUCLEOTIDES) {
                trans[i] += line[i].charAt(j);
              } else {
                trans[i] += "2";
              }
              twos++;
            } else {
              if (USE_NUCLEOTIDES) {
                trans[i] += line[i].charAt(j);
              } else {
                trans[i] += "1";
              }
              ones++;
            }
            trans[i] += j == 0 ? "\t" : "";
          }
        }
        if ((ones > 0 && twos > 0) || INCLUDE_MONOMORPHIC) {
          writer.println(chr + "\t" + markerName + "\t0\t" + position);
          v.add(markerName);
          hash.put(markerName, trans);
        }
      }
      reader.close();
      writer.close();
      indIDs = Array.stringArraySequence(trans.length, "");
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + filename + "\"");
      log.reportException(ioe);
      return;
    }

    try {
      reader = Files.getAppropriateReader(famstruct);
      writer = new PrintWriter(new FileWriter((new File(dir).exists() ? dir : "") + root + ".pre"));
      while (reader.ready()) {
        line = reader.readLine().split("[\\s]+");
        for (int i = 0; i < 5; i++) {
          writer.print(line[i] + "\t");
        }
        writer.print("1");
        index = -2;
        trav = line[6];
        if (famstruct.contains("pedinfo2sample_")) {
          line = line[6].split(":");
          if (line.length != 6) {
            log.reportError(
                "Error - different format than expected for pedinfo2sample_***.txt file (do not alter from what's posted)");
            reader.close();
            writer.close();
            return;
          }
          trav = line[4];
        }
        for (int i = 0; i < indIDs.length; i++) {
          if (indIDs[i].equals(trav)) {
            index = i;
          }
        }
        if (index == -2) {
          log.reportError("Error - Could not find sample " + trav
              + " from the pedigree file in the genotype file");
        }
        for (int i = 0; i < v.size(); i++) {
          line = hash.get(v.elementAt(i));
          writer.print("\t" + line[index]);
        }
        writer.println();
      }
      reader.close();
      writer.close();
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + famstruct + "\"");
      log.reportException(ioe);
      return;
    }

    generateHaploviewBatch(dir, root, true, log);
  }

  public static void plinkMapToHaploviewInfo(String from, String to, Logger log) {
    new SnpMarkerSet(from).writeToFile(to, SnpMarkerSet.HAPLOVIEW_INFO_FORMAT, log);
  }

  public static void splitBedByChromosome(String root) {
    for (int chr = 1; chr <= 23; chr++) {
      System.out.print(".");
      CmdLine.run("plink --noweb --bfile " + root + " --chr " + chr + " --make-bed --out " + root
          + ".chr" + chr, "./");
    }
    System.out.println();
  }

  public static void splitPedByChromosome(String root) {
    for (int chr = 1; chr <= 23; chr++) {
      System.out.print(".");
      CmdLine.run("plink --noweb --file " + root + " --chr " + chr + " --recode --out " + root
          + ".chr" + chr, "./");
    }
    System.out.println();
  }
}
