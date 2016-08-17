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
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.IntVector;
import org.genvisis.common.ext;
import org.genvisis.filesys.SnpMarkerSet;

public class PlinkToLinkage {
  public static void convert(String dir, String root) {
    for (int chr = 1; chr <= 23; chr++) {
      new File(dir + "chrom" + ext.chrome(chr) + ".pre").delete();
      new File(dir + "map" + ext.chrome(chr) + ".dat").delete();
      CmdLine.run("plink --bfile " + root + " --chr " + chr + " --recode12 --out chrom"
                  + ext.chrome(chr), dir);
      if (new File(dir + "chrom" + ext.chrome(chr) + ".map").exists()) {
        new File(dir + "chrom" + ext.chrome(chr)
                 + ".ped").renameTo(new File(dir + "chrom" + ext.chrome(chr) + ".pre"));
        new SnpMarkerSet(dir + "chrom" + ext.chrome(chr) + ".map").createLinkageMap()
                                                                  .createFile(dir + "map"
                                                                              + ext.chrome(chr) + ".dat");
        new File(dir + "chrom" + ext.chrome(chr) + ".ped").delete();
        new File(dir + "chrom" + ext.chrome(chr) + ".map").delete();
      }
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    // String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\";
    // String root = "slim";
    // int slim = 2000;
    String dir = "";
    String root = "slimmedDown";
    int slim = 0;
    boolean create = true;
    String pedigree = "";
    // String pedigree = "fullPedigree.dat";
    // String pedigree = "FullPedigreeWithoutDuplicates.dat";

    String usage =
        "\n" + "link.PlinkToLinkage requires 0-1 arguments\n" + "   (1) directory (i.e. dir=" + dir
                   + " (default))\n" + "   (2) plink root name (i.e. root=" + root + " (default))\n"
                   + "   (3) number of markers to slim down to (i.e. slim=" + slim + " (default))\n"
                   + "   (4) make LinkageFormat files (i.e. -create (" + (create ? "" : "not the ")
                   + " default))\n" + "   (5) pedigree file to update to (i.e. ped=" + pedigree
                   + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("root=")) {
        root = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("slim=")) {
        slim = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("-create")) {
        create = true;
        numArgs--;
      } else if (arg.startsWith("ped=")) {
        pedigree = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      if (slim > 0) {
        root = slim(dir, root, slim);
      }
      if (create) {
        convert(dir, root);
      }
      if (!pedigree.equals("")) {
        updatePedigreeStructure(dir, pedigree);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static String slim(String dir, String root, int slim) {
    PrintWriter writer;
    SnpMarkerSet bim;
    String[] markerNames;
    IntVector xs, autos;
    byte[] chrs;
    int numPres, numWanted, numToDel;

    bim = new SnpMarkerSet(dir + root + ".bim");
    markerNames = bim.getMarkerNames();
    chrs = bim.getChrs();
    xs = new IntVector();
    autos = new IntVector();
    for (int i = 0; i < markerNames.length; i++) {
      if (chrs[i] > 0 && chrs[i] < 23) {
        autos.add(i);
      } else if (chrs[i] == 23) {
        xs.add(i);
      }
    }

    numPres = xs.size();
    numWanted = numPres > slim / 10 ? slim / 10 : numPres;
    numToDel = numPres - numWanted;
    for (int i = numToDel - 1; i >= 0; i--) {
      xs.removeElementAt((int) ((double) i * (double) numPres / numToDel));
    }

    numPres = autos.size();
    numWanted = slim - numWanted;
    numWanted = numPres > numWanted ? numWanted : numPres;
    numToDel = numPres - numWanted;
    for (int i = numToDel - 1; i >= 0; i--) {
      autos.removeElementAt((int) ((double) i * (double) numPres / numToDel));
    }

    try {
      writer = new PrintWriter(new FileWriter(dir + "slimmedList.txt"));
      for (int i = 0; i < autos.size(); i++) {
        writer.println(markerNames[autos.elementAt(i)]);
      }
      for (int i = 0; i < xs.size(); i++) {
        writer.println(markerNames[xs.elementAt(i)]);
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing list");
      e.printStackTrace();
    }

    CmdLine.run("plink --bfile " + root + " --extract slimmedList.txt --make-bed --out slimmedDown",
                dir);

    return "slimmedDown";
  }

  public static void updatePedigreeStructure(String dir, String pedigree) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    Hashtable<String, String[]> hash = new Hashtable<String, String[]>();
    int count;
    String backup;

    for (int chr = 1; chr <= 23; chr++) {
      count = 0;
      try {
        reader = new BufferedReader(new FileReader(dir + "chrom" + ext.chrome(chr) + ".pre"));
        while (reader.ready()) {
          line = reader.readLine().trim().split("[\\s]+");
          if (count == 0) {
            count = line.length - 6;
          } else if (line.length - 6 != count) {
            System.err.println("Error - mismatched line length for " + line[0] + "-" + line[1]);
          }
          hash.put(line[0] + "\t" + line[1], line);
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + dir + "chrom" + ext.chrome(chr) + ".pre"
                           + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + dir + "chrom" + ext.chrome(chr) + ".pre"
                           + "\"");
        System.exit(2);
      }

      try {
        reader = new BufferedReader(new FileReader(dir + pedigree));
        backup = Files.backup("chrom" + ext.chrome(chr) + ".pre", dir, dir, true);
        writer = new PrintWriter(new FileWriter(dir + "chrom" + ext.chrome(chr) + ".pre"));
        while (reader.ready()) {
          line = reader.readLine().trim().split("[\\s]+");
          writer.print(Array.toStr(Array.subArray(line, 0, 6)));
          if (hash.containsKey(line[0] + "\t" + line[1])) {
            line = hash.get(line[0] + "\t" + line[1]);
            for (int i = 6; i < line.length; i++) {
              writer.print("\t" + line[i]);
            }
          } else {
            for (int i = 0; i < count; i++) {
              writer.print("\t0");
            }
          }
          writer.println();
        }
        reader.close();
        writer.close();
        new File(dir + backup).delete();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + dir + pedigree + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + dir + pedigree + "\"");
        System.exit(2);
      }
    }

  }
}
