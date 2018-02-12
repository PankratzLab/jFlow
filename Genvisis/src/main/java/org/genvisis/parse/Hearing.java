package org.genvisis.parse;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.DoubleVector;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.link.LinkageMap;
import com.google.common.primitives.Doubles;

public class Hearing {

  public static boolean AFFECTEDS_ONLY = true;
  public static final String DBSNP_LOCAL = "local_6K_b129.bcp";
  public static final String RUTGERS_GENETIC_MAP = "6K_marker.map.txt";
  public static final String[][] RENAMES = {{"rs860732", "rs449511"}, {"rs1136206", "rs9788"}};

  public static void parse(String dir, String source) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String trav;
    Hashtable<String, String> hash = new Hashtable<String, String>();
    Vector<String> names;
    DoubleVector positions;

    names = new Vector<String>();
    for (int chr = 1; chr <= 22; chr++) {
      try {
        reader = new BufferedReader(new FileReader(dir + source + "chr" + chr + "_SNP.noLD.ped"));
        writer = Files.openAppropriateWriter(dir + "mrkr" + ext.chrome(chr) + ".dat");
        while (reader.ready()) {
          line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
          if (chr == 1) {
            names.add(ArrayUtils.toStr(ArrayUtils.subArray(line, 0, 6)) + "\t" + line[1]
                      + "\t.\t.");
          }
          writer.print(line[0] + "\t" + line[1]);
          for (int i = 6; i < line.length; i++) {
            if (line[i].length() != 3 || line[i].charAt(1) != '/') {
              System.err.println("Error - invalid genotype: '" + line[i] + "'");

            }
            writer.print("\t" + line[i].charAt(0) + "\t" + line[i].charAt(2));
          }
          writer.println();
        }
        reader.close();
        writer.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + dir + source + "chr" + chr + "_SNP.noLD.ped"
                           + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + dir + source + "chr" + chr + "_SNP.noLD.ped"
                           + "\"");
        System.exit(2);
      }
    }

    hash = HashVec.loadFileToHashString(dir + source + "pheno.dat", false);
    try {
      writer = Files.openAppropriateWriter(dir + "example_struct.dat");
      for (int i = 0; i < names.size(); i++) {
        line = names.elementAt(i).split(PSF.Regex.GREEDY_WHITESPACE);
        trav = hash.get(line[1]);
        if (trav == null) {
          System.err.println("No phenotype for individual '" + line[1] + "'");
          line[5] = "0";
        } else {
          line[5] = trav.equals("A") ? "2"
                                     : (trav.equals("U") ? (AFFECTEDS_ONLY ? "0" : "1") : "0");
        }
        writer.println(ArrayUtils.toStr(line));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to example_struct.dat");
      e.printStackTrace();
    }

    for (int chr = 1; chr <= 22; chr++) {
      names = new Vector<String>();
      positions = new DoubleVector();
      try {
        reader = new BufferedReader(new FileReader(dir + source + "chr" + chr + "_SNP.noLD.map"));
        while (reader.ready()) {
          line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
          if (line.length != 3) {
            System.err.println("Error - what am I supposed to do with '" + ArrayUtils.toStr(line)
                               + "'??");
          }
          if (Integer.parseInt(line[0]) != chr) {
            System.err.println("Error - a marker for chromosome " + line[0] + " was found in "
                               + source + "chr" + chr + "_SNP.noLD.map");
          }
          for (String[] element : RENAMES) {
            if (line[1].equalsIgnoreCase(element[0])) {
              line[1] = element[1];
            }
          }
          names.add(line[1]);
          positions.add(Double.parseDouble(line[2]));
        }
        reader.close();

        new LinkageMap(chr, ArrayUtils.toStringArray(names), 2, Doubles.toArray(positions), false,
                       true).createFile(dir + "map" + ext.chrome(chr) + ".dat");
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + dir + source + "chr" + chr + "_SNP.noLD.ped"
                           + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + dir + source + "chr" + chr + "_SNP.noLD.ped"
                           + "\"");
        System.exit(2);
      }
    }
  }

  public static void testBits(String dir, String runtime, String filename) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;

    try {
      new File(dir + runtime).mkdirs();
      reader = new BufferedReader(new FileReader(dir + filename));
      writer = Files.openAppropriateWriter(dir + runtime + "re_chrom01.pre");
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        writer.println(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4]
                       + "\t" + line[5] + "\t" + (int) (Math.random() * 2 + 1) + "\t"
                       + (int) (Math.random() * 2 + 1));
      }
      reader.close();
      writer.close();

      new LinkageMap(1, new String[] {"dummy"}, 2, new double[] {0}, false,
                     false).createFile(dir + runtime + "map01.dat");

      writer = Files.openAppropriateWriter(dir + runtime + "useful.opt");
      writer.println("% Read input in LINKAGE style format:");
      writer.println("	PREFILE re_chrom01.pre");
      writer.println("	DATFILE map01.dat");
      writer.println("");
      writer.println("	% Run multipoint analysis on all [equally weighted] pairs and output to  .txt file");
      writer.println("	MODEL mpt lin all equal chrom01.lin.out chromf01.lin.out");
      writer.println("");
      writer.println("% Other options:");
      writer.println("MAXMEMORY 411");
      writer.close();

    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + filename + "\" not found in current directory");
      fnfe.printStackTrace();
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + filename + "\"");
      ioe.printStackTrace();
      System.exit(2);
    }
  }

  //
  // public static void combineSheets() {
  // PrintWriter writer;
  //
  // try {
  // writer = Files.openAppropriateWriter("C:\\Documents and
  // Settings\\npankrat\\My Documents\\ext\\Cuts\\Eight.ahk");
  // for (int chr = 10; chr >= 1; chr--) {
  // writer.println("Sleep, 500");
  // writer.println("Send, {ENTER}");
  // writer.println("Sleep, 2000");
  // writer.println("MouseClick, right, 164, 974");
  // writer.println("Sleep, 500");
  // writer.println("MouseClick, left, 98, 85");
  // writer.println("Sleep, 500");
  // writer.println("MouseClick, left, 105, 78");
  // writer.println("Sleep, 500");
  // writer.println("MouseClick, left, 92, 109");
  // writer.println("Sleep, 500");
  // writer.println("MouseClick, left, 153, 260");
  // writer.println("Sleep, 500");
  // writer.println("MouseClick, left, 132, 968");
  // writer.println("MouseClick, left, 132, 968");
  // writer.println("Sleep, 500");
  // writer.println("Send, chr"+chr+"{ENTER}");
  // writer.println("Sleep, 500");
  // writer.println("Send, {ALTDOWN}{TAB}{ALTUP}");
  // writer.println("Sleep, 500");
  // writer.println("Send, {UP}");
  // }
  // writer.close();
  // } catch (Exception e) {
  // System.err.println("Error writing script");
  // e.printStackTrace();
  // }
  // }

  public static void main(String[] args) {
    int numArgs = args.length;
    // String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\hearing\\";
    // String dir = "D:/BOSS/LinkageMergedIBC/";
    String dir = "D:/BOSS/LinkageJustIBC/";
    String source = "00src\\raw_genotypes\\";
    boolean parse = false;
    boolean test = true;
    boolean update = false;
    int markerIndex = 2;
    int chrIndex = 1;
    int cmIndex = 4;
    String db = RUTGERS_GENETIC_MAP;

    String usage = "\\n" + "parse.Hearing requires 0-1 arguments\n" + "   (1) directory (i.e. dir="
                   + dir + " (default))\n" + "   (2) source (i.e. src=" + source + " (default))\n"
                   + "   (3) update genetic positions (i.e. -update (" + (update ? "" : "not the ")
                   + "default))\n" + "   (4) database of snp cM positions (i.e. db=" + db
                   + " (default for updating))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("src=")) {
        source = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-parse")) {
        parse = true;
        numArgs--;
      } else if (arg.startsWith("-update")) {
        parse = true;
        numArgs--;
      } else if (arg.startsWith("db=")) {
        db = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      if (update) {
        LinkageMap.updateMaps(dir, db, markerIndex, chrIndex, cmIndex);
      }
      if (parse) {
        parse(dir, source);
      }
      if (test) {
        testBits(dir, "test/", "struct.dat");
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
