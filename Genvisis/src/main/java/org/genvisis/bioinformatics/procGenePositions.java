package org.genvisis.bioinformatics;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

import org.genvisis.common.Aliases;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

public class procGenePositions {
  public static class GenePosition implements Cloneable {
    public static GenePosition[] newArray(int size) {
      GenePosition[] gps = new GenePosition[size];

      for (int i = 0; i < size; i++) {
        gps[i] = new GenePosition();
      }

      return gps;
    }

    public String name;
    public String chr;
    public String start;
    public String stop;
    public String sense;
    public boolean placed;
    public boolean imputed;

    public boolean originally_unplaced;

    public GenePosition() {
      name = chr = start = stop = sense = ".";
      placed = imputed = originally_unplaced = false;
    }

    @Override
    public Object clone() {
      GenePosition gp;

      gp = new GenePosition();
      gp.name = name;
      gp.chr = chr;
      gp.start = start;
      gp.stop = stop;
      gp.sense = sense;
      gp.placed = placed;
      gp.imputed = imputed;
      gp.originally_unplaced = originally_unplaced;

      return gp;
    }
  }

  // public static final String[] SOURCES = {"reference", "Celera"};
  // public static final boolean[] IMPORTANT = {true, true};
  // public static final String[] SOURCES = {"reference", "Celera", "c22_H2", "c5_H2", "c6_COX",
  // "c6_QBL", "DR53", "CRA_TCAGchr7v2"};
  public static final String[] SOURCES = {"GRCh37.p5-Primary Assembly", "HuRef-Primary Assembly"};

  // public static final boolean[] IMPORTANT = {true, true, false, false, false, false, false,
  // false};

  public static final boolean[] IMPORTANT = {true, true};

  // public static final boolean[] IMPORTANT = {true, true, false, false,
  // false, false, false, false};
  public static final int DOMINANT = 0;


  public static final String GROUP_LABEL_TARGET = "GRCh37.p5-Primary Assembly";

  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    String filename = "seq_gene.md";
    String dir = Aliases.getPathToFileInReferenceDirectory(filename, true, new Logger());

    // latest found in:
    // ftp://ftp.ncbi.nih.gov/genomes/MapView/Homo_sapiens/sequence/BUILD.37.3/initial_release/

    String usage = "\n" + "park.procGenePositions requires 0-1 arguments\n"
        + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      dir = ext.parseDirectoryOfFile(filename);
      filename = ext.removeDirectoryInfo(filename);
      runProcGenePositions(dir, filename);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void runProcGenePositions(String dir, String filename) {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] line, geneIDs;
    String geneId;
    String temp;
    Hashtable<String, GenePosition[]> hash = new Hashtable<String, GenePosition[]>();

    int count;
    GenePosition[] gps, trav;
    GenePosition[][] data;
    double[][] poslarlar;
    double[] poslar;
    int[][] orderler;
    int[] order;
    int maxDigits = 0;
    int low, high;

    try {
      reader = new BufferedReader(new FileReader(dir + filename));
      writer = new PrintWriter(new FileWriter("non-gene genes.out"));
      writer.println(reader.readLine());
      while (reader.ready()) {
        temp = reader.readLine();
        line = temp.split("\t");
        if (line[11].equalsIgnoreCase("GENE")) {
          geneId = line[10].substring(7);
          if (hash.containsKey(geneId)) {
            gps = hash.get(geneId);
          } else {
            hash.put(geneId, gps = GenePosition.newArray(SOURCES.length));
          }

          count = 0;
          for (int i = 0; i < SOURCES.length; i++) {
            if (line[12].equalsIgnoreCase(SOURCES[i])) {
              gps[i].name = line[9];
              if (line[1].indexOf("|") > 0) {
                gps[i].chr = line[1].substring(0, line[1].indexOf("|"));
                gps[i].placed = false;
              } else {
                gps[i].chr = line[1];
                gps[i].placed = true;
              }
              gps[i].start = line[2];
              if (line[2].length() > maxDigits) {
                maxDigits = line[2].length();
              }
              gps[i].stop = line[3];
              gps[i].sense = line[4];
              count++;
            }
          }
          if (count == 0) {
            writer.println(temp);
          }
        }
      }
      reader.close();
      writer.close();

      System.err.println("1");

      geneIDs = HashVec.getKeys(hash);
      for (String geneID2 : geneIDs) {
        trav = hash.get(geneID2);
        temp = "";
        count = 0;
        for (int j = 0; j < SOURCES.length; j++) {
          if (trav[j].chr.equals("Un")) {
            trav[j].chr = ".";
            trav[j].originally_unplaced = true;
          }
          if (temp.equals("") && !trav[j].chr.equals(".")) {
            temp = trav[j].chr;
          } else if (!trav[j].chr.equals(temp) && !trav[j].chr.equals(".")) {
            count++;
          }
        }
        if (count > 0) {
          for (int j = 0; j < SOURCES.length; j++) {
            if (!trav[j].chr.equals(".")) {
              gps = GenePosition.newArray(SOURCES.length);
              gps[j] = (GenePosition) trav[j].clone();
              gps[j].name = gps[j].name + "_asPlacedBy_"
                  + ext.replaceWithLinuxSafeCharacters(SOURCES[j], true);
              hash.put(
                  geneID2 + "_asPlacedBy_" + ext.replaceWithLinuxSafeCharacters(SOURCES[j], true),
                  gps);
            }
          }
          hash.remove(geneID2);
        }
      }

      System.err.println("2");

      geneIDs = HashVec.getKeys(hash);
      data = new GenePosition[geneIDs.length][];
      for (int i = 0; i < geneIDs.length; i++) {
        data[i] = hash.get(geneIDs[i]);
      }
      hash.clear();

      // reader = new BufferedReader(new FileReader("Copy (2) of
      // genes.xls"));
      // reader.readLine();
      // while (reader.ready()) {
      // line = reader.readLine().split("[\\s]+");
      // gps = new GenePosition[2];
      // for (int i = 0; i<SOURCES.length; i++) {
      // gps[i] = new GenePosition();
      // gps[i].geneid = line[1+i*7+0];
      // gps[i].chr = line[1+i*7+1];
      // gps[i].start = line[1+i*7+2];
      // gps[i].stop = line[1+i*7+3];
      // gps[i].sense = line[1+i*7+4];
      // gps[i].placed = line[1+i*7+5].equalsIgnoreCase("true");
      // gps[i].imputed = line[1+i*7+6].equalsIgnoreCase("true");
      // }
      // hash.put(line[0], gps);
      // }
      // reader.close();

      // reader = new BufferedReader(new FileReader("Copy (2) of
      // genes.xls"));
      // reader.readLine();
      // genes = new String[33329];
      // data = new GenePosition[genes.length][];
      // for (int j = 0; j<genes.length; j++) {
      // line = reader.readLine().split("[\\s]+");
      // genes[j] = line[0];
      // data[j] = new GenePosition[2];
      // for (int i = 0; i<SOURCES.length; i++) {
      // data[j][i] = new GenePosition();
      // data[j][i].geneid = line[1+i*7+0];
      // data[j][i].chr = line[1+i*7+1];
      // data[j][i].start = line[1+i*7+2];
      // data[j][i].stop = line[1+i*7+3];
      // data[j][i].sense = line[1+i*7+4];
      // data[j][i].placed = line[1+i*7+5].equalsIgnoreCase("true");
      // data[j][i].imputed = line[1+i*7+6].equalsIgnoreCase("true");
      // }
      // }
      // reader.close();

      System.err.println("3");

      poslarlar = new double[SOURCES.length][geneIDs.length];
      for (int i = 0; i < geneIDs.length; i++) {
        for (int j = 0; j < SOURCES.length; j++) {
          try {
            temp = "";
            if (data[i][j].chr.equals("X")) {
              temp = "23";
            } else if (data[i][j].chr.equals("Y")) {
              temp = "24";
            } else if (data[i][j].chr.equals("MT")) {
              temp = "25";
            } else if (data[i][j].chr.equals("Un")) {
              temp = "26";
            } else if (data[i][j].chr.equals(".")) {
              temp = "30";
            } else {
              temp = data[i][j].chr;
            }
            poslarlar[j][i] = Double.parseDouble(temp + "."
                + ext.formNum(data[i][j].start.equals(".") ? "0" : data[i][j].start, maxDigits));
          } catch (Exception e) {
            poslarlar[j][i] = -1;
          }
        }
      }

      System.err.println("4");

      orderler = new int[SOURCES.length][];
      for (int i = 0; i < SOURCES.length; i++) {
        orderler[i] = Sort.quicksort(poslarlar[i]);
        poslarlar[i] = sortArrayLikeThis(poslarlar[i], orderler[i]);
      }

      System.err.println("5");

      for (int i = 0; i < geneIDs.length; i++) {
        for (int j = 0; j < SOURCES.length; j++) {
          if (IMPORTANT[j] && data[i][j].chr.equals(".")) {
            for (int k = 0; k < SOURCES.length; k++) {
              if (j != k && data[i][j].chr.equals(".") && !data[i][j].originally_unplaced
                  && data[i][k].originally_unplaced) {
                data[i][j].name = data[i][k].name;
                data[i][j].chr = data[i][k].chr;
                data[i][j].start = data[i][k].start;
                data[i][j].stop = data[i][k].stop;
                data[i][j].sense = data[i][k].sense;
                data[i][j].placed = false;
                data[i][j].imputed = true;
              } else if (j != k && data[i][j].chr.equals(".") && !data[i][k].chr.equals(".")
                  && !data[i][k].originally_unplaced) {
                low = high = ext.indexOfInt(i, orderler[k]);
                do {
                  low--;
                } while (low >= 0 && data[orderler[k][low]][k].chr.equals(data[i][k].chr)
                    && data[orderler[k][low]][j].chr.equals("."));
                if (low < 0 || !data[orderler[k][low]][k].chr.equals(data[i][k].chr)) {
                  low = -1;
                }
                do {
                  high++;
                } while (high < geneIDs.length
                    && data[orderler[k][high]][k].chr.equals(data[i][k].chr)
                    && data[orderler[k][high]][j].chr.equals("."));
                if (high == geneIDs.length
                    || !data[orderler[k][high]][k].chr.equals(data[i][k].chr)) {
                  high = -1;
                }
                if (low != -1 || high != -1) {
                  data[i][j].name = data[i][k].name;
                  data[i][j].chr = data[i][k].chr;
                  if ((low != -1 && high == -1) || (low != -1 && high != -1
                      && poslarlar[k][i]
                          - poslarlar[k][orderler[k][low]] < poslarlar[k][orderler[k][high]]
                              - poslarlar[k][i])) {
                    data[i][j].start = (Integer.parseInt(data[orderler[k][low]][k].start)
                        + Integer.parseInt(data[i][k].start)
                        - Integer.parseInt(data[orderler[k][low]][k].start)) + "";
                    data[i][j].stop =
                        (Integer.parseInt(data[i][j].start) + Integer.parseInt(data[i][k].stop)
                            - Integer.parseInt(data[i][k].start)) + "";
                  } else {
                    data[i][j].start = (Integer.parseInt(data[orderler[k][high]][k].start)
                        + Integer.parseInt(data[i][k].start)
                        - Integer.parseInt(data[orderler[k][high]][k].start)) + "";
                    data[i][j].stop =
                        (Integer.parseInt(data[i][j].start) + Integer.parseInt(data[i][k].stop)
                            - Integer.parseInt(data[i][k].start)) + "";
                  }
                  data[i][j].sense = data[i][k].sense;
                  data[i][j].placed = false;
                  data[i][j].imputed = true;
                }
              }
            }
          }
        }
      }

      System.err.println("6");

      poslar = new double[geneIDs.length];
      for (int i = 0; i < geneIDs.length; i++) {
        try {
          temp = "";
          if (data[i][DOMINANT].chr.equals("X")) {
            temp = "23";
          } else if (data[i][DOMINANT].chr.equals("Y")) {
            temp = "24";
          } else if (data[i][DOMINANT].chr.equals("MT")) {
            temp = "25";
          } else if (data[i][DOMINANT].chr.equals("Un")) {
            temp = "26";
          } else if (data[i][DOMINANT].chr.equals(".")) {
            temp = "30";
          } else {
            temp = data[i][DOMINANT].chr;
          }
          poslar[i] = Double.parseDouble(temp + "." + ext.formNum(
              data[i][DOMINANT].start.equals(".") ? "0" : data[i][DOMINANT].start, maxDigits));
        } catch (Exception e) {
          poslar[i] = -1;
        }
      }
      order = Sort.quicksort(poslar);

      System.err.println("7");

      try {
        writer = new PrintWriter(new FileWriter(dir + "genes.xln"));
        writer.print("GeneID");
        for (int i = 0; i < SOURCES.length; i++) {
          if (IMPORTANT[i]) {
            writer.print("\t" + SOURCES[i] + "_name\t" + SOURCES[i] + "_chr\t" + SOURCES[i]
                + "_start\t" + SOURCES[i] + "_stop\t" + SOURCES[i] + "_sense\t" + SOURCES[i]
                + "_placed\t" + SOURCES[i] + "_imputed");
          }
        }
        writer.println();
        for (int i = 0; i < geneIDs.length; i++) {
          writer.print(geneIDs[order[i]].indexOf("_asPlacedBy_") > 0
              ? geneIDs[order[i]].substring(0, geneIDs[order[i]].indexOf("_asPlacedBy_"))
              : geneIDs[order[i]]);
          for (int j = 0; j < SOURCES.length; j++) {
            if (IMPORTANT[j]) {
              writer.print("\t" + data[order[i]][j].name + "\t"
                  + (data[order[i]][j].chr.equals(".") ? "Un" : data[order[i]][j].chr) + "\t"
                  + data[order[i]][j].start + "\t" + data[order[i]][j].stop + "\t"
                  + data[order[i]][j].sense + "\t" + data[order[i]][j].placed + "\t"
                  + data[order[i]][j].imputed);
            }
          }
          writer.println();
        }

        writer.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println(
            "Error: file \"" + dir + "genes.xls" + "\" is still in use or otherwise unavailable");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error writing to \"" + dir + "genes.xls" + "\"");
        System.exit(2);
      }

    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + filename + "\"");
      System.exit(2);
    }

    System.err.println("8");

  }

  public static double[] sortArrayLikeThis(double[] unsorted, int[] order) {
    double[] sorted = new double[unsorted.length];

    for (int i = 0; i < unsorted.length; i++) {
      sorted[i] = unsorted[order[i]];
    }

    return sorted;
  }
}
