package org.genvisis.link;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;

import org.genvisis.common.Files;
import org.genvisis.common.ext;

public class sumOSA {
  public class osaChrom {
    private int maxPos;
    private double baseLOD;
    private double maxLOD;
    private final int totalFams;
    private int maxFams;
    private String pval;
    public double cutoff;
    public int repsCompleted;
    public int repsRequested;

    public osaChrom(int numFams) {
      totalFams = numFams;
    }
  }

  public static int NUM_PTS = 10;

  public static String[] SUGGESTIONS = {"struct111-/", "struct100-/", "mm/struct111-/",
      "gh/struct111-/", "mm/struct100-/", "gh/struct100-/"};

  public static void main(String[] args) {
    if (args.length > 1) {
      System.out.println("Expecting maybe 1 argument: LOD score cutoff [optional].");
    } else {
      try {
        new sumOSA();
      } catch (Exception e) {
        e.printStackTrace();
      }
    }
  }

  public sumOSA() throws IOException {
    Vector<String> runs = new Vector<String>();
    int count = 0;

    for (int i = 1; i <= 23; i++) {
      if (new File("chrom" + (i < 10 ? "0" + i : "" + i)).exists()) {
        count++;
      }
    }
    if (count > 10) {
      runs.add("");
    }

    for (String element : SUGGESTIONS) {
      if (new File(element).exists()) {
        runs.add(element);
      }
    }

    new sumOSA(runs);
  }

  public sumOSA(String[] dirs) throws IOException {
    Vector<String> runs = new Vector<String>();

    for (String dir : dirs) {
      runs.add(dir);
    }

    new sumOSA(runs);
  }

  public sumOSA(Vector<String> runDirs) throws IOException {
    PrintWriter writer = null, batch = null;
    String temp;
    String base;
    osaChrom[][] osac;
    String plug = "plug";
    int chr;
    Vector<osaChrom[][]> runs = new Vector<osaChrom[][]>();

    batch = new PrintWriter(new FileWriter("batchPerms"));
    batch.println("sleep 5");
    batch.println();

    plug = (new File(plug)).getAbsolutePath();
    writer = new PrintWriter(new FileWriter(plug));
    writer.close();

    for (int run = 0; run < runDirs.size(); run++) {
      base = runDirs.elementAt(run);

      try {
        if (new File(base + "AscendingKey.dat").exists()) {
          osac = parseMapmaker(base, batch, plug);
        } else {
          osac = parseHauser(base);
        }
        runs.add(osac);
      } catch (Exception e) {
        System.err.println("Error parsing run in " + base);
        e.printStackTrace();
        System.exit(1);
      }
    }

    temp = "";
    writer = new PrintWriter(new FileWriter("osa_summary.xls"));
    for (int run = 0; run < runDirs.size(); run++) {
      writer.print((run == 0 ? "\t" : "") + "\t" + runDirs.elementAt(run) + "\t\t\t\t\t\t\t\t");
      temp += "\t" + (run == 0 ? "Direction" : "")
          + "\tPosition\tmaxLOD\tbaseLOD\tIncrease\tSignificance\tFamsUsed\tTotalFams\tProportion";
    }
    writer.println();
    writer.println(temp);
    for (int i = 0; i < 23; i++) {
      chr = i + 1;
      for (int k = 0; k < 2; k++) {
        writer.print((k == 0 ? "Chr " + chr + "\t" + "Ascending" : "\t" + "Descending"));
        for (int run = 0; run < runs.size(); run++) {
          osac = runs.elementAt(run);
          if (osac[k][i] != null) {
            writer.print("\t" + osac[k][i].maxPos + "\t" + ext.formDeci(osac[k][i].maxLOD, 2, true)
                + "\t" + ext.formDeci(osac[k][i].baseLOD, 2, true) + "\t"
                + ext.formDeci(osac[k][i].maxLOD - osac[k][i].baseLOD, 2, true) + "\t"
                + osac[k][i].pval + "\t" + osac[k][i].maxFams + "\t" + osac[k][i].totalFams + "\t"
                + ext.formDeci(((double) osac[k][i].maxFams / (double) osac[k][i].totalFams), 3,
                    true)
                + "\t");
          } else {
            writer.print("\t\t\t\t\t\t\t\t\t");
          }
        }
        writer.println();
      }
      writer.println();
    }
    writer.close();
    batch.close();
  }

  public osaChrom[][] parseHauser(String base) throws IOException {
    BufferedReader reader = null;
    String[] line;
    String temp, chrome;
    int chr;
    osaChrom[][] osac = new osaChrom[2][23];

    for (int j = 0; j < 23; j++) {
      chr = j + 1;
      chrome = (chr < 10) ? "0" + chr : "" + chr;
      temp = base + "chrom" + chrome + "/osa" + chrome + ".dat.max";
      if (!new File(temp).exists() || !(new BufferedReader(new FileReader(temp))).ready()) {
        continue;
      }

      reader = new BufferedReader(new FileReader(temp));
      do {
        temp = reader.readLine();
      } while (!temp.startsWith("  Variable"));
      reader.readLine();

      for (int i = 0; i < 2; i++) {
        line = reader.readLine().split("[\\s]+");
        osac[i][j] = new osaChrom(Integer.valueOf(line[10]).intValue());
        osac[i][j].maxPos = (int) Double.valueOf(line[3]).doubleValue();
        osac[i][j].baseLOD = Double.valueOf(line[5]).doubleValue();
        osac[i][j].maxLOD = Double.valueOf(line[4]).doubleValue();
        osac[i][j].maxFams = Integer.valueOf(line[9]).intValue();
        osac[i][j].pval = line[8];
        osac[i][j].repsCompleted = Integer.valueOf(line[7]).intValue();
      }

      reader.close();
    }

    return osac;
  }

  public osaChrom[][] parseMapmaker(String base, PrintWriter batch, String plug)
      throws IOException {
    BufferedReader reader = null;
    String[] line;
    String temp, chrome;
    Vector<String> cutLookup, cutoffLookup;
    double lod;
    int chr, count, total, reps;
    osaChrom[][] osac = new osaChrom[2][23];

    for (int i = 0; i < 2; i++) {
      cutLookup = new Vector<String>();
      cutoffLookup = new Vector<String>();
      temp = base + ((i == 0) ? "Ascending" : "Descending") + "Key.dat";
      if (new File(temp).exists()) {
        reader = new BufferedReader(new FileReader(temp));
        while (reader.ready()) {
          line = reader.readLine().split("[\\s]+");
          cutLookup.add(line[1]);
          cutoffLookup.add(line[2]);
        }
        reader.close();
      } else {
        continue;
      }

      for (int j = 0; j < 23; j++) {
        chr = j + 1;
        chrome = (chr < 10) ? "0" + chr : "" + chr;
        temp = base + "chrom" + chrome + ((i == 0) ? "a" : "d") + "-surface.xls";
        if (!new File(temp).exists() || !(new BufferedReader(new FileReader(temp))).ready()) {
          continue;
        }
        osac[i][j] =
            new osaChrom(Integer.valueOf(cutLookup.elementAt(cutLookup.size() - 1)).intValue());

        reader = new BufferedReader(new FileReader(temp));
        reader.readLine();
        for (int cut = 0; cut < cutLookup.size() - 2; cut++) {
          line = reader.readLine().split("[\\s]+");
          for (int k = 1; k < line.length; k++) {
            lod = Double.valueOf(line[k]).doubleValue();
            if (lod > osac[i][j].maxLOD) {
              osac[i][j].maxLOD = lod;
              osac[i][j].maxPos = (k - 1);
              osac[i][j].maxFams = Integer.valueOf(cutLookup.elementAt(cut + 1)).intValue();
            }
          }
        }
        line = reader.readLine().split("[\\s]+");
        osac[i][j].baseLOD = Double.valueOf(line[osac[i][j].maxPos + 1]).doubleValue();
        reader.close();

        count = total = 0;
        if (new File(base + "chrom" + chrome + "/sims/perm" + chrome + ".out").exists()) {
          reader = new BufferedReader(
              new FileReader(base + "chrom" + chrome + "/sims/perm" + chrome + ".out"));
          while (reader.ready()) {
            line = reader.readLine().split("[\\s]+");
            if (Double.valueOf(line[2]).doubleValue() >= osac[i][j].maxLOD) {
              count++;
            }
            total++;
          }
          reader.close();
          if (count >= 20) {
            osac[i][j].pval = "=" + ext.formDeci((double) count / (double) total, 4, true);
          } else if (count == 0) {
            osac[i][j].pval = "<" + ext.formDeci(1.0 / total, 4, true);
            osac[i][j].repsRequested = total * 10 + 20;
          } else {
            osac[i][j].pval = "~" + ext.formDeci((double) count / (double) total, 4, true);
            osac[i][j].repsRequested = total * 20 / count + 20;
          }
        } else {
          osac[i][j].pval = "??";
          osac[i][j].repsRequested = 20;
        }
        osac[i][j].repsCompleted = total;
      }
    }

    for (int j = 0; j < 23; j++) {
      chr = j + 1;
      chrome = (chr < 10) ? "0" + chr : "" + chr;
      reps = (osac[0][j] != null ? osac[0][j].repsRequested : 0)
          + (osac[1][j] != null ? osac[1][j].repsRequested : 0);
      if (reps > 0) {
        if (reps > 200) {
          reps = 100;
        }
        batch.println("cd " + (new File(base).getAbsolutePath()) + "/chrom" + chrome);
        if (!new File(base + "chrom" + chrome + "/sims").exists()) {
          batch.println("mkdir sims");
        }
        batch.println("cd sims");
        batch.println(
            "cp ../../trait.dat ../../re_chrom" + chrome + ".pre ../../map" + chrome + ".dat .");
        batch.println(Files.getRunString() + " permOSA trait.dat " + chr + " " + reps + " d " + plug
            + " >> perm" + chrome + ".out");
      }
    }

    return osac;
  }
}
