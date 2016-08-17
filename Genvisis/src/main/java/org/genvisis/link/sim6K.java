package org.genvisis.link;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.common.ext;

public class sim6K {
  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    String pre = "fam.pre", map = "6Kpanel_map.dat";
    int numReps = 1000;

    String usage =
        "\n" + "park.sim6K requires 0-1 arguments\n" + "   (1) family structure file (i.e. pre="
                   + pre + " (default)\n" + "   (2) map data (i.e. map=" + map + " (default)\n"
                   + "   (3) number of replicates (i.e. reps=" + numReps + " (default)\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("pre=")) {
        pre = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("map=")) {
        map = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("reps=")) {
        numReps = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      new sim6K(pre, map, numReps);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public sim6K(String struct, String map, int reps) throws IOException {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] line;
    String chrome;

    PrintWriter optfile = null;
    String file = null;

    for (int i = 1; i <= 23; i++) {
      chrome = (i < 10) ? "0" + i : "" + i;

      optfile = new PrintWriter(new FileWriter(file + ".opt"));
      optfile.println("% Read input in LINKAGE style format:\n" + "PREFILE " + struct + "\n"
                      + "DATFILE map" + chrome + ".dat\n\n"
                      + "% Simulate stroke reconstruction pedigrees\n" + "SIMULATE dloc:32.0 npre:"
                      + reps + " rep:" + 1 + " err:0.00 yield:1.0 het:1.0\n\n"
                      + "% Other options:\n" + "MAXMEMORY 100");
      optfile.close();

      Process process = null;
      Runtime runtime = Runtime.getRuntime();
      process = runtime.exec("/software/bin/allegro " + file + ".opt");

      try {
        process.waitFor();
      } catch (Exception e) {
        e.printStackTrace();
      }

      (new File(file + ".opt")).delete();

      for (int repNum = 1; repNum <= reps; repNum++) {
        writer = new PrintWriter(new FileWriter("linkage-" + chrome + "" + repNum + ".pre"));
        file = struct + "." + ext.formNum(repNum + "", String.valueOf(reps).length());
        reader = new BufferedReader(new FileReader(file));
        while (reader.ready()) {
          line = reader.readLine().split("[-\\s]+");
          writer.print(1 + ext.formNum(line[0], 5));
          for (int j = 1; j < line.length; j++) {
            writer.print(" " + line[j]);
          }
          writer.println();
        }
        reader.close();
        (new File(file)).delete();
        writer.close();
      }
    }
  }
}
