package org.genvisis.link;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.StringTokenizer;

import org.genvisis.common.ext;

public class simMakeAllegro {
  public static int NUM_REPS = 500;

  public static String[] FILES = {"2sibs", "3sibs", "4sibs"};

  // public static String[] FILES = {"reconstruct"};
  // public static int[] FAMREPS_MASTA = {335, 24, 3}; // plates 1-10
  // public static int[] FAMREPS_MASTA = {385, 24, 3}; // with 400 families
  public static int[] FAMREPS_MASTA = {790, 56, 7}; // to total 1000

  // sibpairs
  // public static int[] FAMREPS = {516, 36, 6};
  // public static int[] FAMREPS = {774, 54, 9};
  // public static int[] FAMREPS = {500};

  public static void main(String[] args) throws IOException {
    if (args.length != 2) {
      System.out.println("Expecting 2 arguments: locus heterogeneity, reps multiplier.");
    } else {
      try {
        new simMakeAllegro(Double.valueOf(args[0]).doubleValue(),
            Integer.valueOf(args[1]).intValue());
        // stdin.readLine();
      } catch (Exception e) {
        e.printStackTrace();
      }
    }
  }

  public simMakeAllegro(double het, int multiplier) throws IOException {
    BufferedReader reader = null;
    PrintWriter writer = null;
    PrintWriter optfile = null;
    String file = null;

    int[] FAMREPS = new int[FAMREPS_MASTA.length];
    for (int i = 0; i < FAMREPS_MASTA.length; i++) {
      FAMREPS[i] = FAMREPS_MASTA[i] * multiplier;
    }

    for (int i = 0; i < FILES.length; i++) {
      file = FILES[i];
      optfile = new PrintWriter(new FileWriter(file + ".opt"));
      optfile.println("% Read input in LINKAGE style format:\n" + "PREFILE " + file + ".pre\n"
          + "DATFILE linkage.dat\n\n" + "% Simulate stroke reconstruction pedigrees\n"
          + "SIMULATE dloc:32.0 npre:" + NUM_REPS + " rep:" + FAMREPS[i]
          + " err:0.00 yield:1.0 het:" + ext.formDeci(het, 2, true) + "\n\n" // change
                                                                             // from
                                                                             // 8 to
          // 32 for denser map
          // wtf? 32 is
          // disease allele
          // position
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
    }

    for (int repNum = 1; repNum <= NUM_REPS; repNum++) {
      writer = new PrintWriter(new FileWriter("linkage-" + repNum + ".pre"));
      for (int i = 0; i < FILES.length; i++) {
        file = FILES[i] + ".pre." + ext.formNum(repNum + "", String.valueOf(NUM_REPS).length());
        reader = new BufferedReader(new FileReader(file));
        while (reader.ready()) {
          writer.println(translate(reader.readLine(), i));
        }
        reader.close();
        (new File(file)).delete();
      }
      writer.close();
    }
  }

  public String translate(String str, int structureNum) throws IOException {
    String translation;
    StringTokenizer st = new StringTokenizer(str, "- ");

    st.nextToken();
    translation = (structureNum + 1) + ext.formNum(st.nextToken(), 5);
    for (int i = 0; i < 21; i++) { // 21 for 8 markers; 69 for 32 markers
      translation += "  " + st.nextToken();
    }
    return translation;
  }
}
