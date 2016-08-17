package org.genvisis.dead;

import java.io.*;

import org.genvisis.common.*;

public class DrunkardsWalk {

  public static void run(double probability, int sampleSize, int replicates) throws IOException {
    int[] finalScore, timesFlipped, finalLead;
    double random;
    boolean zero;
    int prev, prevPrev;

    PrintWriter writer = new PrintWriter(new FileWriter("temp.xln"));
    finalScore = new int[replicates];
    timesFlipped = new int[replicates];
    finalLead = new int[replicates];
    int step = replicates / 100;
    for (int i = 0; i < replicates; i++) {
      prev = -99;
      prevPrev = -999;
      if (i % step == 0) {
        // System.out.println("Replicate: "+i);
      }
      finalScore[i] = 0;
      for (int j = 0; j < sampleSize; j++) {
        random = Math.random();
        zero = (finalScore[i] == 0 && j > 0);
        if (random >= probability) {
          finalScore[i]++;
        } else {
          finalScore[i]--;
        }
        if (j > 1 && zero && finalScore[i] != prevPrev) {
          timesFlipped[i]++;
          finalLead[i] = finalScore[i];
        }
        writer.print((j == 0 ? "" : "\t") + finalScore[i]);
        // System.out.println(i+"\t"+finalScore[i]+"\t"+timesFlipped[i]+"\t"+finalLead[i]);
        // writer.println(i+"\t"+random+"\t"+finalScore[i]+"\t"+timesFlipped[i]+"\t"+finalLead[i]);

        prevPrev = prev;
        prev = finalScore[i];
      }
      writer.println();
    }
    writer.close();

    CountVector cv = new CountVector();
    for (int i = 0; i < replicates; i++) {
      cv.add(timesFlipped[i] + "");
    }

    String[] values;
    int[] counts;
    values = cv.getValues();
    counts = cv.getCounts();

    for (int i = 0; i < values.length; i++) {
      System.out.println(values[i] + "\t" + counts[i]);
    }

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "Proof.dat";

    String usage = "\n" + "dead.DrunkardsWalk requires 0-1 arguments\n"
                   + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

    for (int i = 0; i < args.length; i++) {
      if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h")
          || args[i].equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (args[i].startsWith("file=")) {
        filename = args[i].split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      run(0.50, 20000, 100);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
