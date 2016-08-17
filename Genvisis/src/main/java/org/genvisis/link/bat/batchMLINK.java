package org.genvisis.link.bat;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.StringTokenizer;
import java.util.Vector;

public class batchMLINK {
  public static void main(String[] args) throws IOException {
    if (args.length != 1) {
      System.out.println("Expecting 1 argument: name of file containing dx models (Example: 0.00001 0.03 0.80 0.80)");
    } else {
      try {
        new batchMLINK(args[0]);
      } catch (Exception e) {
        e.printStackTrace();
      }
    }
  }

  @SuppressWarnings("resource")
  public batchMLINK(String models) throws IOException {
    BufferedReader reader = null;
    PrintWriter writer = null, writer1 = null, writer2 = null, writer3 = null, writer4 = null;
    String chrome;
    Vector<double[]> modelParams = new Vector<double[]>();
    StringTokenizer st;
    double[] handle;

    reader = new BufferedReader(new FileReader(models));

    while (reader.ready()) {
      st = new StringTokenizer(reader.readLine());
      handle = new double[4];
      handle[0] = Double.valueOf(st.nextToken()).doubleValue();
      handle[1] = Double.valueOf(st.nextToken()).doubleValue();
      handle[2] = Double.valueOf(st.nextToken()).doubleValue();
      handle[3] = Double.valueOf(st.nextToken()).doubleValue();
      modelParams.add(handle);
    }
    reader.close();

    writer1 = new PrintWriter(new FileWriter("batch.5"));
    writer2 = new PrintWriter(new FileWriter("batch.6"));
    writer3 = new PrintWriter(new FileWriter("batch.7"));
    writer4 = new PrintWriter(new FileWriter("batch.8"));
    writer1.println("#/bin/sh");
    writer1.println();
    writer1.println("sleep 30");
    writer1.println();
    writer2.println("#/bin/sh");
    writer2.println();
    writer2.println("sleep 30");
    writer2.println();
    writer3.println("#/bin/sh");
    writer3.println();
    writer3.println("sleep 30");
    writer3.println();
    writer4.println("#/bin/sh");
    writer4.println();
    writer4.println("sleep 30");
    writer4.println();

    for (int i = 1; i <= 22; i++) {
      chrome = (i < 10) ? "0" + i : "" + i;

      if (i >= 1 && i <= 3) {
        writer = writer1;
      }
      if (i >= 4 && i <= 7) {
        writer = writer2;
      }
      if (i >= 8 && i <= 14) {
        writer = writer3;
      }
      if (i >= 15 && i <= 22) {
        writer = writer4;
      }
      if (writer == null) {
        System.err.println("Error: Not all chromosomes were told which file to be in.");
      }

      for (int j = 0; j < modelParams.size(); j++) {
        handle = modelParams.elementAt(j);
        writer.println("cd chr" + i);
        writer.println("cp chrom" + chrome + ".ped pedin.dat");
        writer.println("cp map" + chrome + ".dat datain.dat");
        writer.println("jcp makeMap4MLINK " + handle[0] + " " + handle[1] + " " + handle[2] + " "
                       + handle[3]);
        writer.println("./pedin > /dev/null");
        writer.println("cp stream.out stream" + (j + 5) + ".out");
        writer.println("cd ..");
        writer.println();
      }

      writer.println();
      writer.println();
    }

    writer1.close();
    writer2.close();
    writer3.close();
    writer4.close();
  }
}
