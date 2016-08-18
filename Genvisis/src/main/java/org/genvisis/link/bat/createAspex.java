package org.genvisis.link.bat;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.StringTokenizer;
import java.util.Vector;

import org.genvisis.common.ext;

public class createAspex {

  public createAspex(int chromosome) throws IOException {
    BufferedReader map = null;
    BufferedReader struct = null;
    BufferedReader genotype = null;
    PrintWriter writer;
    StringTokenizer st;
    String temp, chrome, blank, id;
    boolean peachy = true;
    int numMarkers;
    Vector<String> markers = new Vector<String>();
    double[] markerDists;
    if (chromosome < 10) {
      chrome = "0" + chromosome;
    } else {
      chrome = "" + chromosome;
    }
    try {
      struct = new BufferedReader(new FileReader("struct.dat"));
    } catch (Exception e) {
      System.err.println("Could not open the family structure file: struct.dat");
      System.err.println("Please rectify");
      peachy = false;
    }
    try {
      map = new BufferedReader(new FileReader("map" + chrome + ".dat"));
    } catch (Exception e) {
      try {
        map = new BufferedReader(new FileReader("/home/npankrat/park/00masters/map" + chrome
                                                + ".dat"));
        System.err.println("Could not find map" + chrome + ".dat in the current directory");
        System.err.println("  using the one in /home/npankrat/park/00masters/");
      } catch (Exception e2) {
        System.err.println("Could not find map" + chrome
                           + ".dat in /home/npankrat/park/00masters/ or in the current directory");
        peachy = false;
      }

    }

    try {
      genotype = new BufferedReader(new FileReader("mrkr" + chrome + ".dat"));
    } catch (Exception e) {
      try {
        genotype = new BufferedReader(new FileReader("/home/npankrat/park/00masters/mrkr" + chrome
                                                     + ".dat"));
        System.err.println("Could not find mrkr" + chrome + ".dat in the current directory");
        System.err.println("  using the one in /home/npankrat/park/00masters/");
      } catch (Exception e2) {
        System.err.println("Could not find mrkr" + chrome
                           + ".dat in /home/npankrat/park/00masters/ or in the current directory");
        peachy = false;
      }

    }

    if (!peachy) {
      System.exit(1);
    }

    st = new StringTokenizer(map.readLine());
    numMarkers = Integer.valueOf(st.nextToken()).intValue() - 1;

    for (int i = 0; i < 6; i++) {
      map.readLine();
    }
    if (chromosome == 23) {
      map.readLine();
    }
    for (int i = 0; i < numMarkers; i++) {
      st = new StringTokenizer(map.readLine());
      st.nextToken();
      st.nextToken();
      st.nextToken();
      markers.add(st.nextToken());
      map.readLine();
    }
    map.readLine();
    st = new StringTokenizer(map.readLine());
    markerDists = new double[numMarkers + 1];
    for (int i = 0; i < numMarkers; i++) {
      markerDists[i] = Double.valueOf(st.nextToken()).doubleValue();
    }
    map.close();

    for (int j = 2; j < 5; j++) {
      switch (j) {
        default:
          writer = new PrintWriter(new FileWriter("chrom" + chrome + "h-mpt.param"));
          break;
        case 1:
          writer = new PrintWriter(new FileWriter("chrom" + chrome + "hd-mpt.param"));
          break;
        case 2:
          writer = new PrintWriter(new FileWriter("chrom" + chrome + "k-mpt.param"));
          break;
        case 3:
          writer = new PrintWriter(new FileWriter("chrom" + chrome + "kd-mpt.param"));
          break;
        case 4:
          writer = new PrintWriter(new FileWriter("chrom" + chrome + "k-2pt.param"));
          break;
      }
      writer.println("set nloc " + numMarkers);
      writer.print("set loc { ");
      for (int i = 0; i < numMarkers; i++) {
        writer.print(markers.elementAt(i) + " ");
      }
      writer.println("}");
      writer.println("set blank \"0\"");
      writer.print("set dist { ");
      if (j == 4) {
        for (int i = 0; i < numMarkers; i++) {
          writer.print("10.00 ");
        }
      } else {
        for (int i = 0; i < numMarkers; i++) {
          writer.print(ext.formDeci(markerDists[i] / 100, 4) + " ");
        }
      }
      writer.println("0.0 }");
      writer.println("set risk 2.0");
      writer.println("set z \"[expr 0.25/$risk] 0.50 [expr 0.50-0.25/$risk]\"");
      if (j < 2) {
        writer.println("set mapping Haldane");
      } else {
        writer.println("set mapping Kosambi");
      }
      writer.println("set most_likely true");
      if (j % 2 == 1) {
        writer.println("set no_Dv true");
      } else {
        writer.println("set no_Dv false");
      }
      if (chromosome == 23) {
        writer.println("set sex_linked true");
      }

      // writer.println("set fixed_freq true");
      writer.close();
    }

    writer = new PrintWriter(new FileWriter("chrom" + chrome + ".pre"));

    for (int i = 0; i < numMarkers; i++) {
      writer.print(markers.elementAt(i) + " ");
    }
    writer.println();

    temp = genotype.readLine();
    st = new StringTokenizer(temp);
    st.nextToken();
    st.nextToken();
    blank = "";
    while (st.hasMoreTokens()) {
      blank += "\t0";
      st.nextToken();
    }

    while (struct.ready()) {
      st = new StringTokenizer(struct.readLine());
      id = st.nextToken() + "\t" + st.nextToken();
      writer.print(id + "\t" + st.nextToken() + "\t" + st.nextToken() + "\t" + st.nextToken() + "\t"
                   + st.nextToken());
      if (!st.nextToken().startsWith("N")) {
        while (!temp.startsWith(id)) {
          temp = genotype.readLine();
        }
        st = new StringTokenizer(temp);
        st.nextToken(); // FamID
        st.nextToken(); // IndID
        while (st.hasMoreTokens()) {
          writer.print("\t" + st.nextToken());
        }
        writer.println();
      } else {
        writer.println(" " + blank);
      }
    }

    genotype.close();
    writer.close();
  }

  public static void main(String[] args) throws IOException {
    if (args.length != 1) {
      System.out.println("Expecting 1 argument: chromosome number.");
    } else {
      try {
        new createAspex(Integer.valueOf(args[0]).intValue());
      } catch (Exception e) {
        e.printStackTrace();
        System.err.println("Error in processing chromosome " + args[0]);
      }
    }
  }
}
