package org.genvisis.kaput;

import java.io.*;
import java.util.*;

public class genoByFam {
  public genoByFam() throws IOException {
    BufferedReader reader = null;
    Vector<String> fams = new Vector<String>();
    Hashtable<String, String> hash = new Hashtable<String, String>();
    PrintWriter writer = new PrintWriter(new FileWriter("genoByFam"));
    // PrintWriter writer = new PrintWriter(new FileWriter("correctByFam"));
    String temp, chrome, i;
    StringTokenizer st;

    reader = new BufferedReader(new FileReader("solar_ped"));
    reader.readLine();
    while (reader.ready()) {
      st = new StringTokenizer(reader.readLine(), ",");
      temp = st.nextToken();
      if (!hash.containsKey(temp)) {
        hash.put(temp, ";");
        fams.add(temp);
      }
    }
    reader.close();

    writer.println("#/bin/sh");
    writer.println();
    writer.println("sleep 20");
    writer.println();

    for (int fam = 0; fam < fams.size(); fam++) {
      i = fams.elementAt(fam);
      writer.println("mkdir fam" + i);
      writer.println("java -classpath /home/npankrat/" + org.genvisis.common.PSF.Java.GENVISIS
                     + " park.zeroByFam " + i);
      writer.println("mv solar_marker." + i + ".* fam" + i);
      writer.println("cp solar_freq.* fam" + i);
      writer.println("cp solar_map.* fam" + i);
      writer.println("cp solar_ped fam" + i);
      // writer.println("cp newAAO.dat fam"+i);
      writer.println("cp batch.this fam" + i);
      writer.println("cp correct.aao fam" + i);
      // writer.println("cp correct.this fam"+i);
      writer.println("cd fam" + i);
      for (int chromosome = 2; chromosome <= 2; chromosome++) {
        chrome =
            (Integer.valueOf(chromosome + "").intValue() < 10) ? "0" + chromosome : "" + chromosome;
        writer.println("mv solar_marker." + i + "." + chrome + " solar_marker." + chrome);
      }
      writer.println("chmod +x batch.this");
      writer.println("./batch.this");
      // writer.println("./correct.this");
      writer.println("cd ..");
      writer.println();
    }
    writer.close();

    Runtime.getRuntime().exec("chmod +x correctByFam");
  }

  public static void main(String[] args) throws IOException {
    try {
      new genoByFam();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
