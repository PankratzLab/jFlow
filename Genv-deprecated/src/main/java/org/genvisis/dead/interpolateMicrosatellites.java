package org.genvisis.dead;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.StringTokenizer;
import java.util.Vector;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Sort;

public class interpolateMicrosatellites {

  public interpolateMicrosatellites(String oriname, String finename) throws IOException {
    BufferedReader reader = null;
    PrintWriter writer = null;
    StringTokenizer st;
    String temp, trav, recall, id, interName = "interpolated.dat";
    Hashtable<String, String> hash = new Hashtable<>();
    Vector<String> v;
    int[] keys;

    try {
      reader = new BufferedReader(new FileReader(finename));
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error - could not find the fine mapping '" + finename + "' file.");
      System.exit(1);
    }

    reader.readLine();
    recall = reader.readLine();
    while (reader.ready()) {
      temp = reader.readLine();
      st = new StringTokenizer(temp);
      trav = st.nextToken();

      hash.put(trav, temp);
    }
    reader.close();

    try {
      reader = new BufferedReader(new FileReader(oriname));
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error - could not find the original chromesome file, '" + finename
                         + "'.");
      System.exit(1);
    }

    if (oriname.startsWith("chromosome")) {
      if (oriname.length() >= 12) {
        interName = "chromosome";
        try {
          if (Integer.valueOf(oriname.substring(10, 12)).intValue() <= 23) {
            interName += Integer.valueOf(oriname.substring(10, 12)).intValue() + "intr.dat";
          }
        } catch (NumberFormatException nfe) {}

        if (interName.equals("chromosome")) {
          try {
            interName += Integer.valueOf(oriname.substring(10, 11)).intValue() + "intr.dat";
          } catch (NumberFormatException nfe) {
            interName = oriname + "-interpolated.dat";
          }
        }
      }

    } else {
      interName = oriname + "-interpolated.dat";
    }

    writer = Files.openAppropriateWriter(interName);

    writer.println(reader.readLine());
    writer.print(reader.readLine());
    st = new StringTokenizer(recall);
    st.nextToken();
    st.nextToken();
    st.nextToken();
    while (st.hasMoreTokens()) {
      writer.print("\t" + st.nextToken() + "\t");
    }
    writer.println();

    while (reader.ready()) {
      temp = reader.readLine();
      st = new StringTokenizer(temp);
      trav = st.nextToken();
      id = st.nextToken() + "-" + st.nextToken();

      if (hash.containsKey(trav)) {
        recall = hash.get(trav);
        writer.print(temp);
        st = new StringTokenizer(recall);
        st.nextToken();
        st.nextToken();
        st.nextToken();

        while (st.hasMoreTokens()) {
          writer.print("\t" + st.nextToken());
        }
        writer.println();
        hash.remove(trav);
      } else {
        System.err.println("Warning - individual " + id + " (" + trav
                           + ") does not have finemapping data");
      }
    }
    reader.close();
    writer.close();

    Enumeration<String> enumer = hash.keys();
    v = new Vector<>();
    while (enumer.hasMoreElements()) {
      st = new StringTokenizer(hash.get(enumer.nextElement()));
      trav = st.nextToken();
      id = st.nextToken() + "-" + st.nextToken();

      v.add(id + " (" + trav + ")");
    }
    if (v.size() > 0) {
      keys = Sort.getSortedIndices(ArrayUtils.toStringArray(v));
      System.err.println("Warning - the following IDs did not have an original entry.");
      for (int i = 0; i < v.size(); i++) {
        System.err.println("        - " + v.elementAt(keys[i]));
      }
    }
  }

  public static void main(String[] args) throws IOException {
    String usage = "\n" + "interpolateMicrosatellites requires 2 arguments:\n"
                   + "   (1) the original chromosome file (i.e. 'chromosome5ori.dat')\n"
                   + "   (2) the fine mapping chromosome file (i.e. 'chromosome5fine.dat')\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      }
    }
    if (args.length != 2) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      new interpolateMicrosatellites(args[0], args[1]);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
