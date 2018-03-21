package org.genvisis.one.MR.enhancerPromoter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

public class TallyEnhancerPromoter {

  public static void main(String[] args) {
    String filename = null;
    String outFile = null;

    for (String arg : args) {
      if (arg.startsWith("filename=")) {
        filename = arg.split("=")[1];
      } else if (arg.startsWith("out=")) {
        outFile = arg.split("=")[1];
      }
    }

    if (filename == null || outFile == null) {
      System.err.println("Missing required arguments: filename= , out= . Aborting.");
      System.exit(0);
    }

    tally(filename, outFile);
  }

  private static void tally(String filename, String outFile) {
    try {
      // set up file reader
      BufferedReader reader = new BufferedReader(new FileReader(filename));
      String[] header = reader.readLine().replaceAll("\"", "").split(",");
      int ensembl = indexOf(header, "Ensembl_Regulatory_Build_feature_type");
      int gene = indexOf(header, "SnpEff_ensembl_Gene_name");
      int chrom = indexOf(header, "CHROM");
      int pos = indexOf(header, "POS");

      Map<String, String[]> current = new HashMap<String, String[]>();

      String prevType = ".";
      String[] outHeader = new String[] {"feature_type", "chrom", "start", "stop", "genes",
                                         "count"};

      BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
      writer.write(toString(outHeader, ",") + "\n");

      while (reader.ready()) {
        String[] line = reader.readLine().replaceAll("\"", "").split(",");
        String[] types = line[ensembl].split(";");

        for (String t : types) {
          if (current.get(t) == null && !t.equals(".")) {
            current.put(t, new String[] {t, line[chrom], line[pos], line[pos],
                                         combineGenes("", line[gene]), "1"});
          } else if (current.get(t) != null) {
            String[] temp = current.get(t);
            temp[3] = line[pos];
            temp[5] = Integer.parseInt(temp[5]) + 1 + "";
            // deal with genes list
            temp[4] = combineGenes(temp[4], line[gene]);
            current.put(t, temp);
          }
        }

        String[] prevTypes = prevType.split(";");
        for (String pt : prevTypes) {
          if (indexOf(types, pt) == -1 && current.get(pt) != null) {
            // if this type is not in the new type list, this segment has ended
            // write out the existing segment and reset the map
            writer.write(toString(current.get(pt), ",") + "\n");
            current.put(pt, null);
          }
        }

        prevType = line[ensembl];
      }

      reader.close();
      writer.close();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  private static int indexOf(String[] array, String val) {
    for (int i = 0; i < array.length; i++) {
      if (array[i].equals(val)) return i;
    }
    return -1;
  }

  private static String toString(String[] list, String sep) {
    String out = "";
    for (int i = 0; i < list.length; i++) {
      out += list[i] + (i == list.length - 1 ? "" : sep);
    }
    return out;
  }

  private static String combineGenes(String a, String b) {
    if (a.equals(b) || b.length() == 0) return a;

    String[] aGenes = a.split("\\|");
    String[] bGenes = b.split("\\|");

    String out = a;

    ArrayList<String> seen = new ArrayList<String>();
    Collections.addAll(seen, aGenes);

    for (String gene : bGenes) {
      if (!seen.contains(gene)) {
        out += (out.equals("") ? gene : "|" + gene);
        seen.add(gene);
      }
    }

    return out;
  }
}
