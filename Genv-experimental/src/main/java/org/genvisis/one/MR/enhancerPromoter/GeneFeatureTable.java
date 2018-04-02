package org.genvisis.one.MR.enhancerPromoter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.ext;

public class GeneFeatureTable {

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

    createTable(filename, outFile);
  }

  private static void createTable(String filename, String outFile) {
    try {
      BufferedReader reader = new BufferedReader(new FileReader(filename));
      String[] header = reader.readLine().replaceAll("\"", "").split(",");
      int feature = ext.indexOfStr("feature_type", header);
      int genes = ext.indexOfStr("genes", header);

      Map<String, String[]> out = new HashMap<>();
      String[] outHeader = new String[] {"Gene", "Promoter", "Promoter_Flanking_Region",
                                         "Open_chromatin", "CTCF_Binding_Site", "TF_binding_site",
                                         "Enhancer"};

      while (reader.ready()) {
        String[] line = reader.readLine().replaceAll("\"", "").split(",");

        String[] g = line[genes].split("\\|");
        int featureIndex = ext.indexOfStr(line[feature], outHeader);
        for (String gene : g) {
          if (out.get(gene) == null) {
            out.put(gene, new String[] {gene, "0", "0", "0", "0", "0", "0"});
          }
          String[] temp = out.get(gene);
          temp[featureIndex] = Integer.parseInt(temp[featureIndex]) + 1 + "";
          out.put(gene, temp);
        }
      }

      BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
      writer.write(ArrayUtils.toStr(outHeader, ",") + "\n");

      for (String k : out.keySet()) {
        writer.write(ArrayUtils.toStr(out.get(k), ",") + "\n");
      }

      reader.close();
      writer.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
}
