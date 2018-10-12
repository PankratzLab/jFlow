package org.genvisis.one.link.bat;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import org.pankratzlab.common.Files;

public class gatherOSAhauser {

  public gatherOSAhauser() throws IOException {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String chrome, skipped = "";

    writer = Files.openAppropriateWriter("osaHauserSummary.out");
    for (int chromosome = 1; chromosome <= 23; chromosome++) {
      chrome = (chromosome < 10) ? "0" + chromosome : "" + chromosome;
      try {
        writer.println("Chromosome " + chromosome + ":");
        reader = new BufferedReader(new FileReader("chrom" + chrome + "/osa" + chrome
                                                   + ".dat.max"));
        for (int i = 0; i < 14; i++) {
          reader.readLine();
        }

        while (reader.ready()) {
          writer.println(reader.readLine());
        }
        writer.println();

      } catch (IOException ioe) {
        writer.println("no data");
        skipped += " " + chromosome;
      }
    }
    if (skipped.length() > 0) {
      System.err.println("Skipped chromosomes" + skipped);
    }
    writer.close();

  }

  public static void main(String[] args) {
    if (args.length != 0) {
      System.out.println("Expecting no arguments.");
    } else {
      try {
        new gatherOSAhauser();
      } catch (Exception e) {
        e.printStackTrace();
      }
    }
  }
}
