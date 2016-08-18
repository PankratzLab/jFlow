package org.genvisis.bioinformatics;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.common.Array;
import org.genvisis.common.HashVec;
import org.genvisis.common.Internat;

public class gatherOmim {
  public gatherOmim(String filename) {
    BufferedReader reader = null;
    PrintWriter writer, errors;
    String[] ids;
    String temp;
    int count;

    ids = Array.toStringArray(HashVec.loadFileToVec(filename, false, false, true));

    try {
      writer = new PrintWriter(new FileWriter("omim_names.out"));
      errors = new PrintWriter(new FileWriter("omim_errors.out"));
      for (String id : ids) {
        Internat.downloadFile("http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=OMIM&dopt=Detailed&tmpl=dispomimTemplate&list_uids="
                              + id, "omim.dat", true);
        try {
          reader = new BufferedReader(new FileReader("omim.dat"));
          count = 0;
          while (reader.ready()) {
            temp = reader.readLine();
            if (temp.indexOf("<tr><td align=\"left\" colspan=\"2\" ") >= 0) {
              writer.println(id + "\t"
                             + temp.substring(temp.indexOf("<b>") + 3, temp.indexOf("</b>"))
                             + (count > 0 ? "\t*******" : ""));
              writer.flush();
              count++;
            }
          }
          reader.close();
        } catch (FileNotFoundException fnfe) {
          errors.println(id + "\tfailed to download");
          errors.flush();
        } catch (IOException ioe) {
          errors.println(id + "\tfailed to parse");
          errors.flush();
        }
        try {
          Thread.sleep(200);
        } catch (InterruptedException ex) {
        }
        new File("omim.dat").delete();
      }
      writer.close();
      errors.close();
    } catch (IOException ioe) {
      System.err.println("Global error with report files.");
    }

  }

  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    String filename = "gatherOmim.dat";

    String usage = "\n" + "park.gatherOmim requires 0-1 arguments\n" + "   (1) filename (i.e. file="
                   + filename + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      new gatherOmim(filename);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
