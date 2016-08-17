package org.genvisis.link.bat;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.StringTokenizer;

import org.genvisis.common.Files;

public class solarZeroByInd {
  public static void main(String[] args) throws IOException {
    if (args.length < 1 || args.length > 2) {
      System.err.println(
          "Error: requires 1-2 arguments - individuals to keep (like a .pre file) and a chromosome [optional]");
    } else {
      try {
        if (args.length == 1) {
          new solarZeroByInd(args[0]);
        }
        if (args.length == 2) {
          new solarZeroByInd(args[0], Integer.valueOf(args[1]).intValue(),
              Integer.valueOf(args[1]).intValue());
        }
      } catch (Exception e) {
        e.printStackTrace();
      }
    }
  }

  public solarZeroByInd(String fams) throws IOException {
    new solarZeroByInd(fams, 1, 23);
  }

  public solarZeroByInd(String fams, int start, int stop) throws IOException {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String temp, blank, id;
    StringTokenizer st;
    Hashtable<String, String> hash = new Hashtable<String, String>();

    reader = new BufferedReader(new FileReader(fams));
    while (reader.ready()) {
      st = new StringTokenizer(reader.readLine());
      hash.put(st.nextToken() + "," + st.nextToken(), "");
    }
    reader.close();

    for (int chromosome = start; chromosome <= stop; chromosome++) {
      String bakFilename =
          Files.getBakFilename("solar.gtypes." + chromosome, super.getClass().getName());
      (new File("solar.gtypes." + chromosome)).renameTo(new File(bakFilename));
      reader = new BufferedReader(new FileReader(bakFilename));
      writer = new PrintWriter(new FileWriter("solar.gtypes." + chromosome));
      writer.println(reader.readLine());
      temp = reader.readLine();
      st = new StringTokenizer(temp, "/");
      blank = "0/0";
      for (int i = 0; i < st.countTokens() - 2; i++) {
        blank += ",0/0";
      }
      boolean done = false;
      while (!done) {
        st = new StringTokenizer(temp, ",");
        id = st.nextToken() + "," + st.nextToken();
        if (hash.containsKey(id)) {
          writer.println(temp);
        } else {
          writer.println(id + "," + blank);
        }
        if (reader.ready()) {
          temp = reader.readLine();
        } else {
          done = true;
        }
      }

      reader.close();
      writer.close();

    }
  }
}
