package org.genvisis.link.bat;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.StringTokenizer;

public class solarZeroByFam {
  public solarZeroByFam(int fam) throws IOException {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String temp, chrome, blank, id;
    StringTokenizer st;

    for (int chromosome = 1; chromosome <= 1; chromosome++) {
      chrome = (Integer.valueOf(chromosome + "").intValue() < 10) ? "0" + chromosome
                                                                  : "" + chromosome;
      reader = new BufferedReader(new FileReader("solar_marker." + chrome));
      writer = new PrintWriter(new FileWriter("solar_marker." + fam + "." + chrome));
      writer.println(reader.readLine());
      temp = reader.readLine();
      st = new StringTokenizer(temp, "/");
      blank = " 0/ 0";
      for (int i = 0; i < st.countTokens() - 2; i++) {
        blank += ", 0/ 0";
      }
      while (reader.ready()) {
        st = new StringTokenizer(temp, ",");
        id = st.nextToken();
        if (id.equals(fam + "")) {
          writer.println(temp);
        } else {
          writer.println(id + "," + st.nextToken() + "," + blank);
        }
        temp = reader.readLine();
      }

      reader.close();
      writer.close();

    }
  }

  public static void main(String[] args) throws IOException {
    try {
      new solarZeroByFam(Integer.valueOf(args[0]).intValue());
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
