package org.pankratzlab.internal.utils.widgets;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.PSF;

public class TabVersion {

  public static void make(String filename) {
    BufferedReader reader;
    PrintWriter writer;

    try {
      reader = new BufferedReader(new FileReader(filename));
      writer = Files.openAppropriateWriter(filename + ".xln");
      while (reader.ready()) {
        writer.println(ArrayUtils.toStr(reader.readLine().trim()
                                              .split(PSF.Regex.GREEDY_WHITESPACE)));
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "MakeTabVersion.dat";

    String usage = "\n" + "widgets.TabVersion requires 0-1 arguments\n"
                   + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

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
      make(filename);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}