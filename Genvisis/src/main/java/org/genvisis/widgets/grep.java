package org.genvisis.widgets;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Vector;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipInputStream;

import org.genvisis.common.Array;
import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class grep {
  public static final int DEFAULT_NUM_LINES = 1000;

  public static void filter(String filename, String outputFilename, String[] withs,
                            String[] withouts) throws Elision {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String line;
    InputStreamReader isReader;
    boolean include;

    try {
      writer = Files.getAppropriateWriter(ext.parseDirectoryOfFile(filename)
                                          + (outputFilename.endsWith(".gz") ? outputFilename.substring(0,
                                                                                                       outputFilename.length()
                                                                                                          - 3)
                                                                            : outputFilename));
      isReader = null;
      if (outputFilename.endsWith(".gz")) {
        isReader = new InputStreamReader(new GZIPInputStream(new FileInputStream(filename)));
      } else if (outputFilename.endsWith(".zip")) {
        isReader = new InputStreamReader(new ZipInputStream(new FileInputStream(filename)));
      } else {
        isReader = new FileReader(filename);
      }

      reader = new BufferedReader(isReader);
      while (reader.ready()) {
        line = reader.readLine();
        include = false;
        for (String with : withs) {
          if (line.contains(with)) {
            include = true;
          }
        }
        for (String without : withouts) {
          if (line.contains(without)) {
            include = false;
          }
        }
        if (include) {
          writer.println(line);
        }
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      throw new Elision();
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      throw new Elision();
    }

  }

  public static void fromParameters(String filename, Logger log) throws Elision {
    Vector<String> params;
    String trav;
    String inputFilename = "large_file.txt";
    String outputFilename = "subset_of_lines_from_large_file.txt";
    ArrayList<String> withs, withouts;
    boolean not;

    params =
        Files.parseControlFile(filename, "grep",
                               new String[] {"input=" + inputFilename, "output=" + outputFilename,
                                             "# any line below will be treated as a string of interest, (though any line starting with a has \"#\" will still be considered a comment). If the line starts with an exclamation point (\"!\"), then this will be treated as a NOT, as in no line with this phrase will pass on to the final file",
                                             "first pattern", " second pattern ",
                                             "!exclude this phrase"},
                               log);

    if (params != null) {
      withs = new ArrayList<String>();
      withouts = new ArrayList<String>();
      for (int i = 0; i < params.size(); i++) {
        trav = params.elementAt(i);
        if (trav.startsWith("input=")) {
          inputFilename = ext.parseStringArg(trav, null);
        } else if (trav.startsWith("output=")) {
          outputFilename = ext.parseStringArg(trav, null);
        } else {
          if (trav.startsWith("!")) {
            not = true;
            trav = trav.substring(1);
          } else {
            not = false;
          }
          if (trav.startsWith(" ")) {
            log.report("Warning- \"" + trav
                       + "\" starts with whitespace, and will require an exact match");
          }
          if (trav.endsWith(" ")) {
            log.report("Warning- \"" + trav
                       + "\" ends with whitespace, and will require an exact match");
          }

          if (not) {
            withouts.add(trav);
          } else {
            withs.add(trav);
          }
        }
      }

      filter(inputFilename, outputFilename, Array.toStringArray(withs),
             Array.toStringArray(withouts));
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "input.txt";
    String outputFilename = "output.txt";
    String searchTerm = "pattern_to_find";
    String exclusionCriteria = "pattern_to_exclude";

    String usage =
        "\n" + "widgets.grep requires 0-1 arguments\n" + "   (1) input filename (i.e. file="
                   + filename + " (default))\n" + "   (2) output filename (i.e. out="
                   + outputFilename + " (default))\n" + "   (3) term to require (i.e. include="
                   + searchTerm + " (default))\n" + "   (4) term to exclude (i.e. exclude="
                   + exclusionCriteria + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("out=")) {
        outputFilename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("include=")) {
        searchTerm = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("exclude=")) {
        exclusionCriteria = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      filter(filename, filename + ".grepped", new String[] {searchTerm},
             new String[] {exclusionCriteria});
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
