package org.pankratzlab.internal.gwas;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.PSF;

public class ConvertToMDRformat {

  // public static final String DEFAULT_DIR = "C:\\Documents and
  // Settings\\npankrat\\My Documents\\gwas\\merged\\results\\";
  public static final String DEFAULT_DIR = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\JustY\\";

  // public static final String DEFAULT_FILES = "tophits.recode";
  // public static final String DEFAULT_FILES = "singleAdditive.recode";
  public static final String DEFAULT_FILES = "justY.recode";

  public static void convertToMDRformat(String dir, String prefix) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    Vector<Vector<String>> alleles = new Vector<>();
    Vector<String> markers;
    int count;

    markers = HashVec.loadFileToVec(dir + prefix + ".map", false, new int[] {1}, false);
    for (int i = 0; i < markers.size(); i++) {
      alleles.add(new Vector<String>());
    }

    try {
      reader = new BufferedReader(new FileReader(dir + prefix + ".ped"));
      writer = Files.openAppropriateWriter(dir + prefix + ".mdr");
      writer.println(ArrayUtils.toStr(markers) + "\tClass");
      while (reader.ready()) {
        line = reader.readLine().split(PSF.Regex.GREEDY_WHITESPACE);
        if (!line[5].equals("0")) {
          for (int i = 0; i < markers.size(); i++) {
            if (line[6 + i * 2 + 0].equals("0") || line[6 + i * 2 + 0].equals(".")
                || line[6 + i * 2 + 1].equals("0") || line[6 + i * 2 + 1].equals(".")) {
              writer.print(".\t");
            } else {
              count = 0;
              for (int j = 0; j < 2; j++) {
                if (alleles.elementAt(i).indexOf(line[6 + i * 2 + j]) == -1) {
                  alleles.elementAt(i).add(line[6 + i * 2 + j]);
                }
                count += alleles.elementAt(i).indexOf(line[6 + i * 2 + j]);
              }
              writer.print(count + "\t");
            }
          }
          writer.println(Integer.parseInt(line[5]) - 1);
        }
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + prefix + ".ped"
                         + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + prefix + ".ped" + "\"");
      System.exit(2);
    }

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String dir = DEFAULT_DIR;
    String prefix = DEFAULT_FILES;

    String usage = "\\n" + "park.gwa.ConvertToMDRformat requires 0-1 arguments\n"
                   + "   (1) prefix of ped/map files (i.e. prefix=" + prefix + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("prefix=")) {
        prefix = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      convertToMDRformat(dir, prefix);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}