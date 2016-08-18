// later we can expand this to do multiple swaps at once
package org.genvisis.park;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.StringTokenizer;

import org.genvisis.common.Files;

public class swapDNA {

  public swapDNA(String filename, int start, int stop) throws IOException {
    BufferedReader reader = null;
    String[] line;
    String temp;
    int count;

    try {
      reader = new BufferedReader(new FileReader(filename));
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error - Cannot find " + filename + " in current directory");
      System.exit(1);
    }

    count = 0;
    do {
      count++;
      temp = "Backup before swapDNA" + "(" + count + ")";
    } while (new File(temp).exists());
    new File(temp).mkdir();
    for (int i = start; i <= stop; i++) {
      if (!Files.copyFile("chromosome" + i + ".dat", temp + "/" + "chromosome" + i + ".dat")) {
        System.err.println("Error - could not copy " + "chromosome" + i + ".dat");
      }
    }

    while (reader.ready()) {
      temp = reader.readLine();
      line = temp.split("[\\s]+");
      if (line.length < 2) {
        System.err.println("Error - Could not process the following line (requires two DNA numbers):");
        System.err.println(temp);
        System.exit(3);
      }
      new swapDNA(line[0], line[1], start, stop);
    }

  }

  public swapDNA(String DNA1, String DNA2, int start, int stop) throws IOException {
    BufferedReader reader = null;
    PrintWriter writer = null;
    StringTokenizer st;
    String temp, dna, famInd1 = "", famInd2 = "";
    Hashtable<String, String> hash;
    String[] line;
    // boolean sorunvarmi = false;

    for (int i = start; i <= stop; i++) {
      try {
        reader = new BufferedReader(new FileReader("chromosome" + i + ".dat"));
        hash = new Hashtable<String, String>();
        while (reader.ready()) {
          line = reader.readLine().split("[\\s]+");
          if (line.length > 2) {
            if (hash.containsKey(line[1] + "\t" + line[2])) {
              System.err.println("  " + line[1] + "\t" + line[2] + " has multiple entries");
              // sorunvarmi = true;
            } else {
              hash.put(line[1] + "\t" + line[2], "ha");
            }
          } else {
            if (!line[0].startsWith("place")) {
              System.out.println(line[0]);
            }
          }
        }
        reader.close();
        // if (sorunvarmi) {
        // System.err.println("Program is exiting prematurely; remove
        // duplicate entires for the same individual and rerun");
        // System.exit(9);
        // }

        reader = new BufferedReader(new FileReader("chromosome" + i + ".dat"));
        while (reader.ready()) {
          st = new StringTokenizer(reader.readLine());
          dna = st.nextToken();
          if (dna.equals(DNA1)) {
            famInd1 = st.nextToken() + "\t" + st.nextToken();
          }
          if (dna.startsWith(DNA2)) {
            famInd2 = st.nextToken() + "\t" + st.nextToken();
          }
        }
        reader.close();

        if (famInd1.equals("")) {
          System.err.println("Error: could not find first DNA, " + DNA1 + ", in file \'chromosome"
                             + i + ".dat\'");
        } else if (famInd2.equals("") && !DNA2.equalsIgnoreCase("trash")) {
          System.err.println("Error: could not find second DNA, " + DNA2 + ", in file \'chromosome"
                             + i + ".dat\'");
        } else {
          (new File("chromosome" + i + ".dat")).renameTo(new File("temp"));
          reader = new BufferedReader(new FileReader("temp"));
          writer = new PrintWriter(new FileWriter("chromosome" + i + ".dat"));

          while (reader.ready()) {
            temp = reader.readLine();
            st = new StringTokenizer(temp);
            dna = st.nextToken();

            if (dna.equals(DNA1)) {
              if (!DNA2.equalsIgnoreCase("trash")) {
                writer.print(DNA1 + "\t" + famInd2);
                st.nextToken();
                st.nextToken();
                while (st.hasMoreTokens()) {
                  writer.print("\t" + st.nextToken());
                }
                writer.println();
              }
            } else if (dna.equals(DNA2)) {
              writer.print(DNA2 + "\t" + famInd1);
              st.nextToken();
              st.nextToken();
              while (st.hasMoreTokens()) {
                writer.print("\t" + st.nextToken());
              }
              writer.println();
            } else {
              writer.println(temp);
            }
          }
          reader.close();
          writer.close();
          (new File("temp")).delete();
        }
      } catch (Exception e) {
        e.printStackTrace();
      }
    }
    System.out.println("Make sure you re-sort chromosome files using park.sortFamInd");

  }

  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    String filename = "swapDNAs.dat", DNA1 = null, DNA2 = null;
    int start = 1, stop = 23;
    // int start=5, stop=5;

    String usage = "\n" + "park.swapDNA requires 1-3 arguments:\n"
                   + "   (1) either a file with two columns: DNA#1 DNA#2\n" + "       (i.e. file="
                   + filename + " (default))\n" + "                   OR\n"
                   + "       two arguments for the DNA numbers to switch\n"
                   + "       (i.e. DNA1=2000PD0137 DNA2=2000PD0138)\n"
                   + "       (note: DNA2 can equal 'trash' if you want to destroy)\n"
                   + "   (2/3) chromosome number, if only doing one\n"
                   + "         (i.e. chr=5 (default is all chromosomes))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("DNA1=")) {
        DNA1 = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("DNA2=")) {
        DNA2 = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("chr=")) {
        start = stop = Integer.valueOf(arg.split("=")[1]).intValue();
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    if (args.length == 0) {
      System.err.println("Warning: using defaults (file=" + filename + ", all chromosomes)");
    }
    if ((DNA1 == null && DNA2 != null) || (DNA1 != null && DNA2 == null)) {
      System.err.println("Error - requires 2 DNA numbers, not 1, in order to swap properly");
      System.exit(1);
    }
    try {
      if (DNA1 != null && DNA2 != null) {
        new swapDNA(DNA1, DNA2, start, stop);
      } else {
        new swapDNA(filename, start, stop);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }

  }
}
