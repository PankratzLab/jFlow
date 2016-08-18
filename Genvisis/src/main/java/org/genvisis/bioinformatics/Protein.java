package org.genvisis.bioinformatics;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.common.Matrix;
import org.genvisis.common.ext;

public class Protein {
  public static final String[][] AMINO_ACID_LOOKUP =
      {
       // amino acid, three letter abbreviation, one letter abbreviation, polarity, acidity
       {"alanine", "Ala", "A", "nonpolar", "neutral"},
       {"arginine", "Arg", "R", "polar", "strongly basic"},
       {"asparagine", "Asn", "N", "polar", "neutral"},
       {"aspartic acid", "Asp", "D", "polar", "acidic"},
       {"cysteine", "Cys", "C", "nonpolar", "neutral"},
       {"glutamic acid", "Glu", "E", "polar", "acidic"},
       {"glutamine", "Gln", "Q", "polar", "neutral"},
       {"glycine", "Gly", "G", "nonpolar", "neutral"},
       {"histidine", "His", "H", "polar", "weakly basic"},
       {"isoleucine", "Ile", "I", "nonpolar", "neutral"},
       {"leucine", "Leu", "L", "nonpolar", "neutral"}, {"lysine", "Lys", "K", "polar", "basic"},
       {"methionine", "Met", "M", "nonpolar", "neutral"},
       {"phenylalanine", "Phe", "F", "nonpolar", "neutral"},
       {"proline", "Pro", "P", "nonpolar", "neutral"}, {"serine", "Ser", "S", "polar", "neutral"},
       {"threonine", "Thr", "T", "polar", "neutral"},
       {"tryptophan", "Trp", "W", "nonpolar", "neutral"},
       {"tyrosine", "Tyr", "Y", "polar", "neutral"}, {"valine", "Val", "V", "nonpolar", "neutral"}};

  public static String[] lookup(String letterCode) {
    int index = ext.indexOfStr(letterCode, Matrix.extractColumn(AMINO_ACID_LOOKUP, 2));
    return index == -1 ? null : AMINO_ACID_LOOKUP[index];
  }

  public static void parseAminoAcids(String filename) {
    BufferedReader reader;
    PrintWriter writer;
    String temp;
    int count;

    try {
      reader = new BufferedReader(new FileReader(filename));
      writer = new PrintWriter(new FileWriter(filename + "_aa.xln"));
      count = 0;
      while (reader.ready()) {
        temp = reader.readLine();
        if (temp.startsWith(">")) {
          System.out.println("Ignoring first line: " + temp);
        } else {
          for (int i = 0; i < temp.length(); i++) {
            count++;
            writer.println(count + "\t" + temp.charAt(i));
          }
        }
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
    String filename =
        "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\Reviews\\24 Rouleau random sequencing\\GCK5.txt";

    String usage = "\n" + "bioinformatics.Protein requires 0-1 arguments\n"
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
      parseAminoAcids(filename);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
