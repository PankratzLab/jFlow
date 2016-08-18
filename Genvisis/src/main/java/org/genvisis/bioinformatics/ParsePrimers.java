package org.genvisis.bioinformatics;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

public class ParsePrimers {
  private static void forIlluminaInfiniumDesignScores(String dir, String variants, String sequence,
                                                      String output) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, seqs;
    String trav;
    Hashtable<String, String> hash, used;

    hash = new Hashtable<String, String>();
    try {
      reader = new BufferedReader(new FileReader(dir + sequence));
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        hash.put(line[1], line[2] + "\t" + line[0] + "\t" + line[3]);
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + sequence + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + sequence + "\"");
      System.exit(2);
    }

    used = new Hashtable<String, String>();
    try {
      reader = new BufferedReader(new FileReader(dir + variants));
      writer = new PrintWriter(new FileWriter(dir + variants + "_IlluminaDesign.csv"));
      writer.println("Locus_Name,Target_Type,Sequence,Chromosome,Coordinate,Genome_Build_Version,Source,Source_Version,Sequence_Orientation,Plus_Minus,Force_Infinium_I");
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        if (!line[0].startsWith("chr")) {
          line[0] = "chr" + line[0];
        }
        trav = line[0] + ":" + line[1];

        if (used.containsKey(trav)) {
          writer.print(trav + "_alt" + used.get(trav) + ",SNP,");
          used.put(trav, (Integer.parseInt(used.get(trav)) + 1) + "");
        } else {
          writer.print(trav + ",SNP,");
          used.put(trav, "1");
        }

        if (hash.containsKey(trav)) {
          seqs = hash.get(trav).split("\t", -1);
          if (!seqs[0].equals(line[2])) {
            System.err.println("Error - mismatched reference alleles between input and seq for position "
                               + trav + " (file=" + line[2] + ", seq=" + seqs[0] + ")");
          }
          writer.print(seqs[1] + "[" + line[3] + "/" + line[2] + "]" + seqs[2]);
        } else {
          System.err.println("Error - did not find sequence for variant: " + line[0] + ":"
                             + line[1]);
        }
        writer.println("," + line[0].substring(3) + "," + line[1]
                       + ",37,resequencing,0,Forward,Plus,FALSE");

      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + variants + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + variants + "\"");
      System.exit(2);
    }

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String dir = "";
    String variants = "variants.txt";
    // String sequence = "Galaxy28-[VariantSequence60flanked].tabular";
    String sequence = "Galaxy28-[VariantSequence60flanked]-allVariants.tabular";
    String output = "primers.csv";

    String usage = "\n" + "bioinformatics.ParsePrimers requires 0-1 arguments\n"
                   + "   (0) working directory (i.e. dir=" + dir + " (default))\n"
                   + "   (1) variants list filename (i.e. var=" + variants + " (default))\n"
                   + "   (2) flanking sequence file (i.e. seq=" + sequence + " (default))\n"
                   + "   (3) filename of output (i.e. out=" + output + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("var=")) {
        variants = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("seq=")) {
        sequence = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("out=")) {
        output = arg.split("=")[1];
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
      dir = "D:/Myron/Indian_Diabetes/SequencingPilot/SingaporeSelections_Designed/";
      forIlluminaInfiniumDesignScores(dir, variants, sequence, output);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
