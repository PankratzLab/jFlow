package org.genvisis.gwas;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class Metasoft {
  public static void applyLambdas(String inputFile, String outputFile, double[] lambdas,
      Logger logger) {
    String[] line;

    try {
      BufferedReader reader = new BufferedReader(new FileReader(inputFile));
      PrintWriter writer = new PrintWriter(new FileWriter(outputFile));
      for (int i = 0; i < lambdas.length; i++) {
        if (lambdas[i] < 1) {
          lambdas[i] = 1;
        }
      }
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        if (line.length != lambdas.length * 2 + 1) {
          System.err.println("Error - mismatched number of lambdas and/or columns in '" + inputFile
              + "' (with " + lambdas.length + " lambdas, expecting " + (lambdas.length * 2 + 1)
              + " (1+" + lambdas.length + "*2) columns");
        }
        for (int i = 0; i < lambdas.length; i++) {
          if (!line[1 + i * 2 + 1].equals("NA")) {
            line[1 + i * 2 + 1] =
                ext.formDeci(Double.parseDouble(line[1 + i * 2 + 1]) * lambdas[i], 9);
          }
        }
        writer.println(Array.toStr(line));
      }
      writer.close();
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + inputFile + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + inputFile + "\"");
      System.exit(2);
    }

  }

  public static void generateRawFile(String[] markers, String[] filenamesAndColumnIDs,
      String outputFilename, Logger log) {
    Files.combineWithLessMemory(markers, filenamesAndColumnIDs, null, "MarkerName", ".",
        outputFilename, log, true, false, false, false);
  }

  public static void processRawToInputFile(String rawFilename, String outputFilename, Logger log) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String[] stdAlleles;
    int numStudies;
    int[] numOpposite, numFlipped;

    try {
      reader = new BufferedReader(new FileReader(rawFilename));
      reader.mark(10000);
      line = reader.readLine().trim().split("[\\s]+");
      numStudies = (line.length - 1) / 4;
      if (line.length != numStudies * 4 + 1) {
        System.err.println(
            "Error determining number of studies, are you now parsing more than 4 parameters?");
        reader.close();
        return;
      }
      reader.reset();
      numOpposite = new int[numStudies];
      numFlipped = new int[numStudies];
      writer = new PrintWriter(new FileWriter(outputFilename));
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        stdAlleles = null;
        writer.print(line[0]);
        // for (int i = 5; i < line.length; i += 4) {
        for (int i = 0; i < numStudies; i++) {
          if (line[1 + i * 4 + 0].equals(".") || line[1 + i * 4 + 1].equals(".")) {
            // then ignore
          } else if (stdAlleles == null) {
            stdAlleles = new String[] {line[1 + i * 4 + 0], line[1 + i * 4 + 1]};
          } else {
            switch (Metal.determineStrandConfig(
                new String[] {line[1 + i * 4 + 0], line[1 + i * 4 + 1]}, stdAlleles)) {
              case Metal.STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND:
                numFlipped[i]++;
              case Metal.STRAND_CONFIG_SAME_ORDER_SAME_STRAND:
                break;
              case Metal.STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND:
                numFlipped[i]++;
              case Metal.STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND:
                numOpposite[i]++;
                line[1 + i * 4 + 2] = (Double.parseDouble(line[1 + i * 4 + 2]) * -1) + "";
                numOpposite[i]++;
              case Metal.STRAND_CONFIG_BOTH_NULL:
                break;
              case Metal.STRAND_CONFIG_DIFFERENT_ALLELES:
                System.err.println("Error - for marker " + line[0] + " Study" + (i + 1)
                    + " has different alleles (" + line[1 + i * 4 + 0] + "/" + line[1 + i * 4 + 1]
                    + ") than the rest (" + stdAlleles[0] + "/" + stdAlleles[1] + ")");
                break;
              case Metal.STRAND_CONFIG_SPECIAL_CASE:
                System.err.println("Warning - marker " + line[0]
                    + " has a special case starting with Study" + (i + 1) + ": alleles ("
                    + line[1 + i * 4 + 0] + "/" + line[1 + i * 4 + 1]
                    + ") where previous had only (" + stdAlleles[0] + "/" + stdAlleles[1] + ")");
                break;
              default:
                System.err.println("Error - unknown determineStrandConfig return code");
                break;
            }
          }
          writer.print("\t" + line[1 + i * 4 + 2] + "\t" + line[1 + i * 4 + 3]);
          // writer.print("\t"+line[1+i*4+0]+"\t"+line[1+i*4+1]+"\t"+line[1+i*4+2]+"\t"+line[1+i*4+3]);
        }
        writer.println();
      }
      writer.close();
      reader.close();
      System.out.println("Study\ttimeFlippedStrand\ttimesOppositeOrder");
      for (int i = 0; i < numStudies; i++) {
        System.out.println((i + 1) + "\t" + numFlipped[i] + "\t" + numOpposite[i]);
      }

    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + rawFilename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + rawFilename + "\"");
      System.exit(2);
    }
  }
}
