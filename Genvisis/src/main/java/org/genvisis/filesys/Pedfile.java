package org.genvisis.filesys;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;
import java.util.Hashtable;

import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.Vectors;
import org.genvisis.common.ext;

import com.google.common.primitives.Ints;

public class Pedfile {
  public static final String[] MISSING_VALUES = {"0"};

  public static void main(String[] args) {
    int numArgs = args.length;
    String pheno = "pheno.dat";
    int col = 2;
    String missVal = "x";

    String usage = "\n" + "filesys.Pedfile requires 0-3 arguments\n"
                   + "   (1) phenotype file to use (i.e. pheno=" + pheno + " (default))\n"
                   + "   (2) columns to use (i.e. col=" + col + " (default))\n"
                   + "   (3) value to use for those missing phenotype (i.e. missing=" + missVal
                   + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("pheno=")) {
        pheno = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("col=")) {
        col = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("missing=")) {
        missVal = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      updatePhenotype(pheno, col, missVal);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void updatePhenotype(String phenotypeFile, int phenoCol, String missingValue) {
    // BufferedReader reader;
    PrintWriter writer;
    // String[] line;
    Hashtable<String, String> hash;
    IntVector[] ivs;
    long time;

    hash = HashVec.loadFileToHashString(phenotypeFile, new int[] {0, 1}, new int[] {phenoCol},
                                        false, "\t", false, false, false);

    time = new Date().getTime();
    ivs = Vectors.initializedArray(IntVector.class, 3);
    for (int chr = 1; chr <= 23; chr++) {
      try {
        // reader = new BufferedReader(new
        // FileReader(Files.getBakFilename("re_chrom"+ext.chrome(chr)+".pre",
        // ext.rootOf(phenotypeFile)+"-"+phenoCol, true)));
        // writer = new PrintWriter(new FileWriter("re_chrom"+ext.chrome(chr)+".pre"));
        // while (reader.ready()) {
        // line = reader.readLine().trim().split("[\\s]+");
        // if (hash.containsKey(line[0]+"\t"+line[1])) {
        // line[5] = hash.get(line[0]+"\t"+line[1]);
        // } else {
        // line[5] = missingValue;
        // }
        // writer.println(Array.toStr(line));
        // }
        // reader.close();

        int c;
        String trav;
        char DELIMITER = '\t';
        FileReader in;
        String pheno;

        in = new FileReader(Files.getBakFilename("re_chrom" + ext.chrome(chr) + ".pre",
                                                 ext.rootOf(phenotypeFile) + "-" + phenoCol, true));
        writer = new PrintWriter(new FileWriter("re_chrom" + ext.chrome(chr) + ".pre"));

        c = 0;
        while (c >= 0) {
          trav = "";
          while ((char) (c = in.read()) != DELIMITER && c >= 0) {
            trav += (char) c;
          }
          trav += "\t";
          while ((char) (c = in.read()) != DELIMITER && c >= 0) {
            trav += (char) c;
          }
          if (hash.containsKey(trav)) {
            pheno = hash.get(trav);
          } else {
            pheno = missingValue;
          }
          trav += "\t";
          for (int i = 0; i < 3; i++) {
            while ((char) (c = in.read()) != DELIMITER && c >= 0) {
              trav += (char) c;
            }
            trav += "\t";
          }
          if (c >= 0) {
            writer.print(trav + pheno + "\t");
          }
          while ((char) (c = in.read()) != DELIMITER && c >= 0) {
            ;
          }

          while ((char) (c = in.read()) != '\n' && c >= 0) {
            if ((char) c != '\r') {
              writer.print((char) c);
            }
          }
          if (c >= 0) {
            writer.println();
          }
        }

        in.close();
        writer.close();
        ivs[0].add(chr);
      } catch (FileNotFoundException fnfe) {
        ivs[1].add(chr);
      } catch (IOException ioe) {
        System.err.println(ioe.getMessage());
        ivs[0].add(chr);
      }
    }
    System.out.println("Finished in " + ext.getTimeElapsed(time));

    System.out.println("Successfully updated chromosomes: "
                       + (ivs[0].size() == 0 ? "none" : ext.listRanges(Ints.toArray(ivs[0]))));
    if (ivs[1].size() > 0) {
      System.out.println("Missing files for chromosomes: " + ext.listRanges(Ints.toArray(ivs[1])));
    }
    if (ivs[2].size() > 0) {
      System.out.println("Error parsing files for chromosomes: "
                         + ext.listRanges(Ints.toArray(ivs[2])));
    }
  }

  private FamilyStructure famStruct;
  private int[][] markerCounts;
  private int[] genotypeCounts;

  private char[][] alleles;

  private double[] freqs;

  public Pedfile(String filename) {
    BufferedReader reader;
    String[] line;
    String[][] ids;
    byte[] genders;
    byte[] affections;
    int count, numMarkers;
    char allele;

    try {
      reader = new BufferedReader(new FileReader(filename));
      count = 0;
      numMarkers = -1;
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        if (numMarkers == -1) {
          numMarkers = (line.length - 6) / 2;
        } else if ((line.length - 6) / 2 != numMarkers) {
          System.err.println("Error - different number of markers listed on line " + (count + 1)
                             + " (" + line[0] + "/" + line[1] + ")");
        }
        count++;
      }
      reader.close();

      reader = new BufferedReader(new FileReader(filename));
      ids = new String[count][];
      genders = new byte[count];
      affections = new byte[count];
      markerCounts = new int[count][numMarkers];
      alleles = new char[numMarkers][2];
      freqs = new double[numMarkers];
      genotypeCounts = new int[numMarkers];
      for (int i = 0; i < count; i++) {
        line = reader.readLine().trim().split("[\\s]+");
        ids[i] = new String[] {line[0], line[1], line[2], line[3]};
        genders[i] = Byte.parseByte(line[4]);
        affections[i] = Byte.parseByte(line[5]);
        for (int j = 0; j < numMarkers; j++) {
          for (int k = 0; k < 2; k++) {
            allele = line[6 + 2 * j + k].charAt(0);
            if (ext.indexOfStr(allele + "", MISSING_VALUES) >= 0) {
              markerCounts[i][j] = -1;
            } else {
              genotypeCounts[j]++;
              if (alleles[j][0] == 0 || alleles[j][0] == allele) {
                alleles[j][0] = allele;
              } else if (alleles[j][1] == 0 || alleles[j][1] == allele) {
                alleles[j][1] = allele;
                markerCounts[i][j]++;
                freqs[j]++;
              } else {
                System.err.println("Error - more than 2 alleles for marker " + (i + 1) + "; first "
                                   + alleles[j][0] + " and " + alleles[j][1] + " and now "
                                   + allele);
              }
            }
          }
        }
      }

      for (int i = 0; i < numMarkers; i++) {
        freqs[i] /= genotypeCounts[i];
        if (freqs[i] > 0.50) { // if not minor allele then flip all information
          freqs[i] = 1 - freqs[i];
          allele = alleles[i][0];
          alleles[i][0] = alleles[i][1];
          alleles[i][1] = allele;
          for (int j = 0; j < count; j++) {
            switch (markerCounts[j][i]) {
              case 0:
                markerCounts[j][i] = 2;
                break;
              case 2:
                markerCounts[j][i] = 0;
                break;
              default:
                break;
            }
          }
        }
      }
      reader.close();

      famStruct = new FamilyStructure(ids, genders, affections);
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }
  }

  public char[][] getAlleles() {
    return alleles;
  }

  public FamilyStructure getFamilyStructure() {
    return famStruct;
  }

  public double[] getFrequencies() {
    return freqs;
  }


  public int[] getGenotypeCounts() {
    return genotypeCounts;
  }

  public int[][] getMarkerCounts() {
    return markerCounts;
  }
}
