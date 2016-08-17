// -Xms1024M -Xmx1024M
package org.genvisis.parse;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.SnpMarkerSet;

public class NIEHS_SNPs {
  public static final String DIR =
      "C:\\Documents and Settings\\npankrat\\My Documents\\UMN\\Myron\\ExcisionPathway\\expression\\NIEHS\\";
  public static final String VARIANT_FILE = "EGP.prettybase.txt";
  public static final String RS_LOOKUP_FILE = "rsEGP.txt";
  public static final String SAMPLE_LOOKUP_FILE = "SampleLookup.xln";

  public static void main(String[] args) {
    try {
      parse();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void parse() {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, genes;
    Hashtable<String, String> hash, sampleLookup;
    Hashtable<String, Vector<String>> geneSNPs;
    Vector<String> v;
    int count;
    String[] markerNames, annotation;
    byte[] chrs;
    int[] positions;
    SnpMarkerSet map;
    String[][] alleles;
    // short[][] alleleCounts;
    byte[][][] data;
    String[] samples;
    int indIndex, snpIndex;
    String trav, name;
    Logger log;

    if (!new File(DIR + ext.rootOf(RS_LOOKUP_FILE) + ".ser").exists()) {
      log = new Logger(DIR + ext.rootOf(RS_LOOKUP_FILE) + ".log");
      hash = new Hashtable<String, String>();
      count = 0;
      try {
        reader = new BufferedReader(new FileReader(DIR + VARIANT_FILE));
        while (reader.ready()) {
          line = reader.readLine().trim().split("[\\s]+");
          HashVec.addToHashIfAbsent(hash, line[0], "");
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        log.reportError("Error: file \"" + DIR + VARIANT_FILE
                        + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        log.reportError("Error reading file \"" + DIR + VARIANT_FILE + "\"");
        log.reportException(ioe);
        System.exit(2);
      }

      try {
        reader = new BufferedReader(new FileReader(DIR + RS_LOOKUP_FILE));
        while (reader.ready()) {
          line = reader.readLine().trim().split("[\\s]+");
          if (!line[0].startsWith("#")) {
            if (!line[2].startsWith("chr")) {
              log.reportError("Error - snp '" + line[5] + "' doesn't have a valid chromosome: "
                              + line[2]);
            }
            trav = line[4] + "-" + line[1].split("-")[0] + "-" + line[2].substring(3);
            if (hash.containsKey(trav)) {
              name = hash.get(trav);
              if (name.equals("")) {
                hash.put(trav, line[5]);
              } else if (!name.equals(line[5])) {
                log.reportError("Error - mismatch for '" + trav + "' (" + name + " and " + line[5]
                                + ")");
              }
            }
          }
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        log.reportError("Error: file \"" + DIR + RS_LOOKUP_FILE
                        + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        log.reportError("Error reading file \"" + DIR + RS_LOOKUP_FILE + "\"");
        log.reportException(ioe);
        System.exit(2);
      }

      annotation = HashVec.getKeys(hash);
      chrs = new byte[annotation.length];
      positions = new int[annotation.length];
      markerNames = new String[annotation.length];
      geneSNPs = new Hashtable<String, Vector<String>>();
      for (int i = 0; i < annotation.length; i++) {
        line = annotation[i].split("-");
        markerNames[i] = hash.get(annotation[i]);
        if (markerNames[i].equals("")) {
          log.reportError("Error - no record for " + annotation[i] + " in " + RS_LOOKUP_FILE);
          markerNames[i] = "chr" + line[2] + "_" + line[0];
        } else if (markerNames[i].equals("notAvailable")) {
          markerNames[i] = "chr" + line[2] + "_" + line[0];
        }
        chrs[i] = Positions.chromosomeNumber(line[2]);
        positions[i] = Integer.parseInt(line[0]);
        HashVec.addToHashVec(geneSNPs, line[1], markerNames[count], false);
      }
      new SnpMarkerSet(markerNames, chrs, positions, null, Array.toMatrix(annotation), true,
                       false).serialize(DIR + ext.rootOf(RS_LOOKUP_FILE) + ".ser");

      genes = HashVec.getKeys(geneSNPs);
      new File(DIR + "genes/").mkdir();
      for (String gene : genes) {
        v = geneSNPs.get(gene);
        try {
          writer = new PrintWriter(new FileWriter(DIR + "genes/" + gene + ".txt"));
          for (int j = 0; j < v.size(); j++) {
            writer.println(v.elementAt(j));
          }
          writer.close();
        } catch (Exception e) {
          log.reportError("Error writing to " + DIR + "genes/" + gene + ".txt");
          e.printStackTrace();
        }
      }
    }

    log = new Logger(DIR + ext.rootOf(VARIANT_FILE) + ".log");
    map = SnpMarkerSet.load(DIR + ext.rootOf(RS_LOOKUP_FILE) + ".ser", false);
    annotation = Matrix.extractColumn(map.getAnnotation(), 0);
    hash = new Hashtable<String, String>();
    for (int i = 0; i < annotation.length; i++) {
      if (hash.containsKey(annotation[i])) {
        log.reportError("Error - multiple instances of " + annotation[i]);
      }
      hash.put(annotation[i], i + "");
    }

    sampleLookup = HashVec.loadFileToHashString(DIR + SAMPLE_LOOKUP_FILE, true);
    samples = HashVec.getKeys(sampleLookup);
    for (int i = 0; i < samples.length; i++) {
      hash.put(samples[i], i + "");
    }

    data = Matrix.tripleArrays(samples.length, annotation.length, 2, (byte) -1);
    alleles = new String[annotation.length][2];
    // alleleCounts = new short[annotation.length][2];
    try {
      reader = new BufferedReader(new FileReader(DIR + VARIANT_FILE));
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");

        if (hash.containsKey(line[1])) {
          indIndex = Integer.parseInt(hash.get(line[1]));
          // System.err.println(line[1]+"\t"+indIndex);
          snpIndex = Integer.parseInt(hash.get(line[0]));
          // System.err.println(line[0]+"\t"+snpIndex);
          if (data[indIndex][snpIndex][0] == -1) {
            for (int i = 0; i < 2; i++) {
              if (line[2 + i].equals("N") || line[2 + i].equals("-")) {
                data[indIndex][snpIndex][i] = 0;
              } else {
                for (int j = 0; j < 2; j++) {
                  if (alleles[snpIndex][j] == null || line[2 + i].equals(alleles[snpIndex][j])) {
                    alleles[snpIndex][j] = line[2 + i];
                    // alleleCounts[snpIndex][j]++;
                    data[indIndex][snpIndex][i] = (byte) (j + 1);
                    j = 2;
                  }
                }
                if (data[indIndex][snpIndex][i] == -1) {
                  data[indIndex][snpIndex][i] = 0;
                  log.reportError(Array.toStr(annotation[snpIndex].split("-")) + "\t" + line[2 + i]
                                  + "\tinstead of either\t" + alleles[snpIndex][0] + "\t"
                                  + alleles[snpIndex][1]);
                }
              }
            }
          } else {
            System.err.println("Error - duplicate attempts to fill " + annotation[snpIndex]
                               + " for sample " + samples[indIndex]);
          }
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + DIR + VARIANT_FILE + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + DIR + VARIANT_FILE + "\"");
      log.reportException(ioe);
      System.exit(2);
    }

    map.writeToFile(DIR + "plink.map", SnpMarkerSet.PLINK_MAP_FORMAT, log);
    try {
      writer = new PrintWriter(new FileWriter(DIR + "plink.ped"));
      for (int i = 0; i < samples.length; i++) {
        writer.print(samples[i] + "\t" + sampleLookup.get(samples[i]).substring(0, 7)
                     + "\t0\t0\t1\t1");
        for (int j = 0; j < annotation.length; j++) {
          if (data[i][j][0] < 1 || data[i][j][1] < 1) {
            writer.print("\t0\t0");
          } else {
            writer.print("\t" + data[i][j][0] + "\t" + data[i][j][1]);
          }
        }
        writer.println();
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + DIR + "plink.ped");
      e.printStackTrace();
    }

  }

}
