// -Xms1024M -Xmx1024M
package org.genvisis.parse;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.FamilyStructure;
import org.genvisis.filesys.SnpMarkerSet;

public class SNPlist {
  public static final String[][] FIELDS = {{"SAMPLE_ID"}, {"ASSAY_ID"}, {"GENOTYPE_ID"}};

  public static void main(String[] args) {
    int numArgs = args.length;
    String dir =
        "C:\\Documents and Settings\\npankrat\\My Documents\\UMN\\Myron\\Indian_Diabetes\\00src\\";
    String snps = "W1_W2.dat";
    String pedigree = "secondPedigree.dat";

    String usage = "\n" + "parse.SNPlist requires 0-1 arguments\n" + "   (1) dir (i.e. dir=" + dir
                   + " (default))\n" + "   (2) SNPs filename (i.e. snps=" + snps + " (default))\n"
                   + "   (2) pedigree (with DNA in column 7) filename (i.e. ped=" + pedigree
                   + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("snps=")) {
        snps = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("ped=")) {
        pedigree = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      parse(dir, snps, pedigree, new Logger());
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void parse(String dir, String filename, String pedigree, Logger log) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String trav;
    Hashtable<String, Hashtable<String, String>> hashes;
    Hashtable<String, String> hash;
    Vector<String> snps;
    int[] indices;
    SnpMarkerSet markers;
    String[] markerNames, inds;
    FamilyStructure famStruct;
    int missingInds, missingGenotypes;
    IntVector iv;

    snps = new Vector<String>();
    hashes = new Hashtable<String, Hashtable<String, String>>();
    try {
      reader = new BufferedReader(new FileReader(dir + filename));
      line = reader.readLine().trim().split("\t");
      indices = ext.indexFactors(FIELDS, line, false, true, true, true);
      while (reader.ready()) {
        line = reader.readLine().trim().split("\t", -1);
        HashVec.addToHashHash(hashes, line[indices[0]], line[indices[1]], line[indices[2]]);
        HashVec.addIfAbsent(line[indices[1]], snps);
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + filename + "\"");
      System.exit(2);
    }

    markers = new SnpMarkerSet(Array.toStringArray(snps));
    markers.parseSNPlocations();
    markers.sortMarkers();
    markers.writeToFile(dir + ext.rootOf(filename) + ".map", SnpMarkerSet.PLINK_MAP_FORMAT, log);

    markerNames = markers.getMarkerNames();
    iv = new IntVector();
    famStruct = new FamilyStructure(dir + pedigree, true);
    inds = famStruct.getDnas();
    missingInds = missingGenotypes = 0;
    try {
      writer = new PrintWriter(new FileWriter(dir + ext.rootOf(filename) + ".ped"));
      for (int i = 0; i < inds.length; i++) {
        hash = hashes.get(inds[i]);
        writer.print(famStruct.getIndividualHeader(i, false));
        if (hash == null) {
          System.err.println("Error - '" + inds[i]
                             + "' was listed in the pedigree as a dna, but was not found among the list of genotypes");
          for (String markerName : markerNames) {
            writer.print("\t0\t0");
          }
          missingInds++;
        } else {
          for (String markerName : markerNames) {
            trav = hash.get(markerName);
            if (trav == null) {
              System.err.println("Error - sample '" + inds[i] + "' no genotypes for " + markerName);
              writer.print("\t0\t0");
              missingGenotypes++;
              iv.addIfAbsent(i);
            } else if (trav.equals("null")) {
              writer.print("\t0\t0");
            } else {
              writer.print("\t" + trav.charAt(0) + "\t" + trav.charAt(trav.length() == 1 ? 0 : 1));
            }
          }
        }
        writer.println();
      }
      System.out.println("There were " + (missingInds > 0 ? missingInds : "no") + " indiviudal"
                         + (missingInds == 1 ? "" : "s")
                         + " in the pedigree file that were completely missing data");
      System.out.println("There were " + (missingInds > 0 ? "an additional " : "")
                         + (iv.size() > 0 ? iv.size() : "no") + " indiviudal"
                         + (iv.size() == 1 ? "" : "s") + " with DNA missing "
                         + (missingGenotypes > 0 ? "a total of " + missingGenotypes : "any")
                         + " genotypes");
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + dir + ext.rootOf(filename) + ".ped");
      e.printStackTrace();
    }


  }
}
