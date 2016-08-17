package org.genvisis.link.bat;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.ext;
import org.genvisis.filesys.SnpMarkerSet;
import org.genvisis.link.LinkageFormat;
import org.genvisis.link.LinkageMap;

public class Mendel {
  public static final String[] THETAS =
      {"0.00", "0.01", "0.05", "0.10", "0.15", "0.20", "0.30", "0.40", "0.50"};

  public static void createFiles(String dir, String pedfile, String mapfile, boolean xLinked) {
    BufferedReader reader;
    PrintWriter writer;
    String[] freq, pens;
    LinkageMap map;
    String[] markerNames;
    double[] cM_distances;

    try {
      reader = new BufferedReader(new FileReader(dir + mapfile));
      for (int i = 0; i < 4; i++) {
        reader.readLine();
      }
      freq = reader.readLine().trim().split("[\\s]+");
      reader.readLine();
      pens = reader.readLine().trim().split("[\\s]+");
      reader.close();

      writer = new PrintWriter(new FileWriter(dir + "mendel.def"));
      // if (mapfile.indexOf("23")!=-1) {
      if (xLinked) {
        writer.println("trait,X-Linked,2,2");
      } else {
        writer.println("trait,autosome,2,2");
      }
      writer.println("A," + Double.parseDouble(freq[0]));
      writer.println("a," + Double.parseDouble(freq[1]));
      writer.println("1");
      writer.println("2");

      map = new LinkageMap(dir + mapfile);
      markerNames = map.getMarkerNames();
      cM_distances = map.getDistances(false);

      for (String markerName : markerNames) {
        writer.println(markerName);
      }
      writer.close();

      writer = new PrintWriter(new FileWriter(dir + "mendel.map"));
      writer.println("trait");
      writer.println();
      for (int i = 0; i < markerNames.length; i++) {
        writer.println(markerNames[i]);
        if (i < markerNames.length - 1) {
          writer.println(" " + cM_distances[i + 1]);
        }
      }
      writer.close();

      writer = new PrintWriter(new FileWriter(dir + "Control.in"));
      writer.println("!input files");
      writer.println("DEFINITION_FILE = mendel.def");
      writer.println("MAP_FILE = mendel.map");
      writer.println("PEDIGREE_FILE = " + pedfile);
      writer.println("PEDIGREE_LINKAGE_FORMAT = True");
      writer.println("AFFECTED_LOCUS_OR_FACTOR=trait");
      writer.println("! output files");
      writer.println("PENETRANCE = " + Double.parseDouble(pens[0]) + " :: A\\A");
      writer.println("PENETRANCE = " + Double.parseDouble(pens[1]) + " :: A\\a");
      writer.println("PENETRANCE = " + Double.parseDouble(pens[2]) + " :: a\\a");
      writer.println("SUMMARY_FILE = mendel.sum");
      writer.println("OUTPUT_FILE = mendel.out");
      writer.println("! analysis specific settings");
      writer.println("ANALYSIS_OPTION = Location_scores");
      writer.println("NUMBER_OF_MARKERS_INCLUDED = 1");
      writer.println("TRAVEL=GRID");
      writer.println("STANDARD_GRID=True");
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + mapfile + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + mapfile + "\"");
      System.exit(2);
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String pedfile = "re_chrom16.pre";
    String mapfile = "map16.dat";
    boolean xLinked = false;

    String usage =
        "\\n" + "park.Mendel requires 0-1 arguments\n" + "   (1) pedigree file (i.e. ped=" + pedfile
            + " (default))\n" + "   (2) map file (i.e. map=" + mapfile + " (default))\n"
            + "   (3) all markers are on the X chromosome (i.e. xLinked=" + xLinked
            + " (default))\n" + "  OR\n" + "   (1) chromosome number (i.e. chr=1)\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("ped=")) {
        pedfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("map=")) {
        mapfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("xLinked=")) {
        xLinked = arg.split("=")[1].toLowerCase().equals("true");
        numArgs--;
      } else if (arg.startsWith("chr=")) {
        int chr = Integer.parseInt(arg.split("=")[1]);
        pedfile = "re_chrom" + ext.chrome(chr) + ".pre";
        mapfile = "map" + ext.chrome(chr) + ".dat";
        if (chr == 23) {
          xLinked = true;
        }
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      createFiles("", pedfile, mapfile, xLinked);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static String[][] parseMaxLods(String filename) {
    BufferedReader reader;
    String temp;
    Vector<String[]> v = new Vector<String[]>();
    String[][] results;
    boolean done = false;
    double lod, max;
    int theta;

    try {
      reader = new BufferedReader(new FileReader(filename));
      // for (int i = 0; i < 5; i++) { // version 8.01
      for (int i = 0; i < 7; i++) { // version 9.00
        reader.readLine();
      }
      while (!done) {
        temp = reader.readLine();
        if (temp == null) {
          System.err.println(
              "Error: no results found in mendel.sum; check doublecheck.out to make sure there aren't more errors");
          System.exit(1);
        }
        if (temp.trim().split("[\\s]+").length > 1) {
          max = Double.NEGATIVE_INFINITY;
          theta = -1;
          for (int i = 0; i < 7; i++) {
            lod = Double.parseDouble(temp.substring(17 + i * 7, 17 + i * 7 + 7));
            if (lod > max) {
              max = lod;
              theta = i;
            }
          }
          v.add(new String[] {temp.trim().split("[\\s]+")[0], max + "", THETAS[theta]});
        } else {
          done = true;
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }

    results = new String[v.size()][];
    for (int i = 0; i < v.size(); i++) {
      results[i] = v.elementAt(i);
    }

    return results;
  }

  public static String[][] runModel(String dir, String pedfile, String mapfile,
      SnpMarkerSet markerSet, double disease, double pp, double pq, double qq) {
    LinkageMap map;
    int chr;
    String[][] results, allResults;
    Hashtable<String, String> chrHash;
    Vector<String> autosomalMarkers, xLinkedMarkers;
    String[] markerNames;
    String trav;

    map = new LinkageMap(dir + mapfile);
    markerNames = map.getMarkerNames();

    autosomalMarkers = new Vector<String>();
    xLinkedMarkers = new Vector<String>();
    if (markerSet != null) {
      chrHash = markerSet.getChrHash();
      for (String markerName : markerNames) {
        trav = chrHash.get(markerName);
        if (trav == null) {
          System.err.println("Error - marker '" + markerName
              + "' was not present in the SnpMarkerSet, assuming it is autosomal");
          autosomalMarkers.add(markerName);
        } else {
          chr = Integer.parseInt(trav.split("[\\s]+")[0]);
          if (chr < 23) {
            autosomalMarkers.add(markerName);
          } else if (chr == 23) {
            xLinkedMarkers.add(markerName);
          } else {
            System.err.println("Error - the chromosome listed for '" + markerName + "' (" + chr
                + ") is invalid, marker will not be analyzed");
          }
        }
      }
    } else {
      for (String markerName : markerNames) {
        autosomalMarkers.add(markerName);
      }
    }

    allResults = new String[autosomalMarkers.size() + xLinkedMarkers.size()][];

    if (autosomalMarkers.size() > 0) {
      map = new LinkageMap(dir + mapfile);
      map.alterPenetrance(dir + "autosomal.dat", disease, pp, pq, qq, false);
      LinkageFormat.filterMarkers(dir + pedfile, dir + "autosomal.pre", dir + "autosomal.dat",
          dir + "autosomal.dat", Array.toStringArray(autosomalMarkers), null);
      Mendel.createFiles(dir, "autosomal.pre", "autosomal.dat", false);
      CmdLine.run("mendel", dir);
      results = Mendel.parseMaxLods(dir + "mendel.sum");
      for (int i = 0; i < results.length; i++) {
        allResults[i] = results[i];
      }
    }

    if (xLinkedMarkers.size() > 0) {
      map = new LinkageMap(dir + "map.dat");
      map.setChr(23);
      map.alterPenetrance(dir + "map23.dat", disease, pp, pq, qq, false);
      LinkageFormat.filterMarkers(dir + pedfile, dir + "chr23.pre", dir + "map23.dat",
          dir + "map23.dat", Array.toStringArray(xLinkedMarkers), null);
      Mendel.createFiles(dir, "chr23.pre", "map23.dat", true);
      CmdLine.run("mendel", dir);
      results = Mendel.parseMaxLods(dir + "mendel.sum");
      for (int i = 0; i < results.length; i++) {
        allResults[autosomalMarkers.size() + i] = results[i];
      }
    }

    return allResults;
  }

}
