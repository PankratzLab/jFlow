package org.genvisis.cnv.park;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.ext;

public class CNVtools {
  public static final String WINDOWS_DIRECTORY =
      "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\00src\\";

  // public static final String LINUX_DIRECTORY =
  // "/archive/parkinsons/gwas/Genotype and intensity data files";
  public static final String LINUX_DIRECTORY = "/home/npankrat/cnv/oldClusters/raw";

  // public static final String PLOT_ROOT = "plots/scratch/";
  public static final String PLOT_ROOT = "";

  public static final String ROOT_BASE = "gwas.split";

  public static final String DEFAULT_MARKER_LOOKUP = "lookup/";

  public static final String[] QC_HEADS = {"LRR_mean", "LRR_median", "LRR_SD", "BAF_mean",
                                           "BAF_median", "BAF_SD", "BAF_DRIFT", "WF", "GCWF"};

  public static final String[] ERRORS = {"large SD for LRR", "drifting BAF values",
                                         "waviness factor values", "Small-sized CNV calls"};

  public static void create(String map) {
    BufferedReader mapReader, reader;
    PrintWriter writer;
    String[] line;
    Hashtable<String, String> hash = new Hashtable<String, String>();
    Vector<String> v = new Vector<String>();
    String slot, subdir;
    boolean firstTime = true;
    File file;
    String[] ids = null;

    hash = createMarkerLookup(DEFAULT_MARKER_LOOKUP);

    try {
      mapReader = new BufferedReader(new FileReader(map));
      writer = new PrintWriter(new FileWriter("penn.cnv.xln"));

      while (mapReader.ready()) {
        line = mapReader.readLine().trim().split("[\\s]+");
        if ((map.endsWith(".map") && line.length != 4)
            || (map.endsWith(".bim") && line.length != 6)) {
          System.err.println("Error - it looks like your map file might not be in plink format, which could pose a problem");
        }
        slot = hash.get(line[1]);
        subdir = slot.substring(0, slot.length() - 1) + "/";
        file = new File(PLOT_ROOT + subdir + getFilename(line[1]));
        if (!file.exists()) {
          System.err.println("Error - '" + line[1] + "' was not found in the expected directory (./"
                             + subdir + "/)");
        } else {
          if (firstTime) {
            try {
              reader = new BufferedReader(new FileReader(file));
              reader.readLine();
              while (reader.ready()) {
                v.add(reader.readLine().trim().split("[\\s]+")[0]);
              }
              reader.close();
            } catch (FileNotFoundException fnfe) {
              System.err.println("Error: file \"" + file.getName()
                                 + "\" not found in current directory");
              System.exit(1);
            } catch (IOException ioe) {
              System.err.println("Error reading file \"" + file.getName() + "\"");
              System.exit(2);
            }
            ids = Array.toStringArray(v);
            writer.print("Name\tChr\tPosition");
            for (String id : ids) {
              writer.print("\t" + id + ".GType" + "\t" + id + ".Log R Ratio" + "\t" + id
                           + ".B Allele Freq");
            }
            writer.println();
            firstTime = false;
            System.out.println("There are " + ids.length
                               + " individuals (remember this in case you need to split up to the process in PennCNV)");
          }

          try {
            reader = new BufferedReader(new FileReader(file));
            reader.readLine();
            writer.print(line[1] + "\t" + line[0] + "\t" + line[3]);
            for (int i = 0; i < ids.length; i++) {
              line = reader.readLine().trim().split("[\\s]+");
              if (!line[0].equals(ids[i])) {
                System.err.println("Error - out of sync with first file (expecting '" + ids[i]
                                   + "'; found '" + line[0] + "')");
                System.exit(1);
              }
              switch (Integer.parseInt(line[11])) {
                case 0:
                  writer.print("\tAA");
                  break;
                case 1:
                  writer.print("\tAB");
                  break;
                case 2:
                  writer.print("\tBB");
                  break;
                case -1:
                  writer.print("\tNC");
                  break;
                default:
                  System.err.println("Error - unexpected B allele count");
                  break;
              }
              writer.print("\t" + line[9] + "\t" + line[8]);
            }
            writer.println();
            reader.close();
          } catch (FileNotFoundException fnfe) {
            System.err.println("Error: file \"" + file.getName()
                               + "\" not found in current directory");
            System.exit(1);
          } catch (IOException ioe) {
            System.err.println("Error reading file \"" + file.getName() + "\"");
            System.exit(2);
          }
        }
      }
      mapReader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + map + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + map + "\"");
      System.exit(2);
    }

  }

  public static Hashtable<String, String> createMarkerLookup(String dir) {
    BufferedReader reader;
    Hashtable<String, String> hash = new Hashtable<String, String>();

    File[] files = new File(dir).listFiles(new FilenameFilter() {
      @Override
      public boolean accept(File file, String filename) {
        return filename.endsWith(".txt");
      }
    });

    for (File file : files) {
      try {
        reader = new BufferedReader(new FileReader(file));
        while (reader.ready()) {
          hash.put(reader.readLine(), file.getName().substring(0, file.getName().lastIndexOf(".")));
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + file.getName() + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + file.getName() + "\"");
        System.exit(2);
      }
    }

    return hash;
  }

  public static String getFilename(String marker) {
    return ext.replaceAllWith(marker, ":", ".") + "_plots.xls";
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    // String filename = "SNCA.map";
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\gwas\\CNV\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\CNV\\allMarkers\\";
    // String dir = "C:\\Documents and Settings\\npankrat\\My
    // Documents\\CNV\\oldClusters\\";
    // String dir = "C:\\transfer\\forGawTrip\\CNVs\\final\\";
    String mapfile = "pd_gwas.bim";
    // String warnings = "all.out";
    // String warnings = "allMarkers.log";
    // String warnings = "parsed.log";
    String warnings = "final.log";

    String usage = "\\n" + "cnv.PennCNV requires 0-1 arguments\n" + "   (1) map file (i.e. map="
                   + mapfile + " (default))\n" + "   (2) errors and warnings capture (i.e. err="
                   + warnings + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("map=")) {
        mapfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("err=")) {
        warnings = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      create(mapfile);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
