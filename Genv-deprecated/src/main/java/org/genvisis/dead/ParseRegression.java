// -Xms1024M -Xmx1024M
package org.genvisis.dead;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;
import org.genvisis.filesys.*;

public class ParseRegression {
  // public static final String[] DEFAULT_DIRS = {"C:\\Documents and Settings\\npankrat\\My
  // Documents\\PD-singleton\\08_AOO_analyses\\01_no_covariates\\"};
  // public static final String[] DEFAULT_DIRS = {"C:\\Documents and Settings\\npankrat\\My
  // Documents\\gwas\\firstRun\\"};
  // public static final String[] DEFAULT_DIRS = {"C:\\Documents and Settings\\npankrat\\My
  // Documents\\gwas\\hits\\"};
  // public static final String[] DEFAULT_DIRS = {"C:\\Documents and Settings\\npankrat\\My
  // Documents\\gwas\\jewish\\"};
  // public static final String[] DEFAULT_DIRS = {"C:\\Documents and Settings\\npankrat\\My
  // Documents\\gwas\\results\\"};
  // public static final String[] DEFAULT_DIRS = {"C:\\Documents and Settings\\npankrat\\My
  // Documents\\gwas\\merged\\results\\unadjusted_aoo_results\\"};
  // public static final String[] DEFAULT_DIRS = {"C:\\Documents and Settings\\npankrat\\My
  // Documents\\gwas\\merged\\results\\JustIU\\"};
  // public static final String[] DEFAULT_DIRS = {"C:\\Documents and Settings\\npankrat\\My
  // Documents\\gwas\\merged\\results\\"};
  // public static final String[] DEFAULT_DIRS = {"C:\\Documents and Settings\\npankrat\\My
  // Documents\\gaw\\results\\allMarkers\\"};
  // public static final String[] DEFAULT_DIRS = {"C:\\Documents and Settings\\npankrat\\My
  // Documents\\gaw\\results\\goodMarkers\\"};
  // public static final String[] DEFAULT_DIRS = {"C:\\Documents and Settings\\npankrat\\My
  // Documents\\ADNI\\results\\"};
  // public static final String[] DEFAULT_DIRS = {"C:\\Documents and Settings\\npankrat\\My
  // Documents\\LOAD\\results\\"};
  public static final String[] DEFAULT_DIRS =
      {"C:\\Documents and Settings\\npankrat\\My Documents\\UMN\\Myron\\ExcisionPathway\\ancestry\\cpru\\"};

  // public static final String[] METHODS = {"linear", "logistic"};
  public static final String[] METHODS = {"logistic"};
  // public static final String[] METHODS = {"linear"};

  public static final String[] MODELS = {"Additive", "Dominant", "Recessive", "Genotypic"};
  // public static final String[] MODELS = {"AdditiveWith", "AdditiveWithout", "Dominant",
  // "Recessive"};
  // public static final String[] MODELS = {"Additive", "APOE_Risk", "APOEcounts",
  // "APOEcountsNoAgeNoSex"};


  public static final String[] MODEL_ABBREVS = {"ADD", "DOM", "REC", "GENO_2DF"};
  // public static final String[] MODEL_ABBREVS = {"ADD", "ADD", "DOM", "REC"};
  // public static final String[] MODEL_ABBREVS = {"ADD", "ADD", "ADD", "ADD"};

  public static final int[] MODEL_OFFSETS = {0, 0, 0, 0};

  public static final double CUTOFF = 0.0001;

  public static final String[] MAP_ROOTS = {"plink", "gwas"};

  public static void parseRegression(String[] dirs) {
    BufferedReader reader;
    PrintWriter writer, extras = null;
    String[] line;
    int count;
    Hashtable<String, String> alleles = new Hashtable<String, String>();
    Vector<String> snpNames, topSNPs;
    Hashtable<String, double[][]> allSNPs, used;
    double[][] data;
    SnpMarkerSet map;
    String[] markerNames;
    byte[] chrs;
    int[] positions;

    for (int meth = 0; meth < METHODS.length; meth++) {
      snpNames = new Vector<String>();
      alleles = new Hashtable<String, String>();
      allSNPs = new Hashtable<String, double[][]>();
      used = new Hashtable<String, double[][]>();
      topSNPs = new Vector<String>();
      for (int i = 0; i < dirs.length; i++) {
        for (int j = 0; j < MODELS.length; j++) {
          try {
            reader = new BufferedReader(new FileReader(dirs[i] + MODELS[j].toLowerCase() + ".assoc."
                                                       + METHODS[meth]));
            System.out.println("Parsing " + dirs[i] + MODELS[j].toLowerCase() + ".assoc."
                               + METHODS[meth]);
            reader.readLine();
            count = 0;
            while (reader.ready()) {
              line = reader.readLine().trim().split("[\\s]+");
              if (line.length == 1 && line[0].equals("25")) {
                System.err.println("Error - That weird thing with the last line of " + dirs[i]
                                   + MODELS[j].toLowerCase() + ".assoc." + METHODS[meth]);
                continue;
              }
              if (line[0].equals("23") || line[0].equals("24") || line[0].equals("26")) {
                continue;
              }
              // if
              // (line[4+MODEL_OFFSETS[j]].equals(MODEL_ABBREVS[j])
              // && Integer.parseInt(line[0])<23) {
              if (line[4 + MODEL_OFFSETS[j]].equals(MODEL_ABBREVS[j])
                  && Integer.parseInt(line[0]) < 26) {
                if (j == 0) {
                  snpNames.add(line[1]);
                  alleles.put(line[1], line[3]);
                  allSNPs.put(line[1], data = Matrix.doubleMatrix(MODELS.length, 3, -999));
                } else {
                  if (!line[1].equals(snpNames.elementAt(count))) {
                    System.err.println("Error - out of sync in " + dirs[i] + MODELS[j].toLowerCase()
                                       + ".assoc." + METHODS[meth] + " - expecting "
                                       + snpNames.elementAt(count) + ", found " + line[1]);
                    System.exit(1);
                  }
                  data = allSNPs.get(line[1]);
                }
                data[j][0] =
                    line[6
                         + MODEL_OFFSETS[j]].equals("NA") ? -999
                                                          : Double.parseDouble(line[6
                                                                                    + MODEL_OFFSETS[j]]);
                data[j][1] =
                    line[7
                         + MODEL_OFFSETS[j]].equals("NA") ? -999
                                                          : Double.parseDouble(line[7
                                                                                    + MODEL_OFFSETS[j]]);
                data[j][2] =
                    line[8
                         + MODEL_OFFSETS[j]].equals("NA") ? -999
                                                          : Double.parseDouble(line[8
                                                                                    + MODEL_OFFSETS[j]]);
                if (data[j][2] != -999 && data[j][2] < CUTOFF) {
                  if (!topSNPs.contains(line[1])) {
                    topSNPs.add(line[1]);
                  }
                }
                count++;
              }
            }
            reader.close();
          } catch (FileNotFoundException fnfe) {
            System.err.println("Error: file \"" + dirs[i] + MODELS[j].toLowerCase() + ".assoc."
                               + METHODS[meth] + "\" not found in current directory");
            System.exit(1);
          } catch (IOException e) {
            System.err.println("Error reading file \"" + dirs[i] + MODELS[j].toLowerCase()
                               + ".assoc." + METHODS[meth] + "\"");
            e.printStackTrace();
            System.exit(2);
          }
        }
        System.out.println("Writing table...");

        map = null;
        for (int j = 0; j < MAP_ROOTS.length && map == null; j++) {
          if (new File(dirs[i] + MAP_ROOTS[j] + ".map").exists()) {
            map = new SnpMarkerSet(dirs[i] + MAP_ROOTS[j] + ".map");
          } else if (new File(dirs[i] + MAP_ROOTS[j] + ".bim").exists()) {
            map = new SnpMarkerSet(dirs[i] + MAP_ROOTS[j] + ".bim");
          }
        }
        if (map == null) {
          System.err.println("Error - could not find a .MAP or .BIM file with any of the following roots: "
                             + Array.toStr(MAP_ROOTS, ", "));
          System.exit(1);
        }
        markerNames = map.getMarkerNames();
        chrs = map.getChrs();
        positions = map.getPositions();

        try {
          writer = new PrintWriter(new FileWriter(dirs[i] + METHODS[meth] + ".xls"));
          writer.println("SNP\tChr\tPosition\tAllele\t"
                         + Array.toStr(MODELS, "\tstatistic\tp-value\t") + "\tstatistic\tp-value");
          for (int j = 0; j < markerNames.length; j++) {
            if (allSNPs.containsKey(markerNames[j])) {
              writer.print(markerNames[j] + "\t" + chrs[j] + "\t" + positions[j] + "\t"
                           + alleles.get(markerNames[j]));
              used.put(markerNames[j], data = allSNPs.remove(markerNames[j]));
              for (int k = 0; k < MODELS.length; k++) {
                writer.print("\t" + (data[k][0] == -999 ? "." : ext.formDeci(data[k][0], 3, true))
                             + "\t" + (data[k][1] == -999 ? "." : ext.formDeci(data[k][1], 4, true))
                             + "\t"
                             + (data[k][2] == -999 ? "." : ext.prettyP(data[k][2], 2, 5, 4, true)));
              }

              writer.println();
            } else {
              if (extras == null) {
                extras =
                    new PrintWriter(new FileWriter(dirs[i]
                                                   + "markers in map file but without results.out"));
              }
              extras.println(markerNames[j] + "\t" + chrs[j] + "\t" + positions[j]);
            }
          }
          if (extras != null) {
            extras.close();
          }
          if (allSNPs.size() > 0) {
            extras =
                new PrintWriter(new FileWriter(dirs[i]
                                               + "marke with results file but not in map file.out"));
            line = HashVec.getKeys(allSNPs);
            for (int k = 0; k < line.length; k++) {
              extras.println(line[k]);
            }
            extras.close();
          }
          writer.close();
        } catch (Exception e) {
          System.err.println("Error writing to " + dirs[i] + METHODS[meth] + ".xls");
          e.printStackTrace();
        }

        System.out.println("Writing top hits...");
        try {
          writer = new PrintWriter(new FileWriter(dirs[i] + METHODS[meth] + "_topHits.xls"));
          writer.println("SNP\t" + Array.toStr(MODELS, "\tp-value\t") + "\tp-value");
          for (int j = 0; j < topSNPs.size(); j++) {
            writer.print(topSNPs.elementAt(j));
            data = used.get(topSNPs.elementAt(j));
            for (int k = 0; k < MODELS.length; k++) {
              writer.print("\t" + (data[k][0] == -999 ? "." : ext.formDeci(data[k][0], 3, true))
                           + "\t" + (data[k][2] == -999 ? "." : data[k][2]));
            }
            writer.println();
          }
          writer.close();
        } catch (FileNotFoundException fnfe) {
          System.err.println("Error: file \"" + dirs[i] + METHODS[meth] + "_topHits.xls"
                             + "\" is currently in use");
          System.exit(1);
        } catch (IOException ioe) {
          System.err.println("Error writing to file \"" + dirs[i] + METHODS[meth] + "_topHits.xls"
                             + "\"");
          System.exit(2);
        }
      }
    }
  }

  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    Vector<String> v = new Vector<String>();

    String usage = "\n" + "park.parseLinear requires 0-1 arguments\n"
                   + "   (1) directory (i.e. relevantDirectory/ (default))\n" + "";

    for (int i = 0; i < args.length; i++) {
      if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h")
          || args[i].equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else {
        if (new File(args[i]).exists() && new File(args[i]).isDirectory()) {
          if (!args[i].endsWith("/") && !args[i].endsWith("\\")) {
            args[i] += "/";
          }
          v.add(args[i]);
          numArgs--;
        } else {
          System.err.println("'" + args[i] + "' is not a valid directory");
        }
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      parseRegression(v.size() == 0 ? DEFAULT_DIRS : Array.toStringArray(v));
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
