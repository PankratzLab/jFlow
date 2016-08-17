package org.genvisis.kaput;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;

import com.google.common.primitives.Doubles;

public class dbExport {
  public static boolean INCLUDE_MISSING = false;

  public static boolean SOLAR = false;

  public static String DEFAULT_FILE = "struct.dat";

  public static String DEFAULT_TRAIT = "AOO";

  public static String[] ALT_DIRS =
      {"C:\\Documents and Settings\\npankrat\\My Documents\\1_CRFdb\\",
       "/home/npankrat/park/00masters/crf_db/", "/home/npankrat/reed/CRFdb/"};

  public static String DEFAULT_DB = "crf_db.dat";

  public dbExport(String db_file, String filename, Vector<String> traits,
                  boolean uniqueIDs) throws IOException {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] allTraits, line, data;
    String varName, trav;
    Hashtable<String, String[]> hash = new Hashtable<String, String[]>();
    int[] indices = new int[traits.size()], keys;
    Vector<String> inds = new Vector<String>();
    boolean missing;
    DoubleVector dv;
    double[] dist;
    reader = Files.getReader(db_file, ALT_DIRS);

    allTraits = reader.readLine().split("\t", -1);
    for (int i = 0; i < allTraits.length; i++) {
      if (traits.contains(allTraits[i].toUpperCase())) {
        indices[traits.indexOf(allTraits[i].toUpperCase())] = i;
      }
    }

    varName = traits.firstElement().toUpperCase();
    if (traits.size() > 3) {
      varName += "_" + (traits.size() - 1) + "covars";
    } else {
      for (int i = 1; i < traits.size(); i++) {
        varName += "_" + traits.elementAt(i).toUpperCase();
      }
    }

    while (reader.ready()) {
      line = reader.readLine().split("\t", -1);
      if (line.length != allTraits.length) {
        System.err.println("Error - problem with crf_db.dat file; not enough columns for line starting with '"
                           + line[0] + "'");
        System.exit(5);
      }
      data = new String[traits.size()];
      missing = false;
      for (int i = 0; i < traits.size(); i++) {
        data[i] = line[indices[i]];
        if (data[i].equals(".")) {
          missing = true;
        }
      }
      if (INCLUDE_MISSING || !missing) {
        hash.put(line[0], data);
        inds.add(line[0]);
      }
    }
    reader.close();

    writer = new PrintWriter(new FileWriter(varName + (SOLAR ? ".dat" : ".xls")));
    if (SOLAR) {
      writer.print("FAMID,ID," + varName);
      for (int i = 1; i < traits.size(); i++) {
        writer.print("," + traits.elementAt(i).toUpperCase());
      }
      writer.println();
    }

    dv = new DoubleVector();
    if (filename.equals("")) {
      keys = Sort.quicksort(Array.toStringArray(inds));
      for (int i = 0; i < inds.size(); i++) {
        trav = inds.elementAt(keys[i]);
        writer.println(formatData(uniqueIDs ? trav : trav.substring(0, 5),
                                  uniqueIDs ? trav
                                            : Integer.valueOf(trav.substring(5)).intValue() + "",
                                  hash.get(trav)));
        if (hash.containsKey(trav)) {
          dv.add(Double.parseDouble(hash.get(trav)[0]));
        }
      }
    } else {
      reader = new BufferedReader(new FileReader(filename));
      reader.mark(500);
      if (!reader.readLine().split("[\\s\\,]+")[0].toLowerCase().startsWith("fam")) {
        reader.reset();
      }
      while (reader.ready()) {
        line = reader.readLine().split("[\\s\\,]+");
        trav = uniqueIDs ? line[1] : line[0] + ext.formNum(line[1], 3);
        writer.println(formatData(line[0], line[1],
                                  hash.containsKey(trav) ? hash.get(trav)
                                                         : new String[traits.size()]));
        if (hash.containsKey(trav)) {
          dv.add(Double.parseDouble((hash.get(trav))[0]));
        }
      }
      reader.close();
    }
    writer.close();

    dist = Doubles.toArray(dv);
    writer = new PrintWriter(new FileWriter("punkd.txt")); // why is this
    // more than in
    // the file?????
    writer.println(Array.toStr(dist, 5, 5, "\n"));
    writer.close();

    if (Math.abs(Array.kurtosis(dist)) > 1) {
      System.err.println("Warning!  Kurtosis is " + ext.formDeci(Array.kurtosis(dist), 3)
                         + " which is too high.");
      System.err.println("See if the kurtosis looks any better for the following transformations...");
      // for (int k = 0; k<Transformations.NUM_TRANSFORMATIONS; k++) {
      // System.err.println(Transformations.getLabel(k)+":
      // "+ext.formDeci(Stats.kurtosis(Transformations.transform(dist,
      // k)), 4));
      // }
      System.err.println();
      // if (Array.stdev(bc.getTransform_MaxLL()) < 0.5) {
      // System.err.println("Estimated Trait Standard Deviation is below
      // 0.5");
      // System.err.println(" multiplying transformed trait by a factor of
      // 3.7");
      // mFactor = 3.7;
      // }
      writer = new PrintWriter(new FileWriter(varName + "_BOXCOX" + (SOLAR ? ".dat" : ".xls")));
      if (SOLAR) {
        writer.print("FAMID,ID," + varName + "_BOXCOX");
        for (int i = 1; i < traits.size(); i++) {
          writer.print("," + traits.elementAt(i).toUpperCase());
        }
        writer.println();
      }

      if (filename.equals("")) {
        keys = Sort.quicksort(Array.toStringArray(inds));
        for (int i = 0; i < inds.size(); i++) {
          trav = inds.elementAt(keys[i]);
          data = hash.get(trav);
          // data[0] =
          // ""+(bc.lookUpValue_MaxLL(Double.parseDouble(data[0]))*mFactor);
          writer.println(formatData(uniqueIDs ? trav : trav.substring(0, 5),
                                    uniqueIDs ? trav
                                              : Integer.valueOf(trav.substring(5)).intValue() + "",
                                    data));
        }
      } else {
        reader = new BufferedReader(new FileReader(filename));
        reader.mark(500);
        if (!reader.readLine().split("[\\s\\,]+")[0].toLowerCase().startsWith("fam")) {
          reader.reset();
        }
        while (reader.ready()) {
          line = reader.readLine().split("[\\s\\,]+");
          trav = uniqueIDs ? line[1] : line[0] + ext.formNum(line[1], 3);
          if (hash.containsKey(trav)) {
            data = hash.get(trav);
            // data[0] =
            // ext.formDeci(bc.lookUpValue_MaxLL(Double.parseDouble(data[0]))*mFactor,
            // 5);
          } else {
            data = new String[traits.size()];
          }
          writer.println(formatData(line[0], line[1], data));
        }
        reader.close();
      }
      writer.close();
    }
  }

  public String formatData(String famid, String indid, String[] data) {
    String delimit = SOLAR ? "," : "\t";
    String output = "";

    for (int i = 0; i < data.length; i++) {
      output +=
          delimit + (data[i] == null || data[i].equals(".") || data[i].equals("") ? SOLAR ? "" : "."
                                                                                  : data[i]);
    }

    return famid + delimit + indid + output;
  }

  public static void main(String[] args) throws IOException {
    String db_file = DEFAULT_DB;
    String filename = "default";
    Vector<String> traits = new Vector<String>();
    String[] line;
    String temp, firstBlood = "";
    boolean uniqueIDs = false;

    DEFAULT_TRAIT = "CRP";
    db_file = "database.txt";
    uniqueIDs = true;

    String usage = "\n" + "park.dbExport requires 1-4 arguments\n"
                   + "   (1) the database file (default: db=" + DEFAULT_DB + ")\n"
                   + "   (2) either a struct file (default: file=" + DEFAULT_FILE + ") or -all\n"
                   + "   (3) match based on unique IndID (i.e. -uniqueids (not the default))\n"
                   + "   (4+) the traits to be extracted (default: " + DEFAULT_TRAIT + ")\n"
                   + "   (optional) -solar (will create comma delimited file with header)\n"
                   + "   (optional) -miss (removes data for individuals with any missing values)\n"
                   + "";

    for (int i = 0; i < args.length; i++) {
      if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h")
          || args[i].equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (args[i].startsWith("-uniqueids")) {
        uniqueIDs = true;
      } else if (args[i].startsWith("-all")) {
        filename = "";
      } else if (args[i].startsWith("-solar")) {
        SOLAR = true;
      } else if (args[i].startsWith("-miss")) {
        INCLUDE_MISSING = true;
      } else if (args[i].startsWith("file=")) {
        filename = args[i].split("=")[1];
        if (!new File(filename).exists()) {
          System.err.println("Error - file '" + filename + "' does not exist");
          System.exit(2);
        }
      } else if (args[i].startsWith("db=")) {
        db_file = args[i].split("=")[1];
      } else {
        if (traits.size() == 0) {
          firstBlood = args[i].toUpperCase();
        }
        traits.add(args[i].toUpperCase());
      }
    }

    if (traits.size() == 0) {
      System.out.println("No trait specified; using default (" + DEFAULT_TRAIT + ")");
      traits.add(DEFAULT_TRAIT);
      firstBlood = DEFAULT_TRAIT;
    }
    traits.add("\n");
    line = Files.getReader(db_file, ALT_DIRS).readLine().split("[\\s]+");
    for (int i = 0; i < line.length; i++) {
      if (traits.contains(line[i].toUpperCase())) {
        traits.add(traits.remove(traits.indexOf(line[i].toUpperCase())));
      }
    }
    do {
      temp = traits.remove(0);
      if (!temp.equals("\n")) {
        System.err.println("Trait '" + temp + "' was not found in the header of " + db_file + "\n");
        System.err.println(usage);
        System.exit(4);
      }
    } while (!temp.equals("\n"));

    if (!traits.contains(firstBlood)) {
      System.err.println("Error - main variable '" + firstBlood
                         + "' was not found in the database; aborting");
      System.exit(1);
    }
    traits.insertElementAt(traits.remove(traits.indexOf(firstBlood)), 0);

    if (filename.equals("default")) {
      if (!new File(DEFAULT_FILE).exists()) {
        System.out.println("Warning - default file (" + DEFAULT_FILE
                           + ") does not exist, so exporting info for all individuals");
        filename = "";
      } else {
        System.out.println("Using default file (" + DEFAULT_FILE + ")");
        filename = DEFAULT_FILE;
      }
    }
    try {
      new dbExport(db_file, filename, traits, uniqueIDs);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
