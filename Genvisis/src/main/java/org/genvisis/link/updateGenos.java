package org.genvisis.link;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Files;

public class updateGenos {
  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    String filename = "genoUpdates.dat";

    String usage = "\n" + "park.updateGenos requires 1 argument:\n"
                   + "   (1) a file with three columns - FamID IndID marker_name old_genotype1 old_genotype2 new_genotype1 new_genotype2\n"
                   + "       (i.e. file=" + filename + " (default))\n" + "";

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
    if (args.length == 0) {
      System.err.println("Warning: using defaults (file=" + filename + ")");
    }
    try {
      new updateGenos(filename);
    } catch (Exception e) {
      e.printStackTrace();
    }

  }

  public updateGenos(String filename) throws IOException {
    BufferedReader reader = null;
    PrintWriter writer;
    Hashtable<String, String> allMarkers = new Hashtable<String, String>();
    Hashtable<String, Vector<String>> people;
    Hashtable<String, Hashtable<String, Vector<String>>> chrs =
        new Hashtable<String, Hashtable<String, Vector<String>>>();
    Vector<String> changes, markerNames;
    String temp, bakFilename;
    String[] line, alteration;
    int chr, count;
    boolean madeChange;

    try {
      reader = new BufferedReader(new FileReader("marker.database"));
    } catch (FileNotFoundException fnfe) {
      try {
        reader = new BufferedReader(new FileReader("/home/npankrat/marker.database"));
      } catch (FileNotFoundException fnfe2) {
        System.err.println("Error - could not find marker.database in current or root directory");
        System.exit(1);
      }
    }

    chr = 0;
    while (reader.ready()) {
      temp = reader.readLine();
      if (temp.equals("") || temp.startsWith("Marker")) {
        continue;
      }
      if (temp.startsWith("chromosome")) {
        chr = Integer.valueOf(temp.substring(10, 12)).intValue();
        continue;
      }
      line = temp.split("[ ]+");
      allMarkers.put(line[0].toUpperCase(), chr + "");
    }
    reader.close();

    count = 0;
    if (!new File(filename).exists()) {
      System.err.println("Error - could not find " + filename + " in current directory");
      System.exit(2);
    }
    reader = new BufferedReader(new FileReader(filename));
    while (reader.ready()) {
      temp = reader.readLine();
      line = temp.split("[\\s]+");
      line[2] = line[2].toUpperCase();
      count++;
      if (line.length < 7) {
        System.err.println("Error - incorrect number of clomuns in line " + count + ". Found:");
        System.err.println(temp);
        System.err.println("    Expecting:");
        System.err.println("FamID IndID marker_name old_genotype1 old_genotype2 new_genotype1 new_genotype2");
        System.exit(1);
      }
      if (!allMarkers.containsKey(line[2])) {
        System.err.println("Error - cannot find " + line[2]
                           + " in the marshfield database. Check spelling and/or the order of the file. Expecting:");
        System.err.println("FamID IndID marker_name old_genotype1 old_genotype2 new_genotype1 new_genotype2");
        System.exit(2);
      }

      if (chrs.containsKey(allMarkers.get(line[2]))) {
        people = chrs.get(allMarkers.get(line[2]));
      } else {
        chrs.put(allMarkers.get(line[2]), people = new Hashtable<String, Vector<String>>());
      }

      if (people.containsKey(line[0] + "\t" + line[1])) {
        changes = people.get(line[0] + "\t" + line[1]);
      } else {
        people.put(line[0] + "\t" + line[1], changes = new Vector<String>());
      }
      changes.add(line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6]);
    }
    reader.close();

    for (int i = 1; i <= 23; i++) {
      chr = i;
      if (chrs.containsKey(chr + "")) {
        madeChange = false;
        bakFilename = Files.getBakFilename("chromosome" + chr + ".dat", super.getClass().getName());
        if (!new File("chromosome" + chr + ".dat").exists()) {
          System.err.println("Error - could not find " + "chromosome" + chr + ".dat"
                             + " in current directory");
          System.exit(2);
        }
        new File("chromosome" + chr + ".dat").renameTo(new File(bakFilename));
        reader = new BufferedReader(new FileReader(bakFilename));
        writer = new PrintWriter(new FileWriter("chromosome" + chr + ".dat"));
        writer.println(reader.readLine());
        temp = reader.readLine();
        line = temp.split("[\\s]+");
        markerNames = new Vector<String>();
        for (int j = 3; j < line.length; j++) {
          markerNames.add(line[j].toUpperCase());
        }
        writer.println(temp);
        people = chrs.get(chr + "");
        while (reader.ready()) {
          temp = reader.readLine();
          line = temp.split("[\\s]+");
          if (people.containsKey(line[1] + "\t" + line[2])) {
            // System.out.println("Doing "+line[1]+"\t"+line[2]);
            changes = people.get(line[1] + "\t" + line[2]);
            for (int j = 0; j < changes.size(); j++) {
              alteration = (changes.elementAt(j)).split("[\\s]+");
              if (!markerNames.contains(alteration[0])) {
                System.err.println("Error - marker " + alteration[0]
                                   + " was not found in chromosome" + chr + ".dat");
                System.err.println("        Skipping alteration");
              } else {
                if (line[3 + markerNames.indexOf(alteration[0]) * 2].equals(alteration[1])
                    && line[3 + markerNames.indexOf(alteration[0]) * 2 + 1].equals(alteration[2])) {
                  System.out.println("Congrats you found the right genotype ("
                                     + line[3 + markerNames.indexOf(alteration[0]) * 2] + " & "
                                     + line[3 + markerNames.indexOf(alteration[0]) * 2 + 1] + ")");
                  line[3 + markerNames.indexOf(alteration[0]) * 2] = alteration[3];
                  line[3 + markerNames.indexOf(alteration[0]) * 2 + 1] = alteration[4];
                  madeChange = true;
                } else if (line[3 + markerNames.indexOf(alteration[0]) * 2].equals(alteration[2])
                           && line[3 + markerNames.indexOf(alteration[0]) * 2
                                   + 1].equals(alteration[1])) {
                  System.out.println("Congrats you found the right genotype ("
                                     + line[3 + markerNames.indexOf(alteration[0]) * 2] + " & "
                                     + line[3 + markerNames.indexOf(alteration[0]) * 2 + 1]
                                     + ") (albeit reversed)");
                  line[3 + markerNames.indexOf(alteration[0]) * 2] = alteration[3];
                  line[3 + markerNames.indexOf(alteration[0]) * 2 + 1] = alteration[4];
                  madeChange = true;
                } else if (line[3 + markerNames.indexOf(alteration[0]) * 2].equals(alteration[3])
                           && line[3 + markerNames.indexOf(alteration[0]) * 2
                                   + 1].equals(alteration[4])) {
                  System.out.println("Genotype " + alteration[1] + "/" + alteration[2]
                                     + " has already been changed to " + alteration[3] + "/"
                                     + alteration[4] + " for " + line[1] + "-" + line[2] + " ("
                                     + alteration[0] + ")");
                } else {
                  System.err.println("Error - Wrong genotypes specified for " + line[1] + "-"
                                     + line[2] + " (" + alteration[0] + ")");
                  System.err.println("        Supposed to change " + alteration[1] + "/"
                                     + alteration[2] + " to " + alteration[3] + "/" + alteration[4]
                                     + " but found "
                                     + line[3 + markerNames.indexOf(alteration[0]) * 2] + "/"
                                     + line[3 + markerNames.indexOf(alteration[0]) * 2 + 1]);
                }
              }
            }
            writer.print(line[0]);
            for (int j = 1; j < line.length; j++) {
              writer.print("\t" + line[j]);
            }
            writer.println();
          } else {
            writer.println(temp);
          }
        }
        reader.close();
        writer.close();
        if (!madeChange) {
          new File("chromosome" + chr + ".dat").delete();
          new File(bakFilename).renameTo(new File("chromosome" + chr + ".dat"));

        }

      }
    }

  }
}
