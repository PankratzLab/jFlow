// looks for duplicate individuals (2 DNA numbers or 2 identical FamID-IndID pairs) in a
// chromosome.dat file
package org.genvisis.link.init;

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
import org.genvisis.common.Files;
import org.genvisis.common.Sort;

public class checkForDupes {
  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    boolean fix = true;

    String usage = "\n" + "park.checkForDupes requires 0-1 arguments\n"
                   + "   (1) backup and fix chromosome#.dat files with composite (i.e. -fix (not default)\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("-fix")) {
        fix = true;
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      new checkForDupes(fix);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public checkForDupes(boolean fix) throws IOException {
    BufferedReader reader = null;
    PrintWriter writer, copy;
    String temp, trav;
    Hashtable<String, Vector<String>> hash = new Hashtable<String, Vector<String>>();
    Hashtable<String, String[]> hashGenos = new Hashtable<String, String[]>();
    Vector<String> DNAs, dupeDNAs, dupeTravs = new Vector<String>();
    String[] line;
    String[][] dupeGenos;
    int numAgree, numDisagree;
    int[] keys;
    String bakFilename;
    int count = 0;
    String skipped = "";

    while (count <= 23 && !new File("chromosome" + (++count) + ".dat").exists()) {
      ;
    }
    if (count > 23) {
      System.err.println("Error no valid chromosome#.dat files in current directory");
      System.exit(1);
    }

    reader = new BufferedReader(new FileReader("chromosome" + count + ".dat"));
    reader.readLine();
    reader.readLine();
    do {
      line = reader.readLine().split("[\\s]+");
      trav = line[1] + "-" + line[2];
      if (hash.containsKey(trav)) {
        DNAs = hash.get(trav);
        DNAs.add(line[0]);
        if (!dupeTravs.contains(trav)) {
          dupeTravs.add(trav);
        }
      } else {
        DNAs = new Vector<String>();
        DNAs.add(line[0]);
        hash.put(trav, DNAs);
      }
    } while (reader.ready());
    reader.close();

    dupeDNAs = new Vector<String>();
    for (int i = 0; i < dupeTravs.size(); i++) {
      DNAs = hash.get(dupeTravs.elementAt(i));
      for (int j = 0; j < DNAs.size(); j++) {
        dupeDNAs.add(DNAs.elementAt(j));
      }
    }

    copy = new PrintWriter(new FileWriter("dupes-merged.out"));
    for (int chromosome = 1; chromosome <= 23; chromosome++) {
      try {
        reader = new BufferedReader(new FileReader("chromosome" + chromosome + ".dat"));

        copy.println("Chromosome " + chromosome);
        copy.println();

        hashGenos.clear();
        while (reader.ready()) {
          line = reader.readLine().split("[\\s]+");
          if (dupeDNAs.contains(line[0])) {
            hashGenos.put(line[0], line);
          }
        }
        reader.close();

        for (int i = 0; i < dupeTravs.size(); i++) {
          DNAs = hash.get(dupeTravs.elementAt(i));
          dupeGenos = new String[DNAs.size() + 1][];
          for (int j = 0; j < DNAs.size(); j++) {
            dupeGenos[j] = hashGenos.get(DNAs.elementAt(j));
          }
          dupeGenos[dupeGenos.length - 1] = new String[dupeGenos[0].length];

          numAgree = 0;
          numDisagree = 0;
          for (int j = 3; j < dupeGenos[0].length; j++) {
            temp = "-1";
            for (int k = 0; k < dupeGenos.length - 1; k++) {
              if (dupeGenos[k][j].equals("0")) {

              } else if (temp.equals("-1")) {
                temp = dupeGenos[k][j];
              } else if (temp.equals(dupeGenos[k][j])) {
                numAgree++;
              } else {
                numDisagree++;
                temp = "0";
              }
            }
            dupeGenos[dupeGenos.length - 1][j] = temp.equals("-1") ? "0" : temp;
          }

          copy.println("Genotypes agree " + numAgree + " times in " + dupeGenos[0].length
                       + " markers, and "
                       + (numDisagree == 0 ? "never disagree."
                                           : "disagree " + numDisagree + " time(s)!"));
          if (numDisagree != 0) {
            System.out.println(dupeTravs.elementAt(i) + " genotypes agree " + numAgree
                               + " times in " + dupeGenos[0].length + " markers, and disagree "
                               + numDisagree + " time(s)!");
          }
          for (int j = 0; j < dupeGenos.length; j++) {
            copy.print(((j < dupeGenos.length - 1) ? DNAs.elementAt(j)
                                                   : dupeTravs.elementAt(i) + "/CO")
                       + "\t" + (dupeTravs.elementAt(i)).replace('-', '\t'));
            for (int k = 0; k < dupeGenos[0].length; k++) {
              if (numDisagree > 3) {
                dupeGenos[dupeGenos.length - 1][k] = "0";
              }
              copy.print("\t" + dupeGenos[j][k]);
            }
            copy.println();
            hashGenos.put(dupeTravs.elementAt(i) + "/CO", dupeGenos[dupeGenos.length - 1]);
          }
          copy.println();
        }

        if (fix) {
          if (!new File("chromosome" + chromosome + ".dat").exists()) {
            System.err.println("Error - could not find " + "chromosome" + chromosome + ".dat"
                               + " in current directory");
            System.exit(2);
          }
          bakFilename =
              Files.getBakFilename("chromosome" + chromosome + ".dat", super.getClass().getName());
          (new File("chromosome" + chromosome + ".dat")).renameTo(new File(bakFilename));
          reader = new BufferedReader(new FileReader(bakFilename));
          writer = new PrintWriter(new FileWriter("chromosome" + chromosome + ".dat"));

          writer.println(reader.readLine());
          writer.println(reader.readLine());
          while (reader.ready()) {
            temp = reader.readLine();
            line = temp.split("[\\s]+");
            if (dupeDNAs.contains(line[0])) {
              trav = line[1] + "-" + line[2];
              DNAs = hash.get(trav);
              line = hashGenos.get(trav + "/CO");

              if (line.length > 1) {
                keys = Sort.quicksort(Array.toStringArray(DNAs), Sort.DESCENDING);
                writer.print(DNAs.elementAt(keys[0]) + "\t" + trav.replace('-', '\t'));
                for (int i = 3; i < line.length; i++) {
                  writer.print("\t" + line[i]);
                }
                writer.println();
              }
            } else {
              writer.println(temp);
            }
          }
          reader.close();
          writer.close();
        }
      } catch (FileNotFoundException fnfe) {
        skipped += " " + chromosome;
      } catch (IOException ioe) {
        System.err.println("Error parsing " + "chromosome" + chromosome + ".dat" + "");
        System.exit(3);
      }
    }
    copy.close();
    if (skipped.length() > 0) {
      System.err.println("Skipped chromosomes: " + skipped);
    }
  }
}
