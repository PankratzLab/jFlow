// created a vector of vectors, no expected problems
package org.genvisis.link.bat;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.StringTokenizer;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.DoubleVector;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;
import org.genvisis.link.LinkageMap;
import org.genvisis.link.Markers;

import com.google.common.primitives.Doubles;

public class Filesystem {
  public static void create(int chr) throws IOException {
    create(chr, true);
  }

  public static void create(int chr, boolean reportOutliers) throws IOException {
    String chrome;
    BufferedReader reader, gen = null;
    PrintWriter writer, error;
    StringTokenizer st;
    String[] line;
    Vector<String> orderedMarkers = new Vector<String>();
    DoubleVector distMarkers = new DoubleVector();
    String[] unorderedMarkers;
    int[] key, quickkey, allAlleles, intAllelesForKey;
    int numMarkers;
    String[] markers;
    String temp, IDed, first, second;
    int tempI, total, count = 0, alleleCount;
    Hashtable<String, String> handle, changes;
    Hashtable<String, Hashtable<String, String>> hash =
        new Hashtable<String, Hashtable<String, String>>();
    Enumeration<String> hashKeys;
    Vector<String> alleles;
    Vector<Vector<String>> alleleSizes;
    int[][] problemAlleles;
    int numMissing;
    double[][] alleleFreqs;

    chrome = ext.chrome(chr);
    try {
      gen = new BufferedReader(new FileReader("chromosome" + chr + ".dat"));
    } catch (FileNotFoundException fnfe) {
      System.err.println("chromosome" + chr + ".dat was not found, failed to make mrkr" + chrome
                         + ".dat");
      System.exit(1);
    }
    writer = new PrintWriter(new FileWriter("mrkr" + chrome + ".dat"));
    error = new PrintWriter(new FileWriter("logfile of errors.out", true));

    gen.readLine();
    st = new StringTokenizer(gen.readLine());
    st.nextToken();
    st.nextToken();
    st.nextToken();
    numMarkers = st.countTokens();
    unorderedMarkers = new String[numMarkers];
    alleleSizes = HashVec.newVecVecString(numMarkers);
    for (int i = 0; i < numMarkers; i++) {
      unorderedMarkers[i] = st.nextToken();
    }

    // new orderMarkers(unorderedMarkers);
    Markers.order(unorderedMarkers, true);

    // (new File("markerMap.dat")).renameTo(new
    // File("/home/npankrat/park/00masters/chromosome"+chromosome+"map.txt"));
    // map = new BufferedReader(new
    // FileReader("/home/npankrat/park/00masters/chromosome"+chromosome+"map.txt"));
    new File("chromosome" + chr + "map.txt").delete();
    (new File("markerMap.dat")).renameTo(new File("chromosome" + chr + "map.txt"));
    reader = new BufferedReader(new FileReader("chromosome" + chr + "map.txt"));

    reader.readLine();
    st = new StringTokenizer(reader.readLine());
    orderedMarkers.add(st.nextToken());
    distMarkers.add(10.0);
    while (reader.ready()) {
      st = new StringTokenizer(reader.readLine());
      distMarkers.add(Double.parseDouble(st.nextToken()));
      st = new StringTokenizer(reader.readLine());
      orderedMarkers.add(st.nextToken());
    }
    reader.close();
    numMarkers = orderedMarkers.size();

    key = new int[numMarkers];
    for (int i = 0; i < numMarkers; i++) {
      key[i] = -1;
      for (int j = 0; j < numMarkers; j++) {
        if (orderedMarkers.elementAt(i).equals(unorderedMarkers[j])) {
          key[i] = j;
          break;
        }
      }
      if (key[i] == -1) {
        System.err.println("markers in data file are fubar at marker # " + i + " couldn't find "
                           + orderedMarkers.elementAt(i));
      }
    }

    for (int i = 0; i < numMarkers; i++) {
      hash.put(i + "", new Hashtable<String, String>());
    }
    markers = new String[numMarkers * 2];
    temp = gen.readLine();
    boolean done = false;
    while (!done) {
      count++;
      st = new StringTokenizer(temp);
      st.nextToken();
      IDed = st.nextToken() + "\t" + st.nextToken();
      writer.print(IDed);
      for (int i = 0; i < numMarkers * 2; i++) {
        temp = st.nextToken();
        markers[i] = temp;
      }
      for (int i = 0; i < numMarkers; i++) {
        temp = markers[key[i] * 2];
        handle = hash.get(i + "");
        if (handle.containsKey(temp)) {
          tempI = Integer.valueOf(handle.get(temp)).intValue();
          handle.put(temp, (tempI + 1) + "");
        } else {
          handle.put(temp, "1");
        }
        first = temp;
        temp = markers[key[i] * 2 + 1];
        handle = hash.get(i + "");
        if (handle.containsKey(temp)) {
          tempI = Integer.valueOf(handle.get(temp)).intValue();
          handle.put(temp, (tempI + 1) + "");
        } else {
          handle.put(temp, "1");
        }
        second = temp;

        if (first.equals(".")) {
          System.err.println("FYI: There's a '.' instead of a '0' in chromosome" + chr + ".dat");
          first = "0";
        }
        if (second.equals(".")) {
          System.err.println("FYI: There's a '.' instead of a '0' in chromosome" + chr + ".dat");
          second = "0";
        }
        if (!first.equals("0") && !alleleSizes.elementAt(i).contains(first)) {
          alleleSizes.elementAt(i).add(first);
        }
        if (!second.equals("0") && !alleleSizes.elementAt(i).contains(second)) {
          alleleSizes.elementAt(i).add(second);
        }

        // what do you do in this situation?
        // if (first.equals("0") || second.equals("0")) {
        // writer.print(" 0 0");
        // if (first.equals("0") && second.equals("0")) {} else {
        // error.println(IDed+" has one allele missing
        // ("+ext.formStr(first,4) + ext.formStr(second,4)+") for marker
        // number "+i+" on chromosome "+chromosome);
        // }
        // } else {
        // writer.print(ext.formStr(first,4) + ext.formStr(second,4));
        // }
        if (first.equals("0") || second.equals("0")) {
          if (!first.equals("0")) {
            writer.print(ext.formStr(first, 4) + ext.formStr(first, 4));
          } else {
            writer.print("   0   0");
            if (first.equals("0") && second.equals("0")) {
            } else {
              error.println(IDed + " has one allele missing (" + ext.formStr(first, 4)
                            + ext.formStr(second, 4) + ") for marker number " + i
                            + " on chromosome " + chr);
            }
          }
        } else {
          writer.print(ext.formStr(first, 4) + ext.formStr(second, 4));
        }
      }
      writer.println();
      if (gen.ready()) {
        temp = gen.readLine();
      } else {
        done = true;
      }
      if (temp.startsWith(" ") && !temp.startsWith("\t")) {
        done = true;
      }
    }
    gen.close();
    writer.close();

    problemAlleles = new int[numMarkers][];
    alleleFreqs = new double[numMarkers][];
    for (int i = 0; i < numMarkers; i++) {
      handle = hash.get(i + "");
      total = count * 2;
      if (handle.containsKey("0")) {
        total -= Integer.valueOf(handle.get("0")).intValue();
      }
      hashKeys = handle.keys();
      intAllelesForKey = new int[handle.size()];
      alleles = new Vector<String>();
      while (hashKeys.hasMoreElements()) {
        temp = hashKeys.nextElement();
        intAllelesForKey[alleles.size()] = Integer.parseInt(temp);
        alleles.add(temp);
      }
      quickkey = Sort.quicksort(intAllelesForKey);
      allAlleles = new int[total];
      alleleCount = 0;
      numMissing = alleles.contains("0") ? 1 : 0;

      alleleFreqs[i] = new double[quickkey.length - numMissing];
      for (int j = numMissing; j < quickkey.length; j++) {
        alleleFreqs[i][j - numMissing] =
            (double) Integer.valueOf(handle.get(alleles.elementAt(quickkey[j]))).intValue()
                                         / (double) total;
        for (int k = 0; k < Integer.valueOf(handle.get(alleles.elementAt(quickkey[j])))
                                   .intValue(); k++) {
          allAlleles[alleleCount++] = Integer.valueOf(alleles.elementAt(quickkey[j])).intValue();
        }

      }
      if (alleleCount != total) {
        System.err.println("Error - for some reason the number of genotypes (" + alleleCount
                           + ") did not reach the expected (" + total + ")");
      }

      if (reportOutliers) {
        problemAlleles[i] = spreadCheck(allAlleles, orderedMarkers.elementAt(i));
      }
    }
    error.close();

    new LinkageMap(chr, Array.toStringArray(orderedMarkers), alleleFreqs,
                   Doubles.toArray(distMarkers), false, false).createFile("map" + chrome + ".dat");

    changes = new Hashtable<String, String>();

    if (reportOutliers) {
      gen = new BufferedReader(new FileReader("mrkr" + chrome + ".dat"));
      writer = new PrintWriter(new FileWriter("genosToCheck.dat", true));
      while (gen.ready()) {
        line = gen.readLine().split("[\\s]+");
        for (int i = 0; i < problemAlleles.length; i++) {
          for (int j = 0; j < problemAlleles[i].length; j++) {
            if (line[i * 2 + 2].equals(problemAlleles[i][j] + "")
                || line[i * 2 + 2 + 1].equals(problemAlleles[i][j] + "")) {
              if (!changes.containsKey(line[0] + "\t" + line[1] + "\t"
                                       + orderedMarkers.elementAt(i))) {
                writer.println(line[0] + "\t" + line[1] + "\t" + orderedMarkers.elementAt(i) + "\t"
                               + line[i * 2 + 2] + "\t" + line[i * 2 + 2 + 1] + "\t\t");
                changes.put(line[0] + "\t" + line[1] + "\t" + orderedMarkers.elementAt(i), "");
              }
            }
          }
        }
      }
      writer.close();
    }
  }

  public static void main(String[] args) throws IOException {
    if (args.length != 1) {
      System.out.println("Expecting 1 argument: chromosome number.");
    } else {
      try {
        create(Integer.parseInt(args[0]));
        // create("2");
      } catch (Exception e) {
        e.printStackTrace();
        System.err.println("Error in processing chromosome " + args[0]);
      }
    }
  }

  public static int[] spreadCheck(int[] source, String header) {
    int[] keys = Sort.quicksort(source);
    Vector<String> missedOpportunities = new Vector<String>();
    int[] missedOpps;
    IntVector missedCounts = new IntVector();
    double mean = Array.mean(Array.toDoubleArray(source));
    double sd = Array.stdev(Array.toDoubleArray(source));
    int ub = (int) (mean + sd * 5), lb = (int) (mean - sd * 5);
    int count;

    if (source.length == 0) {
      System.err.println("There is no data for " + header);
      return new int[0];
    }

    count = 0;
    while (source[keys[count++]] < lb) {
      ;
    }
    lb = source[keys[count]];

    count = source.length - 1;
    while (source[keys[count--]] > ub) {
      ;
    }
    ub = source[keys[count]];

    count = 0;
    while (source[keys[count]] + 10 < lb) {
      // String temp = source[keys[count]]+"";
      if (!missedOpportunities.contains(source[keys[count]] + "")) {
        missedOpportunities.add(source[keys[count]] + "");
        missedCounts.add(1);
      } else {
        int index = missedOpportunities.indexOf(source[keys[count]] + "");
        missedCounts.set(index, missedCounts.get(index) + 1);
      }
      count++;
    }
    count = source.length - 1;
    while (source[keys[count]] - 10 > ub) {
      if (!missedOpportunities.contains(source[keys[count]] + "")) {
        missedOpportunities.add(source[keys[count]] + "");
        missedCounts.add(1);
      } else {
        int index = missedOpportunities.indexOf(source[keys[count]] + "");
        missedCounts.set(index, missedCounts.get(index) + 1);
      }
      count--;
    }

    missedOpps = new int[missedOpportunities.size()];
    for (int i = 0; i < missedOpportunities.size(); i++) {
      missedOpps[i] = Integer.valueOf(missedOpportunities.elementAt(i)).intValue();
    }

    for (int i = 0; i < missedOpps.length; i++) {
      System.err.print(missedOpps[i] + " "
                       + (missedCounts.elementAt(i) == 1 ? ""
                                                         : "(" + missedCounts.elementAt(i)
                                                           + "x) "));
    }
    if (missedOpps.length > 0) {
      System.err.println("-- out of the usual range of " + lb + " - " + ub + " for " + header);
    }
    return missedOpps;
  }

  public int[] boxPlotOutliers(int[] source, String header) {
    int[] keys = Sort.quicksort(source);
    int q1, q3, iqr, lf, uf;
    Vector<String> missedOpportunities = new Vector<String>();
    int[] missedOpps;

    if (source.length == 0) {
      System.err.println("There is no data for " + header);
      return new int[0];
    }

    q1 = source[keys[(int) Math.ceil(source.length * 0.10)]];
    q3 = source[keys[(int) Math.floor(source.length * 0.90)]];
    iqr = q3 - q1;
    lf = q1 - iqr * 5;
    uf = q3 + iqr * 5;

    for (int element : source) {
      if (element < lf || element > uf) {
        System.err.println(element + " is out of the expected range of " + lf + " - " + uf + " for "
                           + header);
        missedOpportunities.add(element + "");
      }
    }

    missedOpps = new int[missedOpportunities.size()];
    for (int i = 0; i < missedOpportunities.size(); i++) {
      missedOpps[i] = Integer.valueOf(missedOpportunities.elementAt(i)).intValue();
    }
    return missedOpps;
  }
}
