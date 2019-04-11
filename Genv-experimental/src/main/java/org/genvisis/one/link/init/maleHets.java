package org.genvisis.one.link.init;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.StringTokenizer;
import java.util.Vector;
import org.pankratzlab.common.Files;

public class maleHets {

  public maleHets() throws IOException {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String temp, fam, id, marker, first, second;
    StringTokenizer st;
    int numLeft;
    Vector<String> listOfErrors = new Vector<>();
    Hashtable<String, String> hash = new Hashtable<>();

    try {
      reader = new BufferedReader(new FileReader("chrom23.pre"));
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: could not find the X pre file \"chrom23.pre\"");
      System.exit(1);
    }
    writer = Files.openAppropriateWriter("logfile of errors.out", true);

    System.out.print("Checking for heterozygotes");
    while (reader.ready()) {
      st = new StringTokenizer(reader.readLine());
      fam = st.nextToken();
      id = st.nextToken();
      st.nextToken();
      st.nextToken();
      if (st.nextToken().equals("1")) {
        st.nextToken();
        numLeft = st.countTokens() / 2;
        for (int i = 1; i <= numLeft; i++) {
          first = st.nextToken();
          second = st.nextToken();
          if (!first.equals(second)) {
            System.out.print(".");
            listOfErrors.add(fam + "\t" + id + "\t" + i);
            writer.println("Individual " + fam + "-" + id + " is heterozygous (" + first + "/"
                           + second + ") at marker " + i);
            if (hash.containsKey(fam + "-" + id)) {
              hash.put(fam + "-" + id,
                       (Integer.valueOf(hash.get(fam + "-" + id)).intValue() + 1) + "");
            } else {
              hash.put(fam + "-" + id, 1 + "");
            }
          }
        }
      }
    }
    reader.close();
    writer.close();

    if (listOfErrors.size() == 0) {
      System.out.print("<<no errors>>");
    } else {
      String bakFilename = Files.getBakFilename("chromosome23.dat", super.getClass().getName());
      (new File("chromosome23.dat")).renameTo((new File(bakFilename)));

      reader = new BufferedReader(new FileReader(bakFilename));
      writer = Files.openAppropriateWriter("chromosome23.dat");

      temp = reader.readLine();
      String crap, lastFam = "";
      fam = "";
      for (int i = 0; i < listOfErrors.size(); i++) {
        st = new StringTokenizer(listOfErrors.elementAt(i));
        lastFam = fam;
        fam = st.nextToken();
        id = st.nextToken();
        marker = st.nextToken();

        if (fam.equals(lastFam)) {
          writer.println(temp);
          while (reader.ready()) {
            writer.println(reader.readLine());
          }
          reader.close();
          writer.close();
          (new File("chromosome23.dat")).renameTo((new File("temp")));
          reader = new BufferedReader(new FileReader("temp"));
          writer = Files.openAppropriateWriter("chromosome23.dat");
          temp = reader.readLine();
        }

        st = new StringTokenizer(temp);
        crap = st.nextToken();
        while (!st.nextToken().equals(fam)) {
          writer.println(temp);
          temp = reader.readLine();
          st = new StringTokenizer(temp);
          crap = st.nextToken();
        }
        while (!st.nextToken().equals(id)) {
          writer.println(temp);
          temp = reader.readLine();
          st = new StringTokenizer(temp);
          crap = st.nextToken();
          st.nextToken();
        }
        writer.print(crap + "\t" + fam + "\t" + id);
        for (int j = 0; j < 2 * (Integer.valueOf(marker).intValue() - 1); j++) {
          writer.print("\t" + st.nextToken());
        }
        writer.print("\t0\t0");
        st.nextToken();
        st.nextToken();
        while (st.hasMoreTokens()) {
          writer.print("\t" + st.nextToken());
        }
        writer.println();
        temp = reader.readLine();
      }
      writer.println(temp);
      while (reader.ready()) {
        writer.println(reader.readLine());
      }

      reader.close();
      writer.close();
    }
    System.out.println("done");

    Enumeration<String> enumer = hash.keys();
    Vector<String> male_identity_crises = new Vector<>();
    while (enumer.hasMoreElements()) {
      temp = enumer.nextElement();
      if (Integer.valueOf(hash.get(temp)).intValue() >= 3) {
        male_identity_crises.add(temp);
      }
    }
    Collections.sort(male_identity_crises);
    for (String m : male_identity_crises) {
      System.err.println("Warning - " + m + " is heterozygous for " + hash.get(m)
                         + " markers on the X chromosome");
    }

  }

  public static void main(String[] args) throws IOException {
    if (args.length > 0) {
      System.out.println("Expecting no arguments.");
    }
    try {
      new maleHets();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}