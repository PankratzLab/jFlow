package org.genvisis.link.bat;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.StringTokenizer;
import java.util.Vector;

import org.genvisis.common.Sort;
import org.genvisis.common.ext;

public class osa {

  public osa(String filename, int start, int stop, boolean dv) throws IOException {
    BufferedReader reader = null;
    PrintWriter writer = null, traitData = null;
    StringTokenizer st = null;
    String temp, key, data, chrome, last, next;
    Vector<String> IDs = new Vector<String>();
    double[] covars;
    double trav, prev;
    int[] keys;
    int increment, countPerInc;
    Hashtable<String, String> hash = new Hashtable<String, String>();

    reader = new BufferedReader(new FileReader(filename));
    temp = reader.readLine();
    st = new StringTokenizer(temp);
    if (st.countTokens() == 2) {
      key = st.nextToken();
      data = st.nextToken();
      hash.put(key, data);
      IDs.add(key);
      while (reader.ready()) {
        st = new StringTokenizer(reader.readLine());
        key = st.nextToken();
        data = st.nextToken();
        hash.put(key, data);
        IDs.add(key);
      }
      // } else if (st.countTokens() == 3) { // if used on an individual
      // basis
      // key = st.nextToken() +"-"+ st.nextToken();
      // data = st.nextToken();
      // hash.put(key, data);
      // IDs.add(key);
      // while (reader.ready()) {
      // st = new StringTokenizer(reader.readLine());
      // key = st.nextToken();
      // data = st.nextToken();
      // hash.put(key, data);
      // IDs.add(key);
      // }
    } else {
      System.err.println("Error: File must contain either 2 or 3 (if Individual ID is included) columns. No exceptions.");
      System.exit(1);
    }
    reader.close();

    covars = new double[IDs.size()];
    for (int i = 0; i < IDs.size(); i++) {
      covars[i] = Double.valueOf(hash.get(IDs.elementAt(i))).doubleValue();
    }
    keys = Sort.quicksort(covars);
    hash.clear();

    for (int chromosome = start; chromosome <= stop; chromosome++) {
      chrome = (chromosome < 10) ? "0" + chromosome : "" + chromosome;
      if (!(new File("chrom" + chrome)).exists()) {
        (new File("chrom" + chrome)).mkdir();
      }

      reader = new BufferedReader(new FileReader("re_chrom" + chrome + ".pre"));
      last = "";
      data = "";
      while (reader.ready()) {
        temp = reader.readLine();
        st = new StringTokenizer(temp);
        next = st.nextToken();
        if (!next.equals(last)) {
          if (!last.equals("")) {
            hash.put(last, data);
          }
          data = "";
        }
        data += temp + "\n";
        last = next;
      }
      hash.put(last, data);
      reader.close();

      countPerInc = 0;
      traitData = new PrintWriter(new FileWriter("AscendingKey.dat"));
      traitData.println("Rank\tNum Fams\tvalue");
      increment = 0;
      prev = -99999.777;
      for (int i = 0; i <= IDs.size(); i++) {
        if (i != IDs.size()) {
          trav = covars[keys[i]];
        } else {
          trav = -99999.777;
        }
        if (trav != prev && prev != -99999.777) {
          increment++;
          writer = new PrintWriter(new FileWriter("chrom" + chrome + "/re_chrom" + chrome + "-a"
                                                  + ext.formNum(increment, 4) + ".pre"));
          for (int j = 0; j < i; j++) {
            if (hash.containsKey(IDs.elementAt(keys[j]))) {
              writer.print(hash.get(IDs.elementAt(keys[j])));
            } else {
              System.err.println("Error: could not find family " + IDs.elementAt(keys[j]));
            }
          }
          writer.close();
          traitData.println(increment + "\t" + countPerInc + "\t" + prev);
          // countPerInc = 0;
        }
        prev = trav;
        countPerInc++;
      }
      traitData.close();

      countPerInc = 0;
      traitData = new PrintWriter(new FileWriter("DescendingKey.dat"));
      traitData.println("Rank\tNum Fams\tvalue");
      increment = 0;
      prev = -99999.777;
      for (int i = IDs.size() - 1; i >= -1; i--) {
        if (i != -1) {
          trav = covars[keys[i]];
        } else {
          trav = -99999.777;
        }
        if (trav != prev && prev != -99999.777) {
          increment++;
          writer = new PrintWriter(new FileWriter("chrom" + chrome + "/re_chrom" + chrome + "-d"
                                                  + ext.formNum(increment, 4) + ".pre"));
          for (int j = IDs.size() - 1; j > i; j--) {
            if (hash.containsKey(IDs.elementAt(keys[j]))) {
              writer.print(hash.get(IDs.elementAt(keys[j])));
            } else {
              System.err.println("Error: could not find family " + IDs.elementAt(keys[j]));
            }
          }
          writer.close();
          traitData.println(increment + "\t" + countPerInc + "\t" + prev);
          // countPerInc = 0;
        }
        prev = trav;
        countPerInc++;
      }
      traitData.close();

      writer = new PrintWriter(new FileWriter("chrom" + chrome + "/batch"));
      writer.println("#/bin/sh");
      writer.println();
      writer.println("cp ../map" + chrome + ".dat .");
      writer.println("java -classpath /home/npankrat/" + org.genvisis.common.PSF.Java.GENVISIS
                     + " park.bat.dat2loc map" + chrome + ".dat");
      for (int ascending = 0; ascending < 2; ascending++) {
        temp = ascending == 0 ? "a" : "d";
        for (int i = 1; i <= increment; i++) {
          writer.println("echo -e \"" + ((chromosome == 23) ? "sex on\\n" : "")
                         + "pairs\\n3\\nload map" + chrome + ".loc\\nprep re_chrom" + chrome + "-"
                         + temp + ext.formNum(i, 4) + ".pre\\nn\\nscan\\nestimate\\n"
                         + ((chromosome == 23) ? "" : (dv ? "n" : "y") + "\\n") + "chrom" + chrome
                         + "-" + temp + ext.formNum(i, 4) + ".out\\nchrom" + chrome + "-" + temp
                         + ext.formNum(i, 4) + "-share.ps\\nchrom" + chrome + "-" + temp
                         + ext.formNum(i, 4)
                         + "-mls.ps\\nquit\\n\" | /software/bin/sibs > /dev/null");
          writer.println();
        }
        writer.println("rm *.ps");
      }
      writer.println("rm *.pre");
      writer.close();
      try {
        Runtime.getRuntime().exec("chmod +x chrom" + chrome + "/batch");
      } catch (Exception e) {
        e.printStackTrace();
      }
    }

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "trait.dat";
    int start = 1;
    int stop = 23;
    boolean dv = true;

    String usage = "\n" + "park.bat.osa has 3 optional arguments:\n"
                   + "   (1) a trait file - 2 columns: FamID covar (default: file=" + filename
                   + ")\n" + "   (2) a chromosome number (default: chr=all)\n"
                   + "   (3) whether to use dominance variance (default: dv=true)\n" + "";

    for (int i = 0; i < args.length; i++) {
      if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h")
          || args[i].equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (args[i].startsWith("file=")) {
        filename = args[i].split("=")[1];
        numArgs--;
      } else if (args[i].startsWith("chr=")) {
        if (!args[i].split("=")[1].equals("all")) {
          start = stop = Integer.valueOf(args[i].split("=")[1]).intValue();
        }
        numArgs--;
      } else if (args[i].startsWith("dv=")) {
        if (args[i].split("=")[1].equals("true")) {
          dv = true;
        } else if (args[i].split("=")[1].equals("false")) {
          dv = false;
        } else {
          System.err.println("Error in syntax - '" + args[i].split("=")[1]
                             + "' is not a valid flag for dv (use true/false)");
        }
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    if (args.length == 0) {
      System.err.println("Warning: using defaults (file=" + filename + ", chr=all, dv=true)");
    }

    try {
      new osa(filename, start, stop, dv);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
