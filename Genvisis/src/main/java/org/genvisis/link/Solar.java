// to do - missing first individual
// to do - don't put 99+ in ptypes file
// to do - any reason why h2 is so far off?

package org.genvisis.link;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.StringTokenizer;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;

public class Solar {
  public static int[] MAX_CM = {285, 275, 231, 212, 212, 198, 191, 172, 167, 175, 164, 167, 124,
                                124, 127, 133, 139, 122, 103, 103, 63, 70, 176};

  public static void createFiles(int chr, String trait) throws IOException {
    BufferedReader reader = null;
    PrintWriter writer = null, famtastic = null;
    String temp, famID, indID;
    StringTokenizer st;
    String[] line;
    Vector<String> markerV = new Vector<String>();
    int numMarkers, numAlleles;
    double total;
    Hashtable<String, String> hash;

    reader = new BufferedReader(new FileReader("map" + ext.chrome(chr) + ".dat"));
    writer = new PrintWriter(new FileWriter("solar.freqs." + chr));

    st = new StringTokenizer(reader.readLine());
    numMarkers = Integer.valueOf(st.nextToken()).intValue() - 1;
    for (int i = 0; i < 6; i++) {
      reader.readLine();
    }
    if (chr == 23) {
      reader.readLine();
    }
    for (int i = 0; i < numMarkers; i++) {
      st = new StringTokenizer(reader.readLine());
      st.nextToken();
      numAlleles = Integer.valueOf(st.nextToken()).intValue();
      st.nextToken();
      temp = st.nextToken();
      markerV.add(temp);
      writer.print(ext.formStr(temp, 9, true));
      st = new StringTokenizer(reader.readLine());
      for (int j = 1; j <= numAlleles; j++) {
        writer.print(" " + j + " " + st.nextToken());
      }
      writer.println();
    }
    writer.close();

    reader.readLine();
    st = new StringTokenizer(reader.readLine());
    st.nextToken();

    writer = new PrintWriter(new FileWriter("solar.map." + chr));
    writer.println(chr);
    total = 0;
    writer.println(ext.formStr(markerV.elementAt(0), 9, true) + "   0.00");
    for (int i = 1; i < numMarkers; i++) {
      total += Double.valueOf(st.nextToken()).doubleValue();
      writer.println(ext.formStr(markerV.elementAt(i), 9, true) + " "
                     + ext.formStr(ext.formDeci(total, 2, true), 6));
    }
    reader.close();
    writer.close();

    reader = new BufferedReader(new FileReader("re_chrom" + ext.chrome(chr) + ".pre"));
    writer = new PrintWriter(new FileWriter("solar.gtypes." + chr));
    famtastic = new PrintWriter(new FileWriter("solar.fam"));
    writer.print("FAMID,ID");
    famtastic.println("FAMID,ID,FA,MO,SEX");
    for (int i = 0; i < markerV.size(); i++) {
      writer.print("," + markerV.elementAt(i));
    }
    writer.println();
    while (reader.ready()) {
      st = new StringTokenizer(reader.readLine());
      famID = st.nextToken();
      indID = st.nextToken();
      writer.print(famID + ",");
      writer.print(indID);
      famtastic.println(famID + "," + indID + "," + st.nextToken() + "," + st.nextToken() + ","
                        + st.nextToken());
      st.nextToken();
      do {
        writer.print("," + st.nextToken() + "/" + st.nextToken());
      } while (st.hasMoreElements());
      writer.println();
    }
    reader.close();
    writer.close();
    famtastic.close();

    // hash = tools.pullTraitFromDB(trait);
    hash = HashVec.loadFileToHashString(trait, new int[] {0, 1}, new int[] {2}, trait.endsWith(","),
                                        null, false, false, false);

    reader = new BufferedReader(new FileReader("re_chrom" + ext.chrome(chr) + ".pre"));
    writer = new PrintWriter(new FileWriter("solar.ptypes"));
    // writer.println("FAMID,ID,"+trait);
    writer.println("FAMID,ID,trait");

    while (reader.ready()) {
      line = reader.readLine().split("[\\s]+");
      temp = hash.containsKey(line[0] + "\t" + line[1]) ? hash.get(line[0] + "\t" + line[1]) : "";
      writer.println(line[0] + "," + line[1] + "," + (temp.equals(".") ? "" : temp));
    }
    reader.close();
    writer.close();
  }

  public static void qsub(String filename, String trait) {
    String commands;
    String[][] iterations;

    // commands = ""+
    // "mkdir chrom##\n"+
    //// "java "+classpath+" park.bat.createSolar chr=# trait="+trait+"\n"+
    // "cp solar.fam chrom##\n"+
    // "mv solar.map.# chrom##\n"+
    // "mv solar.freqs.# chrom##\n"+
    // "mv solar.gtypes.# chrom##\n"+
    // "cp "+filename+" chrom##\n"+
    // "cd chrom##\n"+
    // "echo -e \"load pedigree solar.fam\\nload freq solar.freqs.#\\nload marker
    // solar.gtypes.#\\nibddir .\\nverbosity min\\nibd\\nload map solar.map.#\\nibddir .\\nmibddir
    // .\\nmibd 0 [%2] 1\\nmibddir .\\nautomodel "+filename+" "+trait+"\\npolygenic
    // -screen\\nmibddir .\\nchromosome #\\ninterval 1\\nmultipoint -overwrite\\nquit\\n\" | solar >
    // solar.log\n"+
    // "cd ..\n"+
    // "";

    commands = "" + "mkdir chrom[%1]\n" +
    // "java "+classpath+" park.bat.createSolar chr=# trait="+trait+"\n"+
               "cp solar.fam chrom[%1]\n" + "mv solar.map.[%0] chrom[%1]\n"
               + "mv solar.freqs.[%0] chrom[%1]\n" + "mv solar.gtypes.[%0] chrom[%1]\n" + "cp "
               + filename + " chrom[%1]\n" + "cd chrom[%1]\n"
               + "echo -e \"load pedigree solar.fam\\nload freq solar.freqs.[%0]\\nload marker solar.gtypes.[%0]\\nibddir .\\nverbosity min\\nibd\\nload map solar.map.[%0]\\nibddir .\\nmibddir .\\nmibd 0 [%2] 1\\nmibddir .\\nautomodel solar.ptypes trait\\npolygenic -screen\\nmibddir .\\nchromosome [%0]\\ninterval 1\\nmultipoint -overwrite\\nquit\\n\" | /share/apps/bin/solar > solar.log\n"
               + "cd ..\n" + "";

    iterations = new String[22][];
    for (int chr = 1; chr <= 22; chr++) {
      iterations[chr - 1] = new String[] {chr + "", ext.chrome(chr), MAX_CM[chr - 1] + ""};
    }

    Files.qsub("solar", null, 22, commands, iterations, 10000, 24);
  }

  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    int chr = -1;
    // String trait = "AOO";
    String trait = "pheno.dat";
    boolean batch = false;

    // qsub("solar.ptypes", "trait");
    // System.exit(1);

    String usage = "\n" + "link.Solar requires 1-2 arguments\n"
                   + "   (1) chromosome number (i.e. chr=2)\n" + "   (2) trait (i.e. trait=" + trait
                   + " (default)\n" + "   (3) batch (i.e. batch=" + batch + " (default)\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("chr=")) {
        chr = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("trait=")) {
        trait = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("batch=")) {
        batch = ext.parseBooleanArg(arg);
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    // if (chr==-1) {
    // System.err.println("Error - createSolar requires you to specify which chromosome to do");
    // System.exit(2);
    // }
    try {
      if (batch) {
        qsub("solar.ptypes", "trait");
      } else if (chr == -1) {
        for (chr = 1; chr <= 23; chr++) {
          createFiles(chr, trait);
        }
      } else {
        createFiles(chr, trait);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void cleanPhenotypeFile(String filename) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;

    try {
      reader = new BufferedReader(new FileReader(Files.backup(filename, "", "")));
      writer = new PrintWriter(new FileWriter(filename));
      writer.println(reader.readLine());
      while (reader.ready()) {
        line = reader.readLine().split(",", -1);
        for (int i = 2; i < line.length; i++) {
          if (ext.isMissingValue(line[i])) {
            line[i] = "";
          }
        }
        writer.println(Array.toStr(line, ","));
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }
  }
}
