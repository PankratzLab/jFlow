package org.genvisis.link;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.StringTokenizer;

import org.genvisis.common.ext;

public class empiricalPvalues {
  // public static int NUM_REPS = 5000;

  public static void main(String[] args) throws IOException {
    if (args.length == 0 || args.length > 2) {
      System.out.println("Requires 1-2 arguments: chromosome number and [optional] .");
    } else {
      try {
        if (args.length == 1) {
          new empiricalPvalues(Integer.valueOf(args[0]).intValue(), 5000);
        } else {
          new empiricalPvalues(Integer.valueOf(args[0]).intValue(),
                               Integer.valueOf(args[1]).intValue());
        }
      } catch (Exception e) {
        e.printStackTrace();
      }
    }
  }

  public empiricalPvalues(int chromosome, int numReps) throws IOException {
    BufferedReader reader = null;
    BufferedReader reader2 = null;
    PrintWriter writer = null;
    PrintWriter optfile = null;
    String temp, chrome, file = null;
    StringTokenizer st;
    Hashtable<String, boolean[]> hash = new Hashtable<String, boolean[]>();
    int numMarkers;
    boolean[] handle;
    boolean done = false;
    String trav, pos, sibsPos, sibsDvPos;
    double score, sibs, sibsDv;
    boolean isMale;

    chrome = (chromosome < 10) ? "0" + chromosome : "" + chromosome;

    reader = new BufferedReader(new FileReader("re_chrom" + chrome + ".pre"));
    writer = new PrintWriter(new FileWriter("chrom" + chrome + "-template"));

    st = new StringTokenizer(reader.readLine());
    numMarkers = (st.countTokens() - 6) / 2;

    while (!done) {
      if (st.countTokens() > 0) {
        trav = st.nextToken() + "\t" + st.nextToken();
        writer.println(trav + "\t" + st.nextToken() + "\t" + st.nextToken() + "\t" + st.nextToken()
                       + "\t" + st.nextToken() + "\t1\t1");

        handle = new boolean[numMarkers];
        for (int i = 0; i < numMarkers; i++) {
          handle[i] = (Integer.valueOf(st.nextToken()).intValue() != 0);
          st.nextToken();
        }
        hash.put(trav, handle);
      }
      if (reader.ready()) {
        st = new StringTokenizer(reader.readLine());
      } else {
        done = true;
      }
    }
    writer.close();

    optfile = new PrintWriter(new FileWriter("chrom" + chrome + ".opt"));
    optfile.println("% Read input in LINKAGE style format:\n" + "PREFILE chrom" + chrome
                    + "-template\n" + "DATFILE map" + chrome + ".dat\n\n"
                    + "% Simulate stroke reconstruction pedigrees\n" + "SIMULATE dloc:10.0 npre:"
                    + numReps + " rep:1 err:0.00 yield:1.0 het:1\n\n" + "% Other options:\n"
                    + "MAXMEMORY 100");
    optfile.close();

    Process process = null;
    Runtime runtime = Runtime.getRuntime();
    process = runtime.exec("/software/bin/allegro chrom" + chrome + ".opt");

    try {
      process.waitFor();
    } catch (Exception e) {
      e.printStackTrace();
    }

    (new File("chrome" + chrome + ".opt")).delete();

    PrintWriter batch = new PrintWriter(new FileWriter("batch.1"));
    batch.println("cp batch.1 trash1");
    batch.println("cp batch.1 trash2");
    batch.println();
    batch.println("java -classpath /home/npankrat/" + org.genvisis.common.PSF.Java.GENVISIS
                  + " park.bat.dat2loc map" + chrome + ".dat");

    for (int repNum = 1; repNum <= numReps; repNum++) {
      // for (int repNum=NUM_REPS; repNum>=1; repNum--) {
      writer = new PrintWriter(new FileWriter("chrom" + chrome + "-" + repNum + ".pre"));

      file = "chrom" + chrome + "-template."
             + ext.formNum(repNum + "", String.valueOf(numReps).length());
      reader = new BufferedReader(new FileReader(file));
      while (reader.ready()) {
        st = new StringTokenizer(reader.readLine(), "- ");
        trav = st.nextToken() + "\t" + st.nextToken();
        writer.print(trav + "\t" + st.nextToken() + "\t" + st.nextToken());

        isMale = (st.nextToken().equals("1"));

        writer.print("\t" + (isMale ? "1" : "2") + "\t" + st.nextToken());

        handle = hash.get(trav);
        for (int i = 0; i < numMarkers; i++) {
          if (handle[i]) {
            if (chromosome == 23 && isMale) {
              st.nextToken();
              temp = st.nextToken();
              writer.print("\t" + temp + "\t" + temp);
            } else {
              writer.print("\t" + st.nextToken() + "\t" + st.nextToken());
            }
          } else {
            writer.print("\t0\t0");
            st.nextToken();
            st.nextToken();
          }
        }
        writer.println();
      }
      reader.close();
      (new File(file)).delete();

      writer.close();

      if (chromosome != 23) {
        batch.println("echo -e \"pairs\\n3\\nload map" + chrome + ".loc\\nprep chrom" + chrome + "-"
                      + repNum + ".pre\\nn\\nscan\\nestimate\\ny\\nchrom" + chrome + "-" + repNum
                      + "-mls.out\\ntrash1\\ny\\ntrash2\\ny\\nestimate\\nn\\nchrom" + chrome + "-"
                      + repNum
                      + "-Dv-mls.out\\ntrash1\\ny\\ntrash2\\ny\\nquit\\n\" | /software/bin/sibs > chrom"
                      + chrome + "-" + repNum + ".log"); // >
        // /dev/null"
      } else {
        batch.println("echo -e \"sex on\\npairs\\n3\\nload map" + chrome + ".loc\\nprep chrom"
                      + chrome + "-" + repNum + ".pre\\nn\\nscan\\nestimate\\nchrom" + chrome + "-"
                      + repNum
                      + "-mls.out\\ntrash1\\ny\\ntrash2\\ny\\nquit\\n\" | /software/bin/sibs > chrom"
                      + chrome + ".log"); // >
                                          // /dev/null"
      }
      batch.println();
    }

    batch.println("rm trash1");
    batch.println("rm trash2");
    batch.println("rm *.log");
    batch.println("rm chrom" + chrome + "-*.pre");

    batch.close();
    try {
      Process process2 = null;
      Runtime runtime2 = Runtime.getRuntime();
      process2 = runtime2.exec("chmod +x batch.1");
      process2.waitFor();
      Process process3 = null;
      Runtime runtime3 = Runtime.getRuntime();
      process3 = runtime3.exec("./batch.1");
      process3.waitFor();
    } catch (Exception e) {
      e.printStackTrace();
    }

    writer = new PrintWriter(new FileWriter("chrom" + chrome + "-maximumLods.xls"));
    writer.println("rep\tsibsPos\tsibsLOD\tsibsDvPos\tsibsDvLOD");
    for (int i = 1; i <= numReps; i++) {
      try {
        reader = new BufferedReader(new FileReader("chrom" + chrome + "-" + i + "-mls.out"));
        sibs = sibsDv = -1;
        sibsPos = sibsDvPos = "oops";
        if (chromosome < 23) {
          reader2 = new BufferedReader(new FileReader("chrom" + chrome + "-" + i + "-Dv-mls.out"));
          reader.readLine();
          reader2.readLine();
          temp = reader.readLine();
          while (!temp.equals("")) {
            st = new StringTokenizer(temp);
            pos = st.nextToken();
            st.nextToken();
            st.nextToken();
            st.nextToken();
            score = Double.valueOf(st.nextToken()).doubleValue();
            if (score > sibs) {
              sibs = score;
              sibsPos = pos;
            }

            st = new StringTokenizer(reader2.readLine());
            pos = st.nextToken();
            st.nextToken();
            st.nextToken();
            st.nextToken();
            score = Double.valueOf(st.nextToken()).doubleValue();
            if (score > sibsDv) {
              sibsDv = score;
              sibsDvPos = pos;
            }

            temp = reader.readLine();
          }
        } else {
          do {
            temp = reader.readLine();
          } while (!temp.startsWith("Total"));

          temp = reader.readLine();
          while (!temp.equals("")) {
            st = new StringTokenizer(temp);
            pos = st.nextToken();
            score = Double.valueOf(st.nextToken()).doubleValue();
            if (score > sibs) {
              sibs = sibsDv = score;
              sibsPos = sibsDvPos = pos;
            }

            temp = reader.readLine();
          }
        }

        writer.println(i + "\t" + sibsPos + "\t" + sibs + "\t" + sibsDvPos + "\t" + sibsDv);

      } catch (Exception e) {
        writer.println("no Mapmaker data");
        e.printStackTrace();
        System.err.println("Error processing Mapmaker files for chromosome " + i);
      }
    }
    writer.close();

    System.out.println("done!");

  }
}
