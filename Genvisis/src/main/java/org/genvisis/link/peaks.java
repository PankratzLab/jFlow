package org.genvisis.link;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.StringTokenizer;

import org.genvisis.common.Files;
import org.genvisis.common.ext;

public class peaks {

  public static double CUTOFF = 0.995;

  public static void main(String[] args) throws IOException {
    if (args.length > 1) {
      System.out.println("Expecting 0-1 arguments.");
    } else if (args.length == 1) {
      CUTOFF = Double.valueOf(args[0]).doubleValue();
      try {
        new peaks();
      } catch (Exception e) {
        e.printStackTrace();
      }
    } else {
      try {
        new peaks();
      } catch (Exception e) {
        e.printStackTrace();
      }
    }
  }

  public peaks() throws IOException {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String temp, chrome;
    StringTokenizer st;
    double mm, dv, al, ad, adh, ar, arh, ms, me, score;

    st = new StringTokenizer((new File(".")).getAbsolutePath(), Files.isWindows() ? "\\" : "/");
    while (st.countTokens() > 3) {
      st.nextElement();
    }
    writer = new PrintWriter(new FileWriter(st.nextToken() + " " + st.nextToken() + " peaks.out"));
    writer.println(
        "chr\tMapmaker\tDv\t\tMerlin-sibs\tMerlin-ext\tAllegro\tDominant\tAD het\tRecessive\tAR het");

    for (int i = 1; i <= 23; i++) {
      try {
        mm = dv = al = ad = adh = ar = arh = ms = me = -100;
        chrome = (i < 10) ? "0" + i : "" + i;
        reader = new BufferedReader(new FileReader("summary" + chrome + ".out"));

        reader.readLine();
        temp = reader.readLine();
        if (temp.startsWith("no ")) {
          me = -99;
        } else {
          do {
            st = new StringTokenizer(temp);
            st.nextToken();
            score = Double.valueOf(st.nextToken()).doubleValue();
            if (score > mm) {
              mm = score;
            }
            score = Double.valueOf(st.nextToken()).doubleValue();
            if (score > dv) {
              dv = score;
            }
            temp = reader.readLine();
          } while (!temp.equals(""));
        }

        do {
          temp = reader.readLine();
        } while (!temp.startsWith("Merlin-sibpairs"));

        temp = reader.readLine();
        if (temp.startsWith("no ")) {
          ms = -99;
        } else {
          do {
            st = new StringTokenizer(temp);
            st.nextToken();
            score = Double.valueOf(st.nextToken()).doubleValue();
            if (score > ms) {
              ms = score;
            }
            temp = reader.readLine();
          } while (!temp.equals(""));
        }

        do {
          temp = reader.readLine();
        } while (!temp.startsWith("Merlin-extended"));

        temp = reader.readLine();
        if (temp.startsWith("no ")) {
          me = -99;
        } else {
          do {
            st = new StringTokenizer(temp);
            st.nextToken();
            score = Double.valueOf(st.nextToken()).doubleValue();
            if (score > me) {
              me = score;
            }
            temp = reader.readLine();
          } while (!temp.equals(""));
        }

        do {
          temp = reader.readLine();
        } while (!temp.startsWith("Allegro"));

        temp = reader.readLine();
        if (temp.startsWith("no ")) {
          al = -99;
        } else {
          do {
            st = new StringTokenizer(temp);
            st.nextToken();
            score = Double.valueOf(st.nextToken()).doubleValue();
            if (score > al) {
              al = score;
            }
            temp = reader.readLine();
          } while (!temp.equals(""));
        }

        // do {
        // temp = reader.readLine();
        // } while (!temp.startsWith("Aspex2pt"));
        //
        // temp = reader.readLine();
        // if (temp.startsWith("no ")) {
        // as2pt = -99;
        // } else {
        // do {
        // st = new StringTokenizer(temp);
        // st.nextToken();
        // st.nextToken();
        // score = Double.valueOf(st.nextToken()).doubleValue();
        // if (score > as2pt) {
        // as2pt = score;
        // }
        // temp = reader.readLine();
        // } while (!temp.equals(""));
        // }
        //
        // do {
        // temp = reader.readLine();
        // } while (!temp.startsWith("AspexMpt"));
        //
        // temp = reader.readLine();
        // if (temp.startsWith("no ")) {
        // asMpt = -99;
        // } else {
        // do {
        // st = new StringTokenizer(temp);
        // st.nextToken();
        // score = Double.valueOf(st.nextToken()).doubleValue();
        // if (score > asMpt) {
        // asMpt = score;
        // }
        // temp = reader.readLine();
        // } while (!temp.equals(""));
        // }

        do {
          temp = reader.readLine();
        } while (!temp.startsWith("autosomal dominant"));

        temp = reader.readLine();
        if (temp.startsWith("no ")) {
          ad = -99;
        } else {
          do {
            st = new StringTokenizer(temp);
            st.nextToken();
            score = Double.valueOf(st.nextToken()).doubleValue();
            if (score > ad) {
              ad = score;
            }
            temp = reader.readLine();
          } while (!temp.equals(""));
        }

        do {
          temp = reader.readLine();
        } while (!temp.startsWith("autosomal dominant het"));

        temp = reader.readLine();
        if (temp.startsWith("no ")) {
          adh = -99;
        } else {
          do {
            st = new StringTokenizer(temp);
            st.nextToken();
            score = Double.valueOf(st.nextToken()).doubleValue();
            if (score > adh) {
              adh = score;
            }
            temp = reader.readLine();
          } while (!temp.equals(""));
        }

        do {
          temp = reader.readLine();
        } while (!temp.startsWith("autosomal recessive"));

        temp = reader.readLine();
        if (temp.startsWith("no ")) {
          ar = -99;
        } else {
          do {
            st = new StringTokenizer(temp);
            st.nextToken();
            score = Double.valueOf(st.nextToken()).doubleValue();
            if (score > ar) {
              ar = score;
            }
            temp = reader.readLine();
          } while (!temp.equals(""));
        }

        do {
          temp = reader.readLine();
        } while (!temp.startsWith("autosomal recessive het"));

        temp = reader.readLine();
        if (temp.startsWith("no ")) {
          arh = -99;
        } else {
          do {
            st = new StringTokenizer(temp);
            st.nextToken();
            score = Double.valueOf(st.nextToken()).doubleValue();
            if (score > arh) {
              arh = score;
            }
            temp = reader.readLine();
          } while (!temp.equals(""));
        }

        if (i == 23) {
          writer.print("X");
        } else {
          writer.print(i);
        }

        writer.print("\t" + ((mm <= -99) ? "-" : ext.formDeci(mm, 1, true)) + " ");
        for (double j = CUTOFF; j < mm; j += 0.5) {
          writer.print("*");
        }
        if (mm < CUTOFF + 0.5) {
          writer.print("\t");
        }

        writer.print("\t" + ((dv <= -99) ? "-" : ext.formDeci(dv, 1, true)) + " ");
        for (double j = CUTOFF; j < dv; j += 0.5) {
          writer.print("*");
        }
        if (dv < CUTOFF + 0.5) {
          writer.print("\t");
        }

        writer.print("\t" + ((ms <= -99) ? "-" : ext.formDeci(ms, 1, true)) + " ");
        for (double j = CUTOFF; j < ms; j += 0.5) {
          writer.print("*");
        }
        if (ms < CUTOFF + 0.5) {
          writer.print("\t");
        }

        writer.print("\t" + ((me <= -99) ? "-" : ext.formDeci(me, 1, true)) + " ");
        for (double j = CUTOFF; j < me; j += 0.5) {
          writer.print("*");
        }
        if (me < CUTOFF + 0.5) {
          writer.print("\t");
        }

        writer.print("\t" + ((al <= -99) ? "-" : ext.formDeci(al, 1, true)) + " ");
        for (double j = CUTOFF; j < al; j += 0.5) {
          writer.print("*");
        }
        if (al < CUTOFF + 0.5) {
          writer.print("\t");
        }

        // writer.print("\t"+((asMpt<=-99)?"-":ext.formDeci(asMpt, 1,
        // true))+" ");
        // for (double j=CUTOFF; j<asMpt; j+=0.5) {
        // writer.print("*");
        // }
        // if (asMpt < CUTOFF+0.5) {
        // writer.print("\t");
        // }

        writer.print("\t" + ((ad <= -99) ? "-" : ext.formDeci(ad, 1, true)) + " ");
        if (Double.parseDouble(ext.formDeci(ad, 1, true)) > -10) {
          writer.print("\t");
        }

        writer.print("\t" + ((adh <= -99) ? "-" : ext.formDeci(adh, 1, true)) + " ");
        writer.print("\t");

        writer.print("\t" + ((ar <= -99) ? "-" : ext.formDeci(ar, 1, true)) + " ");
        if (Double.parseDouble(ext.formDeci(ar, 1, true)) > -10) {
          writer.print("\t");
        }

        writer.print("\t" + ((arh <= -99) ? "-" : ext.formDeci(arh, 1, true)) + " ");

        writer.println();
      } catch (Exception e) {
        System.err.println("missed chromosome " + i);
      }
    }
    writer.close();
  }
}
