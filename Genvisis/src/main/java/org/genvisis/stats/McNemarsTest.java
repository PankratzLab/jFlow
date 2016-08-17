package org.genvisis.stats;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.HashVec;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;

// non-paramteric test for paired binary data
public class McNemarsTest {
  public static final int THRESHOLD = 25;

  public static void batch(String dir, String pairs, String database, String[] analysisVariables) {
    BufferedReader reader;
    PrintWriter writer;
    String trav;
    Hashtable<String, String> hash;
    Vector<String[]> v;
    String[][] pairings;
    String[][][] data;
    int[][] round;
    int count;
    McNemarsTest mt;

    v = new Vector<String[]>();
    try {
      reader = new BufferedReader(new FileReader(dir + pairs));
      while (reader.ready()) {
        v.add(reader.readLine().trim().split("[\\s]+"));
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + pairs + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + pairs + "\"");
      System.exit(2);
    }
    pairings = Matrix.toStringArrays(v);
    System.out.println("Found " + pairings.length + " pairings in " + pairs);

    hash = HashVec.loadFileToHashString(dir + database, "UniqueID", analysisVariables, "\t");
    data = new String[pairings.length][2][];
    for (int i = 0; i < pairings.length; i++) {
      for (int j = 0; j < 2; j++) {
        trav = hash.get(pairings[i][j]);
        if (trav == null) {
          System.err.println("Error - UniqueID '" + pairings[i][j] + "' was not found in the file '"
              + database + "'");
        } else {
          data[i][j] = trav.split("[\\s]+");
        }
      }
    }
    try {
      writer = new PrintWriter(new FileWriter(dir + ext.rootOf(pairs) + "_results.xln"));
      writer.println(
          "Variable\tn\ta\tb\tc\td\trisk\tchisq\tp-value\tchisq with continuity\tp-value with continuity\tsign-test\tproper_pvalue");
      for (int i = 0; i < analysisVariables.length; i++) {
        System.err.println("Analyzing " + analysisVariables[i]);
        round = new int[pairings.length][2];
        count = 0;
        for (int j = 0; j < pairings.length; j++) {
          for (int k = 0; k < 2; k++) {
            if (data[j][k][i].equals(".")) {
              round[j][k] = -1;
            } else if (data[j][k][i].equals("1")) {
              round[j][k] = 1;
            } else if (data[j][k][i].equals("0")) {
              round[j][k] = 0;
            } else {
              System.err.println("Error - invalid value: " + data[j][k][i]);
              round[j][k] = -1;
            }
          }
          if (round[j][0] != -1 && round[j][1] != -1) {
            count++;
          }
        }
        mt = new McNemarsTest(round);
        mt.setContinuity(false);
        writer.print(analysisVariables[i] + "\t" + count + "\t" + Array.toStr(mt.getABCD()) + "\t"
            + mt.getDirectionOfRisk() + "\t" + mt.getChiSq() + "\t" + mt.getPvalue());
        mt.setContinuity(true);
        writer.print("\t" + mt.getChiSq() + "\t" + mt.getPvalue());
        writer.print("\t" + mt.getPvalueFromBinomialDistribution());
        writer.println("\t" + mt.getPvalue());
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + dir + ext.rootOf(pairs) + "_results.xln");
      e.printStackTrace();
    }

  }

  private int a;
  private int b;
  private int c;
  private int d;
  private double chisq;
  private String variable;

  private boolean continuity;

  public McNemarsTest(int a11, int b10, int c01, int d00) {
    a = a11;
    b = b10;
    c = c01;
    d = d00;
    chisq = compute();
    continuity = true;
  }

  public McNemarsTest(int[][] matrix) {
    int count;

    a = b = c = d = count = 0;
    continuity = true;
    for (int[] element : matrix) {
      if (element[0] == -1 || element[1] == -1) {
        count++;
      } else if (element[0] == 1 && element[1] == 1) {
        a++;
      } else if (element[0] == 1 && element[1] == 0) {
        b++;
      } else if (element[0] == 0 && element[1] == 1) {
        c++;
      } else if (element[0] == 0 && element[1] == 0) {
        d++;
      } else {
        System.err.println(
            "Error - invalid input for McNemar's test (); must be 0, 1, or missing (-1); use the Cochran test for more than two states (not yet implemented)");
      }
    }
    if (count == matrix.length) {
      System.err.println("Error - no valid data");
      a = b = c = d = 0;
      chisq = Double.NaN;
    } else {
      chisq = compute();
    }
  }

  public double compute() {
    if (b + c == 0) {
      System.err.println("Error - Since b and c are both zero, McNemar's test cannot be computed");
      return Double.NaN;
    }
    if (b + c < 20) {
      System.err.println("Error - Since b (" + b + ") + c (" + c + ") is less than 20"
          + (variable == null ? "" : " for variable " + variable)
          + ", a sign test should be used instead of McNemar's test (unfortunately not yet implemented)");
    }
    return Math.pow((double) Math.abs(b - c) - (continuity ? 1 : 0), 2) / (b + c);
  }

  public int[] getABCD() {
    return new int[] {a, b, c, d};
  }

  public double getChiSq() {
    return chisq;
  }

  public String getDirectionOfRisk() {
    return b > c ? "+" : "-";
  }

  public String getDominantDiscordant() {
    return b > c ? "10" : "01";
  }

  public double getPvalue() {
    return (b + c < THRESHOLD ? getPvalueFromBinomialDistribution()
        : getPvalueFromMcNemarTestNoMatterWhat());
  }

  public double getPvalueFromBinomialDistribution() {
    return Math.min(BinomialDistribution.probabilityLTE(Math.min(b, c), b + c, 0.5) * 2, 1);
  }

  public double getPvalueFromMcNemarTestNoMatterWhat() {
    return ProbDist.ChiDist(chisq, 1);
  }

  public void setContinuity(boolean continuity) {
    this.continuity = continuity;
    chisq = compute();
  }

  public void setVariable(String name) {
    variable = name;
  }

}
