package org.genvisis.stats;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;
import java.util.StringTokenizer;
import java.util.Vector;

import org.genvisis.common.ext;

public class Anova {
  public static void main(String[] args) throws IOException {
    Anova a = null;
    if (args.length != 1) {
      System.out.println("Expecting 1 argument: filename with 2 columns - dependent variable independent variable.");
    }
    try {
      a = new Anova(args[0]);
      a.oneway();
      // System.out.println(a.report());
      System.out.println(a.summary());
    } catch (Exception e) {
      e.printStackTrace();
      System.err.println(a.dumpData());
    }
  }

  private final double[][] data; // might later give the option of using double
  // or int, does doubleness cause any problems?

  private double MSE;

  private double[] groupMeans;

  private int[] groupCounts;

  private double[] groupSums;

  private double[] groupSumSquares;

  private Vector<String> groupNames;

  private double SSb;

  private double SSw;

  private int dfN;

  private int dfD;

  private double Fstat;

  private double Fsig;

  private final double Q05_3g = 3.31; // if n > 200 then Q is fixed at 3.31 for 3
  // groups at alpha = 0.05; might create
  // lookup

  private double LSD05;

  private double LSD01;

  private double HSD05;

  private boolean analysis_complete;

  public Anova(double[][] inputData) throws IOException {
    analysis_complete = false;
    data = inputData;
    groupNames = new Vector<String>();
    for (int i = 1; i <= data.length; i++) {
      groupNames.add("Group " + i);
    }
  }

  public Anova(String filename) throws IOException {
    BufferedReader reader = null;
    String depen, indepen;
    StringTokenizer st;
    Hashtable<String, Vector<String>> hash = new Hashtable<String, Vector<String>>();
    Vector<String> v;

    analysis_complete = false;

    try {
      reader = new BufferedReader(new FileReader(filename));

      groupNames = new Vector<String>();
      while (reader.ready()) {
        st = new StringTokenizer(reader.readLine());
        depen = st.nextToken(); // in this case AOO
        indepen = st.nextToken(); // in this case APOE

        if (!hash.containsKey(indepen)) {
          groupNames.add(indepen);
          v = new Vector<String>();
          v.add(depen);
          hash.put(indepen, v);
        } else {
          v = hash.get(indepen);
          v.add(depen);
        }
      }
      reader.close();
    } catch (Exception e) {
      System.err.println("Error: could not process the file \"" + filename + "\".");
      e.printStackTrace();
      System.exit(1);
    }

    data = new double[groupNames.size()][];
    for (int i = 0; i < data.length; i++) {
      v = hash.get(groupNames.elementAt(i));
      data[i] = new double[v.size()];
      for (int j = 0; j < v.size(); j++) {
        data[i][j] = Integer.valueOf(v.elementAt(j)).intValue();
      }
    }
  }

  public void assignGroupNames(Vector<String> inputNames) throws IOException {
    groupNames = inputNames;
  }

  public String dumpData() throws IOException {
    String list = "";

    for (double[] element : data) {
      for (int j = 0; j < element.length; j++) {
        list += element[j] + "\t";
      }
      list += "\n";
    }

    return list;
  }

  public String getCounts() throws IOException {
    String list = "";

    for (int i = 0; i < data.length; i++) {
      list += "\tn" + (i + 1) + "=" + groupSums[i] + "/" + data[i].length;
    }

    return list;
  }

  public String getMeans() throws IOException {
    String list = "";

    if (!analysis_complete) {
      return "Analysis not complete, run Anova.oneway()";
    } else {
      for (int i = 0; i < data.length; i++) {
        list += "\t" + ext.formDeci(groupMeans[i], 2, true);
      }
    }

    return list;
  }

  public double getSig() throws IOException {
    return Fsig;
  }

  public String getStdDevs() throws IOException {
    String list = "";
    double[] groupSDs = new double[groupMeans.length];
    double[] sumSqDiffs = new double[groupMeans.length];

    for (int i = 0; i < data.length; i++) {
      for (int j = 0; j < data[i].length; j++) {
        sumSqDiffs[i] += (data[i][j] - groupMeans[i]) * (data[i][j] - groupMeans[i]);
      }
      groupSDs[i] = Math.sqrt(sumSqDiffs[i] / (data[i].length - 1));
    }

    if (!analysis_complete) {
      return "Analysis not complete, run Anova.oneway()";
    } else {
      for (int i = 0; i < data.length; i++) {
        list += "\t" + ext.formDeci(groupSDs[i], 2, true);
      }
    }

    return list;
  }

  public void oneway() throws IOException {
    int sumOfSums, totalN;
    double calc;

    groupSums = new double[data.length];
    groupSumSquares = new double[data.length];
    groupMeans = new double[data.length];
    groupCounts = new int[data.length];
    SSb = 0;
    SSw = 0;
    sumOfSums = 0;
    totalN = 0;

    for (int i = 0; i < data.length; i++) {
      for (int j = 0; j < data[i].length; j++) {
        groupSums[i] += data[i][j];
        groupSumSquares[i] += data[i][j] * data[i][j];
      }
      groupMeans[i] = (groupSums[i]) / (data[i].length);
      groupCounts[i] = data[i].length;
      calc = (groupSums[i]) * (groupSums[i]) / data[i].length;
      SSb += calc;
      SSw += (groupSumSquares[i]) - calc;
      sumOfSums += groupSums[i];
      totalN += data[i].length;

    }
    SSb -= ((double) sumOfSums) * ((double) sumOfSums) / (totalN);
    dfN = data.length - 1;
    dfD = totalN - data.length;
    Fstat = (SSb / dfN) / (SSw / dfD);
    MSE = SSw / dfD;
    LSD05 =
        ProbDist.TDistReverse(0.05, totalN) * Math.sqrt(2 * MSE) / Math.sqrt(totalN / data.length);
    LSD01 =
        ProbDist.TDistReverse(0.01, totalN) * Math.sqrt(2 * MSE) / Math.sqrt(totalN / data.length);
    HSD05 = Q05_3g * Math.sqrt(2 * MSE) / Math.sqrt(totalN / data.length);
    Fsig = ProbDist.FDist(Fstat, dfN, dfD);

    analysis_complete = true;
  }

  public String report() throws IOException {
    String list = "";

    if (!analysis_complete) {
      return "Analysis not complete, run Anova.oneway()";
    } else {
      for (int i = 0; i < data.length; i++) {
        list += groupNames.elementAt(i) + "\t" + ext.formDeci(groupMeans[i], 2, true) + "\n";
      }
      list += "SSb = " + ext.formDeci(SSb, 2, true) + "\n";
      list += "SSw = " + ext.formDeci(SSw, 2, true) + "\n";
      list += "f = " + ext.formDeci(Fstat, 5, true) + " with " + dfN + " df in numerator, " + dfD
              + " df in denominator\n";
      list += "p-value = " + ProbDist.FDist(Fstat, dfN, dfD) + "\n";
      list += "The LSD minimum mean difference for alpha = 0.05 is " + ext.formDeci(LSD05, 5, true)
              + "\n";
      list += "The LSD minimum mean difference for alpha = 0.01 is " + ext.formDeci(LSD01, 5, true)
              + "\n";
      list += "The HSD minimum mean difference for alpha = 0.05 is " + ext.formDeci(HSD05, 5, true)
              + "\n";
      for (int i = 0; i < data.length; i++) {
        for (int j = i + 1; j < data.length; j++) {
          list += "Mean diff b/w " + groupNames.elementAt(i) + " & " + groupNames.elementAt(j)
                  + " = " + ext.formDeci(groupMeans[i] - groupMeans[j], 2, true)
                  + ((Math.abs(groupMeans[i] - groupMeans[j]) > LSD05) ? "*" : "")
                  + ((Math.abs(groupMeans[i] - groupMeans[j]) > LSD01) ? "*" : "") + "\n";
        }
      }

    }

    return list;
  }

  public String summary() throws IOException {
    String list = "";

    if (!analysis_complete) {
      return "Analysis not complete, run Anova.oneway()";
    } else {
      list += ext.formDeci(Fstat, 5, true) + "\t" + Fsig + "\t" + dfN + "\t" + dfD;
      for (int i = 0; i < data.length; i++) {
        list += "\t" + groupCounts[i];
      }
      for (int i = 0; i < data.length; i++) {
        list += "\t" + ext.formDeci(groupMeans[i], 2, true);
      }
      list += "\t" + ext.formDeci(LSD05, 5, true);
      list += "\t" + ext.formDeci(LSD01, 5, true);
      for (int i = 0; i < data.length; i++) {
        for (int j = i + 1; j < data.length; j++) {
          list += "\t" + ext.formDeci(groupMeans[i] - groupMeans[j], 2, true);
          // list+=((Math.abs(groupMeans[i] - groupMeans[j]) >
          // LSD05)?"*":"")+((Math.abs(groupMeans[i] - groupMeans[j])
          // > LSD01)?"*":"");
        }
      }

    }

    return list;
  }
}
