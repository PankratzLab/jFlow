package org.genvisis.stats;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Serializable;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

// Some Resources used...
// http://ctj.sagepub.com/content/1/6/553.full.pdf
// http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0017238
// http://www.amsciepub.com/doi/pdf/10.2466/pr0.1966.19.1.3
// http://physiolgenomics.physiology.org/content/16/1/99
// http://www.john-uebersax.com/stat/icc.htm
// http://www.biomedcentral.com/content/pdf/1471-2105-15-312.pdf
/**
 * Class to compute the intraclass correlation coefficient. Currently only a "One-way Random" model
 * is used;
 *
 * Issue: not absolutely positive this handles inconsistent group sizes correctly//TODO
 */
public class ICC implements Serializable {
  /**
   *
   */
  private static final long serialVersionUID = 1L;
  private double[] parsedData;
  private final Logger log;
  private String[] response;
  private final String[] maskedResponses;
  private final String[] onlyTheseResponses;
  private ResponseEffect[] rowEffects;
  private double sumTotal, meanTotal, SSTotal, MSTotal, MSWithin, MSBetween, ICC;
  private int nTotal;
  private boolean valid;
  private final boolean verbose;

  /**
   * @param initData double array that will form the values used to compute the ICC
   * @param response must be same length as input data. Represents the groupings to compute the ICC.
   *        For example, to compute ICC between sexes, this array would contain 1's and 2's
   * @param maskedResponses if not null, these responses will be excluded
   * @param onlyTheseResponses if not null, only these responses will be used.
   * @param verbose
   * @param log
   *        <p>
   *        Note : must call {@link ICC#computeICC()} to compute the ICC
   */
  public ICC(double[] initData, String[] response, String[] maskedResponses,
             String[] onlyTheseResponses, boolean verbose, Logger log) {
    super();
    parsedData = initData;
    this.log = log;
    this.response = response;
    rowEffects = null;
    this.maskedResponses = maskedResponses;
    this.onlyTheseResponses = onlyTheseResponses;
    valid = verifyData(initData, response, log);
    ICC = Double.NaN;
    nTotal = 0;
    this.verbose = verbose;
    init(initData);// modifies parsed data
    populateFullStats();
  }

  public void shrink() {
    response = null;
    rowEffects = null;
    parsedData = null;
  }

  /**
   * Computes the ICC
   */
  public void computeICC() {
    computeMSBetween();
    computeMSWithin();
    ICC = (MSBetween - MSWithin) / (MSBetween + MSWithin);
  }

  public ResponseEffect[] getRowEffects() {
    return rowEffects;
  }

  public int getNumEffects(int minGroupSize) {
    int num = 0;
    for (ResponseEffect rowEffect : rowEffects) {
      if (rowEffect.getN() > minGroupSize) {
        num++;
      }
    }
    return num;
  }

  public int getNumEffects() {
    return rowEffects.length;
  }

  public double getICC() {
    return ICC;
  }

  public double getMSTotal() {
    return MSTotal;
  }

  public String[] getRowEffectDataStringFormat() {
    double[][] rowEffectData = getRowEffectData();
    String[] rowEffectDataStringFormat = new String[rowEffectData.length];
    for (int i = 0; i < rowEffectDataStringFormat.length; i++) {
      rowEffectDataStringFormat[i] = Array.toStr(rowEffectData[i]);
    }
    return rowEffectDataStringFormat;
  }

  public double[][] getRowEffectData() {
    double[][] rowEffectData = new double[rowEffects.length][];
    for (int i = 0; i < rowEffectData.length; i++) {
      rowEffectData[i] = rowEffects[i].getData();
    }
    return rowEffectData;
  }

  private void populateFullStats() {
    if (valid) {
      nTotal = parsedData.length;
      sumTotal = Array.sum(parsedData);
      meanTotal = Array.mean(parsedData, sumTotal);
      SSTotal = Array.sumSq(parsedData, meanTotal);
      MSTotal = Array.variance(parsedData, SSTotal);

    } else {
      log.reportError("Error - data is not valid, cannot populate full stats");
    }
  }

  private void computeMSBetween() {
    MSBetween = 0;
    for (ResponseEffect rowEffect : rowEffects) {
      if (rowEffect.isValid()) {
        // System.out.println("Mean Total " + meanTotal + "\tmean\t" + i + "\t" +
        // rowEffects[i].getRowMean() + "\t" + rowEffects[i].getN());
        MSBetween += rowEffect.getN() * Math.pow(rowEffect.getRowMean() - meanTotal, 2);
      }
    }
    if (verbose) {
      System.out.println("Between Rows\tSS: " + MSBetween + "\tMS: "
                         + (MSBetween / (rowEffects.length - 1)));
    }
    MSBetween /= (rowEffects.length - 1);
  }

  public double getBetweenRowsSS() {
    return MSBetween;
  }

  public double getBetweenRowsMS() {
    return MSBetween;
  }

  public double getWithinRowsSS() {
    return MSWithin;
  }

  public double getWithinRowsMS() {
    return MSWithin;
  }

  public int getnTotal() {
    return nTotal;
  }

  public void dump(String file) {
    try {
      PrintWriter writer = new PrintWriter(new FileWriter(file));
      writer.println("Group\tResponse");
      for (int i = 0; i < rowEffects.length; i++) {
        for (int j = 0; j < rowEffects[i].getData().length; j++) {
          if (rowEffects[i].getN() > 1) {
            writer.println("GROUP_" + i + "\t" + rowEffects[i].getData()[j]);
          }
        }
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + file);
      log.reportException(e);
    }
  }

  private void computeMSWithin() {
    MSWithin = 0;
    for (ResponseEffect rowEffect : rowEffects) {
      if (rowEffect.isValid()) {
        MSWithin += rowEffect.getSS();
      }
    }
    if (verbose) {
      System.out.println("Within Rows\tSS: " + MSWithin + "\tMS: "
                         + (MSWithin / (rowEffects.length)));
    }
    MSWithin /= (nTotal - rowEffects.length);
    // MSWithin /= rowEffects.length;

  }

  /**
   * Parses the data into groups defined by the responses
   */
  private void init(double[] data) {
    if (valid) {
      // we first assemble the unique class membership from the rowLabels
      ArrayList<String> uniqueClasses = new ArrayList<String>();
      int dataPointsToUse = 0;
      Hashtable<String, ArrayList<Integer>> track = new Hashtable<String, ArrayList<Integer>>();
      for (int i = 0; i < response.length; i++) {
        if (!isMasked(response[i], maskedResponses, onlyTheseResponses) && !Double.isNaN(data[i])) {
          if (!track.containsKey(response[i])) {
            track.put(response[i], new ArrayList<Integer>());
            uniqueClasses.add(response[i]);
          }
          track.get(response[i]).add(i);
          dataPointsToUse++;
        } else {
          // System.out.println(response[i] + "\t" + data[i]);
        }
      }
      rowEffects = new ResponseEffect[uniqueClasses.size()];
      parsedData = new double[dataPointsToUse];
      int numValidClasses = 0;
      int addIndex = 0;
      // Next we parse the classes to their respective row effects
      for (int i = 0; i < uniqueClasses.size(); i++) {
        String currentLabel = uniqueClasses.get(i);
        ArrayList<Integer> indices = track.get(currentLabel);
        double[] tmpdata = new double[indices.size()];
        for (int j = 0; j < tmpdata.length; j++) {
          tmpdata[j] = data[indices.get(j)];
          parsedData[addIndex] = tmpdata[j];
          addIndex++;
        }
        rowEffects[i] = new ResponseEffect(currentLabel, tmpdata);
        if ((!rowEffects[i].isValid() || rowEffects[i].getN() < 2) && verbose) {
          log.reportError("Error - data for class " + currentLabel + " containing "
                          + rowEffects[i].getN() + " "
                          + (rowEffects[i].getN() == 1 ? "data point is " : "data points are")
                          + " not valid , will not include in the ICC");
        } else {
          numValidClasses++;
        }
      }
      if (numValidClasses < 2) {
        valid = false;
        log.reportError("Error - must have at least two valid classes, cannot compute ICC");

      }
      if (parsedData.length != addIndex) {
        log.reportError("Error - could not add all data");
        valid = false;
      } else {
        if (verbose) {
          log.report("Info - detected " + numValidClasses + " valid classes for ICC computation");
          if (numValidClasses <= 5) {
            log.report(Array.toStr(uniqueClasses.toArray(new String[uniqueClasses.size()])));
          }
        }
      }
    } else {
      log.reportError("Error - data is not valid, cannot compute ICC");
    }

  }

  private static boolean verifyData(double[] data, String[] rowLabels, Logger log) {
    boolean valid = true;
    if (data == null) {
      log.reportError("Error - data were not found");
      valid = false;
    }
    if (rowLabels == null) {
      log.reportError("Error - row labels were not found");
      valid = false;
    }
    if (valid) {
      if (data.length != rowLabels.length) {
        log.reportError("Error - labels and data must have the same length");
        valid = false;
      }
    }
    return valid;
  }

  private static class ResponseEffect implements Serializable {
    /**
     *
     */
    private static final long serialVersionUID = 1L;
    private String label;
    private final double[] data;
    private double sum, rowMean, SS, MS;
    private final int n;
    private final boolean valid;

    public ResponseEffect(String label, double[] data) {
      super();
      this.data = data;
      n = data.length;
      valid = verifyData(data);
      sum = Double.NaN;
      rowMean = Double.NaN;
      SS = Double.NaN;
      MS = Double.NaN;
      popluateRowMetrics();
      get();// TODO
    }

    private void popluateRowMetrics() {
      if (valid) {
        sum = Array.sum(data);
        rowMean = Array.mean(data, sum);
        SS = Array.sumSq(data, rowMean);
        MS = Array.variance(data, SS);
        // System.out.println(SS + "\t" + MS);
      }
    }

    private static boolean verifyData(double[] data) {
      for (double element : data) {
        if (Double.isNaN(element)) {
          return false;
        }
      }
      return data.length >= 2;// skip single data points from computation
    }

    public double[] getData() {
      return data;
    }

    public String getLabel() {
      return label;
    }

    public double getSum() {
      return sum;
    }

    public double getMS() {
      return MS;
    }

    public void get() {
      getLabel();
      getSum();
      getMS();
    }

    public double getRowMean() {
      return rowMean;
    }

    public double getSS() {
      return SS;
    }

    public int getN() {
      return n;
    }

    public boolean isValid() {
      return valid;
    }
  }

  private static boolean isMasked(String response, String[] exclude, String[] include) {
    boolean masked = false;
    if (response == null) {
      masked = true;
    }
    if (!masked && exclude != null) {
      masked = ext.indexOfStr(response, exclude, true, true) >= 0;
    }
    if (!masked && include != null) {
      masked = ext.indexOfStr(response, include, true, true) < 0;
    }
    return masked;
  }

  public static void test2(String fullPathTofile) {
    String[] stuffs = HashVec.loadFileToStringArray(fullPathTofile, false, new int[] {0, 1}, false);
    double[] data = new double[stuffs.length * 2];
    String[] response = new String[stuffs.length * 2];
    int index = 0;
    for (int i = 0; i < stuffs.length; i++) {
      String[] stuff = stuffs[i].trim().split("\t");
      data[index] = Double.parseDouble(stuff[0]);
      response[index] = i + "";
      index++;
      data[index] = Double.parseDouble(stuff[1]);
      response[index] = i + "";
      index++;
    }

    ICC icc = new ICC(data, response, null, null, true, new Logger());
    icc.computeICC();
    System.out.println("ICC=" + icc.getICC());
  }

  public static void test() {
    try {
      BufferedReader in =
                        new BufferedReader(new InputStreamReader(new URL("http://www.uvm.edu/~dhowell/StatPages/More_Stuff/icc/PartnerCorr.dat").openStream()));
      int lines = 0;
      while (in.ready()) {
        lines++;
        in.readLine();
      }
      System.out.println(lines);
      in.close();
      String[] response = new String[lines * 2];
      double[] data = new double[lines * 2];
      in =
         new BufferedReader(new InputStreamReader(new URL("http://www.uvm.edu/~dhowell/StatPages/More_Stuff/icc/PartnerCorr.dat").openStream()));
      int index = 0;
      int count = 1;
      while (in.ready()) {

        String[] line = in.readLine().trim().split("[\\s]+");
        data[index] = Double.parseDouble(line[0]);
        response[index] = count + "";
        index++;
        data[index] = Double.parseDouble(line[1]);
        response[index] = count + "";
        index++;
        count++;
      }
      in.close();
      ICC icc = new ICC(data, response, null, null, true, new Logger());
      icc.computeICC();
      System.out.println("ICC=" + icc.getICC());

    } catch (MalformedURLException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }

  public static void main(String[] args) {
    test();// should be 0.8036863991984644
    // String filename = "C:/data/test/ICC/gedi.exome.matched.unrelated.PC40";
    // test2(filename);// should be 0.35725853348245573

  }
}
