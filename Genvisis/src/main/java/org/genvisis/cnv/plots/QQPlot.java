// -Xms1024M -Xmx1024M
package org.genvisis.cnv.plots;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.GraphicsEnvironment;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;
import org.genvisis.common.Aliases;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF.Colors.BLUES;
import org.genvisis.common.PSF.Colors.GREENS;
import org.genvisis.common.PSF.Colors.ORANGES;
import org.genvisis.common.PSF.Colors.REDS;
import org.genvisis.common.PSF.Colors.VIOLETS;
import org.genvisis.common.PSF.Colors.YELLOWS;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;
import org.genvisis.filesys.SerialFloatArray;

public class QQPlot {

  public static final long serialVersionUID = 1L;

  public static final String[] DEFAULT_FILES = {
                                                // "C:\\Documents and Settings\\npankrat\\My
                                                // Documents\\Downloads\\allelic.txt"
                                                // "C:\\Documents and Settings\\npankrat\\My
                                                // Documents\\LOAD\\QQplots\\dominant.txt"
                                                "C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\QQplots\\no_covariates.txt,0",
                                                // "C:\\Documents and Settings\\npankrat\\My
                                                // Documents\\LOAD\\QQplots\\4PCs.txt",
                                                "C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\QQplots\\E2_E4.txt,0",
      // "C:\\Documents and Settings\\npankrat\\My
      // Documents\\LOAD\\QQplots\\4PCs_E2_E4.txt",
      // "C:\\Documents and Settings\\npankrat\\My
      // Documents\\LOAD\\QQplots\\E4_binary.txt"
  };

  public static final Color[] COLOR_SCHEME = new Color[] {Color.BLACK, Color.GRAY,
                                                          BLUES.MIDNIGHT_EXPRESS, // dark dark
                                                          BLUES.PERSIAN_BLUE, // dark blue
                                                          REDS.VENETIAN_RED, // deep red
                                                          VIOLETS.BLUE_VIOLET, // deep purple
                                                          GREENS.GREEN, // dark green
                                                          BLUES.DODGER_BLUE, // light blue
                                                          BLUES.SLATE_BLUE, // light purple
                                                          GREENS.GREEN_YELLOW, // light green
                                                          VIOLETS.ORCHID, // pink
                                                          BLUES.NAVY, BLUES.CORNFLOWER_BLUE,
                                                          BLUES.DARK_SLATE_BLUE, BLUES.SLATE_BLUE,
                                                          BLUES.MEDIUM_SLATE_BLUE,
                                                          BLUES.LIGHT_SLATE_BLUE, BLUES.MEDIUM_BLUE,
                                                          BLUES.ROYAL_BLUE, Color.BLUE,
                                                          BLUES.DODGER_BLUE, BLUES.DEEP_SKY_BLUE,
                                                          BLUES.LIGHT_SKY_BLUE,
                                                          BLUES.LIGHT_SKY_BLUE, BLUES.STEEL_BLUE,
                                                          BLUES.LIGHT_STEEL_BLUE, BLUES.LIGHT_BLUE,
                                                          BLUES.POWDER_BLUE, BLUES.PALE_TURQUOISE,
                                                          BLUES.DARK_TURQUOISE,
                                                          BLUES.MEDIUM_TURQUOISE, BLUES.TURQUOISE,
                                                          BLUES.AQUA, BLUES.LIGHT_CYAN,
                                                          YELLOWS.AMBER, ORANGES.MANGO_TANGO,};

  private final QQPanel qqPanel;
  private final Logger log;

  /**
   * @param labels labels for each set of pvals
   * @param pvals array of pvals for each label in labels
   * @param log10 use -log10 of p values
   * @param rotated display QQ rotated
   * @param symmetric force symmetric axes
   * @param maxValue maximum -log10 pval to display
   * @param log
   */
  public QQPlot(String[] labels, double[][] pvals, boolean[][] pvalsToUseForLambdaCalc,
                boolean log10, boolean rotated, boolean symmetric, float maxValue, Logger log) {
    this.log = log;
    log.report("Loading data for " + ext.listWithCommas(labels));

    qqPanel = new QQPanel(labels, pvals, pvalsToUseForLambdaCalc, log10, rotated, maxValue,
                          COLOR_SCHEME, log);
    qqPanel.setSymmetricAxes(symmetric);
    if (!log10) {
      qqPanel.setForcePlotXmax(1);
      qqPanel.setForcePlotYmax(1);
    }

    qqPanel.setTrackedSize((int) qqPanel.getSize().getWidth(), (int) qqPanel.getSize().getHeight());
  }

  public void screenCap(String outFile) {
    int descrHeight = qqPanel.getPvals().length * 35;
    qqPanel.setSize(new Dimension(2000, 1440 + descrHeight));
    qqPanel.validate();
    SwingUtilities.invokeLater(new Runnable() {

      @Override
      public void run() {
        qqPanel.screenCapture(outFile);
      }
    });
  }

  protected QQPanel getQqPanel() {
    return qqPanel;
  }

  protected Logger getLog() {
    return log;
  }

  /**
   * @param filenames files to load pvals from
   * @param plotLabel label for the plot
   * @param displayQuantiles display quantiles instead of -log10 pvals
   * @param displayStandardQQ display standard QQ plot with -log10 pvals
   * @param displayRotatedQQ display rotated QQ plot with -log10 pvals
   * @param maxToPlot -log10(p) at which to start truncating
   * @param symmetric force symmetric axes
   * @param maxValue maximum -log10 p-value to plot
   * @param log
   */
  public static QQPlot loadPvals(String[] filenames, String plotLabel, boolean displayQuantiles,
                                 boolean displayStandardQQ, boolean displayRotatedQQ,
                                 double maxToPlot, double mafLwrBnd, boolean symmetric,
                                 float maxValue, Logger log) {
    BufferedReader reader;
    String[] labels;
    double[][] pvals;
    boolean[][] usePvalForLambda;
    int count;
    String fullLine;
    String[] lineParts;
    String trav;
    boolean error, header;
    int[] cols;
    int[] mafCols;
    double minPval;
    String delimiter;
    int invalids;

    if (filenames == null || filenames.length == 0) {
      reportQQError("There are no files selected for viewing in QQ plot; please add the names of the files you want to view after QQ_FILENAMES= in the project's .properties file (also be sure to uncomment the property by removing the hash symbol (\"#\"))",
                    log);
      return null;
    }

    if (maxToPlot < 0) {
      minPval = -1;
    } else {
      minPval = 1 / Math.pow(10, maxToPlot);
    }

    error = false;
    labels = new String[filenames.length];
    pvals = new double[filenames.length][];
    usePvalForLambda = new boolean[filenames.length][];
    cols = ArrayUtils.intArray(filenames.length, 0);
    mafCols = ArrayUtils.intArray(filenames.length, -1);

    for (int i = 0; i < filenames.length; i++) {
      if (filenames[i].indexOf("=") > 0) {
        labels[i] = filenames[i].substring(filenames[i].lastIndexOf("=") + 1);
        filenames[i] = filenames[i].substring(0, filenames[i].lastIndexOf("="));
      } else {
        labels[i] = ext.rootOf(filenames[i]);
      }
      if (filenames[i].indexOf(",") > 0) {
        try {
          String[] inds = filenames[i].substring(filenames[i].indexOf(',') + 1).split(",");
          if (inds.length == 1) {
            cols[i] = Integer.parseInt(inds[0]);
            filenames[i] = filenames[i].substring(0, filenames[i].lastIndexOf(","));
            labels[i] += " '" + Files.getHeaderOfFile(filenames[i], log)[cols[i]] + "'";
            if (mafLwrBnd > 0) {
              String[] v = Files.getHeaderOfFile(filenames[i], log);
              int[] poss = ext.indexFactors(Aliases.ALLELE_FREQS, v, false);
              for (int p : poss) {
                if (p >= 0) {
                  mafCols[i] = p;
                  break;
                }
              }
            }
          } else if (inds.length == 2) {
            cols[i] = Integer.parseInt(inds[0]);
            if (mafLwrBnd > 0) {
              mafCols[i] = Integer.parseInt(inds[1]);
            }
            filenames[i] = filenames[i].substring(0, filenames[i].lastIndexOf(","));
            labels[i] += " '" + Files.getHeaderOfFile(filenames[i], log)[cols[i]] + "'";
          }
        } catch (Exception e) {}
      } else {
        String[] v = Files.getHeaderOfFile(filenames[i], log);
        int[] poss = ext.indexFactors(Aliases.PVALUES, v, false, log, false);
        for (int p : poss) {
          if (p >= 0) {
            cols[i] = p;
            labels[i] += " '" + v[cols[i]] + "'";
            break;
          }
        }
        if (mafLwrBnd > 0) {
          poss = ext.indexFactors(Aliases.ALLELE_FREQS, v, false);
          for (int p : poss) {
            if (p >= 0) {
              mafCols[i] = p;
              break;
            }
          }
        }
      }
      log.report("Loading " + labels[i]);
      if (!Files.exists(filenames[i])) {
        reportQQError("Error - file " + filenames[i] + " does not exist", log);
        return null;
      }
      delimiter = Files.determineDelimiter(filenames[i], log);

      int lineCount = Files.countLines(filenames[i], 0);
      String firstLine = Files.getFirstLineOfFile(filenames[i], log);

      try {
        if (delimiter.equals(",")) {
          trav = ext.splitCommasIntelligently(firstLine, true, log)[cols[i]];
        } else {
          trav = firstLine.trim().split(delimiter, -1)[cols[i]];
        }
      } catch (Exception e) {
        log.reportError("Error - could not parse " + filenames[i] + " completely:");
        log.reportError(firstLine);
        log.reportException(e);
        return null;
      }
      try {
        if (!ext.isMissingValue(trav)) {
          Double.parseDouble(trav); // see if first value is a valid double
          header = false;
        } else {
          header = false;
        }
      } catch (NumberFormatException nfe) {
        header = true;
        lineCount--;
      }

      invalids = 0;
      double[] pv = new double[lineCount];
      boolean[] use = ArrayUtils.booleanArray(lineCount, true);
      boolean commaDelim = delimiter.equals(",");

      try {
        reader = Files.getReader(filenames[i], true, true);
        if (header) {
          reader.readLine();
        }

        count = 0;
        while ((fullLine = reader.readLine()) != null) {
          lineParts = commaDelim ? ext.splitCommasIntelligently(fullLine, true, log)
                                 : fullLine.trim().split(delimiter, -1);
          trav = lineParts[cols[i]];
          if (!ext.isMissingValue(trav)) {
            try {
              pv[count] = Double.parseDouble(trav);
              if (pv[count] > 0) {
                if (minPval > 0 && pv[count] < minPval) {
                  pv[count] = minPval;
                }
                if (mafLwrBnd > 0) {
                  if (mafCols[i] >= 0) {
                    trav = lineParts[mafCols[i]];
                    if (!ext.isMissingValue(trav)) {
                      try {
                        use[count] = Double.parseDouble(trav) > mafLwrBnd;
                      } catch (NumberFormatException nfe) {
                        use[count] = false;
                      }
                    }
                  }
                }
                count++;
              } else {
                if (invalids < 3) {
                  reportQQError("Error - one of the p-values in file " + filenames[i]
                                + " is near zero (" + trav + ") for line:\n"
                                + ext.replaceAllWith(fullLine, delimiter, "  "), log);
                }
                invalids++;
              }
            } catch (NumberFormatException nfe) {
              if (invalids < 3) {
                reportQQError("Error - one of the p-values in file " + filenames[i]
                              + " is not a number (" + trav + ") for line:\n"
                              + ext.replaceAllWith(fullLine, delimiter, "  "), log);
              }
              invalids++;
            }
          }
        }
        reader.close();

        if (invalids > 2) {
          reportQQError("There were " + invalids
                        + " total markers that had an invalid p-value for file " + filenames[i],
                        log);
        }
        log.report("Loaded " + count + " lines of data");

        if (count < lineCount) {
          pvals[i] = new double[count];
          usePvalForLambda[i] = new boolean[count];
          System.arraycopy(pv, 0, pvals[i], 0, count);
          System.arraycopy(use, 0, usePvalForLambda[i], 0, count);
          pv = null;
          use = null;
        } else {
          pvals[i] = pv;
          usePvalForLambda[i] = use;
        }

        int[] inds = Sort.getSortedIndices(pvals[i]);
        pvals[i] = Sort.getOrdered(pvals[i], inds);
        usePvalForLambda[i] = Sort.getOrdered(usePvalForLambda[i], inds);
      } catch (FileNotFoundException fnfe) {
        log.reportError("Error - missing file: \"" + filenames[i] + "\"");
        error = false;
      } catch (IOException ioe) {
        log.reportError("Error reading file \"" + filenames[i] + "\"");
        error = false;
      }
    }

    if (error) {
      return null;
    }

    if (displayQuantiles) {
      return new QQPlot(labels, pvals, usePvalForLambda, false, false, false, maxValue, log);
    }
    if (displayStandardQQ) {
      return new QQPlot(labels, pvals, usePvalForLambda, true, false, symmetric, maxValue, log);
    }
    if (displayRotatedQQ) {
      return new QQPlot(labels, pvals, usePvalForLambda, true, true, false, maxValue, log);
    }
    return null;
  }

  private static void reportQQError(String error, Logger log) {
    if (!GraphicsEnvironment.isHeadless()) {
      JOptionPane.showMessageDialog(null, error, "Error", JOptionPane.ERROR_MESSAGE);
    } else {
      log.reportError(error);
    }
  }

  public static void computeCI(String dir, String prefix, int max, Logger log) {
    PrintWriter writer;
    int count, length;
    float[] array;
    float[][] dists;

    if (max == -1) {
      max = Integer.MAX_VALUE;
    }
    count = 1;
    while (new File(dir + prefix + "." + count + ".results").exists() && count < max) {
      count++;
    }
    count--;
    log.report("Found " + (count == 0 ? "no" : count + "") + " replicates to process");
    dists = new float[0][0];
    length = -1;
    for (int i = 1; i <= count; i++) {
      if (i % 10 == 0) {
        log.report(i + "");
      }
      array = SerialFloatArray.load(dir + prefix + "." + i + ".results").getArray();
      Arrays.sort(array);
      if (i == 1) {
        length = array.length;
        dists = new float[length][count];
      } else if (array.length != length) {
        log.reportError("Error - mismatched number of p-values at rep " + i);
      }
      for (int j = 0; j < array.length; j++) {
        dists[j][i - 1] = array[j];
        if (j == array.length - 1 && Float.isNaN(dists[j][i - 1])) {
          log.reportError("Error - rep " + i + " has NaNs");
        }
      }
    }

    try {
      writer = Files.openAppropriateWriter(dir + "CIs.xln");
      writer.println("2.5%ile\t5%ile\t95%ile\t97.5%ile");
      for (float[] dist : dists) {
        writer.println(ArrayUtils.toStr(ArrayUtils.quants(dist, new double[] {0.025, 0.05, 0.95,
                                                                              0.975})));
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + "CIs.xln");
      log.reportException(e);
    }
  }

  public static void main(String[] args) {
    QQPlotFrame.main(args);
  }
}
