package org.genvisis.cnv.qc;

import java.util.HashSet;
import java.util.Vector;

import org.genvisis.cnv.analysis.ProjectCNVFiltering;
import org.genvisis.cnv.analysis.cnvTrio;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Array;
import org.genvisis.common.CNVFilter;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;

public class CNVTrioFilter extends CNVFilter {
  public static final double DEFAULT_MAX_BEAST_HEIGHT_PARENTS = 0.25;
  public static final double DEFAULT_MIN_BEAST_HEIGHT_DIFFERENCE = 0.5;
  public static final double DEFAULT_MAX_TRIO_LRR_SD = 0.3;
  public static final double DEFAULT_MAX_TRIO_1585_SD = 0.8;
  public static final int DEFAULT_MAX_NUM_CALLS = 100;

  public static final String COMMAND_MAX_BEAST_HEIGHT_PARENTS = "maxBeastHeightParents=";
  public static final String COMMAND_MIN_BEAST_HEIGHT_DIFFERENCE = "minBeastHeightDifference=";
  public static final String COMMAND_MAX_TRIO_LRR_SD = "maxTrioLrrSD=";
  public static final String COMMAND_MAX_TRIO_1585_SD = "maxTrio1585SD=";
  public static final String COMMAND_MAX_NUM_CALLS = "numCalls=";
  public static final String COMMAND_CNV_TRIO_CRF = "filterCNVTrios";
  public static final String COMMAND_CNV_TRIO_CRF_DESCRIPTION =
      " - filter a file of cnv trio results";

  public static final double NO_FILTER_MAX_BEAST_HEIGHT_PARENTS = 0;
  public static final double NO_FILTER_MIN_BEAST_HEIGHT_DIFFERENCE = 0;
  public static final double NO_FILTER_MAX_TRIO_LRR_SD = Double.MAX_VALUE;
  public static final double NO_FILTER_MAX_TRIO_1585_SD = Double.MAX_VALUE;
  public static final int NO_FILTER_NUM_CALLS = Integer.MAX_VALUE;

  private static boolean checkMaxThreshold(double d, double maxThreshold) {
    return d <= maxThreshold;
  }

  private static boolean checkMinThreshold(double d, double minThreshold) {
    return d > minThreshold;
  }

  private static boolean checkTripleMax(double d1, double d2, double d3, double maxThreshold) {
    return checkMaxThreshold(d1, maxThreshold) && checkMaxThreshold(d2, maxThreshold)
           && checkMaxThreshold(d3, maxThreshold);

  }

  public static void fromParametersTrio(String filename, Logger log) {
    Vector<String> params;
    params = Files.parseControlFile(filename, CNVTrioFilter.COMMAND_CNV_TRIO_CRF, getParserParams(),
                                    log);
    if (params != null) {
      main(Array.toStringArray(params));
    }
  }

  public static String[] getDefaultCNVTrioParams() {
    String[] params = new String[11];
    params[0] = "# CNV Trio based Filter";

    params[1] = "# maximum beast height for parents (and also the minimum for an offspring)";
    params[2] = COMMAND_MAX_BEAST_HEIGHT_PARENTS + DEFAULT_MAX_BEAST_HEIGHT_PARENTS;

    params[3] = "# minimum beast height difference between both parents and offspring";
    params[4] = COMMAND_MIN_BEAST_HEIGHT_DIFFERENCE + DEFAULT_MIN_BEAST_HEIGHT_DIFFERENCE;

    params[5] = "# maximum log R ratio standard deviation for all trio members";
    params[6] = COMMAND_MAX_TRIO_LRR_SD + DEFAULT_MAX_TRIO_LRR_SD;

    params[7] = "# maximum BAF 15/85 standard deviation for all trio members";
    params[8] = COMMAND_MAX_TRIO_1585_SD + DEFAULT_MAX_TRIO_1585_SD;

    params[9] = "# maximum number of cnv calls for an offspring";
    params[10] = COMMAND_MAX_NUM_CALLS + DEFAULT_MAX_NUM_CALLS;
    params = Array.concatAll(getDefaultCNVParams(), params);
    return params;
  }

  public static CNVTrioFilter setupCNVTrioFilterFromArgs(Project proj, String[] args,
                                                         boolean defaults, Logger log) {
    CNVTrioFilter filter = new CNVTrioFilter(log);
    ProjectCNVFiltering.setupCNVFilterFromArgs(proj, args, filter, defaults, log);
    if (defaults) {
      filter.setTrioDefualts();
    }
    for (String arg : args) {
      if (arg.startsWith(COMMAND_MAX_BEAST_HEIGHT_PARENTS)) {
        filter.setMaxBeastHeightParents(ext.parseDoubleArg(arg));
        filter.addCommandLineFilter(arg, COMMAND_MAX_SCORE);
      } else if (arg.startsWith(COMMAND_MIN_BEAST_HEIGHT_DIFFERENCE)) {
        filter.setMinBeastHeightDifference(ext.parseDoubleArg(arg));
        filter.addCommandLineFilter(arg, COMMAND_MIN_BEAST_HEIGHT_DIFFERENCE);
      } else if (arg.startsWith(COMMAND_MAX_TRIO_LRR_SD)) {
        filter.setMaxTrioLrrSD(ext.parseDoubleArg(arg));
        filter.addCommandLineFilter(arg, COMMAND_MAX_TRIO_LRR_SD);
      } else if (arg.startsWith(COMMAND_MAX_TRIO_1585_SD)) {
        filter.setMaxTrio1585SD(ext.parseDoubleArg(arg));
        filter.addCommandLineFilter(arg, COMMAND_MAX_TRIO_1585_SD);
      } else if (arg.startsWith(COMMAND_MAX_NUM_CALLS)) {
        filter.setNumCalls(ext.parseIntArg(arg));
        filter.addCommandLineFilter(arg, COMMAND_MAX_NUM_CALLS);
      }
    }
    return filter;
  }

  private double maxBeastHeightParents;

  private double minBeastHeightDifference;

  private double maxTrioLrrSD;

  private double maxTrio1585SD;

  private int numCalls;

  public CNVTrioFilter(double maxBeastHeightParents, double minBeastHeightDifference,
                       double maxTrioLrrSD, double maxTrio1585SD, int numCalls, int minNumMarkers,
                       int maxNumMarkers, int minSize, int maxSize, double minScore,
                       double maxScore, Segment[] problemRegions, Segment[] centromereMidpoints,
                       Segment[] commonReference, int[][] centromereBoundaries,
                       boolean breakupCentromeres, boolean commonIn, HashSet<String> indHash,
                       int build, int CN, Logger log) {
    super(minNumMarkers, maxNumMarkers, minSize, maxSize, minScore, maxScore, problemRegions,
          centromereMidpoints, commonReference, centromereBoundaries, breakupCentromeres, commonIn,
          indHash, build, CN, log);
    this.maxBeastHeightParents = maxBeastHeightParents;
    this.minBeastHeightDifference = minBeastHeightDifference;
    this.maxTrioLrrSD = maxTrioLrrSD;
    this.maxTrio1585SD = maxTrio1585SD;
    this.numCalls = numCalls;

    // TODO Auto-generated constructor stub
  }

  public CNVTrioFilter(double maxBeastHeightParents, double minBeastHeightDifference,
                       double maxTrioLrrSD, double maxTrio1585SD, int numCalls, Logger log) {
    super(log);
    this.maxBeastHeightParents = maxBeastHeightParents;
    this.minBeastHeightDifference = minBeastHeightDifference;
    this.maxTrioLrrSD = maxTrioLrrSD;
    this.maxTrio1585SD = maxTrio1585SD;
    this.numCalls = numCalls;
  }

  public CNVTrioFilter(Logger log) {
    super(log);
  }

  public CNVFilterPass getCNVTrioFilterPass(cnvTrio CNVTrio) {
    CNVFilterPass filterPass = getCNVFilterPass(CNVTrio);
    CNVTrio.computeMinHeightDist(maxBeastHeightParents);
    if (filterPass.passedFilter()) {
      if (!checkMaxThreshold(Math.abs(CNVTrio.getFAHEIGHT()), maxBeastHeightParents)) {
        // filterPass.prepFail();
        // filterPass.addReasonFailing("father's beast height > " + maxBeastHeightParents, ";");
        filterPass.setFailed("father's beast height > " + maxBeastHeightParents, ";");
      }
      if (!checkMaxThreshold(Math.abs(CNVTrio.getMOHEIGHT()), maxBeastHeightParents)) {
        // filterPass.prepFail();
        // filterPass.addReasonFailing("mother's beast height > " + maxBeastHeightParents, ";");
        filterPass.setFailed("mother's beast height > " + maxBeastHeightParents, ";");
      }
      if (!checkTripleMax(CNVTrio.getILRR_SD(), CNVTrio.getFALRR_SD(), CNVTrio.getMOLRR_SD(),
                          maxTrioLrrSD)) {
        // filterPass.prepFail();
        // filterPass.addReasonFailing("trio had LRR SD > " + maxTrioLrrSD, ";");
        filterPass.setFailed("trio had LRR SD > " + maxTrioLrrSD, ";");
      }
      if (!checkTripleMax(CNVTrio.getIBAF1585_SD(), CNVTrio.getFABAF1585_SD(),
                          CNVTrio.getMOBAF1585_SD(), maxTrio1585SD)) {
        // filterPass.prepFail();
        // filterPass.addReasonFailing("trio had BAF 15/85 SD > " + maxTrio1585SD, ";");
        filterPass.setFailed("trio had BAF 15/85 SD > " + maxTrio1585SD, ";");
      }
      if (!checkMinThreshold(CNVTrio.getMinBeastHeightDifference(), minBeastHeightDifference)) {
        // filterPass.prepFail();
        // filterPass.addReasonFailing("minimum beast height (" +
        // CNVTrio.getMinBeastHeightDifference() + ") <" + minBeastHeightDifference, ";");
        filterPass.setFailed("minimum beast height (" + CNVTrio.getMinBeastHeightDifference()
                             + ") <" + minBeastHeightDifference, ";");
      }
      if (CNVTrio.getNumCalls() > numCalls) {
        // filterPass.prepFail();
        // filterPass.addReasonFailing("number of calls > " + numCalls, ";");
        filterPass.setFailed("number of calls > " + numCalls, ";");
      }
    }
    return filterPass;
  }

  public double getMaxBeastHeightParents() {
    return maxBeastHeightParents;
  }

  public double getMaxTrio1585SD() {
    return maxTrio1585SD;
  }

  public double getMaxTrioLrrSD() {
    return maxTrioLrrSD;
  }

  public double getMinBeastHeightDifference() {
    return minBeastHeightDifference;
  }

  public int getNumCalls() {
    return numCalls;
  }

  public void setAllTrioFilters(double maxBeastHeightParents, double minBeastHeightDifference,
                                double maxTrioLrrSD, double maxTrio1585SD, int numCalls) {
    this.maxBeastHeightParents = maxBeastHeightParents;
    this.minBeastHeightDifference = minBeastHeightDifference;
    this.maxTrioLrrSD = maxTrioLrrSD;
    this.maxTrio1585SD = maxTrio1585SD;
    this.numCalls = numCalls;
  }

  public void setMaxBeastHeightParents(double maxBeastHeightParents) {
    this.maxBeastHeightParents = maxBeastHeightParents;
  }

  public void setMaxTrio1585SD(double maxTrio1585SD) {
    this.maxTrio1585SD = maxTrio1585SD;
  }

  public void setMaxTrioLrrSD(double maxTrioLrrSD) {
    this.maxTrioLrrSD = maxTrioLrrSD;
  }

  public void setMinBeastHeightDifference(double minBeastHeightDifference) {
    this.minBeastHeightDifference = minBeastHeightDifference;
  }

  public void setNumCalls(int numCalls) {
    this.numCalls = numCalls;
  }

  public void setTrioDefualts() {
    setAllTrioFilters(DEFAULT_MAX_BEAST_HEIGHT_PARENTS, DEFAULT_MIN_BEAST_HEIGHT_DIFFERENCE,
                      DEFAULT_MAX_TRIO_LRR_SD, DEFAULT_MAX_TRIO_1585_SD, DEFAULT_MAX_NUM_CALLS);
  }

}
