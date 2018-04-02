package org.genvisis.cnv.var;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.LocusSet;

public class MosaicRegion extends CNVariant implements Serializable {

  /**
   *
   */
  private static final long serialVersionUID = 1L;
  public static final String[] ADD_HEADER = new String[] {"bpWeightedScore", "nearestStatScore",
                                                          "pdfScore", "delta", "f", "customF",
                                                          "numFMarkers", "numCustomFMarkers",
                                                          "beastScore", "beastHeight",
                                                          "beastLength"};
  private final double nearestStateScore;
  private final double bpWeightedScore;
  private final double pdfScore;
  private final double delta;
  private double f;
  private final double customF;
  private int numFMarkers;
  private int numCustomFMarkers;
  private int beastLength;
  private double beastHeight;
  private double beastScore;

  public MosaicRegion(CNVariant cnv, double bpWeightedScore, double nearestStateScore,
                      double pdfScore, double delta, double f, double customF) {
    super(cnv);
    this.bpWeightedScore = bpWeightedScore;
    this.nearestStateScore = nearestStateScore;
    this.pdfScore = pdfScore;
    this.delta = delta;
    this.f = f;
    this.customF = customF;
    numFMarkers = 0;
    numCustomFMarkers = 0;
    beastScore = 0;
  }

  public MosaicRegion(CNVariant cnv, MosaicRegion another) {
    super(cnv);
    bpWeightedScore = another.bpWeightedScore;
    nearestStateScore = another.nearestStateScore;
    pdfScore = another.pdfScore;
    delta = another.delta;
    f = another.f;
    customF = another.customF;
  }

  public void setF(double f) {
    this.f = f;
  }

  public void setBeastLength(int beastLength) {
    this.beastLength = beastLength;
  }

  public void setBeastHeight(double beastHeight) {
    this.beastHeight = beastHeight;
  }

  public void setBeastScore(double beastScore) {
    this.beastScore = beastScore;
  }

  public void setNumFMarkers(int numFMarkers) {
    this.numFMarkers = numFMarkers;
  }

  public void setNumCustomFMarkers(int numCustomFMarkers) {
    this.numCustomFMarkers = numCustomFMarkers;
  }

  @Override
  public String toAnalysisString() {
    ArrayList<String> tmp = new ArrayList<>();
    tmp.add(bpWeightedScore + "");
    tmp.add(nearestStateScore + "");
    tmp.add(pdfScore + "");
    tmp.add(delta + "");
    tmp.add(f + "");
    tmp.add(customF + "");
    tmp.add(numFMarkers + "");
    tmp.add(numCustomFMarkers + "");
    tmp.add(beastScore + "");
    tmp.add(beastHeight + "");
    tmp.add(beastLength + "");

    String[] s = ArrayUtils.concatAll(toPlinkFormat().split("\t"), ArrayUtils.toStringArray(tmp));
    return ArrayUtils.toStr(s);
  }

  private static MosaicRegion fromAnalysisString(String[] data, int[] indices) {
    CNVariant tmp = new CNVariant(data);
    double bpWeightedScore = Double.parseDouble(data[indices[CNVariant.PLINK_CNV_HEADER.length]]);
    double nearestStateScore = Double.parseDouble(data[indices[CNVariant.PLINK_CNV_HEADER.length
                                                               + 1]]);
    double pdfScore = Double.parseDouble(data[indices[CNVariant.PLINK_CNV_HEADER.length + 2]]);
    double delta = Double.parseDouble(data[indices[CNVariant.PLINK_CNV_HEADER.length + 3]]);
    double f = Double.parseDouble(data[indices[CNVariant.PLINK_CNV_HEADER.length + 4]]);
    double customF = Double.parseDouble(data[indices[CNVariant.PLINK_CNV_HEADER.length + 5]]);
    int numFMarkers = Integer.parseInt(data[indices[CNVariant.PLINK_CNV_HEADER.length + 6]]);
    int numCustomFMarkers = Integer.parseInt(data[indices[CNVariant.PLINK_CNV_HEADER.length + 7]]);
    double beastScore = Double.parseDouble(data[indices[CNVariant.PLINK_CNV_HEADER.length + 8]]);
    double beastHeight = Double.parseDouble(data[indices[CNVariant.PLINK_CNV_HEADER.length + 9]]);
    int beastLength = Integer.parseInt(data[indices[CNVariant.PLINK_CNV_HEADER.length + 10]]);

    MosaicRegion mosaicRegion = new MosaicRegion(tmp, bpWeightedScore, nearestStateScore, pdfScore,
                                                 delta, f, customF);

    mosaicRegion.f = f;
    mosaicRegion.numFMarkers = numFMarkers;
    mosaicRegion.numCustomFMarkers = numCustomFMarkers;
    mosaicRegion.beastScore = beastScore;
    mosaicRegion.beastHeight = beastHeight;
    mosaicRegion.beastLength = beastLength;
    return mosaicRegion;
  }

  private static boolean validHeader(String[] header) {
    return ArrayUtils.countIf(getIndices(header), -1) == 0;
  }

  private static int[] getIndices(String[] header) {
    // TODO update with new index factors
    List<String> headerL = Arrays.asList(header);
    List<String> tmp = Arrays.asList(getBaseHeader());
    int[] indices = new int[tmp.size()];
    for (int i = 0; i < indices.length; i++) {
      indices[i] = headerL.indexOf(tmp.get(i));
    }
    return indices;
  }

  private static String[] getBaseHeader() {
    return ArrayUtils.concatAll(PLINK_CNV_HEADER, ADD_HEADER);
  }

  @Override
  public String[] getHeader() {
    return getBaseHeader();
  }

  public double getCustomF() {
    return customF;
  }

  public double getPdfScore() {
    return pdfScore;
  }

  public double getNearestStateScore() {
    return nearestStateScore;
  }

  public double getBpWeightedScore() {
    return bpWeightedScore;
  }

  public static LocusSet<MosaicRegion> loadMosLocSet(String file, Logger log) {
    List<MosaicRegion> regions = loadMosFile(file, log);
    return new LocusSet<>(regions.toArray(new MosaicRegion[regions.size()]), true, log);
  }

  private static List<MosaicRegion> loadMosFile(String filename, Logger log) {

    List<MosaicRegion> v = new ArrayList<>();
    try {
      BufferedReader reader = Files.getReader(filename, true, true);
      String[] header = reader.readLine().trim().split("\t");
      if (!validHeader(header)) {
        reader.close();
        throw new IllegalArgumentException("Invalid input file " + filename
                                           + ", expecting header of "
                                           + ArrayUtils.toStr(getBaseHeader()));
      }
      int[] indices = getIndices(header);

      while (reader.ready()) {
        String[] line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        MosaicRegion var = MosaicRegion.fromAnalysisString(line, indices);
        v.add(var);

      }
      reader.close();

      return v;
    } catch (FileNotFoundException fnfe) {
      log.reportError("file \"" + filename + "\" not found in current directory");
      log.reportException(fnfe);
    } catch (IOException ioe) {
      log.reportError("reading file \"" + filename + "\"");
      log.reportException(ioe);
    }

    return null;
  }

}
