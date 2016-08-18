package org.genvisis.cnv.qc;

import com.google.common.primitives.Doubles;

import java.util.ArrayList;

import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.MarkerDataLoader;
import org.genvisis.common.Array;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;

/**
 * Class to compute a QC metric based on B-allele MINOR allele frequency (BMAF) "min(medianBAF,
 * 1-medianBAF)
 *
 */
public class CNVBMAF extends CNVBDeviation {

  public static final double DEFUALT_BAF_HET_PENALTY = 0;
  public static final double DEFUALT_GENO_HET_PENALTY = 0;

  public static final double DEFUALT_BAF_HET_LOWER = .15;
  public static final double DEFUALT_BAF_HET_UPPER = .85;

  /**
   * Metric is divided by all markers
   */
  public static final String SUMMARY_BY_NUM_ALL_MARKERS = "SUMMARY_BY_NUM_ALL_MARKERS";
  /**
   * Metric is divided by all markers, and the average heterzygous bmaf is subtracted from the score
   */
  public static final String SUMMARY_BY_NUM_ALL_MARKERS_AVG_HET_PENALTY =
      "SUMMARY_BY_NUM_ALL_MARKERS_AVG_HET_PENALTY";

  /**
   * Metric is divided by sum of all markers
   */
  public static final String SUMMARY_BY_SUM_BMAF_ALL_MARKERS = "SUMMARY_BY_BMAF_ALL_MARKERS";

  private final ArrayList<Double> bmafAll;// keeping as ArrayList for now in case we want to use
                                          // medians later
  private final ArrayList<Double> bmafNonHet;
  private final ArrayList<Double> bmafHet;

  private int countsHetGeno, countsHetBAF;
  private double bmafMetric, percentHet;

  public CNVBMAF(String[] intensityOnlyFlags, double gcThreshold) {
    super(intensityOnlyFlags, gcThreshold);
    bmafAll = new ArrayList<Double>();
    bmafNonHet = new ArrayList<Double>();
    bmafHet = new ArrayList<Double>();
    countsHetGeno = 0;
    countsHetBAF = 0;
    bmafMetric = 0;// TODO, could set to NaN instead
    percentHet = 0;
  }

  /**
   * 
   * @param markerName
   * @param genotype
   * @param baf
   * @param bmaf is B-allele MINOR allele frequency (BMAF) min(medianBAF, 1-medianBAF)
   * @param gc
   */
  public void add(String markerName, byte genotype, float baf, double bmaf, float gc) {
    if (super.add(markerName, genotype, baf, gc)) {
      // && !isBafHet(baf) //TODO, could add this back in, but seems to give better separation if
      // not there
      if (genotype != 1) {
        bmafNonHet.add(bmaf);
      } else {
        bmafHet.add(bmaf);
        if (genotype == 1) {
          countsHetGeno++;
        }
        if (isBafHet(baf)) {// TODO exclusive counts, or overlapping?
          countsHetBAF++;
        }
      }
      bmafAll.add(bmaf);
    }
  }

  /**
   * @param genoHetPenalty the number of heterozygous genotypes times this number will be subtracted
   *        from the final metric. Set to 0 if no penalty is desired
   * @param bafHetPenalty the number of heterozygous bafs times this number will be subtracted from
   *        the final metric. Set to 0 if no penalty is desired
   * @param summaryType see {@link #SUMMARY_BY_SUM_BMAF_ALL_MARKERS} and
   *        {@link #SUMMARY_BY_NUM_ALL_MARKERS} for options
   * @param log
   */
  public void summarize(double genoHetPenalty, double bafHetPenalty, String summaryType,
                        Logger log) {
    super.summarize();
    if (bmafNonHet.size() > 0) {
      bmafMetric = Array.sum(Doubles.toArray(bmafNonHet));
      if (summaryType.equals(SUMMARY_BY_NUM_ALL_MARKERS)) {
        bmafMetric /= bmafAll.size();
      } else if (summaryType.equals(SUMMARY_BY_NUM_ALL_MARKERS_AVG_HET_PENALTY)) {
        bmafMetric /= bmafAll.size();
        if (bmafHet.size() > 0) {
          bmafMetric -= Array.mean(Doubles.toArray(bmafHet));
        }
      } else if (summaryType.equals(SUMMARY_BY_SUM_BMAF_ALL_MARKERS)) {
        bmafMetric /= Array.sum(Doubles.toArray(bmafAll));
      } else {
        log.reportError("Error - internal error, invalid summary type");
        bmafMetric = Double.NaN;
      }
    }
    if (bmafAll.size() > 0) {// could have been all copy number probes, filtered by gc, etc..
      percentHet = (double) countsHetGeno / bmafAll.size();
    }
    bmafMetric -= countsHetGeno * genoHetPenalty;
    bmafMetric -= countsHetBAF * bafHetPenalty;

  }

  public double getBmafMetric() {
    return bmafMetric;
  }

  public double getPercentHet() {
    return percentHet;
  }

  private static boolean isBafHet(float baf) {
    return baf > DEFUALT_BAF_HET_LOWER && baf < DEFUALT_BAF_HET_UPPER;
  }

  /**
   * Helper class to facilitate a marker data based computation across all samples
   *
   */
  public static class PoplulationBAFs {
    private final CNVBMAF[] cnvbmafs;

    public PoplulationBAFs(int numSamples, String[] intensityOnlyFlags, double gcThreshold) {
      cnvbmafs = initPopulation(numSamples, intensityOnlyFlags, gcThreshold);
    }

    public void add(String markerName, byte[] genotypes, float[] bafs, float[] gcs) {
      float bamf = (float) Array.median(Array.toDoubleArray(Array.removeNaN(bafs)));
      bamf = Math.min(bamf, 1 - bamf);
      for (int i = 0; i < cnvbmafs.length; i++) {
        cnvbmafs[i].add(markerName, genotypes[i], bafs[i], bamf, gcs[i]);
      }
    }

    public void summarize(double genoHetPenalty, double bafHetPenalty, String summaryType,
                          Logger log) {
      for (CNVBMAF cnvbmaf : cnvbmafs) {
        cnvbmaf.summarize(genoHetPenalty, bafHetPenalty, summaryType, log);
      }
    }

    public CNVBMAF[] getCnvbmafs() {
      return cnvbmafs;
    }

    public double[] getCnvbmafsMetrics() {
      double[] cnvbmafsMetrics = new double[cnvbmafs.length];
      for (int i = 0; i < cnvbmafsMetrics.length; i++) {
        cnvbmafsMetrics[i] = cnvbmafs[i].getBmafMetric();
      }
      return cnvbmafsMetrics;
    }

  }

  private static CNVBMAF[] initPopulation(int numSamples, String[] intensityOnlyFlags,
                                          double gcThreshold) {
    CNVBMAF[] cnvbmafs = new CNVBMAF[numSamples];
    for (int i = 0; i < cnvbmafs.length; i++) {
      cnvbmafs[i] = new CNVBMAF(intensityOnlyFlags, gcThreshold);
    }
    return cnvbmafs;
  }

  public static void test() {
    Project proj = new Project(null, false);
    // String display = proj.getFilename(proj.DISPLAY_MARKERS_FILENAME);
    String display = proj.DISPLAY_MARKERS_FILENAMES.getValue()[0];
    String[] markers = HashVec.loadFileToStringArray(display, false, new int[] {0}, true);
    MarkerDataLoader markerDataLoader =
        MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markers);
    PoplulationBAFs poplulationBDeviation =
        new PoplulationBAFs(proj.getSamples().length, DEFAULT_INTENSITY_ONLY_FLAGS,
                            DEFAULT_GC_THRESHOLD);
    for (int i = 0; i < markers.length; i++) {
      MarkerData markerData = markerDataLoader.requestMarkerData(i);
      poplulationBDeviation.add(markers[i], markerData.getAbGenotypes(), markerData.getBAFs(),
                                markerData.getGCs());
      markerDataLoader.releaseIndex(i);

    }
    poplulationBDeviation.summarize(DEFUALT_GENO_HET_PENALTY, DEFUALT_BAF_HET_PENALTY,
                                    SUMMARY_BY_NUM_ALL_MARKERS, new Logger());
  }

  public static void main(String[] args) {
    test();
  }
}
