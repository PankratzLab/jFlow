package org.genvisis.cnv.analysis.pca;

import org.genvisis.cnv.analysis.pca.PrincipalComponentsCompute.PRE_PROCESSING_METHOD;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.stats.StatsCrossTabs;
import org.pankratzlab.common.stats.StatsCrossTabs.STAT_TYPE;

/**
 * @author lane0212 <br>
 *         From "A Multiple Testing Correction Method for Genetic Association Studies Using
 *         Correlated Single Nucleotide Polymorphisms" <br>
 *         http://www.ncbi.nlm.nih.gov/pubmed/18271029 <br>
 *         We extend this method to all double data, and unchunk the blocks <br>
 *         Test data set can be found at
 *         http://sourceforge.net/projects/simplem/files/simpleM_Ex.zip/download : should be 200 in
 *         block mode,191 in single frame mode<br>
 */
public class SimpleM {

  public enum MODE {
    /**
     * We divide the data into blocks and process each block separately
     */
    BLOCK_MODE,
    /**
     * All data is processed at once
     */
    SINGLE_FRAME_MODE;
  }

  private static int BLOCK_SIZE_FROM_PAPER = 133;
  private static double PCA_CUT_OFF_FROM_PAPER = 0.995;

  private final double[][] dataM;// has length of number of variables to test, i.e.
                                 // data[marker][genotypes for marker],
  private String[] dataTitles;
  private final double pcaCutoff;
  private final int blockSize;
  private final MODE mode;
  private final boolean verbose, valid, precorrelated;
  private final Logger log;

  private boolean verify() {
    boolean verify = true;
    if (pcaCutoff < 0 || pcaCutoff > 1) {
      log.reportError("The cutoff must be between 0 and 1 inclusive");
      verify = false;
    }
    return verify;
  }

  public int determineM() {
    if (valid) {
      int m = 0;

      log.reportTimeInfo("Determining effective M using mode " + mode + " and pca cutoff of "
                         + pcaCutoff);
      if (precorrelated) {
        log.reportTimeInfo("Treating data as a pre-correlated matrix");
      }
      switch (mode) {
        case BLOCK_MODE:
          log.reportTimeWarning("This should be parallelized, sorry for now");
          log.reportTimeInfo("Block sizes set to " + blockSize);
          m = scanBlocks();
          break;
        case SINGLE_FRAME_MODE:
          PrincipalComponentsCompute principalComponentsCompute = computePcs(dataM, dataTitles,
                                                                             precorrelated, verbose,
                                                                             log);
          m = getMeff(principalComponentsCompute, pcaCutoff);
          break;
        default:
          log.reportError("Invalid mode " + mode);
          break;
      }
      log.reportTimeInfo("Finished determining effective M using mode " + mode
                         + " and pca cutoff of " + pcaCutoff);
      log.reportTimeInfo("Total variables = " + dataM.length + " , Inferred effective M = " + m);
      return m;
    } else {
      log.reportError("Cannot determine number of tests");
      return -1;
    }
  }

  private int scanBlocks() {
    int index = 0;
    int start = 0;
    int stop = 0;
    int numTotal = dataM.length;
    int m = 0;
    while (stop < numTotal) {
      index++;
      int myDiff = numTotal - stop;
      if (myDiff <= blockSize) {
        break;
      } else {
        stop = start + index * blockSize;
        log.reportTimeInfo("Subsetting dataset from " + start + " -> " + stop);
        double[][] block = ArrayUtils.subArray(dataM, start, stop);
        String[] blockTitles = ArrayUtils.subArray(dataTitles, start, stop);
        start = stop;
        PrincipalComponentsCompute principalComponentsCompute = computePcs(block, blockTitles,
                                                                           precorrelated, verbose,
                                                                           log);
        int mTmp = getMeff(principalComponentsCompute, pcaCutoff);
        m += mTmp;

        log.reportTimeInfo("Index :" + index + " Current m = " + mTmp + " Total m = " + m);
      }
    }
    if (start < dataM.length) {
      double[][] block = ArrayUtils.subArray(dataM, start, numTotal);
      String[] blockTitles = ArrayUtils.subArray(dataTitles, start, numTotal);
      System.out.println("Subsetting from " + start + " -> " + numTotal);
      PrincipalComponentsCompute principalComponentsCompute = computePcs(block, blockTitles,
                                                                         precorrelated, verbose,
                                                                         log);
      int mTmp = getMeff(principalComponentsCompute, pcaCutoff);
      log.reportTimeInfo("Index :" + index + " Current m = " + mTmp + " Total m = " + m);
      m += mTmp;
    }
    return m;
  }

  private static int getMeff(PrincipalComponentsCompute principalComponentsCompute,
                             double pcaCutoff) {
    int m = 0;
    double[] singularValues = principalComponentsCompute.getSingularValues();
    double cut = pcaCutoff * ArrayUtils.sum(singularValues);
    double curSum = 0;
    for (double singularValue : singularValues) {
      if (curSum <= cut) {
        curSum += singularValue;
        m++;
      }
    }
    return m;
  }

  private static PrincipalComponentsCompute computePcs(double[][] data, String[] dataTitles,
                                                       boolean preCorrelated, boolean verbose,
                                                       Logger log) {
    PrincipalComponentsCompute principalComponentsCompute = null;
    if (!preCorrelated) {
      StatsCrossTabs statsCrossTabs = new StatsCrossTabs(data, null, null, dataTitles,
                                                         STAT_TYPE.PEARSON_CORREL, verbose, log);
      statsCrossTabs.computeTable(true);
      principalComponentsCompute = PrincipalComponentsCompute.getPrincipalComponents(data.length
                                                                                     - 1,
                                                                                     PRE_PROCESSING_METHOD.NONE,
                                                                                     statsCrossTabs.getStatisticTable(),
                                                                                     verbose, log);
    } else {
      principalComponentsCompute = PrincipalComponentsCompute.getPrincipalComponents(data.length
                                                                                     - 1,
                                                                                     PRE_PROCESSING_METHOD.NONE,
                                                                                     data, verbose,
                                                                                     log);

    }
    return principalComponentsCompute;
  }

  public static class Builder {

    private double pcaCutoff = PCA_CUT_OFF_FROM_PAPER;
    private int blockSize = BLOCK_SIZE_FROM_PAPER;
    private MODE mode = MODE.SINGLE_FRAME_MODE;
    private boolean verbose = true;
    private final String[] dataTitles = null;
    private boolean precorrelated;

    /**
     * @param pcaCutoff the percent variance that must be explained by the tests
     * @return
     */
    public Builder pcaCutoff(double pcaCutoff) {
      this.pcaCutoff = pcaCutoff;
      return this;
    }

    /**
     * @param blockSize split variables up a block size of this
     * @return
     */
    public Builder blockSize(int blockSize) {
      this.blockSize = blockSize;
      return this;
    }

    /**
     * @param mode see {@link MODE}
     * @return
     */
    public Builder mode(MODE mode) {
      this.mode = mode;
      return this;
    }

    public Builder verbose(boolean verbose) {
      this.verbose = verbose;
      return this;
    }

    /**
     * @param precorrelated if the data is already a correlation matrix, otherwise it will be
     *          computed
     * @return
     */
    public Builder precorrelated(boolean precorrelated) {
      this.precorrelated = precorrelated;
      return this;
    }

    public SimpleM build(double[][] dataM, Logger log) {
      return new SimpleM(this, dataM, log);
    }
  }

  private SimpleM(Builder builder, double[][] dataM, Logger log) {
    this.dataM = dataM;
    if (builder.dataTitles == null) {
      dataTitles = ArrayUtils.stringArraySequence(dataM.length, "Variable_");
    } else {
      dataTitles = builder.dataTitles;
    }
    pcaCutoff = builder.pcaCutoff;
    blockSize = builder.blockSize;
    mode = builder.mode;
    verbose = builder.verbose;
    precorrelated = builder.precorrelated;
    this.log = log;
    valid = verify();
  }

  public static void main(String[] args) {
    test();
  }

  private static void test() {
    String filename = "C:/bin/simpleM_Ex/snpSample.txt";// should be 200 in block mode,191 in single
                                                        // frame mode
    String[][] testData = HashVec.loadFileToStringMatrix(filename, false, null);
    double[][] matrix = ArrayUtils.toDoubleArrays(testData, false);
    String[] titles = new String[matrix.length];
    for (int i = 0; i < titles.length; i++) {
      titles[i] = "var" + i;
    }
    Builder builder = new Builder();
    SimpleM simpleM = builder.build(matrix, new Logger(ext.rootOf(filename, false)));
    simpleM.determineM();
  }
}