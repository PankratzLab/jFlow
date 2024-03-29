package org.genvisis.fcs.sub;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.distribution.MixtureMultivariateNormalDistribution;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.distribution.fitting.MultivariateNormalMixtureExpectationMaximization;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.Pair;
import org.genvisis.fcs.FCSDataLoader;
import org.genvisis.fcs.AbstractPanel2.AXIS_SCALE;
import org.genvisis.fcs.FCSDataLoader.LOAD_STATE;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Logger;

import edu.stanford.facs.logicle.Logicle;

public class EMModel {

  MultivariateNormalMixtureExpectationMaximization mnmem;
  double[][] transformed;
  int[] clusterAssigns;

  public EMModel(MultivariateNormalMixtureExpectationMaximization mnmem, double[][] transformed) {
    this.mnmem = mnmem;
    this.transformed = transformed;
  }

  private static double count(double max) {
    double dec = max;
    double decCnt = 0;
    while (dec > 1) {
      dec = dec / 10d;
      decCnt++;
    }
    return decCnt;
  }

  private static Logicle getBiexScale(double max) {
    double w = 2; // linear decades // W=1 is weird, W=3 -> "scale() didn't
    // converge"
    double a = 0; // negative decades
    double decCnt = count(max);

    return new Logicle(max, Math.min(w, decCnt / 2), decCnt, a);
  }

  public static Logicle getLogicle(double[] data) {
    return getBiexScale(ArrayUtils.max(data));
  }

  private static double[][] transformData(double[][] data, Logicle[] mehs, Logger log) {
    log.report("Transforming data to Logicle scale");

    double[][] transformed = new double[data.length][data[0].length];
    for (int i = 0; i < mehs.length; i++) {
      for (int j = 0; j < data[0].length; j++) {
        transformed[i][j] =
            EMInitializer.DATA_SCALES[i] == AXIS_SCALE.BIEX
                ? mehs[i].scale(data[i][j])
                : data[i][j];
      }
    }
    log.report("Finished transforming data to Logicle scale");
    return transformed;
  }

  private static double[][] process(FCSDataLoader dataLoader, Logger log) {

    double[][] data = new double[EMInitializer.DATA_COLUMNS.length][];
    for (int i = 0; i < EMInitializer.DATA_COLUMNS.length; i++) {
      data[i] = dataLoader.getData(EMInitializer.DATA_COLUMNS[i], true);
    }

    Logicle[] mehs = new Logicle[EMInitializer.DATA_COLUMNS.length];
    for (int i = 0; i < mehs.length; i++) {
      mehs[i] = EMInitializer.DATA_SCALES[i] == AXIS_SCALE.BIEX ? getLogicle(data[i]) : null;
    }
    return transformData(data, mehs, log);
  }

  private static double[][] transpose(double[][] matrixData) {
    RealMatrix m = MatrixUtils.createRealMatrix(matrixData);
    return m.transpose().getData();
  }

  private static double[][] selectMaxVarInitialMix(
      double[][] transposedfullData, double percentage, Logger log) {
    int num = (int) Math.round(transposedfullData.length * percentage);
    log.report("Selecting " + num + " initial points");
    int cnt = transposedfullData[0].length;
    int numEach = num / (2 * cnt);
    RealMatrix m = MatrixUtils.createRealMatrix((2 * cnt) * numEach, cnt);

    for (int i = 0; i < cnt; i++) {
      sort(transposedfullData, i);
      int countMin = 0, countMax = 0;
      while (countMin < numEach) {
        m.setRow(countMin + (numEach * (2 * i)), transposedfullData[countMin]);
        countMin++;
      }
      while (countMax < numEach) {
        m.setRow(
            countMax + (numEach * ((2 * i) + 1)),
            transposedfullData[transposedfullData.length - countMax - 1]);
        countMax++;
      }
    }

    log.report("Finished selecting " + num + " initial points");

    return m.getData();
  }

  private static void sort(double[][] transposedfullData, final int sortInd) {
    Arrays.sort(
        transposedfullData,
        new Comparator<double[]>() {

          @Override
          public int compare(double[] o1, double[] o2) {
            return Double.valueOf(o1[sortInd]).compareTo(o2[sortInd]);
          }
        });
  }

  private static MultivariateNormalMixtureExpectationMaximization getClusters(
      double[][] transposed, double[][] init, int clusters, Logger log) {
    log.report("Detecting clusters, num clusters set to " + clusters);

    MixtureMultivariateNormalDistribution initialMix =
        MultivariateNormalMixtureExpectationMaximization.estimate(init, clusters);

    List<Pair<Double, MultivariateNormalDistribution>> clusts = initialMix.getComponents();
    int numSamples = 1000;
    double[][] newInit = new double[clusters * numSamples][];
    int index = 0;
    Random rand = new Random();
    Logicle l = getBiexScale(10000);
    int noise = 1000;
    int noiseHalf = noise / 2;

    for (Pair<Double, MultivariateNormalDistribution> pair : clusts) {
      for (int i = 0; i < numSamples; i++) {
        if (i < numSamples / 4) {
          newInit[index] =
              new double[] {
                l.scale(0 + rand.nextInt(noise) - noiseHalf),
                l.scale(0 + rand.nextInt(noise) - noiseHalf)
              };
        } else if (i < numSamples / 2) {
          newInit[index] =
              new double[] {
                l.scale(10000 + rand.nextInt(noise) - noiseHalf),
                l.scale(0 + rand.nextInt(noise) - noiseHalf)
              };
        } else if (i > numSamples / 2 && i < (numSamples / 2 + numSamples / 4)) {
          newInit[index] =
              new double[] {
                l.scale(10000 + rand.nextInt(noise) - noiseHalf),
                l.scale(10000 + rand.nextInt(noise) - noiseHalf)
              };
        } else if (i > numSamples / 2) {
          newInit[index] =
              new double[] {
                l.scale(0 + rand.nextInt(noise) - noiseHalf),
                l.scale(10000 + rand.nextInt(noise) - noiseHalf)
              };
        } else {
          newInit[index] = pair.getValue().sample();
        }
        index++;
      }
    }
    initialMix = MultivariateNormalMixtureExpectationMaximization.estimate(newInit, clusters);

    return fit(transposed, initialMix);
  }

  private static MultivariateNormalMixtureExpectationMaximization fit(
      double[][] transposed, MixtureMultivariateNormalDistribution initialMix) {
    MultivariateNormalMixtureExpectationMaximization mle =
        new MultivariateNormalMixtureExpectationMaximization(transposed);

    mle.fit(initialMix, 2000, 1E-5);

    return mle;
  }

  public static EMModel run(FCSDataLoader loader) {
    if (loader.getLoadState() != LOAD_STATE.LOADED) {
      System.err.println("Error - data not loaded yet; please wait a while and try again.");
      return null;
    }
    Logger log = new Logger();
    long t1 = System.currentTimeMillis();
    double[][] transformed = process(loader, log);
    double[][] transposed = transpose(transformed);

    MultivariateNormalMixtureExpectationMaximization mnmem =
        fit(transposed, EMInitializer.load(EMInitializer.DIST_DIR));

    // double[][] initData = selectMaxVarInitialMix(transposed, .1, log);
    //
    // MultivariateNormalMixtureExpectationMaximization mnmem = null;
    //
    // try {
    // MultivariateNormalMixtureExpectationMaximization mnmem1 = getClusters(transposed, initData,
    // 4, log);
    // mnmem = mnmem1;
    // } catch (Exception e) {
    // System.err.println("Error - 4 clusters failed: " + e.getMessage());
    // for (int clust = 6; clust > 1; clust--) {
    // try {
    // MultivariateNormalMixtureExpectationMaximization mnmem1 = getClusters(transposed, initData,
    // clust, log);
    // mnmem = mnmem1;
    // break;
    // } catch (Exception e1) {
    // System.out.println("Clustering " + clust + " -- " + e1.getMessage());
    // }
    // }
    // }
    //
    // new Thread(new Runnable() {
    // @Override
    // public void run() {
    // for (int clust = 2; clust < 4; clust++) {
    // MultivariateNormalMixtureExpectationMaximization mnmem1 = getClusters(transposed, initData,
    // clust, log);
    // System.out.println("LogL - " + clust + ": " + mnmem1.getLogLikelihood());
    // }
    // }
    // }).run();
    //
    // System.out.println("LogL - 4: " + mnmem.getLogLikelihood());
    //
    log.reportTimeElapsed("Clustered in ", t1);

    return new EMModel(mnmem, transformed);
  }

  public int[] getClusterAssigns() {
    if (clusterAssigns == null) {
      MixtureMultivariateNormalDistribution mmnd = mnmem.getFittedModel();

      List<Pair<Double, MultivariateNormalDistribution>> clusts = mmnd.getComponents();

      double[][] means = new double[clusts.size()][];
      double[][] sds = new double[clusts.size()][];

      clusterAssigns = new int[transformed[0].length];
      int clust = 0;
      for (Pair<Double, MultivariateNormalDistribution> pair : clusts) {
        means[clust] = pair.getValue().getMeans();
        sds[clust] = pair.getValue().getStandardDeviations();
        clust++;
      }

      for (int i = 0; i < transformed[0].length; i++) {
        double[] line = new double[transformed.length];
        for (int l = 0; l < transformed.length; l++) {
          line[l] = transformed[l][i];
        }
        double curMax = .01;
        ArrayList<Integer> clustsBelong = new ArrayList<>();
        for (int j = 0; j < means.length; j++) {
          // determine if coords belongs to this cluster
          double tmp = clusts.get(j).getValue().density(line);
          if (tmp > curMax) {
            clustsBelong.add(j);
            curMax = tmp;
          }
        }
        clusterAssigns[i] =
            clustsBelong.isEmpty() ? 0 : (1 + clustsBelong.get(clustsBelong.size() - 1));
      }
    }

    return clusterAssigns;
  }
}
