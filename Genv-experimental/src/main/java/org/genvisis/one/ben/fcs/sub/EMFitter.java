package org.genvisis.one.ben.fcs.sub;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math3.distribution.MixtureMultivariateNormalDistribution;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.distribution.fitting.MultivariateNormalMixtureExpectationMaximization;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.Pair;
import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.one.ben.fcs.FCSDataLoader;
import org.genvisis.one.ben.fcs.FCSDataLoader.DATA_SET;
import org.genvisis.one.ben.fcs.FCSDataLoader.LOAD_STATE;

import edu.stanford.facs.logicle.Logicle;

public class EMFitter {
  
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

  private static Logicle getLogicle(double[] data) {
    return getBiexScale(Array.max(data));
  }

  private static double[][] transformData(double[][] data, Logicle[] mehs, Logger log) {
    log.report("Transforming data to Logicle scale");

    double[][] transformed = new double[data.length][data[0].length];
    for (int i = 0; i < mehs.length; i++) {
      for (int j = 0; j < data[0].length; j++) {
        transformed[i][j] = mehs[i].scale(data[i][j]);
      }
    }
    log.report("Finished transforming data to Logicle scale");
    return transformed;
  }

  private static double[][] process(FCSDataLoader dataLoader, Logger log) {
    ArrayList<String> params = dataLoader.getAllDisplayableNames(DATA_SET.COMPENSATED);
    double[][] data = new double[params.size()][];
    for (int i = 0; i < params.size(); i++) {
      data[i] = dataLoader.getData(params.get(i), true);
    }
    
    Logicle[] mehs = new Logicle[params.size()];
    for (int i = 0; i < mehs.length; i++) {
      mehs[i] = getLogicle(data[i]);
    }
    return transformData(data, mehs, log);
  }
  
  private static double[][] transpose(double[][] matrixData) {
    RealMatrix m = MatrixUtils.createRealMatrix(matrixData);
    return m.transpose().getData();
  }

  private static double[][] selectMaxVarInitialMix(double[][] transposedfullData, double percentage, Logger log) {
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
        m.setRow(countMax + (numEach * ((2 * i) + 1)), transposedfullData[transposedfullData.length - countMax - 1]);
        countMax++;
      }
      
    }
    
    log.report("Finished selecting " + num + " initial points");

    return m.getData();
  }
  

  private static void sort(double[][] transposedfullData, final int sortInd) {
    Arrays.sort(transposedfullData, new Comparator<double[]>() {
      @Override
      public int compare(double[] o1, double[] o2) {
        return Double.valueOf(o1[sortInd]).compareTo(o2[sortInd]);
      }
    });
  }

  private static MultivariateNormalMixtureExpectationMaximization getClusters(double[][] transposed, double[][] init, int clusters, Logger log) {
    log.report("Detecting clusters, num clusters set to " + clusters);
    MixtureMultivariateNormalDistribution initialMix = MultivariateNormalMixtureExpectationMaximization.estimate(init, clusters);

    MultivariateNormalMixtureExpectationMaximization mle = new MultivariateNormalMixtureExpectationMaximization(transposed);

    mle.fit(initialMix, 500, 1E-15);

    return mle;
  }
  
  
  public int[] run(FCSDataLoader loader) {
    if (loader.getLoadState() != LOAD_STATE.LOADED) {
      System.err.println("Error - data not loaded yet; please wait a while and try again.");
      return null;
    }
    Logger log = new Logger();
    long t1 = System.currentTimeMillis();
    double[][] transformed = process(loader, log);
    double[][] transposed = transpose(transformed);
    
    double[][] initData = selectMaxVarInitialMix(transposed, .01, log);
    
    MultivariateNormalMixtureExpectationMaximization mnmem = null;
    for (int clust = 2; clust < 50; clust++) {
      try {
        MultivariateNormalMixtureExpectationMaximization mnmem1 = getClusters(transposed, initData, clust, log);
        mnmem = mnmem1;
        break;
      } catch (Exception e) {
        System.out.println("Clustering " + clust + " -- " + e.getMessage());
      }
    }
    MixtureMultivariateNormalDistribution mmnd = mnmem.getFittedModel();
    
    List<Pair<Double, MultivariateNormalDistribution>> clusts = mmnd.getComponents();

    double[][] means = new double[clusts.size()][];
    double[][] sds = new double[clusts.size()][];
    
    int[] clustAssigns = new int[transformed[0].length];
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
      double curMax = 0.0001;
      ArrayList<Integer> clustsBelong = new ArrayList<>();
      for (int j = 0; j < means.length; j++) {
        // determine if coords belongs to this cluster
        double tmp = clusts.get(j).getValue().density(line);
        if (tmp > curMax) {
          clustsBelong.add(j);
          curMax = tmp;
        }
      }
      
      clustAssigns[i] = clustsBelong.isEmpty() ? 0 : clustsBelong.get(clustsBelong.size() - 1);
      log.reportTimeElapsed("Clustered in ", t1);
    }
    
    return clustAssigns;
  }
  
}
