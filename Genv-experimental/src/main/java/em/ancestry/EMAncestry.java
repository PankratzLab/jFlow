package em.ancestry;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.distribution.MixtureMultivariateNormalDistribution;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.distribution.fitting.MultivariateNormalMixtureExpectationMaximization;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.util.Pair;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.common.matrix.MatrixOperations;
import org.genvisis.common.matrix.MatrixOperations.SCALE_METHOD;

/**
 * Prototype to perform unsupervised clustering of genotype PCA results in n-dimensions using
 * {@link MultivariateNormalMixtureExpectationMaximization} to determine ancestry clusters <br>
 * Currently reports probability of cluster membership (relative to other clusters) and distance to
 * cluster assigned
 */
public class EMAncestry {

  /**
   * @param file file to load
   * @param indices data indices to load
   * @param log
   * @return
   */
  private static double[][] loadData(File file, int[] indices, Logger log) {
    log.reportTimeInfo("Loading data from " + file.getAbsolutePath());
    try {
      List<String> lines = FileUtils.readLines(file, "UTF-8");
      log.reportTimeInfo("Loading "
                         + ArrayUtils.toStr(ArrayUtils.subArray(lines.get(0).trim().split("\t"),
                                                                indices)));
      double[][] data = new double[indices.length][lines.size() - 1];
      int index = 0;
      boolean skip = true;
      for (String line : lines) {

        if (!skip) {
          double[] doubleValues = (double[]) ArrayUtils.toDoubleArray(ArrayUtils.subArray(line.trim()
                                                                                              .split("\t"),
                                                                                          indices));
          for (int i = 0; i < doubleValues.length; i++) {
            data[i][index] = doubleValues[i];
          }
          index++;
        }
        if (skip) {
          skip = false;
        }
      }

      log.reportTimeInfo("Finished loading " + lines.size() + " lines from "
                         + file.getAbsolutePath());

      return data;
    } catch (IOException e) {
      log.reportException(e);
    }
    throw new IllegalStateException("Could not load data");
  }

  /**
   * @param matrixData
   * @return transposed data
   */
  private static double[][] transpose(double[][] matrixData) {
    RealMatrix m = MatrixUtils.createRealMatrix(matrixData);
    return m.transpose().getData();
  }

  /**
   * @param data multidimensional data to clusters
   * @param clusters the number of clusters to detect
   * @param log
   * @return the parameters of the clusters
   */
  private static MultivariateNormalMixtureExpectationMaximization getClusters(double[][] transposed,
                                                                              double[][] init,
                                                                              int clusters,
                                                                              Logger log) {
    log.reportTimeInfo("Detecting clusters, num clusters set to " + clusters);
    MixtureMultivariateNormalDistribution initialMix = MultivariateNormalMixtureExpectationMaximization.estimate(init,
                                                                                                                 clusters);

    MultivariateNormalMixtureExpectationMaximization mle = new MultivariateNormalMixtureExpectationMaximization(transposed);

    mle.fit(initialMix, 10000, 1E-15);

    return mle;
  }

  private static class Cluster {

    private final int clusterAssigned;
    private final List<Double> allPosteriorProbabilities;
    private final double posteriorProbabilityOfAssignedCluster;
    private final double mahalanobisDistance;

    public Cluster(int clusterAssigned, double posteriorProbabilityOfAssignedCluster,
                   List<Double> allPosteriorProbabilities, double mahalanobisDistance) {
      super();
      this.clusterAssigned = clusterAssigned;
      this.posteriorProbabilityOfAssignedCluster = posteriorProbabilityOfAssignedCluster;
      this.allPosteriorProbabilities = allPosteriorProbabilities;
      this.mahalanobisDistance = mahalanobisDistance;
    }

    /*
     * (non-Javadoc)
     * @see java.lang.Object#toString()
     */
    @Override
    public String toString() {
      return clusterAssigned + "\t" + posteriorProbabilityOfAssignedCluster + "\t"
             + mahalanobisDistance + "\t" + ArrayUtils.toStr(allPosteriorProbabilities);
    }
  }

  public static RealMatrix getInvCovarianceMatrix(RealMatrix m) {
    if (m.getRowDimension() != m.getColumnDimension()) {
      throw new IllegalArgumentException("Number of columns does not match number of rows");
    }
    // adapted from https://github.com/apache/mahout/blob/7dff35bc3c61c3e0b95e8e59406742be15203a69/mr/src/main/java/org/apache/mahout/common/distance/MahalanobisDistanceMeasure.java
    // See http://www.mlahanas.de/Math/svd.htm for details,
    // which specifically details the case of covariance matrix inversion
    // Complexity: O(min(nm2,mn2))
    SingularValueDecomposition svd = new SingularValueDecomposition(m);
    RealMatrix sInv = svd.getS();
    // Inverse Diagonal Elems
    for (int i = 0; i < sInv.getRowDimension(); i++) {
      double diagElem = sInv.getEntry(i, i);
      if (diagElem > 0.0) {
        sInv.addToEntry(i, i, 1 / diagElem);
      } else {
        throw new IllegalStateException("Eigen Value equals to 0 found.");
      }
    }
    return svd.getU().multiply(sInv.multiply(svd.getU().transpose()));
  }

  //  https://en.wikipedia.org/wiki/Mahalanobis_distance
  private static double computeMahalanobisDistance(double[] xy, double[] means,
                                                   RealMatrix invCovariances) {

    RealVector v = new ArrayRealVector(xy);
    RealVector meanVector = new ArrayRealVector(means);

    return Math.sqrt(v.subtract(meanVector)
                      .dotProduct(invCovariances.preMultiply(v.subtract(meanVector))));

  }

  private static List<Cluster> getClusters(MixtureMultivariateNormalDistribution mmnd,
                                           double[][] data, Logger log) {

    List<Pair<Double, MultivariateNormalDistribution>> clusts = mmnd.getComponents();
    List<Cluster> clusters = new ArrayList<>();
    List<RealMatrix> invCovariances = new ArrayList<>();
    log.reportTimeInfo("Computing inverse covariance matrices for clusters");
    for (int j = 0; j < clusts.size(); j++) {
      invCovariances.add(getInvCovarianceMatrix(clusts.get(j).getValue().getCovariances()));
    }

    for (int i = 0; i < data[0].length; i++) {
      if (i % 10000 == 0) {
        log.reportTimeInfo("Assigning clusters for " + (i + 1) + " of " + data[0].length
                           + " data points");

      }
      double[] xy = new double[data.length];
      for (int j = 0; j < xy.length; j++) {
        xy[j] = data[j][i];
      }

      //PDF values are all >0
      double curMax = -1;

      List<Double> posteriorProbs = new ArrayList<>();
      int mostLikelyCluster = -1;
      double mahalanobisDistance = -1;
      // pdf of this point in full MixtureMultivariateNormalDistribution 
      double fullPdf = mmnd.density(xy);
      for (int j = 0; j < clusts.size(); j++) {
        // posterior probability given this point and this component (cluster) of the MixtureMultivariateNormalDistribution
        double tmp = clusts.get(j).getValue().density(xy) * clusts.get(j).getFirst().doubleValue();
        clusts.get(j).getValue().getCovariances();
        tmp = tmp / fullPdf;
        posteriorProbs.add(tmp);
        // assign to most likely cluster
        if (tmp > curMax) {
          curMax = tmp;
          mostLikelyCluster = j;
          mahalanobisDistance = computeMahalanobisDistance(xy, clusts.get(j).getValue().getMeans(),
                                                           invCovariances.get(j));
        }
      }
      clusters.add(new Cluster(mostLikelyCluster, curMax, posteriorProbs, mahalanobisDistance));
    }
    log.reportTimeInfo("Finished applying detected clusters");
    return clusters;

  }

  public static void main(String[] args) {
    String pcFile = args[0];
    int numClusters = 4;
    int[] indicesInFile = new int[] {1, 2, 3};

    Logger log = new Logger(ext.parseDirectoryOfFile(pcFile) + "cluster.log");
    log.reportTimeInfo("Loading " + pcFile);
    double[][] datas = loadData(new File(pcFile), indicesInFile, log);
    String[] headerData = ArrayUtils.subArray(Files.getHeaderOfFile(pcFile, log), indicesInFile);
    String[] samples = HashVec.loadFileToStringArray(pcFile, true, new int[] {0}, false);

    log.reportTimeInfo("transposing data " + pcFile);

    double[][] dataT = transpose(datas);
    MultivariateNormalMixtureExpectationMaximization maximization = getClusters(dataT, dataT,
                                                                                numClusters, log);
    MixtureMultivariateNormalDistribution mmnd = maximization.getFittedModel();
    List<Cluster> clusters = getClusters(mmnd, datas, log);
    PrintWriter writer = Files.getAppropriateWriter(ext.addToRoot(pcFile, ".clusters"));
    StringJoiner joiner = new StringJoiner("\t");
    joiner.add("SAMPLE");
    joiner.add("CLUSTER_ASSIGNED");
    joiner.add("PROBALITY_OF_CLUSTER_ASSIGNED");
    joiner.add("MAHALANOBIS_DISTANCE");
    //    mahalanobisDistance
    for (int i = 0; i < numClusters; i++) {
      joiner.add("PROBALITY_OF_CLUSTER_" + i);
    }
    for (String h : headerData) {
      joiner.add(h);
    }
    writer.println(joiner.toString());
    for (int i = 0; i < clusters.size(); i++) {
      StringJoiner joinerSample = new StringJoiner("\t");
      joinerSample.add(samples[i]);
      joinerSample.add(clusters.get(i).toString());
      for (int j = 0; j < datas.length; j++) {
        joinerSample.add(Double.toString(datas[j][i]));
      }
      writer.println(joinerSample.toString());
    }

    writer.close();

    PrintWriter writerDist = Files.getAppropriateWriter(ext.addToRoot(pcFile, ".model"));
    StringJoiner joinerModel = new StringJoiner("\t");
    joinerModel.add("CLUSTER");
    joinerModel.add("WEIGHT");
    for (String h : headerData) {
      joinerModel.add(h + "_MEAN");
      joinerModel.add(h + "_SD");
    }
    writerDist.println(joinerModel.toString());
    for (int i = 0; i < numClusters; i++) {
      StringJoiner joinerModelD = new StringJoiner("\t");

      joinerModelD.add(Integer.toString(i));
      joinerModelD.add(Double.toString(mmnd.getComponents().get(i).getFirst().doubleValue()));

      for (int j = 0; j < datas.length; j++) {
        joinerModelD.add(Double.toString(mmnd.getComponents().get(i).getValue().getMeans()[j]));
        joinerModelD.add(Double.toString(mmnd.getComponents().get(i).getValue()
                                             .getStandardDeviations()[j]));

      }
      writerDist.println(joinerModelD.toString());
    }
    writerDist.close();
  }

}
