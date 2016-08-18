package org.genvisis.mining;

// public int n

import java.util.Random;

public class UnivariateGaussianMixtureModel {
  public int size;
  public double[] weight;
  public int[] pos;
  public PVector[] param;

  /**
   * Class constructor.
   */
  public UnivariateGaussianMixtureModel() {
    size = 0;
  }

  public void setSize(int n) {
    size = n;
    weight = new double[n];
    param = new PVector[n];
    pos = new int[n];
  }

  public UnivariateGaussianMixtureModel(int n) {
    size = n;
    weight = new double[n];
    param = new PVector[n];
    pos = new int[n];
  }

  /**
   * Computes the density value of this mixture value.
   *
   * @param x a point
   * @return density value of this mixture value
   */
  public double density(Point x) {
    double cumul = 0.0d;
    for (int i = 0; i < size; i++) {
      cumul += weight[i] * densityOfGaussian(x, param[i]);
    }
    return cumul;
  }

  /**
   * Normalizes the weights of the mixture models \f$ \alpha_i \f$ such as \f$ \sum_{i=1}^n \alpha_i
   * = 1 \f$.
   */
  public void normalizeWeights() {
    double sum = 0;
    int i;
    for (i = 0; i < size; i++) {
      sum += weight[i];
    }
    for (i = 0; i < size; i++) {
      weight[i] /= sum;
    }
  }

  /**
   * Method toString.
   *
   * @return string describing the mixture model
   */
  @Override
  public String toString() {
    String output = String.format("Mixture containing %d components\n", size);
    for (int i = 0; i < size; i++) {
      output += String.format("  Component %4d: ", pos[i]);
      output += String.format("Weight = %8.6f ", weight[i]);
      output += String.format("Parameters = %s\n", param[i]);
    }
    return output;
  }

  /**
   * Creates and returns a copy of the instance.
   *
   * @return a clone of the instance
   */
  @Override
  public UnivariateGaussianMixtureModel clone() {
    UnivariateGaussianMixtureModel mm = new UnivariateGaussianMixtureModel(size);
    mm.weight = weight.clone();
    mm.pos = pos.clone();
    for (int i = 0; i < size; i++) {
      mm.param[i] = param[i].clone();
    }
    return mm;
  }

  /**
   * Return the dimension of the parameters of the mixture model.
   *
   * @return parameters's dimension
   */
  public int getDimension() {
    return param[0].getDimension();
  }

  /**
   * Draws points from the considered mixture model.
   *
   * @param m number of points to draw
   * @return a point
   */
  public Point[] drawRandomPoints(int m) {

    // Array of points
    Point[] points = new Point[m];

    // Cumulative array
    int n = size;
    double[] t = new double[n];
    double sum = 0;
    for (int i = 0; i < n; i++) {
      sum += weight[i];
      t[i] = sum;
    }

    // Loop
    for (int i = 0; i < m; i++) {

      // Random number between 0 and 1
      double r = Math.random();

      // Find generative class
      int idx = 0;
      while (t[idx] < r && idx < n - 1) {
        idx++;
      }

      // Draw and return the point from the idx-th model
      points[i] = drawRandomPointFromGaussian(param[idx]);
    }
    return points;
  }

  /**
   * Computes the density value \f$ f(x;\mu,\sigma^2) \f$.
   *
   * @param x point
   * @param param parameters (source, natural, or expectation)
   * @return \f$ f(x;\mu,\sigma^2) = \frac{1}{ \sqrt{2\pi \sigma^2} } \exp \left( -
   *         \frac{(x-\mu)^2}{ 2 \sigma^2} \right) \f$
   */
  public static double densityOfGaussian(Point x, PVector param) {
    return Math.exp(-(x.value - param.array[0]) * (x.value - param.array[0])
                    / (2.0d * param.array[1]))
           / (Math.sqrt(2.0d * Math.PI * param.array[1]));
  }

  /**
   * Draws a point from the considered distribution.
   *
   * @param L source parameters \f$ \mathbf{\Lambda} = ( \mu , \sigma^2 )\f$
   * @return a point
   */
  public static Point drawRandomPointFromGaussian(PVector L) {
    double mean = L.array[0];
    double variance = L.array[1];

    // Draw the point
    Random rand = new Random();
    return new Point(mean + rand.nextGaussian() * Math.sqrt(variance));
  }
}
