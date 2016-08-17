package org.genvisis.mining;

public class PredictionEvaluation {
  public static final int MAE = 0; // Mean Absolute Error (MAE)

  public static final int AE = 1; // Average Error (AE)

  public static final int MAPE = 2; // Mean Absolute Percentage Error (MAPE)

  public static final int SSE = 3; // Sum of Squared Errors (SSE)

  public static final int RMSE = 4; // Root Mean Squared Error (RMSE)

  public static double AE(double[] targets, double[] predicteds) {
    double err = 0;

    checkTargetsAndPredicteds(targets, predicteds);
    for (int i = 0; i < targets.length; i++) {
      err += targets[i] - predicteds[i];
    }

    return err / targets.length;
  }

  public static double calculateError(double[] targets, double[] predicteds) {
    return calculateError(targets, predicteds, RMSE);
  }

  public static double calculateError(double[] targets, double[] predicteds, int method) {
    switch (method) {
      case MAE:
        return MAE(targets, predicteds);
      case AE:
        return AE(targets, predicteds);
      case MAPE:
        return MAPE(targets, predicteds);
      case SSE:
        return SSE(targets, predicteds);
      case RMSE:
        return RMSE(targets, predicteds);
      default:
        System.err
            .println("Error - '" + method + "' does not map to an implemented method; using RMSE");
        return RMSE(targets, predicteds);
    }
  }

  public static void checkTargetsAndPredicteds(double[] targets, double[] predicteds) {
    if (targets.length != predicteds.length) {
      System.err.println("Error - targets array and predicteds array were of unequal lengths ("
          + targets.length + " and " + predicteds.length + ", respectively)");
      System.exit(1);
    }
    if (targets.length == 0) {
      System.err.println("Warning - targets array and predicteds array both have lengths of zero");
    }
  }

  public static double MAE(double[] targets, double[] predicteds) {
    double err = 0;

    checkTargetsAndPredicteds(targets, predicteds);
    for (int i = 0; i < targets.length; i++) {
      err += Math.abs(targets[i] - predicteds[i]);
    }

    return err / targets.length;
  }

  public static double MAPE(double[] targets, double[] predicteds) {
    double err = 0;

    checkTargetsAndPredicteds(targets, predicteds);
    for (int i = 0; i < targets.length; i++) {
      err += Math.abs((targets[i] - predicteds[i]) / targets[i]);
    }

    return err / targets.length * 100;
  }

  public static double RMSE(double[] targets, double[] predicteds) {
    return Math.sqrt(SSE(targets, predicteds) / targets.length);
  }

  public static double SSE(double[] targets, double[] predicteds) {
    double err = 0;

    checkTargetsAndPredicteds(targets, predicteds);
    for (int i = 0; i < targets.length; i++) {
      err += Math.pow(targets[i] - predicteds[i], 2);
    }

    return err;
  }
}
