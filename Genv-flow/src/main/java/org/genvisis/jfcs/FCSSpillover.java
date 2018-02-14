package org.genvisis.jfcs;

public class FCSSpillover {

  private String[] parameterNames;
  private double[][] spilloverCoefficients;

  private FCSSpillover() {

  }

  public static FCSSpillover parse(String spilloverKeywordValue) {

    String[] sTmp = spilloverKeywordValue.split(",");
    if (sTmp.length < 2) {
      return null;
    }

    int n = Integer.parseInt(sTmp[0]);

    if (sTmp.length != 1 + n + n * n) return null;

    FCSSpillover spill = new FCSSpillover();

    spill.setParameterNames(new String[n]);
    spill.setCoefficients(new double[n][n]);

    for (int i = 0; i < n; i++) {
      spill.getParameterNames()[i] = sTmp[i + 1];
    }

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        spill.getCoefficients()[i][j] = Double.parseDouble(sTmp[1 + n + n * i + j]);
      }
    }

    return spill;
  }

  public String[] getParameterNames() {
    return parameterNames;
  }

  void setParameterNames(String[] parameterNames) {
    this.parameterNames = parameterNames;
  }

  public double[][] getCoefficients() {
    return spilloverCoefficients;
  }

  void setCoefficients(double[][] spilloverCoefficients) {
    this.spilloverCoefficients = spilloverCoefficients;
  }

}
