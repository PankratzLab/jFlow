package org.genvisis.seq.cnv;

import java.io.Serializable;

public abstract class CNVExtraInfo implements Serializable {
  public enum EXTRA_INFO_TYPE {
                               EXOME_DEPTH;
  }

  private static final long serialVersionUID = 1L;
  protected double dExtra;
  protected String sExtra;
  protected boolean boolExtra;
  protected EXTRA_INFO_TYPE type;

  public CNVExtraInfo(EXTRA_INFO_TYPE type) {
    this.type = type;
  }

  public double getdExtra() {
    return dExtra;
  }

  public String getsExtra() {
    return sExtra;
  }

  public boolean isBoolExtra() {
    return boolExtra;
  }

  public void setBoolExtra(boolean boolExtra) {
    this.boolExtra = boolExtra;
  }

  public void setdExtra(double dExtra) {
    this.dExtra = dExtra;
  }

  public void setsExtra(String sExtra) {
    this.sExtra = sExtra;
  }

}
