package org.genvisis.clinFlow;

import java.io.File;

public class Conversion {

  public Conversion(File dir2, File out2, String linkD, String f) {
    this.dir = dir2;
    this.out = out2;
    this.link = linkD;
    this.fcs = f;
  }

  public final File dir;
  public final File out;
  public final String link;
  public final String fcs;
}
