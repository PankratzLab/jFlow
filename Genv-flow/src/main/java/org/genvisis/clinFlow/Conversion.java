package org.genvisis.clinFlow;

import java.io.File;

public class Conversion {

  public Conversion(File dir2, File out2, String f) {
    this.dir = dir2;
    this.out = out2;
    this.fcs = f;
  }

  final File dir;
  final File out;
  final String fcs;
}
