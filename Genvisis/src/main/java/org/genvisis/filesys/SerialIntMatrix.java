package org.genvisis.filesys;

import java.io.Serializable;

import org.genvisis.common.SerializedFiles;

public class SerialIntMatrix implements Serializable {
  public static final long serialVersionUID = 1L;

  public static SerialIntMatrix load(String filename, boolean jar) {
    return (SerialIntMatrix) SerializedFiles.readSerial(filename, jar, true);
  }

  private final int[][] matrix;

  public SerialIntMatrix(int[][] matrix) {
    this.matrix = matrix;
  }

  public int[][] getMatrix() {
    return matrix;
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename);
  }
}
