package org.genvisis.filesys;

import java.io.Serializable;
import org.genvisis.common.SerializedFiles;

public class SerialIntMatrix implements Serializable {

  public static final long serialVersionUID = 1L;
  private final int[][] matrix;

  public int[][] getMatrix() {
    return matrix;
  }

  public SerialIntMatrix(int[][] matrix) {
    this.matrix = matrix;
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename);
  }

  public static SerialIntMatrix load(String filename) {
    return (SerialIntMatrix) SerializedFiles.readSerial(filename, true);
  }
}
