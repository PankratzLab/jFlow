package org.genvisis.filesys;

import java.io.Serializable;
import org.genvisis.common.SerializedFiles;

public class SerialByteMatrix implements Serializable {

  public static final long serialVersionUID = 1L;
  private final byte[][] matrix;

  public byte[][] getMatrix() {
    return matrix;
  }

  public SerialByteMatrix(byte[][] matrix) {
    this.matrix = matrix;
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename);
  }

  public static SerialByteMatrix load(String filename) {
    return (SerialByteMatrix) SerializedFiles.readSerial(filename, true);
  }
}
