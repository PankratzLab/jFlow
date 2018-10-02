package org.pankratzlab.shared.filesys;

import java.io.Serializable;
import org.pankratzlab.common.SerializedFiles;

public class SerialDoubleArray implements Serializable {

  public static final long serialVersionUID = 1L;
  private final double[] array;

  public SerialDoubleArray(double[] array) {
    this.array = array;
  }

  public double[] getArray() {
    return array;
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename);
  }

  public static SerialDoubleArray load(String filename) {
    return (SerialDoubleArray) SerializedFiles.readSerial(filename, true);
  }
}
