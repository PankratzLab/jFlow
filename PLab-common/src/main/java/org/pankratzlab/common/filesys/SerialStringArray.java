package org.pankratzlab.common.filesys;

import java.io.Serializable;

import org.pankratzlab.common.Files;
import org.pankratzlab.common.SerializedFiles;

public class SerialStringArray implements Serializable {

  public static final long serialVersionUID = 1L;
  private final String[] array;

  public SerialStringArray(String[] array) {
    this.array = array;
  }

  public String[] getArray() {
    return array;
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename);
  }

  public static SerialStringArray load(String filename) {
    return (SerialStringArray) SerializedFiles.readSerial(filename, true);
  }

  public static void dump(String filename) {
    Files.writeArray(load(filename).getArray(), filename + ".xln");
  }
}
