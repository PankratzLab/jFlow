package org.genvisis.filesys;

import java.io.Serializable;

import org.genvisis.common.Files;
import org.genvisis.common.SerializedFiles;

public class SerialStringArray implements Serializable {
  public static final long serialVersionUID = 1L;

  public static void dump(String filename) {
    Files.writeList(load(filename, false).getArray(), filename + ".xln");
  }

  public static SerialStringArray load(String filename, boolean jar) {
    return (SerialStringArray) SerializedFiles.readSerial(filename, jar, true);
  }

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
}
