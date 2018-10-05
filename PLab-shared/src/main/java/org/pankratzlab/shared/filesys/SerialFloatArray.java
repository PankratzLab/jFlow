package org.pankratzlab.shared.filesys;

import java.io.PrintWriter;
import java.io.Serializable;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.SerializedFiles;

public class SerialFloatArray implements Serializable {

  public static final long serialVersionUID = 1L;
  private final float[] array;

  public SerialFloatArray(float[] array) {
    this.array = array;
  }

  public float[] getArray() {
    return array;
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename);
  }

  public static SerialFloatArray load(String filename) {
    return (SerialFloatArray) SerializedFiles.readSerial(filename, true);
  }

  public static void dump(String filename) {
    float[] all = load(filename).getArray();

    try {
      PrintWriter writer = Files.openAppropriateWriter(filename + ".xln");
      for (float element : all) {
        writer.println(element);
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + filename + ".xln");
      e.printStackTrace();
    }
  }
}
