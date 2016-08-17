package org.genvisis.cnv.filesys;

import java.io.Serializable;
import java.util.Hashtable;

import org.genvisis.common.SerializedFiles;

public class OutlierHashtable implements Serializable {
  private static final long serialVersionUID = 1L;
  public static final String[] DATA_ITEMS = new String[] {};
  public static final int MAX_NUM_OF_SAMPLES_BIG_HASH = 1000000;
  public static final int MAX_NUM_OF_SAMPLES_SMALL_HASH = 10000; // Need to make this dynamic

  public static OutlierHashtable load(String filename) {
    return (OutlierHashtable) SerializedFiles.readSerial(filename);
  }

  private final long sampFingerPrintFromProj;
  private final long markFingerPrintFromProj;
  Hashtable<Integer, Float> outlierHashtableSmall;

  Hashtable<Long, Float> outlierHashtableBig;

  public OutlierHashtable(long sampleFingerPrintOfTheProj, long markerFingerPrintOfTheProj) {
    sampFingerPrintFromProj = sampleFingerPrintOfTheProj;
    markFingerPrintFromProj = markerFingerPrintOfTheProj;
    outlierHashtableBig = new Hashtable<Long, Float>();
  }

  public void add(int markerIndexInProj, int sampleIndexInProj, byte dataItem, float value) {
    if (outlierHashtableSmall == null) {
      outlierHashtableBig.put((long) (markerIndexInProj * MAX_NUM_OF_SAMPLES_BIG_HASH
                                      + sampleIndexInProj * 10 + dataItem),
                              value);
    } else {
      outlierHashtableSmall.put(markerIndexInProj * MAX_NUM_OF_SAMPLES_SMALL_HASH
                                + sampleIndexInProj * 10 + dataItem, value);
    }
  }

  public long getMarkFingerPrint() {
    return markFingerPrintFromProj;
  }

  public long getSampFingerPrint() {
    return sampFingerPrintFromProj;
  }

  public float getValue(int markerIndexInProj, int sampleIndexInProj, byte dataItem) {
    if (outlierHashtableSmall == null) {
      return outlierHashtableBig.get((long) markerIndexInProj * MAX_NUM_OF_SAMPLES_BIG_HASH
                                     + sampleIndexInProj * 10 + dataItem);
    } else {
      return outlierHashtableSmall.get(markerIndexInProj * MAX_NUM_OF_SAMPLES_SMALL_HASH
                                       + sampleIndexInProj * 10 + dataItem);
    }
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename);
  }
}
