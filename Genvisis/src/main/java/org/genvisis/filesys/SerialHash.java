package org.genvisis.filesys;

import java.io.Serializable;
import java.util.Hashtable;

import org.genvisis.common.SerializedFiles;

public class SerialHash implements Serializable {
  public static final long serialVersionUID = 1L;

  public static void createSerializedIntHash(String filename, Hashtable<Integer, Integer> intHash) {
    SerialHash sHash = new SerialHash();
    sHash.setIntHash(intHash);
    sHash.serialize(filename);
  }

  public static void createSerializedLongHash(String filename, Hashtable<Long, Long> longHash) {
    SerialHash sHash = new SerialHash();
    sHash.setLongHash(longHash);
    sHash.serialize(filename);
  }

  public static void createSerializedStringArrayHash(String filename,
      Hashtable<String, String[]> stringArrayHash) {
    SerialHash sHash = new SerialHash();
    sHash.setStringArrayHash(stringArrayHash);
    sHash.serialize(filename);
  }

  public static void createSerializedStringHash(String filename,
      Hashtable<String, String> stringHash) {
    SerialHash sHash = new SerialHash();
    sHash.setStringHash(stringHash);
    sHash.serialize(filename);
  }

  public static void createSerializedStringIntArrayHash(String filename,
      Hashtable<String, int[]> stringIntArrayHash) {
    SerialHash sHash = new SerialHash();
    sHash.setStringIntArrayHash(stringIntArrayHash);
    sHash.serialize(filename);
  }

  private static SerialHash load(String filename, boolean jar) {
    return (SerialHash) SerializedFiles.readSerial(filename, jar, true);
  }

  public static Hashtable<Integer, Integer> loadSerializedIntHash(String filename) {
    return load(filename, false).getIntHash();
  }

  public static Hashtable<Long, Long> loadSerializedLongHash(String filename) {
    return load(filename, false).getLongHash();
  }

  public static Hashtable<String, String[]> loadSerializedStringArrayHash(String filename) {
    return load(filename, false).getStringArrayHash();
  }

  public static Hashtable<String, String> loadSerializedStringHash(String filename) {
    return load(filename, false).getStringHash();
  }

  public static Hashtable<String, int[]> loadSerializedStringIntArrayHash(String filename) {
    return load(filename, false).getStringIntArrayHash();
  }

  private Hashtable<String, String> stringHash;

  private Hashtable<String, String[]> stringArrayHash;

  private Hashtable<String, int[]> stringIntArrayHash;

  private Hashtable<Integer, Integer> intHash;

  private Hashtable<Long, Long> longHash;

  public SerialHash() {
    stringHash = null;
    stringArrayHash = null;
    stringIntArrayHash = null;
    intHash = null;
    longHash = null;
  }

  private Hashtable<Integer, Integer> getIntHash() {
    return intHash;
  }

  private Hashtable<Long, Long> getLongHash() {
    return longHash;
  }

  private Hashtable<String, String[]> getStringArrayHash() {
    return stringArrayHash;
  }

  private Hashtable<String, String> getStringHash() {
    return stringHash;
  }

  private Hashtable<String, int[]> getStringIntArrayHash() {
    return stringIntArrayHash;
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename);
  }

  private void setIntHash(Hashtable<Integer, Integer> intHash) {
    this.intHash = intHash;
  }

  private void setLongHash(Hashtable<Long, Long> longHash) {
    this.longHash = longHash;
  }

  private void setStringArrayHash(Hashtable<String, String[]> stringArrayHash) {
    this.stringArrayHash = stringArrayHash;
  }

  private void setStringHash(Hashtable<String, String> stringHash) {
    this.stringHash = stringHash;
  }

  private void setStringIntArrayHash(Hashtable<String, int[]> stringIntArrayHash) {
    this.stringIntArrayHash = stringIntArrayHash;
  }
}
