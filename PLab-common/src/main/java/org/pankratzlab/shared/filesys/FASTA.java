package org.pankratzlab.shared.filesys;

public final class FASTA {

  private static final String INDEX_EXT = ".fai";
  private static final String DICTIONARY_EXT = ".dict";

  private FASTA() {
    // prevent instantiation of utility class
  }

  public static final String getIndex(String fasta) {
    return fasta + INDEX_EXT;
  }

  public static final String getDictionary(String fasta) {
    return fasta + DICTIONARY_EXT;
  }
}
