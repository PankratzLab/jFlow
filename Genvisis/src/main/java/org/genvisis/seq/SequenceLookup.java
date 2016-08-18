package org.genvisis.seq;

import java.util.Hashtable;

import org.genvisis.common.Logger;

/**
 * @author lane0212<br>
 *         Supplies a lookup for Fasta formatted nucleotides, useful when a position can be multiple
 *         nucleotides such as Affy Probe sequences<br>
 *
 *         We do not handle "U" - only handling ACTG combos
 */
public class SequenceLookup {
  private final Hashtable<String, String> nucleotideLookup;
  private final Logger log;

  public SequenceLookup(Logger log) {
    super();
    nucleotideLookup = getNucleotideLookup();
    this.log = log;
  }

  /**
   * @param combo such as "AG" meaning A or G, and "AT" meaning A or T
   * @return the IUPAC representation
   */
  public String getNucForCombo(String combo) {
    if (!nucleotideLookup.containsKey(combo)) {
      log.reportTimeError("Invalid combo " + combo);
      return null;
    } else {
      return nucleotideLookup.get(combo);
    }
  }

  /**
   * 
   * Derived from https://en.wikipedia.org/wiki/FASTA_format
   * 
   * @return
   */
  private static Hashtable<String, String> getNucleotideLookup() {
    Hashtable<String, String> lookup = new Hashtable<String, String>();
    addCombo(new String[] {"A"}, "A", lookup);
    addCombo(new String[] {"C"}, "C", lookup);
    addCombo(new String[] {"G"}, "G", lookup);
    addCombo(new String[] {"T"}, "T", lookup);

    addCombo(new String[] {"A", "G"}, "R", lookup);
    addCombo(new String[] {"C", "T"}, "Y", lookup);
    addCombo(new String[] {"G", "T"}, "K", lookup);
    addCombo(new String[] {"A", "C"}, "M", lookup);
    addCombo(new String[] {"C", "G"}, "S", lookup);
    addCombo(new String[] {"A", "T"}, "W", lookup);

    addCombo(new String[] {"C", "G", "T"}, "B", lookup);
    addCombo(new String[] {"A", "G", "T"}, "D", lookup);
    addCombo(new String[] {"A", "C", "T"}, "H", lookup);
    addCombo(new String[] {"A", "C", "G"}, "V", lookup);

    addCombo(new String[] {"A", "T", "G", "C"}, "N", lookup);

    return lookup;

  }

  private static void addCombo(String[] combo, String value, Hashtable<String, String> lookup) {
    for (int i = 0; i < combo.length; i++) {
      String curLook = combo[i];
      for (int j = 0; j < combo.length; j++) {
        if (i != j) {
          curLook += combo[j];
        }
      }
      System.out.println("putting " + curLook + " to " + value);
      lookup.put(curLook, value);
    }
  }

  public static void main(String[] args) {
    getNucleotideLookup();
  }

}
