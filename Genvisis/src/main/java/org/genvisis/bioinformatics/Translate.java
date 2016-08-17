package org.genvisis.bioinformatics;

import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;

public class Translate {

  public static class Translation {
    public String str;
    public int start;

    public Translation() {}
  }

  public static final String[][] CODE =
      {{"TTT", "Phe", "F"}, {"TTC", "Phe", "F"}, {"TTA", "Leu", "L"}, {"TTG", "Leu", "L"},
          {"CTT", "Leu", "L"}, {"CTC", "Leu", "L"}, {"CTA", "Leu", "L"}, {"CTG", "Leu", "L"},
          {"ATT", "Ile", "I"}, {"ATC", "Ile", "I"}, {"ATA", "Ile", "I"}, {"ATG", "Met", "M"},
          {"GTT", "Val", "V"}, {"GTC", "Val", "V"}, {"GTA", "Val", "V"}, {"GTG", "Val", "V"},
          {"TCT", "Ser", "S"}, {"TCC", "Ser", "S"}, {"TCA", "Ser", "S"}, {"TCG", "Ser", "S"},
          {"CCT", "Pro", "P"}, {"CCC", "Pro", "P"}, {"CCA", "Pro", "P"}, {"CCG", "Pro", "P"},
          {"ACT", "Thr", "T"}, {"ACC", "Thr", "T"}, {"ACA", "Thr", "T"}, {"ACG", "Thr", "T"},
          {"GCT", "Ala", "A"}, {"GCC", "Ala", "A"}, {"GCA", "Ala", "A"}, {"GCG", "Ala", "A"},
          {"TAT", "Tyr", "Y"}, {"TAC", "Tyr", "Y"}, {"TAA", "STOP", "X"}, {"TAG", "STOP", "X"},
          {"CAT", "His", "H"}, {"CAC", "His", "H"}, {"CAA", "Gln", "Q"}, {"CAG", "Gln", "Q"},
          {"AAT", "Asn", "N"}, {"AAC", "Asn", "N"}, {"AAA", "Lys", "K"}, {"AAG", "Lys", "K"},
          {"GAT", "Asp", "D"}, {"GAC", "Asp", "D"}, {"GAA", "Glu", "E"}, {"GAG", "Glu", "E"},
          {"TGT", "Cys", "C"}, {"TGC", "Cys", "C"}, {"TGA", "STOP", "X"}, {"TGG", "Trp", "W"},
          {"CGT", "Arg", "R"}, {"CGC", "Arg", "R"}, {"CGA", "Arg", "R"}, {"CGG", "Arg", "R"},
          {"AGT", "Ser", "S"}, {"AGC", "Ser", "S"}, {"AGA", "Arg", "R"}, {"AGG", "Arg", "R"},
          {"GGT", "Gly", "G"}, {"GGC", "Gly", "G"}, {"GGA", "Gly", "G"}, {"GGG", "Gly", "G"}};


  public static String[] findSearchTerm(String target, String str) {
    Vector<String> v;
    int index;

    v = new Vector<String>();
    do {
      index = str.indexOf(target);
      if (index > -1) {
        str =
            str.substring(0, index) + target.toUpperCase() + str.substring(index + target.length());
        v.add(index + "");
      }
    } while (index > -1);

    System.out.println("Found search sequence " + v.size() + " time(s)");
    return new String[] {str, ext.listWithCommas(Array.toStringArray(v))};
  }

  public static void main(String[] args) {
    String str;
    String output;
    Translation[] translations;
    String[] searchPositions;
    String search;
    Translation longestTranslation;
    String forwardSeq, reverseSeq;
    boolean reverse;

    longestTranslation = new Translation();
    longestTranslation.str = "";
    reverse = false;
    str = ext.getClipboard();
    output = str + "\n\n\n\n";
    if (str.startsWith("find:")) {
      search = str.substring(0, str.indexOf("\n")).trim().split("[\\s]+")[1];
      str = str.substring(str.indexOf("\n") + 1);
    } else {
      search = null;
    }

    forwardSeq = procSequence(str);
    output += "Forward:\n" + forwardSeq + "\n\n";

    translations = translateSequence(forwardSeq);
    for (int i = 0; i < translations.length; i++) {
      output += "Translation #" + (i + 1) + " (starting at " + translations[i].start + "):\n"
          + translations[i].str + "\n\n";
      if (translations[i].str.length() > longestTranslation.str.length()) {
        longestTranslation = translations[i];
      }
    }

    if (search != null) {
      searchPositions = findSearchTerm(search, forwardSeq);
      output += "Forward marked:\n" + searchPositions[0] + "\n\n";
      output += "Forward positions:\t" + searchPositions[1] + "\n\n";
    }

    output += "\n\n\n";

    reverseSeq = reverseSequence(forwardSeq);
    output += "Reverse:\n" + reverseSeq + "\n\n";

    translations = translateSequence(reverseSeq);
    for (int i = 0; i < translations.length; i++) {
      output += "Translation #" + (i + 1) + " (starting at " + translations[i].start + "):\n"
          + translations[i].str + "\n\n";
      if (translations[i].str.length() > longestTranslation.str.length()) {
        longestTranslation = translations[i];
        reverse = true;
      }
    }
    output += "\n\n\n";

    if (search != null) {
      searchPositions = findSearchTerm(search, reverseSeq);
      output += "Reverse marked:\n" + searchPositions[0] + "\n\n";
      output += "Reverse positions:\t" + searchPositions[1] + "\n\n";
    }

    output += "\n\n\n";
    str = reverse ? reverseSeq : forwardSeq;

    for (int i = 0; i < longestTranslation.str.length(); i++) {
      output += i * 3 + longestTranslation.start + "\t" + str.substring(
          longestTranslation.start - 1 + i * 3, longestTranslation.start - 1 + i * 3 + 3);
      output += "\t" + (i + 1) + "\t" + longestTranslation.str.charAt(i) + "\n";
    }

    ext.setClipboard(output);
  }

  public static String procSequence(String str) {
    return ext.replaceAllWith(str,
        new String[][] {{"0", ""}, {"1", ""}, {"2", ""}, {"3", ""}, {"4", ""}, {"5", ""}, {"6", ""},
            {"7", ""}, {"8", ""}, {"9", ""}, {" ", ""}, {"\t", ""}, {"\n", ""}, {"\r", ""}})
        .toLowerCase();
  }

  public static String reverseSequence(String str) {
    return ext.replaceAllWith(str, new String[][] {{"a", "1"}, {"c", "2"}, {"g", "3"}, {"t", "4"},
        {"1", "t"}, {"2", "g"}, {"3", "c"}, {"4", "a"}});
  }

  public static Translation[] translateSequence(String str) {
    Translation[] translations;
    String[] codons;
    String codon;
    int frame, index;

    translations = new Translation[3];
    for (int i = 0; i < translations.length; i++) {
      translations[i] = new Translation();
    }
    codons = Matrix.extractColumn(CODE, 0);

    for (int i = 0; i < str.length() - 2; i++) {
      frame = i % 3;
      codon = str.substring(i, i + 3);
      index = ext.indexOfStr(codon.toUpperCase(), codons);
      if (translations[frame].str == null) {
        if (CODE[index][1].equals("Met")) {
          translations[frame].str = CODE[index][2];
          translations[frame].start = i + 1;
        }
      } else if (translations[frame].str.charAt(translations[frame].str.length() - 1) != 'X') {
        translations[frame].str += CODE[index][2];
      }
    }

    for (int i = 0; i < translations.length; i++) {
      if (translations[i].str != null) {
        translations[i].str = translations[i].str.substring(0, translations[i].str.length() - 1);
      }
    }

    return translations;
  }

}
