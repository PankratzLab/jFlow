package org.genvisis.one.ben;

import java.util.HashMap;

public class IDMasker {

  private static final String ALPH = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

  private static String coded(int val) {
    return String.format("%s%02d", ALPH.charAt(val / 100), val % 100);
  }

  private static String mask(int indi) {
    int bnd1 = 2600;
    int bnd2 = 6760000;
    if (indi < bnd1) {
      return "A00_A00_" + coded(indi);
    } else if (indi < bnd2) {
      int v = indi / bnd1;
      int c = indi - (bnd1 * v);
      return "A00_" + coded(v) + "_" + coded(c);
    } else {
      int v = indi / bnd2;
      int v2 = (indi - (v * bnd2)) / bnd1;
      int v3 = indi - (v * bnd2) - (v2 * bnd1);
      return coded(v) + "_" + coded(v2) + "_" + coded(v3);
    }
  }

  public static String[] maskIDs(String[] ids) {
    String[] result = new String[ids.length];
    for (int i = 0; i < ids.length; i++) {
      result[i] = mask(i);
    }
    return result;
  }

  public static HashMap<String, String> maskIDsMap(String[] ids) {
    HashMap<String, String> result = new HashMap<>();
    for (int i = 0; i < ids.length; i++) {
      result.put(ids[i], mask(i));
    }
    return result;
  }

}
