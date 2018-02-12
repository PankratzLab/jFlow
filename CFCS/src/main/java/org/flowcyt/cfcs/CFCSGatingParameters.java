package org.flowcyt.cfcs;
// CFCSGatingParameters.java

/*
 * ------------------------------------------------------------------------- *\ This software and
 * documentation are provided 'as is' and Tree Star, Inc., its contractors and partners specifically
 * disclaim all other warranties, expressed or implied, including but not limited to implied
 * warranties of merchantability and fitness for a particular purpose, or during any particular date
 * range. By using this software, you are agreeing to these limits of liability, and to hold Tree
 * Star harmless for any information, accurate or erroneous, that might be generated by the program.
 * This software is intended for research use only. Christopher Lane <cdl@best.classes> for Tree
 * Star 1/18/2002 Copyright 2002 \*
 * -------------------------------------------------------------------------
 */

import java.util.HashMap;
import java.util.Map;

public final class CFCSGatingParameters extends CFCSAbstractParameters {

  private final Map<Object, CFCSAbstractParameter> cache = new HashMap<Object, CFCSAbstractParameter>();

  // --------------------------------------------------------------------

  private final static String GATE_PREFIX = "$G";

  private final static String[][] GATE_PROPERTIES = {{"E", "LogDecadesAndOffset", "N"},
                                                     {"F", "Filter", "N"}, {"N", "ShortName", "N"},
                                                     {"P", "EmittedPercent", "N"},
                                                     {"R", "Range", "N"}, {"S", "FullName", "N"},
                                                     {"T", "DetectorType", "N"},
                                                     {"V", "Voltage", "N"},};

  // --------------------------------------------------------------------

  /* friendly */
  CFCSGatingParameters(final CFCSKeywords keywords) {
    super(keywords, CFCSGatingParameter.class);
  }

  // --------------------------------------------------------------------

  protected final String getPrefix() {
    return GATE_PREFIX;
  }

  // --------------------------------------------------------------------

  protected final String[][] getProperties() {
    return GATE_PROPERTIES;
  }

  // --------------------------------------------------------------------

  protected final String getCountKeyword() {
    return CFCSKeywords.GATE_KEYWORD;
  }

  // --------------------------------------------------------------------
  // Does this keyword name look like a gating parameter?  This is similar
  // to the same routine in CFCSParameters but since we can't overide
  // static methods, this can't be done in CFCSAbstractParameters

  /* friendly */
  static boolean isParameter(String keyword) {
    keyword = keyword.toUpperCase();

    if (!keyword.startsWith(GATE_PREFIX)) return false;

    String suffix = null;

    for (int i = 0; i < GATE_PROPERTIES.length; i++) {
      String code = GATE_PROPERTIES[i][PARAMETER_CODE];

      if (keyword.endsWith(code)) {
        suffix = code;
        break;
      }
    }

    if (suffix == null) return false;

    String root = keyword.substring(GATE_PREFIX.length(), keyword.length() - suffix.length());

    try {
      Integer.parseInt(root);
    } catch (NumberFormatException exception) {
      return false;
    }

    return true;
  }

  // --------------------------------------------------------------------

  public final CFCSGatingParameter getParameter(final int index) {
    final Object key = new Integer(index);

    if (cache.containsKey(key)) {
      return (CFCSGatingParameter) ((CFCSGatingParameter) cache.get(key)).copy();
    }

    CFCSGatingParameter parameter = new CFCSGatingParameter();

    super.getParameter(index, parameter);

    cache.put(key, parameter.copy());

    return parameter;
  }

  // --------------------------------------------------------------------

  public final void replaceParameter(final int index, final CFCSGatingParameter parameter) {
    super.replaceParameter(index, parameter);

    cache.remove(new Integer(index));
  }

  // --------------------------------------------------------------------

  public final void deleteParameter(final int index) {
    super.deleteParameter(index);

    cache.remove(new Integer(index));
  }

  // --------------------------------------------------------------------

}
