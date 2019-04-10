package org.flowcyt.cfcs;
// CFCSParameters.java

/*
 * ------------------------------------------------------------------------- *\ This software and
 * documentation are provided 'as is' and Tree Star, Inc., its contractors and partners specifically
 * disclaim all other warranties, expressed or implied, including but not limited to implied
 * warranties of merchantability and fitness for a particular purpose, or during any particular date
 * range. By using this software, you are agreeing to these limits of liability, and to hold Tree
 * Star harmless for any information, accurate or erroneous, that might be generated by the program.
 * This software is intended for research use only. Christopher Lane <cdl@best.classes> for Tree
 * Star 1/14/2002 Copyright 2002 \*
 * -------------------------------------------------------------------------
 */

import java.util.HashMap;
import java.util.Map;

public final class CFCSParameters extends CFCSAbstractParameters {

  private final Map<Object, CFCSAbstractParameter> cache = new HashMap<>();

  // --------------------------------------------------------------------

  private final static String DATA_PREFIX = "$P";

  private final static String[][] DATA_PROPERTIES = {{"S", "FullName", "N"},
                                                     {"N", "ShortName", "N"},
                                                     {"B", "FieldSizeString", "Y"},
                                                     {"R", "Range", "Y"}, {"G", "Gain", "N"},
                                                     {"F", "Filter", "N"},
                                                     {"L", "ExcitationWavelength", "N"},
                                                     {"O", "LaserPower", "N"},
                                                     {"P", "EmittedPercent", "N"},
                                                     {"V", "Voltage", "N"},
                                                     {"T", "DetectorType", "N"},
                                                     {"E", "LogDecadesAndOffset", "Y"},
                                                     {"D", "PreferredDisplay", "N"},
                                                     {"Calibration", "Calibration", "N"},};

  // --------------------------------------------------------------------

  /* friendly */
  CFCSParameters(final CFCSKeywords keywords) {
    super(keywords, CFCSParameter.class);
  }

  // --------------------------------------------------------------------

  protected final String getPrefix() {
    return DATA_PREFIX;
  }

  // --------------------------------------------------------------------

  protected final String[][] getProperties() {
    return DATA_PROPERTIES;
  }

  // --------------------------------------------------------------------

  protected final String getCountKeyword() {
    return CFCSKeywords.PARAMETER_KEYWORD;
  }

  // --------------------------------------------------------------------
  // Does this keyword name look like a data parameter? This is similar
  // to the same routine in CFCSGatingParameters but since we can't overide
  // static methods, this can't be done in CFCSAbstractParameters

  /* friendly */
  static boolean isParameter(String keyword) {
    keyword = keyword.toUpperCase();

    if (!keyword.startsWith(DATA_PREFIX)) return false;

    String suffix = null;

    for (int i = 0; i < DATA_PROPERTIES.length; i++) {
      String code = DATA_PROPERTIES[i][PARAMETER_CODE];

      if (keyword.endsWith(code)) {
        suffix = code;
        break;
      }
    }

    if (suffix == null) return false;

    String root = keyword.substring(DATA_PREFIX.length(), keyword.length() - suffix.length());

    try {
      Integer.parseInt(root);
    } catch (NumberFormatException exception) {
      return false;
    }

    return true;
  }

  // --------------------------------------------------------------------

  public final CFCSParameter getParameter(final int index) {
    final Object key = new Integer(index);

    if (cache.containsKey(key)) {
      return (CFCSParameter) ((CFCSParameter) cache.get(key)).copy();
    }

    CFCSParameter parameter = new CFCSParameter();

    super.getParameter(index, parameter);

    cache.put(key, parameter.copy());

    return parameter;
  }

  // --------------------------------------------------------------------
  public final void addParameter(final CFCSParameter parameter) {
    if (keywords.getDatatype() != CFCSDatatype.ASCII) {
      if (parameter.getRange() > Math.pow(2, parameter.getFieldSize())) {
        throw new CFCSError(CFCSInconsistentAttribute, "Range > pow(2, FieldSize)");
      }
    }
    super.addParameter(parameter);
  }

  // --------------------------------------------------------------------
  public final void replaceParameter(final int index, final CFCSParameter parameter) {
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
