package org.flowcyt.cfcs;
// CFCSKeyword.java

/*
 * ------------------------------------------------------------------------- *\ This software and
 * documentation are provided 'as is' and Tree Star, Inc., its contractors and partners specifically
 * disclaim all other warranties, expressed or implied, including but not limited to implied
 * warranties of merchantability and fitness for a particular purpose, or during any particular date
 * range.
 * 
 * By using this software, you are agreeing to these limits of liability, and to hold Tree Star
 * harmless for any information, accurate or erroneous, that might be generated by the program. This
 * software is intended for research use only.
 * 
 * Christopher Lane <cdl@best.classes> for Tree Star 1/14/2002 Copyright 2002 \*
 * -------------------------------------------------------------------------
 */


public final class CFCSKeyword implements Cloneable, CFCSErrorCodes {

  private static final char FCS_DEFINED_CHAR = '$';

  private String name = null;
  private String value = null;
  private int segIdx = CFCSDataSet.UNDEFINED;

  // --------------------------------------------------------------------

  /* friendly */
  CFCSKeyword(final String name, final String value, final int segIdx) {
    setKeywordName(name);
    setKeywordValue(value);
    setKeywordSource(segIdx);
  }

  // --------------------------------------------------------------------

  /* friendly */
  CFCSKeyword(final String name, final String value) {
    setKeywordName(name);
    setKeywordValue(value);
    setKeywordSource(CFCSDataSet.TEXT);
  }

  // --------------------------------------------------------------------

  /* friendly */
  CFCSKeyword(final String name, final int value) {
    setKeywordName(name);
    setKeywordIntegerValue(value);
    setKeywordSource(CFCSDataSet.TEXT);
  }

  // --------------------------------------------------------------------

  /* friendly */
  CFCSKeyword(final String name, final long value) {
    setKeywordName(name);
    setKeywordLongValue(value);
    setKeywordSource(CFCSDataSet.TEXT);
  }

  // --------------------------------------------------------------------

  // /* friendly */
  static final boolean isSet(final String value) {
    return value != null;
  }

  /* friendly */
  static final boolean isNotSet(final String value) {
    return value == null;
  }

  // --------------------------------------------------------------------

  /* friendly */
  final boolean isEmpty(final String value) {
    return (isNotSet(value) || value.length() == 0);
  }

  // --------------------------------------------------------------------

  /* friendly */
  final CFCSKeyword copy() {
    CFCSKeyword duplicate = null;

    try {
      duplicate = (CFCSKeyword) this.clone();
    } catch (CloneNotSupportedException exception) {
      throw new CFCSError(CFCSSystemError, exception);
    }

    return duplicate;
  }

  // --------------------------------------------------------------------

  /* friendly */
  final boolean nameStartsWith(final String prefix) {
    if (isNotSet(name)) {
      throw new CFCSError(CFCSKeywordUndefined);
    }

    return (name.toUpperCase()).startsWith(prefix.toUpperCase());
  }

  // --------------------------------------------------------------------

  /* friendly */
  final boolean isFCSDefined() {
    if (isNotSet(name)) {
      throw new CFCSError(CFCSKeywordUndefined);
    }

    return (name.charAt(0) == FCS_DEFINED_CHAR);
  }

  // --------------------------------------------------------------------

  public CFCSKeyword() {}

  // --------------------------------------------------------------------
  // Get/set this keyword's name as/to a string.

  public final String getKeywordName() {
    if (isNotSet(name)) {
      throw new CFCSError(CFCSKeywordUndefined);
    }

    return name;
  }

  public final void setKeywordName(final String name) {
    if (isEmpty(name)) { /* 3.2.9 */
      throw new CFCSError(CFCSIllegalName);
    }
    /* else */ this.name = name;
  }

  // --------------------------------------------------------------------
  // Get/set this keyword's value as/to a string.

  public final String getKeywordValue() {
    if (isNotSet(value)) {
      throw new CFCSError(CFCSKeywordUndefined);
    }

    return value;
  }

  public final void setKeywordValue(final String value) {
    if (isEmpty(value)) { /* 3.2.9 */
      throw new CFCSError(CFCSIllegalValue);
    }
    /* else */ this.value = value;
  }

  // --------------------------------------------------------------------
  // Get/set this keyword's value as/to a long.

  public final long getKeywordLongValue() {
    Long value = null;

    try {
      // JS, Feb 24, 2006
      String s = getKeywordValue().trim();
      int n = s.indexOf('.');
      if (n > 0)
        s = s.substring(0, n);
      value = new Long(s);
    } catch (NumberFormatException exception) {
      throw new CFCSError(CFCSBadValueConversion, exception);
    }

    return value.longValue();
  }

  public final void setKeywordLongValue(final long value) {
    this.value = (new Long(value)).toString();
  }

  // --------------------------------------------------------------------
  // Get/set this keyword's value as/to an integer.

  public final int getKeywordIntegerValue() {
    Integer value = null;

    try {
      // JS, Feb 24, 2006
      String s = getKeywordValue().trim();
      int n = s.indexOf('.');
      if (n > 0)
        s = s.substring(0, n);
      value = new Integer(s);
    } catch (NumberFormatException exception) {
      throw new CFCSError(CFCSBadValueConversion, exception);
    }

    return value.intValue();
  }

  public final void setKeywordIntegerValue(final int value) {
    this.value = (new Integer(value)).toString();
  }

  // --------------------------------------------------------------------
  // Get/set this keyword's value as/to a double.

  public final double getKeywordDoubleValue() {
    Double value = null;

    try {
      value = new Double(getKeywordValue().trim());
    } catch (NumberFormatException exception) {
      throw new CFCSError(CFCSBadValueConversion, exception);
    }

    return value.doubleValue();
  }

  public final void setKeywordDoubleValue(final double value) {
    this.value = (new Double(value)).toString();
  }

  // --------------------------------------------------------------------
  // Get/set the source segment of this keyword.

  public final int getKeywordSource() {
    if (segIdx == CFCSDataSet.UNDEFINED) {
      throw new CFCSError(CFCSKeywordUndefined);
    }

    return segIdx;
  }

  public final void setKeywordSource(final int segIdx) {
    if (segIdx < CFCSDataSet.OTHER_START && segIdx != CFCSDataSet.TEXT
        && segIdx != CFCSDataSet.ANALYSIS) {
      throw new CFCSError(CFCSIllegalSegment, segIdx);
    }
    /* else */ this.segIdx = segIdx;
  }

  // --------------------------------------------------------------------

  public final String toString() {
    return getKeywordValue();
  }

  // --------------------------------------------------------------------

}
