package org.genvisis.jfcs;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.nio.ByteOrder;
import java.nio.charset.Charset;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.codec.CharEncoding;

import com.google.common.base.CharMatcher;

public class FCSKeywords {

  byte delimByte = (byte) 12; // default
  Map<String, String> raw;
  private String[] parameterNames;
  private String[] parameterNamesLong;

  protected FCSKeywords(byte[] textBytes) throws UnsupportedEncodingException {
    this.delimByte = textBytes[0];
    parseRaw(textBytes);
    parseParamNames();
  }

  public int getParameterCount() {
    return getKeywordInt(FCS_KEYWORD.PAR.keyword);
  }

  public String getParameterLongName(int index) {
    return parameterNamesLong[index];
  }

  public String getParameterShortName(int index) {
    return parameterNames[index];
  }

  public double[] getParamScaling() {
    double[] scal = new double[getParameterCount()];
    for (int i = 0, count = getParameterCount(); i < count; i++) {
      scal[i] = 1;
      try {
        scal[i] = Double.parseDouble(getKeyword(FCS_KEYWORD.PnG.getKey(i + 1)));
      } catch (Exception e) {
      }
    }
    return scal;
  }

  public int getEventCount() {
    return getKeywordInt(FCS_KEYWORD.TOT.keyword)
        - (hasKeyword(FCS_KEYWORD.ABRT.keyword) ? getKeywordInt(FCS_KEYWORD.ABRT.keyword) : 0)
        - (hasKeyword(FCS_KEYWORD.LOST.keyword) ? getKeywordInt(FCS_KEYWORD.LOST.keyword) : 0);
  }

  public ByteOrder getEndianness() {
    String ord = raw.get(FCS_KEYWORD.BYTEORD.keyword).trim();
    if ("4,3,2,1".equals(ord)) {
      return ByteOrder.BIG_ENDIAN;
    } else if ("1,2,3,4".equals(ord)) {
      return ByteOrder.LITTLE_ENDIAN;
    }
    return null;
  }

  public Map<String, String> getRaw() {
    return raw;
  }

  public String setKeyword(FCS_KEYWORD key, String value) {
    return setKeyword(key.keyword, value);
  }

  public String setKeyword(String key, String value) {
    if (key.length() == 0 || value.length() == 0) {
      throw new IllegalArgumentException("Key/Value cannot be empty: {" + key + "," + value + "}");
    }
    if (!CharMatcher.ascii().matchesAllOf(key)) {
      throw new IllegalArgumentException("Key must be ascii only: " + key);
    }
    if (!Charset.forName(CharEncoding.UTF_8).newEncoder().canEncode(value)) {
      throw new IllegalArgumentException("Value must be UTF-8 only: " + key);
    }
    // TODO check for delimiter byte (need to know delim byte!)
    // TODO check for illegal adds (params, previous names)
    return raw.put(key, value);
  }

  public String getKeyword(FCS_KEYWORD keyword) {
    return getKeyword(keyword.keyword);
  }

  public String getKeyword(String keyword) {
    // TODO error checking / string checking (has $?)
    return raw.get(keyword).trim();
  }

  public boolean hasKeyword(FCS_KEYWORD keyword) {
    return hasKeyword(keyword.keyword);
  }

  public boolean hasKeyword(String keyword) {
    // TODO error checking / string checking (has $?)
    return raw.containsKey(keyword);
  }

  public int getKeywordInt(String keyword) {
    // TODO lots of potential for error checking
    return Integer.parseInt(raw.get(keyword).trim());
  }

  public long getKeywordLong(String keyword) {
    // TODO lots of potential for error checking
    return Long.parseLong(raw.get(keyword).trim());
  }

  public byte[] constructRaw() throws IOException {
    return constructRaw(this.delimByte);
  }

  public byte[] constructRaw(byte delimByte) throws IOException {
    ByteArrayOutputStream out = new ByteArrayOutputStream();
    for (Entry<String, String> r : raw.entrySet()) {
      out.write(delimByte);
      out.write(r.getKey().getBytes());
      out.write(delimByte);
      out.write(r.getValue().getBytes());
    }
    out.write(delimByte);
    return out.toByteArray();
  }

  private void parseRaw(byte[] textBytes) throws UnsupportedEncodingException {
    byte delimByte = textBytes[0];

    raw = new LinkedHashMap<>();

    String v = null;
    byte[] sub;
    boolean isKey = true;
    int s = 1, e = 1;

    while (s < textBytes.length) {
      while (true) {
        e++;
        if (textBytes[e] == delimByte) {
          if (e == textBytes.length - 1 || textBytes[e + 1] != delimByte) {
            // delims can exist in values if doubled
            break;
          }
        }
      }
      sub = new byte[e - s];
      System.arraycopy(textBytes, s, sub, 0, e - s);
      if (!isKey) {
        // values can be UTF-8
        raw.put(v, new String(sub, isKey ? "US-ASCII" : "UTF-8"));
      } else {
        // keys are ascii
        v = new String(sub, "US-ASCII");
      }
      isKey = !isKey;
      s = e + 1;
      e = e + 1;
    }
  }

  private void parseParamNames() {

    Map<String, String> paramNames = new HashMap<>();
    for (String k : raw.keySet()) {
      FCS_KEYWORD v = FCS_KEYWORD.findKey(k);
      if (v == FCS_KEYWORD.PnN) {
        paramNames.put(k, raw.get(k).trim());
      }
      // System.out.println(k + "\t" + (v != null ? v.keyword : "null") + "\t" + raw.get(k));
    }

    setParameterNames(new String[paramNames.size()]);
    parameterNamesLong = new String[paramNames.size()];
    for (int i = 1; i < getParameterNames().length + 1; i++) {
      getParameterNames()[i - 1] = paramNames.get(FCS_KEYWORD.PnN.getKey(i));
      parameterNamesLong[i - 1] =
          hasKeyword(FCS_KEYWORD.PnS.getKey(i))
              ? getKeyword(FCS_KEYWORD.PnS.getKey(i))
              : getParameterNames()[i - 1];
    }
  }

  String[] getParameterNames() {
    return parameterNames;
  }

  void setParameterNames(String[] parameterNames) {
    this.parameterNames = parameterNames;
  }
}
