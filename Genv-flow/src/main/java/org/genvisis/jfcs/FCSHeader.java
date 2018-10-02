package org.genvisis.jfcs;

import java.nio.charset.StandardCharsets;
import org.pankratzlab.common.ArrayUtils;

public class FCSHeader {

  public String getFileFormat() {
    return fileFormat;
  }

  public int getTextStart() {
    return textStart;
  }

  public int getTextStop() {
    return textStop;
  }

  public int getDataStart() {
    return dataStart;
  }

  public int getDataStop() {
    return dataStop;
  }

  public int getAnalysisStart() {
    return analysisStart;
  }

  public int getAnalysisStop() {
    return analysisStop;
  }

  public static final int HEADER_NUM_BYTES = 58;
  private static final int TEXT_START_END_BYTE = 17;
  private static final int TEXT_END_END_BYTE = 25;
  private static final int DATA_START_END_BYTE = 33;
  private static final int DATA_END_END_BYTE = 41;
  private static final int ANALYSIS_START_END_BYTE = 49;
  private static final int ANALYSIS_END_END_BYTE = 57;

  public byte[] getBytes() {
    byte[] byt = ArrayUtils.byteArray(HEADER_NUM_BYTES, (byte) 32);
    System.arraycopy(fileFormat.getBytes(StandardCharsets.US_ASCII), 0, byt, 0, 6);

    byte[] t1 = Integer.toString(textStart).getBytes(StandardCharsets.US_ASCII);
    System.arraycopy(t1, 0, byt, TEXT_START_END_BYTE - t1.length, t1.length);

    byte[] t2 = Integer.toString(textStop).getBytes(StandardCharsets.US_ASCII);
    System.arraycopy(t2, 0, byt, TEXT_END_END_BYTE - t2.length, t2.length);

    byte[] d1 = Integer.toString(dataStart).getBytes(StandardCharsets.US_ASCII);
    System.arraycopy(d1, 0, byt, DATA_START_END_BYTE - d1.length, d1.length);

    byte[] d2 = Integer.toString(dataStop).getBytes(StandardCharsets.US_ASCII);
    System.arraycopy(d2, 0, byt, DATA_END_END_BYTE - d2.length, d2.length);

    byte[] a1 = Integer.toString(analysisStart).getBytes(StandardCharsets.US_ASCII);
    System.arraycopy(a1, 0, byt, ANALYSIS_START_END_BYTE - a1.length, a1.length);

    byte[] a2 = Integer.toString(analysisStop).getBytes(StandardCharsets.US_ASCII);
    System.arraycopy(a2, 0, byt, ANALYSIS_END_END_BYTE - a2.length, a2.length);

    return byt;
  }

  String fileFormat;
  int textStart;
  int textStop;
  int dataStart;
  int dataStop;
  int analysisStart;
  int analysisStop;
}
