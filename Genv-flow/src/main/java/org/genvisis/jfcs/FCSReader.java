package org.genvisis.jfcs;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteOrder;
import java.util.HashSet;
import java.util.Set;
import java.util.regex.Pattern;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.genvisis.bgen.BGENBitMath;
import org.genvisis.common.ext;

public class FCSReader {

  String filePath;
  private RandomAccessFile fileObj;

  private FCSHeader header;
  private FCSKeywords keys;
  private FCSSpillover spill;
  FCSData data;

  private FCSReader() {}

  public static FCSReader open(String file) throws IOException {
    // check if file exists / can be read?
    FCSReader reader = new FCSReader();
    reader.setFile(file);
    reader.openFile();
    reader.readHeader();
    reader.checkHeader();
    reader.readText();
    reader.setSpill();
    reader.initData();
    reader.data.readData(reader);
    reader.closeInternal();
    // reader.readAnalysis();
    return reader;
  }

  private void setFile(String file) {
    this.filePath = file;
  }

  private void openFile() throws FileNotFoundException {
    setSystemFile(new RandomAccessFile(filePath, "r"));
  }

  private void readHeader() throws IOException {
    if (getHeader() != null) return;

    byte[] headBytes = new byte[58];
    int read = getSystemFile().read(headBytes);
    if (read != headBytes.length) {
      throw new IOException("");
    }

    String hdrStr = new String(headBytes, "US-ASCII");
    String[] pts = hdrStr.split("\\s+");

    setHeader(new FCSHeader());
    getHeader().fileFormat = pts[0];
    getHeader().textStart = Integer.parseInt(pts[1]);
    getHeader().textStop = Integer.parseInt(pts[2]);
    if (pts.length >= 7) {
      getHeader().dataStart = Integer.parseInt(pts[3]);
      getHeader().dataStop = Integer.parseInt(pts[4]);
      getHeader().analysisStart = Integer.parseInt(pts[5]);
      getHeader().analysisStop = Integer.parseInt(pts[6]);
    }
  }

  private void checkHeader() throws IOException {
    switch (getHeader().fileFormat) {
      case "FCS3.1":
        break;
      case "FCS3.0":
        break;
      case "FCS2.0":
      case "FCS1.0":
      default:
        throw new IOException("Unsupported FCS File format: " + getHeader().fileFormat);
    }
  }

  private void readText() throws IOException {
    if (getSystemFile().getFilePointer() != getHeader().textStart) {
      getSystemFile().seek(getHeader().textStart);
    }
    byte[] textBytes = new byte[(getHeader().textStop + 1) - getHeader().textStart];
    getSystemFile().read(textBytes);
    this.setKeywords(new FCSKeywords(textBytes));
  }

  private void setSpill() {
    FCS_KEYWORD s = getKeywords().hasKeyword(FCS_KEYWORD.SPILLOVER) ? FCS_KEYWORD.SPILLOVER
                                                                    : getKeywords().hasKeyword(FCS_KEYWORD.SPILL) ? FCS_KEYWORD.SPILL
                                                                                                                  : null;
    this.setSpillover(s == null ? null : FCSSpillover.parse(getKeywords().getKeyword(s)));
  }

  private void initData() {
    char dataType = getKeywords().getKeyword(FCS_KEYWORD.DATATYPE).charAt(0);
    switch (dataType) {
      case 'D':
        data = new FCSDoubleData();
        break;
      case 'F':
        data = new FCSDoubleData();
        break;
      case 'A':
      case 'I':
      default:
        throw new RuntimeException("Data type {" + dataType + "} not supported.");
    }
  }

  private void closeInternal() throws IOException {
    this.fileObj.close();
  }

  public void dispose() {
    header = null;
    keys = null;
    spill = null;
    data.dispose();
    data = null;
  }

  public FCSKeywords getKeywords() {
    return keys;
  }

  private void setKeywords(FCSKeywords keys) {
    this.keys = keys;
  }

  protected RandomAccessFile getSystemFile() {
    return fileObj;
  }

  protected void setSystemFile(RandomAccessFile fileObj) {
    this.fileObj = fileObj;
  }

  public FCSHeader getHeader() {
    return header;
  }

  private void setHeader(FCSHeader header) {
    this.header = header;
  }

  public FCSSpillover getSpillover() {
    return spill;
  }

  private void setSpillover(FCSSpillover spill) {
    this.spill = spill;
  }

  public double[] getParamAsDoubles(String paramName, boolean compensated) {
    return data.getParamData(ext.indexOfStr(paramName,
                                            compensated ? spill.getParameterNames()
                                                        : keys.getParameterNames()),
                             compensated, (double[]) null);
  }

  public double[] getEventAsDoubles(int event, boolean getCompensated) {
    return data.getEventData(event, getCompensated, (double[]) null);
  }

}

interface FCSData {

  void readData(FCSReader reader) throws IOException;

  void dispose();

  boolean hasCompensatedData();

  default public float[] getParamData(int paramIndex, boolean compensated, float[] emptyData) {
    throw new UnsupportedOperationException();
  }

  default public double[] getParamData(int paramIndex, boolean compensated, double[] emptyData) {
    throw new UnsupportedOperationException();
  }

  default public float[] getEventData(int eventIndex, boolean compensated, float[] emptyData) {
    throw new UnsupportedOperationException();
  }

  default public double[] getEventData(int eventIndex, boolean compensated, double[] emptyData) {
    throw new UnsupportedOperationException();
  }

}

class FCSDoubleData implements FCSData {

  double[][] data;
  double[][] compData;

  @Override
  public double[] getParamData(int paramIndex, boolean compensated, double[] emptyData) {
    if (paramIndex >= (compensated ? compData : data).length || paramIndex < 0) {
      throw new ArrayIndexOutOfBoundsException("Invalid index: " + paramIndex);
    }
    if (emptyData != null) {
      if (emptyData.length != (compensated ? compData : data)[paramIndex].length) {
        throw new IllegalArgumentException("Expected data is wrong length; given: "
                                           + emptyData.length + "; expected: "
                                           + (compensated ? compData : data)[paramIndex].length);
      }
    } else {
      emptyData = new double[(compensated ? compData : data)[paramIndex].length];
    }
    System.arraycopy((compensated ? compData : data)[paramIndex], 0, emptyData, 0,
                     emptyData.length);
    return emptyData;
  }

  @Override
  public float[] getParamData(int paramIndex, boolean compensated, float[] emptyData) {
    if (paramIndex >= (compensated ? compData : data).length || paramIndex < 0) {
      throw new ArrayIndexOutOfBoundsException("Invalid index: " + paramIndex);
    }
    if (emptyData != null) {
      if (emptyData.length != (compensated ? compData : data)[paramIndex].length) {
        throw new IllegalArgumentException("Expected data is wrong length; given: "
                                           + emptyData.length + "; expected: "
                                           + (compensated ? compData : data)[paramIndex].length);
      }
    } else {
      emptyData = new float[(compensated ? compData : data)[paramIndex].length];
    }
    for (int e = 0; e < emptyData.length; e++) {
      // LOSING PRECISION!
      emptyData[e] = (float) (compensated ? compData : data)[paramIndex][e];
    }
    return emptyData;
  }

  @Override
  public double[] getEventData(int eventIndex, boolean compensated, double[] emptyData) {
    if (eventIndex >= (compensated ? compData : data)[0].length || eventIndex < 0) {
      throw new ArrayIndexOutOfBoundsException("Invalid index: " + eventIndex);
    }
    if (emptyData != null) {
      if (emptyData.length != (compensated ? compData : data).length) {
        throw new IllegalArgumentException("Expected data is wrong length; given: "
                                           + emptyData.length + "; expected: "
                                           + (compensated ? compData : data).length);
      }
    } else {
      emptyData = new double[(compensated ? compData : data).length];
    }
    for (int p = 0; p < (compensated ? compData : data).length; p++) {
      emptyData[p] = (compensated ? compData : data)[p][eventIndex];
    }
    return emptyData;
  }

  @Override
  public float[] getEventData(int eventIndex, boolean compensated, float[] emptyData) {
    if (eventIndex >= (compensated ? compData : data)[0].length || eventIndex < 0) {
      throw new ArrayIndexOutOfBoundsException("Invalid index: " + eventIndex);
    }
    if (emptyData != null) {
      if (emptyData.length != (compensated ? compData : data).length) {
        throw new IllegalArgumentException("Expected data is wrong length; given: "
                                           + emptyData.length + "; expected: "
                                           + (compensated ? compData : data).length);
      }
    } else {
      emptyData = new float[(compensated ? compData : data).length];
    }
    for (int p = 0; p < (compensated ? compData : data).length; p++) {
      // LOSING PRECISION!
      emptyData[p] = (float) (compensated ? compData : data)[p][eventIndex];
    }
    return emptyData;
  }

  @Override
  public boolean hasCompensatedData() {
    return compData != null;
  }

  @Override
  public void dispose() {
    this.data = null;
    this.compData = null;
  }

  public void readData(FCSReader reader) throws IOException {
    if (this.data != null) return;
    long dataStart = reader.getHeader().dataStart;
    if (dataStart <= 0) {
      dataStart = reader.getKeywords().getKeywordLong(FCS_KEYWORD.BEGINDATA.keyword);
    }
    long dataEnd = reader.getHeader().dataStop;
    if (dataEnd <= 0) {
      dataEnd = reader.getKeywords().getKeywordLong(FCS_KEYWORD.ENDDATA.keyword);
    }
    if (reader.getSystemFile().getFilePointer() != dataStart) {
      reader.getSystemFile().seek(dataStart);
    }
    int events = reader.getKeywords().getEventCount();
    int params = reader.getKeywords().getParameterCount();
    ByteOrder end = reader.getKeywords().getEndianness();

    int[] paramBytes = new int[params];
    int[] paramByteInd = new int[params];
    int byteSum = 0;
    for (int i = 0; i < params; i++) {
      int p = reader.getKeywords()
                    .getKeywordInt(FCS_KEYWORD.PnB.keyword.replace("n", Integer.toString(i + 1)));
      p /= 8; // require discrete bytes
      paramByteInd[i] = byteSum;
      byteSum += p;
      paramBytes[i] = p;
    }

    long dataBytes = events * byteSum;
    long dataLength = (dataEnd - dataStart + 1);
    if (dataBytes != dataLength) {
      // TODO error
    }

    FCSSpillover spill = reader.getSpillover();
    DenseMatrix64F spillMatrix = new DenseMatrix64F(spill.getCoefficients());
    CommonOps.invert(spillMatrix);

    data = new double[params][events];
    compData = new double[spill.getParameterNames().length][events];
    int[] compInds = new int[reader.getKeywords().getParameterNames().length];
    for (int i = 0; i < compInds.length; i++) {
      compInds[i] = ext.indexOfStr(reader.getKeywords().getParameterNames()[i],
                                   spill.getParameterNames());
    }

    double[] paramScaling = reader.getKeywords().getParamScaling();

    byte[] parse = new byte[byteSum];
    for (int e = 0; e < events; e++) {
      reader.getSystemFile().read(parse);
      for (int p = 0; p < params; p++) {
        data[p][e] = (paramBytes[p] == 64 ? BGENBitMath.bytesToDouble(end == ByteOrder.LITTLE_ENDIAN,
                                                                      parse[paramByteInd[p]],
                                                                      parse[paramByteInd[p] + 1],
                                                                      parse[paramByteInd[p] + 2],
                                                                      parse[paramByteInd[p] + 3],
                                                                      parse[paramByteInd[p] + 4],
                                                                      parse[paramByteInd[p] + 5],
                                                                      parse[paramByteInd[p] + 6],
                                                                      parse[paramByteInd[p] + 7])
                                          // double types /should/ be 64 bits, but not always true
                                          : BGENBitMath.bytesToFloat(end == ByteOrder.LITTLE_ENDIAN,
                                                                     parse[paramByteInd[p]],
                                                                     parse[paramByteInd[p] + 1],
                                                                     parse[paramByteInd[p] + 2],
                                                                     parse[paramByteInd[p] + 3]))
                     / paramScaling[p];
      }
      for (int p = 0; p < params; p++) {
        if (compInds[p] == -1) {
          continue;
        }
        float sum = 0;
        for (int p1 = 0; p1 < params; p1++) {
          if (compInds[p1] == -1) {
            continue;
          }
          sum += data[p1][e] * spillMatrix.get(compInds[p1], compInds[p]);
        }
        compData[compInds[p]][e] = sum;
      }
    }
    parse = null;
  }

}

class FCSFloatData implements FCSData {

  float[][] data;
  float[][] compData;

  @Override
  public double[] getParamData(int paramIndex, boolean compensated, double[] emptyData) {
    if (paramIndex >= (compensated ? compData : data).length || paramIndex < 0) {
      throw new ArrayIndexOutOfBoundsException("Invalid index: " + paramIndex);
    }
    if (emptyData != null) {
      if (emptyData.length != (compensated ? compData : data)[paramIndex].length) {
        throw new IllegalArgumentException("Expected data is wrong length; given: "
                                           + emptyData.length + "; expected: "
                                           + (compensated ? compData : data)[paramIndex].length);
      }
    } else {
      emptyData = new double[(compensated ? compData : data)[paramIndex].length];
    }
    for (int e = 0; e < emptyData.length; e++) {
      // GAINING PRECISION!
      emptyData[e] = (float) (compensated ? compData : data)[paramIndex][e];
    }
    return emptyData;
  }

  @Override
  public float[] getParamData(int paramIndex, boolean compensated, float[] emptyData) {
    if (paramIndex >= (compensated ? compData : data).length || paramIndex < 0) {
      throw new ArrayIndexOutOfBoundsException("Invalid index: " + paramIndex);
    }
    if (emptyData != null) {
      if (emptyData.length != (compensated ? compData : data)[paramIndex].length) {
        throw new IllegalArgumentException("Expected data is wrong length; given: "
                                           + emptyData.length + "; expected: "
                                           + (compensated ? compData : data)[paramIndex].length);
      }
    } else {
      emptyData = new float[(compensated ? compData : data)[paramIndex].length];
    }
    System.arraycopy((compensated ? compData : data)[paramIndex], 0, emptyData, 0,
                     emptyData.length);
    return emptyData;
  }

  @Override
  public double[] getEventData(int eventIndex, boolean compensated, double[] emptyData) {
    if (eventIndex >= (compensated ? compData : data)[0].length || eventIndex < 0) {
      throw new ArrayIndexOutOfBoundsException("Invalid index: " + eventIndex);
    }
    if (emptyData != null) {
      if (emptyData.length != (compensated ? compData : data).length) {
        throw new IllegalArgumentException("Expected data is wrong length; given: "
                                           + emptyData.length + "; expected: "
                                           + (compensated ? compData : data).length);
      }
    } else {
      emptyData = new double[(compensated ? compData : data).length];
    }
    for (int p = 0; p < (compensated ? compData : data).length; p++) {
      // GAINING PRECISION!
      emptyData[p] = (compensated ? compData : data)[p][eventIndex];
    }
    return emptyData;
  }

  @Override
  public float[] getEventData(int eventIndex, boolean compensated, float[] emptyData) {
    if (eventIndex >= (compensated ? compData : data)[0].length || eventIndex < 0) {
      throw new ArrayIndexOutOfBoundsException("Invalid index: " + eventIndex);
    }
    if (emptyData != null) {
      if (emptyData.length != (compensated ? compData : data).length) {
        throw new IllegalArgumentException("Expected data is wrong length; given: "
                                           + emptyData.length + "; expected: "
                                           + (compensated ? compData : data).length);
      }
    } else {
      emptyData = new float[(compensated ? compData : data).length];
    }
    for (int p = 0; p < (compensated ? compData : data).length; p++) {
      emptyData[p] = (compensated ? compData : data)[p][eventIndex];
    }
    return emptyData;
  }

  @Override
  public boolean hasCompensatedData() {
    return compData != null;
  }

  @Override
  public void dispose() {
    this.data = null;
    this.compData = null;
  }

  public void readData(FCSReader reader) throws IOException {
    if (this.data != null) return;
    long dataStart = reader.getHeader().dataStart;
    if (dataStart <= 0) {
      dataStart = reader.getKeywords().getKeywordLong(FCS_KEYWORD.BEGINDATA.keyword);
    }
    long dataEnd = reader.getHeader().dataStop;
    if (dataEnd <= 0) {
      dataEnd = reader.getKeywords().getKeywordLong(FCS_KEYWORD.ENDDATA.keyword);
    }
    if (reader.getSystemFile().getFilePointer() != dataStart) {
      reader.getSystemFile().seek(dataStart);
    }
    int events = reader.getKeywords().getEventCount();
    int params = reader.getKeywords().getParameterCount();
    ByteOrder end = reader.getKeywords().getEndianness();

    int[] paramBytes = new int[params];
    int[] paramByteInd = new int[params];
    int byteSum = 0;
    for (int i = 0; i < params; i++) {
      int p = reader.getKeywords()
                    .getKeywordInt(FCS_KEYWORD.PnB.keyword.replace("n", Integer.toString(i + 1)));
      p /= 8; // require discrete bytes
      paramByteInd[i] = byteSum;
      byteSum += p;
      paramBytes[i] = p;
    }

    long dataBytes = events * byteSum;
    long dataLength = (dataEnd - dataStart + 1);
    if (dataBytes != dataLength) {
      // TODO error
    }

    double[] paramScaling = reader.getKeywords().getParamScaling();

    FCSSpillover spill = reader.getSpillover();
    DenseMatrix64F spillMatrix = new DenseMatrix64F(spill.getCoefficients());
    CommonOps.invert(spillMatrix);

    long t = System.nanoTime();
    data = new float[params][events];
    compData = new float[spill.getParameterNames().length][events];
    int[] compInds = new int[reader.getKeywords().getParameterNames().length];
    for (int i = 0; i < compInds.length; i++) {
      compInds[i] = ext.indexOfStr(reader.getKeywords().getParameterNames()[i],
                                   spill.getParameterNames());
    }

    byte[] parse = new byte[byteSum];
    for (int e = 0; e < events; e++) {
      reader.getSystemFile().read(parse);
      for (int p = 0; p < params; p++) {
        data[e][p] = (float) (BGENBitMath.bytesToFloat(end == ByteOrder.LITTLE_ENDIAN,
                                                       parse[paramByteInd[p]],
                                                       parse[paramByteInd[p] + 1],
                                                       parse[paramByteInd[p] + 2],
                                                       parse[paramByteInd[p] + 3])
                              / paramScaling[p]);
      }
      for (int p = 0; p < params; p++) {
        if (compInds[p] == -1) {
          continue;
        }
        float sum = 0;
        for (int p1 = 0; p1 < params; p1++) {
          if (compInds[p1] == -1) {
            continue;
          }
          sum += data[p1][e] * spillMatrix.get(compInds[p], compInds[p1]);
        }
        compData[compInds[p]][e] = sum;
      }
    }
    parse = null;
    System.out.println("Read " + events + " events for " + params + " params in "
                       + ext.getTimeElapsedNanos(t));

  }

}

enum FCS_KEYWORD {
  // required by spec:
  BEGINANALYSIS("$BEGINANALYSIS", Pattern.quote("$BEGINANALYSIS"), true),
  BEGINDATA("$BEGINDATA", Pattern.quote("$BEGINDATA"), true),
  BEGINSTEXT("$BEGINSTEXT", Pattern.quote("$BEGINSTEXT"), true),
  BYTEORD("$BYTEORD", Pattern.quote("$BYTEORD"), true),
  DATATYPE("$DATATYPE", Pattern.quote("$DATATYPE"), true),
  ENDANALYSIS("$ENDANALYSIS", Pattern.quote("$ENDANALYSIS"), true),
  ENDDATA("$ENDDATA", Pattern.quote("$ENDDATA"), true),
  ENDSTEXT("$ENDSTEXT", Pattern.quote("$ENDSTEXT"), true),
  MODE("$MODE", Pattern.quote("$MODE"), true),
  NEXTDATA("$NEXTDATA", Pattern.quote("$NEXTDATA"), true),
  PAR("$PAR", Pattern.quote("$PAR"), true),
  PnB("$PnB", "\\$P\\d+B", true),
  PnE("$PnE", "\\$P\\d+E", true),
  PnN("$PnN", "\\$P\\d+N", true),
  PnR("$PnR", "\\$P\\d+R", true),
  TOT("$TOT", "\\$TOT", true),

  ABRT("$ABRT", "\\$ABRT", false),
  BTIM("$BTIM", "\\$BTIM", false),
  CELLS("$CELLS", "\\$CELLS", false),
  COM("$COM", "\\$COM", false),
  CSMODE("$CSMODE", "\\$CSMODE", false),
  CSVBITS("$CSVBITS", "\\$CSVBITS", false),
  CSVnFLAG("$CSVnFLAG", "\\$CSV\\d+FLAG", false),
  CYT("$CYT", "\\$CYT", false),
  CYTSN("$CYTSN", "\\$CYTSN", false),
  DATE("$DATE", "\\$DATE", false),
  ETIM("$ETIM", "\\$ETIM", false),
  EXP("$EXP", "\\$EXP", false),
  FIL("$FIL", "\\$FIL", false),
  GATE("$GATE", "\\$GATE", false),
  GATING("$GATING", "\\$GATING", false),
  GnE("$GnE", "\\$G\\d+E", false),
  GnF("$GnF", "\\$G\\d+F", false),
  GnN("$GnN", "\\$G\\d+N", false),
  GnP("$GnP", "\\$G\\d+P", false),
  GnR("$GnR", "\\$G\\d+R", false),
  GnS("$GnS", "\\$G\\d+S", false),
  GnT("$GnT", "\\$G\\d+T", false),
  GnV("$GnV", "\\$G\\d+V", false),
  INST("$INST", "\\$INST", false),
  LAST_MODIFIED("$LAST_MODIFIED", "\\$LAST_MODIFIED", false),
  LAST_MODIFIER("$LAST_MODIFIER", "\\$LAST_MODIFIER", false),
  LOST("$LOST", "\\$LOST", false),
  OP("$OP", "\\$OP", false),
  ORIGINALITY("$ORIGINALITY", "\\$ORIGINALITY", false),
  PKn("$PKn", "\\$PK\\d+", false),
  PKNn("$PKNn", "\\$PKN\\d+", false),
  PLATEID("$PLATEID", "\\$PLATEID", false),
  PLATENAME("$PLATENAME", "\\$PLATENAME", false),
  PnCALIBRATION("$PnCALIBRATION", "\\$P\\d+CALIBRATION", false),
  PnD("$PnD", "\\$P\\d+D", false),
  PnF("$PnF", "\\$P\\d+F", false),
  PnG("$PnG", "\\$P\\d+G", false),
  PnL("$PnL", "\\$P\\d+L", false),
  PnO("$PnO", "\\$P\\d+O", false),
  PnP("$PnP", "\\$P\\d+P", false),
  PnS("$PnS", "\\$P\\d+S", false),
  PnT("$PnT", "\\$P\\d+T", false),
  PnV("$PnV", "\\$P\\d+V", false),
  PROJ("$PROJ", "\\$PROJ", false),
  RnI("$RnI", "\\$R\\d+I", false),
  RnW("$RnW", "\\$R\\d+W", false),
  SMNO("$SMNO", "\\$SMNO", false),
  SPILLOVER("$SPILLOVER", "\\$SPILLOVER", false),
  SRC("$SRC", "\\$SRC", false),
  SYS("$SYS", "\\$SYS", false),
  TIMESTEP("$TIMESTEP", "\\$TIMESTEP", false),
  TR("$TR", "\\$TR", false),
  VOL("$VOL", "\\$VOL", false),
  WELLID("$WELLID", "\\$WELLID", false),
  // custom
  SPILL("SPILL", "\\SPILL", false);

  FCS_KEYWORD(String key, String pattern, boolean required) {
    this.required = required;
    this.keyword = key;
    this.patternStr = pattern;
    this.pattern = Pattern.compile(patternStr);
    this.isInd = key.contains("n");
  }

  static Set<FCS_KEYWORD> getRequiredKeywords() {
    HashSet<FCS_KEYWORD> keys = new HashSet<>();
    for (FCS_KEYWORD k : values()) {
      if (k.required) {
        keys.add(k);
      }
    }
    return keys;
  }

  boolean required;
  String keyword;
  String patternStr;
  boolean isInd;
  Pattern pattern;

  String getKey() {
    return keyword;
  }

  String getKey(int index) {
    return isInd ? keyword.replace("n", Integer.toString(index)) : keyword;
  }

  static FCS_KEYWORD findKey(String poss) {
    try {
      // check if valid
      FCS_KEYWORD word = FCS_KEYWORD.valueOf(poss.startsWith("$") ? poss.substring(1) : poss);
      return word;
    } catch (IllegalArgumentException e) {
      for (FCS_KEYWORD word : values()) {
        if (!word.isInd) {
          continue;
        } else if (word.pattern.matcher(poss).matches()) {
          return word;
        }
      }
      return null;
    }
  }

  public static boolean isKey(String poss) {
    return findKey(poss) != null;
  }

}

class FCSHeader {

  String fileFormat;
  int textStart;
  int textStop;
  int dataStart;
  int dataStop;
  int analysisStart;
  int analysisStop;
}
