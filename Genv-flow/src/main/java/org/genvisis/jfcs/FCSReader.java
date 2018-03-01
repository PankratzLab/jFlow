package org.genvisis.jfcs;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteOrder;
import java.util.Arrays;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.genvisis.bgen.BGENBitMath;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.ext;

public class FCSReader {

  String filePath;
  private RandomAccessFile fileObj;

  private String fcsVersion;
  private FCSKeywords keys;
  private FCSSpillover spill;
  FCSData data;

  private FCSReader() {}

  public static FCSKeywords readKeywords(String file) throws IOException {
    FCSReader reader = new FCSReader();
    reader.setFile(file);
    reader.openFile();
    FCSHeader header = reader.readHeader();
    reader.setVersion(header.getFileFormat());
    reader.checkHeader(header);
    reader.readText(header);
    FCSKeywords keys = reader.getKeywords();
    reader.closeInternal();
    reader.dispose();
    return keys;
  }

  public static FCSReader open(String file) throws IOException {
    // check if file exists / can be read?
    FCSReader reader = new FCSReader();
    reader.setFile(file);
    reader.openFile();
    FCSHeader header = reader.readHeader();
    reader.setVersion(header.getFileFormat());
    reader.checkHeader(header);
    reader.readText(header);
    reader.setSpill();
    reader.initData();
    reader.data.readData(reader, header);
    reader.closeInternal();
    // reader.readAnalysis();
    return reader;
  }

  public static void write(FCSReader reader, String newFile) throws IOException {
    RandomAccessFile fO = new RandomAccessFile(newFile, "rw");

    long dataLength = reader.data.getDataByteCount();
    byte[] textBytes;

    int dB = 5000;

    // placeholders
    reader.getKeywords().setKeyword(FCS_KEYWORD.BEGINDATA, "5000");
    reader.getKeywords().setKeyword(FCS_KEYWORD.ENDDATA, Long.toString(5000 + dataLength - 1));
    textBytes = reader.getKeywords().constructRaw();

    while ((textBytes.length + 256 - 1 + 5) >= dB) {
      dB += 1000;
      if (dB > 9999) {
        reader.getKeywords().setKeyword(FCS_KEYWORD.BEGINDATA, Integer.toString(dB));
        reader.getKeywords().setKeyword(FCS_KEYWORD.ENDDATA, Long.toString(dB + dataLength - 1));
        textBytes = reader.getKeywords().constructRaw();
      }
    }

    long dE = dB + dataLength - 1;
    reader.getKeywords().setKeyword(FCS_KEYWORD.BEGINDATA, Integer.toString(dB));
    reader.getKeywords().setKeyword(FCS_KEYWORD.ENDDATA, Long.toString(dE));
    textBytes = reader.getKeywords().constructRaw();

    FCSHeader newHeader = new FCSHeader();
    newHeader.fileFormat = reader.getVersion();

    newHeader.textStart = 256; // apparent default?
    newHeader.textStop = newHeader.textStart + textBytes.length - 1;

    if (newHeader.textStop + 5 + dataLength - 1 > 99999999) {
      newHeader.dataStart = 0;
      newHeader.dataStop = 0;
    } else {
      newHeader.dataStart = newHeader.textStop + 5;
      newHeader.dataStop = (int) (newHeader.dataStart + dataLength - 1);
    }

    // analysis data not supported
    // TODO probably should enforce keywords == 0 also
    newHeader.analysisStart = 0;
    newHeader.analysisStop = 0;

    // write header
    fO.write(newHeader.getBytes());
    // write spaces in between
    byte[] filler = ArrayUtils.byteArray(198, (byte) 32); // 256 - 58 = 198
    fO.write(filler);
    // write text
    fO.write(textBytes);
    // write spaces in between
    filler = ArrayUtils.byteArray(newHeader.dataStart - (newHeader.textStop + 1), (byte) 32);
    fO.write(filler);
    // write data
    reader.data.writeData(reader, fO);
    // write spaces in between
    // write analysis
    fO.close();
  }

  private String getVersion() {
    return fcsVersion;
  }

  private void setVersion(String v) {
    this.fcsVersion = v;
  }

  private void setFile(String file) {
    this.filePath = file;
  }

  private void openFile() throws FileNotFoundException {
    setSystemFile(new RandomAccessFile(filePath, "r"));
  }

  private FCSHeader readHeader() throws IOException {
    FCSHeader header = new FCSHeader();

    byte[] headBytes = new byte[FCSHeader.HEADER_NUM_BYTES];
    // TODO check if system file is open / closed
    // TODO check position

    int read = getSystemFile().read(headBytes);
    if (read != headBytes.length) {
      throw new IOException("Failed to read the proper number of header bytes (expected "
                            + FCSHeader.HEADER_NUM_BYTES + "; found " + read);
    }

    header.fileFormat = new String(Arrays.copyOfRange(headBytes, 0, 6), "US-ASCII").trim();
    header.textStart = Integer.parseInt(new String(Arrays.copyOfRange(headBytes, 10, 18),
                                                   "US-ASCII").trim());
    header.textStop = Integer.parseInt(new String(Arrays.copyOfRange(headBytes, 18, 26),
                                                  "US-ASCII").trim());
    header.dataStart = Integer.parseInt(new String(Arrays.copyOfRange(headBytes, 26, 34),
                                                   "US-ASCII").trim());
    header.dataStop = Integer.parseInt(new String(Arrays.copyOfRange(headBytes, 34, 42),
                                                  "US-ASCII").trim());
    header.analysisStart = Integer.parseInt(new String(Arrays.copyOfRange(headBytes, 42, 50),
                                                       "US-ASCII").trim());
    header.analysisStop = Integer.parseInt(new String(Arrays.copyOfRange(headBytes, 50, 58),
                                                      "US-ASCII").trim());
    return header;
  }

  private void checkHeader(FCSHeader header) throws IOException {
    switch (header.getFileFormat()) {
      case "FCS3.1":
        break;
      case "FCS3.0":
        break;
      case "FCS2.0":
      case "FCS1.0":
      default:
        throw new IOException("Unsupported FCS File format: " + header.getFileFormat());
    }
  }

  private void readText(FCSHeader header) throws IOException {
    if (getSystemFile().getFilePointer() != header.getTextStart()) {
      getSystemFile().seek(header.getTextStart());
    }
    byte[] textBytes = new byte[(header.getTextStop() + 1) - header.getTextStart()];
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
    keys = null;
    spill = null;
    if (data != null) {
      data.dispose();
      data = null;
    }
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

  void readData(FCSReader reader, FCSHeader header) throws IOException;

  long getDataByteCount();

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

  void writeData(FCSReader reader, RandomAccessFile file) throws IOException;

}

class FCSDoubleData implements FCSData {

  int[] paramBytes;
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

  @Override
  public long getDataByteCount() {
    return data[0].length * (long) ArrayUtils.sum(paramBytes);
  }

  public void readData(FCSReader reader, FCSHeader header) throws IOException {
    if (this.data != null) return;
    long dataStart = header.getDataStart();
    if (dataStart <= 0) {
      dataStart = reader.getKeywords().getKeywordLong(FCS_KEYWORD.BEGINDATA.keyword);
    }
    long dataEnd = header.getDataStop();
    if (dataEnd <= 0) {
      dataEnd = reader.getKeywords().getKeywordLong(FCS_KEYWORD.ENDDATA.keyword);
    }
    if (reader.getSystemFile().getFilePointer() != dataStart) {
      reader.getSystemFile().seek(dataStart);
    }
    int events = reader.getKeywords().getEventCount();
    int params = reader.getKeywords().getParameterCount();
    ByteOrder end = reader.getKeywords().getEndianness();

    paramBytes = new int[params];
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
        data[p][e] = (paramBytes[p] == 8 ? BGENBitMath.bytesToDouble(end == ByteOrder.LITTLE_ENDIAN,
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

  public void writeData(FCSReader reader, RandomAccessFile newFile) throws IOException {
    int params = data.length;
    int[] paramBytes = new int[params];
    int[] paramByteInd = new int[params];
    double[] paramScaling = reader.getKeywords().getParamScaling();
    ByteOrder end = reader.getKeywords().getEndianness();
    int byteSum = 0;
    for (int i = 0; i < params; i++) {
      int p = reader.getKeywords()
                    .getKeywordInt(FCS_KEYWORD.PnB.keyword.replace("n", Integer.toString(i + 1)));
      p /= 8; // require discrete bytes
      paramByteInd[i] = byteSum;
      byteSum += p;
      paramBytes[i] = p;
    }

    byte[] rawLine;
    for (int e = 0; e < data[0].length; e++) {
      rawLine = new byte[byteSum];
      for (int p = 0; p < params; p++) {
        if (paramBytes[p] == 8) {
          System.arraycopy(BGENBitMath.doubleToBytes(end == ByteOrder.LITTLE_ENDIAN,
                                                     data[p][e] * paramScaling[p]),
                           0, rawLine, paramByteInd[p], 8);
        } else {
          System.arraycopy(BGENBitMath.floatToBytes(end == ByteOrder.LITTLE_ENDIAN,
                                                    (float) (data[p][e] * paramScaling[p])),
                           0, rawLine, paramByteInd[p], 4);
        }
      }
      newFile.write(rawLine);
    }

  }

}

class FCSFloatData implements FCSData {

  int[] paramBytes;
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

  @Override
  public long getDataByteCount() {
    return data[0].length * (long) ArrayUtils.sum(paramBytes);
  }

  public void readData(FCSReader reader, FCSHeader header) throws IOException {
    if (this.data != null) return;
    long dataStart = header.getDataStart();
    if (dataStart <= 0) {
      dataStart = reader.getKeywords().getKeywordLong(FCS_KEYWORD.BEGINDATA.keyword);
    }
    long dataEnd = header.getDataStop();
    if (dataEnd <= 0) {
      dataEnd = reader.getKeywords().getKeywordLong(FCS_KEYWORD.ENDDATA.keyword);
    }
    if (reader.getSystemFile().getFilePointer() != dataStart) {
      reader.getSystemFile().seek(dataStart);
    }
    int events = reader.getKeywords().getEventCount();
    int params = reader.getKeywords().getParameterCount();
    ByteOrder end = reader.getKeywords().getEndianness();

    paramBytes = new int[params];
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
        data[p][e] = (float) (BGENBitMath.bytesToFloat(end == ByteOrder.LITTLE_ENDIAN,
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

  public void writeData(FCSReader reader, RandomAccessFile newFile) throws IOException {
    int params = data.length;
    int[] paramByteInd = new int[params];
    double[] paramScaling = reader.getKeywords().getParamScaling();
    ByteOrder end = reader.getKeywords().getEndianness();
    int byteSum = 0;
    for (int i = 0; i < params; i++) {
      int p = reader.getKeywords()
                    .getKeywordInt(FCS_KEYWORD.PnB.keyword.replace("n", Integer.toString(i + 1)));
      p /= 8; // require discrete bytes
      paramByteInd[i] = byteSum;
      byteSum += p;
    }

    byte[] rawLine;
    for (int e = 0; e < data[0].length; e++) {
      rawLine = new byte[byteSum];
      for (int p = 0; p < params; p++) {
        System.arraycopy(BGENBitMath.floatToBytes(end == ByteOrder.LITTLE_ENDIAN,
                                                  (float) (data[p][e] * paramScaling[p])),
                         0, rawLine, paramByteInd[p], 4);
      }
      newFile.write(rawLine);
    }

  }

}
