package org.genvisis.one;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import org.genvisis.CLI;
import org.genvisis.common.Files;

public class UKBBExtract {

  private static void extract(String binFile, int markerStart, int numToExtract,
                              String outputFilename) throws IOException {
    System.out.println("Extracting A/B/X/Y data for variants " + markerStart + " to "
                       + (markerStart + numToExtract));
    int nInd = 488377; // number of individuals in UKBB data
    RandomAccessFile binIn = new RandomAccessFile(binFile, "r");
    PrintWriter writer = Files.getAppropriateWriter(outputFilename);
    writer.println("A\tB\tX\tY");
    int binBlockSize = nInd * 8;
    byte[] intensBytes = new byte[binBlockSize];
    float a, b, x, y;
    double log2 = Math.log(2);
    int logEvery = Math.max(1, numToExtract / 5);
    for (int i = markerStart; i < markerStart + numToExtract; i++) {
      if ((i - markerStart) % logEvery == 0) {
        System.out.println("Wrote marker " + i);
      }
      binIn.seek(((long) i) * ((long) binBlockSize));
      binIn.read(intensBytes);
      for (int bitInd = 0, binInd = 0; bitInd < nInd; bitInd++) {
        byte[] intA = {intensBytes[binInd++], intensBytes[binInd++], intensBytes[binInd++],
                       intensBytes[binInd++]};
        byte[] intB = {intensBytes[binInd++], intensBytes[binInd++], intensBytes[binInd++],
                       intensBytes[binInd++]};

        a = ByteBuffer.wrap(intA).order(ByteOrder.LITTLE_ENDIAN).getFloat();
        b = ByteBuffer.wrap(intB).order(ByteOrder.LITTLE_ENDIAN).getFloat();
        if (a <= 0 || b <= 0) {
          x = Float.NaN;
          y = Float.NaN;
        } else {
          x = (float) (Math.log(a / b) / log2);
          y = (float) ((Math.log(a * b) / log2) / 2);
        }
        writer.println(a + "\t" + b + "\t" + x + "\t" + y);
      }
    }
    System.out.println("Done!");
    writer.close();
    binIn.close();
  }

  public static void main(String[] args) {
    CLI cli = new CLI(UKBBExtract.class);

    cli.addArg("file", "Binary Intensity File");
    cli.addArgWithDefault("start", "Index of start marker", 0);
    cli.addArgWithDefault("num", "Number of markers to extract", 100);
    cli.addArgWithDefault("out", "Output file name", "./UKBB_xy_0_100_out.xln");

    cli.parseWithExit(args);

    String file = cli.get("file");
    int start = cli.getI("start");
    int num = cli.getI("num");
    String out = cli.has("out") ? cli.get("out")
                                : "./UKBB_xy_" + start + "_" + (start + num) + "_out.xln";

    try {
      extract(file, start, num, out);
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

}
