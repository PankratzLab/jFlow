package org.genvisis.cnv.filesys;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Date;

import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class Dump {
  public static void dumpMdRaf(String filename, int[] indicesOfMarkersToDump, Logger log) {
    try {
      MarkerData[] mkData;

      if (indicesOfMarkersToDump == null) {
        indicesOfMarkersToDump = new int[] {0};
      }
      mkData = TransposeData.loadFromRAF(filename, indicesOfMarkersToDump);
      for (MarkerData element : mkData) {
        element.dump(null,
                     ext.parseDirectoryOfFile(filename) + ext.rootOf(filename) + "_dump_"
                           + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + ".xln",
                     null, false, log);
      }
    } catch (Exception e) {
      log.reportError("Error dumping data from " + filename + " to a textfile");
      log.reportException(e);
    }
  }

  public static void dumpPlinkBim(String filename) {
    PrintWriter writer;

    try {
      MarkerSet set = MarkerSet.load(filename, false);
      String[] markerNames = set.getMarkerNames();
      byte[] chrs = set.getChrs();
      int[] positions = set.getPositions();

      writer =
          new PrintWriter(new FileWriter(ext.parseDirectoryOfFile(filename) + ext.rootOf(filename)
                                         + "_dump_" + "_" + set.getFingerprint()
                                         + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date()))
                                         + ".xln"));
      writer.println("MarkerName\tChr\tPosition");
      for (int i = 0; i < markerNames.length; i++) {
        writer.println(markerNames[i] + "\t" + chrs[i] + "\t" + positions[i]);
      }
      writer.close();

    } catch (Exception e) {
      System.err.println("Error dumping data from " + filename + " to a textfile");
      e.printStackTrace();
    }
  }

  public static void dumpSampRaf(String filename) {
    PrintWriter writer;

    try {
      Sample samp = Sample.loadFromRandomAccessFile(filename, false);
      float[] gcs = samp.getGCs();
      float[] xs = samp.getXs();
      float[] ys = samp.getYs();
      float[] thetas = samp.getThetas();
      float[] rs = samp.getRs();
      float[] lrrs = samp.getLRRs();
      float[] bafs = samp.getBAFs();
      byte[] forwardGenotypes = samp.getForwardGenotypes();
      byte[] abGenotypes = samp.getAB_Genotypes();

      writer =
          new PrintWriter(new FileWriter(ext.parseDirectoryOfFile(filename) + ext.rootOf(filename)
                                         + "_dump_" + samp.getFingerprint() + ".xln"));
      writer.println((xs == null ? "" : "X") + (ys == null ? "" : "\tY")
                     + (thetas == null ? "" : "\tTheta") + (rs == null ? "" : "\tR")
                     + (bafs == null ? "" : "\tBAF") + (lrrs == null ? "" : "\tLRR")
                     + (gcs == null ? "" : "\tGC_score")
                     + (abGenotypes == null ? "" : "\tAB_Genotypes")
                     + (forwardGenotypes == null ? "" : "\tForward_Genotypes"));
      for (int i = 0; i < samp.getDataLength(); i++) {
        writer.println((xs == null ? "" : xs[i]) + (ys == null ? "" : "\t" + ys[i])
                       + (thetas == null ? "" : "\t" + thetas[i]) + (rs == null ? "" : "\t" + rs[i])
                       + (bafs == null ? "" : "\t" + bafs[i]) + (lrrs == null ? "" : "\t" + lrrs[i])
                       + (gcs == null ? "" : "\t" + gcs[i])
                       + (abGenotypes == null ? ""
                                              : "\t" + (abGenotypes[i] == -1 ? "--"
                                                                             : Sample.AB_PAIRS[abGenotypes[i]]))
                       + (forwardGenotypes == null ? ""
                                                   : "\t"
                                                     + Sample.ALLELE_PAIRS[forwardGenotypes[i]]));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error dumping data from " + filename + " to a textfile");
      e.printStackTrace();
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    // String filename = "sample.fsamp";
    // String filename = "D:/GEDI_exome/transposed/markers.1.mdRAF";
    // String filename = "D:/COGA_exome/samples/100_AA8317_8001548631_47321_D01.sampRAF";
    String filename = "D:/COGA_exome/transposed/markers.0.mdRAF";
    // String filename = "D:/COGA_exome/transposed/outliers.ser";
    // String filename = "D:/GEDI_exome/samples/1002900362.sampRAF";
    // int[] indicesOfMarkersToDump = null;
    int[] indicesOfMarkersToDump = new int[] {195235};
    String[] commandTemp;

    String usage = "\n" + "cnv.filesys.Dump requires 0-1 arguments\n"
                   + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("markerindices=")) {
        commandTemp = arg.split("=")[1].split(",");
        indicesOfMarkersToDump = new int[commandTemp.length];
        for (int j = 0; j < indicesOfMarkersToDump.length; j++) {
          indicesOfMarkersToDump[j] = Integer.parseInt(commandTemp[j]);
        }
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    try {
      if (filename.endsWith(Sample.SAMPLE_FILE_EXTENSION)) {
        dumpSampRaf(filename);
      } else if (filename.endsWith(".bim")) {
        dumpPlinkBim(filename);
      } else if (filename.endsWith(MarkerData.MARKER_DATA_FILE_EXTENSION)) {
        dumpMdRaf(filename, indicesOfMarkersToDump, new Logger());
      } else if (filename.endsWith("outliers.ser")) {
        Sample.dumpOutOfRangeValues(filename, ext.parseDirectoryOfFile(filename)
                                              + ext.rootOf(filename) + "_dump.xln",
                                    false);
      } else {
        System.err.println("Error - unrecongized file type: " + filename);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
