package org.genvisis.one.JL;

import java.io.PrintWriter;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.shared.filesys.Positions;
import org.pankratzlab.shared.filesys.Segment;

public class SegmentContigSubset {

  public static void subset(String segFile, String[] contigs, Logger log) {
    Segment[] q = Segment.loadRegions(segFile, 0, 1, 2, 0, true, true, true, 100);
    String output = ext.addToRoot(segFile, "_" + ArrayUtils.toStr(contigs, "_"));
    try {
      PrintWriter writer = Files.openAppropriateWriter(output);
      for (Segment element : q) {
        String contig = Positions.getChromosomeUCSC(element.getChr(), true);
        if (ext.indexOfStr(contig, contigs) >= 0) {
          writer.println(contig + "\t" + element.getStart() + "\t" + element.getStop());
        }
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + output);
      log.reportException(e);
    }

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String segFile = "Segments.txt";

    // String segFile = "SegmentSubSet.dat";
    String[] contigs = new String[] {"chr1"};

    String usage = "\n" + "one.JL.SegmentSubSet requires 0-1 arguments\n";
    usage += "   (1) a file of segments to subset (i.e. segs=" + segFile + " (default))\n" + "";
    usage += "   (2) a comma delimited list of contigs (i.e. contigs="
             + ArrayUtils.toStr(contigs, ",") + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("segs=")) {
        segFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("contigs=")) {
        contigs = arg.split("=")[1].split(",");
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      subset(segFile, contigs, new Logger());
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
