package org.genvisis.cnv.manage;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;
import org.genvisis.filesys.PlainTextExport;

public class DumpSerializedObject {
  private static void dump(String projectPropertyFile, String filename, String logFile) {
    Project proj = new Project(projectPropertyFile, logFile, false);
    Object object;
    object = SerializedFiles.readSerial(filename, false, proj.getLog(), false, false);

    if (object instanceof TextExport) {
      ((TextExport) object).exportToText(proj, ext.parseDirectoryOfFile(filename)
                                               + ext.rootOf(filename) + "_dump.xln");
      return;
    } else if (object instanceof PlainTextExport) {
      ((PlainTextExport) object).exportToText(ext.parseDirectoryOfFile(filename)
                                              + ext.rootOf(filename) + "_dump.xln", proj.getLog());
    }

    proj.getLog()
        .report("Information on class:" + "\n" + "object.getClass().getName()="
                + object.getClass().getName() + "\n" + "object.getClass()=" + object.getClass()
                + "\n" + "object.toString()=" + object.toString());
    /**
     * Classes we should support now:
     *
     * outliers.ser --* AnnotationCollection --* Centroids --* ClusterFilterCollection --* MarkerSet
     * CNVariant --- would only write one cnvariant to one file... --* SNPMarkerSet --* SegmentList
     * --* SegmentLists
     *
     */

    /**
     * Classes we might consider supporting in the future:
     *
     * BurdenMatrix ChromatinAccessibility GenotypeMatrix LDdatabase SerialStringMatrix SuperNovo
     */

    // if (object instanceof AnnotationCollection) {
    // proj.getLog().report("Detected an AnnotationCollection file");
    //// AnnotationCollection.dump();
    // }


    if (filename.endsWith("outliers.ser")) { // Hashtable<String, Float>
      Sample.dumpOutOfRangeValues(filename, ext.parseDirectoryOfFile(filename)
                                            + ext.rootOf(filename) + "_dump.xln",
                                  false);
    }


  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String proj = null;
    String filename = null;
    String logfile = null;

    String usage = "\n" + "widgets.DumpSerializedObject requires 0-1 arguments\n"
                   + "   (1) filename (e.g. file=outliers.ser (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        proj = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = arg.split("=")[1];
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
      if (proj != null) {
        dump(proj, filename, logfile);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
