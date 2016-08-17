package org.genvisis.one.JL;

import java.io.FileNotFoundException;

import org.genvisis.common.Files;

public class prepARICSRA {

  public static void main(String[] args) {
    int numArgs = args.length;

    String sraRunTable = "/Volumes/Work/data/aric_sra/prep/SraRunTable.txt";

    String usage = "\n" + "this requires 0-1 arguments\n"
        + "   (1) SRA data table (i.e. sraRunTable=" + sraRunTable + " (default))\n";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("sraRunTable=")) {
        sraRunTable = arg.split("=")[1];
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
      prep(sraRunTable);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  private static void prep(String sraRunTable) {
    // String[] extract = new String[]{}
    try {
      Files.getAppropriateReader(sraRunTable);

    } catch (FileNotFoundException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

  }

}
