package org.genvisis.link;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.common.ext;

public class Allegro {
  public static void makeOptFile(int chr, String model, boolean unique) throws IOException {
    String chrome = ext.chrome(chr);
    PrintWriter writer =
                       new PrintWriter(new FileWriter("useful" + (unique ? chrome : "") + ".opt"));

    writer.println("% Read input in LINKAGE style format:");
    writer.println("PREFILE re_chrom" + chrome + ".pre");
    // writer.println("PREFILE trimd.pedin."+chrome);

    if (model.equals("NPL")) {
      writer.println("DATFILE map" + chrome + ".dat");
    } else if (model.equals("DOM")) {
      new LinkageMap(chr).createDominantMap();
      writer.println("DATFILE map" + chrome + ".D.dat");
    } else if (model.equals("REC")) {
      new LinkageMap(chr).createRecessiveMap();
      writer.println("DATFILE map" + chrome + ".R.dat");
    } else {
      System.err.println("Error: model must be either NPL, DOM, or REC, not " + model);
      System.exit(1);
    }

    writer.println("");
    // writer.println("% Won't run with steps higher than 2 for this
    // data.");
    // writer.println("MAXSTEPLENGTH 1.1");
    // writer.println("");
    writer.println("% Run multipoint analysis on all [equally weighted] pairs and output to  .txt file");

    if (model.equals("NPL")) {
      writer.println("MODEL mpt lin all equal chrom" + chrome + ".lin.out chromf" + chrome
                     + ".lin.out");
      writer.println("MODEL mpt exp all equal chrom" + chrome + ".exp.out chromf" + chrome
                     + ".exp.out");
    } else {
      writer.print("MODEL mpt par ");
      if (chr == 23) {
        writer.print("X ");
      }
      writer.print("het ");
      if (model.equals("DOM")) {
        writer.println("chrom" + chrome + ".d.out chromf" + chrome + ".d.out");
      } else if (model.equals("REC")) {
        writer.println("chrom" + chrome + ".r.out chromf" + chrome + ".r.out");
      } else {
        System.err.println("Error - invalid model");
      }
    }
    writer.println("");
    writer.println("% Other options:");
    writer.println("MAXMEMORY 2000");

    writer.close();
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    int chr = 2;
    String model = "NPL";
    boolean unique = true;

    String usage = "\\n" + "link.Allegro requires 0-1 arguments\n" + "   (1) chromosome (i.e. chr="
                   + chr + " (default))\n" + "   (2) model (i.e. model=" + model
                   + " (default) options include: NPL, DOM, REC)\n"
                   + "   (3) unique opt file (i.e. -unique (" + (unique ? "" : "not the ")
                   + "default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("chr")) {
        chr = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("model=")) {
        model = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-unique")) {
        unique = true;
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      makeOptFile(chr, model, unique);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
