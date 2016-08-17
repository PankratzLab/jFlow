package org.genvisis.one;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Hashtable;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

/**
 * Currently specific to spector data
 *
 */
public class PrepSuperNovo {
  private static String[] TRIO_ENDINGS = new String[] {"C", "D", "M"};

  public static void main(String[] args) {
    String dir = "/lustre/lanej0/Project_Spector_Project_014/140516_SN261_0548_AC4GVNACXX/sam/";
    String extension = ".sorted.dedup.realigned.recal.bam";
    String outputDir =
        "/lustre/lanej0/Project_Spector_Project_014/140516_SN261_0548_AC4GVNACXX/superNovo/";

    prepDir(dir, extension, outputDir, new Logger());
  }

  public static void prepDir(String dir, String extension, String outputDir, Logger log) {
    String[] bams = Files.list(dir, extension, false);
    Hashtable<String, Hashtable<String, String>> trios =
        new Hashtable<String, Hashtable<String, String>>();
    for (String bam : bams) {
      String id = bam.split("_")[0];
      String trioId = id.substring(0, id.length() - 1);
      char trioMember = id.charAt(id.length() - 1);

      int index = ext.indexOfStr(trioMember + "", TRIO_ENDINGS);
      if (index < 0) {
        System.out.println("ERROR not trio\t" + trioMember + "\t" + id);
        return;
      } else if (!trios.containsKey(trioId)) {
        trios.put(trioId, new Hashtable<String, String>());
      }
      trios.get(trioId).put(trioMember + "", ext.removeDirectoryInfo(bam));

    }
    String trioListFile = outputDir + "trios.trio";
    new File(outputDir).mkdirs();
    try {
      PrintWriter writer = new PrintWriter(new FileWriter(trioListFile));

      for (String trio : trios.keySet()) {
        String trioentry = trio;
        boolean write = true;
        for (String element : TRIO_ENDINGS) {
          if (trios.get(trio).containsKey(element)) {
            trioentry += "\t" + trios.get(trio).get(element);
          } else {
            write = false;
          }
        }
        if (write) {
          writer.println(trioentry);
        } else {
          log.reportTimeError("Skipping un filled trio for " + trioentry);
        }
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + trioListFile);
      log.reportException(e);
    }

  }
}
