package org.genvisis.seq.analysis;

import java.io.File;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;

public class Ancestry {

  public static void main(String[] args) {
    int numArgs = args.length;
    String vcf = "avcf.vcf.gz";
    String hapMapRsIds = null;
    String g1000RsIds = null;
    String vpopFile = null;
    String outputDirectory = null;

    String usage = "\n" + "seq.analysis.ancestry requires 0-1 arguments\n";
    usage += "   (1) vcf (i.e. vcf=" + vcf + " (default))\n" + "";
    usage += "   (2) file of hapMap rs ids to extract (i.e. hapMap=" + vcf + " (default))\n" + "";
    usage +=
        "   (3) file of 1000 genome rs ids to extract (i.e. g1000=" + vcf + " (default))\n" + "";
    usage += "   (4) vpopFile  (i.e. vpop=" + vpopFile + " (default))\n" + "";
    usage += PSF.Ext.getOutputDirCommand(4, "");

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("vcf=")) {
        vcf = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("hapMap=")) {
        hapMapRsIds = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("g1000=")) {
        g1000RsIds = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("vpop=")) {
        vpopFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith(PSF.Ext.OUTPUT_DIR_COMMAND)) {
        outputDirectory = arg.split("=")[1];
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
      runAncestry(vcf, hapMapRsIds, g1000RsIds, vpopFile, outputDirectory);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void runAncestry(String vcf, String hapMapRsIds, String g1000RsIds, String vpopFile,
                                 String outputDirectory) {
    new File(outputDirectory).mkdirs();

    Logger log = new Logger(outputDirectory + "ancestry.log");
    log.reportTimeInfo("Output directory: " + outputDirectory);

    String curVCF = VCFOps.extractIDs(vcf, g1000RsIds, outputDirectory, true, true, null, null,
                                      false, true, log);
    log.reportTimeInfo("Current VCF: " + curVCF);
    curVCF = VCFOps.extractIDs(curVCF, hapMapRsIds, outputDirectory, true, true, null, null, false,
                               true, log);
    log.reportTimeInfo("Current VCF: " + curVCF);

    String[] filenames =
        VcfPopulation.splitVcfByPopulation(curVCF, vpopFile, false, false, false, log);
    for (String filename : filenames) {
      if (filename.endsWith(".USE.vcf.gz")) {
        curVCF = filename;
        log.reportTimeInfo("Current VCF: " + curVCF);

        curVCF = VCFOps.removeFilteredVariants(curVCF, true, true, log);
        log.reportTimeInfo("Current VCF: " + curVCF);

        String copyVcf =
            ext.parseDirectoryOfFile(curVCF) + ext.rootOf(vpopFile) + "ancestry.vcf.gz";
        Files.copyFileUsingFileChannels(new File(curVCF), new File(copyVcf), log);
        Files.copyFileUsingFileChannels(new File(curVCF + ".tbi"), new File(copyVcf + ".tbi"), log);
        curVCF = copyVcf;
        log.reportTimeInfo("Current VCF: " + curVCF);

        VCFOps.vcfGwasQC(curVCF, log);
        return;
      }
    }
  }

}
