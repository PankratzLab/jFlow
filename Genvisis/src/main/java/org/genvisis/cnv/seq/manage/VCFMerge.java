package org.genvisis.cnv.seq.manage;

import org.genvisis.cnv.Resources;
import org.genvisis.seq.ReferenceGenome;
import org.genvisis.seq.analysis.GATK;
import org.genvisis.seq.manage.VCFOps;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.GenomeBuild;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ext;

/**
 * @author lane0212 Merge a list of vcfs using GATK;
 */
public class VCFMerge {

  public static void merge(String[] vcfs, String mergeOut, String gatkLoc,
                           String referenceGenomeFasta, int numthreads, int memoryInMB,
                           Logger log) {
    GATK gatk = new GATK(gatkLoc, referenceGenomeFasta, null, null, memoryInMB, true, false, log);
    System.out.println("HID1");

    if (vcfs != null) {
      System.out.println("HID2");

      if (mergeOut != null) {
        System.out.println("HID3");

        if (!Files.exists(mergeOut)) {
          System.out.println("HID4");

          if (Files.exists("", vcfs)) {
            if (VCFOps.verifyIndices(vcfs, log)) {
              log.reportTimeInfo("Beginning to merge " + vcfs.length + " vcfs");
              boolean merged = gatk.mergeVCFs(vcfs, mergeOut, numthreads, false, log);
              if (merged) {
                log.reportTimeInfo("Merged " + ArrayUtils.toStr(vcfs, "\n") + " to " + mergeOut);
                VCFOps.extractSamps(mergeOut, log);
              } else {
                log.reportError("Could not merge " + ArrayUtils.toStr(vcfs, "\n") + " to "
                                + mergeOut);
              }
            } else {
              log.reportError("Could not verify all index files ");
            }
          } else {
            log.reportFilesNotFound(vcfs);
          }
        } else {
          log.reportFileExists(mergeOut);
        }
      } else {
        log.reportError("Must provide an output file for the merged results");
      }
    } else {
      log.reportFilesNotFound(vcfs);
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String[] vcfs = null;
    String mergeOut = null;
    String referenceGenome = null;
    GenomeBuild genomeBuild = ReferenceGenome.DEFAULT_BUILD;
    String gatk = GATK.DEFAULT_GATK;
    int numthreads = 1;
    String usage = "\n" + "seq.manage.VCFMerge requires 0-1 arguments\n";
    usage += "   (1) comma delimited full paths to vcf files (i.e. vcfs= (no default))\n" + "";
    usage += "   (2) full path to the merged output (i.e. mergeOut= (no default))\n" + "";
    usage += "   (3) full path to the gatk directory (i.e. gatk=" + gatk + " (default))\n" + "";
    usage += "   (4) full path to the reference genome  (i.e. ref=refGenome.fasta (default is for "
             + genomeBuild + "))\n" + "";
    usage += PSF.Ext.getNumThreadsCommand(5, numthreads);
    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("vcfs=")) {
        vcfs = arg.split("=")[1].split(",");
        numArgs--;
      } else if (arg.startsWith("mergeOut=")) {
        mergeOut = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("gatk=")) {
        gatk = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("ref=")) {
        referenceGenome = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
        numthreads = ext.parseIntArg(arg);
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
      Logger log = new Logger(vcfs == null ? null : ext.addToRoot(vcfs[0], "mergeLog.log"));
      if (referenceGenome == null) {
        referenceGenome = Resources.genome(genomeBuild, log).getFASTA().get();
      }
      merge(vcfs, mergeOut, gatk, referenceGenome, numthreads, PSF.Ext.DEFAULT_MEMORY_MB, log);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
