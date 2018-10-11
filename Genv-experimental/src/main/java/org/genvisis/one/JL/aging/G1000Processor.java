package org.genvisis.one.JL.aging;

import java.io.File;
import java.util.ArrayList;
import org.genvisis.seq.SeqVariables.ASSAY_TYPE;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.qsub.Qsub;
import org.pankratzlab.internal.utils.sra.SRAPipeline;

/**
 *
 */
public class G1000Processor {

  private static final String ASCP_COMMAND_BASE = "/home/pankrat2/lane0212/.aspera/connect/bin/ascp -i /home/pankrat2/lane0212/.aspera/connect/etc/asperaweb_id_dsa.openssh -Tr -Q -l 600M -L- fasp-g1k@fasp.1000genomes.ebi.ac.uk:vol1/ftp/phase3/data/";
  private static final String EXOME_FLAG = "exome_alignment";

  private G1000Processor() {

  }

  private static String generateMtDNACNCommand(String localFile, String jarFile, ASSAY_TYPE aType,
                                               String refGenome, String captureBed) {
    ArrayList<String> mtDNACmd = new ArrayList<>();
    String outDir = ext.parseDirectoryOfFile(localFile) + "mtDNACN/";
    new File(outDir).mkdirs();
    String bamList = localFile + ".bamList.txt";
    Files.write(localFile, bamList);
    mtDNACmd.add("java");
    mtDNACmd.add("-jar");
    mtDNACmd.add(jarFile);
    mtDNACmd.add("seq.analysis.MitoSeqCN");
    mtDNACmd.add("bams=" + bamList);
    mtDNACmd.add("assemblyName=" + ASSEMBLY_NAME.GRCH37.toString());
    mtDNACmd.add("refGenome=" + refGenome);
    mtDNACmd.add("assayType=" + aType.toString());
    if (aType == ASSAY_TYPE.WXS) {
      mtDNACmd.add("captureBed=" + captureBed);
    }
    mtDNACmd.add("outDir=" + outDir);
    mtDNACmd.add("threads=" + 1);
    return ArrayUtils.toStr(mtDNACmd, " ");
  }

  private static String generateQsub(String localFile, String jarFile, ASSAY_TYPE aType,
                                     String refGenome, String captureBed) {
    String qsubFile = localFile + ".qsub";
    ArrayList<String> qsubInfo = new ArrayList<>();
    qsubInfo.add(generateMtDNACNCommand(localFile, jarFile, aType, refGenome, captureBed));
    qsubInfo.add("rm " + localFile);

    Qsub.qsub(qsubFile, ArrayUtils.toStr(qsubInfo, "\n"), 10000, 24, 1);
    return qsubFile;

  }

  private static String getSample(String remoteFile) {
    return ext.removeDirectoryInfo(remoteFile).split("\\.")[0];
  }

  private static String getLocalFile(String remoteFile, String outDir) {
    return outDir + getSample(remoteFile) + "/" + ext.removeDirectoryInfo(remoteFile);
  }

  private static String getDLCommand(String remoteFile, String localFile) {
    String cmd = ASCP_COMMAND_BASE + remoteFile + " " + localFile + "\n";
    cmd = cmd + ASCP_COMMAND_BASE + remoteFile + ".bai" + " " + ext.rootOf(localFile, false)
          + ".bai";

    return cmd;
  }

  private static void setup(String fileList, String capture, String outDir, String jarFile,
                            String refGenome) {
    String[] files = HashVec.loadFileToStringArray(fileList, false, new int[] {0}, true);
    new File(outDir).mkdirs();
    Logger log = new Logger(outDir + "1000g.log");
    ArrayList<String> process = new ArrayList<>();
    for (String file : files) {
      if (file.endsWith("bam")) {
        ASSAY_TYPE aType = ASSAY_TYPE.WGS;
        if (file.contains(EXOME_FLAG)) {
          aType = ASSAY_TYPE.WXS;
        }
        String localFile = getLocalFile(file, outDir);
        new File(ext.parseDirectoryOfFile(localFile)).mkdirs();
        process.add("echo \"start dl to " + localFile + "\"");
        process.add(getDLCommand(file, localFile));
        String qsubFile = generateQsub(localFile, jarFile, aType, refGenome, capture);
        process.add("qsub -q pankratz " + qsubFile);

        process.add("echo \"end dl to " + localFile + "\"");
        log.reportTimeInfo("preparing " + file);
      }
    }

    String runFile = outDir + "1000g.run";
    Files.writeIterable(process, runFile);
  }

  /**
   * @param args
   */
  public static void main(String[] args) {

    CLI c = new CLI(SRAPipeline.class);

    String fileList = "/scratch.global/lanej/1000G/fileList.txt";
    c.addArgWithDefault("fileList", "list of files", fileList);
    c.addArgWithDefault("capture", "capture bed file",
                        "/scratch.global/lanej/1000G/capture/20130108.exome.targets.bed");

    c.addArgWithDefault(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, "/scratch.global/lanej/1000G/");
    c.addArgWithDefault("jarFile", "jar to use", "/home/pankrat2/lane0212/genvisis1000G.jar");
    c.addArgWithDefault(CLI.ARG_REFERENCE_GENOME, CLI.DESC_REFERENCE_GENOME,
                        "/home/pankrat2/public/bin/ref/GRCh37_canon.fa");

    c.parseWithExit(args);
    setup(c.get("fileList"), c.get("capture"), c.get(CLI.ARG_OUTDIR), c.get("jarFile"),
          c.get(CLI.ARG_REFERENCE_GENOME));
  }

}
