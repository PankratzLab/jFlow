package org.genvisis.gwas;

import java.io.File;
import java.util.Date;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;

public class Qc {

  public static final String QC_DIR = "quality_control/";

  public static final String MARKER_QC_DIR = QC_DIR + "marker_qc/";
  public static final String SAMPLE_QC_DIR = QC_DIR + "sample_qc/";
  public static final String LD_PRUNING_DIR = QC_DIR + "ld_pruning/";
  public static final String GENOME_DIR = QC_DIR + "genome/";
  public static final String ANCESTRY_DIR = QC_DIR + "ancestry/";

  /** A rough listing of the Folders created by fullGamut */
  public static String[] FOLDERS_CREATED = {MARKER_QC_DIR, SAMPLE_QC_DIR, LD_PRUNING_DIR,
                                            GENOME_DIR, ANCESTRY_DIR};
  /** A rough listing of the files created, by folder, by fullGamut */
  public static String[][] FILES_CREATED =
                                         {{"plink.bed", "freq.frq", "missing.imiss",
                                           /* "test.missing.missing", *//* not actually necessary */ "hardy.hwe",
                                           "mishap.missing.hap", "gender.assoc", "gender.missing",
                                           "miss_drops.dat"},
                                          {"plink.bed", "missing.imiss"},
                                          {"plink.bed", "plink.prune.in"},
                                          {"plink.bed", "plink.genome", "plink.genome_keep.dat"},
                                          {"plink.bed", "unrelateds.txt"}};

  public static void exportFullGamut(String dir, String plinkPrefix,
                                     boolean keepGenomeInfoForRelatedsOnly) {
    dir = ext.verifyDirFormat(dir) + "quality_control/";
    if (!dir.startsWith("/") && !dir.contains(":")) {
      dir = (new File("./" + dir)).getAbsolutePath();
    }

    String plink = plinkPrefix == null ? "plink" : plinkPrefix;

    StringBuilder cmds = new StringBuilder();
    cmds.append("cd ").append(dir).append("\n");
    cmds.append("mkdir marker_qc/").append("\n");
    cmds.append("cd marker_qc/").append(";\n");
    cmds.append("if [ ! -f ").append(plink).append(".bed ] ; then").append("\n");
    cmds.append("plink2 --bfile ../../").append(plink).append(" --make-bed --noweb --out ./")
        .append(plink).append(";\n");
    cmds.append("fi;\n");
    cmds.append("if [ ! -f ").append(plink).append("_geno20.bed ] ; then").append("\n");
    cmds.append("plink2 --bfile ").append(plink).append(" --geno 0.2 --make-bed --noweb --out ./")
        .append(plink).append("_geno20;\n");
    cmds.append("fi;\n");

    cmds.append("plink2 --bfile ").append(plink)
        .append("_geno20 --mind 0.1 --make-bed --noweb --out ").append(plink).append(";\n");

    cmds.append("if [ ! -f freq.frq ] ; then").append("\n");
    cmds.append("\t").append("plink2 --bfile ").append(plink)
        .append(" --maf 0 --geno 1 --mind 1 --freq --out freq --noweb").append(";\n");
    cmds.append("fi;\n");
    cmds.append("if [ ! -f missing.imiss ] ; then").append("\n");
    cmds.append("\t").append("plink2 --bfile ").append(plink)
        .append(" --maf 0 --geno 1 --mind 1 --missing --out missing --noweb").append(";\n");
    cmds.append("fi;\n");
    cmds.append("if [ ! -f test.missing.missing ] ; then").append("\n");
    cmds.append("\t").append("plink2 --bfile ").append(plink)
        .append(" --maf 0 --geno 1 --mind 1 --test-missing --out test.missing --noweb")
        .append(";\n");
    cmds.append("fi;\n");
    cmds.append("if [ ! -f hardy.hwe ] ; then").append("\n");
    cmds.append("\t").append("plink2 --bfile ").append(plink)
        .append(" --maf 0 --geno 1 --mind 1 --hardy --out hardy --noweb").append(";\n");
    cmds.append("fi;\n");
    cmds.append("if [ ! -f mishap.missing.hap ] ; then").append("\n");
    cmds.append("\t").append("plink2 --bfile ").append(plink)
        .append(" --maf 0 --geno 1 --mind 1 --test-mishap --out mishap --noweb").append(";\n");
    cmds.append("fi;\n");
    cmds.append("if [ ! -f gender.assoc ] ; then").append("\n");
    cmds.append("\t").append("plink2 --bfile ").append(plink)
        .append(" --maf 0 --geno 1 --mind 1 --pheno plink.fam --mpheno 3 --assoc --out gender --noweb")
        .append(";\n");
    cmds.append("fi;\n");
    cmds.append("if [ ! -f gender.missing ] ; then").append("\n");
    cmds.append("\t").append("plink2 --bfile ").append(plink)
        .append(" --maf 0 --geno 1 --mind 1 --pheno plink.fam --mpheno 3 --test-missing --out gender --noweb")
        .append(";\n");
    cmds.append("fi;\n");

    // TODO
    // if (!Files.exists(dir+"marker_qc/miss_drops.dat")) {
    // new File(dir+"marker_qc/miss.crf").delete();
    // int runCode1 = MarkerQC.parseParameters(dir+"marker_qc/miss.crf", log, false); // first
    // generates a file with the defaults
    // if (runCode1 != 0) {
    // // ERROR! TODO not sure if we should quit here; for now, continue;
    // }
    // Files.writeList(Array.addStrToArray("dir="+dir+"marker_qc/",
    // HashVec.loadFileToStringArray(dir+"marker_qc/miss.crf", false, new int[] {0}, false), 1),
    // dir+"marker_qc/miss.crf");
    // int runCode2 = MarkerQC.parseParameters(dir+"marker_qc/miss.crf", log, false); // second runs
    // the filtering
    // if (runCode2 != 0) {
    // // ERROR TODO not sure if we should quit here; for now, continue;
    // }
    // }

    cmds.append("cd ..\n");
    cmds.append("mkdir sample_qc/\n");
    cmds.append("cd sample_qc/\n");
    cmds.append("if [ ! -f ").append(plink).append(".bed ] ; then").append("\n");
    cmds.append("\t").append("plink2 --bfile ../marker_qc/").append(plink)
        .append(" --exclude ../marker_qc/miss_drops.dat --make-bed --noweb\n");
    cmds.append("fi;\n");
    cmds.append("if [ ! -f missing.imiss ] ; then").append("\n");
    cmds.append("\t").append("plink2 --bfile ").append(plink)
        .append(" --maf 0 --geno 1 --mind 1 --missing --out missing --noweb").append(";\n");
    cmds.append("fi;\n");

    cmds.append("cd ..\n");
    cmds.append("mkdir ld_pruning/\n");
    cmds.append("cd ld_pruning/\n");
    cmds.append("if [ ! -f ").append(plink).append(".bed ] ; then").append("\n");
    cmds.append("\t").append("plink2 --bfile ../sample_qc/").append(plink)
        .append(" --mind 0.05 --make-bed --noweb\n");
    cmds.append("fi;\n");
    cmds.append("if [ ! -f ").append(plink).append(".prune.in ] ; then").append("\n");
    cmds.append("\t").append("plink2 --bfile ../sample_qc/").append(plink)
        .append(" --indep-pairwise 50 5 0.3 --noweb\n");
    cmds.append("fi;\n");

    cmds.append("cd ..\n");
    cmds.append("mkdir genome/\n");
    cmds.append("cd genome/\n");
    cmds.append("if [ ! -f ").append(plink).append(".bed ] ; then").append("\n");
    cmds.append("\t").append("plink2 --bfile ../ld_pruning/").append(plink)
        .append(" --extract ../ld_pruning/").append(plink).append(".prune.in --make-bed --noweb\n");
    cmds.append("fi;\n");
    cmds.append("if [ ! -f ").append(plink).append(".genome ] ; then").append("\n");
    cmds.append("\t").append("plink2 --noweb --bfile ").append(plink)
        .append(" --genome" + (keepGenomeInfoForRelatedsOnly ? " --min 0.1" : "")).append("\n");
    cmds.append("fi;\n");
    if (!keepGenomeInfoForRelatedsOnly) {
      cmds.append("if [ ! -f mds20.mds ] ; then").append("\n");
      cmds.append("\tplink --bfile " + plink + " --read-genome " + plink
                  + ".genome --cluster --mds-plot 20 --out mds20 --noweb\n");
      cmds.append("fi;\n");
    }
    if (!Files.exists(dir + "genome/" + plink + ".genome_keep.dat")) {
      String geno = dir + "genome/" + plink + ".genome";
      String fam = dir + "genome/" + plink + ".fam";
      String imiss = dir + "marker_qc/missing.imiss";
      String lrrsd = dir + "genome/lrr_sd.xln";
      if (!Files.exists(lrrsd)) {
        lrrsd = dir + "../../lrr_sd.xln";
      }
      int level = 4;
      // Plink.flagRelateds(geno, fam, imiss, lrrsd, Plink.FLAGS, Plink.THRESHOLDS, level, false);
      cmds.append(Files.getRunString() + " gwas.Plink relate=").append(geno).append(" fam=")
          .append(fam).append(" imiss=").append(imiss).append(" lrr_sd=").append(lrrsd)
          .append(" level=").append(level).append("\n");
    }

    // TODO fill in when ancestry method is finished, currently does nothing but create more files
    // new File(dir+"ancestry/").mkdirs();
    // if (!Files.exists(dir+"ancestry/unrelateds.txt")) {
    // Files.copyFile(dir+"genome/" + plink + ".genome_keep.dat", dir+"ancestry/unrelateds.txt");
    // }
    // if (!Files.exists(dir+"ancestry/" + plink + ".bed")) {
    // CmdLine.runDefaults("plink2 --bfile ../genome/" + plink + " --make-bed --noweb",
    // dir+"ancestry/", log);
    // }
    //
    // ancestry(dir+"ancestry/");

  }

  // TODO CAUTION, MAKE SURE ALL CHANGES TO fullGamut ARE ALSO CHANGED IN exportFullGamut
  public static void fullGamut(String dir, String plinkPrefix,
                               boolean keepGenomeInfoForRelatedsOnly, Logger log) {
    long time;

    time = new Date().getTime();

    dir = ext.verifyDirFormat(dir) + "quality_control/";
    if (!dir.startsWith("/") && !dir.contains(":")) {
      dir = ext.verifyDirFormat((new File("./" + dir)).getAbsolutePath());
    }

    String plink = plinkPrefix == null ? "plink" : plinkPrefix;

    new File(dir + "marker_qc/").mkdirs();
    if (!Files.exists(dir + "marker_qc/" + plink + ".bed")) {
      log.report(ext.getTime() + "]\tRunning initial --make-bed");
      CmdLine.runDefaults("plink2 --bfile ../../" + plink + " --make-bed --noweb --out ./" + plink,
                          dir + "marker_qc/", log);
    }
    if (!Files.exists(dir + "marker_qc/" + plink + ".bed")) {
      System.err.println("Error - CRITICAL ERROR, creating initial PLINK files failed.");
      return;
    }
    if (!Files.exists(dir + "marker_qc/" + plink + "_geno20.bed")) {
      log.report(ext.getTime() + "]\tRunning --geno 0.2");
      CmdLine.runDefaults("plink2 --bfile " + plink + " --geno 0.2 --make-bed --noweb --out ./"
                          + plink + "_geno20", dir + "marker_qc/", log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(dir + "marker_qc/" + plink + ".bed")) {
      log.report(ext.getTime() + "]\tRunning --mind 0.1");
      CmdLine.runDefaults("plink2 --bfile " + plink
                          + "_geno20 --mind 0.1 --make-bed --noweb --out ./" + plink,
                          dir + "marker_qc/", log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(dir + "marker_qc/freq.frq")) {
      log.report(ext.getTime() + "]\tRunning --freq");
      CmdLine.runDefaults("plink2 --bfile " + plink
                          + " --geno 1 --mind 1 --freq --out freq --noweb", dir + "marker_qc/",
                          log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(dir + "marker_qc/missing.imiss")) {
      log.report(ext.getTime() + "]\tRunning --missing");
      CmdLine.runDefaults("plink2 --bfile " + plink
                          + " --geno 1 --mind 1 --missing --out missing --noweb",
                          dir + "marker_qc/", log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(dir + "marker_qc/test.missing.missing")) {
      log.report(ext.getTime() + "]\tRunning --test-missing");
      CmdLine.runDefaults("plink2 --bfile " + plink
                          + " --geno 1 --mind 1 --test-missing --out test.missing --noweb",
                          dir + "marker_qc/", log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(dir + "marker_qc/hardy.hwe")) {
      log.report(ext.getTime() + "]\tRunning --hardy");
      CmdLine.runDefaults("plink2 --bfile " + plink
                          + " --geno 1 --mind 1 --hardy --out hardy --noweb", dir + "marker_qc/",
                          log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(dir + "marker_qc/mishap.missing.hap")) {
      log.report(ext.getTime() + "]\tRunning --test-mishap");
      CmdLine.runDefaults("plink2 --bfile " + plink
                          + " --geno 1 --mind 1 --test-mishap --out mishap --noweb",
                          dir + "marker_qc/", log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(dir + "marker_qc/gender.assoc")) {
      log.report(ext.getTime() + "]\tRunning --assoc gender");
      CmdLine.runDefaults("plink2 --bfile " + plink + " --geno 1 --mind 1 --pheno " + plink
                          + ".fam --mpheno 3 --assoc --out gender --noweb", dir + "marker_qc/",
                          log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(dir + "marker_qc/gender.missing")) {
      log.report(ext.getTime() + "]\tRunning --test-missing gender");
      CmdLine.runDefaults("plink2 --bfile " + plink + " --geno 1 --mind 1 --pheno " + plink
                          + ".fam --mpheno 3 --test-missing --out gender --noweb",
                          dir + "marker_qc/", log);
    }
    PSF.checkInterrupted();

    if (!Files.exists(dir + "marker_qc/miss_drops.dat")) {
      new File(dir + "marker_qc/miss.crf").delete();
      int runCode1 = MarkerQC.parseParameters(dir + "marker_qc/miss.crf", log, false); // first
                                                                                       // generates
                                                                                       // a file
                                                                                       // with the
                                                                                       // defaults
      if (runCode1 != 0) {
        // ERROR! TODO not sure if we should quit here; for now, continue;
      }
      Files.writeList(Array.addStrToArray("dir=" + dir + "marker_qc/",
                                          HashVec.loadFileToStringArray(dir + "marker_qc/miss.crf",
                                                                        false, new int[] {0},
                                                                        false),
                                          1),
                      dir + "marker_qc/miss.crf");
      int runCode2 = MarkerQC.parseParameters(dir + "marker_qc/miss.crf", log, false); // second
                                                                                       // runs the
                                                                                       // filtering
      if (runCode2 != 0) {
        // ERROR TODO not sure if we should quit here; for now, continue;
      }
    }
    PSF.checkInterrupted();

    new File(dir + "sample_qc/").mkdirs();
    if (!Files.exists(dir + "sample_qc/" + plink + ".bed")) {
      log.report(ext.getTime() + "]\tRunning --exclude miss_drops.dat");
      CmdLine.runDefaults("plink2 --bfile ../marker_qc/" + plink
                          + " --exclude ../marker_qc/miss_drops.dat --make-bed --noweb --out "
                          + plink, dir + "sample_qc/", log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(dir + "sample_qc/missing.imiss")) {
      log.report(ext.getTime() + "]\tRunning --missing");
      CmdLine.runDefaults("plink2 --bfile " + plink
                          + " --geno 1 --mind 1 --missing --out missing --noweb",
                          dir + "sample_qc/", log);
    }
    PSF.checkInterrupted();

    new File(dir + "ld_pruning/").mkdirs();
    if (!Files.exists(dir + "ld_pruning/" + plink + ".bed")) {
      log.report(ext.getTime()
                 + "]\tRunning --mind 0.05 (removes samples with callrate <95% for the markers that did pass QC)");
      CmdLine.runDefaults("plink2 --bfile ../sample_qc/" + plink
                          + " --mind 0.05 --make-bed --noweb --out " + plink, dir + "ld_pruning/",
                          log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(dir + "ld_pruning/" + plink + ".prune.in")) {
      log.report(ext.getTime() + "]\tRunning --indep-pairwise 50 5 0.3");
      CmdLine.runDefaults("plink2 --noweb --bfile " + plink + " --indep-pairwise 50 5 0.3 --out "
                          + plink, dir + "ld_pruning/", log);
    }
    PSF.checkInterrupted();

    new File(dir + "genome/").mkdirs();
    if (!Files.exists(dir + "genome/" + plink + ".bed")) {
      log.report(ext.getTime() + "]\tRunning --extract " + plink + ".prune.in");
      CmdLine.runDefaults("plink2 --bfile ../ld_pruning/" + plink + " --extract ../ld_pruning/"
                          + plink + ".prune.in --make-bed --noweb --out " + plink, dir + "genome/",
                          log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(dir + "genome/" + plink + ".genome")) {
      log.report(ext.getTime() + "]\tRunning --genome"
                 + (keepGenomeInfoForRelatedsOnly ? " --min 0.1" : ""));
      CmdLine.runDefaults("plink2 --noweb --bfile " + plink + " --genome"
                          + (keepGenomeInfoForRelatedsOnly ? " --min 0.1" : "") + " --out " + plink,
                          dir + "genome/", log);
    }
    PSF.checkInterrupted();
    if (!keepGenomeInfoForRelatedsOnly && !Files.exists(dir + "genome/mds20.mds")) {
      log.report(ext.getTime() + "]\tRunning --mds-plot 20");
      CmdLine.runDefaults("plink2 --bfile " + plink + " --read-genome " + plink
                          + ".genome --cluster --mds-plot 20 --out mds20 --noweb", dir + "genome/",
                          log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(dir + "genome/" + plink + ".genome_keep.dat")) {
      log.report(ext.getTime() + "]\tRunning flagRelateds");
      String lrrFile = dir + "genome/lrr_sd.xln";
      if (!Files.exists(lrrFile)) {
        lrrFile = dir + "../../lrr_sd.xln";
      }
      Plink.flagRelateds(dir + "genome/" + plink + ".genome", dir + "genome/" + plink + ".fam",
                         dir + "marker_qc/missing.imiss", lrrFile, Plink.FLAGS, Plink.THRESHOLDS, 4,
                         false);
    }
    PSF.checkInterrupted();

    new File(dir + "ancestry/").mkdirs();
    if (!Files.exists(dir + "ancestry/unrelateds.txt")) {
      log.report(ext.getTime() + "]\tCopying genome/" + plink
                 + ".genome_keep.dat to ancestry/unrelateds.txt");
      Files.copyFile(dir + "genome/" + plink + ".genome_keep.dat", dir + "ancestry/unrelateds.txt");
    }
    PSF.checkInterrupted();
    if (!Files.exists(dir + "ancestry/" + plink + ".bed")) {
      log.report(ext.getTime() + "]\tRunning --extract " + plink
                 + ".prune.in (again, this time to ancestry/)");
      CmdLine.runDefaults("plink2 --bfile ../genome/" + plink + " --make-bed --noweb --out "
                          + plink, dir + "ancestry/", log);
    }
    PSF.checkInterrupted();


    System.out.println("Finished this round in " + ext.getTimeElapsed(time));
  }

  public static void fromParameters(String filename, Logger log) {
    Vector<String> params;

    params = Files.parseControlFile(filename, "gwas.Qc",
                                    new String[] {"dir=./",
                                                  "# Make keepGenomeInfoForRelatedsOnly=false if the sample size is small and you want to run MDS plot as well",
                                                  "keepGenomeInfoForRelatedsOnly=true"},
                                    log);

    if (params != null) {
      params.add("log=" + log.getFilename());
      main(Array.toStringArray(params));
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String dir = "./";
    boolean keepGenomeInfoForRelatedsOnly = true;
    String logfile = null;
    Logger log;

    String usage = "\n" + "gwas.Qc requires 0-1 arguments\n"
                   + "   (1) directory with plink.* files (i.e. dir=" + dir + " (default))\n"
                   + "   (2) if no MDS will be run, smaller file (i.e. keepGenomeInfoForRelatedsOnly="
                   + keepGenomeInfoForRelatedsOnly + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = ext.parseStringArg(arg, "./");
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("keepGenomeInfoForRelatedsOnly=")) {
        keepGenomeInfoForRelatedsOnly = ext.parseBooleanArg(arg);
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    if (logfile == null) {
      log = new Logger(dir + "fullGamutOfMarkerAndSampleQC.log");
    } else {
      log = new Logger(logfile);
    }
    try {
      fullGamut(dir, null, keepGenomeInfoForRelatedsOnly, log);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
