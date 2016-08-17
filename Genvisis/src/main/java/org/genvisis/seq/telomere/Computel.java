package org.genvisis.seq.telomere;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Map;
import java.util.Scanner;

import org.apache.commons.cli.Options;
import org.genvisis.CLI;
import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Zip;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.BamOps;

/**
 * Wrapper for the computel telomere length estimator ala
 * https://github.com/lilit-nersisyan/computel
 *
 *
 * Computel is a little annoying to set up, so this will reflect that...
 *
 *
 * NOTE: I had to re-build the samtools version to operate using MSI's libs for C
 *
 * Also, Currently looks like computel takes a perfectly nice ~10 gigabyte bam file and converts it
 * into a massive 250GB monster
 */
public class Computel {

  private static final String SAM_TO_FASTQ_LOC = "SamToFastq.jar";

  private static boolean convertToFasta(String inputBam, String computelLoc, String r1, String r2,
                                        Logger log) {
    String[] inputs = new String[] {inputBam};
    String[] outputs = new String[] {r1, r2};
    ArrayList<String> command = new ArrayList<String>();
    command.add("java");
    command.add("-jar");
    command.add(computelLoc + SAM_TO_FASTQ_LOC);
    command.add("I=");
    command.add(inputBam);
    command.add("F=");
    command.add(r1);
    command.add("F2=");
    command.add(r2);

    return CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", inputs, outputs, true,
                                            false, false, log);
  }

  private static void copyDirectory(File sourceLocation, File targetLocation) throws IOException {

    if (sourceLocation.isDirectory()) {
      if (!targetLocation.exists()) {
        targetLocation.mkdir();
      }

      String[] children = sourceLocation.list();
      for (String element : children) {
        copyDirectory(new File(sourceLocation, element), new File(targetLocation, element));
      }
    } else {

      InputStream in = new FileInputStream(sourceLocation);
      OutputStream out = new FileOutputStream(targetLocation);

      // Copy the bits from instream to outstream
      byte[] buf = new byte[1024];
      int len;
      while ((len = in.read(buf)) > 0) {
        out.write(buf, 0, len);
      }
      in.close();
      out.close();
    }
  }

  public static void main(String[] args) {
    String targetBam = "bam.bam";
    String computelLocation = "computel/";
    String outDir = "out/";

    Options options = CLI.defaultOptions();
    final String bam = "bam";
    CLI.addArg(options, bam, "bam file to analyze", targetBam);

    final String computel = "computel";
    CLI.addArg(options, computel, "full computel directory (as git clone ideally)",
               computelLocation);

    final String outdir = "out";
    CLI.addArg(options, outdir, "the output directory for results", outDir);

    Map<String, String> parsed = CLI.parseWithExit(Computel.class, options, args);

    test(parsed.get(bam), parsed.get(computel), parsed.get(outdir));
  }

  private static String manageConfig(Logger log, String finalOutDirectory, String bowtieSamDir,
                                     String r1, String r2,
                                     int readLength) throws FileNotFoundException {
    Scanner s = new Scanner(new File(finalOutDirectory + "src/examples/config_unix.txt"));

    String config = s.useDelimiter("\\Z").next();
    s.close();
    config = processConfig(finalOutDirectory, bowtieSamDir, config, r1, r2, readLength);
    log.report(config);
    return config;
  }

  /**
   * @param computelOperatingDir Where the config will be setup relative to
   */
  private static String processConfig(String computelOperatingDir, String bowtieSamDir,
                                      String config, String r1, String r2, int readLength) {
    config = config.replaceAll("scripts.dir	./scripts",
                               "scripts.dir	" + computelOperatingDir + "src/scripts");
    String btieBuild = bowtieSamDir + "bowtie2-2.1.0-linux/bowtie2-build";
    Files.chmod(btieBuild);

    config = config.replaceAll("bowtie.build.path	./bowtie2-2.1.0-linux/bowtie2-build",
                               "bowtie.build.path	" + btieBuild);
    String btieAln = bowtieSamDir + "bowtie2-2.1.0-linux/bowtie2-align";
    Files.chmod(btieAln);

    config = config.replaceAll("bowtie.align.path	./bowtie2-2.1.0-linux/bowtie2-align",
                               "bowtie.align.path	" + btieAln);
    String samTools = bowtieSamDir + "samtools-0.1.19-linux/samtools";
    Files.chmod(samTools);
    config = config.replaceAll("samtools.path	./samtools-0.1.19-linux/samtools",
                               "samtools.path	" + samTools);

    String samToFastQ = computelOperatingDir + SAM_TO_FASTQ_LOC;
    Files.chmod(samToFastQ);
    config = config.replaceAll("picard.samtofastq.jar	./SamToFastq.jar",
                               "picard.samtofastq.jar	" + samToFastQ);

    config = config.replaceAll("read.length	76", "read.length	" + readLength);

    config = config.replaceAll("files.with.prefix	T", "files.with.prefix	F");

    config =
        config.replaceAll("fastq	examples/tel_reads.fq, examples/tel_reads1.fq, examples/tel_reads2.fq",
                          "fastq	" + r1 + "," + r2);

    config = config.replaceAll("compute.base.cov	F", "compute.base.cov	T");

    config =
        config.replaceAll("base.index.pathtoprefix	./examples/base.index/base_index",
                          "base.index.pathtoprefix	" + computelOperatingDir + "src/examples/base.index/base_index");

    config = config.replaceAll("output.dir	output",
                               "output.dir	" + computelOperatingDir + "results");
    // compute.base.cov F
    return config;
  }

  private static boolean runComputel(String config, String computelCommandR, Logger log) {
    String[] inputs = new String[] {config, computelCommandR};
    ext.parseDirectoryOfFile(config);
    // TODO, cleanup computels mess, and determine appropriate outputs
    String[] outputs = null;
    // String[] outputs = execDir + "";
    // new String[] { r1, r2 };
    ArrayList<String> command = new ArrayList<String>();
    command.add("Rscript");
    command.add(computelCommandR);
    command.add(config);

    return CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", inputs, outputs, true,
                                            false, false, log);
  }

  private static void runComputel(String inputBam, String outputDir, String computelDirectory,
                                  Logger log) {
    String finalOutDirectory = outputDir + ext.rootOf(inputBam) + "/";
    new File(finalOutDirectory).mkdirs();

    try {
      copyDirectory(new File(computelDirectory), new File(finalOutDirectory));
    } catch (IOException e) {
      log.reportException(e);
    }
    String bowtieSamDir = finalOutDirectory + "bowtie2_samtools_binaries_for_linux/";
    if (!Files.exists(bowtieSamDir)) {
      Zip.unzipFile(finalOutDirectory + "bowtie2_samtools_binaries_for_linux.zip", bowtieSamDir);
    }
    try {

      String r1 = finalOutDirectory + "src/examples/analysis_reads1.fq";
      String r2 = finalOutDirectory + "src/examples/analysis_reads2.fq";
      boolean converted = convertToFasta(inputBam, finalOutDirectory, r1, r2, log);
      if (converted) {
        int readLength = BamOps.estimateReadSize(inputBam, 100000, log);

        String config = manageConfig(log, finalOutDirectory, bowtieSamDir, r1, r2, readLength);
        String configFile = finalOutDirectory + "computelConfig.txt";
        String computelCommand = finalOutDirectory + "src/scripts/computel.cmd.R";
        Files.write(config, configFile);
        runComputel(configFile, computelCommand, log);
      }
    } catch (FileNotFoundException e) {
      log.reportException(e);
    }

  }

  public static void test(String bam, String computelLocation, String outDir) {
    runComputel(bam, outDir, computelLocation, new Logger(outDir + "computel.Log"));
  }

  private Computel() {

  }
}
