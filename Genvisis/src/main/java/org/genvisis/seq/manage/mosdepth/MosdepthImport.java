package org.genvisis.seq.manage.mosdepth;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.cli.ParseException;
import org.genvisis.cnv.Resources;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.CmdLine;
import org.pankratzlab.common.Command;
import org.pankratzlab.common.Elision;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.GenomeBuild;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;

public class MosdepthImport {

  private static final String ARG_BIN_BED = "binBed";
  private static final String DESC_BIN_BED = "Bins-to-use .BED file; determines which bins to import as markers.  This can be generated for the entire genome with "
                                             + FASTAToBedConversion.class.getCanonicalName();
  private static final String ARG_GENO_VCF = "genoVCF";
  private static final String DESC_GENO_VCF = "VCF file with genotype information with the same markers as the selected marker VCF file.";
  private static final String ARG_SELECTED_SNP_VCF = "selectedVCF";
  private static final String DESC_SELECTED_SNP_VCF = "VCF file with selected markers for (at least) the bins in the 'binsBed' file.  These markers should also be present in the file used in the 'genoVCF' file.  This can be generated using "
                                                      + NGSBinSNPSelector.class.getCanonicalName();
  private static final String ARG_MOS_DIR = "mosDir";
  private static final String DESC_MOS_DIR = "Mosdepth results files directory. These should be created by running the mosdepth program using the same .BED file (or a superset file) as is used for the 'binsBed' argument.";
  private static final String ARG_MOS_EXT = "mosExt";
  private static final String DESC_MOS_EXT = "Mosdepth results files extension.";
  private static final String ARG_CRAMCOUNT_DIR = "cramReadsDir";
  private static final String DESC_CRAMCOUNT_DIR = "CRAM read file directory.  These can be generated with "
                                                   + CRAMSnpReader.class.getCanonicalName();
  private static final String ARG_CRAM_DIR = "cramDir";
  private static final String DESC_CRAM_DIR = "Directory with raw CRAM files.";

  private static CLI generateBaseCLI() {
    CLI cli = new CLI(MosdepthPipeline.class);

    cli.addArg("projDir", "Project directory", true);
    cli.addArg("projName", "Project name", true);
    cli.addArg("propDir", "Project property files directory", true);
    cli.addArg("jobID", "Job Identifier");
    cli.addArg(CLI.ARG_THREADS, CLI.DESC_THREADS, false);
    cli.addArg("build", "Genome build, one of " + ArrayUtils.toStr(GenomeBuild.values(), ", "),
               true);

    return cli;
  }

  private static void tryCLI1(String[] args) throws ParseException {
    CLI cli = generateBaseCLI();
    cli.addArg(ARG_BIN_BED, DESC_BIN_BED, true);
    cli.addArg(ARG_SELECTED_SNP_VCF, DESC_SELECTED_SNP_VCF, true);
    cli.addArg(ARG_GENO_VCF, DESC_GENO_VCF, false);
    cli.addArg(ARG_MOS_DIR, DESC_MOS_DIR, true);
    cli.addArg(ARG_MOS_EXT, DESC_MOS_EXT, ".bed", false);
    cli.addArg(ARG_CRAMCOUNT_DIR, DESC_CRAMCOUNT_DIR, true);

    cli.parse(args);

    MosdepthPipeline mi = new MosdepthPipeline();
    mi.setProjectDir(cli.get("projDir"));
    mi.setProjectName(cli.get("projName"));
    mi.setProjectPropertiesDir(cli.get("propDir"));
    mi.setNumThreads(cli.has(CLI.ARG_THREADS) ? cli.getI(CLI.ARG_THREADS)
                                              : Runtime.getRuntime().availableProcessors());

    mi.setBinsToUseBED(cli.get(ARG_BIN_BED));
    if (cli.has(ARG_GENO_VCF)) {
      mi.setGenotypeVCF(cli.get(ARG_GENO_VCF));
    }
    mi.setSelectedMarkerVCF(cli.get(ARG_SELECTED_SNP_VCF));
    mi.setMosdepthDirectory(cli.get(ARG_MOS_DIR),
                            cli.has(ARG_MOS_EXT) ? cli.get(ARG_MOS_EXT) : ".bed");
    mi.setCRAMReadDirectory(cli.get(ARG_CRAMCOUNT_DIR));
    mi.setJobID(cli.has("jobID") ? cli.get("jobID") : null);
    try {
      mi.run();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (Elision e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }

  private static void tryCLI2(String[] args) throws ParseException {
    CLI cli = generateBaseCLI();
    cli.addArg(ARG_BIN_BED, DESC_BIN_BED, true);
    cli.addArg(ARG_SELECTED_SNP_VCF, DESC_SELECTED_SNP_VCF, true);
    cli.addArg(ARG_GENO_VCF, DESC_GENO_VCF, false);
    cli.addArg(ARG_CRAM_DIR, DESC_CRAM_DIR, true);
    cli.addArg(ARG_CRAMCOUNT_DIR, DESC_CRAMCOUNT_DIR, false);
    cli.addArg(ARG_MOS_DIR, DESC_MOS_DIR, false);
    cli.addArg(ARG_MOS_EXT, DESC_MOS_EXT, ".bed", false);

    cli.parse(args);

    String preMos = cli.has(ARG_MOS_DIR) ? null : cli.get(ARG_MOS_DIR);
    String preCnt = cli.has(ARG_CRAMCOUNT_DIR) ? null : cli.get(ARG_CRAMCOUNT_DIR);
    if (preMos == null) {
      new Logger().reportTime("No read depth files specified - pre-processing CRAM files with Mosdepth.  If there are a large number of files, consider performing this step separately and starting again.");
      preMos = runMosdepth(cli.get(ARG_CRAM_DIR), cli.get("projDir"), cli.get(ARG_BIN_BED),
                           GenomeBuild.valueOf(cli.get("build").toUpperCase()),
                           cli.getI(CLI.ARG_THREADS));
    }
    if (preCnt == null) {
      new Logger().reportTime("No allele count files specified - pre-processing CRAM files with "
                              + CRAMSnpReader.class.getName()
                              + ".  If there are a large number of files, consider performing this step separately and starting again.");
      preCnt = runCRAMAlleleCounter(cli.get(ARG_SELECTED_SNP_VCF), cli.get(ARG_CRAM_DIR),
                                    ext.verifyDirFormat(cli.get("projDir")),
                                    GenomeBuild.valueOf(cli.get("build").toUpperCase()));
    }

    MosdepthPipeline mi = new MosdepthPipeline();
    mi.setProjectDir(cli.get("projDir"));
    mi.setProjectName(cli.get("projName"));
    mi.setProjectPropertiesDir(cli.get("propDir"));
    mi.setNumThreads(cli.has(CLI.ARG_THREADS) ? cli.getI(CLI.ARG_THREADS)
                                              : Runtime.getRuntime().availableProcessors());

    mi.setBinsToUseBED(cli.get(ARG_BIN_BED));
    if (cli.has(ARG_GENO_VCF)) {
      mi.setGenotypeVCF(cli.get(ARG_GENO_VCF));
    }
    mi.setSelectedMarkerVCF(cli.get(ARG_SELECTED_SNP_VCF));
    mi.setMosdepthDirectory(preMos, cli.has(ARG_MOS_EXT) ? cli.get(ARG_MOS_EXT) : ".bed");
    mi.setCRAMReadDirectory(preCnt);
    mi.setJobID(cli.has("jobID") ? cli.get("jobID") : null);
    try {
      mi.run();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (Elision e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }

  private static String runMosdepth(String cramDir, String projDir, String bedFile,
                                    GenomeBuild build, int numThreads) {
    Logger log = new Logger();
    String fasta = Resources.genome(build, log).getFASTA().get();
    String[] crams = Files.list(cramDir, "", ".cram", false, false);
    if (crams.length == 0) {
      crams = Files.list(cramDir, "", ".bam", false, false);
    }
    if (crams.length == 0) {
      throw new IllegalArgumentException("No .cram or .bam files found in " + cramDir + ".");
    }
    new File(projDir + "tempMos/").mkdir();
    ExecutorService exec = Executors.newFixedThreadPool(Math.min(numThreads, crams.length));
    for (String cram : crams) {
      exec.submit(new Runnable() {

        public void run() {
          try {
            Command command = Command.builder("mosdepth", "-n", "-b", bedFile, "-f", fasta,
                                              ext.rootOf(cram, true) + ".mos", cramDir + cram)
                                     .dir(projDir + "tempMos/")
                                     .expectedOutputFiles(projDir + "tempMos/"
                                                          + ext.rootOf(cram, true)
                                                          + ".mos.regions.bed.gz")
                                     .build();
            CmdLine.builder(log).build().run(command);
          } catch (Exception e) {
            e.printStackTrace();
          }
        }
      });
    }
    exec.shutdown();
    try {
      exec.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
    } catch (InterruptedException e) {}
    return projDir + "tempMos/";
  }

  private static String runCRAMAlleleCounter(String selectedSNPVCF, String cramDir, String projDir,
                                             GenomeBuild build) {
    System.setProperty("samjdk.reference_fasta",
                       Resources.genome(build, new Logger()).getFASTA().get());
    System.setProperty("reference_fasta", Resources.genome(build, new Logger()).getFASTA().get());
    new File(projDir + "tempSrc/").mkdir();
    new CRAMSnpReader(selectedSNPVCF, cramDir, projDir + "tempSrc/", build, true).run();
    return projDir + "tempSrc/";
  }

  public static void main(String[] args) {
    if (args.length == 0 && Files.isWindows()) {
      MosdepthPipeline mi = new MosdepthPipeline();
      mi.setProjectDir("G:\\bamTesting\\topmed\\project\\");
      mi.setProjectName("TopmedMosdepth");
      mi.setProjectPropertiesDir("D:\\projects\\");
      mi.setNumThreads(Runtime.getRuntime().availableProcessors());

      mi.setBinsToUseBED("G:\\bamTesting\\topmed\\ReferenceGenomeBins_chrs_XYM_hg38.bed");
      // mi.setGenotypeVCF("G:\\bamTesting\\EwingWGS\\ES_recalibrated_snps_indels.vcf.gz");
      mi.setSelectedMarkerVCF("G:\\bamTesting\\topmed\\selected_topmed.vcf");
      mi.setMosdepthDirectory("G:\\bamTesting\\topmed\\00mos\\", ".bed");
      mi.setCRAMReadDirectory("G:\\bamTesting\\topmed\\00cnt\\");
      mi.setJobID(null);
      try {
        mi.run();
      } catch (IOException | Elision e) {
        e.printStackTrace();
      }
    } else {
      try {
        tryCLI1(args);
      } catch (ParseException e) {
        try {
          tryCLI2(args);
        } catch (ParseException e1) {
          // TODO display full error message with all options and various execution paths
        }
      }

    }
  }

}
