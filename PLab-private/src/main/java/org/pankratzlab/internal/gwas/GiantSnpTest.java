package org.pankratzlab.internal.gwas;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;

import org.pankratzlab.common.CLI;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.qsub.Qsub;

public class GiantSnpTest extends SnpTest {

  public void runGiantSnpTest(SnpTest run) throws IOException {
    for (int i = 1; i < 27; i++) {
      String dataFile = run.dataDirectory + "chr" + i + ".dose.vcf.gz";
      if (!org.pankratzlab.common.Files.exists(dataFile)) {
        log.report("No data found for chromosome " + i);
        continue;
      }
      Files.createDirectories(Paths.get("./chr" + i));
      String outFile = "./chr" + i + "/output_chr" + i + ".out";
      String excludes = run.excludesFile;
      String includes = run.includesFile;
      String cmd = new SnpTestCommand(run.snpTestExec, dataFile, run.sampleFile, true, run.vcfField,
                                      run.pheno, run.covars, i, outFile, excludes, includes)
                                                                                            .getCommand();
      try (PrintWriter in = new PrintWriter("./chr" + i + "/input.txt")) {
        in.write(cmd);
      } catch (Exception e) {
        e.printStackTrace();
        run.log.reportError("Writing the SnpTest input broke. Check your arguments.");
      }

      StringBuilder scriptExecCmd = new StringBuilder("cd ").append(ext.pwd() + "chr" + i)
                                                            .append("\n")
                                                            .append(org.pankratzlab.common.Files.getRunString(59392))
                                                            .append(" org.genvisis.one.ScriptExecutor token=finito threads="
                                                                    + Runtime.getRuntime()
                                                                             .availableProcessors());
      Qsub.qsub("./chr" + i + "/runSnpTestChr" + i + ".qsub", scriptExecCmd.toString(), 61440, 96,
                1);
    }
  }

  private static final String ARG_SNPTEST = "snpTest";
  private static final String ARG_DATADIR = "dataDir";
  private static final String ARG_INFOTEMP = "infoTemplate";
  private static final String ARG_SAMPLE = "sample";
  private static final String ARG_PHENO = "pheno";
  private static final String ARG_COVARS = "covars";
  private static final String ARG_REPL = "repl";
  private static final String ARG_EXCLUDES = "exclude";
  private static final String ARG_INCLUDES = "include";

  private static final String DESC_SNPTEST = "SnpTest executable (full path included unless on path)";
  private static final String DESC_DATADIR = "Data file directory";
  private static final String DESC_INFOTEMP = "Info filename template (e.g. chr#.info.gz) - assumes chromosomally-separated data files.";
  private static final String DESC_SAMPLE = "Sample file (SnpTest-format)";
  private static final String DESC_PHENO = "Name of phenotype to test";
  private static final String DESC_COVARS = "Covariates to test, leave blank to use all available (uses the -cov_all flag)";
  private static final String DESC_REPL = "Special character in data file template into which chromosome number will be placed.";
  private static final String DESC_EXCLUDES = "File with a list of samples that should be excluded from analysis. These should match samples from the sample file.";
  private static final String DESC_INCLUDES = "File with a list of samples that should be included. Samples not in this file will be excluded.";

  public static void main(String[] args) throws IOException {
    CLI cli = new CLI(SnpTest.class);

    cli.addArg(ARG_SNPTEST, DESC_SNPTEST);
    cli.addArg(ARG_DATADIR, DESC_DATADIR);
    cli.addArg(ARG_INFOTEMP, DESC_INFOTEMP);
    cli.addArg(ARG_SAMPLE, DESC_SAMPLE);
    cli.addArg(ARG_PHENO, DESC_PHENO);

    cli.addArg(ARG_COVARS, DESC_COVARS, false);
    cli.addArg(ARG_REPL, DESC_REPL, false);
    cli.addArg(ARG_EXCLUDES, DESC_EXCLUDES, false);
    cli.addArg(ARG_INCLUDES, DESC_INCLUDES, false);

    cli.parseWithExit(args);

    String snpTest = cli.get(ARG_SNPTEST);
    String samp = cli.get(ARG_SAMPLE);
    String phen = cli.get(ARG_PHENO);
    String dataDir = cli.get(ARG_DATADIR);

    SnpTest run = new SnpTest();
    GiantSnpTest giant = new GiantSnpTest();
    run.setDataDirectory(dataDir);
    run.snpTestExec = snpTest;
    run.setSampleFile(samp);
    run.setPheno(phen);

    if (cli.has(ARG_COVARS)) {
      run.setCovars(cli.get(ARG_COVARS).split(","));
    }
    if (cli.has(ARG_REPL)) {
      run.setRepl(cli.get(ARG_REPL));
    }
    if (cli.has(ARG_EXCLUDES)) {
      run.setExcludes(cli.get(ARG_EXCLUDES));
    }
    if (cli.has(ARG_INCLUDES)) {
      run.setIncludes(cli.get(ARG_INCLUDES));
    }

    giant.runGiantSnpTest(run);

  }

}
