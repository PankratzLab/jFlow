/**
 * 
 */
package org.pankratzlab.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Bundled;
import org.pankratzlab.common.CmdLine;
import org.pankratzlab.common.Command;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.WorkerHive;
import org.pankratzlab.common.ext;
import org.pankratzlab.core.CLI;

/**
 * Class to calculate power of SKAT tests
 */
public class SKATPowerCalculator {

  /**
   * 
   */
  private static final String ALPHA = "alpha";
  /**
   * 
   */
  private static final String POWER_CALC_TXT = "PowerCalc.txt";
  /**
   * 
   */
  private static final String ANALYSIS_NAME = "analysisName";
  /**
   * 
   */
  private static final String CAUSAL_MAF_CUTOFFS = "causalMAFCutoffs";
  /**
   * 
   */
  private static final String CAUSAL_PERCENTS = "causalPercents";
  /**
   * 
   */
  private static final String PREVALENCE = "prevalence";
  /**
   * 
   */
  private static final String SUB_REGION_LENGTH = "subRegionLength";
  /**
   * 
   */
  private static final String CONTROLS = "controls";
  /**
   * 
   */
  private static final String CASES = "cases";

  private static List<String> loadSkatPowerScript(Logger log) {

    List<String> lines = new ArrayList<>();
    BufferedReader br;
    try {
      br = new BufferedReader(new InputStreamReader(getskatPowerScriptStream()));
      String sCurrentLine;

      while ((sCurrentLine = br.readLine()) != null) {
        lines.add(sCurrentLine.trim());
      }
      br.close();
    } catch (IOException e) {
      log.reportException(e);
    }

    return lines;

  }

  private static InputStream getskatPowerScriptStream() {
    return Bundled.getStream("skatPower.Rscript");
  }

  private static List<String> generateArgs(CLI c, String rootOut, double causalPercent,
                                           double causalMAFCutoffs) {
    String outFile = getPowerCalcOut(rootOut, causalPercent, causalMAFCutoffs) + POWER_CALC_TXT;

    List<String> args = new ArrayList<>();

    args.add("outFile=\"" + outFile + "\"");
    args.add("cases = " + c.getI(CASES));
    args.add("controls = " + c.getI(CONTROLS));
    args.add("Haplotypes = NULL");
    args.add("SNP.Location = NULL");
    args.add("SubRegion.Length  = " + c.getI(SUB_REGION_LENGTH));
    args.add("Prevalence  = " + c.getD(PREVALENCE));
    args.add("Case.Prop  = cases / (cases + controls)");
    args.add("Causal.Percent = " + causalPercent);
    args.add("Causal.MAF.Cutoff = " + causalMAFCutoffs);
    args.add("alpha = " + c.getD(ALPHA));
    args.add("N.Sample.ALL = (cases + controls)");
    args.add("Weight.Param = c(1, 25)");
    args.add("N.Sim = 100");
    args.add("OR.Type = \"Log\"");
    args.add("MaxOR = 2");
    args.add("Negative.Percent = 0");
    return args;
  }

  private static String getPowerCalcOut(String rootOut, double causalPercent,
                                        double causalMAFCutoff) {
    return rootOut + "Causal.Percent_" + causalPercent + "_Causal.MAF.Cutoff_" + causalMAFCutoff;
  }

  private static List<String> generateScript(List<String> baseSkatPowerScript, CLI c,
                                             String rootOut, double causalPercent,
                                             double causalMAFCutoff) {
    List<String> script = new ArrayList<>();
    script.addAll(generateArgs(c, rootOut, causalPercent, causalMAFCutoff));
    script.addAll(baseSkatPowerScript);
    return script;

  }

  private static void run(CLI c) {
    new File(c.get(CLI.ARG_OUTDIR)).mkdirs();
    String rootOut = c.get(CLI.ARG_OUTDIR) + c.get(ANALYSIS_NAME);
    Logger log = new Logger(rootOut + ".log");
    log.reportTimeInfo("Reporting to root out " + rootOut);
    double[] causalPercents = ArrayUtils.toDoubleArray(c.get(CAUSAL_PERCENTS).split(","));
    log.reportTimeInfo(causalPercents.length + " " + CAUSAL_PERCENTS);

    double[] causalMAFCutoffs = ArrayUtils.toDoubleArray(c.get(CAUSAL_MAF_CUTOFFS).split(","));
    log.reportTimeInfo(causalMAFCutoffs.length + " " + CAUSAL_MAF_CUTOFFS);

    List<String> baseSkatPowerScript = loadSkatPowerScript(log);

    WorkerHive<Boolean> hive = new WorkerHive<>(c.getI(CLI.ARG_THREADS), 10, log);
    for (double causalPercent : causalPercents) {
      for (double causalMAFCutoff : causalMAFCutoffs) {
        String currentScript = getPowerCalcOut(rootOut, causalPercent, causalMAFCutoff)
                               + ".Rscript";
        Files.writeIterable(generateScript(baseSkatPowerScript, c, rootOut, causalPercent,
                                           causalMAFCutoff),
                            currentScript);
        final Command command = Command.basic("Rscript", currentScript);
        hive.addCallable(() -> CmdLine.basic(log).run(command));
      }
    }
    hive.execute(true);
    summarize(rootOut, log);
  }

  private static boolean summarize(String rootOut, Logger log) {
    String dir = ext.parseDirectoryOfFile(rootOut);
    String out = rootOut + "summary.txt";
    List<String> script = new ArrayList<>();
    script.add("library(reshape2)");
    script.add("summaries <- do.call(rbind,lapply(list.files(\"" + dir + "\",pattern = \""
               + POWER_CALC_TXT + "$\",full = TRUE,recursive = FALSE),read.delim))");
    script.add("write.table(summaries,file = paste0(\"" + out
               + "\") ,sep = \"\\t\",quote = FALSE, row.names = FALSE)");
    script.add("write.table(acast(summaries, Causal.Percent ~ Causal.MAF.Cutoff, value.var = \"Power\"),file = paste0(\""
               + out + ".mat\") ,sep = \"\\t\",quote = FALSE, row.names = TRUE)");
    String currentScript = out + ".Rscript";
    Files.writeIterable(script, currentScript);
    return CmdLine.basic(log).run(Command.basic("Rscript", currentScript));

  }

  public static void main(String[] args) {
    CLI c = new CLI(SKATPowerCalculator.class);

    c.addArg(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, true);
    c.addArgWithDefault(CLI.ARG_THREADS, CLI.DESC_THREADS, "24");

    c.addArgWithDefault(ANALYSIS_NAME, "name of the analysis", "skat.power");
    c.addArgWithDefault(PREVALENCE, "a value of disease prevalence", "0.10");
    c.addArgWithDefault(CASES, "number of cases", "21250");
    c.addArgWithDefault(CONTROLS, "number of controls", "90000");
    c.addArgWithDefault(CAUSAL_PERCENTS,
                        "comma-delimited list of a value of the percentage of causal SNPs among rare SNPs (MAF < Causal.MAF.Cutoff)",
                        "1,2,5,10,20,30,40,50");
    c.addArgWithDefault(ALPHA, "significance level", "1E-6");
    c.addArgWithDefault(CAUSAL_MAF_CUTOFFS,
                        "comma-delimited lista value of MAF cutoff for the causal SNPs. Only SNPs that have MAFs smaller than this are considered as causal SNPs",
                        "0.005, 0.01, 0.02, 0.03, 0.05, 0.10");

    c.addArgWithDefault(SUB_REGION_LENGTH,
                        "a value of the length of subregions, set to something like 3000 to drastically decrease runtime (-1 for no subregions)",
                        "-1");

    c.parseWithExit(args);
    run(c);

  }

}
