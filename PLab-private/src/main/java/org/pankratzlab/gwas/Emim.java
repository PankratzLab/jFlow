package org.pankratzlab.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.genvisis.cnv.Resources;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ext;
import org.pankratzlab.shared.qsub.Qsub;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.ImmutableSortedSet;
import com.google.common.collect.Sets;
import com.google.common.io.Closeables;

public class Emim {

  public enum EMIM_MODEL {
    GENOTYPIC("Genotypic", 2, true),
    DOMINANT("Dominant", 1, true),
    ADDITIVE("Additive", 1, true),
    IMPRINTING_MATERNAL("Imprinting_M", 1, false),
    IMPRINTING_PATERNAL("Imprinting_P", 1, false);

    private static final Set<EMIM_MODEL> VALUE_SET;
    private static final Set<EMIM_MODEL> OPTIONAL_SET;

    static {
      Set<EMIM_MODEL> valueSet = Sets.newHashSet();
      Set<EMIM_MODEL> optionalSet = Sets.newLinkedHashSet();
      for (EMIM_MODEL model : EMIM_MODEL.values()) {
        valueSet.add(model);
        if (model.optional) {
          optionalSet.add(model);
        }
      }
      VALUE_SET = ImmutableSet.copyOf(valueSet);
      OPTIONAL_SET = ImmutableSortedSet.copyOf(optionalSet);
    }

    private final String name;
    private final int degreesFreedom;
    private final boolean optional;

    EMIM_MODEL(String name, int degreesFreedom, boolean optional) {
      this.name = name;
      this.degreesFreedom = degreesFreedom;
      this.optional = optional;
    }

    @Override
    public String toString() {
      return name;
    }

    public int getDegreesOfFreedom() {
      return degreesFreedom;
    }

    public boolean isOptional() {
      return optional;
    }

    public static Set<EMIM_MODEL> valueSet() {
      return VALUE_SET;
    }

    public static Set<EMIM_MODEL> optionalSet() {
      return OPTIONAL_SET;
    }

  }

  private enum EMIM_PARAM {
    ESTIMATE_R1("   << estimate R1 (0=no, 1=yes)"),
    ESTIMATE_R2("   << estimate R2 (0=no, 1=yes)"),
    R2_EQUALS_R1("   << R2=R1 (0=no, 1=yes)"),
    R2_EQUALS_R1_SQUARED("   << R2=R1squared\t(0=no, 1=yes)"),
    ESTIMATE_S1("   << estimate S1 (0=no, 1=yes)"),
    ESTIMATE_S2("   << estimate S2 (0=no, 1=yes)"),
    S2_EQUALS_S1("   << S2=S1 (0=no, 1=yes)"),
    S2_EQUALS_S1_SQUARED("   << S2=S1squared\t(0=no, 1=yes)"),
    ESTIMATE_IM("   << estimate Im (0=no, 1=yes)"),
    ESTIMATE_IP("   << estimate Ip (0=no, 1=yes)");

    private final String lineSuffix;

    EMIM_PARAM(String lineSuffix) {
      this.lineSuffix = lineSuffix;
    }

    String[] getReplacement(boolean setTo) {
      return setTo ? new String[] {"0" + lineSuffix, "1" + lineSuffix}
                   : new String[] {"1" + lineSuffix, "0" + lineSuffix};
    }
  }

  private static final String EMIM_RISK_SNP_FILENAME = "risksnplist.txt";
  private static final String EMIM_MINOR_SNP_FILENAME = "minorsnplist.txt";

  private static void generateEmimParams(String filenameOriginal, HashSet<EMIM_PARAM> setParams,
                                         Logger log) {
    String[][] replacements = new String[EMIM_PARAM.values().length][];
    for (int i = 0; i < EMIM_PARAM.values().length; i++) {
      EMIM_PARAM param = EMIM_PARAM.values()[i];
      replacements[i] = param.getReplacement(setParams.contains(param));
    }

    replaceLines(filenameOriginal, "emimparams.dat", replacements, log);
  }

  private static void generateRiskAlleles(String riskAlleleFile, Logger log) {
    if (!Files.exists(riskAlleleFile)) {
      log.reportError(riskAlleleFile
                      + " is not in the current directory; aborting risk allele generation");
      return;
    }
    if (!Files.exists(EMIM_MINOR_SNP_FILENAME)) {
      log.reportError(EMIM_MINOR_SNP_FILENAME
                      + " is not in the current directory; aborting risk allele generation");
      return;
    }

    Map<String, String> riskAlleles = HashVec.loadFileToHashString(riskAlleleFile, false);
    BufferedReader reader = null;
    PrintWriter writer = Files.getAppropriateWriter(EMIM_RISK_SNP_FILENAME);

    try {
      reader = Files.getAppropriateReader(EMIM_MINOR_SNP_FILENAME);
      while (reader.ready()) {
        String line = reader.readLine();
        String[] linePieces = line.split(PSF.Regex.GREEDY_WHITESPACE);
        if (linePieces.length != 3) {
          log.reportError("Invalid line in " + EMIM_MINOR_SNP_FILENAME + ": " + line
                          + "; aborting risk allele generation");
          return;
        }
        String markerNum = linePieces[0];
        String markerID = linePieces[1];
        String allele = linePieces[2];
        if (riskAlleles.containsKey(markerID)) {
          allele = riskAlleles.get(markerID);
        }

        writer.println(markerNum + "\t" + markerID + "\t" + allele);
      }
    } catch (FileNotFoundException e) {
      log.reportFileNotFound(EMIM_MINOR_SNP_FILENAME);
    } catch (IOException e) {
      log.reportIOException(EMIM_MINOR_SNP_FILENAME);
    } finally {
      Closeables.closeQuietly(reader);
      writer.close();
    }

  }

  private static void setTo(String runType, EMIM_MODEL model) {
    String filenameOriginal;
    Logger log;

    if (!Files.exists("emimparams.dat")) {
      System.err.println("Error - emimparams.dat is not in the current directory; aborting rewrite");
      return;
    }

    log = new Logger();
    filenameOriginal = Files.backup("emimparams.dat", "./", "./", true);

    HashSet<EMIM_PARAM> setParams = new HashSet<>();
    if (runType.equals("C")) {
      switch (model) {
        case GENOTYPIC:
          setParams.add(EMIM_PARAM.ESTIMATE_R1);
          setParams.add(EMIM_PARAM.ESTIMATE_R2);
          break;
        case DOMINANT:
          setParams.add(EMIM_PARAM.R2_EQUALS_R1);
          break;
        case ADDITIVE:
          setParams.add(EMIM_PARAM.R2_EQUALS_R1_SQUARED);
          break;
        default:
          log.reportError("Invalid EMIM model for " + runType + " run type: " + model.toString());
          return;
      }

    } else if (runType.equals("CM")) {
      switch (model) {
        case GENOTYPIC:
          setParams.add(EMIM_PARAM.ESTIMATE_R1);
          setParams.add(EMIM_PARAM.ESTIMATE_R2);
          setParams.add(EMIM_PARAM.ESTIMATE_S1);
          setParams.add(EMIM_PARAM.ESTIMATE_S2);
          break;
        case DOMINANT:
          setParams.add(EMIM_PARAM.R2_EQUALS_R1);
          setParams.add(EMIM_PARAM.S2_EQUALS_S1);
          break;
        case ADDITIVE:
          setParams.add(EMIM_PARAM.R2_EQUALS_R1_SQUARED);
          setParams.add(EMIM_PARAM.S2_EQUALS_S1_SQUARED);
          break;
        default:
          log.reportError("Invalid EMIM model for " + runType + " run type: " + model.toString());
          return;
      }
    } else if (runType.equals("CIm")) {
      switch (model) {
        case GENOTYPIC:
          setParams.add(EMIM_PARAM.ESTIMATE_R1);
          setParams.add(EMIM_PARAM.ESTIMATE_R2);
          setParams.add(EMIM_PARAM.ESTIMATE_IM);
          break;
        case DOMINANT:
          setParams.add(EMIM_PARAM.R2_EQUALS_R1);
          setParams.add(EMIM_PARAM.ESTIMATE_IM);
          break;
        case ADDITIVE:
          setParams.add(EMIM_PARAM.R2_EQUALS_R1_SQUARED);
          setParams.add(EMIM_PARAM.ESTIMATE_IM);
          break;
        default:
          log.reportError("Invalid EMIM model for " + runType + " run type: " + model.toString());
          return;
      }
    } else if (runType.equals("CIp")) {
      switch (model) {
        case GENOTYPIC:
          setParams.add(EMIM_PARAM.ESTIMATE_R1);
          setParams.add(EMIM_PARAM.ESTIMATE_R2);
          setParams.add(EMIM_PARAM.ESTIMATE_IP);
          break;
        case DOMINANT:
          setParams.add(EMIM_PARAM.R2_EQUALS_R1);
          setParams.add(EMIM_PARAM.ESTIMATE_IP);
          break;
        case ADDITIVE:
          setParams.add(EMIM_PARAM.R2_EQUALS_R1_SQUARED);
          setParams.add(EMIM_PARAM.ESTIMATE_IP);
          break;
        default:
          log.reportError("Invalid EMIM model for " + runType + " run type: " + model.toString());
          return;
      }
    } else if (runType.equals("M")) {
      switch (model) {
        case GENOTYPIC:
          setParams.add(EMIM_PARAM.ESTIMATE_S1);
          setParams.add(EMIM_PARAM.ESTIMATE_S2);
          break;
        case DOMINANT:
          setParams.add(EMIM_PARAM.S2_EQUALS_S1);
          break;
        case ADDITIVE:
          setParams.add(EMIM_PARAM.S2_EQUALS_S1_SQUARED);
          break;
        default:
          log.reportError("Invalid EMIM model for " + runType + " run type: " + model.toString());
          return;
      }
    } else if (runType.equals("POO")) {
      switch (model) {
        case IMPRINTING_MATERNAL:
          setParams.add(EMIM_PARAM.ESTIMATE_IM);
          break;
        case IMPRINTING_PATERNAL:
          setParams.add(EMIM_PARAM.ESTIMATE_IP);
          break;
        default:
          log.reportError("Invalid EMIM model for " + runType + " run type: " + model.toString());
      }
    } else {
      log.reportError("Invalid Run Type: " + runType);
      return;
    }

    generateEmimParams(filenameOriginal, setParams, log);

    System.out.println("emimparams.dat set to " + model.toString() + " model of " + runType
                       + " effect");

  }

  private static void listSexMarkers(String bimFile, String sexFile) throws NumberFormatException,
                                                                     IOException {
    PrintWriter writer;
    BufferedReader reader;
    String line;
    String[] parts;
    int chr;

    writer = Files.getAppropriateWriter(sexFile);
    reader = Files.getAppropriateReader(bimFile);
    line = null;
    while ((line = reader.readLine()) != null) {
      parts = line.trim().split(PSF.Regex.GREEDY_WHITESPACE, -1);
      chr = Integer.parseInt(parts[0]);
      if (chr > 22) {
        writer.println(parts[1]);
      }
    }
    writer.flush();
    writer.close();
    reader.close();
  }

  protected static String scriptAllInDir(String runDir, String plinkDirAndRoot,
                                         String relativePlinkRoot, String excludeFile,
                                         String keepFile, String riskAlleleFile, double pThreshold,
                                         Collection<EMIM_MODEL> models, boolean phaseWithShapeit,
                                         String resultPrefix, Logger log) {
    String commands;
    String currDir = ext.verifyDirFormat(runDir);
    boolean forceRun = false, forceParse = false;
    if (currDir.charAt(0) != '/' && !currDir.contains(":")) {
      currDir = (new File("./" + currDir)).getAbsolutePath() + "/";
    }
    final String premimPBS = currDir + ext.rootOf(plinkDirAndRoot, true) + "_runPremim.pbs";
    final String emimPBS = currDir + ext.rootOf(plinkDirAndRoot, true) + "_runEmim.pbs";
    boolean runPremimPBS = false;

    if (excludeFile.equals("GEN")) {
      excludeFile = "sexChrMarkers.txt";
      if (Files.exists(currDir + excludeFile)) {
        log.report(currDir + excludeFile + " already exists, skipping exclude file generation");
      } else {
        try {
          listSexMarkers(plinkDirAndRoot + ".bim", currDir + excludeFile);
          forceRun = true;
          forceParse = true;
        } catch (Exception e) {
          excludeFile = null;
          e.printStackTrace();
        }
      }
    }

    commands = "";

    if (!forceRun && Files.checkAllFiles(currDir, true, false, log, "plink_prep.log",
                                         "emimPrep.bed", "emimPrep.bim", "emimPrep.fam")) {
      log.report(currDir + "emimPrep PLINK files already exist, skipping PLINK file generation");
    } else {
      forceRun = true;
      forceParse = true;
      commands += "plink2 --noweb --allow-no-sex --bfile " + relativePlinkRoot
                  + (excludeFile != null ? " --exclude " + excludeFile : "")
                  + (keepFile != null ? " --keep " + keepFile : "")
                  + (riskAlleleFile != null ? " --a1-allele " + riskAlleleFile : "")
                  + " --make-bed --out emimPrep\n" + "mv emimPrep.log plink_prep.log\n" + "\n";
    }

    if (!forceRun
        && Files.checkAllFiles(currDir, true, false, log, "plink_mendel.log", "plink.mendel",
                               "plink.lmendel", "plink.fmendel", "plink.imendel")) {
      log.report(currDir
                 + "plink.*mendel files already exist, skipping mendelian error calculation");
    } else {
      forceParse = true;
      commands += "plink2 --noweb --allow-no-sex --bfile emimPrep --keep-allele-order --mendel\n"
                  + "mv plink.log plink_mendel.log\n" + "\n";
    }

    if (!forceRun && Files.checkAllFiles(currDir, true, false, log, "plink_tdt.log", "plink.tdt")) {
      log.report(currDir + "plink.tdt already exists, skipping TDT");
    } else {
      forceParse = true;
      commands += "plink2 --noweb --allow-no-sex --bfile emimPrep --keep-allele-order --tdt --ci 0.95\n"
                  + "mv plink.log plink_tdt.log\n" + "\n";
    }

    if (!forceRun && Files.checkAllFiles(currDir, true, false, log, "plink_hwe.log", "plink.hwe")) {
      log.report(currDir
                 + "plink.hwe already exists, skipping Hardy-Weinberg Equilibrium calculation");
    } else {
      forceParse = true;
      commands += "plink2 --noweb --allow-no-sex --bfile emimPrep --keep-allele-order --hardy\n"
                  + "mv plink.log plink_hwe.log\n" + "\n";
    }

    if (!forceRun
        && Files.checkAllFiles(currDir, true, false, log, "plink_freq.log", "plink.frq")) {
      log.report(currDir + "plink.frq already exists, skipping Minor Allele Frequency calculation");
    } else {
      forceParse = true;
      commands += "plink2 --noweb --allow-no-sex --bfile emimPrep --keep-allele-order --freq\n"
                  + "mv plink.log plink_freq.log\n" + "\n";
    }

    if (!forceRun && Files.checkAllFiles(currDir, true, false, log, EMIM_MINOR_SNP_FILENAME,
                                         EMIM_RISK_SNP_FILENAME)) {
      log.report(EMIM_MINOR_SNP_FILENAME + " and " + EMIM_RISK_SNP_FILENAME
                 + " already exist, skipping risk allele designation");
    } else {
      commands += "\npremim -rout " + EMIM_MINOR_SNP_FILENAME + " emimPrep.bed\n";
      if (riskAlleleFile != null) {
        commands += Files.getRunString() + " gwas.Emim riskAlleles=" + riskAlleleFile + "\n\n";
      } else {
        commands += "cp " + EMIM_MINOR_SNP_FILENAME + " " + EMIM_RISK_SNP_FILENAME + "\n";
      }

    }

    if (!forceRun
        && Files.checkAllFiles(currDir, true, false, log, "premim.log", "emimparams.dat",
                               "emimmarkers.dat", "caseparenttrios.dat", "caseparents.dat",
                               "casemotherduos.dat", "casefatherduos.dat", "casemothers.dat",
                               "casefathers.dat", "cases.dat", "conparents.dat",
                               "conmotherduos.dat", "confatherduos.dat", "cons.dat")) {
      log.report("Outputs of PREMIM in " + currDir + " already exist, skipping PREMIM");
    } else {
      forceRun = true;
      forceParse = true;
      if (phaseWithShapeit) {
        commands += "premim -im -a -ihap -shapeit " + Resources.shapeit(log).getShapeit().get()
                    + " -shapeit-thread 24 -rfile risksnplist.txt emimPrep.bed\n";
        commands += "\nqsub " + emimPBS;
        commands = "cd " + currDir + "\n" + commands;
        Qsub.qsub(premimPBS, commands, 62000, 24, 24);
        commands = "";
        runPremimPBS = true;
      } else {
        commands += "premim -cg -a -rfile risksnplist.txt emimPrep.bed\n";
      }

    }
    String parseCommands = "";
    for (EMIM_MODEL model : models) {
      Set<String> runIDs = Sets.newHashSet();
      if (model.equals(EMIM_MODEL.IMPRINTING_PATERNAL)
          || model.equals(EMIM_MODEL.IMPRINTING_MATERNAL)) {
        runIDs.add("POO");
      } else {
        runIDs.add("C");
        runIDs.add("M");
        runIDs.add("CM");
        runIDs.add("CIm");
        runIDs.add("CIp");
      }
      boolean parseModel = false;
      for (String runID : runIDs) {
        Set<String> outputs = Sets.newHashSet();
        outputs.add("emimsummary_" + runID + "_" + model.toString() + ".out");
        outputs.add("emimparams_" + runID + "_" + model.toString() + ".dat");

        if (!forceRun && Files.checkAllFiles(currDir, outputs, false, false, log)) {
          log.report("Results already exist in " + currDir + " for " + runID + "_"
                     + model.toString() + " model, skipping " + runID + "_" + model.toString()
                     + " EMIM");
        } else {
          parseModel = true;
          commands += "\n" + Files.getRunString() + " gwas.Emim run=" + runID + " model="
                      + model.toString() + "\n" + "emim\n" + "mv emimsummary.out emimsummary_"
                      + runID + "_" + model.toString() + ".out\n" + "rm emimresults.out\n"
                      + "cp emimparams.dat emimparams_" + runID + "_" + model.toString() + ".dat\n";
        }
      }

      if (model.isOptional()
          && (forceParse || parseModel
              || !Files.exists(currDir + formParsedOutFileName(resultPrefix, model)))) {
        parseCommands += Files.getRunString() + " gwas.Emim parse=./" + " hwe=plink.hwe"
                         + " frq=plink.frq" + " pThreshold=" + pThreshold + " model="
                         + model.toString()
                         + (resultPrefix == null ? "" : " resultPrefix=" + resultPrefix) + "\n\n";
      }
    }

    commands += parseCommands;

    if (commands.equals("") && !runPremimPBS) {
      log.report("Warning - No EMIM commands were generated for " + currDir
                 + ", remove directory to re-run EMIM Pipeline from scratch");
      return null;
    }
    commands = "cd " + currDir + "\n" + commands;

    Qsub.qsub(emimPBS, commands, 32000, 12, 1);

    return runPremimPBS ? premimPBS : emimPBS;
  }

  public static String scriptAll(String plinkPrefix, String excludeFile, String keepFile,
                                 double pThreshold) {
    return scriptAllInDir("./", plinkPrefix, plinkPrefix, excludeFile, keepFile, null, pThreshold,
                          EMIM_MODEL.valueSet(), true, null, new Logger());
    // String commands;
    // String currDir;
    //
    // if (excludeFile.equals("GEN")) {
    // try {
    // listSexMarkers(plinkPrefix + ".bim", "sexChrMarkers.txt");
    // excludeFile = "sexChrMarkers.txt";
    // } catch (NumberFormatException | IOException e) {
    // excludeFile = null;
    // e.printStackTrace();
    // }
    // }
    //
    // currDir = (new File("./")).getAbsolutePath();
    // commands = "cd " + currDir + "\n" +
    // "plink2 --noweb --bfile "
    // + plinkPrefix
    // + (excludeFile != null ? " --exclude " + excludeFile : "")
    // + (keepFile != null ? " --keep " + keepFile : "")
    // + " --make-bed --out emimPrep\n"+
    // "plink2 --noweb --bfile emimPrep --mendel\n"+
    // "plink2 --noweb --bfile emimPrep --tdt\n"+
    // "plink2 --noweb --bfile emimPrep --hardy\n"+
    // "premim -cg -a -rout risksnplist.txt emimPrep.bed\n"+
    // "\n"+
    // "jcp gwas.Emim run=C\n"+
    // "emim\n"+
    // "mv emimsummary.out emimsummary_C.out\n"+
    // "mv emimresults.out emimresults_C.out\n"+
    // "cp emimparams.dat emimparams_C.dat\n"+
    // "\n"+
    // "jcp gwas.Emim run=CM\n"+
    // "emim\n"+
    // "mv emimsummary.out emimsummary_CM.out\n"+
    // "mv emimresults.out emimresults_CM.out\n"+
    // "cp emimparams.dat emimparams_CM.dat\n"+
    // "\n"+
    // "jcp gwas.Emim run=M\n"+
    // "emim\n"+
    // "mv emimsummary.out emimsummary_M.out\n"+
    // "mv emimresults.out emimresults_M.out\n"+
    // "cp emimparams.dat emimparams_M.dat\n"+
    // "\n"+
    // "jcp gwas.Emim parse=./ hwe=plink.hwe pThreshold=" + pThreshold;
    //
    // Files.qsub(plinkPrefix+"_runEmim.pbs", commands, 5000, 24, 1);
  }

  private static String formParsedOutFileName(String resultPrefix, EMIM_MODEL model) {
    return (resultPrefix == null ? "" : resultPrefix + "_") + "results_pVals_" + model.toString()
           + ".xln";
  }

  public static void parse(String dir, String resultPrefix, String hweFile, String frqFile,
                           double pValueThreshold, EMIM_MODEL model) {
    if (!model.isOptional()) {
      // Imprinting models are included for comparison with each Genotype model
      return;
    }
    String resultsFileC = dir + "emimsummary_C_" + model.toString() + ".out";
    String resultsFileM = dir + "emimsummary_M_" + model.toString() + ".out";
    String resultsFileCM = dir + "emimsummary_CM_" + model.toString() + ".out";
    String resultsFileCIm = dir + "emimsummary_CIm_" + model.toString() + ".out";
    String resultsFileCIp = dir + "emimsummary_CIp_" + model.toString() + ".out";
    String resultsFileIm = dir + "emimsummary_POO_Imprinting_M.out";
    String resultsFileIp = dir + "emimsummary_POO_Imprinting_P.out";
    String resultsFileTdt = dir + "plink.tdt";
    String mapFile = dir + "emimPrep.bim";
    String mendelErrorFile = dir + "plink.lmendel";
    String outfile = dir + formParsedOutFileName(resultPrefix, model);

    ResultsPackager.parseEmimFormat(resultsFileC, resultsFileM, resultsFileCM, resultsFileCIm,
                                    resultsFileCIp, resultsFileIm, resultsFileIp, resultsFileTdt,
                                    mapFile, mendelErrorFile, hweFile, frqFile, pValueThreshold,
                                    outfile, new Logger("EMIMparser.log"));
  }

  public static void replaceLines(String filenameOriginal, String filenameWithReplacements,
                                  String[][] relacements, Logger log) {
    BufferedReader reader;
    PrintWriter writer;

    try {
      reader = Files.getAppropriateReader(filenameOriginal);
      writer = Files.openAppropriateWriter(filenameWithReplacements);
      while (reader.ready()) {
        writer.println(ext.replaceAllWith(reader.readLine(), relacements));
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + filenameOriginal + "\" not found in current directory");
      return;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + filenameOriginal + "\"");
      return;
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String runType = "C";
    String riskAlleleFile = null;
    String dir = null;
    double pValueThreshold = 1.1;
    String hweFile = null;
    String frqFile = null;
    String plinkPrefix = null;
    EMIM_MODEL model = EMIM_MODEL.DOMINANT;
    String excludeFile = "GEN";
    String keepFile = null;
    String resultPrefix = null;

    String usage = "\n" + "gwas.Emim requires 0-1 arguments\n"
                   + "   (1) run type (either C, CM, or M) (i.e. run=" + runType + " (default))\n"
                   + "   (2) model (" + ArrayUtils.toStr(EMIM_MODEL.values(), ",")
                   + ") (i.e. model=" + model.toString() + " (default))\n" + "  OR\n"
                   + "   (1) desired risk allele file (i.e. riskAlleles=forceRiskAllele.txt (not the default))\n"
                   + "  OR\n"
                   + "   (1) generate script that runs the full process (i.e. script=plinkPrefix (not the default))\n"
                   + "   (2) p-value threshold to filter on (piped to parse method) (i.e. pThreshold="
                   + pValueThreshold + " (default))\n"
                   + "   (3) (optional) a file of markers to exclude - the default will generate a file of sex markers, as Emim won't parse these. (i.e. exclude=drops.dat (not the default))\n"
                   + "   (4) (optional) a file of tab-delimited FID/IID pairs to keep (i.e. keep=completeTrios.dat (not the default))\n"
                   + "  OR\n" + "   (1) directory to parse (i.e. parse=./ (not the default))\n"
                   + "   (2) file prefix for results (i.e. resultPrefix=" + resultPrefix
                   + " (default))\n" + "   (3) p-value threshold to filter on (i.e. pThreshold="
                   + pValueThreshold + " (default))\n" + "   (4) model "
                   + Arrays.toString(EMIM_MODEL.values()) + " (i.e. model=" + model.toString()
                   + " (default))\n"
                   + "   (5) (optional) plink.hwe file to merge with results (i.e. hwe=" + hweFile
                   + " (default))\n"
                   + "   (6) (optional) plink.frq file to merge with results (i.e. frq=" + frqFile
                   + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("run=")) {
        runType = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("plinkPrefix=")) {
        plinkPrefix = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("exclude=")) {
        excludeFile = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("keep=")) {
        keepFile = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("parse=")) {
        dir = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("pThreshold=")) {
        pValueThreshold = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("hwe=")) {
        hweFile = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("frq=")) {
        frqFile = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("resultPrefix=")) {
        resultPrefix = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("riskAlleles=")) {
        riskAlleleFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("model=")) {
        String modelString = ext.parseStringArg(arg, null);
        model = null;
        for (EMIM_MODEL m : EMIM_MODEL.values()) {
          if (modelString.equals(m.toString())) {
            model = m;
            break;
          }
        }
        if (model == null) {
          System.err.println("Invalid model: " + modelString);
        } else {
          numArgs--;
        }
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      if (riskAlleleFile != null) {
        generateRiskAlleles(riskAlleleFile, new Logger("GenerateRiskAlleles.log"));
      } else if (dir != null) {
        parse(dir, resultPrefix, hweFile, frqFile, pValueThreshold, model);
      } else if (plinkPrefix != null) {
        scriptAll(plinkPrefix, excludeFile, keepFile, pValueThreshold);
      } else {
        setTo(runType, model);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
