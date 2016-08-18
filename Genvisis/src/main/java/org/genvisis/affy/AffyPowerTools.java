package org.genvisis.affy;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;

public class AffyPowerTools {
  public static final String[] AFFY_QC_LIB_FILES =
      {"GenomeWideSNP_6.chrXprobes", "GenomeWideSNP_6.chrYprobes", "GenomeWideSNP_6.r2.qca",
       "GenomeWideSNP_6.r2.qcc", "GenomeWideSNP_6.cdf"};
  public static final String[] AFFY_GENOTYPE_SUMMARIZE_LIB_FILES =
      {"GenomeWideSNP_6.cdf", "GenomeWideSNP_6.birdseed-v2.models", "GenomeWideSNP_6.specialSNPs",
       "GenomeWideSNP_6.chrXprobes", "GenomeWideSNP_6.chrYprobes", "hapmap.quant-norm.normalization-target.txt"};
  public static final String[] AFFY_ANALYSIS_OUTPUTS_PREFIXES =
      {"probesGenotype_", "probesSummarize_", "probesSummarize_small"};
  public static final String[] AFFY_ANALYSIS_OUTPUTS_SUFFIXES =
      {".summary.txt", ".report.txt", ".calls.txt", ".confidences.txt", ".log"};
  public static final String[] AFFY_GENVISIS_OUTPUTS_SUFFIXES = {".GENAFFY"};
  public static final String[] AFFY_GENVISIS_OUTPUTS_PREFIXES = {"SNP_", "CN_"};
  public static final String[] AFFY_PROBELIST_HEADER = {"probeset_id"};
  public static final String[] AFFY_CEL_LIST_HEADER = {"cel_files"};
  public static final String[] AFFY_BATCH_NAMES =
      {"AFFYQC", "AFFYGenotype_", "AFFYSummarize_", "GEN_AFFY", "B"};
  public static final String[] AFFY_AUXILIARY_LISTS =
      {"celList.txt", "celList_small.txt", "gw6_probesets.txt"};
  public static final String[] AFFY_ANALYSIS_TYPES =
      {"apt-geno-qc", "apt-probeset-genotype", "apt-probeset-summarize"};
  public static final String AFFY_EXTENSION = ".CEL";
  public static final String[] PENN_CNV_OUTPUTS = {"gw6.genocluster", "gw6.lrr_baf.txt"};
  public static final String[] PENN_BATCH_NAMES = {"PennGenoClust", "PennLRRBAF"};
  public static final String AFFY_TABLES_CLASS = "affy.AffySNP6Tables";
  public static final String DEFAULT_CLASS_PATH = org.genvisis.common.PSF.Java.GENVISIS;
  public static final int DEFUALT_LINE_BUFFER = 500;

  public static void affyQC(String pbsDir, String dataDir, String aptExeDir, String affyLib,
                            String affyResultsDir, int numJobs, int totalMemory,
                            double walltimeRequestedInHours, Logger log) {
    String affyChrX = affyLib + AFFY_QC_LIB_FILES[0];
    String affyChrY = affyLib + AFFY_QC_LIB_FILES[1];
    String affyQCA = affyLib + AFFY_QC_LIB_FILES[2];
    String affyQCC = affyLib + AFFY_QC_LIB_FILES[3];
    String affyCDF = affyLib + AFFY_QC_LIB_FILES[4];
    String[] celFiles = getMatchedFiles(dataDir, log, AFFY_EXTENSION);
    log.report("Found " + celFiles.length + " " + AFFY_EXTENSION + " files");
    int step = writeLists(dataDir, affyResultsDir, celFiles, numJobs, log, true,
                          AFFY_CEL_LIST_HEADER[0], "cel");
    log.report("Which means the step for " + numJobs + " jobs would be " + step + " cel files");
    String qc = aptExeDir + AFFY_ANALYSIS_TYPES[0] + " -c " + affyCDF + " --qca-file " + affyQCA
                + " -qcc-file " + affyQCC + " --chrX-probes " + affyChrX + " --chrY-probes "
                + affyChrY + " --out-file " + affyResultsDir + "[%0].qc --cel-files "
                + affyResultsDir + "qc[%0].txt";
    String[] qcs = getCommands(qc, numJobs, "");
    writeCommands(qcs, pbsDir + AFFY_BATCH_NAMES[0] + ".PBS", numJobs, totalMemory,
                  walltimeRequestedInHours, log);
  }

  private static String[] getFiles(String dataDir) {
    String[] files;
    files = new File(dataDir).list(new FilenameFilter() {
      @Override
      public boolean accept(File file, String filename) {
        return file.length() > 1000;
      }
    });
    return files;
  }

  public static void genotypeAndSummarize(String pbsDir, String dataDir, String aptExeDir,
                                          String affyLib, String affyResultsDir, String celList,
                                          String markerFile, int numSNPJobs, int numCNJobs,
                                          int numSNPBatches, int numCNBatches, int totalMemory,
                                          double walltimeRequestedInHours, Logger log) {
    String affyCDF = affyLib + AFFY_GENOTYPE_SUMMARIZE_LIB_FILES[0];
    String affyBirdseedModel = affyLib + AFFY_GENOTYPE_SUMMARIZE_LIB_FILES[1];
    String affySpecialSnps = affyLib + AFFY_GENOTYPE_SUMMARIZE_LIB_FILES[2];
    String affyChrX = affyLib + AFFY_GENOTYPE_SUMMARIZE_LIB_FILES[3];
    String affyChrY = affyLib + AFFY_GENOTYPE_SUMMARIZE_LIB_FILES[4];
    // Quant target is from pennCNV, should maybe allow to specify;
    String affyHapMapQuant = affyLib + AFFY_GENOTYPE_SUMMARIZE_LIB_FILES[5];
    if (numSNPJobs > 1 || numCNJobs > 1) {
      log.reportError("Warning - copy number specific and SNP probesets are being separated for downstream analysis");
    }
    if (celList == null) {
      celList = affyResultsDir + AFFY_AUXILIARY_LISTS[0];
      String[] celFiles = getMatchedFiles(dataDir, log, AFFY_EXTENSION);
      writeCelList(celFiles, celList, log);
      log.report("Warning - a celfile list was not provided, generating one to use at " + celList
                 + " using " + celFiles.length + " " + AFFY_EXTENSION + " files");
    }
    if (markerFile == null) {
      markerFile = affyResultsDir + AFFY_AUXILIARY_LISTS[2];
      log.report("Warning - a marker file was not provided , generating one to use at " + markerFile
                 + " using all available probesets");
      markerFile =
          getFullProbesetList(dataDir, affyResultsDir, aptExeDir, affyCDF, markerFile, log);
      log.report("Info - a marker file was generated at " + markerFile);

    }

    ArrayList<ArrayList<String>> allProbesets = collectProbesets(markerFile, log);
    String genotypeCommand = aptExeDir + AFFY_ANALYSIS_TYPES[1] + " -c " + affyCDF
                             + " --table-output true -a birdseed-v2 --set-gender-method cn-probe-chrXY-ratio --read-models-birdseed "
                             + affyBirdseedModel + " --special-snps " + affySpecialSnps
                             + " -out-dir " + affyResultsDir + AFFY_ANALYSIS_OUTPUTS_PREFIXES[0]
                             + "[%0]/ --cel-files " + celList + " --chrX-probes " + affyChrX
                             + " --chrY-probes " + affyChrY + " --probeset-ids " + affyResultsDir
                             + "[%0].txt --set-analysis-name " + AFFY_ANALYSIS_OUTPUTS_PREFIXES[0]
                             + "[%0]";
    String summarizeCommand = aptExeDir + AFFY_ANALYSIS_TYPES[2] + " --cdf-file " + affyCDF
                              + " --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --target-sketch "
                              + affyHapMapQuant + " --out-dir " + affyResultsDir
                              + AFFY_ANALYSIS_OUTPUTS_PREFIXES[1] + "[%0]/ --cel-files " + celList
                              + " --probeset-ids " + affyResultsDir
                              + "[%0].txt --set-analysis-name " + AFFY_ANALYSIS_OUTPUTS_PREFIXES[1]
                              + "[%0]";
    writeProbeBatches(affyResultsDir, allProbesets.get(0), numSNPJobs, numSNPBatches, true,
                      AFFY_PROBELIST_HEADER[0], log, "SNP_");
    writeProbeBatches(affyResultsDir, allProbesets.get(1), numCNJobs, numCNBatches, true,
                      AFFY_PROBELIST_HEADER[0], log, "CN_");
    writeCommandBatches(pbsDir + AFFY_BATCH_NAMES[1], genotypeCommand, numSNPJobs, numSNPBatches,
                        AFFY_GENVISIS_OUTPUTS_PREFIXES[0], totalMemory, walltimeRequestedInHours,
                        log);
    writeCommandBatches(pbsDir + AFFY_BATCH_NAMES[2], summarizeCommand, numSNPJobs, numSNPBatches,
                        AFFY_GENVISIS_OUTPUTS_PREFIXES[0], totalMemory, walltimeRequestedInHours,
                        log);
    writeCommandBatches(pbsDir + AFFY_BATCH_NAMES[2], summarizeCommand, numCNJobs, numCNBatches,
                        AFFY_GENVISIS_OUTPUTS_PREFIXES[1], totalMemory, walltimeRequestedInHours,
                        log);

    String help =
        generateGenoSumLogHelp(pbsDir, dataDir, aptExeDir, affyLib, affyResultsDir, celList,
                               markerFile, numSNPJobs, numCNJobs, numSNPBatches, numCNBatches, log);
    log.report(help);
  }

  public static void pennCNVAffyPipeline(String pbsDir, String dataDir, String pennCNVAffybin,
                                         String pennCNVAffyLib, String affyResultsDir,
                                         String pennCnvOutputDir, String sexFile, String locFile,
                                         int numSNPJobs, int numCNJobs, int numSNPBatches,
                                         int numCNBatches, int totalMemory,
                                         double walltimeRequestedInHours, Logger log) {
    String clusterCommand = "perl " + pennCNVAffybin + "generate_affy_geno_cluster.pl "
                            + affyResultsDir + affyResultsDir + AFFY_ANALYSIS_OUTPUTS_PREFIXES[0]
                            + "[%0]" + AFFY_ANALYSIS_OUTPUTS_SUFFIXES[2] + " " + affyResultsDir
                            + AFFY_ANALYSIS_OUTPUTS_PREFIXES[0] + "[%0]"
                            + AFFY_ANALYSIS_OUTPUTS_SUFFIXES[3] + " " + affyResultsDir
                            + AFFY_ANALYSIS_OUTPUTS_PREFIXES[1] + "[%0]"
                            + AFFY_ANALYSIS_OUTPUTS_SUFFIXES[0] + " -locfile " + locFile
                            + " -sexfile " + sexFile + " -out " + pennCnvOutputDir + "[%0]"
                            + PENN_CNV_OUTPUTS[0];
    writeCommandBatches(pbsDir + PENN_BATCH_NAMES[0], clusterCommand, numCNJobs, numCNBatches,
                        AFFY_GENVISIS_OUTPUTS_PREFIXES[1], totalMemory, walltimeRequestedInHours,
                        log);
    writeCommandBatches(pbsDir + PENN_BATCH_NAMES[0], clusterCommand, numSNPJobs, numSNPBatches,
                        AFFY_GENVISIS_OUTPUTS_PREFIXES[0], totalMemory, walltimeRequestedInHours,
                        log);

    String lrrBafCalc = "perl " + pennCNVAffybin + "normalize_affy_geno_cluster.pl "
                        + pennCnvOutputDir + "[%0]" + PENN_CNV_OUTPUTS[0] + " " + affyResultsDir
                        + AFFY_ANALYSIS_OUTPUTS_PREFIXES[1] + "[%0]"
                        + AFFY_ANALYSIS_OUTPUTS_SUFFIXES[0] + " -locfile " + locFile + " -out "
                        + pennCnvOutputDir + "[%0]LRRBAF" + PENN_CNV_OUTPUTS[1];
    writeCommandBatches(pbsDir + PENN_BATCH_NAMES[1], lrrBafCalc, numCNJobs, numCNBatches,
                        AFFY_GENVISIS_OUTPUTS_PREFIXES[1], totalMemory, walltimeRequestedInHours,
                        log);
    writeCommandBatches(pbsDir + PENN_BATCH_NAMES[1], lrrBafCalc, numSNPJobs, numSNPBatches,
                        AFFY_GENVISIS_OUTPUTS_PREFIXES[0], totalMemory, walltimeRequestedInHours,
                        log);
    // String splitSignalFileStart = "perl " + pennCNVbin + "kcolumn.pl " + outDir +
    // affyGenoClusterFolderOut + jobID + "/" + genoLrrBaf + " split 2 -tab -head 3 --name
    // --start_split 5001 -out " + outDir + "Kcol/Kcol" + jobID + "/gw6_split --filenameout " +
    // outDir + "Kcol/Kcol" + jobID + "/" + pennSplitOut;

    // calls=" + affyResultsDir + AFFY_ANALYSIS_OUTPUTS_PREFIXES[0] + "[%0]" +
    // AFFY_ANALYSIS_OUTPUTS_SUFFIXES[2] + " conf=" + affyResultsDir +
    // AFFY_ANALYSIS_OUTPUTS_PREFIXES[0] + "[%0]" + AFFY_ANALYSIS_OUTPUTS_SUFFIXES[3] + " sig=" +
    // affyResultsDir + AFFY_ANALYSIS_OUTPUTS_PREFIXES[1] + "[%0]" +
    // AFFY_ANALYSIS_OUTPUTS_SUFFIXES[0]

  }

  public static void parseAffyTables(String pbsDir, String affyResultsDir, String projectSource,
                                     int numSNPJobs, int numCNJobs, int numSNPBatches,
                                     int numCNBatches, int totalMemory,
                                     double walltimeRequestedInHours, String javaLocation,
                                     String javaCp, int numLinesBuffer, Logger log) {
    if (javaLocation == null) {
      javaLocation = "java";
    }
    if (javaCp == null) {
      javaCp = DEFAULT_CLASS_PATH;
      log.report("Warning - a class path to GENVISIS was not provided, setting to default "
                 + DEFAULT_CLASS_PATH);

    }
    if (numLinesBuffer == 0) {
      numLinesBuffer = DEFUALT_LINE_BUFFER;
      log.report("Warning - a line buffer size was not provided , setting to default of "
                 + DEFUALT_LINE_BUFFER);
    }
    String affySNP6TablesSNPCommand = javaLocation + " -Xmx" + totalMemory + "m -cp " + javaCp + " "
                                      + AFFY_TABLES_CLASS + " calls=" + affyResultsDir
                                      + AFFY_ANALYSIS_OUTPUTS_PREFIXES[0] + "[%0]/"
                                      + AFFY_ANALYSIS_OUTPUTS_PREFIXES[0] + "[%0]"
                                      + AFFY_ANALYSIS_OUTPUTS_SUFFIXES[2] + " conf="
                                      + affyResultsDir + AFFY_ANALYSIS_OUTPUTS_PREFIXES[0] + "[%0]/"
                                      + AFFY_ANALYSIS_OUTPUTS_PREFIXES[0] + "[%0]"
                                      + AFFY_ANALYSIS_OUTPUTS_SUFFIXES[3] + " sig=" + affyResultsDir
                                      + AFFY_ANALYSIS_OUTPUTS_PREFIXES[1] + "[%0]/"
                                      + AFFY_ANALYSIS_OUTPUTS_PREFIXES[1] + "[%0]"
                                      + AFFY_ANALYSIS_OUTPUTS_SUFFIXES[0] + " out=" + projectSource
                                      + "[%0]/" + " numLinesBuffer=" + numLinesBuffer + " -"
                                      + AFFY_GENVISIS_OUTPUTS_PREFIXES[0];
    String affySNP6TablesCNCommand =
        javaLocation + " -Xmx" + totalMemory + "m -cp " + javaCp + " " + AFFY_TABLES_CLASS + " sig="
                                     + affyResultsDir + AFFY_ANALYSIS_OUTPUTS_PREFIXES[1] + "[%0]/"
                                     + AFFY_ANALYSIS_OUTPUTS_PREFIXES[1] + "[%0]"
                                     + AFFY_ANALYSIS_OUTPUTS_SUFFIXES[0] + " out=" + projectSource
                                     + "[%0]/" + " numLinesBuffer=" + numLinesBuffer + " -"
                                     + AFFY_GENVISIS_OUTPUTS_PREFIXES[1];
    writeCommandBatches(pbsDir + AFFY_BATCH_NAMES[3], affySNP6TablesSNPCommand, numSNPJobs,
                        numSNPBatches, AFFY_GENVISIS_OUTPUTS_PREFIXES[0], totalMemory,
                        walltimeRequestedInHours, log);
    writeCommandBatches(pbsDir + AFFY_BATCH_NAMES[3], affySNP6TablesCNCommand, numCNJobs,
                        numCNBatches, AFFY_GENVISIS_OUTPUTS_PREFIXES[1], totalMemory,
                        walltimeRequestedInHours, log);
    String help =
        generateAffyTableLogHelp(pbsDir, affyResultsDir, numSNPJobs, numCNJobs, numSNPBatches,
                                 numCNBatches, totalMemory, walltimeRequestedInHours, javaLocation,
                                 javaCp, numLinesBuffer, log);
    log.report(help);

  }

  public static void mergeFiles(String pbsDir, String affyResultsDir, String projectSource,
                                int numJobs, int totalMemory, double walltimeRequestedInHours,
                                String javaLocation, String javaCp, Logger log) {
    String mergeCommand =
        javaLocation + " -cp " + javaCp + " " + AFFY_TABLES_CLASS + "-merge affyResultsDir="
                          + affyResultsDir + " out=" + projectSource + " threads=" + numJobs;
    String[] mergeCommands = new String[1];
    mergeCommands[0] = mergeCommand;

  }

  private static String generateAffyTableLogHelp(String pbsDir, String affyResultsDir,
                                                 int numSNPJobs, int numCNJobs, int numSNPBatches,
                                                 int numCNBatches, int totalMemory,
                                                 double walltimeRequestedInHours,
                                                 String javaLocation, String javaCp,
                                                 int numLinesBuffer, Logger log) {
    String help = "\n\n\n";
    help += "A few tips to parse the results of genotyping and summarizing in " + affyResultsDir
            + "...\n";
    help += "If you did not split up the analysis you can simply run:\n";
    help += "(1) " + javaLocation + " -cp " + javaCp + " " + AFFY_TABLES_CLASS
            + " calls=YourCallsHere(.calls.txt) conf=YourConfidenceHere(.confidences.txt) sig=YourSignalHere(.summary.txt) out="
            + affyResultsDir + "SNP/ numLinesBuffer=" + numLinesBuffer + " -"
            + AFFY_GENVISIS_OUTPUTS_PREFIXES[0];
    help += "(2) " + javaLocation + "  -cp " + javaCp + " " + AFFY_TABLES_CLASS
            + " sig=YourSignalHere(.summary.txt) out=" + affyResultsDir + "CN/ numLinesBuffer="
            + numLinesBuffer + " -" + AFFY_GENVISIS_OUTPUTS_PREFIXES[1];
    help += "\nIf you have split up the analysis you can either:\n";
    help += "(1) Run each of the commands  found in this log \n";
    help += "(2) Submit the .PBS files to a compute cluster\n";
    help +=
        "Please note that the number of jobs and batches must equal the number of jobs and batches from -gtsum ";
    help += "The final steps for importing into GENVISIS are to:\n";
    help += "(1) create a project file\n";
    help += "(2) run " + javaLocation + " -cp " + javaCp + " " + AFFY_TABLES_CLASS
            + " -merge affyResultsDir=" + affyResultsDir + " proj=Your Project file \n";
    return help;
  }

  private static String generateGenoSumLogHelp(String pbsDir, String dataDir, String aptExeDir,
                                               String affyLib, String affyResultsDir,
                                               String celList, String markerFile, int numSNPJobs,
                                               int numCNJobs, int numSNPBatches, int numCNBatches,
                                               Logger log) {
    String affyCDF = affyLib + AFFY_GENOTYPE_SUMMARIZE_LIB_FILES[0];
    String affyBirdseedModel = affyLib + AFFY_GENOTYPE_SUMMARIZE_LIB_FILES[1];
    String affySpecialSnps = affyLib + AFFY_GENOTYPE_SUMMARIZE_LIB_FILES[2];
    String affyChrX = affyLib + AFFY_GENOTYPE_SUMMARIZE_LIB_FILES[3];
    String affyChrY = affyLib + AFFY_GENOTYPE_SUMMARIZE_LIB_FILES[4];
    String affyHapMapQuant = affyLib + AFFY_GENOTYPE_SUMMARIZE_LIB_FILES[5];

    String help = "\n\n\nA few tips to continue with genotyping and summarizing your .CEL files in "
                  + dataDir + "\n";
    help +=
        "If you would not actually like to break up the analysis, you can run from the command line...\n ";
    help += "\n(1) " + aptExeDir + AFFY_ANALYSIS_TYPES[1] + " -c " + affyCDF
            + " --table-output true -a birdseed-v2 --set-gender-method cn-probe-chrXY-ratio --read-models-birdseed "
            + affyBirdseedModel + " --special-snps " + affySpecialSnps + " -out-dir "
            + affyResultsDir + " --cel-files " + celList + " --chrX-probes " + affyChrX
            + " --chrY-probes " + affyChrY + " --probeset-ids " + markerFile
            + " --set-analysis-name " + AFFY_ANALYSIS_OUTPUTS_PREFIXES[0] + "\n";
    help += "\n(2) " + aptExeDir + AFFY_ANALYSIS_TYPES[2] + " --cdf-file " + affyCDF
            + " --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --target-sketch "
            + affyHapMapQuant + " --out-dir " + affyResultsDir + " --cel-files " + celList
            + " --probeset-ids " + markerFile + " --set-analysis-name "
            + AFFY_ANALYSIS_OUTPUTS_PREFIXES[1] + "\n";
    help += "\n\nIf you would like to split up the analysis you can either:\n";
    help += "(1) Run each of the commands  found in this log \n";
    help += "(2) Submit the .PBS files to a compute cluster\n";
    help += "\nAfter genotyping and summarizing please run the following for more information:";
    help += "java -cp \"" + org.genvisis.common.PSF.Java.GENVISIS
            + "\" affy.AffyPowerTools -parseTables  pbsDir=" + pbsDir + " affyResultsDir="
            + affyResultsDir + " numSNPJobs=" + numSNPJobs + " numCNJobs=" + numCNJobs
            + " numSNPBatches=" + numSNPBatches + " numCNBatches=" + numCNBatches
            + "projectSource=YourDesiredSourceDirectory";
    help +=
        "\n projectSource=YourDesiredSourceDirectory will be the location Genvisis will look to load data ";
    help += "\n*******Please note that the number of jobs and batches must be the same*******";
    return help;

  }

  private static String getFullProbesetList(String dataDir, String affyResultsDir, String aptExeDir,
                                            String affyCDF, String markerFile, Logger log) {
    String smallCelList = affyResultsDir + AFFY_AUXILIARY_LISTS[1];
    String fileToUse = getMatchedFiles(dataDir, log, AFFY_EXTENSION)[0];
    try {
      PrintWriter writer = new PrintWriter(new FileWriter(smallCelList));
      writer.println(AFFY_CEL_LIST_HEADER[0]);
      writer.println(fileToUse);
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to .cel list " + smallCelList);
      e.printStackTrace();
      System.exit(1);
    }
    String psetCommand = aptExeDir + AFFY_ANALYSIS_TYPES[2] + " --cdf-file " + affyCDF
                         + " --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --out-dir "
                         + affyResultsDir + " --cel-files " + smallCelList + " --set-analysis-name "
                         + AFFY_ANALYSIS_OUTPUTS_PREFIXES[2];
    // String psetCommand = aptExeDir + AFFY_ANALYSIS_TYPES[2] + " -a rma --cdf-file " + affyCDF + "
    // -o " + affyResultsDir + " --cel-files " + smallCelList;
    log.report(ext.getTime() + " Info - running a command to extract probeset ids: " + psetCommand);
    CmdLine.run(psetCommand, aptExeDir);
    String toExtractFile =
        getMatchedFiles(affyResultsDir, log,
                        AFFY_ANALYSIS_OUTPUTS_PREFIXES[2] + AFFY_ANALYSIS_OUTPUTS_SUFFIXES[0])[0];
    log.report(ext.getTime() + " Info - extracting probesets from " + toExtractFile);
    extractProbesets(affyResultsDir, markerFile, toExtractFile, log);
    log.report(ext.getTime() + " Info - cleaning up files...");
    deleteFile(smallCelList, log);
    // delete summarize .summary
    deleteFile(toExtractFile, log);
    // delete summarize log
    deleteFile(affyResultsDir + AFFY_ANALYSIS_TYPES[2] + AFFY_ANALYSIS_OUTPUTS_SUFFIXES[4], log);
    // delete summarize .report
    deleteFile(affyResultsDir + AFFY_ANALYSIS_OUTPUTS_PREFIXES[2]
               + AFFY_ANALYSIS_OUTPUTS_SUFFIXES[1], log);
    return markerFile;
  }

  private static void extractProbesets(String affyResultsDir, String markerFile,
                                       String toExtractFile, Logger log) {
    try {
      BufferedReader reader = Files.getAppropriateReader(toExtractFile);
      PrintWriter writer = Files.getAppropriateWriter(markerFile);
      writer.println(AFFY_PROBELIST_HEADER[0]);
      String[] line;
      do {
        line = reader.readLine().trim().split("\t", -1);
      } while (reader.ready() && !line[0].matches(AFFY_PROBELIST_HEADER[0]));
      while (reader.ready()) {
        line = reader.readLine().trim().split("\t", -1);
        if (line[0].startsWith("CN_") || line[0].startsWith("SNP_")
            || line[0].startsWith("AFFX-SNP")) {
          if (line[0].endsWith("-A")) {
            writer.println(line[0].replaceAll("-A", ""));
          } else if (line[0].endsWith("-B")) {
            continue;
          } else {
            writer.println(line[0]);
          }
        }
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error - the file " + toExtractFile + " was not found");
      System.exit(1);
    } catch (Exception e) {
      log.reportError("Error writing to celfile list " + markerFile);
      e.printStackTrace();
      System.exit(1);

    }
  }

  private static void deleteFile(String file, Logger log) {
    File afile = new File(file);
    if (afile.delete()) {
      log.report("Info - cleaning up and deleted file " + file);
    } else {
      log.report("Warning - could not delete file " + file + " ,please manually remove it...");
    }
  }

  private static String[] getMatchedFiles(String dir, Logger log, String Suffix) {
    String[] files = getFiles(dir);
    ArrayList<String> celFiles = new ArrayList<String>();
    for (String file : files) {
      if (file.endsWith(Suffix)) {
        celFiles.add(dir + file);
      }
    }
    if (celFiles.size() == 0) {
      log.reportError("Error - did not find any " + Suffix + " files in " + dir);
      System.exit(1);
    }
    return celFiles.toArray(new String[celFiles.size()]);
  }

  private static void writeCelList(String[] celFiles, String celListFileName, Logger log) {
    try {
      PrintWriter writer = new PrintWriter(new FileWriter(celListFileName));
      writer.println(AFFY_CEL_LIST_HEADER[0]);
      for (String celFile : celFiles) {
        writer.println(celFile);
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to celfile list " + celListFileName);
      e.printStackTrace();
      System.exit(1);

    }
  }

  private static String[][] writeCommandBatches(String batchName, String command, int numJobs,
                                                int numBatches, String probeType, int totalMemory,
                                                double walltimeRequestedInHours, Logger log) {
    String[][] allCommands = new String[numBatches][];
    for (int i = 0; i < numBatches; i++) {
      String[] commands = getCommands(command, numJobs, probeType + AFFY_BATCH_NAMES[4] + i + "_");
      writeCommands(commands, batchName + probeType + "_" + AFFY_BATCH_NAMES[4] + i + ".PBS",
                    numJobs, totalMemory, walltimeRequestedInHours, log);
      allCommands[i] = commands;
    }
    return allCommands;
  }

  private static String[] getCommands(String baseCommand, int numJobs, String replacer) {
    String[] commands = new String[numJobs];
    String trav = baseCommand;
    String[][] iterations = Matrix.toMatrix(Array.stringArraySequence(numJobs, ""));
    for (int i = 0; i < iterations.length; i++) {
      for (int j = 0; j < iterations[i].length; j++) {
        commands[i] = ext.replaceAllWith(trav, "[%" + j + "]", replacer + iterations[i][j]);
      }
    }
    return commands;
  }

  private static void writeProbeBatches(String affyResultsDir, ArrayList<String> probesets,
                                        int numJobs, int numBatches, boolean header, String head,
                                        Logger log, String probeType) {
    int step = (int) Math.ceil((double) probesets.size() / (double) numBatches);
    log.report("Found " + probesets.size() + " " + probeType.replace("_", "") + " probesets");
    log.report("Which means the step for " + numBatches + " batches would be " + step
               + " probesets");
    int tracker = 0;
    int batch = 0;
    int count = 0;
    ArrayList<String> probesetBatch = new ArrayList<String>();
    for (int i = 0; i < probesets.size(); i++) {
      probesetBatch.add(probesets.get(i));
      tracker++;
      if (tracker == step || i == probesets.size() - 1) {
        tracker = 0;
        count += probesetBatch.size();
        log.report("The lists for batch " + batch + " contains " + probesetBatch.size()
                   + " probesets");
        writeLists("", affyResultsDir, probesetBatch.toArray(new String[probesetBatch.size()]),
                   numJobs, log, header, head, probeType + "B" + batch + "_");
        batch++;
        probesetBatch = new ArrayList<String>();
      }
    }
    if (count != probesets.size()) {
      log.reportError("Error - not all " + probeType + "were added to the lists");
      System.exit(1);
    } else {
      log.report("Created " + (numJobs * numBatches) + " " + probeType.replace("_", "") + " lists");
    }

  }

  private static int writeLists(String dataDir, String affyResultsDir, String[] list, int numJobs,
                                Logger log, boolean header, String head, String listName) {
    PrintWriter writer;
    int step = (int) Math.ceil((double) list.length / (double) numJobs);
    new File(affyResultsDir).mkdir();
    for (int i = 0; i < numJobs; i++) {
      try {
        writer = new PrintWriter(new FileWriter(affyResultsDir + listName + (i + 1) + ".txt"));
        if (header) {
          writer.println(head);
        }
        for (int j = i * step; j < Math.min(list.length, (i + 1) * step); j++) {
          writer.println(dataDir + list[j]);
        }
        writer.close();
      } catch (Exception e) {
        log.reportError("Error writing to list" + listName + (i + 1) + ".txt");
        System.exit(1);
        e.printStackTrace();
      }
    }
    return step;
  }

  private static void writeCommands(String[] commands, String batchName, int numJobs,
                                    int totalMemory, double walltimeRequestedInHours, Logger log) {
    Files.qsubMultiple(batchName, commands, numJobs, -1, totalMemory, walltimeRequestedInHours);
    log.report("\n\n***Begin Analysis commands for " + batchName + "***\n\n");
    for (String command : commands) {
      log.report(command);
    }
    log.report("\n\n***End Analysis commands for " + batchName + "***");
  }

  private static ArrayList<ArrayList<String>> collectProbesets(String markerFile, Logger log) {
    ArrayList<ArrayList<String>> allProbesets = new ArrayList<ArrayList<String>>();
    ArrayList<String> cnProbesets = new ArrayList<String>();
    ArrayList<String> snpProbesets = new ArrayList<String>();
    try {
      BufferedReader reader = Files.getAppropriateReader(markerFile);
      int numProbesets = 0;
      String[] line;
      do {
        line = reader.readLine().trim().split("\t", -1);
      } while (reader.ready() && !line[0].matches(AFFY_PROBELIST_HEADER[0]));
      if (!reader.ready()) {
        log.reportError("Error - reached the end of the file without finding a line with the following token: "
                        + Array.toStr(AFFY_PROBELIST_HEADER));
        reader.close();
        System.exit(1);
      }
      while (reader.ready()) {
        line = reader.readLine().trim().split("\t", -1);
        numProbesets++;
        if (line[0].startsWith("CN_")) {
          cnProbesets.add(line[0]);
        } else {
          snpProbesets.add(line[0]);
        }
      }
      reader.close();
      log.report("Info - found a total of " + numProbesets + " probesets in " + markerFile);
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error - the marker position file " + markerFile + " was not found");
      System.exit(1);
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + markerFile + "\"");
      System.exit(1);
    }
    allProbesets.add(snpProbesets);
    allProbesets.add(cnProbesets);
    return allProbesets;
  }

  public static void main(String[] args) {
    System.out.println("affy");
    int numArgs = args.length;
    // affyQC( dataDir, aptExeDir, affyLib, affyResultsDir, numJobs, totalMemory,
    // walltimeRequestedInHours, createLists) {
    String dataDir = "C:/data/aptAutomate/";
    String markerFile = null;
    String aptExeDir = "testexe/";
    String affyLib = "doubletestlib/";
    String affyResultsDir = "C:/data/aptAutomate/testResults/";
    String pbsDir = "C:/data/aptAutomate/testResults/";
    String projectSource = "C:/data/aptAutomate/00src/";
    String filename = "";
    String celList = null;
    String logfile = null;
    String javaLocation = null;
    String javaCP = null;
    int numLinesBuffer = 0;
    int numJobs = 2;
    int numCNJobs = 2;
    int numSNPJobs = 3;
    int numSNPBatches = 3;
    int numCNBatches = 3;
    int totalMemory = 24000;
    int walltimeRequestedInHours = 1;
    boolean genotypeSummarize = false;
    boolean parseTables = false;
    boolean merge = false;

    // TODO
    String usage = "TODO";
    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        return;
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("dataDir=")) {
        dataDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("pbsDir=")) {
        pbsDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("markerFile=")) {
        markerFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("projectSource=")) {
        projectSource = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("aptExeDir=")) {
        aptExeDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("affyLib=")) {
        affyLib = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("affyResultsDir=")) {
        affyResultsDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("celList=")) {
        celList = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("logfile=")) {
        logfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("java=")) {
        javaLocation = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("javaCP=")) {
        javaCP = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("numJobs=")) {
        numJobs = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("numSNPJobs=")) {
        numSNPJobs = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("numCNJobs=")) {
        numCNJobs = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("numSNPBatches=")) {
        numSNPBatches = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("numCNBatches=")) {
        numCNBatches = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("totalMemoryInmb=")) {
        totalMemory = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("walltimeInHours=")) {
        walltimeRequestedInHours = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("numLines=")) {
        numLinesBuffer = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("-gtsum")) {
        genotypeSummarize = true;
        numArgs--;
      } else if (arg.startsWith("-parseTables")) {
        parseTables = true;
        numArgs--;
      } else if (arg.startsWith("-merge")) {
        merge = true;
        numArgs--;
      }

    }
    if (numArgs != 0) {
      System.err.println(usage);
      return;
    }
    try {
      if (logfile == null) {
        logfile = "APTlog.txt";
      }

      Logger log = new Logger(affyResultsDir + logfile);
      if (genotypeSummarize) {
        genotypeAndSummarize(pbsDir, dataDir, aptExeDir, affyLib, affyResultsDir, celList,
                             markerFile, numSNPJobs, numCNJobs, numSNPBatches, numCNBatches,
                             totalMemory, walltimeRequestedInHours, log);
      } else if (parseTables) {
        parseAffyTables(pbsDir, affyResultsDir, projectSource, numSNPJobs, numCNJobs, numSNPBatches,
                        numCNBatches, totalMemory, walltimeRequestedInHours, javaLocation, javaCP,
                        numLinesBuffer, log);
      } else if (merge) {
        if (filename.equals("")) {
          System.err.println("Error - must hava a valid project file to merge affy Results");
        } else {
          Project proj = new Project(filename, false);
          mergeFiles(pbsDir, affyResultsDir, proj.getProperty(proj.SOURCE_DIRECTORY), numJobs,
                     totalMemory, walltimeRequestedInHours, javaLocation, javaCP, log);
        }
      } else {
        affyQC(pbsDir, dataDir, aptExeDir, affyLib, affyResultsDir, numJobs, totalMemory,
               walltimeRequestedInHours, log);
      }
      Files.backup(logfile, affyResultsDir, affyResultsDir);
    } catch (Exception e) {
      e.printStackTrace();
    }

    String base = "/lustre/pankrat2/normDat/";
    // String intitialCelList = base+"lists/bigListFinal.txt";
    // String pbsOutDir = "/home/pankrat2/lanej/PBS/";
    String affyChunk = "cels";
    String batchChunk = "SNP_";
    String genvisisSource = "/lustre/pankrat2/normDat/output/ARICGenvisis/00src/";
    String javaClass = "affy.AffySNP6Tables";
    String project = "/home/pankrat2/lanej/projects/dbGaP_ARIC_11908.properties";

    int lineBuffer = 20;
    int memory = 2999;
    // double wallTime = 30.00;

    String affyQCfolderOut = "QCOut";
    String affyGenofolderOut = "genoTypeOut";
    String affySumfolderOut = "summarOut";
    String affyGenoClusterFolderOut = "clusterOut";
    // String affyDetectFolder = "indCELDetect";
    String affyChunkProbe = "probes";
    String quantNorm = "quant-norm.pm-only.med-polish.expr.summary.txt";

    String finalCelList = base + "lists/bigListFinal.txt";
    String lists = "/scratch/normDat/affyCelDetect/lists/";

    // String snpProbesetList = "AllGenoTypingSnps.txt";
    // String cnProbesetList = "AllCopyNumberProbes.txt";
    String park = "/home/pankrat2/lanej/" + org.genvisis.common.PSF.Java.GENVISIS + " -Xmx"
                  + memory / 1024 + "G ";
    String sexFile = lists + "file_sex.txt";
    String outDir = base + "output/";
    String pennCNVbin = base + "pennCNV/bin/";
    String pennCNVlib = base + "pennCNV/lib/";
    String affyCDF = affyLib + "GenomeWideSNP_6.cdf";
    String affyBirdseedModel = affyLib + "GenomeWideSNP_6.birdseed-v2.models";
    String affySpecialSnps = affyLib + "GenomeWideSNP_6.specialSNPs";
    String affyHapMapQuant = affyLib + "hapmap.quant-norm.normalization-target.txt";
    // String affyRefInput = affyLib + "GenomeWideSNP_6.hapmap270.na33.r1.a5.ref";
    String affyChrX = affyLib + "GenomeWideSNP_6.chrXprobes";
    String affyChrY = affyLib + "GenomeWideSNP_6.chrYprobes";
    String affyQCA = affyLib + "GenomeWideSNP_6.r2.qca";
    String affyQCC = affyLib + "GenomeWideSNP_6.r2.qcc";
    // String annoFile = affyLib + "GenomeWideSNP_6.na33.annot.db";
    String affyGenoQC = affyLib + "apt-geno-qc";
    String affyGenotype = affyLib + "apt-probeset-genotype";
    String affySummarize = affyLib + "apt-probeset-summarize";
    String affyChpToText = affyLib + "apt-chp-to-txt";
    // String qcOut = outDir + "qcOut.txt";
    // String birdseedReport = "birdseed-v2.report.txt";
    String birdseedCalls = "birdseed-v2.calls.txt";
    String birdseedConf = "birdseed-v2.confidences.txt";
    String locFile = pennCNVlib + "affygw6.hg18.pfb";
    String genoCluster = "gw6.genocluster";
    String genoLrrBaf = "gw6.lrr_baf.txt";
    String pennSplitOut = "signalListFile.txt";
    // String pennHmm = pennCNVlib + "affygw6.hmm";

    // String detect_cnv = pennCNVbin + "penncnv/";

    // QC params
    // double callRateCut = .95;

    String[] affyQCJobs = new String[numJobs];
    String[] affyGenoJobs = new String[numJobs];
    String[] affySumJobs = new String[numJobs];
    String[] affyGenoClusterJobs = new String[numJobs];
    // String[] affyCNClusterJobs = new String[numJobs];
    String[] affySNP6Tables = new String[numJobs];
    String[] LRRBAFJobs = new String[numJobs];
    String[] startKColumn = new String[numJobs];
    String[] stopKColumn = new String[numJobs];
    String[] detectCNVs = new String[numJobs];
    String[] chpToTxt = new String[numJobs];

    for (int j = 0; j < numJobs; j++) {
      String batchID = "";
      if (j < 10) {
        batchID = batchChunk + "B0" + j + "_";
      }
      if (j >= 10) {
        batchID = batchChunk + "B" + j + "_";
      }

      // String affyGenoPBS = pbsOutDir + "GT" + batchID + ".pbs";
      // String affySumPBS = pbsOutDir + "summarOut" + batchID + ".pbs";
      // String affyPennGenoClustPBS = pbsOutDir + "PenGenoClust" + batchID + ".pbs";
      // String affyChpToTxtPBS = pbsOutDir + "ChpToTxt" + batchID + ".pbs";
      // String affyPennCNClustPBS = pbsOutDir + "PenCNClust" + batchID + ".pbs";
      // String affyPennGenoLRRBAFPBS = pbsOutDir + "PennLRR_BAF" + batchID + ".pbs";
      // String affySNP6TablesPBS = pbsOutDir + "AS6T" + batchID + ".pbs";
      //
      // String affyQCPBS = pbsOutDir + "QC" + batchID + ".pbs";
      // String affyKColStartPBS = pbsOutDir + "KColStart" + batchID + ".pbs";
      // String affyKColStopPBS = pbsOutDir + "KColStop" + batchID + ".pbs";
      // String affyDetectPBS = pbsOutDir + "JLDetectB0" + j + ".pbs";

      for (int i = 0; i < numJobs; i++) {
        String jobID = "";
        String detectBatch = "";
        if (i < 10) {
          jobID = batchID + "0" + i;
          detectBatch = "0" + i;
        }
        if (i >= 10) {
          jobID = batchID + i;
          detectBatch = "" + i;
        }

        String qc = affyGenoQC + " -c " + affyCDF + " --qca-file " + affyQCA + " -qcc-file "
                    + affyQCC + " --chrX-probes " + affyChrX + " --chrY-probes " + affyChrY
                    + " --out-file " + outDir + affyQCfolderOut + "/" + jobID + ".qc --cel-files "
                    + lists + affyChunk + jobID;
        affyQCJobs[i] = qc;

        String genotypeCommand = affyGenotype + " -c " + affyCDF
                                 + " --cc-chp-output --table-output true -a birdseed-v2 --set-gender-method cn-probe-chrXY-ratio --read-models-birdseed "
                                 + affyBirdseedModel + " --special-snps " + affySpecialSnps
                                 + " -out-dir " + outDir + affyGenofolderOut + jobID
                                 + "/ --cel-files " + finalCelList + " --chrX-probes " + affyChrX
                                 + " --chrY-probes " + affyChrY + " --probeset-ids " + lists
                                 + affyChunkProbe + jobID;
        affyGenoJobs[i] = genotypeCommand;

        String aptChptoTxt =
            affyChpToText + " " + outDir + affyGenofolderOut + jobID + "/cc-chp*/*.chp -o " + outDir
                             + "ARICGenvisis/00src/" + affyGenofolderOut + jobID + "/cc-chp";
        chpToTxt[i] = aptChptoTxt;

        String summarizeCommand = affySummarize + " --cdf-file " + affyCDF
                                  + " --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --target-sketch "
                                  + affyHapMapQuant + " --out-dir " + outDir + affySumfolderOut
                                  + jobID + "/ --cel-files " + finalCelList + " --probeset-ids "
                                  + lists + affyChunkProbe + jobID;

        affySumJobs[i] = summarizeCommand;

        String affySNP6TablesCommand = "/soft/java/jdk1.7.0_45/bin/java -cp " + park + " "
                                       + javaClass + " proj=" + project + " calls=" + outDir
                                       + affyGenofolderOut + jobID + "/" + birdseedCalls + " conf="
                                       + outDir + affyGenofolderOut + jobID + "/" + birdseedConf
                                       + " sig=" + outDir + affySumfolderOut + jobID + "/"
                                       + quantNorm + " out=" + genvisisSource + jobID + "/ -"
                                       + batchChunk + " numLinesBuffer=" + lineBuffer;

        affySNP6Tables[i] = affySNP6TablesCommand;
        // Include sex file later
        String generate_affy_geno_cluster =
            "perl " + pennCNVbin + "generate_affy_geno_cluster.pl " + outDir + affyGenofolderOut
                                            + jobID + "/" + birdseedCalls + " " + outDir
                                            + affyGenofolderOut + jobID + "/" + birdseedConf + " "
                                            + outDir + affySumfolderOut + jobID + "/" + quantNorm
                                            + " " + batchChunk + " -locfile " + locFile
                                            + " -sexfile " + sexFile + " -out " + outDir
                                            + affyGenoClusterFolderOut + jobID + "/" + genoCluster;

        affyGenoClusterJobs[i] = generate_affy_geno_cluster;

        // String generate_affy_CN_cluster = "perl " + pennCNVbin+"generate_affy_CN_cluster.pl "+
        // outDir + affyGenofolderOut + jobID +"/" + birdseedCalls
        // + " " + outDir + affyGenofolderOut + jobID +"/" + birdseedConf + " "
        // + outDir + affySumfolderOut + jobID +"/" +quantNorm + " -locfile " +
        // locFile + " -sexfile "+ sexFile + " -out " + outDir + affyGenoClusterFolderOut + jobID +
        // "/" + genoCluster;
        //
        // affyCNClusterJobs[i] = generate_affy_CN_cluster;

        String lrrBafCalc = "perl " + pennCNVbin + "normalize_affy_geno_cluster.pl " + outDir
                            + affyGenoClusterFolderOut + jobID + "/" + genoCluster + " " + outDir
                            + affySumfolderOut + jobID + "/" + quantNorm + " -locfile " + locFile
                            + " -out " + outDir + affyGenoClusterFolderOut + jobID + "/"
                            + genoLrrBaf;

        LRRBAFJobs[i] = lrrBafCalc;

        String splitSignalFileStart = "perl " + pennCNVbin + "kcolumn.pl " + outDir
                                      + affyGenoClusterFolderOut + jobID + "/" + genoLrrBaf
                                      + " split 2 -tab -head 3 --name --start_split 5001 -out "
                                      + outDir + "Kcol/Kcol" + jobID + "/gw6_split --filenameout "
                                      + outDir + "Kcol/Kcol" + jobID + "/" + pennSplitOut;

        startKColumn[i] = splitSignalFileStart;

        String splitSignalFileStop = "perl " + pennCNVbin + "kcolumn.pl " + outDir
                                     + affyGenoClusterFolderOut + jobID + "/" + genoLrrBaf
                                     + " split 2 -tab -head 3 --name --end_split 5000 -out "
                                     + outDir + "Kcol/Kcol" + jobID + "/gw6_split --filenameout "
                                     + outDir + "Kcol/Kcol" + jobID + "/" + pennSplitOut;
        stopKColumn[i] = splitSignalFileStop;

        // String detectCNV = "perl " + detect_cnv + "detect_cnv.pl --test -hmm " + pennHmm + " -pfb
        // " + locFile + " -conf -gcmodel " + pennCNVlib + "affygw6.hg18.gcmodel -log " + outDir +
        // affyDetectFolder + detectBatch + "/gw6.log -out " + outDir + affyDetectFolder +
        // detectBatch + "/gw6.rawcnv" + " -list " + lists + "indDetect" + detectBatch;
        String detectCNV =
            "perl /home/pankrat2/shared/bin/penncnv/detect_cnv.pl --test -hmm /home/pankrat2/shared/aric_gw6/pennCNVAffy/gw6/lib/affygw6.hmm"
                           + " -pfb /home/pankrat2/shared/aric_gw6/pennCNVAffy/gw6/lib/affygw6.hg18.pfb -conf -gcmodel /home/pankrat2/shared/aric_gw6/pennCNVAffy/gw6/lib/affygw6.hg18.gcmodel -log /scratch/normDat/affyCelDetect/pennResult/"
                           + detectBatch + "gw6.log -out /scratch/normDat/affyCelDetect/pennResult/"
                           + detectBatch + "gw6.rawcnv" + " -list " + lists + "x" + detectBatch;

        detectCNVs[i] = detectCNV;
        // System.out.println(qc);
        // System.out.println(generate_affy_geno_cluster);
        // System.out.println(detectCNV);
      }

      // Files.qsubMultiple(affyQCPBS , affyQCJobs, numJobs, memory ,totalMemory, wallTime);
      // Files.qsubMultiple(affyGenoPBS , affyGenoJobs, numJobs, memory ,totalMemory, wallTime);
      // Files.qsubMultiple(affyChpToTxtPBS , chpToTxt, numJobs, memory ,totalMemory, wallTime);

      // Files.qsubMultiple(affySumPBS , affySumJobs, numJobs, memory ,totalMemory, wallTime);
      // Files.qsubMultiple(affySNP6TablesPBS , affySNP6Tables, numJobs, memory ,totalMemory,
      // wallTime);
      // Files.qsubMultiple(affyPennGenoClustPBS , affyGenoClusterJobs, numJobs, memory
      // ,totalMemory, wallTime);
      // Files.qsubMultiple(affyPennCNClustPBS , affyGenoClusterJobs, numJobs, memory ,totalMemory,
      // wallTime);

      // Files.qsubMultiple(affyPennGenoLRRBAFPBS , LRRBAFJobs, numJobs, memory ,totalMemory,
      // wallTime);

      // Files.qsubMultiple(affyKColStartPBS , startKColumn, numJobs, memory ,totalMemory,
      // wallTime);
      // Files.qsubMultiple(affyKColStopPBS , stopKColumn, numJobs, memory ,totalMemory, wallTime);
      // Files.qsubMultiple(affyDetectPBS, detectCNVs, numJobs, memory, totalMemory, wallTime);
    }
  }

}

// works
// String qcCommand = affyGenoQC + " -c " + affyCDF + " --qca-file " + affyQCA + " -qcc-file " +
// affyQCC + " --chrX-probes " + affyChrX +" --chrY-probes " + affyChrY + " --cel-files "
// + intitialCelList + " --out-file " + qcOut;
// String qcCommand = affyGenoQC + " -c " + affyCDF + " --qca-file " + affyQCA + " -qcc-file " +
// affyQCC + " --chrX-probes " + affyChrX +" --chrY-probes " + affyChrY + " --out-file
// $wdir/qcOut.qc" + " *.CEL";
// //works
// String genotypeCommand = affyGenotype + " --summaries --write-models -c " + affyCDF + " -a
// birdseed-v2 --read-models-birdseed " +
// affyBirdseedModel + " --special-snps " + affySpecialSnps + " -out-dir " +
// outDir + " --cel-files " + intitialCelList + " --chrX-probes " + affyChrX +" --chrY-probes " +
// affyChrY ;
// works
// String summarizeCommand = "apt-probeset-summarize --cdf-file " + affyCDF +
// " --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --target-sketch "+
// affyHapMapQuant + " --out-dir " + outDir + " --cel-files " + finalCelList;
// works
// String copyNumberWorkflow = "apt-copynumber-workflow -v 1 --cdf-file " + affyCDF +
// " --chrX-probes " + affyChrX +" --chrY-probes " + affyChrY + " --special-snps "
// + affySpecialSnps + " --annotation-file " + annoFile + " --reference-input " +
// affyRefInput + " -out-dir " + outDir + " --cel-files " + finalCelList + " --text-output true";
// works
// String generate_affy_geno_cluster = "perl " + affyPennCnv+"generate_affy_geno_cluster.pl "+
// birdseedCalls + " " + birdseedConf + " " + quantNorm + " -locfile " +
// locFile + " -sexfile " + sexFile + " -out " + genoCluster;
// works
// String lrrBafCalc = "perl " + affyPennCnv+"normalize_affy_geno_cluster.pl "+
// genoCluster + " " + quantNorm + " -locfile " + locFile +
// " -out " + genoLrrBaf;
// //works
// String splitSignalFileStart = "perl " + pennCNV+"kcolumn.pl "+
// genoLrrBaf + " split 2 -tab -head 3 --name_by_header -out gw6 --filenameout " +pennSplitOut ;
// //works
// String detectCNV = "perl " + pennCNV+"detect_cnv.pl --test -hmm " + pennHmm +
// " -pfb " + locFile + " -list " + pennSplitOut + " -log gw6.log -out gw6.rawcnv";
//

// System.out.println(qcCommand);
// CmdLine.run(qcCommand , outDir );
//
// genFinalCelList(qcOut , finalCelList , callRateCut);
// System.out.println(genotypeCommand);
// CmdLine.run(genotypeCommand , outDir);
//
// System.out.println(summarizeCommand);
// CmdLine.run(summarizeCommand , outDir);
//
// genSexFile(birdseedReport, sexFile);
// //System.out.println(copyNumberWorkflow);
// //CmdLine.run(copyNumberWorkflow , outDir);
// System.out.println(generate_affy_geno_cluster);
// CmdLine.run(generate_affy_geno_cluster , outDir);

// System.out.println(lrrBafCalc);
// CmdLine.run(lrrBafCalc , outDir);
// run command from output directory

// System.out.println(splitSignalFileStart);
// CmdLine.run(splitSignalFileStart , outDir);

// System.out.println(detectCNV);
// CmdLine.run(detectCNV , outDir);
