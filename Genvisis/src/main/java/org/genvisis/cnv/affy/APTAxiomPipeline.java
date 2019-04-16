package org.genvisis.cnv.affy;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

import org.genvisis.cnv.Resources;
import org.genvisis.cnv.Resources.Resource;
import org.genvisis.cnv.affy.APTAffy6Pipeline.Probesets;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.seq.GenomeBuild;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.CmdLine;
import org.pankratzlab.common.Command;
import org.pankratzlab.common.Elision;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ext;

public class APTAxiomPipeline {

  public static final String CEL_EXTENSION = ".cel";
  public static final String CEL_GZ_EXTENSION = ".cel.gz";
  private static final String CEL_LIST_HEADER = "cel_files";
  private static final String AFFY_PROBELIST_HEADER = "probeset_id";
  private static final String APT_GENOTYPE_AXIOM = "apt-genotype-axiom";

  private final String libraryFilePath; // contains .cdf file etc
  private final String axiomXMLDefFile;
  private final String aptExeDir;// holds "apt-genotype-axiom"
  private final Logger log;

  public APTAxiomPipeline(String libraryFilePath, String aptExeDir, Logger log) {
    this.libraryFilePath = libraryFilePath;
    this.axiomXMLDefFile = Resources.axiomTx(log).getAPTGenotypeAxiomXML().get();
    this.aptExeDir = aptExeDir;
    this.log = log;
    if (!new File(aptExeDir).exists()) {
      log.reportError(aptExeDir + " did not exist (Affy exe directory)");
      throw new IllegalArgumentException();
    }
    validatePreReq();
    log.reportTimeInfo("Validated Affy Power Tools initialization");
  }

  private void validatePreReq() {
    String exeFile = aptExeDir + APT_GENOTYPE_AXIOM;
    if (!Files.exists(exeFile)) {
      log.reportError(exeFile + " did not exist");
      throw new IllegalArgumentException();
    } else {
      Files.chmod(exeFile, false);
    }
  }

  private String generateCelList(String[] celFiles, String outDir, String analysisName) {
    String celListFile = outDir + analysisName + ".celList.txt";
    String[] toWrite = new String[] {CEL_LIST_HEADER};
    toWrite = ArrayUtils.concatAll(toWrite, celFiles);
    Files.writeArray(toWrite, celListFile);
    return celListFile;
  }

  /**
   * run a simple command to extract all probesets, return a formatted file to use for downstream
   * analysis
   */
  private Probesets getAnalysisProbesetList(String celFile, String outDir, String analysisName,
                                            Set<String> markerNames) {
    String smallCelList = outDir + "tmp" + ext.getTimestampForFilename();
    String currentAnalysis = analysisName + "_probelist";
    try {
      PrintWriter writer = Files.openAppropriateWriter(smallCelList);
      writer.println(CEL_LIST_HEADER);
      writer.println(celFile);
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to .cel list " + smallCelList);
      e.printStackTrace();
    }

    ArrayList<String> psetCommand = new ArrayList<>();
    psetCommand.add(aptExeDir + APT_GENOTYPE_AXIOM);
    psetCommand.add("--arg-file");
    psetCommand.add(axiomXMLDefFile);
    psetCommand.add("--analysis-files-path");
    psetCommand.add(libraryFilePath);
    psetCommand.add("--out-dir");
    psetCommand.add(outDir);
    psetCommand.add("--cel-files");
    psetCommand.add(smallCelList);
    psetCommand.add("--analysis-name");
    psetCommand.add(currentAnalysis);

    log.report(ext.getTime() + " Info - running a command to extract probeset ids: " + psetCommand);

    String probeResults = outDir + currentAnalysis + ".calls.txt";

    boolean progress = CmdLine.builder(log).build()
                              .run(Command.builder(psetCommand).necessaryInputFiles(celFile)
                                          .expectedOutputFiles(probeResults).build());

    log.reportTimeInfo("Parsing " + probeResults + " to obtain all probesIds");
    ArrayList<String> probesetIdsAll = new ArrayList<>(1800000);
    ArrayList<String> probesetIdsSNP = new ArrayList<>(925000);

    probesetIdsAll.add(AFFY_PROBELIST_HEADER);
    probesetIdsSNP.add(AFFY_PROBELIST_HEADER);

    ArrayList<String> markersNotUsed = new ArrayList<>();

    try {
      BufferedReader reader = Files.getAppropriateReader(probeResults);
      while (reader.ready()) {
        String[] line = reader.readLine().trim().split("\t");
        String pID = line[0];
        if (markerNames.contains(pID)) {
          probesetIdsAll.add(pID);
          if (!pID.startsWith("CN_")) {
            probesetIdsSNP.add(pID);
          }
        } else {
          markersNotUsed.add(pID);
        }
      }
      reader.close();
    } catch (

    FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + probeResults + "\" not found in current directory");
      return null;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + probeResults + "\"");
      return null;
    }

    log.reportTimeInfo("Detected " + (probesetIdsAll.size() - 1) + " total probesets with "
                       + (probesetIdsSNP.size() - 1) + " SNP_ probesets");
    if (markersNotUsed.size() > 0) {
      log.reportTimeInfo(markersNotUsed.size() + " markers where skipped");
      String pIdSkipFile = outDir + analysisName + ".probesetIdsSkipped.txt";
      Files.writeArray(ArrayUtils.toStringArray(markersNotUsed), pIdSkipFile);

    }
    String pIdAllFile = outDir + analysisName + ".probesetIdsAll.txt";
    String pIdSnpFile = outDir + analysisName + ".probesetIdsSNPS.txt";
    Files.writeArray(ArrayUtils.toStringArray(probesetIdsAll), pIdAllFile);
    Files.writeArray(ArrayUtils.toStringArray(probesetIdsSNP), pIdSnpFile);

    new File(smallCelList).delete();
    Probesets probesets = new Probesets(pIdSnpFile, pIdAllFile);
    probesets.setFail(!progress);
    return probesets;
  }

  private static class APTAxiomResult {

    private final String callFile;
    private final String confFile;
    private final String intensityFile;
    private boolean failed;

    public boolean isFailed() {
      return failed;
    }

    public void setFailed(boolean failed) {
      this.failed = failed;
    }

    public APTAxiomResult(String callFile, String confFile, String intensityFile) {
      super();
      this.callFile = callFile;
      this.confFile = confFile;
      this.intensityFile = intensityFile;
    }

    public String getCallFile() {
      return callFile;
    }

    public String getConfFile() {
      return confFile;
    }

    public String getIntensityFile() {
      return intensityFile;
    }

  }

  private APTAxiomResult runAPTAxiom(String celListFile, String pIDFile, String analysisName,
                                     String outDirRoot) {

    String outCurrent = outDirRoot + analysisName + "_Genotypes/";
    new File(outCurrent).mkdirs();
    ArrayList<String> aptCommand = new ArrayList<>();
    aptCommand.add(aptExeDir + APT_GENOTYPE_AXIOM);
    aptCommand.add("--arg-file");
    aptCommand.add(axiomXMLDefFile);
    aptCommand.add("--analysis-files-path");
    aptCommand.add(libraryFilePath);
    aptCommand.add("--summaries");
    aptCommand.add("--probeset-ids");
    aptCommand.add(pIDFile);
    aptCommand.add("--analysis-name");
    aptCommand.add(analysisName);
    aptCommand.add("-out-dir");
    aptCommand.add(outCurrent);
    aptCommand.add("--cel-files");
    aptCommand.add(celListFile);

    String callFile = outCurrent + analysisName + ".calls.txt";
    String confFile = outCurrent + analysisName + ".confidences.txt";
    String summaryFile = outCurrent + analysisName + ".summary.txt";
    String reportFile = outCurrent + analysisName + ".report.txt";

    APTAxiomResult aptResult = new APTAxiomResult(callFile, confFile, summaryFile);

    boolean progress = CmdLine.builder(log).build()
                              .run(Command.builder(aptCommand)
                                          .expectedOutputFiles(aptResult.getCallFile(),
                                                               aptResult.getConfFile(),
                                                               aptResult.getIntensityFile(),
                                                               reportFile)
                                          .dir(outCurrent).build());

    aptResult.setFailed(!progress);
    return aptResult;
  }

  private static void validateCelSelection(String[] celFiles, Logger log) {
    HashMap<String, String> uniq = new HashMap<>();
    boolean error = false;
    for (String celFile : celFiles) {
      if (uniq.containsKey(ext.removeDirectoryInfo(celFile))) {
        log.reportError(ext.removeDirectoryInfo(celFile)
                        + " was seen multiple times, perhaps across directories");
        error = true;
      } else {
        uniq.put(ext.removeDirectoryInfo(celFile), ext.removeDirectoryInfo(celFile));
      }
      if (celFile.length() != celFile.replaceAll(PSF.Regex.GREEDY_WHITESPACE, "").length()) {
        log.reportError(celFile + " contained whitespace");
      }
    }
    if (error) {
      log.reportError("Invalid cel file list, apt will not run");
      throw new IllegalArgumentException();
    }
  }

  public static void run(String libraryFilePath, Project proj, String aptExeDir,
                         int numThreads) throws Elision {
    Logger log = proj.getLog();
    if (!proj.SOURCE_FILENAME_EXTENSION.getValue().toLowerCase().equals(CEL_EXTENSION)
        && !proj.SOURCE_FILENAME_EXTENSION.getValue().toLowerCase().equals(CEL_GZ_EXTENSION)) {
      log.reportError("Source file extension should be \"" + CEL_EXTENSION + "\" or \""
                      + CEL_GZ_EXTENSION + "\" for Affymetrix projects.  Parsing may fail.");
    }
    String[] celFiles = Files.list(proj.SOURCE_DIRECTORY.getValue(), null,
                                   proj.SOURCE_FILENAME_EXTENSION.getValue(), true, true);
    GenomeBuild build = proj.GENOME_BUILD_VERSION.getValue();

    validateCelSelection(celFiles, log);
    log.reportTimeInfo("Found " + celFiles.length + " " + CEL_EXTENSION + " files to process");

    String markerPositions = proj.MARKER_POSITION_FILENAME.getValue();
    if (markerPositions == null || !Files.exists(markerPositions)) {
      if (markerPositions != null) {
        log.reportError("Could not find marker position file " + markerPositions);
      }
      // TODO FIXME integrate more axiom array subtypes
      Resource markerPos = Resources.axiomTx(log).genome(build).getMarkerPositions();
      if (!markerPos.isAvailable()) {
        throw new IllegalArgumentException("Affymetrix marker positions file for build "
                                           + build.getBuild()
                                           + " is not available.  Parsing cannot continue.");
      } else {
        log.reportTime("No marker positions file found - copying resource from "
                       + ext.parseDirectoryOfFile(markerPos.get()) + " to "
                       + proj.MARKER_POSITION_FILENAME.getValue());
        Files.copyFile(markerPos.get(), proj.MARKER_POSITION_FILENAME.getValue());
        markerPositions = proj.MARKER_POSITION_FILENAME.getValue();
      }
    }
    Set<String> markers = HashVec.loadFileToHashSet(markerPositions, new int[] {0}, "", true);

    APTAxiomPipeline pipeline = new APTAxiomPipeline(libraryFilePath, aptExeDir, log);
    Probesets probeSets = pipeline.getAnalysisProbesetList(celFiles[0],
                                                           proj.PROJECT_DIRECTORY.getValue(),
                                                           proj.PROJECT_NAME.getValue(), markers);
    if (probeSets.isFail()) {
      throw new Elision("Critical Error - Generating probe set failed!");
    }

    String celListFile = pipeline.generateCelList(celFiles, proj.PROJECT_DIRECTORY.getValue(),
                                                  proj.PROJECT_NAME.getValue());
    APTAxiomResult aptResult = pipeline.runAPTAxiom(celListFile, probeSets.getAllFile(),
                                                    proj.PROJECT_NAME.getValue(),
                                                    proj.PROJECT_DIRECTORY.getValue());
    if (aptResult.isFailed()) {
      throw new Elision("Critical Error - Genotyping failed!");
    }

    AffyParsingPipeline parser = new AffyParsingPipeline();
    parser.setProject(proj);
    parser.setGenotypeCallFile(aptResult.getCallFile());
    parser.setConfidencesFile(aptResult.getConfFile());
    parser.setNormIntensitiesFile(aptResult.getIntensityFile());
    parser.run();

  }

  public static final int DEFAULT_MARKER_BUFFER = 100;
  public static final int DEFAULT_MAX_WRITERS = 1000000;

  public static void main(String[] args) {
    CLI cli = new CLI(APTAxiomPipeline.class);

    String DESC_CEL_PATH = "A directory or full path to a file containing " + CEL_EXTENSION
                           + " files for analysis";
    String DESC_LIB_PATH = "A directory with AffyPowerTools executables (should contain apt-genotype-axiom. Available at http://www.affymetrix.com/)";
    String DESC_EXE_PATH = "A directory with Affymetrix Library files (should contain a .cdf file, a .sketch file, etc. Available at http://www.affymetrix.com/)";
    String ARG_CEL_PATH = "cels";
    String ARG_EXE_PATH = "aptExeDir";
    String ARG_LIB_PATH = "libraryFilePath";

    cli.addArg(CLI.ARG_PROJ, CLI.DESC_PROJ, true);
    cli.addArg(ARG_CEL_PATH, DESC_CEL_PATH, true);
    cli.addArg(ARG_EXE_PATH, DESC_EXE_PATH, true);
    cli.addArg(ARG_LIB_PATH, DESC_LIB_PATH, true);
    cli.addArg(CLI.ARG_THREADS, CLI.DESC_THREADS, "1", false);

    cli.parseWithExit(args);

    try {
      Project proj = new Project(cli.get(CLI.ARG_PROJ));
      try {
        run(cli.get(ARG_LIB_PATH), proj, cli.get(ARG_EXE_PATH), cli.getI(CLI.ARG_THREADS));
      } catch (Elision e1) {
        System.err.println(e1.getMessage());
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
