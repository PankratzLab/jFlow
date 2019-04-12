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
import org.genvisis.cnv.affy.AffyPipeline.Probesets;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.seq.GenomeBuild;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CmdLine;
import org.pankratzlab.common.Command;
import org.pankratzlab.common.Elision;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ext;

public class AxiomPipeline {

  public static final String CEL_EXTENSION = ".cel";
  public static final String CEL_GZ_EXTENSION = ".cel.gz";
  private static final String CEL_LIST_HEADER = "cel_files";
  private static final String AFFY_PROBELIST_HEADER = "probeset_id";
  private static final String APT_GENOTYPE_AXIOM = "apt-genotype-axiom";

  private final String libraryFilePath; // contains .cdf file etc
  private final String axiomXMLDefFile;
  private final String aptExeDir;// holds "apt-genotype-axiom"
  private final Logger log;

  public AxiomPipeline(String libraryFilePath, String aptExeDir, Logger log) {
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
    ArrayList<String> genotypeCommand = new ArrayList<>();
    genotypeCommand.add(aptExeDir + APT_GENOTYPE_AXIOM);
    genotypeCommand.add("--arg-file");
    genotypeCommand.add(axiomXMLDefFile);
    genotypeCommand.add("--analysis-files-path");
    genotypeCommand.add(libraryFilePath);
    genotypeCommand.add("--summaries");
    genotypeCommand.add("--probeset-ids");
    genotypeCommand.add(pIDFile);
    genotypeCommand.add("--analysis-name");
    genotypeCommand.add(analysisName);
    genotypeCommand.add("-out-dir");
    genotypeCommand.add(outCurrent);
    genotypeCommand.add("--cel-files");
    genotypeCommand.add(celListFile);

    String callFile = outCurrent + analysisName + ".calls.txt";
    String confFile = outCurrent + analysisName + ".confidences.txt";
    String summaryFile = outCurrent + analysisName + ".summary.txt";
    String reportFile = outCurrent + analysisName + ".report.txt";

    APTAxiomResult genotypeResult = new APTAxiomResult(callFile, confFile, summaryFile);

    boolean progress = CmdLine.builder(log).build()
                              .run(Command.builder(genotypeCommand)
                                          .expectedOutputFiles(genotypeResult.getCallFile(),
                                                               genotypeResult.getConfFile(),
                                                               genotypeResult.getIntensityFile(),
                                                               reportFile)
                                          .dir(outCurrent).build());

    genotypeResult.setFailed(!progress);
    return genotypeResult;
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

  private static void run(String libraryFilePath, Project proj, String aptExeDir,
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

    AxiomPipeline pipeline = new AxiomPipeline(libraryFilePath, aptExeDir, log);
    Probesets probeSets = pipeline.getAnalysisProbesetList(celFiles[0],
                                                           proj.PROJECT_DIRECTORY.getValue(),
                                                           proj.PROJECT_NAME.getValue(), markers);
    if (probeSets.isFail()) {
      throw new Elision("Critical Error - Generating probe set failed!");
    }

    String celListFile = pipeline.generateCelList(celFiles, proj.PROJECT_DIRECTORY.getValue(),
                                                  proj.PROJECT_NAME.getValue());
    APTAxiomResult genotypeResult = pipeline.runAPTAxiom(celListFile, probeSets.getSnpOnlyFile(),
                                                         proj.PROJECT_NAME.getValue(),
                                                         proj.PROJECT_DIRECTORY.getValue());
    if (genotypeResult.isFailed()) {
      throw new Elision("Critical Error - Genotyping failed!");
    }

    AffyParsingPipeline parser = new AffyParsingPipeline();
    parser.setProject(proj);
    parser.setGenotypeCallFile(genotypeResult.getCallFile());
    parser.setConfidencesFile(genotypeResult.getConfFile());
    parser.setNormIntensitiesFile(genotypeResult.getIntensityFile());
    parser.run();

  }

  public static final int DEFAULT_MARKER_BUFFER = 100;
  public static final int DEFAULT_MAX_WRITERS = 1000000;

  public static void main(String[] args) {
    String analysisName = "Genvisis_affy_pipeline";
    String cels = "~/Affy6/cels/";
    String aptExeDir = "~/apt-1.18.0-x86_64-intel-linux/bin/";
    String libraryFilePath = "";
    int numThreads = 1;
    int markerBuffer = DEFAULT_MARKER_BUFFER;
    int maxWritersOpen = DEFAULT_MAX_WRITERS;
    int numArgs = args.length;
    GenomeBuild build = GenomeBuild.HG18;
    String projFile = null;

    String usage = "\n" + "affy.AxiomPipeline requires 0-1 arguments\n"
                   + "   (1) analysis name (i.e. analysisName=" + analysisName + " (default))\n"
                   + "   (2) a directory or full path to a file containing " + CEL_EXTENSION
                   + " files for analysis (i.e. cels=" + cels + " (default))\n"
                   + "   (4) directory with Affy Power Tools executables (should contain apt-probeset-genotype, etc. Available at http://www.affymetrix.com/) (i.e. aptExeDir="
                   + aptExeDir + " (default))\n" + "   (8) optional: number of threads (i.e. "
                   + PSF.Ext.NUM_THREADS_COMMAND + "=" + numThreads + " (default))\n"
                   + "   (9) optional: number of markers to buffer when splitting files (i.e. markerBuffer="
                   + markerBuffer + " (default))\n"
                   + "   (10) optional: maximum number of writers to open, if this is less than the sample size parsing will slow drastically (i.e. maxWritersOpen="
                   + maxWritersOpen + " (default))\n"
                   + "   (11) optional: use the full affymetrix cdf, which contains more mitochondrial probesets (i.e. -full (not the default))\n"
                   + "   (12) specify the genome build to use - ensuring the build matches your marker positions (i.e. build="
                   + build + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("analysisName=")) {
        analysisName = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("cels=")) {
        cels = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("proj=")) {
        projFile = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith("aptExeDir=")) {
        aptExeDir = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("libraryFilePath=")) {
        libraryFilePath = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
        numThreads = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("build=")) {
        try {
          build = GenomeBuild.valueOf(ext.parseStringArg(arg, ""));
          numArgs--;
        } catch (IllegalArgumentException ile) {
          System.err.println("Invalid build " + ext.parseStringArg(arg, ""));
          System.err.println("Options Are: ");
          for (int j = 0; j < GenomeBuild.values().length; j++) {
            System.err.println(GenomeBuild.values()[j]);
          }
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
      Project proj = new Project(projFile);
      try {
        run(libraryFilePath, proj, aptExeDir, numThreads);
      } catch (Elision e1) {
        System.err.println(e1.getMessage());
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
