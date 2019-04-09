package org.genvisis.cnv.affy;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import org.genvisis.cnv.Resources;
import org.genvisis.cnv.Resources.Resource;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.Markers;
import org.genvisis.seq.GenomeBuild;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CmdLine;
import org.pankratzlab.common.Command;
import org.pankratzlab.common.Elision;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ext;

public class AxiomPipeline {

  public static final String CEL_EXTENSION = ".cel";
  public static final String CEL_GZ_EXTENSION = ".cel.gz";
  private static final String CEL_LIST_HEADER = "cel_files";
  private static final String AFFY_PROBELIST_HEADER = "probeset_id";

  private final String axiomXMLDefFile; // e.g. Axiom_tx_v1.r5.apt-probeset-genotype.AxiomSS1.apt2.xml
  private final String cdfFile;
  private final String aptExeDir;// holds "apt-geno-qc", "apt-probeset-genotype",
                                 // "apt-probeset-summarize" etc
  private final Logger log;

  public AxiomPipeline(String xmlDefsFile, String aptExeDir, Logger log) {
    this.axiomXMLDefFile = xmlDefsFile;
    // TODO don't auto-discover this file unless unspecified, and warn if doing so
    this.cdfFile = Files.list(ext.verifyDirFormat(ext.parseDirectoryOfFile(xmlDefsFile)), "",
                              ".cdf", false, true)[0];
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
    if (!Files.exists(axiomXMLDefFile)) {
      log.reportError(axiomXMLDefFile + " could not be found.");
      throw new IllegalArgumentException(axiomXMLDefFile + " could not be found.");
    }

    for (int i = 0; i < AFFY_ANALYSIS_TYPES.values().length; i++) {
      String exeFile = aptExeDir + AFFY_ANALYSIS_TYPES.values()[i].getExe();
      if (!Files.exists(exeFile)) {
        log.reportError(exeFile + " did not exist");
        throw new IllegalArgumentException();
      } else {
        Files.chmod(exeFile, false);
      }
    }
  }

  /**
   * The basic analyses we perform
   */
  private static enum AFFY_ANALYSIS_TYPES {

    GENERATE_PROBE_LIST("apt-probeset-summarize"),
    GENOTYPE("apt-probeset-genotype"),
    NORMALIZE("apt-probeset-summarize");

    private String exe;

    private AFFY_ANALYSIS_TYPES(String exe) {
      this.exe = exe;
    }

    public String getExe() {
      if (Files.isWindows()) {
        return exe + ".exe";
      } else {
        return exe; // *nix systems are much superior
      }
    }

  }

  private String generateCelList(String[] celFiles, String outDir, String analysisName) {
    String celListFile = outDir + analysisName + ".celList.txt";
    String[] toWrite = new String[] {CEL_LIST_HEADER};
    toWrite = ArrayUtils.concatAll(toWrite, celFiles);
    Files.writeArray(toWrite, celListFile);
    return celListFile;
  }

  private static class Probesets {

    private final String snpOnlyFile;
    private final String allFile;
    private boolean fail;

    public Probesets(String snpOnlyFile, String allFile) {
      super();
      this.snpOnlyFile = snpOnlyFile;
      this.allFile = allFile;
    }

    public boolean isFail() {
      return fail;
    }

    public void setFail(boolean fail) {
      this.fail = fail;
    }

    public String getSnpOnlyFile() {
      return snpOnlyFile;
    }

    public String getAllFile() {
      return allFile;
    }

  }

  /**
   * run a simple command to extract all probesets, return a formatted file to use for downstream
   * analysis
   */
  private Probesets getAnalysisProbesetList(String celFile, String outDir, String analysisName,
                                            String markerPositionFile) {
    String smallCelList = outDir + "tmp" + ext.getTimestampForFilename();
    String currentAnalysis = analysisName + "_probelistGenerator";
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
    psetCommand.add(aptExeDir + AFFY_ANALYSIS_TYPES.GENERATE_PROBE_LIST.getExe());
    psetCommand.add("--xml-file");
    psetCommand.add(axiomXMLDefFile);
    psetCommand.add("--cdf-file");
    psetCommand.add(cdfFile);
    psetCommand.add("--analysis");
    psetCommand.add("quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true");
    psetCommand.add("--out-dir");
    psetCommand.add(outDir);
    psetCommand.add("--cel-files");
    psetCommand.add(smallCelList);
    psetCommand.add("--set-analysis-name");
    psetCommand.add(currentAnalysis);

    log.report(ext.getTime() + " Info - running a command to extract probeset ids: " + psetCommand);

    String probeResults = outDir + currentAnalysis + ".summary.txt";

    boolean progress = CmdLine.builder(log).build()
                              .run(Command.builder(psetCommand).necessaryInputFiles(celFile)
                                          .expectedOutputFiles(probeResults).build());

    log.reportTimeInfo("Parsing " + probeResults + " to obtain all probesIds");
    ArrayList<String> probesetIdsAll = new ArrayList<>(1800000);
    ArrayList<String> probesetIdsSNP = new ArrayList<>(925000);

    probesetIdsAll.add(AFFY_PROBELIST_HEADER);
    probesetIdsSNP.add(AFFY_PROBELIST_HEADER);
    String tmpMarkerSet = outDir + analysisName + "tmpMarkerSet.ser";
    Markers.orderMarkers(null, markerPositionFile, tmpMarkerSet, log);
    @SuppressWarnings("deprecation")
    MarkerSetInfo markerSet = org.genvisis.cnv.filesys.MarkerSet.load(tmpMarkerSet);
    String[] names = markerSet.getMarkerNames();

    HashMap<String, String> track = new HashMap<>();
    ArrayList<String> markersNotUsed = new ArrayList<>();

    for (String name : names) {
      track.put(name, name);
    }
    try {
      BufferedReader reader = Files.getAppropriateReader(probeResults);
      while (reader.ready()) {
        String[] line = reader.readLine().trim().split("\t");
        if (line[0].startsWith("CN_") || line[0].startsWith("SNP_")
            || line[0].startsWith("AFFX-SNP")) {
          String pId = line[0];
          if (!pId.endsWith("-B")) {
            String probeParsed = pId.replaceAll("-A", "");
            if (track.containsKey(probeParsed)) {
              probesetIdsAll.add(probeParsed);
              if (!pId.startsWith("CN_")) {
                probesetIdsSNP.add(probeParsed);
              }
            } else {
              markersNotUsed.add(probeParsed);
            }
          }
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
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

  private static class NormalizationResult {

    private final String quantNormFile;
    private boolean fail;

    public NormalizationResult(String quantNormFile) {
      super();
      this.quantNormFile = quantNormFile;
    }

    public boolean isFail() {
      return fail;
    }

    public void setFail(boolean fail) {
      this.fail = fail;
    }

    public String getQuantNormFile() {
      return quantNormFile;
    }

  }

  private NormalizationResult normalize(String celListFile, String pIDFile, String analysisName,
                                        String outDirRoot) {
    String outCurrent = outDirRoot + analysisName + "_Normalization/";
    new File(outCurrent).mkdirs();
    ArrayList<String> normalizeCommand = new ArrayList<>();
    normalizeCommand.add(aptExeDir + AFFY_ANALYSIS_TYPES.NORMALIZE.getExe());
    normalizeCommand.add("--xml-file");
    normalizeCommand.add(axiomXMLDefFile);
    normalizeCommand.add("--cdf-file");
    normalizeCommand.add(cdfFile);
    normalizeCommand.add("--analysis");
    normalizeCommand.add("quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true");
    normalizeCommand.add("--out-dir");
    normalizeCommand.add(outCurrent);
    normalizeCommand.add("--cel-files");
    normalizeCommand.add(celListFile);
    normalizeCommand.add("--set-analysis-name");
    normalizeCommand.add(analysisName);

    String quantNormFile = outCurrent + analysisName + ".summary.txt";
    String report = outCurrent + analysisName + ".report.txt";

    boolean progress = CmdLine.builder(log).build()
                              .run(Command.builder(normalizeCommand).dir(outCurrent)
                                          .expectedOutputFiles(quantNormFile, report).build());

    NormalizationResult normalizationResult = new NormalizationResult(quantNormFile);
    normalizationResult.setFail(!progress);
    return normalizationResult;
  }

  private static class GenotypeResult {

    private final String callFile;
    private final String confFile;
    private boolean failed;

    public boolean isFailed() {
      return failed;
    }

    public void setFailed(boolean failed) {
      this.failed = failed;
    }

    public GenotypeResult(String callFile, String confFile) {
      super();
      this.callFile = callFile;
      this.confFile = confFile;
    }

    public String getCallFile() {
      return callFile;
    }

    public String getConfFile() {
      return confFile;
    }

  }

  private GenotypeResult genotype(String celListFile, String pIDFile, String analysisName,
                                  String outDirRoot) {

    String outCurrent = outDirRoot + analysisName + "_Genotypes/";
    new File(outCurrent).mkdirs();
    ArrayList<String> genotypeCommand = new ArrayList<>();
    genotypeCommand.add(aptExeDir + AFFY_ANALYSIS_TYPES.GENOTYPE.getExe());
    genotypeCommand.add("--xml-file");
    genotypeCommand.add(axiomXMLDefFile);
    genotypeCommand.add("--cdf-file");
    genotypeCommand.add(cdfFile);
    genotypeCommand.add("--table-output");
    genotypeCommand.add("true");
    genotypeCommand.add("--set-gender-method");
    genotypeCommand.add("cn-probe-chrXY-ratio");
    genotypeCommand.add("--set-analysis-name");
    genotypeCommand.add(analysisName);
    genotypeCommand.add("-out-dir");
    genotypeCommand.add(outCurrent);
    genotypeCommand.add("--cel-files");
    genotypeCommand.add(celListFile);

    String callFile = outCurrent + analysisName + ".calls.txt";
    String confFile = outCurrent + analysisName + ".confidences.txt";
    String reportFile = outCurrent + analysisName + ".report.txt";

    GenotypeResult genotypeResult = new GenotypeResult(callFile, confFile);

    boolean progress = CmdLine.builder(log).build()
                              .run(Command.builder(genotypeCommand)
                                          .expectedOutputFiles(genotypeResult.getCallFile(),
                                                               genotypeResult.getConfFile(),
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

  private static void run(String xmlDefsFile, String markerPositions, Project proj,
                          String aptExeDir, int numThreads) throws Elision {
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

    if (markerPositions == null || !Files.exists(markerPositions)) {
      if (markerPositions != null) {
        log.reportError("Could not find marker position file " + markerPositions);
      }
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

    AxiomPipeline pipeline = new AxiomPipeline(xmlDefsFile, aptExeDir, log);
    Probesets probeSets = pipeline.getAnalysisProbesetList(celFiles[0],
                                                           proj.PROJECT_DIRECTORY.getValue(),
                                                           proj.PROJECT_NAME.getValue(),
                                                           markerPositions);
    if (probeSets.isFail()) {
      throw new Elision("Critical Error - Generating probe set failed!");
    }

    String celListFile = pipeline.generateCelList(celFiles, proj.PROJECT_DIRECTORY.getValue(),
                                                  proj.PROJECT_NAME.getValue());
    GenotypeResult genotypeResult = pipeline.genotype(celListFile, probeSets.getSnpOnlyFile(),
                                                      proj.PROJECT_NAME.getValue(),
                                                      proj.PROJECT_DIRECTORY.getValue());
    if (genotypeResult.isFailed()) {
      throw new Elision("Critical Error - Genotyping failed!");
    }

    NormalizationResult normalizationResult = pipeline.normalize(celListFile,
                                                                 probeSets.getAllFile(),
                                                                 proj.PROJECT_NAME.getValue(),
                                                                 proj.PROJECT_DIRECTORY.getValue());
    if (normalizationResult.isFail()) {
      throw new Elision("Critical Error - Normalization failed!");
    }

    AffyParsingPipeline parser = new AffyParsingPipeline();
    parser.setProject(proj);
    parser.setGenotypeCallFile(genotypeResult.getCallFile());
    parser.setConfidencesFile(genotypeResult.getConfFile());
    parser.setNormIntensitiesFile(normalizationResult.getQuantNormFile());
    parser.run();

  }

  public static final int DEFAULT_MARKER_BUFFER = 100;
  public static final int DEFAULT_MAX_WRITERS = 1000000;

  public static void main(String[] args) {
    String analysisName = "Genvisis_affy_pipeline";
    String cels = "~/Affy6/cels/";
    String aptExeDir = "~/apt-1.18.0-x86_64-intel-linux/bin/";
    String xmlDefsFile = "";
    String outDir = "~/Affy6_results/";
    String markerPositions = null;
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
                   + aptExeDir + " (default))\n"
                   + "   (7) full path to a file of marker positions (i.e. markerPositions="
                   + "\"~/resources/hg18.affy6.markerPostions\" (not the default))\n"
                   + "   (8) optional: number of threads (i.e. " + PSF.Ext.NUM_THREADS_COMMAND + "="
                   + numThreads + " (default))\n"
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
      } else if (arg.startsWith("markerPositions=")) {
        markerPositions = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("cels=")) {
        cels = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("proj=")) {
        projFile = ext.parseStringArg(arg);
        numArgs--;
      } else if (arg.startsWith("outDir=")) {
        outDir = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("aptExeDir=")) {
        aptExeDir = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("xmlDefsFile=")) {
        xmlDefsFile = ext.parseStringArg(arg, "");
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
        if (markerPositions == null && Files.exists(proj.MARKER_POSITION_FILENAME.getValue())) {
          proj.getLog()
              .reportTimeWarning("No MarkerPositions argument present; using project file at "
                                 + proj.MARKER_POSITION_FILENAME.getValue());
          markerPositions = proj.MARKER_POSITION_FILENAME.getValue();
        }
        run(xmlDefsFile, markerPositions, proj, aptExeDir, numThreads);
      } catch (Elision e1) {
        System.err.println(e1.getMessage());
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
