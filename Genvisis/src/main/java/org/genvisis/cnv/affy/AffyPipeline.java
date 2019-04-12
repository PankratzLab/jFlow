package org.genvisis.cnv.affy;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import org.genvisis.cnv.LaunchProperties;
import org.genvisis.cnv.LaunchProperties.DefaultLaunchKeys;
import org.genvisis.cnv.Resources;
import org.genvisis.cnv.Resources.Resource;
import org.genvisis.cnv.analysis.CentroidCompute;
import org.genvisis.cnv.analysis.CentroidCompute.CentroidBuilder;
import org.genvisis.cnv.filesys.Centroids;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.filesys.Project.SOURCE_FILE_DELIMITERS;
import org.genvisis.cnv.manage.Markers;
import org.genvisis.cnv.manage.SourceFileParser;
import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
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

/**
 * @author lane0212 Maybe the last affy raw data processor.
 *         <p>
 *         Takes us from .CEL -> to Genvisis. Requires affy power tools and related library files
 *         <p>
 *         Basically follows http://penncnv.openbioinformatics.org/en/latest/user-guide/affy/
 */
public class AffyPipeline {

  public static final String AFFY_CEL_EXTENSION = ".cel";
  public static final String AFFY_CEL_GZ_EXTENSION = ".cel.gz";
  private static final String AFFY_CEL_LIST_HEADER = "cel_files";
  private static final String AFFY_PROBELIST_HEADER = "probeset_id";

  private final String aptExeDir;// holds "apt-geno-qc", "apt-probeset-genotype",
                                 // "apt-probeset-summarize" etc
  private final String aptLibDir;
  private final boolean full;
  private final Logger log;

  public AffyPipeline(String aptExeDir, String aptLibDir, boolean full, Logger log) {
    this.aptExeDir = aptExeDir;
    this.aptLibDir = aptLibDir;
    this.full = full;
    this.log = log;
    if (!new File(aptExeDir).exists()) {
      log.reportError(aptExeDir + " did not exist (Affy exe directory)");
      throw new IllegalArgumentException();
    }
    if (!new File(aptLibDir).exists()) {
      log.reportError(aptLibDir + " did not exist (Affy lib directory)");
      throw new IllegalArgumentException();
    }
    validatePreReq();
    log.reportTimeInfo("Validated Affy Power Tools initialization");
  }

  private void validatePreReq() {

    for (int i = 0; i < AFFY_LIB_FILES.values().length; i++) {
      if (!Files.exists(aptLibDir + AFFY_LIB_FILES.values()[i].getLibFile(full))) {
        log.reportError(aptLibDir + AFFY_LIB_FILES.values()[i].getLibFile(full) + " did not exist");
        throw new IllegalArgumentException();
      }
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
   * Required library files for our pipeline
   */
  private static enum AFFY_LIB_FILES {
    GW6_CDF("GenomeWideSNP_6.cdf"),
    GW6_BIRDSEED_MODELS("GenomeWideSNP_6.birdseed-v2.models"),
    GW6_SPECIAL_SNPS("GenomeWideSNP_6.specialSNPs"),
    GW6_CHRX("GenomeWideSNP_6.chrXprobes"),
    GW6_CHRY("GenomeWideSNP_6.chrYprobes");

    private String libFile;

    private AFFY_LIB_FILES(String libFile) {
      this.libFile = libFile;
    }

    public String getLibFile(boolean full) {
      if (full) {
        switch (this) {
          case GW6_CDF:
          case GW6_SPECIAL_SNPS:
            return ext.addToRoot(libFile, ".Full");
          case GW6_CHRY:
          case GW6_BIRDSEED_MODELS:
          case GW6_CHRX:
          default:
            break;

        }
      }
      return libFile;

    }
  }

  /**
   * The basic analyses we perform
   */
  public static enum AFFY_ANALYSIS_TYPES {

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
    String[] toWrite = new String[] {AFFY_CEL_LIST_HEADER};
    toWrite = ArrayUtils.concatAll(toWrite, celFiles);
    Files.writeArray(toWrite, celListFile);
    return celListFile;
  }

  public static class Probesets {

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
      writer.println(AFFY_CEL_LIST_HEADER);
      writer.println(celFile);
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to .cel list " + smallCelList);
      e.printStackTrace();
    }

    ArrayList<String> psetCommand = new ArrayList<>();
    psetCommand.add(aptExeDir + AFFY_ANALYSIS_TYPES.GENERATE_PROBE_LIST.getExe());
    psetCommand.add("--cdf-file");
    psetCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_CDF.getLibFile(full));
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
                                        String outDirRoot, String targetSketch) {
    String outCurrent = outDirRoot + analysisName + "_Normalization/";
    new File(outCurrent).mkdirs();
    ArrayList<String> normalizeCommand = new ArrayList<>();
    normalizeCommand.add(aptExeDir + AFFY_ANALYSIS_TYPES.NORMALIZE.getExe());
    normalizeCommand.add("--cdf-file");
    normalizeCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_CDF.getLibFile(full));
    normalizeCommand.add("--analysis");
    normalizeCommand.add("quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true");
    normalizeCommand.add("--target-sketch");
    normalizeCommand.add(targetSketch);
    normalizeCommand.add("--out-dir");
    normalizeCommand.add(outCurrent);
    normalizeCommand.add("--cel-files");
    normalizeCommand.add(celListFile);
    normalizeCommand.add("--probeset-ids");
    normalizeCommand.add(pIDFile);
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
    genotypeCommand.add("-c");
    genotypeCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_CDF.getLibFile(full));
    genotypeCommand.add("--table-output");
    genotypeCommand.add("true");
    genotypeCommand.add("-a");
    genotypeCommand.add("birdseed-v2");
    genotypeCommand.add("--set-gender-method");
    genotypeCommand.add("cn-probe-chrXY-ratio");
    genotypeCommand.add("--read-models-birdseed");
    genotypeCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_BIRDSEED_MODELS.getLibFile(full));
    genotypeCommand.add("--special-snps");
    genotypeCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_SPECIAL_SNPS.getLibFile(full));
    genotypeCommand.add("--chrX-probes");
    genotypeCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_CHRX.getLibFile(full));
    genotypeCommand.add("--chrY-probes");
    genotypeCommand.add(aptLibDir + AFFY_LIB_FILES.GW6_CHRY.getLibFile(full));
    genotypeCommand.add("--probeset-ids");
    genotypeCommand.add(pIDFile);
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

  public static void run(Project proj, String aptExeDir, String aptLibDir, String quantNormTarget,
                         boolean full, String markerPositions, int numThreads) throws Elision {
    Logger log = proj.getLog();
    if (!proj.SOURCE_FILENAME_EXTENSION.getValue().toLowerCase().equals(AFFY_CEL_EXTENSION)
        && !proj.SOURCE_FILENAME_EXTENSION.getValue().toLowerCase().equals(AFFY_CEL_GZ_EXTENSION)) {
      log.reportError("Source file extension should be \"" + AFFY_CEL_EXTENSION + "\" or \""
                      + AFFY_CEL_GZ_EXTENSION + "\" for Affymetrix projects.  Parsing may fail.");
    }
    String[] celFiles = Files.list(proj.SOURCE_DIRECTORY.getValue(), null,
                                   proj.SOURCE_FILENAME_EXTENSION.getValue(), true, true);
    GenomeBuild build = proj.GENOME_BUILD_VERSION.getValue();

    runInternal(celFiles, log, quantNormTarget, full, markerPositions, build, proj, aptExeDir,
                aptLibDir, numThreads);

  }

  public static void run(String aptExeDir, String aptLibDir, String cels, String outDir,
                         String quantNormTarget, String analysisName, String markerPositions,
                         boolean full, int numThreads, GenomeBuild build,
                         boolean prepImputation) throws Elision {
    new File(outDir).mkdirs();
    Logger log = new Logger(outDir + "affyPipeline.log");
    String[] celFiles;
    if (Files.isDirectory(cels)) {
      celFiles = Files.list(cels, null, AFFY_CEL_EXTENSION, false, true);
    } else {
      celFiles = HashVec.loadFileToStringArray(cels, false, new int[] {0}, true);
    }
    Project proj = createProject(log, outDir, analysisName, build);

    runInternal(celFiles, log, quantNormTarget, full, markerPositions, build, proj, aptExeDir,
                aptLibDir, numThreads);

    if (proj.getSourceFileHeaders(true) == null || proj.getSourceFileHeaders(true).size() == 0
        || Files.exists(proj.PROJECT_DIRECTORY.getValue() + Project.HEADERS_FILENAME)) {
      new File(proj.PROJECT_DIRECTORY.getValue() + Project.HEADERS_FILENAME).delete();
    }
    proj = new Project(proj.getPropertyFilename());
    SourceFileParser.createFiles(proj, numThreads);
    TransposeData.transposeData(proj, 2000000000, false);
    CentroidCompute.computeAndDumpCentroids(proj, proj.CUSTOM_CENTROIDS_FILENAME.getValue(),
                                            new CentroidBuilder(), numThreads, 2);
    Centroids.recompute(proj, proj.CUSTOM_CENTROIDS_FILENAME.getValue(), false, numThreads);
    TransposeData.transposeData(proj, 2000000000, false);
    SampleData.createMinimalSampleData(proj);
    GenvisisWorkflow.setupImputationDefaults(proj);

  }

  private static Project createProject(Logger log, String outDir, String analysisName,
                                       GenomeBuild build) {
    String propFileDir = LaunchProperties.get(DefaultLaunchKeys.PROJECTS_DIR);
    log.reportTimeInfo("Generating Genvisis project properties file in " + propFileDir);

    String outDirFull = new File(outDir).getAbsolutePath();
    String projectFile = propFileDir + analysisName + ".properties";
    Files.write("PROJECT_DIRECTORY=" + outDirFull, projectFile);
    Project proj = new Project(projectFile);
    proj.PROJECT_NAME.setValue(analysisName);
    proj.PROJECT_DIRECTORY.setValue(outDirFull);
    proj.SOURCE_DIRECTORY.setValue(analysisName + "_00src");
    proj.PROJECT_NAME.setValue(analysisName);
    proj.XY_SCALE_FACTOR.setValue((double) 100);

    proj.ARRAY_TYPE.setValue(ARRAY.AFFY_GW6);

    proj.SOURCE_FILENAME_EXTENSION.setValue(".txt.gz");
    proj.SOURCE_FILE_DELIMITER.setValue(SOURCE_FILE_DELIMITERS.TAB);
    proj.ID_HEADER.setValue("[FILENAME_ROOT]");
    proj.LONG_FORMAT.setValue(false);
    proj.GENOME_BUILD_VERSION.setValue(build);
    proj.MARKER_DATA_DIRECTORY.getValue(true, false);
    return proj;
  }

  private static void runInternal(String[] celFiles, Logger log, String quantNormTarget,
                                  boolean full, String markerPositions, GenomeBuild build,
                                  Project proj, String aptExeDir, String aptLibDir,
                                  int numThreads) throws Elision {

    validateCelSelection(celFiles, log);
    log.reportTimeInfo("Found " + celFiles.length + " " + AFFY_CEL_EXTENSION + " files to process");

    if (quantNormTarget == null || !Files.exists(quantNormTarget)) {
      throw new IllegalArgumentException("A valid target sketch file is required, and available from http://www.openbioinformatics.org/penncnv/download/gw6.tar.gz");
    }
    if (full) {
      log.reportTimeInfo("Running with full affymetrix cdf");
    } else {
      log.reportTimeInfo("Running with default affymetrix cdf, use the \"-full\" command to use the full version");
    }

    if (markerPositions == null || !Files.exists(markerPositions)) {
      if (markerPositions != null) {
        log.reportError("Could not find marker position file " + markerPositions);
      }
      Resource markerPos = Resources.affy(log).genome(build).getMarkerPositions();
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

    AffyPipeline pipeline = new AffyPipeline(aptExeDir, aptLibDir, full, log);
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
                                                                 proj.PROJECT_DIRECTORY.getValue(),
                                                                 quantNormTarget);
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
    String targetSketch = "~/resources/hapmap.quant-norm.normalization-target.txt";
    String aptExeDir = "~/apt-1.18.0-x86_64-intel-linux/bin/";
    String aptLibDir = "~/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/";
    String outDir = "~/Affy6_results/";
    String markerPositions = null;
    int numThreads = 1;
    int markerBuffer = DEFAULT_MARKER_BUFFER;
    int maxWritersOpen = DEFAULT_MAX_WRITERS;
    int numArgs = args.length;
    boolean full = false;
    boolean prepImputation = false;
    GenomeBuild build = GenomeBuild.HG18;
    String projFile = null;

    String usage = "\n" + "affy.AffyPipeline requires 0-1 arguments\n"
                   + "   (1) analysis name (i.e. analysisName=" + analysisName + " (default))\n"
                   + "   (2) a directory or full path to a file containing " + AFFY_CEL_EXTENSION
                   + " files for analysis (i.e. cels=" + cels + " (default))\n"
                   + "   (3) a target sketch file (such as hapmap.quant-norm.normalization-target.txt Available at ) (i.e. sketch="
                   + targetSketch + " (default))\n"
                   + "   (4) directory with Affy Power Tools executables (should contain apt-probeset-genotype, etc. Available at http://www.affymetrix.com/) (i.e. aptExeDir="
                   + aptExeDir + " (default))\n"
                   + "   (5) directory with Affy Power Tools library files (should contain GenomeWideSNP_6.cdf, etc. Available at http://www.affymetrix.com/) (i.e. aptLibDir="
                   + aptLibDir + " (default))\n" + "   (6) output directory  (i.e. outDir=" + outDir
                   + " (default))\n"
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
      } else if (arg.startsWith("sketch=")) {
        targetSketch = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("aptExeDir=")) {
        aptExeDir = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("aptLibDir=")) {
        aptLibDir = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
        numThreads = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("-full")) {
        full = true;
        numArgs--;
      } else if (arg.startsWith("-prepImputation")) {
        prepImputation = true;
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
      if (projFile != null) {
        Project proj = new Project(projFile);
        try {
          if (markerPositions == null && Files.exists(proj.MARKER_POSITION_FILENAME.getValue())) {
            proj.getLog()
                .reportTimeWarning("No MarkerPositions argument present; using project file at "
                                   + proj.MARKER_POSITION_FILENAME.getValue());
            markerPositions = proj.MARKER_POSITION_FILENAME.getValue();
          }
          run(proj, aptExeDir, aptLibDir, targetSketch, full, markerPositions, numThreads);
        } catch (Elision e1) {
          System.err.println(e1.getMessage());
        }
      } else {
        run(aptExeDir, aptLibDir, cels, outDir, targetSketch, analysisName, markerPositions, full,
            numThreads, build, prepImputation);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
