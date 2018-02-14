package org.genvisis.imputation;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import javax.annotation.Nullable;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.PlinkData;
import org.genvisis.cnv.manage.VCFData;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;
import org.genvisis.gwas.FurtherAnalysisQc;
import org.genvisis.gwas.Qc;
import com.google.common.base.Joiner;
import com.google.common.base.Predicates;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Sets;

/**
 * Export genotypes to either PLINK or VCF format and, depending on selected options, generate
 * scripts to run ShapeIt and Minimac programs externally. <br />
 * <br />
 * Imputation TODO's:<br />
 * <ul>
 * <li>Make content aware (i.e., if PLINK files already exist, skip chrs that exist (or overwrite,
 * or fail))<br />
 * </li>
 * <li>Multithread data file export<br />
 * </li>
 * <li>Use ScriptExecutor for ShapeIt/Minimac scripts</li>
 * </ul>
 */
public class ImputationPipeline {

  public static final String PROJ_ARG = "proj=";
  public static final String REF_ARG = "ref=";
  public static final String DROP_SAMPLES_ARG = "dropSamples=";
  public static final String KEEP_SAMPLES_ARG = "keepSamples=";
  public static final String DROP_MARKERS_ARG = "dropMarkers=";
  public static final String KEEP_MARKERS_ARG = "keepMarkers=";
  public static final String PLINK_DIR_ARG = "plinkDir=";
  public static final String PLINK_PREFIX_ARG = "plinkPrefix=";
  public static final String OUT_DIR_AND_ROOT_ARG = "outDirAndRoot=";
  public static final String OUT_DIR_ARG = "outDir=";
  public static final String HAPS_DIR_ARG = "hapsDir=";
  public static final String CHRS_ARG = "chrs=";
  public static final String RUN_TYPE_ARG = "type=";
  public static final String USE_GRC_ARG = "useGRC=";

  private final Project proj;
  private final Set<String> dropMarkers;
  private final Set<String> dropSamples;
  private final Set<String> keepMarkers;
  private final Set<String> keepSamples;
  private final Map<String, Marker> prepMarkers = new HashMap<String, Marker>();

  public ImputationPipeline(Project proj, String referenceFile, KeepDrops keepDrops) {
    this.proj = proj;
    ImputationPrep prep = new ImputationPrep(proj, referenceFile);
    Set<Marker> markers = prep.getMatchingMarkers();
    for (Marker m : markers) {
      prepMarkers.put(m.getName(), m);
    }
    dropSamples = generateSampleSet(proj, keepDrops.getDropSamplesFile());
    keepSamples = generateSampleSet(proj, keepDrops.getKeepSamplesFile());
    dropMarkers = generateMarkerSet(proj, keepDrops.getDropMarkersFile());
    keepMarkers = generateMarkerSet(proj, keepDrops.getKeepMarkersFile());
  }

  private static Set<String> generateSampleSet(Project proj, String fidiidFile) {
    if (!Files.exists(fidiidFile)) {
      if (fidiidFile != null) {
        proj.getLog().reportTimeWarning("File of FID IIDs doesn't exist: " + fidiidFile);
      }
      return null;
    }
    SampleData sampleData = proj.getSampleData(false);
    Set<String> fidiids = HashVec.loadFileToHashSet(fidiidFile, new int[] {0, 1}, "\t", false);
    Set<String> samples = fidiids.stream().map(sampleData::lookupDNA).filter(Predicates.notNull())
                                 .collect(ImmutableSet.toImmutableSet());
    if (samples.size() != fidiids.size()) {
      Set<String> missingIDs = fidiids.stream().filter(fidiid -> !sampleData.lookupContains(fidiid))
                                      .collect(ImmutableSet.toImmutableSet());
      proj.getLog()
          .reportError("Not all Samples in " + fidiidFile + " could be found in SampleData: "
                       + ext.listWithCommas(missingIDs));
    }
    return samples;
  }

  private static Set<String> generateMarkerSet(Project proj, String markerFile) {
    if (!Files.exists(markerFile)) {
      if (markerFile != null) {
        proj.getLog().reportTimeWarning("File of markers doesn't exist: " + markerFile);
      }
      return null;
    }
    Set<String> markers = HashVec.loadFileToHashSet(markerFile, false);
    Set<String> missingIDs = Sets.difference(markers,
                                             proj.getMarkerSet().getMarkerNameMap().keySet());
    if (!missingIDs.isEmpty()) {
      proj.getLog().reportError("Not all Markers in " + markerFile + " could be found in Project: "
                                + ext.listWithCommas(missingIDs));
      markers.removeAll(missingIDs);
    }
    return Collections.unmodifiableSet(markers);
  }

  private Set<String> getSamplesToExport() {
    Set<String> samplesToExport;
    if (keepSamples != null) {
      samplesToExport = Sets.newHashSet(keepSamples);
    } else {
      samplesToExport = Sets.newHashSet(proj.getSamples());
    }
    if (dropSamples != null) {
      samplesToExport.removeAll(dropSamples);
    }
    return Collections.unmodifiableSet(samplesToExport);
  }

  private Set<String> getMarkersToExport() {
    Set<String> markersToExport;
    if (keepMarkers != null && keepMarkers.size() > 0) {
      markersToExport = Sets.newHashSet(keepMarkers);
    } else {
      markersToExport = Sets.newHashSet(proj.getMarkerNames());
    }
    if (dropMarkers != null) {
      markersToExport.removeAll(dropMarkers);
    }
    markersToExport.retainAll(prepMarkers.keySet());
    return markersToExport;
  }

  private Set<String> getChrMarkers(int chr) {
    Set<String> markersToExport = getMarkersToExport();
    Set<String> chrMarkers = new HashSet<String>();
    for (String m : markersToExport) {
      if (prepMarkers.get(m).getChr() == (byte) chr) {
        chrMarkers.add(m);
      }
    }
    return chrMarkers;
  }

  private ArrayList<String> getMarkersSortedNoDupes(int chr) {
    Set<String> chrMarkers = getChrMarkers(chr);

    int[] pos = new int[chrMarkers.size()];
    String[] mkr = chrMarkers.toArray(new String[chrMarkers.size()]);
    for (int i = 0; i < mkr.length; i++) {
      pos[i] = prepMarkers.get(mkr[i]).getPosition();
    }

    int[] indices = Sort.getSortedIndices(pos);

    ArrayList<String> mkrs = new ArrayList<String>();
    for (int i = 0; i < indices.length; i++) {
      if (i == 0 || pos[indices[i]] != pos[indices[i - 1]]) { // skip if prev (in sorted array) was
                                                             // same position
        mkrs.add(mkr[indices[i]]);
      }
    }
    return mkrs;
  }

  public void exportToPlink(String plinkDirAndRoot, int[] chrs) {
    // TODO (??) Only alphanumeric characters in FID/IID
    (new File(ext.parseDirectoryOfFile(plinkDirAndRoot + ".bim"))).mkdirs();
    HashSet<String> toDrop = new HashSet<String>();
    if (keepSamples != null && keepSamples.size() > 0) {
      for (String s : proj.getSamples()) {
        if (keepSamples.contains(s) && !dropSamples.contains(s)) {
          continue;
        }
        toDrop.add(s);
      }
    } else {
      toDrop.addAll(dropSamples);
    }
    String[] writtenDNAs = PlinkData.createFamFile(proj, plinkDirAndRoot, toDrop,
                                                   PlinkData.ExportIDScheme.CONCAT_FID_TO_IID);
    if (writtenDNAs == null) {
      // TODO error
      return;
    }
    int[] indicesOfTargetSamplesInProj = PlinkData.getIndicesOfTargetSamplesInProj(proj,
                                                                                   writtenDNAs,
                                                                                   proj.getLog());

    String clusterFilterFileName = proj.CLUSTER_FILTER_COLLECTION_FILENAME.getValue();

    // TODO multi-thread
    for (int chr : chrs) {
      proj.getLog().report("Exporting chr" + chr);
      ArrayList<String> mkrs = getMarkersSortedNoDupes(chr);

      String[] targetMarkers = mkrs.toArray(new String[mkrs.size()]);
      int[] indicesOfTargetMarkersInProj = new int[targetMarkers.length];
      HashMap<String, Byte> chrsOfTargetMarkers = new HashMap<String, Byte>();
      HashMap<String, Integer> posOfTargetMarkers = new HashMap<String, Integer>();
      PlinkData.getIndicesOfTargetMarkers(proj, targetMarkers, indicesOfTargetMarkersInProj,
                                          chrsOfTargetMarkers, posOfTargetMarkers, proj.getLog());

      String dirAndRoot = plinkDirAndRoot + "_chr" + chr;
      boolean success = PlinkData.createBedFileSnpMajor10KperCycle(proj,
                                                                   ImmutableSet.copyOf(targetMarkers),
                                                                   chrsOfTargetMarkers,
                                                                   posOfTargetMarkers,
                                                                   indicesOfTargetSamplesInProj,
                                                                   proj.GC_THRESHOLD.getValue()
                                                                                    .floatValue(),
                                                                   clusterFilterFileName,
                                                                   dirAndRoot, proj.getLog());

      if (success) {
        PrintWriter refWriter = Files.getAppropriateWriter(dirAndRoot + "_alleles.ref");
        for (String s : targetMarkers) {
          refWriter.println(s + "\t" + prepMarkers.get(s).getRef());
        }
        refWriter.flush();
        refWriter.close();
        Files.copyFile(plinkDirAndRoot + ".fam", dirAndRoot + ".fam");
      }
    }
  }

  public void exportToVCF(String vcfDirAndRoot, int[] chrs, boolean useGRC) {
    String[] samplesToExport = getSamplesToExport().toArray(new String[0]);
    String[] markersToExport = getMarkersToExport().toArray(new String[0]);

    VCFData.exportGenvisisToVCF(proj, samplesToExport, markersToExport, true, useGRC, chrs,
                                vcfDirAndRoot);
  }

  protected static class ImputationPipeRunner {

    public static void runVCF(String projPropFile, int[] chrs, String refFile, KeepDrops keepDrops,
                              String vcfDirAndRoot, boolean useGRC) {
      ImputationPipeline ip = setupPipe(projPropFile, refFile, keepDrops);
      ip.exportToVCF(vcfDirAndRoot, chrs, useGRC);
    }

    public static void runPlink(String projPropFile, int[] chrs, String refFile,
                                KeepDrops keepDrops, String outputDirAndRoot) {
      ImputationPipeline ip = setupPipe(projPropFile, refFile, keepDrops);
      ip.exportToPlink(outputDirAndRoot, chrs);
    }

    public static void runShapeIt(String projPropFile, int[] chrs, String plinkFileDir,
                                  String plinkPrefix, String outDir) {
      Project proj = new Project(projPropFile);
      new ImputationImpl.ShapeIt(proj, chrs, plinkFileDir, plinkPrefix, outDir).createScripts();
    }

    public static void runMinimac(String projPropFile, int[] chrs, String hapsDir, String outDir) {
      Project proj = new Project(projPropFile);
      new ImputationImpl.MiniMac(proj, chrs, hapsDir, outDir).createScripts();
    }

    public static void runPlinkAndShapeIt(String projPropFile, int[] chrs, String refFile,
                                          KeepDrops keepDrops, String outputDir) {
      runPlink(projPropFile, chrs, refFile, keepDrops, outputDir + "/plink/plink");
      runShapeIt(projPropFile, chrs, outputDir + "/plink/", "plink_chr", outputDir + "/haps/");
    }

    public static void runPlinkShapeItAndMinimac(String projPropFile, int[] chrs, String refFile,
                                                 KeepDrops keepDrops, String outputDir) {
      runPlinkAndShapeIt(projPropFile, chrs, refFile, keepDrops, outputDir);
      runMinimac(projPropFile, chrs, outputDir + "/haps/", outputDir);
    }

    private static ImputationPipeline setupPipe(String projPropFile, String refFile,
                                                KeepDrops keepDrops) {
      Project proj = new Project(projPropFile);
      return new ImputationPipeline(proj, refFile, keepDrops);
    }

  }

  public static class KeepDrops {

    private final String dropSamplesFile;
    private final String keepSamplesFile;
    private final String dropMarkersFile;
    private final String keepMarkersFile;

    /**
     * @param dropSamplesFile file of samples to drop
     * @param keepSamplesFile file of samples to keep
     * @param dropMarkersFile file of markers to drop
     * @param keepMarkersFile file of markers to keep
     */
    public KeepDrops(@Nullable String dropSamplesFile, @Nullable String keepSamplesFile,
                     @Nullable String dropMarkersFile, @Nullable String keepMarkersFile) {
      super();
      this.dropSamplesFile = dropSamplesFile;
      this.keepSamplesFile = keepSamplesFile;
      this.dropMarkersFile = dropMarkersFile;
      this.keepMarkersFile = keepMarkersFile;
    }

    public KeepDrops(@Nullable String plinkDir, @Nullable String keepSamplesFile,
                     @Nullable String keepMarkersFile) {
      this(defaultDropSamples(plinkDir), keepSamplesFile, defaultDropMarkers(plinkDir),
           keepMarkersFile);
    }

    public KeepDrops(@Nullable String plinkDir) {
      this(plinkDir, null, null);
    }

    private static String defaultDropMarkers(String plinkDir) {
      if (plinkDir == null) return null;
      String dir = plinkDir + Qc.QC_SUBDIR + FurtherAnalysisQc.FURTHER_ANALYSIS_DIR;
      return dir + FurtherAnalysisQc.MARKER_QC_DROPS;
    }

    private static String defaultDropSamples(String plinkDir) {
      if (plinkDir == null) return null;
      String dir = plinkDir + Qc.QC_SUBDIR + FurtherAnalysisQc.FURTHER_ANALYSIS_DIR;
      return dir + FurtherAnalysisQc.SAMPLE_QC_DROPS;
    }

    public String getDropSamplesFile() {
      return dropSamplesFile;
    }

    public String getKeepSamplesFile() {
      return keepSamplesFile;
    }

    public String getDropMarkersFile() {
      return dropMarkersFile;
    }

    public String getKeepMarkersFile() {
      return keepMarkersFile;
    }

  }

  public enum IMPUTATION_PIPELINE_PATH {
    VCF_ONLY(REF_ARG, DROP_SAMPLES_ARG, KEEP_SAMPLES_ARG, DROP_MARKERS_ARG, KEEP_MARKERS_ARG,
             OUT_DIR_AND_ROOT_ARG, USE_GRC_ARG),
    PLINK_ONLY(REF_ARG, DROP_SAMPLES_ARG, KEEP_SAMPLES_ARG, DROP_MARKERS_ARG, KEEP_MARKERS_ARG,
               OUT_DIR_AND_ROOT_ARG),
    PLINK_SHAPEIT(REF_ARG, DROP_SAMPLES_ARG, KEEP_SAMPLES_ARG, DROP_MARKERS_ARG, KEEP_MARKERS_ARG,
                  OUT_DIR_ARG),
    PLINK_SHAPEIT_MINIMAC(REF_ARG, PLINK_DIR_ARG),
    SHAPEIT(PLINK_DIR_ARG, PLINK_PREFIX_ARG, OUT_DIR_ARG),
    MINIMAC(HAPS_DIR_ARG, OUT_DIR_ARG);

    private final Set<String> reqs;

    IMPUTATION_PIPELINE_PATH(String... reqs) {
      this.reqs = ImmutableSet.copyOf(reqs);
    }

    public Set<String> getReqs() {
      return reqs;
    }

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String projFile = null;
    String refFile = null;
    String dropSamples = null;
    String keepSamples = null;
    String dropMarkers = null;
    String keepMarkers = null;
    String plinkSubdir = null;
    String outDirAndRoot = null;
    String hapsDir = null;
    String outDir = null;
    String plinkPrefix = null;
    int[] chrs = null;
    boolean useGRC = true;
    IMPUTATION_PIPELINE_PATH path = null;

    String usage = "\n" + "org.genvisis.imputation.ImputationPipeline requires 3+ arguments\n"
                   + "   (1) Project properties filename (i.e. " + PROJ_ARG + projFile
                   + " (default))\n"
                   + "   (2) A comma-separated list of chromosomes to export, or null for all (may be omitted) (i.e. "
                   + CHRS_ARG + chrs + " (default))\n" + "   (3) Imputation Pipeline path (i.e. "
                   + RUN_TYPE_ARG + " one of "
                   + ArrayUtils.toStr(IMPUTATION_PIPELINE_PATH.values(), ", ") + "))\n"
                   + "   --------------------- \n"
                   + "   The following arguments may be necessary depending on chosen pipeline:\n"
                   + "   (a) Reference Panel / Site List file, with mkr, chr, pos, ref, and alt columns (i.e. "
                   + REF_ARG + refFile + " (default))\n" + "   (b) File of samples to drop (i.e. "
                   + DROP_SAMPLES_ARG + "dropSamples.txt (no default))\n"
                   + "   (c) File of samples to keep (i.e. " + KEEP_SAMPLES_ARG
                   + "keepSamples.txt (no default))\n" + "   (d) File of markers to drop (i.e. "
                   + DROP_MARKERS_ARG + "dropMarkers.txt (no default))\n"
                   + "   (e) File of markers to keep (i.e. " + KEEP_MARKERS_ARG
                   + "keepMarkers.txt (no default))\n"
                   + "   (f) Subdirectory in which to create PLINK files (i.e. " + PLINK_DIR_ARG
                   + plinkSubdir + " (default))\n" + "   (g) PLINK output prefix (i.e. "
                   + PLINK_PREFIX_ARG + plinkPrefix + " (default))\n"
                   + "   (h) Output directory and fileroot (i.e " + OUT_DIR_AND_ROOT_ARG
                   + outDirAndRoot + " (default))\n" + "   (i) Output directory (i.e " + OUT_DIR_ARG
                   + outDir + " (default))\n"
                   + "   (j) Export contigs as 'chr1' instead of '1' (i.e. " + USE_GRC_ARG + useGRC
                   + " (default))\n" + "   (k) Directory with output from ShapeIt (i.e. "
                   + HAPS_DIR_ARG + hapsDir + " (default))\n" + "   --------------------- \n"
                   + "   Additional pipeline argument requirements are as follows:\n";
    for (IMPUTATION_PIPELINE_PATH pathOption : IMPUTATION_PIPELINE_PATH.values()) {
      usage += "\t" + pathOption.toString() + ":\n";
      usage += "\t" + Joiner.on(", ").join(pathOption.getReqs()) + "\n";
    }

    for (int i = 0; i < args.length; i++) {
      if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h")
          || args[i].equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (args[i].startsWith(PROJ_ARG)) {
        projFile = ext.parseStringArg(args[i]);
        numArgs--;
      } else if (args[i].startsWith(REF_ARG)) {
        refFile = ext.parseStringArg(args[i]);
        numArgs--;
      } else if (args[i].startsWith(DROP_SAMPLES_ARG)) {
        dropSamples = ext.parseStringArg(args[i]);
        numArgs--;
      } else if (args[i].startsWith(KEEP_SAMPLES_ARG)) {
        keepSamples = ext.parseStringArg(args[i]);
        numArgs--;
      } else if (args[i].startsWith(DROP_MARKERS_ARG)) {
        dropMarkers = ext.parseStringArg(args[i]);
        numArgs--;
      } else if (args[i].startsWith(KEEP_MARKERS_ARG)) {
        keepMarkers = ext.parseStringArg(args[i]);
        numArgs--;
      } else if (args[i].startsWith(PLINK_DIR_ARG)) {
        plinkSubdir = ext.parseStringArg(args[i]);
        numArgs--;
      } else if (args[i].startsWith(PLINK_PREFIX_ARG)) {
        plinkPrefix = ext.parseStringArg(args[i]);
        numArgs--;
      } else if (args[i].startsWith(OUT_DIR_AND_ROOT_ARG)) {
        outDirAndRoot = ext.parseStringArg(args[i]);
        new File(ext.parseDirectoryOfFile(outDirAndRoot)).mkdirs();
        numArgs--;
      } else if (args[i].startsWith(OUT_DIR_ARG)) {
        outDir = ext.verifyDirFormat(ext.parseStringArg(args[i]));
        new File(outDir).mkdirs();
        numArgs--;
      } else if (args[i].startsWith(HAPS_DIR_ARG)) {
        hapsDir = ext.parseStringArg(args[i]);
        numArgs--;
      } else if (args[i].startsWith(CHRS_ARG)) {
        chrs = ext.parseIntArrayArg(args[i]);
        numArgs--;
      } else if (args[i].startsWith(RUN_TYPE_ARG)) {
        path = IMPUTATION_PIPELINE_PATH.valueOf(ext.parseStringArg(args[i]));
        numArgs--;
      } else if (args[i].startsWith(USE_GRC_ARG)) {
        useGRC = ext.parseBooleanArg(args[i]);
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + args[i]);
      }
    }

    if (numArgs != 0 || path == null) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      switch (path) {
        case VCF_ONLY:
          ImputationPipeRunner.runVCF(projFile, chrs, refFile,
                                      new KeepDrops(dropSamples, keepSamples, dropMarkers,
                                                    keepMarkers),
                                      outDirAndRoot, useGRC);
          break;
        case PLINK_ONLY:
          ImputationPipeRunner.runPlink(projFile, chrs, refFile,
                                        new KeepDrops(dropSamples, keepSamples, dropMarkers,
                                                      keepMarkers),
                                        outDirAndRoot);
          break;
        case PLINK_SHAPEIT:
          ImputationPipeRunner.runPlinkAndShapeIt(projFile, chrs, refFile,
                                                  new KeepDrops(dropSamples, keepSamples,
                                                                dropMarkers, keepMarkers),
                                                  outDir);
          break;
        case PLINK_SHAPEIT_MINIMAC:
          ImputationPipeRunner.runPlinkShapeItAndMinimac(projFile, chrs, refFile,
                                                         new KeepDrops(dropSamples, keepSamples,
                                                                       dropMarkers, keepMarkers),
                                                         outDir);
          break;
        case MINIMAC:
          ImputationPipeRunner.runMinimac(projFile, chrs, hapsDir, outDir);
          break;
        case SHAPEIT:
          ImputationPipeRunner.runShapeIt(projFile, chrs, plinkSubdir, plinkPrefix, outDir);
          break;
        default:
          System.err.println("Error - unrecognized imputation path: " + path);
          break;
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
