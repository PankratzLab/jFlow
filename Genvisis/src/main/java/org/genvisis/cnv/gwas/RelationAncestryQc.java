package org.genvisis.cnv.gwas;

import java.io.File;
import java.util.Collections;
import java.util.Date;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CmdLine;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ext;
import org.pankratzlab.gwas.Plink;
import org.pankratzlab.common.CLI;
import com.google.common.collect.Maps;

public class RelationAncestryQc extends Qc {

  public static final String MARKER_QC_DIR = "marker_qc/";
  public static final String SAMPLE_QC_DIR = "sample_qc/";
  public static final String LD_PRUNING_DIR = "ld_pruning/";
  public static final String GENOME_DIR = "genome/";
  public static final String ANCESTRY_DIR = "ancestry/";

  public static final String UNRELATEDS_FILENAME = "unrelateds.txt";

  /** A rough listing of the Folders created by fullGamut */
  public static final String[] FOLDERS_CREATED = {Qc.QC_SUBDIR + MARKER_QC_DIR,
                                                  Qc.QC_SUBDIR + SAMPLE_QC_DIR,
                                                  Qc.QC_SUBDIR + LD_PRUNING_DIR,
                                                  Qc.QC_SUBDIR + GENOME_DIR,
                                                  Qc.QC_SUBDIR + ANCESTRY_DIR};
  /** A rough listing of the files created, by folder, by fullGamut */
  // TODO: This does not accommodate cases where the plinkroot is something other than
  // Qc.DEFAULT_PLINKROOT
  // Also ought to be automated...
  public static final String[][] FILES_CREATED = {{Qc.DEFAULT_PLINKROOT + ".bed", "freq.frq",
                                                   "missing.imiss",
                                                   /* "test.missing.missing", *//*
                                                                                 * not actually
                                                                                 * necessary
                                                                                 */ "hardy.hwe",
                                                   "mishap.missing.hap", "gender.assoc",
                                                   "gender.missing", MARKER_QC_DROPS},
                                                  {Qc.DEFAULT_PLINKROOT + ".bed", "missing.imiss"},
                                                  {Qc.DEFAULT_PLINKROOT + ".bed",
                                                   Qc.DEFAULT_PLINKROOT + ".prune.in"},
                                                  {Qc.DEFAULT_PLINKROOT + ".bed",
                                                   Qc.DEFAULT_PLINKROOT + ".genome",
                                                   Qc.DEFAULT_PLINKROOT + ".genome_keep.dat"},
                                                  {Qc.DEFAULT_PLINKROOT + ".bed",
                                                   UNRELATEDS_FILENAME}};
  public static final String ARGS_KEEPGENOME = "keepGenomeInfoForRelatedsOnly";

  public static final Map<QcMetric, String> DEFAULT_QC_METRIC_THRESHOLDS = MarkerQC.DEFAULT_METRIC_THRESHOLDS;
  public static final Set<QcMetric> CUSTOMIZABLE_QC_METRICS = Collections.unmodifiableSet(EnumSet.of(QcMetric.CALLRATE));

  /**
   * @see Qc#Qc(String, String, Map, Logger)
   */
  public RelationAncestryQc(String soureDir, String plinkPrefix,
                            Map<QcMetric, String> markerQCThresholds, Logger log) {
    super(soureDir, plinkPrefix, markerQCThresholds, log);
  }

  public void run(boolean keepGenomeInfoForRelatedsOnly) {
    long time = new Date().getTime();

    if (!markerQc(RelationAncestryQc.MARKER_QC_DIR)) {
      return;
    }

    new File(qcDir + RelationAncestryQc.SAMPLE_QC_DIR).mkdirs();
    if (!Files.exists(qcDir + RelationAncestryQc.SAMPLE_QC_DIR + plinkroot + ".bed")) {
      log.report(ext.getTime() + "]\tRunning --exclude " + MARKER_QC_DROPS);
      CmdLine.runDefaults("plink2 --bfile ../" + RelationAncestryQc.MARKER_QC_DIR + plinkroot
                          + " --exclude ../" + RelationAncestryQc.MARKER_QC_DIR + MARKER_QC_DROPS
                          + " --make-bed --noweb --out " + plinkroot,
                          qcDir + RelationAncestryQc.SAMPLE_QC_DIR, log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(qcDir + RelationAncestryQc.SAMPLE_QC_DIR + "missing.imiss")) {
      log.report(ext.getTime() + "]\tRunning --missing");
      CmdLine.runDefaults("plink2 --bfile " + plinkroot
                          + " --geno 1 --mind 1 --missing --out missing --noweb",
                          qcDir + RelationAncestryQc.SAMPLE_QC_DIR, log);
    }
    PSF.checkInterrupted();

    new File(qcDir + RelationAncestryQc.LD_PRUNING_DIR).mkdirs();
    if (!Files.exists(qcDir + RelationAncestryQc.LD_PRUNING_DIR + plinkroot + ".bed")) {
      log.report(ext.getTime()
                 + "]\tRunning --mind 0.05 (removes samples with callrate <95% for the markers that did pass QC)");
      CmdLine.runDefaults("plink2 --bfile ../" + RelationAncestryQc.SAMPLE_QC_DIR + plinkroot
                          + " --mind 0.05 --make-bed --noweb --out " + plinkroot,
                          qcDir + RelationAncestryQc.LD_PRUNING_DIR, log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(qcDir + RelationAncestryQc.LD_PRUNING_DIR + plinkroot + ".prune.in")) {
      log.report(ext.getTime() + "]\tRunning --indep-pairwise 50 5 0.3");
      CmdLine.runDefaults("plink2 --noweb --bfile " + plinkroot
                          + " --indep-pairwise 50 5 0.3 --out " + plinkroot,
                          qcDir + RelationAncestryQc.LD_PRUNING_DIR, log);
    }
    PSF.checkInterrupted();

    new File(qcDir + RelationAncestryQc.GENOME_DIR).mkdirs();
    if (!Files.exists(qcDir + RelationAncestryQc.GENOME_DIR + plinkroot + ".bed")) {
      log.report(ext.getTime() + "]\tRunning --extract " + plinkroot + ".prune.in");
      CmdLine.runDefaults("plink2 --bfile ../" + RelationAncestryQc.LD_PRUNING_DIR + plinkroot
                          + " --extract ../" + RelationAncestryQc.LD_PRUNING_DIR + plinkroot
                          + ".prune.in --make-bed --noweb --out " + plinkroot,
                          qcDir + RelationAncestryQc.GENOME_DIR, log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(qcDir + RelationAncestryQc.GENOME_DIR + plinkroot + ".genome")) {
      log.report(ext.getTime() + "]\tRunning --genome"
                 + (keepGenomeInfoForRelatedsOnly ? " --min 0.1" : ""));
      CmdLine.runDefaults("plink2 --noweb --bfile " + plinkroot + " --genome"
                          + (keepGenomeInfoForRelatedsOnly ? " --min 0.1" : "") + " --out "
                          + plinkroot, qcDir + RelationAncestryQc.GENOME_DIR, log);
    }
    PSF.checkInterrupted();
    if (!keepGenomeInfoForRelatedsOnly
        && !Files.exists(qcDir + RelationAncestryQc.GENOME_DIR + "mds20.mds")) {
      log.report(ext.getTime() + "]\tRunning --mds-plot 20");
      CmdLine.runDefaults("plink2 --bfile " + plinkroot + " --read-genome " + plinkroot
                          + ".genome --cluster --mds-plot 20 --out mds20 --noweb",
                          qcDir + RelationAncestryQc.GENOME_DIR, log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(qcDir + RelationAncestryQc.GENOME_DIR + plinkroot + ".genome_keep.dat")) {
      log.report(ext.getTime() + "]\tRunning flagRelateds");
      String lrrFile = qcDir + RelationAncestryQc.GENOME_DIR + "lrr_sd.xln";
      if (!Files.exists(lrrFile)) {
        lrrFile = qcDir + "../../lrr_sd.xln";
      }
      Plink.flagRelateds(qcDir + RelationAncestryQc.GENOME_DIR + plinkroot + ".genome",
                         qcDir + RelationAncestryQc.GENOME_DIR + plinkroot + ".fam",
                         qcDir + RelationAncestryQc.MARKER_QC_DIR + "missing.imiss", lrrFile,
                         Plink.FLAGS, Plink.THRESHOLDS, 4, false);
    }
    PSF.checkInterrupted();

    new File(qcDir + RelationAncestryQc.ANCESTRY_DIR).mkdirs();
    if (!Files.exists(qcDir + RelationAncestryQc.ANCESTRY_DIR + UNRELATEDS_FILENAME)) {
      log.report(ext.getTime() + "]\tCopying " + RelationAncestryQc.GENOME_DIR + plinkroot
                 + ".genome_keep.dat to " + RelationAncestryQc.ANCESTRY_DIR + UNRELATEDS_FILENAME);
      Files.copyFile(qcDir + RelationAncestryQc.GENOME_DIR + plinkroot + ".genome_keep.dat",
                     qcDir + RelationAncestryQc.ANCESTRY_DIR + UNRELATEDS_FILENAME);
    }
    PSF.checkInterrupted();
    if (!Files.exists(qcDir + RelationAncestryQc.ANCESTRY_DIR + plinkroot + ".bed")) {
      log.report(ext.getTime() + "]\tRunning --extract " + plinkroot
                 + ".prune.in (again, this time to " + RelationAncestryQc.ANCESTRY_DIR + ")");
      CmdLine.runDefaults("plink2 --bfile ../" + RelationAncestryQc.GENOME_DIR + plinkroot
                          + " --make-bed --noweb --out " + plinkroot,
                          qcDir + RelationAncestryQc.ANCESTRY_DIR, log);
    }
    PSF.checkInterrupted();

    System.out.println("Finished this round in " + ext.getTimeElapsed(time));
  }

  /**
   * @param keepGenomeInfoForRelatedsOnly true to save disk usage if unrelated genome info is not
   *          required
   * @return full path to plinkroot of QC'd plink dataset
   */
  public static void fullGamut(String dir, String plinkPrefix,
                               boolean keepGenomeInfoForRelatedsOnly, Logger log) {
    fullGamut(dir, plinkPrefix, keepGenomeInfoForRelatedsOnly, log, DEFAULT_QC_METRIC_THRESHOLDS);
  }

  /**
   * @param keepGenomeInfoForRelatedsOnly true to save disk usage if unrelated genome info is not
   *          required
   * @param markerQCThresholds TODO
   * @return full path to plinkroot of QC'd plink dataset
   */
  public static void fullGamut(String dir, String plinkPrefix,
                               boolean keepGenomeInfoForRelatedsOnly, Logger log,
                               Map<QcMetric, String> markerQCThresholds) {
    new RelationAncestryQc(dir, plinkPrefix, markerQCThresholds,
                           log).run(keepGenomeInfoForRelatedsOnly);
  }

  public static void fromParameters(String filename, Logger log) {
    List<String> params;

    params = Files.parseControlFile(filename, "gwas.Qc",
                                    new String[] {"dir=./",
                                                  "# Make keepGenomeInfoForRelatedsOnly=false if the sample size is small and you want to run MDS plot as well",
                                                  "keepGenomeInfoForRelatedsOnly=true"},
                                    log);

    if (params != null) {
      params.add("log=" + log.getFilename());
      main(ArrayUtils.toStringArray(params));
    }
  }

  public static void main(String[] args) {
    CLI c = new CLI(RelationAncestryQc.class);
    c.addArgWithDefault(CLI.ARG_INDIR, "directory with binary plink dataset", "./");
    c.addArgWithDefault(CLI.ARG_PLINKROOT, CLI.DESC_PLINKROOT, Qc.DEFAULT_PLINKROOT);
    c.addArgWithDefault(RelationAncestryQc.ARGS_KEEPGENOME, "if no MDS will be run, smaller file",
                        String.valueOf(true));
    Map<QcMetric, String> markerQCThresholds = Maps.newEnumMap(DEFAULT_QC_METRIC_THRESHOLDS);
    for (QcMetric metric : CUSTOMIZABLE_QC_METRICS) {
      String defaultThreshold = markerQCThresholds.get(metric);
      c.addArgWithDefault(metric.getKey(), metric.getCLIDescription(), defaultThreshold);
    }
    c.addArgWithDefault(CLI.ARG_LOG, CLI.DESC_LOG, "fullGamutOfMarkerAndSampleQC.log");

    c.parseWithExit(args);

    String dir = c.get(CLI.ARG_INDIR);
    String inputPlinkroot = c.get(CLI.ARG_PLINKROOT);
    boolean keepGenomeInfoForRelatedsOnly = Boolean.parseBoolean(c.get(RelationAncestryQc.ARGS_KEEPGENOME));
    for (QcMetric metric : CUSTOMIZABLE_QC_METRICS) {
      markerQCThresholds.put(metric, c.get(metric.getKey()));
    }
    Logger log = new Logger(dir + c.get(CLI.ARG_LOG));

    try {
      RelationAncestryQc.fullGamut(dir, inputPlinkroot, keepGenomeInfoForRelatedsOnly, log,
                                   markerQCThresholds);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
