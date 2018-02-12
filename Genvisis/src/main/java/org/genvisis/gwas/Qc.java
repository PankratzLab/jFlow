package org.genvisis.gwas;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.gwas.MarkerQC.QC_METRIC;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;

public abstract class Qc {

  protected static final String DEFAULT_PLINKROOT = "plink";
  public static final String QC_SUBDIR = "quality_control/";
  public static final String MARKER_QC_DROPS = "miss_drops.dat";

  protected final String sourceDir;
  protected final String qcDir;
  protected final String plinkroot;
  protected final Logger log;
  protected final Map<QC_METRIC, String> markerQCThresholds;

  /**
   * @param sourceDir Directory with plink files to run from
   * @param plinkPrefix prefix of plink binaries
   * @param markerQCThresholds thresholds to apply for each desired marker QC metric, {@code null}
   *          for defaults
   * @param log
   */
  protected Qc(String sourceDir, String plinkPrefix, Map<QC_METRIC, String> markerQCThresholds,
               Logger log) {
    super();
    sourceDir = ext.verifyDirFormat(sourceDir);
    if (!sourceDir.startsWith("/") && !sourceDir.contains(":")) {
      sourceDir = ext.verifyDirFormat((new File("./" + sourceDir)).getAbsolutePath());
    }
    this.sourceDir = sourceDir;
    this.qcDir = sourceDir + QC_SUBDIR;
    this.plinkroot = plinkPrefix == null ? DEFAULT_PLINKROOT : plinkPrefix;
    this.markerQCThresholds = markerQCThresholds == null ? MarkerQC.DEFAULT_METRIC_THRESHOLDS
                                                         : markerQCThresholds;
    this.log = log;
  }

  private boolean makeMarkerQCBED(String subDir, String sampleSubsetFile) {
    new File(qcDir + subDir).mkdirs();
    String inRoot = sourceDir + plinkroot;
    String outRoot = qcDir + subDir + plinkroot;
    List<String> commands = Lists.newArrayList("plink2", "--noweb", "--bfile", inRoot, "--make-bed",
                                               "--out", outRoot);
    if (sampleSubsetFile != null) {
      commands.add("--keep");
      commands.add(sampleSubsetFile);
    }
    Set<String> inputs = PSF.Plink.getPlinkBedBimFamSet(inRoot);
    Set<String> outputs = PSF.Plink.getPlinkBedBimFamSet(outRoot);
    return CmdLine.runCommandWithFileChecks(commands, "", inputs, outputs, true, false, true, log);
  }

  private boolean runHWE(String subDir, String mind10, String hweSampleSubsetFile) {
    List<String> commands = Lists.newArrayList("plink2", "--bfile", mind10, "--geno", "1", "--mind",
                                               "1", "--hardy", "--out", "hardy", "--noweb");
    if (hweSampleSubsetFile != null) {
      commands.add("--keep");
      commands.add(hweSampleSubsetFile);
    }
    Set<String> inputs = PSF.Plink.getPlinkBedBimFamSet(mind10);
    Set<String> outputs = ImmutableSet.of("hardy.hwe");
    return CmdLine.runCommandWithFileChecks(commands, qcDir + subDir, inputs, outputs, true, false,
                                            true, log);
  }

  protected boolean markerQc(String subDir) {
    return markerQc(subDir, null, null);
  }

  protected boolean markerQc(String subDir, String sampleSubsetFile, String hweSampleSubsetFile) {
    subDir = ext.verifyDirFormat(subDir);
    if (!makeMarkerQCBED(subDir, sampleSubsetFile)) {
      log.reportError("CRITICAL ERROR, creating initial PLINK files failed.");
      return false;
    }

    String geno20 = plinkroot + "_geno20";
    if (!Files.exists(qcDir + subDir + geno20 + ".bed")) {
      log.report(ext.getTime() + "]\tRunning --geno 0.2");
      CmdLine.runDefaults("plink2 --bfile " + plinkroot + " --geno 0.2 --make-bed --noweb --out ./"
                          + geno20, qcDir + subDir, log);
    }
    PSF.checkInterrupted();
    String geno20mind10 = geno20 + "_mind10";
    if (!Files.exists(qcDir + subDir + geno20mind10 + ".bed")) {
      log.report(ext.getTime() + "]\tRunning --mind 0.1");
      CmdLine.runDefaults("plink2 --bfile " + geno20 + " --mind 0.1 --make-bed --noweb --out ./"
                          + geno20mind10, qcDir + subDir, log);
    }
    PSF.checkInterrupted();
    String mind10 = plinkroot + "_mind10";
    if (!Files.exists(qcDir + subDir + mind10 + ".bed")) {
      log.report(ext.getTime() + "]\tRemoving trimmed samples from --mind 0.1");
      CmdLine.runDefaults("plink2 --bfile " + plinkroot + " --keep " + geno20mind10
                          + ".fam --make-bed --noweb --out ./" + mind10, qcDir + subDir, log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(qcDir + subDir + "freq.frq")) {
      log.report(ext.getTime() + "]\tRunning --freq");
      CmdLine.runDefaults("plink2 --bfile " + mind10
                          + " --geno 1 --mind 1 --freq --out freq --noweb", qcDir + subDir, log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(qcDir + subDir + "missing.imiss")) {
      log.report(ext.getTime() + "]\tRunning --missing");
      CmdLine.runDefaults("plink2 --bfile " + mind10
                          + " --geno 1 --mind 1 --missing --out missing --noweb", qcDir + subDir,
                          log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(qcDir + subDir + "test.missing.missing")) {
      log.report(ext.getTime() + "]\tRunning --test-missing");
      CmdLine.runDefaults("plink2 --bfile " + mind10
                          + " --geno 1 --mind 1 --test-missing --out test.missing --noweb",
                          qcDir + subDir, log);
    }
    PSF.checkInterrupted();
    runHWE(subDir, mind10, hweSampleSubsetFile);
    PSF.checkInterrupted();
    if (!Files.exists(qcDir + subDir + "mishap.missing.hap")) {
      log.report(ext.getTime() + "]\tRunning --test-mishap");
      CmdLine.runDefaults("plink2 --bfile " + geno20mind10
                          + " --geno 1 --mind 1 --test-mishap --out mishap --noweb", qcDir + subDir,
                          log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(qcDir + subDir + "gender.assoc")) {
      log.report(ext.getTime() + "]\tRunning --assoc gender");
      CmdLine.runDefaults("plink2 --bfile " + mind10 + " --geno 1 --mind 1 --pheno " + mind10
                          + ".fam --mpheno 3 --assoc --out gender --noweb", qcDir + subDir, log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(qcDir + subDir + "gender.missing")) {
      log.report(ext.getTime() + "]\tRunning --test-missing gender");
      CmdLine.runDefaults("plink2 --bfile " + mind10 + " --geno 1 --mind 1 --pheno " + mind10
                          + ".fam --mpheno 3 --test-missing --out gender --noweb", qcDir + subDir,
                          log);
    }
    PSF.checkInterrupted();
    if (!Files.exists(qcDir + subDir + MARKER_QC_DROPS)) {
      MarkerQC.generateCRF(qcDir + subDir, qcDir + subDir + "miss.crf", markerQCThresholds);
      int runCode = MarkerQC.parseParameters(qcDir + subDir + "miss.crf", log, false);
      if (runCode != 0) {
        log.reportError("Failed to perform marker QC with " + qcDir + subDir + "miss.crf");
        return false;
      }
    }
    PSF.checkInterrupted();
    return true;
  }

  /**
   * @deprecated Use {@link RelationAncestryQc#fromParameters(String,Logger)} instead
   */
  @Deprecated
  public static void fromParameters(String filename, Logger log) {
    RelationAncestryQc.fromParameters(filename, log);
  }

  /**
   * @deprecated Use {@link RelationAncestryQc#main(String[])} instead
   */
  @Deprecated
  public static void main(String[] args) {
    RelationAncestryQc.main(args);
  }
}
