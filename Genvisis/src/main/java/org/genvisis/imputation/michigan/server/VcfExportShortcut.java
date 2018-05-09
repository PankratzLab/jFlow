package org.genvisis.imputation.michigan.server;

import java.util.HashMap;
import java.util.Map;
import java.util.Optional;
import org.genvisis.CLI;
import org.genvisis.cnv.LaunchProperties;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.filesys.Project.SOURCE_FILE_DELIMITERS;
import org.genvisis.cnv.filesys.SourceFileHeaderData;
import org.genvisis.cnv.gui.ProjectCreationGUI;
import org.genvisis.cnv.manage.SourceFileParser;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.common.collect.MultisetUtils;
import org.genvisis.gwas.FurtherAnalysisQc;
import org.genvisis.gwas.MarkerQC;
import org.genvisis.gwas.MarkerQC.QC_METRIC;
import org.genvisis.gwas.Qc;
import org.genvisis.imputation.ImputationPipeline;
import org.genvisis.imputation.ImputationPipeline.IMPUTATION_PIPELINE_PATH;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;

public class VcfExportShortcut {

  private final Project proj;
  private final Logger log;

  private String manifest = null;

  private double callrateThresh = 0.98; // MarkerQC.DEFAULT_ILLUMINA_CALLRATE_THRESHOLD;
  private double hweThresh = 0.0000001; // 1E-7, FurtherAnalysisQc.BONFERRONI_CORRECTED_P_THRESHOLD
  private double generalQcThresh = 0.0000001; // 1E-7, FurtherAnalysisQc.BONFERRONI_CORRECTED_P_THRESHOLD
  private String sampleDropsFile = null;
  private String markerDropsFile = null;

  private VcfExportShortcut(String projName, String sourceDir, String projDir, String pedigreeFile,
                            Logger log) {
    proj = createOrGetProject(projName, sourceDir, projDir, pedigreeFile, log);
    this.log = log;
  }

  private void setIlluminaManifest(String string) {
    this.manifest = string;
  }

  private void setQCThreshold(double d) {
    this.generalQcThresh = d;
  }

  private void setHWEThreshold(double d) {
    this.hweThresh = d;
  }

  private void setCallrateThreshold(double d) {
    this.callrateThresh = d;
  }

  private static Project createOrGetProject(String projName, String sourceDir, String outDir,
                                            String pedigreeFile, Logger log) {
    final Project proj;
    if (LaunchProperties.projectExists(projName)) {
      String projFile = LaunchProperties.formProjectPropertiesFilename(projName);
      log.reportTime("Found existing project - " + projFile);
      proj = new Project(projFile);
    } else {
      proj = Project.initializeProject(projName);
      proj.SOURCE_DIRECTORY.setValue(sourceDir);
      MultisetUtils.maxCount(ProjectCreationGUI.getSourceExtensionCounts(sourceDir))
                   .map(Multiset.Entry::getElement)
                   .ifPresent(proj.SOURCE_FILENAME_EXTENSION::setValue);
      proj.PROJECT_DIRECTORY.setValue(outDir);
      proj.ARRAY_TYPE.setValue(ARRAY.ILLUMINA);
      Map<String, SourceFileHeaderData> headers = SourceFileHeaderData.validate(sourceDir,
                                                                                proj.SOURCE_FILENAME_EXTENSION.getValue(),
                                                                                false, log,
                                                                                Optional.empty());
      SourceFileHeaderData firstHeader = headers.values().stream().findFirst()
                                                .orElseThrow(() -> new IllegalStateException("No source files parsed"));
      if (firstHeader.getColSampleIdent() >= 0) {
        proj.ID_HEADER.setValue(firstHeader.getCols()[firstHeader.getColSnpIdent()]);
      } else {
        proj.ID_HEADER.setValue(SourceFileParser.FILENAME_AS_ID_OPTION);
      }

      proj.SOURCE_FILE_DELIMITER.setValue(SOURCE_FILE_DELIMITERS.getDelimiter(firstHeader.getSourceFileDelimiter()));
      proj.saveProperties();
      proj.setSourceFileHeaders(headers);
      proj.setLog(log);
      checkAndSetPed(proj, pedigreeFile);
      log.reportTime("Created Genvisis project at " + proj.getPropertyFilename());
    }
    return proj;
  }

  private static void checkAndSetPed(Project proj, String pedigreeFile) {
    proj.getLog()
        .reportTime("Validating Pedigree file, assuming standard pedigree.dat file format (FID, IID, FA, MO, SEX, PHENO, DNA)");
    String[] iids = HashVec.loadFileToStringArray(pedigreeFile, false, new int[] {1}, false);
    String[] uniqs = ArrayUtils.unique(iids);

    if (iids.length != uniqs.length) {
      Multiset<String> idSet = HashMultiset.create();
      for (String i : iids) {
        idSet.add(i);
      }
      idSet.removeIf((s) -> {
        return idSet.count(s) == 0;
      });
      int dupCount = idSet.size();
      StringBuilder dupSB = new StringBuilder(" duplicate IDs -- ");
      for (String s : idSet.elementSet()) {
        dupSB.append(s).append(": ").append(idSet.count(s)).append(", ");
      }
      dupSB.deleteCharAt(dupSB.length() - 1);
      dupSB.deleteCharAt(dupSB.length() - 1);
      throw new IllegalArgumentException("Error - IIDs in Pedigree file must be unique, found "
                                         + dupCount + dupSB.toString());
    }
    proj.PEDIGREE_FILENAME.setValue(pedigreeFile);
  }

  public void runPipeline(int numThreads, String putativeWhtFile, String outputDirAndRoot,
                          String imputationReferenceFile, boolean useGRC) {
    Map<QC_METRIC, String> qcMap = new HashMap<>();
    qcMap.put(QC_METRIC.CALLRATE, "<" + callrateThresh);
    qcMap.put(QC_METRIC.HWE, "<" + hweThresh);
    String generalQc = "<" + generalQcThresh;
    qcMap.put(QC_METRIC.MISHAP_HETERO, generalQc);
    qcMap.put(QC_METRIC.MISHAP_MIN, generalQc);
    qcMap.put(QC_METRIC.P_MISS, generalQc);
    qcMap.put(QC_METRIC.P_GENDER, generalQc);
    qcMap.put(QC_METRIC.P_GENDER_MISS, generalQc);
    StringBuilder outputScript = new StringBuilder(GenvisisWorkflow.setupImputation(proj,
                                                                                    numThreads,
                                                                                    putativeWhtFile,
                                                                                    qcMap, true,
                                                                                    manifest));

    String sampDrop = null;
    if (sampleDropsFile != null) {
      sampDrop = sampleDropsFile;
    } else {
      sampDrop = proj.PROJECT_DIRECTORY.getValue() + GenvisisWorkflow.PLINK_SUBDIR
                 + FurtherAnalysisQc.FURTHER_ANALYSIS_DIR + FurtherAnalysisQc.SAMPLE_QC_DROPS;
    }

    String markDrop = null;
    if (markerDropsFile != null) {
      markDrop = markerDropsFile;
    } else {
      markDrop = proj.PROJECT_DIRECTORY.getValue() + GenvisisWorkflow.PLINK_SUBDIR + Qc.QC_SUBDIR
                 + FurtherAnalysisQc.FURTHER_ANALYSIS_DIR + FurtherAnalysisQc.MARKER_QC_DROPS;
    }

    outputScript.append("## Imputation VCF Export\n");
    outputScript.append("echo \">>>> start Imputation VCF Export at: \" `date`").append("\n");

    outputScript.append(Files.getRunString());
    outputScript.append(" ").append(ImputationPipeline.class.getName());
    outputScript.append(" ").append(ImputationPipeline.PROJ_ARG).append(proj.getPropertyFilename());
    outputScript.append(" ").append(ImputationPipeline.REF_ARG).append(imputationReferenceFile);
    outputScript.append(" ").append(ImputationPipeline.DROP_SAMPLES_ARG).append(sampDrop);
    outputScript.append(" ").append(ImputationPipeline.DROP_MARKERS_ARG).append(markDrop);
    outputScript.append(" ").append(ImputationPipeline.EXPORT_IIDS).append("TRUE");
    outputScript.append(" ").append(ImputationPipeline.USE_GRC_ARG)
                .append(Boolean.toString(useGRC));
    outputScript.append(" ").append(ImputationPipeline.OUT_DIR_AND_ROOT_ARG)
                .append(outputDirAndRoot);
    outputScript.append(" ").append(ImputationPipeline.RUN_TYPE_ARG)
                .append(IMPUTATION_PIPELINE_PATH.VCF_ONLY);
    outputScript.append("\n").append("echo \"<<<< end Imputation VCF Export at: \" `date`")
                .append("\n");
    outputScript.append("\n\n");

    String file = proj.PROJECT_DIRECTORY.getValue() + "ImputationExport.";
    String suggFile = file + ext.getTimestampForFilename() + ".run";

    Files.write(outputScript.toString(), suggFile);

    log.reportTime("Created Imputation Export script at " + suggFile);
  }

  private static final String ARG_PROJ_NAME = "projName";
  private static final String ARG_SRC_DIR = "sourceDir";
  private static final String ARG_PROJ_DIR = "projDir";

  private static final String ARG_ILL_MAN = "manifest";

  private static final String ARG_PED = "pedFile";
  private static final String ARG_PUT_WHT = "putativeWhiteFile";
  private static final String ARG_VCF_OUT = "vcfDirRoot";
  private static final String ARG_IMP_REF = "imputationReference";
  private static final String ARG_GRC_OUT = "useGRC";
  private static final String ARG_CALLRATE = "callrate";
  private static final String ARG_HWE = "hwe";
  private static final String ARG_QC = "qc";

  private static final String DESC_PROJ_NAME = "Project Name";
  private static final String DESC_SRC_DIR = "Source File Directory";
  private static final String DESC_PROJ_DIR = "Directory for project files";

  private static final String DESC_ILL_MAN = "Illumina manifest file (e.g. HumanOmni2.5-4v1_H.csv)";

  private static final String DESC_PED = "Pedigree File";
  private static final String DESC_PUT_WHT = "Putative whites file (file with two columns of FID/IID and no header)";
  private static final String DESC_VCF_OUT = "VCF output directory and file root";
  private static final String DESC_IMP_REF = "Imputation marker reference file";
  private static final String DESC_GRC_OUT = "Export contigs with \"chr\" prepend (defaults to true)";
  private static final String DESC_CALLRATE = "Callrate Threshold";
  private static final String DESC_HWE = "HWE Threshold";
  private static final String DESC_QC = "General QC Threshold";

  public static void main(String[] args) {
    CLI cli = new CLI(VcfExportShortcut.class);

    cli.addArg(ARG_PROJ_NAME, DESC_PROJ_NAME, true);
    cli.addArg(ARG_SRC_DIR, DESC_SRC_DIR, true);
    cli.addArg(ARG_PROJ_DIR, DESC_PROJ_DIR, true);
    cli.addArg(ARG_VCF_OUT, DESC_VCF_OUT, true);
    cli.addArg(ARG_IMP_REF, DESC_IMP_REF, true);
    cli.addArg(ARG_PUT_WHT, DESC_PUT_WHT, true);
    cli.addArg(ARG_PED, DESC_PED, false);

    cli.addArg(ARG_ILL_MAN, DESC_ILL_MAN, true);

    cli.addArg(CLI.ARG_THREADS, CLI.DESC_THREADS,
               Integer.toString(Runtime.getRuntime().availableProcessors()), false);
    cli.addArg(ARG_GRC_OUT, DESC_GRC_OUT, "true", false);
    cli.addArg(ARG_CALLRATE, DESC_CALLRATE,
               MarkerQC.DEFAULT_ILLUMINA_CALLRATE_THRESHOLD.substring(1), false);
    cli.addArg(ARG_HWE, DESC_HWE, FurtherAnalysisQc.BONFERRONI_CORRECTED_P_THRESHOLD.substring(1),
               false);
    cli.addArg(ARG_QC, DESC_QC, FurtherAnalysisQc.BONFERRONI_CORRECTED_P_THRESHOLD.substring(1),
               false);

    cli.parseWithExit(args);

    String projName = cli.get(ARG_PROJ_NAME);
    String srcDir = cli.get(ARG_SRC_DIR);
    String projDir = cli.get(ARG_PROJ_DIR);
    String pedFile = cli.get(ARG_PED);
    String putWht = cli.get(ARG_PUT_WHT);
    String vcfOut = cli.get(ARG_VCF_OUT);
    String impRef = cli.get(ARG_IMP_REF);

    int numThreads = cli.getI(CLI.ARG_THREADS);
    boolean useGRC = Boolean.parseBoolean(cli.get(ARG_GRC_OUT));

    Logger log = new Logger();
    VcfExportShortcut export = new VcfExportShortcut(projName, srcDir, projDir, pedFile, log);

    export.setIlluminaManifest(cli.get(ARG_ILL_MAN));
    export.setCallrateThreshold(cli.getD(ARG_CALLRATE));
    export.setHWEThreshold(cli.getD(ARG_HWE));
    export.setQCThreshold(cli.getD(ARG_QC));

    export.runPipeline(numThreads, putWht, vcfOut, impRef, useGRC);

    log.reportTime("Please check (and, if necessary, adjust) parameters (memory, threads, arguemnts) prior to submitting to queue.");
    log.reportTime("Sample and Marker \"keep\" files can be set in the ImputationExport script - the default files are created automatically and can be edited or used as-is.");

  }

}