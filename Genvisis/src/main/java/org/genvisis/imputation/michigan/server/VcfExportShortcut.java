package org.genvisis.imputation.michigan.server;

import java.util.HashMap;
import java.util.Map;
import java.util.Optional;
import org.genvisis.cnv.LaunchProperties;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.filesys.Project.SOURCE_FILE_DELIMITERS;
import org.genvisis.cnv.filesys.SourceFileHeaderData;
import org.genvisis.cnv.gui.ProjectCreationGUI;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.SourceFileParser;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.common.collect.MultisetUtils;
import org.genvisis.gwas.FurtherAnalysisQc;
import org.genvisis.gwas.MarkerQC.QC_METRIC;
import org.genvisis.gwas.Qc;
import org.genvisis.imputation.ImputationPipeline;
import org.genvisis.imputation.ImputationPipeline.IMPUTATION_PIPELINE_PATH;
import org.genvisis.qsub.Qsub;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;

public class VcfExportShortcut {

  private final Project proj;
  private final Logger log;

  private int numThreads;
  private String putativeWhtFile;
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

  public void runPipeline() {
    // TODO set options / variables
    Map<QC_METRIC, String> qcMap = new HashMap<>();
    qcMap.put(QC_METRIC.CALLRATE, "<" + callrateThresh);
    qcMap.put(QC_METRIC.HWE, "<" + hweThresh);
    String generalQc = "<" + generalQcThresh;
    qcMap.put(QC_METRIC.MISHAP_HETERO, generalQc);
    qcMap.put(QC_METRIC.MISHAP_MIN, generalQc);
    qcMap.put(QC_METRIC.P_MISS, generalQc);
    qcMap.put(QC_METRIC.P_GENDER, generalQc);
    qcMap.put(QC_METRIC.P_GENDER_MISS, generalQc);
    String outputFile = GenvisisWorkflow.setupImputation(proj, numThreads, putativeWhtFile, qcMap);

    log.reportTime("Created Genvisis Workflow script at " + outputFile);
  }

  public void setupVCFExport(String outputDirAndRoot, boolean useGRC) {
    String refFile = useGRC ? Resources.genome(proj.GENOME_BUILD_VERSION.getValue(), proj.getLog())
                                       .getGRCFASTA().getAbsolute()
                            : Resources.genome(proj.GENOME_BUILD_VERSION.getValue(), proj.getLog())
                                       .getFASTA().getAbsolute();

    String sampDrop = null;
    if (sampleDropsFile != null) {
      sampDrop = sampleDropsFile;
    } else {
      sampDrop = proj.PROJECT_DIRECTORY.getValue() + GenvisisWorkflow.PLINK_SUBDIR
                 + FurtherAnalysisQc.FURTHER_ANALYSIS_DIR + FurtherAnalysisQc.SAMPLE_QC_DROPS;
    }

    // TODO WHERE ARE THESE FILES ANYWAY!?!?!  MUST TEST!
    String markDrop = null;
    if (markerDropsFile != null) {
      markDrop = markerDropsFile;
    } else {
      markDrop = proj.PROJECT_DIRECTORY.getValue() + GenvisisWorkflow.PLINK_SUBDIR + Qc.QC_SUBDIR
                 + FurtherAnalysisQc.FURTHER_ANALYSIS_DIR + FurtherAnalysisQc.MARKER_QC_DROPS;
    }

    StringBuilder exportSB = new StringBuilder();
    exportSB.append(Files.getRunString());
    exportSB.append(ImputationPipeline.class.getName());
    exportSB.append(" ").append(ImputationPipeline.PROJ_ARG).append(proj.getPropertyFilename());
    exportSB.append(" ").append(ImputationPipeline.REF_ARG).append(refFile);
    exportSB.append(" ").append(ImputationPipeline.DROP_SAMPLES_ARG).append(sampDrop);
    exportSB.append(" ").append(ImputationPipeline.DROP_MARKERS_ARG).append(markDrop);
    exportSB.append(" ").append(ImputationPipeline.EXPORT_IIDS).append("TRUE");
    exportSB.append(" ").append(ImputationPipeline.USE_GRC_ARG).append(Boolean.toString(useGRC));
    exportSB.append(" ").append(ImputationPipeline.OUT_DIR_AND_ROOT_ARG).append(outputDirAndRoot);
    exportSB.append(" ").append(ImputationPipeline.RUN_TYPE_ARG)
            .append(IMPUTATION_PIPELINE_PATH.VCF_ONLY);

    String file = proj.PROJECT_DIRECTORY.getValue() + "ImputationExport.";
    String suggFile = file + ext.getTimestampForFilename() + ".qsub";

    Qsub.qsubDefaults(suggFile, exportSB.toString());

    log.reportTime("Created Imputation Export script at " + suggFile);
  }

  /*-
   * TODO:
   *  VcFExportShortcut export = new VcfExportShortcut(...);
   *  export.runPipeline();
   *  "Run/Submit Genvisis Workflow script first"
   *  export.setupVCFExport();
   *  "Run/Submit Imputation Export script second"
   *  "Can set or edit sample/marker drops files"
   */

}
