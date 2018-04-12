package org.genvisis.imputation.michigan.server;

import java.util.Map;
import java.util.Optional;
import org.genvisis.cnv.LaunchProperties;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.filesys.Project.SOURCE_FILE_DELIMITERS;
import org.genvisis.cnv.filesys.SourceFileHeaderData;
import org.genvisis.cnv.gui.ProjectCreationGUI;
import org.genvisis.cnv.manage.SourceFileParser;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.common.Logger;
import org.genvisis.common.collect.MultisetUtils;
import org.genvisis.gwas.FurtherAnalysisQc;
import com.google.common.collect.Multiset;

public class VcfExportShortcut {

  private final Project proj;
  private final Logger log;

  private int numThreads;
  private String pedigreeFile;
  private String putativeWhtFile;
  private double callrateThresh = 0.98; // MarkerQC.DEFAULT_ILLUMINA_CALLRATE_THRESHOLD;
  private double hweThresh = 0.0000001; // 1E-7, FurtherAnalysisQc.BONFERRONI_CORRECTED_P_THRESHOLD
  private double generalQcThresh = 0.0000001; // 1E-7, FurtherAnalysisQc.BONFERRONI_CORRECTED_P_THRESHOLD
  private String sampleDropsFile = FurtherAnalysisQc.SAMPLE_QC_DROPS;
  private String markerDropsFile = FurtherAnalysisQc.MARKER_QC_DROPS;

  public VcfExportShortcut(String projName, String sourceDir, String outDir, Logger log) {
    proj = createOrGetProject(projName, sourceDir, outDir, log);
    this.log = log;
  }

  private static Project createOrGetProject(String projName, String sourceDir, String outDir,
                                            Logger log) {
    final Project proj;
    if (LaunchProperties.projectExists(projName)) {
      proj = new Project(LaunchProperties.formProjectPropertiesFilename(projName));
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
    }
    return proj;
  }

  public void runPipeline() {
    // TODO set options / variables
    GenvisisWorkflow.setupImputation(proj);
  }

  //  public void runVCFExport() {
  //    String refFile;
  //    String dropSamp, keepSamp, dropMark, keepMark;
  //    String outputDirAndRoot;
  //    int[] chrs;
  //    boolean useGRC;
  //    KeepDrops kd = new KeepDrops(dropSamp, keepSamp, dropMark, keepMark);
  //    ImputationPipeline exportPipe = new ImputationPipeline(proj, refFile, kd);
  //    exportPipe.exportToVCF(outputDirAndRoot, chrs, useGRC);
  //  }

}
