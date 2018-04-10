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
import com.google.common.collect.Multiset;

public class VcfExportShortcut {

  private final Project proj;
  private final Logger log;

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
    GenvisisWorkflow.setupImputation(proj);
  }

}
