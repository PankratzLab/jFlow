package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;
import org.genvisis.cnv.Resources;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.GcAdjustor;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.Requirement.ResourceRequirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.Variables;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.CLI;

public class GCModelStep extends Step {

  public static final String NAME = "Compute GCMODEL File";
  public static final String DESC = "";

  public static GCModelStep create(Project proj) {
    final ResourceRequirement gcBaseResourceReq = new ResourceRequirement("GC Base file",
                                                                          Resources.genome(proj.GENOME_BUILD_VERSION.getValue(),
                                                                                           proj.getLog())
                                                                                   .getModelBase());
    final Requirement<File> gcModelOutputReq = new Requirement.OutputFileRequirement("GCModel output file must be specified.",
                                                                                     new File(proj.GC_MODEL_FILENAME.getValue()));
    final RequirementSet reqSet = RequirementSetBuilder.and().add(gcBaseResourceReq)
                                                       .add(gcModelOutputReq);
    return new GCModelStep(proj, gcBaseResourceReq, gcModelOutputReq, reqSet);
  }

  final Project proj;
  final Requirement.ResourceRequirement gcBaseResourceReq;
  final Requirement<File> gcModelOutputReq;

  private GCModelStep(Project proj, ResourceRequirement gcBase, Requirement<File> output,
                      RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.noneOf(Requirement.Flag.class));
    this.gcBaseResourceReq = gcBase;
    this.gcModelOutputReq = output;
    this.proj = proj;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    String setGCOutputFile = proj.GC_MODEL_FILENAME.getValue();
    String gcOutputFile = variables.get(gcModelOutputReq).getAbsolutePath();
    if (!ext.verifyDirFormat(setGCOutputFile).equals(gcOutputFile)) {
      proj.GC_MODEL_FILENAME.setValue(gcOutputFile);
    }
  }

  @Override
  public void run(Variables variables) {
    String gcBaseFile = gcBaseResourceReq.getResource().getAbsolute();
    String gcOutputFile = variables.get(gcModelOutputReq).getAbsolutePath();
    org.genvisis.cnv.qc.GcAdjustor.GcModel.gcModel(proj, gcBaseFile, gcOutputFile, 100);
  }

  @Override
  public String getCommandLine(Variables variables) {
    String kvCmd = "";

    String setGCOutputFile = proj.GC_MODEL_FILENAME.getValue();
    String gcOutputFile = variables == null ? null
                                            : variables.get(gcModelOutputReq).getAbsolutePath();
    if (gcOutputFile != null && !ext.verifyDirFormat(setGCOutputFile).equals(gcOutputFile)) {
      kvCmd += " GC_MODEL_FILENAME=" + gcOutputFile;
    }

    String projPropFile = proj.getPropertyFilename();
    StringBuilder cmd = new StringBuilder();
    if (kvCmd.length() > 0) {
      cmd.append(Files.getRunString()).append(GenvisisWorkflow.PROJ_PROP_UPDATE_STR + projPropFile)
         .append(kvCmd).append("\n");
    }
    String gcBaseFile = gcBaseResourceReq.getResource().getAbsolute();
    return cmd.append(Files.getRunString())
              .append(" " + GcAdjustor.class.getName() + " " + CLI.ARG_PROJ + "="
                      + proj.getPropertyFilename() + " " + CLI.ARG_LOG + "="
                      + proj.getLog().getFilename() + " " + GcAdjustor.GC_BASE_FILE + "="
                      + gcBaseFile)
              .toString();
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    return Files.exists(variables.get(gcModelOutputReq));
  }

}
