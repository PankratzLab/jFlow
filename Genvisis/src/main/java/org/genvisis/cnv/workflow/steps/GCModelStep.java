package org.genvisis.cnv.workflow.steps;

import java.util.EnumSet;
import java.util.Map;
import org.genvisis.CLI;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.qc.GcAdjustor;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.Requirement.ResourceRequirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.common.Files;
import org.genvisis.common.ext;

public class GCModelStep extends Step {

  public static final String NAME = "Compute GCMODEL File";
  public static final String DESC = "";

  public static GCModelStep create(Project proj, double priority) {
    final ResourceRequirement gcBaseResourceReq = new ResourceRequirement("GC Base file",
                                                                          Resources.genome(proj.GENOME_BUILD_VERSION.getValue(),
                                                                                           proj.getLog())
                                                                                   .getModelBase());
    final Requirement gcModelOutputReq = new Requirement.OutputFileRequirement("GCModel output file must be specified.",
                                                                               proj.GC_MODEL_FILENAME.getValue());
    final RequirementSet reqSet = RequirementSetBuilder.and().add(gcBaseResourceReq)
                                                       .add(gcModelOutputReq);
    return new GCModelStep(gcBaseResourceReq, gcModelOutputReq, reqSet, priority);
  }

  final Requirement.ResourceRequirement gcBaseResourceReq;
  final Requirement gcModelOutputReq;

  private GCModelStep(ResourceRequirement gcBase, Requirement output, RequirementSet reqSet,
                      double priority) {
    super(NAME, DESC, reqSet, EnumSet.noneOf(Requirement.Flag.class), priority);
    this.gcBaseResourceReq = gcBase;
    this.gcModelOutputReq = output;
  }

  @Override
  public void setNecessaryPreRunProperties(Project proj,
                                           Map<Step, Map<Requirement, String>> variables) {
    String setGCOutputFile = proj.GC_MODEL_FILENAME.getValue();
    String gcOutputFile = variables.get(this).get(gcModelOutputReq);
    if (!ext.verifyDirFormat(setGCOutputFile).equals(gcOutputFile)) {
      proj.GC_MODEL_FILENAME.setValue(gcOutputFile);
    }
  }

  @Override
  public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String gcBaseFile = gcBaseResourceReq.getResource().getAbsolute();
    String gcOutputFile = variables.get(this).get(gcModelOutputReq);
    org.genvisis.cnv.qc.GcAdjustor.GcModel.gcModel(proj, gcBaseFile, gcOutputFile, 100);
  }

  @Override
  public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String kvCmd = "";

    String setGCOutputFile = proj.GC_MODEL_FILENAME.getValue();
    String gcOutputFile = variables == null ? null : variables.get(this).get(gcModelOutputReq);
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
  public boolean checkIfOutputExists(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String gcOutputFile = variables.get(this).get(gcModelOutputReq);
    return Files.exists(gcOutputFile);
  }

}
