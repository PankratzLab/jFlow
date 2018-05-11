package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.MitoPipeline;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.Requirement.FileRequirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.Variables;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Elision;
import org.genvisis.common.Files;

public class SampleDataStep extends Step {

  public static final String REQ_CREATE_MINIMAL = "Create a minimal SampleData.txt file from sample files";
  public static final String REQ_PEDIGREE = "Either a Pedigree.dat file, or any file with a header containing all of the following elements (in any order):  \""
                                            + ArrayUtils.toStr(MitoPipeline.PED_INPUT, ", ") + "\"";

  public static SampleDataStep create(Project proj, Step parseSamplesStep) {
    Requirement<Step> parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    Requirement<Boolean> createMinimalSampleDataReq = new Requirement.BoolRequirement(REQ_CREATE_MINIMAL,
                                                                                      true);
    String pedPreset = proj.PEDIGREE_FILENAME.getValue();

    Requirement<File> pedigreeReq = new FileRequirement(REQ_PEDIGREE, new File(pedPreset));

    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(createMinimalSampleDataReq)
                                                                                 .add(pedigreeReq));
    return new SampleDataStep(proj, createMinimalSampleDataReq, pedigreeReq, reqSet);
  }

  public static final String NAME = "Create SampleData.txt File";
  public static final String DESC = "";

  final Project proj;
  final Requirement<Boolean> createMinimalSampleDataReq;
  final Requirement<File> pedigreeReq;

  public SampleDataStep(Project proj, Requirement<Boolean> createMin, Requirement<File> pedReq,
                        RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.noneOf(Requirement.Flag.class));
    this.proj = proj;
    this.createMinimalSampleDataReq = createMin;
    this.pedigreeReq = pedReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    // Nothing to do
  }

  @Override
  public void run(Variables variables) {
    Boolean minimal = variables.get(createMinimalSampleDataReq);
    File pedFile = variables.get(pedigreeReq);
    String pedFilePath = null;
    if (!minimal && pedFile != null) {
      pedFilePath = pedFile.getAbsolutePath();
    }

    proj.getLog().report("Creating SampleData.txt");
    try {
      int retStat = SampleData.createSampleData(pedFilePath, null, proj);
      if (retStat == -1) {
        throw new RuntimeException("Error during SampleData creation - please check log and try again.");
      }
    } catch (Elision e) {
      throw new RuntimeException(e.getMessage());
    }
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    return Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false));
  }

  @Override
  public String getCommandLine(Variables variables) {
    String projPropFile = proj.getPropertyFilename();
    boolean minimal = variables.get(createMinimalSampleDataReq);
    File pedFile = minimal ? null : variables.get(pedigreeReq);
    StringBuilder cmd = new StringBuilder();
    cmd.append(Files.getRunString()).append(" cnv.var.SampleData proj=").append(projPropFile);
    if (pedFile != null && Files.exists(pedFile)) {
      cmd.append(" ped=").append(pedFile.getAbsolutePath());
    }
    return cmd.toString();
  }

}
