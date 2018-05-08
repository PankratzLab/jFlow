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
  public static final String REQ_SAMPLE_MAP = "A Sample_Map.csv file, with at least two columns having headers \""
                                              + MitoPipeline.SAMPLEMAP_INPUT[1] + "\" and \""
                                              + MitoPipeline.SAMPLEMAP_INPUT[2] + "\"";

  public static SampleDataStep create(Step parseSamplesStep, Project proj) {
    Requirement<Step> parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    Requirement<Boolean> createMinimalSampleDataReq = new Requirement.BoolRequirement(REQ_CREATE_MINIMAL,
                                                                                      true);
    String pedPreset = proj.PEDIGREE_FILENAME.getValue();

    Requirement<File> pedigreeReq = new FileRequirement(REQ_PEDIGREE, new File(pedPreset));

    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(createMinimalSampleDataReq)
                                                                                 .add(pedigreeReq));
    return new SampleDataStep(createMinimalSampleDataReq, pedigreeReq, reqSet);
  }

  public static final String NAME = "Create SampleData.txt File";
  public static final String DESC = "";

  final Requirement<Boolean> createMinimalSampleDataReq;
  final Requirement<File> pedigreeReq;

  public SampleDataStep(Requirement<Boolean> createMin, Requirement<File> pedReq,
                        RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.noneOf(Requirement.Flag.class));
    this.createMinimalSampleDataReq = createMin;
    this.pedigreeReq = pedReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Project proj, Variables variables) {
    // Nothing to do
  }

  @Override
  public void run(Project proj, Variables variables) {
    Boolean minimal = variables.get(createMinimalSampleDataReq);
    File pedFile = minimal ? null : variables.get(pedigreeReq);

    proj.getLog().report("Creating SampleData.txt");
    try {
      int retStat = SampleData.createSampleData(pedFile.getAbsolutePath(), null, proj);
      if (retStat == -1) {
        throw new RuntimeException("SampleData already exists - please delete and try again.");
      }
    } catch (Elision e) {
      throw new RuntimeException(e.getMessage());
    }
  }

  @Override
  public boolean checkIfOutputExists(Project proj, Variables variables) {
    return Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false));
  }

  @Override
  public String getCommandLine(Project proj, Variables variables) {
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
