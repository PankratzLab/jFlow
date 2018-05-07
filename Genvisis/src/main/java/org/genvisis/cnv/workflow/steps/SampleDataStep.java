package org.genvisis.cnv.workflow.steps;

import java.util.EnumSet;
import java.util.Map;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.MitoPipeline;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.Requirement.FileRequirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
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
    Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    Requirement createMinimalSampleDataReq = new Requirement.BoolRequirement(REQ_CREATE_MINIMAL,
                                                                             true);
    String pedPreset = proj.PEDIGREE_FILENAME.getValue();

    Requirement pedigreeReq = new FileRequirement(REQ_PEDIGREE, pedPreset);

    // check for SampleMap only if we haven't found a pedigree
    final String sampMapPreset = Files.exists(pedPreset) ? null
                                                         : GenvisisWorkflow.getLocationOfSampleMap(proj);

    final Requirement sampMapReq = new FileRequirement(REQ_SAMPLE_MAP, sampMapPreset);
    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(createMinimalSampleDataReq)
                                                                                 .add(pedigreeReq)
                                                                                 .add(sampMapReq));
    return new SampleDataStep(createMinimalSampleDataReq, pedigreeReq, sampMapReq, reqSet);
  }

  public static final String NAME = "Create SampleData.txt File";
  public static final String DESC = "";

  final Requirement createMinimalSampleDataReq;
  final Requirement pedigreeReq;
  final Requirement sampMapReq;

  public SampleDataStep(Requirement createMin, Requirement pedReq, Requirement sampMap,
                        RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.noneOf(Requirement.Flag.class));
    this.createMinimalSampleDataReq = createMin;
    this.pedigreeReq = pedReq;
    this.sampMapReq = sampMap;
  }

  @Override
  public void setNecessaryPreRunProperties(Project proj,
                                           Map<Step, Map<Requirement, String>> variables) {
    // Nothing to do
  }

  @Override
  public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
    Boolean minimal = Boolean.parseBoolean(variables.get(this).get(createMinimalSampleDataReq));
    String pedFile = minimal ? null : variables.get(this).get(pedigreeReq);
    String sampleMapCsv = minimal ? null : variables.get(this).get(sampMapReq);

    proj.getLog().report("Creating SampleData.txt");
    try {
      int retStat = SampleData.createSampleData(pedFile, sampleMapCsv, proj);
      if (retStat == -1) {
        throw new RuntimeException("SampleData already exists - please delete and try again.");
      }
    } catch (Elision e) {
      throw new RuntimeException(e.getMessage());
    }
  }

  @Override
  public boolean checkIfOutputExists(Project proj, Map<Step, Map<Requirement, String>> variables) {
    return Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false));
  }

  @Override
  public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String projPropFile = proj.getPropertyFilename();
    Boolean minimal = Boolean.parseBoolean(variables.get(this).get(createMinimalSampleDataReq));
    String pedFile = minimal ? "" : variables.get(this).get(pedigreeReq);
    String sampleMapCsv = minimal ? "" : variables.get(this).get(sampMapReq);
    StringBuilder cmd = new StringBuilder();
    cmd.append(Files.getRunString()).append(" cnv.var.SampleData proj=").append(projPropFile);
    if (!"".equals(pedFile)) {
      cmd.append(" ped=").append(pedFile);
    }
    if (!"".equals(sampleMapCsv)) {
      cmd.append(" sampleMap=").append(sampleMapCsv);
    }
    return cmd.toString();
  }

}
