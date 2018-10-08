package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.hmm.PFB;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.Variables;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.CLI;

public class ComputePFBStep extends Step {

  public static final String NAME = "Compute Population BAF files";
  public static final String DESC = "";

  public static ComputePFBStep create(Project proj, final Step parseSamplesStep) {
    final Requirement<Step> parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    final Requirement<File> sampleSubsetReq = new Requirement.FileRequirement("A Sample subset file must exist.",
                                                                              new File(proj.SAMPLE_SUBSET_FILENAME.getValue()));
    String defaultOutputFile;
    if (Files.exists(proj.SAMPLE_SUBSET_FILENAME.getValue())) {
      defaultOutputFile = ext.rootOf(proj.SAMPLE_SUBSET_FILENAME.getValue()) + ".pfb";
    } else {
      defaultOutputFile = proj.CUSTOM_PFB_FILENAME.getValue();
    }
    final Requirement<File> outputFileReq = new Requirement.OutputFileRequirement("PFB (population BAF) output file must be specified.",
                                                                                  new File(defaultOutputFile));

    final RequirementSet reqSet = RequirementSetBuilder.and()
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(parseSamplesStepReq)
                                                                                 .add(sampleSubsetReq))
                                                       .add(outputFileReq);
    return new ComputePFBStep(proj, sampleSubsetReq, outputFileReq, reqSet);
  }

  final Project proj;
  final Requirement<File> sampleSubsetReq;
  final Requirement<File> outputFileReq;

  private ComputePFBStep(Project proj, Requirement<File> sampleSubReq, Requirement<File> outReq,
                         RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.noneOf(Requirement.Flag.class));
    this.proj = proj;
    this.sampleSubsetReq = sampleSubReq;
    this.outputFileReq = outReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    String setSubSampFile = proj.SAMPLE_SUBSET_FILENAME.getValue();
    String subSampFile = variables.get(sampleSubsetReq).getAbsolutePath();
    String setPFBFile = proj.CUSTOM_PFB_FILENAME.getValue();
    String pfbOutputFile = variables.get(outputFileReq).getAbsolutePath();

    if (!ext.verifyDirFormat(setSubSampFile).equals(subSampFile)) {
      proj.SAMPLE_SUBSET_FILENAME.setValue(subSampFile);
    }
    if (!ext.verifyDirFormat(setPFBFile).equals(pfbOutputFile)) {
      proj.CUSTOM_PFB_FILENAME.setValue(pfbOutputFile);
    }
  }

  @Override
  public void run(Variables variables) {
    org.genvisis.cnv.hmm.PFB.populationBAF(proj);
  }

  @Override
  public String getCommandLine(Variables variables) {
    String kvCmd = "";

    String setSubSampFile = proj.SAMPLE_SUBSET_FILENAME.getValue();
    String subSampFile = variables == null ? null
                                           : variables.get(sampleSubsetReq).getAbsolutePath();
    String setPFBFile = proj.CUSTOM_PFB_FILENAME.getValue();
    String pfbOutputFile = variables == null ? null
                                             : variables.get(outputFileReq).getAbsolutePath();

    if (subSampFile != null && !ext.verifyDirFormat(setSubSampFile).equals(subSampFile)) {
      kvCmd += " SAMPLE_SUBSET_FILENAME=" + subSampFile;
    }
    if (pfbOutputFile != null && !ext.verifyDirFormat(setPFBFile).equals(pfbOutputFile)) {
      kvCmd += " CUSTOM_PFB_FILENAME=" + pfbOutputFile;
    }

    String projPropFile = proj.getPropertyFilename();
    StringBuilder cmd = new StringBuilder();
    if (kvCmd.length() > 0) {
      cmd.append(Files.getRunString()).append(GenvisisWorkflow.PROJ_PROP_UPDATE_STR + projPropFile)
         .append(kvCmd).append("\n");
    }
    return cmd.append(Files.getRunString()).append(" ").append(PFB.class.getName()).append(" ")
              .append(CLI.ARG_PROJ).append("=").append(proj.getPropertyFilename()).append(" ")
              .append(CLI.ARG_LOG).append("=").append(proj.getLog().getFilename()).toString();
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    String subSampFile = variables.get(sampleSubsetReq).getAbsolutePath();
    String pfbOutputFile = variables.get(outputFileReq).getAbsolutePath();
    return Files.exists(pfbOutputFile) || Files.exists(ext.rootOf(subSampFile) + ".pfb");
  }
}
