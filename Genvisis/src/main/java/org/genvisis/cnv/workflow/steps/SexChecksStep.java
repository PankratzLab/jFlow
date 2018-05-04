package org.genvisis.cnv.workflow.steps;

import java.util.EnumSet;
import java.util.Map;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.MarkerBlastQC;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.common.Files;

public class SexChecksStep extends Step {

  public static final String NAME = "";
  public static final String DESC = "";

  public static SexChecksStep create(final Project proj, final Step parseSamplesStep,
                                     final Step markerBlastStep, final Step sampleDataStep,
                                     final Step transposeStep, final Step sampleQCStep,
                                     final double priority) {
    final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    final Requirement sampleDataStepReq = new Requirement.StepRequirement(sampleDataStep);
    final Requirement transposeStepReq = new Requirement.StepRequirement(transposeStep);
    final Requirement sampleQCStepReq = new Requirement.StepRequirement(sampleQCStep);
    final Requirement addToSampleDataReq = new Requirement.OptionalBoolRequirement("Add Estimated Sex to Sample Data",
                                                                                   true);

    final Requirement oneHittersReq = new Requirement.FileRequirement("List of markers that do not cross hybridize",
                                                                      MarkerBlastQC.defaultOneHitWondersFilename(proj.BLAST_ANNOTATION_FILENAME.getValue()));
    final Requirement markerBlastStepReq = new Requirement.StepRequirement(markerBlastStep);
    final Requirement noCrossHybeReq = new Requirement.BoolRequirement("Use only X and Y chromosome R values to identify sex discriminating markers",
                                                                       false);

    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(sampleDataStepReq).add(transposeStepReq)
                                                       .add(sampleQCStepReq)
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(oneHittersReq)
                                                                                 .add(markerBlastStepReq)
                                                                                 .add(noCrossHybeReq));
    return new SexChecksStep(addToSampleDataReq, noCrossHybeReq, oneHittersReq, reqSet, priority);
  }

  final Requirement addToSampleDataReq;
  final Requirement oneHittersReq;
  final Requirement noCrossHybeReq;

  public SexChecksStep(Requirement addToSD, Requirement noCross, Requirement oneHit,
                       RequirementSet reqSet, double priority) {
    super(NAME, DESC, reqSet, EnumSet.noneOf(Requirement.Flag.class), priority);
    this.addToSampleDataReq = addToSD;
    this.noCrossHybeReq = noCross;
    this.oneHittersReq = oneHit;
  }

  @Override
  public void setNecessaryPreRunProperties(Project proj,
                                           Map<Step, Map<Requirement, String>> variables) {
    // Nothing to do here
  }

  @Override
  public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
    proj.getLog().report("Running SexCheck");
    boolean addToSampleData = Boolean.parseBoolean(variables.get(this).get(addToSampleDataReq));
    String discriminatingMarkersFile;
    if (Boolean.parseBoolean(variables.get(this).get(noCrossHybeReq))) {
      discriminatingMarkersFile = null;
    } else {
      discriminatingMarkersFile = variables.get(this).get(oneHittersReq);
      if (!Files.exists(discriminatingMarkersFile)) {
        MarkerBlastQC.getOneHitWonders(proj, proj.BLAST_ANNOTATION_FILENAME.getValue(),
                                       discriminatingMarkersFile, 0.8, proj.getLog());
      }
    }
    org.genvisis.cnv.qc.SexChecks.sexCheck(proj, addToSampleData, discriminatingMarkersFile);
  }

  @Override
  public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String projPropFile = proj.getPropertyFilename();
    StringBuilder cmd = new StringBuilder();
    boolean addToSampleData = Boolean.parseBoolean(variables.get(this).get(addToSampleDataReq));
    String discriminatingMarkersFile;
    if (Boolean.parseBoolean(variables.get(this).get(noCrossHybeReq))) {
      discriminatingMarkersFile = null;
    } else {
      discriminatingMarkersFile = variables.get(this).get(oneHittersReq);
      if (!Files.exists(discriminatingMarkersFile)) {
        cmd.append(Files.getRunString())
           .append(" cnv.qc.MarkerBlastQC proj=" + projPropFile + " blastVCF="
                   + proj.BLAST_ANNOTATION_FILENAME.getValue())
           .append("\n");
      }
    }
    return cmd.append(Files.getRunString()).append(" cnv.qc.SexChecks -check proj=" + projPropFile)
              .toString()
           + (discriminatingMarkersFile == null ? "" : " useMarkers=" + discriminatingMarkersFile)
           + (addToSampleData ? "" : " -skipSampleData");
  }

  @Override
  public boolean checkIfOutputExists(Project proj, Map<Step, Map<Requirement, String>> variables) {
    return Files.exists(proj.SEXCHECK_RESULTS_FILENAME.getValue());
  }

}
