package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.MarkerBlastQC;
import org.genvisis.cnv.qc.SexChecks;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.Variables;
import org.pankratzlab.common.Files;

public class SexChecksStep extends Step {

  public static final String NAME = "Run Sex Checks";
  public static final String DESC = "";
  public static final String NO_CROSS_HYBE_REQUIREMENT = "Use only X and Y chromosome R values to identify sex discriminating markers";
  public static final String ADD_ESTSEX_TO_SAMPDATA_REQUIREMENT = "Add Estimated Sex to Sample Data";

  public static SexChecksStep create(final Project proj, final Step markerBlastStep,
                                     final Step sampleDataStep, final Step transposeStep,
                                     final Step sampleQCStep) {
    final Requirement<Step> sampleDataStepReq = new Requirement.StepRequirement(sampleDataStep);
    final Requirement<Step> transposeStepReq = new Requirement.StepRequirement(transposeStep);
    final Requirement<Step> sampleQCStepReq = new Requirement.StepRequirement(sampleQCStep);
    final Requirement<Boolean> addToSampleDataReq = new Requirement.OptionalBoolRequirement(ADD_ESTSEX_TO_SAMPDATA_REQUIREMENT,
                                                                                            true);

    final Requirement<File> oneHittersReq = new Requirement.FileRequirement("List of markers that do not cross hybridize",
                                                                            new File(MarkerBlastQC.defaultOneHitWondersFilename(proj.BLAST_ANNOTATION_FILENAME.getValue())));
    final Requirement<Step> markerBlastStepReq = new Requirement.StepRequirement(markerBlastStep);
    final Requirement<Boolean> noCrossHybeReq = new Requirement.BoolRequirement(NO_CROSS_HYBE_REQUIREMENT,
                                                                                false);

    final RequirementSet reqSet = RequirementSetBuilder.and().add(sampleDataStepReq)
                                                       .add(transposeStepReq).add(sampleQCStepReq)
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(oneHittersReq)
                                                                                 .add(markerBlastStepReq)
                                                                                 .add(noCrossHybeReq));
    return new SexChecksStep(proj, addToSampleDataReq, noCrossHybeReq, oneHittersReq, reqSet);
  }

  final Project proj;
  final Requirement<Boolean> addToSampleDataReq;
  final Requirement<File> oneHittersReq;
  final Requirement<Boolean> noCrossHybeReq;

  public SexChecksStep(Project proj, Requirement<Boolean> addToSD, Requirement<Boolean> noCross,
                       Requirement<File> oneHit, RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.noneOf(Requirement.Flag.class));
    this.proj = proj;
    this.addToSampleDataReq = addToSD;
    this.noCrossHybeReq = noCross;
    this.oneHittersReq = oneHit;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    // Nothing to do here
  }

  @Override
  public void run(Variables variables) {
    proj.getLog().report("Running SexCheck");
    boolean addToSampleData = variables.get(addToSampleDataReq);
    String discriminatingMarkersFile;
    if (variables.get(noCrossHybeReq)) {
      discriminatingMarkersFile = null;
    } else {
      discriminatingMarkersFile = variables.get(oneHittersReq).getAbsolutePath();
      if (!Files.exists(discriminatingMarkersFile)) {
        MarkerBlastQC.getOneHitWonders(proj, proj.BLAST_ANNOTATION_FILENAME.getValue(),
                                       discriminatingMarkersFile, 0.8, proj.getLog());
      }
    }
    org.genvisis.cnv.qc.SexChecks.sexCheck(proj, addToSampleData, discriminatingMarkersFile);
  }

  @Override
  public String getCommandLine(Variables variables) {
    String projPropFile = proj.getPropertyFilename();
    StringBuilder cmd = new StringBuilder();
    boolean addToSampleData = variables.get(addToSampleDataReq);
    File discriminatingMarkersFile;
    if (variables.get(noCrossHybeReq)) {
      discriminatingMarkersFile = null;
    } else {
      discriminatingMarkersFile = variables.get(oneHittersReq);
      if (!Files.exists(discriminatingMarkersFile)) {
        discriminatingMarkersFile = new File(MarkerBlastQC.defaultOneHitWondersFilename(proj.BLAST_ANNOTATION_FILENAME.getValue()));
        cmd.append(Files.getRunString())
           .append(" " + MarkerBlastQC.class.getName() + " proj=" + projPropFile + " blastVCF="
                   + proj.BLAST_ANNOTATION_FILENAME.getValue())
           .append("\n");
      }
    }
    return cmd.append(Files.getRunString())
              .append(" " + SexChecks.class.getName() + " -check proj=" + projPropFile).toString()
           + (discriminatingMarkersFile == null ? "" : " useMarkers=" + discriminatingMarkersFile)
           + (addToSampleData ? "" : " -skipSampleData");
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    return Files.exists(proj.SEXCHECK_RESULTS_FILENAME.getValue());
  }

}
