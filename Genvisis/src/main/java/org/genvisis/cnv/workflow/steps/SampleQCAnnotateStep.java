package org.genvisis.cnv.workflow.steps;

import java.util.EnumSet;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.SampleQC;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.cnv.workflow.Variables;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.ext;

public class SampleQCAnnotateStep extends Step {

  public static final String NAME = "Identify Excluded Samples";
  public static final String DESC = "";

  public static SampleQCAnnotateStep create(Project proj, final Step sampleQCStep) {
    final Requirement<Step> sampleQCStepReq = new Requirement.StepRequirement(sampleQCStep);
    final Requirement<Double> lrrSdThresholdReq = new Requirement.DoubleRequirement("lrrSDThreshold",
                                                                                    "LRR SD Threshold",
                                                                                    proj.LRRSD_CUTOFF.getValue(),
                                                                                    proj.LRRSD_CUTOFF.getMinValue(),
                                                                                    proj.LRRSD_CUTOFF.getMaxValue());

    final Requirement<Double> callrateThresholdReq = new Requirement.DoubleRequirement("callrateThreshold",
                                                                                       "Callrate Threshold",
                                                                                       proj.SAMPLE_CALLRATE_THRESHOLD.getValue(),
                                                                                       proj.SAMPLE_CALLRATE_THRESHOLD.getMinValue(),
                                                                                       proj.SAMPLE_CALLRATE_THRESHOLD.getMaxValue());

    return new SampleQCAnnotateStep(proj, sampleQCStepReq, lrrSdThresholdReq, callrateThresholdReq);
  }

  final Project proj;
  final Requirement<Double> callrateThresholdReq;
  final Requirement<Double> lrrSdThresholdReq;

  public SampleQCAnnotateStep(Project proj, Requirement<Step> sampleQCStepReq,
                              Requirement<Double> lrrSdReq, Requirement<Double> callrateReq) {
    super(NAME, DESC,
          RequirementSetBuilder.and().add(sampleQCStepReq).add(lrrSdReq).add(callrateReq),
          EnumSet.noneOf(Requirement.Flag.class));
    this.proj = proj;
    this.lrrSdThresholdReq = lrrSdReq;
    this.callrateThresholdReq = callrateReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    double projLrrSdThreshold = proj.LRRSD_CUTOFF.getValue();
    double lrrSdThreshold = variables.get(lrrSdThresholdReq);
    double projCallrateThreshold = proj.SAMPLE_CALLRATE_THRESHOLD.getValue();
    double callrateThreshold = variables.get(callrateThresholdReq);

    if (projLrrSdThreshold != lrrSdThreshold) {
      proj.LRRSD_CUTOFF.setValue(lrrSdThreshold);
    }
    if (projCallrateThreshold != callrateThreshold) {
      proj.SAMPLE_CALLRATE_THRESHOLD.setValue(callrateThreshold);
    }
  }

  @Override
  public void run(Variables variables) {
    SampleQC.parseExcludes(proj);
  }

  @Override
  public String getCommandLine(Variables variables) {
    return getStepCommandLine(proj, variables);
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    String sampleDataFile = proj.SAMPLE_DATA_FILENAME.getValue();
    if (!Files.exists(sampleDataFile)) {
      return false;
    }
    String[] header = Files.getHeaderOfFile(sampleDataFile, proj.getLog());

    // These columns should always added by SampleQC
    String[] baseHeader = {SampleQC.EXCLUDE_HEADER, "ExcludeNote"};

    int[] facts = ext.indexFactors(baseHeader, header, false);
    for (int i : facts) {
      if (i == -1) {
        return false;
      }
    }

    return true;
  }

  public static void main(String[] args) {
    Project proj = Step.parseProject(args);
    StepBuilder sb = new StepBuilder(proj);
    Step samplesStep = sb.generateSamplesParsingStep();
    SampleQCAnnotateStep step = SampleQCAnnotateStep.create(proj, samplesStep);
    Variables variables = step.parseArguments(args);
    Step.run(proj, step, variables);
  }

}
