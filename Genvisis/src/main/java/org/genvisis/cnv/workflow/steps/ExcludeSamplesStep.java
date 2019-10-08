package org.genvisis.cnv.workflow.steps;

import java.awt.Font;
import java.util.EnumSet;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.SwingConstants;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.LrrSd;
import org.genvisis.cnv.qc.SampleQC;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepAssist;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.cnv.workflow.Variables;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.ext;

public class ExcludeSamplesStep extends Step implements StepAssist {

  public static final String NAME = "Identify Excluded Samples";
  public static final String DESC = "";

  public static ExcludeSamplesStep create(Project proj, final Step sampleQCStep) {
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

    return new ExcludeSamplesStep(proj, sampleQCStepReq, lrrSdThresholdReq, callrateThresholdReq);
  }

  final Project proj;
  final Requirement<Double> callrateThresholdReq;
  final Requirement<Double> lrrSdThresholdReq;

  public ExcludeSamplesStep(Project proj, Requirement<Step> sampleQCStepReq,
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
    ExcludeSamplesStep step = ExcludeSamplesStep.create(proj, samplesStep);
    Variables variables = step.parseArguments(args);
    Step.run(proj, step, variables);
  }

  private JLabel assistLabel;

  @Override
  public JComponent getStepAssistComponent() {
    if (assistLabel == null) {
      assistLabel = new JLabel();
      assistLabel.setFont(assistLabel.getFont().deriveFont(Font.PLAIN, 12));
      assistLabel.setHorizontalAlignment(SwingConstants.CENTER);
      assistLabel.setHorizontalTextPosition(SwingConstants.CENTER);
    }
    return assistLabel;
  }

  @Override
  public void updateStepAssistComponent(Variables variables) {
    double projLrrSdThreshold = proj.LRRSD_CUTOFF.getValue();
    double projCallrateThreshold = proj.SAMPLE_CALLRATE_THRESHOLD.getValue();
    setNecessaryPreRunProperties(variables);
    SampleQC sampleQC = SampleQC.loadSampleQCWithoutSideEffects(proj, LrrSd.SAMPLE_COLUMN,
                                                                LrrSd.NUMERIC_COLUMNS, false,
                                                                false);
    int numExcluded = sampleQC.addExcludes();
    StringBuilder msg = new StringBuilder();
    msg.append("<html>These options will result in <b>").append(numExcluded).append("</b> ")
       .append(numExcluded > 1 ? "samples" : "sample").append(" being excluded.</html>");
    assistLabel.setText(msg.toString());
    proj.LRRSD_CUTOFF.setValue(projLrrSdThreshold);
    proj.SAMPLE_CALLRATE_THRESHOLD.setValue(projCallrateThreshold);
  }

}
