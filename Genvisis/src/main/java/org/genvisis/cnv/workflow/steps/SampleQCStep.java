package org.genvisis.cnv.workflow.steps;

import java.util.EnumSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.LrrSd;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.Requirement.StepRequirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.cnv.workflow.Variables;
import org.genvisis.common.Files;
import org.genvisis.common.PSF;

public class SampleQCStep extends Step {

  public static final String NAME = "Run Sample QC Metrics";
  public static final String DESC = "";

  private Requirement<Integer> numThreadsReq;

  public static SampleQCStep create(Step parseSamplesStep, Requirement<Integer> numThreadsReq) {
    RequirementSet reqSet = RequirementSetBuilder.and().add(new StepRequirement(parseSamplesStep))
                                                 .add(numThreadsReq);
    return new SampleQCStep(reqSet);
  }

  private SampleQCStep(RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.of(Requirement.Flag.MULTITHREADED));
  }

  @Override
  public void setNecessaryPreRunProperties(Project proj, Variables variables) {
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
    GenvisisWorkflow.maybeSetProjNumThreads(proj, numThreads);
  }

  @Override
  public void run(Project proj, Variables variables) {
    proj.getLog().report("Running LrrSd");
    int numThreads = proj.NUM_THREADS.getValue();
    LrrSd.init(proj, null, null, numThreads, false);
  }

  @Override
  public String getCommandLine(Project proj, Variables variables) {
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
    String projPropFile = proj.getPropertyFilename();
    StringBuilder cmd = new StringBuilder();
    cmd.append(Files.getRunString()).append(" cnv.qc.LrrSd").append(" proj=").append(projPropFile)
       .append(" ").append(PSF.Ext.NUM_THREADS_COMMAND).append(numThreads)
       .append(" projectMarkers=TRUE");
    return cmd.toString();
  }

  @Override
  public boolean checkIfOutputExists(Project proj, Variables variables) {
    return Files.exists(proj.SAMPLE_QC_FILENAME.getValue(false, false));
  }
}
