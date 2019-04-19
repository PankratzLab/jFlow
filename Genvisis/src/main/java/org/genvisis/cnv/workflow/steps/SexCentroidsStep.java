package org.genvisis.cnv.workflow.steps;

import java.util.EnumSet;

import org.genvisis.cnv.filesys.Centroids;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.cnv.workflow.Variables;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.PSF;

public class SexCentroidsStep extends Step {

  public static final String NAME = "Create Sex-Specific Centroids; Filter PFB file";
  public static final String DESC = "";

  public static SexCentroidsStep create(Project proj, Step pfbStep,
                                        Requirement<Integer> numThreadsReq) {
    return new SexCentroidsStep(proj, pfbStep, numThreadsReq);
  }

  private final Project proj;
  private final Requirement<Integer> numThreadsReq;

  private SexCentroidsStep(Project proj, Step pfbStep, Requirement<Integer> numThreadsReq) {
    super(NAME, DESC,
          RequirementSetBuilder.and().add(numThreadsReq)
                               .add(new Requirement.StepRequirement(pfbStep)),
          EnumSet.of(Requirement.Flag.RUNTIME, Requirement.Flag.MULTITHREADED));
    this.numThreadsReq = numThreadsReq;
    this.proj = proj;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
    GenvisisWorkflow.maybeSetProjNumThreads(proj, numThreads);
  }

  @Override
  public void run(Variables variables) {
    String malePFB;
    String femalePFB;
    String centFilePathM;
    String centFilePathF;
    String outputDir = proj.DATA_DIRECTORY.getValue();
    malePFB = outputDir + "males.pfb";
    femalePFB = outputDir + "females.pfb";
    centFilePathM = outputDir + "sexSpecific_Male.cent";
    centFilePathF = outputDir + "sexSpecific_Female.cent";

    int numThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
    Centroids.computeSexSpecificCentroids(proj, new String[] {malePFB, femalePFB},
                                          new String[] {centFilePathM, centFilePathF}, numThreads);

  }

  @Override
  public String getCommandLine(Variables variables) {
    int numThreads = StepBuilder.resolveThreads(proj,
                                                variables == null ? -1
                                                                  : variables.get(numThreadsReq));
    String mainCmd = Files.getRunString() + " cnv.filesys.Centroids proj="
                     + proj.getPropertyFilename() + " -sexSpecific " + PSF.Ext.NUM_THREADS_COMMAND
                     + numThreads;
    return mainCmd;
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    String malePFB;
    String femalePFB;
    String centFilePathM;
    String centFilePathF;
    String outputDir = proj.DATA_DIRECTORY.getValue();
    malePFB = outputDir + "males.pfb";
    femalePFB = outputDir + "females.pfb";
    centFilePathM = outputDir + "sexSpecific_Male.cent";
    centFilePathF = outputDir + "sexSpecific_Female.cent";
    boolean exists = Files.exists(malePFB);
    exists = exists && Files.exists(femalePFB);
    exists = exists && Files.exists(centFilePathM);
    exists = exists && Files.exists(centFilePathF);
    return exists;
  }

}
