package org.genvisis.cnv.workflow.steps;

import java.util.EnumSet;
import java.util.Map;
import org.genvisis.CLI;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.IlluminaMarkerBlast;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.common.Files;
import org.genvisis.common.ext;
import com.google.common.collect.ImmutableMap;

public class IlluminaMarkerBlastStep extends Step {

  public static final String NAME = "Run Marker BLAST Annotation";
  public static final String DESC = "";

  public static IlluminaMarkerBlastStep create(Project proj, final Step parseSamplesStep,
                                               Requirement numThreadsReq) {
    final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    final Requirement manifestFileReq = new Requirement.FileRequirement(ext.capitalizeFirst(IlluminaMarkerBlast.DESC_MANIFEST),
                                                                        IlluminaMarkerBlast.EXAMPLE_MANIFEST);

    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(manifestFileReq).add(numThreadsReq);
    return new IlluminaMarkerBlastStep(manifestFileReq, numThreadsReq, reqSet);
  }

  private IlluminaMarkerBlastStep(Requirement manifestFileReq, Requirement numThreadsReq,
                                  RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.of(Requirement.Flag.MEMORY, Requirement.Flag.RUNTIME,
                                         Requirement.Flag.MULTITHREADED));
    this.manifestFileReq = manifestFileReq;
    this.numThreadsReq = numThreadsReq;
  }

  final Requirement manifestFileReq;
  final Requirement numThreadsReq;

  @Override
  public void setNecessaryPreRunProperties(Project proj,
                                           Map<Step, Map<Requirement, String>> variables) {
    // Not necessary for this step

  }

  @Override
  public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String manifestFile = variables.get(this).get(manifestFileReq);
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(this).get(numThreadsReq));
    new IlluminaMarkerBlast(proj, numThreads, manifestFile).blastEm();
  }

  @Override
  public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String manifestFile = variables.get(this).get(manifestFileReq);
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(this).get(numThreadsReq));
    ImmutableMap.Builder<String, String> argsBuilder = ImmutableMap.builder();
    argsBuilder.put(CLI.ARG_PROJ, proj.getPropertyFilename());
    argsBuilder.put(IlluminaMarkerBlast.ARG_MANIFEST, manifestFile);
    argsBuilder.put(CLI.ARG_THREADS, String.valueOf(numThreads));
    return Files.getRunString() + " "
           + CLI.formCmdLine(IlluminaMarkerBlast.class, argsBuilder.build());
  }

  @Override
  public boolean checkIfOutputExists(Project proj, Map<Step, Map<Requirement, String>> variables) {
    return Files.exists(proj.BLAST_ANNOTATION_FILENAME.getValue());
  }
}
