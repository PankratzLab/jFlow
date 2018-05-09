package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;
import org.genvisis.CLI;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.IlluminaMarkerBlast;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.cnv.workflow.Variables;
import org.genvisis.common.Files;
import org.genvisis.common.ext;
import com.google.common.collect.ImmutableMap;

public class IlluminaMarkerBlastStep extends Step {

  public static final String NAME = "Run Marker BLAST Annotation";
  public static final String DESC = "";

  public static IlluminaMarkerBlastStep create(Project proj,
                                               final ParseSamplesStep parseSamplesStep,
                                               Requirement<Integer> numThreadsReq) {
    final Requirement<Step> parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    final Requirement<File> manifestFileReq = new Requirement.FileRequirement(ext.capitalizeFirst(IlluminaMarkerBlast.DESC_MANIFEST)
                                                                              + "  (e.g. "
                                                                              + IlluminaMarkerBlast.EXAMPLE_MANIFEST
                                                                              + ")", new File(""));

    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(manifestFileReq).add(numThreadsReq);
    return new IlluminaMarkerBlastStep(proj, manifestFileReq, numThreadsReq, reqSet);
  }

  private IlluminaMarkerBlastStep(Project proj, Requirement<File> manifestFileReq,
                                  Requirement<Integer> numThreadsReq, RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.of(Requirement.Flag.MEMORY, Requirement.Flag.RUNTIME,
                                         Requirement.Flag.MULTITHREADED));
    this.proj = proj;
    this.manifestFileReq = manifestFileReq;
    this.numThreadsReq = numThreadsReq;
  }

  final Project proj;
  final Requirement<File> manifestFileReq;
  final Requirement<Integer> numThreadsReq;

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    // Not necessary for this step

  }

  @Override
  public void run(Variables variables) {
    String manifestFile = variables.get(manifestFileReq).getAbsolutePath();
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
    new IlluminaMarkerBlast(proj, numThreads, manifestFile).blastEm();
  }

  @Override
  public String getCommandLine(Variables variables) {
    String manifestFile = variables.get(manifestFileReq).getAbsolutePath();
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
    ImmutableMap.Builder<String, String> argsBuilder = ImmutableMap.builder();
    argsBuilder.put(CLI.ARG_PROJ, proj.getPropertyFilename());
    argsBuilder.put(IlluminaMarkerBlast.ARG_MANIFEST, manifestFile);
    argsBuilder.put(CLI.ARG_THREADS, String.valueOf(numThreads));
    return Files.getRunString() + " "
           + CLI.formCmdLine(IlluminaMarkerBlast.class, argsBuilder.build());
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    return Files.exists(proj.BLAST_ANNOTATION_FILENAME.getValue());
  }
}
