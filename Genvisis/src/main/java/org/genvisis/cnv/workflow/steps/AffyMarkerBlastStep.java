package org.genvisis.cnv.workflow.steps;

import java.util.EnumSet;
import java.util.Map;
import org.genvisis.CLI;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.AffyMarkerBlast;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.common.Files;
import org.genvisis.common.ext;
import com.google.common.collect.ImmutableMap;

public class AffyMarkerBlastStep extends Step {

  public static final String NAME = "Run Marker BLAST Annotation";
  public static final String DESC = "";

  public static AffyMarkerBlastStep create(final Step parseSamplesStep, Requirement numThreadsReq) {
    final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);

    final Requirement probeFileReq = new Requirement.FileRequirement(ext.capitalizeFirst(AffyMarkerBlast.DESC_PROBE_FILE),
                                                                     AffyMarkerBlast.EXAMPLE_PROBE_FILE);
    final Requirement annotFileReq = new Requirement.FileRequirement(ext.capitalizeFirst(AffyMarkerBlast.DESC_ANNOT_FILE),
                                                                     AffyMarkerBlast.EXAMPLE_ANNOT_FILE);

    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(probeFileReq).add(annotFileReq)
                                                       .add(numThreadsReq);
    return new AffyMarkerBlastStep(probeFileReq, annotFileReq, numThreadsReq, reqSet);
  }

  final Requirement probeFileReq;
  final Requirement annotFileReq;
  final Requirement numThreadsReq;

  private AffyMarkerBlastStep(Requirement probeFileReq, Requirement annotFileReq,
                              Requirement numThreadsReq, RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.of(Requirement.Flag.MEMORY, Requirement.Flag.RUNTIME,
                                         Requirement.Flag.MULTITHREADED));
    this.probeFileReq = probeFileReq;
    this.annotFileReq = annotFileReq;
    this.numThreadsReq = numThreadsReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Project proj, Map<Requirement, String> variables) {
    // Not necessary for this step

  }

  @Override
  public void run(Project proj, Map<Requirement, String> variables) {
    String annotFile = variables.get(annotFileReq);
    String probeFile = variables.get(probeFileReq);
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
    new AffyMarkerBlast(proj, numThreads, probeFile, annotFile).blastEm();
  }

  @Override
  public String getCommandLine(Project proj, Map<Requirement, String> variables) {
    String annotFile = variables.get(annotFileReq);
    String probeFile = variables.get(probeFileReq);
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
    ImmutableMap.Builder<String, String> argsBuilder = ImmutableMap.builder();
    argsBuilder.put(CLI.ARG_PROJ, proj.getPropertyFilename());
    argsBuilder.put(AffyMarkerBlast.ARG_PROBE_FILE, probeFile);
    argsBuilder.put(AffyMarkerBlast.ARG_ANNOT_FILE, annotFile);
    argsBuilder.put(CLI.ARG_THREADS, String.valueOf(numThreads));
    return Files.getRunString() + " " + CLI.formCmdLine(AffyMarkerBlast.class, argsBuilder.build());
  }

  @Override
  public boolean checkIfOutputExists(Project proj, Map<Requirement, String> variables) {
    return Files.exists(proj.BLAST_ANNOTATION_FILENAME.getValue());
  }

}
