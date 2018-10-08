package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.AffyMarkerBlast;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.cnv.workflow.Variables;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.CLI;
import com.google.common.collect.ImmutableMap;

public class AffyMarkerBlastStep extends Step {

  public static final String NAME = "Run Marker BLAST Annotation";
  public static final String DESC = "";

  public static AffyMarkerBlastStep create(Project proj, final Step parseSamplesStep,
                                           Requirement<Integer> numThreadsReq) {
    final Requirement<Step> parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);

    final Requirement<File> probeFileReq = new Requirement.FileRequirement(ext.capitalizeFirst(AffyMarkerBlast.DESC_PROBE_FILE),
                                                                           new File(AffyMarkerBlast.EXAMPLE_PROBE_FILE));
    final Requirement<File> annotFileReq = new Requirement.FileRequirement(ext.capitalizeFirst(AffyMarkerBlast.DESC_ANNOT_FILE),
                                                                           new File(AffyMarkerBlast.EXAMPLE_ANNOT_FILE));

    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(probeFileReq).add(annotFileReq)
                                                       .add(numThreadsReq);
    return new AffyMarkerBlastStep(proj, probeFileReq, annotFileReq, numThreadsReq, reqSet);
  }

  final Project proj;
  final Requirement<File> probeFileReq;
  final Requirement<File> annotFileReq;
  final Requirement<Integer> numThreadsReq;

  private AffyMarkerBlastStep(Project proj, Requirement<File> probeFileReq,
                              Requirement<File> annotFileReq, Requirement<Integer> numThreadsReq,
                              RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.of(Requirement.Flag.MEMORY, Requirement.Flag.RUNTIME,
                                         Requirement.Flag.MULTITHREADED));
    this.proj = proj;
    this.probeFileReq = probeFileReq;
    this.annotFileReq = annotFileReq;
    this.numThreadsReq = numThreadsReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    // Not necessary for this step

  }

  @Override
  public void run(Variables variables) {
    String annotFile = variables.get(annotFileReq).getPath();
    String probeFile = variables.get(probeFileReq).getPath();
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
    new AffyMarkerBlast(proj, numThreads, probeFile, annotFile).blastEm();
  }

  @Override
  public String getCommandLine(Variables variables) {
    String annotFile = variables.get(annotFileReq).getPath();
    String probeFile = variables.get(probeFileReq).getPath();
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
    ImmutableMap.Builder<String, String> argsBuilder = ImmutableMap.builder();
    argsBuilder.put(CLI.ARG_PROJ, proj.getPropertyFilename());
    argsBuilder.put(AffyMarkerBlast.ARG_PROBE_FILE, probeFile);
    argsBuilder.put(AffyMarkerBlast.ARG_ANNOT_FILE, annotFile);
    argsBuilder.put(CLI.ARG_THREADS, String.valueOf(numThreads));
    return Files.getRunString() + " " + CLI.formCmdLine(AffyMarkerBlast.class, argsBuilder.build());
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    return Files.exists(proj.BLAST_ANNOTATION_FILENAME.getValue());
  }

}
