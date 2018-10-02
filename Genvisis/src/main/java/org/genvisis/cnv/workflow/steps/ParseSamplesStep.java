package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ext;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.cnv.workflow.Variables;

public class ParseSamplesStep extends Step {

  public static final String NAME = "Parse Sample Files";
  public static final String DESC = "";

  public static ParseSamplesStep create(Project proj, final Step markerPositionsStep,
                                        Requirement<Integer> numThreadsReq) {
    final Requirement<File> markerPositionsReq = new Requirement.FileRequirement("Marker Positions file must already exist.",
                                                                                 new File(proj.MARKER_POSITION_FILENAME.getValue(false,
                                                                                                                                 false)));
    final RequirementSet reqSet = RequirementSetBuilder.and();
    if (markerPositionsStep == null) {
      reqSet.add(markerPositionsReq).add(numThreadsReq);
    } else {
      final Requirement<Step> markerPositionsStepReq = new Requirement.StepRequirement(markerPositionsStep);
      reqSet.add(RequirementSetBuilder.or().add(markerPositionsReq).add(markerPositionsStepReq))
            .add(numThreadsReq);
    }
    return new ParseSamplesStep(proj, markerPositionsReq, numThreadsReq, reqSet);
  }

  final Project proj;
  final Requirement<File> markerPositionsReq;
  final Requirement<Integer> numThreadsReq;

  private ParseSamplesStep(Project proj, Requirement<File> markerPosReq,
                           Requirement<Integer> numThreadsReq, RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.of(Requirement.Flag.MEMORY, Requirement.Flag.RUNTIME,
                                         Requirement.Flag.MULTITHREADED));
    this.proj = proj;
    this.markerPositionsReq = markerPosReq;
    this.numThreadsReq = numThreadsReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    String projFile = proj.MARKER_POSITION_FILENAME.getValue(false, false);
    String mkrFile = variables.get(markerPositionsReq).getAbsolutePath();
    mkrFile = ext.verifyDirFormat(mkrFile);
    mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
    if (!mkrFile.equals(projFile)) {
      proj.MARKER_POSITION_FILENAME.setValue(mkrFile);
    }
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
    GenvisisWorkflow.maybeSetProjNumThreads(proj, numThreads);
  }

  @Override
  public void run(Variables variables) {
    int numThreads = proj.NUM_THREADS.getValue();
    proj.getLog().report("Parsing sample files");
    int retCode = org.genvisis.cnv.manage.SourceFileParser.createFiles(proj, numThreads);
    switch (retCode) {
      case 0:
        throw new RuntimeException("Operation failure, please check log for more information.");
      case 1:
      case 6:
      default:
        break;
    }
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    String sampleDirectory = proj.SAMPLE_DIRECTORY.getValue(false, false);
    boolean mkrSetFile = proj.MARKERSET_FILENAME.exists();
    boolean returnValue = mkrSetFile;
    returnValue = returnValue && proj.getSampleList() != null;
    returnValue = returnValue && Files.exists(sampleDirectory);

    int numSamples = returnValue ? proj.getSampleList().getSamples().length : 0;
    returnValue = returnValue && numSamples > 0;
    returnValue = returnValue
                  && Files.countFiles(sampleDirectory, Sample.SAMPLE_FILE_EXTENSION) == numSamples;
    // checking the validity / completeness of each sample would be a Good Thing, but too
    // costly time-wise for larger projects
    return returnValue;
  }

  @Override
  public String getCommandLine(Variables variables) {
    String projPropFile = proj.getPropertyFilename();
    StringBuilder kvCmd = new StringBuilder(Files.getRunString()).append(GenvisisWorkflow.PROJ_PROP_UPDATE_STR)
                                                                 .append(projPropFile);
    StringBuilder kvPairs = new StringBuilder();
    String projFile = proj.MARKER_POSITION_FILENAME.getValue(false, false);
    String mkrFile = variables.get(markerPositionsReq).getAbsolutePath();
    mkrFile = ext.verifyDirFormat(mkrFile);
    mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
    if (!mkrFile.equals(projFile)) {
      kvPairs.append(" MARKER_POSITION_FILENAME=").append(mkrFile);
    }
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
    if (numThreads != proj.NUM_THREADS.getValue()) {
      kvPairs.append(" ").append(proj.NUM_THREADS.getName()).append("=").append(numThreads);
    }
    StringBuilder command = new StringBuilder();
    if (kvPairs.length() != 0) {
      command.append(kvCmd).append(kvPairs).append("\n");
    }
    command.append(Files.getRunString()).append(" cnv.manage.SourceFileParser proj=")
           .append(projPropFile).append(" ").append(PSF.Ext.NUM_THREADS_COMMAND).append(numThreads);
    return command.toString();
  }

}
