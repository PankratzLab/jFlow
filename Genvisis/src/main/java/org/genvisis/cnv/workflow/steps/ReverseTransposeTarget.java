package org.genvisis.cnv.workflow.steps;

import java.io.IOException;
import java.util.EnumSet;
import org.genvisis.CLI;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.manage.TempFileTranspose;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.Variables;
import org.genvisis.common.Files;

public class ReverseTransposeTarget extends Step {

  public static final String NAME = "Reverse-Transpose Marker Files to Sample Files";
  public static final String DESC = "";

  public static ReverseTransposeTarget create(Project proj, Step affyCelParsingStep) {
    return new ReverseTransposeTarget(proj, affyCelParsingStep);
  }

  final Project proj;

  private ReverseTransposeTarget(Project proj, Step affyCelParsingStep) {
    super(NAME, DESC,
          RequirementSetBuilder.and().add(new Requirement.StepRequirement(affyCelParsingStep)),
          EnumSet.of(Requirement.Flag.MEMORY, Requirement.Flag.MULTITHREADED,
                     Requirement.Flag.RUNTIME));
    this.proj = proj;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    // Nothing to do here
  }

  @Override
  public void run(Variables variables) {
    String temp = proj.PROJECT_DIRECTORY.getValue() + "temp/";
    TempFileTranspose tft = new TempFileTranspose(proj, temp, "");
    tft.setupMarkerListFile();
    try {
      tft.runFirst();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
    tft.setupSampleListFile();
    try {
      tft.runSecond();
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  @Override
  public String getCommandLine(Variables variables) {
    String projPropFile = proj.getPropertyFilename();
    StringBuilder cmd1 = new StringBuilder();
    cmd1.append(Files.getRunString());
    cmd1.append(" ").append(TempFileTranspose.class.getName());
    cmd1.append(CLI.ARG_PROJ).append("=").append(projPropFile);
    cmd1.append(" type=M jobID=$PBS_JOBID -setup");
    cmd1.append("\n");
    cmd1.append(Files.getRunString());
    cmd1.append(" ").append(TempFileTranspose.class.getName());
    cmd1.append(CLI.ARG_PROJ).append("=").append(projPropFile);
    cmd1.append(" type=S jobID=$PBS_JOBID -setup");

    return cmd1.toString();
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    return Files.exists(proj.SAMPLE_DIRECTORY.getValue())
           && Files.countFiles(proj.SAMPLE_DIRECTORY.getValue(false, false),
                               Sample.SAMPLE_FILE_EXTENSION) == proj.getSamples().length;
  }

}
