package org.genvisis.cnv.workflow.steps;

import java.util.EnumSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.Variables;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;

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
    if (!Files.exists(proj.MARKER_DATA_DIRECTORY.getValue() + TransposeData.TEMP_SAMPLES_FILE)) {
      String[] samples = proj.getSamples();
      if (samples == null) {
        samples = HashVec.loadFileToStringArray(proj.PROJECT_DIRECTORY.getValue()
                                                + "ListOfSamples.txt", false, null, false);
      }
      Files.writeArray(samples,
                       proj.MARKER_DATA_DIRECTORY.getValue() + TransposeData.TEMP_SAMPLES_FILE);
    }
    TransposeData.reverseTransposeStreaming(proj,
                                            ext.replaceWithLinuxSafeCharacters(proj.PROJECT_NAME.getValue()));
  }

  @Override
  public String getCommandLine(Variables variables) {
    String projPropFile = proj.getPropertyFilename();
    StringBuilder cmd = new StringBuilder();
    return cmd.append(Files.getRunString())
              .append(" " + TransposeData.class.getName() + " -reverseStream proj=" + projPropFile)
              .toString();
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    return Files.exists(proj.SAMPLE_DIRECTORY.getValue())
           && Files.countFiles(proj.SAMPLE_DIRECTORY.getValue(false, false),
                               Sample.SAMPLE_FILE_EXTENSION) == proj.getSamples().length;
  }

}
