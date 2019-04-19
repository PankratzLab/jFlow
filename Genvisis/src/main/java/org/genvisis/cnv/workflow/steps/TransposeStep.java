package org.genvisis.cnv.workflow.steps;

import java.util.EnumSet;

import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.Variables;
import org.pankratzlab.common.Files;

public class TransposeStep extends Step {

  public static final String NAME = "Transpose Data into Marker-Dominant Files";
  public static final String DESC = "";

  public static TransposeStep create(Project proj, Step parseSamplesStep) {
    return new TransposeStep(proj, parseSamplesStep);
  }

  final Project proj;

  private TransposeStep(Project proj, Step parseSamplesStep) {
    super(NAME, DESC,
          RequirementSetBuilder.and().add(new Requirement.StepRequirement(parseSamplesStep)),
          EnumSet.of(Requirement.Flag.MEMORY));
    this.proj = proj;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    // Nothing to do here
  }

  @Override
  public void run(Variables variables) {
    proj.getLog().report("Transposing data");
    TransposeData.transposeData(proj, 2000000000, false); // compact if no LRR was provided
  }

  @Override
  public String getCommandLine(Variables variables) {
    String projPropFile = proj.getPropertyFilename();
    StringBuilder cmd = new StringBuilder();
    return cmd.append(Files.getRunString()).append(" cnv.manage.TransposeData -transpose proj="
                                                   + projPropFile + " max=" + 2000000000)
              .toString();
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    return Files.countFiles(proj.MARKER_DATA_DIRECTORY.getValue(false, false),
                            MarkerData.MARKER_DATA_FILE_EXTENSION) > 0;
  }

}
