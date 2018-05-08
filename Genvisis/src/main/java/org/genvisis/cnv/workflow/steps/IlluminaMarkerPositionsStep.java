package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.qc.MarkerBlast;
import org.genvisis.cnv.qc.MarkerBlast.FILE_SEQUENCE_TYPE;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.Variables;
import org.genvisis.common.Files;

public class IlluminaMarkerPositionsStep extends Step {

  public static final String NAME = "Create Marker Positions (if not already exists)";
  public static final String DESC = "";

  public static IlluminaMarkerPositionsStep create(Project proj) {
    final Requirement<File> manifestReq = new Requirement.FileRequirement("An Illumina Manifest file.",
                                                                          new File(proj.getLocationOfSNP_Map(false)));

    return new IlluminaMarkerPositionsStep(manifestReq);
  }

  final Requirement<File> manifestReq;

  private IlluminaMarkerPositionsStep(Requirement<File> manifestReq) {
    super(NAME, DESC, RequirementSetBuilder.and().add(manifestReq),
          EnumSet.noneOf(Requirement.Flag.class));
    this.manifestReq = manifestReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Project proj, Variables variables) {
    // not needed for step
  }

  @Override
  public void run(Project proj, Variables variables) {
    proj.getLog().report("Generating marker positions file");
    String manifest = variables.get(manifestReq).getAbsolutePath();
    if (Files.exists(manifest)) {
      MarkerBlast.extractMarkerPositionsFromManifest(manifest, ARRAY.ILLUMINA,
                                                     FILE_SEQUENCE_TYPE.MANIFEST_FILE,
                                                     proj.MARKER_POSITION_FILENAME.getValue(false,
                                                                                            false),
                                                     Files.determineDelimiter(manifest,
                                                                              proj.getLog()),
                                                     proj.getLog());
    }
  }

  @Override
  public boolean checkIfOutputExists(Project proj, Variables variables) {
    return Files.exists(proj.MARKER_POSITION_FILENAME.getValue(false, false));
  }

  @Override
  public String getCommandLine(Project proj, Variables variables) {
    String projFile = proj.getPropertyFilename();
    File manifest = variables.get(manifestReq);
    String baseCommand = Files.getRunString() + " cnv.manage.Markers proj=" + projFile;
    return baseCommand + " snps=" + manifest + " -manifest";
  }

}
