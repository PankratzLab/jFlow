package org.genvisis.cnv.workflow.steps;

import java.util.EnumSet;
import java.util.Map;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.qc.MarkerBlast;
import org.genvisis.cnv.qc.MarkerBlast.FILE_SEQUENCE_TYPE;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.common.Files;

public class IlluminaMarkerPositionsStep extends Step {

  public static final String NAME = "Create Marker Positions (if not already exists)";
  public static final String DESC = "";

  public static IlluminaMarkerPositionsStep create(Project proj, double priority) {
    final Requirement snpMapReq = new Requirement.FileRequirement("An Illumina SNP_map file.",
                                                                  proj.getLocationOfSNP_Map(false));
    final Requirement manifestReq = new Requirement.FileRequirement("An Illumina Manifest file.",
                                                                    proj.getLocationOfSNP_Map(false));

    final RequirementSet reqSet = RequirementSetBuilder.and()
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(snpMapReq)
                                                                                 .add(manifestReq));
    return new IlluminaMarkerPositionsStep(snpMapReq, manifestReq, reqSet, priority);
  }

  final Requirement snpMapReq;
  final Requirement manifestReq;

  private IlluminaMarkerPositionsStep(Requirement snpMapReq, Requirement manifestReq,
                                      RequirementSet reqSet, double priority) {
    super(NAME, DESC, reqSet, EnumSet.noneOf(Requirement.Flag.class), priority);
    this.snpMapReq = snpMapReq;
    this.manifestReq = manifestReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Project proj,
                                           Map<Step, Map<Requirement, String>> variables) {
    // not needed for step
  }

  @Override
  public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
    proj.getLog().report("Generating marker positions file");
    String snpMap = variables.get(this).get(snpMapReq);
    String manifest = variables.get(this).get(manifestReq);
    if (Files.exists(snpMap)) {
      org.genvisis.cnv.manage.Markers.generateMarkerPositions(proj, snpMap);
    } else if (Files.exists(manifest)) {
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
  public boolean checkIfOutputExists(Project proj, Map<Step, Map<Requirement, String>> variables) {
    return Files.exists(proj.MARKER_POSITION_FILENAME.getValue(false, false));
  }

  @Override
  public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String projFile = proj.getPropertyFilename();
    String snpMap = variables.get(this).get(snpMapReq);
    String manifest = variables.get(this).get(manifestReq);
    String baseCommand = Files.getRunString() + " cnv.manage.Markers proj=" + projFile;
    if (Files.exists(snpMap)) {
      return baseCommand + " snps=" + snpMap;
    } else {
      return baseCommand + " snps=" + manifest + " -manifest";
    }
  }

}
