package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;

import org.genvisis.cnv.Resources;
import org.genvisis.cnv.affy.APTAxiomPipeline;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.AffyMarkerBlast;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.Requirement.BoolRequirement;
import org.genvisis.cnv.workflow.Requirement.DirRequirement;
import org.genvisis.cnv.workflow.Requirement.FileRequirement;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.Variables;
import org.pankratzlab.common.Elision;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ext;

public class AxiomCELProcessingStep extends Step {

  public static final String NAME = "Process Affymetrix Axiom CEL files";
  public static final String DESC = "";

  private final Project proj;

  private final Requirement<Integer> numThreadsReq;
  private final DirRequirement aptExeReq;
  private final DirRequirement aptLibReq;
  private final FileRequirement annotFileReq;
  private final BoolRequirement skipGenoReq;

  public static final String DESC_APT_EXT = "Directory with AffyPowerTools executables (should contain apt-genotype axiom, etc. Available at http://www.affymetrix.com/)";
  public static final String DESC_APT_LIB = "Directory with AffyPowerTools library files (should contain a .cdf file, a .sketch file, etc. Available at http://www.affymetrix.com/)";

  public static AxiomCELProcessingStep create(Project proj, Requirement<Integer> numThreadsReq) {
    DirRequirement aptExtReq = new DirRequirement(DESC_APT_EXT, new File(""));
    DirRequirement aptLibReq = new DirRequirement(DESC_APT_LIB, new File(""));
    String annotFile = proj.SNP_DATA_FILE.getValue().equals("") ? AffyMarkerBlast.EXAMPLE_ANNOT_FILE
                                                                : proj.SNP_DATA_FILE.getValue();
    FileRequirement annotFileReq = new Requirement.FileRequirement(ext.capitalizeFirst(AffyMarkerBlast.DESC_ANNOT_FILE),
                                                                   new File(annotFile));
    BoolRequirement skipGenoReq = new BoolRequirement("Do not import forward genotypes.", false);

    return new AxiomCELProcessingStep(proj, aptExtReq, aptLibReq, annotFileReq, skipGenoReq,
                                      numThreadsReq);
  }

  private AxiomCELProcessingStep(Project proj, DirRequirement aptExtReq, DirRequirement aptLibReq,
                                 FileRequirement annotFileReq, BoolRequirement skipGenoReq,
                                 Requirement<Integer> numThreadsReq) {
    super(NAME, DESC,
          RequirementSetBuilder.and().add(aptExtReq).add(aptLibReq)
                               .add(RequirementSetBuilder.or().add(annotFileReq).add(skipGenoReq))
                               .add(numThreadsReq),
          EnumSet.of(Requirement.Flag.MULTITHREADED, Requirement.Flag.MEMORY,
                     Requirement.Flag.RUNTIME));
    this.proj = proj;
    this.aptExeReq = aptExtReq;
    this.aptLibReq = aptLibReq;
    this.numThreadsReq = numThreadsReq;
    this.annotFileReq = annotFileReq;
    this.skipGenoReq = skipGenoReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    // no-op
  }

  @Override
  public void run(Variables variables) {
    String aptExeDir = ext.verifyDirFormat(variables.get(aptExeReq).getPath());
    String aptLibDir = ext.verifyDirFormat(variables.get(aptLibReq).getPath());
    String annotFile = variables.get(this.skipGenoReq) ? null
                                                       : variables.get(annotFileReq).getPath();
    int numThreads = variables.get(numThreadsReq);

    try {
      APTAxiomPipeline.run(aptLibDir, proj, aptExeDir, annotFile, numThreads);
    } catch (Elision e) {
      throw new RuntimeException(e);
    }
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    return Files.countFiles(proj.MARKER_DATA_DIRECTORY.getValue(false, false),
                            MarkerData.MARKER_DATA_FILE_EXTENSION) > 0;
  }

  @Override
  public String getCommandLine(Variables variables) {
    Resources.axiomTx(proj.getLog()).genome(proj.GENOME_BUILD_VERSION.getValue())
             .getMarkerPositions().get();
    StringBuilder cmd = new StringBuilder();
    cmd.append(Files.getRunString());
    cmd.append(" ").append(APTAxiomPipeline.class.getName());
    cmd.append(" proj=").append(proj.getPropertyFilename());
    cmd.append(" aptExeDir=").append(ext.verifyDirFormat(variables.get(aptExeReq).getPath()));
    cmd.append(" libraryFilePath=").append(ext.verifyDirFormat(variables.get(aptLibReq).getPath()));
    if (variables.get(this.skipGenoReq)) {
      cmd.append(" -skipGenotypes");
    } else {
      cmd.append(" annotFile=").append(variables.get(annotFileReq).getPath());
    }
    cmd.append(" ").append(PSF.Ext.NUM_THREADS_COMMAND).append(variables.get(numThreadsReq));
    return cmd.toString();
  }

}
