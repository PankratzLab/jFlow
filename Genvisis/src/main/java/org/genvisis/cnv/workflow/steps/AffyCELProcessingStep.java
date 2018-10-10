package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;
import org.genvisis.cnv.Resources;
import org.genvisis.cnv.affy.AffyPipeline;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.Requirement.DirRequirement;
import org.genvisis.cnv.workflow.Requirement.FileRequirement;
import org.genvisis.cnv.workflow.Requirement.OptionalBoolRequirement;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.seq.GenomeBuild;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.Variables;
import org.pankratzlab.common.Elision;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ext;

public class AffyCELProcessingStep extends Step {

  public static final String NAME = "Process Affymetrix CEL files";
  public static final String DESC = "";

  private final Project proj;

  private final OptionalBoolRequirement fullReq;
  private final Requirement<Integer> numThreadsReq;
  private final DirRequirement aptExeReq;
  private final DirRequirement aptLibReq;
  private final FileRequirement sketchReq;

  public static final String DESC_FULL = "Use the full affymetrix cdf, which contains more mitochondrial probesets.";
  public static final String DESC_APT_EXT = "Directory with Affy Power Tools executables (should contain apt-probeset-genotype, etc. Available at http://www.affymetrix.com/)";
  public static final String DESC_APT_LIB = "Directory with Affy Power Tools library files (should contain GenomeWideSNP_6.cdf, etc. Available at http://www.affymetrix.com/)";
  public static final String DESC_SKETCH = "A target sketch file (such as hapmap.quant-norm.normalization-target.txt)";

  public static AffyCELProcessingStep create(Project proj, Requirement<Integer> numThreadsReq) {

    OptionalBoolRequirement fullReq = new OptionalBoolRequirement(DESC_FULL, false);
    DirRequirement aptExtReq = new DirRequirement(DESC_APT_EXT, new File(""));
    DirRequirement aptLibReq = new DirRequirement(DESC_APT_LIB, new File(""));
    FileRequirement sketchReq = new FileRequirement(DESC_SKETCH, new File(""));

    return new AffyCELProcessingStep(proj, aptExtReq, aptLibReq, sketchReq, fullReq, numThreadsReq);
  }

  private AffyCELProcessingStep(Project proj, DirRequirement aptExtReq, DirRequirement aptLibReq,
                                FileRequirement sketchReq, OptionalBoolRequirement fullReq,
                                Requirement<Integer> numThreadsReq) {
    super(NAME, DESC,
          RequirementSetBuilder.and().add(aptExtReq).add(aptLibReq).add(sketchReq).add(fullReq)
                               .add(numThreadsReq),
          EnumSet.of(Requirement.Flag.MULTITHREADED, Requirement.Flag.MEMORY,
                     Requirement.Flag.RUNTIME));
    this.proj = proj;
    this.aptExeReq = aptExtReq;
    this.aptLibReq = aptLibReq;
    this.sketchReq = sketchReq;
    this.fullReq = fullReq;
    this.numThreadsReq = numThreadsReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    // no-op
  }

  @Override
  public void run(Variables variables) {
    String aptExeDir = ext.verifyDirFormat(variables.get(aptExeReq).getPath());
    String aptLibDir = ext.verifyDirFormat(variables.get(aptLibReq).getPath());
    String quantNormTarget = variables.get(sketchReq).getPath();
    boolean full = variables.get(fullReq);
    int numThreads = variables.get(numThreadsReq);

    try {
      AffyPipeline.run(proj, aptExeDir, aptLibDir, quantNormTarget, full, numThreads);
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
    Resources.affy(proj.getLog()).genome(GenomeBuild.HG19).getMarkerPositions().get(); // download if necessary

    StringBuilder cmd = new StringBuilder();
    cmd.append(Files.getRunString());
    cmd.append(" ").append(AffyPipeline.class.getName());
    cmd.append(" proj=").append(proj.getPropertyFilename());
    cmd.append(" aptExeDir=").append(ext.verifyDirFormat(variables.get(aptExeReq).getPath()));
    cmd.append(" aptLibDir=").append(ext.verifyDirFormat(variables.get(aptLibReq).getPath()));
    cmd.append(" sketch=").append(variables.get(sketchReq).getPath());
    cmd.append(" ").append(PSF.Ext.NUM_THREADS_COMMAND).append(variables.get(numThreadsReq));
    if (variables.get(fullReq)) {
      cmd.append(" -full");
    }
    return cmd.toString();
  }

}
