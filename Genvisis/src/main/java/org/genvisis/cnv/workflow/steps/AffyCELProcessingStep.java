package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;
import org.genvisis.affy.AffyPipeline;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.SOURCE_FILE_DELIMITERS;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.Requirement.FileRequirement;
import org.genvisis.cnv.workflow.Requirement.OptionalBoolRequirement;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.Variables;
import org.genvisis.common.Files;
import org.genvisis.common.PSF;

public class AffyCELProcessingStep extends Step {

  public static final String NAME = "Process Affymetrix CEL files";
  public static final String DESC = "";

  Project proj;

  OptionalBoolRequirement fullReq;
  Requirement<Integer> numThreadsReq;
  FileRequirement aptExeReq;
  FileRequirement aptLibReq;
  FileRequirement sketchReq;

  public static final String DESC_MKR_BUFF = "Number of markers to buffer when splitting files.";
  public static final String DESC_MAX_WRIT = "Maximum number of writers to open, if this is less than the sample size parsing will slow drastically.";
  public static final String DESC_FULL = "Use the full affymetrix cdf, which contains more mitochondrial probesets.";
  public static final String DESC_APT_EXT = "Directory with Affy Power Tools executables (should contain apt-probeset-genotype, etc. Available at http://www.affymetrix.com/)";
  public static final String DESC_APT_LIB = "Directory with Affy Power Tools library files (should contain GenomeWideSNP_6.cdf, etc. Available at http://www.affymetrix.com/)";
  public static final String DESC_SKETCH = "A target sketch file (such as hapmap.quant-norm.normalization-target.txt)";

  public static AffyCELProcessingStep create(Project proj, Requirement<Integer> numThreadsReq) {

    OptionalBoolRequirement fullReq = new OptionalBoolRequirement(DESC_FULL, false);
    FileRequirement aptExtReq = new FileRequirement(DESC_APT_EXT, new File(""));
    FileRequirement aptLibReq = new FileRequirement(DESC_APT_LIB, new File(""));
    FileRequirement sketchReq = new FileRequirement(DESC_SKETCH, new File(""));

    return new AffyCELProcessingStep(proj, aptExtReq, aptLibReq, sketchReq, fullReq, numThreadsReq);
  }

  private AffyCELProcessingStep(Project proj, FileRequirement aptExtReq, FileRequirement aptLibReq,
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
    String aptExeDir = variables.get(aptExeReq).getPath();
    String aptLibDir = variables.get(aptLibReq).getPath();
    String quantNormTarget = variables.get(sketchReq).getPath();
    int markerBuffer = AffyPipeline.DEFAULT_MARKER_BUFFER;
    int maxWritersOpen = AffyPipeline.DEFAULT_MAX_WRITERS;
    boolean full = variables.get(fullReq);
    int numThreads = variables.get(numThreadsReq);

    AffyPipeline.run(proj, aptExeDir, aptLibDir, quantNormTarget, markerBuffer, maxWritersOpen,
                     full, numThreads);
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    // four measures that indicate success, from AffyPipeline
    return proj.SOURCE_DIRECTORY.getValue().equals("00src_CEL/")
           && proj.SOURCE_FILENAME_EXTENSION.getValue().equals(".txt.gz")
           && proj.SOURCE_FILE_DELIMITER.getValue() == SOURCE_FILE_DELIMITERS.TAB
           && proj.ID_HEADER.getValue().equals("[FILENAME_ROOT]");
  }

  @Override
  public String getCommandLine(Variables variables) {
    StringBuilder cmd = new StringBuilder();
    cmd.append(Files.getRunString());
    cmd.append(" ").append(AffyPipeline.class.getName());
    cmd.append(" proj=").append(proj.getPropertyFilename());
    cmd.append(" aptExeDir=").append(variables.get(aptExeReq).getPath());
    cmd.append(" aptLibDir=").append(variables.get(aptLibReq).getPath());
    cmd.append(" ").append(PSF.Ext.NUM_THREADS_COMMAND).append(variables.get(numThreadsReq));
    if (variables.get(fullReq)) {
      cmd.append(" -full");
    }
    return cmd.toString();
  }

}
