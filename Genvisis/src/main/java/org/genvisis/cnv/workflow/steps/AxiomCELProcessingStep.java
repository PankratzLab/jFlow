package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;

import org.genvisis.cnv.affy.APTAxiomPipeline;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.Requirement.BoolRequirement;
import org.genvisis.cnv.workflow.Requirement.DirRequirement;
import org.genvisis.cnv.workflow.Requirement.FileRequirement;
import org.genvisis.cnv.workflow.Requirement.StepRequirement;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.cnv.workflow.Variables;
import org.pankratzlab.common.Elision;
import org.pankratzlab.common.Files;
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
  private final FileRequirement xmlSchemaReq;

  public static final String DESC_APT_EXT = "Directory with AffyPowerTools executables (should contain apt-genotype axiom, etc. Available at http://www.affymetrix.com/)";
  public static final String DESC_APT_LIB = "Directory with AffyPowerTools library files (should contain a .cdf file, a .sketch file, etc. Available at http://www.affymetrix.com/)";

  public static AxiomCELProcessingStep create(Project proj, Step axiomManifestStep,
                                              Requirement<Integer> numThreadsReq) {
    StepRequirement manifestReq = new StepRequirement(axiomManifestStep);

    DirRequirement aptExtReq = new DirRequirement("aptExe", DESC_APT_EXT, new File(""));
    DirRequirement aptLibReq = new DirRequirement("aptLib", DESC_APT_LIB, new File(""));

    FileRequirement xmlSchemaReq = new FileRequirement("xmlSchema",
                                                       "An apt-genotype-axiom XML schema file (likely in the same directory as the library files, e.g. Axiom_tx_v1.r5.apt-genotype-axiom.AxiomCN_GT1.apt2.custom.xml)",
                                                       new File(""));

    String annotFile = proj.SNP_DATA_FILE.getValue()
                                         .equals("") ? AxiomManifestParsingStep.AXIOM_EXAMPLE_MANIFEST
                                                     : proj.SNP_DATA_FILE.getValue();
    FileRequirement annotFileReq = new Requirement.FileRequirement("annotFile",
                                                                   AxiomManifestParsingStep.AXIOM_MANIFEST_DESC,
                                                                   new File(annotFile));
    BoolRequirement skipGenoReq = new BoolRequirement("skipGenotypes",
                                                      "Do not import forward genotypes.", false);

    return new AxiomCELProcessingStep(proj, manifestReq, aptExtReq, aptLibReq, xmlSchemaReq,
                                      annotFileReq, skipGenoReq, numThreadsReq);
  }

  private AxiomCELProcessingStep(Project proj, StepRequirement manifestReq,
                                 DirRequirement aptExtReq, DirRequirement aptLibReq,
                                 FileRequirement xmlSchemaReq, FileRequirement annotFileReq,
                                 BoolRequirement skipGenoReq, Requirement<Integer> numThreadsReq) {
    super(NAME, DESC,
          RequirementSetBuilder.and().add(manifestReq).add(aptExtReq).add(aptLibReq)
                               .add(xmlSchemaReq)
                               .add(RequirementSetBuilder.or().add(annotFileReq).add(skipGenoReq))
                               .add(numThreadsReq),
          EnumSet.of(Requirement.Flag.MULTITHREADED, Requirement.Flag.MEMORY,
                     Requirement.Flag.RUNTIME));
    this.proj = proj;
    this.xmlSchemaReq = xmlSchemaReq;
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
    String xmlFile = variables.get(xmlSchemaReq).getPath();
    String aptExeDir = ext.verifyDirFormat(variables.get(aptExeReq).getPath());
    String aptLibDir = ext.verifyDirFormat(variables.get(aptLibReq).getPath());
    String annotFile = variables.get(this.skipGenoReq) ? null
                                                       : variables.get(annotFileReq).getPath();
    int numThreads = variables.get(numThreadsReq);

    try {
      APTAxiomPipeline.run(xmlFile, aptLibDir, proj, aptExeDir, annotFile, numThreads);
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
    return getStepCommandLine(proj, variables);
  }

  public static void main(String[] args) {
    Project proj = Step.parseProject(args);
    StepBuilder sb = new StepBuilder(proj);
    AxiomManifestParsingStep manifestStep = sb.generateAxiomManifestParsingStep();
    Step step = sb.generateAxiomCELProcessingStep(manifestStep);
    Variables variables = step.parseArguments(args);
    Step.run(proj, step, variables);
  }

}
