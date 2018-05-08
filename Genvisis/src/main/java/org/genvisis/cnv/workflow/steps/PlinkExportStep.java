package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;
import java.util.StringJoiner;
import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.PlinkData;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.Variables;
import org.genvisis.common.Files;
import org.genvisis.common.PSF;

public class PlinkExportStep extends Step {

  public static PlinkExportStep create(Project proj, final Step parseSamplesStep) {
    final Requirement<Step> parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    final Requirement<File> pedigreeRequirement = new Requirement.FileRequirement("A pedigree.dat file must exist.",
                                                                                  new File(proj.PEDIGREE_FILENAME.getValue(false,
                                                                                                                           false)));
    final Requirement<Boolean> createPedigreeRequirement = new Requirement.BoolRequirement("Create a minimal pedigree.dat file [will pull information from SexChecks step results].",
                                                                                           false);

    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(pedigreeRequirement)
                                                                                 .add(createPedigreeRequirement));
    return new PlinkExportStep(proj, pedigreeRequirement, createPedigreeRequirement, reqSet);
  }

  final Requirement<File> pedigreeRequirement;
  final Requirement<Boolean> createPedigreeRequirement;
  final Project proj;

  public static final String NAME = "";
  public static final String DESC = "";

  private PlinkExportStep(Project proj, Requirement<File> pedigreeRequirement,
                          Requirement<Boolean> createPedigreeRequirement, RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.of(Requirement.Flag.MEMORY));
    this.proj = proj;
    this.pedigreeRequirement = pedigreeRequirement;
    this.createPedigreeRequirement = createPedigreeRequirement;
  }

  @Override
  public Variables getDefaultRequirementValues() {
    Variables varMap = super.getDefaultRequirementValues();
    if (!Files.exists(pedigreeRequirement.getDefaultValue())) {
      // if no pedigree, default to creating a minimal one
      varMap.put(createPedigreeRequirement, Boolean.TRUE);
    }
    return varMap;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    if (!variables.get(createPedigreeRequirement)) {
      String projPedFile = proj.PEDIGREE_FILENAME.getValue(false, false);
      File pedFile = variables.get(pedigreeRequirement);
      if (!pedFile.getAbsolutePath().equals(projPedFile)) {
        proj.PEDIGREE_FILENAME.setValue(pedFile.getAbsolutePath());
      }
    }
  }

  @Override
  public void run(Variables variables) {
    if (variables.get(createPedigreeRequirement)) {
      proj.getLog().report("Creating Pedigree File");
      Pedigree.build(proj, null, null, false);
    }
    if (!Files.exists(proj.PEDIGREE_FILENAME.getValue())) {
      throw new RuntimeException("Creation of Pedigree file in [Create/Run PLINK Files] step failed.");
    }

    proj.getLog().report("Running PLINK");

    boolean create = PlinkData.saveGenvisisToPlinkBedSet(proj,
                                                         GenvisisWorkflow.PLINK_SUBDIR
                                                               + GenvisisWorkflow.PLINKROOT,
                                                         null, null,
                                                         PlinkData.ExportIDScheme.DNA_DNA);
    if (!create) {
      throw new RuntimeException("Creation of initial PLINK files failed.");
    }
    proj.PLINK_DIR_FILEROOTS.addValue(proj.PROJECT_DIRECTORY.getValue()
                                      + GenvisisWorkflow.PLINK_SUBDIR + GenvisisWorkflow.PLINKROOT);
  }

  @Override
  public String getCommandLine(Variables variables) {
    String kvCmd = "";

    if (!variables.get(createPedigreeRequirement)) {
      String projPedFile = proj.PEDIGREE_FILENAME.getValue(false, false);
      File pedFile = variables.get(pedigreeRequirement);
      if (!pedFile.getAbsolutePath().equals(projPedFile)) {
        kvCmd += " PEDIGREE_FILENAME=" + pedFile.getAbsolutePath();
      }
    }

    String projPropFile = proj.getPropertyFilename();
    StringBuilder cmd = new StringBuilder();
    if (kvCmd.length() > 0) {
      cmd.append(Files.getRunString()).append(GenvisisWorkflow.PROJ_PROP_UPDATE_STR)
         .append(projPropFile).append(kvCmd).append("\n");
    }
    if (variables.get(createPedigreeRequirement)) {
      cmd.append(Files.getRunString()).append(" cnv.filesys.Pedigree proj=").append(projPropFile)
         .append("\n");
    }
    cmd.append(new StringJoiner(" ").add(Files.getRunString()).add(PlinkData.class.getName())
                                    .add("-genvisisToBed")
                                    .add("plinkdata=" + GenvisisWorkflow.PLINK_SUBDIR
                                         + GenvisisWorkflow.PLINKROOT)
                                    .add("proj=" + proj.getPropertyFilename())
                                    .add(PlinkData.ARG_EXPORT_ID_SCHEME
                                         + PlinkData.ExportIDScheme.DNA_DNA));
    return cmd.toString();
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    boolean plinkFilesExist = Files.checkAllFiles(GenvisisWorkflow.getPlinkDir(proj),
                                                  PSF.Plink.getPlinkBedBimFamSet(GenvisisWorkflow.PLINKROOT),
                                                  false, proj.getLog());
    boolean pedGenerated = variables.get(createPedigreeRequirement);
    boolean pedCheck = pedGenerated ? Files.exists(proj.PEDIGREE_FILENAME.getValue()) : true;
    return plinkFilesExist && pedCheck;
  }

}
