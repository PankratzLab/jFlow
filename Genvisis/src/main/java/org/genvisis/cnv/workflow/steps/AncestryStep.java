package org.genvisis.cnv.workflow.steps;

import java.util.EnumSet;
import java.util.Map;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.Requirement.FileRequirement;
import org.genvisis.cnv.workflow.Requirement.ResourceRequirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.gwas.Ancestry;

public class AncestryStep extends Step {

  public static final String NAME = "Run Ancestry Checks";
  public static final String DESC = "";

  public static AncestryStep create(Project proj, Step gwasQCStep, double priority) {
    final Requirement gwasQCStepReq = new Requirement.StepRequirement(gwasQCStep);
    final Requirement.ResourceRequirement hapMapFoundersReq = new Requirement.ResourceRequirement("PLINK root of HapMap founders",
                                                                                                  Resources.hapMap(proj.getLog())
                                                                                                           .getUnambiguousHapMapFounders());

    final RequirementSet reqSet = RequirementSetBuilder.and().add(gwasQCStepReq)
                                                       .add(putativeWhitesReq)
                                                       .add(hapMapFoundersReq);
    final Requirement.ResourceRequirement hapMapAncestryReq = new Requirement.ResourceRequirement("HapMap Samples Ancestry File",
                                                                                                  Resources.hapMap(proj.getLog())
                                                                                                           .getHapMapAncestries());
    return new AncestryStep(hapMapFoundersReq, hapMapAncestryReq, reqSet, priority);
  }

  private final static Requirement putativeWhitesReq = new FileRequirement("File with FID/IID pairs of putative white samples",
                                                                           "");

  final ResourceRequirement hapMapFoundersReq;
  final ResourceRequirement hapMapAncestryReq;

  private AncestryStep(ResourceRequirement hapFound, ResourceRequirement hapAnc,
                       RequirementSet reqSet, double priority) {
    super(NAME, DESC, reqSet, EnumSet.noneOf(Requirement.Flag.class), priority);
    this.hapMapFoundersReq = hapFound;
    this.hapMapAncestryReq = hapAnc;
  }

  @Override
  public void setNecessaryPreRunProperties(Project proj,
                                           Map<Step, Map<Requirement, String>> variables) {
    // not needed for step
  }

  @Override
  public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String putativeWhites = variables.get(this).get(putativeWhitesReq);
    String hapMapPlinkRoot = hapMapFoundersReq.getResource().getAbsolute();
    hapMapAncestryReq.getResource().get();
    String ancestryDir = GenvisisWorkflow.getAncestryDir(proj);
    Ancestry.runPipeline(ancestryDir, putativeWhites, hapMapPlinkRoot, proj,
                         new Logger(ancestryDir + "ancestry.log"));
  }

  @Override
  public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String putativeWhites = variables.get(this).get(putativeWhitesReq);
    String hapMapPlinkRoot = hapMapFoundersReq.getResource().getAbsolute();
    hapMapAncestryReq.getResource().get();
    String ancestryDir = GenvisisWorkflow.getAncestryDir(proj);
    String command = Files.getRunString() + " gwas.Ancestry -runPipeline dir=" + ancestryDir;
    command += " putativeWhites=" + putativeWhites;
    command += " proj=" + proj.getPropertyFilename();
    command += " hapMapPlinkRoot=" + hapMapPlinkRoot;
    command += " log=" + ancestryDir + "ancestry.log";
    return command;
  }

  @Override
  public boolean checkIfOutputExists(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String ancestryDir = GenvisisWorkflow.getAncestryDir(proj);
    return Files.exists(ancestryDir + Ancestry.RACE_FREQS_FILENAME)
           && Files.exists(ancestryDir + Ancestry.RACE_IMPUTATIONAS_FILENAME);
  }

}
