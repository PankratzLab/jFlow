package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;

import org.genvisis.cnv.Resources;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.gwas.Ancestry;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.Requirement.FileRequirement;
import org.genvisis.cnv.workflow.Requirement.OptionalFileRequirement;
import org.genvisis.cnv.workflow.Requirement.ResourceRequirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.Variables;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;

public class AncestryStep extends Step {

  public static final String NAME = "Run Ancestry Checks";
  public static final String DESC = "";

  public static AncestryStep create(Project proj, Step gwasQCStep) {
    final Requirement<Step> gwasQCStepReq = new Requirement.StepRequirement(gwasQCStep);
    final Requirement.ResourceRequirement hapMapFoundersReq = new Requirement.ResourceRequirement("PLINK root of HapMap founders",
                                                                                                  Resources.hapMap(proj.getLog())
                                                                                                           .getUnambiguousHapMapFounders());

    String defaultSnpFile = null;
    if (proj.ARRAY_TYPE.getValue() == ARRAY.AFFY_GW6
        || proj.ARRAY_TYPE.getValue() == ARRAY.AFFY_GW6_CN) {
      defaultSnpFile = Resources.affy(proj.getLog()).getRSIDLookup().get();
    }
    if (defaultSnpFile == null) {
      defaultSnpFile = "";
    }

    final Requirement.OptionalFileRequirement snpIDLookupFileReq = new OptionalFileRequirement("snpNameLookupFile",
                                                                                               "A SNP name replacement file with two columns, the first being the original SNP name and the second containing the replacement name.",
                                                                                               new File(defaultSnpFile));
    final RequirementSet reqSet = RequirementSetBuilder.and().add(gwasQCStepReq)
                                                       .add(putativeWhitesReq)
                                                       .add(hapMapFoundersReq)
                                                       .add(snpIDLookupFileReq);
    final Requirement.ResourceRequirement hapMapAncestryReq = new Requirement.ResourceRequirement("HapMap Samples Ancestry File",
                                                                                                  Resources.hapMap(proj.getLog())
                                                                                                           .getHapMapAncestries());
    return new AncestryStep(proj, snpIDLookupFileReq, hapMapFoundersReq, hapMapAncestryReq, reqSet);
  }

  private final static Requirement<File> putativeWhitesReq = new FileRequirement("putativeWhitesFile",
                                                                                 "File with FID/IID pairs of putative white samples",
                                                                                 new File(""));

  final Project proj;
  final OptionalFileRequirement snpIDLookupReq;
  final ResourceRequirement hapMapFoundersReq;
  final ResourceRequirement hapMapAncestryReq;

  private AncestryStep(Project proj, OptionalFileRequirement snpIDLookupReq,
                       ResourceRequirement hapFound, ResourceRequirement hapAnc,
                       RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.noneOf(Requirement.Flag.class));
    this.proj = proj;
    this.snpIDLookupReq = snpIDLookupReq;
    this.hapMapFoundersReq = hapFound;
    this.hapMapAncestryReq = hapAnc;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    // not needed for step
  }

  @Override
  public void run(Variables variables) {
    String putativeWhites = variables.get(putativeWhitesReq) == null ? null
                                                                     : variables.get(putativeWhitesReq)
                                                                                .getAbsolutePath();
    String hapMapPlinkRoot = hapMapFoundersReq.getResource().getAbsolute();
    hapMapAncestryReq.getResource().get();
    String ancestryDir = GenvisisWorkflow.getAncestryDir(proj);
    File f = variables.get(snpIDLookupReq);
    String snpIDFile = f == null || f.getPath().equals("") ? null : f.getAbsolutePath();
    Ancestry.runPipeline(ancestryDir, putativeWhites, hapMapPlinkRoot, snpIDFile, proj,
                         new Logger(ancestryDir + "ancestry.log"));
  }

  @Override
  public String getCommandLine(Variables variables) {
    String putativeWhites = variables.get(putativeWhitesReq) == null ? null
                                                                     : variables.get(putativeWhitesReq)
                                                                                .getAbsolutePath();
    String hapMapPlinkRoot = hapMapFoundersReq.getResource().getAbsolute();
    hapMapAncestryReq.getResource().get();
    File f = variables.get(snpIDLookupReq);
    String snpIDFile = f == null || f.getPath().equals("") ? null : f.getAbsolutePath();
    String ancestryDir = GenvisisWorkflow.getAncestryDir(proj);
    String command = Files.getRunString() + " " + Ancestry.class.getName() + " -runPipeline dir="
                     + ancestryDir;
    command += " putativeWhites=" + putativeWhites;
    command += " proj=" + proj.getPropertyFilename();
    command += " hapMapPlinkRoot=" + hapMapPlinkRoot;
    if (snpIDFile != null) {
      if (Files.exists(snpIDFile)) {
        command += " snpLookup=" + snpIDFile;
      } else {
        proj.getLog()
            .reportTimeWarning("Specified SNP name lookup file couldn't be found: " + snpIDFile);
      }
    }
    command += " log=" + ancestryDir + "ancestry.log";
    return command;
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    String ancestryDir = GenvisisWorkflow.getAncestryDir(proj);
    return Files.exists(ancestryDir + Ancestry.RACE_FREQS_FILENAME)
           && Files.exists(ancestryDir + Ancestry.RACE_IMPUTATIONS_FILENAME);
  }

}
