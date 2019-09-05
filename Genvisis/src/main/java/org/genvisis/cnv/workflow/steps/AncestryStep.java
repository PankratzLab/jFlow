package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.io.PrintWriter;
import java.util.EnumSet;

import org.genvisis.cnv.Resources;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.gwas.Ancestry;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.ProgramRequirement;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.Requirement.BoolRequirement;
import org.genvisis.cnv.workflow.Requirement.FileRequirement;
import org.genvisis.cnv.workflow.Requirement.OptionalFileRequirement;
import org.genvisis.cnv.workflow.Requirement.ResourceRequirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.cnv.workflow.Variables;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.ext;

public class AncestryStep extends Step {

  public static final String NAME = "Run Ancestry Checks";
  public static final String DESC = "";
  public static final String PUTATIVE_ALL_FILENAME = "putativeWhites_AllSamples.txt";

  public static AncestryStep create(Project proj, Step gwasQCStep) {
    final Requirement<Step> gwasQCStepReq = new Requirement.StepRequirement(gwasQCStep);
    final Requirement.ResourceRequirement hapMapFoundersReq = new Requirement.ResourceRequirement("PLINK root of HapMap founders",
                                                                                                  Resources.hapMap(proj.getLog())
                                                                                                           .getUnambiguousHapMapFounders());

    String defaultSnpFile = null;
    if (proj.ARRAY_TYPE.getValue() == ARRAY.AFFY_GW6
        || proj.ARRAY_TYPE.getValue() == ARRAY.AFFY_GW6_CN) {
      defaultSnpFile = Resources.affy(proj.getLog()).getRSIDLookup().get();
    } else if (proj.ARRAY_TYPE.getValue() == ARRAY.AFFY_AXIOM) {
      String[] lookups = Files.list(proj.PROJECT_DIRECTORY.getValue(),
                                    AxiomManifestParsingStep.RSLOOKUP_SUFFIX);
      defaultSnpFile = lookups.length > 0 ? lookups[0] : "";
    }
    if (defaultSnpFile == null) {
      defaultSnpFile = "";
    }

    final Requirement.OptionalFileRequirement snpIDLookupFileReq = new OptionalFileRequirement("snpNameLookupFile",
                                                                                               "A SNP name replacement file with two columns, the first being the original SNP name and the second containing the replacement name.",
                                                                                               new File(defaultSnpFile));
    final ProgramRequirement plinkExeReq = new ProgramRequirement(CLI.ARG_PLINK_EXE,
                                                                  ext.capitalizeFirst(CLI.DESC_PLINK_EXE),
                                                                  CLI.DEF_PLINK_EXE);
    FileRequirement putativeWhitesReq = new FileRequirement("putativeWhitesFile",
                                                            "File with FID/IID pairs of putative white samples",
                                                            new File(""));
    BoolRequirement useAllReq = new BoolRequirement("useAll",
                                                    "Use all samples (as putative whites)", false);

    final RequirementSet reqSet = RequirementSetBuilder.and().add(gwasQCStepReq)
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(putativeWhitesReq)
                                                                                 .add(useAllReq))
                                                       .add(hapMapFoundersReq).add(plinkExeReq)
                                                       .add(snpIDLookupFileReq);
    final Requirement.ResourceRequirement hapMapAncestryReq = new Requirement.ResourceRequirement("HapMap Samples Ancestry File",
                                                                                                  Resources.hapMap(proj.getLog())
                                                                                                           .getHapMapAncestries());
    return new AncestryStep(proj, snpIDLookupFileReq, hapMapFoundersReq, hapMapAncestryReq,
                            putativeWhitesReq, useAllReq, reqSet, plinkExeReq);
  }

  final Project proj;
  final OptionalFileRequirement snpIDLookupReq;
  final ResourceRequirement hapMapFoundersReq;
  final ResourceRequirement hapMapAncestryReq;
  final ProgramRequirement plinkExeReq;
  final FileRequirement putativeWhitesReq;
  final BoolRequirement useAllSamplesReq;

  private AncestryStep(Project proj, OptionalFileRequirement snpIDLookupReq,
                       ResourceRequirement hapFound, ResourceRequirement hapAnc,
                       FileRequirement putativeWhitesReq, BoolRequirement useAllReq,
                       RequirementSet reqSet, ProgramRequirement plinkExeReq) {
    super(NAME, DESC, reqSet, EnumSet.noneOf(Requirement.Flag.class));
    this.proj = proj;
    this.snpIDLookupReq = snpIDLookupReq;
    this.hapMapFoundersReq = hapFound;
    this.hapMapAncestryReq = hapAnc;
    this.plinkExeReq = plinkExeReq;
    this.putativeWhitesReq = putativeWhitesReq;
    this.useAllSamplesReq = useAllReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    // not needed for step
  }

  @Override
  public void run(Variables variables) {
    String putativeWhites = null;
    if (variables.get(useAllSamplesReq).booleanValue()) {
      putativeWhites = proj.PROJECT_DIRECTORY.getValue() + PUTATIVE_ALL_FILENAME;
      if (!Files.exists(putativeWhites)
          || Files.countLines(putativeWhites, 0) != proj.getSamples().length) {
        PrintWriter writer = Files.getAppropriateWriter(putativeWhites);
        for (String samp : proj.getSamples()) {
          writer.println(samp + "\t" + samp);
        }
        writer.close();
      }
    } else if (variables.get(putativeWhitesReq) != null) {
      putativeWhites = variables.get(putativeWhitesReq).getAbsolutePath();
    }

    String hapMapPlinkRoot = hapMapFoundersReq.getResource().getAbsolute();
    hapMapAncestryReq.getResource().get();
    String ancestryDir = GenvisisWorkflow.getAncestryDir(proj);
    File f = variables.get(snpIDLookupReq);
    String snpIDFile = f == null || f.getPath().equals("") ? null : f.getAbsolutePath();
    new Ancestry(ancestryDir, proj, variables.get(plinkExeReq)).runPipeline(putativeWhites,
                                                                            hapMapPlinkRoot,
                                                                            snpIDFile);
  }

  @Override
  public String getCommandLine(Variables variables) {
    return getStepCommandLine(proj, variables);
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    String ancestryDir = GenvisisWorkflow.getAncestryDir(proj);
    return Files.exists(ancestryDir + Ancestry.RACE_FREQS_FILENAME)
           && Files.exists(ancestryDir + Ancestry.RACE_IMPUTATIONS_FILENAME);
  }

  public static void main(String[] args) {
    Project proj = Step.parseProject(args);
    StepBuilder sb = new StepBuilder(proj);
    Step samplesStep = sb.generateSamplesParsingStep();
    PlinkExportStep plinkExportStep = sb.generatePlinkExportStep(samplesStep);
    GwasQCStep gwasQCStep = sb.generateGwasQCStep(plinkExportStep);
    AncestryStep ancestryStep = sb.generateAncestryStep(gwasQCStep);
    Variables variables = ancestryStep.parseArguments(args);
    Step.run(proj, ancestryStep, variables);
  }

}
