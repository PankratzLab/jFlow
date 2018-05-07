package org.genvisis.cnv.workflow.steps;

import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import org.genvisis.CLI;
import org.genvisis.cnv.analysis.pca.PCImputeRace;
import org.genvisis.cnv.analysis.pca.PCImputeRace.RACE;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.common.Files;
import org.genvisis.common.PSF;
import org.genvisis.gwas.Ancestry;
import org.genvisis.gwas.FurtherAnalysisQc;
import org.genvisis.gwas.MarkerQC.QC_METRIC;
import org.genvisis.gwas.Qc;
import org.genvisis.gwas.RelationAncestryQc;
import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class FurtherAnalysisQCStep extends Step {

  public static final String NAME = "Run Further Analysis QC";
  public static final String DESC = "";

  final Map<QC_METRIC, Requirement> metricRequirements;

  public static FurtherAnalysisQCStep create(Project proj, Step plinkExportStep, Step gwasQCStep,
                                             Step ancestryStep) {
    final Requirement plinkExportStepReq = new Requirement.StepRequirement(plinkExportStep);
    final Requirement gwasQCStepReq = new Requirement.StepRequirement(gwasQCStep);
    final Requirement ancestryStepReq = new Requirement.StepRequirement(ancestryStep);
    final Requirement unrelatedsFileReq = new Requirement.FileRequirement("File with list of unrelated FID/IID pairs to use for marker QC",
                                                                          "");
    final Requirement europeansFilesReq = new Requirement.FileRequirement("File with list of European samples to use for Hardy-Weinberg equilibrium tests",
                                                                          "");
    final RequirementSet reqSet = RequirementSetBuilder.and().add(plinkExportStepReq)
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(gwasQCStepReq)
                                                                                 .add(unrelatedsFileReq))
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(ancestryStepReq)
                                                                                 .add(europeansFilesReq));

    Map<QC_METRIC, Requirement> metricRequirements = Maps.newEnumMap(QC_METRIC.class);
    for (QC_METRIC metric : QC_METRIC.values()) {
      Map<QC_METRIC, String> defaultThresholds = FurtherAnalysisQc.getDefaultMarkerQCThresholds(proj.getArrayType());
      String defaultVal = defaultThresholds.get(metric);
      final Requirement metricReq = new Requirement.ThresholdRequirement(metric.getUserDescription(),
                                                                         defaultVal);
      reqSet.add(metricReq);
      metricRequirements.put(metric, metricReq);
    }

    return new FurtherAnalysisQCStep(metricRequirements, unrelatedsFileReq, europeansFilesReq,
                                     reqSet);
  }

  final Requirement unrelatedsFileReq;
  final Requirement europeansFilesReq;

  private FurtherAnalysisQCStep(Map<QC_METRIC, Requirement> metricReqs, Requirement unrelReq,
                                Requirement euroReq, RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.noneOf(Requirement.Flag.class));
    this.unrelatedsFileReq = unrelReq;
    this.europeansFilesReq = euroReq;
    this.metricRequirements = metricReqs;
  }

  @Override
  public void setNecessaryPreRunProperties(Project proj,
                                           Map<Step, Map<Requirement, String>> variables) {
    // not needed for step
  }

  @Override
  public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
    Map<Requirement, String> stepVars = variables.get(this);

    String unrelatedsFile = resolveUnrelatedsFile(proj, stepVars);

    String europeansFile = resolveEuropeansFile(proj, stepVars);

    Map<QC_METRIC, String> markerQCThresholds = Maps.newEnumMap(QC_METRIC.class);
    for (QC_METRIC metric : QC_METRIC.values()) {
      Requirement req = metricRequirements.get(metric);
      markerQCThresholds.put(metric, stepVars.get(req));
    }
    new FurtherAnalysisQc(GenvisisWorkflow.getPlinkDir(proj), GenvisisWorkflow.PLINKROOT,
                          markerQCThresholds, unrelatedsFile, europeansFile, proj.getLog())
                                                                                           .runFurtherAnalysisQC();
  }

  @Override
  public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
    Map<Requirement, String> stepVars = variables.get(this);

    String unrelatedsFile = resolveUnrelatedsFile(proj, stepVars);

    String europeansFile = resolveEuropeansFile(proj, stepVars);

    List<String> commandChunks = Lists.newArrayList();
    commandChunks.add(Files.getRunString());
    commandChunks.add(FurtherAnalysisQc.class.getName());
    commandChunks.add(CLI.formCmdLineArg(FurtherAnalysisQc.ARG_UNRELATEDS, unrelatedsFile));
    commandChunks.add(CLI.formCmdLineArg(FurtherAnalysisQc.ARG_EUROPEANS, europeansFile));
    commandChunks.add(CLI.formCmdLineArg(CLI.ARG_INDIR, GenvisisWorkflow.getPlinkDir(proj)));
    commandChunks.add(CLI.formCmdLineArg(CLI.ARG_PLINKROOT, GenvisisWorkflow.PLINKROOT));
    for (QC_METRIC metric : QC_METRIC.values()) {
      Requirement req = metricRequirements.get(metric);
      commandChunks.add(CLI.formCmdLineArg(metric.getKey(), stepVars.get(req)));
    }
    return Joiner.on(' ').join(commandChunks);
  }

  @Override
  public boolean checkIfOutputExists(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String dir = GenvisisWorkflow.getPlinkDir(proj) + Qc.QC_SUBDIR
                 + FurtherAnalysisQc.FURTHER_ANALYSIS_DIR;
    String qcdPlinkroot = GenvisisWorkflow.PLINKROOT
                          + FurtherAnalysisQc.FURTHER_ANALYSIS_QC_PLINK_SUFFIX;
    return PSF.Plink.bedBimFamExist(dir + qcdPlinkroot)
           && Files.exists(dir + FurtherAnalysisQc.SAMPLE_QC_DROPS, false)
           && Files.exists(dir + FurtherAnalysisQc.MARKER_QC_DROPS, false);
  }

  private String resolveUnrelatedsFile(Project proj, Map<Requirement, String> stepVars) {
    String unrelatedsFile = stepVars.get(unrelatedsFileReq);
    if (!Files.exists(unrelatedsFile)) {
      unrelatedsFile = GenvisisWorkflow.getAncestryDir(proj)
                       + RelationAncestryQc.UNRELATEDS_FILENAME;
    }
    return unrelatedsFile;
  }

  private String resolveEuropeansFile(Project proj, Map<Requirement, String> stepVars) {
    String europeansFile = stepVars.get(europeansFilesReq);
    if (europeansFile == null || "".equals(europeansFile)) {
      String raceImputationFilename = GenvisisWorkflow.getAncestryDir(proj)
                                      + Ancestry.RACE_IMPUTATIONAS_FILENAME;
      europeansFile = PCImputeRace.formRaceListFilename(RACE.WHITE, raceImputationFilename);
    }
    return europeansFile;
  }

}
