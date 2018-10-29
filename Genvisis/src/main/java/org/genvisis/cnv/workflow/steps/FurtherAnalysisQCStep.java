package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import org.genvisis.cnv.analysis.pca.PCImputeRace;
import org.genvisis.cnv.analysis.pca.PCImputeRace.RACE;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.gwas.Ancestry;
import org.genvisis.cnv.gwas.utils.FurtherAnalysisQc;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.Variables;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.PSF;
import org.pankratzlab.utils.gwas.Qc;
import org.pankratzlab.utils.gwas.QcMetric;
import org.pankratzlab.utils.gwas.RelationAncestryQc;
import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class FurtherAnalysisQCStep extends Step {

  public static final String NAME = "Run Further Analysis QC";
  public static final String DESC = "";

  final Map<QcMetric, Requirement<String>> metricRequirements;

  public static FurtherAnalysisQCStep create(Project proj, Step plinkExportStep, Step gwasQCStep,
                                             Step ancestryStep) {
    final Requirement<Step> plinkExportStepReq = new Requirement.StepRequirement(plinkExportStep);
    final Requirement<Step> gwasQCStepReq = new Requirement.StepRequirement(gwasQCStep);
    final Requirement<Step> ancestryStepReq = new Requirement.StepRequirement(ancestryStep);
    final Requirement<File> unrelatedsFileReq = new Requirement.FileRequirement("File with list of unrelated FID/IID pairs to use for marker QC",
                                                                                new File(""));
    final Requirement<File> europeansFilesReq = new Requirement.FileRequirement("File with list of European samples to use for Hardy-Weinberg equilibrium tests",
                                                                                new File(""));
    final RequirementSet reqSet = RequirementSetBuilder.and().add(plinkExportStepReq)
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(gwasQCStepReq)
                                                                                 .add(unrelatedsFileReq))
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(ancestryStepReq)
                                                                                 .add(europeansFilesReq));

    Map<QcMetric, Requirement<String>> metricRequirements = Maps.newEnumMap(QcMetric.class);
    for (QcMetric metric : QcMetric.values()) {
      Map<QcMetric, String> defaultThresholds = FurtherAnalysisQc.getDefaultMarkerQCThresholds(proj.getArrayType());
      String defaultVal = defaultThresholds.get(metric);
      final Requirement<String> metricReq = new Requirement.ThresholdRequirement(metric.getUserDescription(),
                                                                                 defaultVal);
      reqSet.add(metricReq);
      metricRequirements.put(metric, metricReq);
    }

    return new FurtherAnalysisQCStep(proj, metricRequirements, unrelatedsFileReq, europeansFilesReq,
                                     reqSet);
  }

  final Project proj;
  final Requirement<File> unrelatedsFileReq;
  final Requirement<File> europeansFilesReq;

  private FurtherAnalysisQCStep(Project proj, Map<QcMetric, Requirement<String>> metricReqs,
                                Requirement<File> unrelReq, Requirement<File> euroReq,
                                RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.noneOf(Requirement.Flag.class));
    this.proj = proj;
    this.unrelatedsFileReq = unrelReq;
    this.europeansFilesReq = euroReq;
    this.metricRequirements = metricReqs;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    // not needed for step
  }

  @Override
  public void run(Variables variables) {

    String unrelatedsFile = resolveUnrelatedsFile(variables);

    String europeansFile = resolveEuropeansFile(variables);

    Map<QcMetric, String> markerQCThresholds = Maps.newEnumMap(QcMetric.class);
    for (QcMetric metric : QcMetric.values()) {
      Requirement<String> req = metricRequirements.get(metric);
      markerQCThresholds.put(metric, variables.get(req));
    }
    new FurtherAnalysisQc(GenvisisWorkflow.getPlinkDir(proj), GenvisisWorkflow.PLINKROOT,
                          markerQCThresholds, unrelatedsFile, europeansFile, proj.getLog())
                                                                                           .runFurtherAnalysisQC();
  }

  @Override
  public String getCommandLine(Variables variables) {
    String unrelatedsFile = resolveUnrelatedsFile(variables);

    String europeansFile = resolveEuropeansFile(variables);

    List<String> commandChunks = Lists.newArrayList();
    commandChunks.add(Files.getRunString());
    commandChunks.add(FurtherAnalysisQc.class.getName());
    commandChunks.add(CLI.formCmdLineArg(FurtherAnalysisQc.ARG_UNRELATEDS, unrelatedsFile));
    commandChunks.add(CLI.formCmdLineArg(FurtherAnalysisQc.ARG_EUROPEANS, europeansFile));
    commandChunks.add(CLI.formCmdLineArg(CLI.ARG_INDIR, GenvisisWorkflow.getPlinkDir(proj)));
    commandChunks.add(CLI.formCmdLineArg(CLI.ARG_PLINKROOT, GenvisisWorkflow.PLINKROOT));
    for (QcMetric metric : QcMetric.values()) {
      Requirement<String> req = metricRequirements.get(metric);
      commandChunks.add(CLI.formCmdLineArg(metric.getKey(), variables.get(req)));
    }
    return Joiner.on(' ').join(commandChunks);
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    String dir = GenvisisWorkflow.getPlinkDir(proj) + Qc.QC_SUBDIR
                 + FurtherAnalysisQc.FURTHER_ANALYSIS_DIR;
    String qcdPlinkroot = GenvisisWorkflow.PLINKROOT
                          + FurtherAnalysisQc.FURTHER_ANALYSIS_QC_PLINK_SUFFIX;
    return PSF.Plink.bedBimFamExist(dir + qcdPlinkroot)
           && Files.exists(dir + FurtherAnalysisQc.SAMPLE_QC_DROPS, false)
           && Files.exists(dir + FurtherAnalysisQc.MARKER_QC_DROPS, false);
  }

  private String resolveUnrelatedsFile(Variables stepVars) {
    String unrelatedsFile = stepVars.get(unrelatedsFileReq).getPath();
    if (unrelatedsFile == null || unrelatedsFile.equals("") || !Files.exists(unrelatedsFile)) {
      unrelatedsFile = GenvisisWorkflow.getAncestryDir(proj)
                       + RelationAncestryQc.UNRELATEDS_FILENAME;
    }
    return unrelatedsFile;
  }

  private String resolveEuropeansFile(Variables stepVars) {
    String europeansFile = stepVars.get(europeansFilesReq).getPath();
    if (europeansFile == null || "".equals(europeansFile) || !Files.exists(europeansFile)) {
      String raceImputationFilename = GenvisisWorkflow.getAncestryDir(proj)
                                      + Ancestry.RACE_IMPUTATIONS_FILENAME;
      europeansFile = PCImputeRace.formRaceListFilename(RACE.WHITE, raceImputationFilename);
    }
    return europeansFile;
  }

}