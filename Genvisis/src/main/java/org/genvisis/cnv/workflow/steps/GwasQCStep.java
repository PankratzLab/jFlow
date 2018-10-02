package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.Variables;
import org.pankratzlab.common.Files;
import org.pankratzlab.core.CLI;
import org.pankratzlab.gwas.MarkerQC;
import org.pankratzlab.gwas.PlinkMendelianChecker;
import org.pankratzlab.gwas.RelationAncestryQc;
import org.pankratzlab.gwas.MarkerQC.QC_METRIC;
import org.pankratzlab.shared.gwas.Qc;
import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class GwasQCStep extends Step {

  public static final String NAME = "Run GWAS QC";
  public static final String DESC = "";

  public static GwasQCStep create(Project proj, Step plinkExportStep) {
    final Requirement<Step> plinkExportStepReq = new Requirement.StepRequirement(plinkExportStep);
    String defaultCallrate;
    switch (proj.getArrayType()) {
      case AFFY_GW6:
      case AFFY_GW6_CN:
        defaultCallrate = MarkerQC.DEFAULT_AFFY_CALLRATE_THRESHOLD;
        break;
      case AFFY_AXIOM:
      case ILLUMINA:
        defaultCallrate = MarkerQC.DEFAULT_ILLUMINA_CALLRATE_THRESHOLD;
        break;
      default:
        throw new IllegalArgumentException("Invalid " + proj.getArrayType().getClass().getName()
                                           + ": " + proj.getArrayType().toString());
    }
    final Requirement<String> callrateReq = new Requirement.ThresholdRequirement(QC_METRIC.CALLRATE.getUserDescription(),
                                                                                 defaultCallrate);
    final RequirementSet reqSet = RequirementSetBuilder.and().add(plinkExportStepReq)
                                                       .add(callrateReq);

    return new GwasQCStep(proj, callrateReq, reqSet);
  }

  final Project proj;
  final Requirement<String> callrateReq;

  private GwasQCStep(Project proj, Requirement<String> callrateReq, RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.noneOf(Requirement.Flag.class));
    this.proj = proj;
    this.callrateReq = callrateReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    // not needed for step
  }

  @Override
  public void run(Variables variables) {
    String dir = GenvisisWorkflow.getPlinkDir(proj);
    Map<QC_METRIC, String> markerQCThresholds = Maps.newEnumMap(RelationAncestryQc.DEFAULT_QC_METRIC_THRESHOLDS);
    markerQCThresholds.put(QC_METRIC.CALLRATE, variables.get(callrateReq));
    new RelationAncestryQc(dir, GenvisisWorkflow.PLINKROOT, markerQCThresholds,
                           proj.getLog()).run(false);
    if (new File(dir + Qc.QC_SUBDIR + RelationAncestryQc.GENOME_DIR + GenvisisWorkflow.PLINKROOT
                 + ".genome").exists()) {
      proj.GENOME_CLUSTER_FILENAME.setValue(dir + Qc.QC_SUBDIR + RelationAncestryQc.GENOME_DIR
                                            + GenvisisWorkflow.PLINKROOT + ".genome");
      proj.saveProperties();
    }
    new PlinkMendelianChecker(proj).run();
  }

  @Override
  public String getCommandLine(Variables variables) {
    String dir = GenvisisWorkflow.getPlinkDir(proj);
    Map<QC_METRIC, String> markerQCThresholds = Maps.newEnumMap(RelationAncestryQc.DEFAULT_QC_METRIC_THRESHOLDS);
    markerQCThresholds.put(QC_METRIC.CALLRATE, variables.get(callrateReq));

    List<String> commandChunks = Lists.newArrayList();
    commandChunks.add(Files.getRunString());
    commandChunks.add(RelationAncestryQc.class.getName());
    commandChunks.add(CLI.formCmdLineArg(CLI.ARG_INDIR, GenvisisWorkflow.getPlinkDir(proj)));
    commandChunks.add(CLI.formCmdLineArg(CLI.ARG_PLINKROOT, GenvisisWorkflow.PLINKROOT));
    commandChunks.add(CLI.formCmdLineArg(RelationAncestryQc.ARGS_KEEPGENOME, "false"));
    commandChunks.add(CLI.formCmdLineArg(QC_METRIC.CALLRATE.getKey(), variables.get(callrateReq)));
    commandChunks.add("\n" + Files.getRunString());
    commandChunks.add(GenvisisWorkflow.PROJ_PROP_UPDATE_STR + proj.getPropertyFilename());
    commandChunks.add(proj.GENOME_CLUSTER_FILENAME.getName() + "=" + dir + Qc.QC_SUBDIR
                      + RelationAncestryQc.GENOME_DIR + GenvisisWorkflow.PLINKROOT + ".genome");
    commandChunks.add("\n" + Files.getRunString());
    commandChunks.add(PlinkMendelianChecker.class.getName());
    commandChunks.add("proj=" + proj.getPropertyFilename());
    return Joiner.on(' ').join(commandChunks);
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    String dir = GenvisisWorkflow.getPlinkDir(proj);
    for (int i = 0; i < org.pankratzlab.gwas.RelationAncestryQc.FOLDERS_CREATED.length; i++) {
      for (int j = 0; j < org.pankratzlab.gwas.RelationAncestryQc.FILES_CREATED[i].length; j++) {
        if (!Files.exists(dir + org.pankratzlab.gwas.RelationAncestryQc.FOLDERS_CREATED[i]
                          + org.pankratzlab.gwas.RelationAncestryQc.FILES_CREATED[i][j])) {
          return false;
        }
      }
    }
    return Files.checkAllFiles(PlinkMendelianChecker.parseOutputDirectory(proj),
                               PlinkMendelianChecker.OUTPUTS, false, proj.getLog());
  }

}
