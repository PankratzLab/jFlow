package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import org.genvisis.CLI;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.common.Files;
import org.genvisis.gwas.MarkerQC;
import org.genvisis.gwas.MarkerQC.QC_METRIC;
import org.genvisis.gwas.PlinkMendelianChecker;
import org.genvisis.gwas.Qc;
import org.genvisis.gwas.RelationAncestryQc;
import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class GwasQCStep extends Step {

  public static final String NAME = "Run GWAS QC";
  public static final String DESC = "";

  public static GwasQCStep create(Project proj, Step plinkExportStep) {
    final Requirement plinkExportStepReq = new Requirement.StepRequirement(plinkExportStep);
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
    final Requirement callrateReq = new Requirement.ThresholdRequirement(QC_METRIC.CALLRATE.getUserDescription(),
                                                                         defaultCallrate);
    final RequirementSet reqSet = RequirementSetBuilder.and().add(plinkExportStepReq)
                                                       .add(callrateReq);

    return new GwasQCStep(callrateReq, reqSet);
  }

  final Requirement callrateReq;

  private GwasQCStep(Requirement callrateReq, RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.noneOf(Requirement.Flag.class));
    this.callrateReq = callrateReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Project proj,
                                           Map<Step, Map<Requirement, String>> variables) {
    // not needed for step
  }

  @Override
  public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String dir = GenvisisWorkflow.getPlinkDir(proj);
    Map<QC_METRIC, String> markerQCThresholds = Maps.newEnumMap(RelationAncestryQc.DEFAULT_QC_METRIC_THRESHOLDS);
    markerQCThresholds.put(QC_METRIC.CALLRATE, variables.get(this).get(callrateReq));
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
  public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
    Map<Requirement, String> stepVars = variables.get(this);

    String dir = GenvisisWorkflow.getPlinkDir(proj);
    Map<QC_METRIC, String> markerQCThresholds = Maps.newEnumMap(RelationAncestryQc.DEFAULT_QC_METRIC_THRESHOLDS);
    markerQCThresholds.put(QC_METRIC.CALLRATE, variables.get(this).get(callrateReq));

    List<String> commandChunks = Lists.newArrayList();
    commandChunks.add(Files.getRunString());
    commandChunks.add(RelationAncestryQc.class.getName());
    commandChunks.add(CLI.formCmdLineArg(CLI.ARG_INDIR, GenvisisWorkflow.getPlinkDir(proj)));
    commandChunks.add(CLI.formCmdLineArg(CLI.ARG_PLINKROOT, GenvisisWorkflow.PLINKROOT));
    commandChunks.add(CLI.formCmdLineArg(RelationAncestryQc.ARGS_KEEPGENOME, "false"));
    commandChunks.add(CLI.formCmdLineArg(QC_METRIC.CALLRATE.getKey(), stepVars.get(callrateReq)));
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
  public boolean checkIfOutputExists(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String dir = GenvisisWorkflow.getPlinkDir(proj);
    for (int i = 0; i < org.genvisis.gwas.RelationAncestryQc.FOLDERS_CREATED.length; i++) {
      for (int j = 0; j < org.genvisis.gwas.RelationAncestryQc.FILES_CREATED[i].length; j++) {
        if (!Files.exists(dir + org.genvisis.gwas.RelationAncestryQc.FOLDERS_CREATED[i]
                          + org.genvisis.gwas.RelationAncestryQc.FILES_CREATED[i][j])) {
          return false;
        }
      }
    }
    return Files.checkAllFiles(PlinkMendelianChecker.parseOutputDirectory(proj),
                               PlinkMendelianChecker.OUTPUTS, false, proj.getLog());
  }

}
