package org.genvisis.cnv.workflow.steps;

import java.util.EnumSet;
import java.util.Map;
import java.util.Set;
import org.genvisis.cnv.analysis.pca.PCA;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity.CHROMOSOME_X_STRATEGY;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity.CORRECTION_TYPE;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.MitoPipeline;
import org.genvisis.cnv.manage.PRoCtOR;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.common.Files;

public class PCCorrectionStep extends Step {

  public static final String NAME = "Create PC-Corrected Project";
  public static final String DESC = "";

  final Requirement numPCsReq;
  final Requirement outputBaseReq;
  final Requirement callrateReq;
  final Requirement recomputeLrrReq;
  final Requirement tempDirReq;
  final Requirement correctionStrategyReq;
  final Requirement sexChromosomeStrategyReq;
  final Requirement setupCNVCalling;
  final Requirement numThreadsReq;

  public static PCCorrectionStep create(Step parseSamplesStep, Requirement numThreadsReq,
                                        double priority) {
    final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    final Requirement numPCsReq = new Requirement.PosIntRequirement("Number of principal components for correction.",
                                                                    MitoPipeline.DEFAULT_NUM_COMPONENTS);
    final Requirement outputBaseReq = new Requirement.OutputFileRequirement("Output file path (relative to project directory) and baseName for principal components correction files",
                                                                            MitoPipeline.FILE_BASE) {

      @Override
      public boolean checkRequirement(Project proj, String arg, Set<Step> stepSelections,
                                      Map<Step, Map<Requirement, String>> variables) {
        String outputBase = proj.PROJECT_DIRECTORY.getValue() + arg;
        String finalReport = outputBase + PCA.FILE_EXTs[0];
        return super.checkRequirement(proj, finalReport, stepSelections, variables);
      }
    };
    final Requirement callrateReq = new Requirement.DoubleRequirement("Call-rate filter for determining high-quality markers",
                                                                      MitoPipeline.DEFAULT_MKR_CALLRATE_FILTER,
                                                                      0.0, 1.0);
    final Requirement recomputeLrrReq = new Requirement.OptionalBoolRequirement("Re-compute Log-R Ratio values? (usually false if LRRs already exist)",
                                                                                false);
    final Requirement tempDirReq = new Requirement.OptionalFileRequirement("Temporary directory for intermediate files (which tend to be very large)",
                                                                           "");
    final Requirement correctionStrategyReq = new Requirement.EnumRequirement("Correction Type",
                                                                              CORRECTION_TYPE.XY);
    final Requirement sexChromosomeStrategyReq = new Requirement.EnumRequirement("Sex Chromosome Strategy",
                                                                                 CHROMOSOME_X_STRATEGY.BIOLOGICAL);
    final Requirement setupCNVCalling = new Requirement.OptionalBoolRequirement("Create script with steps to process corrected data and call CNVs?",
                                                                                false);

    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(numPCsReq).add(outputBaseReq)
                                                       .add(callrateReq).add(recomputeLrrReq)
                                                       .add(tempDirReq).add(correctionStrategyReq)
                                                       .add(sexChromosomeStrategyReq)
                                                       .add(numThreadsReq).add(setupCNVCalling);

    return new PCCorrectionStep(numPCsReq, outputBaseReq, callrateReq, recomputeLrrReq, tempDirReq,
                                correctionStrategyReq, sexChromosomeStrategyReq, setupCNVCalling,
                                numThreadsReq, reqSet, priority);

  }

  private PCCorrectionStep(Requirement numPCsReq, Requirement outputBaseReq,
                           Requirement callrateReq, Requirement recomputeLrrReq,
                           Requirement tempDirReq, Requirement correctionStrategyReq,
                           Requirement sexChromosomeStrategyReq, Requirement setupCNVCalling,
                           Requirement numThreadsReq, RequirementSet reqSet, double priority) {
    super(NAME, DESC, reqSet, EnumSet.of(Requirement.Flag.MEMORY, Requirement.Flag.RUNTIME,
                                         Requirement.Flag.MULTITHREADED),
          priority);
    this.numPCsReq = numPCsReq;
    this.outputBaseReq = outputBaseReq;
    this.callrateReq = callrateReq;
    this.recomputeLrrReq = recomputeLrrReq;
    this.tempDirReq = tempDirReq;
    this.correctionStrategyReq = correctionStrategyReq;
    this.sexChromosomeStrategyReq = sexChromosomeStrategyReq;
    this.setupCNVCalling = setupCNVCalling;
    this.numThreadsReq = numThreadsReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Project proj,
                                           Map<Step, Map<Requirement, String>> variables) {
    // not needed for step
  }

  @Override
  public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
    int numComponents = Integer.parseInt(variables.get(this).get(numPCsReq));
    String outputBase = variables.get(this).get(outputBaseReq);
    double markerCallRateFilter = Double.parseDouble(variables.get(this).get(callrateReq));
    boolean recomputeLRRPCs = Boolean.parseBoolean(variables.get(this).get(recomputeLrrReq));
    String tmpDir = variables.get(this).get(tempDirReq);
    if ("".equals(tmpDir.trim())) {
      tmpDir = null;
    }
    CORRECTION_TYPE type = CORRECTION_TYPE.valueOf(variables.get(this).get(correctionStrategyReq));
    CHROMOSOME_X_STRATEGY strategy = CHROMOSOME_X_STRATEGY.valueOf(variables.get(this)
                                                                            .get(sexChromosomeStrategyReq));

    int totalThreads = StepBuilder.resolveThreads(proj, variables.get(this).get(numThreadsReq));
    boolean cnvCalling = Boolean.parseBoolean(variables.get(this).get(setupCNVCalling));
    String retMsg = PRoCtOR.shadow(proj, tmpDir, outputBase, markerCallRateFilter, recomputeLRRPCs,
                                   type, strategy, numComponents, totalThreads, cnvCalling);
    if (retMsg != null && !"".equals(retMsg)) {
      throw new RuntimeException(retMsg);
    }
  }

  @Override
  public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
    int numComponents = Integer.parseInt(variables.get(this).get(numPCsReq));
    String outputBase = variables.get(this).get(outputBaseReq);
    double markerCallRateFilter = Double.parseDouble(variables.get(this).get(callrateReq));
    boolean recomputeLRRPCs = Boolean.parseBoolean(variables.get(this).get(recomputeLrrReq));
    String tmpDir = variables.get(this).get(tempDirReq);
    if ("".equals(tmpDir.trim())) {
      tmpDir = null;
    }
    String correctionType = variables.get(this).get(correctionStrategyReq);
    String strategy = variables.get(this).get(sexChromosomeStrategyReq);

    int totalThreads = StepBuilder.resolveThreads(proj, variables.get(this).get(numThreadsReq));

    boolean cnvCalling = Boolean.parseBoolean(variables.get(this).get(setupCNVCalling));
    String projPropFile = proj.getPropertyFilename();
    StringBuilder cmd = new StringBuilder();
    cmd.append(Files.getRunString()).append(" org.genvisis.cnv.manage.PRoCtOR").append(" proj=")
       .append(projPropFile).append(" numComponents=").append(numComponents).append(" outputBase=")
       .append(outputBase).append(" callrate=").append(markerCallRateFilter)
       .append(" recomputeLRR=").append(recomputeLRRPCs).append(" type=").append(correctionType)
       .append(" sexStrategy=").append(strategy).append(" numThreads=").append(totalThreads);
    if (tmpDir != null) {
      cmd.append(" tmp=").append(tmpDir);
    }
    if (cnvCalling) {
      cmd.append(" -callCNVs");
    }

    return cmd.toString();
  }

  @Override
  public boolean checkIfOutputExists(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String outputBase = proj.PROJECT_DIRECTORY.getValue() + variables.get(this).get(outputBaseReq);
    String finalReport = outputBase + PCA.FILE_EXTs[0];
    return Files.exists(finalReport);
  }

}
