package org.genvisis.cnv.workflow.steps;

import java.io.File;
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
import org.genvisis.cnv.workflow.Requirement.OptionalFileRequirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.cnv.workflow.Variables;
import org.genvisis.common.Files;

public class PCCorrectionStep extends Step {

  public static final String NAME = "Create PC-Corrected Project";
  public static final String DESC = "";

  final Project proj;
  final Requirement<Integer> numPCsReq;
  final Requirement<File> outputBaseReq;
  final Requirement<Double> callrateReq;
  final Requirement<Boolean> recomputeLrrReq;
  final Requirement<File> tempDirReq;
  final Requirement<CORRECTION_TYPE> correctionStrategyReq;
  final Requirement<CHROMOSOME_X_STRATEGY> sexChromosomeStrategyReq;
  final Requirement<Boolean> setupCNVCalling;
  final Requirement<Integer> numThreadsReq;

  public static PCCorrectionStep create(Project proj, Step parseSamplesStep,
                                        Requirement<Integer> numThreadsReq) {
    final Requirement<Step> parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    final Requirement<Integer> numPCsReq = new Requirement.PosIntRequirement("Number of principal components for correction.",
                                                                             MitoPipeline.DEFAULT_NUM_COMPONENTS);
    final Requirement<File> outputBaseReq = new OptionalFileRequirement("Output file path (relative to project directory) and baseName for principal components correction files",
                                                                        new File(proj.PROJECT_DIRECTORY.getValue()
                                                                                 + MitoPipeline.FILE_BASE)) {

      @Override
      public boolean checkRequirement(String arg, Set<Step> stepSelections,
                                      Map<Step, Variables> variables) {
        String outputBase = proj.PROJECT_DIRECTORY.getValue() + arg;
        String finalReport = outputBase + PCA.FILE_EXTs[0];
        return super.checkRequirement(finalReport, stepSelections, variables);
      }
    };
    final Requirement<Double> callrateReq = new Requirement.DoubleRequirement("Call-rate filter for determining high-quality markers",
                                                                              MitoPipeline.DEFAULT_MKR_CALLRATE_FILTER,
                                                                              0.0, 1.0);
    final Requirement<Boolean> recomputeLrrReq = new Requirement.OptionalBoolRequirement("Re-compute Log-R Ratio values? (usually false if LRRs already exist)",
                                                                                         false);
    final Requirement<File> tempDirReq = new Requirement.OptionalFileRequirement("Temporary directory for intermediate files (which tend to be very large)",
                                                                                 new File(""));
    final Requirement<CORRECTION_TYPE> correctionStrategyReq = new Requirement.EnumRequirement<CORRECTION_TYPE>("Correction Type",
                                                                                                                CORRECTION_TYPE.XY);
    final Requirement<CHROMOSOME_X_STRATEGY> sexChromosomeStrategyReq = new Requirement.EnumRequirement<CHROMOSOME_X_STRATEGY>("Sex Chromosome Strategy",
                                                                                                                               CHROMOSOME_X_STRATEGY.BIOLOGICAL);
    final Requirement<Boolean> setupCNVCalling = new Requirement.OptionalBoolRequirement("Create script with steps to process corrected data and call CNVs?",
                                                                                         false);

    return new PCCorrectionStep(proj, parseSamplesStepReq, numPCsReq, outputBaseReq, callrateReq,
                                recomputeLrrReq, tempDirReq, correctionStrategyReq,
                                sexChromosomeStrategyReq, setupCNVCalling, numThreadsReq);
  }

  private static RequirementSet createReqSet(Requirement<Step> parseSamplesStepReq,
                                             Requirement<Integer> numPCsReq,
                                             Requirement<File> outputBaseReq,
                                             Requirement<Double> callrateReq,
                                             Requirement<Boolean> recomputeLrrReq,
                                             Requirement<File> tempDirReq,
                                             Requirement<CORRECTION_TYPE> correctionStrategyReq,
                                             Requirement<CHROMOSOME_X_STRATEGY> sexChromosomeStrategyReq,
                                             Requirement<Boolean> setupCNVCalling,
                                             Requirement<Integer> numThreadsReq) {
    return RequirementSetBuilder.and().add(parseSamplesStepReq).add(numPCsReq).add(outputBaseReq)
                                .add(callrateReq).add(recomputeLrrReq).add(tempDirReq)
                                .add(correctionStrategyReq).add(sexChromosomeStrategyReq)
                                .add(numThreadsReq).add(setupCNVCalling);
  }

  private PCCorrectionStep(Project proj, Requirement<Step> parseSamplesStepReq,
                           Requirement<Integer> numPCsReq, Requirement<File> outputBaseReq,
                           Requirement<Double> callrateReq, Requirement<Boolean> recomputeLrrReq,
                           Requirement<File> tempDirReq,
                           Requirement<CORRECTION_TYPE> correctionStrategyReq,
                           Requirement<CHROMOSOME_X_STRATEGY> sexChromosomeStrategyReq,
                           Requirement<Boolean> setupCNVCalling,
                           Requirement<Integer> numThreadsReq) {
    super(NAME, DESC,
          createReqSet(parseSamplesStepReq, numPCsReq, outputBaseReq, callrateReq, recomputeLrrReq,
                       tempDirReq, correctionStrategyReq, sexChromosomeStrategyReq, setupCNVCalling,
                       numThreadsReq),
          EnumSet.of(Requirement.Flag.MEMORY, Requirement.Flag.RUNTIME,
                     Requirement.Flag.MULTITHREADED));
    this.proj = proj;
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
  public void setNecessaryPreRunProperties(Variables variables) {
    // not needed for step
  }

  @Override
  public void run(Variables variables) {
    int numComponents = variables.get(numPCsReq);
    File outputBase = variables.get(outputBaseReq);
    double markerCallRateFilter = variables.get(callrateReq);
    boolean recomputeLRRPCs = variables.get(recomputeLrrReq);
    String tmpDir = variables.hasValid(tempDirReq) ? variables.get(tempDirReq).getPath() : null;
    CORRECTION_TYPE type = variables.get(correctionStrategyReq);
    CHROMOSOME_X_STRATEGY strategy = variables.get(sexChromosomeStrategyReq);

    int totalThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
    boolean cnvCalling = variables.get(setupCNVCalling);
    String retMsg = PRoCtOR.shadow(proj, tmpDir, outputBase.getPath(), markerCallRateFilter,
                                   recomputeLRRPCs, type, strategy, numComponents, totalThreads,
                                   cnvCalling);
    if (retMsg != null && !"".equals(retMsg)) {
      throw new RuntimeException(retMsg);
    }
  }

  @Override
  public String getCommandLine(Variables variables) {
    int numComponents = variables.get(numPCsReq);
    String outputBase = variables.get(outputBaseReq).getPath();
    double markerCallRateFilter = variables.get(callrateReq);
    boolean recomputeLRRPCs = variables.get(recomputeLrrReq);
    String tmpDir = variables.hasValid(tempDirReq) ? variables.get(tempDirReq).getAbsolutePath()
                                                   : null;
    CORRECTION_TYPE correctionType = (CORRECTION_TYPE) variables.get(correctionStrategyReq);
    CHROMOSOME_X_STRATEGY strategy = (CHROMOSOME_X_STRATEGY) variables.get(sexChromosomeStrategyReq);

    int totalThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));

    boolean cnvCalling = variables.get(setupCNVCalling);
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
  public boolean checkIfOutputExists(Variables variables) {
    String outputBase = proj.PROJECT_DIRECTORY.getValue() + variables.get(outputBaseReq).getPath();
    String finalReport = outputBase + PCA.FILE_EXTs[0];
    return Files.exists(finalReport);
  }

}
