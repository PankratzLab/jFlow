package org.genvisis.cnv.workflow.steps;

import java.util.EnumSet;
import java.util.Map;
import org.genvisis.cnv.analysis.pca.PCA;
import org.genvisis.cnv.analysis.pca.PCAPrep;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsCompute.PRE_PROCESSING_METHOD;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.MitoPipeline;
import org.genvisis.cnv.qc.GcAdjustor;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.common.Files;
import org.genvisis.common.PSF;

public class MitoCNEstimateStep extends Step {

  public static final String NAME = "Create Mitochondrial Copy-Number Estimates File";
  public static final String DESC = "";

  public static MitoCNEstimateStep create(Project proj, Step transposeStep,
                                          Requirement numThreadsReq, double priority) {
    // FIXME http://genvisis.org/MitoPipeline/#illumina_marker_lists has illumina markers.. this
    // should be linked to, or
    // these steps split or something...
    final Requirement transposeStepReq = new Requirement.StepRequirement(transposeStep);
    final Requirement medianMarkersReq = new Requirement.FileRequirement("MedianMarkers file must exist.",
                                                                         "");
    final Requirement lrrSdThresholdReq = new Requirement.DoubleRequirement("LRR SD threshold to filter samples.",
                                                                            proj.LRRSD_CUTOFF.getValue(),
                                                                            proj.LRRSD_CUTOFF.getMinValue(),
                                                                            proj.LRRSD_CUTOFF.getMaxValue());
    final Requirement callrateThresholdReq = new Requirement.DoubleRequirement("Call rate threshold to filter markers.",
                                                                               MitoPipeline.DEFAULT_MKR_CALLRATE_FILTER,
                                                                               0.0, 1.0);
    final Requirement qcPassingOnlyReq = new Requirement.OptionalBoolRequirement("Compute PCs with samples passing QC only",
                                                                                 true);
    final Requirement imputeNaNs = new Requirement.OptionalBoolRequirement("Impute mean value for NaN",
                                                                           true);
    final Requirement recomputeLrrPCMarkersReq = new Requirement.OptionalBoolRequirement("Should recompute Log-R ratio for PC markers?",
                                                                                         true);
    final Requirement recomputeLrrMedianMarkersReq = new Requirement.OptionalBoolRequirement("Should recompute Log-R ratio for median markers?",
                                                                                             true);
    final Requirement homozygousOnlyReq = new Requirement.OptionalBoolRequirement("Homozygous only?",
                                                                                  true);
    final Requirement gcRegressionDistanceReq = new Requirement.PosIntRequirement("Regression distance for the GC adjustment",
                                                                                  GcAdjustor.DEFAULT_REGRESSION_DISTANCE[0]);
    final Requirement pcSelectionSamplesReq = new Requirement.OptionalFileRequirement("A file listing a subset of samples (DNA ID) to use for determining optimal PC selection, typically a list of unrelated and single race samples. If a list is not provided, only samples passing sample qc thresholds will be used.",
                                                                                      "");
    final Requirement externalBetaFileReq = new Requirement.OptionalFileRequirement("An external beta file to optimize PC selection.",
                                                                                    "");

    final RequirementSet reqSet = RequirementSetBuilder.and().add(transposeStepReq)
                                                       .add(medianMarkersReq).add(lrrSdThresholdReq)
                                                       .add(callrateThresholdReq)
                                                       .add(qcPassingOnlyReq).add(imputeNaNs)
                                                       .add(recomputeLrrPCMarkersReq)
                                                       .add(recomputeLrrMedianMarkersReq)
                                                       .add(homozygousOnlyReq)
                                                       .add(gcRegressionDistanceReq)
                                                       .add(numThreadsReq)
                                                       .add(pcSelectionSamplesReq)
                                                       .add(externalBetaFileReq);
    return new MitoCNEstimateStep(medianMarkersReq, lrrSdThresholdReq, callrateThresholdReq,
                                  qcPassingOnlyReq, imputeNaNs, recomputeLrrPCMarkersReq,
                                  recomputeLrrMedianMarkersReq, homozygousOnlyReq,
                                  gcRegressionDistanceReq, pcSelectionSamplesReq,
                                  externalBetaFileReq, numThreadsReq, reqSet, priority);
  }

  final Requirement medianMarkersReq;
  final Requirement lrrSdThresholdReq;
  final Requirement callrateThresholdReq;
  final Requirement qcPassingOnlyReq;
  final Requirement imputeNaNs;
  final Requirement recomputeLrrPCMarkersReq;
  final Requirement recomputeLrrMedianMarkersReq;
  final Requirement homozygousOnlyReq;
  final Requirement gcRegressionDistanceReq;
  final Requirement pcSelectionSamplesReq;
  final Requirement externalBetaFileReq;
  final Requirement numThreadsReq;

  private MitoCNEstimateStep(Requirement medianMarkersReq, Requirement lrrSdThresholdReq,
                             Requirement callrateThresholdReq, Requirement qcPassingOnlyReq,
                             Requirement imputeNaNs, Requirement recomputeLrrPCMarkersReq,
                             Requirement recomputeLrrMedianMarkersReq,
                             Requirement homozygousOnlyReq, Requirement gcRegressionDistanceReq,
                             Requirement pcSelectionSamplesReq, Requirement externalBetaFileReq,
                             Requirement numThreadsReq, RequirementSet reqSet, double priority) {
    super(NAME, DESC, reqSet, EnumSet.of(Requirement.Flag.MULTITHREADED), priority);
    this.medianMarkersReq = medianMarkersReq;
    this.lrrSdThresholdReq = lrrSdThresholdReq;
    this.callrateThresholdReq = callrateThresholdReq;
    this.qcPassingOnlyReq = qcPassingOnlyReq;
    this.imputeNaNs = imputeNaNs;
    this.recomputeLrrPCMarkersReq = recomputeLrrPCMarkersReq;
    this.recomputeLrrMedianMarkersReq = recomputeLrrMedianMarkersReq;
    this.homozygousOnlyReq = homozygousOnlyReq;
    this.gcRegressionDistanceReq = gcRegressionDistanceReq;
    this.pcSelectionSamplesReq = pcSelectionSamplesReq;
    this.externalBetaFileReq = externalBetaFileReq;
    this.numThreadsReq = numThreadsReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Project proj,
                                           Map<Step, Map<Requirement, String>> variables) {
    double sampleLRRSdFilter = Double.parseDouble(variables.get(this).get(lrrSdThresholdReq));
    if (sampleLRRSdFilter < 0) {
      switch (proj.ARRAY_TYPE.getValue()) {
        case AFFY_GW6:
        case AFFY_GW6_CN:
          proj.LRRSD_CUTOFF.setValue(0.35);
          proj.getLog()
              .reportTimeInfo("Setting " + proj.LRRSD_CUTOFF.getName()
                              + " to default 0.35 for array " + proj.ARRAY_TYPE.getValue());
          break;
        case ILLUMINA:
          proj.LRRSD_CUTOFF.setValue(0.30);
          proj.getLog()
              .reportTimeInfo("Setting " + proj.LRRSD_CUTOFF.getName()
                              + " to default 0.30 for array " + proj.ARRAY_TYPE.getValue());
          break;
        default:
          throw new IllegalArgumentException("Invalid Array type");
      }
    } else {
      proj.LRRSD_CUTOFF.setValue(sampleLRRSdFilter);
    }
  }

  @Override
  public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String medianMarkers = variables.get(this).get(medianMarkersReq);
    double markerCallRateFilter = Double.parseDouble(variables.get(this).get(callrateThresholdReq));
    // FIXME: This gcCorrect assignment was carried over from the old indexed version but
    // appears incorrect
    boolean gcCorrect = Boolean.parseBoolean(variables.get(this).get(qcPassingOnlyReq));
    boolean imputeMeanForNaN = Boolean.parseBoolean(variables.get(this).get(imputeNaNs));
    boolean recomputeLRRPCs = Boolean.parseBoolean(variables.get(this)
                                                            .get(recomputeLrrPCMarkersReq));
    boolean recomputeLRRMedian = Boolean.parseBoolean(variables.get(this)
                                                               .get(recomputeLrrMedianMarkersReq));
    boolean homozygousOnly = Boolean.parseBoolean(variables.get(this).get(homozygousOnlyReq));
    int bpGcModel = GcAdjustor.GcModel.DEFAULT_GC_MODEL_BIN_FASTA;
    int regressionDistance = Integer.parseInt(variables.get(this).get(gcRegressionDistanceReq));
    int numComponents = MitoPipeline.DEFAULT_NUM_COMPONENTS;
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(this).get(numThreadsReq));
    String outputBase = MitoPipeline.FILE_BASE;

    String betaOptFile = variables.get(this).get(pcSelectionSamplesReq);
    String betaFile = variables.get(this).get(externalBetaFileReq);

    boolean markerQC = true;
    double[] pvalOpt = MitoPipeline.DEFAULT_PVAL_OPTS;
    String pedFile = null;
    String useFile = null;
    boolean sampLrr = true;
    boolean plot = false;
    int retCode = PCAPrep.prepPCA(proj, numThreads, outputBase, markerQC, markerCallRateFilter,
                                  useFile, proj.getSampleList(), proj.getLog());
    if (retCode == PCAPrep.SUCCESS_CODE) {
      MitoPipeline.estimateMtDNACN(proj, numThreads, medianMarkers, numComponents, outputBase,
                                   homozygousOnly, markerCallRateFilter, betaOptFile, pedFile,
                                   recomputeLRRPCs, recomputeLRRMedian, sampLrr, imputeMeanForNaN,
                                   gcCorrect, bpGcModel, regressionDistance,
                                   proj.GENOME_BUILD_VERSION.getValue(), pvalOpt, betaFile, plot,
                                   false, PRE_PROCESSING_METHOD.NONE, proj.getLog());
    } else {
      throw new RuntimeException(PCAPrep.errorMessage(retCode));
    }
  }

  @Override
  public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String medianMarkers = variables.get(this).get(medianMarkersReq);
    double lrrSD = Double.parseDouble(variables.get(this).get(lrrSdThresholdReq));
    double markerCallRateFilter = Double.parseDouble(variables.get(this).get(callrateThresholdReq));
    // FIXME: This gcCorrect assignment was carried over from the old indexed version but
    // appears incorrect
    boolean gcCorrect = Boolean.parseBoolean(variables.get(this).get(qcPassingOnlyReq));
    boolean imputeMeanForNaN = Boolean.parseBoolean(variables.get(this).get(imputeNaNs));
    boolean recomputeLRRPCs = Boolean.parseBoolean(variables.get(this)
                                                            .get(recomputeLrrPCMarkersReq));
    boolean recomputeLRRMedian = Boolean.parseBoolean(variables.get(this)
                                                               .get(recomputeLrrMedianMarkersReq));
    boolean homozygousOnly = Boolean.parseBoolean(variables.get(this).get(homozygousOnlyReq));
    int bpGcModel = GcAdjustor.GcModel.DEFAULT_GC_MODEL_BIN_FASTA;
    int regressionDistance = Integer.parseInt(variables.get(this).get(gcRegressionDistanceReq));
    int numComponents = MitoPipeline.DEFAULT_NUM_COMPONENTS;
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(this).get(numThreadsReq));
    String outputBase = MitoPipeline.FILE_BASE;

    String betaOptFile = variables.get(this).get(pcSelectionSamplesReq);
    String betaFile = variables.get(this).get(externalBetaFileReq);
    boolean sampLrr = true;

    String projPropFile = proj.getPropertyFilename();
    StringBuilder cmd = new StringBuilder();
    cmd.append(Files.getRunString()).append(" org.genvisis.cnv.manage.MitoPipeline")
       .append(" proj=").append(projPropFile).append(" mitochondrialMarkers=").append(medianMarkers)
       .append(" numComponents=").append(numComponents).append(" imputeMeanForNaN=")
       .append(imputeMeanForNaN).append(" recomputeLRR_PCs=").append(recomputeLRRPCs)
       .append(" recomputeLRR_Median=").append(recomputeLRRMedian).append(" gcCorrect=")
       .append(gcCorrect).append(" bpGcModel=").append(bpGcModel).append(" LRRSD=").append(lrrSD)
       .append(" markerCallRate=").append(markerCallRateFilter).append(" regressionDistance=")
       .append(regressionDistance).append(" sampLRR=").append(sampLrr).append(" ")
       .append(PSF.Ext.NUM_THREADS_COMMAND).append(numThreads).append(" log=")
       .append(proj.getLog().getFilename()).append(" output=").append(outputBase);
    if (!"".equals(betaOptFile)) {
      cmd.append(" ").append(MitoPipeline.PC_OPT_FILE).append("=").append(betaOptFile);
    }
    if (!"".equals(betaFile)) {
      cmd.append(" betas=").append(betaFile);
    }
    if (!homozygousOnly) {
      cmd.append(" -allCalls ");
    }

    cmd.append(" -SkipProjectCreationWithLongUndocumentedFlag ");

    return cmd.toString();
  }

  @Override
  public boolean checkIfOutputExists(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String outputBase = proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.FILE_BASE;
    String finalReport = outputBase + PCA.FILE_EXTs[0];
    return Files.exists(finalReport);
  }

}
