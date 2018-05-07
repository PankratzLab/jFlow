package org.genvisis.cnv.workflow;

import org.apache.commons.collections4.SortedBidiMap;
import org.apache.commons.collections4.bidimap.DualTreeBidiMap;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.workflow.steps.ABLookupStep;
import org.genvisis.cnv.workflow.steps.AffyMarkerBlastStep;
import org.genvisis.cnv.workflow.steps.AncestryStep;
import org.genvisis.cnv.workflow.steps.AnnotateSampleDataStep;
import org.genvisis.cnv.workflow.steps.CallCNVsStep;
import org.genvisis.cnv.workflow.steps.ComputePFBStep;
import org.genvisis.cnv.workflow.steps.FurtherAnalysisQCStep;
import org.genvisis.cnv.workflow.steps.GCModelStep;
import org.genvisis.cnv.workflow.steps.GwasQCStep;
import org.genvisis.cnv.workflow.steps.IlluminaMarkerBlastStep;
import org.genvisis.cnv.workflow.steps.IlluminaMarkerPositionsStep;
import org.genvisis.cnv.workflow.steps.MarkerQCStep;
import org.genvisis.cnv.workflow.steps.MitoCNEstimateStep;
import org.genvisis.cnv.workflow.steps.MosaicArmsStep;
import org.genvisis.cnv.workflow.steps.PCCorrectionStep;
import org.genvisis.cnv.workflow.steps.ParseSamplesStep;
import org.genvisis.cnv.workflow.steps.PlinkExportStep;
import org.genvisis.cnv.workflow.steps.SampleDataStep;
import org.genvisis.cnv.workflow.steps.SampleQCStep;
import org.genvisis.cnv.workflow.steps.SexCentroidsStep;
import org.genvisis.cnv.workflow.steps.SexChecksStep;
import org.genvisis.cnv.workflow.steps.TransposeStep;

/**
 * Helper class to minimize manual bookkeeping when instantiating steps. Each
 * {@code generateXXXXStep} method should use the {@link #priority()} method to get its priority,
 * and call {@link #register(Step)} on the constructed step.
 * <p>
 * TODO: to reduce the risk of coding mistakes, convert the priority and register methods to
 * intrinsic functions of the steps themselves
 * </p>
 */
public class StepBuilder {

  private SortedBidiMap<Double, Step> buildSteps;
  private double p;
  private final Requirement numThreadsReq;
  public static final String NUM_THREADS_DESC = "Number of threads";
  public static final String PUTATIVE_WHITE_FILE_DESCRIPTION = "File with FID/IID pairs of putative white samples";

  public StepBuilder(Project proj) {
    numThreadsReq = new Requirement.PosIntRequirement(NUM_THREADS_DESC,
                                                      proj.NUM_THREADS.getValue());
    buildSteps = new DualTreeBidiMap<>();
    p = 0.0;
  }

  public Requirement getNumThreadsReq() {
    return numThreadsReq;
  }

  /**
   * @return All steps {@link #register(Step)}ed by this step builder thus far
   */
  public SortedBidiMap<Double, Step> getSortedSteps() {
    return buildSteps;
  }

  /**
   * @return The next step priority
   */
  private double priority() {
    return ++p;
  }

  /**
   * Register the given step in the list returned by {@link #getSortedSteps()}
   */
  Step register(Step s) {
    buildSteps.put(priority(), s);
    return s;
  }

  Step generateIlluminaMarkerPositionsStep(Project proj) {
    return register(IlluminaMarkerPositionsStep.create(proj));
  }

  Step generateIlluminaMarkerBlastAnnotationStep(Project proj, final Step parseSamplesStep) {
    return register(IlluminaMarkerBlastStep.create(proj, parseSamplesStep, numThreadsReq));
  }

  Step generateAffyMarkerBlastAnnotationStep(final Step parseSamplesStep) {
    return register(AffyMarkerBlastStep.create(parseSamplesStep, numThreadsReq));
  }

  Step generateParseSamplesStep(Project proj) {
    return generateParseSamplesStep(proj, null);
  }

  Step generateParseSamplesStep(Project proj, final Step markerPositionsStep) {
    return register(ParseSamplesStep.create(proj, markerPositionsStep, numThreadsReq));
  }

  Step generateCreateSampleDataStep(Project proj, final Step parseSamplesStep) {
    return register(SampleDataStep.create(parseSamplesStep, proj));
  }

  Step generateTransposeStep(Project proj, final Step parseSamplesStep) {
    return register(TransposeStep.create(parseSamplesStep));
  }

  Step generateGCModelStep(Project proj) {
    return register(GCModelStep.create(proj));
  }

  Step generateSampleQCStep(Project proj, final Step parseSamplesStep) {
    return register(SampleQCStep.create(parseSamplesStep, numThreadsReq));
  }

  Step generateMarkerQCStep(Project proj, final Step parseSamplesStep) {
    return register(MarkerQCStep.create(proj, parseSamplesStep, numThreadsReq));
  }

  Step generateSexChecksStep(Project proj, final Step parseSamplesStep, final Step markerBlastStep,
                             final Step sampleDataStep, final Step transposeStep,
                             final Step sampleQCStep) {
    return register(SexChecksStep.create(proj, parseSamplesStep, markerBlastStep, sampleDataStep,
                                         transposeStep, sampleQCStep));
  }

  Step generatePlinkExportStep(Project proj, final Step parseSamplesStep) {
    return register(PlinkExportStep.create(proj, parseSamplesStep));
  }

  Step generateGwasQCStep(Project proj, Step plinkExportStep) {
    return register(GwasQCStep.create(proj, plinkExportStep));
  }

  Step generateAncestryStep(Project proj, final Step gwasQCStep) {
    return register(AncestryStep.create(proj, gwasQCStep));
  }

  Step generateFurtherAnalysisQCStep(Project proj, Step plinkExportStep, Step gwasQCStep,
                                     Step ancestryStep) {
    return register(FurtherAnalysisQCStep.create(proj, plinkExportStep, gwasQCStep, ancestryStep));
  }

  Step generateMosaicArmsStep(Project proj, final Step parseSamplesStep) {
    return register(MosaicArmsStep.create(proj, parseSamplesStep, numThreadsReq));
  }

  Step generateAnnotateSampleDataStep(Project proj, final Step sampleQCStep,
                                      final Step createSampleDataStep, final Step gwasQCStep) {
    return register(AnnotateSampleDataStep.create(proj, sampleQCStep, createSampleDataStep,
                                                  gwasQCStep));
  }

  Step generateMitoCNEstimateStep(Project proj, Step transposeStep) {
    return register(MitoCNEstimateStep.create(proj, transposeStep, numThreadsReq));
  }

  Step generatePFBStep(Project proj, final Step parseSamplesStep) {
    return register(ComputePFBStep.create(proj, parseSamplesStep));
  }

  Step generateSexCentroidsStep() {
    return register(SexCentroidsStep.create(numThreadsReq));
  }

  Step generateCNVStep(Project proj, Step pfbStep, Step gcModelStep) {
    return register(CallCNVsStep.create(proj, pfbStep, gcModelStep, numThreadsReq));
  }

  Step generatePCCorrectedProjectStep(Project proj, final Step parseSamplesStep) {
    return register(PCCorrectionStep.create(parseSamplesStep, numThreadsReq));
  }

  Step generateABLookupStep(final Step parseSamplesStep) {
    return register(ABLookupStep.create(parseSamplesStep));
  }

  public static int resolveThreads(Project proj, String arg) {
    int numThreads = Requirement.checkIntArgOrNeg1(arg);
    if (numThreads <= 0) {
      numThreads = proj.NUM_THREADS.getValue();
    }
    return numThreads;
  }

}
