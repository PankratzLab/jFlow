package org.genvisis.cnv.workflow;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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

  private Map<Step, Double> priorityMap;
  private List<Step> buildSteps;
  private double p;
  private final Requirement<Integer> numThreadsReq;
  public static final String NUM_THREADS_DESC = "Number of threads";
  public static final String PUTATIVE_WHITE_FILE_DESCRIPTION = "File with FID/IID pairs of putative white samples";

  public StepBuilder(Project proj) {
    numThreadsReq = new Requirement.PosIntRequirement(NUM_THREADS_DESC,
                                                      proj.NUM_THREADS.getValue());
    buildSteps = new ArrayList<>();
    priorityMap = new HashMap<>();
    p = 0.0;
  }

  public Requirement<Integer> getNumThreadsReq() {
    return numThreadsReq;
  }

  /**
   * @return All steps {@link #register(Step)}ed by this step builder thus far
   */
  public List<Step> getSortedSteps() {
    buildSteps.sort(new Comparator<Step>() {

      @Override
      public int compare(Step o1, Step o2) {
        return Double.compare(priorityMap.get(o1), priorityMap.get(o2));
      }
    });
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
    priorityMap.put(s, priority());
    buildSteps.add(s);
    return s;
  }

  Step generateIlluminaMarkerPositionsStep(Project proj) {
    return register(IlluminaMarkerPositionsStep.create(proj));
  }

  Step generateIlluminaMarkerBlastAnnotationStep(Project proj, final Step parseSamplesStep) {
    return register(IlluminaMarkerBlastStep.create(proj, parseSamplesStep, numThreadsReq));
  }

  Step generateAffyMarkerBlastAnnotationStep(final Project proj, final Step parseSamplesStep) {
    return register(AffyMarkerBlastStep.create(proj, parseSamplesStep, numThreadsReq));
  }

  Step generateParseSamplesStep(Project proj) {
    return generateParseSamplesStep(proj, null);
  }

  Step generateParseSamplesStep(Project proj, final Step markerPositionsStep) {
    return register(ParseSamplesStep.create(proj, markerPositionsStep, numThreadsReq));
  }

  Step generateCreateSampleDataStep(Project proj, final Step parseSamplesStep) {
    return register(SampleDataStep.create(proj, parseSamplesStep));
  }

  Step generateTransposeStep(Project proj, final Step parseSamplesStep) {
    return register(TransposeStep.create(proj, parseSamplesStep));
  }

  Step generateGCModelStep(Project proj) {
    return register(GCModelStep.create(proj));
  }

  Step generateSampleQCStep(Project proj, final Step parseSamplesStep) {
    return register(SampleQCStep.create(proj, parseSamplesStep, numThreadsReq));
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

  Step generateSexCentroidsStep(Project proj, final Step pfbStep) {
    return register(SexCentroidsStep.create(proj, pfbStep, numThreadsReq));
  }

  Step generateCNVStep(Project proj, Step pfbStep, Step gcModelStep) {
    return register(CallCNVsStep.create(proj, pfbStep, gcModelStep, numThreadsReq));
  }

  Step generatePCCorrectedProjectStep(Project proj, final Step parseSamplesStep) {
    return register(PCCorrectionStep.create(proj, parseSamplesStep, numThreadsReq));
  }

  Step generateABLookupStep(Project proj, final Step parseSamplesStep) {
    return register(ABLookupStep.create(proj, parseSamplesStep));
  }

  public static int resolveThreads(Project proj, int numThreads) {
    if (numThreads <= 0) {
      numThreads = proj.NUM_THREADS.getValue();
    }
    return numThreads;
  }

}