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
  <T extends Step> T register(T s) {
    priorityMap.put(s, priority());
    buildSteps.add(s);
    return s;
  }

  IlluminaMarkerPositionsStep generateIlluminaMarkerPositionsStep(Project proj) {
    return register(IlluminaMarkerPositionsStep.create(proj));
  }

  IlluminaMarkerBlastStep generateIlluminaMarkerBlastAnnotationStep(Project proj,
                                                                    ParseSamplesStep parseSamplesStep) {
    return register(IlluminaMarkerBlastStep.create(proj, parseSamplesStep, numThreadsReq));
  }

  AffyMarkerBlastStep generateAffyMarkerBlastAnnotationStep(final Project proj,
                                                            ParseSamplesStep parseSamplesStep) {
    return register(AffyMarkerBlastStep.create(proj, parseSamplesStep, numThreadsReq));
  }

  ParseSamplesStep generateParseSamplesStep(Project proj) {
    return generateParseSamplesStep(proj, null);
  }

  ParseSamplesStep generateParseSamplesStep(Project proj,
                                            IlluminaMarkerPositionsStep markerPositionsStep) {
    return register(ParseSamplesStep.create(proj, markerPositionsStep, numThreadsReq));
  }

  SampleDataStep generateCreateSampleDataStep(Project proj, ParseSamplesStep parseSamplesStep) {
    return register(SampleDataStep.create(proj, parseSamplesStep));
  }

  TransposeStep generateTransposeStep(Project proj, ParseSamplesStep parseSamplesStep) {
    return register(TransposeStep.create(proj, parseSamplesStep));
  }

  GCModelStep generateGCModelStep(Project proj) {
    return register(GCModelStep.create(proj));
  }

  SampleQCStep generateSampleQCStep(Project proj, ParseSamplesStep parseSamplesStep) {
    return register(SampleQCStep.create(proj, parseSamplesStep, numThreadsReq));
  }

  MarkerQCStep generateMarkerQCStep(Project proj, ParseSamplesStep parseSamplesStep) {
    return register(MarkerQCStep.create(proj, parseSamplesStep, numThreadsReq));
  }

  SexChecksStep generateSexChecksStep(Project proj, ParseSamplesStep parseSamplesStep,
                                      IlluminaMarkerBlastStep markerBlastStep,
                                      SampleDataStep sampleDataStep, TransposeStep transposeStep,
                                      SampleQCStep sampleQCStep) {
    return register(SexChecksStep.create(proj, parseSamplesStep, markerBlastStep, sampleDataStep,
                                         transposeStep, sampleQCStep));
  }

  SexChecksStep generateSexChecksStep(Project proj, ParseSamplesStep parseSamplesStep,
                                      AffyMarkerBlastStep markerBlastStep,
                                      SampleDataStep sampleDataStep, TransposeStep transposeStep,
                                      SampleQCStep sampleQCStep) {
    return register(SexChecksStep.create(proj, parseSamplesStep, markerBlastStep, sampleDataStep,
                                         transposeStep, sampleQCStep));
  }

  PlinkExportStep generatePlinkExportStep(Project proj, ParseSamplesStep parseSamplesStep) {
    return register(PlinkExportStep.create(proj, parseSamplesStep));
  }

  GwasQCStep generateGwasQCStep(Project proj, PlinkExportStep plinkExportStep) {
    return register(GwasQCStep.create(proj, plinkExportStep));
  }

  AncestryStep generateAncestryStep(Project proj, GwasQCStep gwasQCStep) {
    return register(AncestryStep.create(proj, gwasQCStep));
  }

  FurtherAnalysisQCStep generateFurtherAnalysisQCStep(Project proj, PlinkExportStep plinkExportStep,
                                                      GwasQCStep gwasQCStep,
                                                      AncestryStep ancestryStep) {
    return register(FurtherAnalysisQCStep.create(proj, plinkExportStep, gwasQCStep, ancestryStep));
  }

  MosaicArmsStep generateMosaicArmsStep(Project proj, ParseSamplesStep parseSamplesStep) {
    return register(MosaicArmsStep.create(proj, parseSamplesStep, numThreadsReq));
  }

  AnnotateSampleDataStep generateAnnotateSampleDataStep(Project proj, SampleQCStep sampleQCStep,
                                                        SampleDataStep createSampleDataStep,
                                                        GwasQCStep gwasQCStep) {
    return register(AnnotateSampleDataStep.create(proj, sampleQCStep, createSampleDataStep,
                                                  gwasQCStep));
  }

  MitoCNEstimateStep generateMitoCNEstimateStep(Project proj, TransposeStep transposeStep) {
    return register(MitoCNEstimateStep.create(proj, transposeStep, numThreadsReq));
  }

  ComputePFBStep generatePFBStep(Project proj, ParseSamplesStep parseSamplesStep) {
    return register(ComputePFBStep.create(proj, parseSamplesStep));
  }

  SexCentroidsStep generateSexCentroidsStep(Project proj, ComputePFBStep pfbStep) {
    return register(SexCentroidsStep.create(proj, pfbStep, numThreadsReq));
  }

  CallCNVsStep generateCNVStep(Project proj, Step pfbStep, GCModelStep gcModelStep) {
    return register(CallCNVsStep.create(proj, pfbStep, gcModelStep, numThreadsReq));
  }

  PCCorrectionStep generatePCCorrectedProjectStep(Project proj, ParseSamplesStep parseSamplesStep) {
    return register(PCCorrectionStep.create(proj, parseSamplesStep, numThreadsReq));
  }

  ABLookupStep generateABLookupStep(Project proj, ParseSamplesStep parseSamplesStep) {
    return register(ABLookupStep.create(proj, parseSamplesStep));
  }

  public static int resolveThreads(Project proj, int numThreads) {
    if (numThreads <= 0) {
      numThreads = proj.NUM_THREADS.getValue();
    }
    return numThreads;
  }

}
