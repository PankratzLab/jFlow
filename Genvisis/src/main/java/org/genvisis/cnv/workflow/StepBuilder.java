package org.genvisis.cnv.workflow;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.workflow.steps.ABLookupStep;
import org.genvisis.cnv.workflow.steps.AffyCELProcessingStep;
import org.genvisis.cnv.workflow.steps.AffyMarkerBlastStep;
import org.genvisis.cnv.workflow.steps.AncestryStep;
import org.genvisis.cnv.workflow.steps.AnnotateSampleDataStep;
import org.genvisis.cnv.workflow.steps.AxiomCELProcessingStep;
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
import org.genvisis.cnv.workflow.steps.ReverseTransposeTarget;
import org.genvisis.cnv.workflow.steps.SampleDataStep;
import org.genvisis.cnv.workflow.steps.SampleQCStep;
import org.genvisis.cnv.workflow.steps.SexCentroidsStep;
import org.genvisis.cnv.workflow.steps.SexChecksStep;
import org.genvisis.cnv.workflow.steps.TransposeStep;
import org.pankratzlab.common.CLI;

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

  private Project proj;
  private Map<Step, Double> priorityMap;
  private Map<Class<? extends Step>, Step> stepInstanceMap;
  private List<Step> buildSteps;
  private double p;
  private final Requirement<Integer> numThreadsReq;
  public static final String NUM_THREADS_DESC = "Number of threads";
  public static final String PUTATIVE_WHITE_FILE_DESCRIPTION = "File with FID/IID pairs of putative white samples";

  public StepBuilder(Project proj) {
    this.proj = proj;
    numThreadsReq = createNumThreadsReq(proj);
    buildSteps = new ArrayList<>();
    priorityMap = new HashMap<>();
    stepInstanceMap = new HashMap<>();
    p = 0.0;
  }

  public static Requirement<Integer> createNumThreadsReq(Project proj) {
    return new Requirement.PosIntRequirement(CLI.ARG_THREADS, NUM_THREADS_DESC,
                                             proj.NUM_THREADS.getValue());
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
    stepInstanceMap.put(s.getClass(), s);
    priorityMap.put(s, priority());
    buildSteps.add(s);
    return s;
  }

  IlluminaMarkerPositionsStep generateIlluminaMarkerPositionsStep() {
    return stepInstanceMap.containsKey(IlluminaMarkerPositionsStep.class) ? (IlluminaMarkerPositionsStep) stepInstanceMap.get(IlluminaMarkerPositionsStep.class)
                                                                          : register(IlluminaMarkerPositionsStep.create(proj));
  }

  IlluminaMarkerBlastStep generateIlluminaMarkerBlastAnnotationStep(Step parseSamplesStep) {
    return register(IlluminaMarkerBlastStep.create(proj, parseSamplesStep, numThreadsReq));
  }

  AxiomCELProcessingStep generateAxiomCELProcessingStep() {
    return register(AxiomCELProcessingStep.create(proj, numThreadsReq));
  }

  AffyCELProcessingStep generateAffyCELProcessingStep() {
    return register(AffyCELProcessingStep.create(proj, numThreadsReq));
  }

  AffyMarkerBlastStep generateAffyMarkerBlastAnnotationStep(ReverseTransposeTarget parseSamplesStep) {
    return register(AffyMarkerBlastStep.create(proj, parseSamplesStep, numThreadsReq));
  }

  ParseSamplesStep generateParseSamplesStep() {
    return generateParseSamplesStep(null);
  }

  ParseSamplesStep generateParseSamplesStep(IlluminaMarkerPositionsStep markerPositionsStep) {
    return register(ParseSamplesStep.create(proj, markerPositionsStep, numThreadsReq));
  }

  SampleDataStep generateCreateSampleDataStep(ReverseTransposeTarget reverseTransposeTarget) {
    return register(SampleDataStep.create(proj, reverseTransposeTarget));
  }

  SampleDataStep generateCreateSampleDataStep(ParseSamplesStep parseSamplesStep) {
    return register(SampleDataStep.create(proj, parseSamplesStep));
  }

  ReverseTransposeTarget generateReverseTransposeStep(Step parseAffyCELs) {
    return register(ReverseTransposeTarget.create(proj, parseAffyCELs));
  }

  TransposeStep generateTransposeStep(Step parseSamplesStep) {
    return register(TransposeStep.create(proj, parseSamplesStep));
  }

  GCModelStep generateGCModelStep() {
    return register(GCModelStep.create(proj));
  }

  SampleQCStep generateSampleQCStep(Step parseSamplesStep) {
    return register(SampleQCStep.create(proj, parseSamplesStep, numThreadsReq));
  }

  MarkerQCStep generateMarkerQCStep(Step parseSamplesStep) {
    return register(MarkerQCStep.create(proj, parseSamplesStep, numThreadsReq));
  }

  public SexChecksStep generateSexChecksStep(Step markerBlastStep, SampleDataStep sampleDataStep,
                                             Step transposeStep, SampleQCStep sampleQCStep) {
    return register(SexChecksStep.create(proj, markerBlastStep, sampleDataStep, transposeStep,
                                         sampleQCStep));
  }

  PlinkExportStep generatePlinkExportStep(ParseSamplesStep parseSamplesStep) {
    return register(PlinkExportStep.create(proj, parseSamplesStep));
  }

  PlinkExportStep generatePlinkExportStep(ReverseTransposeTarget reverseTransposeStep) {
    return register(PlinkExportStep.create(proj, reverseTransposeStep));
  }

  GwasQCStep generateGwasQCStep(PlinkExportStep plinkExportStep) {
    return register(GwasQCStep.create(proj, plinkExportStep));
  }

  AncestryStep generateAncestryStep(GwasQCStep gwasQCStep) {
    return register(AncestryStep.create(proj, gwasQCStep));
  }

  FurtherAnalysisQCStep generateFurtherAnalysisQCStep(PlinkExportStep plinkExportStep,
                                                      GwasQCStep gwasQCStep,
                                                      AncestryStep ancestryStep) {
    return register(FurtherAnalysisQCStep.create(proj, plinkExportStep, gwasQCStep, ancestryStep));
  }

  MosaicArmsStep generateMosaicArmsStep(ParseSamplesStep parseSamplesStep) {
    return register(MosaicArmsStep.create(proj, parseSamplesStep, numThreadsReq));
  }

  MosaicArmsStep generateMosaicArmsStep(ReverseTransposeTarget reverseTransposeStep) {
    return register(MosaicArmsStep.create(proj, reverseTransposeStep, numThreadsReq));
  }

  AnnotateSampleDataStep generateAnnotateSampleDataStep(SampleQCStep sampleQCStep,
                                                        SampleDataStep createSampleDataStep,
                                                        GwasQCStep gwasQCStep) {
    return register(AnnotateSampleDataStep.create(proj, sampleQCStep, createSampleDataStep,
                                                  gwasQCStep));
  }

  MitoCNEstimateStep generateMitoCNEstimateStep(TransposeStep transposeStep) {
    return register(MitoCNEstimateStep.create(proj, transposeStep, numThreadsReq));
  }

  MitoCNEstimateStep generateMitoCNEstimateStep(ReverseTransposeTarget reverseTransposeStep) {
    return register(MitoCNEstimateStep.create(proj, reverseTransposeStep, numThreadsReq));
  }

  ComputePFBStep generatePFBStep(ParseSamplesStep parseSamplesStep) {
    return register(ComputePFBStep.create(proj, parseSamplesStep));
  }

  ComputePFBStep generatePFBStep(ReverseTransposeTarget reverseTransposeStep) {
    return register(ComputePFBStep.create(proj, reverseTransposeStep));
  }

  SexCentroidsStep generateSexCentroidsStep(ComputePFBStep pfbStep) {
    return register(SexCentroidsStep.create(proj, pfbStep, numThreadsReq));
  }

  CallCNVsStep generateCNVStep(Step pfbStep, GCModelStep gcModelStep) {
    return register(CallCNVsStep.create(proj, pfbStep, gcModelStep, numThreadsReq));
  }

  public PCCorrectionStep generatePCCorrectedProjectStep(Step parseSamplesStep) {
    return register(PCCorrectionStep.create(proj, parseSamplesStep, numThreadsReq));
  }

  ABLookupStep generateABLookupStep(ParseSamplesStep parseSamplesStep) {
    return register(ABLookupStep.create(proj, parseSamplesStep));
  }

  ABLookupStep generateABLookupStep(ReverseTransposeTarget reverseTransposeStep) {
    return register(ABLookupStep.create(proj, reverseTransposeStep));
  }

  public Step generateSamplesParsingStep() {
    switch (proj.getArrayType()) {
      case AFFY_AXIOM:
        AxiomCELProcessingStep parseAxiomCELs = generateAxiomCELProcessingStep();
        return generateReverseTransposeStep(parseAxiomCELs);
      case AFFY_GW6:
      case AFFY_GW6_CN:
        Step parseAffyCELs = generateAffyCELProcessingStep();
        return generateReverseTransposeStep(parseAffyCELs);
      case ILLUMINA:
        IlluminaMarkerPositionsStep markerPositions = generateIlluminaMarkerPositionsStep();
        return generateParseSamplesStep(markerPositions);
      case NGS:
      default:
        throw new UnsupportedOperationException("GenvisisWorkflow does not currently support arrays of type "
                                                + proj.getArrayType() + ".");
    }
  }

  public Step generateMarkerBlastStep() {
    switch (proj.getArrayType()) {
      case AFFY_AXIOM:
        AxiomCELProcessingStep parseAxiomCELs = generateAxiomCELProcessingStep();
        ReverseTransposeTarget reverseTranspose = generateReverseTransposeStep(parseAxiomCELs);
        return generateAffyMarkerBlastAnnotationStep(reverseTranspose);
      case AFFY_GW6:
      case AFFY_GW6_CN:
        Step parseAffyCELs = generateAffyCELProcessingStep();
        ReverseTransposeTarget reverseTranspose1 = generateReverseTransposeStep(parseAffyCELs);
        return generateAffyMarkerBlastAnnotationStep(reverseTranspose1);
      case ILLUMINA:
        IlluminaMarkerPositionsStep markerPositions = generateIlluminaMarkerPositionsStep();
        ParseSamplesStep parseSamplesStep = generateParseSamplesStep(markerPositions);
        return generateIlluminaMarkerBlastAnnotationStep(parseSamplesStep);
      case NGS:
      default:
        throw new UnsupportedOperationException("GenvisisWorkflow does not currently support arrays of type "
                                                + proj.getArrayType() + ".");
    }
  }

  public static int resolveThreads(Project proj, int numThreads) {
    if (numThreads <= 0) {
      numThreads = proj.NUM_THREADS.getValue();
    }
    return numThreads;
  }

}
