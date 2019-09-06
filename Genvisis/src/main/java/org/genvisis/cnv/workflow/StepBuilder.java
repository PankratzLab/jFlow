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
import org.genvisis.cnv.workflow.steps.AxiomManifestParsingStep;
import org.genvisis.cnv.workflow.steps.CallCNVsStep;
import org.genvisis.cnv.workflow.steps.ComputePFBStep;
import org.genvisis.cnv.workflow.steps.FurtherAnalysisQCStep;
import org.genvisis.cnv.workflow.steps.GCModelStep;
import org.genvisis.cnv.workflow.steps.GwasQCStep;
import org.genvisis.cnv.workflow.steps.IdentifyProblemMarkersStep;
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
import org.genvisis.cnv.workflow.steps.SampleQCAnnotateStep;
import org.genvisis.cnv.workflow.steps.SampleQCStep;
import org.genvisis.cnv.workflow.steps.SexCentroidsStep;
import org.genvisis.cnv.workflow.steps.SexChecksStep;
import org.genvisis.cnv.workflow.steps.TransposeStep;
import org.pankratzlab.common.CLI;

import com.google.common.collect.ClassToInstanceMap;
import com.google.common.collect.MutableClassToInstanceMap;

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
  private ClassToInstanceMap<Step> stepInstanceMap;
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
    stepInstanceMap = MutableClassToInstanceMap.create();
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

  public AxiomManifestParsingStep generateAxiomManifestParsingStep() {
    return stepInstanceMap.containsKey(AxiomManifestParsingStep.class) ? stepInstanceMap.getInstance(AxiomManifestParsingStep.class)
                                                                       : register(AxiomManifestParsingStep.create(proj));
  }

  IlluminaMarkerPositionsStep generateIlluminaMarkerPositionsStep() {
    return stepInstanceMap.containsKey(IlluminaMarkerPositionsStep.class) ? stepInstanceMap.getInstance(IlluminaMarkerPositionsStep.class)
                                                                          : register(IlluminaMarkerPositionsStep.create(proj));
  }

  IlluminaMarkerBlastStep generateIlluminaMarkerBlastAnnotationStep(Step parseSamplesStep) {
    return stepInstanceMap.containsKey(IlluminaMarkerBlastStep.class) ? stepInstanceMap.getInstance(IlluminaMarkerBlastStep.class)
                                                                      : register(IlluminaMarkerBlastStep.create(proj,
                                                                                                                parseSamplesStep,
                                                                                                                numThreadsReq));
  }

  public AxiomCELProcessingStep generateAxiomCELProcessingStep(AxiomManifestParsingStep axiomManifestParsingStep) {
    return stepInstanceMap.containsKey(AxiomCELProcessingStep.class) ? stepInstanceMap.getInstance(AxiomCELProcessingStep.class)
                                                                     : register(AxiomCELProcessingStep.create(proj,
                                                                                                              axiomManifestParsingStep,
                                                                                                              numThreadsReq));
  }

  public AffyCELProcessingStep generateAffyCELProcessingStep() {
    return stepInstanceMap.containsKey(AffyCELProcessingStep.class) ? stepInstanceMap.getInstance(AffyCELProcessingStep.class)
                                                                    : register(AffyCELProcessingStep.create(proj,
                                                                                                            numThreadsReq));
  }

  AffyMarkerBlastStep generateAffyMarkerBlastAnnotationStep(ReverseTransposeTarget parseSamplesStep) {
    return stepInstanceMap.containsKey(AffyMarkerBlastStep.class) ? stepInstanceMap.getInstance(AffyMarkerBlastStep.class)
                                                                  : register(AffyMarkerBlastStep.create(proj,
                                                                                                        parseSamplesStep,
                                                                                                        numThreadsReq));
  }

  ParseSamplesStep generateParseSamplesStep(IlluminaMarkerPositionsStep markerPositionsStep) {
    return stepInstanceMap.containsKey(ParseSamplesStep.class) ? stepInstanceMap.getInstance(ParseSamplesStep.class)
                                                               : register(ParseSamplesStep.create(proj,
                                                                                                  markerPositionsStep,
                                                                                                  numThreadsReq));
  }

  public SampleDataStep generateCreateSampleDataStep(Step samplesParsingStep) {
    return stepInstanceMap.containsKey(SampleDataStep.class) ? stepInstanceMap.getInstance(SampleDataStep.class)
                                                             : register(SampleDataStep.create(proj,
                                                                                              samplesParsingStep));
  }

  public ReverseTransposeTarget generateReverseTransposeStep(Step parseAffyCELs) {
    return stepInstanceMap.containsKey(ReverseTransposeTarget.class) ? stepInstanceMap.getInstance(ReverseTransposeTarget.class)
                                                                     : register(ReverseTransposeTarget.create(proj,
                                                                                                              parseAffyCELs));
  }

  public TransposeStep generateTransposeStep(Step parseSamplesStep) {
    return stepInstanceMap.containsKey(TransposeStep.class) ? stepInstanceMap.getInstance(TransposeStep.class)
                                                            : register(TransposeStep.create(proj,
                                                                                            parseSamplesStep));
  }

  GCModelStep generateGCModelStep() {
    return stepInstanceMap.containsKey(GCModelStep.class) ? stepInstanceMap.getInstance(GCModelStep.class)
                                                          : register(GCModelStep.create(proj));
  }

  public SampleQCStep generateSampleQCStep(Step parseSamplesStep) {
    return stepInstanceMap.containsKey(SampleQCStep.class) ? stepInstanceMap.getInstance(SampleQCStep.class)
                                                           : register(SampleQCStep.create(proj,
                                                                                          parseSamplesStep,
                                                                                          numThreadsReq));
  }

  public SampleQCAnnotateStep generateSampleQCAnnotationStep(Step sampleQCStep) {
    return stepInstanceMap.containsKey(SampleQCAnnotateStep.class) ? stepInstanceMap.getInstance(SampleQCAnnotateStep.class)
                                                                   : register(SampleQCAnnotateStep.create(proj,
                                                                                                          sampleQCStep));
  }

  public MarkerQCStep generateMarkerQCStep(Step parseSamplesStep) {
    return stepInstanceMap.containsKey(MarkerQCStep.class) ? stepInstanceMap.getInstance(MarkerQCStep.class)
                                                           : register(MarkerQCStep.create(proj,
                                                                                          parseSamplesStep,
                                                                                          numThreadsReq));
  }

  public SexChecksStep generateSexChecksStep(Step markerBlastStep, SampleDataStep sampleDataStep,
                                             Step transposeStep, SampleQCStep sampleQCStep) {
    return stepInstanceMap.containsKey(SexChecksStep.class) ? stepInstanceMap.getInstance(SexChecksStep.class)
                                                            : register(SexChecksStep.create(proj,
                                                                                            markerBlastStep,
                                                                                            sampleDataStep,
                                                                                            transposeStep,
                                                                                            sampleQCStep));
  }

  public PlinkExportStep generatePlinkExportStep(Step parseSamplesStep) {
    return stepInstanceMap.containsKey(PlinkExportStep.class) ? stepInstanceMap.getInstance(PlinkExportStep.class)
                                                              : register(PlinkExportStep.create(proj,
                                                                                                parseSamplesStep));
  }

  public GwasQCStep generateGwasQCStep(PlinkExportStep plinkExportStep) {
    return stepInstanceMap.containsKey(GwasQCStep.class) ? stepInstanceMap.getInstance(GwasQCStep.class)
                                                         : register(GwasQCStep.create(proj,
                                                                                      plinkExportStep));
  }

  public AncestryStep generateAncestryStep(GwasQCStep gwasQCStep) {
    return stepInstanceMap.containsKey(AncestryStep.class) ? stepInstanceMap.getInstance(AncestryStep.class)
                                                           : register(AncestryStep.create(proj,
                                                                                          gwasQCStep));
  }

  public FurtherAnalysisQCStep generateFurtherAnalysisQCStep(PlinkExportStep plinkExportStep,
                                                             GwasQCStep gwasQCStep,
                                                             AncestryStep ancestryStep) {
    return stepInstanceMap.containsKey(FurtherAnalysisQCStep.class) ? stepInstanceMap.getInstance(FurtherAnalysisQCStep.class)
                                                                    : register(FurtherAnalysisQCStep.create(proj,
                                                                                                            plinkExportStep,
                                                                                                            gwasQCStep,
                                                                                                            ancestryStep));
  }

  IdentifyProblemMarkersStep generateIdentifyProblemMarkersStep(MarkerQCStep markerQCStep,
                                                                FurtherAnalysisQCStep faqcStep) {
    return stepInstanceMap.containsKey(IdentifyProblemMarkersStep.class) ? stepInstanceMap.getInstance(IdentifyProblemMarkersStep.class)
                                                                         : register(IdentifyProblemMarkersStep.create(proj,
                                                                                                                      markerQCStep,
                                                                                                                      faqcStep));
  }

  MosaicArmsStep generateMosaicArmsStep(ParseSamplesStep parseSamplesStep) {
    return stepInstanceMap.containsKey(MosaicArmsStep.class) ? stepInstanceMap.getInstance(MosaicArmsStep.class)
                                                             : register(MosaicArmsStep.create(proj,
                                                                                              parseSamplesStep,
                                                                                              numThreadsReq));
  }

  MosaicArmsStep generateMosaicArmsStep(ReverseTransposeTarget reverseTransposeStep) {
    return stepInstanceMap.containsKey(MosaicArmsStep.class) ? stepInstanceMap.getInstance(MosaicArmsStep.class)
                                                             : register(MosaicArmsStep.create(proj,
                                                                                              reverseTransposeStep,
                                                                                              numThreadsReq));
  }

  AnnotateSampleDataStep generateAnnotateSampleDataStep(SampleQCStep sampleQCStep,
                                                        SampleDataStep createSampleDataStep,
                                                        GwasQCStep gwasQCStep) {
    return stepInstanceMap.containsKey(AnnotateSampleDataStep.class) ? stepInstanceMap.getInstance(AnnotateSampleDataStep.class)
                                                                     : register(AnnotateSampleDataStep.create(proj,
                                                                                                              sampleQCStep,
                                                                                                              createSampleDataStep,
                                                                                                              gwasQCStep));
  }

  public MitoCNEstimateStep generateMitoCNEstimateStep(Step markersParsingStep) {
    return stepInstanceMap.containsKey(MitoCNEstimateStep.class) ? stepInstanceMap.getInstance(MitoCNEstimateStep.class)
                                                                 : register(MitoCNEstimateStep.create(proj,
                                                                                                      markersParsingStep,
                                                                                                      numThreadsReq));
  }

  ComputePFBStep generatePFBStep(ParseSamplesStep parseSamplesStep) {
    return stepInstanceMap.containsKey(ComputePFBStep.class) ? stepInstanceMap.getInstance(ComputePFBStep.class)
                                                             : register(ComputePFBStep.create(proj,
                                                                                              parseSamplesStep));
  }

  ComputePFBStep generatePFBStep(ReverseTransposeTarget reverseTransposeStep) {
    return stepInstanceMap.containsKey(ComputePFBStep.class) ? stepInstanceMap.getInstance(ComputePFBStep.class)
                                                             : register(ComputePFBStep.create(proj,
                                                                                              reverseTransposeStep));
  }

  SexCentroidsStep generateSexCentroidsStep(ComputePFBStep pfbStep) {
    return stepInstanceMap.containsKey(SexCentroidsStep.class) ? stepInstanceMap.getInstance(SexCentroidsStep.class)
                                                               : register(SexCentroidsStep.create(proj,
                                                                                                  pfbStep,
                                                                                                  numThreadsReq));
  }

  CallCNVsStep generateCNVStep(Step pfbStep, GCModelStep gcModelStep) {
    return stepInstanceMap.containsKey(CallCNVsStep.class) ? stepInstanceMap.getInstance(CallCNVsStep.class)
                                                           : register(CallCNVsStep.create(proj,
                                                                                          pfbStep,
                                                                                          gcModelStep,
                                                                                          numThreadsReq));
  }

  public PCCorrectionStep generatePCCorrectedProjectStep(Step parseSamplesStep,
                                                         SexChecksStep sexChecksStep) {
    return stepInstanceMap.containsKey(PCCorrectionStep.class) ? stepInstanceMap.getInstance(PCCorrectionStep.class)
                                                               : register(PCCorrectionStep.create(proj,
                                                                                                  parseSamplesStep,
                                                                                                  sexChecksStep,
                                                                                                  numThreadsReq));
  }

  ABLookupStep generateABLookupStep(ParseSamplesStep parseSamplesStep) {
    return stepInstanceMap.containsKey(ABLookupStep.class) ? stepInstanceMap.getInstance(ABLookupStep.class)
                                                           : register(ABLookupStep.create(proj,
                                                                                          parseSamplesStep));
  }

  ABLookupStep generateABLookupStep(ReverseTransposeTarget reverseTransposeStep) {
    return stepInstanceMap.containsKey(ABLookupStep.class) ? stepInstanceMap.getInstance(ABLookupStep.class)
                                                           : register(ABLookupStep.create(proj,
                                                                                          reverseTransposeStep));
  }

  public Step generateMarkersParsingStep() {
    switch (proj.getArrayType()) {
      case AFFY_AXIOM:
        AxiomManifestParsingStep manifestStep = generateAxiomManifestParsingStep();
        return generateAxiomCELProcessingStep(manifestStep);
      case AFFY_GW6:
      case AFFY_GW6_CN:
        return generateAffyCELProcessingStep();
      case ILLUMINA:
        IlluminaMarkerPositionsStep markerPositions = generateIlluminaMarkerPositionsStep();
        ParseSamplesStep parseSamplesStep = generateParseSamplesStep(markerPositions);
        return generateTransposeStep(parseSamplesStep);
      case NGS:
      default:
        throw new UnsupportedOperationException("GenvisisWorkflow does not currently support arrays of type "
                                                + proj.getArrayType() + ".");
    }

  }

  public Step generateSamplesParsingStep() {
    switch (proj.getArrayType()) {
      case AFFY_AXIOM:
        AxiomManifestParsingStep manifestStep = generateAxiomManifestParsingStep();
        AxiomCELProcessingStep parseAxiomCELs = generateAxiomCELProcessingStep(manifestStep);
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
        AxiomManifestParsingStep manifestStep = generateAxiomManifestParsingStep();
        AxiomCELProcessingStep parseAxiomCELs = generateAxiomCELProcessingStep(manifestStep);
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
