package org.genvisis.cnv.workflow;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.genvisis.cnv.Launch;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.gui.GenvisisWorkflowGUI;
import org.genvisis.cnv.hmm.CNVCaller;
import org.genvisis.cnv.qc.IlluminaMarkerBlast;
import org.genvisis.cnv.workflow.steps.AffyCELProcessingStep;
import org.genvisis.cnv.workflow.steps.AffyMarkerBlastStep;
import org.genvisis.cnv.workflow.steps.AffymetrixManifestParsingStep;
import org.genvisis.cnv.workflow.steps.AncestryStep;
import org.genvisis.cnv.workflow.steps.AxiomCELProcessingStep;
import org.genvisis.cnv.workflow.steps.ComputePFBStep;
import org.genvisis.cnv.workflow.steps.FurtherAnalysisQCStep;
import org.genvisis.cnv.workflow.steps.GCModelStep;
import org.genvisis.cnv.workflow.steps.GwasQCStep;
import org.genvisis.cnv.workflow.steps.IlluminaMarkerBlastStep;
import org.genvisis.cnv.workflow.steps.IlluminaMarkerPositionsStep;
import org.genvisis.cnv.workflow.steps.MarkerQCStep;
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
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.qsub.Qsub;
import org.pankratzlab.utils.gwas.Qc;
import org.pankratzlab.utils.gwas.QcMetric;
import org.pankratzlab.utils.gwas.RelationAncestryQc;

public class GenvisisWorkflow {

  private static final String NUM_THREADS_DESC = "Number of Threads to Use";
  public static final String PROJ_PROP_UPDATE_STR = " " + Project.class.getName() + " "
                                                    + CLI.ARG_PROJ + "=";
  public static final String PLINK_SUBDIR = "plink/";
  public static final String PLINKROOT = "plink";
  final Project proj;
  private final List<Step> steps;
  Logger log;
  private final Launch launch;

  public GenvisisWorkflow(Project project, Launch launch) {
    if (project == null) {
      throw new IllegalArgumentException(this.getClass().getName()
                                         + " cannot be constructed with a null "
                                         + Project.class.getName());
    }
    proj = project;
    log = project.getLog();
    this.launch = launch;

    steps = Collections.unmodifiableList(generateSteps(project.IS_PC_CORRECTED_PROJECT.getValue()));
  }

  public void showDialogAndRun() {
    GenvisisWorkflowGUI gui;
    gui = new GenvisisWorkflowGUI(proj, launch, steps);
    if (!gui.getCancelled()) {
      gui.setModal(true);
      gui.setVisible(true);

      if (gui.getCancelled()) {
        return;
      }
    }
  }

  private List<Step> generateSteps(boolean isPCCorrectedProject) {
    switch (proj.getArrayType()) {
      case AFFY_AXIOM:
        return generateAxiomSteps(isPCCorrectedProject);
      case AFFY_GW6:
      case AFFY_GW6_CN:
        return generateAffySteps(isPCCorrectedProject);
      case ILLUMINA:
        return isPCCorrectedProject ? generatePCCorrectedIlluminaSteps() : generateIlluminaSteps();
      case NGS_WGS:
        return generateNGSWGSSteps(isPCCorrectedProject);
      case NGS_WES:
      default:
        throw new UnsupportedOperationException("GenvisisWorkflow does not currently support arrays of type "
                                                + proj.getArrayType() + ".");
    }
  }

  private List<Step> generateNGSWGSSteps(boolean isPCCorrectedProject) {
    StepBuilder sb = new StepBuilder(proj);

    Step mdRAFParsingStep = sb.generateMosdepthImportStep();

    ReverseTransposeTarget reverseTransposeStep = sb.generateReverseTransposeStep(mdRAFParsingStep);
    SampleDataStep createSampleDataStep = sb.generateCreateSampleDataStep(reverseTransposeStep);
    GCModelStep gcModelStep = sb.generateGCModelStep();
    SampleQCStep sampleQCStep = sb.generateSampleQCStep(reverseTransposeStep);
    sb.generateSampleQCAnnotationStep(sampleQCStep);
    MarkerQCStep markerQCStep = sb.generateMarkerQCStep(reverseTransposeStep);
    SexChecksStep sexChecksStep = sb.generateSexChecksStep(null, createSampleDataStep,
                                                           reverseTransposeStep, sampleQCStep);
    sb.generateABLookupStep(reverseTransposeStep);
    PlinkExportStep plinkExportStep = sb.generatePlinkExportStep(reverseTransposeStep);
    GwasQCStep gwasQCStep = sb.generateGwasQCStep(plinkExportStep);
    AncestryStep ancestryStep = sb.generateAncestryStep(gwasQCStep);
    FurtherAnalysisQCStep faqcStep = sb.generateFurtherAnalysisQCStep(plinkExportStep, gwasQCStep,
                                                                      ancestryStep);
    sb.generateIdentifyProblemMarkersStep(markerQCStep, faqcStep);
    sb.generateMosaicArmsStep(reverseTransposeStep);
    sb.generateAnnotateSampleDataStep(sampleQCStep, createSampleDataStep, gwasQCStep);
    sb.generateMitoCNEstimateStep(mdRAFParsingStep);
    ComputePFBStep pfbStep = sb.generatePFBStep(reverseTransposeStep);
    sb.generateSexCentroidsStep(pfbStep);
    sb.generateCNVStep(pfbStep, gcModelStep);
    if (!isPCCorrectedProject) {
      sb.generatePCCorrectedProjectStep(reverseTransposeStep, sexChecksStep);
    }

    return sb.getSortedSteps();
  }

  private List<Step> generateAxiomSteps(boolean isPCCorrectedProject) {
    StepBuilder sb = new StepBuilder(proj);
    AffymetrixManifestParsingStep manifestStep;
    AxiomCELProcessingStep parseAxiomCELs;
    AffyMarkerBlastStep affyMarkerBlastStep = null;
    ReverseTransposeTarget reverseTransposeStep = null;
    SampleDataStep createSampleDataStep;
    GCModelStep gcModelStep;
    SampleQCStep sampleQCStep;
    manifestStep = sb.generateAffymetrixManifestParsingStep();
    parseAxiomCELs = sb.generateAxiomCELProcessingStep(manifestStep);
    reverseTransposeStep = sb.generateReverseTransposeStep(parseAxiomCELs);
    affyMarkerBlastStep = sb.generateAffyMarkerBlastAnnotationStep(reverseTransposeStep);
    createSampleDataStep = sb.generateCreateSampleDataStep(reverseTransposeStep);
    gcModelStep = sb.generateGCModelStep();
    sampleQCStep = sb.generateSampleQCStep(reverseTransposeStep);
    sb.generateSampleQCAnnotationStep(sampleQCStep);
    MarkerQCStep markerQCStep = sb.generateMarkerQCStep(reverseTransposeStep);
    SexChecksStep sexChecksStep = sb.generateSexChecksStep(affyMarkerBlastStep,
                                                           createSampleDataStep,
                                                           reverseTransposeStep, sampleQCStep);
    sb.generateABLookupStep(reverseTransposeStep);
    PlinkExportStep plinkExportStep = sb.generatePlinkExportStep(reverseTransposeStep);
    GwasQCStep gwasQCStep = sb.generateGwasQCStep(plinkExportStep);
    AncestryStep ancestryStep = sb.generateAncestryStep(gwasQCStep);
    FurtherAnalysisQCStep faqcStep = sb.generateFurtherAnalysisQCStep(plinkExportStep, gwasQCStep,
                                                                      ancestryStep);
    sb.generateIdentifyProblemMarkersStep(markerQCStep, faqcStep);
    sb.generateMosaicArmsStep(reverseTransposeStep);
    sb.generateAnnotateSampleDataStep(sampleQCStep, createSampleDataStep, gwasQCStep);
    sb.generateMitoCNEstimateStep(parseAxiomCELs);
    ComputePFBStep pfbStep = sb.generatePFBStep(reverseTransposeStep);
    sb.generateSexCentroidsStep(pfbStep);
    sb.generateCNVStep(pfbStep, gcModelStep);
    if (!isPCCorrectedProject) {
      sb.generatePCCorrectedProjectStep(reverseTransposeStep, sexChecksStep);
    }

    return sb.getSortedSteps();
  }

  private List<Step> generateAffySteps(boolean isPCCorrectedProject) {
    StepBuilder sb = new StepBuilder(proj);

    AffymetrixManifestParsingStep manifestStep;
    AffyCELProcessingStep parseAffyCELs;
    AffyMarkerBlastStep affyMarkerBlastStep = null;
    ReverseTransposeTarget reverseTransposeStep = null;
    SampleDataStep createSampleDataStep;
    GCModelStep gcModelStep;
    SampleQCStep sampleQCStep;
    manifestStep = sb.generateAffymetrixManifestParsingStep();
    parseAffyCELs = sb.generateAffyCELProcessingStep(manifestStep);
    reverseTransposeStep = sb.generateReverseTransposeStep(parseAffyCELs);
    affyMarkerBlastStep = sb.generateAffyMarkerBlastAnnotationStep(reverseTransposeStep);
    createSampleDataStep = sb.generateCreateSampleDataStep(reverseTransposeStep);
    gcModelStep = sb.generateGCModelStep();
    sampleQCStep = sb.generateSampleQCStep(reverseTransposeStep);
    sb.generateSampleQCAnnotationStep(sampleQCStep);
    MarkerQCStep markerQCStep = sb.generateMarkerQCStep(reverseTransposeStep);
    SexChecksStep sexChecksStep = sb.generateSexChecksStep(affyMarkerBlastStep,
                                                           createSampleDataStep,
                                                           reverseTransposeStep, sampleQCStep);
    sb.generateABLookupStep(reverseTransposeStep);
    PlinkExportStep plinkExportStep = sb.generatePlinkExportStep(reverseTransposeStep);
    GwasQCStep gwasQCStep = sb.generateGwasQCStep(plinkExportStep);
    AncestryStep ancestryStep = sb.generateAncestryStep(gwasQCStep);
    FurtherAnalysisQCStep faqcStep = sb.generateFurtherAnalysisQCStep(plinkExportStep, gwasQCStep,
                                                                      ancestryStep);
    sb.generateIdentifyProblemMarkersStep(markerQCStep, faqcStep);
    sb.generateMosaicArmsStep(reverseTransposeStep);
    sb.generateAnnotateSampleDataStep(sampleQCStep, createSampleDataStep, gwasQCStep);
    sb.generateMitoCNEstimateStep(parseAffyCELs);
    ComputePFBStep pfbStep = sb.generatePFBStep(reverseTransposeStep);
    sb.generateSexCentroidsStep(pfbStep);
    sb.generateCNVStep(pfbStep, gcModelStep);
    if (!isPCCorrectedProject) {
      sb.generatePCCorrectedProjectStep(reverseTransposeStep, sexChecksStep);
    }

    return sb.getSortedSteps();
  }

  private List<Step> generateIlluminaSteps() {
    StepBuilder sb = new StepBuilder(proj);
    ParseSamplesStep parseSamplesStep;
    IlluminaMarkerBlastStep illumMarkerBlastStep = null;
    SampleDataStep createSampleDataStep;
    GCModelStep gcModelStep;
    SampleQCStep sampleQCStep;
    TransposeStep transposeStep = null;
    IlluminaMarkerPositionsStep markerPositions = sb.generateIlluminaMarkerPositionsStep();
    parseSamplesStep = sb.generateParseSamplesStep(markerPositions);
    illumMarkerBlastStep = sb.generateIlluminaMarkerBlastAnnotationStep(parseSamplesStep);
    createSampleDataStep = sb.generateCreateSampleDataStep(parseSamplesStep);
    transposeStep = sb.generateTransposeStep(parseSamplesStep);
    gcModelStep = sb.generateGCModelStep();
    sampleQCStep = sb.generateSampleQCStep(parseSamplesStep);
    sb.generateSampleQCAnnotationStep(sampleQCStep);
    MarkerQCStep markerQCStep = sb.generateMarkerQCStep(parseSamplesStep);
    SexChecksStep sexChecksStep = sb.generateSexChecksStep(illumMarkerBlastStep,
                                                           createSampleDataStep, transposeStep,
                                                           sampleQCStep);
    sb.generateABLookupStep(parseSamplesStep);
    PlinkExportStep plinkExportStep = sb.generatePlinkExportStep(parseSamplesStep);
    GwasQCStep gwasQCStep = sb.generateGwasQCStep(plinkExportStep);
    AncestryStep ancestryStep = sb.generateAncestryStep(gwasQCStep);
    FurtherAnalysisQCStep faqcStep = sb.generateFurtherAnalysisQCStep(plinkExportStep, gwasQCStep,
                                                                      ancestryStep);
    sb.generateIdentifyProblemMarkersStep(markerQCStep, faqcStep);
    sb.generateMosaicArmsStep(parseSamplesStep);
    sb.generateAnnotateSampleDataStep(sampleQCStep, createSampleDataStep, gwasQCStep);
    sb.generateMitoCNEstimateStep(transposeStep);
    ComputePFBStep pfbStep = sb.generatePFBStep(parseSamplesStep);
    sb.generateSexCentroidsStep(pfbStep);
    sb.generateCNVStep(pfbStep, gcModelStep);
    sb.generatePCCorrectedProjectStep(parseSamplesStep, sexChecksStep);

    return sb.getSortedSteps();
  }

  private List<Step> generatePCCorrectedIlluminaSteps() {
    StepBuilder sb = new StepBuilder(proj);
    ReverseTransposeTarget reverseTransposeStep = null;
    IlluminaMarkerBlastStep illumMarkerBlastStep = null;
    SampleDataStep createSampleDataStep;
    GCModelStep gcModelStep;
    SampleQCStep sampleQCStep;
    reverseTransposeStep = sb.generateReverseTransposeStep(Step.EMPTY_STEP);
    illumMarkerBlastStep = sb.generateIlluminaMarkerBlastAnnotationStep(reverseTransposeStep);
    createSampleDataStep = sb.generateCreateSampleDataStep(reverseTransposeStep);
    gcModelStep = sb.generateGCModelStep();
    sampleQCStep = sb.generateSampleQCStep(reverseTransposeStep);
    sb.generateSampleQCAnnotationStep(sampleQCStep);
    MarkerQCStep markerQCStep = sb.generateMarkerQCStep(reverseTransposeStep);
    sb.generateSexChecksStep(illumMarkerBlastStep, createSampleDataStep, reverseTransposeStep,
                             sampleQCStep);
    sb.generateABLookupStep(reverseTransposeStep);
    PlinkExportStep plinkExportStep = sb.generatePlinkExportStep(reverseTransposeStep);
    GwasQCStep gwasQCStep = sb.generateGwasQCStep(plinkExportStep);
    AncestryStep ancestryStep = sb.generateAncestryStep(gwasQCStep);
    FurtherAnalysisQCStep faqcStep = sb.generateFurtherAnalysisQCStep(plinkExportStep, gwasQCStep,
                                                                      ancestryStep);
    sb.generateIdentifyProblemMarkersStep(markerQCStep, faqcStep);
    sb.generateMosaicArmsStep(reverseTransposeStep);
    sb.generateAnnotateSampleDataStep(sampleQCStep, createSampleDataStep, gwasQCStep);
    sb.generateMitoCNEstimateStep(reverseTransposeStep);
    ComputePFBStep pfbStep = sb.generatePFBStep(reverseTransposeStep);
    sb.generateSexCentroidsStep(pfbStep);
    sb.generateCNVStep(pfbStep, gcModelStep);

    return sb.getSortedSteps();
  }

  public static String getLocationOfSampleMap(Project proj) {
    String filename;

    String projDir = proj.PROJECT_DIRECTORY.getValue();
    String snpMap = "Sample_Map.csv";
    String snpMapGz = "Sample_Map.csv.gz";
    if (Files.exists(projDir + snpMap)) {
      filename = projDir + snpMap;
    } else if (Files.exists(projDir + snpMapGz)) {
      filename = projDir + snpMapGz;
    } else {
      String srcDir = proj.SOURCE_DIRECTORY.getValue();
      if (Files.exists(srcDir + snpMap)) {
        filename = srcDir + snpMap;
      } else if (Files.exists(srcDir + snpMapGz)) {
        filename = srcDir + snpMapGz;
      } else {
        return null;
      }
    }
    return filename;
  }

  public static String getPlinkDir(Project proj) {
    return ext.verifyDirFormat(proj.PROJECT_DIRECTORY.getValue() + PLINK_SUBDIR);
  }

  public static String getAncestryDir(Project proj) {
    return getPlinkDir(proj) + Qc.QC_SUBDIR + RelationAncestryQc.ANCESTRY_DIR;
  }

  public static void maybeSetProjNumThreads(Project proj, int numThreads) {
    if (numThreads != proj.NUM_THREADS.getValue()) {
      proj.NUM_THREADS.setValue(numThreads);
    }
  }

  public static String setupImputationDefaults(Project proj) {
    return setupIlluminaImputation(proj, Runtime.getRuntime().availableProcessors(), null, null,
                                   false, null, true);
  }

  public static String setupAffyImputation(Project proj, int numThreads, String putativeWhitesFile,
                                           Map<QcMetric, String> faqcThreshs, String aptExeDir,
                                           String aptLibDir, String sketch, boolean includeQC) {
    StepBuilder sb = new StepBuilder(proj);

    Requirement<Integer> numThreadsReq = sb.getNumThreadsReq();
    AffymetrixManifestParsingStep manifestStep = sb.generateAffymetrixManifestParsingStep();
    AffyCELProcessingStep parseCELFiles = sb.generateAffyCELProcessingStep(manifestStep);
    ReverseTransposeTarget reverseTranspose = sb.generateReverseTransposeStep(parseCELFiles);
    SampleDataStep sampleData = sb.generateCreateSampleDataStep(reverseTranspose);
    AffyMarkerBlastStep blast = sb.generateAffyMarkerBlastAnnotationStep(reverseTranspose);
    SampleQCStep sampleQc = null;
    SampleQCAnnotateStep sampleExcludes = null;
    MarkerQCStep markerQc = null;
    SexChecksStep sexChecks = null;
    PlinkExportStep exportPlink = null;
    GwasQCStep gwasQc = null;
    AncestryStep ancestry = null;
    FurtherAnalysisQCStep faqcStep = null;
    if (includeQC) {
      sampleQc = sb.generateSampleQCStep(reverseTranspose);
      sampleExcludes = sb.generateSampleQCAnnotationStep(sampleQc);
      markerQc = sb.generateMarkerQCStep(reverseTranspose);
      sexChecks = sb.generateSexChecksStep(blast, sampleData, reverseTranspose, sampleQc);
      exportPlink = sb.generatePlinkExportStep(reverseTranspose);
      gwasQc = sb.generateGwasQCStep(exportPlink);
      ancestry = sb.generateAncestryStep(gwasQc);
      faqcStep = sb.generateFurtherAnalysisQCStep(exportPlink, gwasQc, ancestry);
    }

    Variables stepReqs;
    Map<Step, Variables> varMap = new HashMap<>();

    stepReqs = parseCELFiles.getDefaultRequirementValues();
    for (Requirement<?> r1 : stepReqs.keys()) {
      if (r1.getDescription().equals(AffyCELProcessingStep.DESC_SKETCH)) {
        stepReqs.parseOrFail(r1, sketch);
      } else if (r1.getDescription().equals(AffyCELProcessingStep.DESC_APT_EXT)) {
        stepReqs.parseOrFail(r1, aptExeDir);
      } else if (r1.getDescription().equals(AffyCELProcessingStep.DESC_APT_LIB)) {
        stepReqs.parseOrFail(r1, aptLibDir);
      }
    }
    varMap.put(parseCELFiles, stepReqs);
    varMap.put(reverseTranspose, reverseTranspose.getDefaultRequirementValues());

    stepReqs = sampleData.getDefaultRequirementValues();
    for (Requirement<?> r1 : stepReqs.keys()) {
      if (r1.getDescription().equalsIgnoreCase(SampleDataStep.REQ_CREATE_MINIMAL)) {
        stepReqs.parseOrFail(r1, "false");
        break;
      }
    }
    varMap.put(sampleData, stepReqs);

    if (includeQC) {
      varMap.put(sampleQc, sampleQc.getDefaultRequirementValues());
      varMap.put(sampleExcludes, sampleExcludes.getDefaultRequirementValues());
      varMap.put(markerQc, markerQc.getDefaultRequirementValues());

      stepReqs = sexChecks.getDefaultRequirementValues();
      for (Requirement<?> r1 : stepReqs.keys()) {
        if (r1.getDescription().equals(SexChecksStep.NO_CROSS_HYBE_REQUIREMENT)) {
          stepReqs.parseOrFail(r1, "false");
        } else if (r1.getDescription().equals(SexChecksStep.ADD_ESTSEX_TO_SAMPDATA_REQUIREMENT)) {
          stepReqs.parseOrFail(r1, "true");
        }
      }

      varMap.put(sexChecks, stepReqs);
      varMap.put(exportPlink, exportPlink.getDefaultRequirementValues());

      stepReqs = gwasQc.getDefaultRequirementValues();
      if (faqcThreshs != null && !faqcThreshs.isEmpty()) {
        fixQCThreshs(stepReqs, faqcThreshs);
      }
      varMap.put(gwasQc, stepReqs);

      stepReqs = ancestry.getDefaultRequirementValues();
      if (putativeWhitesFile != null) {
        Requirement<?> r = null;
        for (Requirement<?> r1 : stepReqs.keys()) {
          if (r1.getDescription().equals(StepBuilder.PUTATIVE_WHITE_FILE_DESCRIPTION)) {
            r = r1;
            break;
          }
        }
        if (r == null) {
          throw new IllegalStateException();
        }
        stepReqs.parseOrFail(r, putativeWhitesFile);
      }
      varMap.put(ancestry, stepReqs);

      stepReqs = faqcStep.getDefaultRequirementValues();
      if (faqcThreshs != null && !faqcThreshs.isEmpty()) {
        fixQCThreshs(stepReqs, faqcThreshs);
      }
      varMap.put(faqcStep, stepReqs);
    }

    // override threads defaults
    if (numThreads > 0) {
      for (Step s : varMap.keySet()) {
        if (varMap.get(s).hasValid(numThreadsReq)) {
          varMap.get(s).put(numThreadsReq, numThreads);
        }
      }
    }

    StringBuilder output = new StringBuilder("## Genvisis Project Pipeline - Stepwise Commands\n\n");

    addStepInfo(output, parseCELFiles, parseCELFiles.getCommandLine(varMap.get(parseCELFiles)));
    addStepInfo(output, reverseTranspose,
                reverseTranspose.getCommandLine(varMap.get(reverseTranspose)));
    addStepInfo(output, sampleData, sampleData.getCommandLine(varMap.get(sampleData)));
    if (includeQC) {
      addStepInfo(output, sampleQc, sampleQc.getCommandLine(varMap.get(sampleQc)));
      addStepInfo(output, sampleExcludes,
                  sampleExcludes.getCommandLine(varMap.get(sampleExcludes)));
      addStepInfo(output, markerQc, markerQc.getCommandLine(varMap.get(markerQc)));
      addStepInfo(output, sexChecks, sexChecks.getCommandLine(varMap.get(sexChecks)));
      addStepInfo(output, exportPlink, exportPlink.getCommandLine(varMap.get(exportPlink)));
      addStepInfo(output, gwasQc, gwasQc.getCommandLine(varMap.get(gwasQc)));
      addStepInfo(output, ancestry, ancestry.getCommandLine(varMap.get(ancestry)));
      addStepInfo(output, faqcStep, faqcStep.getCommandLine(varMap.get(faqcStep)));
    }

    return output.toString();
  }

  public static String setupIlluminaImputation(Project proj, int numThreads,
                                               String putativeWhitesFile,
                                               Map<QcMetric, String> faqcThreshs,
                                               boolean parseSource, String man, boolean includeQC) {
    StepBuilder sb = new StepBuilder(proj);

    Requirement<Integer> numThreadsReq = sb.getNumThreadsReq();

    IlluminaMarkerPositionsStep createMkrPos = man != null ? sb.generateIlluminaMarkerPositionsStep()
                                                           : null;
    ParseSamplesStep parseSamples = sb.generateParseSamplesStep(createMkrPos);
    TransposeStep transpose = sb.generateTransposeStep(parseSamples);
    SampleDataStep sampleData = sb.generateCreateSampleDataStep(parseSamples);
    IlluminaMarkerBlastStep blast = sb.generateIlluminaMarkerBlastAnnotationStep(parseSamples);
    SampleQCStep sampleQc = null;
    SampleQCAnnotateStep sampleExcludes = null;
    MarkerQCStep markerQc = null;
    SexChecksStep sexChecks = null;
    PlinkExportStep exportPlink = null;
    GwasQCStep gwasQc = null;
    AncestryStep ancestry = null;
    FurtherAnalysisQCStep faqcStep = null;
    if (includeQC) {
      sampleQc = sb.generateSampleQCStep(parseSamples);
      sampleExcludes = sb.generateSampleQCAnnotationStep(sampleQc);
      markerQc = sb.generateMarkerQCStep(parseSamples);
      sexChecks = sb.generateSexChecksStep(blast, sampleData, transpose, sampleQc);
      exportPlink = sb.generatePlinkExportStep(parseSamples);
      gwasQc = sb.generateGwasQCStep(exportPlink);
      ancestry = sb.generateAncestryStep(gwasQc);
      faqcStep = sb.generateFurtherAnalysisQCStep(exportPlink, gwasQc, ancestry);
    }

    Variables stepReqs;
    Map<Step, Variables> varMap = new HashMap<>();

    if (createMkrPos != null) {
      stepReqs = createMkrPos.getDefaultRequirementValues();
      for (Requirement<?> r1 : stepReqs.keys()) {
        if (r1.getDescription().contains("Manifest")) {
          stepReqs.parseOrFail(r1, man);
        }
      }
      varMap.put(createMkrPos, stepReqs);
    }
    if (parseSource) {
      varMap.put(parseSamples, parseSamples.getDefaultRequirementValues());
      varMap.put(transpose, transpose.getDefaultRequirementValues());
      stepReqs = blast.getDefaultRequirementValues();
      for (Requirement<?> r1 : stepReqs.keys()) {
        if (r1.getDescription().contains(IlluminaMarkerBlast.DESC_MANIFEST.substring(1))) {
          stepReqs.parseOrFail(r1, man);
          break;
        }
      }
      varMap.put(blast, stepReqs);

      stepReqs = sampleData.getDefaultRequirementValues();
      for (Requirement<?> r1 : stepReqs.keys()) {
        if (r1.getDescription().equalsIgnoreCase(SampleDataStep.REQ_CREATE_MINIMAL)) {
          stepReqs.parseOrFail(r1, "false");
          break;
        }
      }
      varMap.put(sampleData, stepReqs);
    }

    if (includeQC) {
      varMap.put(sampleQc, sampleQc.getDefaultRequirementValues());
      varMap.put(sampleExcludes, sampleExcludes.getDefaultRequirementValues());
      varMap.put(markerQc, markerQc.getDefaultRequirementValues());

      stepReqs = sexChecks.getDefaultRequirementValues();
      for (Requirement<?> r1 : stepReqs.keys()) {
        if (r1.getDescription().equals(SexChecksStep.NO_CROSS_HYBE_REQUIREMENT)) {
          stepReqs.parseOrFail(r1, "false");
        } else if (r1.getDescription().equals(SexChecksStep.ADD_ESTSEX_TO_SAMPDATA_REQUIREMENT)) {
          stepReqs.parseOrFail(r1, "true");
        }
      }

      varMap.put(sexChecks, stepReqs);
      varMap.put(exportPlink, exportPlink.getDefaultRequirementValues());

      stepReqs = gwasQc.getDefaultRequirementValues();
      if (faqcThreshs != null && !faqcThreshs.isEmpty()) {
        fixQCThreshs(stepReqs, faqcThreshs);
      }
      varMap.put(gwasQc, stepReqs);

      stepReqs = ancestry.getDefaultRequirementValues();
      if (putativeWhitesFile != null) {
        Requirement<?> r = null;
        for (Requirement<?> r1 : stepReqs.keys()) {
          if (r1.getDescription().equals(StepBuilder.PUTATIVE_WHITE_FILE_DESCRIPTION)) {
            r = r1;
            break;
          }
        }
        if (r == null) {
          throw new IllegalStateException();
        }
        stepReqs.parseOrFail(r, putativeWhitesFile);
      }
      varMap.put(ancestry, stepReqs);

      stepReqs = faqcStep.getDefaultRequirementValues();
      if (faqcThreshs != null && !faqcThreshs.isEmpty()) {
        fixQCThreshs(stepReqs, faqcThreshs);
      }
      varMap.put(faqcStep, stepReqs);
    }

    // override threads defaults
    if (numThreads > 0) {
      for (Step s : varMap.keySet()) {
        if (varMap.get(s).hasValid(numThreadsReq)) {
          varMap.get(s).put(numThreadsReq, numThreads);
        }
      }
    }

    StringBuilder output = new StringBuilder("## Genvisis Project Pipeline - Stepwise Commands\n\n");

    if (createMkrPos != null) {
      addStepInfo(output, createMkrPos, createMkrPos.getCommandLine(varMap.get(createMkrPos)));
    }
    if (parseSource) {
      addStepInfo(output, parseSamples, parseSamples.getCommandLine(varMap.get(parseSamples)));
      addStepInfo(output, transpose, transpose.getCommandLine(varMap.get(transpose)));
      addStepInfo(output, sampleData, sampleData.getCommandLine(varMap.get(sampleData)));
      addStepInfo(output, blast, blast.getCommandLine(varMap.get(blast)));
    }
    if (includeQC) {
      addStepInfo(output, sampleQc, sampleQc.getCommandLine(varMap.get(sampleQc)));
      addStepInfo(output, sampleExcludes,
                  sampleExcludes.getCommandLine(varMap.get(sampleExcludes)));
      addStepInfo(output, markerQc, markerQc.getCommandLine(varMap.get(markerQc)));
      addStepInfo(output, sexChecks, sexChecks.getCommandLine(varMap.get(sexChecks)));
      addStepInfo(output, exportPlink, exportPlink.getCommandLine(varMap.get(exportPlink)));
      addStepInfo(output, gwasQc, gwasQc.getCommandLine(varMap.get(gwasQc)));
      addStepInfo(output, ancestry, ancestry.getCommandLine(varMap.get(ancestry)));
      addStepInfo(output, faqcStep, faqcStep.getCommandLine(varMap.get(faqcStep)));
    }

    return output.toString();
  }

  private static void fixQCThreshs(Variables reqMap, Map<QcMetric, String> threshMap) {
    Map<String, String> qc = new HashMap<>();
    for (QcMetric met : threshMap.keySet()) {
      qc.put(met.getUserDescription(), threshMap.get(met));
    }
    for (Requirement<?> r1 : reqMap.keys()) {
      if (qc.containsKey(r1.getDescription())) {
        reqMap.parseOrFail(r1, qc.get(r1.getDescription()));
      }
    }
  }

  public static void setupCNVCalling(String projectProperties) {
    Project pcProj = new Project(projectProperties);
    StepBuilder sb = new StepBuilder(pcProj);
    Step transpose = sb.generateTransposeStep(null);
    // Create new sample data, run sex checks?
    GCModelStep gc = sb.generateGCModelStep();
    ComputePFBStep pfb = sb.generatePFBStep((ParseSamplesStep) null);
    SexCentroidsStep cent = sb.generateSexCentroidsStep(pfb);
    Step cnv = sb.generateCNVStep(pfb, gc);
    Variables cnvOpts = new Variables();
    List<Requirement<?>> reqs = cnv.getRequirements().getFlatRequirementsList();
    for (Requirement<?> req : reqs) {
      if (req.getDescription().equals(CNVCaller.CNV_SCOPE_DESC)) {
        cnvOpts.parseOrFail(req, CNVCaller.CALLING_SCOPE.BOTH.toString());
      } else if (req.getDescription().equals(NUM_THREADS_DESC)) {
        cnvOpts.parseOrFail(req, "" + (Runtime.getRuntime().availableProcessors() - 1));
      }
    }

    String s1 = transpose.getCommandLine(null);
    String s2 = gc.getCommandLine(null);
    String s3 = pfb.getCommandLine(null);
    String s4 = cent.getCommandLine(null);
    String s5 = cnv.getCommandLine(cnvOpts);

    String file = pcProj.PROJECT_DIRECTORY.getValue() + "CNVCallingPipeline.";
    String suggFile = file + ext.getTimestampForFilename() + ".pbs";
    String runFile = file + ext.getTimestampForFilename() + ".run";

    StringBuilder output = new StringBuilder("## Genvisis Project Pipeline - Stepwise Commands\n\n");

    addStepInfo(output, transpose, s1);
    addStepInfo(output, gc, s2);
    addStepInfo(output, pfb, s3);
    addStepInfo(output, cent, s4);
    addStepInfo(output, cnv, s5);

    Qsub.qsubDefaults(suggFile, output.toString());
    Files.write(output.toString(), runFile);
  }

  public static void addStepInfo(StringBuilder output, Step step, String stepCmd) {
    output.append("## ").append(step.getName()).append("\n");
    output.append("echo \">>>> start ").append(step.getName()).append(" at: \" `date`")
          .append("\n");
    output.append(stepCmd).append("\n");
    output.append("echo \"<<<< end ").append(step.getName()).append(" at: \" `date`").append("\n");
    output.append("\n\n");
  }

}
