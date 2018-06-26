package org.genvisis.cnv.workflow;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.genvisis.CLI;
import org.genvisis.cnv.Launch;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.gui.GenvisisWorkflowGUI;
import org.genvisis.cnv.hmm.CNVCaller;
import org.genvisis.cnv.qc.IlluminaMarkerBlast;
import org.genvisis.cnv.workflow.steps.AffyCELProcessingStep;
import org.genvisis.cnv.workflow.steps.AffyMarkerBlastStep;
import org.genvisis.cnv.workflow.steps.AncestryStep;
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
import org.genvisis.cnv.workflow.steps.SampleQCStep;
import org.genvisis.cnv.workflow.steps.SexCentroidsStep;
import org.genvisis.cnv.workflow.steps.SexChecksStep;
import org.genvisis.cnv.workflow.steps.TransposeStep;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.gwas.MarkerQC.QC_METRIC;
import org.genvisis.gwas.Qc;
import org.genvisis.gwas.RelationAncestryQc;
import org.genvisis.qsub.Qsub;

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

    steps = Collections.unmodifiableList(generateSteps(!project.IS_PC_CORRECTED_PROJECT.getValue()));
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

  private List<Step> generateSteps(boolean allowCorrectionStep) {
    boolean isAffy = proj.getArrayType() == ARRAY.AFFY_GW6
                     || proj.getArrayType() == ARRAY.AFFY_GW6_CN
                     || proj.getArrayType() == ARRAY.AFFY_AXIOM;
    if (isAffy) {
      return generateAffySteps(allowCorrectionStep);
    } else {
      return generateIlluminaSteps(allowCorrectionStep);
    }
  }

  private List<Step> generateAffySteps(boolean allowCorrectionStep) {
    StepBuilder sb = new StepBuilder(proj);
    AffyCELProcessingStep parseAffyCELs;
    AffyMarkerBlastStep affyMarkerBlastStep = null;
    ReverseTransposeTarget reverseTransposeStep = null;
    SampleDataStep createSampleDataStep;
    GCModelStep gcModelStep;
    SampleQCStep sampleQCStep;
    parseAffyCELs = sb.generateAffyCELProcessingStep(proj);
    reverseTransposeStep = sb.generateReverseTransposeStep(proj, parseAffyCELs);
    affyMarkerBlastStep = sb.generateAffyMarkerBlastAnnotationStep(proj, reverseTransposeStep);
    createSampleDataStep = sb.generateCreateSampleDataStep(proj, reverseTransposeStep);
    gcModelStep = sb.generateGCModelStep(proj);
    sampleQCStep = sb.generateSampleQCStep(proj, reverseTransposeStep);
    sb.generateMarkerQCStep(proj, reverseTransposeStep);
    sb.generateSexChecksStep(proj, reverseTransposeStep, affyMarkerBlastStep, createSampleDataStep,
                             sampleQCStep);
    sb.generateABLookupStep(proj, reverseTransposeStep);
    PlinkExportStep plinkExportStep = sb.generatePlinkExportStep(proj, reverseTransposeStep);
    GwasQCStep gwasQCStep = sb.generateGwasQCStep(proj, plinkExportStep);
    AncestryStep ancestryStep = sb.generateAncestryStep(proj, gwasQCStep);
    sb.generateFurtherAnalysisQCStep(proj, plinkExportStep, gwasQCStep, ancestryStep);
    sb.generateMosaicArmsStep(proj, reverseTransposeStep);
    sb.generateAnnotateSampleDataStep(proj, sampleQCStep, createSampleDataStep, gwasQCStep);
    sb.generateMitoCNEstimateStep(proj, reverseTransposeStep);
    ComputePFBStep pfbStep = sb.generatePFBStep(proj, reverseTransposeStep);
    sb.generateSexCentroidsStep(proj, pfbStep);
    sb.generateCNVStep(proj, pfbStep, gcModelStep);
    if (allowCorrectionStep) {
      sb.generatePCCorrectedProjectStep(proj, reverseTransposeStep);
    }

    return sb.getSortedSteps();
  }

  private List<Step> generateIlluminaSteps(boolean allowCorrectionStep) {
    StepBuilder sb = new StepBuilder(proj);
    ParseSamplesStep parseSamplesStep;
    IlluminaMarkerBlastStep illumMarkerBlastStep = null;
    SampleDataStep createSampleDataStep;
    GCModelStep gcModelStep;
    SampleQCStep sampleQCStep;
    TransposeStep transposeStep = null;
    IlluminaMarkerPositionsStep markerPositions = sb.generateIlluminaMarkerPositionsStep(proj);
    parseSamplesStep = sb.generateParseSamplesStep(proj, markerPositions);
    illumMarkerBlastStep = sb.generateIlluminaMarkerBlastAnnotationStep(proj, parseSamplesStep);
    createSampleDataStep = sb.generateCreateSampleDataStep(proj, parseSamplesStep);
    transposeStep = sb.generateTransposeStep(proj, parseSamplesStep);
    gcModelStep = sb.generateGCModelStep(proj);
    sampleQCStep = sb.generateSampleQCStep(proj, parseSamplesStep);
    sb.generateMarkerQCStep(proj, parseSamplesStep);
    sb.generateSexChecksStep(proj, parseSamplesStep, illumMarkerBlastStep, createSampleDataStep,
                             transposeStep, sampleQCStep);
    sb.generateABLookupStep(proj, parseSamplesStep);
    PlinkExportStep plinkExportStep = sb.generatePlinkExportStep(proj, parseSamplesStep);
    GwasQCStep gwasQCStep = sb.generateGwasQCStep(proj, plinkExportStep);
    AncestryStep ancestryStep = sb.generateAncestryStep(proj, gwasQCStep);
    sb.generateFurtherAnalysisQCStep(proj, plinkExportStep, gwasQCStep, ancestryStep);
    sb.generateMosaicArmsStep(proj, parseSamplesStep);
    sb.generateAnnotateSampleDataStep(proj, sampleQCStep, createSampleDataStep, gwasQCStep);
    sb.generateMitoCNEstimateStep(proj, transposeStep);
    ComputePFBStep pfbStep = sb.generatePFBStep(proj, parseSamplesStep);
    sb.generateSexCentroidsStep(proj, pfbStep);
    sb.generateCNVStep(proj, pfbStep, gcModelStep);
    if (allowCorrectionStep) {
      sb.generatePCCorrectedProjectStep(proj, parseSamplesStep);
    }

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
                                           Map<QC_METRIC, String> faqcThreshs, String aptExeDir,
                                           String aptLibDir, String sketch, boolean includeQC) {
    StepBuilder sb = new StepBuilder(proj);

    Requirement<Integer> numThreadsReq = sb.getNumThreadsReq();
    AffyCELProcessingStep parseCELFiles = sb.generateAffyCELProcessingStep(proj);
    ReverseTransposeTarget reverseTranspose = sb.generateReverseTransposeStep(proj, parseCELFiles);
    SampleDataStep sampleData = sb.generateCreateSampleDataStep(proj, reverseTranspose);
    AffyMarkerBlastStep blast = sb.generateAffyMarkerBlastAnnotationStep(proj, reverseTranspose);
    SampleQCStep sampleQc = null;
    MarkerQCStep markerQc = null;
    SexChecksStep sexChecks = null;
    PlinkExportStep exportPlink = null;
    GwasQCStep gwasQc = null;
    AncestryStep ancestry = null;
    FurtherAnalysisQCStep faqcStep = null;
    if (includeQC) {
      sampleQc = sb.generateSampleQCStep(proj, reverseTranspose);
      markerQc = sb.generateMarkerQCStep(proj, reverseTranspose);
      sexChecks = sb.generateSexChecksStep(proj, reverseTranspose, blast, sampleData, sampleQc);
      exportPlink = sb.generatePlinkExportStep(proj, reverseTranspose);
      gwasQc = sb.generateGwasQCStep(proj, exportPlink);
      ancestry = sb.generateAncestryStep(proj, gwasQc);
      faqcStep = sb.generateFurtherAnalysisQCStep(proj, exportPlink, gwasQc, ancestry);
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
                                               Map<QC_METRIC, String> faqcThreshs,
                                               boolean parseSource, String man, boolean includeQC) {
    StepBuilder sb = new StepBuilder(proj);

    Requirement<Integer> numThreadsReq = sb.getNumThreadsReq();

    IlluminaMarkerPositionsStep createMkrPos = man != null ? sb.generateIlluminaMarkerPositionsStep(proj)
                                                           : null;
    ParseSamplesStep parseSamples = sb.generateParseSamplesStep(proj, createMkrPos);
    TransposeStep transpose = sb.generateTransposeStep(proj, parseSamples);
    SampleDataStep sampleData = sb.generateCreateSampleDataStep(proj, parseSamples);
    IlluminaMarkerBlastStep blast = sb.generateIlluminaMarkerBlastAnnotationStep(proj,
                                                                                 parseSamples);
    SampleQCStep sampleQc = null;
    MarkerQCStep markerQc = null;
    SexChecksStep sexChecks = null;
    PlinkExportStep exportPlink = null;
    GwasQCStep gwasQc = null;
    AncestryStep ancestry = null;
    FurtherAnalysisQCStep faqcStep = null;
    if (includeQC) {
      sampleQc = sb.generateSampleQCStep(proj, parseSamples);
      markerQc = sb.generateMarkerQCStep(proj, parseSamples);
      sexChecks = sb.generateSexChecksStep(proj, parseSamples, blast, sampleData, transpose,
                                           sampleQc);
      exportPlink = sb.generatePlinkExportStep(proj, parseSamples);
      gwasQc = sb.generateGwasQCStep(proj, exportPlink);
      ancestry = sb.generateAncestryStep(proj, gwasQc);
      faqcStep = sb.generateFurtherAnalysisQCStep(proj, exportPlink, gwasQc, ancestry);
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
        if (r1.getDescription().equalsIgnoreCase(IlluminaMarkerBlast.DESC_MANIFEST)) {
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
      addStepInfo(output, markerQc, markerQc.getCommandLine(varMap.get(markerQc)));
      addStepInfo(output, sexChecks, sexChecks.getCommandLine(varMap.get(sexChecks)));
      addStepInfo(output, exportPlink, exportPlink.getCommandLine(varMap.get(exportPlink)));
      addStepInfo(output, gwasQc, gwasQc.getCommandLine(varMap.get(gwasQc)));
      addStepInfo(output, ancestry, ancestry.getCommandLine(varMap.get(ancestry)));
      addStepInfo(output, faqcStep, faqcStep.getCommandLine(varMap.get(faqcStep)));
    }

    return output.toString();
  }

  private static void fixQCThreshs(Variables reqMap, Map<QC_METRIC, String> threshMap) {
    Map<String, String> qc = new HashMap<>();
    for (QC_METRIC met : threshMap.keySet()) {
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
    Step transpose = sb.generateTransposeStep(pcProj, null);
    // Create new sample data, run sex checks?
    GCModelStep gc = sb.generateGCModelStep(pcProj);
    ComputePFBStep pfb = sb.generatePFBStep(pcProj, (ParseSamplesStep) null);
    SexCentroidsStep cent = sb.generateSexCentroidsStep(pcProj, pfb);
    Step cnv = sb.generateCNVStep(pcProj, pfb, gc);
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
