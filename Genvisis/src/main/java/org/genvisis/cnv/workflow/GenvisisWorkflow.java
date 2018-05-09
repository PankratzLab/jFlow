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
    StepBuilder sb = new StepBuilder(proj);

    ParseSamplesStep parseSamplesStep;
    IlluminaMarkerBlastStep illumMarkerBlastStep = null;
    AffyMarkerBlastStep affyMarkerBlastStep = null;
    if (proj.getArrayType() == ARRAY.AFFY_GW6 || proj.getArrayType() == ARRAY.AFFY_GW6_CN
        || proj.getArrayType() == ARRAY.AFFY_AXIOM) {
      parseSamplesStep = sb.generateParseSamplesStep(proj);
      affyMarkerBlastStep = sb.generateAffyMarkerBlastAnnotationStep(proj, parseSamplesStep);
    } else {
      IlluminaMarkerPositionsStep markerPositions = sb.generateIlluminaMarkerPositionsStep(proj);
      parseSamplesStep = sb.generateParseSamplesStep(proj, markerPositions);
      illumMarkerBlastStep = sb.generateIlluminaMarkerBlastAnnotationStep(proj, parseSamplesStep);
    }
    SampleDataStep createSampleDataStep = sb.generateCreateSampleDataStep(proj, parseSamplesStep);
    TransposeStep transposeStep = sb.generateTransposeStep(proj, parseSamplesStep);
    GCModelStep gcModelStep = sb.generateGCModelStep(proj);
    SampleQCStep sampleQCStep = sb.generateSampleQCStep(proj, parseSamplesStep);
    sb.generateMarkerQCStep(proj, parseSamplesStep);
    if (proj.getArrayType() == ARRAY.AFFY_GW6 || proj.getArrayType() == ARRAY.AFFY_GW6_CN
        || proj.getArrayType() == ARRAY.AFFY_AXIOM) {
      sb.generateSexChecksStep(proj, parseSamplesStep, affyMarkerBlastStep, createSampleDataStep,
                               transposeStep, sampleQCStep);
    } else {
      sb.generateSexChecksStep(proj, parseSamplesStep, illumMarkerBlastStep, createSampleDataStep,
                               transposeStep, sampleQCStep);
    }
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
    return setupImputation(proj, Runtime.getRuntime().availableProcessors(), null, null, false,
                           null);
  }

  public static String setupImputation(Project proj, int numThreads, String putativeWhitesFile,
                                       Map<QC_METRIC, String> faqcThreshs, boolean parseSource,
                                       String man) {
    StepBuilder sb = new StepBuilder(proj);

    Requirement<Integer> numThreadsReq = sb.getNumThreadsReq();
    IlluminaMarkerPositionsStep createMkrPos = man != null ? sb.generateIlluminaMarkerPositionsStep(proj)
                                                           : null;
    ParseSamplesStep parseSamples = sb.generateParseSamplesStep(proj, createMkrPos);
    TransposeStep transpose = sb.generateTransposeStep(proj, parseSamples);
    SampleDataStep sampleData = sb.generateCreateSampleDataStep(proj, parseSamples);
    IlluminaMarkerBlastStep blast = sb.generateIlluminaMarkerBlastAnnotationStep(proj,
                                                                                 parseSamples);
    SampleQCStep sampleQc = sb.generateSampleQCStep(proj, parseSamples);
    MarkerQCStep markerQc = sb.generateMarkerQCStep(proj, parseSamples);
    SexChecksStep sexChecks = sb.generateSexChecksStep(proj, parseSamples, blast, sampleData,
                                                       transpose, sampleQc);
    PlinkExportStep exportPlink = sb.generatePlinkExportStep(proj, parseSamples);
    GwasQCStep gwasQc = sb.generateGwasQCStep(proj, exportPlink);
    AncestryStep ancestry = sb.generateAncestryStep(proj, gwasQc);
    FurtherAnalysisQCStep faqcStep = sb.generateFurtherAnalysisQCStep(proj, exportPlink, gwasQc,
                                                                      ancestry);

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
    varMap.put(sampleQc, sampleQc.getDefaultRequirementValues());
    varMap.put(markerQc, markerQc.getDefaultRequirementValues());

    stepReqs = sexChecks.getDefaultRequirementValues();
    for (Requirement<?> r1 : stepReqs.keys()) {
      if (r1.getDescription().equals(SexChecksStep.NO_CROSS_HYBE_REQUIREMENT)) {
        stepReqs.parseOrFail(r1, "true");
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

    // override threads defaults
    if (numThreads > 0) {
      for (Step s : varMap.keySet()) {
        if (varMap.get(s).has(numThreadsReq)) {
          varMap.get(s).put(numThreadsReq, numThreads);
        }
      }
    }

    String s0 = createMkrPos == null ? "" : createMkrPos.getCommandLine(varMap.get(createMkrPos));
    String s1 = parseSamples.getCommandLine(varMap.get(parseSamples));
    String s2 = transpose.getCommandLine(varMap.get(transpose));
    String s3 = sampleData.getCommandLine(varMap.get(sampleData));
    String s4 = blast.getCommandLine(varMap.get(blast));
    String s5 = sampleQc.getCommandLine(varMap.get(sampleQc));
    String s6 = markerQc.getCommandLine(varMap.get(markerQc));
    String s7 = sexChecks.getCommandLine(varMap.get(sexChecks));
    String s8 = exportPlink.getCommandLine(varMap.get(exportPlink));
    String s9 = gwasQc.getCommandLine(varMap.get(gwasQc));
    String s10 = ancestry.getCommandLine(varMap.get(ancestry));
    String s11 = faqcStep.getCommandLine(varMap.get(faqcStep));

    StringBuilder output = new StringBuilder("## Genvisis Project Pipeline - Stepwise Commands\n\n");

    if (createMkrPos != null) {
      addStepInfo(output, createMkrPos, s0);
    }
    if (parseSource) {
      addStepInfo(output, parseSamples, s1);
      addStepInfo(output, transpose, s2);
      addStepInfo(output, sampleData, s3);
      addStepInfo(output, blast, s4);
    }
    addStepInfo(output, sampleQc, s5);
    addStepInfo(output, markerQc, s6);
    addStepInfo(output, sexChecks, s7);
    addStepInfo(output, exportPlink, s8);
    addStepInfo(output, gwasQc, s9);
    addStepInfo(output, ancestry, s10);
    addStepInfo(output, faqcStep, s11);

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
    ComputePFBStep pfb = sb.generatePFBStep(pcProj, null);
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
