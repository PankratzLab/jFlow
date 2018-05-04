package org.genvisis.cnv.workflow;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import org.genvisis.cnv.Launch;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.gui.GenvisisWorkflowGUI;
import org.genvisis.cnv.hmm.CNVCaller;<<<<<<<Upstream,based on origin/master
import org.genvisis.cnv.hmm.CNVCaller.CALLING_SCOPE;
import org.genvisis.cnv.hmm.CNVCaller.PFB_MANAGEMENT_TYPE;
import org.genvisis.cnv.hmm.PFB;
import org.genvisis.cnv.manage.MitoPipeline;
import org.genvisis.cnv.manage.PRoCtOR;
import org.genvisis.cnv.manage.PlinkData;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.cnv.prop.Property;
import org.genvisis.cnv.qc.AffyMarkerBlast;
import org.genvisis.cnv.qc.GcAdjustor;=======>>>>>>>7d 6642e  Begin extracting steps into classes
import org.genvisis.cnv.qc.IlluminaMarkerBlast;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.gwas.MarkerQC.QC_METRIC;
import org.genvisis.gwas.Qc;
import org.genvisis.gwas.RelationAncestryQc;
import org.genvisis.qsub.Qsub;

public class GenvisisWorkflow {

  private static final String NUM_THREADS_DESC = "Number of Threads to Use";
  static final String PROJ_PROP_UPDATE_STR = " org.genvisis.cnv.filesys.Project proj=";
  public static final String PLINK_SUBDIR = "plink/";
  static final String PLINKROOT = "plink";
  final Project proj;
  private final SortedSet<Step> steps;
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

    steps = Collections.unmodifiableSortedSet(generateSteps(!project.IS_PC_CORRECTED_PROJECT.getValue()));
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

  private SortedSet<Step> generateSteps(boolean allowCorrectionStep) {
    StepBuilder sb = new StepBuilder();

    Step parseSamplesStep;
    Step markerBlastStep;
    if (proj.getArrayType() == ARRAY.AFFY_GW6 || proj.getArrayType() == ARRAY.AFFY_GW6_CN
        || proj.getArrayType() == ARRAY.AFFY_AXIOM) {
      parseSamplesStep = sb.generateParseSamplesStep();
      markerBlastStep = sb.generateAffyMarkerBlastAnnotationStep(parseSamplesStep);
    } else {
      Step markerPositions = sb.generateIlluminaMarkerPositionsStep();
      parseSamplesStep = sb.generateParseSamplesStep(markerPositions);
      markerBlastStep = sb.generateIlluminaMarkerBlastAnnotationStep(parseSamplesStep);
    }
    Step createSampleDataStep = sb.generateCreateSampleDataStep(parseSamplesStep);
    Step transposeStep = sb.generateTransposeStep(parseSamplesStep);
    Step gcModelStep = sb.generateGCModelStep();
    Step sampleQCStep = sb.generateSampleQCStep(parseSamplesStep);
    sb.generateMarkerQCStep(parseSamplesStep);
    sb.generateSexChecksStep(parseSamplesStep, markerBlastStep, createSampleDataStep, transposeStep,
                             sampleQCStep);
    sb.generateABLookupStep(parseSamplesStep);
    Step plinkExportStep = sb.generatePlinkExportStep(parseSamplesStep);
    Step gwasQCStep = sb.generateGwasQCStep(plinkExportStep);
    Step ancestryStep = sb.generateAncestryStep(gwasQCStep);
    sb.generateFurtherAnalysisQCStep(plinkExportStep, gwasQCStep, ancestryStep);
    sb.generateMosaicArmsStep(parseSamplesStep);
    sb.generateAnnotateSampleDataStep(sampleQCStep, createSampleDataStep, gwasQCStep);
    sb.generateMitoCNEstimateStep(transposeStep);
    Step pfbStep = sb.generatePFBStep(parseSamplesStep);
    sb.generateSexCentroidsStep();
    sb.generateCNVStep(pfbStep, gcModelStep);
    if (allowCorrectionStep) {
      sb.generatePCCorrectedProjectStep(parseSamplesStep);
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

  public static int resolveThreads(Project proj, String arg) {
    int numThreads = Requirement.checkIntArgOrNeg1(arg);
    if (numThreads <= 0) {
      numThreads = proj.NUM_THREADS.getValue();
    }
    return numThreads;
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
    GenvisisWorkflow workflow = new GenvisisWorkflow(proj, null);
    StepBuilder sb = new StepBuilder(workflow);

    Requirement threadsReq = workflow.getNumThreadsReq();
    Step createMkrPos = man != null ? sb.generateIlluminaMarkerPositionsStep() : null;
    Step parseSamples = sb.generateParseSamplesStep(createMkrPos);
    Step transpose = sb.generateTransposeStep(parseSamples);
    Step sampleData = sb.generateCreateSampleDataStep(parseSamples);
    Step blast = sb.generateIlluminaMarkerBlastAnnotationStep(parseSamples);
    Step sampleQc = sb.generateSampleQCStep(parseSamples);
    Step markerQc = sb.generateMarkerQCStep(parseSamples);
    Step sexChecks = sb.generateSexChecksStep(parseSamples, blast, sampleData, transpose, sampleQc);
    Step exportPlink = sb.generatePlinkExportStep(parseSamples);
    Step gwasQc = sb.generateGwasQCStep(exportPlink);
    Step ancestry = sb.generateAncestryStep(gwasQc);
    Step faqcStep = sb.generateFurtherAnalysisQCStep(exportPlink, gwasQc, ancestry);

    Map<Requirement, String> stepReqs;
    Map<Step, Map<Requirement, String>> varMap = new HashMap<>();

    if (createMkrPos != null) {
      stepReqs = createMkrPos.getDefaultRequirementValues();
      for (Requirement r1 : stepReqs.keySet()) {
        if (r1.getDescription().contains("Manifest")) {
          stepReqs.put(r1, man);
        }
      }
      varMap.put(createMkrPos, stepReqs);
    }
    if (parseSource) {
      varMap.put(parseSamples, parseSamples.getDefaultRequirementValues());
      varMap.put(transpose, transpose.getDefaultRequirementValues());
      stepReqs = blast.getDefaultRequirementValues();
      for (Requirement r1 : stepReqs.keySet()) {
        if (r1.getDescription().equalsIgnoreCase(IlluminaMarkerBlast.DESC_MANIFEST)) {
          stepReqs.put(r1, man);
          break;
        }
      }
      varMap.put(blast, stepReqs);

      stepReqs = sampleData.getDefaultRequirementValues();
      for (Requirement r1 : stepReqs.keySet()) {
        if (r1.getDescription().equalsIgnoreCase(StepBuilder.CREATE_MINIMAL_SAMPDATA_REQUIREMENT)) {
          stepReqs.put(r1, "false");
          break;
        }
      }
      varMap.put(sampleData, stepReqs);
    }
    varMap.put(sampleQc, sampleQc.getDefaultRequirementValues());
    varMap.put(markerQc, markerQc.getDefaultRequirementValues());

    stepReqs = sexChecks.getDefaultRequirementValues();
    for (Requirement r1 : stepReqs.keySet()) {
      if (r1.getDescription().equals(StepBuilder.NO_CROSS_HYBE_REQUIREMENT)) {
        stepReqs.put(r1, "true");
      } else if (r1.getDescription().equals(StepBuilder.ADD_ESTSEX_TO_SAMPDATA_REQUIREMENT)) {
        stepReqs.put(r1, "true");
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
      Requirement r = null;
      for (Requirement r1 : stepReqs.keySet()) {
        if (r1.getDescription().equals(StepBuilder.PUTATIVE_WHITE_FILE_DESCRIPTION)) {
          r = r1;
          break;
        }
      }
      if (r == null) {
        throw new IllegalStateException();
      }
      stepReqs.put(r, putativeWhitesFile);
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
        if (varMap.get(s).containsKey(threadsReq)) {
          varMap.get(s).put(threadsReq, Integer.toString(numThreads));
        }
      }
    }

    String s0 = createMkrPos == null ? "" : createMkrPos.getCommandLine(proj, varMap);
    String s1 = parseSamples.getCommandLine(proj, varMap);
    String s2 = transpose.getCommandLine(proj, varMap);
    String s3 = sampleData.getCommandLine(proj, varMap);
    String s4 = blast.getCommandLine(proj, varMap);
    String s5 = sampleQc.getCommandLine(proj, varMap);
    String s6 = markerQc.getCommandLine(proj, varMap);
    String s7 = sexChecks.getCommandLine(proj, varMap);
    String s8 = exportPlink.getCommandLine(proj, varMap);
    String s9 = gwasQc.getCommandLine(proj, varMap);
    String s10 = ancestry.getCommandLine(proj, varMap);
    String s11 = faqcStep.getCommandLine(proj, varMap);

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

  private static void fixQCThreshs(Map<Requirement, String> reqMap,
                                   Map<QC_METRIC, String> threshMap) {
    Map<String, String> qc = new HashMap<>();
    for (QC_METRIC met : threshMap.keySet()) {
      qc.put(met.getUserDescription(), threshMap.get(met));
    }
    for (Requirement r1 : reqMap.keySet()) {
      if (qc.containsKey(r1.getDescription())) {
        reqMap.put(r1, qc.get(r1.getDescription()));
      }
    }
  }

  public static void setupCNVCalling(String projectProperties) {
    Project pcProj = new Project(projectProperties);
    StepBuilder sb = new StepBuilder();
    Step transpose = sb.generateTransposeStep(null);
    // Create new sample data, run sex checks?
    Step gc = sb.generateGCModelStep();
    Step pfb = sb.generatePFBStep(null);
    Step cent = sb.generateSexCentroidsStep();
    Step cnv = sb.generateCNVStep(pfb, gc);
    Map<Step, Map<Requirement, String>> stepOpts = new HashMap<>();
    HashMap<Requirement, String> cnvOpts = new HashMap<>();
    List<Requirement> reqs = cnv.getRequirements().getFlatRequirementsList();
    for (Requirement req : reqs) {
      if (req.getDescription().equals(CNVCaller.CNV_SCOPE_DESC)) {
        cnvOpts.put(req, CNVCaller.CALLING_SCOPE.BOTH.toString());
      } else if (req.getDescription().equals(NUM_THREADS_DESC)) {
        cnvOpts.put(req, "" + (Runtime.getRuntime().availableProcessors() - 1));
      }
    }
    stepOpts.put(cnv, cnvOpts);

    String s1 = transpose.getCommandLine(pcProj, null);
    String s2 = gc.getCommandLine(pcProj, null);
    String s3 = pfb.getCommandLine(pcProj, null);
    String s4 = cent.getCommandLine(pcProj, null);
    String s5 = cnv.getCommandLine(pcProj, stepOpts);

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
