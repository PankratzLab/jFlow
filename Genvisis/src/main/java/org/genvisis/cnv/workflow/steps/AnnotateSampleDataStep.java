package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;
import java.util.Map;
import java.util.Set;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.SampleQC;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.Variables;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.ext;
import org.pankratzlab.gwas.RelationAncestryQc;
import org.pankratzlab.utils.gwas.Qc;

public class AnnotateSampleDataStep extends Step {

  public static final String NAME = "Annotate Sample Data File";
  public static final String DESC = "";

  public static AnnotateSampleDataStep create(Project proj, final Step sampleQCStep,
                                              final Step createSampleDataStep,
                                              final Step gwasQCStep) {
    final Requirement<Step> sampleQCStepReq = new Requirement.StepRequirement(sampleQCStep);
    final Requirement<Step> createSampleDataStepReq = new Requirement.StepRequirement(createSampleDataStep);
    final Requirement<Boolean> skipIDingDuplicatesReq = new Requirement.BoolRequirement("Skip identifying duplicates",
                                                                                        false);
    final Requirement<Step> gwasQCStepReq = new Requirement.StepRequirement(gwasQCStep);
    final Requirement<Boolean> notGcCorrectedLrrSdReq = new Requirement.BoolRequirement("Do not use GC corrected LRR SD?",
                                                                                        false);
    final Requirement<String> gcCorrectedLrrSdReq = new Requirement<String>("GC Corrected LRR SD must exist in Sample QC File",
                                                                            Requirement.RequirementInputType.NONE) {

      @Override
      public boolean checkRequirement(String arg, Set<Step> stepSelections,
                                      Map<Step, Variables> variables) {
        String sampleQCFile = proj.SAMPLE_QC_FILENAME.getValue();
        return Files.exists(sampleQCFile)
               && ext.indexOfStr("LRR_SD_Post_Correction",
                                 Files.getHeaderOfFile(sampleQCFile, proj.getLog())) != -1;
      }

      @Override
      public String parseValue(String raw) {
        return raw;
      }

    };
    final Requirement<Double> lrrSdThresholdReq = new Requirement.DoubleRequirement("LRR SD Threshold",
                                                                                    proj.LRRSD_CUTOFF.getValue(),
                                                                                    proj.LRRSD_CUTOFF.getMinValue(),
                                                                                    proj.LRRSD_CUTOFF.getMaxValue());

    final Requirement<Double> callrateThresholdReq = new Requirement.DoubleRequirement("Callrate Threshold",
                                                                                       proj.SAMPLE_CALLRATE_THRESHOLD.getValue(),
                                                                                       proj.SAMPLE_CALLRATE_THRESHOLD.getMinValue(),
                                                                                       proj.SAMPLE_CALLRATE_THRESHOLD.getMaxValue());
    final Requirement<Integer> numQReq = new Requirement.PosIntRequirement("Number of Quantiles to Generate",
                                                                           10);
    final Requirement<Boolean> replaceFIDIIDReq = new Requirement.OptionalBoolRequirement("Replace FID and IID with data from Pedigree",
                                                                                          false);

    final RequirementSet reqSet = RequirementSetBuilder.and().add(sampleQCStepReq)
                                                       .add(createSampleDataStepReq)
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(skipIDingDuplicatesReq)
                                                                                 .add(gwasQCStepReq))
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(notGcCorrectedLrrSdReq)
                                                                                 .add(gcCorrectedLrrSdReq))
                                                       .add(lrrSdThresholdReq)
                                                       .add(callrateThresholdReq).add(numQReq)
                                                       .add(replaceFIDIIDReq);
    return new AnnotateSampleDataStep(proj, replaceFIDIIDReq, numQReq, notGcCorrectedLrrSdReq,
                                      skipIDingDuplicatesReq, lrrSdThresholdReq,
                                      callrateThresholdReq, reqSet);
  }

  final Project proj;
  final Requirement<Boolean> skipIDingDuplicatesReq;
  final Requirement<Double> callrateThresholdReq;
  final Requirement<Double> lrrSdThresholdReq;
  final Requirement<Boolean> replaceFIDIIDReq;
  final Requirement<Integer> numQReq;
  final Requirement<Boolean> notGcCorrectedLrrSdReq;

  public AnnotateSampleDataStep(Project proj, Requirement<Boolean> replaceIDReq,
                                Requirement<Integer> numQReq, Requirement<Boolean> notGCReq,
                                Requirement<Boolean> skipIDingDupReq, Requirement<Double> lrrSdReq,
                                Requirement<Double> callrateReq, RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.noneOf(Requirement.Flag.class));
    this.proj = proj;
    this.skipIDingDuplicatesReq = skipIDingDupReq;
    this.lrrSdThresholdReq = lrrSdReq;
    this.callrateThresholdReq = callrateReq;
    this.replaceFIDIIDReq = replaceIDReq;
    this.numQReq = numQReq;
    this.notGcCorrectedLrrSdReq = notGCReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    double projLrrSdThreshold = proj.LRRSD_CUTOFF.getValue();
    double lrrSdThreshold = variables.get(lrrSdThresholdReq);
    double projCallrateThreshold = proj.SAMPLE_CALLRATE_THRESHOLD.getValue();
    double callrateThreshold = variables.get(callrateThresholdReq);

    if (projLrrSdThreshold != lrrSdThreshold) {
      proj.LRRSD_CUTOFF.setValue(lrrSdThreshold);
    }
    if (projCallrateThreshold != callrateThreshold) {
      proj.SAMPLE_CALLRATE_THRESHOLD.setValue(callrateThreshold);
    }
  }

  @Override
  public void run(Variables variables) {
    boolean checkDuplicates = !variables.get(skipIDingDuplicatesReq).booleanValue();
    String duplicatesSetFile = null;
    if (checkDuplicates) {
      duplicatesSetFile = GenvisisWorkflow.getPlinkDir(proj) + Qc.QC_SUBDIR
                          + RelationAncestryQc.GENOME_DIR + GenvisisWorkflow.PLINKROOT
                          + ".genome_duplicatesSet.dat";
    }
    boolean gcCorrectedLrrSd = !variables.get(notGcCorrectedLrrSdReq).booleanValue();
    int numQ = variables.get(numQReq);
    boolean correctFidIids = variables.get(replaceFIDIIDReq);
    SampleQC.parseAndAddToSampleData(proj, numQ, 0, false, gcCorrectedLrrSd, duplicatesSetFile,
                                     correctFidIids);
  }

  @Override
  public String getCommandLine(Variables variables) {

    double projLrrSdThreshold = proj.LRRSD_CUTOFF.getValue();
    double lrrSdThreshold = variables.get(lrrSdThresholdReq);
    double projCallrateThreshold = proj.SAMPLE_CALLRATE_THRESHOLD.getValue();
    double callrateThreshold = variables.get(callrateThresholdReq);

    String projPropFile = proj.getPropertyFilename();

    boolean checkDuplicates = !variables.get(skipIDingDuplicatesReq).booleanValue();
    String duplicatesSetFile = null;
    if (checkDuplicates) {
      duplicatesSetFile = GenvisisWorkflow.getPlinkDir(proj) + Qc.QC_SUBDIR
                          + RelationAncestryQc.GENOME_DIR + GenvisisWorkflow.PLINKROOT
                          + ".genome_duplicatesSet.dat";
    }
    boolean gcCorrectedLrrSd = !variables.get(notGcCorrectedLrrSdReq).booleanValue();
    int numQ = variables.get(numQReq);
    boolean correctFidIids = variables.get(replaceFIDIIDReq);

    String kvCmd = "";

    if (projLrrSdThreshold != lrrSdThreshold) {
      kvCmd += " LRRSD_CUTOFF=" + lrrSdThreshold;
    }
    if (projCallrateThreshold != callrateThreshold) {
      kvCmd += " SAMPLE_CALLRATE_THRESHOLD=" + callrateThreshold;
    }

    StringBuilder cmd = new StringBuilder();
    if (kvCmd.length() > 0) {
      cmd.append(Files.getRunString()).append(GenvisisWorkflow.PROJ_PROP_UPDATE_STR + projPropFile)
         .append(kvCmd).append("\n");
    }
    cmd.append(Files.getRunString())
       .append(" cnv.qc.SampleQC proj=" + projPropFile + " numQ=" + numQ + " justQuantiles=false"
               + " gcCorrectedLrrSd=" + gcCorrectedLrrSd + " duplicatesSetFile=" + duplicatesSetFile
               + " correctFidIids=" + correctFidIids);
    return cmd.toString();
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    String sampleDataFile = proj.SAMPLE_DATA_FILENAME.getValue();
    if (!Files.exists(sampleDataFile)) {
      return false;
    }
    String[] header = Files.getHeaderOfFile(sampleDataFile, proj.getLog());

    // These columns should always added by SampleQC
    String[] baseHeader = {SampleQC.EXCLUDE_HEADER, "ExcludeNote"};

    int[] facts = ext.indexFactors(baseHeader, header, false);
    for (int i : facts) {
      if (i == -1) {
        return false;
      }
    }

    boolean checkDuplicates = !variables.get(skipIDingDuplicatesReq).booleanValue();
    File duplicateSetFile = new File(GenvisisWorkflow.getPlinkDir(proj) + Qc.QC_SUBDIR
                                     + RelationAncestryQc.GENOME_DIR + GenvisisWorkflow.PLINKROOT
                                     + ".genome_duplicatesSet.dat");
    // If there was no duplicate file than these headers won't be added
    if (checkDuplicates && duplicateSetFile.exists()) {
      // These columns are only added if checkDuplicates occurred
      String[] dupHeader = {SampleQC.DUPLICATE_ID_HEADER, "Use", "UseNote", "Use_cnv",
                            "Use_cnvNote"};
      facts = ext.indexFactors(dupHeader, header, false);
      for (int i : facts) {
        if (i == -1) {
          return false;
        }
      }
    }
    return true;
  }

}
