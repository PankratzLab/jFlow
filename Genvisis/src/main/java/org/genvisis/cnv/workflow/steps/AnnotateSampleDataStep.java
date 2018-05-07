package org.genvisis.cnv.workflow.steps;

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
import org.genvisis.common.Files;
import org.genvisis.common.ext;
import org.genvisis.gwas.Qc;
import org.genvisis.gwas.RelationAncestryQc;

public class AnnotateSampleDataStep extends Step {

  public static final String NAME = "";
  public static final String DESC = "";

  public static AnnotateSampleDataStep create(Project proj, final Step sampleQCStep,
                                              final Step createSampleDataStep,
                                              final Step gwasQCStep, double priority) {
    final Requirement sampleQCStepReq = new Requirement.StepRequirement(sampleQCStep);
    final Requirement createSampleDataStepReq = new Requirement.StepRequirement(createSampleDataStep);
    final Requirement skipIDingDuplicatesReq = new Requirement.BoolRequirement("Skip identifying duplicates",
                                                                               false);
    final Requirement gwasQCStepReq = new Requirement.StepRequirement(gwasQCStep);
    final Requirement notGcCorrectedLrrSdReq = new Requirement.BoolRequirement("Do not use GC corrected LRR SD?",
                                                                               false);
    final Requirement gcCorrectedLrrSdReq = new Requirement("GC Corrected LRR SD must exist in Sample QC File",
                                                            Requirement.RequirementInputType.NONE) {

      @Override
      public boolean checkRequirement(Project proj, String arg, Set<Step> stepSelections,
                                      Map<Step, Map<Requirement, String>> variables) {
        String sampleQCFile = proj.SAMPLE_QC_FILENAME.getValue();
        return Files.exists(sampleQCFile)
               && ext.indexOfStr("LRR_SD_Post_Correction",
                                 Files.getHeaderOfFile(sampleQCFile, proj.getLog())) != -1;
      }

    };
    final Requirement lrrSdThresholdReq = new Requirement.DoubleRequirement("LRR SD Threshold",
                                                                            proj.LRRSD_CUTOFF.getValue(),
                                                                            proj.LRRSD_CUTOFF.getMinValue(),
                                                                            proj.LRRSD_CUTOFF.getMaxValue());

    final Requirement callrateThresholdReq = new Requirement.DoubleRequirement("Callrate Threshold",
                                                                               proj.SAMPLE_CALLRATE_THRESHOLD.getValue(),
                                                                               proj.SAMPLE_CALLRATE_THRESHOLD.getMinValue(),
                                                                               proj.SAMPLE_CALLRATE_THRESHOLD.getMaxValue());
    final Requirement numQReq = new Requirement.PosIntRequirement("Number of Quantiles to Generate",
                                                                  10);
    final Requirement replaceFIDIIDReq = new Requirement.OptionalBoolRequirement("Replace FID and IID with data from Pedigree",
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
    return new AnnotateSampleDataStep(replaceFIDIIDReq, numQReq, notGcCorrectedLrrSdReq,
                                      skipIDingDuplicatesReq, lrrSdThresholdReq,
                                      callrateThresholdReq, reqSet, priority);
  }

  final Requirement skipIDingDuplicatesReq;
  final Requirement callrateThresholdReq;
  final Requirement lrrSdThresholdReq;
  final Requirement replaceFIDIIDReq;
  final Requirement numQReq;
  final Requirement notGcCorrectedLrrSdReq;

  public AnnotateSampleDataStep(Requirement replaceIDReq, Requirement numQReq, Requirement notGCReq,
                                Requirement skipIDingDupReq, Requirement lrrSdReq,
                                Requirement callrateReq, RequirementSet reqSet, double priority) {
    super(NAME, DESC, reqSet, EnumSet.noneOf(Requirement.Flag.class), priority);
    this.skipIDingDuplicatesReq = skipIDingDupReq;
    this.lrrSdThresholdReq = lrrSdReq;
    this.callrateThresholdReq = callrateReq;
    this.replaceFIDIIDReq = replaceIDReq;
    this.numQReq = numQReq;
    this.notGcCorrectedLrrSdReq = notGCReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Project proj,
                                           Map<Step, Map<Requirement, String>> variables) {
    double projLrrSdThreshold = proj.LRRSD_CUTOFF.getValue();
    double lrrSdThreshold = Double.parseDouble(variables.get(this).get(lrrSdThresholdReq));
    double projCallrateThreshold = proj.SAMPLE_CALLRATE_THRESHOLD.getValue();
    double callrateThreshold = Double.parseDouble(variables.get(this).get(callrateThresholdReq));

    if (projLrrSdThreshold != lrrSdThreshold) {
      proj.LRRSD_CUTOFF.setValue(lrrSdThreshold);
    }
    if (projCallrateThreshold != callrateThreshold) {
      proj.SAMPLE_CALLRATE_THRESHOLD.setValue(callrateThreshold);
    }
  }

  @Override
  public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
    boolean checkDuplicates = !Boolean.parseBoolean(variables.get(this)
                                                             .get(skipIDingDuplicatesReq));
    String duplicatesSetFile = null;
    if (checkDuplicates) {
      duplicatesSetFile = GenvisisWorkflow.getPlinkDir(proj) + Qc.QC_SUBDIR
                          + RelationAncestryQc.GENOME_DIR + GenvisisWorkflow.PLINKROOT
                          + ".genome_duplicatesSet.dat";
    }
    boolean gcCorrectedLrrSd = !Boolean.parseBoolean(variables.get(this)
                                                              .get(notGcCorrectedLrrSdReq));
    int numQ = Integer.parseInt(variables.get(this).get(numQReq));
    boolean correctFidIids = Boolean.parseBoolean(variables.get(this).get(replaceFIDIIDReq));
    SampleQC.parseAndAddToSampleData(proj, numQ, 0, false, gcCorrectedLrrSd, duplicatesSetFile,
                                     correctFidIids);
  }

  @Override
  public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {

    double projLrrSdThreshold = proj.LRRSD_CUTOFF.getValue();
    double lrrSdThreshold = Double.parseDouble(variables.get(this).get(lrrSdThresholdReq));
    double projCallrateThreshold = proj.SAMPLE_CALLRATE_THRESHOLD.getValue();
    double callrateThreshold = Double.parseDouble(variables.get(this).get(callrateThresholdReq));

    String projPropFile = proj.getPropertyFilename();

    boolean checkDuplicates = !Boolean.parseBoolean(variables.get(this)
                                                             .get(skipIDingDuplicatesReq));
    String duplicatesSetFile = null;
    if (checkDuplicates) {
      duplicatesSetFile = GenvisisWorkflow.getPlinkDir(proj) + Qc.QC_SUBDIR
                          + RelationAncestryQc.GENOME_DIR + GenvisisWorkflow.PLINKROOT
                          + ".genome_duplicatesSet.dat";
    }
    boolean gcCorrectedLrrSd = !Boolean.parseBoolean(variables.get(this)
                                                              .get(notGcCorrectedLrrSdReq));
    int numQ = Integer.parseInt(variables.get(this).get(numQReq));
    boolean correctFidIids = Boolean.parseBoolean(variables.get(this).get(replaceFIDIIDReq));

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
  public boolean checkIfOutputExists(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String sampleDataFile = proj.SAMPLE_DATA_FILENAME.getValue();
    if (!Files.exists(sampleDataFile)) {
      return false;
    }
    boolean checkDuplicates = !Boolean.parseBoolean(variables.get(this)
                                                             .get(skipIDingDuplicatesReq));
    String[] header = Files.getHeaderOfFile(sampleDataFile, proj.getLog());
    if (checkDuplicates
        && ext.indexOfStr(SampleQC.DUPLICATE_ID_HEADER, header, false, true) == -1) {
      return false;
    }
    String[] reqHdr = {SampleQC.EXCLUDE_HEADER, "ExcludeNote", "Use", "UseNote", "Use_cnv",
                       "Use_cnvNote"};
    int[] facts = ext.indexFactors(reqHdr, header, false);
    for (int i : facts) {
      if (i == -1) {
        return false;
      }
    }
    return true;
  }

}
