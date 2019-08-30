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
import org.pankratzlab.utils.gwas.Qc;
import org.pankratzlab.utils.gwas.RelationAncestryQc;

public class AnnotateSampleDataStep extends Step {

  public static final String NAME = "Annotate Sample Data File";
  public static final String DESC = "";

  public static AnnotateSampleDataStep create(Project proj, final Step sampleQCStep,
                                              final Step createSampleDataStep,
                                              final Step gwasQCStep) {
    final Requirement<Step> sampleQCStepReq = new Requirement.StepRequirement(sampleQCStep);
    final Requirement<Step> createSampleDataStepReq = new Requirement.StepRequirement(createSampleDataStep);
    final Requirement<Boolean> skipIDingDuplicatesReq = new Requirement.BoolRequirement("skipDuplicates",
                                                                                        "Skip identifying duplicates",
                                                                                        false);
    final Requirement<Step> gwasQCStepReq = new Requirement.StepRequirement(gwasQCStep);
    final Requirement<Boolean> notGcCorrectedLrrSdReq = new Requirement.BoolRequirement("dontUseGCLRRSD",
                                                                                        "Do not use GC corrected LRR SD?",
                                                                                        false);
    final Requirement<String> gcCorrectedLrrSdReq = new Requirement<String>("GC Corrected LRR SD must exist in Sample QC File") {

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
    final Requirement<Integer> numQReq = new Requirement.PosIntRequirement("numQuantiles",
                                                                           "Number of Quantiles to Generate",
                                                                           10);
    final Requirement<Boolean> replaceFIDIIDReq = new Requirement.OptionalBoolRequirement("replaceIDs",
                                                                                          "Replace FID and IID with data from Pedigree",
                                                                                          false);

    final RequirementSet reqSet = RequirementSetBuilder.and().add(sampleQCStepReq)
                                                       .add(createSampleDataStepReq)
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(skipIDingDuplicatesReq)
                                                                                 .add(gwasQCStepReq))
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(notGcCorrectedLrrSdReq)
                                                                                 .add(gcCorrectedLrrSdReq))
                                                       .add(numQReq).add(replaceFIDIIDReq);
    return new AnnotateSampleDataStep(proj, replaceFIDIIDReq, numQReq, notGcCorrectedLrrSdReq,
                                      skipIDingDuplicatesReq, reqSet);
  }

  final Project proj;
  final Requirement<Boolean> skipIDingDuplicatesReq;
  final Requirement<Boolean> replaceFIDIIDReq;
  final Requirement<Integer> numQReq;
  final Requirement<Boolean> notGcCorrectedLrrSdReq;

  public AnnotateSampleDataStep(Project proj, Requirement<Boolean> replaceIDReq,
                                Requirement<Integer> numQReq, Requirement<Boolean> notGCReq,
                                Requirement<Boolean> skipIDingDupReq, RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.noneOf(Requirement.Flag.class));
    this.proj = proj;
    this.skipIDingDuplicatesReq = skipIDingDupReq;
    this.replaceFIDIIDReq = replaceIDReq;
    this.numQReq = numQReq;
    this.notGcCorrectedLrrSdReq = notGCReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    // nothing to do
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
    SampleQC.parseAndAddToSampleDataWithoutExcludes(proj, numQ, 0, false, gcCorrectedLrrSd,
                                                    duplicatesSetFile, correctFidIids);
  }

  @Override
  public String getCommandLine(Variables variables) {
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

    StringBuilder cmd = new StringBuilder();
    cmd.append(Files.getRunString()).append(" ").append(SampleQC.class.getName())
       .append(" proj=" + projPropFile + " numQ=" + numQ + " justQuantiles=false"
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

    int[] facts;

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
