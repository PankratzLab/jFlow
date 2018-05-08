package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;
import java.util.Map;
import java.util.Set;
import org.genvisis.cnv.filesys.Centroids;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.hmm.CNVCaller;
import org.genvisis.cnv.hmm.CNVCaller.CALLING_SCOPE;
import org.genvisis.cnv.hmm.CNVCaller.PFB_MANAGEMENT_TYPE;
import org.genvisis.cnv.prop.Property;
import org.genvisis.cnv.workflow.GenvisisWorkflow;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.cnv.workflow.Variables;
import org.genvisis.common.Files;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;

public class CallCNVsStep extends Step {

  public static final String NAME = "Call CNVs";
  public static final String DESC = "";

  public static CallCNVsStep create(Project proj, Step pfbStep, Step gcModelStep,
                                    Requirement<Integer> numThreadsReq) {
    final Requirement<File> hmmFile = new Requirement.FileRequirement("Hidden Markov Model File Must Exist",
                                                                      new File(proj.HMM_FILENAME.getValue()));
    final Requirement<Step> pfbStepReq = new Requirement.StepRequirement(pfbStep);
    final Requirement<File> pfbFileReq = new Requirement.FileRequirement("PFB File Must Exist",
                                                                         new File(proj.CUSTOM_PFB_FILENAME.getValue()));
    final Requirement<Step> gcModelStepReq = new Requirement.StepRequirement(gcModelStep);
    final Requirement<File> gcModelFileReq = new Requirement.FileRequirement("GCMODEL File Must Exist",
                                                                             new File(proj.GC_MODEL_FILENAME.getValue()));
    final Requirement<CALLING_SCOPE> callingTypeReq = new Requirement.EnumRequirement<CALLING_SCOPE>(CNVCaller.CNV_SCOPE_DESC,
                                                                                                     CNVCaller.CALLING_SCOPE.AUTOSOMAL);
    final Requirement<Boolean> useCentroidsReq = new Requirement.OptionalBoolRequirement("If calling chromosomal CNVs, use sex-specific centroids to recalculate LRR/BAF values?",
                                                                                         true);
    final Requirement<File> outputFileReq = new Requirement.OutputFileRequirement("Output filename.",
                                                                                  new File("cnvs/genvisis.cnv")) {

      @Override
      public boolean checkRequirement(Project proj, String arg, Set<Step> stepSelections,
                                      Map<Step, Variables> variables) {
        return super.checkRequirement(proj, proj.PROJECT_DIRECTORY.getValue() + arg, stepSelections,
                                      variables);
      }
    };

    final RequirementSet reqSet = RequirementSetBuilder.and().add(hmmFile)
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(pfbStepReq)
                                                                                 .add(pfbFileReq))
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(gcModelStepReq)
                                                                                 .add(gcModelFileReq))
                                                       .add(callingTypeReq).add(useCentroidsReq)
                                                       .add(numThreadsReq).add(outputFileReq);
    return new CallCNVsStep(hmmFile, pfbFileReq, gcModelFileReq, callingTypeReq, useCentroidsReq,
                            outputFileReq, numThreadsReq, reqSet);
  }

  final Requirement<File> hmmFile;
  final Requirement<File> pfbFileReq;
  final Requirement<File> gcModelFileReq;
  final Requirement<CALLING_SCOPE> callingTypeReq;
  final Requirement<Boolean> useCentroidsReq;
  final Requirement<File> outputFileReq;
  final Requirement<Integer> numThreadsReq;

  private CallCNVsStep(Requirement<File> hmmFile, Requirement<File> pfbFileReq,
                       Requirement<File> gcModelFileReq, Requirement<CALLING_SCOPE> callingTypeReq,
                       Requirement<Boolean> useCentroidsReq, Requirement<File> outputFileReq,
                       Requirement<Integer> numThreadsReq, RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.of(Requirement.Flag.MEMORY, Requirement.Flag.MULTITHREADED));
    this.hmmFile = hmmFile;
    this.pfbFileReq = pfbFileReq;
    this.gcModelFileReq = gcModelFileReq;
    this.callingTypeReq = callingTypeReq;
    this.useCentroidsReq = useCentroidsReq;
    this.outputFileReq = outputFileReq;
    this.numThreadsReq = numThreadsReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Project proj, Variables variables) {
    String hmmP = proj.HMM_FILENAME.getValue();
    String hmmG = variables.get(hmmFile).getAbsolutePath();
    if (!hmmP.equals(hmmG)) {
      proj.HMM_FILENAME.setValue(hmmG);
    }
    String pfbP = proj.CUSTOM_PFB_FILENAME.getValue();
    String pfbG = variables.get(pfbFileReq).getAbsolutePath();
    if (!pfbP.equals(pfbG)) {
      proj.CUSTOM_PFB_FILENAME.setValue(pfbG);
    }
    String gcmP = proj.GC_MODEL_FILENAME.getValue();
    String gcmG = variables.get(gcModelFileReq).getAbsolutePath();
    if (!gcmP.equals(gcmG)) {
      proj.GC_MODEL_FILENAME.setValue(gcmG);
    }
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
    GenvisisWorkflow.maybeSetProjNumThreads(proj, numThreads);
  }

  @Override
  public void run(Project proj, Variables variables) {
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
    GenvisisWorkflow.maybeSetProjNumThreads(proj, numThreads);
    String output = variables.get(outputFileReq).getPath(); // gets PROJ_DIR prepended, so NOT ABSOLUTE

    CALLING_SCOPE scope = variables.get(callingTypeReq);

    (new File(ext.parseDirectoryOfFile(proj.PROJECT_DIRECTORY.getValue() + output))).mkdirs();

    String[] samples = proj.getSamples();
    boolean useCentroids = Boolean.valueOf(variables.get(useCentroidsReq));
    Centroids[] cents = new Centroids[] {null, null};
    if (useCentroids) {
      if (Files.exists(proj.SEX_CENTROIDS_FEMALE_FILENAME.getValue())
          && Files.exists(proj.SEX_CENTROIDS_MALE_FILENAME.getValue())) {
        cents[0] = Centroids.load(proj.SEX_CENTROIDS_MALE_FILENAME.getValue());
        cents[1] = Centroids.load(proj.SEX_CENTROIDS_FEMALE_FILENAME.getValue());
      }
    }

    if (scope != CALLING_SCOPE.CHROMOSOMAL) {
      CNVCaller.callAutosomalCNVs(proj, output, samples, null, null, null,
                                  CNVCaller.DEFAULT_MIN_SITES, CNVCaller.DEFAULT_MIN_CONF,
                                  PFB_MANAGEMENT_TYPE.PENNCNV_DEFAULT, numThreads, 1);
      String file = proj.PROJECT_DIRECTORY.getValue() + output;
      if (Files.exists(file)) {
        proj.CNV_FILENAMES.addValue(file);
      }
    }
    if (scope != CALLING_SCOPE.AUTOSOMAL) {
      CNVCaller.callGenomeCnvs(proj, output, cents, null, CNVCaller.DEFAULT_MIN_SITES,
                               CNVCaller.DEFAULT_MIN_CONF, PFB_MANAGEMENT_TYPE.PENNCNV_DEFAULT,
                               numThreads, 1);
      String[] files = {proj.PROJECT_DIRECTORY.getValue() + output + "_23M.cnv",
                        proj.PROJECT_DIRECTORY.getValue() + output + "_23F.cnv",
                        proj.PROJECT_DIRECTORY.getValue() + output + "_24M.cnv"};
      for (String f : files) {
        if (Files.exists(f)) {
          proj.CNV_FILENAMES.addValue(f);
        }
      }
    }

    proj.saveProperties(new Property[] {proj.CNV_FILENAMES});
  }

  @Override
  public String getCommandLine(Project proj, Variables variables) {
    String kvCmd = "";

    String hmmP = proj.HMM_FILENAME.getValue();
    String hmmG = variables.get(hmmFile).getAbsolutePath();
    if (hmmG != null && !hmmP.equals(hmmG)) {
      kvCmd += " HMM_FILENAME=" + hmmG;
    }
    String pfbP = proj.CUSTOM_PFB_FILENAME.getValue();
    String pfbG = variables.get(pfbFileReq).getAbsolutePath();
    if (pfbG != null && !pfbP.equals(pfbG)) {
      kvCmd += " CUSTOM_PFB_FILENAME=" + pfbG;
    }
    String gcmP = proj.GC_MODEL_FILENAME.getValue();
    String gcmG = variables.get(gcModelFileReq).getAbsolutePath();
    if (gcmG != null && !gcmP.equals(gcmG)) {
      kvCmd += " GC_MODEL_FILENAME=" + gcmG;
    }

    int numThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
    if (numThreads != proj.NUM_THREADS.getValue()) {
      kvCmd += " " + proj.NUM_THREADS.getName() + "=" + numThreads;
    }
    String projPropFile = proj.getPropertyFilename();
    StringBuilder cmd = new StringBuilder();
    if (kvCmd.length() > 0) {
      cmd.append(Files.getRunString()).append(GenvisisWorkflow.PROJ_PROP_UPDATE_STR + projPropFile)
         .append(kvCmd).append("\n");
    }

    boolean useCentroids = Boolean.valueOf(variables.get(useCentroidsReq));
    CALLING_SCOPE scope = variables.get(callingTypeReq);

    String autoCmd = cmd.append(Files.getRunString())
                        .append(" cnv.hmm.CNVCaller proj=" + projPropFile)
                        .append(" out=" + variables.get(outputFileReq).getPath()).append(" ")
                        .append(PSF.Ext.NUM_THREADS_COMMAND).append(numThreads).toString();
    String genomeCmd = autoCmd + " -genome";
    if (!useCentroids) {
      genomeCmd += " -noCentroids";
    }

    switch (scope) {
      case AUTOSOMAL:
        return autoCmd;
      case CHROMOSOMAL:
        return genomeCmd;
      case BOTH:
        return autoCmd + "\n" + genomeCmd;
    }
    return autoCmd;
  }

  @Override
  public boolean checkIfOutputExists(Project proj, Variables variables) {
    String output = variables.get(outputFileReq).getPath();
    return Files.exists(proj.PROJECT_DIRECTORY.getValue() + output);
  }

}
