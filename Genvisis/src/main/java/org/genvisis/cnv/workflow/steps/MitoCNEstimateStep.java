package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.genvisis.cnv.analysis.pca.PCA;
import org.genvisis.cnv.analysis.pca.PCAPrep;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsCompute.PRE_PROCESSING_METHOD;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.MitoPipeline;
import org.genvisis.cnv.qc.GcAdjustor;
import org.genvisis.cnv.qc.IntensityMarkers;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.Requirement.BoolRequirement;
import org.genvisis.cnv.workflow.RequirementSet;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.cnv.workflow.Variables;
import org.pankratzlab.common.Files;

public class MitoCNEstimateStep extends Step {

  public static final String NAME = "Create Mitochondrial Copy-Number Estimates File";
  public static final String DESC = "";
  public static final String ALL_MITO_MARKERS_FILENAME = "AllMitoMarkers.txt";

  public static MitoCNEstimateStep create(Project proj, Step markersParsingStep,
                                          Requirement<Integer> numThreadsReq) {
    // FIXME http://genvisis.org/MitoPipeline/#illumina_marker_lists has illumina markers.. this
    // should be linked to, or
    // these steps split or something...
    final Requirement<Step> transposeStepReq = new Requirement.StepRequirement(markersParsingStep);

    final Requirement<Boolean> hasAnyMitoMarkersReq = new Requirement<Boolean>("Project must have mitochondrial markers",
                                                                               "", false) {

      @Override
      public boolean checkRequirement(String arg, Set<Step> stepSelections,
                                      Map<Step, Variables> variables) {
        return proj.getMarkerSet() != null && proj.getMarkerSet().getChrMap().containsKey((byte) 26)
               && proj.getMarkerSet().getChrMap().get((byte) 26).size() != 0;
      }

      @Override
      public Boolean parseValue(String raw) {
        return Boolean.valueOf(raw);
      }

    };
    final Requirement<File> medianMarkersReq = new Requirement.FileRequirement("cleanMitoMarkers",
                                                                               "A file containing mitochondrial markers with clean clusters must exist (see http://genvisis.org/MitoPipeline/scatter.html for more details).",
                                                                               new File(""));
    final BoolRequirement useAllMitoMarkersReq = new BoolRequirement("useAllMitoMarkers",
                                                                     "Use all mitochondrial markers",
                                                                     false);

    final Requirement<Double> lrrSdThresholdReq = new Requirement.DoubleRequirement("lrrSDThreshold",
                                                                                    "LRR SD threshold to filter samples.",
                                                                                    proj.LRRSD_CUTOFF.getValue(),
                                                                                    proj.LRRSD_CUTOFF.getMinValue(),
                                                                                    proj.LRRSD_CUTOFF.getMaxValue());
    final Requirement<Double> callrateThresholdReq = new Requirement.DoubleRequirement("callrateThreshold",
                                                                                       "Call rate threshold to filter markers.",
                                                                                       MitoPipeline.DEFAULT_MKR_CALLRATE_FILTER,
                                                                                       0.0, 1.0);
    final Requirement<Boolean> qcPassingOnlyReq = new Requirement.OptionalBoolRequirement("computePCsWithPassingSamps",
                                                                                          "Compute PCs with samples passing QC only",
                                                                                          true);
    final Requirement<Boolean> imputeNaNs = new Requirement.OptionalBoolRequirement("imputeNaNMean",
                                                                                    "Impute mean value for NaN",
                                                                                    true);
    final Requirement<Boolean> recomputeLrrPCMarkersReq = new Requirement.OptionalBoolRequirement("recomputePC",
                                                                                                  "Should recompute Log-R ratio for PC markers?",
                                                                                                  true);
    final Requirement<Boolean> recomputeLrrMedianMarkersReq = new Requirement.OptionalBoolRequirement("recomputeMedian",
                                                                                                      "Should recompute Log-R ratio for median markers?",
                                                                                                      true);
    final Requirement<Boolean> homozygousOnlyReq = new Requirement.OptionalBoolRequirement("homozygousOnly",
                                                                                           "Homozygous only?",
                                                                                           true);
    final Requirement<Integer> gcRegressionDistanceReq = new Requirement.PosIntRequirement("regressionDistance",
                                                                                           "Regression distance for the GC adjustment",
                                                                                           GcAdjustor.DEFAULT_REGRESSION_DISTANCE[0]);
    final Requirement<File> pcSelectionSamplesReq = new Requirement.OptionalFileRequirement("subsetFile",
                                                                                            "A file listing a subset of samples (DNA ID) to use for determining optimal PC selection, typically a list of unrelated and single race samples. If a list is not provided, only samples passing sample qc thresholds will be used.",
                                                                                            new File(""));
    final Requirement<File> externalBetaFileReq = new Requirement.OptionalFileRequirement("betaFile",
                                                                                          "An external beta file to optimize PC selection.",
                                                                                          new File(""));

    final RequirementSet reqSet = RequirementSetBuilder.and().add(transposeStepReq)
                                                       .add(hasAnyMitoMarkersReq)
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(medianMarkersReq)
                                                                                 .add(useAllMitoMarkersReq))
                                                       .add(lrrSdThresholdReq)
                                                       .add(callrateThresholdReq)
                                                       .add(qcPassingOnlyReq).add(imputeNaNs)
                                                       .add(recomputeLrrPCMarkersReq)
                                                       .add(recomputeLrrMedianMarkersReq)
                                                       .add(homozygousOnlyReq)
                                                       .add(gcRegressionDistanceReq)
                                                       .add(numThreadsReq)
                                                       .add(pcSelectionSamplesReq)
                                                       .add(externalBetaFileReq);
    return new MitoCNEstimateStep(proj, useAllMitoMarkersReq, medianMarkersReq, lrrSdThresholdReq,
                                  callrateThresholdReq, qcPassingOnlyReq, imputeNaNs,
                                  recomputeLrrPCMarkersReq, recomputeLrrMedianMarkersReq,
                                  homozygousOnlyReq, gcRegressionDistanceReq, pcSelectionSamplesReq,
                                  externalBetaFileReq, numThreadsReq, reqSet);
  }

  final Project proj;
  final Requirement<Boolean> useAllMitoMarkersReq;
  final Requirement<File> cleanMitoMarkersReq;
  final Requirement<Double> lrrSdThresholdReq;
  final Requirement<Double> callrateThresholdReq;
  final Requirement<Boolean> qcPassingOnlyReq;
  final Requirement<Boolean> imputeNaNs;
  final Requirement<Boolean> recomputeLrrPCMarkersReq;
  final Requirement<Boolean> recomputeLrrMedianMarkersReq;
  final Requirement<Boolean> homozygousOnlyReq;
  final Requirement<Integer> gcRegressionDistanceReq;
  final Requirement<File> pcSelectionSamplesReq;
  final Requirement<File> externalBetaFileReq;
  final Requirement<Integer> numThreadsReq;

  private MitoCNEstimateStep(Project proj, Requirement<Boolean> useAllMitoMarkersReq,
                             Requirement<File> cleanMitoMarkersReq,
                             Requirement<Double> lrrSdThresholdReq,
                             Requirement<Double> callrateThresholdReq,
                             Requirement<Boolean> qcPassingOnlyReq, Requirement<Boolean> imputeNaNs,
                             Requirement<Boolean> recomputeLrrPCMarkersReq,
                             Requirement<Boolean> recomputeLrrMedianMarkersReq,
                             Requirement<Boolean> homozygousOnlyReq,
                             Requirement<Integer> gcRegressionDistanceReq,
                             Requirement<File> pcSelectionSamplesReq,
                             Requirement<File> externalBetaFileReq,
                             Requirement<Integer> numThreadsReq, RequirementSet reqSet) {
    super(NAME, DESC, reqSet, EnumSet.of(Requirement.Flag.MULTITHREADED));
    this.proj = proj;
    this.useAllMitoMarkersReq = useAllMitoMarkersReq;
    this.cleanMitoMarkersReq = cleanMitoMarkersReq;
    this.lrrSdThresholdReq = lrrSdThresholdReq;
    this.callrateThresholdReq = callrateThresholdReq;
    this.qcPassingOnlyReq = qcPassingOnlyReq;
    this.imputeNaNs = imputeNaNs;
    this.recomputeLrrPCMarkersReq = recomputeLrrPCMarkersReq;
    this.recomputeLrrMedianMarkersReq = recomputeLrrMedianMarkersReq;
    this.homozygousOnlyReq = homozygousOnlyReq;
    this.gcRegressionDistanceReq = gcRegressionDistanceReq;
    this.pcSelectionSamplesReq = pcSelectionSamplesReq;
    this.externalBetaFileReq = externalBetaFileReq;
    this.numThreadsReq = numThreadsReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    double sampleLRRSdFilter = variables.get(lrrSdThresholdReq);
    if (sampleLRRSdFilter < 0) {
      switch (proj.ARRAY_TYPE.getValue()) {
        case AFFY_GW6:
        case AFFY_GW6_CN:
          proj.LRRSD_CUTOFF.setValue(0.35);
          proj.getLog()
              .reportTimeInfo("Setting " + proj.LRRSD_CUTOFF.getName()
                              + " to default 0.35 for array " + proj.ARRAY_TYPE.getValue());
          break;
        case ILLUMINA:
          proj.LRRSD_CUTOFF.setValue(0.30);
          proj.getLog()
              .reportTimeInfo("Setting " + proj.LRRSD_CUTOFF.getName()
                              + " to default 0.30 for array " + proj.ARRAY_TYPE.getValue());
          break;
        default:
          throw new IllegalArgumentException("Invalid Array type");
      }
    } else {
      proj.LRRSD_CUTOFF.setValue(sampleLRRSdFilter);
    }

    if (!Files.exists(proj.INTENSITY_PC_MARKERS_FILENAME.getValue())) {
      Files.writeIterable(IntensityMarkers.getIntensityMarkers(proj),
                          proj.INTENSITY_PC_MARKERS_FILENAME.getValue());
    }

  }

  @Override
  public void run(Variables variables) {
    String medianMarkers = null;
    if (variables.get(useAllMitoMarkersReq)) {
      medianMarkers = proj.PROJECT_DIRECTORY.getValue() + ALL_MITO_MARKERS_FILENAME;
      if (!Files.exists(medianMarkers)) {
        List<String> mitoNames = proj.getMarkerSet().getChrMap().get((byte) 26).stream()
                                     .map(Marker::getName).collect(Collectors.toList());
        Files.writeIterable(mitoNames, medianMarkers);
      }
    } else if (variables.hasValid(cleanMitoMarkersReq)) {
      medianMarkers = variables.get(cleanMitoMarkersReq).getPath();
    }
    double markerCallRateFilter = variables.get(callrateThresholdReq);
    // FIXME: This gcCorrect assignment was carried over from the old indexed version but
    // appears incorrect
    boolean gcCorrect = variables.get(qcPassingOnlyReq);
    boolean imputeMeanForNaN = variables.get(imputeNaNs);
    boolean recomputeLRRPCs = variables.get(recomputeLrrPCMarkersReq);
    boolean recomputeLRRMedian = variables.get(recomputeLrrMedianMarkersReq);
    boolean homozygousOnly = variables.get(homozygousOnlyReq);
    int bpGcModel = GcAdjustor.GcModel.DEFAULT_GC_MODEL_BIN_FASTA;
    int regressionDistance = variables.get(gcRegressionDistanceReq);
    int numComponents = MitoPipeline.DEFAULT_NUM_COMPONENTS;
    int numThreads = StepBuilder.resolveThreads(proj, variables.get(numThreadsReq));
    String outputBase = MitoPipeline.FILE_BASE;

    String betaOptFile = variables.get(pcSelectionSamplesReq).getAbsolutePath();
    String betaFile = variables.get(externalBetaFileReq).getAbsolutePath();

    boolean markerQC = true;
    double[] pvalOpt = MitoPipeline.DEFAULT_PVAL_OPTS;
    String pedFile = null;
    String useFile = null;
    boolean sampLrr = true;
    boolean plot = false;

    int retCode = PCAPrep.prepPCA(proj, numThreads, outputBase, markerQC, markerCallRateFilter,
                                  useFile, proj.getSampleList(), proj.getLog());
    if (retCode == PCAPrep.SUCCESS_CODE) {
      MitoPipeline.estimateMtDNACN(proj, numThreads, medianMarkers, numComponents, outputBase,
                                   homozygousOnly, markerCallRateFilter, betaOptFile, pedFile,
                                   recomputeLRRPCs, recomputeLRRMedian, sampLrr, imputeMeanForNaN,
                                   gcCorrect, bpGcModel, regressionDistance,
                                   proj.GENOME_BUILD_VERSION.getValue(), pvalOpt, betaFile, plot,
                                   false, PRE_PROCESSING_METHOD.NONE, proj.getLog());
    } else {
      throw new RuntimeException(PCAPrep.errorMessage(retCode));
    }
  }

  @Override
  public String getCommandLine(Variables variables) {
    return getStepCommandLine(proj, variables);
  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    String outputBase = proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.FILE_BASE;
    String finalReport = outputBase + PCA.FILE_EXTs[0];
    return Files.exists(finalReport);
  }

  public static void main(String[] args) {
    Project proj = Step.parseProject(args);
    StepBuilder sb = new StepBuilder(proj);
    Step markersParsingStep = sb.generateMarkersParsingStep();
    Step step = sb.generateMitoCNEstimateStep(markersParsingStep);
    Variables variables = step.parseArguments(args);
    Step.run(proj, step, variables);
  }

}
