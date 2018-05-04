package org.genvisis.cnv.workflow;

import java.io.File;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.StringJoiner;
import org.genvisis.CLI;
import org.genvisis.cnv.analysis.Mosaicism;
import org.genvisis.cnv.analysis.pca.PCA;
import org.genvisis.cnv.analysis.pca.PCAPrep;
import org.genvisis.cnv.analysis.pca.PCImputeRace;
import org.genvisis.cnv.analysis.pca.PCImputeRace.RACE;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsCompute.PRE_PROCESSING_METHOD;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity.CHROMOSOME_X_STRATEGY;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity.CORRECTION_TYPE;
import org.genvisis.cnv.filesys.ABLookup;
import org.genvisis.cnv.filesys.ABLookup.ABSource;
import org.genvisis.cnv.filesys.Centroids;
import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.hmm.CNVCaller;
import org.genvisis.cnv.hmm.CNVCaller.CALLING_SCOPE;
import org.genvisis.cnv.hmm.CNVCaller.PFB_MANAGEMENT_TYPE;
import org.genvisis.cnv.hmm.PFB;
import org.genvisis.cnv.manage.MitoPipeline;
import org.genvisis.cnv.manage.PRoCtOR;
import org.genvisis.cnv.manage.PlinkData;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.prop.Property;
import org.genvisis.cnv.qc.AffyMarkerBlast;
import org.genvisis.cnv.qc.GcAdjustor;
import org.genvisis.cnv.qc.IlluminaMarkerBlast;
import org.genvisis.cnv.qc.LrrSd;
import org.genvisis.cnv.qc.MarkerBlast;
import org.genvisis.cnv.qc.MarkerBlast.FILE_SEQUENCE_TYPE;
import org.genvisis.cnv.qc.MarkerMetrics;
import org.genvisis.cnv.qc.SampleQC;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.steps.SampleDataStep;
import org.genvisis.cnv.workflow.steps.SexChecksStep;
import org.genvisis.cnv.workflow.steps.TransposeStep;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.gwas.Ancestry;
import org.genvisis.gwas.FurtherAnalysisQc;
import org.genvisis.gwas.MarkerQC;
import org.genvisis.gwas.MarkerQC.QC_METRIC;
import org.genvisis.gwas.PlinkMendelianChecker;
import org.genvisis.gwas.Qc;
import org.genvisis.gwas.RelationAncestryQc;
import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

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

  private SortedSet<Step> buildSteps;
  private double p;
  private static final Requirement NUM_THREADS_REQ = new Requirement.PosIntRequirement(NUM_THREADS_DESC,
                                                                                       1);

  public StepBuilder() {
    buildSteps = Sets.newTreeSet();
    p = 0.0;
  }

  /**
   * @return All steps {@link #register(Step)}ed by this step builder thus far
   */
  public SortedSet<Step> getSortedSteps() {
    return buildSteps;
  }

  /**
   * @return The next step priority
   */
  double priority() {
    return ++p;
  }

  /**
   * Register the given step in the list returned by {@link #getSortedSteps()}
   */
  Step register(Step s) {
    buildSteps.add(s);
    return s;
  }

  Step generateIlluminaMarkerPositionsStep(Project proj) {
    final Requirement snpMapReq = new Requirement.FileRequirement("An Illumina SNP_map file.",
                                                                  proj.getLocationOfSNP_Map(false));
    final Requirement manifestReq = new Requirement.FileRequirement("An Illumina Manifest file.",
                                                                    proj.getLocationOfSNP_Map(false));

    final RequirementSet reqSet = RequirementSetBuilder.and()
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(snpMapReq)
                                                                                 .add(manifestReq));
    return register(new Step("Create Marker Positions (if not already exists)", "", reqSet,
                             EnumSet.noneOf(Requirement.Flag.class), priority()) {

      @Override
      public void setNecessaryPreRunProperties(Project proj,
                                               Map<Step, Map<Requirement, String>> variables) {
        // not needed for step
      }

      @Override
      public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
        proj.getLog().report("Generating marker positions file");
        String snpMap = variables.get(this).get(snpMapReq);
        String manifest = variables.get(this).get(manifestReq);
        if (Files.exists(snpMap)) {
          org.genvisis.cnv.manage.Markers.generateMarkerPositions(proj, snpMap);
        } else if (Files.exists(manifest)) {
          MarkerBlast.extractMarkerPositionsFromManifest(manifest, ARRAY.ILLUMINA,
                                                         FILE_SEQUENCE_TYPE.MANIFEST_FILE,
                                                         proj.MARKER_POSITION_FILENAME.getValue(false,
                                                                                                false),
                                                         Files.determineDelimiter(manifest,
                                                                                  proj.getLog()),
                                                         proj.getLog());
        }
      }

      @Override
      public boolean checkIfOutputExists(Project proj,
                                         Map<Step, Map<Requirement, String>> variables) {
        return Files.exists(proj.MARKER_POSITION_FILENAME.getValue(false, false));
      }

      @Override
      public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
        String projFile = proj.getPropertyFilename();
        String snpMap = variables.get(this).get(snpMapReq);
        String manifest = variables.get(this).get(manifestReq);
        String baseCommand = Files.getRunString() + " cnv.manage.Markers proj=" + projFile;
        if (Files.exists(snpMap)) {
          return baseCommand + " snps=" + snpMap;
        } else {
          return baseCommand + " snps=" + manifest + " -manifest";
        }
      }

    });
  }

  Step generateIlluminaMarkerBlastAnnotationStep(Project proj, final Step parseSamplesStep) {
    final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    final Requirement manifestFileReq = new Requirement.FileRequirement(ext.capitalizeFirst(IlluminaMarkerBlast.DESC_MANIFEST),
                                                                        IlluminaMarkerBlast.EXAMPLE_MANIFEST);

    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(manifestFileReq)
                                                       .add(GenvisisWorkflow.getNumThreadsReq());

    return register(new Step("Run Marker BLAST Annotation", "", reqSet,
                             EnumSet.of(Requirement.Flag.MEMORY, Requirement.Flag.RUNTIME,
                                        Requirement.Flag.MULTITHREADED),
                             priority()) {

      @Override
      public void setNecessaryPreRunProperties(Project proj,
                                               Map<Step, Map<Requirement, String>> variables) {
        // Not necessary for this step

      }

      @Override
      public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
        String manifestFile = variables.get(this).get(manifestFileReq);
        int numThreads = StepBuilder.resolveThreads(proj,
                                                    variables.get(this)
                                                             .get(GenvisisWorkflow.getNumThreadsReq()));
        new IlluminaMarkerBlast(proj, numThreads, manifestFile).blastEm();
      }

      @Override
      public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
        String manifestFile = variables.get(this).get(manifestFileReq);
        int numThreads = StepBuilder.resolveThreads(proj,
                                                    variables.get(this)
                                                             .get(GenvisisWorkflow.getNumThreadsReq()));
        ImmutableMap.Builder<String, String> argsBuilder = ImmutableMap.builder();
        argsBuilder.put(CLI.ARG_PROJ, proj.getPropertyFilename());
        argsBuilder.put(IlluminaMarkerBlast.ARG_MANIFEST, manifestFile);
        argsBuilder.put(CLI.ARG_THREADS, String.valueOf(numThreads));
        return Files.getRunString() + " "
               + CLI.formCmdLine(IlluminaMarkerBlast.class, argsBuilder.build());
      }

      @Override
      public boolean checkIfOutputExists(Project proj,
                                         Map<Step, Map<Requirement, String>> variables) {
        return Files.exists(proj.BLAST_ANNOTATION_FILENAME.getValue());
      }
    });
  }

  Step generateAffyMarkerBlastAnnotationStep(final Step parseSamplesStep) {
    final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);

    final Requirement probeFileReq = new Requirement.FileRequirement(ext.capitalizeFirst(AffyMarkerBlast.DESC_PROBE_FILE),
                                                                     AffyMarkerBlast.EXAMPLE_PROBE_FILE);
    final Requirement annotFileReq = new Requirement.FileRequirement(ext.capitalizeFirst(AffyMarkerBlast.DESC_ANNOT_FILE),
                                                                     AffyMarkerBlast.EXAMPLE_ANNOT_FILE);

    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(probeFileReq).add(annotFileReq)
                                                       .add(GenvisisWorkflow.getNumThreadsReq());

    return register(new Step("Run Marker BLAST Annotation", "", reqSet,
                             EnumSet.of(Requirement.Flag.MEMORY, Requirement.Flag.RUNTIME,
                                        Requirement.Flag.MULTITHREADED),
                             priority()) {

      @Override
      public void setNecessaryPreRunProperties(Project proj,
                                               Map<Step, Map<Requirement, String>> variables) {
        // Not necessary for this step

      }

      @Override
      public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
        String annotFile = variables.get(this).get(annotFileReq);
        String probeFile = variables.get(this).get(probeFileReq);
        int numThreads = StepBuilder.resolveThreads(proj,
                                                    variables.get(this)
                                                             .get(GenvisisWorkflow.getNumThreadsReq()));
        new AffyMarkerBlast(proj, numThreads, probeFile, annotFile).blastEm();
      }

      @Override
      public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
        String annotFile = variables.get(this).get(annotFileReq);
        String probeFile = variables.get(this).get(probeFileReq);
        int numThreads = StepBuilder.resolveThreads(proj,
                                                    variables.get(this)
                                                             .get(GenvisisWorkflow.getNumThreadsReq()));
        ImmutableMap.Builder<String, String> argsBuilder = ImmutableMap.builder();
        argsBuilder.put(CLI.ARG_PROJ, proj.getPropertyFilename());
        argsBuilder.put(AffyMarkerBlast.ARG_PROBE_FILE, probeFile);
        argsBuilder.put(AffyMarkerBlast.ARG_ANNOT_FILE, annotFile);
        argsBuilder.put(CLI.ARG_THREADS, String.valueOf(numThreads));
        return Files.getRunString() + " "
               + CLI.formCmdLine(AffyMarkerBlast.class, argsBuilder.build());
      }

      @Override
      public boolean checkIfOutputExists(Project proj,
                                         Map<Step, Map<Requirement, String>> variables) {
        return Files.exists(proj.BLAST_ANNOTATION_FILENAME.getValue());
      }
    });
  }

  Step generateParseSamplesStep(Project proj) {
    return generateParseSamplesStep(proj, null);
  }

  Step generateParseSamplesStep(Project proj, final Step markerPositionsStep) {

    final Requirement markerPositionsReq = new Requirement.FileRequirement("Marker Positions file must already exist.",
                                                                           proj.MARKER_POSITION_FILENAME.getValue(false,
                                                                                                                  false));

    final RequirementSet reqSet = RequirementSetBuilder.and();
    if (markerPositionsStep == null) {
      reqSet.add(markerPositionsReq).add(GenvisisWorkflow.getNumThreadsReq());
    } else {
      final Requirement markerPositionsStepReq = new Requirement.StepRequirement(markerPositionsStep);
      reqSet.add(RequirementSetBuilder.or().add(markerPositionsReq).add(markerPositionsStepReq))
            .add(GenvisisWorkflow.getNumThreadsReq());
    }

    return register(new Step("Parse Sample Files", "", reqSet,
                             EnumSet.of(Requirement.Flag.MEMORY, Requirement.Flag.RUNTIME,
                                        Requirement.Flag.MULTITHREADED),
                             priority()) {

      @Override
      public void setNecessaryPreRunProperties(Project proj,
                                               Map<Step, Map<Requirement, String>> variables) {
        String projFile = proj.MARKER_POSITION_FILENAME.getValue(false, false);
        String mkrFile = variables.get(this).get(markerPositionsReq);
        mkrFile = ext.verifyDirFormat(mkrFile);
        mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
        if (!mkrFile.equals(projFile)) {
          proj.MARKER_POSITION_FILENAME.setValue(mkrFile);
        }
        int numThreads = StepBuilder.resolveThreads(proj,
                                                    variables.get(this)
                                                             .get(GenvisisWorkflow.getNumThreadsReq()));
        GenvisisWorkflow.maybeSetProjNumThreads(proj, numThreads);
      }

      @Override
      public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
        int numThreads = proj.NUM_THREADS.getValue();
        proj.getLog().report("Parsing sample files");
        int retCode = org.genvisis.cnv.manage.SourceFileParser.createFiles(proj, numThreads);
        switch (retCode) {
          case 0:
            throw new RuntimeException("Operation failure, please check log for more information.");
          case 1:
          case 6:
          default:
            break;
        }
      }

      @Override
      public boolean checkIfOutputExists(Project proj,
                                         Map<Step, Map<Requirement, String>> variables) {
        String sampleDirectory = proj.SAMPLE_DIRECTORY.getValue(false, false);
        boolean mkrSetFile = Files.exists(proj.MARKERSET_FILENAME.getValue(false, false));
        boolean returnValue = mkrSetFile;
        returnValue = returnValue && proj.getSampleList() != null;
        returnValue = returnValue && Files.exists(sampleDirectory);

        int numSamples = returnValue ? proj.getSampleList().getSamples().length : 0;
        returnValue = returnValue && numSamples > 0;
        returnValue = returnValue && Files.countFiles(sampleDirectory,
                                                      Sample.SAMPLE_FILE_EXTENSION) == numSamples;
        // checking the validity / completeness of each sample would be a Good Thing, but too
        // costly time-wise for larger projects
        return returnValue;
      }

      @Override
      public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
        String projPropFile = proj.getPropertyFilename();
        StringBuilder kvCmd = new StringBuilder(Files.getRunString()).append(GenvisisWorkflow.PROJ_PROP_UPDATE_STR)
                                                                     .append(projPropFile);
        StringBuilder kvPairs = new StringBuilder();
        String projFile = proj.MARKER_POSITION_FILENAME.getValue(false, false);
        String mkrFile = variables.get(this).get(markerPositionsReq);
        mkrFile = ext.verifyDirFormat(mkrFile);
        mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
        if (!mkrFile.equals(projFile)) {
          kvPairs.append(" MARKER_POSITION_FILENAME=").append(mkrFile);
        }
        int numThreads = StepBuilder.resolveThreads(proj,
                                                    variables.get(this)
                                                             .get(GenvisisWorkflow.getNumThreadsReq()));
        if (numThreads != proj.NUM_THREADS.getValue()) {
          kvPairs.append(" ").append(proj.NUM_THREADS.getName()).append("=").append(numThreads);
        }
        StringBuilder command = new StringBuilder();
        if (kvPairs.length() != 0) {
          command.append(kvCmd).append(kvPairs).append("\n");
        }
        command.append(Files.getRunString()).append(" cnv.manage.SourceFileParser proj=")
               .append(projPropFile).append(" ").append(PSF.Ext.NUM_THREADS_COMMAND)
               .append(numThreads);
        return command.toString();
      }

    });
  }

  Step generateCreateSampleDataStep(Project proj, final Step parseSamplesStep) {
    return register(SampleDataStep.create(parseSamplesStep, proj, priority()));
  }

  Step generateTransposeStep(Project proj, final Step parseSamplesStep) {
    return register(TransposeStep.create(parseSamplesStep, priority()));
  }

  Step generateGCModelStep(Project proj) {
    final Requirement.ResourceRequirement gcBaseResourceReq = new Requirement.ResourceRequirement("GC Base file",
                                                                                                  Resources.genome(proj.GENOME_BUILD_VERSION.getValue(),
                                                                                                                   proj.getLog())
                                                                                                           .getModelBase());
    final Requirement gcModelOutputReq = new Requirement.OutputFileRequirement("GCModel output file must be specified.",
                                                                               proj.GC_MODEL_FILENAME.getValue());

    final RequirementSet reqSet = RequirementSetBuilder.and().add(gcBaseResourceReq)
                                                       .add(gcModelOutputReq);
    return register(new Step("Compute GCMODEL File", "", reqSet,
                             EnumSet.noneOf(Requirement.Flag.class), priority()) {

      @Override
      public void setNecessaryPreRunProperties(Project proj,
                                               Map<Step, Map<Requirement, String>> variables) {
        String setGCOutputFile = proj.GC_MODEL_FILENAME.getValue();
        String gcOutputFile = variables.get(this).get(gcModelOutputReq);
        if (!ext.verifyDirFormat(setGCOutputFile).equals(gcOutputFile)) {
          proj.GC_MODEL_FILENAME.setValue(gcOutputFile);
        }
      }

      @Override
      public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
        String gcBaseFile = gcBaseResourceReq.getResource().getAbsolute();
        String gcOutputFile = variables.get(this).get(gcModelOutputReq);
        org.genvisis.cnv.qc.GcAdjustor.GcModel.gcModel(proj, gcBaseFile, gcOutputFile, 100);
      }

      @Override
      public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
        String kvCmd = "";

        String setGCOutputFile = proj.GC_MODEL_FILENAME.getValue();
        String gcOutputFile = variables == null ? null : variables.get(this).get(gcModelOutputReq);
        if (gcOutputFile != null && !ext.verifyDirFormat(setGCOutputFile).equals(gcOutputFile)) {
          kvCmd += " GC_MODEL_FILENAME=" + gcOutputFile;
        }

        String projPropFile = proj.getPropertyFilename();
        StringBuilder cmd = new StringBuilder();
        if (kvCmd.length() > 0) {
          cmd.append(Files.getRunString())
             .append(GenvisisWorkflow.PROJ_PROP_UPDATE_STR + projPropFile).append(kvCmd)
             .append("\n");
        }
        String gcBaseFile = gcBaseResourceReq.getResource().getAbsolute();
        return cmd.append(Files.getRunString())
                  .append(" " + GcAdjustor.class.getName() + " " + CLI.ARG_PROJ + "="
                          + proj.getPropertyFilename() + " " + CLI.ARG_LOG + "="
                          + proj.getLog().getFilename() + " " + GcAdjustor.GC_BASE_FILE + "="
                          + gcBaseFile)
                  .toString();
      }

      @Override
      public boolean checkIfOutputExists(Project proj,
                                         Map<Step, Map<Requirement, String>> variables) {
        String gcOutputFile = variables.get(this).get(gcModelOutputReq);
        return Files.exists(gcOutputFile);
      }

    });
  }

  Step generateSampleQCStep(Project proj, final Step parseSamplesStep) {
    final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(GenvisisWorkflow.getNumThreadsReq());
    return register(new Step("Run Sample QC Metrics", "", reqSet,
                             EnumSet.of(Requirement.Flag.MULTITHREADED), priority()) {

      @Override
      public void setNecessaryPreRunProperties(Project proj,
                                               Map<Step, Map<Requirement, String>> variables) {
        int numThreads = StepBuilder.resolveThreads(proj,
                                                    variables.get(this)
                                                             .get(GenvisisWorkflow.getNumThreadsReq()));
        GenvisisWorkflow.maybeSetProjNumThreads(proj, numThreads);
      }

      @Override
      public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
        proj.getLog().report("Running LrrSd");
        int numThreads = proj.NUM_THREADS.getValue();
        LrrSd.init(proj, null, null, numThreads, false);
      }

      @Override
      public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
        int numThreads = StepBuilder.resolveThreads(proj,
                                                    variables.get(this)
                                                             .get(GenvisisWorkflow.getNumThreadsReq()));
        String projPropFile = proj.getPropertyFilename();
        StringBuilder cmd = new StringBuilder();
        cmd.append(Files.getRunString()).append(" cnv.qc.LrrSd").append(" proj=")
           .append(projPropFile).append(" ").append(PSF.Ext.NUM_THREADS_COMMAND).append(numThreads)
           .append(" projectMarkers=TRUE");
        return cmd.toString();
      }

      @Override
      public boolean checkIfOutputExists(Project proj,
                                         Map<Step, Map<Requirement, String>> variables) {
        return Files.exists(proj.SAMPLE_QC_FILENAME.getValue(false, false));
      }
    });
  }

  Step generateMarkerQCStep(Project proj, final Step parseSamplesStep) {
    final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);

    final Requirement exportAllReq = new Requirement.OptionalBoolRequirement("Export all markers in project.",
                                                                             true);

    String[] tgtMkrFiles = proj.TARGET_MARKERS_FILENAMES.getValue();
    final Requirement targetMarkersReq = new Requirement.FileRequirement("A targetMarkers files listing the markers to QC.",
                                                                         tgtMkrFiles != null && tgtMkrFiles.length >= 1 ? tgtMkrFiles[0]
                                                                                                                        : "");
    final Set<String> sampleDataHeaders;
    if (Files.exists(proj.SAMPLE_DATA_FILENAME.getValue()) && proj.getSampleData(false) != null) {
      sampleDataHeaders = proj.getSampleData(false).getMetaHeaders();
    } else {
      sampleDataHeaders = Sets.newHashSet();
    }
    final Set<String> defaultBatchHeaders = Sets.intersection(sampleDataHeaders,
                                                              MarkerMetrics.DEFAULT_SAMPLE_DATA_BATCH_HEADERS);
    final Requirement.ListSelectionRequirement batchHeadersReq = new Requirement.ListSelectionRequirement("SampleData column headers to use as batches for batch effects calculations",
                                                                                                          sampleDataHeaders,
                                                                                                          defaultBatchHeaders,
                                                                                                          true);

    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(exportAllReq)
                                                                                 .add(targetMarkersReq))
                                                       .add(batchHeadersReq)
                                                       .add(GenvisisWorkflow.getNumThreadsReq());

    return register(new Step("Run Marker QC Metrics", "", reqSet,
                             EnumSet.of(Requirement.Flag.MULTITHREADED), priority()) {

      @Override
      public void setNecessaryPreRunProperties(Project proj,
                                               Map<Step, Map<Requirement, String>> variables) {
        // Nothing to do here
      }

      @Override
      public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
        boolean allMarkers = Boolean.parseBoolean(variables.get(this).get(exportAllReq));
        String tgtFile = allMarkers ? null : variables.get(this).get(targetMarkersReq);
        boolean[] samplesToExclude = proj.getSamplesToExclude();
        int numThreads = StepBuilder.resolveThreads(proj,
                                                    variables.get(this)
                                                             .get(GenvisisWorkflow.getNumThreadsReq()));
        Set<String> batchHeaders = ImmutableSet.copyOf(Requirement.ListSelectionRequirement.parseArgValString(variables.get(this)
                                                                                                                       .get(batchHeadersReq)));
        MarkerMetrics.fullQC(proj, samplesToExclude, tgtFile, true, batchHeaders, numThreads);
      }

      @Override
      public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
        boolean allMarkers = Boolean.parseBoolean(variables.get(this).get(exportAllReq));
        String tgtFile = variables.get(this).get(targetMarkersReq);
        int numThreads = StepBuilder.resolveThreads(proj,
                                                    variables.get(this)
                                                             .get(GenvisisWorkflow.getNumThreadsReq()));
        List<String> batchHeaders = Requirement.ListSelectionRequirement.parseArgValString(variables.get(this)
                                                                                                    .get(batchHeadersReq));
        StringJoiner args = new StringJoiner(" ");
        args.add(Files.getRunString());
        args.add(MarkerMetrics.class.getCanonicalName());
        args.add("-fullQC");
        args.add("proj=" + proj.getPropertyFilename());
        if (!allMarkers) args.add("markers=" + tgtFile);
        String batchHeadersArg = String.join(",", batchHeaders);
        args.add(MarkerMetrics.BATCH_HEADERS_ARG + "=" + batchHeadersArg);
        args.add(PSF.Ext.NUM_THREADS_COMMAND + numThreads);
        return args.toString();
      }

      @Override
      public boolean checkIfOutputExists(Project proj,
                                         Map<Step, Map<Requirement, String>> variables) {
        String markerMetricsFile = proj.MARKER_METRICS_FILENAME.getValue();
        return Files.exists(markerMetricsFile);
      }

    });
  }

  Step generateSexChecksStep(Project proj, final Step parseSamplesStep, final Step markerBlastStep,
                             final Step sampleDataStep, final Step transposeStep,
                             final Step sampleQCStep) {
    return register(SexChecksStep.create(proj, parseSamplesStep, markerBlastStep, sampleDataStep,
                                         transposeStep, sampleQCStep, priority()));
  }

  Step generatePlinkExportStep(Project proj, final Step parseSamplesStep) {
    final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    final Requirement pedigreeRequirement = new Requirement.FileRequirement("A pedigree.dat file must exist.",
                                                                            proj.PEDIGREE_FILENAME.getValue(false,
                                                                                                            false));
    final Requirement createPedigreeRequirement = new Requirement.BoolRequirement("Create a minimal pedigree.dat file [will pull information from SexChecks step results].",
                                                                                  false);

    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(pedigreeRequirement)
                                                                                 .add(createPedigreeRequirement));

    return register(new Step("Create PLINK Files", "", reqSet, EnumSet.of(Requirement.Flag.MEMORY),
                             priority()) {

      @Override
      public Map<Requirement, String> getDefaultRequirementValues() {
        Map<Requirement, String> varMap = super.getDefaultRequirementValues();
        if (!Files.exists(proj.PEDIGREE_FILENAME.getValue())) {
          // if no pedigree, default to creating a minimal one
          varMap.put(createPedigreeRequirement, Boolean.TRUE.toString());
        }
        return varMap;
      }

      @Override
      public void setNecessaryPreRunProperties(Project proj,
                                               Map<Step, Map<Requirement, String>> variables) {
        if (!Boolean.parseBoolean(variables.get(this).get(createPedigreeRequirement))) {
          String projPedFile = proj.PEDIGREE_FILENAME.getValue(false, false);
          String pedFile = variables.get(this).get(pedigreeRequirement);
          if (!pedFile.equals(projPedFile)) {
            proj.PEDIGREE_FILENAME.setValue(pedFile);
          }
        }
      }

      @Override
      public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
        if (Boolean.parseBoolean(variables.get(this).get(createPedigreeRequirement))) {
          proj.getLog().report("Creating Pedigree File");
          Pedigree.build(proj, null, null, false);
        }
        if (!Files.exists(proj.PEDIGREE_FILENAME.getValue())) {
          throw new RuntimeException("Creation of Pedigree file in [Create/Run PLINK Files] step failed.");
        }

        proj.getLog().report("Running PLINK");

        boolean create = PlinkData.saveGenvisisToPlinkBedSet(proj,
                                                             GenvisisWorkflow.PLINK_SUBDIR
                                                                   + GenvisisWorkflow.PLINKROOT,
                                                             null, null,
                                                             PlinkData.ExportIDScheme.DNA_DNA);
        if (!create) {
          throw new RuntimeException("Creation of initial PLINK files failed.");
        }
        proj.PLINK_DIR_FILEROOTS.addValue(proj.PROJECT_DIRECTORY.getValue()
                                          + GenvisisWorkflow.PLINK_SUBDIR
                                          + GenvisisWorkflow.PLINKROOT);
      }

      @Override
      public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
        String kvCmd = "";

        if (!Boolean.parseBoolean(variables.get(this).get(createPedigreeRequirement))) {
          String projPedFile = proj.PEDIGREE_FILENAME.getValue(false, false);
          String pedFile = variables.get(this).get(pedigreeRequirement);
          if (!pedFile.equals(projPedFile)) {
            kvCmd += " PEDIGREE_FILENAME=" + pedFile;
          }
        }

        String projPropFile = proj.getPropertyFilename();
        StringBuilder cmd = new StringBuilder();
        if (kvCmd.length() > 0) {
          cmd.append(Files.getRunString()).append(GenvisisWorkflow.PROJ_PROP_UPDATE_STR)
             .append(projPropFile).append(kvCmd).append("\n");
        }
        if (Boolean.parseBoolean(variables.get(this).get(createPedigreeRequirement))) {
          cmd.append(Files.getRunString()).append(" cnv.filesys.Pedigree proj=")
             .append(projPropFile).append("\n");
        }
        cmd.append(new StringJoiner(" ").add(Files.getRunString()).add(PlinkData.class.getName())
                                        .add("-genvisisToBed")
                                        .add("plinkdata=" + GenvisisWorkflow.PLINK_SUBDIR
                                             + GenvisisWorkflow.PLINKROOT)
                                        .add("proj=" + proj.getPropertyFilename())
                                        .add(PlinkData.ARG_EXPORT_ID_SCHEME
                                             + PlinkData.ExportIDScheme.DNA_DNA));
        return cmd.toString();
      }

      @Override
      public boolean checkIfOutputExists(Project proj,
                                         Map<Step, Map<Requirement, String>> variables) {
        boolean plinkFilesExist = Files.checkAllFiles(GenvisisWorkflow.getPlinkDir(proj),
                                                      PSF.Plink.getPlinkBedBimFamSet(GenvisisWorkflow.PLINKROOT),
                                                      false, proj.getLog());
        boolean pedGenerated = Boolean.parseBoolean(variables.get(this)
                                                             .get(createPedigreeRequirement));
        boolean pedCheck = pedGenerated ? Files.exists(proj.PEDIGREE_FILENAME.getValue()) : true;
        return plinkFilesExist && pedCheck;
      }

    });
  }

  Step generateGwasQCStep(Project proj, Step plinkExportStep) {
    final Requirement plinkExportStepReq = new Requirement.StepRequirement(plinkExportStep);
    String defaultCallrate;
    switch (proj.getArrayType()) {
      case AFFY_GW6:
      case AFFY_GW6_CN:
        defaultCallrate = MarkerQC.DEFAULT_AFFY_CALLRATE_THRESHOLD;
        break;
      case AFFY_AXIOM:
      case ILLUMINA:
        defaultCallrate = MarkerQC.DEFAULT_ILLUMINA_CALLRATE_THRESHOLD;
        break;
      default:
        throw new IllegalArgumentException("Invalid " + proj.getArrayType().getClass().getName()
                                           + ": " + proj.getArrayType().toString());
    }
    final Requirement callrateReq = new Requirement.ThresholdRequirement(QC_METRIC.CALLRATE.getUserDescription(),
                                                                         defaultCallrate);
    final RequirementSet reqSet = RequirementSetBuilder.and().add(plinkExportStepReq)
                                                       .add(callrateReq);

    return register(new Step("Run GWAS QC", "", reqSet, EnumSet.noneOf(Requirement.Flag.class),
                             priority()) {

      @Override
      public void setNecessaryPreRunProperties(Project proj,
                                               Map<Step, Map<Requirement, String>> variables) {
        // not needed for step
      }

      @Override
      public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
        String dir = GenvisisWorkflow.getPlinkDir(proj);
        Map<QC_METRIC, String> markerQCThresholds = Maps.newEnumMap(RelationAncestryQc.DEFAULT_QC_METRIC_THRESHOLDS);
        markerQCThresholds.put(QC_METRIC.CALLRATE, variables.get(this).get(callrateReq));
        new RelationAncestryQc(dir, GenvisisWorkflow.PLINKROOT, markerQCThresholds,
                               proj.getLog()).run(false);
        if (new File(dir + Qc.QC_SUBDIR + RelationAncestryQc.GENOME_DIR + GenvisisWorkflow.PLINKROOT
                     + ".genome").exists()) {
          proj.GENOME_CLUSTER_FILENAME.setValue(dir + Qc.QC_SUBDIR + RelationAncestryQc.GENOME_DIR
                                                + GenvisisWorkflow.PLINKROOT + ".genome");
          proj.saveProperties();
        }
        new PlinkMendelianChecker(proj).run();
      }

      @Override
      public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
        Map<Requirement, String> stepVars = variables.get(this);

        String dir = GenvisisWorkflow.getPlinkDir(proj);
        Map<QC_METRIC, String> markerQCThresholds = Maps.newEnumMap(RelationAncestryQc.DEFAULT_QC_METRIC_THRESHOLDS);
        markerQCThresholds.put(QC_METRIC.CALLRATE, variables.get(this).get(callrateReq));

        List<String> commandChunks = Lists.newArrayList();
        commandChunks.add(Files.getRunString());
        commandChunks.add(RelationAncestryQc.class.getName());
        commandChunks.add(CLI.formCmdLineArg(CLI.ARG_INDIR, GenvisisWorkflow.getPlinkDir(proj)));
        commandChunks.add(CLI.formCmdLineArg(CLI.ARG_PLINKROOT, GenvisisWorkflow.PLINKROOT));
        commandChunks.add(CLI.formCmdLineArg(RelationAncestryQc.ARGS_KEEPGENOME, "false"));
        commandChunks.add(CLI.formCmdLineArg(QC_METRIC.CALLRATE.getKey(),
                                             stepVars.get(callrateReq)));
        commandChunks.add("\n" + Files.getRunString());
        commandChunks.add(GenvisisWorkflow.PROJ_PROP_UPDATE_STR + proj.getPropertyFilename());
        commandChunks.add(proj.GENOME_CLUSTER_FILENAME.getName() + "=" + dir + Qc.QC_SUBDIR
                          + RelationAncestryQc.GENOME_DIR + GenvisisWorkflow.PLINKROOT + ".genome");
        commandChunks.add("\n" + Files.getRunString());
        commandChunks.add(PlinkMendelianChecker.class.getName());
        commandChunks.add("proj=" + proj.getPropertyFilename());
        return Joiner.on(' ').join(commandChunks);
      }

      @Override
      public boolean checkIfOutputExists(Project proj,
                                         Map<Step, Map<Requirement, String>> variables) {
        String dir = GenvisisWorkflow.getPlinkDir(proj);
        for (int i = 0; i < org.genvisis.gwas.RelationAncestryQc.FOLDERS_CREATED.length; i++) {
          for (int j = 0; j < org.genvisis.gwas.RelationAncestryQc.FILES_CREATED[i].length; j++) {
            if (!Files.exists(dir + org.genvisis.gwas.RelationAncestryQc.FOLDERS_CREATED[i]
                              + org.genvisis.gwas.RelationAncestryQc.FILES_CREATED[i][j])) {
              return false;
            }
          }
        }
        return Files.checkAllFiles(PlinkMendelianChecker.parseOutputDirectory(proj),
                                   PlinkMendelianChecker.OUTPUTS, false, proj.getLog());
      }
    });
  }

  Step generateAncestryStep(Project proj, final Step gwasQCStep) {
    final Requirement gwasQCStepReq = new Requirement.StepRequirement(gwasQCStep);
    final Requirement putativeWhitesReq = new Requirement.FileRequirement("File with FID/IID pairs of putative white samples",
                                                                          "");
    final Requirement.ResourceRequirement hapMapFoundersReq = new Requirement.ResourceRequirement("PLINK root of HapMap founders",
                                                                                                  Resources.hapMap(proj.getLog())
                                                                                                           .getUnambiguousHapMapFounders());

    final RequirementSet reqSet = RequirementSetBuilder.and().add(gwasQCStepReq)
                                                       .add(putativeWhitesReq)
                                                       .add(hapMapFoundersReq);
    final Requirement.ResourceRequirement hapMapAncestryReq = new Requirement.ResourceRequirement("HapMap Samples Ancestry File",
                                                                                                  Resources.hapMap(proj.getLog())
                                                                                                           .getHapMapAncestries());

    return register(new Step("Run Ancestry Checks", "", reqSet,
                             EnumSet.noneOf(Requirement.Flag.class), priority()) {

      @Override
      public void setNecessaryPreRunProperties(Project proj,
                                               Map<Step, Map<Requirement, String>> variables) {
        // not needed for step
      }

      @Override
      public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
        String putativeWhites = variables.get(this).get(putativeWhitesReq);
        String hapMapPlinkRoot = hapMapFoundersReq.getResource().getAbsolute();
        hapMapAncestryReq.getResource().get();
        String ancestryDir = GenvisisWorkflow.getAncestryDir(proj);
        Ancestry.runPipeline(ancestryDir, putativeWhites, hapMapPlinkRoot, proj,
                             new Logger(ancestryDir + "ancestry.log"));
      }

      @Override
      public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
        String putativeWhites = variables.get(this).get(putativeWhitesReq);
        String hapMapPlinkRoot = hapMapFoundersReq.getResource().getAbsolute();
        hapMapAncestryReq.getResource().get();
        String ancestryDir = GenvisisWorkflow.getAncestryDir(proj);
        String command = Files.getRunString() + " gwas.Ancestry -runPipeline dir=" + ancestryDir;
        command += " putativeWhites=" + putativeWhites;
        command += " proj=" + proj.getPropertyFilename();
        command += " hapMapPlinkRoot=" + hapMapPlinkRoot;
        command += " log=" + ancestryDir + "ancestry.log";
        return command;
      }

      @Override
      public boolean checkIfOutputExists(Project proj,
                                         Map<Step, Map<Requirement, String>> variables) {
        String ancestryDir = GenvisisWorkflow.getAncestryDir(proj);
        return Files.exists(ancestryDir + Ancestry.RACE_FREQS_FILENAME)
               && Files.exists(ancestryDir + Ancestry.RACE_IMPUTATIONAS_FILENAME);
      }
    });
  }

  Step generateFurtherAnalysisQCStep(Project proj, Step plinkExportStep, Step gwasQCStep,
                                     Step ancestryStep) {
    final Requirement plinkExportStepReq = new Requirement.StepRequirement(plinkExportStep);
    final Requirement gwasQCStepReq = new Requirement.StepRequirement(gwasQCStep);
    final Requirement ancestryStepReq = new Requirement.StepRequirement(ancestryStep);
    final Requirement unrelatedsFileReq = new Requirement.FileRequirement("File with list of unrelated FID/IID pairs to use for marker QC",
                                                                          "");
    final Requirement europeansFilesReq = new Requirement.FileRequirement("File with list of European samples to use for Hardy-Weinberg equilibrium tests",
                                                                          "");

    final RequirementSet reqSet = RequirementSetBuilder.and().add(plinkExportStepReq)
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(gwasQCStepReq)
                                                                                 .add(unrelatedsFileReq))
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(ancestryStepReq)
                                                                                 .add(europeansFilesReq));
    final Map<QC_METRIC, Requirement> metricRequirements = Maps.newEnumMap(QC_METRIC.class);
    for (QC_METRIC metric : QC_METRIC.values()) {
      Map<QC_METRIC, String> defaultThresholds = FurtherAnalysisQc.getDefaultMarkerQCThresholds(proj.getArrayType());
      String defaultVal = defaultThresholds.get(metric);
      final Requirement metricReq = new Requirement.ThresholdRequirement(metric.getUserDescription(),
                                                                         defaultVal);
      reqSet.add(metricReq);
      metricRequirements.put(metric, metricReq);
    }

    return register(new Step("Run Further Analysis QC", "", reqSet,
                             EnumSet.noneOf(Requirement.Flag.class), priority()) {

      @Override
      public void setNecessaryPreRunProperties(Project proj,
                                               Map<Step, Map<Requirement, String>> variables) {
        // not needed for step
      }

      @Override
      public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
        Map<Requirement, String> stepVars = variables.get(this);

        String unrelatedsFile = resolveUnrelatedsFile(stepVars);

        String europeansFile = resolveEuropeansFile(stepVars);

        Map<QC_METRIC, String> markerQCThresholds = Maps.newEnumMap(QC_METRIC.class);
        for (QC_METRIC metric : QC_METRIC.values()) {
          Requirement req = metricRequirements.get(metric);
          markerQCThresholds.put(metric, stepVars.get(req));
        }
        new FurtherAnalysisQc(GenvisisWorkflow.getPlinkDir(proj), GenvisisWorkflow.PLINKROOT,
                              markerQCThresholds, unrelatedsFile, europeansFile, proj.getLog())
                                                                                               .runFurtherAnalysisQC();
      }

      @Override
      public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
        Map<Requirement, String> stepVars = variables.get(this);

        String unrelatedsFile = resolveUnrelatedsFile(stepVars);

        String europeansFile = resolveEuropeansFile(stepVars);

        List<String> commandChunks = Lists.newArrayList();
        commandChunks.add(Files.getRunString());
        commandChunks.add(FurtherAnalysisQc.class.getName());
        commandChunks.add(CLI.formCmdLineArg(FurtherAnalysisQc.ARG_UNRELATEDS, unrelatedsFile));
        commandChunks.add(CLI.formCmdLineArg(FurtherAnalysisQc.ARG_EUROPEANS, europeansFile));
        commandChunks.add(CLI.formCmdLineArg(CLI.ARG_INDIR, GenvisisWorkflow.getPlinkDir(proj)));
        commandChunks.add(CLI.formCmdLineArg(CLI.ARG_PLINKROOT, GenvisisWorkflow.PLINKROOT));
        for (QC_METRIC metric : QC_METRIC.values()) {
          Requirement req = metricRequirements.get(metric);
          commandChunks.add(CLI.formCmdLineArg(metric.getKey(), stepVars.get(req)));
        }
        return Joiner.on(' ').join(commandChunks);
      }

      @Override
      public boolean checkIfOutputExists(Project proj,
                                         Map<Step, Map<Requirement, String>> variables) {
        String dir = GenvisisWorkflow.getPlinkDir(proj) + Qc.QC_SUBDIR
                     + FurtherAnalysisQc.FURTHER_ANALYSIS_DIR;
        String qcdPlinkroot = GenvisisWorkflow.PLINKROOT
                              + FurtherAnalysisQc.FURTHER_ANALYSIS_QC_PLINK_SUFFIX;
        return PSF.Plink.bedBimFamExist(dir + qcdPlinkroot)
               && Files.exists(dir + FurtherAnalysisQc.SAMPLE_QC_DROPS, false)
               && Files.exists(dir + FurtherAnalysisQc.MARKER_QC_DROPS, false);
      }

      private String resolveUnrelatedsFile(Map<Requirement, String> stepVars) {
        String unrelatedsFile = stepVars.get(unrelatedsFileReq);
        if (!Files.exists(unrelatedsFile)) {
          unrelatedsFile = GenvisisWorkflow.getAncestryDir(proj)
                           + RelationAncestryQc.UNRELATEDS_FILENAME;
        }
        return unrelatedsFile;
      }

      private String resolveEuropeansFile(Map<Requirement, String> stepVars) {
        String europeansFile = stepVars.get(europeansFilesReq);
        if (europeansFile == null || "".equals(europeansFile)) {
          String raceImputationFilename = GenvisisWorkflow.getAncestryDir(proj)
                                          + Ancestry.RACE_IMPUTATIONAS_FILENAME;
          europeansFile = PCImputeRace.formRaceListFilename(RACE.WHITE, raceImputationFilename);
        }
        return europeansFile;
      }

    });
  }

  Step generateMosaicArmsStep(Project proj, final Step parseSamplesStep) {
    final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(GenvisisWorkflow.getNumThreadsReq());
    return register(new Step("Create Mosaic Arms File", "", reqSet,
                             EnumSet.of(Requirement.Flag.MULTITHREADED), priority()) {

      @Override
      public void setNecessaryPreRunProperties(Project proj,
                                               Map<Step, Map<Requirement, String>> variables) {

        int numThreads = StepBuilder.resolveThreads(proj,
                                                    variables.get(this)
                                                             .get(GenvisisWorkflow.getNumThreadsReq()));
        GenvisisWorkflow.maybeSetProjNumThreads(proj, numThreads);
      }

      @Override
      public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
        Mosaicism.findOutliers(proj);
      }

      @Override
      public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
        String kvCmd = "";

        int numThreads = StepBuilder.resolveThreads(proj,
                                                    variables.get(this)
                                                             .get(GenvisisWorkflow.getNumThreadsReq()));
        if (numThreads != proj.NUM_THREADS.getValue()) {
          kvCmd += " " + proj.NUM_THREADS.getName() + "=" + numThreads;
        }

        String projPropFile = proj.getPropertyFilename();
        StringBuilder cmd = new StringBuilder();
        if (kvCmd.length() > 0) {
          cmd.append(Files.getRunString())
             .append(GenvisisWorkflow.PROJ_PROP_UPDATE_STR + projPropFile).append(kvCmd)
             .append("\n");
        }
        return cmd.append(Files.getRunString())
                  .append(" cnv.analysis.Mosaicism proj=" + proj.getPropertyFilename()).toString();
      }

      @Override
      public boolean checkIfOutputExists(Project proj,
                                         Map<Step, Map<Requirement, String>> variables) {
        return Files.exists(proj.RESULTS_DIRECTORY.getValue(false, false) + "Mosaicism.xln");
      }
    });
  }

  Step generateAnnotateSampleDataStep(Project proj, final Step sampleQCStep,
                                      final Step createSampleDataStep, final Step gwasQCStep) {
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

    return register(new Step("Annotate Sample Data File", "", reqSet,
                             EnumSet.noneOf(Requirement.Flag.class), priority()) {

      @Override
      public void setNecessaryPreRunProperties(Project proj,
                                               Map<Step, Map<Requirement, String>> variables) {
        double projLrrSdThreshold = proj.LRRSD_CUTOFF.getValue();
        double lrrSdThreshold = Double.parseDouble(variables.get(this).get(lrrSdThresholdReq));
        double projCallrateThreshold = proj.SAMPLE_CALLRATE_THRESHOLD.getValue();
        double callrateThreshold = Double.parseDouble(variables.get(this)
                                                               .get(callrateThresholdReq));

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
        double callrateThreshold = Double.parseDouble(variables.get(this)
                                                               .get(callrateThresholdReq));

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
          cmd.append(Files.getRunString())
             .append(GenvisisWorkflow.PROJ_PROP_UPDATE_STR + projPropFile).append(kvCmd)
             .append("\n");
        }
        cmd.append(Files.getRunString())
           .append(" cnv.qc.SampleQC proj=" + projPropFile + " numQ=" + numQ
                   + " justQuantiles=false" + " gcCorrectedLrrSd=" + gcCorrectedLrrSd
                   + " duplicatesSetFile=" + duplicatesSetFile + " correctFidIids="
                   + correctFidIids);
        return cmd.toString();
      }

      @Override
      public boolean checkIfOutputExists(Project proj,
                                         Map<Step, Map<Requirement, String>> variables) {
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

    });
  }

  Step generateMitoCNEstimateStep(Project proj, Step transposeStep) {
    // FIXME http://genvisis.org/MitoPipeline/#illumina_marker_lists has illumina markers.. this
    // should be linked to, or
    // these steps split or something...
    final Requirement transposeStepReq = new Requirement.StepRequirement(transposeStep);
    final Requirement medianMarkersReq = new Requirement.FileRequirement("MedianMarkers file must exist.",
                                                                         "");
    final Requirement lrrSdThresholdReq = new Requirement.DoubleRequirement("LRR SD threshold to filter samples.",
                                                                            proj.LRRSD_CUTOFF.getValue(),
                                                                            proj.LRRSD_CUTOFF.getMinValue(),
                                                                            proj.LRRSD_CUTOFF.getMaxValue());
    final Requirement callrateThresholdReq = new Requirement.DoubleRequirement("Call rate threshold to filter markers.",
                                                                               MitoPipeline.DEFAULT_MKR_CALLRATE_FILTER,
                                                                               0.0, 1.0);
    final Requirement qcPassingOnlyReq = new Requirement.OptionalBoolRequirement("Compute PCs with samples passing QC only",
                                                                                 true);
    final Requirement imputeNaNs = new Requirement.OptionalBoolRequirement("Impute mean value for NaN",
                                                                           true);
    final Requirement recomputeLrrPCMarkersReq = new Requirement.OptionalBoolRequirement("Should recompute Log-R ratio for PC markers?",
                                                                                         true);
    final Requirement recomputeLrrMedianMarkersReq = new Requirement.OptionalBoolRequirement("Should recompute Log-R ratio for median markers?",
                                                                                             true);
    final Requirement homozygousOnlyReq = new Requirement.OptionalBoolRequirement("Homozygous only?",
                                                                                  true);
    final Requirement gcRegressionDistanceReq = new Requirement.PosIntRequirement("Regression distance for the GC adjustment",
                                                                                  GcAdjustor.DEFAULT_REGRESSION_DISTANCE[0]);
    final Requirement pcSelectionSamplesReq = new Requirement.OptionalFileRequirement("A file listing a subset of samples (DNA ID) to use for determining optimal PC selection, typically a list of unrelated and single race samples. If a list is not provided, only samples passing sample qc thresholds will be used.",
                                                                                      "");
    final Requirement externalBetaFileReq = new Requirement.OptionalFileRequirement("An external beta file to optimize PC selection.",
                                                                                    "");

    final RequirementSet reqSet = RequirementSetBuilder.and().add(transposeStepReq)
                                                       .add(medianMarkersReq).add(lrrSdThresholdReq)
                                                       .add(callrateThresholdReq)
                                                       .add(qcPassingOnlyReq).add(imputeNaNs)
                                                       .add(recomputeLrrPCMarkersReq)
                                                       .add(recomputeLrrMedianMarkersReq)
                                                       .add(homozygousOnlyReq)
                                                       .add(gcRegressionDistanceReq)
                                                       .add(GenvisisWorkflow.getNumThreadsReq())
                                                       .add(pcSelectionSamplesReq)
                                                       .add(externalBetaFileReq);

    return register(new Step("Create Mitochondrial Copy-Number Estimates File", "", reqSet,
                             EnumSet.of(Requirement.Flag.MULTITHREADED), priority()) {

      @Override
      public void setNecessaryPreRunProperties(Project proj,
                                               Map<Step, Map<Requirement, String>> variables) {
        double sampleLRRSdFilter = Double.parseDouble(variables.get(this).get(lrrSdThresholdReq));
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
      }

      @Override
      public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
        String medianMarkers = variables.get(this).get(medianMarkersReq);
        double markerCallRateFilter = Double.parseDouble(variables.get(this)
                                                                  .get(callrateThresholdReq));
        // FIXME: This gcCorrect assignment was carried over from the old indexed version but
        // appears incorrect
        boolean gcCorrect = Boolean.parseBoolean(variables.get(this).get(qcPassingOnlyReq));
        boolean imputeMeanForNaN = Boolean.parseBoolean(variables.get(this).get(imputeNaNs));
        boolean recomputeLRRPCs = Boolean.parseBoolean(variables.get(this)
                                                                .get(recomputeLrrPCMarkersReq));
        boolean recomputeLRRMedian = Boolean.parseBoolean(variables.get(this)
                                                                   .get(recomputeLrrMedianMarkersReq));
        boolean homozygousOnly = Boolean.parseBoolean(variables.get(this).get(homozygousOnlyReq));
        int bpGcModel = GcAdjustor.GcModel.DEFAULT_GC_MODEL_BIN_FASTA;
        int regressionDistance = Integer.parseInt(variables.get(this).get(gcRegressionDistanceReq));
        int numComponents = MitoPipeline.DEFAULT_NUM_COMPONENTS;
        int numThreads = StepBuilder.resolveThreads(proj,
                                                    variables.get(this)
                                                             .get(GenvisisWorkflow.getNumThreadsReq()));
        String outputBase = MitoPipeline.FILE_BASE;

        String betaOptFile = variables.get(this).get(pcSelectionSamplesReq);
        String betaFile = variables.get(this).get(externalBetaFileReq);

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
                                       recomputeLRRPCs, recomputeLRRMedian, sampLrr,
                                       imputeMeanForNaN, gcCorrect, bpGcModel, regressionDistance,
                                       proj.GENOME_BUILD_VERSION.getValue(), pvalOpt, betaFile,
                                       plot, false, PRE_PROCESSING_METHOD.NONE, proj.getLog());
        } else {
          throw new RuntimeException(PCAPrep.errorMessage(retCode));
        }
      }

      @Override
      public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
        String medianMarkers = variables.get(this).get(medianMarkersReq);
        double lrrSD = Double.parseDouble(variables.get(this).get(lrrSdThresholdReq));
        double markerCallRateFilter = Double.parseDouble(variables.get(this)
                                                                  .get(callrateThresholdReq));
        // FIXME: This gcCorrect assignment was carried over from the old indexed version but
        // appears incorrect
        boolean gcCorrect = Boolean.parseBoolean(variables.get(this).get(qcPassingOnlyReq));
        boolean imputeMeanForNaN = Boolean.parseBoolean(variables.get(this).get(imputeNaNs));
        boolean recomputeLRRPCs = Boolean.parseBoolean(variables.get(this)
                                                                .get(recomputeLrrPCMarkersReq));
        boolean recomputeLRRMedian = Boolean.parseBoolean(variables.get(this)
                                                                   .get(recomputeLrrMedianMarkersReq));
        boolean homozygousOnly = Boolean.parseBoolean(variables.get(this).get(homozygousOnlyReq));
        int bpGcModel = GcAdjustor.GcModel.DEFAULT_GC_MODEL_BIN_FASTA;
        int regressionDistance = Integer.parseInt(variables.get(this).get(gcRegressionDistanceReq));
        int numComponents = MitoPipeline.DEFAULT_NUM_COMPONENTS;
        int numThreads = StepBuilder.resolveThreads(proj,
                                                    variables.get(this)
                                                             .get(GenvisisWorkflow.getNumThreadsReq()));
        String outputBase = MitoPipeline.FILE_BASE;

        String betaOptFile = variables.get(this).get(pcSelectionSamplesReq);
        String betaFile = variables.get(this).get(externalBetaFileReq);
        boolean sampLrr = true;

        String projPropFile = proj.getPropertyFilename();
        StringBuilder cmd = new StringBuilder();
        cmd.append(Files.getRunString()).append(" org.genvisis.cnv.manage.MitoPipeline")
           .append(" proj=").append(projPropFile).append(" mitochondrialMarkers=")
           .append(medianMarkers).append(" numComponents=").append(numComponents)
           .append(" imputeMeanForNaN=").append(imputeMeanForNaN).append(" recomputeLRR_PCs=")
           .append(recomputeLRRPCs).append(" recomputeLRR_Median=").append(recomputeLRRMedian)
           .append(" gcCorrect=").append(gcCorrect).append(" bpGcModel=").append(bpGcModel)
           .append(" LRRSD=").append(lrrSD).append(" markerCallRate=").append(markerCallRateFilter)
           .append(" regressionDistance=").append(regressionDistance).append(" sampLRR=")
           .append(sampLrr).append(" ").append(PSF.Ext.NUM_THREADS_COMMAND).append(numThreads)
           .append(" log=").append(proj.getLog().getFilename()).append(" output=")
           .append(outputBase);
        if (!"".equals(betaOptFile)) {
          cmd.append(" ").append(MitoPipeline.PC_OPT_FILE).append("=").append(betaOptFile);
        }
        if (!"".equals(betaFile)) {
          cmd.append(" betas=").append(betaFile);
        }
        if (!homozygousOnly) {
          cmd.append(" -allCalls ");
        }

        cmd.append(" -SkipProjectCreationWithLongUndocumentedFlag ");

        return cmd.toString();
      }

      @Override
      public boolean checkIfOutputExists(Project proj,
                                         Map<Step, Map<Requirement, String>> variables) {
        String outputBase = proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.FILE_BASE;
        String finalReport = outputBase + PCA.FILE_EXTs[0];// PrincipalComponentsResiduals.MT_REPORT_EXT[0];
        // boolean mkrFiles = true;
        // for (String file :
        // PrincipalComponentsResiduals.MT_REPORT_MARKERS_USED) {
        // if (!Files.exists(outputBase + file)) {
        // mkrFiles = false;
        // break;
        // }
        // }
        return Files.exists(finalReport) /* && mkrFiles */;
      }
    });
  }

  Step generatePFBStep(Project proj, final Step parseSamplesStep) {
    final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    final Requirement sampleSubsetReq = new Requirement.FileRequirement("A Sample subset file must exist.",
                                                                        proj.SAMPLE_SUBSET_FILENAME.getValue());
    String defaultOutputFile;
    if (Files.exists(proj.SAMPLE_SUBSET_FILENAME.getValue())) {
      defaultOutputFile = ext.rootOf(proj.SAMPLE_SUBSET_FILENAME.getValue()) + ".pfb";
    } else {
      defaultOutputFile = proj.CUSTOM_PFB_FILENAME.getValue();
    }
    final Requirement outputFileReq = new Requirement.OutputFileRequirement("PFB (population BAF) output file must be specified.",
                                                                            defaultOutputFile);

    final RequirementSet reqSet = RequirementSetBuilder.and()
                                                       .add(RequirementSetBuilder.or()
                                                                                 .add(parseSamplesStepReq)
                                                                                 .add(sampleSubsetReq))
                                                       .add(outputFileReq);

    return register(new Step("Compute Population BAF files", "", reqSet,
                             EnumSet.noneOf(Requirement.Flag.class), priority()) {

      @Override
      public void setNecessaryPreRunProperties(Project proj,
                                               Map<Step, Map<Requirement, String>> variables) {
        String setSubSampFile = proj.SAMPLE_SUBSET_FILENAME.getValue();
        String subSampFile = variables.get(this).get(sampleSubsetReq);
        String setPFBFile = proj.CUSTOM_PFB_FILENAME.getValue();
        String pfbOutputFile = variables.get(this).get(outputFileReq);

        if (!ext.verifyDirFormat(setSubSampFile).equals(subSampFile)) {
          proj.SAMPLE_SUBSET_FILENAME.setValue(subSampFile);
        }
        if (!ext.verifyDirFormat(setPFBFile).equals(pfbOutputFile)) {
          proj.CUSTOM_PFB_FILENAME.setValue(pfbOutputFile);
        }
      }

      @Override
      public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
        org.genvisis.cnv.hmm.PFB.populationBAF(proj);
      }

      @Override
      public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
        String kvCmd = "";

        String setSubSampFile = proj.SAMPLE_SUBSET_FILENAME.getValue();
        String subSampFile = variables == null ? null : variables.get(this).get(sampleSubsetReq);
        String setPFBFile = proj.CUSTOM_PFB_FILENAME.getValue();
        String pfbOutputFile = variables == null ? null : variables.get(this).get(outputFileReq);

        if (subSampFile != null && !ext.verifyDirFormat(setSubSampFile).equals(subSampFile)) {
          kvCmd += " SAMPLE_SUBSET_FILENAME=" + subSampFile;
        }
        if (pfbOutputFile != null && !ext.verifyDirFormat(setPFBFile).equals(pfbOutputFile)) {
          kvCmd += " CUSTOM_PFB_FILENAME=" + pfbOutputFile;
        }

        String projPropFile = proj.getPropertyFilename();
        StringBuilder cmd = new StringBuilder();
        if (kvCmd.length() > 0) {
          cmd.append(Files.getRunString())
             .append(GenvisisWorkflow.PROJ_PROP_UPDATE_STR + projPropFile).append(kvCmd)
             .append("\n");
        }
        return cmd.append(Files.getRunString()).append(" ").append(PFB.class.getName()).append(" ")
                  .append(CLI.ARG_PROJ).append("=").append(proj.getPropertyFilename()).append(" ")
                  .append(CLI.ARG_LOG).append(proj.getLog().getFilename()).toString();
      }

      @Override
      public boolean checkIfOutputExists(Project proj,
                                         Map<Step, Map<Requirement, String>> variables) {
        String subSampFile = variables.get(this).get(sampleSubsetReq);
        String pfbOutputFile = variables.get(this).get(outputFileReq);
        return Files.exists(pfbOutputFile) || Files.exists(ext.rootOf(subSampFile) + ".pfb");
      }
    });
  }

  Step generateSexCentroidsStep() {
    final RequirementSet reqSet = RequirementSetBuilder.and()
                                                       .add(GenvisisWorkflow.getNumThreadsReq());
    return register(new Step("Create Sex-Specific Centroids; Filter PFB file", "", reqSet,
                             EnumSet.of(Requirement.Flag.RUNTIME, Requirement.Flag.MULTITHREADED),
                             priority()) {

      @Override
      public void setNecessaryPreRunProperties(Project proj,
                                               Map<Step, Map<Requirement, String>> variables) {

        int numThreads = StepBuilder.resolveThreads(proj,
                                                    variables.get(this)
                                                             .get(GenvisisWorkflow.getNumThreadsReq()));
        GenvisisWorkflow.maybeSetProjNumThreads(proj, numThreads);
      }

      @Override
      public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
        String malePFB;
        String femalePFB;
        String centFilePathM;
        String centFilePathF;
        String outputDir = proj.DATA_DIRECTORY.getValue();
        malePFB = outputDir + "males.pfb";
        femalePFB = outputDir + "females.pfb";
        centFilePathM = outputDir + "sexSpecific_Male.cent";
        centFilePathF = outputDir + "sexSpecific_Female.cent";

        int numThreads = StepBuilder.resolveThreads(proj,
                                                    variables.get(this)
                                                             .get(GenvisisWorkflow.getNumThreadsReq()));
        Centroids.computeSexSpecificCentroids(proj, new String[] {malePFB, femalePFB},
                                              new String[] {centFilePathM, centFilePathF},
                                              numThreads);

      }

      @Override
      public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
        int numThreads = StepBuilder.resolveThreads(proj,
                                                    variables == null ? "-1"
                                                                      : variables.get(this)
                                                                                 .get(GenvisisWorkflow.getNumThreadsReq()));
        String mainCmd = Files.getRunString() + " cnv.filesys.Centroids proj="
                         + proj.getPropertyFilename() + " -sexSpecific "
                         + PSF.Ext.NUM_THREADS_COMMAND + numThreads;
        return mainCmd;
      }

      @Override
      public boolean checkIfOutputExists(Project proj,
                                         Map<Step, Map<Requirement, String>> variables) {
        String malePFB;
        String femalePFB;
        String centFilePathM;
        String centFilePathF;
        String outputDir = proj.DATA_DIRECTORY.getValue();
        malePFB = outputDir + "males.pfb";
        femalePFB = outputDir + "females.pfb";
        centFilePathM = outputDir + "sexSpecific_Male.cent";
        centFilePathF = outputDir + "sexSpecific_Female.cent";
        boolean exists = Files.exists(malePFB);
        exists = exists && Files.exists(femalePFB);
        exists = exists && Files.exists(centFilePathM);
        exists = exists && Files.exists(centFilePathF);
        return exists;
      }

    });
  }

  Step generateCNVStep(Project proj, Step pfbStep, Step gcModelStep) {
    final Requirement hmmFile = new Requirement.FileRequirement("Hidden Markov Model File Must Exist",
                                                                proj.HMM_FILENAME.getValue());
    final Requirement pfbStepReq = new Requirement.StepRequirement(pfbStep);
    final Requirement pfbFileReq = new Requirement.FileRequirement("PFB File Must Exist",
                                                                   proj.CUSTOM_PFB_FILENAME.getValue());
    final Requirement gcModelStepReq = new Requirement.StepRequirement(gcModelStep);
    final Requirement gcModelFileReq = new Requirement.FileRequirement("GCMODEL File Must Exist",
                                                                       proj.GC_MODEL_FILENAME.getValue());
    final Requirement callingTypeReq = new Requirement.EnumRequirement(CNVCaller.CNV_SCOPE_DESC,
                                                                       CNVCaller.CALLING_SCOPE.AUTOSOMAL);
    final Requirement useCentroidsReq = new Requirement.OptionalBoolRequirement("If calling chromosomal CNVs, use sex-specific centroids to recalculate LRR/BAF values?",
                                                                                true);
    final Requirement outputFileReq = new Requirement.OutputFileRequirement("Output filename.",
                                                                            "cnvs/genvisis.cnv") {

      @Override
      public boolean checkRequirement(Project proj, String arg, Set<Step> stepSelections,
                                      Map<Step, Map<Requirement, String>> variables) {
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
                                                       .add(GenvisisWorkflow.getNumThreadsReq())
                                                       .add(outputFileReq);

    return register(new Step("Call CNVs", "", reqSet,
                             EnumSet.of(Requirement.Flag.MEMORY, Requirement.Flag.MULTITHREADED),
                             priority()) {

      @Override
      public void setNecessaryPreRunProperties(Project proj,
                                               Map<Step, Map<Requirement, String>> variables) {
        String hmmP = proj.HMM_FILENAME.getValue();
        String hmmG = variables.get(this).get(hmmFile);
        if (!hmmP.equals(hmmG)) {
          proj.HMM_FILENAME.setValue(hmmG);
        }
        String pfbP = proj.CUSTOM_PFB_FILENAME.getValue();
        String pfbG = variables.get(this).get(pfbFileReq);
        if (!pfbP.equals(pfbG)) {
          proj.CUSTOM_PFB_FILENAME.setValue(pfbG);
        }
        String gcmP = proj.GC_MODEL_FILENAME.getValue();
        String gcmG = variables.get(this).get(gcModelFileReq);
        if (!gcmP.equals(gcmG)) {
          proj.GC_MODEL_FILENAME.setValue(gcmG);
        }
        int numThreads = StepBuilder.resolveThreads(proj,
                                                    variables.get(this)
                                                             .get(GenvisisWorkflow.getNumThreadsReq()));
        GenvisisWorkflow.maybeSetProjNumThreads(proj, numThreads);
      }

      @Override
      public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
        int numThreads = StepBuilder.resolveThreads(proj,
                                                    variables.get(this)
                                                             .get(GenvisisWorkflow.getNumThreadsReq()));
        GenvisisWorkflow.maybeSetProjNumThreads(proj, numThreads);
        String output = variables.get(this).get(outputFileReq); // gets PROJ_DIR prepended, so NOT
                                                               // ABSOLUTE

        CALLING_SCOPE scope = CALLING_SCOPE.valueOf(variables.get(this).get(callingTypeReq));

        (new File(ext.parseDirectoryOfFile(proj.PROJECT_DIRECTORY.getValue() + output))).mkdirs();

        String[] samples = proj.getSamples();
        boolean useCentroids = Boolean.valueOf(variables.get(this).get(useCentroidsReq));
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
      public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
        String kvCmd = "";

        String hmmP = proj.HMM_FILENAME.getValue();
        String hmmG = variables.get(this).get(hmmFile);
        if (hmmG != null && !hmmP.equals(hmmG)) {
          kvCmd += " HMM_FILENAME=" + hmmG;
        }
        String pfbP = proj.CUSTOM_PFB_FILENAME.getValue();
        String pfbG = variables.get(this).get(pfbFileReq);
        if (pfbG != null && !pfbP.equals(pfbG)) {
          kvCmd += " CUSTOM_PFB_FILENAME=" + pfbG;
        }
        String gcmP = proj.GC_MODEL_FILENAME.getValue();
        String gcmG = variables.get(this).get(gcModelFileReq);
        if (gcmG != null && !gcmP.equals(gcmG)) {
          kvCmd += " GC_MODEL_FILENAME=" + gcmG;
        }

        int numThreads = StepBuilder.resolveThreads(proj,
                                                    variables.get(this)
                                                             .get(GenvisisWorkflow.getNumThreadsReq()));
        if (numThreads != proj.NUM_THREADS.getValue()) {
          kvCmd += " " + proj.NUM_THREADS.getName() + "=" + numThreads;
        }
        String projPropFile = proj.getPropertyFilename();
        StringBuilder cmd = new StringBuilder();
        if (kvCmd.length() > 0) {
          cmd.append(Files.getRunString())
             .append(GenvisisWorkflow.PROJ_PROP_UPDATE_STR + projPropFile).append(kvCmd)
             .append("\n");
        }

        boolean useCentroids = Boolean.valueOf(variables.get(this).get(useCentroidsReq));
        CALLING_SCOPE scope = CALLING_SCOPE.valueOf(variables.get(this).get(callingTypeReq));

        String autoCmd = cmd.append(Files.getRunString())
                            .append(" cnv.hmm.CNVCaller proj=" + projPropFile)
                            .append(" out=" + variables.get(this).get(outputFileReq)).append(" ")
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
      public boolean checkIfOutputExists(Project proj,
                                         Map<Step, Map<Requirement, String>> variables) {
        String output = variables.get(this).get(outputFileReq);
        return Files.exists(proj.PROJECT_DIRECTORY.getValue() + output);
      }

    });
  }

  Step generatePCCorrectedProjectStep(Project proj, final Step parseSamplesStep) {
    final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    final Requirement numPCsReq = new Requirement.PosIntRequirement("Number of principal components for correction.",
                                                                    MitoPipeline.DEFAULT_NUM_COMPONENTS);
    final Requirement outputBaseReq = new Requirement.OutputFileRequirement("Output file path (relative to project directory) and baseName for principal components correction files",
                                                                            MitoPipeline.FILE_BASE) {

      @Override
      public boolean checkRequirement(Project proj, String arg, Set<Step> stepSelections,
                                      Map<Step, Map<Requirement, String>> variables) {
        String outputBase = proj.PROJECT_DIRECTORY.getValue() + arg;
        String finalReport = outputBase + PCA.FILE_EXTs[0];
        return super.checkRequirement(proj, finalReport, stepSelections, variables);
      }
    };
    final Requirement callrateReq = new Requirement.DoubleRequirement("Call-rate filter for determining high-quality markers",
                                                                      MitoPipeline.DEFAULT_MKR_CALLRATE_FILTER,
                                                                      0.0, 1.0);
    final Requirement recomputeLrrReq = new Requirement.OptionalBoolRequirement("Re-compute Log-R Ratio values? (usually false if LRRs already exist)",
                                                                                false);
    final Requirement tempDirReq = new Requirement.OptionalFileRequirement("Temporary directory for intermediate files (which tend to be very large)",
                                                                           "");
    final Requirement correctionStrategyReq = new Requirement.EnumRequirement("Correction Type",
                                                                              CORRECTION_TYPE.XY);
    final Requirement sexChromosomeStrategyReq = new Requirement.EnumRequirement("Sex Chromosome Strategy",
                                                                                 CHROMOSOME_X_STRATEGY.BIOLOGICAL);
    final Requirement setupCNVCalling = new Requirement.OptionalBoolRequirement("Create script with steps to process corrected data and call CNVs?",
                                                                                false);

    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
                                                       .add(numPCsReq).add(outputBaseReq)
                                                       .add(callrateReq).add(recomputeLrrReq)
                                                       .add(tempDirReq).add(correctionStrategyReq)
                                                       .add(sexChromosomeStrategyReq)
                                                       .add(GenvisisWorkflow.getNumThreadsReq())
                                                       .add(setupCNVCalling);

    return register(new Step("Create PC-Corrected Project", "", reqSet,
                             EnumSet.of(Requirement.Flag.MEMORY, Requirement.Flag.RUNTIME,
                                        Requirement.Flag.MULTITHREADED),
                             priority()) {

      @Override
      public void setNecessaryPreRunProperties(Project proj,
                                               Map<Step, Map<Requirement, String>> variables) {
        // not needed for step
      }

      @Override
      public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
        int numComponents = Integer.parseInt(variables.get(this).get(numPCsReq));
        String outputBase = variables.get(this).get(outputBaseReq);
        double markerCallRateFilter = Double.parseDouble(variables.get(this).get(callrateReq));
        boolean recomputeLRRPCs = Boolean.parseBoolean(variables.get(this).get(recomputeLrrReq));
        String tmpDir = variables.get(this).get(tempDirReq);
        if ("".equals(tmpDir.trim())) {
          tmpDir = null;
        }
        CORRECTION_TYPE type = CORRECTION_TYPE.valueOf(variables.get(this)
                                                                .get(correctionStrategyReq));
        CHROMOSOME_X_STRATEGY strategy = CHROMOSOME_X_STRATEGY.valueOf(variables.get(this)
                                                                                .get(sexChromosomeStrategyReq));

        int totalThreads = StepBuilder.resolveThreads(proj,
                                                      variables.get(this)
                                                               .get(GenvisisWorkflow.getNumThreadsReq()));
        boolean cnvCalling = Boolean.parseBoolean(variables.get(this).get(setupCNVCalling));
        String retMsg = PRoCtOR.shadow(proj, tmpDir, outputBase, markerCallRateFilter,
                                       recomputeLRRPCs, type, strategy, numComponents, totalThreads,
                                       cnvCalling);
        if (retMsg != null && !"".equals(retMsg)) {
          throw new RuntimeException(retMsg);
        }
      }

      @Override
      public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
        int numComponents = Integer.parseInt(variables.get(this).get(numPCsReq));
        String outputBase = variables.get(this).get(outputBaseReq);
        double markerCallRateFilter = Double.parseDouble(variables.get(this).get(callrateReq));
        boolean recomputeLRRPCs = Boolean.parseBoolean(variables.get(this).get(recomputeLrrReq));
        String tmpDir = variables.get(this).get(tempDirReq);
        if ("".equals(tmpDir.trim())) {
          tmpDir = null;
        }
        String correctionType = variables.get(this).get(correctionStrategyReq);
        String strategy = variables.get(this).get(sexChromosomeStrategyReq);

        int totalThreads = StepBuilder.resolveThreads(proj,
                                                      variables.get(this)
                                                               .get(GenvisisWorkflow.getNumThreadsReq()));

        boolean cnvCalling = Boolean.parseBoolean(variables.get(this).get(setupCNVCalling));
        String projPropFile = proj.getPropertyFilename();
        StringBuilder cmd = new StringBuilder();
        cmd.append(Files.getRunString()).append(" org.genvisis.cnv.manage.PRoCtOR").append(" proj=")
           .append(projPropFile).append(" numComponents=").append(numComponents)
           .append(" outputBase=").append(outputBase).append(" callrate=")
           .append(markerCallRateFilter).append(" recomputeLRR=").append(recomputeLRRPCs)
           .append(" type=").append(correctionType).append(" sexStrategy=").append(strategy)
           .append(" numThreads=").append(totalThreads);
        if (tmpDir != null) {
          cmd.append(" tmp=").append(tmpDir);
        }
        if (cnvCalling) {
          cmd.append(" -callCNVs");
        }

        return cmd.toString();
      }

      @Override
      public boolean checkIfOutputExists(Project proj,
                                         Map<Step, Map<Requirement, String>> variables) {
        String outputBase = proj.PROJECT_DIRECTORY.getValue()
                            + variables.get(this).get(outputBaseReq);
        String finalReport = outputBase + PCA.FILE_EXTs[0];
        return Files.exists(finalReport);
      }

    });

  }

  Step generateABLookupStep(final Step parseSamplesStep) {
    final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq);
    return register(new Step("Generate AB Lookup File", "", reqSet,
                             EnumSet.of(Requirement.Flag.RUNTIME), priority()) {

      @Override
      public void setNecessaryPreRunProperties(Project proj,
                                               Map<Step, Map<Requirement, String>> variables) {
        // Nothing to do here
      }

      @Override
      public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
        String filename = proj.PROJECT_DIRECTORY.getValue()
                          + ext.addToRoot(ABLookup.DEFAULT_AB_FILE, "_parsed");
        ABLookup.parseABLookup(proj, ABSource.VCF, filename);

        if (ABLookup.fillInMissingAlleles(proj, filename, proj.getLocationOfSNP_Map(true), false)) {
          ABLookup.applyABLookupToFullSampleFiles(proj, filename);
        } else {
          throw new RuntimeException("Failed to fill in missing alleles - please check log for more info.");
        }
      }

      @Override
      public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
        String filename = proj.PROJECT_DIRECTORY.getValue()
                          + ext.addToRoot(ABLookup.DEFAULT_AB_FILE, "_parsed");
        String projFile = proj.getPropertyFilename();
        String mapFile = proj.getLocationOfSNP_Map(true);

        List<String> baseCommand = ImmutableList.of(Files.getRunString(), ABLookup.class.getName(),
                                                    CLI.formCmdLineArg(CLI.ARG_PROJ, projFile));
        List<String> commandVcf = Lists.newArrayList(baseCommand);
        commandVcf.add(CLI.formCmdLineArg(CLI.ARG_OUTFILE, filename));
        commandVcf.add(CLI.formCmdLineFlag(ABLookup.FLAGS_VCF));

        List<String> commandPartial = Lists.newArrayList(baseCommand);
        commandPartial.add(CLI.formCmdLineArg(ABLookup.ARGS_PARTAB, filename));
        commandPartial.add(CLI.formCmdLineArg(ABLookup.ARGS_MAP, mapFile));

        List<String> commandProp = Lists.newArrayList(ImmutableList.of(Files.getRunString(),
                                                                       Project.class.getName(),
                                                                       CLI.formCmdLineArg(CLI.ARG_PROJ,
                                                                                          projFile)));
        commandProp.add(CLI.formCmdLineArg(proj.AB_LOOKUP_FILENAME.getName(), filename));

        List<String> commandApply = Lists.newArrayList(baseCommand);
        commandApply.add(CLI.formCmdLineFlag(ABLookup.FLAGS_APPLYAB));

        StringBuilder cmd = new StringBuilder();

        cmd.append(Joiner.on(" ").join(commandVcf)).append("\n");
        cmd.append(Joiner.on(" ").join(commandPartial)).append("\n");
        cmd.append(Joiner.on(" ").join(commandProp)).append("\n");
        cmd.append(Joiner.on(" ").join(commandApply));

        return cmd.toString();
      }

      @Override
      public boolean checkIfOutputExists(Project proj,
                                         Map<Step, Map<Requirement, String>> variables) {
        return Files.exists(proj.AB_LOOKUP_FILENAME.getValue(false, false));
      }
    });
  }

  static int resolveThreads(Project proj, String arg) {
    int numThreads = Requirement.checkIntArgOrNeg1(arg);
    if (numThreads <= 0) {
      numThreads = proj.NUM_THREADS.getValue();
    }
    return numThreads;
  }

}
