package org.genvisis.cnv.analysis.pca;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.ExtProjectDataParser;
import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.cnv.qc.GcAdjustorParameter.GcAdjustorParameters;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.ArrayUtils.BooleanClassifier;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.SerializedFiles;
import org.pankratzlab.common.WorkerTrain;
import org.pankratzlab.common.WorkerTrain.AbstractProducer;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.stats.ICC;
import org.pankratzlab.common.stats.LeastSquares.LS_TYPE;
import org.pankratzlab.common.stats.Rscript.COLUMNS_MULTIPLOT;
import org.pankratzlab.common.stats.Rscript.ErrorBars;
import org.pankratzlab.common.stats.Rscript.GEOM_POINT_SIZE;
import org.pankratzlab.common.stats.Rscript.GeomText;
import org.pankratzlab.common.stats.Rscript.PLOT_DEVICE;
import org.pankratzlab.common.stats.Rscript.RScatter;
import org.pankratzlab.common.stats.Rscript.RScatters;
import org.pankratzlab.common.stats.Rscript.SCATTER_TYPE;
import org.pankratzlab.common.stats.StatsCrossTabs.STAT_TYPE;
import org.pankratzlab.common.stats.StatsCrossTabs.StatsCrossTabRank;
import org.pankratzlab.common.stats.StatsCrossTabs.VALUE_TYPE;

import com.google.common.primitives.Ints;

public class CorrectionIterator implements Serializable {

  /**
   *
   */
  private static final long serialVersionUID = 1L;
  private static final String BOX_Y = "mtDNA CN estimate";
  public static final String FINAL_EST_DIR = "typed/";
  private final Project proj;
  private final String markesToEvaluate;
  private final String samplesToBuildModels;
  private final ITERATION_TYPE iType;
  private final ORDER_TYPE oType;
  private final MODEL_BUILDER_TYPE bType;
  private final String outputDir;
  private final LS_TYPE lType;
  private final int numthreads;
  private IterationResult iterationResult;
  private final boolean recomputeLRR;
  private final double pcPercent;

  public CorrectionIterator(Project proj, String markesToEvaluate, String samplesToBuildModels,
                            ITERATION_TYPE iType, ORDER_TYPE oType, MODEL_BUILDER_TYPE bType,
                            String outputDir, LS_TYPE lType, boolean recomputeLRR, double pcPercent,
                            int numthreads) {
    super();
    this.proj = proj;
    this.markesToEvaluate = markesToEvaluate;
    this.samplesToBuildModels = samplesToBuildModels;
    this.iType = iType;
    this.oType = oType;
    this.bType = bType;
    this.outputDir = outputDir;
    this.lType = lType;
    this.numthreads = numthreads;
    this.pcPercent = pcPercent;
    this.recomputeLRR = recomputeLRR;
    // this.lrrSdCut = lrrSdCut;
    // this.callRateCut = callRateCut;

  }

  public void run() {
    iterationResult = run(proj, markesToEvaluate, samplesToBuildModels, iType, oType, bType,
                          outputDir, lType, recomputeLRR, pcPercent, numthreads);
  }

  public IterationResult getIterationResult() {
    return iterationResult;
  }

  public enum ITERATION_TYPE {
    /**
     * The evaluation happens with the addition of other independent variables in addition to PCs
     */
    WITHOUT_INDEPS,
    /**
     * other independent variables are not added
     */
    // WITH_INDEPS;

  }

  public enum MODEL_BUILDER_TYPE {
    /**
     * We load the additional sample (union with not excluded) file for model building
     */
    WITH_QC_BUILDERS,
    /**
     * We build models with everyone, except excluded individuals
     */
    // WITHOUT_QC_BUILDERS;

  }

  public enum ORDER_TYPE {
    /**
     * PCS are ranked by the amount of variance explained in the estimate of interest
     */
    // RANK_R2,
    /**
     * PCs are regressed by their natural ordering (PC1,2,3)
     */
    NATURAL,
    /**
     * PCS are ranked by the amount of variance explained within a stepwise regression
     */
    STEPWISE_RANK_R2,

    /**
     * PCS are ranked by spearman abs r
     */
    // RANK_R,
    /**
     * PCs are filtered for an association with known QC metrics from {@link LrrSd}
     */
    // QC_ASSOCIATION;

  }

  private IterationResult run(Project proj, String markesToEvaluate, String samplesToBuildModels,
                              ITERATION_TYPE iType, ORDER_TYPE oType, MODEL_BUILDER_TYPE bType,
                              String outputDir, LS_TYPE lType, boolean recomputeLRR,
                              double pcPercent, int numthreads) {
    Logger log = proj.getLog();
    CorrectionEvaluator cEvaluator = null;
    String output = outputDir + "correctionEval_" + iType + "_" + oType + "_" + bType;
    IterationResult iterationResult = new IterationResult(output, iType, oType, bType, pcPercent);
    EvaluationResult[] precomputed = null;
    if (Files.exists(iterationResult.getOutputSer())) {
      try {
        precomputed = EvaluationResult.readSerial(iterationResult.getOutputSer(), log);
        log.reportTimeInfo("Loading precomputed results from " + iterationResult.getOutputSer());

      } catch (Exception e) {// to be safe
        precomputed = null;
      }
    }

    proj.getLog().reportTimeInfo("Loading " + proj.INTENSITY_PC_FILENAME.getValue());

    new File(outputDir).mkdirs();
    log.reportTimeInfo("Beginning iteration evaluation:");
    log.reportTimeInfo("PC file: " + proj.INTENSITY_PC_FILENAME.getValue());
    log.reportTimeInfo("Iteration type : " + iType);
    log.reportTimeInfo("Order type : " + oType);
    log.reportTimeInfo("Model building type: " + bType);

    boolean[] samplesForModels = null;
    boolean valid = true;
    switch (bType) {
      // case WITHOUT_QC_BUILDERS:
      //
      // samplesForModels = Array.booleanArray(proj.getNumberOfParsedSamples(), true);
      // break;
      case WITH_QC_BUILDERS:
        if (!Files.exists(samplesToBuildModels)) {
          log.reportError("Model building type was set to " + bType + " but the sample file "
                          + samplesToBuildModels + " did not exist");
          valid = false;
        } else {
          log.reportTimeInfo("Loading model builders from " + samplesToBuildModels);
          String[] sampsForMods = HashVec.loadFileToStringArray(samplesToBuildModels, false,
                                                                new int[] {0}, false, true, "\t");
          log.reportTimeInfo("Loaded " + sampsForMods.length + " model builders from "
                             + samplesToBuildModels);

          int[] indices = ext.indexLargeFactors(sampsForMods, proj.getSamples(), true,
                                                proj.getLog(), true);
          samplesForModels = ArrayUtils.booleanArray(proj.getNumberOfParsedSamples(), false);
          for (int i = 0; i < indices.length; i++) {
            samplesForModels[indices[i]] = true;
          }
        }
        log.reportTimeInfo("Loaded " + ArrayUtils.booleanArraySum(samplesForModels)
                           + " model builders from " + samplesToBuildModels);

        log.reportTimeInfo(ArrayUtils.booleanArraySum(samplesForModels)
                           + " model builders from after QC filtering");

        break;
      default:
        break;

    }

    switch (iType) {
      case WITHOUT_INDEPS:
        log.reportTimeInfo("Evaluating with " + ArrayUtils.booleanArraySum(samplesForModels)
                           + " samples, no additional independent variables");
        break;

      default:
        return null;

    }
    // } else {
    // log.reportTimeWarning("Loading precomputed sample preparation from " +
    // iterationResult.getBasePrep());
    // BasicPrep basicPrep = BasicPrep.readSerial(iterationResult.getBasePrep(), log);
    // iterationResult.setBasicPrep(basicPrep);
    // samplesForModels = basicPrep.getSamplesForModels();
    //
    // }
    PrincipalComponentsResiduals pcResiduals = proj.loadPcResids();
    int totalNumPCs = pcResiduals.getTotalNumComponents();
    int numSamples = ArrayUtils.booleanArraySum(samplesForModels);
    log.reportTimeInfo("Detected " + pcResiduals.getTotalNumComponents() + "available PCs in "
                       + proj.INTENSITY_PC_FILENAME.getValue());
    log.reportTimeInfo("Using a total of " + numSamples + " samples");
    int numPCsToEvaluate = Math.round((float) pcPercent * numSamples);
    log.reportTimeInfo("PC percent set to " + pcPercent + " giving " + numPCsToEvaluate
                       + " total pcs to use out of " + totalNumPCs);
    numPCsToEvaluate = Math.min(numPCsToEvaluate, totalNumPCs);
    log.reportTimeInfo("Setting number of pcs to " + numPCsToEvaluate);
    proj.INTENSITY_PC_NUM_COMPONENTS.setValue(numPCsToEvaluate);
    pcResiduals = proj.loadPcResids();
    pcResiduals.fillInMissing();
    pcResiduals.setMarkersToAssessFile(markesToEvaluate);
    Files.writeArray(ArrayUtils.subArray(proj.getSamples(), samplesForModels),
                     outputDir + iType + "_" + oType + "_" + bType + "_samplesForModels.txt");
    pcResiduals.setHomozygousOnly(true);
    proj.getLog().reportTimeWarning("In gc-correction mode now, using "
                                    + proj.GC_CORRECTION_PARAMETERS_FILENAMES.getValue()[0]);
    GcAdjustorParameters params = GcAdjustorParameters.readSerial(proj.GC_CORRECTION_PARAMETERS_FILENAMES.getValue()[0],
                                                                  proj.getLog());
    if (params.getCentroids() == null && recomputeLRR) {
      throw new IllegalArgumentException("Must have centroids");
    } else {
      if (!recomputeLRR) {
        proj.getLog().reportTimeInfo("Using gc correction (no lrr recomp) from "
                                     + proj.GC_CORRECTION_PARAMETERS_FILENAMES.getValue()[0]);

      } else {
        proj.getLog().reportTimeInfo("Using centroids and gc correction from "
                                     + proj.GC_CORRECTION_PARAMETERS_FILENAMES.getValue()[0]);
      }
    }
    // GcAdjustorParameters params =
    pcResiduals.setParams(params);
    pcResiduals.computeAssessmentDataMedians();
    cEvaluator = new CorrectionEvaluator(proj, pcResiduals, null, null, null, null, lType);
    int[] order = null;
    double[][] extraIndeps = null;
    StatsCrossTabRank sTabRank = null;

    iterationResult.setValid(valid);
    if (valid) {
      // sTabRank = pcResiduals.getStatRankFor(pcResiduals.getMedians(), extraIndeps,
      // samplesForModels, "RAW_MEDIANS", STAT_TYPE.LIN_REGRESSION, VALUE_TYPE.STAT, false,
      // numthreads, proj.getLog());
      // sTabRank.dump(iterationResult.getOutputRank(), oType != ORDER_TYPE.NATURAL , log);

      if (precomputed == null) {
        sTabRank = pcResiduals.getStatRankFor(pcResiduals.getMedians(), extraIndeps,
                                              samplesForModels, "RAW_MEDIANS",
                                              STAT_TYPE.LIN_REGRESSION, VALUE_TYPE.STAT,
                                              oType == ORDER_TYPE.STEPWISE_RANK_R2, numthreads,
                                              proj.getLog());

        sTabRank.dump(iterationResult.getOutputRank(),
                      oType != ORDER_TYPE.NATURAL && oType != ORDER_TYPE.STEPWISE_RANK_R2, log);

        switch (oType) {
          case NATURAL:
            order = null;
            break;
          // case RANK_R2:
          // order = new int[sTabRank.getOrder().length];
          // for (int i = 0; i < sTabRank.getOrder().length; i++) {
          // order[i] = sTabRank.getOrder()[i] + 1;// one based for pcs
          // }
          // break;
          case STEPWISE_RANK_R2:
            order = new int[sTabRank.getOrder().length];
            for (int i = 0; i < sTabRank.getOrder().length; i++) {
              order[i] = sTabRank.getOrder()[i] + 1;// one based for pcs
            }
            // log.reportTimeInfo("PC steArray.toStr(order));
            break;

          default:
            break;
        }
      } else {
        log.reportTimeWarning("Skipping PC selection, relying on pre-computed results");
      }
      // boolean[] samplesToEvaluate = proj.getSamplesToInclude(null);
      boolean[] samplesToEvaluate = samplesForModels;

      if (order != null && oType == ORDER_TYPE.STEPWISE_RANK_R2) {
        String out = outputDir + iType + "_" + oType + "_" + bType + "_PC_" + oType
                     + "_Selection.txt";
        ArrayList<String> outOrder = new ArrayList<>();
        outOrder.add(oType + "_PC_RANK\tOriginal_PC_Rank");
        for (int i = 0; i < order.length; i++) {
          outOrder.add((i + 1) + "\tPC" + order[i]);
        }
        Files.writeIterable(outOrder, out);
      }
      Files.writeArray(ArrayUtils.subArray(proj.getSamples(), samplesToEvaluate),
                       outputDir + iType + "_" + oType + "_" + bType + "_samplesForEval.txt");

      cEvaluator = new CorrectionEvaluator(proj, pcResiduals, precomputed, order,
                                           new boolean[][] {samplesForModels, samplesToEvaluate},
                                           extraIndeps, lType);
      BasicPrep basicPrep = new BasicPrep(cEvaluator.getParser().getNumericData(),
                                          cEvaluator.getParser().getNumericDataTitles(),
                                          samplesToEvaluate, samplesForModels, null);
      BasicPrep.serialize(basicPrep, iterationResult.getBasePrep());
      iterationResult.setBasicPrep(basicPrep);
      log.reportTimeInfo(ArrayUtils.booleanArraySum(samplesForModels) + " samples for models");
      log.reportTimeInfo(ArrayUtils.booleanArraySum(samplesToEvaluate) + " samples for evaluation");

      ArrayList<EvaluationResult> store = new ArrayList<>();

      try (PrintWriter writer = Files.openAppropriateWriter(iterationResult.getOutputSummary());
           WorkerTrain<EvaluationResult> train = new WorkerTrain<>(cEvaluator, numthreads,
                                                                   numthreads, proj.getLog())) {
        int index = 0;
        while (train.hasNext()) {
          EvaluationResult result = train.next();
          result.setItType(iType);
          result.setOrType(oType);
          result.setbType(bType);
          if (index == 0) {
            writer.println(ArrayUtils.toStr(result.getHeader()));
          }
          writer.println(ArrayUtils.toStr(result.getData()));
          index++;
          result.shrink();
          store.add(result);
        }
      } catch (Exception e) {
        proj.getLog().reportError("Error writing to " + iterationResult.getOutputSummary());
        proj.getLog().reportException(e);
      }
      EvaluationResult.serialize(store.toArray(new EvaluationResult[store.size()]),
                                 iterationResult.getOutputSer());
    }

    return iterationResult;
  }

  private static class BasicPrep implements Serializable {

    /**
     *
     */
    private static final long serialVersionUID = 1L;
    private final double[][] numericData;
    private final boolean[] samplesForModels;
    private final boolean[] samplesToEvaluate;
    private final String[] numericTitles;

    // private StatsCrossTabRank sTabRank;

    public BasicPrep(double[][] numericData, String[] numericTitles, boolean[] samplesToEvaluate,
                     boolean[] samplesForModels, StatsCrossTabRank sTabRank) {
      super();
      this.numericData = numericData;
      this.samplesToEvaluate = samplesToEvaluate;
      this.samplesForModels = samplesForModels;
      this.numericTitles = numericTitles;
      // this.sTabRank = sTabRank;
    }

    public boolean[] getSamplesForModels() {
      return samplesForModels;
    }

    public boolean[] getSamplesToEvaluate() {
      return samplesToEvaluate;
    }

    public double[][] getNumericData() {
      return numericData;
    }

    public String[] getNumericTitles() {
      return numericTitles;
    }

    public static void serialize(BasicPrep basicPrep, String fileName) {
      SerializedFiles.writeSerial(basicPrep, fileName, true);
    }

  }

  private static class IterationResult implements Serializable {

    /**
     *
     */
    private static final long serialVersionUID = 1L;
    private static final String[] HERIT_ADDITIONS = new String[] {"NUM_PC",
                                                                  "MERLIN_PROPORTION_HERITABLITY",
                                                                  "SOLAR_PROPORTION_HERITABLITY",
                                                                  "SOLAR_PVAL", "SOLAR_ST_ERROR",
                                                                  "SOLAR_KURT",
                                                                  "SOLAR_KURT_WARNING"};
    private String outputRoot;
    private String outputSer;
    private final String outputRank;
    private final String outputSummary;
    private final String rankRscript;
    private final String evalRscript;
    private String heritRscript;
    private String heritSummary;
    private BasicPrep basicPrep;
    private final ITERATION_TYPE iType;
    private final MODEL_BUILDER_TYPE bType;
    private final ORDER_TYPE oType;
    private final String rankplot;
    private final String evalPlot;
    private String heritPlot;
    private final String basePrep;
    private final String boxPlot;

    public IterationResult(String outputRoot, ITERATION_TYPE iType, ORDER_TYPE oType,
                           MODEL_BUILDER_TYPE bType, double pcPercent) {
      super();
      this.outputRoot = outputRoot;
      outputSer = outputRoot + ".summary.ser";
      outputSummary = outputRoot + ".summary.txt";
      outputRank = outputRoot + ".rank.txt";
      rankRscript = outputRoot + ".rank.rscript";
      rankplot = outputRoot + ".rank.Plot.pdf";
      evalRscript = outputRoot + ".eval.rscript";
      evalPlot = outputRoot + ".eval.Plot.pdf";
      basePrep = outputRoot + "basePrep.ser";
      boxPlot = outputRoot + ".boxPlot.txt";
      this.iType = iType;
      this.oType = oType;
      this.bType = bType;

    }

    // public double getPcPercent() {
    // return pcPercent;
    // }
    //
    // public boolean isValid() {
    // return valid;
    // }

    public String getBasePrep() {
      return basePrep;
    }

    public void setBasicPrep(BasicPrep basicPrep) {
      this.basicPrep = basicPrep;
    }

    public BasicPrep getBasicPrep() {
      return basicPrep;
    }

    public RScatter plotRank(Logger log) {
      RScatter rScatter = new RScatter(outputRank, rankRscript, ext.rootOf(outputRank), rankplot,
                                       "OriginalOrder", new String[] {"Stat"}, SCATTER_TYPE.POINT,
                                       log);
      rScatter.setyRange(new double[] {0, 1});
      rScatter.setxLabel("PC (" + oType + " - sorted)");
      rScatter.setyLabel("Rsq");
      rScatter.setTitle(iType + " " + bType);
      rScatter.setgPoint_SIZE(GEOM_POINT_SIZE.GEOM_POINT);
      rScatter.execute();
      return rScatter;
    }

    public BooleanClassifier[] estimateBoxPlots(Project proj, String[] sampleDataStratCats,
                                                String[] numericStratCats, List<RScatter> scatters,
                                                Logger log) {
      String boxDir = ext.parseDirectoryOfFile(boxPlot) + "boxPlots/";
      String quantDir = ext.parseDirectoryOfFile(boxPlot) + "quants/";
      new File(boxDir).mkdirs();
      new File(quantDir).mkdirs();

      ExtProjectDataParser.ProjectDataParserBuilder builder = new ExtProjectDataParser.ProjectDataParserBuilder();
      builder.sampleBased(true);
      builder.treatAllNumeric(false);
      builder.requireAll(true);
      builder.verbose(true);
      builder.dataKeyColumnName("DNA");
      builder.stringDataTitles(sampleDataStratCats);
      if (numericStratCats != null) {
        int[] indices = ext.indexFactors(numericStratCats,
                                         Files.getHeaderOfFile(proj.SAMPLE_DATA_FILENAME.getValue(),
                                                               proj.getLog()),
                                         true);
        ArrayList<Integer> use = new ArrayList<>();
        for (int i = 0; i < indices.length; i++) {
          if (indices[i] >= 0) {
            use.add(i);
          }
        }
        builder.numericDataTitles(ArrayUtils.subArray(numericStratCats, Ints.toArray(use)));
      }
      BooleanClassifier[] classifiers = new BooleanClassifier[sampleDataStratCats.length];
      ArrayList<RScatter> rScatters = new ArrayList<>();
      try {
        log.reportTimeInfo("Loading " + proj.SAMPLE_DATA_FILENAME.getValue());
        ExtProjectDataParser parser = builder.build(proj, proj.SAMPLE_DATA_FILENAME.getValue());
        parser.determineIndicesFromTitles();
        parser.loadData();

        log.reportTimeInfo("Finished loading " + proj.SAMPLE_DATA_FILENAME.getValue());
        for (int i = 0; i < sampleDataStratCats.length; i++) {
          BooleanClassifier booleanClassifier = ArrayUtils.classifyStringsToBoolean(parser.getStringDataForTitle(sampleDataStratCats[i]),
                                                                                    new String[] {"NaN"});
          classifiers[i] = booleanClassifier;
          EvaluationResult[] evaluationResults = EvaluationResult.readSerial(outputSer, log);
          try {
            String outputBox = ext.addToRoot(boxPlot, sampleDataStratCats[i]);
            PrintWriter writer = Files.openAppropriateWriter(outputBox);
            writer.print(sampleDataStratCats[i]);
            ArrayList<String> pcYs = new ArrayList<>();
            for (int j = 0; j < evaluationResults.length; j++) {
              String title = "PC" + j;
              writer.print("\t" + title);
              pcYs.add(title);
            }
            writer.println();
            for (int sampleIndex = 0; sampleIndex < evaluationResults[0].getEstimateData().length; sampleIndex++) {
              for (int groupIndex = 0; groupIndex < booleanClassifier.getClassified().length; groupIndex++) {
                int n = 0;
                for (int j = 0; j < getBasicPrep().getSamplesToEvaluate().length; j++) {
                  if (getBasicPrep().getSamplesToEvaluate()[j]
                      && booleanClassifier.getClassified()[groupIndex][j]) {
                    n++;
                  }
                }
                if (getBasicPrep().getSamplesToEvaluate()[sampleIndex]
                    && booleanClassifier.getClassified()[groupIndex][sampleIndex]) {
                  writer.print(booleanClassifier.getTitles()[groupIndex] + "_n_" + n);
                  for (EvaluationResult evaluationResult : evaluationResults) {
                    writer.print("\t" + evaluationResult.getEstimateData()[sampleIndex]);
                  }
                  writer.println();
                }

              }
            }
            writer.close();

            plotFirst(sampleDataStratCats, log, boxDir, rScatters, i, evaluationResults, outputBox,
                      pcYs);
            plotSkips(sampleDataStratCats, log, boxDir, rScatters, i, evaluationResults, outputBox,
                      pcYs);
            RScatter rScatterFirstLast = plotFirstLast(sampleDataStratCats, log, boxDir, i,
                                                       outputBox, pcYs);
            rScatters.add(rScatterFirstLast);

            if (numericStratCats != null) {

              for (String numericStratCat : numericStratCats) {

                String outputBoxSub = ext.addToRoot(outputBox,
                                                    sampleDataStratCats[i] + "_" + numericStratCat);
                double[] data = parser.getNumericDataForTitle(numericStratCat);
                if (data != null) {
                  PrintWriter writerSub = Files.openAppropriateWriter(outputBoxSub);
                  writerSub.print(sampleDataStratCats[i] + "\t" + numericStratCat);
                  ArrayList<String> pcYsub = new ArrayList<>();

                  for (int PC = 0; PC < evaluationResults.length; PC++) {
                    String title = "PC" + PC;
                    writerSub.print("\t" + title);
                    pcYsub.add(title);
                  }
                  writerSub.println();
                  for (int sampleIndex = 0; sampleIndex < evaluationResults[0].getEstimateData().length; sampleIndex++) {
                    for (int groupIndex = 0; groupIndex < booleanClassifier.getClassified().length; groupIndex++) {
                      int n = 0;
                      for (int k = 0; k < getBasicPrep().getSamplesToEvaluate().length; k++) {
                        if (getBasicPrep().getSamplesToEvaluate()[k]
                            && booleanClassifier.getClassified()[groupIndex][k]
                            && !Double.isNaN(data[k])) {
                          n++;
                        }
                      }

                      if (getBasicPrep().getSamplesToEvaluate()[sampleIndex]
                          && booleanClassifier.getClassified()[groupIndex][sampleIndex]
                          && !Double.isNaN(data[sampleIndex])) {
                        writerSub.print(booleanClassifier.getTitles()[groupIndex] + "_n_" + n + "\t"
                                        + data[sampleIndex]);
                        for (EvaluationResult evaluationResult : evaluationResults) {
                          writerSub.print("\t" + evaluationResult.getEstimateData()[sampleIndex]);
                        }
                        writerSub.println();
                      }

                    }
                  }
                  writerSub.close();
                  RScatter rScatterNumeric = new RScatter(outputBoxSub,
                                                          ext.addToRoot(outputBoxSub,
                                                                        "sub.rscript"),
                                                          ext.removeDirectoryInfo(outputBoxSub) + "sub",
                                                          boxDir + ext.removeDirectoryInfo(outputBoxSub) + "sub.pdf",
                                                          sampleDataStratCats[i],
                                                          new String[] {numericStratCat,
                                                                        pcYsub.get(0),
                                                                        pcYsub.get(pcYsub.size()
                                                                                   - 1)},
                                                          SCATTER_TYPE.BOX, log);
                  rScatterNumeric.setOverWriteExisting(true);
                  rScatterNumeric.setFontsize(12);

                  rScatterNumeric.setxLabel("PC (" + oType + " - sorted)");
                  rScatterNumeric.setTitle(iType + " " + bType);
                  rScatterNumeric.execute();
                  rScatters.add(rScatterNumeric);
                }
              }
            }

          } catch (Exception e) {
            log.reportError("Error writing to " + boxPlot);
            log.reportException(e);
          }
        }
      } catch (FileNotFoundException e) {
        e.printStackTrace();
      }
      // String finalBoxPlot = boxDir + ext.removeDirectoryInfo(outputRoot) + "finalBoxPlot";
      // RScatters rScattersAll = new RScatters(rScatters.toArray(new RScatter[rScatters.size()]),
      // finalBoxPlot + ".rscript", finalBoxPlot + ".pdf", COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1,
      // PLOT_DEVICE.PDF, log);
      scatters.addAll(rScatters);
      // rScattersAll.execute();
      return classifiers;
      // RScatter
    }

    private void plotFirst(String[] sampleDataStratCats, Logger log, String dir,
                           List<RScatter> rScatters, int i, EvaluationResult[] evaluationResults,
                           String outputBox, ArrayList<String> pcYs) {
      RScatter rScatterSubset = new RScatter(outputBox, ext.addToRoot(outputBox, "sub.rscript"),
                                             ext.removeDirectoryInfo(outputBox) + "sub",
                                             dir + ext.removeDirectoryInfo(outputBox) + "sub.pdf",
                                             sampleDataStratCats[i],
                                             ArrayUtils.subArray(ArrayUtils.toStringArray(pcYs), 0,
                                                                 Math.min(20,
                                                                          evaluationResults.length)),
                                             SCATTER_TYPE.BOX, log);
      rScatterSubset.setOverWriteExisting(false);
      rScatterSubset.setyLabel(BOX_Y);
      rScatterSubset.setxLabel("PC (" + oType + " - sorted)");
      rScatterSubset.setTitle(iType + " " + bType);
      rScatterSubset.execute();
      rScatterSubset.setFontsize(12);

      rScatters.add(rScatterSubset);
    }

    private void plotSkips(String[] sampleDataStratCats, Logger log, String dir,
                           List<RScatter> rScatters, int i, EvaluationResult[] evaluationResults,
                           String outputBox, ArrayList<String> pcYs) {
      ArrayList<String> skips = new ArrayList<>();
      int skip = 0;
      while (skip < evaluationResults.length) {
        skips.add(pcYs.get(skip));
        skip += 10;
      }

      RScatter rScatterSkip = new RScatter(outputBox, ext.addToRoot(outputBox, "skips.rscript"),
                                           ext.removeDirectoryInfo(outputBox) + "skips",
                                           dir + ext.removeDirectoryInfo(outputBox) + "skips.pdf",
                                           sampleDataStratCats[i], ArrayUtils.toStringArray(skips),
                                           SCATTER_TYPE.BOX, log);
      rScatterSkip.setOverWriteExisting(false);
      rScatterSkip.setyLabel(BOX_Y);
      rScatterSkip.setFontsize(12);

      rScatterSkip.setxLabel("PC (" + oType + " - sorted)");
      rScatterSkip.setTitle(iType + " " + bType);
      rScatterSkip.execute();
      rScatters.add(rScatterSkip);
    }

    private RScatter plotFirstLast(String[] sampleDataStratCats, Logger log, String dir, int i,
                                   String outputBox, List<String> pcYs) {
      RScatter rScatterFirstLast = new RScatter(outputBox,
                                                ext.addToRoot(outputBox, "subFirstlast.rscript"),
                                                ext.removeDirectoryInfo(outputBox) + "subFirstlast",
                                                dir + ext.removeDirectoryInfo(outputBox) + "subFirstlast.pdf",
                                                sampleDataStratCats[i],
                                                new String[] {pcYs.get(0),
                                                              pcYs.get(pcYs.size() - 1)},
                                                SCATTER_TYPE.BOX, log);
      rScatterFirstLast.setOverWriteExisting(false);
      rScatterFirstLast.setFontsize(12);
      rScatterFirstLast.setxLabel("PC (" + oType + " - sorted)");
      rScatterFirstLast.setTitle(iType + " " + bType);
      rScatterFirstLast.execute();
      return rScatterFirstLast;
    }

    public ITERATION_TYPE getiType() {
      return iType;
    }

    public MODEL_BUILDER_TYPE getbType() {
      return bType;
    }

    public ORDER_TYPE getoType() {
      return oType;
    }

    public RScatter plotSummary(String[] dataColumns, String index,
                                EvaluationResult evaluationResult, int modelN, double[] trimRangeX,
                                double[] trimRangeY, Logger log) {
      RScatter rScatter = new RScatter(outputSummary, evalRscript,
                                       ext.rootOf(outputSummary) + "_" + index,
                                       ext.rootOf(ext.addToRoot(evalPlot, index + ""), false)
                                                                                + ".jpeg",
                                       "Evaluated", dataColumns, SCATTER_TYPE.POINT, log);
      rScatter.setyRange(new double[] {-1, 1});
      rScatter.setxLabel("PC (" + oType + " - sorted)");
      rScatter.setFontsize(12);
      rScatter.setTitle(iType + " " + bType);
      String[] availableYs = rScatter.getrSafeYColumns();

      String[] altYs = new String[availableYs.length];
      String[] titles = ArrayUtils.toStringArray(evaluationResult.getCorrelTitles());
      for (int i = 0; i < availableYs.length; i++) {
        if (availableYs[i].equals("Rsquare_correction")) {
          altYs[i] = "n_" + modelN;
        } else {
          String toSearch = availableYs[i].replaceAll("PEARSON_CORREL_", "")
                                          .replaceAll("SPEARMAN_CORREL_", "");
          int dataIndex = ext.indexOfStr(toSearch, titles);
          if (dataIndex >= 0) {
            String tmp = availableYs[i];
            tmp = tmp.replaceAll("SPEARMAN_CORREL_", "").replaceAll("NO_STRAT", "ALL")
                     .replaceAll("EVAL_DATA_", "");
            System.out.println(evaluationResult.getNumIndsCorrel().size());

            tmp = tmp + "_n_" + evaluationResult.getNumIndsCorrel().get(dataIndex);
            altYs[i] = tmp;
          } else {
            altYs[i] = availableYs[i];
          }
        }
      }
      rScatter.setrSafeAltYColumnNames(altYs);

      rScatter.setxRange(null);
      rScatter.setOverWriteExisting(true);
      rScatter.execute();

      RScatter rScatterTrim = new RScatter(outputSummary, evalRscript,
                                           ext.rootOf(outputSummary) + "_" + index + "_" + index,
                                           ext.rootOf(ext.addToRoot(evalPlot,
                                                                    index + "" + "_" + index),
                                                      false) + ".jpeg",
                                           "Evaluated", dataColumns, SCATTER_TYPE.POINT, log);
      rScatterTrim.setxLabel("PC (" + oType + " - sorted)");
      rScatterTrim.setTitle(iType + " " + bType);
      rScatterTrim.setrSafeAltYColumnNames(altYs);
      if (availableYs.length == 1 && availableYs[0].equals("Rsquare_correction")) {
        rScatterTrim.setyLabel("R SQUARED");
      } else {
        rScatterTrim.setyLabel("R (SPEARMAN)");
      }
      rScatterTrim.setFontsize(12);

      rScatterTrim.setxRange(trimRangeX);
      rScatterTrim.setyRange(trimRangeY);

      rScatterTrim.setOverWriteExisting(true);
      rScatterTrim.execute();
      return rScatter;
    }

    public HeritPlot plotHeritability(Project proj, String pedFile, boolean[] samplesToEvaluate,
                                      double[] otherData, String otherDataTitle, Logger log) {
      heritRscript = outputRoot + ".summary.heritability.rscript";
      heritPlot = outputRoot + ".summary.heritability.pdf";
      heritSummary = outputRoot + ".summary.heritability_summary.parsed.xln";
      String tmpHerit = outputRoot + ".summary.heritability_summary.xln";

      if (pedFile != null) {

        if (!Files.exists(heritSummary)
            || Files.countLines(heritSummary, 0) != Files.countLines(outputSummary, 0)) {
          if (!Files.exists(tmpHerit)
              || Files.countLines(tmpHerit, 0) != Files.countLines(outputSummary, 0)) {

            EvaluationResult.prepareHeritability(proj, pedFile, samplesToEvaluate, outputSer,
                                                 otherData, otherDataTitle);
          }
          try {
            BufferedReader reader = Files.getAppropriateReader(tmpHerit);
            String[] toExtract = new String[] {"Merlin_est.", "Solar_est.", "Solar_p",
                                               "Solar_StdError", "Solar_Kurt", "Solar_KurtWarning"};

            System.out.println(tmpHerit);
            int[] indices = ext.indexFactors(toExtract, Files.getHeaderOfFile(tmpHerit, log), true);
            if (ArrayUtils.countIf(indices, -1) > 0) {
              log.reportError("Could not find " + ArrayUtils.toStr(toExtract) + " in " + tmpHerit);
              return null;
            }
            PrintWriter writer = Files.openAppropriateWriter(heritSummary);
            writer.println(ArrayUtils.toStr(HERIT_ADDITIONS) + "\t" + reader.readLine().trim());
            int index = 0;
            while (reader.ready()) {
              String[] line = reader.readLine().trim().split("\t");
              try {

                double merlinEst = Double.parseDouble(line[indices[0]].replaceAll("%", "")) / 100;
                double solareEst = Double.parseDouble(line[indices[1]]);
                double solareP = Double.parseDouble(line[indices[2]]);
                double solareStError = Double.parseDouble(line[indices[3]]);
                double solarKurt = Double.parseDouble(line[indices[4]]);
                String kurtWarning = line[indices[5]];
                writer.println(index + "\t" + merlinEst + "\t" + solareEst + "\t" + solareP + "\t"
                               + solareStError + "\t" + solarKurt + "\t" + kurtWarning + "\t"
                               + ArrayUtils.toStr(line));
              } catch (NumberFormatException nfe) {
                log.reportTimeWarning("Skipping line " + ArrayUtils.toStr(line)
                                      + " , invalid estimate");
              }
              index++;
            }
            reader.close();
            writer.close();
          } catch (FileNotFoundException fnfe) {
            log.reportError("Error: file \"" + tmpHerit + "\" not found in current directory");
            return null;
          } catch (IOException ioe) {
            log.reportError("Error reading file \"" + tmpHerit + "\"");
            return null;
          } catch (Exception e) {
            log.reportError("Error writing to " + heritSummary);
            log.reportException(e);
          }
        } else {
          log.reportFileExists(heritSummary);
        }
        int numFamIndex = ext.indexOfStr("n_Families_size>1", Files.getHeaderOfFile(tmpHerit, log));
        int numSampsIndex = ext.indexOfStr("n_Samples", Files.getHeaderOfFile(tmpHerit, log));
        int numFam = 0;
        int numSamps = 0;
        double[] estimates = new double[] {-1};

        try {
          String[] line = Files.getHeaderOfFile(tmpHerit, "\t", new String[] {"Model"}, log);
          numFam = Integer.parseInt(line[numFamIndex]);
          numSamps = Integer.parseInt(line[numSampsIndex]);
          double[] tmp = ArrayUtils.toDoubleArray(HashVec.loadFileToStringArray(heritSummary, true,
                                                                                new int[] {1},
                                                                                false));
          if (tmp.length > 0) {
            estimates = tmp;
          }
        } catch (Exception e) {
          log.reportException(e);
        }
        double[] pvals = ArrayUtils.toDoubleArray(HashVec.loadFileToStringArray(heritSummary, true,
                                                                                new int[] {3},
                                                                                false));
        double[] solarHerit = ArrayUtils.toDoubleArray(HashVec.loadFileToStringArray(heritSummary,
                                                                                     true,
                                                                                     new int[] {2},
                                                                                     false));
        double[] solarStError = ArrayUtils.toDoubleArray(HashVec.loadFileToStringArray(heritSummary,
                                                                                       true,
                                                                                       new int[] {4},
                                                                                       false));
        String[] solarWarnKurt = HashVec.loadFileToStringArray(heritSummary, true, new int[] {6},
                                                               false);

        ArrayList<GeomText> sigGTexts = new ArrayList<>();
        for (int i = 0; i < pvals.length; i++) {
          if (pvals[i] < 0.05) {
            GeomText geomText = new GeomText(i, Math.min(solarHerit[i] + solarStError[i], 1), 0,
                                             "*", 12);
            sigGTexts.add(geomText);
          }
          if (solarWarnKurt[i].equals("NORMAL")) {
            GeomText geomText = new GeomText(i, Math.max(solarHerit[i] - solarStError[i], 0), 0,
                                             "N", 12);
            sigGTexts.add(geomText);
          }
        }

        RScatter rScatter = new RScatter(heritSummary, heritRscript, ext.rootOf(heritSummary),
                                         heritPlot, HERIT_ADDITIONS[0],
                                         new String[] {HERIT_ADDITIONS[1], HERIT_ADDITIONS[2]},
                                         SCATTER_TYPE.POINT, log);
        ErrorBars errorBars = new ErrorBars(new String[] {HERIT_ADDITIONS[2]},
                                            new String[] {HERIT_ADDITIONS[4]});
        rScatter.setErrorBars(errorBars);
        rScatter.setyRange(new double[] {0, 1});
        rScatter.setxLabel("PC (" + oType + " - sorted)");
        rScatter.setTitle(iType + " " + bType + "; nInd=" + numSamps + ", nFam=" + numFam);
        rScatter.setyLabel("Proportion Heritability");
        rScatter.setFontsize(12);
        rScatter.setOverWriteExisting(true);
        rScatter.setgPoint_SIZE(GEOM_POINT_SIZE.GEOM_POINT);
        rScatter.setOverWriteExisting(true);
        rScatter.setgTexts(sigGTexts.toArray(new GeomText[sigGTexts.size()]));
        rScatter.execute();
        return new HeritPlot(rScatter, numSamps, numFam, estimates);
      }

      return null;
    }

    public void setValid(boolean valid) {}

    public String getOutputSer() {
      return outputSer;
    }

    public void setOutputSer(String outputSer) {
      this.outputSer = outputSer;
    }

    public void setOutputRoot(String outputRoot) {
      this.outputRoot = outputRoot;
    }

    public String getOutputRank() {
      return outputRank;
    }

    public String getOutputSummary() {
      return outputSummary;
    }

  }

  private static CorrectionIterator[] getIterations(Project proj, String markesToEvaluate,
                                                    String samplesToBuildModels, String outputDir,
                                                    LS_TYPE lType, boolean recomputeLRR,
                                                    double pcPercent, int numthreads) {
    ArrayList<CorrectionIterator> cIterators = new ArrayList<>();
    // System.out.println("JDOFJSDF remember the pcs");
    for (int i = 0; i < ITERATION_TYPE.values().length; i++) {
      for (int j = 0; j < ORDER_TYPE.values().length; j++) {
        for (int j2 = 0; j2 < MODEL_BUILDER_TYPE.values().length; j2++) {
          cIterators.add(new CorrectionIterator(proj, markesToEvaluate, samplesToBuildModels,
                                                ITERATION_TYPE.values()[i], ORDER_TYPE.values()[j],
                                                MODEL_BUILDER_TYPE.values()[j2], outputDir, lType,
                                                recomputeLRR, pcPercent, numthreads));
        }
      }
    }
    return cIterators.toArray(new CorrectionIterator[cIterators.size()]);
  }

  private static class HeritPlot {

    private final RScatter rScatter;
    private final int numSamps;
    private final int numFam;
    private final double[] estimates;

    public HeritPlot(RScatter rScatter, int numSamps, int numFam, double[] estimates) {
      super();
      this.rScatter = rScatter;
      this.numSamps = numSamps;
      this.numFam = numFam;
      this.estimates = estimates;
    }

    public RScatter getrScatter() {
      return rScatter;
    }

    public double[] getEstimates() {
      return estimates;
    }

    public int getNumSamps() {
      return numSamps;
    }

    public int getNumFam() {
      return numFam;
    }

  }

  public static CorrectionIterator[] runAll(Project proj, String markesToEvaluate,
                                            String samplesToBuildModels, String outputDir,
                                            String pcFile, String pedFile, LS_TYPE lType,
                                            boolean recomputeLRR, double pcPercent, boolean plot,
                                            int numthreads) {
    if (pcFile != null) {
      proj.INTENSITY_PC_FILENAME.setValue(pcFile);
    } else {}
    if (outputDir == null) {
      outputDir = proj.PROJECT_DIRECTORY.getValue()
                  + ext.rootOf(proj.INTENSITY_PC_FILENAME.getValue()) + "_eval/";
    }
    if (!Files.exists(proj.MARKER_DATA_DIRECTORY.getValue() + "markers.0.mdRAF")) {
      proj.getLog().reportTimeWarning("Did not see " + proj.MARKER_DATA_DIRECTORY.getValue()
                                      + "markers.0.mdRAF, attempting to transpose now");
      TransposeData.transposeData(proj, 2000000000, false);
    }
    new File(outputDir).mkdirs();
    String lrrSdCurrent = outputDir + "lrrSD.txt";
    proj.SAMPLE_QC_FILENAME.setValue(lrrSdCurrent);

    CorrectionIterator[] cIterators = getIterations(proj, markesToEvaluate, samplesToBuildModels,
                                                    outputDir, lType, recomputeLRR, pcPercent,
                                                    numthreads);
    ArrayList<RScatter> rScatters = new ArrayList<>();

    for (CorrectionIterator cIterator : cIterators) {
      cIterator.run();

    }

    String[] customPlotFiles = Files.list(proj.PROJECT_DIRECTORY.getValue(), null,
                                          ".pc_evaluation.titles.txt", false, true);

    ArrayList<String[]> plotters = new ArrayList<>();
    for (String customPlotFile : customPlotFiles) {
      String[] groups = HashVec.loadFileToStringArray(customPlotFile, true, new int[] {0}, false);
      String[] names = HashVec.loadFileToStringArray(customPlotFile, true, new int[] {1}, false);
      List<String>[] arrayStringList = new List[ArrayUtils.unique(groups).length];
      for (int j = 0; j < arrayStringList.length; j++) {
        arrayStringList[j] = new ArrayList<>(10);
      }
      for (int j = 0; j < groups.length; j++) {
        arrayStringList[Integer.parseInt(groups[j]) - 1].add("SPEARMAN_CORREL_EVAL_DATA_" + names[j]
                                                             + "_CUSTOM_PHENO_TAG_NO_STRAT");
      }
      for (List<String> element : arrayStringList) {
        plotters.add(ArrayUtils.toStringArray(element));
      }
    }

    // For making our basic plots
    // String[] plotTitlesForMain = new String[] { "Rsquare_correction",
    // "ICC_EVAL_CLASS_DUPLICATE_ALL", "ICC_EVAL_CLASS_DUPLICATE_SAME_VISIT", "ICC_EVAL_CLASS_FC",
    // "SPEARMAN_CORREL_AGE", "SPEARMAN_CORREL_EVAL_DATA_SEX",
    // "SPEARMAN_CORREL_EVAL_DATA_resid.mtDNaN.qPCR.MT001",
    // "SPEARMAN_CORREL_EVAL_DATA_resid.mtDNA.qPCR",
    // "SPEARMAN_CORREL_EVAL_DATA_Mt_DNA_relative_copy_number",
    // "SPEARMAN_CORREL_EVAL_DATA_Ratio.ND1", "SPEARMAN_CORREL_EVAL_DATA_qpcr.qnorm.exprs" };
    // String[] plotTitlesForMito = new String[] { "Rsquare_correction",
    // "ICC_EVAL_CLASS_FC_NO_STRAT", "SPEARMAN_CORREL_AGE_NO_STRATs",
    // "SPEARMAN_CORREL_EVAL_DATA_Mt_DNA_relative_copy_number_NO_STRAT",
    // "SPEARMAN_CORREL_EVAL_DATA_SEX_NO_STRATs" };
    //
    // String[] plotTitlesForMitoFC = new String[] {
    // "SPEARMAN_CORREL_EVAL_DATA_Mt_DNA_relative_copy_number_NO_STRAT",
    // "SPEARMAN_CORREL_EVAL_DATA_Mt_DNA_relative_copy_number_PT",
    // "SPEARMAN_CORREL_EVAL_DATA_Mt_DNA_relative_copy_number_BU",
    // "SPEARMAN_CORREL_EVAL_DATA_Mt_DNA_relative_copy_number_NY",
    // "SPEARMAN_CORREL_EVAL_DATA_Mt_DNA_relative_copy_number_DK" };
    // String[] exomeMito = new String[] { "SPEARMAN_CORREL_EVAL_DATA_qpcr.qnorm.exprs_NO_STRAT",
    // "SPEARMAN_CORREL_EVAL_DATA_qpcr.qnorm.exprs_M",
    // "SPEARMAN_CORREL_EVAL_DATA_qpcr.qnorm.exprs_W",
    // "SPEARMAN_CORREL_EVAL_DATA_qpcr.qnorm.exprs_F" };
    // plotTitlesForMitoFC = Array.concatAll(plotTitlesForMitoFC, exomeMito);
    //
    // String[] plotTitlesForMitoAge = new String[] { "SPEARMAN_CORREL_AGE_NO_STRAT",
    // "SPEARMAN_CORREL_AGE_PT", "SPEARMAN_CORREL_AGE_BU", "SPEARMAN_CORREL_AGE_NY",
    // "SPEARMAN_CORREL_AGE_DK", "SPEARMAN_CORREL_AGE_M", "SPEARMAN_CORREL_AGE_W",
    // "SPEARMAN_CORREL_AGE_F", "SPEARMAN_CORREL_AGE_J" };
    // String[] plotTitlesForMitoSex = new String[] { "SPEARMAN_CORREL_EVAL_DATA_SEX_NO_STRAT",
    // "SPEARMAN_CORREL_EVAL_DATA_SEX_PT", "SPEARMAN_CORREL_EVAL_DATA_SEX_BU",
    // "SPEARMAN_CORREL_EVAL_DATA_SEX_NY", "SPEARMAN_CORREL_EVAL_DATA_SEX_DK",
    // "SPEARMAN_CORREL_EVAL_DATA_SEX_M", "SPEARMAN_CORREL_EVAL_DATA_SEX_W",
    // "SPEARMAN_CORREL_EVAL_DATA_SEX_F", "SPEARMAN_CORREL_EVAL_DATA_SEX_J" };
    // String[] quickTitles = new String[] { "SPEARMAN_CORREL_AGE_NO_STRAT",
    // "SPEARMAN_CORREL_EVAL_DATA_SEX_NO_STRAT",
    // "SPEARMAN_CORREL_EVAL_DATA_qpcr.qnorm.exprs_NO_STRAT",
    // "SPEARMAN_CORREL_EVAL_DATA_Mt_DNA_relative_copy_number_NO_STRAT" };

    // String[][] plotTitlesForSummary = new String[][] { plotTitlesForMain, plotTitlesForMito,
    // plotTitlesForMitoFC, plotTitlesForMitoAge, plotTitlesForMitoSex, quickTitles };

    String[][] plotTitlesForSummary = new String[0][];

    String[][] fileLoad = plotters.toArray(new String[plotters.size()][]);
    String[][] tmp = new String[plotTitlesForSummary.length + fileLoad.length][];

    for (int i = 0; i < plotTitlesForSummary.length; i++) {
      tmp[i] = plotTitlesForSummary[i];
    }
    for (int i = 0; i < fileLoad.length; i++) {
      tmp[i + plotTitlesForSummary.length] = fileLoad[i];
    }
    plotTitlesForSummary = tmp;

    // String[] subsetDataHeritability = new String[] { "EVAL_DATA_Mt_DNA_relative_copy_number" };
    //
    // String[] numericStratCats = new String[] { "EVAL_DATA_Mt_DNA_relative_copy_number",
    // "EVAL_DATA_qpcr.qnorm.exprs" };
    // IterSummaryProducer producer = new IterSummaryProducer(proj, cIterators,
    // plotTitlesForSummary, CorrectionEvaluator.INDEPS_CATS, numericStratCats, pedFile, new
    // double[] { 0, 150 }, new double[] { -.5, 1 }, plot);

    IterSummaryProducer producer = new IterSummaryProducer(proj, cIterators, plotTitlesForSummary,
                                                           CorrectionEvaluator.INDEPS_CATS, null,
                                                           pedFile, new double[] {0, 150},
                                                           new double[] {-.5, 1}, plot);

    // if (pedFile != null) {
    // producer.setSubsetDataHeritability(subsetDataHeritability);
    // }

    try (WorkerTrain<IterSummary> summaryTrain = new WorkerTrain<>(producer, numthreads, numthreads,
                                                                   proj.getLog())) {
      while (summaryTrain.hasNext()) {
        IterSummary iterSummary = summaryTrain.next();
        RScatter[] rScattersTmp = iterSummary.getrScatters();
        for (RScatter element : rScattersTmp) {
          rScatters.add(element);
        }
      }
    }

    String outputRoot = outputDir + "finalSummary";

    Logger log = proj.getLog();
    String finalCompFile = outputDir + "/typed/" + "evals.finalComp.txt";
    try {
      PrintWriter writer = Files.openAppropriateWriter(finalCompFile);
      ArrayList<String> header = new ArrayList<>();
      header.add("MODEL_BUILDER_TYPE");
      header.add("ORDER_TYPE");
      header.add("ITERATION_TYPE");
      header.add("NUM_SAMPLES_MODEL");
      header.add("NUM_SAMPLES_EVALUATED");
      header.add("NUM_TOTAL_PCS");

      header.add("STAT");
      header.add("EVALUATED");
      header.add("MIN_STAT_PC");
      header.add("MIN_STAT");
      header.add("MAX_STAT_PC");
      header.add("MAX_STAT");
      header.add("MIN_PVAL_PC");
      header.add("MIN_PVAL");
      header.add("MAX_PVAL_PC");
      header.add("MAX_PVAL");
      header.add("AVERAGE_STAT");
      header.add("MEDIAN_STAT");
      header.add("STD_STAT");
      header.add("PC_0_STAT");
      header.add("PC_0_PVAL");
      header.add("PC_15_STAT");
      header.add("PC_15_PVAL");
      header.add("PC_50_STAT");
      header.add("PC_50_PVAL");
      header.add("PC_150_STAT");
      header.add("PC_150_PVAL");
      header.add("PC_MAX_STAT");
      header.add("PC_MAX_PVAL");
      writer.println(ArrayUtils.toStr(ArrayUtils.toStringArray(header)));
      for (CorrectionIterator correctionIterator : cIterators) {
        EvaluationResult[] evaluationResults = EvaluationResult.readSerial(correctionIterator.getIterationResult()
                                                                                             .getOutputSer(),
                                                                           log);
        if (evaluationResults.length > 0) {
          MODEL_BUILDER_TYPE mBuilder_TYPE = correctionIterator.getIterationResult().getbType();
          ORDER_TYPE oType = correctionIterator.getIterationResult().getoType();
          ITERATION_TYPE iType = correctionIterator.getIterationResult().getiType();
          ArrayList<String> names = evaluationResults[0].getCorrelTitles();
          ArrayList<double[]> statsSpear = evaluationResults[0].getSpearmanCorrel();
          printStatSummary(writer, correctionIterator, "SPEARMAN", evaluationResults, mBuilder_TYPE,
                           oType, iType, statsSpear, null, names);
          ArrayList<double[]> statsPear = evaluationResults[0].getPearsonCorrels();
          printStatSummary(writer, correctionIterator, "PEARSON", evaluationResults, mBuilder_TYPE,
                           oType, iType, statsPear, null, names);
          printStatSummary(writer, correctionIterator, "ICC", evaluationResults, mBuilder_TYPE,
                           oType, iType, statsPear, evaluationResults[0].getIccs(),
                           evaluationResults[0].getIccTitles());
        }
      }

      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + finalCompFile);
      log.reportException(e);
    }
    if (rScatters.size() > 0) {
      RScatters finalScatters = new RScatters(rScatters.toArray(new RScatter[rScatters.size()]),
                                              outputRoot + ".rscript", outputRoot + ".pdf",
                                              COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_2,
                                              PLOT_DEVICE.PDF, proj.getLog());

      finalScatters.execute();
    }
    return cIterators;
  }

  private static void printStatSummary(PrintWriter writer, CorrectionIterator correctionIterator,
                                       String type, EvaluationResult[] evaluationResults,
                                       MODEL_BUILDER_TYPE mBuilder_TYPE, ORDER_TYPE oType,
                                       ITERATION_TYPE iType, List<double[]> statsAL,
                                       ArrayList<ICC> statsICC, ArrayList<String> name) {

    for (int i = 0; i < (statsAL == null ? statsICC.size() : statsAL.size()); i++) {
      double[] stats = new double[evaluationResults.length];
      double[] pval = new double[evaluationResults.length];
      int numSamps = evaluationResults[0].getNumIndsCorrel().get(i);
      double stat0 = Double.NaN;
      double pval0 = Double.NaN;
      double stat50 = Double.NaN;
      double pval50 = Double.NaN;
      double stat15 = Double.NaN;
      double pval15 = Double.NaN;
      double stat150 = Double.NaN;
      double pval150 = Double.NaN;
      double statMax = Double.NaN;
      double pvalMax = Double.NaN;

      for (int j = 0; j < evaluationResults.length; j++) {
        if (type.equals("SPEARMAN")) {
          stats[j] = evaluationResults[j].getSpearmanCorrel().get(i)[0];
          pval[j] = evaluationResults[j].getSpearmanCorrel().get(i)[1];
          if (j == 0) {
            stat0 = stats[j];
            pval0 = pval[j];
          }
          if (j == 15) {
            stat15 = stats[j];
            pval15 = pval[j];
          }
          if (j == 50) {
            stat50 = stats[j];
            pval50 = pval[j];
          }
          if (j == 150) {
            stat150 = stats[j];
            pval150 = pval[j];
          }
          if (j == evaluationResults.length - 1) {
            statMax = stats[j];
            pvalMax = pval[j];
          }
        } else if (type.equals("PEARSON")) {
          stats[j] = evaluationResults[j].getPearsonCorrels().get(i)[0];
          pval[j] = evaluationResults[j].getPearsonCorrels().get(i)[1];
          if (j == 0) {
            stat0 = stats[j];
            pval0 = pval[j];
          }
          if (j == 15) {
            stat15 = stats[j];
            pval15 = pval[j];
          }
          if (j == 50) {
            stat50 = stats[j];
            pval50 = pval[j];
          }
          if (j == 150) {
            stat150 = stats[j];
            pval150 = pval[j];
          }
          if (j == evaluationResults.length - 1) {
            statMax = stats[j];
            pvalMax = pval[j];
          }
        } else if (type.equals("ICC")) {
          stats[j] = evaluationResults[j].getIccs().get(i).getICC();
          pval[j] = Double.NaN;
          if (j == 0) {
            stat0 = stats[j];
            pval0 = pval[j];
          }
          if (j == 15) {
            stat15 = stats[j];
            pval15 = pval[j];
          }
          if (j == 50) {
            stat50 = stats[j];
            pval50 = pval[j];
          }
          if (j == 150) {
            stat150 = stats[j];
            pval150 = pval[j];
          }
          if (j == evaluationResults.length - 1) {
            statMax = stats[j];
            pvalMax = pval[j];
          }
        }
        if (evaluationResults[j].getNumIndsCorrel().get(i) != numSamps) {
          System.err.println("Error mismatched number of samples");
          System.exit(1);
        }
      }

      String evalName = name.get(i);
      // if ((evalName.contains("qpcr.qnorm.exprs") || evalName.toLowerCase().contains("age") ||
      // evalName.toLowerCase().contains("center") || evalName.toLowerCase().contains("dupli") ||
      // evalName.toLowerCase().contains("sex"))) {
      double maxStat = ArrayUtils.max(stats);
      double minStat = ArrayUtils.min(stats);
      int maxStatPC = ArrayUtils.maxIndex(stats);
      int minStatPC = ArrayUtils.minIndex(stats);

      double maxPval = ArrayUtils.max(pval);
      double minPval = ArrayUtils.min(pval);
      int maxPvalPC = ArrayUtils.maxIndex(pval);
      int minPvalPC = ArrayUtils.minIndex(pval);

      double avgStat = ArrayUtils.mean(stats);
      double medianStat = ArrayUtils.mean(stats);
      double stdvStat = ArrayUtils.stdev(stats);

      int numSamplesModel = ArrayUtils.booleanArraySum(correctionIterator.getIterationResult()
                                                                         .getBasicPrep()
                                                                         .getSamplesForModels());
      ArrayList<String> result = new ArrayList<>();
      result.add(mBuilder_TYPE.toString());
      result.add(oType.toString());
      result.add(iType.toString());
      result.add(numSamplesModel + "");
      result.add(numSamps + "");
      result.add(evaluationResults.length - 1 + "");
      result.add(type);

      result.add(evalName);
      result.add(minStatPC + "");
      result.add(minStat + "");
      result.add(maxStatPC + "");
      result.add(maxStat + "");
      result.add(minPvalPC + "");
      result.add(minPval + "");
      result.add(maxPvalPC + "");
      result.add(maxPval + "");
      result.add(avgStat + "");
      result.add(medianStat + "");
      result.add(stdvStat + "");
      result.add(stat0 + "");
      result.add(pval0 + "");
      result.add(stat15 + "");
      result.add(pval15 + "");
      result.add(stat50 + "");
      result.add(pval50 + "");
      result.add(stat150 + "");
      result.add(pval150 + "");
      result.add(statMax + "");
      result.add(pvalMax + "");

      writer.println(ArrayUtils.toStr(ArrayUtils.toStringArray(result)));
    }
    // }
  }

  private static class IterSummary {

    private final RScatter[] rScatters;

    public RScatter[] getrScatters() {
      return rScatters;
    }

    public IterSummary(RScatter[] rScatters, MODEL_BUILDER_TYPE mBuilder_TYPE, ITERATION_TYPE iType,
                       ORDER_TYPE otType) {
      super();
      this.rScatters = rScatters;
    }

  }

  public static class IterSummaryProducer extends AbstractProducer<IterSummary> {

    private final Project proj;
    private final CorrectionIterator[] cIterators;
    private final String[][] plotTitlesForSummary;
    private String[] subsetDataHeritability;
    private final String[] stratCats, numericStratCats;
    private final double[] trimRangeX;
    private final double[] trimRangeY;
    private final String pedFile;
    private final boolean plot;
    private int index;

    public IterSummaryProducer(Project proj, CorrectionIterator[] cIterators,
                               String[][] plotTitlesForSummary, String[] stratCats,
                               String[] numericStratCats, String pedFile, double[] trimRangeX,
                               double[] trimRangeY, boolean plot) {
      super();
      this.proj = proj;
      this.cIterators = cIterators;
      this.pedFile = pedFile;
      this.plotTitlesForSummary = plotTitlesForSummary;
      index = 0;
      this.stratCats = stratCats;
      this.numericStratCats = numericStratCats;
      this.trimRangeX = trimRangeX;
      this.trimRangeY = trimRangeY;
      this.plot = plot;
    }

    public void setSubsetDataHeritability(String[] subsetDataHeritability) {
      this.subsetDataHeritability = subsetDataHeritability;
    }

    @Override
    public boolean hasNext() {
      // TODO Auto-generated method stub
      return index < cIterators.length;
    }

    @Override
    public Callable<IterSummary> next() {
      final CorrectionIterator tmp = cIterators[index];
      final Logger log = proj.getLog();
      Callable<IterSummary> callable = new Callable<IterSummary>() {

        @Override
        public IterSummary call() throws Exception {
          IterationResult iterationResult = tmp.getIterationResult();

          String originalSer = iterationResult.getOutputSer();
          BooleanClassifier[] classifiers = null;
          String typedOut = ext.parseDirectoryOfFile(originalSer) + FINAL_EST_DIR;
          String outputRoot = typedOut + iterationResult.getbType() + "_"
                              + iterationResult.getoType() + "_" + iterationResult.getiType()
                              + "_finalSummary";
          new File(typedOut).mkdirs();
          EvaluationResult[] evaluationResults = EvaluationResult.readSerial(originalSer,
                                                                             proj.getLog());
          ArrayList<RScatter> scatters = new ArrayList<>();
          if (plot) {
            if (stratCats != null && numericStratCats != null) {
              try {
                classifiers = iterationResult.estimateBoxPlots(proj, stratCats, numericStratCats,
                                                               scatters, proj.getLog());
              } catch (Exception e) {

              }

            }

            int modelN = ArrayUtils.booleanArraySum(iterationResult.getBasicPrep()
                                                                   .getSamplesForModels());
            for (int i = 0; i < plotTitlesForSummary.length; i++) {
              RScatter rScatterSummary = iterationResult.plotSummary(plotTitlesForSummary[i],
                                                                     i + "", evaluationResults[0],
                                                                     modelN, trimRangeX, trimRangeY,
                                                                     proj.getLog());
              rScatterSummary.setSeriesLabeler(null);
              rScatterSummary.setgTexts(null);
              rScatterSummary.execute();

              scatters.add(rScatterSummary);
            }

            String[] fullHeader = Files.getHeaderOfFile(iterationResult.getOutputSummary(),
                                                        proj.getLog());
            for (String element : fullHeader) {
              if (element.startsWith("SPEARMAN_CORREL_EVAL_DATA_")
                  && element.endsWith("NO_STRAT")) {
                String title = element.replaceAll("SPEARMAN_CORREL_EVAL_DATA_", "")
                                      .replaceAll("NO_STRAT", "");
                RScatter rScatterSummary = iterationResult.plotSummary(new String[] {element},
                                                                       title, evaluationResults[0],
                                                                       modelN, trimRangeX,
                                                                       trimRangeY, proj.getLog());
                rScatterSummary.setSeriesLabeler(null);
                rScatterSummary.setgTexts(null);
                rScatterSummary.execute();
              }
            }
            if (pedFile != null) {
              classifyHerit(tmp, log, iterationResult, classifiers, scatters, null, null, null);
              if (subsetDataHeritability == null) {

                scatters.add(iterationResult.plotHeritability(proj, pedFile,
                                                              iterationResult.getBasicPrep()
                                                                             .getSamplesToEvaluate(),
                                                              null, null, log)
                                            .getrScatter());
                scatters.add(iterationResult.plotRank(log));
              } else {
                if (iterationResult.getBasicPrep() == null) {
                  log.reportError("must have data basic prep object");
                  return null;
                } else {

                  scatters.add(iterationResult.plotHeritability(proj, pedFile,
                                                                iterationResult.getBasicPrep()
                                                                               .getSamplesToEvaluate(),
                                                                null, null, log)
                                              .getrScatter());

                  for (int i = 0; i < subsetDataHeritability.length; i++) {
                    if (!subsetDataHeritability[i].toLowerCase().contains("sex")) {

                      IterationResult tmpR = tmp.getIterationResult();
                      tmpR.setOutputSer(originalSer);
                      int index = ext.indexOfStr(subsetDataHeritability[i],
                                                 iterationResult.getBasicPrep().getNumericTitles());
                      if (index >= 0) {

                        double[] data = iterationResult.getBasicPrep().getNumericData()[index];
                        boolean[] currentModel = new boolean[iterationResult.getBasicPrep()
                                                                            .getSamplesToEvaluate().length];
                        for (int j = 0; j < currentModel.length; j++) {
                          if (iterationResult.getBasicPrep().getSamplesToEvaluate()[j]
                              && !Double.isNaN(data[j])) {
                            currentModel[j] = true;
                          } else {
                            currentModel[j] = false;
                          }
                        }
                        classifyHerit(tmp, log, iterationResult, classifiers, scatters,
                                      subsetDataHeritability[i], currentModel, data);

                        String newDir = ext.parseDirectoryOfFile(originalSer)
                                        + subsetDataHeritability[i] + "_herit/";
                        new File(newDir).mkdirs();
                        String newSer = newDir
                                        + ext.removeDirectoryInfo(ext.addToRoot(originalSer,
                                                                                "." + subsetDataHeritability[i]
                                                                                             + ".summary"));
                        Files.copyFileUsingFileChannels(new File(originalSer), new File(newSer),
                                                        log);
                        tmpR.setOutputRoot(ext.rootOf(ext.rootOf(newSer, false), false));
                        tmpR.setOutputSer(newSer);
                        RScatter tmpherit = tmpR.plotHeritability(proj, pedFile, currentModel, null,
                                                                  null, log)
                                                .getrScatter();
                        tmpherit.setTitle(subsetDataHeritability[i] + " - matched samples n="
                                          + ArrayUtils.booleanArraySum(currentModel));
                        scatters.add(tmpR.plotHeritability(proj, pedFile, currentModel, null, null,
                                                           log)
                                         .getrScatter());
                      } else {
                        log.reportTimeWarning("Skipping " + subsetDataHeritability[i]);
                      }
                    }
                  }
                  // TODO this whole thing needs to be refactored, what a garble
                  scatters.add(iterationResult.plotRank(log));
                }
              }

            } else {
              scatters.add(iterationResult.plotRank(log));
            }
            RScatters finalScattersTyped = new RScatters(scatters.toArray(new RScatter[scatters.size()]),
                                                         outputRoot + ".rscript",
                                                         outputRoot + ".pdf",
                                                         COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_2,
                                                         PLOT_DEVICE.PDF, proj.getLog());

            finalScattersTyped.execute();
          }

          String outAllEstimates = outputRoot + ".estimates.txt.gz";

          // System.out.println(outAllEstimates);
          // System.out.println(evaluationResults.length);

          String[] samples = proj.getSamples();
          try {
            PrintWriter writer = Files.getAppropriateWriter(outAllEstimates);
            writer.print("DNA");
            for (int j = 0; j < evaluationResults.length; j++) {
              writer.print("\tPC" + j);
            }
            writer.println();
            for (int j = 0; j < samples.length; j++) {
              writer.print(samples[j]);
              for (EvaluationResult evaluationResult : evaluationResults) {
                writer.print("\t" + evaluationResult.getEstimateData()[j]);
              }
              writer.println();
            }
            writer.close();
          } catch (Exception e) {
            log.reportError("Error writing to " + outAllEstimates);
            log.reportException(e);
          }
          return new IterSummary(scatters.toArray(new RScatter[scatters.size()]),
                                 iterationResult.getbType(), iterationResult.getiType(),
                                 iterationResult.getoType());
        }

        private void classifyHerit(final CorrectionIterator tmp, final Logger log,
                                   IterationResult iterationResult, BooleanClassifier[] classifiers,
                                   List<RScatter> scatters, String subsetDataHeritability,
                                   boolean[] subsetMask, double[] otherData) {
          if (classifiers != null) {
            String originalSer = iterationResult.getOutputSer();
            for (int classifyIndex = 0; classifyIndex < classifiers.length; classifyIndex++) {
              if (classifiers[classifyIndex].getTitles().length < 10
                  && !stratCats[classifyIndex].contains("sex")) {

                try {
                  Thread.sleep(1000);
                } catch (InterruptedException ie) {}
                ArrayList<String> currentClass = new ArrayList<>();
                ArrayList<HeritPlot> hPlots = new ArrayList<>();
                ArrayList<GeomText> otherDataGeomTexts = new ArrayList<>();
                boolean[] allModel = ArrayUtils.booleanArray(iterationResult.getBasicPrep()
                                                                            .getSamplesToEvaluate().length,
                                                             false);
                for (int titleIndex = 0; titleIndex < classifiers[classifyIndex].getTitles().length; titleIndex++) {
                  IterationResult tmpR = tmp.getIterationResult();
                  String newDir = ext.parseDirectoryOfFile(originalSer) + stratCats[classifyIndex]
                                  + "_herit"
                                  + (subsetDataHeritability == null ? "" : subsetDataHeritability)
                                  + "/";
                  new File(newDir).mkdirs();
                  String newSer = newDir
                                  + ext.removeDirectoryInfo(ext.addToRoot(originalSer,
                                                                          "." + classifiers[classifyIndex].getTitles()[titleIndex]
                                                                                       + (subsetDataHeritability == null ? ""
                                                                                                                         : subsetDataHeritability)
                                                                                       + ".summary"));
                  Files.copyFileUsingFileChannels(new File(iterationResult.getOutputSer()),
                                                  new File(newSer), log);
                  tmpR.setOutputRoot(ext.rootOf(ext.rootOf(newSer, false), false));
                  tmpR.setOutputSer(newSer);

                  boolean[] currentModel = new boolean[iterationResult.getBasicPrep()
                                                                      .getSamplesToEvaluate().length];
                  for (int j = 0; j < currentModel.length; j++) {
                    if (iterationResult.getBasicPrep().getSamplesToEvaluate()[j]
                        && classifiers[classifyIndex].getClassified()[titleIndex][j]
                        && (subsetMask == null || subsetMask[j])) {
                      currentModel[j] = true;
                      allModel[j] = true;
                    } else {
                      currentModel[j] = false;
                    }
                  }
                  HeritPlot heritPlot = tmpR.plotHeritability(proj, pedFile, currentModel, null,
                                                              null, log);
                  RScatter tmpherit = heritPlot.getrScatter();
                  tmpherit.setTitle(classifiers[classifyIndex].getTitles()[titleIndex]
                                    + " - matched samples n="
                                    + ArrayUtils.booleanArraySum(currentModel));
                  tmpherit.execute();
                  GeomText[] gTexts = tmpherit.getgTexts();
                  for (GeomText gText : gTexts) {
                    otherDataGeomTexts.add(gText);
                  }
                  // scatters.add(tmpR.plotHeritability(proj, pedFile, currentModel, log));
                  hPlots.add(heritPlot);
                  currentClass.add(tmpherit.getDataFile());

                  if (otherData != null) {
                    newSer = newDir
                             + ext.removeDirectoryInfo(ext.addToRoot(originalSer,
                                                                     ".raw" + classifiers[classifyIndex].getTitles()[titleIndex]
                                                                                  + (subsetDataHeritability == null ? ""
                                                                                                                    : subsetDataHeritability)
                                                                                  + ".summary"));
                    Files.copyFileUsingFileChannels(new File(iterationResult.getOutputSer()),
                                                    new File(newSer), log);
                    tmpR.setOutputRoot(ext.rootOf(ext.rootOf(newSer, false), false));
                    tmpR.setOutputSer(newSer);

                    HeritPlot heritPlotRaw = tmpR.plotHeritability(proj, pedFile, currentModel,
                                                                   otherData,
                                                                   subsetDataHeritability, log);

                    GeomText geomTextRaw = new GeomText(0, heritPlotRaw.getEstimates()[0], 0,
                                                        stratCats[classifyIndex] + classifiers[classifyIndex].getTitles()[titleIndex]
                                                                                              + "_qpcr",
                                                        3);
                    otherDataGeomTexts.add(geomTextRaw);
                    // System.out.println(geomTextRaw.getCommand());
                  }

                }

                if (otherData != null) {
                  String newDir = ext.parseDirectoryOfFile(originalSer) + stratCats[classifyIndex]
                                  + "_herit"
                                  + (subsetDataHeritability == null ? "" : subsetDataHeritability)
                                  + "/";
                  IterationResult tmpFullModelRSub = tmp.getIterationResult();

                  String newSer = newDir
                                  + ext.removeDirectoryInfo(ext.addToRoot(originalSer,
                                                                          ".rawAll" + (subsetDataHeritability == null ? ""
                                                                                                                      : subsetDataHeritability)
                                                                                       + ".summary"));
                  Files.copyFileUsingFileChannels(new File(iterationResult.getOutputSer()),
                                                  new File(newSer), log);
                  tmpFullModelRSub.setOutputRoot(ext.rootOf(ext.rootOf(newSer, false), false));
                  tmpFullModelRSub.setOutputSer(newSer);

                  HeritPlot heritPlotRaw = tmpFullModelRSub.plotHeritability(proj, pedFile,
                                                                             allModel, otherData,
                                                                             subsetDataHeritability,
                                                                             log);

                  GeomText geomTextRaw = new GeomText(0, heritPlotRaw.getEstimates()[0], 0,
                                                      stratCats[classifyIndex] + "_qpcr", 3);
                  otherDataGeomTexts.add(geomTextRaw);
                  // System.out.println(geomTextRaw.getCommand());
                }

                String stratCatHerit = stratCats[classifyIndex]
                                       + (subsetDataHeritability == null ? ""
                                                                         : subsetDataHeritability);
                String comboHerit = ext.rootOf(originalSer, false) + stratCatHerit
                                    + ".heritability.summary";
                String[] adds = new String[hPlots.size() + 1];

                IterationResult tmpFullModelR = tmp.getIterationResult();
                String newDir = ext.parseDirectoryOfFile(originalSer) + stratCats[classifyIndex]
                                + "_herit"
                                + (subsetDataHeritability == null ? "" : subsetDataHeritability)
                                + "/";
                new File(newDir).mkdirs();
                String newSer = newDir
                                + ext.removeDirectoryInfo(ext.addToRoot(originalSer,
                                                                        "." + stratCatHerit
                                                                                     + ".summary"));

                Files.copyFileUsingFileChannels(new File(iterationResult.getOutputSer()),
                                                new File(newSer), log);

                tmpFullModelR.setOutputRoot(ext.rootOf(ext.rootOf(newSer, false), false));
                tmpFullModelR.setOutputSer(newSer);
                HeritPlot all = tmpFullModelR.plotHeritability(proj, pedFile, allModel, null, null,
                                                               log);

                adds[adds.length - 1] = stratCatHerit + "_ALL_nInd_" + all.getNumSamps()
                                        + "_numFam_" + all.getNumFam();
                currentClass.add(all.getrScatter().getDataFile());
                for (int i = 0; i < adds.length - 1; i++) {
                  String add = classifiers[classifyIndex].getTitles()[i] + "_nInd_"
                               + hPlots.get(i).getNumSamps() + "_numFam_"
                               + hPlots.get(i).getNumFam();
                  adds[i] = add;
                }
                String[] heritColumns = new String[currentClass.size()];
                String[] errorColumns = new String[currentClass.size()];
                String[][] columns = Files.paste(ArrayUtils.toStringArray(currentClass), comboHerit,
                                                 new int[] {0, 1, 2, 3, 4, 5, 6}, 0, adds, null,
                                                 log);
                for (int j = 0; j < columns.length; j++) {
                  heritColumns[j] = columns[j][2];
                  errorColumns[j] = columns[j][4];

                }
                ErrorBars errorBars = new ErrorBars(heritColumns, errorColumns);
                RScatter classHerit = new RScatter(comboHerit, comboHerit + ".rscript",
                                                   ext.removeDirectoryInfo(comboHerit),
                                                   comboHerit + ".pdf", columns[0][0], heritColumns,
                                                   SCATTER_TYPE.POINT, log);
                classHerit.setTitle(all.getrScatter().getTitle() + " \nHeritability by "
                                    + stratCats[classifyIndex]
                                    + (subsetDataHeritability == null ? ""
                                                                      : " and "
                                                                        + subsetDataHeritability));
                classHerit.setxLabel(all.getrScatter().getxLabel());
                classHerit.setyLabel(all.getrScatter().getyLabel());

                classHerit.setErrorBars(errorBars);
                classHerit.setFontsize(12);
                String[] availY = classHerit.getrSafeYColumns();
                String[] alt = new String[availY.length];
                for (int i = 0; i < alt.length; i++) {
                  alt[i] = availY[i].replaceAll("SOLAR_PROPORTION_HERITABLITY_", "");
                }
                classHerit.getErrorBars().setrSafeDataColumns(alt);
                classHerit.setrSafeAltYColumnNames(alt);
                classHerit.setOverWriteExisting(true);
                classHerit.setgTexts(otherDataGeomTexts.toArray(new GeomText[otherDataGeomTexts.size()]));
                classHerit.execute();

                if (subsetDataHeritability != null) {
                  // System.out.println(comboHerit);
                }
                scatters.add(classHerit);
              } else {
                log.reportTimeWarning("Skipping " + stratCats[classifyIndex]
                                      + " due to size of categories");
              }
            }
          }
        }
      };
      index++;
      return callable;
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String proj = null;
    String markers = null;
    String defaultDir = null;
    LS_TYPE lstype = LS_TYPE.SVD;
    int numThreads = 3;
    String samplesToBuildModels = null;
    String pcFile = null;
    String pedFile = null;
    double pcPercent = 0.05;
    boolean recomputeLRR = true;

    String usage = "\n" + "cnv.analysis.pca.CorrectionEvaluator requires 0-1 arguments\n";
    usage += "   (1) project filename (i.e. proj=" + proj + " (default))\n" + "";
    usage += "   (2) markers to Evaluate (i.e. markers=" + markers + " (default))\n" + "";
    usage += PSF.Ext.getOutputDirCommand(3, defaultDir);
    usage += PSF.Ext.getNumThreadsCommand(4, numThreads);
    usage += "   (5) type of linear regression (i.e. lstype=" + lstype + " ( default, options are "
             + java.util.Arrays.asList(LS_TYPE.values()) + "\n" + "";
    usage += "   (6) samples to generate models (i.e.samples= (no default))\n" + "";
    usage += "   (7) ped file to generate heritability (i.e.ped= (no default))\n" + "";
    usage += "   (8) alternate pc file to use (i.e.pcFile= (no default))\n" + "";
    usage += "   (9) percent of pcs to use base of the samples to build the models(i.e.pcPercent="
             + pcPercent + "(no default))\n" + "";
    usage += "   (10) skip on the fly centroids, only use gc-correction (i.e.recomputeLRR="
             + recomputeLRR + "(no default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        proj = arg.split("=")[1];
        numArgs--;
      }

      else if (arg.startsWith("markers=")) {
        markers = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("recomputeLRR=")) {
        recomputeLRR = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("samples=")) {
        samplesToBuildModels = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("ped=")) {
        pedFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("lstype=")) {
        lstype = LS_TYPE.valueOf(arg);
        numArgs--;
      } else if (arg.startsWith("pcFile=")) {
        pcFile = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("pcPercent=")) {
        pcPercent = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith(PSF.Ext.OUTPUT_DIR_COMMAND)) {
        defaultDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
        numThreads = ext.parseIntArg(arg);
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      Project project = new Project(proj);

      runAll(project, markers, samplesToBuildModels, defaultDir, pcFile, pedFile, lstype,
             recomputeLRR, pcPercent, true, numThreads);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
// if (Double.isNaN(lrrSdCut)) {
// project.getLog().reportTimeInfo("updating array specific lrr-sd");
// if (project.getArrayType() == ARRAY.ILLUMINA) {
// lrrSdCut = 0.30;
// } else {
// lrrSdCut = 0.35;
// }
// }
// else if (args[i].startsWith("lrrSdCut=")) {
// lrrSdCut = ext.parseDoubleArg(args[i]);
// numArgs--;
// } else if (args[i].startsWith("callRateCut=")) {
// callRateCut = ext.parseDoubleArg(args[i]);
// numArgs--;
// }

// if (!Files.exists(proj.SAMPLE_QC_FILENAME.getValue())) {
// String markerQCFile = outputDir + "markersForSampleQC.txt";
// ArrayList<String> toEvaluate = new ArrayList<String>();
// String[] names = proj.getMarkerNames();
// for (int i = 0; i < names.length; i++) {
// if (!proj.getArrayType().isCNOnly(names[i])) {
// toEvaluate.add(names[i]);
// }
// }
// Files.writeList(Array.toStringArray(toEvaluate), markerQCFile);
// LrrSd.init(proj, null, markerQCFile, markerQCFile, numthreads, null, false);
// }

// RScatter rScatterEvery = new RScatter(outputBox, ext.addToRoot(outputBox, ".rscript"),
// ext.removeDirectoryInfo(outputBox), dir + ext.removeDirectoryInfo(outputBox) + ".pdf",
// sampleDataStratCats[i], Array.toStringArray(pcYs), SCATTER_TYPE.BOX, log);
// rScatterEvery.setOverWriteExisting(true);
//
// rScatterEvery.setWidth(100);
// rScatterEvery.setHeight(60);
// rScatterEvery.setxLabel("PC (" + oType + " - sorted)");
// rScatterEvery.setTitle(iType + " " + bType);
// rScatterEvery.execute();
// rScatters.add(rScatterEvery);

// IterationResult iterationResult = cIterators[i].getIterationResult();
// RScatter rScatter = iterationResult.plotSummary(new String[] { "Rsquare_correction",
// "ICC_EVAL_CLASS_DUPLICATE_ALL", "ICC_EVAL_CLASS_DUPLICATE_SAME_VISIT", "ICC_EVAL_CLASS_FC",
// "SPEARMAN_CORREL_AGE", "SPEARMAN_CORREL_EVAL_DATA_SEX",
// "SPEARMAN_CORREL_EVAL_DATA_resid.mtDNaN.qPCR.MT001",
// "SPEARMAN_CORREL_EVAL_DATA_resid.mtDNA.qPCR",
// "SPEARMAN_CORREL_EVAL_DATA_Mt_DNA_relative_copy_number" }, proj.getLog());
// rScatters.add(rScatter);
// rScatters.add(iterationResult.plotRank(proj.getLog()));
// if (pedFile != null) {
// proj.getLog().reportTimeInfo("Since a ped file " + pedFile + " was provided, we will generate
// heritability estimates");
// EvaluationResult.prepareHeritability(proj, pedFile, iterationResult.getOutputSer());
// }

// System.exit(1);

// writer.println("PC\tGROUP\tVALUE");
// for (int PC = 0; PC < evaluationResults.length; PC++) {
// for (int j = 0; j < evaluationResults[PC].getEstimateData().length; j++) {
// for (int j2 = 0; j2 < booleanClassifier.getClassified().length; j2++) {
// if (booleanClassifier.getClassified()[j2][j]) {
// writer.println(PC + "\t" + booleanClassifier.getTitles()[j2] + "\t" +
// evaluationResults[PC].getEstimateData()[j]);
// }
// }
// }
// }
// writer.close();
// TODO alt column names, ICC center for affy

//
// RScatter rScatterSummaryMax = iterationResult.plotSummary(plotTitlesForSummary[i], i,
// proj.getLog());
// rScatterSummaryMax.setDirectLableGtexts(true);
// rScatterSummaryMax.setOnlyMaxMin(true);
// rScatterSummaryMax.setPlotVar(rScatterSummary.getPlotVar() + "max");
// rScatterSummaryMax.setFontsize(font);
// rScatterSummaryMax.getSeriesLabeler().setLabelMin(false);
//
// RScatter rScatterSummaryMin = iterationResult.plotSummary(plotTitlesForSummary[i], i,
// proj.getLog());
// rScatterSummaryMin.setDirectLableGtexts(true);
// rScatterSummaryMin.setOnlyMaxMin(true);
// rScatterSummaryMin.setPlotVar(rScatterSummary.getPlotVar() + "min");
// rScatterSummaryMin.setFontsize(font);
// rScatterSummaryMin.getSeriesLabeler().setLabelMax(false);

// private static double[][] loadIndeps(CorrectionEvaluator cEvaluator, String[] indepHeaders,
// double[][] indepMasks, String[] catHeaders, String[] catHeaderMask, Logger log) {
// ExtProjectDataParser parser = cEvaluator.getParser();
//
// double[][] extraIndeps = new
// double[Array.booleanArraySum(parser.getDataPresent())][indepHeaders.length];
//
// for (int i = 0; i < indepHeaders.length; i++) {
// int curSum = 0;
// if (parser.getNumericDataForTitle(indepHeaders[i]).length <= 0) {
// log.reportTimeError("Did not find " + indepHeaders[i] + " to load for independent variable");
// return null;
// } else {
// double[] data = parser.getNumericDataForTitle(indepHeaders[i]);
// for (int j = 0; j < data.length; j++) {
// boolean mask = false;
// for (int k = 0; k < indepMasks[i].length; k++) {
// if (data[j] == indepMasks[i][k]) {
// mask = true;
// }
// }
// if (mask) {
// extraIndeps[j][i] = Double.NaN;
// } else {
// extraIndeps[j][i] = data[j];
// curSum++;
// }
// }
// }
// log.reportTimeInfo(curSum + " samples for independant variable " + indepHeaders[i]);
// }
// for (int i = 0; i < catHeaders.length; i++) {
//
// if (parser.getStringDataForTitle(catHeaders[i]).length <= 0) {
// log.reportTimeError("Did not find " + catHeaders[i] + " to load for independent variable");
// return null;
// } else {
// String[] data = parser.getStringDataForTitle(catHeaders[i]);
// CategoricalPredictor predictor = new CategoricalPredictor(data, catHeaderMask, log);
// DummyCoding dummyCoding = predictor.createDummyCoding(false);
// for (int j = 0; j < dummyCoding.getDummyBoolean().length; j++) {
// double[] dummyData = Array.doubleArray(dummyCoding.getTitles().length, Double.NaN);
// if (dummyCoding.getDummyBoolean()[j]) {
// for (int k = 0; k < dummyData.length; k++) {
// dummyData[k] = dummyCoding.getDummyData()[k][j];
// }
// }
// extraIndeps[j] = Array.concatDubs(extraIndeps[j], dummyData);
// }
// }
// }
//
// return extraIndeps;
// }
// private void plotQuants(EvaluationResult[] evaluationResults, String outputQuant, double[] data,
// int numq, String sampleDataStratCats, Logger log) throws IOException {
// String realOut = ext.addToRoot(outputQuant, ".numQ_" + numq);
// PrintWriter writerQuant = Files.openAppropriateWriter(realOut);
// ArrayList<String> pcYsubQ = new ArrayList<String>();
// String dataTitle = sampleDataStratCats + ".quant_" + numq;
// pcYsubQ.add(dataTitle);
// writerQuant.print(dataTitle);
// System.out.println(realOut);
//
// for (int PC = 0; PC < evaluationResults.length; PC++) {
//
// String title = "PC" + PC;
// String part = title + "quant.matched.dist_" + numq;
// String full = title + "quant.full.dist_" + numq;
// writerQuant.print("\t" + part + "\t" + full);
// pcYsubQ.add(full);
// pcYsubQ.add(part);
// }
// // Quantiles[] quantiles = Quantiles.qetQuantilesFor(numQ, variableDomMatrix, titles,
// proj.getLog());
//
// double[][] toQuant = new double[evaluationResults.length * 2 + 1][data.length];
// toQuant[0] = data;
// int index = 1;
// for (int pcIndex = 0; pcIndex < evaluationResults.length; pcIndex++) {
// double[] estimateMatched = new double[getBasicPrep().getSamplesToEvaluate().length];
// double[] estimateFull = new double[getBasicPrep().getSamplesToEvaluate().length];
// for (int k = 0; k < getBasicPrep().getSamplesToEvaluate().length; k++) {
// estimateMatched[k] = Double.NaN;
// estimateFull[k] = Double.NaN;
// if (getBasicPrep().getSamplesToEvaluate()[k]) {
// estimateFull[k] = evaluationResults[pcIndex].getEstimateData()[k];
// if (!Double.isNaN(data[k])) {
// estimateMatched[k] = evaluationResults[pcIndex].getEstimateData()[k];
// }
// }
// }
// toQuant[index] = estimateMatched;
// index++;
// toQuant[index] = estimateFull;
// index++;
// }
// Quantiles[] quantiles = null;
// try {
// quantiles = Quantiles.qetQuantilesFor(numq, toQuant, Array.toStringArray(pcYsubQ), log);
// } catch (Exception e) {
// e.printStackTrace();
// System.exit(1);
// }
// // System.out.println(quantiles.length);
//
// for (int sampleQuantileIndex = 0; sampleQuantileIndex <
// quantiles[0].getQuantileMembership().length; sampleQuantileIndex++) {
// if (quantiles[0].getQuantileMembership()[sampleQuantileIndex] > 0) {
// for (int i = 0; i < quantiles.length; i++) {
//
// if (i > 0) {
// writerQuant.print("\t");
// }
// writerQuant.print(quantiles[i].getQuantileMembership()[sampleQuantileIndex]);
// }
// writerQuant.println();
// }
// }
//
// writerQuant.close();
// // System.out.println("DSF");
// // System.exit(1);
// }

// case QC_ASSOCIATION:
//
// SelectionResult result = PCSelector.select(proj, 0.05, STAT_TYPE.SPEARMAN_CORREL,
// SELECTION_TYPE.EFFECTIVE_M_CORRECTED);
// order = result.getOrder();
// if (result == null || order.length < 1) {
// log.reportTimeError("Could not select PCs from QC metrics, trying again");
// Files.copyFile(proj.SAMPLE_QC_FILENAME.getValue(), proj.SAMPLE_QC_FILENAME.getValue() +
// ext.getTimestampForFilename());
// new File(proj.SAMPLE_QC_FILENAME.getValue()).delete();
// LrrSd.init(proj, null, null, numthreads);
// result = PCSelector.select(proj, 0.05, STAT_TYPE.SPEARMAN_CORREL,
// SELECTION_TYPE.EFFECTIVE_M_CORRECTED);
// order = result.getOrder();
// if (order.length < 1) {
// return null;
// }
// }
// case WITH_INDEPS:
// extraIndeps = loadIndeps(cEvaluator, CorrectionEvaluator.INDEPS, new double[][] { { 0, 3, 4, 5,
// 6, 7, 8, 9, 10, Double.NaN }, { -1, Double.NaN } }, CorrectionEvaluator.INDEPS_CATS, new String[]
// { "NaN" }, log);
//
// if (extraIndeps == null) {
// log.reportTimeError("type = " + iType + " and were missing some of the following " +
// Array.toStr(CorrectionEvaluator.INDEPS));
// log.reportTimeError("Available = " + Array.toStr(cEvaluator.getParser().getNumericDataTitles()));
// valid = false;
// } else {
// boolean[] tmpInclude = new boolean[samplesForModels.length];
// Arrays.fill(tmpInclude, false);
// for (int i = 0; i < extraIndeps.length; i++) {
//
// if (samplesForModels[i]) {
// boolean hasNan = false;
// for (int j = 0; j < extraIndeps[i].length; j++) {
// if (Double.isNaN(extraIndeps[i][j])) {
// hasNan = true;
// }
// }
// if (!hasNan) {
// tmpInclude[i] = true;
// }
// }
// }
// log.reportTimeInfo("Original number of samples: " + Array.booleanArraySum(samplesForModels));
// log.reportTimeInfo("Number of samples with valid independant variables, final evaluation set: " +
// Array.booleanArraySum(tmpInclude));
// samplesForModels = tmpInclude;
// }
// break;
