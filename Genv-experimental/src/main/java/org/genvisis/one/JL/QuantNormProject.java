package org.genvisis.one.JL;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.concurrent.Callable;

import org.genvisis.cnv.analysis.CentroidCompute;
import org.genvisis.cnv.filesys.Centroids;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.cnv.qc.LrrSd;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.SerializedFiles;
import org.pankratzlab.common.WorkerTrain;
import org.pankratzlab.common.WorkerTrain.AbstractProducer;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.stats.Histogram.DynamicHistogram;
import org.pankratzlab.common.stats.QuantileNormalization;
import org.pankratzlab.common.stats.Rscript.COLUMNS_MULTIPLOT;
import org.pankratzlab.common.stats.Rscript.PLOT_DEVICE;
import org.pankratzlab.common.stats.Rscript.RScatter;
import org.pankratzlab.common.stats.Rscript.RScatters;
import org.pankratzlab.common.stats.Rscript.SCATTER_TYPE;
import org.pankratzlab.common.stats.Rscript.VertLine;

public class QuantNormProject {

  private static final double THRESHOLD_FACTOR = 1.5;

  private static void quantNormProject(Project proj, int numthreads) {
    numthreads = 7;
    Project projNorm = Project.prepareNewProject(proj, "quantNorm");
    QNormProducer producer = new QNormProducer(proj, projNorm, proj.getSamples());
    try (WorkerTrain<Hashtable<String, Float>> train = new WorkerTrain<>(producer, numthreads, 2,
                                                                         proj.getLog())) {
      Hashtable<String, Float> outliers = new Hashtable<>();
      int index = 0;
      long time = System.currentTimeMillis();

      while (train.hasNext()) {

        if (index % numthreads == 0) {
          proj.getLog()
              .reportTimeInfo("Normalized " + index + " samples in " + ext.getTimeElapsed(time));
          time = System.currentTimeMillis();

        }
        Hashtable<String, Float> tmp = train.next();
        outliers.putAll(tmp);
        proj.getLog().reportTimeInfo(index + "");
        index++;
      }

      if (!Files.exists(projNorm.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser")) {
        SerializedFiles.writeSerial(outliers, projNorm.SAMPLE_DIRECTORY.getValue(true, true)
                                              + "outliers.ser");
        TransposeData.transposeData(projNorm, 2000000000, false);
        CentroidCompute.computeAndDumpCentroids(projNorm);
        Centroids.recompute(projNorm, projNorm.CUSTOM_CENTROIDS_FILENAME.getValue());
      }
    }

    summarizeMetrics(proj, projNorm, numthreads);

  }

  private static void summarizeMetrics(Project projOriginal, Project projCorrected,
                                       int numThreads) {
    String prior = ext.addToRoot(projCorrected.SAMPLE_QC_FILENAME.getValue(), ".priorToCorrection");
    projOriginal.SAMPLE_QC_FILENAME.setValue(prior);
    if (!Files.exists(projCorrected.SAMPLE_QC_FILENAME.getValue())) {
      LrrSd.init(projCorrected, null, null, null, null, numThreads);
    }
    if (!Files.exists(projOriginal.SAMPLE_QC_FILENAME.getValue())) {
      LrrSd.init(projOriginal, null, null, null, null, numThreads);
    }
    String priorNewCent = ext.addToRoot(projCorrected.SAMPLE_QC_FILENAME.getValue(), ".newCents");
    if (!Files.exists(priorNewCent)) {
      projOriginal.SAMPLE_QC_FILENAME.setValue(priorNewCent);
      String cent = ext.addToRoot(projOriginal.CUSTOM_CENTROIDS_FILENAME.getValue(), "newCent");
      projOriginal.CUSTOM_CENTROIDS_FILENAME.setValue(cent);
      CentroidCompute.computeAndDumpCentroids(projOriginal);
      LrrSd.init(projOriginal, null, null, null, cent, numThreads);
    }

    String comboQC = ext.addToRoot(projCorrected.SAMPLE_QC_FILENAME.getValue(), ".combinedQC");
    String[] orginalFiles = new String[] {prior, priorNewCent,
                                          projCorrected.SAMPLE_QC_FILENAME.getValue()};
    String[] titles = new String[] {"ORIGINAL", "ORIGINAL_NEW_CENT", "QuantNorm"};
    String[] fullHeader = ArrayUtils.concatAll(new String[] {LrrSd.SAMPLE_COLUMN},
                                               LrrSd.NUMERIC_COLUMNS);
    int[] indices = ext.indexFactors(fullHeader,
                                     Files.getHeaderOfFile(orginalFiles[0], projOriginal.getLog()),
                                     true);
    String[][] newColums = Files.paste(orginalFiles, comboQC, indices, 0, titles,
                                       new String[] {LrrSd.SAMPLE_COLUMN}, projOriginal.getLog());
    ArrayList<Integer> lrrIndices = new ArrayList<>();
    for (int i = 0; i < LrrSd.NUMERIC_COLUMNS.length; i++) {
      if (LrrSd.NUMERIC_COLUMNS[i].startsWith("LRR_SD")) {
        lrrIndices.add(i);
      }
    }
    DynamicHistogram[][] histograms = new DynamicHistogram[orginalFiles.length][lrrIndices.size()];
    for (int i = 0; i < histograms.length; i++) {
      histograms[i] = DynamicHistogram.initHistograms(histograms[i].length, 0, 1, 2);
    }

    String status = "STATUS";
    String comboBox = comboQC + ".box";
    try {
      PrintWriter writer = Files.openAppropriateWriter(comboBox);
      writer.println(status + "\t" + ArrayUtils.toStr(fullHeader));
      for (int i = 0; i < orginalFiles.length; i++) {

        BufferedReader reader = Files.getAppropriateReader(orginalFiles[i]);
        reader.readLine();
        while (reader.ready()) {
          String[] line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
          String[] data = ArrayUtils.subArray(line, indices);
          if (!line[indices[0]].equals(LrrSd.SAMPLE_COLUMN)) {
            writer.println(titles[i] + "\t" + ArrayUtils.toStr(data));
            for (int j = 0; j < lrrIndices.size(); j++) {
              histograms[i][j].addDataPointToHistogram(Double.parseDouble(data[lrrIndices.get(j)
                                                                               + 1]));
            }
          }
        }
        reader.close();

      }
      writer.close();

    } catch (FileNotFoundException fnfe) {
      projOriginal.getLog().reportError("Error: file(s) \"" + ArrayUtils.toStr(orginalFiles)
                                        + "\" not found in current directory");
      return;
    } catch (IOException ioe) {
      projOriginal.getLog()
                  .reportError("Error reading file(s) \"" + ArrayUtils.toStr(orginalFiles) + "\"");
      return;
    } catch (Exception e) {
      projOriginal.getLog().reportError("Error writing to " + comboBox);
      projOriginal.getLog().reportException(e);
    }
    String gcLookDir = projCorrected.PROJECT_DIRECTORY.getValue() + "gc_analysis/";
    new File(gcLookDir).mkdirs();
    ArrayList<RScatter> rscatterHist = new ArrayList<>();
    ArrayList<DynamicHistogram> all = new ArrayList<>();
    ArrayList<String> allTitles = new ArrayList<>();

    for (int i = 0; i < histograms[0].length; i++) {
      ArrayList<DynamicHistogram> tmps = new ArrayList<>();
      String out = gcLookDir + "gc_" + LrrSd.NUMERIC_COLUMNS[lrrIndices.get(i)] + ".hist";
      String curQc = LrrSd.NUMERIC_COLUMNS[lrrIndices.get(i)];
      String[] titleTmp = ArrayUtils.tagOn(titles, curQc, null);
      for (String element : titleTmp) {
        allTitles.add(element);
      }
      for (DynamicHistogram[] histogram : histograms) {
        tmps.add(histogram[i]);
        all.add(histogram[i]);
      }
      DynamicHistogram.dumpToSameFile(tmps.toArray(new DynamicHistogram[tmps.size()]), titleTmp,
                                      out, false, projCorrected.getLog());
      RScatter rtmp = new RScatter(out, out + ".rscript", ext.removeDirectoryInfo(out),
                                   out + ".jpeg", "Bin", titleTmp, SCATTER_TYPE.POINT,
                                   projCorrected.getLog());
      rtmp.setyLabel("Count " + curQc);
      rtmp.setxLabel("Bin " + curQc);
      rtmp.setTitle("Original V News " + curQc);

      RScatter rtmpCumulative = new RScatter(out, out + "cu.rscript",
                                             ext.removeDirectoryInfo(out) + "cu", out + "cu.jpeg",
                                             "Bin", ArrayUtils.tagOn(titleTmp, "Cumulative_", null),
                                             SCATTER_TYPE.POINT, projCorrected.getLog());
      rtmpCumulative.setyLabel("Cumulative Count " + curQc);
      rtmpCumulative.setxLabel("Bin " + curQc);
      rtmpCumulative.setTitle("Original V News " + curQc);
      rscatterHist.add(rtmp);
      rscatterHist.add(rtmpCumulative);

    }
    String allHist = gcLookDir + "Gc_summaryHistAll.txt";
    DynamicHistogram.dumpToSameFile(all.toArray(new DynamicHistogram[all.size()]),
                                    ArrayUtils.toStringArray(allTitles), allHist, false,
                                    projCorrected.getLog());
    VertLine[] vertLines = new VertLine[] {new VertLine(0.30)};
    RScatter rtmpCumulativeAll = new RScatter(allHist, allHist + "cu.rscript",
                                              ext.removeDirectoryInfo(allHist) + "cu",
                                              allHist + "cu.jpeg", "Bin",
                                              ArrayUtils.tagOn(ArrayUtils.toStringArray(allTitles),
                                                               "Cumulative_", null),
                                              SCATTER_TYPE.POINT, projCorrected.getLog());
    rtmpCumulativeAll.setyLabel("Cumulative Count All");
    rtmpCumulativeAll.setxLabel("Bin ");
    rtmpCumulativeAll.setTitle("Original V News ");
    rtmpCumulativeAll.setVertLines(vertLines);
    rscatterHist.add(rtmpCumulativeAll);

    String finalHistRoot = gcLookDir + "Gc_summaryHist";

    RScatters rScattersHists = new RScatters(rscatterHist.toArray(new RScatter[rscatterHist.size()]),
                                             finalHistRoot + ".rscript", finalHistRoot + ".pdf",
                                             COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1, PLOT_DEVICE.PDF,
                                             projCorrected.getLog());
    rScattersHists.execute();

    ArrayList<RScatter> rScatters = new ArrayList<>();

    for (int i = 1; i < newColums[0].length; i++) {

      String curFile = gcLookDir + "gc_" + LrrSd.NUMERIC_COLUMNS[i - 1];
      RScatter rScatter = new RScatter(comboQC, curFile + ".rscript",
                                       ext.removeDirectoryInfo(curFile), curFile + ".pdf",
                                       newColums[0][i],
                                       new String[] {newColums[1][i], newColums[2][i]},
                                       SCATTER_TYPE.POINT, projCorrected.getLog());
      rScatter.setTitle("Original V News " + LrrSd.NUMERIC_COLUMNS[i - 1]);
      rScatter.setxLabel("Original " + LrrSd.NUMERIC_COLUMNS[i - 1]);
      rScatter.setyLabel("News " + LrrSd.NUMERIC_COLUMNS[i - 1]);
      rScatter.setRegLines(true);
      rScatters.add(rScatter);

      String curFileBox = gcLookDir + "gc_" + LrrSd.NUMERIC_COLUMNS[i - 1] + "_box";
      RScatter rScatterBox = new RScatter(comboBox, curFileBox + ".rscript",
                                          ext.removeDirectoryInfo(curFileBox), curFileBox + ".pdf",
                                          status, new String[] {LrrSd.NUMERIC_COLUMNS[i - 1]},
                                          SCATTER_TYPE.BOX, projCorrected.getLog());
      rScatterBox.setTitle("Original V News " + LrrSd.NUMERIC_COLUMNS[i - 1]);
      rScatterBox.setyLabel(LrrSd.NUMERIC_COLUMNS[i - 1]);
      rScatters.add(rScatterBox);
    }
    String finalRoot = gcLookDir + "Gc_summary";
    RScatters rScatters2 = new RScatters(rScatters.toArray(new RScatter[rScatters.size()]),
                                         finalRoot + ".rscript", finalRoot + ".pdf",
                                         COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1, PLOT_DEVICE.PDF,
                                         projCorrected.getLog());
    rScatters2.execute();
  }

  private static class QNormProducer extends AbstractProducer<Hashtable<String, Float>> {

    private final Project projOriginal;
    private final Project projNorm;
    private final String[] samples;
    private int index;

    public QNormProducer(Project projOriginal, Project projNorm, String[] samples) {
      super();
      this.projOriginal = projOriginal;
      this.projNorm = projNorm;
      projNorm.SAMPLE_DIRECTORY.getValue(true, false);
      this.samples = samples;
    }

    @Override
    public boolean hasNext() {
      return index < samples.length;
    }

    @Override
    public Callable<Hashtable<String, Float>> next() {
      final String currentSamp = samples[index];
      Callable<Hashtable<String, Float>> callable = new Callable<Hashtable<String, Float>>() {

        @Override
        public Hashtable<String, Float> call() throws Exception {
          // TODO Auto-generated method stub
          return quantNormSample(projOriginal, projNorm, currentSamp);
        }
      };
      index++;
      return callable;
    }
  }

  private static Hashtable<String, Float> quantNormSample(Project projOriginal, Project projNorm,
                                                          String sample) {
    Hashtable<String, Float> outliers = new Hashtable<>();
    String output = projNorm.SAMPLE_DIRECTORY.getValue() + sample + Sample.SAMPLE_FILE_EXTENSION;

    if (!Files.exists(output)) {
      Sample samp = projOriginal.getFullSampleFromRandomAccessFile(sample);
      double[][] xy = new double[][] {ArrayUtils.toDoubleArray(samp.getXs()),
                                      ArrayUtils.toDoubleArray(samp.getYs())};
      QuantileNormalization qn = new QuantileNormalization(xy, projNorm.getLog());
      qn.setThresholdFactor(THRESHOLD_FACTOR);
      projNorm.getLog().reportTimeInfo("Beginning normalizing " + sample);
      long time = System.currentTimeMillis();
      qn.normalize();
      projNorm.getLog()
              .reportTimeInfo("Finished normalizing " + sample + " in " + ext.getTimeElapsed(time));
      double[][] xyNorm = qn.getNormData();
      time = System.currentTimeMillis();
      Sample sampNorm = new Sample(sample, samp.getFingerprint(), samp.getGCs(),
                                   ArrayUtils.toFloatArray(xyNorm[0]),
                                   ArrayUtils.toFloatArray(xyNorm[1]), samp.getBAFs(),
                                   samp.getLRRs(), samp.getForwardGenotypes(),
                                   samp.getAB_Genotypes(), samp.getCanXYBeNegative());
      sampNorm.saveToRandomAccessFile(output, outliers, sample);
      projNorm.getLog()
              .reportTimeInfo("Saved sample " + sample + " in " + ext.getTimeElapsed(time));
    }

    return outliers;
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = null;

    String usage = "\n" + "one.JL.QuantNormProject requires 0-1 arguments\n"
                   + "   (1) filename (i.e. proj=" + filename + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
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
      Project proj = new Project(filename);
      System.out.println(filename);
      proj.setPropertyFilename(filename);
      proj.PROJECT_PROPERTIES_FILENAME.setValue(filename);
      quantNormProject(proj, proj.NUM_THREADS.getValue());
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
