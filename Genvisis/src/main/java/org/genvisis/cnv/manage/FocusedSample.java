package org.genvisis.cnv.manage;

import java.util.Hashtable;
import java.util.concurrent.Callable;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;

/**
 * Subsets samples
 *
 */
public class FocusedSample {

  private static class WorkerSubset implements Callable<FocusedSample> {
    private final String sample;
    private final Project original;
    private final Project newFocus;
    private final int[] focusedIndices;
    private final long newFingerPrint;
    private final boolean overwriteExisting;

    public WorkerSubset(String sample, Project original, Project newFocus, int[] focusedIndices,
                        long newFingerPrint, boolean overwriteExisting) {
      super();
      this.sample = sample;
      this.original = original;
      this.newFocus = newFocus;
      this.focusedIndices = focusedIndices;
      this.newFingerPrint = newFingerPrint;
      this.overwriteExisting = overwriteExisting;
    }

    @Override
    public FocusedSample call() throws Exception {
      Sample samp = original.getFullSampleFromRandomAccessFile(sample);
      FocusedSample focusedSample = null;
      if (samp != null) {
        String newSampleFileName =
            newFocus.SAMPLE_DIRECTORY.getValue(true, true) + sample + Sample.SAMPLE_FILE_EXTENSION;
        focusedSample = new FocusedSample(focusedIndices, samp, newFingerPrint, overwriteExisting);
        focusedSample.saveToRandomAccessFile(newSampleFileName, samp.getSampleName(),
                                             original.getLog());
      } else {
        original.getLog().reportTimeError("Could not load sample " + sample);
      }
      // TODO Auto-generated method stub
      return focusedSample;
    }

  }

  public static boolean focusAProject(Project original, Project newFocus, String[] markersToUse,
                                      boolean[] samplesToUse, long newFingerPrint, int numThreads,
                                      boolean overwriteExisting, Logger log) {
    boolean focused = true;
    if (original.PROJECT_DIRECTORY.getValue().equals(newFocus.PROJECT_DIRECTORY.getValue())) {
      log.reportTimeError("The focused project must have a different project directory than the original, halting");
      focused = false;
    } else if (original.SAMPLE_DIRECTORY.getValue(false, true)
                                        .equals(newFocus.SAMPLE_DIRECTORY.getValue(false, true))) {
      log.reportTimeError("The focused project must have a different sample directory than the original, halting");
      focused = false;
    } else if (markersToUse == null || samplesToUse == null) {
      log.reportTimeError("Please provide subsets... if you want they can be all markers or samples");
      focused = false;
    } else {
      log.reportTimeInfo("Markers to export = " + markersToUse.length);
      log.reportTimeInfo("Samples to export = " + Array.booleanArraySum(samplesToUse));
      int[] markerToUseIndices =
          ext.indexLargeFactors(markersToUse, original.getMarkerNames(), true, log, true, false);
      WorkerHive<FocusedSample> hive = new WorkerHive<FocusedSample>(numThreads, 10, log);
      WorkerSubset[] workerSubsets = getWorkers(original, newFocus, markerToUseIndices,
                                                samplesToUse, newFingerPrint, overwriteExisting);
      hive.addCallables(workerSubsets);
      hive.setReportEvery(100);
      hive.execute(true);
      FocusedSample[] focusedSamples =
          hive.getResults().toArray(new FocusedSample[hive.getResults().size()]);
      writeOutliers(newFocus, focusedSamples);
      return focused;
    }
    return focused;
  }

  private static WorkerSubset[] getWorkers(Project original, Project newFocus, int[] focusedIndices,
                                           boolean[] samplesToUse, long newFingerPrint,
                                           boolean overwriteExisting) {
    String[] samples = samplesToUse == null ? original.getSamples()
                                            : Array.subArray(original.getSamples(), samplesToUse);
    WorkerSubset[] workerSubsets = new WorkerSubset[samples.length];
    for (int i = 0; i < samples.length; i++) {
      workerSubsets[i] = new WorkerSubset(samples[i], original, newFocus, focusedIndices,
                                          newFingerPrint, overwriteExisting);
    }
    return workerSubsets;
  }

  // private boolean fail;

  public static void main(String[] args) {
    test();
  }

  public static void test() {
    Project proj = new Project(null, false);
    Project proj2 = new Project(null, false);
    proj2.setProperty(proj.PROJECT_DIRECTORY, proj.PROJECT_DIRECTORY.getValue() + "subSample/");
    String[] markers = proj.getTargetMarkers();
    focusAProject(proj, proj2, markers, null, proj.getMarkerSet().getFingerprint(), 8, true,
                  proj.getLog());
    proj2.getMarkerLookup();
  }

  private static void writeOutliers(Project newFocus, FocusedSample[] focusedSamples) {
    Hashtable<String, Float> outliers = new Hashtable<String, Float>();
    for (FocusedSample focusedSample2 : focusedSamples) {
      outliers.putAll(focusedSample2.getOutliers());
    }
    String outlierFileName = newFocus.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser";
    newFocus.getLog().reportTimeInfo("Writing outliers to " + outlierFileName);
    SerializedFiles.writeSerial(outliers, outlierFileName);

  }

  private final Sample focusedSample;

  private final boolean overwriteExisting;

  private final Hashtable<String, Float> outliers;

  public FocusedSample(int[] focusedIndices, Sample sample, long newFingerPrint,
                       boolean overwriteExisting) {
    super();
    this.overwriteExisting = overwriteExisting;
    focusedSample = getFocusedSample(sample, focusedIndices, newFingerPrint);
    outliers = new Hashtable<String, Float>();
    // this.fail =false;
  }

  private Sample getFocusedSample(Sample sample, int[] focusedIndices, long newFingerPrint) {
    byte[] abGenotypes =
        sample.getAB_Genotypes() == null ? null
                                         : Array.subArray(sample.getAB_Genotypes(), focusedIndices);
    byte[] forwardGenotypes = sample.getForwardGenotypes() == null ? null
                                                                   : Array.subArray(sample.getForwardGenotypes(),
                                                                                    focusedIndices);

    float[] xs = sample.getXs() == null ? null : Array.subArray(sample.getXs(), focusedIndices);
    float[] ys = sample.getYs() == null ? null : Array.subArray(sample.getYs(), focusedIndices);

    float[] lrrs =
        sample.getLRRs() == null ? null : Array.subArray(sample.getLRRs(), focusedIndices);
    float[] bafs =
        sample.getBAFs() == null ? null : Array.subArray(sample.getBAFs(), focusedIndices);

    float[] gcs = sample.getGCs() == null ? null : Array.subArray(sample.getGCs(), focusedIndices);
    boolean canXYBeNegative = sample.getCanXYBeNegative();
    Sample focusedSample = new Sample(sample.getSampleName(), newFingerPrint, gcs, xs, ys, bafs,
                                      lrrs, forwardGenotypes, abGenotypes, canXYBeNegative);
    return focusedSample;
  }

  public Hashtable<String, Float> getOutliers() {
    return outliers;
  }

  public boolean saveToRandomAccessFile(String filename, String sampleName, Logger log) {
    boolean saved = true;
    if (!overwriteExisting && Files.exists(filename)) {
      log.reportTimeWarning(filename + " already exists, will not overwrite");
      saved = false;
    } else {
      focusedSample.saveToRandomAccessFile(filename, outliers, focusedSample.getSampleName());
    }
    return saved;
  }

}
