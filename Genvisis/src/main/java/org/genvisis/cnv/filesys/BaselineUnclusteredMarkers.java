package org.genvisis.cnv.filesys;

import java.io.Serializable;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Set;

import org.pankratzlab.common.Files;
import org.pankratzlab.common.ProgressMonitor;
import org.pankratzlab.common.SerializedFiles;

public class BaselineUnclusteredMarkers implements Serializable {

  private static final long serialVersionUID = 1L;
  private static final String BASELINE_UNCLUSTERED_MARKERS_FILENAME = "baselineUnclusteredMarkers";

  public static final float MARKER_CALLRATE_THRESHOLD = 0.5f;

  private final Hashtable<String, Float> markerCallrate;

  private BaselineUnclusteredMarkers() {
    markerCallrate = new Hashtable<>();
  }

  private void serialize(Project proj) {
    SerializedFiles.writeSerial(this, getBaselineUnclusteredMarkersFile(proj));
  }

  public boolean markerUnclustered(String marker) {
    return markerCallrate.get(marker) != null;
  }

  public float markerCallrate(String marker) {
    Float callrate = markerCallrate.get(marker);
    return callrate == null ? -1f : callrate;
  }

  /**
   * @return a {@link Set} of the names of all markers for which {@link #markerUnclustered(String)}
   *         would return true
   */
  public Set<String> unclusteredMarkers() {
    return Collections.unmodifiableSet(markerCallrate.keySet());
  }

  private static String getBaselineUnclusteredMarkersFile(Project proj) {
    return proj.DATA_DIRECTORY.getValue() + BASELINE_UNCLUSTERED_MARKERS_FILENAME
           + Files.SERIALIZED_FILE_EXTENSION;
  }

  public static boolean createBaselineUnclusteredMarkersFileFromSamples(Project proj) {
    boolean success = true;
    String taskName = "createBaselineUnclusteredMarkersFileFromSamples_sampleCalls";
    proj.getProgressMonitor().beginDeterminateTask(taskName, "Summing marker calls for each sample",
                                                   proj.getNumberOfParsedSamples(),
                                                   ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
    int[] autosomalMarkerIndices = proj.getAutosomalMarkerIndices();
    int[] autosomalMarkerCalls = new int[autosomalMarkerIndices.length];
    for (String sampleID : proj.getSamples()) {
      proj.getProgressMonitor().updateTask(taskName);
      Sample sample = proj.getPartialSampleFromRandomAccessFile(sampleID, false, false, false,
                                                                false, true);
      if (sample.getAB_Genotypes() != null && sample.getAB_Genotypes().length != 0) {
        for (int i = 0; i < autosomalMarkerIndices.length; i++) {
          if (sample.getAB_Genotypes()[autosomalMarkerIndices[i]] != -1) {
            autosomalMarkerCalls[i]++;
          }
        }
      } else if (sample.getForwardGenotypes() == null || sample.getForwardGenotypes().length == 0) {
        for (int i = 0; i < autosomalMarkerIndices.length; i++) {
          if (sample.getForwardGenotypes()[autosomalMarkerIndices[i]] != 0) {
            autosomalMarkerCalls[i]++;
          }
        }
      } else {
        proj.getLog()
            .reportError("Could not generate Baseline Unclustered Markers File, no genotypes for samples");
        success = false;
        break;
      }
    }
    proj.getProgressMonitor().endTask(taskName);
    if (success) {
      String[] autosomalMarkers = proj.getAutosomalMarkers();
      taskName = "createBaselineUnclusteredMarkersFileFromSamples_callrates";
      proj.getProgressMonitor().beginDeterminateTask(taskName, "Checking marker callrates",
                                                     autosomalMarkers.length,
                                                     ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
      BaselineUnclusteredMarkers bum = new BaselineUnclusteredMarkers();
      for (int i = 0; i < autosomalMarkers.length; i++) {
        proj.getProgressMonitor().updateTask(taskName);
        float callrate = (float) autosomalMarkerCalls[i] / proj.getNumberOfParsedSamples();
        if (callrate < MARKER_CALLRATE_THRESHOLD) {
          bum.markerCallrate.put(autosomalMarkers[i], callrate);
        }
      }
      proj.getProgressMonitor().endTask(taskName);
      bum.serialize(proj);
    }
    return success;
  }

  public static boolean baselineUnclusteredMarkersFileExists(Project proj) {
    return Files.exists(getBaselineUnclusteredMarkersFile(proj));
  }

  public static BaselineUnclusteredMarkers getProjBaselineUnclusteredMarkers(Project proj) {
    if (baselineUnclusteredMarkersFileExists(proj)) {
      return (BaselineUnclusteredMarkers) SerializedFiles.readSerial(getBaselineUnclusteredMarkersFile(proj));
    }
    return null;
  }

}
