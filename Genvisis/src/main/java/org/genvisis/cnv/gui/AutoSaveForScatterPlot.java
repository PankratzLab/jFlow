package org.genvisis.cnv.gui;

import java.io.File;

import org.genvisis.cnv.filesys.AnnotationCollection;
import org.genvisis.cnv.filesys.ClusterFilterCollection;

public class AutoSaveForScatterPlot implements Runnable {

  private ClusterFilterCollection clusterFilters;
  private String clusterFilterFilename;
  private AnnotationCollection annotations;
  private String annotationFilename;
  private final int period;
  private boolean killed;
  private boolean isClusterFiltersUpdated;
  private boolean isAnnotationsUpdated;

  public AutoSaveForScatterPlot(ClusterFilterCollection clusterFilterCollection,
                                String clusterFilterFilename,
                                AnnotationCollection AnnotationCollection,
                                String annotationFilename, int periodInSeconds) {
    clusterFilters = clusterFilterCollection;
    this.clusterFilterFilename = clusterFilterFilename;
    annotations = AnnotationCollection;
    this.annotationFilename = annotationFilename;
    period = periodInSeconds;
    killed = false;
    isClusterFiltersUpdated = false;
    isAnnotationsUpdated = false;
  }

  public void addToAutoSave(AnnotationCollection collection, String annotationFilename) {
    annotations = collection;
    this.annotationFilename = annotationFilename;
  }

  public void addToAutoSave(ClusterFilterCollection collection, String clusterFilterFilename) {
    clusterFilters = collection;
    this.clusterFilterFilename = clusterFilterFilename;
  }

  public boolean isAnnotationNull() {
    return annotations == null;
  }

  public boolean isClusterFilterNull() {
    return clusterFilters == null;
  }

  public void kill() {
    if (clusterFilters != null) {
      new File(clusterFilterFilename).delete();
    }
    if (annotations != null) {
      new File(annotationFilename).delete();
    }

    killed = true;
  }

  @Override
  public void run() {
    while (!killed) {
      saveNow();
      try {
        Thread.sleep(period * 1000);
      } catch (InterruptedException ie) {
      }
    }
  }

  public void saveNow() {
    if (isClusterFiltersUpdated) {
      clusterFilters.serialize(clusterFilterFilename);
      isClusterFiltersUpdated = false;
    }

    if (isAnnotationsUpdated) {
      annotations.serialize(annotationFilename);
      isAnnotationsUpdated = false;
    }
  }

  public void setAnnotationUpdated(boolean isUpdated) {
    isAnnotationsUpdated = isUpdated;
  }

  public void setClusterFilterUpdated(boolean isUpdated) {
    isClusterFiltersUpdated = isUpdated;
  }
}
