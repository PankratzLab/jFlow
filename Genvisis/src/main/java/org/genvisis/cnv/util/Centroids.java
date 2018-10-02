package org.genvisis.cnv.util;

import org.genvisis.cnv.analysis.CentroidCompute;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity;
import org.genvisis.cnv.filesys.ClusterFilterCollection;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.common.Logger;

public class Centroids {

  public static CentroidCompute prepareProperCentroid(ARRAY array, MarkerData markerData,
                                                       int[] sampleSex,
                                                       boolean[] samplesToUseCluster,
                                                       double missingnessThreshold,
                                                       double confThreshold,
                                                       ClusterFilterCollection clusterFilterCollection,
                                                       boolean medianCenter, Logger log) {
    CentroidCompute centroid;
    if (PrincipalComponentsIntensity.isAffyIntensityOnly(array, markerData)) {
      // centroid = markerData.getCentroid(sampleSex, samplesToUseCluster, true,
      // missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, log);
      centroid = markerData.getCentroid(sampleSex, samplesToUseCluster, false, missingnessThreshold,
                                        confThreshold, clusterFilterCollection, medianCenter, log);
      CentroidCompute.setFakeAB(markerData, centroid, clusterFilterCollection, .1f, log);
    } else {
      centroid = markerData.getCentroid(sampleSex, samplesToUseCluster, false, missingnessThreshold,
                                        confThreshold, clusterFilterCollection, medianCenter, log);
      if (array.isCNOnly(markerData.getMarkerName())) {
        CentroidCompute.setFakeAB(markerData, centroid, clusterFilterCollection, 0, log);
      }
    }
    centroid.computeCentroid();
    return centroid;
  }

}
