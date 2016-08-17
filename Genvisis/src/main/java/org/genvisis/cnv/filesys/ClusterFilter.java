package org.genvisis.cnv.filesys;

import java.io.Serializable;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.IntVector;

import com.google.common.primitives.Ints;

/**
 * This is a data structure to hold a single filter that can be used to screen the data points.
 * 
 * @author npankrat and zxu
 *
 */
public class ClusterFilter implements Serializable {
  private static final long serialVersionUID = 1L;

  private byte plotType;
  private float rawXMin;
  private float rawYMin;
  private float rawXmax;
  private float rawYmax;
  private byte newGenotype;

  public ClusterFilter() {}

  public ClusterFilter(byte plotType, float rawXMin, float rawYMin, float rawXmax, float rawYmax) {
    this(plotType, rawXMin, rawYMin, rawXmax, rawYmax, (byte) 0);
  }

  public ClusterFilter(byte plotType, float rawXMin, float rawYMin, float rawXmax, float rawYmax,
      byte newGenotype) {
    this.plotType = plotType;
    this.rawXMin = rawXMin;
    this.rawYMin = rawYMin;
    this.rawXmax = rawXmax;
    this.rawYmax = rawYmax;
    this.newGenotype = newGenotype;
  }

  public ClusterFilter(byte plotType, float rawXMin, float rawYMin, float rawXmax, float rawYmax,
      MarkerData markerData) {
    this.plotType = plotType;
    this.rawXMin = rawXMin;
    this.rawYMin = rawYMin;
    this.rawXmax = rawXmax;
    this.rawYmax = rawYmax;
    // this.newGenotype = suggestedNewGenoTypeByNearbyCentroid(markerData);
    newGenotype = suggestedNewGenoTypeByNearbyPoint(markerData);
  }

  public byte getCluterGenotype() {
    return newGenotype;
  }

  public byte getPlotType() {
    return plotType;
  }

  public float getXMax() {
    return rawXmax;
  }

  public float getXMin() {
    return rawXMin;
  }

  public float getYMax() {
    return rawYmax;
  }

  public float getYMin() {
    return rawYMin;
  }

  public void setClusterGenotype(byte newGenotype) {
    this.newGenotype = newGenotype;
  }

  public void setPlotType(byte plotType) {
    this.plotType = plotType;
  }

  public byte suggestedNewGenoTypeByNearbyCentroid(MarkerData markerData) {
    float[] realX;
    float[] realY;
    byte result = -1;
    float xSum, ySum;
    Hashtable<String, IntVector> hash;
    String cluster;
    float distance;
    float distancetemp = 0;
    byte[] genotypes;
    IntVector iv;
    float[][] clusterCenters;
    int[] genotypeCount;
    byte oldGenotype;
    int genotypeFrequencyCount;
    int[] genotypeIndices;

    switch (getPlotType()) {
      case 0:
        // realX = markerData.getX_Raws();
        // realY = markerData.getY_Raws();
        // break;
        // case 1:
        realX = markerData.getXs();
        realY = markerData.getYs();
        break;
      case 1:
        realX = markerData.getThetas();
        realY = markerData.getRs();
        break;
      case 2:
        realX = markerData.getBAFs();
        realY = markerData.getLRRs();
        break;
      default:
        realX = markerData.getXs();
        realY = markerData.getYs();
    }

    hash = new Hashtable<String, IntVector>();
    genotypes = markerData.getAbGenotypes();
    // iterate through all samples
    for (int i = 0; i < genotypes.length; i++) {
      if (realX[i] >= rawXMin && realY[i] >= rawYMin && realX[i] <= rawXmax
          && realY[i] <= rawYmax) {
        cluster = "3"; // "3" identifies the data points within the cluster filter
      } else {
        cluster = genotypes[i] + "";
      }
      if (hash.containsKey(cluster)) {
        iv = hash.get(cluster);
      } else {
        hash.put(cluster, iv = new IntVector());
      }
      iv.add(i);
    }

    // Find the existing genotype of the data points within the cluster filter
    iv = hash.get("3");
    if (iv != null) {
      /*
       * Search for the genotype of majority of the data points If there is any data point with
       * missing value found, the majority's genotype is set to Missing Value.
       */
      genotypeCount = new int[] {0, 0, 0};
      oldGenotype = -2;
      for (int i = 0; i < 3; i++) {
        genotypeIndices = Ints.toArray(iv);
        for (int genotypeIndice : genotypeIndices) {
          if (genotypes[genotypeIndice] == i) {
            genotypeCount[i]++;
          } else if (genotypes[genotypeIndice] == -1) {
            oldGenotype = -1;
          }
        }
      }
      // Search for the maximum in int[] genotypeCount
      if (oldGenotype == -2) {
        genotypeFrequencyCount = -1;
        for (byte i = 0; i < 3; i++) {
          if (genotypeCount[i] > genotypeFrequencyCount) {
            genotypeFrequencyCount = genotypeCount[i];
            oldGenotype = i;
          }
        }
      }

      // Find the genotype to suggest, by the closest distance
      // keys = HashVec.getKeys(hash);
      // clusterCenters = new float[keys.length][2];
      clusterCenters = new float[4][2];
      distance = Float.MAX_VALUE;

      for (byte i = 3; i >= 0; i--) {
        iv = hash.get(i + "");
        if (iv != null) {
          xSum = 0;
          ySum = 0;
          genotypeIndices = Ints.toArray(iv);
          for (int j = 0; j < iv.size(); j++) {
            xSum = xSum + realX[genotypeIndices[j]];
            ySum = ySum + realY[genotypeIndices[j]];
          }
          clusterCenters[i] = new float[] {xSum / iv.size(), ySum / iv.size()};
          if (i != 3) {
            distancetemp =
                (float) Math.sqrt(Math.pow(clusterCenters[i][0] - clusterCenters[3][0], 2)
                    + Math.pow(clusterCenters[i][1] - clusterCenters[3][1], 2));
            if (distancetemp < distance) {
              distance = distancetemp;
              result = i;
            }
          }
        } else if (i != 3 && genotypeCount[i] > 0 && genotypeCount[i] == Math.max(genotypeCount[0],
            Math.max(genotypeCount[1], genotypeCount[2]))) {
          result = i;
          break;
        }
      }

      if (oldGenotype == result) {
        result = -1;
      }
    }

    return result;
  }

  public byte suggestedNewGenoTypeByNearbyPoint(MarkerData markerData) {
    float[] realX;
    float[] realY;
    float xSum, ySum;
    Vector<Integer> indexOfPointsOutsideTheCluster;
    float distancetemp = 0;
    byte[] genotypes;
    int[] counterGenotypeInsideTheCluster;
    int[] counterGenotypesOutsideTheCluster;
    float minDist;
    int indexOfNearbyPoint;
    int maxCount;
    byte genotypeOld;
    byte genotypeNew;
    // byte nGenotypesTotal;
    byte nGenotypesOutsideTheCluster;

    switch (getPlotType()) {
      case 0:
        // realX = markerData.getX_Raws();
        // realY = markerData.getY_Raws();
        // break;
        // case 1:
        realX = markerData.getXs();
        realY = markerData.getYs();
        break;
      case 1:
        realX = markerData.getThetas();
        realY = markerData.getRs();
        break;
      case 2:
        realX = markerData.getBAFs();
        realY = markerData.getLRRs();
        break;
      default:
        realX = markerData.getXs();
        realY = markerData.getYs();
    }

    genotypes = markerData.getAbGenotypes();
    // genotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerName,
    // gcThreshold);
    indexOfPointsOutsideTheCluster = new Vector<Integer>(genotypes.length);
    xSum = 0;
    ySum = 0;
    counterGenotypeInsideTheCluster = new int[] {0, 0, 0, 0};
    counterGenotypesOutsideTheCluster = new int[] {0, 0, 0, 0};
    genotypeOld = -2;
    for (int i = 0; i < genotypes.length; i++) {
      if (realX[i] >= rawXMin && realY[i] >= rawYMin && realX[i] <= rawXmax
          && realY[i] <= rawYmax) {
        xSum += realX[i];
        ySum += realY[i];
        if (genotypes[i] < 0) {
          genotypeOld = -1;
        } else {
          counterGenotypeInsideTheCluster[genotypes[i] + 1]++;
        }
      } else {
        indexOfPointsOutsideTheCluster.add(i);
        counterGenotypesOutsideTheCluster[genotypes[i] + 1]++;
      }
    }
    xSum = xSum / (genotypes.length - indexOfPointsOutsideTheCluster.size());
    ySum = ySum / (genotypes.length - indexOfPointsOutsideTheCluster.size());
    nGenotypesOutsideTheCluster = 0;
    // nGenotypesTotal = 0;
    maxCount = 0;
    if (genotypeOld == -2) {
      for (byte i = 1; i < counterGenotypeInsideTheCluster.length; i++) {
        if (counterGenotypeInsideTheCluster[i] > maxCount) {
          maxCount = counterGenotypeInsideTheCluster[i];
          genotypeOld = i;
        }
        if (counterGenotypesOutsideTheCluster[i] > 0) {
          // nGenotypesTotal ++;
          nGenotypesOutsideTheCluster++;
        } else if (counterGenotypeInsideTheCluster[i] > 0) {
          // nGenotypesTotal ++;
        }
      }
    }

    minDist = Float.MAX_VALUE;
    indexOfNearbyPoint = -1;
    for (int i = 0; i < indexOfPointsOutsideTheCluster.size(); i++) {
      distancetemp =
          (float) Math.sqrt(Math.pow(realX[indexOfPointsOutsideTheCluster.elementAt(i)] - xSum, 2)
              + Math.pow(realY[indexOfPointsOutsideTheCluster.elementAt(i)] - ySum, 2));
      if (distancetemp < minDist && genotypes[indexOfPointsOutsideTheCluster.elementAt(i)] != -1) {
        minDist = distancetemp;
        indexOfNearbyPoint = indexOfPointsOutsideTheCluster.elementAt(i);
      }
    }

    if (indexOfNearbyPoint == -1) {
      genotypeNew = -1;
    } else if (nGenotypesOutsideTheCluster == 3) {
      genotypeNew = genotypes[indexOfNearbyPoint];
      if (genotypeNew == genotypeOld) {
        genotypeNew = -1;
      }
    } else {
      genotypeNew = 1;
    }

    return genotypeNew;
  }

}
