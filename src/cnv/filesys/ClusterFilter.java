package cnv.filesys;

import java.io.Serializable;
import java.util.Hashtable;
import java.util.Vector;

import common.IntVector;

/**
 * This is a data structure to hold a single filter that can be used to screen the data points. 
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
	
	public ClusterFilter () {
	}
	
	public ClusterFilter (byte plotType, float rawXMin, float rawYMin, float rawXmax, float rawYmax) {
		this(plotType, rawXMin, rawYMin, rawXmax, rawYmax, (byte) 0);
	}
	
	public ClusterFilter (byte plotType, float rawXMin, float rawYMin, float rawXmax, float rawYmax, byte newGenotype) {
		this.plotType = plotType;
		this.rawXMin = rawXMin;
		this.rawYMin = rawYMin;
		this.rawXmax = rawXmax;
		this.rawYmax = rawYmax;
		this.newGenotype = newGenotype;
	}
	
	public ClusterFilter (byte plotType, float rawXMin, float rawYMin, float rawXmax, float rawYmax, MarkerData markerData) {
		this.plotType = plotType;
		this.rawXMin = rawXMin;
		this.rawYMin = rawYMin;
		this.rawXmax = rawXmax;
		this.rawYmax = rawYmax;
//		this.newGenotype = suggestedNewGenoTypeByNearbyCentroid(markerData);
		this.newGenotype = suggestedNewGenoTypeByNearbyPoint(markerData);
	}
	
	public byte getPlotType () {
		return this.plotType;
	}
	
	public void setClusterGenotype (byte newGenotype) {
		this.newGenotype=newGenotype;
	}
	
	public byte getCluterGenotype () {
		return this.newGenotype;
	}
	
	public float getXMin () {
		return rawXMin;
	}

	public float getYMin () {
		return rawYMin;
	}

	public float getXMax () {
		return rawXmax;
	}

	public float getYMax () {
		return rawYmax;
	}

	public byte suggestedNewGenoTypeByNearbyCentroid (MarkerData markerData) {
		float[] realX;
		float[] realY;
		byte result=-1;
		float xSum, ySum;
		Hashtable<String,IntVector> hash;
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
		
		switch(getPlotType()) {
		case 0:
			realX = markerData.getX_Raws();
			realY = markerData.getY_Raws();
			break;
		case 1:
			realX = markerData.getXs();
			realY = markerData.getYs();
			break;
		case 2:
			realX = markerData.getThetas();
			realY = markerData.getRs();
			break;
		case 3:
			realX = markerData.getBAFs();
			realY = markerData.getLRRs();
			break;
		default:
			realX = markerData.getXs();
			realY = markerData.getYs();
		}

		hash = new Hashtable<String, IntVector>();
		genotypes = markerData.getAB_Genotypes();
		// iterate through all samples
		for (int i=0; i<genotypes.length; i++) {
			if (       realX[i] >= rawXMin
					&& realY[i] >= rawYMin
					&& realX[i] <= rawXmax
					&& realY[i] <= rawYmax) {
				cluster = "3";	//"3" identifies the data points within the cluster filter
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
		if (iv!=null) {
			/*
			 * Search for the genotype of majority of the data points
			 * If there is any data point with missing value found, the majority's genotype is set to Missing Value. 
			 */
			genotypeCount = new int[] {0,0,0};
			oldGenotype=-2;
			for (int i=0; i<3; i++){
				genotypeIndices = iv.toArray();
				for (int j=0; j<genotypeIndices.length; j++) {
					if (genotypes[genotypeIndices[j]]==i) {
						genotypeCount[i]++;
					} else if (genotypes[genotypeIndices[j]]==-1) {
						oldGenotype=-1;
					}
				}
			}
			// Search for the maximum in int[] genotypeCount 
			if (oldGenotype==-2) {
				genotypeFrequencyCount=-1;
				for (byte i=0; i<3; i++){
					if (genotypeCount[i]>genotypeFrequencyCount) {
						genotypeFrequencyCount=genotypeCount[i];
						oldGenotype=i;
					}
				}
			}

			// Find the genotype to suggest, by the closest distance  
//			keys = HashVec.getKeys(hash);
//			clusterCenters = new float[keys.length][2];
			clusterCenters = new float[4][2];
			distance = Float.MAX_VALUE;

			for (byte i = 3; i>=0; i--) {
				iv = hash.get(i+"");
				if (iv!=null) {
					xSum = 0;
					ySum = 0;
					genotypeIndices = iv.toArray();
					for (int j=0; j<iv.size(); j++) {
						xSum = xSum + realX[genotypeIndices[j]];
						ySum = ySum + realY[genotypeIndices[j]];
					}
					clusterCenters[i] = new float[] {xSum/iv.size(), ySum/iv.size()};
					if (i!=3) {
						distancetemp = (float) Math.sqrt(Math.pow(clusterCenters[i][0]-clusterCenters[3][0],2)+Math.pow(clusterCenters[i][1]-clusterCenters[3][1],2));
						if (distancetemp < distance) {
							distance = distancetemp;
							result = i;
						}
					}
				} else if (i!=3 && genotypeCount[i]>0 && genotypeCount[i]==Math.max(genotypeCount[0], Math.max(genotypeCount[1], genotypeCount[2]))) {
					result = i;
					break;
				}
			}

			if (oldGenotype == result) {
				result=-1;
			}
		}
		
		return result;
	}

	public byte suggestedNewGenoTypeByNearbyPoint (MarkerData markerData) {
		float[] realX;
		float[] realY;
		float xSum, ySum;
		Vector<Integer> indexOfPointsOutsideTheCluster;
		float distancetemp = 0;
		byte[] genotypes;
		int[] genotypeCount;
		float minDist;
		int indexOfNearbyPoint;
		int maxCount;
		byte oldGenotype;
		byte newGenotype;
		
		switch(getPlotType()) {
		case 0:
			realX = markerData.getX_Raws();
			realY = markerData.getY_Raws();
			break;
		case 1:
			realX = markerData.getXs();
			realY = markerData.getYs();
			break;
		case 2:
			realX = markerData.getThetas();
			realY = markerData.getRs();
			break;
		case 3:
			realX = markerData.getBAFs();
			realY = markerData.getLRRs();
			break;
		default:
			realX = markerData.getXs();
			realY = markerData.getYs();
		}

		indexOfPointsOutsideTheCluster = new Vector<Integer>();
		genotypes = markerData.getAB_Genotypes();
		xSum = 0;
		ySum = 0;
		genotypeCount = new int[] {0,0,0};
		oldGenotype = -2;
		for (int i=0; i<genotypes.length; i++) {
			if (realX[i] >= rawXMin && realY[i] >= rawYMin && realX[i] <= rawXmax && realY[i] <= rawYmax) {
				xSum += realX[i];
				ySum += realY[i];
				if (genotypes[i] < 0) {
					oldGenotype = -1;
				} else {
					genotypeCount[genotypes[i]] ++;
				}
			} else {
				indexOfPointsOutsideTheCluster.add(i);
			}
		}
		xSum = xSum / (genotypes.length - indexOfPointsOutsideTheCluster.size());
		ySum = ySum / (genotypes.length - indexOfPointsOutsideTheCluster.size());
		maxCount = 0;
		if (oldGenotype == -2) {
			for (byte i = 0; i < 3; i++){
				if (genotypeCount[i] > maxCount) {
					maxCount = genotypeCount[i];
					oldGenotype = i;
				}
			}
		}

		minDist = Float.MAX_VALUE;
		indexOfNearbyPoint = -1;
		for (int i = 0; i < indexOfPointsOutsideTheCluster.size(); i++) {
			distancetemp = (float) Math.sqrt(Math.pow(realX[indexOfPointsOutsideTheCluster.elementAt(i)] - xSum, 2) + Math.pow(realY[indexOfPointsOutsideTheCluster.elementAt(i)] - ySum, 2));
			if (distancetemp < minDist && genotypes[indexOfPointsOutsideTheCluster.elementAt(i)] != -1) {
				minDist = distancetemp;
				indexOfNearbyPoint = indexOfPointsOutsideTheCluster.elementAt(i);
			}
		}

		if (indexOfNearbyPoint == -1) {
			newGenotype = -1;
		} else {
			newGenotype = genotypes[indexOfNearbyPoint];
			if (newGenotype == oldGenotype) {
				newGenotype = -1;
			}
		}

		return newGenotype;
	}
}
