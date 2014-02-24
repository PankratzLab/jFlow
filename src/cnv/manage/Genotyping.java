package cnv.manage;

// TODO this is an experimental concept that should not be released
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import cnv.filesys.MarkerData;
import cnv.filesys.Project;

import common.Elision;
import common.Logger;

public class Genotyping {
//	private String sampleId;
//	private float[] x;
//	private float[] y;
//	private byte[] genotypeOld;
	private byte[] genotypeNew;
	private byte nClusters;
	private double[][] centroids;
	private double[][] marginsMin;
	private double[][] clusterSizes;
	private String markerAnnotations;

//	public GenotypeClustering (float[] x, float[] y, byte[] genotypeOld, byte[] genotypeNew, byte nClusters, double[][] centroids, double[][] marginsMin, double[][] clusterSizes, String markerAnnotations) {
	public Genotyping (byte[] genotypeNew, byte nClusters, double[][] centroids, double[][] marginsMin, double[][] clusterSizes, String markerAnnotations) {
//		this.sampleId = sampleId;
//		this.x = x;
//		this.y = y;
//		this.genotypeOld = genotypeOld;
		this.genotypeNew = genotypeNew;
		this.nClusters = nClusters;
		this.centroids = centroids;
		this.marginsMin = marginsMin;
		this.clusterSizes = clusterSizes;
		this.markerAnnotations = markerAnnotations;
	}

	public byte[] getGenotypeNew() {
		return genotypeNew;
	}

	public byte getClusters() {
		return nClusters;
	}

	public double[][] getCentroids() {
		return centroids;
	}

	public double[][] getMarginMin() {
		return marginsMin;
	}

	public double[][] getClusterSizes() {
		return clusterSizes;
	}

	public String getMarkerAnnotations() {
		return markerAnnotations;
	}

	public static Genotyping clusteringAMarker(float[] x, float[] y) {
		int MAX_N_CLUSTERS = new String[] {"-/-", "A/-", "B/-", "A/A", "A/B", "B/B", "AAA", "AAB", "ABB", "BBB"}.length;
		int N_FOR_DISPERSE_TEST = 20;

//		byte[] clusterLabels;
//		byte nClusters;
//		double[][] clusterCentroids;
//		double[][] clusterMargins;
//		double[][] clusterSizes;
//		String annotations;
		double[] distOfMins;
		int[] indexOfMins;
		int[] sortedIndex;
		double distCurrent;
		int nLowerMeans;
		double lowerMean;
//		double midMean;

		try {
			if (x.length != y.length) {
				throw new Elision("Error - Clustering: the lengths of the two arrays have to be equal.");
			}
		} catch (Elision e) {
			e.printStackTrace();
		}

		distOfMins = new double[x.length];
		indexOfMins = new int[x.length];
		for (int i = 0; i < x.length; i ++) {
			distOfMins[i] = Double.MAX_VALUE;
			indexOfMins[i] = -1;
			for (int j = i + 1; j < y.length; j ++) {
				distCurrent = Math.sqrt(Math.pow(x[i] - x[j], 2) + Math.pow(y[i] - y[j], 2));
				if (distCurrent < distOfMins[i]) {
					distOfMins[i] = distCurrent;
					indexOfMins[i] = j;
				}
			}
		}

		distCurrent = Double.MAX_VALUE;
		sortedIndex = getSortedIndex(distOfMins);

		lowerMean = 0;
		nLowerMeans = x.length - MAX_N_CLUSTERS - N_FOR_DISPERSE_TEST;
		for (int i = 0; i < nLowerMeans; i ++) {
			lowerMean += distOfMins[sortedIndex[i]];
		}
		lowerMean = lowerMean / nLowerMeans;
		
//		clusterLabels = new byte[x.length];
		
//		return new Genotyping(clusterLabels, nClusters, clusterCentroids, clusterMargins, clusterSizes, annotations);
		return null;
	}
	
	public static int[] getSortedIndex(double[] data) {
		int[] result;
		double[] dataTmp;
		double swap;
		
		dataTmp = new double[data.length];
		for (int i = 0; i < dataTmp.length; i++) {
			dataTmp[i] = data[i];
		}
		
		result = new int[data.length];
		for (int i = 0; i < result.length; i++) {
			result[i] = i;
		}
		
		for (int i = 0; i < dataTmp.length; i++) {
			for (int j = i; j < dataTmp.length; j++) {
				if (dataTmp[j] < dataTmp[i]) {
					swap = dataTmp[i];
					dataTmp[i] = dataTmp[j];
					dataTmp[j] = swap;
					result[i] = j;
				}
			}
		}
		return result; 
	}

	public static byte[] putIntoClusters(int[] indexOfMins, int[] sortedIndex, int startingIndex, int endingIndex) {
		byte[] clusterLabels;
		byte nextAvailableLabel;
//		Vector<Vector<Integer>> result;

		clusterLabels = new byte[indexOfMins.length];
		for (int i = 0; i < clusterLabels.length; i++) {
			clusterLabels[i] = (byte) -1;
		}
		if (startingIndex < 0) {
			startingIndex = 0;
		}
		if (endingIndex < 0) {
			endingIndex = indexOfMins.length;
		}

		nextAvailableLabel = 0;
		for (int i = startingIndex; i < endingIndex; i++) {
			if (clusterLabels[sortedIndex[i]] == -1 && clusterLabels[indexOfMins[sortedIndex[i]]] == -1) {
				clusterLabels[sortedIndex[i]] = nextAvailableLabel;
				clusterLabels[indexOfMins[sortedIndex[i]]] = nextAvailableLabel;
				nextAvailableLabel ++;
			} else if (clusterLabels[sortedIndex[i]] == -1) {
				clusterLabels[sortedIndex[i]] = clusterLabels[indexOfMins[sortedIndex[i]]];
			} else {
				clusterLabels[indexOfMins[sortedIndex[i]]] = clusterLabels[sortedIndex[i]];
			}
		}
		
		return clusterLabels;
	}
	
	public static void test (String projFilename, String markerName, String outFileName) {
		Project proj;
		MarkerData markerData;
		byte[] genotypeNew;
		float[] x, y;
		byte[] genotypeOld;
		PrintWriter out;

		proj = new Project(projFilename, false);
		markerData = new MarkerDataLoader(proj, new String[] {markerName}, 1, new Logger()).getMarkerData(0);
		x = markerData.getXs();
		y = markerData.getYs();
		genotypeOld = markerData.getAB_Genotypes();
		genotypeNew = clusteringAMarker(x, y).getGenotypeNew();
		
		try {
			out = new PrintWriter(new FileOutputStream(proj.getDir(Project.DATA_DIRECTORY) + outFileName));
			for (int i = 0; i < genotypeOld.length; i++) {
				out.println(i + "\t" + x[i] + "\t" + y[i] + "\t" + genotypeOld[i] + "\t" + genotypeNew[i]);
			}
			out.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	public static void main (String args[]) {
		String projFilename = "C:/workspace/Genvisis/projects/GEDI_exome.properties";
		String markerName = "exm221098";
		String outFileName = "newGenotype_" + markerName + ".txt";

		test(projFilename, markerName, outFileName);
	}
}
