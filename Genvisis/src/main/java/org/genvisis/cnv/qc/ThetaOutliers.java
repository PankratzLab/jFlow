package org.genvisis.cnv.qc;

import java.io.*;
import java.util.*;

import org.genvisis.cnv.filesys.ClusterFilterCollection;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.MarkerDataLoader;
import org.genvisis.common.*;
import org.genvisis.mining.Distance;

public class ThetaOutliers {
	public static final String[] QC_HEADS = {"LRR_mean", "LRR_median", "LRR_SD", "BAF_mean", "BAF_median", "BAF_SD", "BAF_DRIFT", "WF", "GCWF"};
	public static final String[] ERRORS = {"large SD for LRR", "drifting BAF values", "waviness factor values", "Small-sized CNV calls"};
	public static final String QC_SUMMARY_FILE = "Sample_QC.xln";


	/**
	 * Detect whether a null data point should be reclustered with a known genotype.
	 * @param proj The project where the data is.
	 * 
	 */
	public static void loadData(Project proj, boolean useClusterFilters, byte stdDev) {
		PrintWriter writer;
		String[] sampleList;
		String filename, output;
		ClusterFilterCollection clusterFilterCollection;
		int[] result;
		MarkerDataLoader markerDataLoader;
		MarkerData markerData;
		String[] markerNames;
		long time;
		Logger log;
		
		log = proj.getLog();
		output = proj.PROJECT_DIRECTORY.getValue()+"resultOfReclusterByTheta_sd"+stdDev+".txt";

		// load data: sample list
		filename = proj.SAMPLE_SUBSET_FILENAME.getValue(true, false);
		if (ext.rootOf(filename) == null || ext.rootOf(filename).equals("")) {
			sampleList = proj.getSampleList().getSamples();
		} else if (Files.exists(filename, proj.JAR_STATUS.getValue())) {
			log.report("filename: "+filename);
			sampleList = HashVec.loadFileToStringArray(filename, false, new int[] {0}, false);
		} else {
			proj.message("Failed to load \""+filename+"\"");
			return;
		}

		try {
			writer = new PrintWriter(new FileWriter(output));
			writer.println("Name\tChr\tPosition\tindividualID");

	        time = new Date().getTime();
	        markerNames = proj.getMarkerNames();
			markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames);
	        clusterFilterCollection = proj.getClusterFilterCollection();

	        for (int i = 0; i < markerNames.length; i++) {
	        	markerData = markerDataLoader.requestMarkerData(i);
	        	if (i % 100 == 0) {
	        		log.report(ext.getTime()+"\tMarker "+i+" of "+markerNames.length);
	        	}

	        	if (useClusterFilters) {
	        		result = reclusterNullGenotypeByTheta(markerData, clusterFilterCollection, 20, log);
	        	} else {
	        		result = reclusterNullGenotypeByTheta(markerData, null, stdDev, log);
	        	}
				for (int k=0; result!=null && k<result.length; k++) {
					writer.println(markerData.getMarkerName()+"\t"+markerData.getChr()+"\t"+markerData.getPosition()+"\t"+sampleList[k]);
				}
				markerDataLoader.releaseIndex(i);
			}
			writer.close();
			log.report("Reclusterable Null genotypes file is now ready at: "+output);
			log.report("Finished searching for ThetaOutliers in " + ext.getTimeElapsed(time));
		} catch (Exception e) {
			log.reportError("Error writing to '" + output + "'");
			log.reportException(e);
		}
	}

	public static int[] reclusterNullGenotypeByTheta(MarkerData markerData, ClusterFilterCollection clusterFilterCollection, int numberOfStdDev, Logger log) {
		byte[] genotypes;
		double[] rs, thetas;
		DoubleVector[] rsByGenotype, thetasByGenotype;
		double[] meanR, meanTheta, sdTheta;
		ByteVector nonMissingGenotype;
		double[] distance;
		byte n;
		IntVector result;

		if (clusterFilterCollection != null) {
			genotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerData.getMarkerName(), 0, log);
		} else {
			genotypes = markerData.getAbGenotypes();
		}
//		for (int i=0; i<genotypes.length; i++) {
//			if (genotypes[i]>-1) {
//				System.out.println("\t"+genotypes.length);
//			}
//		}
//		rs = Array.normalize(Array.removeNaN(Array.toDoubleArray(markerData.getRs())));
//		thetas = Array.normalize(Array.removeNaN(Array.toDoubleArray(markerData.getThetas())));
//		rs = Array.removeNaN(Array.toDoubleArray(markerData.getRs()));
//		thetas = Array.removeNaN(Array.toDoubleArray(markerData.getThetas()));
		rs = Array.toDoubleArray(markerData.getRs());
		thetas = Array.toDoubleArray(markerData.getThetas());
		rsByGenotype = new DoubleVector[4];
		thetasByGenotype = new DoubleVector[4];
		for (int i=0; i<4; i++) {
			rsByGenotype[i] = new DoubleVector();
			thetasByGenotype[i] = new DoubleVector();
		}
//		if (genotypes.length!=rs.length) {
//			System.out.println(markerData.getMarkerName());
//		} else {
			for (int i=0; i<genotypes.length; i++) {
				rsByGenotype[genotypes[i]+1].add(rs[i]);
				thetasByGenotype[genotypes[i]+1].add(thetas[i]);
			}
//		}
		nonMissingGenotype = new ByteVector();
		meanR = new double[4];
		meanTheta = new double[4];
		sdTheta = new double[4];
		distance = new double[4];
		for (int i=1; i<4; i++) {
			if (rsByGenotype[i].size() >= 5) {
				nonMissingGenotype.add((byte)i);
				meanR[i] = Array.mean(rsByGenotype[i].toArray());
				meanTheta[i] = Array.mean(thetasByGenotype[i].toArray());
				sdTheta[i] = Array.stdev(thetasByGenotype[i].toArray());
//				if (markerData.getMarkerName().equals("rs2139063") || markerData.getMarkerName().equals("rs35687686")) {
//					System.out.println(markerData.getMarkerName()+"\t"+meanR[i]);
//				}
			}
		}
		if (nonMissingGenotype.size()<3 && rsByGenotype[0].size()>0) {
			result = new IntVector();
			for (int j=0; nonMissingGenotype.size()>0 && j<rsByGenotype[0].size(); j++) {
				if (nonMissingGenotype.size()==1) {
					n = nonMissingGenotype.elementAt((byte)0);
				} else {
					for (byte i=1; i<4; i++) {
						if (rsByGenotype[i].size()!=0) {
							//TODO actually need to standardize the R and Theta
							distance[i] = Distance.euclidean(new double[] {rsByGenotype[0].elementAt(j), thetasByGenotype[0].elementAt(j)}, new double[] {meanR[i], meanTheta[i]});
						}
					}
					if (distance[nonMissingGenotype.elementAt((byte)0)] < distance[nonMissingGenotype.elementAt((byte)1)] ) {
						n = nonMissingGenotype.elementAt((byte)0);
					} else {
						n = nonMissingGenotype.elementAt((byte)1);
					}
				}
				if (rsByGenotype[0].elementAt(j)>0.5*meanR[n] && (thetasByGenotype[0].elementAt(j)>=(meanTheta[n]+numberOfStdDev*sdTheta[n]) || thetasByGenotype[0].elementAt(j)<=(meanTheta[n]-numberOfStdDev*sdTheta[n]))) {
//				if (   rsByGenotype[0].elementAt(j)>0.5*meanR[n]
//					&& (   (thetasByGenotype[0].elementAt(j)<(meanTheta[n]+15*sdTheta[n]) && thetasByGenotype[0].elementAt(j)>=(meanTheta[n]+10*sdTheta[n]))
//						|| (thetasByGenotype[0].elementAt(j)<(meanTheta[n]-15*sdTheta[n]) && thetasByGenotype[0].elementAt(j)>=(meanTheta[n]-10*sdTheta[n]))
//					   )
//					) {
					//TODO assign the new genotype. Currently just highlighting
					result.add(j);
//					System.out.println(rsByGenotype[0].elementAt(j)+","+thetasByGenotype[0].elementAt(j)+" where mean R is "+meanR[n]+" and theta is "+meanTheta[n]+" +/- "+sdTheta[n]+" (n="+rsByGenotype[n].size()+")");
				}
			}
			return result.toArray();
		} else {
			return null;
		}
	}
	
	public static void main(String[] args) {
		String filename;
		Project proj;

		filename = "D:/home/npankrat/projects/GEDI_exomeRAF.properties";
//		Project proj = new Project("C:/workspace/Genvisis/projects/GEDI_exome.properties", false);
		proj = new Project(filename, false);
		loadData(proj, true, (byte) 12);

//		MarkerData[] markers = MarkerSet.loadFromList(proj, new String[] {"rs17080321", "rs7898873", "rs17080321", "rs7898873"});
//		MarkerData[] markers = MarkerSet.loadFromList(proj, new String[] {"rs17246013", "rs34771052", "rs17080321", "rs7898873", "rs17231443", "rs2227433", "rs9907972", "rs34148246", "rs11572080", "rs34942735", "rs4646168"});
//		for (int i = 0; i < markers.length; i++) {
//			int[] results = reclusterNullGenotypeByTheta(markers[i], null);
//			System.out.println(markers[i].getMarkerName()+"\t("+results.length+")\t"+Array.toStr(results));
//		}
//		
	}
}
