package cnv.qc;

import java.io.*;
import java.util.*;

import javax.swing.JOptionPane;

import mining.Distance;
//import cnv.analysis.FilterCalls;
import cnv.filesys.ClusterFilterCollection;
import cnv.filesys.Sample;
import cnv.filesys.MarkerData;
import cnv.filesys.MarkerDataCollection;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.var.CNVariant;
import cnv.var.SampleData;
import common.*;

public class ThetaOutliers {
	public static final String[] QC_HEADS = {"LRR_mean", "LRR_median", "LRR_SD", "BAF_mean", "BAF_median", "BAF_SD", "BAF_DRIFT", "WF", "GCWF"};
	public static final String[] ERRORS = {"large SD for LRR", "drifting BAF values", "waviness factor values", "Small-sized CNV calls"};
	public static final String QC_SUMMARY_FILE = "Sample_QC.xln";


	/**
	 * Detect whether a null data point should be reclustered with a known genotype.
	 * @param proj The project where the data is.
	 * 
	 */
	public static void loadData(Project proj, String markerListfilename, byte stdDev) {
		PrintWriter writer;
		String[] sampleList, markerList = null, line;
		String filename, output;
		ClusterFilterCollection clusterFilterCollection;
		MarkerData[] markerData = null;
		int[] result;
		
		output = proj.getProjectDir()+"resultOfReclusterByTheta_sd"+stdDev+".txt";

		// load data: sample list
		filename = proj.getFilename(Project.SAMPLE_SUBSET_FILENAME, true, false);
		if (ext.rootOf(filename) == null || ext.rootOf(filename).equals("")) {
			sampleList = proj.getSampleList().getSamples();
		} else if (Files.exists(filename, proj.getJarStatus())) {
			System.out.print("filename: "+filename);
			sampleList = HashVec.loadFileToStringArray(filename, false, new int[] {0}, false);
		} else {
			JOptionPane.showMessageDialog(null, "Failed to load \""+filename+"\"", "Error", JOptionPane.ERROR_MESSAGE);
			return;
		}

		// load data: marker list and makerData
		if (markerListfilename == null) {
			//TODO !!! some users might move the .scat files, and thus markerList will not match the genotypes derived from markerData 
			markerList = proj.getMarkerSet().getMarkerNames();
	        String[] files = Files.list(proj.getDir(Project.PLOT_DIRECTORY), ".scat", false);
	        if ((new File(proj.getDir(Project.DATA_DIRECTORY)+"clusterFilters.ser")).exists()) {
	        	clusterFilterCollection = ClusterFilterCollection.load(proj.getDir(Project.DATA_DIRECTORY)+"clusterFilters.ser", proj.getJarStatus());
	        } else {
	        	clusterFilterCollection = null;
	        }
			try {
				writer = new PrintWriter(new FileWriter(output));
				writer.println("Name\tChr\tPosition\tindividualID");
				for (int i=0; i<files.length; i++) {
					markerData = MarkerDataCollection.load(proj.getDir(Project.PLOT_DIRECTORY)+files[i], proj.getJarStatus()).getCollection();
					for (int j=0; j<markerData.length; j++) {
						//TODO how to get
//						System.out.print("\n"+markerData[j].getMarkerName());
//						System.out.print("\t"+markerData[j].getChr());
//						System.out.print("\t"+markerData[j].getPosition());
//						result = reclusterNullGenotypeByTheta(markerData[j], clusterFilterCollection, 20);
						result = reclusterNullGenotypeByTheta(markerData[j], null, stdDev);
						for (int k=0; result!=null && k<result.length; k++) {
//							System.out.println(markerData[j].getMarkerName());
//							System.out.println(markerData[j].getChr());
//							System.out.println(markerData[j].getPosition());
//							System.out.println(result.length);
//							System.out.println(k);
//							System.out.println(sampleList.length);
							writer.println(markerData[j].getMarkerName()+"\t"+markerData[j].getChr()+"\t"+markerData[j].getPosition()+"\t"+sampleList[k]);
						}
					}
				}
				writer.close();
				System.out.println("Reclusterable Null genotypes file is now ready at: "+output);
			} catch (Exception e) {
				System.err.println("Error writing to '" + output + "'");
				e.printStackTrace();
			}
		} else if (Files.exists(markerListfilename, proj.getJarStatus())) {
			Vector<String> markerNames = new Vector<String>();
			try {
				BufferedReader reader = new BufferedReader(new FileReader(markerListfilename));
				while (reader.ready()) {
					line = reader.readLine().trim().split("\t", -1);
					if (!markerNames.contains(line[0])) {
						markerNames.add(line[0]);
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				JOptionPane.showMessageDialog(null, "Error - could not find \""+filename+"\"", "Error", JOptionPane.ERROR_MESSAGE);
			} catch (Exception e) {
				System.err.println("Error reading file \""+filename+"\"");
				System.exit(2);
			}
			System.out.println("sampleList.length="+sampleList.length+"\t"+"markerList.length="+markerList.length);
			markerList = Array.toStringArray(markerNames);
			markerData = MarkerSet.loadFromList(proj, markerList);
		} else {
			JOptionPane.showMessageDialog(null, "Failed to load \""+markerListfilename+"\"", "Error", JOptionPane.ERROR_MESSAGE);
			return;
		}
	}

	public static int[] reclusterNullGenotypeByTheta(MarkerData markerData, ClusterFilterCollection clusterFilterCollection, int numberOfStdDev) {
		byte[] genotypes;
		double[] rs, thetas;
		DoubleVector[] rsByGenotype, thetasByGenotype;
		double[] meanR, meanTheta, sdTheta;
		ByteVector nonMissingGenotype;
		double[] distance;
		byte n;
		IntVector result;

		if (clusterFilterCollection != null) {
			genotypes = markerData.getAB_GenotypesAfterFilters(clusterFilterCollection, markerData.getMarkerName(), 0);
		} else {
			genotypes = markerData.getAB_Genotypes();
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
		Project proj = new Project("C:/workspace/Genvisis/projects/GEDI_exome.properties", false);
		loadData(proj, null, (byte) 12);

//		MarkerData[] markers = MarkerSet.loadFromList(proj, new String[] {"rs17080321", "rs7898873", "rs17080321", "rs7898873"});
//		MarkerData[] markers = MarkerSet.loadFromList(proj, new String[] {"rs17246013", "rs34771052", "rs17080321", "rs7898873", "rs17231443", "rs2227433", "rs9907972", "rs34148246", "rs11572080", "rs34942735", "rs4646168"});
//		for (int i = 0; i < markers.length; i++) {
//			int[] results = reclusterNullGenotypeByTheta(markers[i], null);
//			System.out.println(markers[i].getMarkerName()+"\t("+results.length+")\t"+Array.toStr(results));
//		}
//		
	}
}
