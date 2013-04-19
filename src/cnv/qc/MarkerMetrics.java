package cnv.qc;

import java.io.*;
import java.util.*;

import stats.Ttest;
import cnv.filesys.*;
import cnv.manage.MarkerDataLoaderRunnable;
import cnv.var.SampleData;
import common.*;
import db.FilterDB;

public class MarkerMetrics {
	public static final String[] FULL_QC_HEADER = {"MarkerName", "Chr", "CallRate", "meanTheta_AA", "meanTheta_AB", "meanTheta_BB", "diffTheta_AB-AA", "diffTheta_BB-AB", "sdTheta_AA", "sdTheta_AB", "sdTheta_BB", "meanR_AA", "meanR_AB", "meanR_BB", "num_AA", "num_AB", "num_BB", "pct_AA", "pct_AB", "pct_BB", "MAF", "HetEx"};
	
	public static void fullQC(Project proj, String fileWithListOfSamplesToUse, String subsetOfMarkers, Logger log) {
		Hashtable<String, String> hash;
        String[] samples;
        boolean[] samplesToUse;

        if (fileWithListOfSamplesToUse != null) {
        	hash = HashVec.loadFileToHashString(fileWithListOfSamplesToUse, false);
	        samples = proj.getSamples();
	        samplesToUse = new boolean[samples.length];
	        for (int i = 0; i < samples.length; i++) {
	        	samplesToUse[i] = hash.containsKey(samples[i]);
			}
        } else {
        	samplesToUse = null;
        }
		
        fullQC(proj, samplesToUse, subsetOfMarkers, log);
	}
	
	public static void fullQC(Project proj, boolean[] samplesToExclude, String markersToInclude, Logger log) {
		PrintWriter writer;
		String[] samples;
		float[] thetas, rs;
		MarkerData markerData;
        byte[] abGenotypes;
        String markerName;
        ClusterFilterCollection clusterFilterCollection;
        float gcThreshold;
        long time;
        MarkerDataLoaderRunnable markerDataLoader;
        String[] markerList;
        String line, eol;
//        int countAA, countAB, countBB, countNull;
        int[] count;
        double[] sumTheta, sumR, meanTheta, sdTheta;
        double temp;
        
        if (System.getProperty("os.name").startsWith("Windows")) {
        	eol = "\r\n";
		} else {
			eol = "\n";
		}
        
        samples = proj.getSamples();
 		
        if ((new File(proj.getDir(Project.DATA_DIRECTORY)+"clusterFilters.ser")).exists()) {
        	clusterFilterCollection = ClusterFilterCollection.load(proj.getDir(Project.DATA_DIRECTORY)+"clusterFilters.ser", proj.getJarStatus());
        } else {
        	clusterFilterCollection = null;
        }
		
        gcThreshold = Float.parseFloat(proj.getProperty(Project.GC_THRESHOLD));

		try {
			writer = new PrintWriter(new FileWriter(proj.getDir(Project.RESULTS_DIRECTORY)+"markerGenderChecks.xln"));
			writer.println(Array.toStr(FULL_QC_HEADER));
			
			if (markersToInclude != null) {
				markerList = HashVec.loadFileToStringArray(proj.getDir(Project.PROJECT_DIRECTORY)+markersToInclude, false, new int[] {0}, false);
			} else {
				markerList = proj.getMarkerNames();
			}
			markerDataLoader = MarkerDataLoaderRunnable.loadMarkerDataFromList(proj, markerList);
			line = "";
			time = new Date().getTime();
			for (int i = 0; i < markerList.length; i++) {
				markerData = markerDataLoader.requestMarkerData(i);

				markerName = markerData.getMarkerName();
				thetas = markerData.getThetas();
				rs = markerData.getRs();
				abGenotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerName, gcThreshold);
				
				
				// this is where all code will go
				count = new int[4];
				sumTheta = new double[count.length];
				sumR = new double[count.length];
				for (int j = 0; j < samples.length; j++) {
					if (samplesToExclude==null || samplesToExclude[j]) {
						count[abGenotypes[j] + 1] ++;
						sumTheta[abGenotypes[j] + 1] += thetas[j];
						sumR[abGenotypes[j] + 1] += rs[j];
					}
				}

				meanTheta = new double[count.length];
				for (int j = 1; j < meanTheta.length; j++) {
					meanTheta[j] = sumTheta[j] / count[j];
				}
				sdTheta = new double[count.length];
				for (int j = 0; j < samples.length; j++) {
					if (samplesToExclude==null || samplesToExclude[j]) {
						temp = (thetas[j] - meanTheta[ abGenotypes[j] + 1 ]);
						sdTheta[ abGenotypes[j] + 1 ] +=  temp * temp;
					}
				}
				for (int j = 1; j < sdTheta.length; j++) {
					if (count[j] == 0) {
						sdTheta[j] = Double.NaN;	
					} else {
						sdTheta[j] = Math.sqrt(sdTheta[j] / (count[j]-1));
					}
				}

				line += (markerName
						+ "\t" + markerData.getChr()
						+ "\t" + (1- ((float)count[0] / (count[0] + count[1] + count[2] + count[3])))
						+ "\t" + meanTheta[1]
						+ "\t" + meanTheta[2]
						+ "\t" + meanTheta[3]
						+ "\t" + (meanTheta[2] - meanTheta[1])
						+ "\t" + (meanTheta[3] - meanTheta[2])
						+ "\t" + sdTheta[1]
						+ "\t" + sdTheta[2]
						+ "\t" + sdTheta[3]
						+ "\t" + (sumR[1] / count[1])
						+ "\t" + (sumR[2] / count[2])
						+ "\t" + (sumR[3] / count[3])
						+ "\t" + count[1]
						+ "\t" + count[2]
						+ "\t" + count[3]
						+ "\t" + ((float) count[1] / (count[0] + count[1] + count[2] + count[3]))
						+ "\t" + ((float) count[2] / (count[0] + count[1] + count[2] + count[3]))
						+ "\t" + ((float) count[3] / (count[0] + count[1] + count[2] + count[3]))
						+ "\t" + (float) (count[1]<count[3]? (count[1] + count[2]) : (count[2] + count[3])) / (count[0] + count[1] + 2 * count[2] + count[3])
						+ "\t" + AlleleFreq.HetExcess(count[1], count[2], count[3])[0]
						+ eol
						);
				
				if (line.length() > 25000) {
					writer.print(line);
					writer.flush();
					line = "";
				}
			}
			writer.print(line);
			writer.close();
			log.report("Finished analyzing "+markerList.length+" in "+ext.getTimeElapsed(time));

		} catch (Exception e) {
			log.reportError("Error writing results");
			e.printStackTrace();
		}
	}
	
	private static void separationOfSexes(Project proj, String subset, Logger log) {
		PrintWriter writer;
		String[] samples;
		float[] xs, ys;
		MarkerData markerData;
        SampleData sampleData;
        int[] sexes;
        byte[] abGenotypes;
        String markerName;
        ClusterFilterCollection clusterFilterCollection;
        float gcThreshold;
        long time;
        DoubleVector[][][] values; // x/y, genotype, sex
        double femaleComp;
        double[][] zScores, tVals, counts;
        double[] maleValues, femaleValues;
        double zMean, tMean, zMin, tMin, count;
        double maf, countBs;
        MarkerDataLoaderRunnable markerDataLoader;
        String[] markerList;
        String line, eol;
        
        if (System.getProperty("os.name").startsWith("Windows")) {
        	eol = "\r\n";
		} else {
			eol = "\n";
		}
        
        sampleData = proj.getSampleData(2, false);
        samples = proj.getSamples();
        sexes = new int[samples.length];
        for (int i = 0; i < samples.length; i++) {
        	sexes[i] = Math.max(0, sampleData.getSexForIndividual(samples[i]));
		}
 		
        if ((new File(proj.getDir(Project.DATA_DIRECTORY)+"clusterFilters.ser")).exists()) {
        	clusterFilterCollection = ClusterFilterCollection.load(proj.getDir(Project.DATA_DIRECTORY)+"clusterFilters.ser", proj.getJarStatus());
        } else {
        	clusterFilterCollection = null;
        }
		
        gcThreshold = Float.parseFloat(proj.getProperty(Project.GC_THRESHOLD));

		try {
			writer = new PrintWriter(new FileWriter(proj.getDir(Project.RESULTS_DIRECTORY)+"markerGenderChecks.xln"));
			writer.println("SNP\tX_zAA\tY_zBB\tX_tAA\tY_tBB\tMean_Zs\tMean_Ts\tMin_Z\tMin_T\tMAF");
			
			if (subset != null) {
				markerList = HashVec.loadFileToStringArray(proj.getDir(Project.PROJECT_DIRECTORY)+subset, false, new int[] {0}, false);
			} else {
				markerList = proj.getMarkerNames();
			}
			markerDataLoader = MarkerDataLoaderRunnable.loadMarkerDataFromList(proj, markerList);
			time = new Date().getTime();
			line = "";
			for (int i = 0; i < markerList.length; i++) {
				markerData = markerDataLoader.requestMarkerData(i);

				markerName = markerData.getMarkerName();
				xs = markerData.getXs();
				ys = markerData.getYs();
				abGenotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerName, gcThreshold);
				
				values = new DoubleVector[3][4][3]; // x/y/r, genotype, sex
				for (int s = 0; s < samples.length; s++) {
					if (ext.isValidDouble(xs[s]+"")) {
						if (values[0][abGenotypes[s]+1][sexes[s]] == null) {
							values[0][abGenotypes[s]+1][sexes[s]] = new DoubleVector();
							values[1][abGenotypes[s]+1][sexes[s]] = new DoubleVector();
							values[2][abGenotypes[s]+1][sexes[s]] = new DoubleVector();
						}
						values[0][abGenotypes[s]+1][sexes[s]].add(xs[s]);
						values[1][abGenotypes[s]+1][sexes[s]].add(ys[s]);
						values[2][abGenotypes[s]+1][sexes[s]].add(xs[s]+ys[s]);
					}
				}
				
				zScores = new double[3][3];
				tVals = new double[3][3];
				counts = new double[3][3];
				for (int k = 0; k < values.length; k++) {
					for (int ab = 1; ab < values[k].length; ab++) {
						if (values[k][ab][1] != null) {
							maleValues = values[k][ab][1].toArray();
							counts[k][ab-1] += maleValues.length;
//							System.out.println(markerName+"\t"+FullSample.AB_PAIRS[ab-1]+"/"+(k==0?"X":"Y")+"\t"+maleValues.length+" males\t"+Array.mean(maleValues));
							
							// use all the females with a missing genotype for comparison, if any females have this genotype, only use that group if they outnumber those with a missing genotype 
							femaleValues = null;
							if (values[k][ab][2] != null) {
								femaleValues = values[k][ab][2].toArray();
								counts[k][ab-1] += femaleValues.length;
							}
							if (values[k][0][2] != null && (femaleValues == null || values[k][0][2].size() > femaleValues.length)) {
								femaleValues = values[k][0][2].toArray();
							}

							if (femaleValues == null) {
								femaleComp = 0;
								tVals[k][ab-1] = Double.NaN;
							} else {
//								System.out.println(markerName+"\t"+FullSample.AB_PAIRS[ab-1]+"/"+(k==0?"X":"Y")+"\t"+femaleValues.length+" females\t"+Array.mean(femaleValues));
								femaleComp = Array.mean(femaleValues);
								tVals[k][ab-1] = Math.abs(new Ttest(maleValues, femaleValues).getT());											
							}
							zScores[k][ab-1] = (Array.mean(maleValues) - femaleComp)/Array.stdev(maleValues, false);
						} else {
							zScores[k][ab-1] = Double.NaN;
							tVals[k][ab-1] = Double.NaN;
						}
					}
				}
				
//				int[][] interestedPairs = new int[][] {{0,0}, {1,2}, {2,1}};
				int[][] interestedPairs = new int[][] {{0,0}, {1,2}};
				zMean = count = 0;
				zMin = Double.MAX_VALUE;
				line += markerName;
				for (int k = 0; k < interestedPairs.length; k++) {
					if (Double.isNaN(zScores[interestedPairs[k][0]][interestedPairs[k][1]])) {
						line += "\t.";
					} else {
						line += "\t"+zScores[interestedPairs[k][0]][interestedPairs[k][1]];
						zMean += zScores[interestedPairs[k][0]][interestedPairs[k][1]]*counts[interestedPairs[k][0]][interestedPairs[k][1]];
						zMin = Math.min(zMin, zScores[interestedPairs[k][0]][interestedPairs[k][1]]);
						count += counts[interestedPairs[k][0]][interestedPairs[k][1]];
					}
				}
				zMean /= count;

				
				tMean = count = 0;
				tMin = Double.MAX_VALUE;
				for (int k = 0; k < interestedPairs.length; k++) {
					if (Double.isNaN(tVals[interestedPairs[k][0]][interestedPairs[k][1]])) {
						line += "\t.";
					} else {
						line += "\t"+tVals[interestedPairs[k][0]][interestedPairs[k][1]];
						tMean += tVals[interestedPairs[k][0]][interestedPairs[k][1]]*counts[interestedPairs[k][0]][interestedPairs[k][1]];
						tMin = Math.min(tMin, tVals[interestedPairs[k][0]][interestedPairs[k][1]]);
						count += counts[interestedPairs[k][0]][interestedPairs[k][1]];
					}
				}
				tMean /= count;
				
				line += "\t"+zMean+"\t"+tMean+"\t"+zMin+"\t"+tMin;
				
				count = countBs = 0;
				for (int k = 0; k < abGenotypes.length; k++) {
					if (abGenotypes[k] >= 0) {
						countBs += abGenotypes[k];
						count++;
					}
				}
				maf = countBs/(2*count);
				if (maf > 0.5) {
					maf = 1-maf;
				}
				line += "\t"+maf+eol;
				if (line.length() > 25000) {
					writer.print(line);
					writer.flush();
					line = "";
				}
			}
			writer.print(line);
			System.out.println("Finished analyzing "+markerList.length+" in "+ext.getTimeElapsed(time));

			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing results");
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String markersSubset = null;
		String samplesExclude = null;
		String logfile = null;
		Logger log;
		Project proj;
		String filename = Project.DEFAULT_PROJECT;

//		subset = "data/test.txt";
		filename = "C:/workspace/Genvisis/projects/gedi_exome.properties";

		String usage = "\n" + 
				"cnv.qc.MarkerMetrics requires 0-1 arguments\n" + 
				"   (1) project file (i.e. proj="+filename+" (default))\n"+
				"   (2) filename of subset of markers to include / otherwise all markers (i.e. markersubset=" + markersSubset + " (default))\n" + 
				"   (3) filename of subset of samples to exclude / otherwise all markers (i.e. samplesexclude=" + samplesExclude + " (default))\n" + 
//				"   (4) look for separation between males and females (i.e. file=" + markersSubset + " (default))\n" + 
				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("subset=")) {
				markersSubset = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("sampleexclude=")) {
				samplesExclude = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		
		try {
			proj = new Project(filename, false);
			log = new Logger(logfile);
			
//			separationOfSexes(proj, subset, log);
			
			fullQC(proj, samplesExclude, markersSubset, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
