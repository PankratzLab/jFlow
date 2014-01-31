package cnv.qc;

import java.io.*;
import java.util.*;

import stats.Ttest;
import cnv.filesys.*;
import cnv.manage.MarkerDataLoader;
import cnv.var.SampleData;
import common.*;
import db.FilterDB;

public class MarkerMetrics {
	public static final String[] FULL_QC_HEADER = {"MarkerName", "Chr", "CallRate", "meanTheta_AA", "meanTheta_AB", "meanTheta_BB", "diffTheta_AB-AA", "diffTheta_BB-AB", "sdTheta_AA", "sdTheta_AB", "sdTheta_BB", "meanR_AA", "meanR_AB", "meanR_BB", "num_AA", "num_AB", "num_BB", "pct_AA", "pct_AB", "pct_BB", "MAF", "HetEx"};
	public static final String[] LRR_VARIANCE_HEADER = {"MarkerName", "Chr", "Position", "SD_LRR", "MeanAbsLRR", "SD_BAF1585", "MeanAbsBAF1585"};
	
	public static final String DEFAULT_REVIEW_CRITERIA = "cnv/qc/default_review.criteria";
	public static final String DEFAULT_EXCLUSION_CRITERIA = "cnv/qc/default_exclusion.criteria";
	
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
        MarkerDataLoader markerDataLoader;
        String[] markerNames;
        String line, eol;
        int[] counts;
        double[] sumTheta, sumR, meanTheta, sdTheta;
        double temp;

        if (System.getProperty("os.name").startsWith("Windows")) {
        	eol = "\r\n";
		} else {
			eol = "\n";
		}
        
        samples = proj.getSamples();
        clusterFilterCollection = proj.getClusterFilterCollection();
        gcThreshold = Float.parseFloat(proj.getProperty(Project.GC_THRESHOLD));

		try {
			writer = new PrintWriter(new FileWriter(proj.getFilename(Project.MARKER_METRICS_FILENAME, true, false)));
			writer.println(Array.toStr(FULL_QC_HEADER));
			
			if (markersToInclude != null) {
				markerNames = HashVec.loadFileToStringArray(proj.getDir(Project.PROJECT_DIRECTORY)+markersToInclude, false, new int[] {0}, false);
			} else {
				markerNames = proj.getMarkerNames();
			}
			markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames);
			line = "";
			time = new Date().getTime();
			System.out.println("hi!");
			for (int i = 0; i < markerNames.length; i++) {
				markerData = markerDataLoader.requestMarkerData(i);
				if (i % 1000 == 0) {
					System.out.println(ext.getTime()+"\tMarker "+i+" of "+markerNames.length);
				}

				markerName = markerData.getMarkerName();
				thetas = markerData.getThetas();
				rs = markerData.getRs();
				abGenotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerName, gcThreshold);
				
				counts = new int[4];
				sumTheta = new double[counts.length];
				sumR = new double[counts.length];
				for (int j = 0; j < samples.length; j++) {
					if (samplesToExclude==null || !samplesToExclude[j]) {
						counts[abGenotypes[j] + 1] ++;
						sumTheta[abGenotypes[j] + 1] += thetas[j];
						sumR[abGenotypes[j] + 1] += rs[j];
					}
				}

				meanTheta = new double[counts.length];
				for (int j = 1; j < meanTheta.length; j++) {
					meanTheta[j] = sumTheta[j] / counts[j];
				}
				sdTheta = new double[counts.length];
				for (int j = 0; j < samples.length; j++) {
					if (samplesToExclude==null || !samplesToExclude[j]) {
						temp = (thetas[j] - meanTheta[ abGenotypes[j] + 1 ]);
						sdTheta[ abGenotypes[j] + 1 ] +=  temp * temp;
					}
				}
				for (int j = 1; j < sdTheta.length; j++) {
					if (counts[j] == 0) {
						sdTheta[j] = Double.NaN;	
					} else {
						sdTheta[j] = Math.sqrt(sdTheta[j] / (counts[j]-1));
					}
				}

				line += markerName
						+ "\t" + markerData.getChr()
						+ "\t" + (1- ((float)counts[0] / (counts[0] + counts[1] + counts[2] + counts[3])))
						+ "\t" + meanTheta[1]
						+ "\t" + meanTheta[2]
						+ "\t" + meanTheta[3]
						+ "\t" + (meanTheta[2] - meanTheta[1])
						+ "\t" + (meanTheta[3] - meanTheta[2])
						+ "\t" + sdTheta[1]
						+ "\t" + sdTheta[2]
						+ "\t" + sdTheta[3]
						+ "\t" + (sumR[1] / counts[1])
						+ "\t" + (sumR[2] / counts[2])
						+ "\t" + (sumR[3] / counts[3])
						+ "\t" + counts[1]
						+ "\t" + counts[2]
						+ "\t" + counts[3]
						+ "\t" + ((float) counts[1] / (counts[0] + counts[1] + counts[2] + counts[3]))
						+ "\t" + ((float) counts[2] / (counts[0] + counts[1] + counts[2] + counts[3]))
						+ "\t" + ((float) counts[3] / (counts[0] + counts[1] + counts[2] + counts[3]))
						+ "\t" + (float) (counts[1]<counts[3]? (counts[1] + counts[2]) : (counts[2] + counts[3])) / (counts[0] + counts[1] + 2 * counts[2] + counts[3])
						+ "\t" + AlleleFreq.HetExcess(counts[1], counts[2], counts[3])[0]
						+ eol;
				
				if (line.length() > 25000) {
					writer.print(line);
					writer.flush();
					line = "";
				}
				markerDataLoader.releaseIndex(i);
			}
			writer.print(line);
			writer.close();
			log.report("Finished analyzing "+markerNames.length+" in "+ext.getTimeElapsed(time));
		} catch (Exception e) {
			log.reportError("Error writing marker metrics to "+proj.getFilename(Project.MARKER_METRICS_FILENAME, false, false));
			e.printStackTrace();
		}
	}
	
	public static void lrrVariance(Project proj, boolean[] samplesToInclude, String markersToInclude, Logger log) {
		PrintWriter writer;
		float[] lrrs, bafs;
		boolean[] useBAFs;
		MarkerData markerData;
        String markerName;
        long time;
        MarkerDataLoader markerDataLoader;
        String[] markerNames;
        String line, eol;

        if (System.getProperty("os.name").startsWith("Windows")) {
        	eol = "\r\n";
		} else {
			eol = "\n";
		}
        
		try {
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+"lrrVariance.xln"));
			writer.println(Array.toStr(LRR_VARIANCE_HEADER));
			
			if (markersToInclude != null) {
				markerNames = HashVec.loadFileToStringArray(proj.getDir(Project.PROJECT_DIRECTORY)+markersToInclude, false, new int[] {0}, false);
			} else {
				markerNames = proj.getMarkerNames();
			}
			markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames);
			line = "";
			time = new Date().getTime();
			System.out.println("hi!");
			for (int i = 0; i < markerNames.length; i++) {
				markerData = markerDataLoader.requestMarkerData(i);
				if (i % 1000 == 0) {
					System.out.println(ext.getTime()+"\tMarker "+i+" of "+markerNames.length);
				}

				markerName = markerData.getMarkerName();
				lrrs = markerData.getLRRs();
				bafs = markerData.getBAFs();

				line += markerName+"\t"+markerData.getChr()+"\t"+markerData.getPosition();

				if (lrrs == null) {
					System.err.println("Error - null lrr array for marker "+markerName);
				}
				lrrs = Array.removeNaN(Array.subArray(lrrs, samplesToInclude));
				line += "\t"+Array.stdev(lrrs, true);
				for (int j = 0; j < lrrs.length; j++) {
					lrrs[j] = Math.abs(lrrs[j]);
				}
				line += "\t"+Array.mean(lrrs);
				
				useBAFs = Array.clone(samplesToInclude);
				for (int j = 0; j < bafs.length; j++) {
					if (bafs[j] < 0.15 || bafs[j] > 0.85) {
						useBAFs[j] = false;
					}
				}
				bafs = Array.removeNaN(Array.subArray(bafs, useBAFs));
				line += "\t"+Array.stdev(bafs, true);
				for (int j = 0; j < bafs.length; j++) {
					bafs[j] = Math.abs(bafs[j] - 0.50f);
				}
				line += "\t"+Array.mean(bafs);

				line += eol;
				
				if (line.length() > 25000) {
					writer.print(line);
					writer.flush();
					line = "";
				}
				markerDataLoader.releaseIndex(i);
			}
			writer.print(line);
			writer.close();
			log.report("Finished analyzing "+markerNames.length+" in "+ext.getTimeElapsed(time));
		} catch (Exception e) {
			log.reportError("Error writing results");
			e.printStackTrace();
		}
	}
	
	public static void separationOfSexes(Project proj, String subset, Logger log) {
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
        MarkerDataLoader markerDataLoader;
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
 		
        clusterFilterCollection = proj.getClusterFilterCollection();
        gcThreshold = Float.parseFloat(proj.getProperty(Project.GC_THRESHOLD));

		try {
			writer = new PrintWriter(new FileWriter(proj.getDir(Project.RESULTS_DIRECTORY)+"markerGenderChecks.xln"));
			writer.println("SNP\tX_zAA\tY_zBB\tX_tAA\tY_tBB\tMean_Zs\tMean_Ts\tMin_Z\tMin_T\tMAF");
			
			if (subset != null) {
				markerList = HashVec.loadFileToStringArray(proj.getDir(Project.PROJECT_DIRECTORY)+subset, false, new int[] {0}, false);
			} else {
				markerList = proj.getMarkerNames();
			}
			markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerList);
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

	public static void filterMetrics(Project proj, Logger log) {
		String markerMetricsFilename, reviewCriteriaFilename, exclusionCriteriaFilename;

		markerMetricsFilename = proj.getFilename(Project.MARKER_METRICS_FILENAME, false, false);
		if (!Files.exists(markerMetricsFilename)) {
			log.reportError("Error - marker metrics file not found at "+markerMetricsFilename);
			return;
		}
		
		reviewCriteriaFilename = proj.getFilename(Project.MARKER_REVIEW_CRITERIA_FILENAME, false, false);
		if (Files.exists(reviewCriteriaFilename)) {
			log.report("Using "+reviewCriteriaFilename+" for the review criteria");
		} else {
			log.report("Could not find "+reviewCriteriaFilename+", so generating from default parameters");
			Files.copyFileFromJar(DEFAULT_REVIEW_CRITERIA, reviewCriteriaFilename);
		}

		exclusionCriteriaFilename = proj.getFilename(Project.MARKER_EXCLUSION_CRITERIA_FILENAME, false, false);
		if (Files.exists(exclusionCriteriaFilename)) {
			log.report("Using "+exclusionCriteriaFilename+" for the review criteria");
		} else {
			log.report("Could not find "+reviewCriteriaFilename+", so generating from default parameters");
			Files.copyFileFromJar(DEFAULT_EXCLUSION_CRITERIA, exclusionCriteriaFilename);
		}
		
		FilterDB.filter(markerMetricsFilename, reviewCriteriaFilename, proj.getDir(Project.RESULTS_DIRECTORY)+"markersToReview.out", log);
		FilterDB.filter(markerMetricsFilename, exclusionCriteriaFilename, proj.getDir(Project.RESULTS_DIRECTORY)+"markersToExclude.out", log);
	}
	
	public static void tallyFlaggedReviewedChangedAndDropped(Project proj, Logger log) {
		BufferedReader reader;
		PrintWriter writer, writerMissed;
		String[] line;
		Hashtable<String, Vector<String>> warningHash;
		Hashtable<String, String> flaggedMarkers, droppedMarkers, annotatedMarkers, reclusteredMarkers, allOtherMarkers;
		Vector<String> v;
		String[] filenames;
		String dir;
		String[] header, expectedHeader, warnings, markerNames;
		MarkerDataLoader markerDataLoader;
		MarkerData markerData;
		ClusterFilterCollection clusterFilterCollection;
		AnnotationCollection annotationCollection;
		byte[] genotypes;
		boolean zeroedOut;
		float gcThreshold;
		String[] keys;
		int numAnnotated, numReclustered, numDropped;
		String warning;
		boolean problem;
		String missedOutputFile;
		boolean annotated;
		
		dir = proj.getProjectDir();
		gcThreshold = proj.getFloat(Project.GC_THRESHOLD);

		filenames = new String[] {"results/markersToExclude.out", "results/markersToReview.out"};
		filenames = new String[] {"results/markersToBoth.out"};
		filenames = new String[] {"results/markersToReviewOrFlaggedByQC.out"};
		filenames = new String[] {"results/markersToReview.out", "results/markersFlaggedByQC.out", "results/markersToReviewOrFlaggedByQC.out"};
		expectedHeader = new String[] {"Unit", "Reason flagged", "numFlags"};
		
		problem = false;
		for (int i = 0; i < filenames.length; i++) {
			if (Files.exists(dir+filenames[i])) {
				header = Files.getHeaderOfFile(dir+filenames[i], log);
				if (!ext.checkHeader(header, expectedHeader, false)) {
					problem = true; 
				}
			} else {
				log.reportError("Error - could not find "+dir+filenames[i]);
				problem = true;
			}
		}
		if (problem) {
			System.exit(1);
		}
		
        clusterFilterCollection = proj.getClusterFilterCollection();
        annotationCollection = proj.getAnnotationCollection();

        flaggedMarkers = new Hashtable<String, String>();
		for (int i = 0; i < filenames.length; i++) {
			v = new Vector<String>();
			warningHash = new Hashtable<String, Vector<String>>();
			try {
				reader = Files.getAppropriateReader(dir+filenames[i]);
				if (reader == null) {
					log.reportError("Error - could not find "+dir+filenames[i]);
				}
				header = reader.readLine().trim().split("\t");
				ext.checkHeader(header, expectedHeader, true);
				while (reader.ready()) {
					line = reader.readLine().trim().split("\t");
					if (!flaggedMarkers.containsKey(line[0])) {
						flaggedMarkers.put(line[0], "flaggedButNotAnnotated");
					}
					flaggedMarkers.put(line[0], flaggedMarkers.get(line[0])+";"+line[1]);
					warnings = line[1].split(";");
					for (int j = 0; j < warnings.length; j++) {
						warning = warnings[j].substring(0, warnings[j].indexOf(" (")>0?warnings[j].indexOf(" ("):warnings[j].length());
						if (!warningHash.containsKey(warning)) {
							warningHash.put(warning, new Vector<String>());
						}
						warningHash.get(warning).add(line[0]);
						HashVec.addIfAbsent(line[0], v);
					}
					
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + dir+filenames[i] + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + dir+filenames[i] + "\"");
				System.exit(2);
			}
			
			markerNames = Array.toStringArray(v);
			log.report(markerNames.length+" markers met at least one criterion in "+filenames[i]);
			markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames);
			annotatedMarkers = new Hashtable<String, String>();
			reclusteredMarkers = new Hashtable<String, String>();
			droppedMarkers = new Hashtable<String, String>();
			for (int j = 0; j < markerNames.length; j++) {
				if (j % 1000 == 0) {
					System.out.println((j+1)+" of "+markerNames.length);
				}
				do {
					markerData = markerDataLoader.getMarkerData(j);
					if (markerData == null) {
						System.out.println("waiting...");
						try {
							Thread.sleep(1000);
						} catch (InterruptedException ie) {}
					}
				} while (markerData == null);
				zeroedOut = true;
				genotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerNames[j], gcThreshold);
				for (int k = 0; k < genotypes.length; k++) {
					if (genotypes[k] != -1) {
						zeroedOut = false;
					}
				}
				if (zeroedOut) {
					droppedMarkers.put(markerNames[j], "");
				}
				if (annotationCollection != null && annotationCollection.markerHasAnyAnnotation(markerNames[j])) {
					annotatedMarkers.put(markerNames[j], "");
				}
				if (clusterFilterCollection != null && clusterFilterCollection.getClusterFilters(markerNames[j]) != null) {
					reclusteredMarkers.put(markerNames[j], "");
				}
				markerDataLoader.releaseIndex(j);
			}
			
			keys = HashVec.getKeys(warningHash);
			try {
				writer = new PrintWriter(new FileWriter(dir+filenames[i]+".out"));
				writer.println("WarningCriterion\tnumMarkers\tnumAnnotated\tnumReclustered\tnumDropped");
				for (int j = 0; j < keys.length; j++) {
					v = warningHash.get(keys[j]);
					numReclustered = numAnnotated = numDropped = 0;
					for (int k = 0; k < v.size(); k++) {
						if (reclusteredMarkers.containsKey(v.elementAt(k))) {
							numReclustered++;
						}
						if (annotatedMarkers.containsKey(v.elementAt(k))) {
							numAnnotated++;
						}
						if (droppedMarkers.containsKey(v.elementAt(k))) {
							numDropped++;
						}
					}
					writer.println(keys[j]+"\t"+v.size()+"\t"+numAnnotated+"\t"+numReclustered+"\t"+numDropped);
				}
				numReclustered = numAnnotated = numDropped = 0;
				for (int j = 0; j < markerNames.length; j++) {
					if (reclusteredMarkers.containsKey(markerNames[j])) {
						numReclustered++;
					}
					if (annotatedMarkers.containsKey(markerNames[j])) {
						numAnnotated++;
					}
					if (droppedMarkers.containsKey(markerNames[j])) {
						numDropped++;
					}
				}
				writer.println("Any criteria\t"+markerNames.length+"\t"+numAnnotated+"\t"+numReclustered+"\t"+numDropped);

				missedOutputFile = dir+filenames[i]+"_missed.out";
				writerMissed = new PrintWriter(new FileWriter(missedOutputFile));
				allOtherMarkers = HashVec.loadToHashNull(proj.getMarkerNames());
				for (int j = 0; j < markerNames.length; j++) {
					if (!annotatedMarkers.containsKey(markerNames[j])) {
						writerMissed.println(markerNames[j]+"\t"+flaggedMarkers.get(markerNames[j]));
					}
					allOtherMarkers.remove(markerNames[j]);
				}
				numReclustered = numAnnotated = numDropped = 0;
				markerNames = HashVec.getKeys(allOtherMarkers, false, false);
				writer.print("Everything else\t"+markerNames.length);

				for (int j = 0; j < markerNames.length; j++) {
					if (annotationCollection != null && annotationCollection.markerHasAnyAnnotation(markerNames[j])) {
						numAnnotated++;
					}
					if (clusterFilterCollection != null && clusterFilterCollection.getClusterFilters(markerNames[j]) != null) {
						numReclustered++;
					} else {
						allOtherMarkers.remove(markerNames[j]);
					}
				}
				markerNames = HashVec.getKeys(allOtherMarkers, false, false);
				markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames);
				for (int j = 0; j < markerNames.length; j++) {
					do {
						markerData = markerDataLoader.getMarkerData(j);
						if (markerData == null) {
							System.out.println("waiting...");
							try {
								Thread.sleep(500);
							} catch (InterruptedException ie) {}
						}
					} while (markerData == null);
					zeroedOut = true;
					genotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerNames[j], gcThreshold);
					for (int k = 0; k < genotypes.length; k++) {
						if (genotypes[k] != -1) {
							zeroedOut = false;
						}
					}
					if (annotationCollection != null && annotationCollection.markerHasAnyAnnotation(markerNames[j])) {
						annotated = true;
					} else {
						annotated = false;
					}
					
					if (zeroedOut) {
						if (!annotated) {
							writerMissed.println(markerNames[j]+"\tdroppedButUnannotated");
						}
						numDropped++;
					} else {
						if (!annotated) {
							writerMissed.println(markerNames[j]+"\treclusteredButUnannotated");
						}
					}
					markerDataLoader.releaseIndex(j);
				}
				writer.println("\t"+numAnnotated+"\t"+numReclustered+"\t"+numDropped);
				writer.close();
				writerMissed.close();
				if (new File(missedOutputFile).length() == 0) {
					new File(missedOutputFile).delete();
				}
			} catch (Exception e) {
				System.err.println("Error writing to " + dir+filenames[i]+".out");
				e.printStackTrace();
			}
		}
		
		annotationReport(proj, dir+"results/annotationReport.xln", log);
	}
	
	public static void annotationReport(Project proj, String outputFile, Logger log) {
		PrintWriter writer, writerMissed;
		Hashtable<String, String> reclusteredMarkers, droppedMarkers, allOtherMarkers;
		String[] markerNames, markersWithAnnotation;
		MarkerDataLoader markerDataLoader;
		MarkerData markerData;
		ClusterFilterCollection clusterFilterCollection;
		AnnotationCollection annotationCollection;
		byte[] genotypes;
		boolean zeroedOut;
		float gcThreshold;
		int numReclustered, numDropped;
		char[] annotationKeys;
		String annotation;
		String missedOutputFile;
		
		gcThreshold = proj.getFloat(Project.GC_THRESHOLD);

        clusterFilterCollection = proj.getClusterFilterCollection();
        annotationCollection = proj.getAnnotationCollection();

		markerNames = annotationCollection.getMarkerLists();
		log.report(markerNames.length+" markers have an annotation");
		markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames);
		reclusteredMarkers = new Hashtable<String, String>();
		droppedMarkers = new Hashtable<String, String>();
		for (int j = 0; j < markerNames.length; j++) {
			if (j % 1000 == 0) {
				System.out.println((j+1)+" of "+markerNames.length+" annotated markers");
			}
			do {
				markerData = markerDataLoader.getMarkerData(j);
				if (markerData == null) {
					System.out.println("waiting...");
					try {
						Thread.sleep(500);
					} catch (InterruptedException ie) {}
				}
			} while (markerData == null);
			zeroedOut = true;
			genotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerNames[j], gcThreshold);
			for (int k = 0; k < genotypes.length; k++) {
				if (genotypes[k] != -1) {
					zeroedOut = false;
				}
			}
			if (zeroedOut) {
				droppedMarkers.put(markerNames[j], "");
			}
			if (clusterFilterCollection != null && clusterFilterCollection.getClusterFilters(markerNames[j]) != null) {
				reclusteredMarkers.put(markerNames[j], "");
			}
			markerDataLoader.releaseIndex(j);
		}
		
		annotationKeys = annotationCollection.getKeys();
		try {
			writer = new PrintWriter(new FileWriter(outputFile));
			writer.println("Annotation\tnumMarkers\tnumReclustered\tnumDropped");
			for (int j = 0; j < annotationKeys.length; j++) {
				annotation = annotationCollection.getDescriptionForComment(annotationKeys[j], false, false);
				markersWithAnnotation = annotationCollection.getMarkerLists(annotationKeys[j]);
				numReclustered = numDropped = 0;
				for (int k = 0; k < markersWithAnnotation.length; k++) {
					if (reclusteredMarkers.containsKey(markersWithAnnotation[k])) {
						numReclustered++;
					}
					if (droppedMarkers.containsKey(markersWithAnnotation[k])) {
						numDropped++;
					}
				}
				writer.println(annotation+"\t"+markersWithAnnotation.length+"\t"+numReclustered+"\t"+numDropped);
			}
			writer.println("Any annotation\t"+markerNames.length+"\t"+reclusteredMarkers.size()+"\t"+droppedMarkers.size());
			
			allOtherMarkers = HashVec.loadToHashNull(proj.getMarkerNames());
			for (int j = 0; j < markerNames.length; j++) {
				allOtherMarkers.remove(markerNames[j]);
			}
			numReclustered = numDropped = 0;
			markerNames = HashVec.getKeys(allOtherMarkers, false, false);
			writer.print("Everything else\t"+markerNames.length);

			for (int j = 0; j < markerNames.length; j++) {
				if (clusterFilterCollection != null && clusterFilterCollection.getClusterFilters(markerNames[j]) != null) {
					numReclustered++;
				} else {
					allOtherMarkers.remove(markerNames[j]);
				}
			}
			markerNames = HashVec.getKeys(allOtherMarkers, false, false);
			markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames);
			missedOutputFile = ext.addToRoot(outputFile, "_missed");
			writerMissed = new PrintWriter(new FileWriter(missedOutputFile));
			for (int j = 0; j < markerNames.length; j++) {
				do {
					markerData = markerDataLoader.getMarkerData(j);
					if (markerData == null) {
						System.out.println("waiting...");
						try {
							Thread.sleep(500);
						} catch (InterruptedException ie) {}
					}
				} while (markerData == null);
				zeroedOut = true;
				genotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerNames[j], gcThreshold);
				for (int k = 0; k < genotypes.length; k++) {
					if (genotypes[k] != -1) {
						zeroedOut = false;
					}
				}
				if (zeroedOut) {
					writerMissed.println(markerNames[j]+"\tdropped");
					numDropped++;
				} else {
					writerMissed.println(markerNames[j]+"\treclustered");
				}
				markerDataLoader.releaseIndex(j);
			}
			writer.println("\t"+numReclustered+"\t"+numDropped);
			writer.close();
			writerMissed.close();
			if (new File(missedOutputFile).length() == 0) {
				new File(missedOutputFile).delete();
			}
		} catch (Exception e) {
			System.err.println("Error writing to " + outputFile);
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String markersSubset = null;
		String samples = null;
		String logfile = null;
		Logger log;
		Project proj;
		String filename = Project.DEFAULT_PROJECT;
		boolean sexSeparation = false;
		boolean fullQC = false;
		boolean filter = false;
		boolean lrrVariance = false;
		boolean tally = false;

		String usage = "\n" + 
				"cnv.qc.MarkerMetrics requires 0-1 arguments\n" + 
				"   (1) project file (i.e. proj="+filename+" (default))\n"+
				"   (2) filename of subset of markers to include / otherwise all markers (i.e. markers=" + markersSubset + " (default))\n" + 
				"   (3) filename of subset of samples to include / otherwise all samples (i.e. samples=" + samples + " (default))\n" +
				"  AND\n" + 
				"   (4) look at intensity separation between males and females (i.e. -separation (not the default))\n" + 
				"  OR\n" + 
				"   (4) perform full list of checks on marker quality (i.e. -fullQC (not the default))\n" + 
				"  OR\n" + 
				"   (4) filter markers based on filter criteria (i.e. -filter (not the default))\n" + 
				"  OR\n" + 
				"   (4) check variance of LRR to help determine regions of instability (i.e. -lrrVar (not the default))\n" + 
				"  OR\n" + 
				"   (4) tally the number of reviewed markers that were changed or dropped (i.e. -tally (not the default))\n" + 
				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("markers=")) {
				markersSubset = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("samples=")) {
				samples = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-separation")) {
				sexSeparation = true;
				numArgs--;
			} else if (args[i].startsWith("-fullQC")) {
				fullQC = true;
				numArgs--;
			} else if (args[i].startsWith("-lrrVar")) {
				lrrVariance = true;
				numArgs--;
			} else if (args[i].startsWith("-filter")) {
				filter = true;
				numArgs--;
			} else if (args[i].startsWith("-tally")) {
				tally = true;
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
//			subset = "data/test.txt";
//			filename = "C:/workspace/Genvisis/projects/gedi_exome.properties";
//			filename = "/home/npankrat/projects/GEDI_exomeRAF.properties";
//			filename = "/home/npankrat/projects/exome_chip_win.properties";
//			fullQC = true;
//			filter = true;

//			filename = "C:/workspace/Genvisis/projects/SingaporeReplication.properties";
//			fullQC = true;
//			filter = true;

//			markersSubset = "";
			
//			filename = "/home/npankrat/projects/SingaporeReplication.properties";
			filename = "/home/npankrat/projects/GEDI_exomeRAF.properties";
			tally = true;
			
			proj = new Project(filename, false);
			log = new Logger(logfile);
			
			if (sexSeparation) {
				separationOfSexes(proj, markersSubset, log);
			} 
			if (fullQC) {			
				fullQC(proj, proj.getSamplesToExclude(samples, log), markersSubset, log);
			}
			if (lrrVariance) {
				lrrVariance(proj, proj.getSamplesToInclude(samples, log), markersSubset, log);
			}
			if (filter) {			
				filterMetrics(proj, log);
			}
			if (tally) {
				tallyFlaggedReviewedChangedAndDropped(proj, log);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
