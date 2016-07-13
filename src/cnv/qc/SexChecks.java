// -Xms1024M -Xmx1024M
package cnv.qc;

import java.io.*;
import java.util.*;

import org.apache.commons.math3.stat.inference.TTest;

import cnv.analysis.MosaicismDetect;
import cnv.analysis.MosaicismDetect.MosaicBuilder;
import cnv.filesys.*;
import cnv.manage.MDL;
import cnv.manage.MarkerDataLoader;
import cnv.var.MosaicRegion;
import cnv.var.SampleData;
import common.*;
import filesys.LocusSet;
import filesys.Segment;
import stats.*;

public class SexChecks {
	public static final String EST_SEX_HEADER = "Estimated Sex;1=Male;2=Female;3=Klinefelter;4=UPD Klinefelter;5=Mosaic Klinefelter;6=Triple X;7=Mosaic Triple X;8=Turner;9=Mosaic Turner";
	public static final int[] EST_SEX_MAPPING = {0, 1, 2, 1, 1, 1, 2, 2, 2, 2};
	public static final String[] SEX_HEADER = {"Sample", "FID", "IID", "Sex", EST_SEX_HEADER, "Note", "Check", "Median X R", "Median Y R", "R Ratio Y:X", "Number of X BAFs 10-90%", "Median X LRR", "Median Y LRR"};
	public static final String[] KARYOTYPES = {"", "XY", "XX", "XXY", "XXY", "XXY", "XXX", "XXX", "X", "X"};
	
	private static final float XY_R_RATIO_MIN_SEED_MALE = 0.8f;
	private static final float XY_R_RATIO_MAX_SEED_MALE = 1.2f;
	private static final float XY_R_RATIO_MAX_SEED_FEMALE = 0.2f;
	private static final float NUM_SD_FOR_BAF_10_90_OUTLIERS = 2.0f;
	private static final float NUM_SD_FOR_MALE_X_OUTLIERS = 4.0f;
	private static final float NUM_SD_FOR_FEMALE_X_OUTLIERS = 1.5f;
	private static final float NUM_SD_FOR_FEMALE_X_FULL_ANEUPLOIDY = 4.0f;
	private static final float MAX_SD_FOR_Y_OUTLIERS = 5.0f;
	private static final double SEX_DISCRIMINATING_BASE_P_THRESHOLD = 0.001; // This will be bonferroni corrected for number of markers checked
	private static final double MOSAIC_F_CERTAINTY_THRESHOLD = 0.2;
	private static final double MOSAIC_COVERAGE_CERTAINTY_THRESHOLD = 0.8;
	private static final double MOSAIC_COVERAGE_ABSOLUTE_THRESHOLD = 0.5;

	private Project proj;
	private Logger log;
	private MarkerSet markerSet;
	private String[] sampleNames;
	
	private int[][] indicesByChr;
	
	private boolean[] xKeeps;
	private boolean[] yKeeps;

	private String[] xMarkers;
	private String[] yMarkers;
	
	private boolean[] xUseMarkers;
	private boolean[] yUseMarkers;
	
	private float[] rMedX;
	private float[] rMedY;

	private float[][] lrrsX;
	private float[][] lrrsY;
	
	private float[][] bafsX;
	
	private boolean[] seedMales;
	private boolean[] seedFemales;
	
	private float[] lrrMedX;
	private float[] lrrMedY;
	
	private float[] pctXBaf10_90;
	
	private int[] sexes;
	private boolean[] uncertains;
	private String[] notes;


	
	private SexChecks(Project proj, boolean appendToSampleData, String nonCrossHybridizingMarkersFile) {
		this.proj = proj;
		this.log = proj.getLog();
		
		markerSet = proj.getMarkerSet();
		sampleNames = proj.getSamples();
	    if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }
	    
	    log.report("Finding Sex Chromosome Markers...");
	    generateMarkerLists(nonCrossHybridizingMarkersFile);
	    
		log.report("Loading Sex Chromosome Marker Data...");
		gatherMarkerStats();
		
		log.report("Determining Samples to seed Sex Checks...");
		generateSeedSexLists();
		log.report("Found " + Array.booleanArraySum(seedMales) + " obvious males");
		log.report("Found " + Array.booleanArraySum(seedFemales) + " obvious females");
		log.report("Seeding sex checks using these " + (Array.booleanArraySum(seedMales)
														+ Array.booleanArraySum(seedFemales)) + " samples (out of " + sampleNames.length + " total samples)");
		
		
		log.report("Scanning for markers that express differently by sex");
		
		xUseMarkers = sexDiscriminatingXMarkers();
		yUseMarkers = sexDiscriminatingYMarkers();
		
		log.report("Found " + Array.booleanArraySum(xUseMarkers) + " sex differentiating markers out of " + indicesByChr[23].length + " X chromosome markers");
		log.report("Found " + Array.booleanArraySum(yUseMarkers) + " sex differentiating markers out of " + indicesByChr[24].length + " Y chromosome markers");

		log.report("Calculating median sample LRR for identified X and Y chromosome markers");
		lrrMedX = calcMedianLRRs(lrrsX, Array.booleanArrayToIndices(xUseMarkers));
		lrrMedY = calcMedianLRRs(lrrsY, Array.booleanArrayToIndices(yUseMarkers));
		
		log.report("Calculating sample counts of heterozygote calls for identified X chromosome markers");
		pctXBaf10_90 = calcPctBaf10_90(bafsX, Array.booleanArrayToIndices(xUseMarkers));
		
		log.report("Estimating sex for each sample");
		estimateSexes();
		
		log.report("Writing outputs");
		writeToFile(appendToSampleData);	
	}
	
	private void generateMarkerLists(String nonCrossHybridizingMarkersFile) {
		HashSet<String> nonCrossHybridizingMarkers;
		String[] markerNames = markerSet.getMarkerNames();
		if (nonCrossHybridizingMarkersFile == null) {
			log.reportTimeError("No file of markers that do not cross hybridize was provided, all X and Y chromosome markers will be used to determine sex baselines");
			nonCrossHybridizingMarkers = new HashSet<String>(Arrays.asList(markerNames));
		} else {
			log.report("Using " + nonCrossHybridizingMarkersFile + " to identify markers that do not cross hybridize");
			nonCrossHybridizingMarkers = HashVec.loadFileToHashSet(nonCrossHybridizingMarkersFile, false);
		}
		indicesByChr = markerSet.getIndicesByChr();
		byte[] chrs = markerSet.getChrs();
	    
		xKeeps = new boolean[chrs.length];
		yKeeps = new boolean[chrs.length];
		for (int i = 0; i < chrs.length; i++) {
			switch (chrs[i]) {
				case 23:
					xKeeps[i] = nonCrossHybridizingMarkers.contains(markerNames[i]);
					yKeeps[i] = false;
					break;
				case 24:
					xKeeps[i] = false;
					yKeeps[i] = nonCrossHybridizingMarkers.contains(markerNames[i]);
					break;
				default:
					xKeeps[i] = false;
					yKeeps[i] = false;
					break;
			}
		}
		xMarkers = Array.subArray(markerSet.getMarkerNames(), xKeeps);
		yMarkers = Array.subArray(markerSet.getMarkerNames(), yKeeps);
	}
	
	private void gatherMarkerStats() {
		MDL mdl = new MDL(proj, markerSet, xMarkers, Math.max(proj.NUM_THREADS.getValue() - 1, 1), 100);
		
		float[][] rs;
		
		lrrsX = new float[xMarkers.length][sampleNames.length];
		bafsX = new float[xMarkers.length][sampleNames.length];
		rs = new float[sampleNames.length][xMarkers.length];
		for (int m = 0; mdl.hasNext(); m++) {
			MarkerData markerData = mdl.next();
			float[] xs = markerData.getXs();
			float[] ys = markerData.getYs();
			lrrsX[m] = markerData.getLRRs();
			bafsX[m] = markerData.getBAFs();
			for (int s = 0; s < sampleNames.length; s++) {
				rs[s][m] = Centroids.calcR(xs[s], ys[s]);
				
			}
		}
		mdl.shutdown();
		rMedX = new float[sampleNames.length];
		for (int i = 0; i < sampleNames.length; i++) {
			rMedX[i] = Array.median(Array.removeNonFinites(rs[i]));
		}
		

		mdl = new MDL(proj, markerSet, yMarkers, 1, 100);
		
		lrrsY = new float[yMarkers.length][sampleNames.length];
		rs = new float[sampleNames.length][yMarkers.length];
		for (int m = 0; mdl.hasNext(); m++) {
			MarkerData markerData = mdl.next();
			float[] xs = markerData.getXs();
			float[] ys = markerData.getYs();
			lrrsY[m] = markerData.getLRRs();
			for (int s = 0; s < sampleNames.length; s++) {
				rs[s][m] = Centroids.calcR(xs[s], ys[s]);
			}
		}
		mdl.shutdown();
		rMedY = new float[sampleNames.length];
		for (int i = 0; i < sampleNames.length; i++) {
			rMedY[i] = Array.median(Array.removeNonFinites(rs[i]));
		}
		
		
	}
	
	private void generateSeedSexLists() {
		
		seedMales = new boolean[sampleNames.length];
		seedFemales = new boolean[sampleNames.length];
		
		for (int i = 0; i < sampleNames.length; i++) {
			float rRatio = rMedY[i] / rMedX[i];
			
			if (rRatio > XY_R_RATIO_MIN_SEED_MALE && rRatio < XY_R_RATIO_MAX_SEED_MALE) {
				seedMales[i] = true;
			} else if (rRatio < XY_R_RATIO_MAX_SEED_FEMALE) {
				seedFemales[i] = true;
			}
		}
		
	}
	
	private void estimateSexes() {
		float[] maleMedLRRsX = Array.subArray(lrrMedX, seedMales);
		float[] femaleMedLRRsX = Array.subArray(lrrMedX, seedFemales);
		
		float maleMeanX = Array.mean(maleMedLRRsX, true);
		float maleStdDevX = Array.stdev(maleMedLRRsX, true);
		float femaleMeanX = Array.mean(femaleMedLRRsX, true);
		float femaleStdDevX = Array.stdev(femaleMedLRRsX, true);
		
		float maleMeanBafX10_90 = Array.mean(Array.subArray(pctXBaf10_90, seedMales), true);
		float maleStdDevBafX10_90 = Array.stdev(Array.subArray(pctXBaf10_90, seedMales), true);
		float femaleMeanBafX10_90 = Array.mean(Array.subArray(pctXBaf10_90, seedFemales), true);
		float femaleStdDevBafX10_90 = Array.stdev(Array.subArray(pctXBaf10_90, seedFemales), true);
		
		float[] maleMedLRRsY = Array.subArray(lrrMedY, seedMales);
		float[] femaleMedLRRsY = Array.subArray(lrrMedY, seedFemales);
		
		float maleMeanY = Array.mean(maleMedLRRsY, true);
		float maleStdDevY = Array.stdev(maleMedLRRsY, true);
		float femaleMeanY = Array.mean(femaleMedLRRsY, true);
		float femaleStdDevY = Array.stdev(femaleMedLRRsY, true);
		
		sexes = new int[sampleNames.length];
		uncertains = Array.booleanArray(sampleNames.length, false);
		notes = Array.stringArray(sampleNames.length, "");
		
		boolean[] mosaicismCheckUse = mosaicismUse();
		
		int sdForYOutliers = 0;
		while (sdForYOutliers < MAX_SD_FOR_Y_OUTLIERS) {
			if ((maleMeanY - (sdForYOutliers + 1) * maleStdDevY) > (femaleMeanY + (sdForYOutliers + 1) * femaleStdDevY)) {
				 sdForYOutliers++;
			} else {
				break;
			}
		}
		log.report("Using " + sdForYOutliers + " standard deviations from mean male and female Y LRRs to define sex clusters");
		float maleFloorY = maleMeanY - sdForYOutliers * maleStdDevY;
		log.report("Male mean Y LRR:    " + maleMeanY);
		log.report("Male Std Dev Y LRR: " + maleStdDevY);
		log.report("Male Y LRR Floor:   " + maleFloorY);
		float femaleCeilingY = femaleMeanY + sdForYOutliers * femaleStdDevY;
		log.report("Female mean Y LRR:    " + femaleMeanY);
		log.report("Female SD Y LRR:      " + femaleStdDevY);
		log.report("Female Y LRR Ceiling: " + femaleCeilingY);
		
		for (int i = 0; i < sampleNames.length; i++) {
			boolean male = false;
			boolean female = false;
			if (lrrMedY[i] > maleFloorY) {
				male = true;
			} else if (lrrMedY[i] < femaleCeilingY) {
				female = true;
			} else {
				uncertains[i] = true;
				notes[i] += "Median Y LRR (" + lrrMedY[i] + ") is outside of both male and female acceptance intervals; ";
				if (seedMales[i]) {
					male = true;
				} else if (seedFemales[i]) {
					female = true;
				}
			}
			
			if (male) {
				if (seedFemales[i]) {
					uncertains[i] = true;
					notes[i] += "Ratio of Median X R to Median Y R indicated female; ";
				} else if (!seedMales[i]) {
					notes[i] += "Ratio of Median X R to Median Y R outlier; ";
				}
				if (pctXBaf10_90[i] > (maleMeanBafX10_90 + NUM_SD_FOR_BAF_10_90_OUTLIERS * maleStdDevBafX10_90)) {
					if (checkXMosaicism(i, mosaicismCheckUse)) {
						sexes[i] = 5; // Mosaic Klinefelter
					} else {
						sexes[i] = 3; // Full Klinefelter
					}
				} else if (lrrMedX[i] > (maleMeanX + NUM_SD_FOR_MALE_X_OUTLIERS * maleStdDevX)) {
					sexes[i] = 4; // UPD Klinefelter
				} else {
					sexes[i] = 1; // Male
				}
			} else if (female) {
				if (seedMales[i]) {
					uncertains[i] = true;
					notes[i] += "Ratio of Median X R to Median Y R indicated male; ";
				} else if (!seedFemales[i]) {
					notes[i] += "Ratio of Median X R to Median Y R outlier; ";
				}
				
				if (lrrMedX[i] > (femaleMeanX + NUM_SD_FOR_FEMALE_X_OUTLIERS * femaleStdDevX) && checkXMosaicism(i, mosaicismCheckUse)) {
					if (lrrMedX[i] > (femaleMeanX + NUM_SD_FOR_FEMALE_X_FULL_ANEUPLOIDY * femaleStdDevX)) {
						sexes[i] = 6; // Full Triple X
					} else {
						sexes[i] = 7; // Mosaic Triple X
					}
				} else if (lrrMedX[i] < (femaleMeanX - NUM_SD_FOR_FEMALE_X_OUTLIERS * femaleStdDevX) && checkXMosaicism(i, mosaicismCheckUse)) {
					if (lrrMedX[i] < (femaleMeanX - NUM_SD_FOR_FEMALE_X_FULL_ANEUPLOIDY * femaleStdDevX) && pctXBaf10_90[i] < (femaleMeanBafX10_90 - NUM_SD_FOR_BAF_10_90_OUTLIERS * femaleStdDevBafX10_90)) {
						sexes[i] = 8; // Full Turner
					} else {
						sexes[i] = 9; // Mosaic Turner
					}
				} else {
					sexes[i] = 2; // Female
				}
			} else {
				sexes[i] = 0; // Missing
				uncertains[i] = true;
			}
		}
	}
	
	/*
	 * Generates a boolean array where every marker is true except for X chromosome markers that were excluded.
	 * Autosomal markers are needed for the mosaicism checker but we want to only check mosoacism on "good" X chromosome markers
	 */
	private boolean[] mosaicismUse() {
		boolean[] use = Array.booleanArray(markerSet.getPositions().length, true);
		int[] xIndices = Array.booleanArrayToIndices(xKeeps);
		HashSet<Integer> xInclude = new HashSet<Integer>();
		for (int i = 0; i < xIndices.length; i++) {
			if (xUseMarkers[i]) {
				xInclude.add(xIndices[i]);
			}
		}
		for (int index : indicesByChr[23]) {
			if (!xInclude.contains(index)){
				use[index] = false;
			}
		}
		return use;
	}
	
	private boolean checkXMosaicism(int sample, boolean[] use) {
		float[] bafs = proj.getPartialSampleFromRandomAccessFile(sampleNames[sample], false, false, true, false, false).getBAFs();
		MosaicBuilder mosaicBuilder = new MosaicBuilder();
		mosaicBuilder.use(use);
		MosaicismDetect mosaicismDetect = mosaicBuilder.build(proj, sampleNames[sample], markerSet, Array.toDoubleArray(bafs));
		int xStart = markerSet.getPositions()[indicesByChr[23][0]];
		int xStop = markerSet.getPositions()[indicesByChr[23][indicesByChr[23].length - 1]];
		Segment xSegment = new Segment((byte) 23, xStart, xStop);
		LocusSet<MosaicRegion> xMosaic = mosaicismDetect.callMosaic(xSegment, false);
		int numRegions = xMosaic.getLoci().length;
		if (numRegions == 0) {
			return false;
		}
		String notesAdd = (numRegions == 1 ? "Region" : (numRegions + " regions")) + " of X chromosome mocaicism identified: ";
		double totalCoverage = 0;
		for (MosaicRegion mr : xMosaic.getLoci()) {
			if (mr.getCustomF() < MOSAIC_F_CERTAINTY_THRESHOLD) {
				uncertains[sample] = true;
			}
			double regionCoverage = (double)mr.getSize() / xSegment.getSize();
			totalCoverage += regionCoverage;
			notesAdd += "F=" + mr.getCustomF() + ", " + ext.formDeci(regionCoverage * 100, 4, true) + "% coverage (" + mr.getStart() + " - " + mr.getStop() + "); ";
		}
		notesAdd += "Total Mosaic Coverage: " + ext.formDeci(totalCoverage * 100, 4, true) + "%; ";
		if (totalCoverage < MOSAIC_COVERAGE_CERTAINTY_THRESHOLD) {
			if (totalCoverage > MOSAIC_COVERAGE_ABSOLUTE_THRESHOLD) uncertains[sample] = true;
			else {
				uncertains[sample] = false;
				return false;
			}
		}
		notes[sample] += notesAdd;
		return true;
	}
	
	private void writeToFile (boolean appendToSampleData) {
		
		PrintWriter writer;
		String[] lookup;
		String famIndPair;
		
		SampleData sampleData = proj.getSampleData(0, false);

		try {
			writer = new PrintWriter(new FileWriter(proj.SEXCHECK_RESULTS_FILENAME.getValue(true, false)));
			writer.println(Array.toStr(SEX_HEADER));
			for (int i = 0; i<sampleNames.length; i++) {
			    lookup = sampleData.lookup(sampleNames[i]);
				famIndPair = lookup == null ? null : lookup[1];
				if (famIndPair == null) {
					log.reportError("Error - no data for sample '"+sampleNames[i]+"'");
					writer.print(sampleNames[i]+"\t"+".\t.\t-9");
				} else {
					writer.print(sampleNames[i]+"\t"+famIndPair+"\t"+sampleData.getSexForIndividual(sampleNames[i]));
				}
				writer.println("\t" + sexes[i] + 
							   "\t" + (notes[i].equals("") ? "." : notes[i]) +
							   "\t" + (uncertains[i] ? "1" : "0") + 
							   "\t" + rMedX[i] + 
							   "\t" + rMedY[i] + 
							   "\t" + (rMedY[i] / rMedX[i]) + 
							   "\t" + pctXBaf10_90[i] +
							   "\t" + lrrMedX[i] + 
							   "\t" + lrrMedY[i]);
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + proj.SEXCHECK_RESULTS_FILENAME.getValue());
			log.reportException(e);
		}
		if (appendToSampleData) {
			Hashtable<String, String> linkData = new Hashtable<String, String>();
			for (int i = 0; i < sampleNames.length; i++) {
				linkData.put(sampleNames[i], Integer.toString(sexes[i]));
			}
			if (!sampleData.addData(linkData, "DNA", new String[] {"CLASS=" + EST_SEX_HEADER}, ".", "", log)) {
				log.reportError("Error - failed to write Estimated Sex to sample data file");
			}
		}
	}
	
	private float[] calcMedianLRRs(float[][] lrrs, int[] useMarkers) {
		float[][] lrrsBySample = new float[sampleNames.length][useMarkers.length];
		for (int m = 0; m < useMarkers.length; m++) {
			for (int s = 0; s < sampleNames.length; s++) {
				lrrsBySample[s][m] = lrrs[useMarkers[m]][s];
			}
		}
		
		float[] medianLRRs = new float[sampleNames.length];
		for (int i = 0; i < sampleNames.length; i++) {
			medianLRRs[i] = Array.median(Array.removeNonFinites(lrrsBySample[i]));
		}
		
		return medianLRRs;
	}
	
	private float[] calcPctBaf10_90(float[][] bafs, int[] useMarkers) {
		int[] baf10_90_counts = Array.intArray(sampleNames.length, 0);
		int[] baf_counts = Array.intArray(sampleNames.length, 0);
		for (int m = 0; m < useMarkers.length; m++) {
			for (int s = 0; s < sampleNames.length; s++) {
				baf_counts[s]++;
				if (bafs[useMarkers[m]][s] > 0.1f && bafs[useMarkers[m]][s] < 0.9f) {
					baf10_90_counts[s]++;
				}
			}
		}
		float[] pctBaf19_90 = new float[sampleNames.length];
		for (int s = 0; s < sampleNames.length; s++) {
			if (baf_counts[s] == 0) {
				pctBaf19_90[s] = 0.0f;
			} else {
				pctBaf19_90[s] = (float) baf10_90_counts[s] / baf_counts[s];
			}
		}
		return pctBaf19_90;
	}

	public static int mapEstimatedSexToSex(String estCode) {
	    String[] estCodes = EST_SEX_HEADER.split(";");
	    for (int i = 0; i < estCodes.length; i++) {
	        if (estCodes[i].startsWith(estCode)) {
	            return EST_SEX_MAPPING[i];
	        }
	    }
	    return 0;
	}

	public static void sexCheck(Project proj, boolean appendToSampleData) {
		sexCheck(proj, appendToSampleData, null);
	}
	
	public static void sexCheck(Project proj, boolean appendToSampleData, String nonCrossHybridizingMarkersFile) {
		new SexChecks(proj, appendToSampleData, nonCrossHybridizingMarkersFile);
	}
	
	private boolean[] sexDiscriminatingXMarkers() {
		boolean[] discriminatingMarkers = new boolean[xMarkers.length];
		TTest tTest = new TTest();
		for (int i = 0; i < xMarkers.length; i++) {
			double[] markerLrrs = Array.toDoubleArray(lrrsX[i]);
			double[] maleLrrs = Array.removeNonFinites(Array.subArray(markerLrrs, seedMales));
			double[] femaleLrrs = Array.removeNonFinites(Array.subArray(markerLrrs, seedFemales));
			if (maleLrrs.length < 2 || femaleLrrs.length < 2 || Array.mean(femaleLrrs) <= Array.mean(maleLrrs)) {
				discriminatingMarkers[i] = false;
				continue;
			}
			double pVal = tTest.tTest(maleLrrs, femaleLrrs);
			discriminatingMarkers[i] = pVal < SEX_DISCRIMINATING_BASE_P_THRESHOLD / xMarkers.length;
		}
		return discriminatingMarkers;
	}

	private boolean[] sexDiscriminatingYMarkers() {
		boolean[] discriminatingMarkers = new boolean[yMarkers.length];
		TTest tTest = new TTest();
		for (int i = 0; i < yMarkers.length; i++) {
			double[] markerLrrs = Array.toDoubleArray(lrrsY[i]);
			double[] maleLrrs = Array.removeNonFinites(Array.subArray(markerLrrs, seedMales));
			double[] femaleLrrs = Array.removeNonFinites(Array.subArray(markerLrrs, seedFemales));
			if (maleLrrs.length < 2 || femaleLrrs.length < 2 || Array.mean(maleLrrs) <= Array.mean(femaleLrrs)) {
				discriminatingMarkers[i] = false;
				continue;
			}
			double pVal = tTest.tTest(maleLrrs, femaleLrrs);
			discriminatingMarkers[i] = pVal < SEX_DISCRIMINATING_BASE_P_THRESHOLD / yMarkers.length;
			
		}
		return discriminatingMarkers;
	}
	
	public static void markerByMarker(Project proj) {
		PrintWriter writer;
		String[] samples;
		Vector<double[]> xys, baflrrs; 
		float[] xs, ys, lrrs, bafs;
		Vector<String> intensityDeps;
		LogisticRegression lr;
		String output;
	    SampleData sampleData;
	    int[] sexes;
		MarkerDataLoader markerDataLoader;
		MarkerData markerData;
		String[] markerNames;
		long time;
		Logger log;
	    
		log = proj.getLog();
	    sampleData = proj.getSampleData(2, false);
	    samples = proj.getSamples();
	    sexes = new int[samples.length];
	    for (int i = 0; i < samples.length; i++) {
	    	sexes[i] = sampleData.getSexForIndividual(samples[i]);
		}
		
	
		try {
			writer = new PrintWriter(new FileWriter(proj.RESULTS_DIRECTORY.getValue(false, true)+"markerGenderChecks.xln"));
			writer.println("SNP\tX abs(T)\tY abs(T)\tBAF abs(T)\tLRR abs(T)\tX p\tY p\tXY r2\tBAF p\tLRR p\tBAF/LRR r2");
			
	        time = new Date().getTime();
	        markerNames = proj.getMarkerNames();
			markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames);
	        for (int i = 0; i < markerNames.length; i++) {
	        	markerData = markerDataLoader.requestMarkerData(i);
	        	if (i % 100 == 0) {
	        		log.report(ext.getTime()+"\tMarker "+i+" of "+markerNames.length);
	        	}
	
				output = markerData.getMarkerName();
				
				xs = markerData.getXs();
				ys = markerData.getYs();
				bafs = markerData.getBAFs();
				lrrs = markerData.getLRRs();
				
				intensityDeps = new Vector<String>();
				xys = new Vector<double[]>();
				baflrrs = new Vector<double[]>();
				for (int s = 0; s < samples.length; s++) {
					if (ext.isValidDouble(lrrs[s]+"")) {
						intensityDeps.add(sexes[s]+"");
						xys.add(new double[] {xs[s], ys[s]});
						baflrrs.add(new double[] {bafs[s], lrrs[s]});
					}
				}
				
				
				if (intensityDeps.size()==0) {
					log.reportError("Warning - no data for marker "+markerData.getMarkerName());
					output += "\t.\t.\t.\t.\t.\t.\t.\t.";
				} else {
					output += "\t"+Math.abs(new Ttest(Array.toIntArray(Array.toStringArray(intensityDeps)), Matrix.extractColumn(Matrix.toDoubleArrays(xys), 0)).getPvalue());
					output += "\t"+Math.abs(new Ttest(Array.toIntArray(Array.toStringArray(intensityDeps)), Matrix.extractColumn(Matrix.toDoubleArrays(xys), 1)).getPvalue());
					output += "\t"+Math.abs(new Ttest(Array.toIntArray(Array.toStringArray(intensityDeps)), Matrix.extractColumn(Matrix.toDoubleArrays(baflrrs), 0)).getPvalue());
					output += "\t"+Math.abs(new Ttest(Array.toIntArray(Array.toStringArray(intensityDeps)), Matrix.extractColumn(Matrix.toDoubleArrays(baflrrs), 1)).getPvalue());
				}
	
				lr = null;
				try {
					lr = new LogisticRegression(intensityDeps, xys);
					output += "\t"+lr.getSigs()[1]+"\t"+lr.getSigs()[2]+"\t"+(lr.getRsquare()<0?".":lr.getRsquare());
				} catch (Exception e) {
					output += "\t.\t.\t.";
				}
				try {
					lr = new LogisticRegression(intensityDeps, baflrrs);
					output += "\t"+lr.getSigs()[1]+"\t"+lr.getSigs()[2]+"\t"+(lr.getRsquare()<0?".":lr.getRsquare());
				} catch (Exception e) {
					output += "\t.\t.\t.";
				}
	
				writer.println(output);
				writer.flush();
				markerDataLoader.releaseIndex(i);
			}
	        log.reportError("Finished in " + ext.getTimeElapsed(time));
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing results");
			log.reportException(e);
		}
	
	}

	public static void dropMarkers(String allMarkers, String markersToDrop) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		HashSet<String> hashSet;
	
		hashSet = HashVec.loadFileToHashSet(markersToDrop, false);
	
		try {
			reader = new BufferedReader(new FileReader(allMarkers));
			writer = new PrintWriter(new FileWriter(ext.rootOf(allMarkers)+"_dropped.out"));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (!hashSet.contains(line[0])) {
					writer.println(Array.toStr(line));
				}
			}
			writer.close();
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+allMarkers+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+allMarkers+"\"");
			System.exit(2);
		}
	}

	public static void identifyPseudoautosomalBreakpoints(Project proj) {
			PrintWriter writer;
			String[] samples;
			float[] lrrs;
			MarkerData markerData;
	        SampleData sampleData;
	        int[] sexes;
	        byte[] abGenotypes;
	        String markerName;
	        ClusterFilterCollection clusterFilterCollection;
	        float gcThreshold;
	        long time;
	        DoubleVector[] values; // sex
	        MarkerDataLoader markerDataLoader;
	        String[] markerList;
	        String line, eol;
			MarkerSet markerSet;
	        String[] markerNames;
	        boolean[] sexChrs;
	        byte[] chrs;
	        int[][] genotypeCounts;
	        boolean[] samplesToExclude;
	        Logger log;
	        
	        if (Files.isWindows()) {
	        	eol = "\r\n";
			} else {
				eol = "\n";
			}
	        
	        log = proj.getLog();
	        sampleData = proj.getSampleData(2, false);
	        samplesToExclude = proj.getSamplesToExclude();
	        samples = proj.getSamples();
	        sexes = new int[samples.length];
	        for (int i = 0; i < samples.length; i++) {
	        	sexes[i] = Math.max(0, sampleData.getSexForIndividual(samples[i]));
			}
	        
	        markerSet = proj.getMarkerSet();
	        markerNames = markerSet.getMarkerNames();
	        chrs = markerSet.getChrs();
	        sexChrs = new boolean[chrs.length];
	        for (int i = 0; i < chrs.length; i++) {
	        	sexChrs[i] = chrs[i]>=23;
			}
	        markerList = Array.subArray(markerNames, sexChrs);
	 		
	        clusterFilterCollection = proj.getClusterFilterCollection();
	//        gcThreshold = Float.parseFloat(proj.getProperty(Project.GC_THRESHOLD));
	        gcThreshold = proj.getProperty(proj.GC_THRESHOLD).floatValue();
	
	        try {
				writer = new PrintWriter(new FileWriter(proj.RESULTS_DIRECTORY.getValue(true, false)+"pseudoautosomalSearch.xln"));
				writer.println("SNP\tChr\tPosition\tmLRR_M\tmLRR_F\thet_M\thet_F\tmiss_M\tmiss_F");
				
				markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerList);
				time = new Date().getTime();
				line = "";
				for (int i = 0; i < markerList.length; i++) {
					markerData = markerDataLoader.requestMarkerData(i);
	
					markerName = markerData.getMarkerName();
					lrrs = markerData.getLRRs();
					abGenotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerName, gcThreshold, log);
					
					genotypeCounts = new int[2][4]; // sex, genotype
					values = new DoubleVector[2]; // sex
					values[0] = new DoubleVector();
					values[1] = new DoubleVector();
					for (int s = 0; s < samples.length; s++) {
						if (ext.isValidDouble(lrrs[s]+"") && !samplesToExclude[s]) {
							if (sexes[s] == 1 || sexes[s] == 2) {
								values[sexes[s]-1].add(lrrs[s]);
								genotypeCounts[sexes[s]-1][abGenotypes[s]+1]++;
							}
						}
					}
	
					line += markerName +"\t"+ markerData.getChr() +"\t"+ markerData.getPosition();
					if (values[0].size() > 0) {
						line += "\t"+Array.mean(values[0].toArray());
					} else {
						line += "\t.";
					}
					if (values[1].size() > 0) {
						line += "\t"+Array.mean(values[1].toArray());
					} else {
						line += "\t.";
					}
					if (genotypeCounts[0][1]+genotypeCounts[0][2]+genotypeCounts[0][3] > 0) {
						line += "\t"+(double)genotypeCounts[0][2]/(double)(genotypeCounts[0][1]+genotypeCounts[0][2]+genotypeCounts[0][3]);
					} else {
						line += "\t.";
					}
					if (genotypeCounts[1][1]+genotypeCounts[1][2]+genotypeCounts[1][3] > 0) {
						line += "\t"+(double)genotypeCounts[1][2]/(double)(genotypeCounts[1][1]+genotypeCounts[1][2]+genotypeCounts[1][3]);
					} else {
						line += "\t.";
					}
					line += "\t"+genotypeCounts[0][0];
					line += "\t"+genotypeCounts[1][0];
					line += eol;
					
					if (line.length() > 25000) {
						writer.print(line);
						writer.flush();
						line = "";
					}
					markerDataLoader.releaseIndex(i);
				}
				writer.print(line);
				log.report("Identified pseudo-autosomal breakpoints from "+markerList.length+" markers in "+ext.getTimeElapsed(time));
	
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing results");
				log.reportException(e);
			}
			
		}

	public static void main(String[] args) {
		int numArgs = args.length;
		boolean check = false;
		boolean skipSampleData = false;
		String useMarkers = null;
		String markersToDrop = "data/drops.dat";
		String allMarkers = "data/markerListWithIndices.dat";
		boolean drop = false;
		Project proj;
		String filename = null;
		boolean par = false;

		String usage = "\\n"+
		"qc.SexChecks requires 0-1 arguments\n"+
		"   (1) project properties filename (i.e. proj="+cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n"+
		" AND\n"+
		"   (2) check sex of indiviudals (i.e. -check (not the default))\n"+
		"   (3) skip adding estimated sex to Sample Data (i.e. -skipSampleData (not the default))\n"+
		"   (4) filename of list of markers that do not cross hybridize to use for sex determination (i.e. useMarkers=oneHitWonders.txt (not the default))\n"+
		" OR\n"+
		"   (2) drop markers (i.e. -drop (not the default))\n"+
		"   (3) file with all markers (i.e. all="+allMarkers+" (default file))\n"+
		"   (4) list of bad markers (i.e. drop="+markersToDrop+" (default file))\n"+
		" OR\n"+
		"   (2) check sex chromosomes for pseudoautosomal regions (i.e. -PARcheck (not the default))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-check")) {
				check = true;
				numArgs--;
			} else if (args[i].startsWith("-skipSampleData")) {
				skipSampleData = true;
				numArgs--;
			} else if (args[i].startsWith("useMarkers=")) {
				useMarkers = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-drop")) {
				drop = true ;
				numArgs--;
			} else if (args[i].startsWith("all=")) {
				allMarkers = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("drop=")) {
				markersToDrop = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("-PARcheck")) {
				par = true ;
				numArgs--;
			}
		}

		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}

//		check = true;
//		par = true;
//		filename = "D:/home/npankrat/projects/GEDI_exomeRAF.properties";
		try {
			proj = new Project(filename, false);
			
			if (check) {
				sexCheck(proj, !skipSampleData, useMarkers);
			} else if (par) {
				identifyPseudoautosomalBreakpoints(proj);
			} else if (drop) {
				dropMarkers(allMarkers, markersToDrop);
			} else {
				markerByMarker(proj);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
