// quantile normalize

package cnv.analysis;

import java.io.*;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import cnv.filesys.*;
import cnv.manage.MarkerDataLoader;
import cnv.qc.SexChecks;
import cnv.var.SampleData;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.Matrix;
import common.ext;
import filesys.Segment;

public class AnalysisFormats implements Runnable {
	public static final String[] PROGRAM_OPTIONS = {"PennCNV", "QuantiSNP"};
	public static final String[] OUTPUT_DIRECTORIES = {"PennCNV/", "QuantiSNP/"};
	public static final int PENN_CNV = 1;
	public static final int QUANTISNP = 2;
	// public static final int EM_ITERATIONS = 10;
	public static final int EM_ITERATIONS = 25;
	private Project proj;
	private String[] samples;
	private int program;
	private HashSet<String> hash;
	private int threadCount;

	public AnalysisFormats(Project proj, String[] samples, int program, HashSet<String> hash, int threadCount) {
		this.proj = proj;
		this.samples = samples;
		this.program = program;
		this.hash = hash;
		this.threadCount = threadCount;
	}

	public void run() {
		switch (program) {
		case PENN_CNV:
			penncnv(proj, samples, hash, null, threadCount);
			break;
		case QUANTISNP:
			quantisnp(proj, samples, hash);
			break;
		default:
			System.err.println("Error - invalid program option: "+program);
			break;
		}

	}

	public static void penncnv(Project proj, final String[] samples, final HashSet<String> markersToWrite, String subDir, int threadCount) {
		final String[] markerNames = proj.getMarkerNames();
		final boolean jar;
		final boolean gzip;
		final String dir;
		final String sampleDir;
		final Logger log = proj.getLog();
		
		dir = proj.PENNCNV_DATA_DIRECTORY.getValue(false, false) + (subDir == null ? "" : subDir);
		sampleDir = proj.SAMPLE_DIRECTORY.getValue(false, true);
		new File(dir).mkdirs();
		jar = proj.JAR_STATUS.getValue();
//		gzip = proj.getBoolean(proj.PENNCNV_GZIP_YESNO);
		gzip = proj.PENNCNV_GZIP_YESNO.getValue();
		
//		int threadCount = Runtime.getRuntime().availableProcessors();

		final ConcurrentLinkedQueue<Integer>[] sampleIndexQueues = new ConcurrentLinkedQueue[threadCount];
		for (int i = 0; i < threadCount; i++) {
			sampleIndexQueues[i] = new ConcurrentLinkedQueue<Integer>();
		}
		for (int i = 0; i < samples.length; i++) {
			sampleIndexQueues[i % threadCount].add(i);	
		}

		ExecutorService computeHub = Executors.newFixedThreadPool(threadCount);
		for (int threadI = 0; threadI < threadCount; threadI++) {
			final int myIndex = threadI;
			final long myStartTime = System.currentTimeMillis();
			computeHub.execute(new Runnable() {
				@Override
				public void run() {
					PrintWriter writer;
					int mySampleCount = 0;
					String sampleName;
					float[] lrrs, bafs;
					byte[] genotypes; 
					Sample mySample;
            		int skippedExports = 0;

					while(!sampleIndexQueues[myIndex].isEmpty()) {
						int sampleIndex = sampleIndexQueues[myIndex].poll();
						sampleName = samples[sampleIndex];
						String exportFileName = dir + sampleName + (gzip ? ".gz" : "");
						if (!Files.exists(exportFileName)) {
							log.report(ext.getTime() + "\tExporting " + (sampleIndex + 1) + " of " + samples.length+"\t"+sampleName);
							if (Files.exists(sampleDir + sampleName + Sample.SAMPLE_DATA_FILE_EXTENSION, jar)) {
								mySample = Sample.loadFromRandomAccessFile(sampleDir + sampleName + Sample.SAMPLE_DATA_FILE_EXTENSION, false, false, true, true, true, jar);
							} else {
								log.reportError("Error - the " + sampleName + Sample.SAMPLE_DATA_FILE_EXTENSION + " is not found.");
								return;
							}
							lrrs = mySample.getLRRs();
							bafs = mySample.getBAFs();
							genotypes = mySample.getAB_Genotypes();
						
    						try {
    							writer = Files.getAppropriateWriter(exportFileName);
    							writer.println("Name\t" + sampleName + ".GType\t" + sampleName + ".Log R Ratio\t" + sampleName + ".B Allele Freq");
    							for (int j = 0; j < markerNames.length; j++) {
    								if (markersToWrite == null || markersToWrite.contains(markerNames[j])) {
    									writer.println(markerNames[j] + "\t" + (genotypes[j] == -1 ? "NC" : Sample.AB_PAIRS[genotypes[j]]) + "\t" + lrrs[j] + "\t" + bafs[j]);
    								}
    							}
    							writer.close();
    						} catch (Exception e) {
    							log.reportError("Error writing PennCNV data for " + sampleName);
    							log.reportException(e);
    						}
						} else {
                    		skippedExports++;
                        }
						
						mySampleCount++;
					}

					log.report("Thread " + myIndex + " processed " + mySampleCount + " samples in " + ext.getTimeElapsed(myStartTime)+(skippedExports>0?"; skipped "+skippedExports+" samples that had been exported previously":""));
				}
			});
		}
		
		computeHub.shutdown();
		try {
			computeHub.awaitTermination(Long.MAX_VALUE, java.util.concurrent.TimeUnit.NANOSECONDS);
		} catch (InterruptedException e) {
			log.report("Sample export was interrupted - exported sample files may not be complete or correct.");
		}
		computeHub = null;

	}
	
	private static Centroids[] computeCentroids(final Project proj, final boolean[] includeList, String[] pfbFiles, String[] centFiles, final boolean shiftPFBForSex, int threadCount) {
		PrintWriter writerM, writerF;
		MarkerSet markerSet;
		String sampleDataFile;
		String[] allMarkers, samples, header;//, markersToUse;
		byte[] markerChrs;
		final boolean[] inclSampAll;
		final boolean[] inclSampFemales;
		final boolean[] inclSampMales;
		final int markerCount = Array.booleanArraySum(includeList);
		int[] sampleSex;
		final float[][][] rawCentroidsFemale;
		final float[][][] rawCentroidsMale;
		Vector<String>[] markerLists;
		final Logger log = proj.getLog();
		ExecutorService computeHub;
		final ConcurrentLinkedQueue<Integer>[] markerIndexQueues;
		final Hashtable<Integer, String[][]> pfbInfo;
		final Hashtable<Integer, Integer>[] fullToTruncMarkerIndices;
		Hashtable<String, Vector<String>> sexData;
		final MarkerDataLoader[] markerDataLoaders;
		SampleData sampleData;
		
		markerSet = proj.getMarkerSet();
		sampleData = proj.getSampleData(0, false);


		inclSampAll = proj.getSamplesToInclude(null);
		if (!sampleData.hasExcludedIndividuals()) {
			log.report("Warning - there is no 'Exclude' column in SampleData.txt; centroids will be determined using all samples.");
		}
		samples = proj.getSamples();
		sampleDataFile = proj.SAMPLE_DATA_FILENAME.getValue(false, false);
		header = Files.getHeaderOfFile(sampleDataFile, proj.getLog());
		int sexInd = -1;
		for (int i = 0; i < header.length; i++) {
			if (("CLASS=" + SexChecks.EST_SEX_HEADER).toUpperCase().equals(header[i].toUpperCase())) {
				sexInd = i;
				break;
			}
		}
		if (sexInd == -1) {
			log.reportError("Error - no estimated sex found in sample data file - please run SexChecks with -check argument to generate the required data");
			return null;
		}
		sexData = HashVec.loadFileToHashVec(sampleDataFile, 0, new int[] { sexInd }, "\t", true, false);
		
		inclSampMales = Array.clone(inclSampAll);
		inclSampFemales = Array.clone(inclSampAll);
		sampleSex = new int[inclSampAll.length];
		for (int i = 0; i < samples.length; i++) {
			int sex = sampleData.getSexForIndividual(samples[i]);
			if (sex == -1) {
				sex = Integer.parseInt(sexData.get(samples[i].toUpperCase()).get(0));
			}
			sampleSex[i] = sex;
			if (sex == 1) {
				inclSampFemales[i] = false;
			} else if (sex == 2) {
				inclSampMales[i] = false;
			} else {
				// Leave these for now, but when computing LRRs and BAFs, will need to be crafty....
			}
		}

		if (threadCount == -1) {
			threadCount = Runtime.getRuntime().availableProcessors();
		}
		markerIndexQueues = new ConcurrentLinkedQueue[threadCount];
		markerLists = new Vector[threadCount];
		fullToTruncMarkerIndices = new Hashtable[threadCount];
		markerDataLoaders = new MarkerDataLoader[threadCount];
		for (int i = 0; i < threadCount; i++) {
			markerLists[i] = new Vector<String>();
			markerIndexQueues[i] = new ConcurrentLinkedQueue<Integer>();
			fullToTruncMarkerIndices[i] = new Hashtable<Integer, Integer>();
		}
		
		allMarkers = markerSet.getMarkerNames();
		if (includeList.length != allMarkers.length) {
			log.reportError("Error - mismatched lists; list of markers to include in centroid computation must be the same size as the marker list for the given project"); 
		}
		markerChrs = markerSet.getChrs();
		
		int qInd = 0;
		for (int i = 0; i < markerChrs.length; i++) {
			if (includeList[i]) {
				markerLists[qInd].add(allMarkers[i]);
				fullToTruncMarkerIndices[qInd].put(i, markerLists[qInd].size() - 1);
				markerIndexQueues[qInd].add(Integer.valueOf(i));
				qInd = (qInd + 1) % threadCount;
			}
		}
		
		for (int i = 0; i < threadCount; i++) {
			markerDataLoaders[i] = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, Array.toStringArray(markerLists[i]));
		}
		
		rawCentroidsMale = new float[allMarkers.length][][];
		rawCentroidsFemale = new float[allMarkers.length][][];

		pfbInfo = new Hashtable<Integer, String[][]>();
		
		log.report("Computing sex-specific centroids for " + markerCount + " sex-specific markers on " + threadCount + " thread(s).");
		
		computeHub = Executors.newFixedThreadPool(threadCount);
		for (int i = 0; i < threadCount; i++) {
			final int myIndex = i;
			final long myStartTime = System.currentTimeMillis();
			computeHub.execute(new Runnable() {
				@Override
				public void run() {
					int myMarkerCount = 0;
					while(!markerIndexQueues[myIndex].isEmpty()) {
						Integer indexInt = markerIndexQueues[myIndex].poll();
						if (indexInt == null) continue;
						int index = indexInt.intValue();
						
						if (!includeList[index]) {
							rawCentroidsMale[index] = new float[][]{{Float.NaN, Float.NaN}, {Float.NaN, Float.NaN}, {Float.NaN, Float.NaN}};
							rawCentroidsFemale[index] = new float[][]{{Float.NaN, Float.NaN}, {Float.NaN, Float.NaN}, {Float.NaN, Float.NaN}};
							continue;
						}
						
						int markerIndex = fullToTruncMarkerIndices[myIndex].get(index);
						MarkerData markerData = markerDataLoaders[myIndex].requestMarkerData(markerIndex);
						CentroidCompute centCompM = new CentroidCompute(markerData, 
													null, 
													inclSampMales, 
													false, // NOT intensity only 
													1, // no filtering
													0,  // no filtering
													null,  // no filtering
													true,  // median, not mean
													proj.getLog());
						
						CentroidCompute centCompF = new CentroidCompute(markerData, 
													null, 
													inclSampFemales, 
													false, // NOT intensity only 
													1, // no filtering
													0,  // no filtering
													null,  // no filtering
													true,  // median, not mean
													proj.getLog());
						
						
						centCompM.computeCentroid(true);
						centCompF.computeCentroid(true);
						
						rawCentroidsMale[index] = centCompM.getCentroid();
						rawCentroidsFemale[index] = centCompF.getCentroid();
						
						float[] bafCnt = new float[]{0, 0};
						float[] bafSum = new float[]{0, 0};
						float[] genCnt = new float[]{0, 0};
						float[] bafM = centCompM.getRecomputedBAF();
						byte[] genM = centCompM.getClustGenotypes();
						float[] bafF = centCompF.getRecomputedBAF();
						byte[] genF = centCompF.getClustGenotypes();
						for (int s = 0; s < inclSampAll.length; s++) {
							if (inclSampMales[s]) {
								if (!Float.isNaN(bafM[s])) {
									bafSum[0] += bafM[s];
									bafCnt[0]++;
									if (genM[s] >= 0) {
										genCnt[0]++;
									}
								}
							}
							if (inclSampFemales[s]) {
								if (!Float.isNaN(bafF[s])) {
									bafSum[1] += bafF[s];
									bafCnt[1]++;
									if (genF[s] >= 0) {
										genCnt[1]++;
									}
								}
							}
						}
						
						pfbInfo.put(index, new String[][]{
								{markerData.getMarkerName(), "" + (shiftPFBForSex ? markerData.getChr() - 22 : markerData.getChr()), "" + markerData.getPosition(), "" + (genCnt[0] > 0 ? (bafSum[0] / bafCnt[0]) : 2)},
								{markerData.getMarkerName(), "" + (shiftPFBForSex ? markerData.getChr() - 22 : markerData.getChr()), "" + markerData.getPosition(), "" + (genCnt[1] > 0 ? (bafSum[1] / bafCnt[1]) : 2)}
						});
						if (markerIndex > 0 && markerIndex % 10000 == 0) {
							log.report(ext.getTime() + "\t...sex centroids computed up to marker " + (markerCount - markerIndex) + " of " + markerCount);
						}
						
						markerDataLoaders[myIndex].releaseIndex(markerIndex);
						centCompM = null;
						centCompF = null;
						
						myMarkerCount++;
					}
					
					System.out.println("Thread " + myIndex + " processed " + myMarkerCount + " markers in " + ext.getTimeElapsed(myStartTime));
				}
			});
		}
		
		computeHub.shutdown();
		try {
			computeHub.awaitTermination(Long.MAX_VALUE, java.util.concurrent.TimeUnit.NANOSECONDS);
		} catch (InterruptedException e) {
			log.report("Centroid computation was interrupted - .pfb and .cent files may not be complete or correct.");
		}
		computeHub = null;
		
		int nullCnt = 0;
		for (int i = 0; i < rawCentroidsFemale.length; i++) {
			if (rawCentroidsFemale[i] == null) {
				nullCnt++;
			}
		}
		System.out.println(nullCnt + " null cent entries");
		
		if (pfbFiles != null) {
			log.report("Writing sex-specific PFB files");
				
			try {
				writerM = new PrintWriter(new FileWriter(pfbFiles[0]));
				writerF = new PrintWriter(new FileWriter(pfbFiles[1]));
				
				writerM.println("Name\tChr\tPosition\tPFB");
				writerF.println("Name\tChr\tPosition\tPFB");
				
				for (int i = 0; i < allMarkers.length; i++) {
					if (!includeList[i]) continue;
					String[][] pfbEntry = pfbInfo.get(Integer.valueOf(i));
					writerM.println(pfbEntry[0][0] + "\t" + pfbEntry[0][1] + "\t" + pfbEntry[0][2] + "\t" + pfbEntry[0][3]);
					writerF.println(pfbEntry[1][0] + "\t" + pfbEntry[1][1] + "\t" + pfbEntry[1][2] + "\t" + pfbEntry[1][3]);
				}
				
				writerM.flush();
				writerF.flush();

				writerM.close();
				writerF.close();
			} catch (IOException e1) {
				log.reportError("Error - problem occured when writing to new sex-specific .pfb files");
				log.reportException(e1);
			}
			
			writerM = null;
			writerF = null;
		}
		
		Centroids[] centroids = new Centroids[2]; 
		centroids[0] = new Centroids(rawCentroidsMale, markerSet.getFingerprint());
		centroids[1] = new Centroids(rawCentroidsFemale, markerSet.getFingerprint());
		
		if (centFiles != null) {
			log.report("Writing sex-specific Centroid files");
			
			centroids[0].serialize(centFiles[0]);
			Centroids.exportToText(proj, centFiles[0], centFiles[0] + ".txt", allMarkers);
			
			centroids[1].serialize(centFiles[1]);
			Centroids.exportToText(proj, centFiles[1], centFiles[1] + ".txt", allMarkers);
		}
		
		return centroids;
	}
	
	public static String[] pennCNVSexHackMultiThreaded(Project proj, String gcModelFile, boolean useExcluded, int threadCount) {
		String sampleDataFile;
		final String sampleDir;
		String sexDir, pennDir, pennData;
		final String maleDir;
		final String femaleDir;
		String malePFBFile, femalePFBFile, newGCFile, centFilePathM, centFilePathF;
		final String[] allMarkers;
		final String[] allSamples;
		String[] header;
		final SampleData sampleData;
		MarkerSet ms;
		byte[] markerChrs;
		final boolean jar;
		final boolean gzip;
		final boolean[] includeMarkersList;
		boolean[] includeSamplesList;
		final Hashtable<String, Vector<String>> sexData;
		Centroids[] centroids;
		
		final Logger log = proj.getLog();
		
		pennDir = proj.getProperty(proj.PENNCNV_RESULTS_DIRECTORY);
		pennData = proj.getProperty(proj.PENNCNV_DATA_DIRECTORY);
		sexDir = pennData + "sexSpecific/";
		
		maleDir = sexDir + "male/";
		femaleDir = sexDir + "female/";
		
		new File(sexDir).mkdirs();
		new File(maleDir).mkdir();
		new File(femaleDir).mkdir();
		
		malePFBFile = sexDir + "males.pfb";
		femalePFBFile = sexDir + "females.pfb";
		newGCFile = sexDir + "sexSpecific.gcModel";
		
		centFilePathM = pennData + "sexSpecific/sexSpecific_Male.cent";
		centFilePathF = pennData + "sexSpecific/sexSpecific_Female.cent";

		ms = proj.getMarkerSet();
		sampleData = proj.getSampleData(0, false);
		
		allMarkers = ms.getMarkerNames();
		markerChrs = ms.getChrs();
		Vector<String> markerList = new Vector<String>();
		final Hashtable<String, Integer> markersToIndices = new Hashtable<String, Integer>();
		includeMarkersList = new boolean[allMarkers.length];
		
		for (int i = 0; i < markerChrs.length; i++) {
			switch(markerChrs[i]) {
				case 23:
				case 24:
				case 25:
				case 26:
					includeMarkersList[i] = true;
					markerList.add(allMarkers[i]);
					markersToIndices.put(allMarkers[i], i);
					break;
				default:
					includeMarkersList[i] = false;
					break;
			}
		}
		
		if (Files.exists(centFilePathM) && Files.exists(centFilePathF)) {
			centroids = new Centroids[] {Centroids.load(centFilePathM, proj.JAR_STATUS.getValue()), Centroids.load(centFilePathM, proj.JAR_STATUS.getValue())};
		} else {
			centroids = computeCentroids(proj, includeMarkersList, new String[]{malePFBFile, femalePFBFile}, new String[]{centFilePathM, centFilePathF}, true, threadCount);
		}
		final float[][][] rawCentroidsMale;
		final float[][][] rawCentroidsFemale;
		rawCentroidsMale = centroids[0].getCentroids();
		rawCentroidsFemale = centroids[1].getCentroids();
		
		log.report("Exporting sex-specific sample data");
		
		sampleDir = proj.SAMPLE_DIRECTORY.getValue(false, true);
		jar = proj.JAR_STATUS.getValue();
//		gzip = proj.getBoolean(proj.PENNCNV_GZIP_YESNO);
		gzip = proj.PENNCNV_GZIP_YESNO.getValue();
		
		if (!useExcluded) {
			includeSamplesList = proj.getSamplesToInclude(null);
			if (!sampleData.hasExcludedIndividuals()) {
				log.report("Warning - there is no 'Exclude' column in SampleData.txt; centroids will be determined using all samples.");
			}
			allSamples = Array.subArray(proj.getSamples(), includeSamplesList);
		} else {
			allSamples = proj.getSamples();
		}
		
		sampleDataFile = proj.SAMPLE_DATA_FILENAME.getValue(false, false);
		header = Files.getHeaderOfFile(sampleDataFile, proj.getLog());
		int sexInd = -1;
		for (int i = 0; i < header.length; i++) {
			if (("CLASS=" + SexChecks.EST_SEX_HEADER).toUpperCase().equals(header[i].toUpperCase())) {
				sexInd = i;
				break;
			}
		}
		if (sexInd == -1) {
			log.reportError("Error - no estimated sex found in sample data file - please run SexChecks with -check argument to generate the required data");
			return null;
		}
		sexData = HashVec.loadFileToHashVec(sampleDataFile, 0, new int[] { sexInd }, "\t", true, false);
		
		if (threadCount == -1) {
			threadCount = Runtime.getRuntime().availableProcessors();
		}

		final ConcurrentLinkedQueue<Integer>[] sampleIndexQueues = new ConcurrentLinkedQueue[threadCount];
		for (int i = 0; i < threadCount; i++) {
			sampleIndexQueues[i] = new ConcurrentLinkedQueue<Integer>();
		}
		for (int i = 0; i < allSamples.length; i++) {
			sampleIndexQueues[i % threadCount].add(i);
		}
		

		ExecutorService computeHub = Executors.newFixedThreadPool(threadCount);
		for (int threadI = 0; threadI < threadCount; threadI++) {
			final int myIndex = threadI;
			final long myStartTime = System.currentTimeMillis();
			computeHub.execute(new Runnable() {
				@Override
				public void run() {
					PrintWriter writer;
					int mySampleCount = 0;
					String sampleName;
					float[] thetas, rs;
					byte[] genotypes; 
					int skippedExports = 0;
					
					while(!sampleIndexQueues[myIndex].isEmpty()) {
						Sample mySample;
						int sampleIndex = sampleIndexQueues[myIndex].poll();
						sampleName = allSamples[sampleIndex];
						int sex = sampleData.getSexForIndividual(sampleName);
						if (sex == -1) {
							sex = Integer.parseInt(sexData.get(sampleName.toUpperCase()).get(0));
						}
						boolean compFemale = SexChecks.KARYOTYPES[sex].contains("XX");
						String exportFileName = (compFemale ? femaleDir : maleDir) + sampleName + (gzip ? ".gz" : "");
						if (!Files.exists(exportFileName)) {
							log.report(ext.getTime() + "\tExporting " + (sampleIndex + 1) + " of " + allSamples.length);
							if (Files.exists(sampleDir + sampleName + Sample.SAMPLE_DATA_FILE_EXTENSION, jar)) {
								mySample = Sample.loadFromRandomAccessFile(sampleDir + sampleName + Sample.SAMPLE_DATA_FILE_EXTENSION, false, true, false, false, true, jar);
							} else {
								log.reportError("Error - the " + sampleName + Sample.SAMPLE_DATA_FILE_EXTENSION + " is not found.");
								// TODO okay to just skip this sample instead of halting entirely?
								continue;
							}
							
							
							thetas = mySample.getThetas();
							rs = mySample.getRs();
							genotypes = mySample.getAB_Genotypes();
						
    						try {
    							writer = Files.getAppropriateWriter(exportFileName);
    							writer.println("Name\t" + sampleName + ".GType\t" + sampleName + ".Log R Ratio\t" + sampleName + ".B Allele Freq");
    							for (int j = 0; j < allMarkers.length; j++) {
    								if (!includeMarkersList[j] || null == (compFemale ? rawCentroidsFemale[j] : rawCentroidsMale[j])) continue;
    								
    								float lrr = Centroids.calcLRR(thetas[j], rs[j], (compFemale ? rawCentroidsFemale[j] : rawCentroidsMale[j]));
    								float baf = Centroids.calcBAF(thetas[j], (compFemale ? rawCentroidsFemale[j] : rawCentroidsMale[j]));
    								
    								writer.println(allMarkers[j] + "\t" + (genotypes[j] == -1 ? "NC" : Sample.AB_PAIRS[genotypes[j]]) + "\t" + lrr + "\t" + baf);
    							}
    							writer.close();
    						} catch (Exception e) {
    							log.reportError("Error writing sex-specific ("+ (compFemale ? "female" : "male") +") PennCNV data for " + sampleName);
    							log.reportException(e);
    						}
						} else {
						    skippedExports++;
						}
						mySampleCount++;
					}

					log.report("Thread " + myIndex + " processed " + mySampleCount + " samples in " + ext.getTimeElapsed(myStartTime)+(skippedExports>0?"; skipped "+skippedExports+" samples that had been exported previously":""));
				}
			});
		}
		
		computeHub.shutdown();
		try {
			computeHub.awaitTermination(Long.MAX_VALUE, java.util.concurrent.TimeUnit.NANOSECONDS);
		} catch (InterruptedException e) {
			log.report("Sample export was interrupted - exported sample files may not be complete or correct.");
		}
		computeHub = null;
		
		filterSexSpecificGCModel(proj, gcModelFile, newGCFile, new String[]{"23", "X", "24", "Y", "25", "XY", "26", "M"});
		
		return new String[]{malePFBFile, femalePFBFile, newGCFile};
	}
	
	public static String[] pennCNVSexHackSingleThreaded(Project proj, String gcModelFile) {
		// exports data for chr23-chr26 and recodes them as chr1-chr4 in a new subdirectory ~/penndata/sexSpecific/
		boolean jar, gzip, writeNewPFBs, writeCentroids, writeGCFile;
		boolean[] inclSampAll, inclSampMales, inclSampFemales;
		int[] sampleSex;
		byte[] markerChrs, genM, genF, genotypes;
		float[] bafCnt, bafSum, genCnt, bafM, bafF, thetas, rs;
		String[] allMarkers, sexMarkers, samples, header, centFilePathM, centFilePathF;
		Logger log;
		MarkerSet ms;
		MarkerDataLoader markerDataLoader;
		SampleData sampleData;
		PrintWriter writer;
		Sample samp;
		String sampleDataFile, sampleDir, sexDir, pennDir, pennData, maleDir, femaleDir, malePFBFile, femalePFBFile, newGCFile;
		float[][][] rawCentroidsMale, rawCentroidsFemale;

		Vector<String[]> malePFBs, femalePFBs;
		Vector<String> markerList;
		Hashtable<String, Vector<String>> sexData;
		Hashtable<String, Integer> sexMarkerToIndex = new Hashtable<String, Integer>();
		
		Centroids centroidsMale, centroidsFemale;
		
		log = proj.getLog();

		pennDir = proj.getProperty(proj.PENNCNV_RESULTS_DIRECTORY);
		pennData = proj.getProperty(proj.PENNCNV_DATA_DIRECTORY);
		sexDir = proj.PROJECT_DIRECTORY.getValue() + pennDir + pennData + "sexSpecific/";
		
		maleDir = sexDir + "male/";
		femaleDir = sexDir + "female/";
		
		new File(sexDir).mkdirs();
		new File(maleDir).mkdir();
		new File(femaleDir).mkdir();
		
		malePFBFile = sexDir + "males.pfb";
		femalePFBFile = sexDir + "females.pfb";
		newGCFile = sexDir + "sexSpecific.gcModel";

//		centFilePathM = sexDir + "sexSpecific_Male.cent";
//		centFilePathF = sexDir + "sexSpecific_Female.cent";
		
		centFilePathM = new String[]{pennDir + pennData + "sexSpecific/sexSpecific_Male.cent", ""};
		centFilePathF = new String[]{pennDir + pennData + "sexSpecific/sexSpecific_Female.cent", ""};
		centFilePathM[1] = proj.PROJECT_DIRECTORY.getValue() + centFilePathM[0];
		centFilePathF[1] = proj.PROJECT_DIRECTORY.getValue() + centFilePathF[0];
		
		
		writeCentroids = !Files.exists(centFilePathM[1]) || !Files.exists(centFilePathF[1]);
		writeNewPFBs = !Files.exists(malePFBFile) && !Files.exists(femalePFBFile);
		writeGCFile = !Files.exists(newGCFile);
		
		if (!writeCentroids && !writeNewPFBs && !writeGCFile) {
			return new String[]{malePFBFile, femalePFBFile, newGCFile};
		}
		
		ms = proj.getMarkerSet();
		sampleData = proj.getSampleData(0, false);
		
		allMarkers = ms.getMarkerNames();
		markerChrs = ms.getChrs();
		markerList = new Vector<String>();
		
		for (int i = 0; i < markerChrs.length; i++) {
			switch(markerChrs[i]) {
				case 23:
				case 24:
				case 25:
				case 26:
					markerList.add(allMarkers[i]);
					sexMarkerToIndex.put(allMarkers[i], Integer.valueOf(i));
					break;
			}
		}
		sexMarkers = Array.toStringArray(markerList);
		
		markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, sexMarkers);
		
		inclSampAll = proj.getSamplesToInclude(null);
		if (!sampleData.hasExcludedIndividuals()) {
			log.report("Warning - there is no 'Exclude' column in SampleData.txt; centroids will be determined using all samples.");
		}
		samples = proj.getSamples();//Array.subArray(proj.getSamples(), inclSampAll);
		sampleDataFile = proj.SAMPLE_DATA_FILENAME.getValue(false, false);
		header = Files.getHeaderOfFile(sampleDataFile, proj.getLog());
		int sexInd = -1;
		for (int i = 0; i < header.length; i++) {
			if (("CLASS=" + SexChecks.EST_SEX_HEADER).toUpperCase().equals(header[i].toUpperCase())) {
				sexInd = i;
				break;
			}
		}
		if (sexInd == -1) {
			log.reportError("Error - no estimated sex found in sample data file - please run SexChecks with -check argument to generate the required data");
			return null;
		}
		sexData = HashVec.loadFileToHashVec(sampleDataFile, 0, new int[] { sexInd }, "\t", true, false);
		
		inclSampMales = new boolean[inclSampAll.length];
		inclSampFemales = new boolean[inclSampAll.length];
		sampleSex = new int[inclSampAll.length];
		for (int i = 0; i < inclSampAll.length; i++) {
			int sex = sampleData.getSexForIndividual(samples[i]);
			if (sex == -1) {
				sex = Integer.parseInt(sexData.get(samples[i].toUpperCase()).get(0));
			}
			sampleSex[i] = sex;
			if (sex == 1) {
				inclSampFemales[i] = false;
			} else if (sex == 2) {
				inclSampMales[i] = false;
			} else {
				// Leave these for now, but when computing LRRs and BAFs, will need to be crafty....
			}
		}
		
		// TODO should we write these out as they're computed instead of storing and writing later?
		malePFBs = new Vector<String[]>();
		femalePFBs = new Vector<String[]>();
		
		rawCentroidsMale = new float[sexMarkers.length][][];
		rawCentroidsFemale = new float[sexMarkers.length][][];
		
		log.report("Computing sex-specific centroids for " + sexMarkers.length + " sex-specific markers on one thread.");
		CentroidCompute centCompM;
		CentroidCompute centCompF;
		for (int i = 0; i < sexMarkers.length; i++) {
			MarkerData markerData = markerDataLoader.requestMarkerData(i);
			
			centCompM = new CentroidCompute(markerData, 
											null, 
											inclSampMales, 
											false, // NOT intensity only 
											1, // no filtering
											0,  // no filtering
											null,  // no filtering
											true,  // median, not mean
											proj.getLog());
			
			centCompF = new CentroidCompute(markerData, 
											null, 
											inclSampFemales, 
											false, // NOT intensity only 
											1, // no filtering
											0,  // no filtering
											null,  // no filtering
											true,  // median, not mean
											proj.getLog());
			
			
			centCompM.computeCentroid(true);
			centCompF.computeCentroid(true);
			
			rawCentroidsMale[i] = centCompM.getCentroid();
			rawCentroidsFemale[i] = centCompF.getCentroid();
			
			bafCnt = new float[]{0, 0};
			bafSum = new float[]{0, 0};
			genCnt = new float[]{0, 0};
			bafM = centCompM.getRecomputedBAF();
			genM = centCompM.getClustGenotypes();
			bafF = centCompF.getRecomputedBAF();
			genF = centCompF.getClustGenotypes();
			for (int s = 0; s < inclSampAll.length; s++) {
				if (inclSampMales[s]) {
					if (!Float.isNaN(bafM[s])) {
						bafSum[0] += bafM[s];
						bafCnt[0]++;
						if (genM[s] >= 0) {
							genCnt[0]++;
						}
					}
				}
				if (inclSampFemales[s]) {
					if (!Float.isNaN(bafF[s])) {
						bafSum[1] += bafF[s];
						bafCnt[1]++;
						if (genF[s] >= 0) {
							genCnt[1]++;
						}
					}
				}
			}
			
			malePFBs.add(new String[]{markerData.getMarkerName(), "" + (markerData.getChr() - 22), "" + markerData.getPosition(), "" + (genCnt[0] > 0 ? (bafSum[0] / bafCnt[0]) : 2)});
			femalePFBs.add(new String[]{markerData.getMarkerName(), "" + (markerData.getChr() - 22), "" + markerData.getPosition(), "" + (genCnt[1] > 0 ? (bafSum[1] / bafCnt[1]) : 2)});
			if (i > 0 && i % 10000 == 0) {
				log.report(ext.getTime() + "\t...sex centroids computed up to marker " + i + " of " + sexMarkers.length);
			}
			
			markerDataLoader.releaseIndex(i);
			centCompM = null;
			centCompF = null;
		}
		
		if (writeNewPFBs) {
			log.report("Writing sex-specific PFB files");
			
			try {
				writer = new PrintWriter(new FileWriter(malePFBFile));
				writer.println("Name\tChr\tPosition\tPFB");
				for (String[] male : malePFBs) {
					writer.println(male[0] + "\t" + male[1] + "\t" + male[2] + "\t" + male[3]);
				}
				writer.close();
			} catch (IOException e1) {
				log.reportError("Error - problem occured when writing to new male-only .pfb file");
				log.reportException(e1);
			}
			
			try {
				writer = new PrintWriter(new FileWriter(femalePFBFile));
				writer.println("Name\tChr\tPosition\tPFB");
				for (String[] female : femalePFBs) {
					writer.println(female[0] + "\t" + female[1] + "\t" + female[2] + "\t" + female[3]);
				}
				writer.close();
			} catch (IOException e1) {
				log.reportError("Error - problem occured when writing to new female-only .pfb file");
				log.reportException(e1);
			}
		}
		malePFBs = null;
		femalePFBs = null;
		
		if (writeCentroids) {
			log.report("Writing sex-specific Centroid files");
			
			centroidsMale = new Centroids(rawCentroidsMale, MarkerSet.fingerprint(sexMarkers));
			centroidsMale.serialize(centFilePathM[1]);
			Centroids.exportToText(proj, centFilePathM[0], centFilePathM[0] + ".txt", sexMarkers);
			
			centroidsFemale = new Centroids(rawCentroidsFemale, MarkerSet.fingerprint(sexMarkers));
			centroidsFemale.serialize(centFilePathF[1]);
			Centroids.exportToText(proj, centFilePathF[0], centFilePathF[0] + ".txt", sexMarkers);
		}
		centroidsMale = null;
		centroidsFemale = null;
		
		log.report("Exporting sex-specific sample data");
		
		sampleDir = proj.SAMPLE_DIRECTORY.getValue(false, true);
		jar = proj.JAR_STATUS.getValue();
//		gzip = proj.getBoolean(proj.PENNCNV_GZIP_YESNO);
		gzip = proj.PENNCNV_GZIP_YESNO.getValue();
		
		int skippedExports = 0;
		for (int i = 0; i < samples.length; i++) {
			int sex = sampleData.getSexForIndividual(samples[i]);
			if (sex == -1) {
				sex = Integer.parseInt(sexData.get(samples[i].toUpperCase()).get(0));
			}
			boolean compFemale = SexChecks.KARYOTYPES[sex].contains("XX");

			String exportFileName = (compFemale ? femaleDir : maleDir) + samples[i] + (gzip ? ".gz" : "");
			if (!Files.exists(exportFileName)) {
				log.report(ext.getTime() + "\tTransforming " + (i + 1) + " of " + samples.length);
				if (Files.exists(sampleDir + samples[i] + Sample.SAMPLE_DATA_FILE_EXTENSION, jar)) {
					samp = Sample.loadFromRandomAccessFile(sampleDir + samples[i] + Sample.SAMPLE_DATA_FILE_EXTENSION, false, true, false, false, true, jar);
				} else {
					log.reportError("Error - the " + samples[i] + Sample.SAMPLE_DATA_FILE_EXTENSION + " is not found.");
					// TODO okay to just skip this sample instead of halting entirely?
					continue;
				}
				
				thetas = samp.getThetas();
				rs = samp.getRs();
				genotypes = samp.getAB_Genotypes();

				try {
    				writer = Files.getAppropriateWriter(exportFileName);
    				writer.println("Name\t" + samples[i] + ".GType\t" + samples[i] + ".Log R Ratio\t" + samples[i] + ".B Allele Freq");
    				for (int j = 0; j < sexMarkers.length; j++) {
    					int markerIndex = sexMarkerToIndex.get(sexMarkers[j]).intValue();
    					
    					float lrr = Centroids.calcLRR(thetas[markerIndex], rs[markerIndex], (compFemale ? rawCentroidsFemale[j] : rawCentroidsMale[j]));
    					float baf = Centroids.calcBAF(thetas[markerIndex], (compFemale ? rawCentroidsFemale[j] : rawCentroidsMale[j]));
    					
    					writer.println(sexMarkers[j] + "\t" + (genotypes[markerIndex] == -1 ? "NC" : Sample.AB_PAIRS[genotypes[markerIndex]]) + "\t" + lrr + "\t" + baf);
    				}
    				writer.close();
    			} catch (Exception e) {
    				log.reportError("Error writing sex-specific ("+ (compFemale ? "female" : "male") +") PennCNV data for " + samples[i]);
    				log.reportException(e);
    			}
			} else {
				skippedExports++;
            }
		}

		log.report(skippedExports>0?"Skipped "+skippedExports+" of "+samples.length+" samples that had been exported previously":"");
		
		if (writeGCFile) {
			filterSexSpecificGCModel(proj, gcModelFile, newGCFile, new String[]{"23", "X", "24", "Y", "25", "XY", "26", "M"});
		}
		
		return new String[]{malePFBFile, femalePFBFile, newGCFile};
	}
	
	private static String filterSexSpecificGCModel(Project proj, String gcModelFile, String newGCFile, String[] chrs) {
		// TODO combine method with filter methods in PennCNV - only difference is changing chr #
		BufferedReader reader = null;
		PrintWriter writer = null;
		
		try {
			reader = new BufferedReader(new FileReader(gcModelFile));
			writer = new PrintWriter(new FileWriter(newGCFile));

			String header;
			String temp;
			String[] line;
			if (reader.ready()) {
				header = reader.readLine();
				writer.println(header);
			}
			while(reader.ready()) {
				temp = reader.readLine();
				line = temp.trim().split("[\\s]+");
				for (String chr : chrs) {
					if (line[1].equals(chr)) {
						
						byte chrVal = 0;
						if ("23".equals(chr) || "X".equals(chr)) {
							chrVal = 23;
						} else if ("24".equals(chr) || "Y".equals(chr)) {
							chrVal = 25;
						} else if ("25".equals(chr) || "XY".equals(chr)) {
							chrVal = 26;
						} else if ("26".equals(chr) || "M".equals(chr)) {
							chrVal = 27;
						}
						
						writer.println(line[0] + "\t" + (chrVal - 22) + "\t" + line[2] + "\t" + line[3]);
					}
				}
			}
		} catch (IOException e) {
			proj.getLog().reportError("Error - transforming sex-specific gcModel failed");
			proj.getLog().reportException(e);
			return gcModelFile;
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					proj.getLog().reportError("Error - couldn't properly close file reader for " + gcModelFile);
					proj.getLog().reportException(e);
				}
				reader = null;
			}
			if (writer != null) {
				writer.close();
				writer = null;
			}
		}
		
		return newGCFile;
	}
	
	public static void quantisnp(Project proj, String[] samples, HashSet<String> hash) {
		PrintWriter writer;
		MarkerSet set;
		String[] markerNames = proj.getMarkerNames();
		Sample samp;
		float[] lrrs, bafs;
		byte[] chrs;
		int[] positions;
		
		set = proj.getMarkerSet();
		chrs = set.getChrs();
		positions = set.getPositions();

		new File(proj.PROJECT_DIRECTORY.getValue()+"quanti_data/").mkdirs();
		for (int i = 0; i<samples.length; i++) {
			System.out.println(ext.getTime()+"\tTransforming "+(i+1)+" of "+samples.length);
			samp = proj.getPartialSampleFromRandomAccessFile(samples[i]);
			set.checkFingerprint(samp);
			lrrs = samp.getLRRs();
			bafs = samp.getBAFs();
			
			try {
				writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()+"quanti_data/"+samples[i]));
				writer.println("Name\tChr\tPosition\t"+samples[i]+".Log R Ratio\t"+samples[i]+".B Allele Freq");
				for (int j = 0; j<markerNames.length; j++) {
					if (hash == null || hash.contains(markerNames[j])) {
						writer.println(markerNames[j]+"\t"+chrs[j]+"\t"+positions[j]+"\t"+lrrs[j]+"\t"+bafs[j]);
					}
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing QuantiSNP data for "+samples[i]);
				e.printStackTrace();
			}
		}
	}
	
	public static void batchQuantiSNP(Project proj, int numBatches) {
		Vector<String[]> v = new Vector<String[]>();
		Hashtable<String,String> genders;
		String[] inputs, outputs;
		String commands, gender;

//		genders = HashVec.loadFileToHashString(proj.getFilename(proj.SAMPLE_DATA_FILENAME), "DNA", new String[] {"CLASS=Gender"}, "");
		genders = HashVec.loadFileToHashString(proj.SAMPLE_DATA_FILENAME.getValue(), "DNA", new String[] {"CLASS=Gender"}, "");
		
		inputs = new File(proj.PROJECT_DIRECTORY.getValue()+"quanti_data/").list(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(".qs");
			}
		});
		
		if (inputs == null) {
			System.err.println("Error - QuantiSNP inputs files have not yet been created");
		}

		outputs = new File(proj.RESULTS_DIRECTORY.getValue(false, true)+"QuantiSNP/").list(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith("_output.out");
			}
		});
		
		if (outputs == null) {
			System.out.println("Found "+inputs.length+" samples; creating output directory for QuantiSNP in "+proj.RESULTS_DIRECTORY.getValue(false, true)+"QuantiSNP/");
			outputs = new String[0];
		} else {
			System.out.println("Found "+inputs.length+" samples, as well as results for "+outputs.length+" that have been done (not necessarily the same ones)");
		}
		
		for (int i = 0; i<inputs.length; i++) {
			if (ext.indexOfStr(ext.rootOf(inputs[i])+"_output.out", outputs)==-1) {
				if (genders.containsKey(ext.rootOf(inputs[i]))) {
					gender = genders.get(ext.rootOf(inputs[i]));
					if (gender.equals("M") || gender.equals("1")) {
						gender = "male";
					} else if (gender.equals("F") || gender.equals("2")) {
						gender = "female";
					} else {
						System.err.println("Error - '"+gender+"' is not a valid gender (expecting M/F or 1/2)");
					}
				} else {
					System.err.println("Error - no gender found for subject '"+ext.rootOf(inputs[i])+"'");
					gender = null;
				}					

				v.add(new String[] {ext.rootOf(inputs[i]), gender});
			}
		}

		System.out.println("Made "+numBatches+" batch files that will take care of the "+v.size()+" files yet to parse");

//		commands = "quantisnp.exe --config ../windows/config.dat  --emiters "+EM_ITERATIONS+" --Lsetting 2000000 --maxcopy 3 --printRS --doGCcorrect --gcdir ../gc/b36/ --output "+OUTPUT_DIRECTORIES[1]+"[%0].out --gender [%1]--input-files ../source/[%0].qs 300\n\n";
		commands = "quantisnp --output "+proj.RESULTS_DIRECTORY.getValue(false, true)+OUTPUT_DIRECTORIES[1]+"[%0].out --gender [%1] --input-files ../source/[%0].qs 300\n\n";
		Files.batchIt("batch", null, numBatches, commands, Matrix.toStringArrays(v));
	}	

	// TODO convert this to an Executor
	public static void launch(Project proj, int program, String markers, int numThreads) {
		Vector<Vector<String>> sampleLists = new Vector<Vector<String>>();
		String[] samples = proj.getSamples();
		Thread[] threads;
		HashSet<String> hash;

		for (int i = 0; i<numThreads; i++) {
			sampleLists.add(new Vector<String>());
		}
		for (int i = 0; i<samples.length; i++) {
			sampleLists.elementAt(i%numThreads).add(samples[i]);
		}
		if (markers == null) {
			hash = null;
		} else {
			hash = HashVec.loadFileToHashSet(proj.PROJECT_DIRECTORY.getValue()+markers, false);
		}
		threads = new Thread[numThreads];
		for (int i = 0; i<numThreads; i++) {
			threads[i] = new Thread(new AnalysisFormats(proj, Array.toStringArray(sampleLists.elementAt(i)), program, hash, 1));
			threads[i].start();
			try {
				Thread.sleep(100L);
			} catch (InterruptedException ex) {}
		}
	}
	
	public static void filter(Project proj, String regions, String list, String outfile) {
		PrintWriter writer;
		MarkerSet markers;
		Segment[] segs;
		String[] markerNames;
		byte[] chrs;
		int[] positions;
		int countFromList, countInRegions, countOverlap;
		HashSet<String> hash;
		
		if (outfile == null) {
			System.err.println("Error - outfile is defined as null; need to provide a filename before results can be filtered");
			return;
		}
		
		if (regions.equals("")) {
			segs = new Segment[0];
		} else {
			segs = Segment.loadUCSCregions(proj.PROJECT_DIRECTORY.getValue()+regions, false);
		}

		if (list.equals("")) {
			hash = new HashSet<String>();
		} else {
			hash = HashVec.loadFileToHashSet(proj.PROJECT_DIRECTORY.getValue()+list, false);
		}

		markers = proj.getMarkerSet();
		markerNames = markers.getMarkerNames();
		chrs = markers.getChrs();
		positions = markers.getPositions();
		
		try {
	        writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()+outfile));
	        countFromList = countInRegions = countOverlap = 0;
			for (int i = 0; i<markerNames.length; i++) {
				if (segs.length > 0 && Segment.overlapsAny(new Segment(chrs[i], positions[i], positions[i]), segs)) {
					countInRegions++;
					if (hash.contains(markerNames[i])) {
						countFromList++;
						countOverlap++;
					}
				} else if (hash.contains(markerNames[i])) {
					countFromList++;
				} else {
					writer.println(markerNames[i]);
				}
	        }
			System.out.println("Started off with "+chrs.length+" markers in the dataset");
			System.out.println("   "+countFromList+" of "+hash.size()+" markers on the list were removed ("+ext.formDeci((countFromList-countOverlap)/(double)chrs.length*100, 2, true)+"% of total)");
			System.out.println("   "+countInRegions+" were found within the list of regions and removed ("+ext.formDeci((countInRegions-countOverlap)/(double)chrs.length*100, 2, true)+"% of total)");
			System.out.println("   "+countOverlap+" overlap in filtering criteria ("+ext.formDeci(countOverlap/(double)chrs.length*100, 2, true)+"% of total)");
			System.out.println("Leaving behind "+(chrs.length-countFromList-countInRegions+countOverlap)+" in final marker list ("+ext.formDeci((chrs.length-countFromList-countInRegions+countOverlap)/(double)chrs.length*100, 2, true)+"% of total)");
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+outfile);
	        e.printStackTrace();
        }		
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		int numThreads = 6;
		int program = PENN_CNV;
		String filterRegions = "";
		String filterList = "";
		String markers = null;
		Project proj;

		String usage = "\n"+
		"filesys.AnalysisFormats requires 0-1 arguments\n"+
		"   (1) project properties filename (i.e. proj="+cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n"+
		"   (2) number of threads to use (i.e. threads="+numThreads+" (default))\n"+
		"   (3) filter markers out within specified regions (i.e. filterRegions=problematicRegions.dat (not the default))\n"+
		"   (4) filter markers out from list (i.e. filterList=drops.dat (not the default))\n"+
		"   (5) input/output file of final list of markers to use (all markers if null) (i.e. markers="+markers+" (default))\n"+
		"   (6) program option (i.e. program="+program+" (default))\n";
		for (int i = 0; i<PROGRAM_OPTIONS.length; i++) {
			usage += "           "+(i+1)+" = "+PROGRAM_OPTIONS[i]+"\n";
		}

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("threads=")) {
				numThreads = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("program=")) {
				program = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("filterRegions=")) {
				filterRegions = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("filterList=")) {
				filterList = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("markers=")) {
				markers = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}

//		filterRegions = "data/problematicRegions.dat";
//		filterList = "data/drops.dat";
//		markers = "finalMarkerList.dat";
		try {
			proj = new Project(filename, false);
			if (!filterRegions.equals("") || !filterList.equals("")) {
				filter(proj, filterRegions, filterList, markers);
			} else {
				launch(proj, program, markers, numThreads);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
