package cnv.qc;

import java.io.*;
import java.util.Arrays;

import cnv.filesys.*;
//import java.util.*;
import common.*;

public class LrrSd extends Parallelizable {
	private Project proj;
	private String[] samples;
	private String centroidsFile;
	private int threadNumber;
	private int numThreads;
	private boolean[] markersForCallrate, markersForEverythingElse;

	public LrrSd(Project proj, String[] samples, boolean[] markersForCallrate, boolean[] markersForEverythingElse, String centroidsFile, int threadNumber, int numThreads) {
		this.proj = proj;
		this.samples = samples;
		this.centroidsFile = centroidsFile;
		this.threadNumber = threadNumber;
		this.numThreads = numThreads;
		this.markersForCallrate = markersForCallrate;
		this.markersForEverythingElse = markersForEverythingElse;
	}
	
	public void run() {
		PrintWriter writer;
		Sample fsamp;
		float[][][] cents;
		byte[] chrs, abGenotypes, forwardGenotypes;
		float[] lrrs, bafs, bafsWide;
		double abCallRate, forwardCallRate;
		int[] bafBinCounts;
		boolean multimodal;
		int subIndex = -1;
		Logger log;
		
		log = proj.getLog();
		try {
			if (centroidsFile==null) {
				cents = null;
			} else {
				cents = Centroids.load(centroidsFile, false).getCentroids(); 
			}
			chrs = proj.getMarkerSet().getChrs();
			subIndex = Array.indexOfFirstMaxByte(chrs, (byte) 23);// index is the first byte >= 23, chrs.length if all are less, -1 if none are less, 0 if all are greater!
			if (subIndex <= 0) {
				proj.getLog().reportError("Error - was not able to detect any autosomal markers for sample QC in " + proj.getFilename(Project.MARKERSET_FILENAME));
				return;
			}
			if (chrs[subIndex] != 23) {
				proj.getLog().report("Info - did not detect chromosome 23 in " + proj.getFilename(Project.MARKERSET_FILENAME));
			}
			if (markersForEverythingElse != null) {
				markersForEverythingElse = (subIndex == markersForEverythingElse.length ? markersForEverythingElse : Array.subArray(markersForEverythingElse, 0, subIndex));
			}
			
			
			int numAb = (markersForCallrate == null ? chrs.length : Array.booleanArraySum(markersForCallrate));
			int numAllElse = (markersForEverythingElse == null ? subIndex : Array.booleanArraySum(markersForEverythingElse));
			if (threadNumber == 1) {// we can just show this once
				proj.getLog().report("Info - using " + numAb + " markers for sample call rate qc");
				proj.getLog().report("Info - using " + numAllElse + " autosomal markers for all other sample qc metrics");
			}
			if (numAb == 0 || numAllElse == 0) {
				if (numAb == 0) {
					proj.getLog().report("Error - cannot compute sample call rate with 0 markers, halting");
				}
				if (numAllElse == 0) {
					proj.getLog().report("Error - cannot compute sample qc metrics with 0 markers, halting");
				}
				return;
			}
			if (numAb < 1000) {
				proj.getLog().report("Warning - using " + numAb + (numAb == 1 ? " marker" : " markers") + " for sample call rate may result in inaccurate sample qc, please consider using more");
			}
			if (numAllElse < 1000) {
				proj.getLog().report("Warning - using " + numAllElse + (numAllElse == 1 ? " marker" : " markers") + " for other qc metrics may result in inaccurate sample qc, please consider using more");
			}

			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+"lrr_sd."+threadNumber));
			writer.println("Sample\tLRR_AVG\tLRR_SD\tBAF1585_SD\tAB_callrate\tForward_callrate");

			
			for (int i = 0; i<samples.length; i++) {
	        	log.report((i+1)+" of "+samples.length);
				fsamp = proj.getFullSampleFromRandomAccessFile(samples[i]);
				if (fsamp == null) {
					log.reportError("Error - "+samples[i]+Sample.SAMPLE_DATA_FILE_EXTENSION+" not found in samples directory");
				} else {
					lrrs = cents == null ? fsamp.getLRRs() : fsamp.getLRRs(cents);
					lrrs = (subIndex == lrrs.length ? lrrs : Array.subArray(lrrs, 0, subIndex));
					bafs = fsamp.getBAFs();
					bafs = (subIndex == bafs.length ? bafs : Array.subArray(bafs, 0, subIndex));
					bafsWide = bafs;

					if (markersForEverythingElse != null) {
						lrrs = Array.subArray(lrrs, markersForEverythingElse);
						bafs = Array.subArray(bafs, markersForEverythingElse);
						bafsWide = Array.subArray(bafsWide, markersForEverythingElse);
					}
					
					abGenotypes = fsamp.getAB_Genotypes();
					forwardGenotypes = fsamp.getForwardGenotypes();
					
					if (markersForCallrate != null) {//we do not need autosomal only markers here...
						abGenotypes = (abGenotypes == null ? abGenotypes : Array.subArray(abGenotypes, markersForCallrate));
						forwardGenotypes = (forwardGenotypes == null ? forwardGenotypes : Array.subArray(forwardGenotypes, markersForCallrate));
					}

					bafBinCounts = new int[101];
					for (int j = 0; j < bafs.length; j++) {
						if (!Float.isNaN(bafs[j])) {
							bafBinCounts[(int)Math.floor(bafs[j]*100)]++;
						}
						if (bafs[j] < 0.15 || bafs[j] > 0.85) {
							bafs[j] = Float.NaN;
						}
						if (bafsWide[j] < 0.03 || bafsWide[j] > 0.97) {
							bafsWide[j] = Float.NaN;
						}
					}
					abCallRate = 0;
					if (abGenotypes != null) {
						for (int j = 0; j < abGenotypes.length; j++) {
							if (abGenotypes[j] >= 0) {
								abCallRate++;
							}
						}
						abCallRate /= abGenotypes.length;
						
					}
					forwardCallRate = 0;
					if (forwardGenotypes != null) {
						for (int j = 0; j < forwardGenotypes.length; j++) {
							if (forwardGenotypes[j] > 0) {
								forwardCallRate++;
							}
						}
						forwardCallRate /= forwardGenotypes.length;
					}
					multimodal = Array.isMultimodal(Array.toDoubleArray(Array.removeNaN(bafsWide)), 0.1, 0.5, 0.01);
					writer.println(samples[i] + "\t" + Array.mean(lrrs, true) + "\t" + Array.stdev(lrrs, true) + "\t" + Array.stdev(bafs, true) + "\t" + abCallRate + "\t" + forwardCallRate + "\t" + multimodal + "\t" + Array.toStr(bafBinCounts));
					writer.flush();
				}
			}
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void finalAction() {
		String[] files; 
		
		files = Array.stringArraySequence(numThreads, proj.getProjectDir()+"lrr_sd.");
		Files.cat(files, proj.getProjectDir()+"lrr_sd.xln", Array.intArray(files.length, 0), proj.getLog());
		for (int i = 0; i<files.length; i++) {
			new File(files[i]).delete();
        }
	}
	
	private static boolean[] getMarkerSubset(Project proj, String[] subMarkers) {
		String[] markers = proj.getMarkerNames();
		boolean[] markerSubset = new boolean[markers.length];
		if (subMarkers == null) {
			Arrays.fill(markerSubset, true);
		} else {
			Arrays.fill(markerSubset, false);
			int[] indicesToUse = ext.indexLargeFactors(subMarkers, markers, true, proj.getLog(), true, false);
			for (int i = 0; i < indicesToUse.length; i++) {
				if (indicesToUse[i] < 0) {
					return null;
				} else {
					markerSubset[indicesToUse[i]] = true;
				}
			}
		}
		return markerSubset;
	}
	
	public static void init(Project proj, String customSampleFileList, String centroidsFile, int numThreads) {
		init(proj, customSampleFileList, null, null, centroidsFile, numThreads);
	}

	public static void init(Project proj, String customSampleFileList, String markersForCallrateFile, String markersForEverythingElseFile, String centroidsFile, int numThreads) {
		String[] samples, subsamples;
		String[][] threadSeeds;
		LrrSd[] runables;
		boolean error;
		boolean[] markersForCallrate, markersForEverythingElse;
		Logger log;
		
		error = false;
		log = proj.getLog();
		samples = proj.getSamples();
		if (customSampleFileList != null) {
			subsamples = HashVec.loadFileToStringArray(customSampleFileList, false, new int[] {0}, false);
			for (int i = 0; i < subsamples.length; i++) {
				if (ext.indexOfStr(subsamples[i], samples) == -1) {
					log.reportError("Error - subsample '"+subsamples[i]+"' was not found in the list of samples of project '"+proj.getNameOfProject()+"'");
					error = true;
				}
			}
			if (error) {
				log.reportError("Error - missing some samples, QC will not be performed");
				return;
			} else {
				samples = subsamples;
			}
		}
		
		markersForCallrate = null;
		markersForEverythingElse = null;
		if (markersForCallrateFile != null) {
			markersForCallrate = getMarkerSubset(proj, HashVec.loadFileToStringArray(markersForCallrateFile, false, new int[] { 0 }, false));
			if (markersForCallrate == null) {
				log.reportError("Error - Some markers listed in " + markersForCallrateFile + " were not found in the current project, or were duplicates");
				return;
			}
		}
		if (markersForEverythingElseFile != null) {
			markersForEverythingElse = getMarkerSubset(proj, HashVec.loadFileToStringArray(markersForEverythingElseFile, false, new int[] { 0 }, false));
			if (markersForCallrate == null) {
				log.reportError("Error - Some markers listed in " + markersForEverythingElseFile + " were not found in the current project, or were duplicates");
				return;
			}
		}

		threadSeeds = Parallelizable.splitList(samples, numThreads, false);
		runables = new LrrSd[numThreads];
		for (int i = 0; i<numThreads; i++) {
			runables[i] = new LrrSd(proj, threadSeeds[i], markersForCallrate, markersForEverythingElse, centroidsFile, i + 1, numThreads);
        }
		
		Parallelizable.launch(runables, log);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String centroids = null;
		String filenameOfListOfSamples = null;
		String markersForCallrateFile = null;
		String markersForEverythingElseFile = null;
		int numThreads = 1;
		Project proj;

		String usage = "\n"+
		"cnv.qc.LrrSd requires 0-6 arguments\n"+
		"   (1) project properties filename (i.e. proj="+cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n"+
		"   (2) centroids with which to compute LRRs (i.e. cents=genotype.cent (not the default; to be found in data/ directory))\n"+
		"   (3) number of threads to use (i.e. threads="+numThreads+" (default))\n"+
		"   (4) optional: if you only want to look at a subset of the samples, filename of sample list (i.e. subsample=these.txt (not the default))\n"+
		"   (5) optional: if you only want to compute AB_callrate and Forward_callrate from a subset of the markers, filename of marker list (i.e. callRateMarkers=those.txt (not the default))\n"+
		"   (6) optional: if you only want to compute the other qc metrics (excluding AB_callrate and Forward_callrate) from a subset of the markers, filename of marker list (i.e. otherMarkers=this.txt (not the default))\n"+

		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("subsample=")) {
				filenameOfListOfSamples = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("callRateMarkers=")) {
				markersForCallrateFile = args[i].split("=")[1];
				numArgs--;
			}else if (args[i].startsWith("otherMarkers=")) {
				markersForEverythingElseFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("cents=")) {
				centroids = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("threads=")) {
				numThreads = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
//			filename = "/home/npankrat/projects/GEDI.properties";
//			filenameOfListOfSamples = "D:/data/GEDI/plate51.txt";			
//			
			proj = new Project(filename, false);
			init(proj, filenameOfListOfSamples, markersForCallrateFile, markersForEverythingElseFile, centroids, numThreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
