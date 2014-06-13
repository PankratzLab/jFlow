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
	private boolean[] markersForCallrate;

	public LrrSd(Project proj, String[] samples, boolean[] markersForCallrate, String centroidsFile, int threadNumber, int numThreads) {
		this.proj = proj;
		this.samples = samples;
		this.centroidsFile = centroidsFile;
		this.threadNumber = threadNumber;
		this.numThreads = numThreads;
		this.markersForCallrate = markersForCallrate;
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
		int numCallMarkers;

		try {
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+"lrr_sd."+threadNumber));
			writer.println("Sample\tLRR_SD\tBAF1585_SD\tAB_callrate\tForward_callrate");

			if (centroidsFile==null) {
				cents = null;
			} else {
				cents = Centroids.load(centroidsFile, false).getCentroids(); 
			}
			
			for (int i = 0; i<samples.length; i++) {
	        	System.out.println((i+1)+" of "+samples.length);
				fsamp = proj.getFullSampleFromRandomAccessFile(samples[i]);
				chrs = proj.getMarkerSet().getChrs();
				if (fsamp == null) {
					System.err.println("Error - "+samples[i]+Sample.SAMPLE_DATA_FILE_EXTENSION+" not found in samples directory");
				} else {
					lrrs = cents==null?fsamp.getLRRs():fsamp.getLRRs(cents);
					lrrs = Array.subArray(lrrs, 0, Array.indexOfByte(chrs, (byte)23));
					bafs = fsamp.getBAFs();
					bafs = Array.subArray(bafs, 0, Array.indexOfByte(chrs, (byte)23));
					bafsWide = Array.subArray(bafs, 0, Array.indexOfByte(chrs, (byte)23));
					abGenotypes = fsamp.getAB_Genotypes();
					forwardGenotypes = fsamp.getForwardGenotypes();
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
					numCallMarkers = 0;
					if (abGenotypes != null) {
						for (int j = 0; j < abGenotypes.length; j++) {
							if (markersForCallrate == null || markersForCallrate[j]) {
								numCallMarkers++;
								if (abGenotypes[j] >= 0) {
									abCallRate++;
								}
							}
						}
						if (abCallRate > 0) {
							abCallRate /= numCallMarkers;
						}
					}
					forwardCallRate = 0;
					numCallMarkers = 0;
					if (forwardGenotypes != null) {
						for (int j = 0; j < forwardGenotypes.length; j++) {
							if (markersForCallrate == null || markersForCallrate[j]) {
								numCallMarkers++;
								if (forwardGenotypes[j] > 0) {
									forwardCallRate++;
								}
							}
						}
						if (forwardCallRate > 0) {
							forwardCallRate /= numCallMarkers;
						}
					}
					multimodal = Array.isMultimodal(Array.toDoubleArray(Array.removeNaN(bafsWide)), 0.1, 0.5, 0.01);
					writer.println(samples[i]+"\t"+Array.stdev(lrrs, true)+"\t"+Array.stdev(bafs, true)+"\t"+abCallRate+"\t"+forwardCallRate+"\t"+multimodal+"\t"+Array.toStr(bafBinCounts));
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
		Files.cat(files, proj.getProjectDir()+"lrr_sd.xln", Array.intArray(files.length, 0), new Logger(null));
		for (int i = 0; i<files.length; i++) {
			new File(files[i]).delete();
        }
	}
	
	public static boolean[] getMarkerSubset(Project proj, String markersSubsetFile) {
		String[] subMarkers = null;
		if (markersSubsetFile != null) {
			subMarkers = HashVec.loadFileToStringArray(markersSubsetFile, false, new int[] { 0 }, false);
		}
		return getMarkerSubset(proj, subMarkers);
	}

	private static boolean[] getMarkerSubset(Project proj, String[] subMarkers) {
		String[] markers = proj.getMarkerNames();
		boolean[] markersForCallrate = new boolean[markers.length];
		if (subMarkers == null) {
			Arrays.fill(markersForCallrate, true);
		} else {
			Arrays.fill(markersForCallrate, false);
			int[] indicesToUse = ext.indexLargeFactors(subMarkers, markers, true, new Logger(null), true, false);
			for (int i = 0; i < indicesToUse.length; i++) {
				if (indicesToUse[i] < 0) {
					return null;
				} else {
					markersForCallrate[indicesToUse[i]] = true;
				}
			}
		}
		return markersForCallrate;
	}
	
	public static void init(Project proj, String customSampleFileList, String centroidsFile, int numThreads) {
		init(proj, customSampleFileList, null, centroidsFile, numThreads);
	}

	public static void init(Project proj, String customSampleFileList, String markersForCallrateFile, String centroidsFile, int numThreads) {
		String[] samples, subsamples;
		String[][] threadSeeds;
		LrrSd[] runables;
		boolean error;
		boolean[] markersForCallrate;
		error = false;
		samples = proj.getSamples();
		if (customSampleFileList != null) {
			subsamples = HashVec.loadFileToStringArray(customSampleFileList, false, new int[] {0}, false);
			for (int i = 0; i < subsamples.length; i++) {
				if (ext.indexOfStr(subsamples[i], samples) == -1) {
					System.err.println("Error - subsample '"+subsamples[i]+"' was not found in the list of samples of project '"+proj.getNameOfProject()+"'");
					error = true;
				}
			}
			if (error) {
				System.err.println("Error - missing some samples, QC will not be performed");
				return;
			} else {
				samples = subsamples;
			}
		}
		
		markersForCallrate = getMarkerSubset(proj, markersForCallrateFile);
		if (markersForCallrate == null) {
			System.err.println("Error - missing markers from " + markersForCallrateFile + ", QC will not be performed");
			return;
		}
		
		threadSeeds = Parallelizable.splitList(samples, numThreads, false);
		runables = new LrrSd[numThreads];
		for (int i = 0; i<numThreads; i++) {
			runables[i] = new LrrSd(proj, threadSeeds[i], markersForCallrate, centroidsFile, i + 1, numThreads);
        }
		
		Parallelizable.launch(runables);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = Project.DEFAULT_PROJECT;
		String centroids = null;
		String filenameOfListOfSamples = null;
		String markersForCallrateFile = null;
		int numThreads = 1;
		Project proj;

		String usage = "\n"+
		"cnv.qc.LrrSd requires 0-4 arguments\n"+
		"   (1) project file (i.e. proj="+filename+" (default))\n"+
		"   (2) centroids with which to compute LRRs (i.e. cents=genotype.cent (not the default; to be found in data/ directory))\n"+
		"   (3) number of threads to use (i.e. threads="+numThreads+" (default))\n"+
		"   (4) optional: if you only want to look at a subset of the samples, filename of sample list (i.e. subsample=these.txt (not the default))\n"+
		"   (5) optional: if you only want to compute call rates from a subset of the markers, filename of marker list (i.e. callRateMarkers=those.txt (not the default))\n"+
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

			init(proj, filenameOfListOfSamples, markersForCallrateFile, centroids, numThreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
