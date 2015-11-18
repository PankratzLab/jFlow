package cnv.qc;

import java.io.*;
import java.util.Arrays;

import cnv.filesys.*;
import cnv.hmm.CNVCaller;
import cnv.qc.GcAdjustor.GC_CORRECTION_METHOD;
import cnv.qc.GcAdjustor.GcModel;
//import java.util.*;
import common.*;

public class LrrSd extends Parallelizable {
	private static final String BOUND_SD = "LRR_SD_" + CNVCaller.MIN_LRR_MEDIAN_ADJUST + "_" + CNVCaller.MAX_LRR_MEDIAN_ADJUST;
	private static final String BOUND_SD_CORRECTED = "LRR_SD_Post_Correction_" + CNVCaller.MIN_LRR_MEDIAN_ADJUST + "_" + CNVCaller.MAX_LRR_MEDIAN_ADJUST;

	public static final String[] NUMERIC_COLUMNS = { "LRR_AVG", "LRR_SD", BOUND_SD, "BAF1585_SD", "Genotype_callrate", "Genotype_heterozygosity", "WF_Prior_Correction", "GCWF_Prior_Correction", "WF_Post_Correction", "GCWF_Post_Correction", "LRR_SD_Post_Correction", BOUND_SD_CORRECTED };
	public static final String SAMPLE_COLUMN = "Sample";
	private Project proj;
	private String[] samples;
	private String centroidsFile;
	private int threadNumber;
	private int numThreads;
	private boolean[] markersForCallrate, markersForEverythingElse;
	private GcModel gcModel;

	public LrrSd(Project proj, String[] samples, boolean[] markersForCallrate, boolean[] markersForEverythingElse, String centroidsFile, GcModel gcModel, int threadNumber, int numThreads) {
		this.proj = proj;
		this.samples = samples;
		this.centroidsFile = centroidsFile;
		this.threadNumber = threadNumber;
		this.numThreads = numThreads;
		this.markersForCallrate = markersForCallrate;
		this.markersForEverythingElse = markersForEverythingElse;
		this.gcModel = gcModel;
	}
	
	public void run() {
		PrintWriter writer;
		Sample fsamp;
        float[][][] cents;
        byte[] chrs/*, abGenotypes, forwardGenotypes*/;
//        float[] lrrs, bafs, bafsWide;
//        double abCallRate, forwardCallRate, abHetRate, forwardHetRate, wfPrior, gcwfPrior, wfPost, gcwfPost, lrrsdBound, lrrsdPost, lrrsdPostBound;
//        int[] bafBinCounts;
//        boolean multimodal;
        int subIndex = -1;
        Logger log;
        
		String PROG_KEY = "LRRSTDEV_" + threadNumber;
		String progDesc = "Compute Log-R Ratio Std.Dev. in Thread " + threadNumber;
		
		ProgressMonitor progMon = proj.getProgressMonitor();
		if (progMon != null) {
			System.out.println(progMon);
			progMon.beginDeterminateTask(PROG_KEY, progDesc, samples.length + 1, ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
		}
		
		log = proj.getLog();
		try {
			if (centroidsFile==null) {
				cents = null;
			} else {
				cents = Centroids.load(centroidsFile, false).getCentroids(); 
			}
			
			if (progMon != null) {
				proj.getProgressMonitor().updateTask(PROG_KEY);
			}
			
			chrs = proj.getMarkerSet().getChrs();
			subIndex = Array.indexOfFirstMaxByte(chrs, (byte) 23);// index is the first byte >= 23, chrs.length if all are less, -1 if none are less, 0 if all are greater!
			if (subIndex <= 0) {
//				proj.getLog().reportError("Error - was not able to detect any autosomal markers for sample QC in " + proj.getFilename(proj.MARKERSET_FILENAME));
				proj.getLog().reportError("Error - was not able to detect any autosomal markers for sample QC in " + proj.MARKERSET_FILENAME.getValue());
				return;
			}
			if (chrs[subIndex] != 23) {
//				proj.getLog().report("Info - did not detect chromosome 23 in " + proj.getFilename(proj.MARKERSET_FILENAME));
				proj.getLog().report("Info - did not detect chromosome 23 in " + proj.MARKERSET_FILENAME.getValue());
			}
			if (markersForEverythingElse != null) {
				for (int i = subIndex; i < markersForEverythingElse.length; i++) {
					markersForEverythingElse[i] = false;
				}
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
			
//			writer = new PrintWriter(new FileWriter(ext.rootOf(proj.getFilename(proj.SAMPLE_QC_FILENAME), false) + "." + threadNumber));
			writer = new PrintWriter(new FileWriter(ext.rootOf(proj.SAMPLE_QC_FILENAME.getValue(), false) + "." + threadNumber));
			writer.println(SAMPLE_COLUMN + "\t" + Array.toStr(NUMERIC_COLUMNS));
			
			for (int i = 0; i<samples.length; i++) {
	        	log.report((i+1)+" of "+samples.length);
				fsamp = proj.getFullSampleFromRandomAccessFile(samples[i]);
				if (fsamp == null) {
					log.reportError("Error - "+samples[i]+Sample.SAMPLE_DATA_FILE_EXTENSION+" not found in samples directory");
				} else {
					writer.println(Array.toStr(LrrSdPerSample(proj, samples[i], fsamp, cents, markersForCallrate, markersForEverythingElse, gcModel, GC_CORRECTION_METHOD.GENVISIS_GC, log), "\t"));
		            writer.flush();
				}
				if (progMon != null) {
					proj.getProgressMonitor().updateTask(PROG_KEY);
				}
			}
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		if (progMon != null) {
			proj.getProgressMonitor().endTask(PROG_KEY);
		}
	}
	
	/**
	 * Returns:<br/>
                sampleID,<br/> 
                Array.mean(lrrs, true)<br/> 
                Array.stdev(lrrs, true)<br/>
                lrrsdBound<br/>
                Array.stdev(bafs, true)<br/> 
                (abCallRate > 0 ? abCallRate : forwardCallRate)<br/>
                (abCallRate > 0 ? abHetRate : forwardHetRate)<br/>
                wfPrior<br/>
                gcwfPrior<br/>
                wfPost<br/>
                gcwfPost <br/>
                lrrsdPost<br/>
                lrrsdPostBound<br/>
                multimodal <br/>
                Array.toStr(bafBinCounts)<br/>
	 * 
	 * 
	 * @param proj
	 * @param sampleID
	 * @param fsamp
	 * @param cents
	 * @param markersForCallrate
	 * @param markersForEverythingElse
	 * @param gcModel
	 * @param log
	 * @return
	 */
	public static String[] LrrSdPerSample(Project proj, String sampleID, Sample fsamp, float[][][] cents, boolean[] markersForCallrate, boolean[] markersForEverythingElse, GcModel gcModel, GC_CORRECTION_METHOD correctionMethod, Logger log) {
        byte[] abGenotypes, forwardGenotypes;
        float[] lrrs, bafs, bafsWide;
        double abCallRate, forwardCallRate, abHetRate, forwardHetRate, wfPrior, gcwfPrior, wfPost, gcwfPost, lrrsdBound, lrrsdPost, lrrsdPostBound;
        int[] bafBinCounts;
        boolean multimodal;
        
        lrrs = cents == null ? fsamp.getLRRs() : fsamp.getLRRs(cents);
        bafs = cents == null ? fsamp.getBAFs() : fsamp.getBAFs(cents);
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
        abCallRate = abHetRate = 0;
        if (abGenotypes != null) {
            for (int j = 0; j < abGenotypes.length; j++) {
                if (abGenotypes[j] >= 0) {
                    abCallRate++;
                }
                if (abGenotypes[j] == 1) {
                    abHetRate++;
                }
                
            }
            abHetRate /= abCallRate;
            abCallRate /= abGenotypes.length;
        }
        forwardCallRate = forwardHetRate = 0;
        if (forwardGenotypes != null) {
            for (int j = 0; j < forwardGenotypes.length; j++) {
                if (forwardGenotypes[j] > 0) {
                    forwardCallRate++;
                }
                if (forwardGenotypes[j] == 1) {
                    forwardHetRate++;
                }
            }
            forwardHetRate /= forwardCallRate;
            forwardCallRate /= forwardGenotypes.length;
        }
        wfPrior = Double.NaN;
        gcwfPrior = Double.NaN;
        wfPost = Double.NaN;
        gcwfPost = Double.NaN;
        lrrsdPost = Double.NaN;
        lrrsdPostBound = Double.NaN;
        if (gcModel != null) {
        	//TODO, pass markerset to avoid loading
			GcAdjustor gcAdjustor = GcAdjustor.getComputedAdjustor(proj, null, lrrs, gcModel, correctionMethod, true, true, false);
            if (!gcAdjustor.isFail()) {
                wfPrior = gcAdjustor.getWfPrior();
                gcwfPrior = gcAdjustor.getGcwfPrior();
                wfPost = gcAdjustor.getWfPost();
                gcwfPost = gcAdjustor.getGcwfPost();
                double[] tmp;
                if (markersForEverythingElse == null) {
                    lrrsdPost = Array.stdev(gcAdjustor.getCorrectedIntensities(), true);
                    tmp = CNVCaller.adjustLrr(gcAdjustor.getCorrectedIntensities(), CNVCaller.MIN_LRR_MEDIAN_ADJUST, CNVCaller.MAX_LRR_MEDIAN_ADJUST, false, log);
                } else {
                    double[] subLrr = Array.subArray(gcAdjustor.getCorrectedIntensities(), markersForEverythingElse);
                    lrrsdPost = Array.stdev(subLrr, true);
                    tmp = CNVCaller.adjustLrr(subLrr, CNVCaller.MIN_LRR_MEDIAN_ADJUST, CNVCaller.MAX_LRR_MEDIAN_ADJUST, false, log);
                }
                lrrsdPostBound = Array.stdev(Array.getValuesBetween(tmp, CNVCaller.MIN_LRR_MEDIAN_ADJUST, CNVCaller.MAX_LRR_MEDIAN_ADJUST, false), true);
            }
        }
        
        multimodal = Array.isMultimodal(Array.toDoubleArray(Array.removeNaN(bafsWide)), 0.1, 0.5, 0.01);
        double[] tmp = CNVCaller.adjustLrr(Array.toDoubleArray(lrrs), CNVCaller.MIN_LRR_MEDIAN_ADJUST, CNVCaller.MAX_LRR_MEDIAN_ADJUST, false, proj.getLog());
        lrrsdBound = Array.stdev(Array.getValuesBetween(tmp, CNVCaller.MIN_LRR_MEDIAN_ADJUST, CNVCaller.MAX_LRR_MEDIAN_ADJUST, false), true);
        String[] retVals = new String[]{
                sampleID, 
                Array.mean(lrrs, true) + "", 
                Array.stdev(lrrs, true) + "",
                lrrsdBound + "", 
                Array.stdev(bafs, true) + "", 
                (abCallRate > 0 ? abCallRate : forwardCallRate) + "",
                (abCallRate > 0 ? abHetRate : forwardHetRate) + "",
                wfPrior + "",
                gcwfPrior + "",
                wfPost  + "",
                gcwfPost + "",
                lrrsdPost  + "",
                lrrsdPostBound  + "",
                multimodal + "",
                Array.toStr(bafBinCounts),
        };
        return retVals;
	}
	
	
	public void finalAction() {
		String[] files; 
		
//		files = Array.stringArraySequence(numThreads, ext.rootOf(proj.getFilename(proj.SAMPLE_QC_FILENAME), false) + ".");
//		Files.cat(files, proj.getFilename(proj.SAMPLE_QC_FILENAME), Array.intArray(files.length, 0), proj.getLog());
		files = Array.stringArraySequence(numThreads, ext.rootOf(proj.SAMPLE_QC_FILENAME.getValue(), false) + ".");
		Files.cat(files, proj.SAMPLE_QC_FILENAME.getValue(), Array.intArray(files.length, 0), proj.getLog());
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
		GcModel gcModel;
		Logger log;
		
		
		error = false;
		log = proj.getLog();
		samples = proj.getSamples();
		if (customSampleFileList != null) {
			subsamples = HashVec.loadFileToStringArray(customSampleFileList, false, new int[] {0}, false);
			for (int i = 0; i < subsamples.length; i++) {
				if (ext.indexOfStr(subsamples[i], samples) == -1) {
					log.reportError("Error - subsample '"+subsamples[i]+"' was not found in the list of samples of project '"+proj.PROJECT_NAME.getValue()+"'");
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
		gcModel = null;
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
		if (Files.exists(proj.GC_MODEL_FILENAME.getValue(false, false))) {
			gcModel = GcAdjustor.GcModel.populateFromFile(proj.GC_MODEL_FILENAME.getValue(false, false), false, log);
			if (gcModel == null) {
				log.reportError("Error - detected the gc model defined by " + proj.GC_MODEL_FILENAME + " as " + proj.GC_MODEL_FILENAME.getValue(false, false) + " in property file " + proj.getPropertyFilename() + " exists, but an error occurred while loading the file");
				log.reportError("	   - If you would like to skip WF and GCWF qc metrics, either change the " + proj.GC_MODEL_FILENAME + " property to a filename that does not exist, or change the name of " + proj.GC_MODEL_FILENAME.getValue(false, false));
				return;
			}
		} else {
			log.report("Info - did not find gc model file " + proj.GC_MODEL_FILENAME.getValue(false, false) + ", skipping gc correction and related qc");
		}

		threadSeeds = Parallelizable.splitList(samples, numThreads, false);
		runables = new LrrSd[numThreads];
		for (int i = 0; i<numThreads; i++) {
			runables[i] = new LrrSd(proj, threadSeeds[i], markersForCallrate, markersForEverythingElse, centroidsFile, gcModel, i + 1, numThreads);
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

		"   Note: if a gc model is available as defined by the \"GC_MODEL_FILENAME\" property in the project properties file, WF and GCFW (after adjusting for GC content) will be reported\n" +
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
