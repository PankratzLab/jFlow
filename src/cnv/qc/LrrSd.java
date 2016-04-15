package cnv.qc;

import java.io.*;
import java.util.Arrays;
import java.util.Hashtable;

import cnv.filesys.*;
import cnv.filesys.MarkerSet.PreparedMarkerSet;
import cnv.hmm.CNVCaller;
import cnv.manage.MitoPipeline;
import cnv.qc.GcAdjustor.GC_CORRECTION_METHOD;
import cnv.qc.GcAdjustor.GcModel;
import cnv.var.SampleData;
//import java.util.*;
import common.*;

public class LrrSd extends Parallelizable {
	private static final String BOUND_SD = "LRR_SD_" + CNVCaller.MIN_LRR_MEDIAN_ADJUST + "_" + CNVCaller.MAX_LRR_MEDIAN_ADJUST;
	private static final String BOUND_SD_CORRECTED = "LRR_SD_Post_Correction_" + CNVCaller.MIN_LRR_MEDIAN_ADJUST + "_" + CNVCaller.MAX_LRR_MEDIAN_ADJUST;
	private static final String BOUND_MAD = "LRR_MAD_" + CNVCaller.MIN_LRR_MEDIAN_ADJUST + "_" + CNVCaller.MAX_LRR_MEDIAN_ADJUST;
	private static final String BOUND_MAD_CORRECTED = "LRR_MAD_Post_Correction_" + CNVCaller.MIN_LRR_MEDIAN_ADJUST + "_" + CNVCaller.MAX_LRR_MEDIAN_ADJUST;
	
	public static final String[] NUMERIC_COLUMNS = { "LRR_AVG", "LRR_SD", BOUND_SD, "LRR_MAD", BOUND_MAD, "BAF1585_SD", "Genotype_callrate", "Genotype_heterozygosity", "WF_Prior_Correction", "GCWF_Prior_Correction", "WF_Post_Correction", "GCWF_Post_Correction", "LRR_SD_Post_Correction", BOUND_SD_CORRECTED, "LRR_MAD_Post_Correction", BOUND_MAD_CORRECTED };
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
			PreparedMarkerSet markerSet = PreparedMarkerSet.getPreparedMarkerSet(proj.getMarkerSet());
			for (int i = 0; i<samples.length; i++) {
//	        	log.report((i+1)+" of "+samples.length);
				fsamp = proj.getFullSampleFromRandomAccessFile(samples[i]);
				if (fsamp == null) {
					log.reportError("Error - "+samples[i]+Sample.SAMPLE_FILE_EXTENSION+" not found in samples directory");
				} else {
					writer.println(Array.toStr(LrrSdPerSample(proj, markerSet, samples[i], fsamp, cents, markersForCallrate, markersForEverythingElse, gcModel, GC_CORRECTION_METHOD.GENVISIS_GC, log), "\t"));
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
                Array.mad(Array.removeNaN(lrrs))<br/>
                lrrMadBound<br/>
                Array.stdev(bafs, true)<br/> 
                (abCallRate > 0 ? abCallRate : forwardCallRate)<br/>
                (abCallRate > 0 ? abHetRate : forwardHetRate)<br/>
                wfPrior<br/>
                gcwfPrior<br/>
                wfPost<br/>
                gcwfPost <br/>
                lrrsdPost<br/>
                lrrsdPostBound<br/>
                lrrMadPost<br/>
                lrrMadBoundPost<br/>
                multimodal <br/>
                Array.toStr(bafBinCounts)<br/>
	 * 
	 * 
	 * @param proj
	 * @param pMarkerSet
	 * @param sampleID
	 * @param fsamp
	 * @param cents
	 * @param markersForCallrate
	 * @param markersForEverythingElse
	 * @param gcModel
	 * @param log
	 * @return
	 */
	public static String[] LrrSdPerSample(Project proj, PreparedMarkerSet pMarkerSet, String sampleID, Sample fsamp, float[][][] cents, boolean[] markersForCallrate, boolean[] markersForEverythingElse, GcModel gcModel, GC_CORRECTION_METHOD correctionMethod, Logger log) {
        byte[] abGenotypes, forwardGenotypes;
        float[] lrrs, bafs, bafsWide;
		double abCallRate, forwardCallRate, abHetRate, forwardHetRate, wfPrior, gcwfPrior, wfPost, gcwfPost, lrrsdBound, lrrsdPost, lrrsdPostBound, lrrMadBound, lrrMadPost, lrrMadBoundPost;
        int[] bafBinCounts;
        boolean multimodal;
        
        lrrs = cents == null ? fsamp.getLRRs() : fsamp.getLRRs(cents);
        bafs = cents == null ? fsamp.getBAFs() : fsamp.getBAFs(cents);
		bafsWide = bafs;
		markersForEverythingElse = markersForEverythingElse == null ? proj.getAutosomalMarkerBoolean() : markersForEverythingElse;
		if (markersForEverythingElse == null || markersForEverythingElse.length == 0) {
			proj.getLog().reportTimeWarning("Could not determine appropriate marker subset, lrr_sd.xln data for sample " + fsamp.getSampleName() + " will be based on all markers, not just autosomal markers");
			markersForEverythingElse = null;
		} else {
			lrrs = Array.subArray(lrrs, markersForEverythingElse);
			bafs = Array.subArray(bafs, markersForEverythingElse);
			bafsWide = Array.subArray(bafsWide, markersForEverythingElse);
		}
        abGenotypes = fsamp.getAB_Genotypes();
        forwardGenotypes = fsamp.getForwardGenotypes();
        //TODO, remove cnv only probes using proj Array type if markersForCallrate is not provided...
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
		lrrMadPost = Double.NaN;
		lrrMadBoundPost = Double.NaN;
        if (gcModel != null) {
			GcAdjustor gcAdjustor = GcAdjustor.getComputedAdjustor(proj, pMarkerSet, cents == null ? fsamp.getLRRs() : fsamp.getLRRs(cents), gcModel, correctionMethod, true, true, false);
            if (!gcAdjustor.isFail()) {
                wfPrior = gcAdjustor.getWfPrior();
                gcwfPrior = gcAdjustor.getGcwfPrior();
                wfPost = gcAdjustor.getWfPost();
                gcwfPost = gcAdjustor.getGcwfPost();
                double[] tmp;
                if (markersForEverythingElse == null) {
                	
                    lrrsdPost = Array.stdev(gcAdjustor.getCorrectedIntensities(), true);
                    lrrMadPost = Array.mad(Array.removeNaN(gcAdjustor.getCorrectedIntensities()));
                    tmp = CNVCaller.adjustLrr(gcAdjustor.getCorrectedIntensities(), CNVCaller.MIN_LRR_MEDIAN_ADJUST, CNVCaller.MAX_LRR_MEDIAN_ADJUST, false, log);
				} else {
					double[] subLrr = Array.subArray(gcAdjustor.getCorrectedIntensities(), markersForEverythingElse);
					lrrsdPost = Array.stdev(subLrr, true);
					lrrMadPost = Array.mad(Array.removeNaN(subLrr));
					tmp = CNVCaller.adjustLrr(subLrr, CNVCaller.MIN_LRR_MEDIAN_ADJUST, CNVCaller.MAX_LRR_MEDIAN_ADJUST, false, log);
				}
				tmp = Array.removeNaN(Array.getValuesBetween(tmp, CNVCaller.MIN_LRR_MEDIAN_ADJUST, CNVCaller.MAX_LRR_MEDIAN_ADJUST, false));
				lrrsdPostBound = Array.stdev(tmp, true);
				lrrMadBoundPost = Array.mad(tmp);
			}
		}
        
		multimodal = Array.isMultimodal(Array.toDoubleArray(Array.removeNaN(bafsWide)), 0.1, 0.5, 0.01);
	    double[] dlrrs = Array.toDoubleArray(lrrs);
		double[] tmp = CNVCaller.adjustLrr(dlrrs, CNVCaller.MIN_LRR_MEDIAN_ADJUST, CNVCaller.MAX_LRR_MEDIAN_ADJUST, false, proj.getLog());
		tmp = Array.removeNaN(Array.getValuesBetween(tmp, CNVCaller.MIN_LRR_MEDIAN_ADJUST, CNVCaller.MAX_LRR_MEDIAN_ADJUST, false));
		lrrsdBound = Array.stdev(tmp, true);
		lrrMadBound = Array.mad(tmp);
		
        String[] retVals = new String[]{
                sampleID, 
                Array.mean(lrrs, true) + "", 
                Array.stdev(lrrs, true) + "",
                lrrsdBound + "", 
				Array.mad(Array.removeNaN(dlrrs)) + "",
				lrrMadBound + "",
                Array.stdev(bafs, true) + "", 
                (abCallRate > 0 ? abCallRate : forwardCallRate) + "",
                (abCallRate > 0 ? abHetRate : forwardHetRate) + "",
                wfPrior + "",
                gcwfPrior + "",
                wfPost  + "",
                gcwfPost + "",
                lrrsdPost  + "",
                lrrsdPostBound  + "",
				lrrMadPost + "",
				lrrMadBoundPost + "",
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
	

    /**
     * If the useFile is not null, we return a hash with the subset of individuals. Can return null if useFile does not exist or does not contain any individuals
     * 
     * @param useFile
     * @param log
     * @return
     */
    private static Hashtable<String, String> checkSubset(String useFile, Logger log) {
        Hashtable<String, String> subset = new Hashtable<String, String>();
        if (useFile != null) {
            if (Files.exists(useFile)) {
                subset = HashVec.loadFileToHashString(useFile, 0, new int[] { 0 }, null, false, false);
                if (subset.size() == 0) {
                    log.reportError("Error - did not find any samples in the subset file " + useFile);
                    return null;
                } else {
                    log.report("Analysis will be performed starting with the subset of " + subset.size() + " samples found in " + useFile);
                }
            } else {
                log.reportError("Error - a file list of samples to use was provided, but the file " + useFile + " does not exist");
                return null;
            }
        } else {
            log.report("Info - A subset of samples was not provided with the \"useFile=\" argument, using all parsed samples as input to the SVD");
        }
        return subset;
    }

    /**
     * Check to make sure that sample data has DNA header, and that the QC has not already been added
     */
    private static boolean checkSampleData(Project proj, SampleData sampleData) {
        // This should not happen, but if it is the case we will not attempt to add qc metrics
        boolean addToSampleData = true;
        Logger log = proj.getLog();

        if (!sampleData.containsDNA()) {
            addToSampleData = false;
            log.reportError("Error - sample data did not contain column with header \"DNA\", not adding sample qc summaries to sample data");
        }
        if (qcAdded(proj)) {
            addToSampleData = false;
            log.reportError("Detected that sample data QC metrics have been added already, will not add these again");
            // log.reportError("If new thresholds were used, please remove the columns [" + ext.listWithCommas(SAMPLE_DATA_ADDITION_HEADERS, true) + "] in " + proj.getFilename(proj.SAMPLE_DATA_FILENAME));
            log.reportError("If new thresholds were used, please remove the columns [" + ext.listWithCommas(MitoPipeline.SAMPLE_DATA_ADDITION_HEADERS, true) + "] in " + proj.SAMPLE_DATA_FILENAME.getValue());
        }
        return addToSampleData;
    }

    /**
     * Check all indices for -1 status
     */
    private static boolean checkIndices(Project proj, int[] indices) {
        boolean allGood = true;
        for (int i = 0; i < indices.length; i++) {
            if (indices[i] == -1) {
                allGood = false;
                proj.getLog().reportError("Error - The sample QC file " +proj.SAMPLE_QC_FILENAME.getValue() + " did not contain the proper headings, this should not happen");
            }
        }
        return allGood;
    }

    /**
     * Check the header of the sample data file to see if the sample data qc headers are present
     */
    private static boolean qcAdded(Project proj) {
        boolean added = true;
        // String[] header = Files.getHeaderOfFile(proj.getFilename(proj.SAMPLE_DATA_FILENAME), proj.getLog());
        String[] header = Files.getHeaderOfFile(proj.SAMPLE_DATA_FILENAME.getValue(), proj.getLog());
        int[] indices = ext.indexFactors(MitoPipeline.SAMPLE_DATA_ADDITION_HEADERS, header, true, proj.getLog(), false, false);
        for (int i = 0; i < indices.length; i++) {
            if (indices[i] < 0) {
                added = false;
            }
        }
        return added;
    }
	
	/**
     * @param proj
     *            current project
     * @param outputBase
     * @param markersForABCallRate
     * @param markersForEverythingElse
     * @param numThreads
     *            threads for LRR_SD
     * @param sampleCallRateFilter
     *            filter samples by this call rate (LRR_SD filter is set in project)
     * @param computeQC
     * @param useFile
     *            a further filter of samples that will be used
     * @param log
     */
	public static int[] filterSamples(Project proj, String outputBase, String markersForABCallRate, String markersForEverythingElse, int numThreads, String useFile, boolean gcMetrics) {
        Hashtable<String, String> sampDataQC = new Hashtable<String, String>();
        int[] indices;
        String[] line;
        int numPassing, count;
        Logger log = proj.getLog();
        String delim = "\t";

        SampleData sampleData = proj.getSampleData(0, false);
        // double lrrSdFilter = Double.parseDouble(proj.getProperty(proj.LRRSD_CUTOFF));
        double lrrSdFilter = proj.LRRSD_CUTOFF.getValue();
        double callRateFilter = proj.SAMPLE_CALLRATE_THRESHOLD.getValue();
        boolean addToSampleData = checkSampleData(proj, sampleData);
        Hashtable<String, String> subset = checkSubset(useFile, log);

        if (Files.exists(proj.SAMPLE_QC_FILENAME.getValue())) {
            log.report("The sample qc file " + proj.SAMPLE_QC_FILENAME.getValue() + " already exists");
            log.report("Skipping qc computation, filtering on existing qc file " + proj.SAMPLE_QC_FILENAME.getValue());
        } else {
            log.report("Computing sample QC for all samples...");
            log.report("Will be reporting sample qc to " + proj.SAMPLE_QC_FILENAME.getValue());
			cnv.qc.LrrSd.init(proj, null, markersForABCallRate, markersForEverythingElse, numThreads, null, gcMetrics);
            if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }
        }

        count = 0;
        numPassing = 0;
        try {
            if (Thread.currentThread().isInterrupted()) { throw new RuntimeException(new InterruptedException()); }
            BufferedReader reader = Files.getReader(proj.SAMPLE_QC_FILENAME.getValue(), false, true, false);
            PrintWriter writerUse = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + outputBase + MitoPipeline.PCA_SAMPLES));
            PrintWriter writerSummary = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + outputBase + MitoPipeline.PCA_SAMPLES_SUMMARY));

            writerSummary.println(Array.toStr(MitoPipeline.SAMPLE_QC_SUMMARY));
            if (!reader.ready()) {
                writerUse.close();
                writerSummary.close();
                reader.close();
                log.reportError("Error - QC file (" + proj.SAMPLE_QC_FILENAME.getValue() + ") was empty");
                return new int[] { numPassing, count };
            }
            line = reader.readLine().trim().split(delim);
            indices = ext.indexFactors(MitoPipeline.QC_COLUMNS, line, true, log, true, false);

            if (!checkIndices(proj, indices)) {
                writerUse.close();
                writerSummary.close();
                reader.close();
                log.reportError("Error - could not detect proper header in QC file (" + proj.SAMPLE_QC_FILENAME.getValue() + ")");
                return null;
            }

            while (reader.ready()) {
                line = reader.readLine().trim().split(delim);
                // skip any headers as a result of concatenating the qc results
                if (!line[indices[0]].equals(MitoPipeline.QC_COLUMNS[0])) {
                    // if passes qc
                    if (Double.parseDouble(line[indices[1]]) < lrrSdFilter && Double.parseDouble(line[indices[2]]) > callRateFilter) {
                        sampDataQC.put(line[indices[0]], line[indices[1]] + "\t" + line[indices[2]] + "\t" + "0");
                        // check the subset
                        if (subset.size() == 0 || subset.containsKey(line[indices[0]])) {
                            writerUse.println(line[indices[0]]);
                            writerSummary.println(line[indices[0]] + "\t" + line[indices[1]] + "\t" + line[indices[2]] + "\t" + "TRUE");
                            numPassing++;
                        } else {
                            // sampDataQC.put(line[indices[0]], line[indices[1]] + "\t" + line[indices[2]] + "\t" + "1");
                            writerSummary.println(line[indices[0]] + "\t" + line[indices[1]] + "\t" + line[indices[2]] + "\t" + "FALSE");
                        }
                    } else {
                        sampDataQC.put(line[indices[0]], line[indices[1]] + "\t" + line[indices[2]] + "\t" + "1");
                        writerSummary.println(line[indices[0]] + "\t" + line[indices[1]] + "\t" + line[indices[2]] + "\t" + "FALSE");
                    }
                    count++;
                }
            }

            reader.close();
            writerUse.close();
            writerSummary.close();

            if (numPassing == 0) {
                log.reportError("Error - all Samples were filtered out by the QC step");
                log.reportError("If there are a large number of cnv-only probes on the array, try lowering the call rate threshold for samples or use the \"-markerQC\" option to only select markers with high quality call rates");
                return new int[] { numPassing, count };
            }
            log.report("Info - " + numPassing + " " + (numPassing == 1 ? "sample" : "samples") + " passed the QC threshold" + (subset.size() > 0 ? " and were present in the subset file " + useFile : ""));
        } catch (FileNotFoundException fnfe) {
            log.reportError("Error: file \"" + proj.SAMPLE_QC_FILENAME.getValue() + "\" not found in current directory");
        } catch (IOException ioe) {
			log.reportError("Error reading file \"" + proj.SAMPLE_QC_FILENAME.getValue() + "\"");
		}

		if (addToSampleData) {
			sampleData.addData(sampDataQC, MitoPipeline.DNA_LINKER, Array.tagOn(MitoPipeline.SAMPLE_DATA_ADDITION_HEADERS, outputBase, null), ext.MISSING_VALUES[1], delim, log);
		}
		return new int[] { numPassing, count };
	}

	private static boolean[] getMarkerSubset(Project proj, String[] subMarkers, String[] markers) {
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

	public static void init(Project proj, String customSampleFileList, String centroidsFile, int numThreads, boolean useAllMarkers) {
	    boolean[] callRate = null;
	    boolean[] theRest = null;
	    
	    if (!useAllMarkers) {
	        callRate = Array.booleanNegative(proj.getCNMarkers());
	        theRest = proj.getAutosomalMarkerBoolean();
	        if(callRate.length != theRest.length) {
	            proj.getLog().reportError("Error -- array lengths differ between proj.getCNMarkers() and proj.getAutosomalMarkerBoolean().  Please report or fix this.");
	            return;
	        }
	        for (int i = 0; i < callRate.length; i++) {
	            if (!callRate[i] && theRest[i]) {
	                theRest[i] = false;
	            }
	        }
	    }
	    
	    init(proj, customSampleFileList, callRate, theRest, centroidsFile, numThreads);
	}
	
	public static void init(Project proj, String customSampleFileList, String markersForCallrateFile, String markersForEverythingElseFile, int numThreads, String centroidsFile, boolean gcMetrics) {
    	String[] markers;
    	boolean[] markersForCallrate, markersForEverythingElse;
    	BaselineUnclusteredMarkers bum;
    	Logger log;
    	
    	log = proj.getLog();
    	
    	if(!BaselineUnclusteredMarkers.baselineUnclusteredMarkersFileExists(proj)){
    		log.report("Baseline Unclustered Markers file does not exist and will be created now");
    		if(!BaselineUnclusteredMarkers.createBaselineUnclusteredMarkersFileFromSamples(proj)){
    			log.reportError("Error - Baseline Unclustered Markers file could not be created");
    			return;
    		}
    	}
    	
    	markers = proj.getMarkerNames();
    	if (markersForCallrateFile != null) {
    		markersForCallrate = getMarkerSubset(proj,HashVec.loadFileToStringArray(markersForCallrateFile, false, new int[] { 0 }, false), markers);
    		if (markersForCallrate == null) {
    			log.reportError("Error - Some markers listed in " + markersForCallrateFile + " were not found in the current project, or were duplicates");
    			return;
    		}
    	} else {
    		markersForCallrate = new boolean[markers.length];
    		Arrays.fill(markersForCallrate, true);
    	}
    	bum = BaselineUnclusteredMarkers.getProjBaselineUnclusteredMarkers(proj);
    	for (int i = 0; i < markersForCallrate.length; i++){
    		if(markersForCallrate[i] && bum.markerUnclustered(markers[i])) {
    			markersForCallrate[i] = false;
    		}
    	}
    	markersForEverythingElse = null;
    	if (markersForEverythingElseFile != null) {
    		markersForEverythingElse = getMarkerSubset(proj, HashVec.loadFileToStringArray(markersForEverythingElseFile, false, new int[] { 0 }, false), markers);
    		if (markersForEverythingElse == null) {
    			log.reportError("Error - Some markers listed in " + markersForEverythingElseFile + " were not found in the current project, or were duplicates");
    			return;
    		}
    	}
    	
		init(proj, customSampleFileList, markersForCallrate, markersForEverythingElse, centroidsFile, gcMetrics, numThreads);
    }
	public static void init(Project proj, String customSampleFileList, boolean[] markersForCallrate, boolean[] markersForEverythingElse, String centroidsFile, int numThreads) {
		init(proj, customSampleFileList, markersForCallrate, markersForEverythingElse, centroidsFile, true, numThreads);
	}

	public static void init(Project proj, String customSampleFileList, boolean[] markersForCallrate, boolean[] markersForEverythingElse, String centroidsFile, boolean gcMetrics, int numThreads) {
	    String[] samples, subsamples;
        String[][] threadSeeds;
        LrrSd[] runables;
        boolean error;
        String[] markers;
        BaselineUnclusteredMarkers bum;
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
        
        if(!BaselineUnclusteredMarkers.baselineUnclusteredMarkersFileExists(proj)){
            log.report("Baseline Unclustered Markers file does not exist and will be created now");
            if(!BaselineUnclusteredMarkers.createBaselineUnclusteredMarkersFileFromSamples(proj)){
                log.reportError("Error - Baseline Unclustered Markers file could not be created");
                return;
            }
        }
        
        markers = proj.getMarkerNames();
        if (markersForCallrate == null) {
            markersForCallrate = new boolean[markers.length];
            Arrays.fill(markersForCallrate, true);
        }
        bum = BaselineUnclusteredMarkers.getProjBaselineUnclusteredMarkers(proj);
        for (int i = 0; i < markersForCallrate.length; i++){
            if(markersForCallrate[i] && bum.markerUnclustered(markers[i])) {
                markersForCallrate[i] = false;
            }
        }

        gcModel = null;
		if (gcMetrics && Files.exists(proj.GC_MODEL_FILENAME.getValue(false, false))) {
            gcModel = GcAdjustor.GcModel.populateFromFile(proj.GC_MODEL_FILENAME.getValue(false, false), false, log);
            if (gcModel == null) {
                log.reportError("Error - detected the gc model defined by " + proj.GC_MODEL_FILENAME + " as " + proj.GC_MODEL_FILENAME.getValue(false, false) + " in property file " + proj.getPropertyFilename() + " exists, but an error occurred while loading the file");
                log.reportError("      - If you would like to skip WF and GCWF qc metrics, either change the " + proj.GC_MODEL_FILENAME + " property to a filename that does not exist, or change the name of " + proj.GC_MODEL_FILENAME.getValue(false, false));
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
		boolean filter = false;
		boolean projectDefinedMarkers = false;
		String sampleCallRateFilter = null;
		String outputBase = null;

		String usage = "\n"+
		"cnv.qc.LrrSd requires 0-6 arguments\n"+
		"   (1) project properties filename (i.e. proj="+cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n"+
		"   (2) centroids with which to compute LRRs (i.e. cents=genotype.cent (not the default; to be found in data/ directory))\n"+
		"   (3) number of threads to use (i.e. threads="+numThreads+" (default))\n"+
		"   (4) optional: if you only want to look at a subset of the samples, filename of sample list (i.e. subsample=these.txt (not the default))\n"+
		"   (5) optional: if you only want to compute AB_callrate and Forward_callrate from a subset of the markers, filename of marker list (i.e. callRateMarkers=those.txt (not the default))\n"+
		"   (6) optional: if you only want to compute the other qc metrics (excluding AB_callrate and Forward_callrate) from a subset of the markers, filename of marker list (i.e. otherMarkers=this.txt (not the default))\n"+
		"   (7) optional: if you only want to compute metrics based on autosomal and non-CN markers (i.e. projectMarkers=TRUE (not the default))\n" + 
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
			} else if (args[i].startsWith("-filter")) {
			    filter = true;
			    numArgs--;
			} else if (args[i].startsWith("outBase=")) {
			    outputBase = args[i].split("=")[1];
			    numArgs--;
			} else if (args[i].startsWith("callRate=")) {
			    sampleCallRateFilter = args[i].split("=")[1];
			    numArgs--;
			} else if (args[i].startsWith("projectMarkers=")) {
			    projectDefinedMarkers = ext.parseBooleanArg(args[i]);
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
		    proj = new Project(filename, false);
			if (filter) {
				proj.SAMPLE_CALLRATE_THRESHOLD.setValue(Double.parseDouble(sampleCallRateFilter));
				filterSamples(proj, outputBase, markersForCallrateFile, markersForEverythingElseFile, numThreads, filenameOfListOfSamples, true);
			} else {
			    if (projectDefinedMarkers) {
			        init(proj, filenameOfListOfSamples, centroids, numThreads, false);
			    } else {
					init(proj, filenameOfListOfSamples, markersForCallrateFile, markersForEverythingElseFile, numThreads, centroids, true);
			    }
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
