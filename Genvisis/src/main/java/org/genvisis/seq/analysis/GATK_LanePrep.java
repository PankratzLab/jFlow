package org.genvisis.seq.analysis;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.seq.analysis.GATK.BaseRecalibration;

/**
 * Prepping for the GATK is done on a lane by lane basis as reflected here
 *
 */
public class GATK_LanePrep extends BWA_Analysis {
	private static final String PICARD_METRICS_SUMMARY = "picard_metrics_summary.txt";
	private final Picard picard;
	private Picard.Picard_Analysis[] picard_Analysis;
	private final GATK gatk;
	private GATK.BaseRecalibration[] gRecalibrations;
	private MergeBam.BamMerger[] mBamMergers;
	private final MergeBam mergeBam;

	public GATK_LanePrep(	String rootInputDir, String rootOutputDir, String referenceGenomeFasta,
												boolean verbose, int numWithinSampleThreads, int numBetweenSampleThreads,
												BWA bwa, Picard picard, GATK gatk, MergeBam mergeBam, Logger log) {
		super(rootInputDir, (rootOutputDir == null ? rootInputDir : rootOutputDir),
					referenceGenomeFasta, verbose, numWithinSampleThreads, numBetweenSampleThreads, bwa, log);
		this.picard = picard;
		this.gatk = gatk;
		this.mergeBam = mergeBam;
	}

	public MergeBam.BamMerger[] getmBamMergers() {
		return mBamMergers;
	}

	public void runBWA(String fileOfSamplePairs) {
		init(fileOfSamplePairs);
		if (!isFail()) {
			analyzeBWA_MEM();
		}
	}

	public void resetBwAnalysisIndividuals(MergeBam.BamMerger[] mBamMergers) {
		boolean[] mergeMask = getMergeMask(mBamMergers, getLog());
		if (isVerbose()) {
			getLog().report(ext.getTime()	+ " Info - " + Array.booleanArraySum(mergeMask) + " of "
											+ mBamMergers.length
											+ " sample(s) will go through another round of de-duping and realigining");
		}
		BWA_AnalysisIndividual[] bwAnalysisIndividuals = new BWA_AnalysisIndividual[Array.booleanArraySum(mergeMask)];
		int index = 0;
		for (int i = 0; i < mBamMergers.length; i++) {
			if (mergeMask[i]) {
				bwAnalysisIndividuals[index] = new BWA_AnalysisIndividual(null, getRootOutputDir(),
																																	mBamMergers[i].getBaseId(), null,
																																	null, null, null, getLog());
				bwAnalysisIndividuals[index].setOutput(mBamMergers[i].getOutputBam());
				index++;
			}
		}
		setBwAnalysisIndividuals(bwAnalysisIndividuals);
	}

	public void runPicard() {
		if (!isFail()) {
			BWA_AnalysisIndividual[] bwAnalysisIndividuals = getBwAnalysisIndividuals();
			if (bwAnalysisIndividuals != null) {
				picard_Analysis = new Picard.Picard_Analysis[bwAnalysisIndividuals.length];
				double memoryRatio = (double) 1 / getNumBetweenSampleThreads();// added this because with
																																				// more than 4 samples,

				memoryRatio -= .01;
				if (memoryRatio > Picard.DEFAULT_SORTING_COLLECTION_SIZE_RATIO) {
					memoryRatio = Picard.DEFAULT_SORTING_COLLECTION_SIZE_RATIO;
				} else {
					getLog().report(ext.getTime()	+ " Info - adjusting Picard's memory ratio to "
													+ memoryRatio + " since there are more than 3 samples...");
				}
				ExecutorService executor = Executors.newFixedThreadPool(getNumBetweenSampleThreads());
				Hashtable<String, Future<Picard.Picard_Analysis>> tmpResults = new Hashtable<String, Future<Picard.Picard_Analysis>>();
				for (int i = 0; i < bwAnalysisIndividuals.length; i++) {
					Logger altLog = new Logger(ext.rootOf(getLog().getFilename(), false)	+ "_Picard_ID_"
																			+ bwAnalysisIndividuals[i].getID() + "_Lane_"
																			+ bwAnalysisIndividuals[i].getLane() + ".log");
					tmpResults.put(i	+ "",
													executor.submit(new WorkerPicard(	picard, bwAnalysisIndividuals[i].getID(),
																														bwAnalysisIndividuals[i].getOutput(),
																														memoryRatio, altLog)));
				}
				for (int i = 0; i < bwAnalysisIndividuals.length; i++) {
					try {
						getLog().memoryPercentFree();
						picard_Analysis[i] = tmpResults.get(i + "").get();
						getLog().report(ext.getTime()	+ "Info - retrieving picard results for "
														+ picard_Analysis[i].getFullPathToSamFile());
						if (picard_Analysis[i].isFail() && !isFail()) {
							getLog().reportError("Error - failed picard for "
																		+ picard_Analysis[i].getFullPathToSamFile());
							setFail(true);
						}
					} catch (InterruptedException e) {
						getLog().reportError("Could not complete running Picard on internal index " + i);
						getLog().reportException(e);
						setFail(true);
					} catch (ExecutionException e) {
						getLog().reportError("Could not complete running Picard on internal index " + i);
						getLog().reportException(e);
						setFail(true);
					}
				}
				executor.shutdown();
				try {
					executor.awaitTermination(10, TimeUnit.DAYS);
				} catch (InterruptedException e) {
					getLog().reportException(e);
				}
				String[] picardFiles = new String[picard_Analysis.length];
				for (int i = 0; i < picard_Analysis.length; i++) {
					if (picard_Analysis[i].isAllThere()) {
						picardFiles[i] = picard_Analysis[i].getFullPathToMetricsTxt();
					} else {
						getLog().reportError("Error - could not find picard metrics file for root input files:\n"
																	+ bwAnalysisIndividuals[i].getAvailableFiles("\n"));
					}
				}
				getLog().report(ext.getTime() + "Info - parsing picard metrics files");

				Picard.PicardMetricsParser pMetricsParser = new Picard.PicardMetricsParser(	picardFiles,
																																										getLog());
				pMetricsParser.parse(getRootOutputDir() + PICARD_METRICS_SUMMARY);
				getLog().report(ext.getTime() + "Info - finished parsing picard metrics files");

			} else {
				// TODO better check
			}
		}
	}

  public String[] getCurrentFiles(boolean dedupStage) {
    if (!dedupStage) {
      getLog().reportError("Error - should not be requesting other files with un implemented request");
      return null;
    } else {
      String[] currentFiles = new String[picard_Analysis.length];
      for (int i = 0; i < currentFiles.length; i++) {
        currentFiles[i] = picard_Analysis[i].getFullPathToSortedDeDuppedBamFile();
      }
      return currentFiles;
    }
  }

	public String[] getFilesNotMerged() {
		boolean[] merged = getMergeMask(mBamMergers, getLog());
		ArrayList<String> filesNotMerged = new ArrayList<String>();
		for (int i = 0; i < merged.length; i++) {
			if (!merged[i] && !mBamMergers[i].shouldMerge()) {
				filesNotMerged.add(mBamMergers[i].getInputBams()[0]);// only one file
			}
		}
		if (filesNotMerged.size() > 0) {
			return filesNotMerged.toArray(new String[filesNotMerged.size()]);
		} else {
			return null;
		}
	}

  // public void runIndelRealign() {
  // if (!isFail()) {
  // if (picard_Analysis != null) {
  // gIndelPreps = new GATK.IndelPrep[picard_Analysis.length];
  // ExecutorService executor = Executors.newFixedThreadPool(getNumBetweenSampleThreads());
  // Hashtable<String, Future<GATK.IndelPrep>> tmpResults =
  // new Hashtable<String, Future<GATK.IndelPrep>>();
  // for (int i = 0; i < picard_Analysis.length; i++) {
  // getLog().report(ext.getTime() + "Info - beginning realignment for "
  // + picard_Analysis[i].getFullPathToSortedDeDuppedBamFile());
  // Logger altLog = new Logger(ext.rootOf(getLog().getFilename(), false) + "_IndelPrep_ID_"
  // + getBwAnalysisIndividuals()[i].getID() + "_Lane_"
  // + getBwAnalysisIndividuals()[i].getLane() + ".log");
  // tmpResults.put(i + "",
  // executor.submit(new WorkerIndel(gatk, picard_Analysis[i].getBaseID(),
  // picard_Analysis[i].getFullPathToSortedDeDuppedBamFile(),
  // altLog)));
  // }
  // for (int i = 0; i < picard_Analysis.length; i++) {
  // try {
  // gIndelPreps[i] = tmpResults.get(i + "").get();
  // if (gIndelPreps[i].isFail() && !isFail()) {
  // getLog().reportError("Error - failed indel re-alignment for "
  // + gIndelPreps[i].getDedup_reads_bam());
  // setFail(true);
  // }
  // } catch (InterruptedException e) {
  // getLog().reportError("Error - could running GATK indel Prep on internal index " + i);
  // getLog().reportException(e);
  // setFail(true);
  // } catch (ExecutionException e) {
  // getLog().reportError("Error - could running GATK indel Prep on internal index " + i);
  // getLog().reportException(e);
  // setFail(true);
  // }
  // }
  // executor.shutdown();
  // try {
  // executor.awaitTermination(10, TimeUnit.DAYS);
  // } catch (InterruptedException e) {
  // getLog().reportException(e);
  // }
  // } else {
  // // TODO better check
  // }
  // }
  // }

  public void runBaseRecal() {
    if (!isFail()) {
      if (picard_Analysis != null) {
        gRecalibrations = new GATK.BaseRecalibration[picard_Analysis.length];
        ExecutorService executor = Executors.newFixedThreadPool(getNumBetweenSampleThreads());
        Hashtable<String, Future<GATK.BaseRecalibration>> tmpResults =
                                                                     new Hashtable<String, Future<GATK.BaseRecalibration>>();
        for (int i = 0; i < picard_Analysis.length; i++) {
          Logger altLog = new Logger(ext.rootOf(getLog().getFilename(), false)
                                     + "_BaseRecalibration_ID_"
                                     + getBwAnalysisIndividuals()[i].getID() + "_Lane_"
                                     + getBwAnalysisIndividuals()[i].getLane() + ".log");
          tmpResults.put(i + "",
                         executor.submit(new WorkerRecalibration(gatk,
                                                                 picard_Analysis[i].getBaseID(),
                                                                 picard_Analysis[i].getFullPathToSortedDeDuppedBamFile(),
                                                                 altLog)));
        }
        for (int i = 0; i < picard_Analysis.length; i++) {
          try {
            gRecalibrations[i] = tmpResults.get(i + "").get();
            if (gRecalibrations[i].isFail() && !isFail()) {
              getLog().reportError("Error - failed recalibration for "
                                   + gRecalibrations[i].getDedup_reads_bam());
              setFail(true);
            }
          } catch (InterruptedException e) {
            getLog().reportError("Error - when running GATK Base recalibration on internal index "
                                 + i);
            getLog().reportException(e);
            setFail(true);
          } catch (ExecutionException e) {
            getLog().reportError("Error - when running GATK Base recalibration on internal index "
                                 + i);
            getLog().reportException(e);
            setFail(true);
          }
        }
        executor.shutdown();
        try {
          executor.awaitTermination(10, TimeUnit.DAYS);
        } catch (InterruptedException e) {
          getLog().reportException(e);
        }
      } else {
        // TODO better check
      }
    }
  }

	public void runBamMerge() {
		System.out.println(isFail() + "\t" + mergeBam.isFail());

		// TODO get uniq files
		if (!isFail() && !mergeBam.isFail()) {
			if (gRecalibrations != null) {
				GATK.BaseRecalibration[][] gRecalibrationsToMerge = getCalibrationsToMerge(	gRecalibrations,
																																										mergeBam,
																																										getLog());
				mBamMergers = new MergeBam.BamMerger[gRecalibrationsToMerge.length];
				ExecutorService executor = Executors.newFixedThreadPool(getNumBetweenSampleThreads());
				Hashtable<String, Future<MergeBam.BamMerger>> tmpResults = new Hashtable<String, Future<MergeBam.BamMerger>>();
				for (int i = 0; i < mBamMergers.length; i++) {
					Logger altLog = new Logger(ext.rootOf(getLog().getFilename(), false)
																				+ "_BaseRecalibration_ID_"
																			+ getBwAnalysisIndividuals()[i].getID() + "_Lane_"
																			+ getBwAnalysisIndividuals()[i].getLane() + ".log");
					tmpResults.put(i	+ "",
													executor.submit(new WorkerBamMerger(mergeBam,
																															getAllBaseIds(gRecalibrationsToMerge[i]),
																															gRecalibrationsToMerge[i][0].getNewBaseId(),
																															getRootOutputDir(),
																															MergeBam.BamMerger.MERGE_STAGES[0],
																															getInputFilesFrom(gRecalibrationsToMerge[i]),
																															altLog)));
				}
				for (int i = 0; i < mBamMergers.length; i++) {
					try {
						mBamMergers[i] = tmpResults.get(i + "").get();
						if ((mergeBam.isFail() || mBamMergers[i].isFail()) && !isFail()) {
							getLog().reportError("Error - failed merging for "
																		+ Array.toStr(mBamMergers[i].getInputBams(), "\n"));
							setFail(true);
						}
					} catch (InterruptedException e) {
						getLog().reportError("Error - when running GATK Base recalibraion on internal index "
																	+ i);
						getLog().reportException(e);
						setFail(true);
					} catch (ExecutionException e) {
						getLog().reportError("Error - when running GATK Base recalibraion on internal index "
																	+ i);
						getLog().reportException(e);
						setFail(true);
					}
				}
				executor.shutdown();
				try {
					executor.awaitTermination(10, TimeUnit.DAYS);
				} catch (InterruptedException e) {
					getLog().reportException(e);
				}
			} else {
				// TODO better check
			}
		}
	}

	private static String[] getAllBaseIds(GATK.BaseRecalibration[] gRecalibrationsToMerge) {
		String[] uniq = new String[gRecalibrationsToMerge.length];
		for (int i = 0; i < uniq.length; i++) {
			uniq[i] = gRecalibrationsToMerge[i].getBaseId();
		}
		return Array.unique(uniq);
	}

	public Picard getPicard() {
		return picard;
	}

	public GATK getGatk() {
		return gatk;
	}

	public MergeBam getMergeBam() {
		return mergeBam;
	}

	public int getNumOtherThreads() {
		return getNumBetweenSampleThreads();
	}

	@Override
	public void batch(int numBatches, int memoryInMB, int wallTimeInHours, String baseName) {
		String[] batchesByLane = BWA_AnalysisIndividual.getBatchesByLane(getBwAnalysisIndividuals());// we
																																																	// force
																																																	// different
																																																	// lanes
																																																	// from
																																																	// the
																																																	// same
																																																	// sample
																																																	// to
																																																	// be
																																																	// in
																																																	// the
																																																	// same
																																																	// batch...for
																																																	// downstream
																																																	// merging
		String[][] batchedMatchedFiles = Array.splitUpStringArray(batchesByLane, numBatches, getLog());
		String[][] batches = new String[batchedMatchedFiles.length][1];
		for (int i = 0; i < batches.length; i++) {
			batches[i][0] = "batch_" + i + "_" + baseName;
			Files.writeArray(batchedMatchedFiles[i], getRootOutputDir() + batches[i][0] + ".txt");
		}
		// TODO, change classpath
		String command = Array.toStr(PSF.Load.getAllModules(), "\n");
		command += "\njava -Xmx"	+ memoryInMB + "m -jar ~/genvisisGATK.jar seq.analysis.GATK_LanePrep "
								+ ROOT_INPUT_COMMAND + getRootInputDir() + SPACE + ROOT_OUTPUT_COMMAND
								+ getRootOutputDir() + SPACE;
		command += REFERENCE_GENOME_COMMAND	+ getReferenceGenomeFasta() + SPACE + BWA_LOCATION_COMMAND
								+ getBwa().getBwaLocation() + SPACE;
		command += NUM_BETWEEN_THREADS_COMMAND	+ getNumBetweenSampleThreads() + SPACE
								+ FILE_OF_SAMPLE_PAIRS_COMMAND + getRootOutputDir() + "[%0].txt" + SPACE
								+ NUM_WITHIN_THREADS_COMMAND + getNumWithinSampleThreads() + SPACE;
		command += MergeBam.SAMTOOLS_COMMAND + getMergeBam().getSamtoolsLocation() + SPACE;
		command += Picard.PICARD_LOCATION_COMMAND + getPicard().getPicardLocation() + SPACE;
		command += GATK.GATK_LOCATION_COMMAND + getGatk().getGATKLocation() + SPACE;
		command += GATK.KNOWN_SITES_SNP_LOCATION_COMMAND
								+ Array.toStr(getGatk().getKnownSitesSnpFile(), GATK.SPLIT) + SPACE;
		command += GATK.KNOWN_SITES_INDEL_LOCATION_COMMAND
								+ Array.toStr(getGatk().getKnownSitesIndelFile(), GATK.SPLIT);
		Files.qsub("GATK_Lane_Prep"	+ baseName, command, batches, memoryInMB, wallTimeInHours,
								getNumWithinSampleThreads() * getNumBetweenSampleThreads());
	}

	private static class WorkerPicard implements Callable<Picard.Picard_Analysis> {
		private final Picard picard;
		private final String fullPathToSamFile, baseId;
		private final Logger altLog;
		private final double memoryRatio;

		public WorkerPicard(Picard picard, String baseId, String fullPathToSamFile, double memoryRatio,
												Logger altLog) {
			super();
			this.picard = picard;
			this.baseId = baseId;
			this.fullPathToSamFile = fullPathToSamFile;
			this.memoryRatio = memoryRatio;
			this.altLog = altLog;
		}

		@Override
		public Picard.Picard_Analysis call() {// acts like run
			return picard.picardASam(baseId, fullPathToSamFile, memoryRatio, altLog);
		}
	}

	private static class WorkerIndel implements Callable<GATK.IndelPrep> {
		private final GATK GATK;
		private final String fullPathToDedupReadsBam, baseId;
		private final Logger altLog;

		public WorkerIndel(GATK gATK, String baseId, String fullPathToDedupReadsBam, Logger altLog) {
			super();
			GATK = gATK;
			this.baseId = baseId;
			this.fullPathToDedupReadsBam = fullPathToDedupReadsBam;
			this.altLog = altLog;
		}

		@Override
		public GATK.IndelPrep call() {// acts like run
			return GATK.realignABam(baseId, fullPathToDedupReadsBam, altLog);
		}
	}
	
  private static class WorkerRecalibration implements Callable<GATK.BaseRecalibration> {
    private final GATK GATK;
    private final String fullPathToDedupReadsBam, baseId;
    private final Logger altLog;

    public WorkerRecalibration(GATK gATK, String baseId, String fullPathToDedupReadsBam,
                               Logger altLog) {
      super();
      GATK = gATK;
      this.baseId = baseId;
      this.fullPathToDedupReadsBam = fullPathToDedupReadsBam;
      this.altLog = altLog;
    }

    @Override
    public GATK.BaseRecalibration call() {// acts like run
      return GATK.recalibrateABam(baseId, fullPathToDedupReadsBam, altLog);
    }
  }

	private static class WorkerBamMerger implements Callable<MergeBam.BamMerger> {
		private final MergeBam mergeBam;
		private final String newBaseId, outputDir, mergeStage;
		private final String[] inputBams, baseIds;
		private final Logger altLog;

		public WorkerBamMerger(	MergeBam mergeBam, String[] baseIds, String newBaseId, String outputDir,
														String mergeStage, String[] inputBams, Logger altLog) {
			super();
			this.mergeBam = mergeBam;
			this.baseIds = baseIds;
			this.newBaseId = newBaseId;
			this.outputDir = outputDir;
			this.mergeStage = mergeStage;
			this.inputBams = inputBams;
			this.altLog = altLog;
		}

		@Override
		public MergeBam.BamMerger call() {
			return mergeBam.mergeABam(baseIds, newBaseId, inputBams, outputDir, mergeStage, altLog);
		}
	}

	public static boolean runPrep(String rootInputDir, String rootOutputDir, String fileOfSamplePairs,
																String bwaLocation, String picardLocation, String gATKLocation,
																String samtoolsLocation, String referenceGenomeFasta,
																String[] knownSitesSnpFile, String[] knownSitesIndelFile,
																boolean overwriteExisting, boolean verbose, int numSampleThreads,
																int numOtherThreads, int memoryInMB, int wallTimeInHours,
																boolean batch, int numBatches, Logger log) {
		BWA bwa = new BWA(bwaLocation, overwriteExisting, verbose, log);
		Picard picard = new Picard(picardLocation, null, overwriteExisting, verbose, log);
		GATK gatk = new GATK(	gATKLocation, referenceGenomeFasta, null, knownSitesSnpFile,
													knownSitesIndelFile, verbose, overwriteExisting, log);
		MergeBam mergeBam = new MergeBam(samtoolsLocation, overwriteExisting, verbose, log);
		GATK_LanePrep gLanePrep = new GATK_LanePrep(rootInputDir, rootOutputDir, referenceGenomeFasta,
																								verbose, numSampleThreads, numOtherThreads, bwa,
																								picard, gatk, mergeBam, log);
		if (batch) {
			gLanePrep.init(fileOfSamplePairs);
			gLanePrep.batch(numBatches, memoryInMB, wallTimeInHours, "Batch");
		} else {
			// the following are on a per lane basis

      gLanePrep.runBWA(fileOfSamplePairs);// Initializes all samples to be processed in this run of
                                          // the pipeline
      gLanePrep.runPicard();
      gLanePrep.runBaseRecal();
      gLanePrep.runBamMerge();// should skip if only one lane
      // now on to a per sample basis, if needed
      finalize(gatk, gLanePrep, verbose, log);
    }
    return true;
  }

	private static void finalize(GATK gatk, GATK_LanePrep gLanePrep, boolean verbose, Logger log) {
		String[] mergedFiles = null;// to genotype after merging
		String[] filesNotMerged = null;// to genotype directly


    if (gLanePrep.getmBamMergers() != null
        && Array.booleanArraySum(getMergeMask(gLanePrep.getmBamMergers(), log)) > 0) {
      gLanePrep.resetBwAnalysisIndividuals(gLanePrep.getmBamMergers());// step right before picard,
                                                                       // set output of merge to
                                                                       // input of picard and go
                                                                       // again
      gLanePrep.runPicard();
      if (!gLanePrep.isFail()) {
        mergedFiles = gLanePrep.getCurrentFiles(true);
        if (mergedFiles == null) {
          log.reportError("Error - could not retrieve merged filenames");
          gLanePrep.setFail(true);
        }
      }
    }
    if (!gLanePrep.isFail()) {
      String[] bamsToGenotype = null;
      filesNotMerged = gLanePrep.getFilesNotMerged();
      if (mergedFiles != null) {
        bamsToGenotype = mergedFiles;
      }
      if (filesNotMerged != null && bamsToGenotype != null) {
        bamsToGenotype = Array.concatAll(bamsToGenotype, filesNotMerged);
      } else if (filesNotMerged != null) {
        bamsToGenotype = filesNotMerged;
      }
      if (bamsToGenotype == null) {
        log.reportError("Error - could not find any files to genotype, this should not happen");
      } else {
        GATK_Genotyper gatk_Genotyper = new GATK_Genotyper(gatk, null, null, null,
                                                           gLanePrep.getNumBetweenSampleThreads(),
                                                           gLanePrep.getNumWithinSampleThreads(),
                                                           verbose, log);
        gatk_Genotyper.runSingleSampleAllSites(bamsToGenotype);
      }
    }
  }

	private static GATK.BaseRecalibration[][] getCalibrationsToMerge(	GATK.BaseRecalibration[] gRecalibrations,
																																		MergeBam mergeBam, Logger log) {
		// log.report("Warning - assuming that unique sample Ids are the first two \"_\"-delimited
		// fields of the input fastaq files, and barcodes are the third");
		Hashtable<String, ArrayList<GATK.BaseRecalibration>> track = new Hashtable<String, ArrayList<GATK.BaseRecalibration>>();
		ArrayList<String> unique = new ArrayList<String>();
		for (BaseRecalibration gRecalibration : gRecalibrations) {
			String baseId = parseBaseId(gRecalibration.getBaseId());

			if (!track.containsKey(baseId)) {
				track.put(baseId, new ArrayList<GATK.BaseRecalibration>());
				unique.add(baseId);
			}
			track.get(baseId).add(gRecalibration);
		}
		GATK.BaseRecalibration[][] calibrationsToMerge = new GATK.BaseRecalibration[unique.size()][];
		for (int i = 0; i < unique.size(); i++) {
			ArrayList<GATK.BaseRecalibration> current = track.get(unique.get(i));
			calibrationsToMerge[i] = current.toArray(new GATK.BaseRecalibration[current.size()]);
			String barcode = calibrationsToMerge[i][0].getBarcode();
			ArrayList<String> barcodesToAdd = new ArrayList<String>();
			for (int j = 0; j < calibrationsToMerge[i].length; j++) {
				if (!calibrationsToMerge[i][j].getBarcode().equals(barcode)) {
					log.report(ext.getTime()	+ " Info - since " + calibrationsToMerge[i][0].getBaseId()
											+ " and " + calibrationsToMerge[i][j].getBaseId()
											+ " appear to be the same sample with different barcodes, they will be merged");
					barcodesToAdd.add(calibrationsToMerge[i][j].getBarcode());
				}
			}
			if (barcodesToAdd.size() > 0) {// we will re-header the file here
				String barcodesAdded = Array.toStr(	Array.unique(barcodesToAdd.toArray(new String[barcodesToAdd.size()])),
																						FileNameParser.SPLIT);
				String newSampleId = calibrationsToMerge[i][0].getBaseId()	+ FileNameParser.SPLIT
															+ barcodesAdded;
				for (int j = 0; j < calibrationsToMerge[i].length; j++) {
					calibrationsToMerge[i][j].setNewBaseId(newSampleId);
					// ReHeader reHeader =
					// mergeBam.reHeaderBamFilePriorToMerge(calibrationsToMerge[i][j].getRrd_bam(),
					// calibrationsToMerge[i][j].getBaseId(), newSampleId, log);
					// if (!reHeader.isFail()) {
					// // calibrationsToMerge[i][j].setRrd_bam(reHeader.getReHeaderBam());
					// } else {
					// log.report("Error - could not re header file " +
					// calibrationsToMerge[i][j].getRrd_bam());
					// calibrationsToMerge[i][j].setFail(reHeader.isFail());
					// }
				}
				log.report(ext.getTime() + " Info - new ID is  " + calibrationsToMerge[i][0].getBaseId());
			} else {
				for (int j = 0; j < calibrationsToMerge[i].length; j++) {
					calibrationsToMerge[i][j].setNewBaseId(calibrationsToMerge[i][j].getBaseId());// no
																																												// modification
				}
			}
		}

		return calibrationsToMerge;
	}

	private static String parseBaseId(String baseId) {
		int len = baseId.split(FileNameParser.SPLIT).length;
		if (len < 3) {
			return baseId;
		} else {
			return Array.toStr(Array.subArray(baseId.split(FileNameParser.SPLIT), 0, 2));
		}
	}

	private static String[] getInputFilesFrom(GATK.BaseRecalibration[] gRecalibrations) {
		String[] inputBams = new String[gRecalibrations.length];
		for (int i = 0; i < gRecalibrations.length; i++) {
			inputBams[i] = gRecalibrations[i].getRrd_bam();
		}
		return inputBams;
	}

	private static boolean[] getMergeMask(MergeBam.BamMerger[] mBamMergers, Logger log) {
		boolean[] mergeMask = new boolean[mBamMergers.length];
		for (int i = 0; i < mergeMask.length; i++) {
			mergeMask[i] = mBamMergers[i].shouldMerge();
		}
		return mergeMask;
	}

	public static void main(String[] args) {

		int numArgs = args.length;
		String rootInputDir = null;
		String rootOutputDir = null;
		String referenceGenomeFasta = "";
		String bwaLocation = "";
		String picardLocation = "";
		String gATKLocation = "";
		String samtoolsLocation = "";

		String[] knownSitesSnpFile = new String[] {"NA"};
		String[] knownSitesIndelFile = new String[] {"NA"};

		String fileOfSamplePairs = null;
		boolean verbose = true;
		boolean batch = false;
		int numBatches = 5;
		int memoryInMB = 23000;
		int wallTimeInHours = 24;
		int numWithinSampleThreads = 1;
		int numBetweenSampleThreads = 1;
		boolean overwriteExisting = false;

		String logFile = "GATK_PREP.log";


    String usage = "\n" + "seq.BWA_Analysis requires 2 argument\n";
    usage += "   (1) root input directory (i.e. " + ROOT_INPUT_COMMAND + rootInputDir
             + " (no default))\n" + "";
    usage += "   (2) root output directory (i.e. " + ROOT_OUTPUT_COMMAND + rootOutputDir
             + " (no default))\n" + "";
    usage += "   (3) tab-delimited file with no header of paired .fastq (i.e. "
             + FILE_OF_SAMPLE_PAIRS_COMMAND + fileOfSamplePairs + " (optional, no default))\n" + "";
    usage += "   (4) the full path to a  reference genome in fasta format (i.e."
             + REFERENCE_GENOME_COMMAND + referenceGenomeFasta + " (no default))\n" + "";
    usage += "   (5) the full path to the bwa executable (i.e. " + BWA_LOCATION_COMMAND
             + bwaLocation + " (defualts to systems path))\n" + "";
    usage += "   (6) the full path to the picard directory containing the .jar (i.e. "
             + Picard.PICARD_LOCATION_COMMAND + picardLocation + " (default))\n" + "";
    usage += "   (7) the full path to the GATK executable (i.e. " + GATK.GATK_LOCATION_COMMAND
             + gATKLocation + " (defualts to systems path))\n" + "";
    usage += "   (7) the full path to the samtools directory containing .jar files (i.e. "
             + MergeBam.SAMTOOLS_COMMAND + samtoolsLocation + " (defualts to systems path))\n" + "";
    usage += "   (8) the full path to reference indel files (comma delimited if multiple) (i.e. "
             + GATK.KNOWN_SITES_INDEL_LOCATION_COMMAND + Array.toStr(knownSitesIndelFile, ",")
             + " (default))\n" + "";
    usage += "   (9) the full path to reference snp files (comma delimited if multiple) (i.e. "
             + GATK.KNOWN_SITES_SNP_LOCATION_COMMAND + Array.toStr(knownSitesSnpFile, ",")
             + " (default))\n" + "";

		usage += "   (10) run in quiet mode (i.e. " + QUIET_COMMAND + " (not tbe default))\n" + "";
		usage += "   (11) number of threads per sample for bwa mem (i.e."	+ NUM_BETWEEN_THREADS_COMMAND
							+ numBetweenSampleThreads + " (default))\n" + "";
		usage += "   (12) number of sample threads for bwa mem (i.e."	+ NUM_WITHIN_THREADS_COMMAND
							+ numWithinSampleThreads + " (default))\n" + "";

		usage +=
					"   (13) filename for a log (i.e. " + LOG_FILE_COMMAND + logFile + " (default))\n" + "";
		usage += "   (14) set up a batch analysis for the root input directory for a log (i.e. "
							+ BATCH_COMMAND + " (not the default))\n" + "";
		usage += "   (15) number of batches for a batched analysis (i.e. "	+ NUMBATCHES_COMMAND
							+ numBatches + " (the default))\n" + "";
		usage += "   (16) over-write exsiting files (i.e. "	+ OVERWRITE_EXISTING_COMMAND
							+ " (not the default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith(ROOT_INPUT_COMMAND)) {
				rootInputDir = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(ROOT_OUTPUT_COMMAND)) {
				rootOutputDir = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(FILE_OF_SAMPLE_PAIRS_COMMAND)) {
				fileOfSamplePairs = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(REFERENCE_GENOME_COMMAND)) {
				referenceGenomeFasta = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(BWA_LOCATION_COMMAND)) {
				bwaLocation = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(MergeBam.SAMTOOLS_COMMAND)) {
				samtoolsLocation = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(LOG_FILE_COMMAND)) {
				logFile = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(NUM_BETWEEN_THREADS_COMMAND)) {
				numBetweenSampleThreads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(NUM_WITHIN_THREADS_COMMAND)) {
				numWithinSampleThreads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(NUMBATCHES_COMMAND)) {
				numBatches = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("memoryInMB=")) {
				memoryInMB = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("wallTimeInHours=")) {
				wallTimeInHours = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(QUIET_COMMAND)) {
				verbose = false;
				numArgs--;
			} else if (arg.startsWith(BATCH_COMMAND)) {
				batch = true;
				numArgs--;
			} else if (arg.startsWith(OVERWRITE_EXISTING_COMMAND)) {
				batch = true;
				numArgs--;
			} else if (arg.startsWith(Picard.PICARD_LOCATION_COMMAND)) {
				picardLocation = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(GATK.GATK_LOCATION_COMMAND)) {
				gATKLocation = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(GATK.KNOWN_SITES_SNP_LOCATION_COMMAND)) {
				knownSitesSnpFile = ext.parseStringArg(arg, "").split(GATK.SPLIT);
				numArgs--;
			} else if (arg.startsWith(GATK.KNOWN_SITES_INDEL_LOCATION_COMMAND)) {
				knownSitesIndelFile = ext.parseStringArg(arg, "").split(GATK.SPLIT);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}

		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		Logger log = new Logger((rootOutputDir == null ? rootInputDir : rootOutputDir)
														+ "GATK_PREP.log");
		runPrep(rootInputDir, rootOutputDir, fileOfSamplePairs, bwaLocation, picardLocation,
						gATKLocation, samtoolsLocation, referenceGenomeFasta, knownSitesSnpFile,
						knownSitesIndelFile, overwriteExisting, verbose, numWithinSampleThreads,
						numBetweenSampleThreads, memoryInMB, wallTimeInHours, batch, numBatches, log);
	}
}
