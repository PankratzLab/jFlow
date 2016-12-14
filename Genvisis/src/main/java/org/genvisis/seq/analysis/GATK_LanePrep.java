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

/**
 * Prepping for the GATK is done on a lane by lane basis as reflected here
 *
 */
public class GATK_LanePrep extends BWA_Analysis {
	private static final String PICARD_METRICS_SUMMARY = "picard_metrics_summary.txt";
	private final Picard picard;
  private Picard.PicardAnalysis[] picardAnalyses;
  private Picard.PicardMergeDedupe[] picardMergeDedupes;
	private final GATK gatk;
	private GATK.BaseRecalibration[] gRecalibrations;

	public GATK_LanePrep(	String rootInputDir, String rootOutputDir, String referenceGenomeFasta,
												boolean verbose, int numWithinSampleThreads, int numBetweenSampleThreads,
                       BWA bwa, Picard picard, GATK gatk, Logger log) {
		super(rootInputDir, (rootOutputDir == null ? rootInputDir : rootOutputDir),
					referenceGenomeFasta, verbose, numWithinSampleThreads, numBetweenSampleThreads, bwa, log);
		this.picard = picard;
		this.gatk = gatk;
	}

	public void runBWA(String fileOfSamplePairs) {
		init(fileOfSamplePairs);
		if (!isFail()) {
			analyzeBWA_MEM();
		}
	}

	public void runPicard() {
		if (!isFail()) {
			BWA_AnalysisIndividual[] bwAnalysisIndividuals = getBwAnalysisIndividuals();
			if (bwAnalysisIndividuals != null) {
        picardAnalyses = new Picard.PicardAnalysis[bwAnalysisIndividuals.length];
        double memoryRatio = calcMemoryRatio(bwAnalysisIndividuals.length, getLog());
				ExecutorService executor = Executors.newFixedThreadPool(getNumBetweenSampleThreads());
				Hashtable<String, Future<Picard.PicardAnalysis>> tmpResults =
                                                                    new Hashtable<String, Future<Picard.PicardAnalysis>>();
				for (int i = 0; i < bwAnalysisIndividuals.length; i++) {
					Logger altLog = new Logger(ext.rootOf(getLog().getFilename(), false)	+ "_Picard_ID_"
																			+ bwAnalysisIndividuals[i].getID() + "_Lane_"
																			+ bwAnalysisIndividuals[i].getLane() + ".log");
					tmpResults.put(Integer.toString(i),
                         executor.submit(new WorkerPicard(picard, bwAnalysisIndividuals[i].getID(),
                                                          bwAnalysisIndividuals[i].getOutput(),
                                                          memoryRatio, altLog)));
				}
				for (int i = 0; i < bwAnalysisIndividuals.length; i++) {
					try {
						getLog().memoryPercentFree();
            picardAnalyses[i] = tmpResults.get(Integer.toString(i)).get();
						getLog().report(ext.getTime()	+ "Info - retrieving picard results for "
                            + picardAnalyses[i].getFullPathToSamFile());
            if (picardAnalyses[i].isFail() && !isFail()) {
              getLog().reportError("Failed picard for "
                                       + picardAnalyses[i].getFullPathToSamFile());
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
        String[] picardFiles = new String[picardAnalyses.length];
        for (int i = 0; i < picardAnalyses.length; i++) {
          if (picardAnalyses[i].isAllThere()) {
            picardFiles[i] = picardAnalyses[i].getFullPathToMetricsTxt();
					} else {
            getLog().reportError("Could not find picard metrics file for root input files:\n"
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
      getLog().reportError("Should not be requesting other files with un implemented request");
      return null;
    } else {
      String[] currentFiles = new String[picardAnalyses.length];
      for (int i = 0; i < currentFiles.length; i++) {
        currentFiles[i] = picardAnalyses[i].getFullPathToSortedDeDuppedBamFile();
      }
      return currentFiles;
    }
  }

  private void runMergeDedupe() {
    if (!isFail()) {
      if (picardAnalyses != null) {
        Picard.PicardAnalysis[][] picardAnalysesToMerge = getPicardAnalysesToMerge(picardAnalyses,
                                                                                   getLog());
        picardMergeDedupes = new Picard.PicardMergeDedupe[picardAnalysesToMerge.length];
        double memoryRatio = calcMemoryRatio(picardAnalysesToMerge.length, getLog());
        ExecutorService executor = Executors.newFixedThreadPool(getNumBetweenSampleThreads());
        Hashtable<String, Future<Picard.PicardMergeDedupe>> tmpResults =
                                                                       new Hashtable<String, Future<Picard.PicardMergeDedupe>>();
        for (int i = 0; i < picardMergeDedupes.length; i++) {
          Logger altLog = new Logger(ext.rootOf(getLog().getFilename(), false)
                                     + "_Picard.PicardMergeDedupe_ID_"
                                     + picardAnalysesToMerge[i][0].getBaseID() + ".log");
          String[] inputBams = new String[picardAnalysesToMerge[i].length];
          String[] inputBamIndices = new String[inputBams.length];
          for (int j = 0; j < picardAnalysesToMerge[i].length; j++) {
            inputBams[j] = picardAnalysesToMerge[i][j].getFullPathToSortedDeDuppedBamFile();
            inputBamIndices[j] =
                               picardAnalysesToMerge[i][j].getFullPathToSortedDeDuppedBamFileIndex();
          }
          tmpResults.put(Integer.toString(i),
                         executor.submit(new WorkerMergeDedupe(picard,
                                                               picardAnalysesToMerge[i][0].getBaseID(),
                                                               inputBams, inputBamIndices,
                                                               getRootOutputDir(), memoryRatio,
                                                               altLog)));
        }
        for (int i = 0; i < picardMergeDedupes.length; i++) {
          try {
            picardMergeDedupes[i] = tmpResults.get(Integer.toString(i)).get();
            if ((picardMergeDedupes[i].isFail()) && !isFail()) {
              getLog().reportError("Failed merging and deduping for "
                                       + Array.toStr(picardMergeDedupes[i].getFullPathsToInputBams(),
                                                     "\n"));
              setFail(true);
            }
          } catch (InterruptedException e) {
            getLog().reportError("When running GATK Merge/Dedupe on internal index " + i);
            getLog().reportException(e);
            setFail(true);
          } catch (ExecutionException e) {
            getLog().reportError("When running GATK Merge/Dedupe on internal index " + i);
            getLog().reportException(e);
            setFail(true);
          }
        }
        executor.shutdown();
        BWA_AnalysisIndividual[] bwAnalysisIndividuals =
                                                       new BWA_AnalysisIndividual[picardMergeDedupes.length];
        for (int i = 0; i < picardMergeDedupes.length; i++) {
          bwAnalysisIndividuals[i] = new BWA_AnalysisIndividual(null, getRootOutputDir(),
                                                                picardMergeDedupes[i].getBaseID(),
                                                                null, null, null, null, getLog());
          bwAnalysisIndividuals[i].setOutput(picardMergeDedupes[i].getFullPathToMergedDedupedBam());
        }
        setBwAnalysisIndividuals(bwAnalysisIndividuals);
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

  private void runBaseRecal() {
    if (!isFail()) {
      if (picardMergeDedupes != null) {
        gRecalibrations = new GATK.BaseRecalibration[picardMergeDedupes.length];
        ExecutorService executor = Executors.newFixedThreadPool(getNumBetweenSampleThreads());
        Hashtable<String, Future<GATK.BaseRecalibration>> tmpResults =
                                                                     new Hashtable<String, Future<GATK.BaseRecalibration>>();
        for (int i = 0; i < picardMergeDedupes.length; i++) {
          Logger altLog = new Logger(ext.rootOf(getLog().getFilename(), false)
                                     + "_BaseRecalibration_ID_"
                                     + getBwAnalysisIndividuals()[i].getID() + ".log");
          tmpResults.put(Integer.toString(i),
                         executor.submit(new WorkerRecalibration(gatk,
                                                                 picardMergeDedupes[i].getBaseID(),
                                                                 picardMergeDedupes[i].getFullPathToMergedDedupedBam(),
                                                                 altLog)));
        }
        for (int i = 0; i < picardMergeDedupes.length; i++) {
          try {
            gRecalibrations[i] = tmpResults.get(Integer.toString(i)).get();
            if (gRecalibrations[i].isFail() && !isFail()) {
              getLog().reportError("Failed recalibration for "
                                       + gRecalibrations[i].getDedup_reads_bam());
              setFail(true);
            }
          } catch (InterruptedException e) {
            getLog().reportError("When running GATK Base recalibration on internal index "
                                 + i);
            getLog().reportException(e);
            setFail(true);
          } catch (ExecutionException e) {
            getLog().reportError("When running GATK Base recalibration on internal index "
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

  private void genotype() {
    GATK_Genotyper gatk_Genotyper = new GATK_Genotyper(gatk, null, null, null,
                                                       getNumBetweenSampleThreads(),
                                                       getNumWithinSampleThreads(), isVerbose(),
                                                       getLog());
    String[] genotypeBams = new String[gRecalibrations.length];
    for (int i = 0; i < genotypeBams.length; i++) {
      genotypeBams[i] = gRecalibrations[i].getRrd_bam();
    }
    gatk_Genotyper.runSingleSampleAllSites(genotypeBams);
  }

	public Picard getPicard() {
		return picard;
	}

	public GATK getGatk() {
		return gatk;
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
		String command = Array.toStr(PSF.Load.getAllModules(), "\n");
    command += "\njava -Xmx" + memoryInMB + "m -jar ~/genvisisGATK3.6.jar " + this.getClass().getName()
								+ ROOT_INPUT_COMMAND + getRootInputDir() + SPACE + ROOT_OUTPUT_COMMAND
								+ getRootOutputDir() + SPACE;
		command += REFERENCE_GENOME_COMMAND	+ getReferenceGenomeFasta() + SPACE + BWA_LOCATION_COMMAND
								+ getBwa().getBwaLocation() + SPACE;
		command += NUM_BETWEEN_THREADS_COMMAND	+ getNumBetweenSampleThreads() + SPACE
								+ FILE_OF_SAMPLE_PAIRS_COMMAND + getRootOutputDir() + "[%0].txt" + SPACE
								+ NUM_WITHIN_THREADS_COMMAND + getNumWithinSampleThreads() + SPACE;
		command += Picard.PICARD_LOCATION_COMMAND + getPicard().getPicardLocation() + SPACE;
		command += GATK.GATK_LOCATION_COMMAND + getGatk().getGATKLocation() + SPACE;
		command += GATK.KNOWN_SITES_SNP_LOCATION_COMMAND
								+ Array.toStr(getGatk().getKnownSitesSnpFile(), GATK.SPLIT) + SPACE;
		command += GATK.KNOWN_SITES_INDEL_LOCATION_COMMAND
								+ Array.toStr(getGatk().getKnownSitesIndelFile(), GATK.SPLIT);
		Files.qsub("GATK_Lane_Prep"	+ baseName, command, batches, memoryInMB, wallTimeInHours,
								getNumWithinSampleThreads() * getNumBetweenSampleThreads());
	}

  private static class WorkerPicard implements Callable<Picard.PicardAnalysis> {
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
    public Picard.PicardAnalysis call() {// acts like run
			return picard.picardASam(baseId, fullPathToSamFile, memoryRatio, altLog);
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

  private static class WorkerMergeDedupe implements Callable<Picard.PicardMergeDedupe> {
    private final Picard picard;
    private final String[] fullPathsToInputBams;
    private final String[] fullPathsToInputBamIndices;
    private final String outputDir;
    private final String baseID;
		private final Logger altLog;
    private final double memoryRatio;

    public WorkerMergeDedupe(Picard picard, String baseId, String[] fullPathsToInputBams,
                             String[] fullPathsToInputBamIndices, String outputDir,
                             double memoryRatio, Logger altLog) {
			super();
      this.picard = picard;
      this.baseID = baseId;
      this.fullPathsToInputBams = fullPathsToInputBams;
      this.fullPathsToInputBamIndices = fullPathsToInputBamIndices;
			this.outputDir = outputDir;
      this.memoryRatio = memoryRatio;
			this.altLog = altLog;
		}

		@Override
    public Picard.PicardMergeDedupe call() {// acts like run
      return picard.mergeDedupeBams(baseID, fullPathsToInputBams, fullPathsToInputBamIndices,
                                    outputDir, memoryRatio, altLog);
		}
	}

	public static boolean runPrep(String rootInputDir, String rootOutputDir, String fileOfSamplePairs,
																String bwaLocation, String picardLocation, String gATKLocation,
                                String referenceGenomeFasta, String[] knownSitesSnpFile,
                                String[] knownSitesIndelFile, boolean overwriteExisting,
                                boolean verbose, int numSampleThreads, int numOtherThreads,
                                int memoryInMB, int wallTimeInHours, boolean batch, int numBatches,
                                Logger log) {
		BWA bwa = new BWA(bwaLocation, overwriteExisting, verbose, log);
		Picard picard = new Picard(picardLocation, null, overwriteExisting, verbose, log);
		GATK gatk = new GATK(	gATKLocation, referenceGenomeFasta, null, knownSitesSnpFile,
													knownSitesIndelFile, verbose, overwriteExisting, log);
		GATK_LanePrep gLanePrep = new GATK_LanePrep(rootInputDir, rootOutputDir, referenceGenomeFasta,
																								verbose, numSampleThreads, numOtherThreads, bwa,
                                                picard, gatk, log);
		if (batch) {
			gLanePrep.init(fileOfSamplePairs);
			gLanePrep.batch(numBatches, memoryInMB, wallTimeInHours, "Batch");
		} else {
			// the following are on a per lane basis

      gLanePrep.runBWA(fileOfSamplePairs);// Initializes all samples to be processed in this run of
                                          // the pipeline
      gLanePrep.runPicard();

      gLanePrep.runMergeDedupe();
      gLanePrep.runBaseRecal();
      // gLanePrep.runBamMerge();// should skip if only one lane
      // now on to a per sample basis, if needed
      gLanePrep.genotype();
    }
    return true;
  }


  private static double calcMemoryRatio(int numSamples, Logger log) {
    double memoryRatio = (double) 1 / numSamples;
    memoryRatio -= .01;
    if (memoryRatio > Picard.DEFAULT_SORTING_COLLECTION_SIZE_RATIO) {
      memoryRatio = Picard.DEFAULT_SORTING_COLLECTION_SIZE_RATIO;
      } else {
      log.report(ext.getTime() + " Info - adjusting Picard's memory ratio to " + memoryRatio
                 + " since there are more than 3 samples...");
      }
    return memoryRatio;
    }

  private static Picard.PicardAnalysis[][] getPicardAnalysesToMerge(Picard.PicardAnalysis[] picardAnalyses,
                                                                    Logger log) {
		// log.report("Warning - assuming that unique sample Ids are the first two \"_\"-delimited
		// fields of the input fastaq files, and barcodes are the third");
    Hashtable<String, ArrayList<Picard.PicardAnalysis>> track =
                                                              new Hashtable<String, ArrayList<Picard.PicardAnalysis>>();
		ArrayList<String> unique = new ArrayList<String>();
    for (Picard.PicardAnalysis picardAnalysis : picardAnalyses) {
      String baseId = parseBaseId(picardAnalysis.getBaseID());

			if (!track.containsKey(baseId)) {
        track.put(baseId, new ArrayList<Picard.PicardAnalysis>());
				unique.add(baseId);
			}
      track.get(baseId).add(picardAnalysis);
		}
    Picard.PicardAnalysis[][] analysesToMerge = new Picard.PicardAnalysis[unique.size()][];
		for (int i = 0; i < unique.size(); i++) {
      ArrayList<Picard.PicardAnalysis> current = track.get(unique.get(i));
      analysesToMerge[i] = current.toArray(new Picard.PicardAnalysis[current.size()]);
      String barcode = analysesToMerge[i][0].getBarcode();
			ArrayList<String> barcodesToAdd = new ArrayList<String>();
      for (int j = 0; j < analysesToMerge[i].length; j++) {
        if (!analysesToMerge[i][j].getBarcode().equals(barcode)) {
          log.report(ext.getTime() + " Info - since " + analysesToMerge[i][0].getBaseID() + " and "
                     + analysesToMerge[i][j].getBaseID()
											+ " appear to be the same sample with different barcodes, they will be merged");
          barcodesToAdd.add(analysesToMerge[i][j].getBarcode());
				}
			}
			if (barcodesToAdd.size() > 0) {// we will re-header the file here
				String barcodesAdded = Array.toStr(	Array.unique(barcodesToAdd.toArray(new String[barcodesToAdd.size()])),
																						FileNameParser.SPLIT);
        String newSampleId =
                           analysesToMerge[i][0].getBaseID() + FileNameParser.SPLIT + barcodesAdded;
        for (int j = 0; j < analysesToMerge[i].length; j++) {
          analysesToMerge[i][j].setNewBaseID(newSampleId);
					// ReHeader reHeader =
					// mergeBam.reHeaderBamFilePriorToMerge(calibrationsToMerge[i][j].getRrd_bam(),
					// calibrationsToMerge[i][j].getBaseId(), newSampleId, log);
					// if (!reHeader.isFail()) {
					// // calibrationsToMerge[i][j].setRrd_bam(reHeader.getReHeaderBam());
					// } else {
					// log.report("could not re header file " +
					// calibrationsToMerge[i][j].getRrd_bam());
					// calibrationsToMerge[i][j].setFail(reHeader.isFail());
					// }
				}
        log.report(ext.getTime() + " Info - new ID is  " + analysesToMerge[i][0].getBaseID());
			} else {
        for (int j = 0; j < analysesToMerge[i].length; j++) {
          analysesToMerge[i][j].setNewBaseID(analysesToMerge[i][j].getBaseID());// no modification
				}
			}
		}

    return analysesToMerge;
	}

	private static String parseBaseId(String baseId) {
		int len = baseId.split(FileNameParser.SPLIT).length;
		if (len < 3) {
			return baseId;
		} else {
			return Array.toStr(Array.subArray(baseId.split(FileNameParser.SPLIT), 0, 2));
		}
	}

	public static void main(String[] args) {

		int numArgs = args.length;
		String rootInputDir = null;
		String rootOutputDir = null;
		String referenceGenomeFasta = "";
		String bwaLocation = "";
		String picardLocation = "";
		String gATKLocation = "";

		String[] knownSitesSnpFile = new String[] {"NA"};
		String[] knownSitesIndelFile = new String[] {"NA"};

		String fileOfSamplePairs = null;
		boolean verbose = true;
		boolean batch = false;
		int numBatches = 5;
    int memoryInMB = 61000;
    int wallTimeInHours = 96;
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
    usage += "   (6) the full path to the picard (2.6.0) directory containing the .jar (i.e. "
             + Picard.PICARD_LOCATION_COMMAND + picardLocation + " (default))\n" + "";
    usage += "   (7) the full path to the GATK (3.6) executable (i.e. " + GATK.GATK_LOCATION_COMMAND
             + gATKLocation + " (defualts to systems path))\n" + "";
    usage += "   (8) the full path to reference indel files (comma delimited if multiple) (i.e. "
             + GATK.KNOWN_SITES_INDEL_LOCATION_COMMAND + Array.toStr(knownSitesIndelFile, ",")
             + " (default))\n" + "";
    usage += "   (9) the full path to reference snp files (comma delimited if multiple) (i.e. "
             + GATK.KNOWN_SITES_SNP_LOCATION_COMMAND 
             + Array.toStr(new String[] {"site1", "site2", "site3"}, ",") + " (not the default))\n" + "";

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
            gATKLocation, referenceGenomeFasta, knownSitesSnpFile, knownSitesIndelFile,
            overwriteExisting, verbose, numWithinSampleThreads, numBetweenSampleThreads, memoryInMB,
            wallTimeInHours, batch, numBatches, log);
	}
}
