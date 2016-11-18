package org.genvisis.cnv.hmm;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.concurrent.Callable;

import org.genvisis.cnv.analysis.PennCNV;
import org.genvisis.cnv.filesys.Centroids;
import org.genvisis.cnv.filesys.MarkerSet.PreparedMarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.hmm.PennHmm.ViterbiResult;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.Resources.Resource;
import org.genvisis.cnv.qc.GcAdjustor;
import org.genvisis.cnv.qc.GcAdjustor.GCAdjustorBuilder;
import org.genvisis.cnv.qc.GcAdjustor.GC_CORRECTION_METHOD;
import org.genvisis.cnv.qc.GcAdjustor.GcModel;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.Positions;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.LocusSet;
import org.genvisis.seq.manage.BamImport.NGS_MARKER_TYPE;

import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

/**
 * @author lane0212
 *
 *         Handles the data preparation and finalizations for calling cnvs via the PennCNV methods
 */
public class CNVCaller {
	public static final double MIN_LRR_MEDIAN_ADJUST = -2;
	public static final double MAX_LRR_MEDIAN_ADJUST = 2;
	private static final double MIN_BAF_MEDIAN_ADJUST = .25;
	private static final double MAX_BAF_MEDIAN_ADJUST = .75;
	private static final int PENN_CNV_SIG_FIGS = 4;
	public static final int DEFUALT_MIN_SITES = 3;
	public static final int DEFUALT_MIN_CONF = 3;

	private static final int MIN_MARKERS_PER_CHROMOSOME = 10;

	private final Project proj;
	private final String dna;
	private PennHmm pennHmm;
	private final GcModel gcModel;
	private final PreparedMarkerSet markerSet;
	private boolean[] markersToUse;
	private boolean[] copyNumberDef;
	private double[] analysisLrrs, analysisBafs, analysisPfbs;
	private int[] analysisProjectIndices, analysisPositions;// linker of what is being analyzed to the
																													// order in the project
	private final DATA_ADJUSTMENTS[] dataAdjustments;
	private final PFB_MANAGEMENT_TYPE pManagementType;
	private final boolean debugMode;

	/**
	 * A couple ways that a PFB file can be used;
	 *
	 */
	public enum PFB_MANAGEMENT_TYPE {
																		/**
																		 * Less than 0 pfbs are removed from the analysis
																		 */
																		PENNCNV_DEFAULT,
																		/**
																		 * Less than 0 pfbs are treated as CN only. In genvisis -1 pfbs
																		 * represent snps without genotypes (uncalled, see
																		 * {@link PennCNV#populationBAF(Project)})
																		 */
																		LESS_THAN_0_GO_CN;
	}

	/**
	 * @param proj
	 * @param dna the samples dna ID
	 * @param pennHmm
	 * @param gcModel
	 * @param pfb
	 * @param dataAdjustments the order of operations for adjusting data (such as whether to scale to
	 *        median and then gc correct, or the other way around)
	 * @param markerSet
	 * @param markersToUse the markers to use to make CNV calls
	 * @param copyNumberDef boolean array defining copy number only (non snp) markers... which get
	 *        their own hmm params
	 * @param lrrs intense!city
	 * @param bafs
	 * @param debugMode
	 */
	public CNVCaller(	Project proj, String dna, PennHmm pennHmm, GcModel gcModel, PFB pfb,
										DATA_ADJUSTMENTS[] dataAdjustments, PreparedMarkerSet markerSet,
										boolean[] markersToUse, boolean[] copyNumberDef, double[] lrrs, double[] bafs,
										PFB_MANAGEMENT_TYPE pManagementType, boolean debugMode) {
		super();
		this.proj = proj;
		this.pennHmm = pennHmm;
		this.dna = dna;
		this.gcModel = gcModel;
		this.markerSet = markerSet;
		this.markersToUse = markersToUse;
		analysisLrrs = lrrs;
		analysisBafs = bafs;
		analysisPfbs = pfb.getPfbs();
		this.pManagementType = pManagementType;
		managePfbs();
		this.dataAdjustments = dataAdjustments;
		this.copyNumberDef = copyNumberDef;
		analysisProjectIndices = Array.arrayOfIndices(markerSet.getMarkerNames().length);
		analysisPositions = markerSet.getPositions();
		this.debugMode = debugMode;

		if (lrrs.length != bafs.length || markerSet.getMarkerNames().length != lrrs.length) {
			String error = "BUG: must supply entire lrr and baf data set for the project, consider using the markersToUse array to subset the analysis";
			proj.getLog().reportError(error);
			throw new IllegalArgumentException(error);
		}
		if (debugMode) {
			proj.getLog().reportTimeInfo("Sample: " + dna);
		}
	}

	private void managePfbs() {
		if (markerSet.getMarkerNames().length != analysisPfbs.length) {
			throw new IllegalArgumentException("Method seems to be called at the wrong time, call before subsetting");
		}
		switch (pManagementType) {
			case LESS_THAN_0_GO_CN:
				for (int i = 0; i < analysisPfbs.length; i++) {
					if (analysisPfbs[i] < 1) {
						analysisPfbs[i] = 2;// The switch for CN only analysis
					}
				}
				break;
			case PENNCNV_DEFAULT:
				int numRemoved = 0;
				if (markersToUse == null) {
					markersToUse = Array.booleanArray(markerSet.getMarkerNames().length, true);
				}
				for (int i = 0; i < analysisPfbs.length; i++) {
					if (analysisPfbs[i] < 0) {
						markersToUse[i] = false;
						numRemoved++;
					}
				}
				if (numRemoved > 0) {
					proj.getLog()
							.reportTimeInfo(numRemoved
															+ " markers were removed, defined by pfb having a value less than 0");
				}
				break;
			default:
				throw new IllegalArgumentException("Invalid management type " + pManagementType);
		}
	}

	private static double[] roundToPennCNVSigFigs(double[] array) {
		return Array.round(array, PENN_CNV_SIG_FIGS);
	}

	/**
	 * Warning: only the default PennCNV order of adjustments has been tested
	 */
	public void adjustData() {

		if (dataAdjustments == null) {
			proj.getLog().reportTimeWarning("No data adjustments supplied");
		} else {
			double lrrSd = Array.stdev(Array.getValuesBetween(analysisLrrs, MIN_LRR_MEDIAN_ADJUST,
																												MAX_LRR_MEDIAN_ADJUST));
			for (DATA_ADJUSTMENTS type : dataAdjustments) {
				switch (type) {
					case ROUND_TO_PENNCNV_SIG_FIGS:// May not be necessary
						analysisLrrs = roundToPennCNVSigFigs(analysisLrrs);
						analysisBafs = roundToPennCNVSigFigs(analysisBafs);
						analysisPfbs = roundToPennCNVSigFigs(analysisPfbs);
						// lrrSd = Array.stdev(getValuesBetween(analysisLrrs, MIN_LRR_MEDIAN_ADJUST,
						// MAX_LRR_MEDIAN_ADJUST));
						break;
					case ADJUST_HMM_SD:
						pennHmm = PennHmm.adjustBSD(pennHmm, lrrSd, proj.getLog());
						break;
					case GC_ADJUST:
						if (gcModel == null) {
							String error = "gc model cannot be null if adjustment " + type + " is flagged";
							proj.getLog().reportError(error);
							throw new IllegalArgumentException();
						} else {
							double[] dataToCorrect = analysisLrrs;
							if (analysisProjectIndices.length != markerSet.getMarkerNames().length) {// only
																																												// current
																																												// indicies
																																												// will be
																																												// used
								dataToCorrect = new double[markerSet.getMarkerNames().length];
								for (int j = 0; j < analysisProjectIndices.length; j++) {
									dataToCorrect[analysisProjectIndices[j]] = analysisLrrs[j];
								}
							}
							GCAdjustorBuilder builder = new GCAdjustorBuilder();
							builder.correctionMethod(GC_CORRECTION_METHOD.PENNCNV_GC);
							builder.verbose(debugMode);
							GcAdjustor gcAdjustor = builder.build(proj, markerSet, gcModel, dataToCorrect);
							gcAdjustor.correctIntensities();
							gcAdjustor.computeQCMetrics(true, true);
							analysisLrrs = Array.subArray(gcAdjustor.getCorrectedIntensities(),
																						analysisProjectIndices);
							proj.getLog().reportTimeInfo(gcAdjustor.getAnnotatedQCString());
							lrrSd = Array.stdev(Array.getValuesBetween(	analysisLrrs, MIN_LRR_MEDIAN_ADJUST,
																													MAX_LRR_MEDIAN_ADJUST));

						}
						break;
					case HANDLE_NAN:
						int nanCount = 0;
						for (int j = 0; j < analysisBafs.length; j++) {
							if (Double.isNaN(analysisBafs[j]) || Double.isNaN(analysisLrrs[j])) {
								analysisBafs[j] = 0;
								analysisLrrs[j] = 0;
								nanCount++;
							}
						}
						if (nanCount > 0) {
							proj.getLog().reportTimeInfo("Set "	+ nanCount
																						+ " observations with NaN data to 0 for lrr and baf");
						}
						lrrSd = Array.stdev(Array.getValuesBetween(	analysisLrrs, MIN_LRR_MEDIAN_ADJUST,
																												MAX_LRR_MEDIAN_ADJUST));

						break;
					case MEDIAN_ADJUST:
						analysisLrrs = adjustLrr(	analysisLrrs, MIN_LRR_MEDIAN_ADJUST, MAX_LRR_MEDIAN_ADJUST,
																			debugMode, proj.getLog());
						analysisBafs = adjustBaf(	analysisBafs, MIN_BAF_MEDIAN_ADJUST, MAX_BAF_MEDIAN_ADJUST,
																			debugMode, proj.getLog());
						// TODO, update lrrSd later?
						// lrrSd = Array.stdev(getValuesBetween(analysisLrrs, MIN_LRR_MEDIAN_ADJUST,
						// MAX_LRR_MEDIAN_ADJUST)); PennCNV does not update lrr sd here so we wont either
						break;
					case SUBSET_TO_ANALYSIS_MARKERS:
						if (markersToUse != null) {
							boolean[] tmpExclude = markersToUse;
							if (analysisProjectIndices.length != markerSet.getMarkerNames().length) {
								tmpExclude = Array.subArray(markersToUse, analysisProjectIndices);
							}
							analysisLrrs = Array.subArray(analysisLrrs, tmpExclude);
							analysisBafs = Array.subArray(analysisBafs, tmpExclude);
							analysisPfbs = Array.subArray(analysisPfbs, tmpExclude);
							analysisPositions = Array.subArray(analysisPositions, tmpExclude);
							analysisProjectIndices = Array.subArray(analysisProjectIndices, tmpExclude);
							copyNumberDef = Array.subArray(copyNumberDef, analysisProjectIndices);
							lrrSd = Array.stdev(Array.getValuesBetween(	analysisLrrs, MIN_LRR_MEDIAN_ADJUST,
																													MAX_LRR_MEDIAN_ADJUST));
						}
						break;

					default:
						proj.getLog().reportError("Invalid adjustment type " + type);
						analysisBafs = null;
						analysisLrrs = null;
						break;
				}
			}
		}
	}

	private CNVCallResult callCNVS(	int[] chrsToCall, boolean callReverse, int minNumMarkers,
																	double minConf, int numThreads) {
		WorkerHive<CNVCallResult> hive = new WorkerHive<CNVCallResult>(numThreads, 10, proj.getLog());
		boolean[] finalAnalysisSet = Array.booleanArray(markerSet.getMarkerNames().length, false);
		HashMap<String, ArrayList<Integer>> chrIndices = new HashMap<String, ArrayList<Integer>>();
		byte[] chrs = markerSet.getChrs();

		for (int i = 0; i < analysisProjectIndices.length; i++) {
			String chr = Positions.getChromosomeUCSC(chrs[analysisProjectIndices[i]], true);
			if (!chrIndices.containsKey(chr)) {
				chrIndices.put(chr, new ArrayList<Integer>());
			}
			chrIndices.get(chr).add(i);
			finalAnalysisSet[analysisProjectIndices[i]] = true;
		}
		int[][] snpDists = getSNPDist(proj, markerSet, false, finalAnalysisSet);

		if (chrsToCall == null) {
			chrsToCall = Array.arrayOfIndices(snpDists.length);
		}
		for (int curChr : chrsToCall) {
			String chr = Positions.getChromosomeUCSC(curChr, true);
			if (snpDists[curChr].length > MIN_MARKERS_PER_CHROMOSOME && chrIndices.containsKey(chr)) {
				int[] currentChrIndices = Ints.toArray(chrIndices.get(chr));
				int[] currentChrPositions = Array.subArray(analysisPositions, currentChrIndices);
				double[] currentChrLrr = Array.subArray(analysisLrrs, currentChrIndices);
				double[] currentChrBaf = Array.subArray(analysisBafs, currentChrIndices);
				double[] currentChrPfbs = Array.subArray(analysisPfbs, currentChrIndices);
				boolean[] currentChrCnDef = Array.subArray(copyNumberDef, currentChrIndices);
				String[] currentNames = Array.subArray(
																								Array.subArray(	markerSet.getMarkerNames(),
																																analysisProjectIndices),
																								currentChrIndices);
				CNVCallerWorker worker = new CNVCallerWorker(	proj, dna, (byte) curChr, currentChrPositions,
																											currentNames, pennHmm, currentChrLrr,
																											currentChrBaf, currentChrPfbs,
																											snpDists[curChr], currentChrCnDef,
																											callReverse, debugMode);
				hive.addCallable(worker);
			} else {
				if (debugMode) {
					proj.getLog()
							.reportTimeWarning("There were fewer than "	+ MIN_MARKERS_PER_CHROMOSOME
																	+ " analysis markers on chromosome " + chr
																	+ " in the final call set, skipping");
				}
			}
		}
		hive.execute(true);
		ArrayList<CNVCallResult> results = hive.getResults();
		ArrayList<CNVariant> allCNVs = new ArrayList<CNVariant>();
		ArrayList<CNVariant> allReverse = new ArrayList<CNVariant>();
		ArrayList<CNVariant> allReverseConsensus = new ArrayList<CNVariant>();

		for (int i = 0; i < results.size(); i++) {
			for (int j = 0; j < results.get(i).getChrCNVs().getLoci().length; j++) {
				if (results.get(i).getChrCNVs().getLoci()[j].getNumMarkers() >= minNumMarkers
							&& !Double.isNaN(results.get(i).getChrCNVs().getLoci()[j].getScore())
						&& results.get(i).getChrCNVs().getLoci()[j].getScore() > minConf) {
					allCNVs.add(results.get(i).getChrCNVs().getLoci()[j]);
				}
			}
			if (callReverse) {
				throw new IllegalArgumentException("Call reverse is no longer active since it gives identical results. If you really want to use it, go back in git history");
			}
		}

		LocusSet<CNVariant> allLocusSet =
																		new LocusSet<CNVariant>(allCNVs.toArray(new CNVariant[allCNVs.size()]),
																														true, proj.getLog()) {

																			/**
																			 * 
																			 */
																			private static final long serialVersionUID = 1L;

																		};

		LocusSet<CNVariant> allLocusSetReverse = new LocusSet<CNVariant>(	allReverse.toArray(new CNVariant[allReverse.size()]),
																																			true, proj.getLog()) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};
		LocusSet<CNVariant> allLocusSetReverseConsensus =
																										new LocusSet<CNVariant>(allReverseConsensus.toArray(new CNVariant[allReverseConsensus.size()]),
																																						true, proj.getLog()) {

																											/**
																											 * 
																											 */
																											private static final long serialVersionUID = 1L;

																										};
		return new CNVCallResult(allLocusSet, allLocusSetReverse, allLocusSetReverseConsensus);
	}

	/**
	 * The bare bones worker that should operate on any data passed to it
	 *
	 */
	private static class CNVCallerWorker implements Callable<CNVCallResult> {
		private final Project proj;
		private final String dna;
		private final byte currentChr;
		private final int[] positions;
		private final PennHmm pennHmm;
		private final double[] lrrs;
		private final double[] bafs;
		private final double[] pfbs;
		private final int[] snipDists;

		private final boolean[] cnDef;
		private final boolean verbose;
		private final boolean callReverse;
		private final String[] names;

		private CNVCallerWorker(Project proj, String dna, byte currentChr, int[] positions,
														String[] names, PennHmm pennHmm, double[] lrrs, double[] bafs,
														double[] pfbs, int[] snipDists, boolean[] cnDef, boolean callReverse,
														boolean verbose) {
			super();
			this.proj = proj;
			this.dna = dna;
			this.currentChr = currentChr;
			this.positions = positions;
			this.pennHmm = pennHmm;
			this.lrrs = lrrs;
			this.bafs = bafs;
			this.pfbs = pfbs;
			this.snipDists = snipDists;
			this.cnDef = cnDef;
			this.names = names;
			this.callReverse = callReverse;
			this.verbose = verbose;
		}

		@Override
		public CNVCallResult call() throws Exception {

			try {
				ViterbiResult viterbiResult = PennHmm.ViterbiLogNP_CHMM(pennHmm, lrrs, bafs, pfbs,
																																snipDists, cnDef);
				LocusSet<CNVariant> chrCnvs = viterbiResult.analyzeStateSequence(	proj, dna, dna, currentChr,
																																					positions, names, 2,
																																					false, verbose);
				LocusSet<CNVariant> chrCnvsReverse = null;
				LocusSet<CNVariant> chrCNVsReverseConsensus = null;

				chrCnvs = PennHmm.scoreCNVsSameChr(	pennHmm, chrCnvs, positions, lrrs, bafs, pfbs, cnDef,
																						viterbiResult.getQ(), 2, proj.getLog());
				if (callReverse) {
					throw new IllegalArgumentException("Call reverse is no longer active since it gives identical results.");

				}
				CNVCallResult callResult = new CNVCallResult(	chrCnvs, chrCnvsReverse,
																											chrCNVsReverseConsensus);
				return callResult;
			} catch (Exception e) {
				proj.getLog().reportError("Could not call cnvs for sample "	+ dna + " on chromosome "
																			+ currentChr);
				proj.getLog().reportException(e);
				throw new IllegalStateException("Could not call cnvs for sample "	+ dna + " on chromosome "
																				+ currentChr);
			}

		}
	}

	public static class CNVCallResult {
		private final LocusSet<CNVariant> chrCNVs;
		private final LocusSet<CNVariant> chrCNVsReverseConsensus;
		private final LocusSet<CNVariant> chrCNVsReverse;

		public CNVCallResult(	LocusSet<CNVariant> chrCNVs, LocusSet<CNVariant> chrCNVsReverse,
													LocusSet<CNVariant> chrCNVsReverseConsensus) {
			super();
			this.chrCNVs = chrCNVs;
			this.chrCNVsReverse = chrCNVsReverse;
			this.chrCNVsReverseConsensus = chrCNVsReverseConsensus;

		}

		public LocusSet<CNVariant> getChrCNVs() {
			return chrCNVs;
		}

		public LocusSet<CNVariant> getChrCNVsReverseConsensus() {
			return chrCNVsReverseConsensus;
		}

		public LocusSet<CNVariant> getChrCNVsReverse() {
			return chrCNVsReverse;
		}

	}

	/**
	 * Flags for the PennCNV data adjustments, can be used to assign order of operations as well;
	 *
	 */
	private enum DATA_ADJUSTMENTS {
																	/**
																	 * PennCNV uses 4 sigFigs
																	 */
																	ROUND_TO_PENNCNV_SIG_FIGS,

																	/**
																	 * Subset to the analysis set, such as autosomal only, etc
																	 */
																	SUBSET_TO_ANALYSIS_MARKERS,
																	/**
																	 * Replaces NaN entries in either baf or lrr with 0 for both baf
																	 * and lrr
																	 */
																	HANDLE_NAN,
																	/**
																	 * Does the PennCNV gc adjustment on the current data
																	 */
																	GC_ADJUST,
																	/**
																	 * Does the PennCNV median lrr and baf adjustment on the current
																	 * data
																	 */
																	MEDIAN_ADJUST,

																	/**
																	 * Adjust the hmm model by the standard deviation of the current
																	 * data
																	 */
																	ADJUST_HMM_SD;
	}

	/**
	 * @param proj
	 * @return the plus one index physical distance of each marker (of the analysis set) on each
	 *         chromosome
	 */
	private static int[][] getSNPDist(Project proj, PreparedMarkerSet markerSet, boolean reverse,
																		boolean[] projectIndicesToUse) {
		int[][] chrPos = markerSet.getPositionsByChr();
		int[][] tmp = new int[chrPos.length][];
		int projectIndex = 0;
		for (int i = 0; i < chrPos.length; i++) {
			ArrayList<Integer> updated = new ArrayList<Integer>();
			for (int j = 0; j < chrPos[i].length; j++) {
				if (projectIndicesToUse[projectIndex]) {
					updated.add(chrPos[i][j]);
				}
				projectIndex++;
			}
			tmp[i] = Ints.toArray(updated);
		}
		int[][] snpDists = new int[tmp.length][];
		for (int i = 0; i < tmp.length; i++) {
			int[] distsTmp = new int[tmp[i].length];
			int[] posTmp = reverse ? Array.reverse(tmp[i]) : tmp[i];
			if (distsTmp.length > 0) {
				distsTmp[posTmp.length - 1] = 0;
				for (int j = 0; j < distsTmp.length - 1; j++) {
					int dist = Math.abs(posTmp[j + 1] - posTmp[j]);// abs in case of reverese
					distsTmp[j] = dist > 0 ? dist : 1;
				}
			}
			snpDists[i] = distsTmp;

		}
		return snpDists;

	}

	/**
	 * Median adjust these lrr values like PennCNV
	 */
	public static double[] adjustLrr(	double[] lrrs, double minLrr, double maxLrr, boolean verbose,
																		Logger log) {
		double[] adjusted = new double[lrrs.length];
		double median = Array.median(Array.getValuesBetween(lrrs, minLrr, maxLrr));
		if (verbose) {
			log.reportTimeInfo("Median adjusting lrr values by " + median + "\t" + Array.median(lrrs));
		}
		for (int i = 0; i < adjusted.length; i++) {
			adjusted[i] = lrrs[i] - median;
		}
		return adjusted;
	}

	/**
	 * Median adjust these baf values like PennCNV
	 */
	public static double[] adjustBaf(	double[] bafs, double minBaf, double maxBaf, boolean debugMode,
																		Logger log) {
		double[] adjusted = new double[bafs.length];
		ArrayList<Double> bafsToMedian = new ArrayList<Double>();
		for (int i = 0; i < bafs.length; i++) {
			if (!Double.isNaN(bafs[i]) && bafs[i] > minBaf && bafs[i] < maxBaf) {
				bafsToMedian.add(bafs[i]);
			}
		}
		double median = Array.median(Doubles.toArray(bafsToMedian));
		double factor = median - 0.5;
		if (debugMode) {
			log.reportTimeInfo("Median adjusting baf measures by " + factor);
		}
		for (int i = 0; i < adjusted.length; i++) {
			if (!Double.isNaN(bafs[i]) && bafs[i] > minBaf && bafs[i] < maxBaf) {
				adjusted[i] = bafs[i] - factor;
			} else {
				adjusted[i] = bafs[i];
			}

		}
		return adjusted;
	}

	/**
	 * @return the {@link DATA_ADJUSTMENTS} enums in the order that should replicate PennCNV with
	 *         gcCorrection
	 */
	private static DATA_ADJUSTMENTS[] getPennCNVGCProcessingOrder() {
		ArrayList<DATA_ADJUSTMENTS> da = new ArrayList<CNVCaller.DATA_ADJUSTMENTS>();
		da.add(DATA_ADJUSTMENTS.HANDLE_NAN);
		da.add(DATA_ADJUSTMENTS.GC_ADJUST);
		da.add(DATA_ADJUSTMENTS.SUBSET_TO_ANALYSIS_MARKERS);
		da.add(DATA_ADJUSTMENTS.MEDIAN_ADJUST);
		da.add(DATA_ADJUSTMENTS.ADJUST_HMM_SD);

		return da.toArray(new DATA_ADJUSTMENTS[da.size()]);
	}

	/**
	 * @return the {@link DATA_ADJUSTMENTS} enums in the order that should replicate PennCNV without
	 *         gcCorrection
	 */
	private static DATA_ADJUSTMENTS[] getPennCNVProcessingOrder() {
		ArrayList<DATA_ADJUSTMENTS> da = new ArrayList<CNVCaller.DATA_ADJUSTMENTS>();
		da.add(DATA_ADJUSTMENTS.HANDLE_NAN);
		da.add(DATA_ADJUSTMENTS.SUBSET_TO_ANALYSIS_MARKERS);
		da.add(DATA_ADJUSTMENTS.MEDIAN_ADJUST);
		da.add(DATA_ADJUSTMENTS.ADJUST_HMM_SD);

		return da.toArray(new DATA_ADJUSTMENTS[da.size()]);
	}

	private static CNVCallResult callCNVsFor(	Project proj, PennHmm pennHmm, String sampleName,
																						double[] sampLrrs, double[] sampBafs, GcModel gcModel,
																						PFB pfb, PreparedMarkerSet markerSet,
																						boolean[] markersToUse, boolean[] copyNumberDef,
																						int[] chrsToCall, boolean callReverse,
																						int minNumMarkers, double minConf, int numThreads,
																						PFB_MANAGEMENT_TYPE pManagementType,
																						boolean debugMode) {
		DATA_ADJUSTMENTS[] dAdjustments;
		if (gcModel == null) {
			dAdjustments = getPennCNVProcessingOrder();
			proj.getLog().reportTimeWarning("No gc model was provided, calling on cnvs on raw data");

		} else {
			dAdjustments = getPennCNVGCProcessingOrder();
		}

		CNVCaller caller = new CNVCaller(	proj, sampleName, pennHmm, gcModel, pfb, dAdjustments,
																			markerSet, markersToUse, copyNumberDef, sampLrrs, sampBafs,
																			pManagementType, debugMode);
		caller.adjustData();

		return caller.callCNVS(chrsToCall, callReverse, minNumMarkers, minConf, numThreads);

	}

	/**
	 * @param proj
	 * @param pennHmm required {@link PennHmm}
	 * @param sample sample to call cnvs on
	 * @param gcModel optional {@link GcModel}
	 * @param pfb required {@link PFB}
	 * @param markerSet
	 * @param numThreads analysis is split up by chromosome
	 * @param debugMode report more values and warnings to compare with regular penncnv calls
	 * @return
	 */
	public static CNVCallResult callCNVsFor(Project proj, final PennHmm pennHmm, String sample,
																					double[] lrrs, double[] bafs, final GcModel gcModel,
																					final PFB pfb, final PreparedMarkerSet markerSet,
																					int[] chrsToCall, boolean[] markersToUse,
																					boolean callReverse, int minNumMarkers, double minConf,
																					PFB_MANAGEMENT_TYPE pManagementType, int numThreads,
																					boolean debugMode) {
		String[] markerNames = markerSet.getMarkerNames();
		boolean[] copyNumberDef = Array.booleanArray(markerNames.length, false);
		ARRAY array = proj.getArrayType();
		if (debugMode) {
			proj.getLog()
					.reportTimeInfo("Assigning copy number probes according to "	+ array.toString()
													+ " using the following " + Array.toStr(array.getCnFlags(), ","));
			proj.getLog()
					.reportTimeInfo("BAF values greater than 1 will also be set to copy number only");

		}
		for (int i = 0; i < copyNumberDef.length; i++) {
			copyNumberDef[i] = array.isCNOnly(markerNames[i]) || bafs[i] > 1;
		}
		if (debugMode) {
			proj.getLog().reportTimeInfo("Found "	+ Array.booleanArraySum(copyNumberDef)
																		+ " copy number only markers");
		}
		if (markersToUse == null) {
			int[] autosomalMarkers = proj.getAutosomalMarkerIndices();
			markersToUse = Array.booleanArray(markerSet.getMarkerNames().length, false);
			for (int i = 0; i < autosomalMarkers.length; i++) {
				markersToUse[autosomalMarkers[i]] = true;
				if (array == ARRAY.NGS) {
					String name = markerNames[autosomalMarkers[i]];
					markersToUse[autosomalMarkers[i]] =
																						NGS_MARKER_TYPE.getType(name) != NGS_MARKER_TYPE.VARIANT_SITE;
				}
			}
		}
		if (debugMode) {
			proj.getLog()
					.reportTimeInfo("Found " + Array.booleanArraySum(markersToUse) + " total markers to use");

		}
		CNVCallResult cnvs = callCNVsFor(	proj, pennHmm, sample, lrrs, bafs, gcModel, pfb, markerSet,
																			markersToUse, copyNumberDef, chrsToCall, callReverse,
																			minNumMarkers, minConf, numThreads, pManagementType,
																			debugMode);
		System.gc();
		return cnvs;
	}

	/**
	 * @author lane0212 Useful when calling cnvs across many samples
	 */
	private static class CNVProducer extends AbstractProducer<CNVCallResult> {
		private final Project proj;
		private final PennHmm pennHmm;
		private final GcModel gcModel;
		private final PFB pfb;
		private final String[] samples;
		private final int[] chrsToCall;
		private final boolean[] markersToUse;
		private final Centroids centroids;
		private final int numSampleThreads;
		private int index;
		private final boolean debugMode, callReverse;
		private final PreparedMarkerSet markerSet;
		private final int minNumMarkers;
		private final double minConf;
		private final PFB_MANAGEMENT_TYPE pManagementType;

		/**
		 * @param proj
		 * @param pennHmm
		 * @param gcModel
		 * @param pfb
		 * @param samples call cnvs on these samples
		 * @param chrsToCall call cnvs on these chromosomes only
		 * @param centroids use these centroids for deriving lrr/baf
		 * @param numSampleThreads
		 * @param debugMode
		 */
		private CNVProducer(Project proj, PreparedMarkerSet markerSet, PennHmm pennHmm, GcModel gcModel,
												PFB pfb, String[] samples, int[] chrsToCall, boolean[] markersToUse,
												Centroids centroids, int minNumMarkers, double minConf,
												int numSampleThreads, boolean callReverse,
												PFB_MANAGEMENT_TYPE pManagementType, boolean debugMode) {
			super();
			this.proj = proj;
			this.markerSet = markerSet;
			this.pennHmm = pennHmm;
			this.gcModel = gcModel;
			this.pfb = pfb;
			this.samples = samples;
			this.chrsToCall = chrsToCall;
			this.markersToUse = markersToUse;
			this.centroids = centroids;
			this.numSampleThreads = numSampleThreads;
			index = 0;
			this.callReverse = callReverse;
			this.debugMode = debugMode;
			this.minNumMarkers = minNumMarkers;
			this.minConf = minConf;
			this.pManagementType = pManagementType;
		}

		@Override
		public boolean hasNext() {
			return index < samples.length;
		}

		@Override
		public Callable<CNVCallResult> next() {
			final String sample = samples[index];
			final PennHmm pennHmmTmp = new PennHmm(pennHmm);
			final GcModel gcModelTmp = gcModel == null ? null : new GcModel(gcModel);
			final PFB pfbTmp = new PFB(pfb);
			Callable<CNVCallResult> callable = new Callable<CNVCallResult>() {

				@Override
				public CNVCallResult call() throws Exception {
					Sample curSample = proj.getFullSampleFromRandomAccessFile(sample);
					if (curSample.getFingerprint() != markerSet.getFingerprint()) {
						throw new IllegalStateException("Mismatched fingerprint for sample and marker set "
																							+ curSample.getLRRs().length + "\t"
																						+ markerSet.getMarkerNames().length);
					}
					float[] lrrs = centroids == null	? curSample.getLRRs()
																						: curSample.getLRRs(centroids.getCentroids());
					float[] bafs = centroids == null	? curSample.getBAFs()
																						: curSample.getBAFs(centroids.getCentroids());
					CNVCallResult cnvs = callCNVsFor(	proj, pennHmmTmp, curSample.getSampleName(),
																						Array.toDoubleArray(lrrs), Array.toDoubleArray(bafs),
																						gcModelTmp, pfbTmp, markerSet, chrsToCall, markersToUse,
																						callReverse, minNumMarkers, minConf, pManagementType,
																						numSampleThreads, debugMode);
					System.gc();
					return cnvs;
				}

			};
			index++;
			return callable;
		}
	}

	public static class CNVCallerIterator implements Iterator<CNVCallResult> {
		private final WorkerTrain<CNVCallResult> train;

		private CNVCallerIterator(WorkerTrain<CNVCallResult> train) {
			super();
			this.train = train;
		}

		@Override
		public boolean hasNext() {
			return train.hasNext();
		}

		@Override
		public CNVCallResult next() {
			return train.next();
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException();
		}
	}

	/**
	 * Possible autosomal, and sex chrs cnv calling method
	 */
	public static void callGenomeCnvs(Project proj, String outputFile, String[] maleSamples,
																		String[] femaleSamples, Centroids[] centroids,
																		int minNumMarkers, double minConf,
																		PFB_MANAGEMENT_TYPE pManagementType, int numSampleThreads,
																		int numChrThreads) {
		PreparedMarkerSet markerSet = PreparedMarkerSet.getPreparedMarkerSet(proj.getMarkerSet());
		String output = proj.PROJECT_DIRECTORY.getValue() + outputFile;
		proj.getLog().reportTimeInfo("CNVS will be reported to " + output);
		new File(ext.parseDirectoryOfFile(output)).mkdirs();
		CNVCallerIterator callerIterator = getCallerIterator(proj, markerSet, Array.concatAll(maleSamples, femaleSamples), null, null,
											centroids[0], minNumMarkers, minConf, pManagementType, numSampleThreads,
											numChrThreads);
		writeOutput(callerIterator, output, proj.getLog());
		// will passing null to chrsToCall result in calling on 23/24 also?
		boolean[] chr23 = Array.booleanArray(markerSet.getMarkerNames().length, false);
		int[][] indicesByChr = markerSet.getIndicesByChr();
		for (int i = 0; i < indicesByChr[23].length; i++) {
			chr23[indicesByChr[23][i]] = true;
		}

		boolean[] chr24 = Array.booleanArray(markerSet.getMarkerNames().length, false);
		for (int i = 0; i < indicesByChr[24].length; i++) {
			chr24[indicesByChr[24][i]] = true;
		}
		callerIterator = getCallerIterator(proj, markerSet, maleSamples, new int[] {23}, chr23, centroids[1],
											minNumMarkers, minConf, pManagementType, numSampleThreads, numChrThreads);
		output = proj.PROJECT_DIRECTORY.getValue() + "23_M_" + outputFile;
		writeOutput(callerIterator, output, proj.getLog());
		callerIterator = getCallerIterator(proj, markerSet, femaleSamples, new int[] {23}, chr23, centroids[2],
											minNumMarkers, minConf, pManagementType, numSampleThreads, numChrThreads);
		output = proj.PROJECT_DIRECTORY.getValue() + "23_F_" + outputFile;
		writeOutput(callerIterator, output, proj.getLog());
		callerIterator = getCallerIterator(proj, markerSet, maleSamples, new int[] {24}, chr24, centroids[1],
											minNumMarkers, minConf, pManagementType, numSampleThreads, numChrThreads);
		output = proj.PROJECT_DIRECTORY.getValue() + "24_M_" + outputFile;
		writeOutput(callerIterator, output, proj.getLog());

	}

	/**
	 * @param proj
	 * @param call on these samples only
	 * @param chrsToCall call cnvs on these chromosomes only
	 * @param markersToUse use only these markers when calling cnvs...
	 * @param centroids use these centroids to call cnvs
	 * @param numSampleThreads number of samples analyzed at once.
	 * @param numChrThreads number of chromosomes in each sample analyzed at once
	 *
	 *        NOTE: total thread usage is numSampleThreads*numChrThreads
	 */
	public static CNVCallerIterator getCallerIterator(Project proj, PreparedMarkerSet markerSet,
																										String[] samples, int[] chrsToCall,
																										boolean[] markersToUse, Centroids centroids,
																										int minNumMarkers, double minConf,
																										PFB_MANAGEMENT_TYPE pManagementType,
																										int numSampleThreads, int numChrThreads) {

		PennHmm pennHmmOriginal = PennHmm.loadPennHmm(proj.HMM_FILENAME.getValue(), new Logger());
		if (!Files.exists(proj.CUSTOM_PFB_FILENAME.getValue())) {
			proj.getLog().reportTimeInfo("Did not find "	+ proj.CUSTOM_PFB_FILENAME.getValue()
																		+ ", attempting to generate it now");
			PennCNV.populationBAF(proj);
		}
		PFB pfb = PFB.loadPFB(proj, proj.CUSTOM_PFB_FILENAME.getValue());
		if (!Files.exists(proj.GC_MODEL_FILENAME.getValue(false, false))) {
			Resource gmodelBase = Resources	.genome(proj.GENOME_BUILD_VERSION.getValue(), proj.getLog())
																			.getModelBase();
			if (!Files.exists(proj.GC_MODEL_FILENAME.getValue()) && gmodelBase.isAvailable()) {
				proj.getLog()
						.reportTimeWarning("Generating gcModel for "	+ proj.GENOME_BUILD_VERSION.getValue()
																+ " at " + proj.GC_MODEL_FILENAME.getValue() + " from "
																+ gmodelBase.get());
				proj.getLog().setLevel(3);
				PennCNV.gcModel(proj, gmodelBase.get(), proj.GC_MODEL_FILENAME.getValue(), 100);
			}
		}
		GcModel gcModel = GcAdjustor.GcModel.populateFromFile(
																													proj.GC_MODEL_FILENAME.getValue(false,
																																													false),
																													false, proj.getLog());
		if (gcModel == null) {
			proj.getLog().reportTimeWarning("Calling cnvs without gc correction");
		}
		CNVProducer producer = new CNVProducer(	proj, markerSet, pennHmmOriginal, gcModel, pfb, samples,
																						chrsToCall, markersToUse, centroids, minNumMarkers,
																						minConf, numChrThreads, false, pManagementType, true);
		WorkerTrain<CNVCallResult> train = new WorkerTrain<CNVCallResult>(producer, numSampleThreads, 2,
																																			proj.getLog());
		return new CNVCallerIterator(train);
	}

	/**
	 * This method calls autosomal cnvs only
	 *
	 * @param output relative to the project directory
	 */
	public static void callAutosomalCNVs(	Project proj, String output, String[] samples,
																				int[] chrsToCall, Centroids centroids, int minNumMarkers,
																				double minConf, PFB_MANAGEMENT_TYPE pManagementType,
																				int numSampleThreads, int numChrThreads) {
		PreparedMarkerSet markerSet = PreparedMarkerSet.getPreparedMarkerSet(proj.getMarkerSet());
		output = proj.PROJECT_DIRECTORY.getValue() + output;
		proj.getLog().reportTimeInfo("CNVS will be reported to " + output);
		new File(ext.parseDirectoryOfFile(output)).mkdirs();
		CNVCallerIterator callerIterator = getCallerIterator(	proj, markerSet, samples, chrsToCall, null,
																													centroids, minNumMarkers, minConf,
																													pManagementType, numSampleThreads,
																													numChrThreads);
		writeOutput(callerIterator, output, proj.getLog());
	}
	
	private static void writeOutput(CNVCallerIterator callerIterator, String output, Logger log) {
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
			int sum = 0;
			int cnt = 0;
			while (callerIterator.hasNext()) {
				cnt++;
				LocusSet<CNVariant> cnvs = callerIterator.next().getChrCNVs();
				sum += cnvs.getLoci().length;
				for (int i = 0; i < cnvs.getLoci().length; i++) {
					writer.println(cnvs.getLoci()[i].toPlinkFormat());
				}
			}
			log.reportTimeInfo("Wrote " + sum + " CNVs for " + cnt + " samples to " + output);
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + output);
			log.reportException(e);
		}
	}


	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String output = "genvisis.cnvs";
		int numThreads = 24;
		int minNumMarkers = DEFUALT_MIN_SITES;
		double minConf = DEFUALT_MIN_CONF;

		String usage = "\n" + "cnv.hmm.CNVCaller requires 0-1 arguments\n";
		usage += "   (1) proj (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) output file (relative to project directory) (i.e. out="	+ filename
							+ " (default))\n" + "";
		usage += "   (3) minimum number of markers to report a cnv  (i.e. minMarkers="	+ minNumMarkers
							+ " (default))\n" + "";
		usage +=
					"   (4) minimum confidence report a cnv  (i.e. minConf=" + minConf + " (default))\n" + "";

		usage += PSF.Ext.getNumThreadsCommand(4, numThreads);

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("proj=")) {
				filename = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("out=")) {
				output = ext.parseStringArg(arg, "");

				numArgs--;
			} else if (arg.startsWith("minMarkers=")) {
				minNumMarkers = ext.parseIntArg(arg);

				numArgs--;
			} else if (arg.startsWith("minConf=")) {
				minConf = ext.parseDoubleArg(arg);

				numArgs--;
			} else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(arg);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			Project proj = new Project(filename, false);
			callAutosomalCNVs(proj, output, proj.getSamples(), null, null, minNumMarkers, minConf,
												PFB_MANAGEMENT_TYPE.PENNCNV_DEFAULT, numThreads, 1);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
