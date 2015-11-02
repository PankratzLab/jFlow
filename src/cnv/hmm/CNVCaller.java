package cnv.hmm;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.concurrent.Callable;

import common.Array;
import common.Logger;
import common.PSF;
import common.Positions;
import common.WorkerHive;
import common.WorkerTrain;
import common.ext;
import common.WorkerTrain.Producer;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Project.ARRAY;
import cnv.filesys.Sample;
import cnv.hmm.PennHmm.ViterbiResult;
import cnv.qc.GcAdjustor;
import cnv.qc.GcAdjustor.GcModel;
import cnv.var.CNVariant;
import cnv.var.CNVariant.CNVBuilder;
import cnv.var.LocusSet;
import cnv.var.LocusSet.TO_STRING_TYPE;
import filesys.Segment;

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

	private static final int MIN_MARKERS_PER_CHROMOSOME = 10;

	private Project proj;
	private String dna;
	private PennHmm pennHmm;
	private GcModel gcModel;
	private MarkerSet markerSet;
	private boolean[] markersToUse, copyNumberDef;
	private double[] analysisLrrs, analysisBafs, analysisPfbs;
	private int[] analysisProjectIndices, analysisPositions;// linker of what is being analyzed to the order in the project
	private DATA_ADJUSTMENTS[] dataAdjustments;
	private boolean debugMode;

	public CNVCaller(Project proj, String dna, PennHmm pennHmm, GcModel gcModel, PFB pfb, DATA_ADJUSTMENTS[] dataAdjustments, MarkerSet markerSet, boolean[] markersToUse, boolean[] copyNumberDef, double[] lrrs, double[] bafs, boolean debugMode) {
		super();
		this.proj = proj;
		this.pennHmm = pennHmm;
		this.dna = dna;
		this.gcModel = gcModel;
		this.markerSet = markerSet;
		this.markersToUse = markersToUse;
		this.analysisLrrs = lrrs;
		this.analysisBafs = bafs;
		this.analysisPfbs = pfb.getPfbs();
		this.dataAdjustments = dataAdjustments;
		this.copyNumberDef = copyNumberDef;
		this.analysisProjectIndices = Array.intArray(markerSet.getMarkerNames().length);
		this.analysisPositions = markerSet.getPositions();
		this.debugMode = debugMode;
		if (lrrs.length != bafs.length || markerSet.getMarkerNames().length != lrrs.length) {
			String error = "BUG: must supply entire lrr and baf data set for the project, consider using the markersToUse array to subset the analysis";
			proj.getLog().reportTimeError(error);
			throw new IllegalArgumentException(error);
		}
		if (debugMode) {
			proj.getLog().reportTimeInfo("Sample: " + dna);
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
			double lrrSd = Array.stdev(Array.getValuesBetween(analysisLrrs, MIN_LRR_MEDIAN_ADJUST, MAX_LRR_MEDIAN_ADJUST));
			for (int i = 0; i < dataAdjustments.length; i++) {
				DATA_ADJUSTMENTS type = dataAdjustments[i];
				switch (type) {
				case ROUND_TO_PENNCNV_SIG_FIGS:// May not be necessary
					analysisLrrs = roundToPennCNVSigFigs(analysisLrrs);
					analysisBafs = roundToPennCNVSigFigs(analysisBafs);
					analysisPfbs = roundToPennCNVSigFigs(analysisPfbs);
					// lrrSd = Array.stdev(getValuesBetween(analysisLrrs, MIN_LRR_MEDIAN_ADJUST, MAX_LRR_MEDIAN_ADJUST));
					break;
				case ADJUST_HMM_SD:
					pennHmm = PennHmm.adjustBSD(pennHmm, lrrSd, proj.getLog());
					break;
				case GC_ADJUST:
					if (gcModel == null) {
						String error = "gc model cannot be null if adjustment " + type + " is flagged";
						proj.getLog().reportTimeError(error);
						throw new IllegalArgumentException();
					} else {
						double[] dataToCorrect = analysisLrrs;
						if (analysisProjectIndices.length != markerSet.getMarkerNames().length) {// only current indicies will be used
							dataToCorrect = new double[markerSet.getMarkerNames().length];
							for (int j = 0; j < analysisProjectIndices.length; j++) {
								dataToCorrect[analysisProjectIndices[j]] = analysisLrrs[j];
							}
						}
						GcAdjustor gcAdjustor = new GcAdjustor(proj, gcModel, dataToCorrect, null, true, debugMode);
						gcAdjustor.correctIntensities();
						gcAdjustor.computeQCMetrics(true, true);
						analysisLrrs = Array.subArray(gcAdjustor.getCorrectedIntensities(), analysisProjectIndices);
						proj.getLog().reportTimeInfo(gcAdjustor.getAnnotatedQCString());
						lrrSd = Array.stdev(Array.getValuesBetween(analysisLrrs, MIN_LRR_MEDIAN_ADJUST, MAX_LRR_MEDIAN_ADJUST));

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
						proj.getLog().reportTimeInfo("Set " + nanCount + " observations with NaN data to 0 for lrr and baf");
					}
					lrrSd = Array.stdev(Array.getValuesBetween(analysisLrrs, MIN_LRR_MEDIAN_ADJUST, MAX_LRR_MEDIAN_ADJUST));

					break;
				case MEDIAN_ADJUST:
					analysisLrrs = adjustLrr(analysisLrrs, MIN_LRR_MEDIAN_ADJUST, MAX_LRR_MEDIAN_ADJUST, debugMode, proj.getLog());
					analysisBafs = adjustBaf(analysisBafs, MIN_BAF_MEDIAN_ADJUST, MAX_BAF_MEDIAN_ADJUST,debugMode, proj.getLog());
					// TODO, update lrrSd later?
					// lrrSd = Array.stdev(getValuesBetween(analysisLrrs, MIN_LRR_MEDIAN_ADJUST, MAX_LRR_MEDIAN_ADJUST)); PennCNV does not update lrr sd here so we wont either
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
						lrrSd = Array.stdev(Array.getValuesBetween(analysisLrrs, MIN_LRR_MEDIAN_ADJUST, MAX_LRR_MEDIAN_ADJUST));
					}
					break;

				default:
					proj.getLog().reportTimeError("Invalid adjustment type " + type);
					analysisBafs = null;
					analysisLrrs = null;
					break;
				}
			}
		}
	}

	private CNVCallResult callCNVS(int[] chrsToCall, boolean callReverse, int numThreads) {
		WorkerHive<CNVCallResult> hive = new WorkerHive<CNVCallResult>(numThreads, 10, proj.getLog());
		boolean[] finalAnalysisSet = Array.booleanArray(markerSet.getMarkerNames().length, false);
		Hashtable<String, ArrayList<Integer>> chrIndices = new Hashtable<String, ArrayList<Integer>>();
		byte[] chrs = markerSet.getChrs();

		for (int i = 0; i < analysisProjectIndices.length; i++) {
			String chr = Positions.getChromosomeUCSC(chrs[analysisProjectIndices[i]], true);
			if (!chrIndices.containsKey(chr)) {
				chrIndices.put(chr, new ArrayList<Integer>());
			}
			chrIndices.get(chr).add(i);
			finalAnalysisSet[analysisProjectIndices[i]] = true;
		}
		int[][] snpDists = getSNPDist(proj, false, finalAnalysisSet);
		int[][] snpDistsReverse = getSNPDist(proj, true, finalAnalysisSet);

		if (chrsToCall == null) {
			chrsToCall = Array.intArray(snpDists.length);
		}
		for (int i = 0; i < chrsToCall.length; i++) {
			int curChr = chrsToCall[i];

			String chr = Positions.getChromosomeUCSC(curChr, true);
			// if (i == 16) {
			if (snpDists[curChr].length > MIN_MARKERS_PER_CHROMOSOME && chrIndices.containsKey(chr)) {
				int[] currentChrIndices = Array.toIntArray(chrIndices.get(chr));
				int[] currentChrPositions = Array.subArray(analysisPositions, currentChrIndices);
				double[] currentChrLrr = Array.subArray(analysisLrrs, currentChrIndices);
				double[] currentChrBaf = Array.subArray(analysisBafs, currentChrIndices);
				double[] currentChrPfbs = Array.subArray(analysisPfbs, currentChrIndices);
				boolean[] currentChrCnDef = Array.subArray(copyNumberDef, currentChrIndices);
				String[] currentNames = Array.subArray(Array.subArray(markerSet.getMarkerNames(), analysisProjectIndices), currentChrIndices);
				CNVCallerWorker worker = new CNVCallerWorker(proj, dna, (byte) curChr, currentChrPositions, currentNames, pennHmm, currentChrLrr, currentChrBaf, currentChrPfbs, snpDists[curChr], snpDistsReverse[curChr], currentChrCnDef, callReverse, debugMode);
				hive.addCallable(worker);
			} else {
				if (debugMode) {
					proj.getLog().reportTimeWarning("There were fewer than " + MIN_MARKERS_PER_CHROMOSOME + " analysis markers on chromosome " + chr + " in the final call set, skipping");
				}
			}
			// } else {
			// proj.getLog().reportTimeError("REMVEMEBDE CHR!^");
			// }
		}
		hive.execute(true);
		ArrayList<CNVCallResult> results = hive.getResults();
		ArrayList<CNVariant> allCNVs = new ArrayList<CNVariant>();
		ArrayList<CNVariant> allReverse = new ArrayList<CNVariant>();
		ArrayList<CNVariant> allReverseConsensus = new ArrayList<CNVariant>();

		for (int i = 0; i < results.size(); i++) {
			for (int j = 0; j < results.get(i).getChrCNVs().getLoci().length; j++) {
				allCNVs.add(results.get(i).getChrCNVs().getLoci()[j]);
			}
			for (int j = 0; j < results.get(i).getChrCNVsReverse().getLoci().length; j++) {
				allReverse.add(results.get(i).getChrCNVsReverse().getLoci()[j]);
			}
			for (int j = 0; j < results.get(i).getChrCNVsReverseConsensus().getLoci().length; j++) {
				allReverseConsensus.add(results.get(i).getChrCNVsReverseConsensus().getLoci()[j]);
			}
		}

		LocusSet<CNVariant> allLocusSet = new LocusSet<CNVariant>(allCNVs.toArray(new CNVariant[allCNVs.size()]), true, proj.getLog()) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};

		LocusSet<CNVariant> allLocusSetReverse = new LocusSet<CNVariant>(allReverse.toArray(new CNVariant[allReverse.size()]), true, proj.getLog()) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};
		LocusSet<CNVariant> allLocusSetReverseConsensus = new LocusSet<CNVariant>(allReverseConsensus.toArray(new CNVariant[allReverseConsensus.size()]), true, proj.getLog()) {

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
		private Project proj;
		private String dna;
		private byte currentChr;
		private int[] positions;
		private PennHmm pennHmm;
		private double[] lrrs;
		private double[] bafs;
		private double[] pfbs;
		private int[] snipDists;
		private int[] snipDistsReverse;

		private boolean[] cnDef;
		private boolean verbose, callReverse;
		private String[] names;

		private CNVCallerWorker(Project proj, String dna, byte currentChr, int[] positions, String[] names, PennHmm pennHmm, double[] lrrs, double[] bafs, double[] pfbs, int[] snipDists, int[] snipDistsReverse, boolean[] cnDef, boolean callReverse, boolean verbose) {
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
			this.snipDistsReverse = snipDistsReverse;
			this.cnDef = cnDef;
			this.names = names;
			this.callReverse = callReverse;
			this.verbose = verbose;
		}

		@Override
		public CNVCallResult call() throws Exception {
			ViterbiResult viterbiResult = PennHmm.ViterbiLogNP_CHMM(pennHmm, lrrs, bafs, pfbs, snipDists, cnDef);
			LocusSet<CNVariant> chrCnvs = viterbiResult.analyzeStateSequence(proj, dna, dna, currentChr, positions, names, 2, false, verbose);
			LocusSet<CNVariant> chrCnvsReverse = null;
			LocusSet<CNVariant> chrCNVsReverseConsensus = null;

			chrCnvs = PennHmm.scoreCNVsSameChr(pennHmm, chrCnvs, positions, lrrs, bafs, pfbs, cnDef, proj.getLog());
			if (callReverse) {
				ViterbiResult viterbiResultReverse = PennHmm.ViterbiLogNP_CHMM(pennHmm, Array.reverse(lrrs), Array.reverse(bafs), Array.reverse(pfbs), snipDistsReverse, Array.reverse(cnDef));
				chrCnvsReverse = viterbiResultReverse.analyzeStateSequence(proj, dna, dna, currentChr, positions, names, 2, true, verbose);
				chrCnvsReverse = PennHmm.scoreCNVsSameChr(pennHmm, chrCnvs, positions, lrrs, bafs, pfbs, cnDef, proj.getLog());
				chrCNVsReverseConsensus = developConsensus(chrCnvs, chrCnvsReverse, positions, proj.getLog());
				chrCNVsReverseConsensus = PennHmm.scoreCNVsSameChr(pennHmm, chrCNVsReverseConsensus, positions, lrrs, bafs, pfbs, cnDef, proj.getLog());
			}
			CNVCallResult callResult = new CNVCallResult(chrCnvs, chrCnvsReverse, chrCNVsReverseConsensus);
			return callResult;
		}
	}

	private static LocusSet<CNVariant> developConsensus(LocusSet<CNVariant> chrCnvs, LocusSet<CNVariant> chrCnvsReverse, int[] positions, Logger log) {
		ArrayList<CNVariant> consensus = new ArrayList<CNVariant>();
		if (chrCnvsReverse != null) {
			for (int i = 0; i < chrCnvs.getLoci().length; i++) {
				CNVariant tmp = chrCnvs.getLoci()[i];

				CNVariant[] revereseOverlaps = chrCnvsReverse.getOverLappingLoci(tmp);
				if (revereseOverlaps != null && revereseOverlaps.length > 0) {
					for (int j = 0; j < revereseOverlaps.length; j++) {
						if (tmp.getCN() == revereseOverlaps[j].getCN()) {
							Segment intersect = tmp.getIntersection(revereseOverlaps[j], log);
							if (intersect.getSize() != tmp.getSize()) {
								System.err.println("CONSENSUS FLAG");
							}
							CNVBuilder builder = new CNVBuilder(tmp);
							builder.start(intersect.getStart());
							builder.stop(intersect.getStop());
							builder.numMarkers(Array.getValuesBetween(positions, intersect.getStart(), intersect.getStop(), true).length);
							consensus.add(builder.build());
						}
					}
				}
			}
		}

		LocusSet<CNVariant> consensusSet = new LocusSet<CNVariant>(consensus.toArray(new CNVariant[consensus.size()]), true, log) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};
		return consensusSet;
	}

	public static class CNVCallResult {
		private LocusSet<CNVariant> chrCNVs;
		private LocusSet<CNVariant> chrCNVsReverseConsensus;
		private LocusSet<CNVariant> chrCNVsReverse;

		public CNVCallResult(LocusSet<CNVariant> chrCNVs, LocusSet<CNVariant> chrCNVsReverse, LocusSet<CNVariant> chrCNVsReverseConsensus) {
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
		 * Replaces NaN entries in either baf or lrr with 0 for both baf and lrr
		 */
		HANDLE_NAN, /**
		 * Does the PennCNV gc adjustment on the current data
		 */
		GC_ADJUST, /**
		 * Does the PennCNV median lrr and baf adjustment on the current data
		 */
		MEDIAN_ADJUST,

		/**
		 * Adjust the hmm model by the standard deviation of the current data
		 */
		ADJUST_HMM_SD;
	}

	/**
	 * @param proj
	 * @return the plus one index physical distance of each marker (of the analysis set) on each chromosome
	 */
	private static int[][] getSNPDist(Project proj, boolean reverse, boolean[] projectIndicesToUse) {

		MarkerSet markerSet = proj.getMarkerSet();
		int[][] chrPos = markerSet.getPositionsByChr();
		int projectIndex = 0;
		for (int i = 0; i < chrPos.length; i++) {
			ArrayList<Integer> updated = new ArrayList<Integer>();
			for (int j = 0; j < chrPos[i].length; j++) {
				if (projectIndicesToUse[projectIndex]) {
					updated.add(chrPos[i][j]);
				}
				projectIndex++;
			}
			chrPos[i] = Array.toIntArray(updated);
		}
		int[][] snpDists = new int[chrPos.length][];
		for (int i = 0; i < chrPos.length; i++) {
			int[] distsTmp = new int[chrPos[i].length];
			int[] posTmp = reverse ? Array.reverse(chrPos[i]) : chrPos[i];
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
	public static double[] adjustLrr(double[] lrrs, double minLrr, double maxLrr, boolean verbose, Logger log) {
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
	public static double[] adjustBaf(double[] bafs, double minBaf, double maxBaf, boolean debugMode, Logger log) {
		double[] adjusted = new double[bafs.length];
		ArrayList<Double> bafsToMedian = new ArrayList<Double>();
		for (int i = 0; i < bafs.length; i++) {
			if (!Double.isNaN(bafs[i]) && bafs[i] > minBaf && bafs[i] < maxBaf) {
				bafsToMedian.add(bafs[i]);
			}
		}
		double median = Array.median(Array.toDoubleArray(bafsToMedian));
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
	 * @return the {@link DATA_ADJUSTMENTS} enums in the order that should replicate PennCNV with gcCorrection
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
	 * @return the {@link DATA_ADJUSTMENTS} enums in the order that should replicate PennCNV without gcCorrection
	 */
	private static DATA_ADJUSTMENTS[] getPennCNVProcessingOrder() {
		ArrayList<DATA_ADJUSTMENTS> da = new ArrayList<CNVCaller.DATA_ADJUSTMENTS>();
		da.add(DATA_ADJUSTMENTS.HANDLE_NAN);
		da.add(DATA_ADJUSTMENTS.SUBSET_TO_ANALYSIS_MARKERS);
		da.add(DATA_ADJUSTMENTS.MEDIAN_ADJUST);
		da.add(DATA_ADJUSTMENTS.ADJUST_HMM_SD);

		return da.toArray(new DATA_ADJUSTMENTS[da.size()]);
	}

	private static CNVCallResult callCNVsFor(Project proj, PennHmm pennHmm, String sampleName, double[] sampLrrs, double[] sampBafs, GcModel gcModel, PFB pfb, MarkerSet markerSet, boolean[] markersToUse, boolean[] copyNumberDef, int[] chrsToCall, boolean callReverse, int numThreads, boolean debugMode) {
		DATA_ADJUSTMENTS[] dAdjustments = null;
		if (gcModel == null) {
			dAdjustments = getPennCNVProcessingOrder();
			proj.getLog().reportTimeWarning("No gc model was provided, calling on cnvs on raw data");

		} else {
			dAdjustments = getPennCNVGCProcessingOrder();
		}

		CNVCaller caller = new CNVCaller(proj, sampleName, pennHmm, gcModel, pfb, dAdjustments, markerSet, markersToUse, copyNumberDef, sampLrrs, sampBafs, debugMode);
		caller.adjustData();

		CNVCallResult cnvs = caller.callCNVS(chrsToCall, callReverse, numThreads);
		return cnvs;
	}

	/**
	 * @param proj
	 * @param pennHmm
	 *            required {@link PennHmm}
	 * @param sample
	 *            sample to call cnvs on
	 * @param gcModel
	 *            optional {@link GcModel}
	 * @param pfb
	 *            required {@link PFB}
	 * @param markerSet
	 * @param numThreads
	 *            analysis is split up by chromosome
	 * @param debugMode
	 *            report more values and warnings to compare with regular penncnv calls
	 * @return
	 */
	public static CNVCallResult callCNVsFor(Project proj, final PennHmm pennHmm, String sample, double[] lrrs, double[] bafs, final GcModel gcModel, final PFB pfb, final MarkerSet markerSet, int[] chrsToCall, boolean callReverse, int numThreads, boolean debugMode) {
		String[] markerNames = markerSet.getMarkerNames();
		boolean[] copyNumberDef = Array.booleanArray(markerNames.length, false);
		ARRAY array = proj.getArrayType();
		if (debugMode) {
			proj.getLog().reportTimeInfo("Assigning copy number probes according to " + array.toString() + " using the following " + Array.toStr(array.getCnFlags(), ","));
			proj.getLog().reportTimeInfo("BAF values greater than 1 will also be set to copy number only");

		}
		for (int i = 0; i < copyNumberDef.length; i++) {
			copyNumberDef[i] = array.isCNOnly(markerNames[i]) || bafs[i] > 1;
		}
		if (debugMode) {
			proj.getLog().reportTimeInfo("Found " + Array.booleanArraySum(copyNumberDef) + " copy number only markers");
		}
		int[] autosomalMarkers = proj.getAutosomalMarkerIndices();
		boolean[] markersToUse = Array.booleanArray(markerNames.length, false);
		for (int i = 0; i < autosomalMarkers.length; i++) {
			markersToUse[autosomalMarkers[i]] = true;
		}
		CNVCallResult cnvs = callCNVsFor(proj, pennHmm, sample, lrrs, bafs, gcModel, pfb, markerSet, markersToUse, copyNumberDef, chrsToCall, callReverse, numThreads, debugMode);
		return cnvs;
	}

	/**
	 * @author lane0212 Useful when calling cnvs across many samples
	 */
	private static class CNVProducer implements Producer<CNVCallResult> {
		private Project proj;
		private PennHmm pennHmm;
		private GcModel gcModel;
		private PFB pfb;
		private String[] samples;
		private int[] chrsToCall;
		private int numSampleThreads, index;
		private boolean debugMode, callReverse;

		/**
		 * @param proj
		 * @param pennHmm
		 * @param gcModel
		 * @param pfb
		 * @param samples
		 *            call cnvs on these samples
		 * @param numSampleThreads
		 * @param debugMode
		 */
		private CNVProducer(Project proj, PennHmm pennHmm, GcModel gcModel, PFB pfb, String[] samples, int[] chrsToCall, int numSampleThreads, boolean callReverse, boolean debugMode) {
			super();
			this.proj = proj;
			this.pennHmm = pennHmm;
			this.gcModel = gcModel;
			this.pfb = pfb;
			this.samples = samples;
			this.chrsToCall = chrsToCall;
			this.numSampleThreads = numSampleThreads;
			this.index = 0;
			this.callReverse = callReverse;
			this.debugMode = debugMode;
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
			final MarkerSet markerSet = proj.getMarkerSet();
			Callable<CNVCallResult> callable = new Callable<CNVCallResult>() {

				@Override
				public CNVCallResult call() throws Exception {
					Sample curSample = proj.getFullSampleFromRandomAccessFile(sample);
					CNVCallResult cnvs = callCNVsFor(proj, pennHmmTmp, curSample.getSampleName(), Array.toDoubleArray(curSample.getLRRs()), Array.toDoubleArray(curSample.getBAFs()), gcModelTmp, pfbTmp, markerSet, chrsToCall, callReverse, numSampleThreads, debugMode);
					return cnvs;
				}

			};
			index++;
			return callable;
		}

		@Override
		public void shutdown() {

		}
	}

	/**
	 * @param proj
	 * @param output
	 * @param numSampleThreads
	 *            number of samples analyzed at once.
	 * @param numChrThreads
	 *            number of chromosomes in each sample analyzed at once
	 * 
	 *            NOTE: total thread usage is numSampleThreads*numChrThreads
	 */
	public static void callCNVs(Project proj, String output, int numSampleThreads, int numChrThreads) {
		output = proj.PROJECT_DIRECTORY.getValue() + output;
		proj.getLog().reportTimeInfo("CNVS will be reported to " + output);
		PennHmm pennHmmOriginal = PennHmm.loadPennHmm(proj.HMM_FILENAME.getValue(), new Logger());
		PFB pfb = PFB.loadPFB(proj, proj.CUSTOM_PFB_FILENAME.getValue());
		GcModel gcModel = GcAdjustor.GcModel.populateFromFile(proj.GC_MODEL_FILENAME.getValue(false, false), false, proj.getLog());
		if (gcModel == null) {
			proj.getLog().reportTimeWarning("Calling cnvs without gc correction");
		}
		String[] samples = proj.getSamples();
		CNVProducer producer = new CNVProducer(proj, pennHmmOriginal, gcModel, pfb, samples, null, numChrThreads, false, true);
		WorkerTrain<CNVCallResult> train = new WorkerTrain<CNVCallResult>(producer, numSampleThreads, 2, proj.getLog());
		ArrayList<CNVariant> allCNVs = new ArrayList<CNVariant>();

		int index = 0;
		while (train.hasNext()) {
			index++;
			try {
				train.next().getChrCNVs().addAll(allCNVs);
			} catch (Exception e) {
				proj.getLog().reportTimeError("encountered problems calling cnvs for sample " + index + "\t" + samples[index]);
			}
			proj.getLog().reportTimeInfo("Called CNVs for" + index + " samples");

		}
		LocusSet<CNVariant> finalSet = new LocusSet<CNVariant>(allCNVs.toArray(new CNVariant[allCNVs.size()]), true, proj.getLog()) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};

		finalSet.writeRegions(proj.PROJECT_DIRECTORY.getValue() + output, TO_STRING_TYPE.REGULAR, true, proj.getLog());

	}

	public static void test() {
		Project proj = new Project("/home/pankrat2/lanej/projects/OSv2_hg19.properties", false);
		String hmm = proj.HMM_FILENAME.getValue();
		PennHmm pennHmmOriginal = PennHmm.loadPennHmm(hmm, new Logger());
		PFB pfb = PFB.loadPFB(proj, proj.CUSTOM_PFB_FILENAME.getValue());
		GcModel gcModel = GcAdjustor.GcModel.populateFromFile(proj.GC_MODEL_FILENAME.getValue(false, false), false, proj.getLog());
		int numThreads = 1;
		long trailerTime = System.currentTimeMillis();
		String[] samples = proj.getSamples();
		ArrayList<CNVariant> allCNVs = new ArrayList<CNVariant>();
		String[] sampTmp = new String[] { samples[ext.indexOfStr("7165764002_R06C02", samples)] };
		// String[] sampTmp= samples;
		CNVProducer producer = new CNVProducer(proj, pennHmmOriginal, null, pfb, sampTmp, null, 1, true, true);
		WorkerTrain<CNVCallResult> train = new WorkerTrain<CNVCallResult>(producer, numThreads, 2, proj.getLog());
		int index = 0;
		while (train.hasNext()) {
			index++;
			try {
				train.next().getChrCNVs().addAll(allCNVs);
			} catch (Exception e) {
				System.out.println("Error for sample " + index + "\t" + samples[index]);
				System.exit(1);
			}
			proj.getLog().reportTimeInfo("Called CNVs for" + index + " samples");

		}
		LocusSet<CNVariant> finalSet = new LocusSet<CNVariant>(allCNVs, true, proj.getLog()) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};
		finalSet.writeRegions(proj.PROJECT_DIRECTORY.getValue() + "cnvCaller.cnvs", TO_STRING_TYPE.REGULAR, true, proj.getLog());

		proj.getLog().reportTimeElapsed(trailerTime);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String output = "genvisis.cnvs";
		int numThreads = 24;

		String usage = "\n" + "cnv.hmm.CNVCaller requires 0-1 arguments\n";
		usage += "   (1) proj (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) output file (relative to project directory) (i.e. out=" + filename + " (default))\n" + "";
		usage += PSF.Ext.getNumThreadsCommand(3, numThreads);

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				output = ext.parseStringArg(args[i], "");

				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(args[i]);
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
			Project proj = new Project(filename, false);
			callCNVs(proj, output, numThreads, 1);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
