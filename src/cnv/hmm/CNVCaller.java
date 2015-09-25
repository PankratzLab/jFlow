package cnv.hmm;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.concurrent.Callable;

import common.Array;
import common.Logger;
import common.Positions;
import common.WorkerHive;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.hmm.PennHmm.ViterbiResult;
import cnv.qc.GcAdjustor;
import cnv.qc.GcAdjustor.GcModel;
import cnv.var.CNVariant;
import cnv.var.LocusSet;

/**
 * @author lane0212
 *
 *         Handles the data preparation and finalizations for calling cnvs via the PennCNV methods
 */
public class CNVCaller {
	private static final double MIN_LRR_MEDIAN_ADJUST = -2;
	private static final double MAX_LRR_MEDIAN_ADJUST = 2;
	private static final double MIN_BAF_MEDIAN_ADJUST = .25;
	private static final double MAX_BAF_MEDIAN_ADJUST = .75;
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
	}

	/**
	 * Warning: only the default PennCNV order of adjustments has been tested
	 */
	public void adjustData() {
		if (dataAdjustments == null) {
			proj.getLog().reportTimeWarning("No data adjustments supplied");
		} else {
			for (int i = 0; i < dataAdjustments.length; i++) {
				DATA_ADJUSTMENTS type = dataAdjustments[i];
				switch (type) {
				case ADJUST_HMM_SD:
					pennHmm = PennHmm.adjustBSD(pennHmm, Array.stdev(getValuesBetween(analysisLrrs, MIN_LRR_MEDIAN_ADJUST, MAX_LRR_MEDIAN_ADJUST)), proj.getLog());
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
					break;
				case MEDIAN_ADJUST:
					analysisLrrs = adjustLrr(analysisLrrs, MIN_LRR_MEDIAN_ADJUST, MAX_LRR_MEDIAN_ADJUST, proj.getLog());
					analysisBafs = adjustBaf(analysisBafs, MIN_BAF_MEDIAN_ADJUST, MAX_BAF_MEDIAN_ADJUST, proj.getLog());
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

	private LocusSet<CNVariant> callCNVS(int numThreads) {
		WorkerHive<LocusSet<CNVariant>> hive = new WorkerHive<LocusSet<CNVariant>>(numThreads, 10, proj.getLog());
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
		int[][] snpDists = getSNPDist(proj, finalAnalysisSet);
		for (int i = 0; i < snpDists.length; i++) {
			String chr = Positions.getChromosomeUCSC(i, true);
			if (snpDists[i].length > MIN_MARKERS_PER_CHROMOSOME && chrIndices.containsKey(chr)) {
				int[] currentChrIndices = Array.toIntArray(chrIndices.get(chr));
				int[] currentChrPositions = Array.subArray(analysisPositions, currentChrIndices);
				double[] currentChrLrr = Array.subArray(analysisLrrs, currentChrIndices);
				double[] currentChrBaf = Array.subArray(analysisBafs, currentChrIndices);
				double[] currentChrPfbs = Array.subArray(analysisPfbs, currentChrIndices);
				boolean[] currentChrCnDef = Array.subArray(copyNumberDef, currentChrIndices);
				CNVCallerWorker worker = new CNVCallerWorker(proj, dna, (byte) i, currentChrPositions, pennHmm, currentChrLrr, currentChrBaf, currentChrPfbs, snpDists[i], currentChrCnDef);
				hive.addCallable(worker);
			}
		}
		hive.execute(true);
		ArrayList<LocusSet<CNVariant>> results = hive.getResults();
		ArrayList<CNVariant> allCNVs = new ArrayList<CNVariant>();
		for (int i = 0; i < results.size(); i++) {
			for (int j = 0; j < results.get(i).getLoci().length; j++) {
				allCNVs.add(results.get(i).getLoci()[j]);
			}
		}

		LocusSet<CNVariant> allLocusSet = new LocusSet<CNVariant>(allCNVs, true, proj.getLog()) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};
		return allLocusSet;
	}

	/**
	 * The bare bones worker that should operate on any data passed to it
	 *
	 */
	private static class CNVCallerWorker implements Callable<LocusSet<CNVariant>> {
		private Project proj;
		private String dna;
		private byte currentChr;
		private int[] positions;
		private PennHmm pennHmm;
		private double[] lrrs;
		private double[] bafs;
		private double[] pfbs;
		private int[] snipDists;
		private boolean[] cnDef;

		private CNVCallerWorker(Project proj, String dna, byte currentChr, int[] positions, PennHmm pennHmm, double[] lrrs, double[] bafs, double[] pfbs, int[] snipDists, boolean[] cnDef) {
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
		}

		@Override
		public LocusSet<CNVariant> call() throws Exception {
			ViterbiResult viterbiResult = PennHmm.ViterbiLogNP_CHMM(pennHmm, lrrs, bafs, pfbs, snipDists, cnDef);
			LocusSet<CNVariant> chrCnvs = viterbiResult.analyzeStateSequence(proj, dna, dna, currentChr, positions, 2);
			if (chrCnvs.getLoci().length > 0) {
				chrCnvs = PennHmm.scoreCNVsSameChr(pennHmm, chrCnvs, positions, lrrs, bafs, pfbs, cnDef, proj.getLog());
			}
			return chrCnvs;
		}

	}

	/**
	 * Flags for the PennCNV data adjustments, can be used to assign order of operations as well;
	 *
	 */
	private enum DATA_ADJUSTMENTS {
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
	private static int[][] getSNPDist(Project proj, boolean[] projectIndicesToUse) {

		MarkerSet markerSet = proj.getMarkerSet();
		int[][] chrPos = markerSet.getPositionsByChr();
		int projectIndex = 0;
		for (int i = 0; i < chrPos.length; i++) {
			ArrayList<Integer> updated = new ArrayList<Integer>();
			for (int j = 0; j < chrPos[i].length; j++) {
				if (projectIndicesToUse[projectIndex]) {
					updated.add(chrPos[i][j]);
					projectIndex++;
				}
			}
			chrPos[i] = Array.toIntArray(updated);
		}
		int[][] snpDists = new int[chrPos.length][];
		for (int i = 0; i < chrPos.length; i++) {
			int[] distsTmp = new int[chrPos[i].length];
			if (distsTmp.length > 0) {
				distsTmp[chrPos[i].length - 1] = chrPos[i][chrPos[i].length - 1];
				for (int j = 0; j < distsTmp.length - 1; j++) {
					int dist = chrPos[i][j + 1] - chrPos[i][j];
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
	private static double[] adjustLrr(double[] lrrs, double minLrr, double maxLrr, Logger log) {
		double[] adjusted = new double[lrrs.length];
		double median = Array.median(getValuesBetween(lrrs, minLrr, maxLrr));
		log.reportTimeInfo("Median adjusting lrr values by " + median + "\t" + Array.median(Array.removeNaN(lrrs)));
		for (int i = 0; i < adjusted.length; i++) {
			adjusted[i] = lrrs[i] - median;
		}
		return adjusted;
	}

	/**
	 * Median adjust these baf values like PennCNV
	 */
	private static double[] adjustBaf(double[] bafs, double minBaf, double maxBaf, Logger log) {
		double[] adjusted = new double[bafs.length];
		ArrayList<Double> bafsToMedian = new ArrayList<Double>();
		for (int i = 0; i < bafs.length; i++) {
			if (!Double.isNaN(bafs[i]) && bafs[i] > minBaf && bafs[i] < maxBaf) {
				bafsToMedian.add(bafs[i]);
			}
		}
		double median = Array.median(Array.toDoubleArray(bafsToMedian));
		double factor = median - 0.5;
		log.reportTimeInfo("Median adjusting baf measures by " + factor);
		for (int i = 0; i < adjusted.length; i++) {
			if (!Double.isNaN(bafs[i]) && bafs[i] > minBaf && bafs[i] < maxBaf) {
				adjusted[i] = bafs[i] - factor;
			} else {
				adjusted[i] = bafs[i];
			}

		}
		return adjusted;
	}

	private static double[] getValuesBetween(double[] array, double min, double max) {
		ArrayList<Double> tmp = new ArrayList<Double>();
		for (int i = 0; i < array.length; i++) {
			if (!Double.isNaN(array[i]) && array[i] > min && array[i] < max) {
				tmp.add(array[i]);
			}
		}
		return Array.toDoubleArray(tmp);
	}

	/**
	 * @return the {@link DATA_ADJUSTMENTS} enums in the order that should replicate PennCNV
	 */
	private static DATA_ADJUSTMENTS[] getPennCNVProcessingOrder() {
		ArrayList<DATA_ADJUSTMENTS> da = new ArrayList<CNVCaller.DATA_ADJUSTMENTS>();
		da.add(DATA_ADJUSTMENTS.HANDLE_NAN);
		da.add(DATA_ADJUSTMENTS.GC_ADJUST);
		da.add(DATA_ADJUSTMENTS.SUBSET_TO_ANALYSIS_MARKERS);
		da.add(DATA_ADJUSTMENTS.MEDIAN_ADJUST);
		da.add(DATA_ADJUSTMENTS.ADJUST_HMM_SD);

		return da.toArray(new DATA_ADJUSTMENTS[da.size()]);
	}

	// public CNVCaller(Project proj, String dna, PennHmm pennHmm, GcModel gcModel, PFB pfb, DATA_ADJUSTMENTS[] dataAdjustments, MarkerSet markerSet, boolean[] markersToUse, boolean[] copyNumberDef, double[] lrrs, double[] bafs, boolean debugMode) {

	public static LocusSet<CNVariant> callCNVsFor(Project proj, PennHmm pennHmm, String sampleName, double[] sampLrrs, double sampBafs, GcModel gcModel, PFB pfb, MarkerSet markerSet, boolean[] markersToUse, boolean[] copyNumberDef, boolean debugMode) {
	//	CNVCaller caller = new CNVCaller(proj, sampleName, pennHmm, gcModel, pfb, getPennCNVProcessingOrder(), markerSet, markersToUse, copyNumberDef, sampLrrs, sampBafs, debugMode);
		return null;
	}

}
