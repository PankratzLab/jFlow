package cnv.qc;

import java.util.ArrayList;

import cnv.filesys.MarkerData;
import cnv.filesys.Project;
import cnv.manage.MarkerDataLoader;
import common.Array;
import common.HashVec;

/**
 * Class to compute a QC metric based on B-allele MINOR allele frequency (BMAF) "min(meanBAF, 1-meanBAF)
 *
 */
public class CNVBMAF extends CNVBDeviation {

	public static final double DEFUALT_BAF_PENALTY = 0;
	public static final double DEFUALT_BAF_HET_LOWER = .15;
	public static final double DEFUALT_BAF_HET_UPPER = .85;

	/**
	 * Metric is divided by all markers
	 */
	public static final String SUMMARY_BY_NUM_ALL_MARKERS = "SUMMARY_BY_NUM_ALL_MARKERS";

	/**
	 * Metric is divided by sum of all markers
	 */
	public static final String SUMMARY_BY_SUM_BMAF_ALL_MARKERS = "SUMMARY_BY_BMAF_ALL_MARKERS";

	private ArrayList<Double> bmafAll;// keeping as ArrayList for now in case we want to use medians later
	private ArrayList<Double> bmafNonHet;
	private int countsHet;
	private double bmafMetric, percentHet;

	public CNVBMAF(String[] intensityOnlyFlags, double gcThreshold) {
		super(intensityOnlyFlags, gcThreshold);
		this.bmafAll = new ArrayList<Double>();
		this.bmafNonHet = new ArrayList<Double>();
		this.countsHet = 0;
		this.bmafMetric = 0;// TODO, could set to NaN instead
		this.percentHet = 1;
	}

	public void add(String markerName, byte genotype, float baf, double bmaf, float gc) {
		if (super.add(markerName, genotype, baf, gc)) {
			// checkBafForHet(baf)
			if (genotype != 1 && checkBafForHet(baf)) {
				bmafNonHet.add(bmaf);
			} else {
				countsHet++;
			}
			bmafAll.add(bmaf);
		}
	}

	/**
	 * @param bafPenalty
	 *            the number of heterozygous calls times this number will be subracted from the final metric. Set to 0 if no penalty is desired
	 * @param summaryType
	 *            see {@link #SUMMARY_BY_SUM_BMAF_ALL_MARKERS} and {@link #SUMMARY_BY_NUM_ALL_MARKERS} for options
	 */
	public void summarize(double bafPenalty, String summaryType) {
		super.summarize();
		if (bmafNonHet.size() > 0) {
			bmafMetric = Array.sum(Array.toDoubleArray(bmafNonHet));

			if (summaryType.equals(SUMMARY_BY_NUM_ALL_MARKERS)) {
				bmafMetric /= bmafAll.size();
			} else if (summaryType.equals(SUMMARY_BY_SUM_BMAF_ALL_MARKERS)) {
				bmafMetric /= Array.sum(Array.toDoubleArray(bmafAll));
			}
			percentHet = (double) countsHet / bmafAll.size();
		}
		bmafMetric -= (double) countsHet * bafPenalty;
	}

	public double getBmafMetric() {
		return bmafMetric;
	}

	public double getPercentHet() {
		return percentHet;
	}

	private static boolean checkBafForHet(float baf) {
		return baf < DEFUALT_BAF_HET_LOWER || baf > DEFUALT_BAF_HET_UPPER;
	}

	/**
	 * Helper class to facilitate a marker data based computation across all samples
	 *
	 */
	public static class PoplulationBAFs {
		private CNVBMAF[] cnvbmafs;

		public PoplulationBAFs(int numSamples, String[] intensityOnlyFlags, double gcThreshold) {
			this.cnvbmafs = initPopulation(numSamples, intensityOnlyFlags, gcThreshold);
		}

		public void add(String markerName, byte[] genotypes, float[] bafs, float[] gcs) {
			float bamf = (float) Array.median(Array.toDoubleArray(Array.removeNaN(bafs)));
			bamf = Math.min(bamf, 1 - bamf);
			for (int i = 0; i < cnvbmafs.length; i++) {
				cnvbmafs[i].add(markerName, genotypes[i], bafs[i], bamf, gcs[i]);
			}
		}

		public void summarize() {
			for (int i = 0; i < cnvbmafs.length; i++) {
				cnvbmafs[i].summarize(DEFUALT_BAF_PENALTY, SUMMARY_BY_NUM_ALL_MARKERS);
			}
		}

		public CNVBMAF[] getCnvbmafs() {
			return cnvbmafs;
		}

	}

	private static CNVBMAF[] initPopulation(int numSamples, String[] intensityOnlyFlags, double gcThreshold) {
		CNVBMAF[] cnvbmafs = new CNVBMAF[numSamples];
		for (int i = 0; i < cnvbmafs.length; i++) {
			cnvbmafs[i] = new CNVBMAF(intensityOnlyFlags, gcThreshold);
		}
		return cnvbmafs;
	}

	public static void test() {
		Project proj = new Project(null, false);
		String display = proj.getFilename(Project.DISPLAY_MARKERS_FILENAME);
		String[] markers = HashVec.loadFileToStringArray(display, false, new int[] { 0 }, true);
		MarkerDataLoader markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markers);
		PoplulationBAFs poplulationBDeviation = new PoplulationBAFs(proj.getSamples().length, DEFAULT_INTENSITY_ONLY_FLAGS, DEFAULT_GC_THRESHOLD);
		for (int i = 0; i < markers.length; i++) {
			MarkerData markerData = markerDataLoader.requestMarkerData(i);
			poplulationBDeviation.add(markers[i], markerData.getAbGenotypes(), markerData.getBAFs(), markerData.getGCs());

			markerDataLoader.releaseIndex(i);

		}
		poplulationBDeviation.summarize();
	}

}
