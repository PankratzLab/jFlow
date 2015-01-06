package cnv.qc;

import java.util.ArrayList;

import cnv.filesys.MarkerData;
import cnv.filesys.Project;
import cnv.manage.MarkerDataLoader;
import common.Array;
import common.HashVec;
import common.ext;

/**
 * Class for computing the B-Deviation for markers in a cnv call, similar to http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2963811/pdf/1469.pdf <br>
 * Currently this is designed to compute markers on the fly...and does not maintain a copy of the input data <br>
 * Note: Since affy copy number probes do not have a meaningful baf, we deviate from the above method by skipping all intensity only flags if defined<br>
 * Also, we take the median b-deviation as opposed to the mean
 */
public class CNVBDeviation {

	/**
	 * These are special, CN_ is Affy and cnv is Illumina
	 */
	public static final String[] DEFAULT_INTENSITY_ONLY_FLAGS = { "CN_", "cnv" };
	public static final float DEFAULT_GC_THRESHOLD = 0;

	private String[] intensityOnlyFlags;
	private double gcThreshold, medianBDeviationAll, medianBDeviationHet;
	private ArrayList<Double> bDeviationsAll;// all markers
	private ArrayList<Double> bDeviationsHet;// heterozygous markers

	/**
	 * @param intensityOnlyFlags
	 *            markers that contain these flags will be skipped
	 * @param gcThreshold
	 *            markers with gcs < than this will be skipped
	 */
	public CNVBDeviation(String[] intensityOnlyFlags, double gcThreshold) {
		this.intensityOnlyFlags = intensityOnlyFlags;
		this.gcThreshold = gcThreshold;
		this.medianBDeviationAll = 0;
		this.medianBDeviationHet = 0;
		this.bDeviationsAll = new ArrayList<Double>();
		this.bDeviationsHet = new ArrayList<Double>();
	}

	/**
	 * @param markerName
	 *            used for intensity only determination
	 * @param genotype
	 * @param baf
	 * @param gc
	 *            to filter calls if desired
	 * @return true if the marker was added to the b-deviation container
	 */
	public boolean add(String markerName, byte genotype, float baf, float gc) {
		boolean added = true;
		double dBaf = baf;

		if (Float.isNaN(baf) || ext.indexOfStartsWith(markerName, intensityOnlyFlags, false) > 0 || gc < gcThreshold) {// skip intensity only, and gc filtered
			added = false;
		} else {
			double currentBDeviation = Double.NaN;

			if (genotype == 1) {
				currentBDeviation = Math.abs(dBaf - 0.5);
				bDeviationsHet.add(currentBDeviation);
			} else if (genotype == 0 || genotype == 2) {
				currentBDeviation = Math.min(baf, 1 - dBaf);
			} else {
				currentBDeviation = Math.min(Math.min(baf, 1 - baf), Math.abs(baf - 0.5));
			}
			currentBDeviation = Math.sqrt(currentBDeviation);
			if (!Double.isNaN(currentBDeviation)) {
				bDeviationsAll.add(currentBDeviation);
			} else {
				added = false;
			}
		}
		return added;
	}

	public void summarize() {
		if (bDeviationsAll.size() > 0) {
			double[] bDeviationsDAll = Array.toDoubleArray(bDeviationsAll);
			medianBDeviationAll = Array.median(bDeviationsDAll);
		}
		if (bDeviationsHet.size() > 0) {
			double[] bDeviationsDHet = Array.toDoubleArray(bDeviationsHet);
			medianBDeviationHet = Array.median(bDeviationsDHet);

		}
	}

	public double getMedianBDeviationAll() {
		return medianBDeviationAll;
	}

	public double getMedianBDeviationHet() {
		return medianBDeviationHet;
	}

	/**
	 * Helper class to facilitate a marker data based computation across all samples
	 *
	 */
	public static class PoplulationBDeviation {
		private CNVBDeviation[] cnvbDeviations;

		public PoplulationBDeviation(int numSamples, String[] intensityOnlyFlags, double gcThreshold) {
			this.cnvbDeviations = initPopulation(numSamples, intensityOnlyFlags, gcThreshold);

		}

		public void add(String markerName, byte[] genotypes, float[] bafs, float[] gcs) {
			for (int i = 0; i < cnvbDeviations.length; i++) {
				cnvbDeviations[i].add(markerName, genotypes[i], bafs[i], gcs[i]);
			}
		}

		public void summarize() {
			for (int i = 0; i < cnvbDeviations.length; i++) {
				cnvbDeviations[i].summarize();
			}
		}

		public CNVBDeviation[] getCnvbDeviations() {
			return cnvbDeviations;
		}

	}

	private static CNVBDeviation[] initPopulation(int numSamples, String[] intensityOnlyFlags, double gcThreshold) {
		CNVBDeviation[] cnvbDeviations = new CNVBDeviation[numSamples];
		for (int i = 0; i < cnvbDeviations.length; i++) {
			cnvbDeviations[i] = new CNVBDeviation(intensityOnlyFlags, gcThreshold);
		}
		return cnvbDeviations;
	}

	public static void test() {
		Project proj = new Project(null, false);
		String display = proj.getFilename(Project.DISPLAY_MARKERS_FILENAME);
		String[] markers = HashVec.loadFileToStringArray(display, false, new int[] { 0 }, true);
		MarkerDataLoader markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markers);
		PoplulationBDeviation poplulationBDeviation = new PoplulationBDeviation(proj.getSamples().length, DEFAULT_INTENSITY_ONLY_FLAGS, DEFAULT_GC_THRESHOLD);
		for (int i = 0; i < markers.length; i++) {
			MarkerData markerData = markerDataLoader.requestMarkerData(i);
			poplulationBDeviation.add(markers[i], markerData.getAbGenotypes(), markerData.getBAFs(), markerData.getGCs());

			markerDataLoader.releaseIndex(i);

		}
		poplulationBDeviation.summarize();

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		test();
	}

}
