package org.genvisis.cnv.qc;

import java.util.ArrayList;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.ext;

import com.google.common.primitives.Doubles;

/**
 * Class for computing the B-Deviation for markers in a cnv call, similar to
 * http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2963811/pdf/1469.pdf <br>
 * Currently this is designed to compute markers on the fly...and does not maintain a copy of the
 * input data <br>
 * Note: Since affy copy number probes do not have a meaningful baf, we deviate from the above
 * method by skipping all intensity only flags if defined<br>
 * Also, we take the median b-deviation as opposed to the mean
 *
 */
public class CNVBDeviation {

	/**
	 * These are special, CN_ is Affy and cnv is Illumina
	 */
	public static final String[] DEFAULT_INTENSITY_ONLY_FLAGS = {"CN_", "cnv"};
	public static final float DEFAULT_GC_THRESHOLD = 0;

	private final String[] intensityOnlyFlags;
	private final double gcThreshold;
	private double medianBDeviationAll;
	private double medianBDeviationHet;
	private final ArrayList<Double> bDeviationsAll;// all markers b-deviation
	private final ArrayList<Double> bDeviationsHet;// heterozygous markers b-deviation

	/**
	 * @param intensityOnlyFlags markers that contain these flags will be skipped
	 * @param gcThreshold markers with gcs < than this will be skipped
	 */
	public CNVBDeviation(String[] intensityOnlyFlags, double gcThreshold) {
		this.intensityOnlyFlags = intensityOnlyFlags;
		this.gcThreshold = gcThreshold;
		// TODO, set the 0s to NaN?
		medianBDeviationAll = 0;
		medianBDeviationHet = 0;
		bDeviationsAll = new ArrayList<Double>();
		bDeviationsHet = new ArrayList<Double>();
	}

	/**
	 * @param markerName used for intensity only determination
	 * @param genotype
	 * @param baf
	 * @param gc to filter calls if desired
	 * @return true if the marker was added to the b-deviation container
	 */
	public boolean add(String markerName, byte genotype, float baf, float gc) {
		boolean added = true;
		double dBaf = baf;

		if (Float.isNaN(baf)	|| ext.indexOfStartsWith(markerName, intensityOnlyFlags, false) > 0
				|| gc < gcThreshold) {// skip intensity only, and gc filtered
			added = false;
		} else {
			double currentBDeviation = Double.NaN;// metric from paper

			if (genotype == 1) {
				currentBDeviation = Math.abs(dBaf - 0.5);
				bDeviationsHet.add(currentBDeviation);
			} else if (genotype == 0 || genotype == 2) {
				currentBDeviation = Math.min(baf, 1 - dBaf);

			} else {// no call
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
			double[] bDeviationsDAll = Doubles.toArray(bDeviationsAll);
			medianBDeviationAll = ArrayUtils.median(bDeviationsDAll);
		}
		if (bDeviationsHet.size() > 0) {
			double[] bDeviationsDHet = Doubles.toArray(bDeviationsHet);
			medianBDeviationHet = ArrayUtils.median(bDeviationsDHet);
		}
	}

	public double getMedianBDeviationAll() {
		return medianBDeviationAll;
	}

	public double getMedianBDeviationHet() {
		return medianBDeviationHet;
	}

}
