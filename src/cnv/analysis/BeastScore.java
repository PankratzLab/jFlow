package cnv.analysis;

import java.util.ArrayList;

import cnv.manage.Transforms;
import cnv.var.CNVariant;
import common.Array;
import common.Logger;
import filesys.Segment;

/**
 * Class to compute beast scores similar to the beast algorithm http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0086272.
 * Values are always inverse normalized with 5df, and scaled by the empirically derived SCALE_FACTOR_MAD. NaN data is discourage, but is ignored;
 * Note: Since the input data is always inverse normalized, the MAD for a region (such as a cnv) is the median of the (absolute value) inverse normalized data across that region (scaled by the MAD of the superset i.e the chromosome)
 * Note: one difference to the score scheme is that a min/max length parameter has not been added. i.e(length<=min -> score=-, length>=max -> height goes to min height. This may make long cnv calls have inflated scores. 
 */
public class BeastScore {
	/**
	 * The scale factor was empirically derived and scales all MAD values for consistent comparisons across samples
	 */
	public static final double SCALE_FACTOR_MAD = 0.134894516;
	public static final float DEFAULT_ALPHA = 0.5f;
	private float[] inputData, beastHeights, beastScores, invTransScaledStdev;
	private int[][] indicesToChunk;
	private int[][] indicesForScores;
	private int[] beastLengths;
	private Logger log;

	/**
	 * @param inputData 
	 * 			  Most often LRR values for a single sample
	 * @param indicesToChunk
	 *            Most often the indices of chromosomes
	 * @param indicesForScores
	 *            Most often the indices of markers in cnvs for a particular sample
	 * @param log 
	 * 			  You know, a log
	 */
	public BeastScore(float[] inputData, int[][] indicesToChunk, int[][] indicesForScores, Logger log) {
		super();
		this.inputData = inputData;
		this.indicesToChunk = indicesToChunk;
		this.indicesForScores = indicesForScores;
		this.beastHeights = new float[indicesForScores.length];
		this.beastScores = new float[indicesForScores.length];
		this.beastLengths = new int[indicesForScores.length];
		this.log = log;
	}

	public float[] getBeastScores() {
		return beastScores;
	}

	public float[] getBeastHeights() {
		return beastHeights;
	}

	public int[] getBeastLengths() {
		return beastLengths;
	}

	public float[] getInvTransScaledStdevChunks() {
		return invTransScaledStdev;
	}

	/**
	 * Computes beast score using the default alpha
	 */
	public void computeBeastScores() {
		computeBeastScores(DEFAULT_ALPHA);
	}

	/**
	 * Inverse transforms data, scales chunks (usually chromosomes) by MAD scale factor, uses scaled chromosome values to obtain scaled MAD values for input data
	 * @param alpha
	 * @param computeSTDevRaw
	 *            compute standard deviation of each chunk corresponding to the indicesForScores from raw data
	 * @param computeSTDevtransformed
	 *            compute standard deviation of each chunk corresponding to the indicesForScores from inverseTransformed and scaled data
	 * 
	 */
	public void computeBeastScores(float alpha) {
		float[] inverseTransformedData = transformData(inputData, indicesToChunk, log);
		float[] indicesMADScaled = getscaleMADIndices(indicesToChunk, inverseTransformedData, log);
		float[] inverseTransformedDataScaleMAD = getscaleMADData(inverseTransformedData, indicesToChunk, indicesMADScaled, log);
		this.beastHeights = getBeastHeights(inverseTransformedDataScaleMAD, indicesForScores, log);
		this.beastLengths = getBeastLengths(inverseTransformedDataScaleMAD, indicesForScores, log);
		this.beastScores = getBeastScores(beastHeights, beastLengths, alpha, log);
	}

	/**
	 * 
	 * @param inputData
	 * @param indicesToTransform
	 * @return inputData inverse Transformed with 5 df according to the indicesToChunk
	 */
	public static float[] transformData(float[] inputData, int[][] indicesToChunk, Logger log) {
		return Transforms.transform(inputData, 2, indicesToChunk, Array.booleanArray(indicesToChunk.length, true));
	}

	/**
	 * 
	 * Computes the median across indicesToScale of inverseTransformedData and scales by SCALE_FACTOR_MAD
	 * 
	 * @param indicesToScale
	 * @param inverseTransformedData
	 * @param log
	 * @return
	 */
	public static float[] getscaleMADIndices(int[][] indicesToScale, float[] inverseTransformedData, Logger log) {
		float[] indicesMADScaled = new float[indicesToScale.length];
		for (int i = 0; i < indicesToScale.length; i++) {
			if (indicesToScale[i] != null && indicesToScale[i].length > 0) {
				ArrayList<Float> medianIndices = new ArrayList<Float>();
				for (int j = 0; j < indicesToScale[i].length; j++) {
					if (!Double.isNaN(inverseTransformedData[j])) {
						medianIndices.add(Math.abs(inverseTransformedData[j]));
					}
				}
				indicesMADScaled[i] = (float) (Array.median(Array.toDoubleArray(Array.toFloatArray(medianIndices))) / SCALE_FACTOR_MAD);
			}
		}
		return indicesMADScaled;
	}

	/**
	 * 
	 * scales inverseTransformedData according to indicesToScale and by the factors in scaleMAD
	 * 
	 * @param inverseTransformedData
	 * @param indicesToScale
	 * @param scaleMAD
	 * @param log
	 * @return the inverse transformed data scaled by the MAD according to indicesToScale
	 */
	public static float[] getscaleMADData(float[] inverseTransformedData, int[][] indicesToScale, float[] scaleMAD, Logger log) {
		if (scaleMAD.length != indicesToScale.length) {
			log.reportError("Error - the indices to scale and the factors to scale by must have the same length");
			return null;
		}
		float[] inverseTransformedDataScaleMAD = new float[inverseTransformedData.length];
		for (int i = 0; i < indicesToScale.length; i++) {
			if (indicesToScale != null && indicesToScale[i].length > 0) {
				for (int j = 0; j < indicesToScale[i].length; j++) {
					inverseTransformedDataScaleMAD[indicesToScale[i][j]] = (float) (inverseTransformedData[indicesToScale[i][j]] / (scaleMAD[i]));
				}
			} else {
				log.reportError("Warning - the index " + i + " was missing data");
			}
		}
		return inverseTransformedDataScaleMAD;
	}

	/**
	 * @param inverseTransformedDataScaleMAD
	 * @param indicesForHeights
	 * @return the median of inverseTransformedDataScaleMAD according to the indices in indicesForHeights
	 */
	public static float[] getBeastHeights(float[] inverseTransformedDataScaleMAD, int[][] indicesForHeights, Logger log) {
		float[] beastHeights = new float[indicesForHeights.length];
		for (int i = 0; i < indicesForHeights.length; i++) {
			ArrayList<Float> medianHeightIndices = new ArrayList<Float>();
			for (int j = 0; j < indicesForHeights[i].length; j++) {
				if (!Float.isNaN(inverseTransformedDataScaleMAD[indicesForHeights[i][j]])) {
					medianHeightIndices.add(inverseTransformedDataScaleMAD[indicesForHeights[i][j]]);
				}
			}
			if (medianHeightIndices.size() > 0) {
				beastHeights[i] = (float) (Array.median(Array.toDoubleArray(Array.toFloatArray(medianHeightIndices))));
			} else {
				beastHeights[i] = Float.NaN;
			}
		}
		return beastHeights;
	}

	/**
	 * @param inverseTransformedDataScaleMAD
	 * @param indicesForLengths
	 * @return the number of non - NaN inputs in inverseTransformedDataScaleMAD according to the indices in indicesForLengths
	 */
	public static int[] getBeastLengths(float[] inverseTransformedDataScaleMAD, int[][] indicesForLengths, Logger log) {
		int[] beastLengths = new int[indicesForLengths.length];
		for (int i = 0; i < indicesForLengths.length; i++) {
			beastLengths[i] = 0;
			for (int j = 0; j < indicesForLengths[i].length; j++) {
				if (!Float.isNaN(inverseTransformedDataScaleMAD[indicesForLengths[i][j]])) {
					beastLengths[i]++;
				}
			}
		}
		return beastLengths;
	}

	/**
	 * @param beastHeights
	 * @param beastLengths
	 * @param alpha
	 * @param log
	 * @return the beast scores for each index of beastHeights and beastLengths
	 */
	public static float[] getBeastScores(float[] beastHeights, int[] beastLengths, float alpha, Logger log) {
		float[] beastScores = new float[beastHeights.length];
		if (beastHeights.length != beastLengths.length) {
			log.reportError("Error - heights and lengths must contain the same number of inputs to compute a Beast Score");
			return null;
		} else {
			for (int i = 0; i < beastHeights.length; i++) {
				beastScores[i] = scoreBeast(beastLengths[i], alpha, beastHeights[i], log);
			}
		}
		return beastScores;
	}

	/**
	 * @param length
	 * @param alpha
	 * @param height
	 * @param log
	 * @return the beast score using using length^alpha * height
	 * We do not implement a max/min height
	 */
	public static float scoreBeast(int length, float alpha, float height, Logger log) {
		return (float) Math.abs(Math.pow(length, alpha) * height);
	}

	/**
	 * Helper function to extract indices of markers contained in a CNVariant
	 * 
	 * @param chr
	 * @param positions
	 * @param cnVariant
	 * @param log
	 * @return
	 */
	public static int[] getCNVMarkerIndices(byte[] chr, int[] positions, CNVariant cnVariant, Logger log) {
		int[] indices = new int[cnVariant.getNumMarkers()];
		int count = 0;
		if (chr.length != positions.length) {
			log.reportError("Error - the chromosome and position arrays must be the same length");
			return null;
		}
		for (int i = 0; i < chr.length; i++) {
			if (cnVariant.overlaps(new Segment(chr[i], positions[i], positions[i]))) {
				if (count >= indices.length) {
					log.reportError("Error - found too many markers for cnVariant " + cnVariant.toString() + "\n Perhaps some markers were filtered out prior to calling?");
				}
				indices[count] = i;
				count++;
			}
		}
		if (indices.length < count) {
			log.reportError("Error - not enough markers were found for cnVariant" + cnVariant.toString() + " using the current chr/position data");
			return null;
		}
		return indices;
	}
}


//public static float[] getStandardDeviationForTargets(int[][] indicesToStDev, int[][] indicesForStdev, float[] data, Logger log) {
//	float[] stDev = getSTDevIndices(indicesToStDev, data, log);
//	float[] stDevData = assignStdevToDataAllIndices(data, indicesToStDev, stDev, log);
//	return getStDevByTarget(stDevData, indicesForStdev, log);
//}
//
///**
// * Compute the standard deviation of the data across the indices according to indicesToStDev (usually across chromosomes)
// * 
// * @param indicesToStDev
// * @param data
// * @param log
// * @return
// */
//public static float[] getSTDevIndices(int[][] indicesToStDev, float[] data, Logger log) {
//	float[] stDev = new float[indicesToStDev.length];
//	for (int i = 0; i < indicesToStDev.length; i++) {
//		if (indicesToStDev[i] != null && indicesToStDev[i].length > 0) {
//			ArrayList<Float> stDevIndices = new ArrayList<Float>();
//			for (int j = 0; j < indicesToStDev[i].length; j++) {
//				if (!Double.isNaN(data[j])) {
//					stDevIndices.add(data[j]);
//				}
//			}
//			stDev[i] = Array.stdev(Array.toFloatArray(stDevIndices), true);
//		}
//	}
//	return stDev;
//}
//
///**
// * 
// * Assigns standard standard deviations computed from the global regions (usually chromosomes) to the indices of the full data array as specified by indicesToScale
// * 
// * @param data
// * @param indicesToScale
// * @param stDev
// * @param log
// * @return
// */
//public static float[] assignStdevToDataAllIndices(float[] data, int[][] indicesToScale, float[] stDev, Logger log) {
//	if (stDev.length != indicesToScale.length) {
//		log.reportError("Error - the indices to scale and the factor to scale by must have the same length");
//		return null;
//	}
//	float[] stDevData = new float[data.length];
//	for (int i = 0; i < indicesToScale.length; i++) {
//		if (indicesToScale != null && indicesToScale[i].length > 0) {
//			for (int j = 0; j < indicesToScale[i].length; j++) {
//				stDevData[indicesToScale[i][j]] = stDev[i];
//			}
//		}
//	}
//	return stDevData;
//}
//
///**
// * Assigns standard deviations of the global region (usually chromosome) to the target regions (usually cnvs)
// * 
// * @param stDevAllData
// * @param indicesForStdev
// * @param log
// * @return
// */
//public static float[] getStDevByTarget(float[] stDevAllData, int[][] indicesForStdev, Logger log) {
//	float[] stDev = new float[indicesForStdev.length];
//	for (int i = 0; i < indicesForStdev.length; i++) {
//		stDev[i] = stDevAllData[indicesForStdev[i][0]];
//	}
//	return stDev;
//}
