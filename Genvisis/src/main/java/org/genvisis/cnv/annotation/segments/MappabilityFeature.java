/**
 * 
 */
package org.genvisis.cnv.annotation.segments;

import org.genvisis.common.Logger;
import org.genvisis.seq.manage.BEDFileReader.BEDFeatureSeg;

import htsjdk.tribble.bed.BEDFeature;

/**
 * @author Kitty
 *
 */
public class MappabilityFeature extends BEDFeatureSeg {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private BEDFeatureSeg bFeatureSeg;

	/**
	 * @param bedFeature constructs a {@link BEDFeatureSeg}
	 * @param log
	 */
	public MappabilityFeature(BEDFeature bedFeature, Logger log) {
		super(bedFeature, log);
	}

	/**
	 * @param log
	 * @return mappability, as stored and returned by {@link BEDFeature#getName()}, confusing, yes -
	 *         but seem to be the map track format
	 */
	public double getMappability(Logger log) {
		double mapScore = Double.NaN;
		String tmpScore = bFeatureSeg.getBedFeature().getName();

		try {
			mapScore = Double.parseDouble(tmpScore);

		} catch (NumberFormatException nfe) {
			String error = "Could not convert " + tmpScore + " to mappability score";
			log.reportError(error);
			throw new IllegalArgumentException(error);
		}
		if (mapScore < 0) {
			String error = "Could not convert " + tmpScore + " to mappability score, negative mapScore";
			log.reportError(error);
			throw new IllegalArgumentException(error);
		}
		return mapScore;

	}

}
