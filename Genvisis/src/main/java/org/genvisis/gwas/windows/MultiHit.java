package org.genvisis.gwas.windows;

import java.util.Map;

/**
 * {@link Hittable} implementation that holds multiple p-values for a single position
 */
public class MultiHit extends AbstractHittable {
	private final Map<String, Double> pVals;

	public MultiHit(String name, byte chr, int pos, Map<String, Double> pVals) {
		super(name, chr, pos);
		this.pVals = pVals;
	}

	/**
	 * return The p-value for the given p-value column label
	 */
	public double pVal(String pVal) {
		return pVals.get(pVal);
	}

	@Override
	public double getPval() {
		throw new UnsupportedOperationException(getClass().toString() + " has multiple p-values");
	}
}
