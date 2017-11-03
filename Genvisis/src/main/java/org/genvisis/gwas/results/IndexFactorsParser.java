package org.genvisis.gwas.results;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.ext;

/**
 * Intermediate implementation of a {@code ResultFormatParser} that uses
 * {@code ext#indexFactors(String[], String[], boolean)} to determine where desired data is in the
 * input file header.
 */
public abstract class IndexFactorsParser extends ResultFormatParser {

	/**
	 * Aliases of desired columns
	 */
	String[][] queries;
	/**
	 * Indices of desired columns
	 */
	int[] indices;
	boolean dieIfMissing = true;

	public IndexFactorsParser(String outputDelim, String[][] desiredHeaders) {
		super(outputDelim);
		this.queries = desiredHeaders;
	}

	public IndexFactorsParser(String outputDelim, String[][] desiredHeaders, boolean dieIfAnyMissing) {
		this(outputDelim, desiredHeaders);
		dieIfMissing = dieIfAnyMissing;
	}

	@Override
	public void setInputHeader(String[] inputHeader) {
		super.setInputHeader(inputHeader);
		// search for the desired column names
		this.indices = ext.indexFactors(queries, inputHeader, false, true, true);
		if (dieIfMissing) {
			// throw an exception if any are missing
			for (int i = 0; i < indices.length; i++) {
				if (indices[i] == -1) {
					throw new IllegalStateException("Couldn't find header column " + i + ", alises: {"
																					+ ArrayUtils.toStr(queries[i], ",") + "}.");
				}
			}
		}
	}

	public int[] getIndicesOfInputs() {
		return indices;
	}

}
