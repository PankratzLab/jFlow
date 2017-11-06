package org.genvisis.gwas.results;

import org.genvisis.common.ArrayUtils;

/**
 * Basic implementation of {@code IndexFactorsParser} that subsets the input data line to the
 * desired columns.
 */
public class SimpleResultsParser extends IndexFactorsParser {
	String[] outputHeader;

	public SimpleResultsParser(String outputDelim, String[][] desiredHeaders, String[] outputHeader) {
		super(outputDelim, desiredHeaders);
		this.outputHeader = outputHeader;
	}

	@Override
	public String[] getOutputHeader() {
		return outputHeader;
	}

	@Override
	public String[] parseInputLine(String[] inputLine) {
		return ArrayUtils.subArray(inputLine, getIndicesOfInputs());
	}

}
