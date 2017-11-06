package org.genvisis.gwas.results;

public abstract class ResultFormatParser {
	/**
	 * Output file delimiter
	 */
	String outputDelim;
	/**
	 * Input file header
	 */
	String[] inputHeader;

	public ResultFormatParser(String outDelim) {
		this.outputDelim = outDelim;
	}

	/**
	 * Used by {@code ResultsPackager2}.
	 * 
	 * @param inputHeader Input header values
	 */
	public void setInputHeader(String[] inputHeader) {
		this.inputHeader = inputHeader;
	}

	public String getOutputDelim() {
		return outputDelim;
	};

	/**
	 * @return the output header, which can be dynamically-created based on the input header values
	 */
	public abstract String[] getOutputHeader();

	/**
	 * Parse a line of the input file and return a line of output data
	 * 
	 * @param inputLine
	 * @return
	 */
	public abstract String[] parseInputLine(String[] inputLine);

}
