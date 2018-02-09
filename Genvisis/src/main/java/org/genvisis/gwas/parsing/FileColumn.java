package org.genvisis.gwas.parsing;

import java.util.Map;

public interface FileColumn<T> {
	/**
	 * Get the name of the column. Not expected to be unique. Used if writing output and to link
	 * columns between files.
	 * 
	 * @return The column name
	 */
	String getName();

	/**
	 * Initialize this FileColumn. Repeated calls to this method should be no-ops.
	 * 
	 * @param headerMap Header of the file mapped from String value to index. May be {@code null} if
	 *        the file does not have a header - if this {@link FileColumn} is expecting a header and
	 *        {@code null} is given, an exception is expected to occur.
	 */
	void initialize(Map<String, Integer> headerMap);

	/**
	 * Get this FileColumn's value from a line of data.
	 * 
	 * @param line {@code String[]} of line data
	 * @return {@code <T>}
	 * @throws ParseFailureException thrown if there is a
	 */
	T getValue(String[] line) throws ParseFailureException;

	public boolean dieOnParseFailure();

}
