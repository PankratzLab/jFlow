package org.genvisis.gwas.results.files;

import java.io.FileNotFoundException;

import org.genvisis.common.Files;

/**
 * Class for creating and configuring {@link FileParser} objects.
 */
public abstract class AbstractFileParserFactory {

	protected FileParser parser;

	protected AbstractFileParserFactory(String inputFile, String inputDelim) {
		if (!Files.exists(inputFile)) {
			throw new IllegalArgumentException(new FileNotFoundException("File missing: " + inputFile));
		}
		this.parser = new FileParser(inputFile, inputDelim);
	}

	/**
	 * Add a filter to drop lines, specifying whether or not to throw an exception if a value fails.
	 * 
	 * @param filter {@link ColumnFilter}
	 * @param dieIfFail
	 * @return
	 */
	public AbstractFileParserFactory filter(ColumnFilter filter, boolean dieIfFail) {
		this.parser.addFilter(filter, dieIfFail);
		return this;
	}

	/**
	 * Add a filter to drop lines, without dying if filtering rejects a line.
	 * 
	 * @param filter {@link ColumnFilter}
	 * @return this
	 */
	public AbstractFileParserFactory filter(ColumnFilter filter) {
		return filter(filter, false);
	}

	/**
	 * Skip a specified number of lines at the beginning of the file.<br />
	 * Multiple calls to this method will overwrite the previous value.
	 * 
	 * @param num Number of lines to skip
	 * @return this
	 */
	public AbstractFileParserFactory skipLines(int num) {
		this.parser.setSkippedLines(num);
		return this;
	}

	/**
	 * Skip all lines that start with the specified prefix.<br />
	 * Can be called multiple times with different prefices.
	 * 
	 * @param prefix Line skip prefix
	 * @return this
	 */
	public AbstractFileParserFactory skipPrefix(String prefix) {
		this.parser.addSkippedPrefix(prefix);
		return this;
	}

	/**
	 * Ignore lines with an incorrect number of columns. Default behavior is to throw an exception.
	 * 
	 * @return this
	 */
	public AbstractFileParserFactory skipLinesWithIncorrectNumberOfColumns() {
		this.parser.skipLineWithIncorrectNumberOfColumns();
		return this;
	}

	/**
	 * Throw an exception if blank lines are encountered while reading the file. Default behavior is
	 * to skip blank lines.
	 * 
	 * @return this
	 */
	public AbstractFileParserFactory failOnBlankLines() {
		this.parser.failOnBlankLines();
		return this;
	}

	/**
	 * Flag the parser to not expect a header. May cause failure if data columns expect headers.
	 * 
	 * @return this
	 */
	public AbstractFileParserFactory noHeader() {
		this.parser.setNoHeader();
		return this;
	}

	/**
	 * Flag the parser to not trim input lines.
	 * 
	 * @return this
	 */
	public AbstractFileParserFactory noTrim() {
		this.parser.setNoTrim();
		return this;
	}

	/**
	 * Set the value to be used if parsing a value fails.
	 * 
	 * @param value String
	 * @return this
	 */
	public AbstractFileParserFactory parseFailValue(String value) {
		this.parser.setParseFailValue(value);
		return this;
	}

	/**
	 * Check that the key columns in the file link share a name with the data loaded in the primary
	 * parser.
	 * 
	 * @param link Link to check
	 */
	private void checkLink(FileLink link) {
		for (FileColumn<?> fc : link.getKeys()) {
			boolean found = false;
			for (FileColumn<?> fc1 : parser.getDataColumnsInOrder()) {
				if (fc1.getName().equals(fc.getName())) {
					found = true;
					break;
				}
			}
			for (FileColumn<?> fc1 : parser.getAddlDataColumns()) {
				if (fc1.getName().equals(fc.getName())) {
					found = true;
					break;
				}
			}
			if (!found) {
				throw new IllegalArgumentException("Source data must load linked column " + fc.getName());
			}
		}
	}

	/**
	 * Link another file to this one with a {@link FileLink}. Key columns, specified by
	 * {@link FileLink#keys(FileColumn...)}, must share a name (though not necessarily the same header
	 * value) with a column expected to be loaded (either directly or by a filter) in this file.
	 * 
	 * @param link {@link FileLink}
	 * @return this
	 */
	public AbstractFileParserFactory link(FileLink link) {
		checkLink(link);
		this.parser.addLink(link);
		return this;
	}

	public FileParser build() {
		return parser;
	}

}
