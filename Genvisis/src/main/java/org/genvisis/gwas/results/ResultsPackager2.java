package org.genvisis.gwas.results;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public final class ResultsPackager2 {

	/**
	 * Input file reader
	 */
	BufferedReader reader;
	/**
	 * Output file writer
	 */
	PrintWriter writer;
	/**
	 * Output header column names
	 */
	String[] header;
	/**
	 * Input delimiter, decided internally with {@code ext#determineDelimiter(String)}
	 */
	String inputDelim;
	/**
	 * "Black box" formatter
	 */
	ResultFormatParser resultLineParser;
	/**
	 * Logger object
	 */
	private Logger log;
	/**
	 * Call {@code String#trim()} on input lines after reading?
	 */
	private boolean trim;
	/**
	 * Prefix of lines to skip
	 */
	private String skipPrefix = null;

	public ResultsPackager2(ResultFormatParser lineParser) {
		this(lineParser, new Logger());
	}

	public ResultsPackager2(ResultFormatParser lineParser, Logger log) {
		resultLineParser = lineParser;
		this.log = log;
	}

	/**
	 * Opens the input file reader, scans to the header, determines the delimiter, and sends the split
	 * header to the ResultLineParser. <br />
	 * <br />
	 * Gives options for skipping lines based on a prefix, a given count, or both. If the header line
	 * itself starts with the same prefix, it will be skipped - therefore it is recommended to use
	 * line count skipping instead.
	 * 
	 * @param input Input filename (with full path)
	 * @param lineSkipPrefix Prefix of lines to skip (e.g. '#')
	 * @param initialLinesToSkip Number of lines to skip (not including lines starting with the
	 *        specified prefix)
	 * @param log Logger object
	 * @throws IOException thrown if there is a problem when opening the file or reading to the header
	 */
	public void openInput(String input, boolean trimInputLines, String lineSkipPrefix,
												int initialLinesToSkip) throws IOException {
		reader = Files.getAppropriateReader(input);
		trim = trimInputLines;
		skipPrefix = lineSkipPrefix;
		String line = null;
		int skip = initialLinesToSkip;
		while ((line = reader.readLine()) != null) {
			if (trim) {
				line = line.trim();
			}
			if (skipPrefix != null && line.startsWith(skipPrefix)) {
				continue;
			}
			if (skip > 0) {
				skip--;
				continue;
			}
			inputDelim = ext.determineDelimiter(line);
			header = line.split(inputDelim);
			break;
		}
		resultLineParser.setInputHeader(header);
	}

	/**
	 * Open the output file writer and print the header line if desired.
	 * 
	 * @param output Output file name (with full path)
	 * @param append Append to file if true
	 * @param writeHeader Write the header line to the file if true (will likely be the inverse of
	 *        append, i.e. if append is true, writeHeader will likely be false)
	 */
	public void openOutput(String output, boolean append, boolean writeHeader) {
		writer = Files.getAppropriateWriter(output, append);
		if (writeHeader) {
			if (append) {
				log.reportTimeWarning("Appending to output file but also writing header line.");
			}
			writer.println(ArrayUtils.toStr(resultLineParser.getOutputHeader(),
																			resultLineParser.getOutputDelim()));
		}
	}

	/**
	 * Process the input file line-by-line using the given {@code ResultFormatParser}.
	 * 
	 * @throws IOException thrown if there's a problem when reading input lines.
	 */
	public void parse() throws IOException {
		String line = null;
		long t1 = System.nanoTime();
		while ((line = reader.readLine()) != null) {
			if (trim) {
				line = line.trim();
			}
			if (skipPrefix != null && line.startsWith(skipPrefix)) {
				continue;
			}
			String[] result = resultLineParser.parseInputLine(line.split(inputDelim));
			if (result != null) {
				writer.println(ArrayUtils.toStr(result, resultLineParser.getOutputDelim()));
			}
		}
		writer.flush();
		log.reportTimeElapsed("Parsed file in ", t1);
	}

	/**
	 * Close the internal reader and writer.
	 * 
	 * @throws IOException if there's a problem when closing the reader.
	 */
	public void close() throws IOException {
		reader.close();
		writer.close();
	}

	/**
	 * Static utility method to do everything (largely intended as an example of usage).
	 * 
	 * Input parameters match those required by instance methods.
	 */
	public static void runPackager(String input, boolean trimInputLines, String lineSkipPrefix,
																 int initialLinesToSkip, String output, boolean append,
																 boolean writeHeaderToOutput, ResultFormatParser lineParser,
																 Logger log) throws IOException {
		ResultsPackager2 packager = new ResultsPackager2(lineParser);
		packager.openInput(input, trimInputLines, lineSkipPrefix, initialLinesToSkip);
		packager.openOutput(output, append, writeHeaderToOutput);
		packager.parse();
		packager.close();
	}

}
