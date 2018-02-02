package org.genvisis.gwas.parsing;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.genvisis.common.Files;
import org.genvisis.common.ext;

public class FileParser implements Iterable<DataLine> {

	private String inputFile;
	private String inputFileDelim;

	private int skippedLines = 0;
	private Set<String> skipPrefices;

	private boolean failBlank = false;
	private boolean failCount = true;

	private List<ColumnFilter> filters;
	private Map<ColumnFilter, Boolean> filterDeath;

	private boolean noHeader = false;
	private boolean trimInput = true;

	private String parseFailValue = ".";

	private List<FileLink> linkedParsers;
	private Map<FileLink, Map<FileColumn<?>, FileColumn<?>>> linkedFileColumns;
	protected List<FileColumn<?>> dataInOrder;
	protected List<FileColumn<?>> addlDataToLoad;
	protected List<FileColumn<?>> optlDataToLoad;
	protected List<FileColumn<?>> optlDataFound;

	private BufferedReader reader;
	private boolean opened = false;
	private long lineCount = 0;
	private String[] header;

	/**
	 * Don't use - create FileParsers with {@link FileParserFactory#setup(String, FileColumn...)}
	 * instead.
	 * 
	 * @param inputFile String full path to file
	 * @param inputDelim Input delimiter, or null if the delimiter should be determined from the
	 *        header line
	 */
	FileParser(String inputFile, String inputDelim) {
		this.inputFile = inputFile;
		this.inputFileDelim = inputDelim;
		skipPrefices = new HashSet<>();
		filters = new ArrayList<>();
		filterDeath = new HashMap<>();
		linkedParsers = new ArrayList<>();
		dataInOrder = new ArrayList<>();
		addlDataToLoad = new ArrayList<>();
		optlDataToLoad = new ArrayList<>();
		linkedFileColumns = new HashMap<>();
	}

	/**
	 * Add a filter, specifying whether to die or not if a line fails to pass the filter.
	 * 
	 * @param filter {@link ColumnFilter}
	 * @param dieIfFail boolean
	 */
	void addFilter(ColumnFilter filter, boolean dieIfFail) {
		if (!this.filters.contains(filter)) {
			this.filters.add(filter);
			this.filterDeath.put(filter, dieIfFail);
			for (FileColumn<?> filterCol : filter.getFilterColumns()) {
				if (!dataInOrder.contains(filterCol)) {
					addlDataToLoad.add(filterCol);
				}
			}
		}
	}

	/**
	 * See {@link AbstractFileParserFactory#skipLines(int)}
	 */
	void setSkippedLines(int num) {
		this.skippedLines = num;
	}

	/**
	 * See {@link AbstractFileParserFactory#skipPrefix(String)}
	 */
	void addSkippedPrefix(String prefix) {
		this.skipPrefices.add(prefix);
	}

	/**
	 * See {@link AbstractFileParserFactory#parseFailValue(String)}
	 */
	void setParseFailValue(String pfv) {
		this.parseFailValue = pfv;
	}

	/**
	 * See {@link AbstractFileParserFactory#skipLinesWithIncorrectNumberOfColumns()}
	 */
	void skipLineWithIncorrectNumberOfColumns() {
		this.failCount = false;
	}

	/**
	 * See {@link AbstractFileParserFactory#failOnBlankLines()}
	 */
	void failOnBlankLines() {
		this.failBlank = true;
	}

	/**
	 * See {@link AbstractFileParserFactory#noHeader()}
	 */
	void setNoHeader() {
		this.noHeader = true;
	}

	/**
	 * See {@link AbstractFileParserFactory#noTrim()}
	 */
	void setNoTrim() {
		this.trimInput = false;
	}

	/**
	 * See {@link AbstractFileParserFactory#link(FileLink)}
	 * 
	 * @param link {@link FileLink}
	 */
	void addLink(FileLink link) {
		linkedParsers.add(link);
		linkedFileColumns.put(link, new HashMap<>());
		Set<FileColumn<?>> all = new HashSet<>();
		all.addAll(getDataColumnsInOrder());
		all.addAll(getAddlDataColumns());
		for (FileColumn<?> fc : link.getKeys()) {
			for (FileColumn<?> fc1 : all) {
				if (fc.getName().equals(fc1.getName())) {
					linkedFileColumns.get(link).put(fc, fc1);
					break;
				}
			}
		}
	}

	/**
	 * See {@link AbstractFileParserFactory#optionalColumns(FileColumn...)}
	 * 
	 * @param columns {link {@link FileColumn}...
	 */
	void addOptionalColumns(FileColumn<?>... columns) {
		for (FileColumn<?> fc : columns) {
			this.optlDataToLoad.add(fc);
		}
	}

	private List<FileColumn<?>> getOutputColumnsInOrder() {
		List<FileColumn<?>> cols = new ArrayList<>(dataInOrder);
		for (FileLink fl : linkedParsers) {
			for (FileColumn<?> fc : fl.getValues()) {
				if (!cols.contains(fc)) {
					cols.add(fc);
				}
			}
		}
		return cols;
	}

	private void open() throws IOException {
		if (opened) {
			return;
		}
		reader = Files.getAppropriateReader(inputFile);
		String line = null;
		int skip = skippedLines;
		header = null;
		lineCount = 0;
		Map<String, Integer> hdrMap = null;
		while ((line = reader.readLine()) != null) {
			lineCount++;
			if (trimInput) {
				line = line.trim();
			}
			if (skipPrefices.size() > 0) {
				for (String skp : skipPrefices) {
					if (line.startsWith(skp)) {
						continue;
					}
				}
			}
			if (skip > 0) {
				skip--;
				continue;
			}
			if (inputFileDelim == null) {
				inputFileDelim = ext.determineDelimiter(line);
			}
			if (!noHeader) {
				header = line.split(inputFileDelim);
				hdrMap = new HashMap<String, Integer>();
				for (int i = 0; i < header.length; i++) {
					hdrMap.put(header[i], i);
				}
			}
			for (FileColumn<?> fc : dataInOrder) {
				fc.initialize(hdrMap);
			}
			for (FileColumn<?> fc : addlDataToLoad) {
				fc.initialize(hdrMap);
			}
			optlDataFound = new ArrayList<>();
			for (FileColumn<?> fc : optlDataToLoad) {
				try {
					fc.initialize(hdrMap);
					optlDataFound.add(fc);
				} catch (IllegalStateException e) {
					// not in file
				}
			}
			break;
		}

		opened = true;
	}

	@SuppressWarnings("unchecked")
	private DataLine readLine() throws IOException {
		if (reader == null)
			return null;

		DataLine lineData;
		String[] parts;

		boolean skip = false;
		do {
			skip = false;
			parts = null;
			lineData = null;

			String line = reader.readLine();
			lineCount++;
			if (line == null)
				break;
			if (trimInput) {
				line = line.trim();
			}
			if (failBlank && line.equals("")) {
				throw new IllegalStateException("Line "
																				+ lineCount + " was blank!");
			}
			if (skipPrefices.size() > 0) {
				for (String skp : skipPrefices) {
					if (line.startsWith(skp)) {
						skip = true;
						break;
					}
				}
			}
			if (skip) {
				continue;
			}
			lineData = new DataLine(parseFailValue);
			parts = line.split(inputFileDelim, -1);
			if (failCount && parts.length != header.length) {
				throw new IllegalStateException("Line "
																				+ lineCount + " was the wrong length; expected "
																				+ header.length + ", found " + parts.length);
			}
			for (@SuppressWarnings("rawtypes")
			FileColumn fc : dataInOrder) {
				try {
					// unchecked
					lineData.put(fc, fc.getValue(parts));
				} catch (ParseFailureException e) {
					lineData.fail(fc);
				}
			}
			for (@SuppressWarnings("rawtypes")
			FileColumn fc : addlDataToLoad) {
				try {
					// unchecked
					lineData.put(fc, fc.getValue(parts));
				} catch (ParseFailureException e) {
					lineData.fail(fc);
				}
			}
			for (@SuppressWarnings("rawtypes")
			FileColumn fc : optlDataFound) {
				try {
					// unchecked
					lineData.put(fc, fc.getValue(parts));
				} catch (ParseFailureException e) {
					lineData.fail(fc);
				}
			}
			for (ColumnFilter fc : filters) {
				if (!fc.filter(lineData)) {
					if (filterDeath.get(fc)) {
						String colNames = "";
						List<FileColumn<?>> cols = fc.getFilterColumns();
						for (int f = 0; f < cols.size(); f++) {
							colNames += cols.get(f).getName();
							if (f < cols.size() - 1) {
								colNames += ", ";
							}
						}
						throw new IllegalStateException("Filter on columns " + colNames + " failed on line "
																						+ lineCount);
					}
					skip = true;
					break;
				}
			}
		} while (skip);

		if (lineData != null) {
			for (FileLink fl : linkedParsers) {
				String keyVal = "";
				for (FileColumn<?> fc : fl.getKeys()) {
					keyVal += lineData.getString(linkedFileColumns.get(fl).get(fc));
				}
				DataLine linkedLine = fl.get(keyVal);
				if (linkedLine == null) {
					// missing data
					for (FileColumn<?> valueCol : fl.getValues()) {
						// check if data already exists
						if (!lineData.has(valueCol)) {
							lineData.fail(valueCol);
						}
					}
				} else {
					for (@SuppressWarnings("rawtypes")
					FileColumn valueCol : fl.getValues()) {
						// check if data already exists or was a failure
						if (!lineData.hasValid(valueCol)) {
							// unchecked
							lineData.put(valueCol, linkedLine.getUnsafe(valueCol));
						} else if (lineData.hasValid(valueCol)
											 && linkedLine.hasValid(valueCol)
											 // unchecked
											 && !linkedLine.getUnsafe(valueCol).equals(lineData.getUnsafe(valueCol))) {
							// if a value already exists and isn't the same nor a parse failure
							throw new IllegalStateException("Different value found for column "
																							+ valueCol.getName() + " with linked key; key="
																							+ keyVal + ", value1=" + lineData.getString(valueCol)
																							+ ", value2=" + linkedLine.getString(valueCol));
						}
					}
				}
			}
		}

		return lineData;
	}

	/**
	 * @return An unmodifiable view of the data columns that will be present in the loaded data.
	 */
	public List<FileColumn<?>> getDataColumnsInOrder() {
		return Collections.unmodifiableList(dataInOrder);
	}

	/**
	 * @return An unmodifiable view of the additional columns of data that will be loaded (e.g. for
	 *         {@link ColumnFilter} operations.
	 */
	public List<FileColumn<?>> getAddlDataColumns() {
		return Collections.unmodifiableList(addlDataToLoad);
	}

	/**
	 * @return The input file path
	 */
	public String getInputFile() {
		return inputFile;
	}

	/**
	 * Close the internal {@link BufferedReader}.
	 * 
	 * @throws IOException
	 */
	public void close() throws IOException {
		reader.close();
	}

	public Iterator<DataLine> iterator() {
		try {
			open();
		} catch (IOException e1) {
			throw new RuntimeException(e1);
		}
		return new Iterator<DataLine>() {
			boolean started = false;
			DataLine currentLine = null;

			@Override
			public DataLine next() {
				DataLine line = null;
				try {
					if (!started) {
						// if we're calling next() without having called hasNext() first
						currentLine = readLine();
						started = true;
					}
					line = currentLine;
					// read the next line, which will return null when we're done
					currentLine = readLine();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
				return line;
			}

			@Override
			public boolean hasNext() {
				if (!started) {
					/*
					 * if we haven't started reading, we don't know if we have or will have any data (due to
					 * filters / skipped lines), so load the first line of data.
					 */
					try {
						currentLine = readLine();
					} catch (IOException e) {
						throw new RuntimeException(e);
					}
					started = true;
				}
				return currentLine != null;
			}
		};
	}

	public List<DataLine> load() throws IOException {
		List<DataLine> data = new ArrayList<>();
		Iterator<DataLine> iter = iterator();
		while (iter.hasNext()) {
			data.add(iter.next());
		}
		close();
		return data;
	}

	/**
	 * 
	 * @param dropIfKeyFail Set to true to drop any lines where the key value(s) throw a
	 *        ParseFailureException; set to false to throw a RuntimeException if invalid keys are
	 *        found.
	 * @param keyCols Vararg of FileColumns whose String values are to be concatenated together and
	 *        used as the key
	 * @return
	 * @throws IOException
	 */
	public Map<String, DataLine> load(boolean dropIfKeyFail,
																		FileColumn<?>... keyCols) throws IOException {
		Map<String, DataLine> data = new HashMap<>();
		Iterator<DataLine> iter = iterator();
		while (iter.hasNext()) {
			DataLine line = iter.next();
			try {
				String keyVal = "";
				for (FileColumn<?> fc : keyCols) {
					keyVal += line.get(fc);
				}
				data.put(keyVal, line);
			} catch (ParseFailureException e) {
				if (!dropIfKeyFail) {
					throw new RuntimeException(e);
				}
			}
		}
		close();
		return data;
	}

	public List<DataLine> parseToFileAndLoad(String outputFile, String outDelim) throws IOException {
		return parseToFileAndLoad(outputFile, outDelim, null);
	}

	public List<DataLine> parseToFileAndLoad(String outputFile, String outDelim,
																					 List<FileColumn<?>> outputOrder) throws IOException {
		List<DataLine> data = new ArrayList<>();
		Iterator<DataLine> iter = iterator();

		PrintWriter writer = Files.getAppropriateWriter(outputFile);
		StringBuilder lineOut = new StringBuilder();
		List<FileColumn<?>> outputColumns = outputOrder != null ? outputOrder
																														: getOutputColumnsInOrder();
		for (int i = 0, count = outputColumns.size(); i < count; i++) {
			lineOut.append(outputColumns.get(i).getName());
			if (i < count - 1) {
				lineOut.append(outDelim);
			}
		}
		writer.println(lineOut.toString());

		while (iter.hasNext()) {
			DataLine line = iter.next();
			data.add(line);
			lineOut = new StringBuilder();
			for (int i = 0, count = outputColumns.size(); i < count; i++) {
				lineOut.append(line.getString(outputColumns.get(i)));
				if (i < count - 1) {
					lineOut.append(outDelim);
				}
			}
			writer.println(lineOut.toString());
		}
		writer.close();
		close();
		return data;
	}

	public void parseToFile(String outputFile, String outDelim) throws IOException {
		this.parseToFile(outputFile, outDelim, null);
	}

	public void parseToFile(String outputFile, String outDelim,
													List<FileColumn<?>> outputOrder) throws IOException {
		this.parseToFile(outputFile, outDelim, true, false, outputOrder);
	}

	public void parseToFile(String outputFile, String outDelim, boolean writeHeader,
													boolean append) throws IOException {
		parseToFile(outputFile, outDelim, writeHeader, append, null);
	}

	public void parseToFile(String outputFile, String outDelim, boolean writeHeader,
													boolean append, List<FileColumn<?>> outputOrder) throws IOException {
		Iterator<DataLine> iter = iterator();

		PrintWriter writer = Files.getAppropriateWriter(outputFile, append);
		List<FileColumn<?>> outputColumns = outputOrder == null ? getOutputColumnsInOrder()
																														: outputOrder;
		if (outputOrder != null) {
			List<FileColumn<?>> allOutput = getOutputColumnsInOrder();
			if (!outputColumns.containsAll(allOutput)) {
				for (FileColumn<?> fc : allOutput) {
					if (!outputColumns.contains(fc)) {
						outputColumns.add(fc);
					}
				}
			}
		}

		StringBuilder lineOut = new StringBuilder();
		if (writeHeader) {
			for (int i = 0, count = outputColumns.size(); i < count; i++) {
				lineOut.append(outputColumns.get(i).getName());
				if (i < count - 1) {
					lineOut.append(outDelim);
				}
			}
			writer.println(lineOut.toString());
		}

		while (iter.hasNext()) {
			DataLine line = iter.next();
			lineOut = new StringBuilder();
			for (int i = 0, count = outputColumns.size(); i < count; i++) {
				lineOut.append(line.getString(outputColumns.get(i)));
				if (i < count - 1) {
					lineOut.append(outDelim);
				}
			}
			writer.println(lineOut.toString());
		}
		writer.close();
		close();
	}


}
