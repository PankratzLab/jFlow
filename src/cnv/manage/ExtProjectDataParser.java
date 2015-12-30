package cnv.manage;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.Hashtable;

import common.Array;
import common.Files;
import common.TypedFileParser;
import common.ext;
import common.TypedFileParser.TypedFileLine;
import cnv.filesys.Project;

/**
 * Trying to have a single framework to load and match sample or marker data to a project from external files This has not been heavily tested, so please verify and test before use...
 * 
 *
 */
public class ExtProjectDataParser {
	private static final String DEFUALT_SAMPLE_NAME_COLUMN = "DNA";
	private static final int DEFUALT_SAMPLE_NAME_COLUMN_INDEX = 0;
	private static final String DEFUALT_MISSING_STRING = ".";

	// private static final String DEFUALT_MARKER_NAME_COLUMN = "marker";
	// private static final int DEFUALT_MARKER_NAME_COLUMN_INDEX = 0;
	private Project proj;
	private TypedFileParser typedFileParser;
	private double[][] numericData;
	private String[][] stringData;
	private String[] numericDataTitles, stringDataTitles;
	private int dataKeyColumnIndex;
	private String dataKeyColumnName;
	private String[] dataToLoad;
	private boolean[] dataPresent;
	private boolean requireAll, verbose;
	private boolean skipUnFound;
	private boolean hasHeader;
	private boolean treatAllNumeric;
	private String fullPathToDataFile;
	private String missingString;
	private String commentString;
	private String[] headerFlags;
	private boolean sampleBased;
	private boolean concatMultipleStringEntries;
	private boolean setInvalidNumericToNaN;
	private boolean firstEntryOnly;

	public void loadData() {
		Hashtable<String, Integer> markerIndices = proj.getMarkerIndices();
		proj.getLog().reportTimeInfo("Attempting to load " + (sampleBased ? "sample based" : "marker based") + " data from " + fullPathToDataFile);
		if (treatAllNumeric) {
			treatAllNumeric();
		}
		init();
		try {
			if (hasHeader) {
				scanToData();
			}
			int numInvalidNumerics = 0;
			while (typedFileParser.ready()) {
				TypedFileLine typedFileLine = typedFileParser.readTypedLine();
				if (typedFileLine.isValidLine()) {
					String data = typedFileLine.getStringData()[0][0];
					int dataIndex = -1;
					if (sampleBased) {
						dataIndex = ext.indexOfStr(data, dataToLoad);
					} else {
						if (markerIndices.containsKey(data)) {
							dataIndex = markerIndices.get(data);
						}
					}
					if (dataIndex < 0) {
						if (skipUnFound) {
							if (verbose) {
								proj.getLog().reportTimeWarning("Did not find " + data + " in project's " + (sampleBased ? "samples" : "markers") + ", skipping");
							}
						} else {
							proj.getLog().reportTimeError("Did not find " + data + " in project's " + (sampleBased ? "samples" : "markers") + ", halting");
							typedFileParser.close();
							return;
						}
					} else {
						boolean concat = dataPresent[dataIndex];
						if (!concat || !firstEntryOnly) {
							if (concat && (!concatMultipleStringEntries || typedFileLine.hasNumericData())) {
								if (typedFileLine.hasNumericData()) {
									proj.getLog().reportTimeError("Multiple entries for the same " + (sampleBased ? "sample" : "marker") + " " + data + " were found and numeric indices were requested, cannot load...");

								} else {
									proj.getLog().reportTimeError("Multiple entries for the same " + (sampleBased ? "sample" : "marker") + " " + data + " were found, consider setting the concatMultipleStringEntries option if only string data is needed");
								}
								return;
							}
							dataPresent[dataIndex] = true;
							if (typedFileLine.hasNumericData() && numericData != null) {
								for (int i = 0; i < numericData.length; i++) {
									if (typedFileLine.getNumInvalidNumerics() > 0 && !setInvalidNumericToNaN) {
										proj.getLog().reportTimeError("failed to load " + numericDataTitles[i]);
										proj.getLog().reportTimeInfo("Set the setInvalidNumericToNaN flag to avoid this shutdown");
										return;
									}
									numInvalidNumerics += typedFileLine.getNumInvalidNumerics();
									numericData[i][dataIndex] = typedFileLine.getNumericData()[0][i];
								}
							}
							if (typedFileLine.hasStringData() && stringData != null) {
								for (int i = 0; i < stringData.length; i++) {
									stringData[i][dataIndex] = (concat ? stringData[i][dataIndex] + typedFileParser.getSeparator() + typedFileLine.getStringData()[1][i] : typedFileLine.getStringData()[1][i]);// 0 is reserved for the samples
								}
							}
						}
					}
				}
			}
			typedFileParser.close();
			int dataLoaded = Array.booleanArraySum(dataPresent);
			if (requireAll && dataLoaded != dataToLoad.length) {
				proj.getLog().reportTimeError("All data were required to be present, but did not see all of them");
				stringData = null;
				numericData = null;
			} else {
				proj.getLog().reportTimeInfo("Loaded data for " + dataLoaded + " " + (sampleBased ? "sample(s)" : "marker(s) ") + " of " + dataPresent.length);
				if (numInvalidNumerics > 0) {
					proj.getLog().reportTimeWarning(numInvalidNumerics + "numeric entries had invalid doubles and were set to NaN");
				}
			}
		} catch (IOException e) {
			proj.getLog().reportTimeError("Error reading file " + fullPathToDataFile);
			e.printStackTrace();
		}

	}

	/**
	 * Reads until the data begins
	 * 
	 * @throws IOException
	 */
	private void scanToData() throws IOException {
		if (commentString != null) {
			proj.getLog().reportTimeWarning("The has header option was flagged and comments were set, skipping header and comments");
			String line;
			do {
				line = typedFileParser.readLine();
			} while (line != null && line.startsWith(commentString));
		} else if (headerFlags != null) {
			String line;
			proj.getLog().reportTimeWarning("The has header option was flagged and and header flags were set, skipping lines preceding header");
			do {
				line = typedFileParser.readLine();
				if (verbose) {
					proj.getLog().reportTimeInfo("Skipping line " + line);
				}
			} while (line != null && Array.countIf(ext.indexFactors(headerFlags, line.trim().split(typedFileParser.getSeparator()), true, proj.getLog(), false, false), -1) > 0);
		} else {

			proj.getLog().reportTimeWarning("The has header option was flagged and comments were not set, skipping the first line");
			typedFileParser.readLine();
		}
	}

	public String[] getStringDataAt(int index, boolean subset) {
		return subset ? Array.subArray(stringData[index], dataPresent) : stringData[index];
	}

	public boolean determineIndicesFromTitles() {
		boolean determined = true;
		if (hasHeader) {
			String[] header;
			if (headerFlags == null) {
				header = Files.getHeaderOfFile(fullPathToDataFile, typedFileParser.getSeparator(), new String[] { commentString }, proj.getLog());
			} else {
				header = Files.getLineContaining(fullPathToDataFile, typedFileParser.getSeparator(), headerFlags, proj.getLog());
			}
			if (header != null) {
				if (numericDataTitles == null && stringDataTitles == null) {
					proj.getLog().reportTimeError("Must provide valid arrays to determine indices");
					determined = false;
				}
				if (dataKeyColumnName == null) {
					proj.getLog().reportTimeError("Must provide a data name column to determine indices");
					determined = false;
				}
				determined = determineSampleIndex(determined, header);
				if (determined && stringDataTitles != null) {
					if (typedFileParser.getStringColumns() != null) {
						proj.getLog().reportTimeWarning("String columns were already provided, skipping string column assignment");
					} else {
						int[] tmp = ext.indexFactors(stringDataTitles, header, true, false);
						determined = Array.min(tmp) >= 0;
						if (!determined) {
							proj.getLog().reportTimeError("Could not find all string indices in header " + Array.toStr(header));
						} else {
							typedFileParser.setStringColumns(new int[][] { tmp });
						}
					}
				}
				if (determined && numericDataTitles != null) {
					if (typedFileParser.getNumericColumns() != null) {
						proj.getLog().reportTimeWarning("Numeric columns were already provided, skipping numeric column assignment");
					} else {
						int[] tmp = ext.indexFactors(numericDataTitles, header, true, false);
						determined = Array.min(tmp) >= 0;
						if (!determined) {
							proj.getLog().reportTimeError("Could not find all numeric indices in header " + Array.toStr(header));
						} else {
							typedFileParser.setNumericColumns(new int[][] { tmp });
						}
					}

				}
			} else {
				determined = false;
				proj.getLog().reportTimeError("Could not obtain header with current options");
			}
		} else {
			determined = false;
			proj.getLog().reportTimeError("Header option must be flagged to use this method");
		}

		return determined;
	}

	public double[][] getNumericData() {
		return numericData;
	}

	public String[][] getStringData() {
		return stringData;
	}

	public boolean hasStringDataForTitle(String title) {
		int index = stringDataTitles == null ? -1 : ext.indexOfStr(title, stringDataTitles);
		return index >= 0 && stringData != null;
	}

	public String[] getStringDataForTitle(String title) {
		int index = stringDataTitles == null ? -1 : ext.indexOfStr(title, stringDataTitles);
		String[] data = new String[] {};
		if (index < 0 || stringData == null) {
			proj.getLog().reportTimeError("Data for " + title + " was not found");
			proj.getLog().reportTimeError(stringDataTitles == null ? "No String data Titles available" : Array.toStr(stringDataTitles));
			if (stringData == null) {
				proj.getLog().reportTimeError("No string data available");
			}
		} else {
			data = stringData[index];
		}
		return data;
	}

	public double[] getNumericDataForTitle(String title) {
		int index = numericData == null ? -1 : ext.indexOfStr(title, numericDataTitles);
		double[] data = new double[] {};
		if (index < 0 || numericData == null) {
			proj.getLog().reportTimeError("Data for " + title + " was not found");
			proj.getLog().reportTimeError(numericData == null ? "No Numeric data Titles available" : "Titles available: " + Array.toStr(numericDataTitles));
			if (numericData == null) {
				proj.getLog().reportTimeError("No Numeric data available");
			}
		} else {
			data = numericData[index];
		}
		return data;
	}

	public String[] getNumericDataTitles() {
		return numericDataTitles;
	}

	public String[] getStringDataTitles() {
		return stringDataTitles;
	}

	public boolean[] getDataPresent() {
		return dataPresent;
	}

	public String[] getDataToLoad() {
		return dataToLoad;
	}

	private boolean determineSampleIndex(boolean determined, String[] header) {
		if (determined && !dataKeyColumnName.equals(DEFUALT_SAMPLE_NAME_COLUMN)) {
			dataKeyColumnIndex = ext.indexOfStr(dataKeyColumnName, header);
			if (dataKeyColumnIndex < 0) {
				proj.getLog().reportTimeError("Header must contain sample name column " + dataKeyColumnName);
				determined = false;
				proj.getLog().reportError("Header found was" + Array.toStr(header));
			}
		}
		return determined;
	}

	private void init() {
		String[] header = Files.getHeaderOfFile(fullPathToDataFile, typedFileParser.getSeparator(), proj.getLog());
		this.dataToLoad = sampleBased ? proj.getSamples() : proj.getMarkerNames();
		if (sampleBased) {
			proj.getLog().reportTimeInfo("Sample based loader:");
		} else {
			proj.getLog().reportTimeInfo("Marker based loader:");

		}
		this.dataPresent = new boolean[dataToLoad.length];
		Arrays.fill(dataPresent, false);

		if (typedFileParser.getNumericColumns() != null && typedFileParser.getNumericColumns().length > 1) {
			proj.getLog().reportTimeError("Currently this method only handles one dimensionional numeric extraction, even though its an int[][]");
			return;
		}
		if (typedFileParser.getStringColumns() != null && typedFileParser.getStringColumns().length > 1) {
			proj.getLog().reportTimeError("Currently this method only handles one dimensionional string extraction, even though its an int[][]");
			return;
		}

		typedFileParser.addStringColumns(new int[] { dataKeyColumnIndex });// samples will be in first array position

		if (typedFileParser.getNumericColumns() != null) {
			this.numericData = new double[typedFileParser.getNumericColumns()[0].length][dataToLoad.length];
			if (numericDataTitles == null) {
				proj.getLog().reportTimeInfo("Numeric titles were not provided, assigning...");
				if (hasHeader) {
					numericDataTitles = Array.subArray(header, typedFileParser.getNumericColumns()[0]);
				} else {
					numericDataTitles = Array.toStringArray(typedFileParser.getNumericColumns()[0]);
				}
			}
			for (int i = 0; i < numericData.length; i++) {
				Arrays.fill(numericData[i], Double.NaN);
			}
		}
		if (typedFileParser.getStringColumns() != null && typedFileParser.getStringColumns().length > 1) {// more than just the sample column we added above
			this.stringData = new String[typedFileParser.getStringColumns()[1].length][dataToLoad.length];
			if (stringDataTitles == null) {
				proj.getLog().reportTimeInfo("String titles were not provided, assigning...");
				if (hasHeader) {
					stringDataTitles = Array.subArray(header, typedFileParser.getStringColumns()[1]);
				} else {
					stringDataTitles = Array.toStringArray(typedFileParser.getStringColumns()[1]);
				}
			}

			for (int i = 0; i < stringData.length; i++) {
				Arrays.fill(stringData[i], missingString);
			}
		}
	}

	private void treatAllNumeric() {
		String[] header = Files.getHeaderOfFile(fullPathToDataFile, proj.getLog());
		if (!hasHeader || determineSampleIndex(true, header)) {
			proj.getLog().reportTimeInfo("Treating all columns besides index " + dataKeyColumnIndex + " (" + dataKeyColumnName + ") as numeric column(s)");
			proj.getLog().reportTimeInfo("If data titles are present, they are assumed to be in order of the column header");

			if (!hasHeader && numericDataTitles == null) {
				proj.getLog().reportTimeInfo("Since no header or data titles are present, numeric titles will be the relative indices");
			} else if (numericDataTitles != null && numericDataTitles.length != header.length - 1) {
				proj.getLog().reportTimeError("The length of the data title (" + numericDataTitles.length + ") does not match the data in the file's length (" + (header.length - 1) + ")");
				return;
			}
			int[] tmpNumeric = new int[header.length - 1];
			String[] tmpString = new String[header.length - 1];
			int index = 0;
			for (int i = 0; i < header.length; i++) {
				if (i != dataKeyColumnIndex) {
					tmpNumeric[index] = i;
					if (hasHeader && numericDataTitles == null) {
						tmpString[index] = header[i];
					} else if (numericDataTitles == null) {
						tmpString[index] = "Data" + index;
					} else {
						tmpString[index] = numericDataTitles[index];
					}
					index++;
				}
			}
			numericDataTitles = tmpString;
			typedFileParser.setNumericColumns(new int[][] { tmpNumeric });
			proj.getLog().reportTimeInfo("Using " + numericDataTitles.length + " numeric data columns");
			if (numericDataTitles.length < 10) {
				proj.getLog().reportTimeInfo("Data columns are: " + Array.toStr(numericDataTitles));

			}
		} else {
			proj.getLog().reportTimeError("Could not determine proper sample column index for treating all data as numeric");
		}
	}

	/**
	 * Builder for {@link ExtProjectDataParser}: the default options should handle a file with a header containing DNA, and all other columns (with header) containing numeric data
	 *
	 */
	public static class ProjectDataParserBuilder {
		TypedFileParser.Builder typeBuilder = new TypedFileParser.Builder();
		private String[] numericDataTitles;
		private String[] stringDataTitles;
		private int dataKeyColumnIndex = DEFUALT_SAMPLE_NAME_COLUMN_INDEX;
		private String dataKeyColumnName = DEFUALT_SAMPLE_NAME_COLUMN;
		private String[] dataToLoad;
		private boolean[] dataPresent;
		private boolean requireAll = false;
		private boolean skipUnFound = true;
		private boolean hasHeader = true;
		private boolean treatAllNumeric = true;
		private String missingString = DEFUALT_MISSING_STRING;
		private boolean sampleBased = true;
		private String commentString = null;
		private boolean verbose = true;
		private String[] headerFlags = null;
		private boolean concatMultipleStringEntries = false;
		private boolean setInvalidNumericToNaN = false;
		private boolean firstEntryOnly = false;

		/**
		 * @param separator
		 *            data separator for the file being loaded
		 * @return
		 */
		public ProjectDataParserBuilder separator(String separator) {
			typeBuilder.separator(separator);
			return this;
		}

		/**
		 * 
		 * @param numericDataIndices
		 *            indices to extract for numeric data...all data goes to double[]
		 * @return
		 */
		public ProjectDataParserBuilder numericDataIndices(int[][] numericDataIndices) {
			typeBuilder.numericDataIndices(numericDataIndices);
			return this;
		}

		/**
		 * @param stringDataIndices
		 *            indices to extract for string data...
		 * @return
		 */
		public ProjectDataParserBuilder stringDataIndices(int[][] stringDataIndices) {
			typeBuilder.stringDataIndices(stringDataIndices);
			return this;
		}

		/**
		 * @param numericDataTitles
		 *            optional titles of the numeric data...usually if a header is not present
		 * @return
		 */
		public ProjectDataParserBuilder numericDataTitles(String[] numericDataTitles) {
			this.numericDataTitles = numericDataTitles;
			return this;
		}

		/**
		 * @param numericDataTitles
		 *            optional titles of the string data...usually if a header is not present
		 * @return
		 */
		public ProjectDataParserBuilder stringDataTitles(String[] stringDataTitles) {
			this.stringDataTitles = stringDataTitles;
			return this;
		}

		/**
		 * @param numericDataTitles
		 *            optional titles of the string data...usually if a header is not present
		 * @return
		 */
		public ProjectDataParserBuilder setInvalidNumericToNaN(boolean setInvalidNumericToNaN) {
			this.setInvalidNumericToNaN = setInvalidNumericToNaN;
			typeBuilder.setInvalidNumericToNaN(setInvalidNumericToNaN);
			return this;
		}

		/**
		 * @param dataKeyColumnIndex
		 *            this column will be matched to either the projects samples or markers
		 * @return
		 */
		public ProjectDataParserBuilder dataKeyColumnIndex(int dataKeyColumnIndex) {
			this.dataKeyColumnIndex = dataKeyColumnIndex;
			return this;
		}

		/**
		 * @param firstEntryOnly
		 *            if true, only the first entry of a duplicate will be loaded
		 * @return
		 */
		public ProjectDataParserBuilder firstEntryOnly(boolean firstEntryOnly) {
			this.firstEntryOnly = firstEntryOnly;
			return this;
		}

		/**
		 * @param dataKeyColumnName
		 *            name of the data to match column
		 * @return
		 */
		public ProjectDataParserBuilder dataKeyColumnName(String dataKeyColumnName) {
			this.dataKeyColumnName = dataKeyColumnName;
			return this;
		}

		/**
		 * @param requireAll
		 *            require all samples or markers to be present in the file
		 * @return
		 */
		public ProjectDataParserBuilder requireAll(boolean requireAll) {
			this.requireAll = requireAll;
			return this;
		}

		/**
		 * @param skipUnFound
		 *            if a sample or marker is not found, keep going and just skip it
		 * @return
		 */
		public ProjectDataParserBuilder skipUnFound(boolean skipUnFound) {
			this.skipUnFound = skipUnFound;
			return this;
		}

		/**
		 * @param hasHeader
		 *            the file has a header
		 * @return
		 */
		public ProjectDataParserBuilder hasHeader(boolean hasHeader) {
			this.hasHeader = hasHeader;
			return this;
		}

		/**
		 * @param treatAllNumeric
		 *            all columns will be treated as numeric except for the column to be matched
		 * @return
		 */
		public ProjectDataParserBuilder treatAllNumeric(boolean treatAllNumeric) {
			this.treatAllNumeric = treatAllNumeric;
			return this;
		}

		/**
		 * @param missingString
		 *            for string data, set strings to this if a sample/marker is not found in the data file
		 * @return
		 */
		public ProjectDataParserBuilder missingString(String missingString) {
			this.missingString = missingString;
			return this;
		}

		/**
		 * @param sampleBased
		 *            the file will be matched to samples
		 * @return
		 */
		public ProjectDataParserBuilder sampleBased(boolean sampleBased) {
			this.sampleBased = sampleBased;
			return this;
		}

		/**
		 * @param commentString
		 *            lines starting with this string will be skipped when parsing the header
		 * @return
		 */
		public ProjectDataParserBuilder commentString(String commentString) {
			this.commentString = commentString;
			return this;
		}

		/**
		 * @param headerFlags
		 *            the header contains these lines, and others will be skipped
		 * @return
		 */
		public ProjectDataParserBuilder headerFlags(String[] headerFlags) {
			this.headerFlags = headerFlags;
			return this;
		}

		/**
		 * @param verbose
		 *            verbose reporting
		 * @return
		 */
		public ProjectDataParserBuilder verbose(boolean verbose) {
			this.verbose = verbose;
			return this;
		}

		/**
		 * @param concatMultipleStringEntries
		 *            if multiple entries are found for either a marker or sample, the results will be concatenated
		 * @return
		 */
		public ProjectDataParserBuilder concatMultipleStringEntries(boolean concatMultipleStringEntries) {
			this.concatMultipleStringEntries = concatMultipleStringEntries;
			return this;
		}

		/**
		 * Construct the parser with the options set by the builder
		 * 
		 * @param proj
		 *            the project that will match the data
		 * @param fullPathToDataFile
		 *            data file to load
		 * @return
		 * @throws FileNotFoundException
		 */
		public ExtProjectDataParser build(Project proj, String fullPathToDataFile) throws FileNotFoundException {
			return new ExtProjectDataParser(this, proj, typeBuilder, fullPathToDataFile);
		}
	}

	private ExtProjectDataParser(ProjectDataParserBuilder builder, Project proj, TypedFileParser.Builder typeBuilder, String fullPathToDataFile) throws FileNotFoundException {
		this.proj = proj;
		this.numericDataTitles = builder.numericDataTitles;
		this.stringDataTitles = builder.stringDataTitles;
		this.dataKeyColumnIndex = builder.dataKeyColumnIndex;
		this.dataKeyColumnName = builder.dataKeyColumnName;
		this.dataToLoad = builder.dataToLoad;
		this.dataPresent = builder.dataPresent;
		this.requireAll = builder.requireAll;
		this.skipUnFound = builder.skipUnFound;
		this.hasHeader = builder.hasHeader;
		this.treatAllNumeric = builder.treatAllNumeric;
		this.fullPathToDataFile = fullPathToDataFile;
		this.missingString = builder.missingString;
		this.sampleBased = builder.sampleBased;
		this.commentString = builder.commentString;
		this.verbose = builder.verbose;
		this.headerFlags = builder.headerFlags;
		this.concatMultipleStringEntries = builder.concatMultipleStringEntries;
		this.typedFileParser = typeBuilder.build(Files.getAppropriateReader(fullPathToDataFile));
		this.setInvalidNumericToNaN = builder.setInvalidNumericToNaN;
		this.firstEntryOnly = builder.firstEntryOnly;
	}
}
