package org.genvisis.common;

import java.io.IOException;
import java.io.Reader;

/**
 * a parser that will extract typed data for each line based on indices provided
 *
 */
public class TypedFileParser extends FileParser {
	private int[][] numericColumns;
	private int[][] stringColumns;
	private final boolean setInvalidNumericToNaN;

	public TypedFileLine readTypedLine() throws IOException {
		String[] line = readLineArray();
		double[][] numericData = null;
		String[][] stringData = null;
		boolean validLine = true;
		int numInvalidNumerics = 0;
		// TODO, check length of line
		try {
			if (numericColumns != null) {
				numericData = new double[numericColumns.length][];
				for (int i = 0; i < numericColumns.length; i++) {
					String[] tmp = Array.subArray(line, numericColumns[i]);
					numericData[i] = new double[tmp.length];
					for (int j = 0; j < tmp.length; j++) {
						try {
							numericData[i][j] = Double.parseDouble(tmp[j]);
						} catch (NumberFormatException nfe) {
							if (!setInvalidNumericToNaN) {
								log.reportTimeError("Found invalid number on line " + Array.toStr(line));
								validLine = false;
							} else {
								numInvalidNumerics++;
								numericData[i][j] = Double.NaN;
							}
						}
					}
				}
			}
			if (stringColumns != null) {
				stringData = new String[stringColumns.length][];
				for (int i = 0; i < stringColumns.length; i++) {
					stringData[i] = Array.subArray(line, stringColumns[i]);
				}
			}
			return new TypedFileLine(numericData, stringData, validLine, numInvalidNumerics);
		} catch (ArrayIndexOutOfBoundsException aiob) {
			log.reportTimeError("Could not extract indices from line " + Array.toStr(line));
			if (numericColumns != null) {
				String numeric = "";
				for (int[] numericColumn : numericColumns) {
					numeric += "\t" + Array.toStr(numericColumn);
				}
				log.reportTimeError("Indices for numeric data: " + numeric);
			}
			if (stringColumns != null) {
				String string = "";
				for (int[] stringColumn : stringColumns) {
					string += "\t" + Array.toStr(stringColumn);
				}
				log.reportTimeError("Indices for numeric data: " + string);
			}
			// log.reportException(aiob);
			// aiob.printStackTrace();
			return new TypedFileLine(null, null, false, numInvalidNumerics);
		}
	}

	public int[][] getNumericColumns() {
		return numericColumns;
	}

	public void setNumericColumns(int[][] numericColumns) {
		this.numericColumns = numericColumns;
	}

	public void setStringColumns(int[][] stringColumns) {
		this.stringColumns = stringColumns;
	}

	public void addStringColumns(int[] col) {
		if (stringColumns == null) {
			stringColumns = new int[][] {col};
		} else {
			int[][] tmp = stringColumns;
			stringColumns = new int[stringColumns.length + 1][];
			stringColumns[0] = col;
			for (int i = 1; i < stringColumns.length; i++) {
				stringColumns[i] = tmp[i - 1];
			}
		}
	}

	public int[][] getStringColumns() {
		return stringColumns;
	}

	public static class Builder extends FileParser.Builder<Builder> {
		private int[][] numericDataIndices;
		private int[][] stringDataIndices;
		private boolean setInvalidNumericToNaN;

		public Builder numericDataIndices(int[][] numericDataIndices) {
			this.numericDataIndices = numericDataIndices;
			return this;
		}

		public Builder stringDataIndices(int[][] stringDataIndices) {
			this.stringDataIndices = stringDataIndices;
			return this;
		}

		public Builder setInvalidNumericToNaN(boolean setInvalidNumericToNaN) {
			this.setInvalidNumericToNaN = setInvalidNumericToNaN;
			return this;
		}

		@Override
		public TypedFileParser build(Reader in) {
			return new TypedFileParser(this, in);
		}
	}

	private TypedFileParser(Builder builder, Reader in) {
		super(builder, in);
		numericColumns = builder.numericDataIndices;
		stringColumns = builder.stringDataIndices;
		setInvalidNumericToNaN = builder.setInvalidNumericToNaN;
	}

	public static class TypedFileLine {
		private final double[][] numericData;
		private final String[][] stringData;
		private final boolean hasNumericData;
		private final boolean hasStringData;
		private final boolean validLine;
		private final int numInvalidNumerics;

		public TypedFileLine(	double[][] numericData, String[][] stringData, boolean validLine,
													int numInvalidNumerics) {
			super();
			this.numericData = numericData;
			this.stringData = stringData;
			hasNumericData = numericData != null;
			hasStringData = stringData != null;
			this.validLine = validLine;
			this.numInvalidNumerics = numInvalidNumerics;
		}

		public boolean hasNumericData() {
			return hasNumericData;
		}

		public int getNumInvalidNumerics() {
			return numInvalidNumerics;
		}

		public boolean hasStringData() {
			return hasStringData;
		}

		public double[][] getNumericData() {
			return numericData;
		}

		public String[][] getStringData() {
			return stringData;
		}

		public boolean isValidLine() {
			return validLine;
		}

	}
}
