package parse;

import java.io.IOException;
import java.io.Reader;

import common.Array;

/**
 * a parser that will extract typed data for each line based on indices provided
 *
 */
public class TypedFileParser extends FileParser {
	private int[][] numericColumns;
	private int[][] stringColumns;

	public TypedFileLine readTypedLine() throws IOException {
		String[] line = readLineArray();
		double[][] numericData = null;
		String[][] stringData = null;
		boolean validLine = true;
		// TODO, check length of line
		if (numericColumns != null) {
			numericData = new double[numericColumns.length][];
			try {
				for (int i = 0; i < numericColumns.length; i++) {
					numericData[i] = Array.toDoubleArray(Array.subArray(line, numericColumns[i]));
				}
			} catch (NumberFormatException e) {
				log.reportTimeError("Found invalid number on line " + Array.toStr(line));
				validLine = false;
			}
		}
		if (stringColumns != null) {
			stringData = new String[stringColumns.length][];
			for (int i = 0; i < stringColumns.length; i++) {
				stringData[i] = Array.subArray(line, stringColumns[i]);
			}
		}
		return new TypedFileLine(numericData, stringData, validLine);
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
			stringColumns = new int[][] { col };
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

		public Builder numericDataIndices(int[][] numericDataIndices) {
			this.numericDataIndices = numericDataIndices;
			return this;
		}

		public Builder stringDataIndices(int[][] stringDataIndices) {
			this.stringDataIndices = stringDataIndices;
			return this;
		}

		public TypedFileParser build(Reader in) {
			return new TypedFileParser(this, in);
		}
	}

	private TypedFileParser(Builder builder, Reader in) {
		super(builder, in);
		this.numericColumns = builder.numericDataIndices;
		this.stringColumns = builder.stringDataIndices;
	}

	public static class TypedFileLine {
		private double[][] numericData;
		private String[][] stringData;
		private boolean hasNumericData;
		private boolean hasStringData;
		private boolean validLine;

		public TypedFileLine(double[][] numericData, String[][] stringData, boolean validLine) {
			super();
			this.numericData = numericData;
			this.stringData = stringData;
			this.hasNumericData = numericData != null;
			this.hasStringData = stringData != null;
			this.validLine = validLine;
		}

		public boolean hasNumericData() {
			return hasNumericData;
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
