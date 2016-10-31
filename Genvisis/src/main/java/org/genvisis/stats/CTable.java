package org.genvisis.stats;

import java.util.Collections;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.HashVec;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;

public class CTable {
	// TODO merge with list in ext
	public static String[] DEFAULT_NULL_VALUES = {"", ".", "NA", "na", "N/A", "n/a", "NaN", "Missing",
																								"missing"};

	private Hashtable<String, Integer> cells;
	private String[][] rowLabels, columnLabels;
	private String[] nullValues;
	private boolean showRowProportions;

	public CTable() {}

	public CTable(int[] array1, int[] array2) {
		this(Array.toStringArray(array1), Array.toStringArray(array2));
	}

	public CTable(String[][] matrix) {
		this(matrix[0], matrix[1]);
	}


	// This is a finished copy for the labels format, but haven't been used.
	public CTable(String[] array1, String[] array2) {
		String index;
		Vector<String> rowLabelVector, columnLabelVector;

		if (array1.length != array2.length) {
			System.err.println("Error - arrays are not of the same length");
		}

		// Create the Hashtable cells, ArrayList rowLabels, and ArrayList columnLabels
		cells = new Hashtable<String, Integer>();
		rowLabelVector = new Vector<String>();
		columnLabelVector = new Vector<String>();
		for (int i = 0; i < array1.length; i++) {
			index = array1[i] + " x " + array2[i];
			if (cells.containsKey(index)) {
				cells.put(index, cells.get(index) + 1);
			} else {
				cells.put(index, 1);
			}

			HashVec.addIfAbsent(array1[i], rowLabelVector);
			HashVec.addIfAbsent(array2[i], columnLabelVector);
		}

		Collections.sort(rowLabelVector);
		Collections.sort(columnLabelVector);
		rowLabels = new String[rowLabelVector.size()][2];
		columnLabels = new String[columnLabelVector.size()][2];
		// for(int i=0; i<array1.length; i++) {
		for (int i = 0; i < rowLabelVector.size(); i++) {
			rowLabels[i][0] = rowLabels[i][1] = rowLabelVector.get(i);
		}
		for (int i = 0; i < columnLabelVector.size(); i++) {
			columnLabels[i][0] = columnLabels[i][1] = columnLabelVector.get(i);
		}
		// Need to add missing value processing here

		nullValues = DEFAULT_NULL_VALUES;
		showRowProportions = false;
	}

	public void setCustomNullValues(String[] newNullValues) {
		nullValues = newNullValues;
	}

	public String[] getFinalRowLabels() {
		return Matrix.extractColumn(rowLabels, 1);
	}

	public String[] getFinalColunLabels() {
		return Matrix.extractColumn(columnLabels, 1);
	}


	public void setCustomLabelsAndOrder(String[][] customRowLabelsAndOrder,
																			String[][] customColumnLabelsAndOrder) {
		if (customRowLabelsAndOrder != null) {
			rowLabels = customRowLabelsAndOrder;
		}

		if (customColumnLabelsAndOrder != null) {
			columnLabels = customColumnLabelsAndOrder;
		}
	}

	public String getCTableInHtml() {
		String output;
		int[][] contingencyTable;
		double[][] rowProportions;

		output = "<html><table border=\"1\">";// <tr><td></td>" + columnLabels.get(1) + "</tr>";

		contingencyTable = getContingencyTable();
		rowProportions = Matrix.computeRowProportions(contingencyTable);
		for (int i = -1; i < rowLabels.length; i++) {
			for (int j = -1; j < columnLabels.length; j++) {
				if (i == -1) {
					if (j == -1) {
						output = output + "<tr><td></td>";
					} else {
						output = output + "<td>" + columnLabels[j][1] + "</td>";
					}
				} else if (j == -1) {
					output = output + "<tr><td align=\"center\">" + rowLabels[i][1] + "</td>";
				} else {
					output = output	+ "<td align=\"center\">" + contingencyTable[i][j]
										+ (showRowProportions	? " (" + ext.formDeci(rowProportions[i][j], 3, true) + ")"
																					: "")
										+ "</td>";
				}
			}
			output = output + "</tr>";
		}
		output = output + "</table>" + "</html>";
		return output;
	}

	public static String[][] extrapolateCounts(String[] rows, int[] genotype) {
		int maxValue;
		String[][] extrapolate;

		maxValue = Array.max(genotype);
		extrapolate = new String[2][genotype.length * maxValue];
		for (int i = 0; i < genotype.length; i++) {
			for (int j = 0; j < maxValue; j++) {
				extrapolate[0][i * maxValue + j] = rows[i];
				if (genotype[i] < 0) {
					extrapolate[1][i * maxValue + j] = ".";
				} else if (genotype[i] > j) {
					extrapolate[1][i * maxValue + j] = "A";
				} else {
					extrapolate[1][i * maxValue + j] = "B";
				}
			}
		}

		return extrapolate;
	}

	public int[][] getContingencyTable() {
		String[] possibleRowValues, possibleColumnValues;
		int[][] result;

		result = new int[rowLabels.length][columnLabels.length];
		for (int i = 0; i < rowLabels.length; i++) {
			for (int j = 0; j < columnLabels.length; j++) {
				result[i][j] = 0;
				if (rowLabels[i][0] == null) {
					possibleRowValues = nullValues;
				} else {
					possibleRowValues = new String[] {rowLabels[i][0]};
				}
				if (columnLabels[j][0] == null) {
					possibleColumnValues = nullValues;
				} else {
					possibleColumnValues = new String[] {columnLabels[j][0]};
				}
				for (String possibleRowValue : possibleRowValues) {
					for (String possibleColumnValue : possibleColumnValues) {
						if (cells.containsKey(possibleRowValue + " x " + possibleColumnValue)) {
							result[i][j] += cells.get(possibleRowValue + " x " + possibleColumnValue);
						}
					}
				}
			}
		}
		return result;
	}

}
