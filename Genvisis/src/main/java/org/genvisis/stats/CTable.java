package org.genvisis.stats;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.ext;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Table;
import com.google.common.collect.TreeBasedTable;

public class CTable {
	private static final Set<String> DEFAULT_NULL_VALUES = ImmutableSet.of("", ".", "NA", "na", "N/A",
																																				 "n/a",
																																				 "NaN", "Missing",
																																				 "missing");
	private static final String MISSING = "Missing";

	private Table<String, String, Integer> table;
	private Map<String, String> rowLabels;
	private Map<String, String> colLabels;
	private Map<String, Integer> rowCounts;
	private Map<String, Integer> colCounts;
	private boolean showRowProportions;

	public CTable() {}

	public CTable(int[] array1, int[] array2) {
		this(ArrayUtils.toStringArray(array1), ArrayUtils.toStringArray(array2));
	}

	public CTable(String[][] matrix) {
		this(matrix[0], matrix[1]);
	}

	/**
	 * Construct a new contingency table from a set of parallel arrays. Each index {@code i} from [0,
	 * data.length] adds 1 to the count of the row/column combination {@code rowData[i]} and
	 * {@code colData[i]}.
	 *
	 * @param rowData Array of row data
	 * @param colData Array of column data
	 */
	public CTable(String[] rowData, String[] colData) {
		// These are parallel arrays and thus must be the same size
		if (rowData.length != colData.length) {
			System.err.println("Error - arrays are not of the same length");
		}

		// Parse the data points to build the table
		table = TreeBasedTable.create();
		for (int i = 0; i < rowData.length; i++) {
			String rowKey = sanitizeNull(rowData[i]);
			String colKey = sanitizeNull(colData[i]);
			int count = 1 + safeGet(rowKey, colKey);
			table.put(rowKey, colKey, count);
		}

		rowLabels = buildIdentityMap(table.rowKeySet());
		colLabels = buildIdentityMap(table.columnKeySet());

		// Get total counts for rows/columns
		rowCounts = new HashMap<>();
		for (String rowKey : table.rowKeySet()) {
			rowCounts.put(rowKey, sum(table.row(rowKey)));
		}

		colCounts = new HashMap<>();
		for (String colKey : table.columnKeySet()) {
			colCounts.put(colKey, sum(table.column(colKey)));
		}

		showRowProportions = false;
	}

	public String[] getFinalRowLabels() {
		return table.rowKeySet().toArray(new String[table.rowKeySet().size()]);
	}

	public String[] getFinalColunLabels() {
		return table.columnKeySet().toArray(new String[table.columnKeySet().size()]);
	}


	/**
	 * Each element in these arrays be a length 2 array. Index 0 is the original row or column label,
	 * Index 1 is the new label to use.
	 * 
	 * @param customRowLabelsAndOrder Array mapping row labels
	 * @param customColumnLabelsAndOrder Array mapping column labels
	 */
	public void setCustomLabelsAndOrder(String[][] customRowLabelsAndOrder,
																			String[][] customColumnLabelsAndOrder) {
		if (customRowLabelsAndOrder != null && customRowLabelsAndOrder.length > 0) {
			rowLabels = buildMap(customRowLabelsAndOrder);
		}

		if (customColumnLabelsAndOrder != null && customColumnLabelsAndOrder.length > 0) {
			colLabels = buildMap(customColumnLabelsAndOrder);
		}
	}

	public String getCTableInHtml() {
		StringBuilder output = new StringBuilder("<html><table border=\"1\"><tr><td></td>");

		for (String col : colLabels.keySet()) {
			output.append("<td>");
			output.append(colLabels.get(col));
			output.append("</td>");
		}
		output.append("</tr>");

		for (String row : rowLabels.keySet()) {
			output.append("<tr><td align=\"center\">");
			output.append(rowLabels.get(row));
			output.append("</td>");
			for (String col : colLabels.keySet()) {
				output.append("<td align=\"center\">");
				int count = safeGet(row, col);
				output.append(count);
				if (showRowProportions) {
					double proportion = count / rowCounts.get(row).doubleValue();
					output.append(" (" + ext.formDeci(proportion, 3, true) + ")");
				}
				output.append("</td>");
			}
			output.append("</tr>");
		}
		output.append("</table></html>");

		return output.toString();
	}

	public static String[][] extrapolateCounts(String[] rows, int[] genotype) {
		int maxValue;
		String[][] extrapolate;

		maxValue = ArrayUtils.max(genotype);
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
		// Convert the table to an int[][] of counts
		int[][] result = new int[table.rowKeySet().size()][table.columnKeySet().size()];
		int i = 0;
		for (String row : rowLabels.keySet()) {
			int j = 0;
			for (String col : colLabels.keySet()) {
				result[i][j++] = safeGet(row, col);
			}
			i++;
		}
		return result;
	}

	/**
	 * Helper method to return the value at a given row/column combination. Handles missing values
	 * appropriately.
	 */
	private int safeGet(String row, String col) {
		if (table.contains(row, col)) {
			return table.get(row, col);
		}

		// Missing value
		return 0;
	}

	/**
	 * Helper method to unify all {@link #DEFAULT_NULL_VALUES}
	 *
	 * @return The input string if not {@code null} and not in {@link #DEFAULT_NULL_VALUES}. Otherwise
	 *         {@link #MISSING}
	 */
	private String sanitizeNull(String key) {
		if (key == null || DEFAULT_NULL_VALUES.contains(key)) {
			return MISSING;
		}
		return key;
	}

	/**
	 * Helper method to convert a 2-D string array to a map.
	 */
	private Map<String, String> buildMap(String[][] labels) {
		Map<String, String> map = new LinkedHashMap<>();

		for (int i = 0; i < labels.length; i++) {
			String[] pair = labels[i];
			if (pair.length != 2) {
				throw new IllegalArgumentException("When setting custom row/column labels, must provide arrays of [original label, new label]");
			}
			map.put(pair[0], pair[1]);
		}
		return map;
	}

	/**
	 * Helper method to build a map of strings to themselves
	 */
	private Map<String, String> buildIdentityMap(Set<String> values) {
		// Use a LinkedHashMap to provide consistent ordering.
		Map<String, String> identity = new LinkedHashMap<>();
		for (String v : values) {
			identity.put(v, v);
		}
		return identity;
	}

	/**
	 * Helper method to sum the counts in a row or column
	 */
	private Integer sum(Map<String, Integer> rowOrColumn) {
		int c = 0;
		for (Integer i : rowOrColumn.values()) {
			c += i;
		}
		return c;
	}

}
