package org.pankratzlab.shared.stats;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.ext;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Table;
import com.google.common.collect.TreeBasedTable;

public class CTable {

  private static final Set<String> DEFAULT_NULL_VALUES = ImmutableSet.of("", ".", "NA", "na", "N/A",
                                                                         "n/a", "NaN", "Missing",
                                                                         "missing");
  public static final String MISSING = "Missing";

  private Table<String, String, Integer> table;
  private Map<String, String> rowLabels;
  private Map<String, String> colLabels;
  private Map<String, Integer> rowCounts;
  private Map<String, Integer> colCounts;
  private boolean showRowProportions;

  /**
   * Construct a new contingency table from a set of parallel arrays. Each index {@code i} from [0,
   * data.length] adds 1 to the count of the row/column combination {@code rowData[i]} and
   * {@code colData[i]}.
   *
   * @param rowData Array of row data
   * @param colData Array of column data
   * @param axis Logic to use in adding "missing value" axis. Default: {@link AddMissing#NATURAL}
   */
  public CTable(String[][] rowLabels, String[][] colLabels, String[] rowData, String[] colData) {
    // These are parallel arrays and thus must be the same size
    if (rowData.length != colData.length || rowLabels.length != colLabels.length) {
      System.err.println("Error - arrays are not of the same length");
    }

    // NB: label arrays can be empty, which means "use actual values"
    // need to add "missing" category somewhere? if size == 1 or empty array

    this.rowLabels = buildMap(rowLabels);
    this.colLabels = buildMap(colLabels);

    // Parse the data points to build the table
    table = TreeBasedTable.create();
    for (int i = 0; i < rowData.length; i++) {
      String rowKey = sanitizeData(rowData[i], this.rowLabels);
      String colKey = sanitizeData(colData[i], this.colLabels);
      int count = 1 + safeGet(rowKey, colKey);
      table.put(rowKey, colKey, count);
    }

    validateSize(this.rowLabels);
    validateSize(this.colLabels);

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

  /**
   * Ensure the map has at least at least 2 entries. If it only has one entry, we can attempt to
   * auto-correct this by adding a {@link #MISSING} identity mapping
   *
   * @param labelMap
   * @throws IllegalStateException if we can not get the map to 2 valid entries
   */
  private void validateSize(Map<String, String> labelMap) {
    if (labelMap.size() == 0 || (labelMap.size() == 1 && labelMap.containsKey(MISSING))) {
      throw new IllegalStateException("No values found for contingency table");
    }
    if (labelMap.size() == 1) {
      labelMap.put(MISSING, MISSING);
    }
  }

  public String[] getFinalRowLabels() {
    return table.rowKeySet().toArray(new String[table.rowKeySet().size()]);
  }

  public String[] getFinalColunLabels() {
    return table.columnKeySet().toArray(new String[table.columnKeySet().size()]);
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

  public static CTable extrapolateCounts(String[][] rowLabels, String[][] colLabels, String[] rows,
                                         int[] genotype) {
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

    return new CTable(rowLabels, colLabels, extrapolate[0], extrapolate[1]);
  }

  public int[][] getContingencyTable() {
    // Convert the table to an int[][] of counts
    int[][] result = new int[rowLabels.keySet().size()][colLabels.keySet().size()];
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
   * @param expectedValues set of acceptable key values
   * @return The input string if not {@code null} and not in {@link #DEFAULT_NULL_VALUES}. Otherwise
   *         {@link #MISSING}
   */
  private String sanitizeData(String key, Map<String, String> expectedValues) {
    if (key == null || DEFAULT_NULL_VALUES.contains(key)) {
      key = MISSING;
    }
    if (!expectedValues.containsKey(key)) {
      expectedValues.put(key, key);
    }
    return key;
  }

  /**
   * Helper method to convert a 2-D string array to a map.
   * 
   * @param labels each sub-array is an expected value with an optional display label
   */
  private Map<String, String> buildMap(String[][] labels) {
    if (labels == null) {
      throw new IllegalArgumentException("Error: labels cannot be null.");
    }
    Map<String, String> map = new LinkedHashMap<>();

    for (int i = 0; i < labels.length; i++) {
      String[] pair = labels[i];

      if (pair == null || pair.length == 0 || pair.length > 2) {
        throw new IllegalArgumentException("Error: label mappings must be in the form [observed value, (optional) display value].");
      }
      if (pair[0].equals(MISSING)) {
        map.put(MISSING, "missing value");
      } else if (pair.length == 1) {
        map.put(pair[0], pair[0]);
      } else {
        map.put(pair[0], pair[1]);
      }

    }
    return map;
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
