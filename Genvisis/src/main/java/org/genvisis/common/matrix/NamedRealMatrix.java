/**
 * 
 */
package org.genvisis.common.matrix;

import java.io.PrintWriter;
import java.io.Serializable;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.StringJoiner;
import org.apache.commons.math3.linear.RealMatrix;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;

/**
 * Attaches column and row names to a {@link RealMatrix}
 */
public class NamedRealMatrix implements Serializable {

  /**
   * 
   */
  private static final long serialVersionUID = 1L;
  private final Map<String, Integer> rowNameMap;
  private final Map<String, Integer> columnNameMap;
  private final Map<Integer, String> indexRowMap;
  private final Map<Integer, String> indexColumnMap;
  protected final RealMatrix m;

  /**
   * @param rowNameMap maps row names to row indices
   * @param columnNameMap maps column names to column indices
   * @param m {@link RealMatrix} underlying data
   */
  public NamedRealMatrix(Map<String, Integer> rowNameMap, Map<String, Integer> columnNameMap,
                         RealMatrix m) {
    super();
    validateMap(rowNameMap, m.getRowDimension());
    validateMap(columnNameMap, m.getColumnDimension());
    this.rowNameMap = rowNameMap;
    this.columnNameMap = columnNameMap;
    this.indexRowMap = generateIndexMap(rowNameMap);
    this.indexColumnMap = generateIndexMap(columnNameMap);
    this.m = m;
  }

  private static Map<Integer, String> generateIndexMap(Map<String, Integer> map) {
    Map<Integer, String> iMap = new HashMap<>();
    for (String key : map.keySet()) {
      iMap.put(map.get(key), key);
    }
    return iMap;

  }

  void validateMap(Map<String, Integer> map, int numEntries) {
    if (map == null || map.size() != numEntries) {
      throw new IllegalArgumentException("Mismatched mapped names length, map size=" + map.size()
                                         + " should be " + numEntries);
    }
    HashSet<Integer> indices = new HashSet<>();
    HashSet<String> keys = new HashSet<>();

    for (String key : map.keySet()) {
      if (map.get(key) >= numEntries) {
        throw new IllegalArgumentException("Invalid named index");
      } else {
        indices.add(map.get(key));
        keys.add(key);
      }
    }
    if (indices.size() != map.size()) {
      throw new IllegalArgumentException("Duplicate indices provided");
    }
    if (keys.size() != map.size()) {
      throw new IllegalArgumentException("Duplicate keys provided");
    }
  }

  /**
   * @return the indexRowMap
   */
  public Map<Integer, String> getIndexRowMap() {
    return indexRowMap;
  }

  /**
   * @return the indexColumnMap
   */
  public Map<Integer, String> getIndexColumnMap() {
    return indexColumnMap;
  }

  /**
   * @return the rowNameMap
   */
  public Map<String, Integer> getRowNameMap() {
    return rowNameMap;
  }

  /**
   * @return the columnNameMap
   */
  public Map<String, Integer> getColumnNameMap() {
    return columnNameMap;
  }

  /**
   * @return the m
   */
  public RealMatrix getM() {
    return m;
  }

  /**
   * @param out the text file to write to
   * @param columnOneTitle the first column will have this entry in the header (e.g "PCs", "SAMPLE")
   * @param log
   */
  public void dumpToText(String out, String columnOneTitle, Logger log) {
    PrintWriter writer = Files.getAppropriateWriter(out);
    log.reportTimeInfo("Writing matrix to " + out);
    StringJoiner joiner = new StringJoiner("\t");
    joiner.add(columnOneTitle);

    for (int column = 0; column < m.getColumnDimension(); column++) {
      joiner.add(getIndexColumnMap().get(column));
    }
    writer.println(joiner.toString());

    for (int row = 0; row < m.getRowDimension(); row++) {
      StringJoiner rows = new StringJoiner("\t");
      rows.add(getIndexRowMap().get(row));

      for (int column = 0; column < m.getColumnDimension(); column++) {
        rows.add(Double.toString(m.getEntry(row, column)));
      }
      writer.println(rows.toString());
    }
    writer.close();
  }

}
