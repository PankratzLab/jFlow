/**
 * 
 */
package org.genvisis.common.matrix;

import java.io.Serializable;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import org.apache.commons.math3.linear.RealMatrix;

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
  private final RealMatrix m;

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
  protected Map<Integer, String> getIndexRowMap() {
    return indexRowMap;
  }

  /**
   * @return the indexColumnMap
   */
  protected Map<Integer, String> getIndexColumnMap() {
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

}
