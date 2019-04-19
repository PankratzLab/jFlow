/**
 * 
 */
package org.genvisis.cnv.ejml.matrix;

import java.io.PrintWriter;
import java.io.Serializable;
import java.util.Map;
import java.util.StringJoiner;

import org.ejml.data.DenseMatrix64F;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;

import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;

/**
 * Attaches column and row names<-> index mapping to a {@link DenseMatrix64F} using {@link BiMap}s
 */
public class NamedRealMatrix implements Serializable {

  /**
   * 
   */
  private static final long serialVersionUID = 1L;

  private final BiMap<String, Integer> rowMap;
  private final BiMap<String, Integer> columnMap;
  protected final DenseMatrix64F m;

  /**
   * @param rowNameMap maps row names to row indices
   * @param columnNameMap maps column names to column indices
   * @param m {@link DenseMatrix64F} underlying data
   */
  public NamedRealMatrix(Map<String, Integer> rowNameMap, Map<String, Integer> columnNameMap,
                         DenseMatrix64F m) {
    super();
    this.rowMap = generateBiMap(rowNameMap, m.getNumRows());
    this.columnMap = generateBiMap(columnNameMap, m.getNumCols());
    this.m = m;
  }

  /**
   * Generates and validates the dimension map
   * 
   * @param map these key value pairs will form the {@link BiMap}
   * @param expectedSize
   * @return {@link BiMap}
   */
  private static BiMap<String, Integer> generateBiMap(Map<String, Integer> map, int expectedSize) {
    if (map == null || map.size() != expectedSize) {
      throw new IllegalArgumentException("Mismatched mapped names length, map size="
                                         + (map == null ? 0 : map.size()) + " should be "
                                         + expectedSize);
    }
    BiMap<String, Integer> biMap = HashBiMap.create(map.keySet().size());
    for (Map.Entry<String, Integer> entry : map.entrySet()) {
      if (entry.getValue() >= expectedSize) {
        throw new IllegalArgumentException("Mapped values greater than data dimension");

      }
      biMap.put(entry.getKey(), entry.getValue());
    }
    return biMap;
  }

  /**
   * @param rowName name to get the index of
   * @return the index
   */
  public int getRowIndexFor(String rowName) {
    if (!rowMap.containsKey(rowName)) {
      throw new IllegalArgumentException("invalid row name " + rowName);
    }
    return rowMap.get(rowName);
  }

  /**
   * @param columName name to get the index of
   * @return the index
   */
  public int getColumnIndexFor(String columName) {
    if (!columnMap.containsKey(columName)) {
      throw new IllegalArgumentException("invalid column name " + columName);
    }
    return columnMap.get(columName);
  }

  /**
   * @param row to get the name for
   * @return the name
   */
  public String getNameForRowIndex(int row) {
    if (!rowMap.inverse().containsKey(row)) {
      throw new IllegalArgumentException("invalid row index " + row);
    }
    return rowMap.inverse().get(row);
  }

  /**
   * @param column to get the name for
   * @return the name
   */
  public String getNameForColumnIndex(int column) {
    if (!columnMap.inverse().containsKey(column)) {
      throw new IllegalArgumentException("invalid column index " + column);
    }
    return columnMap.inverse().get(column);
  }

  /**
   * @return the {@link DenseMatrix64F} m
   */
  public DenseMatrix64F getDenseMatrix() {
    return m;
  }

  /**
   * @return the rowMap
   */
  public BiMap<String, Integer> getRowMap() {
    return rowMap;
  }

  /**
   * @return the columnMap
   */
  public BiMap<String, Integer> getColumnMap() {
    return columnMap;
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

    for (int column = 0; column < m.numCols; column++) {
      joiner.add(getNameForColumnIndex(column));
    }
    writer.println(joiner.toString());

    for (int row = 0; row < m.numRows; row++) {
      StringJoiner rows = new StringJoiner("\t");
      rows.add(getNameForRowIndex(row));

      for (int column = 0; column < m.numCols; column++) {
        rows.add(Double.toString(m.get(row, column)));
      }
      writer.println(rows.toString());
    }
    writer.close();
  }

}
