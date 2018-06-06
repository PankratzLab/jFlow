/**
 * 
 */
package org.genvisis.common.matrix;

import java.io.File;
import java.io.PrintWriter;
import java.util.Map;
import java.util.StringJoiner;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

/**
 * 
 *
 */
public class SVD extends NamedRealMatrix {

  /**
   * 
   */
  private static final long serialVersionUID = 1L;
  private final RealMatrix v;
  private final DiagonalMatrix w;

  /**
   * @param rowNameMap see {@link NamedRealMatrix}
   * @param columnNameMap see {@link NamedRealMatrix}
   * @param m see {@link NamedRealMatrix}
   * @param v the v {@link RealMatrix}, (PCs)
   * @param w the w {@link RealMatrix}, eigenvector
   */

  public SVD(Map<String, Integer> rowNameMap, Map<String, Integer> columnNameMap, RealMatrix m,
             RealMatrix v, DiagonalMatrix w) {
    super(rowNameMap, columnNameMap, m);
    this.v = v;
    this.w = w;
  }

  /**
   * @param outputRoot full path - ".pcs" will be appended to the output root
   */
  public void dumpPCsToText(String outputRoot, Logger log) {
    new File(ext.parseDirectoryOfFile(outputRoot)).mkdirs();
    String out = outputRoot + ".pcs";
    PrintWriter writer = Files.getAppropriateWriter(out);
    log.reportTimeInfo("Writing PCs to " + out);
    StringJoiner joiner = new StringJoiner("\t");
    joiner.add("SAMPLE");
    for (int i = 0; i < v.getColumnDimension(); i++) {
      joiner.add("PC" + (i + 1));
    }
    writer.println(joiner.toString());

    for (int i = 0; i < v.getRowDimension(); i++) {
      StringJoiner sample = new StringJoiner("\t");
      sample.add(getIndexColumnMap().get(i));
      for (int j = 0; j < v.getColumnDimension(); j++) {
        sample.add(Double.toString(v.getEntry(j, i)));
      }
      writer.println(sample.toString());

    }
    writer.close();
  }

}
