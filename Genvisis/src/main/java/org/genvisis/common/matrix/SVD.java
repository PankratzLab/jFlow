/**
 * 
 */
package org.genvisis.common.matrix;

import java.io.File;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.StringJoiner;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.ejml.alg.dense.decomposition.svd.SvdImplicitQrDecompose_D64;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.SingularOps;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;

/**
 * Stores the V,W and input data for Singular value decompositions
 */
public class SVD extends NamedRealMatrix {

  /**
   * 
   */
  private static final long serialVersionUID = 1L;
  private final RealMatrix v;
  private final DiagonalMatrix w;
  private final int numComponents;

  /**
   * @param rowNameMap see {@link NamedRealMatrix}
   * @param columnNameMap see {@link NamedRealMatrix}
   * @param m see {@link NamedRealMatrix}
   * @param v the v {@link RealMatrix}, (PCs)
   * @param w the w {@link RealMatrix}, eigenvector
   */

  public SVD(Map<String, Integer> rowNameMap, Map<String, Integer> columnNameMap, RealMatrix m,
             RealMatrix v, DiagonalMatrix w, int numComponents) {
    super(rowNameMap, columnNameMap, m);
    this.v = v;
    this.w = w;
    this.numComponents = numComponents;
    if (numComponents != w.getColumnDimension()) {
      throw new IllegalArgumentException("Invalid number of singular values ("
                                         + w.getColumnDimension() + ") for " + numComponents
                                         + " components");
    }

    if (numComponents != v.getRowDimension()) {
      throw new IllegalArgumentException("Invalid diminsion for V  (" + v.getRowDimension()
                                         + ") for " + numComponents + " components");
    }
  }

  /**
   * @param outputRoot full path - ".pcs" will be appended to the output root
   * @param columnOneTitle the first column will have this entry in the header (e.g "PCs", "SAMPLE")
   * @param log {@link Logger}
   */
  public void dumpPCsToText(String outputRoot, String columnOneTitle, Logger log) {
    new File(ext.parseDirectoryOfFile(outputRoot)).mkdirs();
    String out = outputRoot + ".pcs.gz";
    PrintWriter writer = Files.getAppropriateWriter(out);
    log.reportTimeInfo("Writing PCs to " + out);
    StringJoiner joiner = new StringJoiner("\t");
    joiner.add("SAMPLE");
    for (int component = 0; component < numComponents; component++) {
      joiner.add("PC" + (component + 1));
    }
    writer.println(joiner.toString());

    for (int outputRow = 0; outputRow < v.getColumnDimension(); outputRow++) {
      StringJoiner sample = new StringJoiner("\t");
      sample.add(getIndexColumnMap().get(outputRow));
      for (int component = 0; component < numComponents; component++) {//
        sample.add(Double.toString(v.getEntry(component, outputRow)));
      }
      writer.println(sample.toString());

    }
    writer.close();
  }

  /**
   * @param outputRoot full path - ".singularValues.gz" will be appended to the output root
   * @param columnOneTitle the first column will have this entry in the header (e.g
   *          "SingularValues")
   * @param log {@link Logger}
   */
  public void dumpSingularValuesToText(String outputRoot, String columnOneTitle, Logger log) {
    new File(ext.parseDirectoryOfFile(outputRoot)).mkdirs();
    String out = outputRoot + ".singularValues.gz";
    log.reportTimeInfo("Writing singular values to " + out);
    StringJoiner joiner = new StringJoiner("\n");
    joiner.add(columnOneTitle);

    for (int component = 0; component < numComponents; component++) {
      joiner.add(Double.toString(w.getEntry(component, component)));
    }
  }

  /**
   * @param outputRoot full path - ".loadings.gz" will be appended to the output root
   * @param columnOneTitle the first column will have this entry in the header (e.g "MARKER")
   * @param log {@link Logger}
   */
  public void dumpLoadingsToText(String outputRoot, String columnOneTitle, Logger log) {
    new File(ext.parseDirectoryOfFile(outputRoot)).mkdirs();
    String out = outputRoot + ".loadings.gz";
    PrintWriter writer = Files.getAppropriateWriter(out);
    log.reportTimeInfo("Writing loadings to " + out);
    StringJoiner joiner = new StringJoiner("\t");
    joiner.add(columnOneTitle);

    for (int component = 0; component < numComponents; component++) {
      //   v is transposed releative to m
      joiner.add("LOADING_" + component);
    }
    writer.println(joiner.toString());

    for (int row = 0; row < m.getRowDimension(); row++) {
      StringJoiner loadings = new StringJoiner("\t");
      loadings.add(getIndexRowMap().get(row));
      double[] rowData = m.getRow(row);

      for (int component = 0; component < numComponents; component++) {
        double loading = getLoading(w.getEntry(component, component), rowData, v.getRow(component));

        loadings.add(Double.toString(loading));
      }
      writer.println(loadings.toString());
    }
    writer.close();
  }

  private NamedRealMatrix computeLoadings() {
    //    Will have all markers, but not all "PCs" all the time
    RealMatrix loadings = MatrixUtils.createRealMatrix(getM().getRowDimension(), numComponents);

    Map<String, Integer> rowMap = new HashMap<>();
    for (int row = 0; row < m.getRowDimension(); row++) {
      rowMap.put(getIndexRowMap().get(row), row);
    }

    Map<String, Integer> loadingMap = new HashMap<>();
    for (int component = 0; component < numComponents; component++) {
      loadingMap.put("LOADING" + component, component);
    }
    for (int row = 0; row < m.getRowDimension(); row++) {

      double[] rowData = m.getRow(row);
      for (int component = 0; component < numComponents; component++) {
        double loading = getLoading(w.getEntry(component, component), rowData, v.getRow(component));
        loadings.addToEntry(row, component, loading);
      }
    }
    return new NamedRealMatrix(rowMap, loadingMap, loadings);
  }

  public NamedRealMatrix getExtraploatedPCs(NamedRealMatrix other, Logger log) {
    if (!other.getRowNameMap().keySet().containsAll(getRowNameMap().keySet())
        || getRowNameMap().keySet().size() != other.getRowNameMap().keySet().size()) {
      throw new IllegalArgumentException("All rows from data to be extrapolated must be present");
    }
    log.reportTimeInfo("Extrapolating PCs");
    NamedRealMatrix loadings = computeLoadings();

    Map<String, Integer> pcMap = new HashMap<>();
    for (int component = 0; component < numComponents; component++) {
      pcMap.put("PC" + (component + 1), component);
    }
    RealMatrix pcs = MatrixUtils.createRealMatrix(getM().getColumnDimension(), numComponents);
    for (int row = 0; row < other.getM().getRowDimension(); row++) {
      for (int component = 0; component < numComponents; component++) {
        for (int column = 0; column < other.getM().getColumnDimension(); column++) {
          pcs.addToEntry(column, component, other.getM().getEntry(row, column)
                                            * loadings.getM().getEntry(row, component));
        }
      }
    }
    for (int pc = 0; pc < pcs.getColumnDimension(); pc++) {
      for (int row = 0; row < pcs.getRowDimension(); row++) {
        double current = pcs.getEntry(row, pc);
        pcs.setEntry(row, pc, current / w.getEntry(pc, pc));
      }
    }

    return new NamedRealMatrix(other.getColumnNameMap(), pcMap, pcs);
  }

  private static double getLoading(double singularValue, double[] data, double[] basis) {
    double loading = 0;
    for (int i = 0; i < basis.length; i++) {
      loading += data[i] * basis[i];
    }
    return loading / singularValue;
  }

  /**
   * @param m use {@link SvdImplicitQrDecompose_D64 } to perform SVD on this {@link NamedRealMatrix}
   * @param numComponents number of components that will be stored
   * @param log
   * @return
   */
  public static SVD computeSVD(NamedRealMatrix m, int numComponents, Logger log) {
    log.reportTimeInfo("Computing EJML PCs");
    SvdImplicitQrDecompose_D64 svd = new SvdImplicitQrDecompose_D64(false, false, true, false);
    svd.decompose(MatrixOperations.toDenseMatrix64F(m.getM()));
    log.memoryPercentTotalFree();

    log.reportTimeInfo("Finished Computing EJML PCs");
    DenseMatrix64F tv = svd.getV(null, true);

    DenseMatrix64F tmpW = svd.getW(null);
    SingularOps.descendingOrder(null, false, tmpW, tv, true);
    int numSingular = Math.min(numComponents, Math.min(tmpW.numRows, tmpW.numCols));
    double[] singularValues = new double[numSingular];
    for (int i = 0; i < numSingular; i++) {
      singularValues[i] = tmpW.get(i, i);
    }
    tv.reshape(numComponents, m.getM().getColumnDimension(), true);
    DiagonalMatrix w = new DiagonalMatrix(singularValues);
    RealMatrix v = MatrixUtils.createRealMatrix(tv.numRows, tv.numCols);
    for (int row = 0; row < tv.numRows; row++) {
      for (int col = 0; col < tv.numCols; col++) {
        v.addToEntry(row, col, tv.get(row, col));
      }
    }
    return new SVD(m.getRowNameMap(), m.getColumnNameMap(), m.getM(), v, w, numComponents);
  }

  public void writeSerial(String serFile) {
    SerializedFiles.writeSerial(this, serFile, true);
  }

  public static SVD readSerial(String serFile) {
    return (SVD) SerializedFiles.readSerial(serFile);
  }
}
