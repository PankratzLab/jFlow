/**
 * 
 */
package org.genvisis.common.matrix;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.ejml.alg.dense.decomposition.svd.SvdImplicitQrDecompose_D64;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
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
  private NamedRealMatrix loadings;
  private final int numComponents;

  /**
   * @param rowNameMap see {@link NamedRealMatrix}
   * @param columnNameMap see {@link NamedRealMatrix}
   * @param m see {@link NamedRealMatrix}
   * @param v the v {@link RealMatrix}, (PCs)
   * @param w the w {@link RealMatrix}, eigenvector
   */

  public SVD(Map<String, Integer> rowNameMap, Map<String, Integer> columnNameMap, DenseMatrix64F m,
             RealMatrix v, DiagonalMatrix w) {
    super(rowNameMap, columnNameMap, m);
    this.v = v;
    this.w = w;
    this.numComponents = w.getColumnDimension();

    if (w.getColumnDimension() != v.getRowDimension()) {
      throw new IllegalArgumentException("Diminsion mismatch for  V (num rows="
                                         + v.getRowDimension() + ") and W (num columns = "
                                         + w.getColumnDimension() + ")");
    }
  }

  /**
   * @return the {@link RealMatrix} v, in our case we usually call em "PCs"
   */
  public RealMatrix getPCs() {
    return v;
  }

  /**
   * @return the {@link NamedRealMatrix} loadings
   */
  public NamedRealMatrix getLoadings() {
    return loadings;
  }

  /**
   * @param outputRoot full path - ".pcs.gz" will be appended to the output root
   * @param columnOneTitle the first column will have this entry in the header (e.g "PCs", "SAMPLE")
   * @param log {@link Logger}
   */
  public void dumpPCsToText(String outputRoot, String columnOneTitle, Logger log) {
    new File(ext.parseDirectoryOfFile(outputRoot)).mkdirs();
    String out = outputRoot + ".pcs.gz";
    PrintWriter writer = Files.getAppropriateWriter(out);
    log.reportTimeInfo("Writing PCs to " + out);
    StringJoiner joiner = new StringJoiner("\t");
    joiner.add(columnOneTitle);
    List<String> namedComponents = getNamedComponents(numComponents);
    for (String namedComponent : namedComponents) {
      joiner.add(namedComponent);
    }

    for (int component = 0; component < numComponents; component++) {}
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

    for (int row = 0; row < m.getNumRows(); row++) {
      StringJoiner loadingString = new StringJoiner("\t");
      loadingString.add(getIndexRowMap().get(row));

      for (int component = 0; component < numComponents; component++) {
        loadingString.add(Double.toString(loadings.m.get(row, component)));
      }
      writer.println(loadingString.toString());
    }
    writer.close();
  }

  //TODO - could potentially read marker and dump loading line by line
  private void computeLoadings() {
    //    Will have all markers, but not all "PCs" all the time
    DenseMatrix64F loadingData = new DenseMatrix64F(getM().numRows, numComponents);

    Map<String, Integer> rowMap = new HashMap<>();
    for (int row = 0; row < m.numRows; row++) {
      rowMap.put(getIndexRowMap().get(row), row);
    }

    Map<String, Integer> loadingMap = new HashMap<>();
    for (int component = 0; component < numComponents; component++) {
      loadingMap.put("LOADING" + component, component);
    }
    for (int row = 0; row < m.numRows; row++) {

      DenseMatrix64F rowData = new DenseMatrix64F(1, m.numCols);
      CommonOps.extractRow(m, row, rowData);

      for (int component = 0; component < numComponents; component++) {
        double loading = getLoading(w.getEntry(component, component), rowData.data,
                                    v.getRow(component));
        loadingData.add(row, component, loading);
      }
    }
    this.loadings = new NamedRealMatrix(rowMap, loadingMap, loadingData);
  }

  // TODO, could refactor to static method that just takes singular values "w" and loadings "u"

  public NamedRealMatrix getExtraploatedPCs(NamedRealMatrix other, Logger log) {
    if (!getRowNameMap().keySet().equals(other.getRowNameMap().keySet())) {
      throw new IllegalArgumentException("All rows from data to be extrapolated must be present");
    }
    log.reportTimeInfo("Extrapolating PCs");

    Map<String, Integer> pcMap = new HashMap<>();
    List<String> namedComponents = getNamedComponents(numComponents);
    for (int i = 0; i < namedComponents.size(); i++) {
      pcMap.put(namedComponents.get(i), i);
    }

    DenseMatrix64F pcs = new DenseMatrix64F(other.getM().numCols, numComponents);
    for (int row = 0; row < other.getM().numRows; row++) {
      for (int component = 0; component < numComponents; component++) {
        for (int column = 0; column < other.getM().numCols; column++) {
          pcs.add(column, component,
                  other.getM().get(row, column) * loadings.getM().get(row, component));
        }
      }
    }
    for (int pc = 0; pc < pcs.numCols; pc++) {
      for (int row = 0; row < pcs.numRows; row++) {
        double current = pcs.get(row, pc);
        pcs.set(row, pc, current / w.getEntry(pc, pc));
      }
    }

    return new NamedRealMatrix(other.getColumnNameMap(), pcMap, pcs);
  }

  private static List<String> getNamedComponents(int numComponents) {
    List<String> namedComponents = new ArrayList<String>();
    for (int component = 0; component < numComponents; component++) {
      namedComponents.add("PC" + (component + 1));
    }
    return namedComponents;
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
    svd.decompose(m.getM());
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
    tv.reshape(numComponents, m.getM().numCols, true);
    DiagonalMatrix w = new DiagonalMatrix(singularValues);
    RealMatrix v = MatrixUtils.createRealMatrix(tv.numRows, tv.numCols);
    for (int row = 0; row < tv.numRows; row++) {
      for (int col = 0; col < tv.numCols; col++) {
        v.addToEntry(row, col, tv.get(row, col));
      }
    }

    SVD svdResult = new SVD(m.getRowNameMap(), m.getColumnNameMap(), m.getM(), v, w);
    svdResult.computeLoadings();
    return svdResult;
  }

  public void writeSerial(String serFile) {
    SerializedFiles.writeSerial(this, serFile, true);
  }

  /**
   * @param serFile
   * @return
   */
  public static SVD readSerial(String serFile) {
    return (SVD) SerializedFiles.readSerial(serFile);
  }
}
