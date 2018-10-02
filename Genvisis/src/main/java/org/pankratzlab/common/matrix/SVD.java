/**
 * 
 */
package org.pankratzlab.common.matrix;

import java.io.File;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.ejml.alg.dense.decomposition.svd.SvdImplicitQrDecompose_D64;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ejml.ops.SingularOps;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.SerializedFiles;
import org.pankratzlab.common.ext;

/**
 * Stores the U,V,W results from Singular value decompositions, note U is computed by us
 */
public class SVD implements Serializable {

  /**
   * 
   */
  private static final long serialVersionUID = 1L;
  private final NamedRealMatrix v;
  private final DiagonalMatrix w;
  private NamedRealMatrix loadings;
  private final int numComponents;

  /**
   * @param v the v {@link NamedRealMatrix}, (PCs)
   * @param w the w {@link DiagonalMatrix}, eigenvector
   */

  private SVD(NamedRealMatrix v, DiagonalMatrix w) {
    this.v = v;
    this.w = w;
    this.numComponents = w.getColumnDimension();

    if (w.getColumnDimension() != v.getDenseMatrix().getNumRows()) {
      throw new IllegalArgumentException("Diminsion mismatch for  V (num rows="
                                         + v.getDenseMatrix().getNumRows()
                                         + ") and W (num columns = " + w.getColumnDimension()
                                         + ")");
    }
  }

  /**
   * @return the {@link RealMatrix} v, in our case we usually call em "PCs"
   */
  public NamedRealMatrix getPCs() {
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
    writer.println(joiner.toString());

    for (int outputRow = 0; outputRow < v.getDenseMatrix().getNumCols(); outputRow++) {
      StringJoiner sample = new StringJoiner("\t");
      sample.add(v.getNameForColumnIndex(outputRow));
      for (int component = 0; component < numComponents; component++) {//
        sample.add(Double.toString(v.getDenseMatrix().get(component, outputRow)));
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

    for (int row = 0; row < loadings.getDenseMatrix().getNumRows(); row++) {
      StringJoiner loadingString = new StringJoiner("\t");
      loadingString.add(loadings.getNameForRowIndex(row));

      for (int component = 0; component < numComponents; component++) {
        loadingString.add(Double.toString(loadings.m.get(row, component)));
      }
      writer.println(loadingString.toString());
    }
    writer.close();
  }

  //TODO - could potentially read marker and dump loading line by line
  private void computeLoadings(NamedRealMatrix m) {
    //    Will have all markers, but not all "PCs" all the time
    DenseMatrix64F loadingData = new DenseMatrix64F(m.getDenseMatrix().numRows, numComponents);

    Map<String, Integer> rowMap = new HashMap<>();
    for (int row = 0; row < m.getDenseMatrix().numRows; row++) {
      rowMap.put(m.getNameForRowIndex(row), row);
    }

    Map<String, Integer> loadingMap = new HashMap<>();
    for (int component = 0; component < numComponents; component++) {
      loadingMap.put("LOADING" + component, component);
    }
    for (int row = 0; row < m.getDenseMatrix().numRows; row++) {

      DenseMatrix64F rowData = new DenseMatrix64F(1, m.getDenseMatrix().numCols);
      CommonOps.extractRow(m.getDenseMatrix(), row, rowData);

      for (int component = 0; component < numComponents; component++) {
        DenseMatrix64F componentData = new DenseMatrix64F(1, v.getDenseMatrix().numCols);

        double loading = getLoading(w.getEntry(component, component), rowData.data,
                                    CommonOps.extractRow(v.getDenseMatrix(), component,
                                                         componentData).data);
        loadingData.add(row, component, loading);
      }
    }
    this.loadings = new NamedRealMatrix(rowMap, loadingMap, loadingData);
  }

  // TODO, could refactor to static method that just takes singular values "w" and loadings "u"

  public NamedRealMatrix getExtraploatedPCs(NamedRealMatrix other, Logger log) {
    //we want to make sure the other has all data available, but we do not care if it has more
    if (!other.getRowMap().keySet().containsAll(loadings.getRowMap().keySet())) {
      throw new IllegalArgumentException("All rows from data to be extrapolated must be present");
    }
    log.reportTimeInfo("Extrapolating PCs");

    DenseMatrix64F pcs = new DenseMatrix64F(other.getDenseMatrix().numCols, numComponents);
    for (int row = 0; row < loadings.getDenseMatrix().numRows; row++) {
      int otherRow = other.getRowIndexFor(loadings.getNameForRowIndex(row));
      for (int component = 0; component < numComponents; component++) {
        for (int column = 0; column < other.getDenseMatrix().numCols; column++) {
          pcs.add(column, component, other.getDenseMatrix().get(otherRow, column)
                                     * loadings.getDenseMatrix().get(row, component));
        }
      }
    }
    for (int pc = 0; pc < pcs.numCols; pc++) {
      for (int row = 0; row < pcs.numRows; row++) {
        double current = pcs.get(row, pc);
        pcs.set(row, pc, current / w.getEntry(pc, pc));
      }
    }

    return new NamedRealMatrix(other.getColumnMap(), getNamedComponentsMap(numComponents), pcs);
  }

  private static Map<String, Integer> getNamedComponentsMap(int numComponents) {
    Map<String, Integer> pcMap = new HashMap<>();
    List<String> namedComponents = getNamedComponents(numComponents);
    for (int i = 0; i < namedComponents.size(); i++) {
      pcMap.put(namedComponents.get(i), i);
    }
    return pcMap;
  }

  private static List<String> getNamedComponents(int numComponents) {
    List<String> namedComponents = new ArrayList<>();
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
    svd.decompose(m.getDenseMatrix());
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
    tv.reshape(numComponents, m.getDenseMatrix().numCols, true);
    NamedRealMatrix vNamedRealMatrix = new NamedRealMatrix(getNamedComponentsMap(numComponents),
                                                           m.getColumnMap(), tv);
    DiagonalMatrix w = new DiagonalMatrix(singularValues);

    SVD svdResult = new SVD(vNamedRealMatrix, w);
    svdResult.computeLoadings(m);
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
