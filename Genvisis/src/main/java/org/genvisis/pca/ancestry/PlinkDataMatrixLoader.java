/**
 * 
 */
package org.genvisis.pca.ancestry;

import java.util.HashMap;
import java.util.Map;
import org.apache.commons.math3.linear.RealMatrix;
import org.ejml.data.DenseMatrix64F;
import org.genvisis.common.Logger;
import org.genvisis.common.matrix.MatrixDataLoading;
import org.genvisis.common.matrix.NamedRealMatrix;
import org.genvisis.filesys.DosageData;

/**
 * Implements {@link MatrixDataLoading} to load plink files to a {@link RealMatrix}, using -1,0,1,2
 */
public class PlinkDataMatrixLoader implements MatrixDataLoading {

  private final String dir;
  private final String plinkRoot;
  private final Logger log;

  public PlinkDataMatrixLoader(String dir, String plinkRoot, Logger log) {

    this.dir = dir;
    this.plinkRoot = plinkRoot;
    this.log = log;

  }

  /*
   * (non-Javadoc)
   * @see org.genvisis.cnv.analysis.pca.ancestry.AncestryDataLoading#loadDataToMatrix()
   */
  @Override
  public NamedRealMatrix getData() {
    DosageData d = DosageData.loadPlinkBinary(dir, null, null, plinkRoot, null, true, true);
    String[] samples = new String[d.getIds().length];
    for (int i = 0; i < samples.length; i++) {
      samples[i] = d.getIds()[i][0] + "\t" + d.getIds()[i][1];
    }

    String[] markers = d.getMarkerSet().getMarkerNames();
    Map<String, Integer> markerMap = new HashMap<>();
    for (int i = 0; i < markers.length; i++) {
      markerMap.put(markers[i], i);
    }

    Map<String, Integer> sampleMap = new HashMap<>();
    for (int i = 0; i < samples.length; i++) {
      sampleMap.put(samples[i], i);
    }
    log.reportTimeInfo("Preparing matrix for " + markers.length + " markers and " + samples.length
                       + " samples");
    DenseMatrix64F m = new DenseMatrix64F(markers.length, samples.length);
    for (int column = 0; column < m.numCols; column++) {
      for (int row = 0; row < m.numRows; row++) {
        m.set(row, column, d.getDosageValues()[row][column]);
      }
    }
    log.memoryPercentTotalFree();
    return new NamedRealMatrix(markerMap, sampleMap, m);
  }

}
