/**
 * 
 */
package org.genvisis.pca.ancestry;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.genvisis.common.Logger;
import org.genvisis.common.matrix.MatrixDataLoading;
import org.genvisis.common.matrix.NamedRealMatrix;
import org.genvisis.filesys.DosageData;

/**
 * Implements {@link MatrixDataLoading} to load plink files to a {@link RealMatrix}, using -1,0,1,2
 */
public class PlinkDataLoader implements MatrixDataLoading {

  private final String dir;
  private final String plinkRoot;
  private final Logger log;

  public PlinkDataLoader(String dir, String plinkRoot, Logger log) {

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
      samples[i] = d.getIds()[0] + "\t" + d.getIds()[1];
    }

    String[] markers = d.getMarkerSet().getMarkerNames();

    log.reportTimeInfo("Preparing matrix for " + markers.length + " markers and " + samples.length
                       + " samples");
    RealMatrix m = MatrixUtils.createRealMatrix(markers.length, samples.length);
    for (int column = 0; column < m.getColumnDimension(); column++) {
      for (int row = 0; row < m.getRowDimension(); row++) {
        m.setEntry(row, column, d.getDosageValues()[row][column]);
      }
    }
    return new NamedRealMatrix(markers, samples, m);
  }

}
