/**
 * 
 */
package org.genvisis.cnv.analysis.pca.ancestry;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.genvisis.common.Logger;
import org.genvisis.common.RealMatrixUtils.NamedRealMatrix;
import org.genvisis.filesys.DosageData;

/**
 * Implements {@link AncestryDataLoading} to load plink files to a {@link RealMatrix}
 */
public class PlinkDataLoader implements AncestryDataLoading {

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

    System.out.println(d.getDosageValues().length);
    System.out.println(d.getDosageValues()[0].length);

    String[] samples = new String[d.getIds().length];
    for (int i = 0; i < samples.length; i++) {
      samples[i] = d.getIds()[0] + "\t" + d.getIds()[1];
    }

    String[] markers = d.getMarkerSet().getMarkerNames();

    RealMatrix m = MatrixUtils.createRealMatrix(markers.length, samples.length);
    for (int column = 0; column < m.getColumnDimension(); column++) {
      for (int row = 0; row < m.getRowDimension(); row++) {
        m.setEntry(row, column, d.getDosageValues()[row][column]);
      }
    }

    return new NamedRealMatrix(markers, samples, m);
  }

}
