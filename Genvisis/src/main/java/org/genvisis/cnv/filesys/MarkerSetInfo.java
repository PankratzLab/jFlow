package org.genvisis.cnv.filesys;

import org.genvisis.common.Logger;
import org.genvisis.filesys.Segment;

public interface MarkerSetInfo {

  @Deprecated
  String[] getMarkerNames();

  @Deprecated
  byte[] getChrs();

  @Deprecated
  int[] getPositions();

  @Deprecated
  int[][] getPositionsByChr();

  @Deprecated
  int[][] getIndicesByChr();

  @Deprecated
  int[] getIndicesOfMarkersIn(Segment seg, int[][] indicesByChr, Logger log);

  @Deprecated
  String[] getMarkersIn(Segment seg, int[][] indicesByChr);

  boolean checkFingerprint(Sample samp);

  long getFingerprint();

}
