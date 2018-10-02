package org.genvisis.cnv.filesys;

import org.pankratzlab.common.Logger;
import org.pankratzlab.shared.filesys.Segment;

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
