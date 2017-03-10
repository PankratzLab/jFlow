package org.genvisis.cnv.filesys;

import org.genvisis.common.Logger;
import org.genvisis.filesys.Segment;


public interface MarkerSetInfo {

	String[] getMarkerNames();

	byte[] getChrs();

	int[] getPositions();

	int[][] getPositionsByChr();

	int[][] getIndicesByChr();

	int[] getIndicesOfMarkersIn(Segment seg, int[][] indicesByChr, Logger log);

	String[] getMarkersIn(Segment seg, int[][] indicesByChr);

	boolean checkFingerprint(Sample samp);

	long getFingerprint();

}
