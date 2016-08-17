package org.genvisis.one.JL;

import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class DumpMito {

  public static void main(String[] args) {
    MarkerSet markerSet = MarkerSet.load("/Users/Kitty/git/affySnp6/markers.ser", false);
    // int[][] indices = markerSet.getIndicesByChr();
    String[] mitoCurrent = HashVec.loadFileToStringArray("/Users/Kitty/git/affySnp6/MT.119.txt",
                                                         false, new int[] {0}, true);
    int[] indicesM = ext.indexLargeFactors(mitoCurrent, markerSet.getMarkerNames(), true,
                                           new Logger(), true, true);
    String[] mitos = new String[indicesM.length];
    System.out.print(mitos.length);
    for (int i = 0; i < indicesM.length; i++) {
      mitos[i] = markerSet.getMarkerNames()[indicesM[i]] + "\t" + markerSet.getChrs()[indicesM[i]]
                 + "\t" + markerSet.getPositions()[indicesM[i]];

    }
    Files.writeList(mitos, "/Users/Kitty/git/affySnp6/MT.119.txt");// 110 anno, 9 set to 0
  }

}
