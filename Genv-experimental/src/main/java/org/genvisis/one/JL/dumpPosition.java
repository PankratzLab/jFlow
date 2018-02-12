package org.genvisis.one.JL;

import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.Project;

public class dumpPosition {

  public static void main(String[] args) {
    Project proj = new Project("/home/pankrat2/lanej/projects/FHS.properties");
    MarkerDetailSet markerSet = proj.getMarkerSet();
    markerSet.exportToText(proj, "/home/pankrat2/lanej/test/FHSmarkerPositions.txt");
  }
}
