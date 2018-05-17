/**
 * 
 */
package org.genvisis.cnv.analysis.collapse;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.genvisis.cnv.analysis.BeastScore.BeastVariant;
import org.genvisis.cnv.analysis.collapse.LRRForceCaller.LRRRegion;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.var.MosaicRegion;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Logger;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import com.google.common.collect.ImmutableSet;

/**
 * Methods to force call a region of the genome for copy number status
 */
class CNVForceCaller {

  static final List<METHODS> VALID_METHODS = Arrays.asList(new METHODS[] {METHODS.BEAST,
                                                                          METHODS.MOSAIC,
                                                                          METHODS.LRR});

  private LocusSet<MosaicRegion> mosaicResults;
  private LocusSet<BeastVariant> beastResults;
  private LocusSet<LRRRegion> lrrResults;
  private String dna;

  // private
  private CNVForceCaller() {
    //
  }

  void dumpToFile(PrintWriter writer, boolean writeHeader, Logger log) {
    if (writeHeader) {
      String[] header = ArrayUtils.concatAll(mosaicResults.getLoci()[0].getHeader(),
                                             beastResults.getLoci()[0].getHeader(),
                                             lrrResults.getLoci()[0].getHeader());

      writer.println(ArrayUtils.toStr(header));
    }
    if (mosaicResults.getLoci().length == beastResults.getLoci().length
        && beastResults.getLoci().length == lrrResults.getLoci().length) {
      for (int i = 0; i < mosaicResults.getLoci().length; i++) {
        writer.println(mosaicResults.getLoci()[i].toAnalysisString() + "\t"
                       + beastResults.getLoci()[i].toAnalysisString() + "\t"
                       + lrrResults.getLoci()[i].toAnalysisString());
      }
    } else {
      throw new IllegalArgumentException("Invalid loci length");
    }
  }

  enum METHODS {
    BEAST, PENNCNV, MOSAIC, LRR;
  }

  static <T extends Segment> CNVForceCaller call(Project proj, Sample sample, LocusSet<T> regions,
                                                 List<METHODS> methods, boolean[] markersToUse,
                                                 LocusSet<Segment> distributionalMosExcludes,
                                                 int[][] indicesByChr) {

    if (!VALID_METHODS.containsAll(methods)) {
      List<METHODS> invalidMethods = new ArrayList<>(methods);
      invalidMethods.removeAll(VALID_METHODS);

    }

    int[][] autosomalIndices = new int[][] {proj.getAutosomalMarkerIndices()};
    CNVForceCaller called = new CNVForceCaller();
    called.dna = sample.getSampleName();
    for (METHODS method : methods) {
      proj.getLog()
          .reportTimeInfo("Running method " + method + " for sample " + sample.getSampleName());
      switch (method) {
        case BEAST:
          boolean[] beastMarkers = ArrayUtils.booleanArray(markersToUse.length, false);
          proj.getLog().reportTimeInfo("Starting with " + ArrayUtils.booleanArraySum(markersToUse)
                                       + " markers to use");
          if (proj.ARRAY_TYPE.getValue() == ARRAY.NGS) {
            proj.getLog().reportTimeInfo("subsetting to CN only for " + proj.ARRAY_TYPE.getValue());
            boolean[] tmp = proj.getCNMarkers();
            for (int i = 0; i < tmp.length; i++) {
              beastMarkers[i] = tmp[i] && markersToUse[i];
            }
          }
          proj.getLog().reportTimeInfo(+ArrayUtils.booleanArraySum(beastMarkers)
                                       + " CN only  markers to use");
          BeastForceCaller beastForceCaller = new BeastForceCaller(proj, sample.getSampleName(),
                                                                   autosomalIndices,
                                                                   ArrayUtils.toDoubleArray(sample.getLRRs()),
                                                                   beastMarkers);
          beastForceCaller.forceCall(regions);
          called.beastResults = beastForceCaller.getResults();
          break;
        case MOSAIC:
          ImmutableSet.Builder<Marker> useMarkers = ImmutableSet.builder();
          List<Marker> projOrderMarkers = proj.getMarkerSet().getMarkers();
          for (int i = 0; i < markersToUse.length; i++) {
            if (markersToUse[i]) useMarkers.add(projOrderMarkers.get(i));
          }
          MosaicForceCaller mosaicForceCaller = new MosaicForceCaller(proj, sample.getSampleName(),
                                                                      sample.markerBAFMap(proj.getMarkerSet()),
                                                                      useMarkers.build());
          mosaicForceCaller.setDistributionalExcludes(distributionalMosExcludes);
          mosaicForceCaller.forceCall(regions);
          called.mosaicResults = mosaicForceCaller.getResults();
          break;
        case LRR:
          LRRForceCaller lrrForceCaller = new LRRForceCaller(proj, sample.getSampleName(),
                                                             indicesByChr,
                                                             ArrayUtils.toDoubleArray(sample.getLRRs()),
                                                             markersToUse);

          lrrForceCaller.forceCall(regions);
          called.lrrResults = lrrForceCaller.getResults();
          break;
        case PENNCNV:

        default:
          throw new IllegalArgumentException("Invalid methods provided " + methods);

      }
    }
    return called;
  }
}
