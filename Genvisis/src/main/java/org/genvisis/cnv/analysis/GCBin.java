package org.genvisis.cnv.analysis;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.genvisis.cnv.filesys.MarkerSet.PreparedMarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.hmm.PennHmm.ViterbiResult;
import org.genvisis.cnv.qc.GcAdjustor.GcModel;
import org.genvisis.common.Array;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;

import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;

public class GCBin extends Segment {

  /**
   * 
   */
  private static final long serialVersionUID = 1L;
  private final double avgGcContent;
  private final int numMarkers;
  private int binNumber;

  public GCBin(Segment segment, double avgGcContent, int numMarkers, int binNumber) {
    super(segment.getChr(), segment.getStart(), segment.getStop());
    this.avgGcContent = avgGcContent;
    this.numMarkers = numMarkers;
  }

  public double getAvgGcContent() {
    return avgGcContent;
  }

  public int getNumMarkers() {
    return numMarkers;
  }

  public int getBinNumber() {
    return binNumber;
  }

  public static LocusSet<GCBin> bin(Project proj, GcModel gcModel) {// TODO, split and thread by chr

    double[] gcs = gcModel.getGcs();
    new OpdfGaussian(Array.mean(gcs, true), Math.pow(Array.stdev(gcs, true), 2));
    int numStates = 50;
    int zeroState = 5;
    NormalDistribution nd =
        new org.apache.commons.math3.distribution.NormalDistribution(Array.mean(gcs, true),
                                                                     Array.stdev(gcs, true));
    int[] stateSequence = new int[gcs.length];
    for (int i = 0; i < gcs.length; i++) {
      double cdf = nd.cumulativeProbability(gcs[i]);
      int state = (int) Math.round(cdf * numStates) - 1;
      stateSequence[i] = state + zeroState;
    }
    PreparedMarkerSet markerSet = PreparedMarkerSet.getPreparedMarkerSet(proj.getMarkerSet());
    int[][] indices = markerSet.getIndicesByChr();
    ArrayList<GCBin> allGc = new ArrayList<GCBin>();

    for (int i = 0; i < indices.length; i++) {
      if (indices[i].length > 0) {
        int[] currentStates = Array.subArray(stateSequence, indices[i]);
        int[] positions = Array.subArray(markerSet.getPositions(), indices[i]);
        byte chr = (byte) i;
        ViterbiResult vtr = new ViterbiResult(currentStates, null);
        LocusSet<CNVariant> tmp = vtr.analyzeStateSequence(proj, "GC_CONTENT", "chr" + i, chr,
                                                           positions, null, 0, false, true);
        ArrayList<int[]> indicesStates = vtr.getIndexStateChange();
        for (int j = 0; j < tmp.getLoci().length; j++) {
          CNVariant gcv = tmp.getLoci()[j];
          int[] startStop = indicesStates.get(j);
          int[] stateIndices = Arrays.copyOfRange(indices[i], startStop[0], startStop[1] + 1);
          double gc = Array.mean(Array.subArray(gcs, stateIndices), true);
          allGc.add(new GCBin(gcv, gc, gcv.getNumMarkers(), gcv.getCN()));
        }
      }
    }

    LocusSet<GCBin> binset =
        new LocusSet<GCBin>(allGc.toArray(new GCBin[allGc.size()]), true, proj.getLog()) {

          /**
           * 
           */
          private static final long serialVersionUID = 1L;

        };

    return binset;
    // binset.writeRegions(outFile, TO_STRING_TYPE.REGULAR, true, proj.getLog());
  }

}
