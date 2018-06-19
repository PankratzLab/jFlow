
package org.genvisis.cnv.annotator;

import java.util.ArrayList;
import java.util.List;
import org.genvisis.cnv.analysis.BeastScore;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.filesys.CNVariant;
import org.genvisis.seq.manage.BamImport.NGS_MARKER_TYPE;
import com.google.common.primitives.Ints;

/**
 * {@link Annotator} that computes {@link BeastScore} annotations
 */
public class BeastScoreAnnotator extends AbstractAnnotator<CNVariant> {

  private static final String BEAST_SCORE_ANNO = "BEAST_SCORE";
  private static final String BEAST_HEIGHT_ANNO = "BEAST_HEIGHT";
  private static final String BEAST_LENGTH_ANNO = "BEAST_LENGTH";

  private final SampleData sd;
  private final MarkerDetailSet md;
  private final Project proj;
  private final int[][] indicesByChr;

  public BeastScoreAnnotator(Project proj) {
    // The project is used to load lrr data
    this.proj = proj;

    // We use the SampleData to look up sample IDs
    sd = proj.getSampleData(false);

    // And we need to use the markers
    md = proj.getMarkerSet();

    // Cache marker indices which will be needed for computing beast score
    indicesByChr = md.getIndicesByChr();
  }

  @Override
  protected void calcAnnotations(CNVariant segmentToAnnotate, int start, int stop, String suffix) {
    String dna = sd.lookupDNA(segmentToAnnotate.getFamilyID() + "\t"
                              + segmentToAnnotate.getIndividualID());
    // TODO cache samples?
    Sample samp = proj.getFullSampleFromRandomAccessFile(dna);

    Iterable<Marker> namesIn = md.viewMarkersInSeg(segmentToAnnotate.getChr(), start, stop);
    List<Integer> nonVariant = new ArrayList<>();

    // TODO Markers should know their NGS type?
    for (Marker element : namesIn) {
      if (!ARRAY.NGS.equals(proj.ARRAY_TYPE.getValue())
          || !NGS_MARKER_TYPE.VARIANT_SITE.equals(NGS_MARKER_TYPE.getType(element))) {
        nonVariant.add(md.getMarkerIndexMap().get(element));
      }
    }

    // TODO allow calc of beastScore for single segment
    // TODO beastScore should be rewritten to take a MarkerDetailSet
    int[][] cnvIndices = {Ints.toArray(nonVariant)};
    BeastScore beastScore = new BeastScore(samp.getLRRs(), indicesByChr, cnvIndices, proj.getLog());

    beastScore.computeBeastScores();

    // Record the annotations
    addAnnotation(segmentToAnnotate, BEAST_SCORE_ANNO,
                  Float.valueOf(beastScore.getBeastScores()[0]), suffix);
    addAnnotation(segmentToAnnotate, BEAST_HEIGHT_ANNO,
                  Float.valueOf(beastScore.getBeastHeights()[0]), suffix);
    addAnnotation(segmentToAnnotate, BEAST_LENGTH_ANNO,
                  Float.valueOf(beastScore.getBeastLengths()[0]), suffix);
  }

}
