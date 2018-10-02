/**
 * 
 */
package org.genvisis.cnv.annotation.segments;

import java.util.ArrayList;
import java.util.List;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.pankratzlab.common.Logger;
import org.pankratzlab.shared.filesys.GeneData;
import org.pankratzlab.shared.filesys.GeneTrack;
import org.pankratzlab.shared.filesys.LocusSet;
import org.pankratzlab.shared.filesys.Segment;

/**
 * @author Kitty
 */
public class GeneAnnotator implements SegmentAnnotator {

  private final LocusSet<GeneData> geneSet;
  private final Logger log;

  /**
   * @param geneSet holding available {@link #GeneData} objects to annotate with
   * @param log
   */
  private GeneAnnotator(LocusSet<GeneData> geneSet, Logger log) {
    super();
    this.geneSet = geneSet;
    this.log = log;
  }

  public LocusSet<GeneData> getGeneSet() {
    return geneSet;
  }

  /**
   * @param build the genome build to annotate with
   * @param log
   * @return
   */
  public static GeneAnnotator getDefaultAnnotator(GENOME_BUILD build, Logger log) {
    String gtrackFile = Resources.genome(build, log).getGTrack().get();

    LocusSet<GeneData> geneSet = GeneTrack.load(gtrackFile).convertToLocusSet(log);

    return new GeneAnnotator(geneSet, log);
  }

  /*
   * (non-Javadoc)
   * @see
   * org.genvisis.cnv.annotation.segments.SegmentAnnotator#annotate(org.genvisis.filesys.Segment,
   * org.genvisis.cnv.annotation.segments.SegmentAnotation)
   */
  @Override
  public void annotate(Segment segment, SegmentAnotation segmentAnotation) {
    if (segmentAnotation.getAttributes().containsKey(SegmentAnnotationKeys.GENE.toString())) {
      String error = "This annotation has already been used for Gene lookup";
      log.reportError(error);
      throw new IllegalArgumentException(error);
    }

    List<String> values = new ArrayList<>();
    GeneData[] geneDatas = geneSet.getOverLappingLoci(segment);

    if (geneDatas == null || geneDatas.length == 0) {
      values.add(SegmentAnnotationKeys.GENE.getMissingValue());
    } else {
      for (GeneData geneData : geneDatas) {
        values.add(geneData.getGeneName());
      }
    }
    segmentAnotation.getAttributes().put(SegmentAnnotationKeys.GENE.toString(), values);
  }

}
