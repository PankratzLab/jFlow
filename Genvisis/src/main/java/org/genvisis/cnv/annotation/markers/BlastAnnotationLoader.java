package org.genvisis.cnv.annotation.markers;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.BLAST_ANNOTATION_TYPES;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.BlastAnnotation;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.PROBE_TAG;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerSet;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.ArraySpecialList.ArrayBlastAnnotationList;
import org.pankratzlab.shared.filesys.Segment;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * @author lane0212 Loads summarized blast results in annotation format
 */
public class BlastAnnotationLoader extends AnnotationFileLoader {

  private final byte[] chrs;
  private final int[] pos;
  private final Map<String, Integer> markerIndices;

  /**
   * @param annotationFilename the file created by a {@link BlastAnnotationWriter} run
   * @param markerSet {@link MarkerSet} to load
   * @param indexRequired , should always be set to true
   * @param log
   */
  public BlastAnnotationLoader(String annotationFilename, MarkerDetailSet markerSet,
                               boolean indexRequired, Logger log) {
    super(null, BlastAnnotationTypes.getBaseAnnotations(), annotationFilename, indexRequired, log);
    chrs = markerSet.getChrs();
    pos = markerSet.getPositions();
    markerIndices = markerSet.getMarkerIndices();
  }

  /**
   * @param markers
   * @param otherQueries these queries can be null, but if not the appropriate annotations will be
   *          parsed
   * @return
   */
  public MarkerBlastResult[] loadBlastAnnotationsFor(String[] markers,
                                                     AnnotationParser[]... otherQueries) {

    if (ArrayUtils.unique(markers).length != markers.length) {
      String error = "Internal error, markers for blast annotation retrieval must be unique";
      log.reportError(error);
      throw new IllegalArgumentException(error);
    }

    MarkerBlastResult[] blastAnnotations = initResults(markers);

    Segment[] segs = getSegmentsForMarkers(markers);

    AnnotationQuery annotationQuery = getAnnotationQuery(segs);

    boolean[] found = ArrayUtils.booleanArray(markers.length, false);

    int count = 0;
    while (annotationQuery.hasNext()) {
      count++;
      if (count % 10000 == 0) {
        log.reportTimeInfo("Loaded " + count + " annotations");
      }
      VariantContext vc = annotationQuery.next();
      if (otherQueries != null) {
        for (AnnotationParser[] annotationParsers : otherQueries) {
          for (AnnotationParser annotationParser : annotationParsers) {
            if (annotationParser.shouldAnnotateWith(vc, log)) {
              annotationParser.parseAnnotation(vc, log);
            }
          }
        }
      }
      String id = vc.getID();
      int annoIndex = ext.indexOfStr(id, markers);
      if (annoIndex < 0) {
        log.reportTimeWarning("Query has returned un-desired marker " + id + ", ignoring");
      } else {
        found[annoIndex] = true;
        blastAnnotations[annoIndex].parseAnnotation(vc, log);
      }
    }
    if (ArrayUtils.booleanArraySum(found) != markers.length) {
      String error = markers.length + " markers were expected to be loaded, but only "
                     + ArrayUtils.booleanArraySum(found) + " markers were found";
      for (int i = 0; i < found.length; i++) {
        if (!found[i]) {
          error += "\nMissing " + markers[i];
        }
      }
      log.reportError(error);
      log.reportError(ArrayUtils.toStr(markers));
      // throw new IllegalStateException(error);
    } else {
      log.reportTimeInfo("Loaded " + markers.length + " marker annotations");
    }
    return blastAnnotations;
  }

  private MarkerBlastResult[] initResults(String[] markers) {
    MarkerBlastResult[] blastAnnotations = new MarkerBlastResult[markers.length];
    for (int i = 0; i < blastAnnotations.length; i++) {
      blastAnnotations[i] = new MarkerBlastResult(markers[i], BLAST_ANNOTATION_TYPES.values(), 100);
    }
    return blastAnnotations;
  }

  private Segment[] getSegmentsForMarkers(final String[] markers) {
    Segment[] segs = new Segment[markers.length];
    for (int i = 0; i < segs.length; i++) {
      int markerIndex = markerIndices.get(markers[i]);
      Segment markerSeg = new Segment(chrs[markerIndex], pos[markerIndex], pos[markerIndex]);
      segs[i] = markerSeg;
    }
    return segs;
  }

  /**
   * @author lane0212 Probably could be its own class
   */
  public static class MarkerBlastResult implements AnnotationParser {

    private final BLAST_ANNOTATION_TYPES[] bTypes;
    private final ArrayBlastAnnotationList[] annotationLists;
    private final String markerName;
    private boolean found;

    public MarkerBlastResult(String markerName, BLAST_ANNOTATION_TYPES[] bTypes,
                             int initialCapacity) {
      super();
      this.markerName = markerName;
      this.bTypes = bTypes;
      annotationLists = new ArrayBlastAnnotationList[bTypes.length];
      for (int i = 0; i < annotationLists.length; i++) {
        annotationLists[i] = new ArrayBlastAnnotationList(initialCapacity);
      }
    }

    public boolean hasPerfectMatch(Logger log) {
      return getAnnotationsFor(BLAST_ANNOTATION_TYPES.PERFECT_MATCH, log).size() > 0;
    }

    public int getNumOffTarget(Logger log) {
      return getAnnotationsFor(BLAST_ANNOTATION_TYPES.OFF_T_ALIGNMENTS, log).size();
    }

    public ArrayList<BlastAnnotation> getAnnotationsFor(BLAST_ANNOTATION_TYPES btype, Logger log) {
      int testIndex = getAnnotationIndexFor(btype, log);
      return annotationLists[testIndex];
    }

    public BLAST_ANNOTATION_TYPES[] getbTypes() {
      return bTypes;
    }

    public ArrayBlastAnnotationList[] getAnnotationLists() {
      return annotationLists;
    }

    private int getAnnotationIndexFor(BLAST_ANNOTATION_TYPES bType, Logger log) {
      int index = -1;
      for (int i = 0; i < bTypes.length; i++) {
        if (bTypes[i] == bType) {
          index = i;
          break;
        }
      }
      if (index < 0) {
        String error = "Internal error: Annotation does not contain " + bType;
        log.reportError(error);
        throw new IllegalStateException(error);
      }
      return index;
    }

    @Override
    public void parseAnnotation(VariantContext vc, Logger log) {
      for (int i = 0; i < BLAST_ANNOTATION_TYPES.values().length; i++) {// each annotation type has
                                                                        // a separate key in the
                                                                        // file
        String info = vc.getCommonInfo()
                        .getAttributeAsString(BLAST_ANNOTATION_TYPES.values()[i].toString(),
                                              BLAST_ANNOTATION_TYPES.values()[i].getDefaultValue());
        if (!info.equals(BLAST_ANNOTATION_TYPES.values()[i].getDefaultValue())) {
          List<String> groups = Arrays.asList(info.replaceAll("\\[", "").replaceAll("\\]", "")
                                                  .split("\\s*,\\s*"));
          for (String group : groups) {
            String[] segCigarStrand = group.split("/");
            annotationLists[i].add(new BlastAnnotation(TextCigarCodec.decode(segCigarStrand[0]),
                                                       new Segment(segCigarStrand[1]),
                                                       Strand.valueOf(segCigarStrand[2]),
                                                       PROBE_TAG.valueOf(segCigarStrand[3]),
                                                       Double.parseDouble(segCigarStrand[4])));
          }
        }
      }
    }

    @Override
    public boolean shouldAnnotateWith(VariantContext vc, Logger log) {
      return markerName.equals(vc.getID());
    }

    @Override
    public void setFound(boolean found) {
      this.found = found;
    }

    @Override
    public boolean isFound() {
      return found;
    }

  }
}
