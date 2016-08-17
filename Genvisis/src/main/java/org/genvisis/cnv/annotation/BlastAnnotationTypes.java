package org.genvisis.cnv.annotation;

import java.util.ArrayList;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.util.CNVHelper;
import org.genvisis.common.Logger;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.analysis.Blast.BlastResults;

import htsjdk.samtools.Cigar;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;

public class BlastAnnotationTypes {

  /**
   * @author lane0212 Define a key value entry method in the vcf
   */
  public enum BLAST_ANNOTATION_TYPES {
                                      // /**
                                      // * alignment without gaps
                                      // */
                                      // OFF_T_ALIGNMENT_GAPS_AND_MISMATCHES(VCFHeaderLineType.String,
                                      // "OFF_T_ALIGNMENT_GAPS_AND_MISMATCHES", "Off-target
                                      // Alignments with mismatches AND gaps (" +
                                      // ANNO_DELIM + "-delimited)", ANNO_DELIM, DEFAULT_VALUE),
                                      // /**
                                      // *
                                      // * alignment without gaps
                                      // */
                                      // OFF_T_ALIGNMENT_NO_GAPS(VCFHeaderLineType.String,
                                      // "OFF_T_ALIGNMENT_NO_GAPS", "Off-target
                                      // Alignments with mismatches, but without gaps (" +
                                      // ANNO_DELIM + "-delimited)", ANNO_DELIM,
                                      // DEFAULT_VALUE), /**
                                      // * alignment without mismatches
                                      // */
                                      // OFF_T_ALIGNMENT_NO_MISMATCHES(VCFHeaderLineType.String,
                                      // "OFF_T_ALIGNMENT_NO_MISMATCHES",
                                      // "Off-target Alignments with gaps, but without mismatches ("
                                      // + ANNO_DELIM + "-delimited)",
                                      // ANNO_DELIM, DEFAULT_VALUE), /**
                                      // * alignment without mismatches or gaps
                                      // */
                                      // OFF_T_ALIGNMENT_NO_MISMATCHES_NO_GAPS(VCFHeaderLineType.String,
                                      // "OFF_T_ALIGNMENT_NO_MISMATCHES_NO_GAPS", "Off-target
                                      // Alignments without gaps AND without
                                      // mismatches (" + ANNO_DELIM + "-delimited)", ANNO_DELIM,
                                      // DEFAULT_VALUE),

                                      OFF_T_ALIGNMENTS(VCFHeaderLineType.String, DEFAULT_COUNT,
                                                       DEFAULT_NUMBER,
                                                       "Off-target alignments (Cigar-positive orientation, Segment - positive orientation, Strand -original orientation) sorted by continous reference alignment length ("
                                                                       + ANNO_DELIM + "-delimited)",
                                                       ANNO_DELIM, DEFAULT_VALUE),

                                      ON_T_ALIGNMENTS_NON_PERFECT(VCFHeaderLineType.String,
                                                                  DEFAULT_COUNT, DEFAULT_NUMBER,
                                                                  "On-target alignments (but not perfectly matched)(Cigar-positive orientation, Segment - positive orientation, Strand -original orientation) sorted by continous reference alignment length ("
                                                                                                 + ANNO_DELIM
                                                                                                 + "-delimited)",
                                                                  ANNO_DELIM, DEFAULT_VALUE),

                                      /**
                                       * 
                                       */
                                      PERFECT_MATCH(VCFHeaderLineType.String, DEFAULT_COUNT, DEFAULT_NUMBER, "Alignments (Cigar-positive orientation, Segment - positive orientation, Strand -original orientation) that perfectly cover the region", ANNO_DELIM, DEFAULT_VALUE),
    /**
    * 
    */
    // NON_PERFECT_MATCH_ON_TARGET(VCFHeaderLineType.String, "HAS_PERFECT_MATCH", "There is a
    // perfect match Alignment", ANNO_DELIM, DEFAULT_VALUE_PM),

    ;

    private final String name;
    private final int number;
    private final String Description;
    private final VCFHeaderLineType vType;
    private final VCFHeaderLineCount count;
    private final String sep;
    private final String defaultValue;

    /**
     * @param vType this is not really used on our end
     * @param Description
     * @param sep
     * @param defaultValue
     */
    private BLAST_ANNOTATION_TYPES(VCFHeaderLineType vType, VCFHeaderLineCount count, int number,
                                   String Description, String sep, String defaultValue) {
      this.vType = vType;
      this.count = count;
      this.number = number;
      name = toString();
      this.Description = Description;
      this.sep = sep;
      this.defaultValue = defaultValue;
    }

    public VCFHeaderLineCount getCount() {
      return count;
    }

    public String getDefaultValue() {
      return defaultValue;
    }

    public String getDescription() {
      return Description;
    }

    public String getName() {
      return name;
    }

    public int getNumber() {
      return number;
    }

    public String getSep() {
      return sep;
    }

    public VCFHeaderLineType getvType() {
      return vType;
    }

  }
  /**
   * @author lane0212 Stores the {@link Cigar} string for the blast hit, and the {@link Segment} on
   *         the reference genome
   */
  public static class BlastAnnotation {
    public static String[] toAnnotationString(BlastAnnotation[] blastAnnotations) {
      String[] annotations = new String[blastAnnotations.length];
      for (int i = 0; i < annotations.length; i++) {
        BlastAnnotation tmp = blastAnnotations[i];
        annotations[i] = tmp.getCigar().toString() + "/" + tmp.getRefLoc().getUCSClocation() + "/"
                         + CNVHelper.decode(tmp.getStrand()) + "/" + tmp.getTag() + "/"
                         + tmp.geteValue();
      }
      return annotations;
    }

    private final Cigar cigar; // positive strand
    private final Segment refLoc;// positvie strand
    private final Strand strand;// original alignment strand
    private final PROBE_TAG tag;// A probe, B probe, or Both

    private final double eValue;

    public BlastAnnotation(Cigar cigar, Segment refLoc, Strand strand, PROBE_TAG tag,
                           double eValue) {
      super();
      this.cigar = cigar;
      this.refLoc = refLoc;
      this.strand = strand;
      this.tag = tag;
      this.eValue = eValue;
    }

    public Cigar getCigar() {
      return cigar;
    }

    public double geteValue() {
      return eValue;
    }

    public Segment getRefLoc() {
      return refLoc;
    }

    public Strand getStrand() {
      return strand;
    }

    public PROBE_TAG getTag() {
      return tag;
    }
  }
  public enum PROBE_TAG {

                         A("_A"), B("_B"),
                         /**
                          * When both the A and B probes have identical sequences
                          */
                         BOTH("_AB");
    public static PROBE_TAG parseMarkerTag(String markerWithTag, Logger log) {
      return markerWithTag.endsWith(PROBE_TAG.BOTH.getTag()) ? PROBE_TAG.toTag(markerWithTag.substring(markerWithTag.length()
                                                                                                       - PROBE_TAG.BOTH.getTag()
                                                                                                                       .length()),
                                                                               log)
                                                             : PROBE_TAG.toTag(markerWithTag.substring(markerWithTag.length()
                                                                                                       - 2),
                                                                               log);
    }

    public static PROBE_TAG toTag(String tag, Logger log) {
      if (tag.equals(A.getTag())) {
        return A;
      } else if (tag.equals(B.getTag())) {
        return B;
      } else if (tag.equals(BOTH.getTag())) {
        return BOTH;
      } else {
        String error = "Invalid tag instance " + tag;
        log.reportTimeError(error);
        throw new IllegalArgumentException(error);
      }
    }

    private final String tag;

    private PROBE_TAG(String tag) {
      this.tag = tag;
    }

    public String getTag() {
      return tag;
    }
  }
  public enum TOP_BOT {
                       TOP, BOT, NA,
                       /**
                        * For cnvi
                        */
                       PLUS,
                       /**
                        * For cnvi and indels
                        */
                       MINUS;
  }

  private static final String ANNO_DELIM = ",";

  private static final String DEFAULT_VALUE = ".";

  private static final int DEFAULT_NUMBER = 1;

  private static final VCFHeaderLineCount DEFAULT_COUNT = VCFHeaderLineCount.UNBOUNDED; // unknown
                                                                                        // number of
                                                                                        // entries

  public static AnnotationData[] getAnnotationDatas() {
    ArrayList<AnnotationData> annotations = new ArrayList<AnnotationData>();
    for (int i = 0; i < BLAST_ANNOTATION_TYPES.values().length; i++) {
      BLAST_ANNOTATION_TYPES btype = BLAST_ANNOTATION_TYPES.values()[i];
      annotations.add(new AnnotationData(btype.getvType(), btype.getCount(), btype.getNumber(),
                                         btype.getName(), btype.getDescription(),
                                         btype.getDefaultValue(), btype.getDefaultValue()) {});
    }
    return annotations.toArray(new AnnotationData[annotations.size()]);
  }

  public static Annotation[] getBaseAnnotations() {

    ArrayList<Annotation> annotations = new ArrayList<Annotation>();
    for (int i = 0; i < BLAST_ANNOTATION_TYPES.values().length; i++) {
      BLAST_ANNOTATION_TYPES btype = BLAST_ANNOTATION_TYPES.values()[i];
      annotations.add(new Annotation(btype.getvType(), btype.getCount(), btype.getNumber(),
                                     btype.getName(), btype.getDescription(),
                                     btype.getDefaultValue()) {});
    }
    return annotations.toArray(new Annotation[annotations.size()]);
  }

  public static boolean shouldBeAnnotatedAs(Project proj, BlastResults blastResults,
                                            BLAST_ANNOTATION_TYPES type, Segment seg, Logger log) {
    boolean annotatedAs = false;
    switch (type) {
      case OFF_T_ALIGNMENTS:
        annotatedAs = !blastResults.getSegment().overlaps(seg);
        break;
      case ON_T_ALIGNMENTS_NON_PERFECT:
        if (blastResults.getSegment().overlaps(seg)) {
          if (blastResults.getMismatches() > 0 || blastResults.getGapOpens() > 0
              || blastResults.getAlignmentLength() != proj.getArrayType().getProbeLength()) {
            annotatedAs = true;
          }
        }
        break;
      case PERFECT_MATCH:
        annotatedAs = blastResults.getMismatches() == 0 && blastResults.getGapOpens() == 0
                      && blastResults.getAlignmentLength() == proj.getArrayType().getProbeLength()
                      && blastResults.getSegment().overlaps(seg);
        break;
      default:
        log.reportTimeError("Invalid annotation test " + type);
        break;
    }
    return annotatedAs;
  }
}
