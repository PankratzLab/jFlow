package cnv.annotation;

import filesys.Segment;
import htsjdk.samtools.Cigar;
import htsjdk.variant.vcf.VCFHeaderLineType;

import java.util.ArrayList;

import cnv.filesys.Project;
import common.Logger;
import seq.analysis.Blast.BlastResults;

public class BlastAnnotationTypes {

	private static final String ANNO_DELIM = ",";
	private static final String DEFAULT_VALUE = ".";
	private static final String DEFAULT_VALUE_PM = ".";

	/**
	 * @author lane0212
	 * Define a key value entry method in the vcf
	 */
	public enum BLAST_ANNOTATION_TYPES {
		// /**
		// * alignment without gaps
		// */
		// OFF_T_ALIGNMENT_GAPS_AND_MISMATCHES(VCFHeaderLineType.String, "OFF_T_ALIGNMENT_GAPS_AND_MISMATCHES", "Off-target Alignments with mismatches AND gaps (" + ANNO_DELIM + "-delimited)", ANNO_DELIM, DEFAULT_VALUE),
		// /**
		// *
		// * alignment without gaps
		// */
		// OFF_T_ALIGNMENT_NO_GAPS(VCFHeaderLineType.String, "OFF_T_ALIGNMENT_NO_GAPS", "Off-target Alignments with mismatches, but without gaps (" + ANNO_DELIM + "-delimited)", ANNO_DELIM, DEFAULT_VALUE), /**
		// * alignment without mismatches
		// */
		// OFF_T_ALIGNMENT_NO_MISMATCHES(VCFHeaderLineType.String, "OFF_T_ALIGNMENT_NO_MISMATCHES", "Off-target Alignments with gaps, but without mismatches (" + ANNO_DELIM + "-delimited)", ANNO_DELIM, DEFAULT_VALUE), /**
		// * alignment without mismatches or gaps
		// */
		// OFF_T_ALIGNMENT_NO_MISMATCHES_NO_GAPS(VCFHeaderLineType.String, "OFF_T_ALIGNMENT_NO_MISMATCHES_NO_GAPS", "Off-target Alignments without gaps AND without mismatches (" + ANNO_DELIM + "-delimited)", ANNO_DELIM, DEFAULT_VALUE),

		OFF_T_ALIGNMENTS(VCFHeaderLineType.String, "Off-target alignments sorted by continous reference alignment length (" + ANNO_DELIM + "-delimited)", ANNO_DELIM, DEFAULT_VALUE),

		ON_T_ALIGNMENTS_NON_PERFECT(VCFHeaderLineType.String, "On-target alignments (but not perfectly matched) sorted by continous reference alignment length (" + ANNO_DELIM + "-delimited)", ANNO_DELIM, DEFAULT_VALUE),

		/**
		 * 
		 */
		PERFECT_MATCH(VCFHeaderLineType.String, "There is a perfect match Alignment", ANNO_DELIM, DEFAULT_VALUE_PM), /**
		 * 
		 */
		// NON_PERFECT_MATCH_ON_TARGET(VCFHeaderLineType.String, "HAS_PERFECT_MATCH", "There is a perfect match Alignment", ANNO_DELIM, DEFAULT_VALUE_PM),

		;

		private String name;
		private String Description;
		private VCFHeaderLineType vType;
		private String sep;
		private String defaultValue;

		/**
		 * @param vType this is not really used on our end
		 * @param Description
		 * @param sep
		 * @param defaultValue
		 */
		private BLAST_ANNOTATION_TYPES(VCFHeaderLineType vType, String Description, String sep, String defaultValue) {
			this.vType = vType;
			this.name = this.toString();
			this.Description = Description;
			this.sep = sep;
			this.defaultValue = defaultValue;
		}

		public String getName() {
			return name;
		}

		public String getDescription() {
			return Description;
		}

		public VCFHeaderLineType getvType() {
			return vType;
		}

		public String getDefaultValue() {
			return defaultValue;
		}

		public String getSep() {
			return sep;
		}
	}

	public static boolean shouldBeAnnotatedAs(Project proj, BlastResults blastResults, BLAST_ANNOTATION_TYPES type, Segment seg, Logger log) {
		boolean annotatedAs = false;
		switch (type) {
		case OFF_T_ALIGNMENTS:
			annotatedAs = !blastResults.getSegment().overlaps(seg);
			break;
		case ON_T_ALIGNMENTS_NON_PERFECT:
			if (blastResults.getSegment().overlaps(seg)) {
				if (blastResults.getMismatches() > 0 || blastResults.getGapOpens() > 0 || blastResults.getAlignmentLength() != proj.getArrayType().getProbeLength()) {
					annotatedAs = true;
				}
			}
			break;
		case PERFECT_MATCH:
			annotatedAs = blastResults.getMismatches() == 0 && blastResults.getGapOpens() == 0 && blastResults.getAlignmentLength() == proj.getArrayType().getProbeLength() && blastResults.getSegment().overlaps(seg);
			break;
		default:
			log.reportTimeError("Invalid annotation test " + type);
			break;
		}
		return annotatedAs;
	}

	public static Annotation[] getBaseAnnotations() {

		ArrayList<Annotation> annotations = new ArrayList<Annotation>();
		for (int i = 0; i < BLAST_ANNOTATION_TYPES.values().length; i++) {
			annotations.add(new Annotation(BLAST_ANNOTATION_TYPES.values()[i].getvType(), BLAST_ANNOTATION_TYPES.values()[i].getName(), BLAST_ANNOTATION_TYPES.values()[i].getDescription()) {
			});
		}
		return annotations.toArray(new Annotation[annotations.size()]);
	}

	public static AnnotationData[] getAnnotationDatas() {
		ArrayList<AnnotationData> annotations = new ArrayList<AnnotationData>();
		for (int i = 0; i < BLAST_ANNOTATION_TYPES.values().length; i++) {
			annotations.add(new AnnotationData(BLAST_ANNOTATION_TYPES.values()[i].getvType(), BLAST_ANNOTATION_TYPES.values()[i].getName(), BLAST_ANNOTATION_TYPES.values()[i].getDescription(), BLAST_ANNOTATION_TYPES.values()[i].getDefaultValue()) {
			});
		}
		return annotations.toArray(new AnnotationData[annotations.size()]);
	}

	/**
	 * @author lane0212 Stores the {@link Cigar} string for the blast hit, and the {@link Segment} on the reference genome
	 */
	public static class BlastAnnotation {
		private Cigar cigar;
		private Segment refLoc;

		public BlastAnnotation(Cigar cigar, Segment refLoc) {
			super();
			this.cigar = cigar;
			this.refLoc = refLoc;
		}

		public Cigar getCigar() {
			return cigar;
		}

		public Segment getRefLoc() {
			return refLoc;
		}

		public static String[] toAnnotationString(BlastAnnotation[] blastAnnotations) {
			String[] annotations = new String[blastAnnotations.length];
			for (int i = 0; i < annotations.length; i++) {
				annotations[i] = blastAnnotations[i].getCigar().toString() + "/" + blastAnnotations[i].getRefLoc().getUCSClocation();
			}
			return annotations;
		}
	}
}
