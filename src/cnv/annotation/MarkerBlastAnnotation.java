package cnv.annotation;

import filesys.Segment;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import cnv.annotation.BlastAnnotationTypes.BLAST_ANNOTATION_TYPES;
import cnv.annotation.BlastAnnotationTypes.BlastAnnotation;
import cnv.filesys.Project;
import common.Logger;
import common.ArraySpecialList.ArrayBlastAnnotationList;

public class MarkerBlastAnnotation implements AnnotationParser {

	private BLAST_ANNOTATION_TYPES[] bTypes;
	private ArrayBlastAnnotationList[] annotationLists;
	private MarkerBlastHistogramAnnotation blastAlignmentHistogram;
	private MarkerSeqAnnotation markerSeqAnnotation;
	private String markerName;
	private boolean found;

	public MarkerBlastAnnotation(String markerName, BLAST_ANNOTATION_TYPES[] bTypes, int initialCapacity) {
		super();
		this.markerName = markerName;
		this.bTypes = bTypes;
		this.annotationLists = new ArrayBlastAnnotationList[bTypes.length];
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
			log.reportTimeError(error);
			throw new IllegalStateException(error);
		}
		return index;
	}

	public int[] getAlignmentHistogram(Project proj) {
		return blastAlignmentHistogram.formatCountsForProject(proj);
	}

	public MarkerSeqAnnotation getMarkerSeqAnnotation() {
		return markerSeqAnnotation;
	}

	@Override
	public void parseAnnotation(VariantContext vc, Logger log) {
		for (int i = 0; i < BLAST_ANNOTATION_TYPES.values().length; i++) {// each annotation type has a separate key in the file
			String info = vc.getCommonInfo().getAttributeAsString(BLAST_ANNOTATION_TYPES.values()[i].toString(), BLAST_ANNOTATION_TYPES.values()[i].getDefaultValue());
			if (!info.equals(BLAST_ANNOTATION_TYPES.values()[i].getDefaultValue())) {
				List<String> groups = Arrays.asList(info.replaceAll("\\[", "").replaceAll("\\]", "").split("\\s*,\\s*"));
				for (String group : groups) {
					String[] segCigarStrand = group.split("/");
					annotationLists[i].add(new BlastAnnotation(TextCigarCodec.decode(segCigarStrand[0]), new Segment(segCigarStrand[1]), Strand.toStrand(segCigarStrand[2])));
				}
			}
		}
		this.blastAlignmentHistogram = new MarkerBlastHistogramAnnotation(MarkerBlastHistogramAnnotation.DEFAULT_NAME, MarkerBlastHistogramAnnotation.DEFAULT_DESCRIPTION, null);
		blastAlignmentHistogram.parseAnnotation(vc, log);
		this.markerSeqAnnotation = new MarkerSeqAnnotation();
		markerSeqAnnotation.parseAnnotation(vc, log);
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

	public static MarkerBlastAnnotation[] initForMarkers(final String[] markers) {
		MarkerBlastAnnotation[] blastAnnotations = new MarkerBlastAnnotation[markers.length];
		for (int i = 0; i < blastAnnotations.length; i++) {
			blastAnnotations[i] = new MarkerBlastAnnotation(markers[i], BLAST_ANNOTATION_TYPES.values(), 100);
		}
		return blastAnnotations;
	}

}
