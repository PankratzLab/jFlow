package org.genvisis.cnv.annotation.markers;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.BLAST_ANNOTATION_TYPES;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.BlastAnnotation;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.PROBE_TAG;
import org.genvisis.cnv.annotation.markers.MarkerEvalueHistogramAnnotation.EvalueHistogram;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.ArraySpecialList.ArrayBlastAnnotationList;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.filesys.Segment;

import com.google.common.collect.Maps;

import htsjdk.samtools.TextCigarCodec;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.VariantContext;

public class MarkerBlastAnnotation implements AnnotationParser {

	private final Map<BLAST_ANNOTATION_TYPES, ArrayBlastAnnotationList> annotationLists;
	private MarkerBlastHistogramAnnotation blastAlignmentHistogram;
	private MarkerEvalueHistogramAnnotation markerEvalueHistogramAnnotation;
	private MarkerSeqAnnotation markerSeqAnnotation;
	private final String markerName;
	private boolean found;

	public MarkerBlastAnnotation(String markerName) {
		this(markerName, BLAST_ANNOTATION_TYPES.values(), 100);
	}

	public MarkerBlastAnnotation(String markerName, BLAST_ANNOTATION_TYPES[] bTypes,
															 int initialCapacity) {
		super();
		this.markerName = markerName;
		annotationLists = Maps.newEnumMap(BLAST_ANNOTATION_TYPES.class);
		for (BLAST_ANNOTATION_TYPES bType : bTypes) {
			annotationLists.put(bType, new ArrayBlastAnnotationList(initialCapacity));
		}
	}

	public boolean hasPerfectMatch(Logger log) {
		return !getAnnotationsFor(BLAST_ANNOTATION_TYPES.PERFECT_MATCH, log).isEmpty();
	}

	public int getNumOffTarget(Logger log) {
		return getAnnotationsFor(BLAST_ANNOTATION_TYPES.OFF_T_ALIGNMENTS, log).size();
	}

	public int getNumOnTargetNonPerfect(Logger log) {
		return getAnnotationsFor(BLAST_ANNOTATION_TYPES.ON_T_ALIGNMENTS_NON_PERFECT, log).size();
	}

	public List<BlastAnnotation> getAnnotationsFor(BLAST_ANNOTATION_TYPES bType, Logger log) {
		List<BlastAnnotation> annotations = annotationLists.get(bType);
		if (annotations == null) {
			String error = "Internal error: Annotation does not contain " + bType;
			log.reportError(error);
			throw new IllegalStateException(error);
		}
		return annotations;
	}

	public Set<BLAST_ANNOTATION_TYPES> getbTypes() {
		return annotationLists.keySet();
	}

	public Map<BLAST_ANNOTATION_TYPES, ArrayBlastAnnotationList> getAnnotationLists() {
		return annotationLists;
	}

	public int[] getAlignmentHistogram(Project proj) {
		return blastAlignmentHistogram.formatCountsForProject(proj);
	}

	public MarkerSeqAnnotation getMarkerSeqAnnotation() {
		return markerSeqAnnotation;
	}

	public EvalueHistogram getEvalueHistogram() {
		return markerEvalueHistogramAnnotation.formatHistogram();
	}

	@Override
	public void parseAnnotation(VariantContext vc, Logger log) {
		for (Entry<BLAST_ANNOTATION_TYPES, ArrayBlastAnnotationList> bTypeEntry : annotationLists.entrySet()) {// each
																																																					 // annotation
																																																					 // type
																																																					 // has
																																																					 // a
			// separate key in the file
			BLAST_ANNOTATION_TYPES bType = bTypeEntry.getKey();
			ArrayBlastAnnotationList annotationList = bTypeEntry.getValue();
			String info = vc.getCommonInfo()
											.getAttributeAsString(bType.toString(),
																						bType.getDefaultValue());
			if (!info.equals(bType.getDefaultValue())) {
				List<String> groups = Arrays.asList(info.replaceAll("\\[", "").replaceAll("\\]", "")
																								.split("\\s*,\\s*"));
				for (String group : groups) {
					String[] segCigarStrand = group.split("/");
					int[] location = Positions.parseUCSClocation(segCigarStrand[1]);
					if (location != null) {
						Segment segment = new Segment((byte) location[0], location[1], location[2]);
						annotationList.add(new BlastAnnotation(TextCigarCodec.decode(segCigarStrand[0]),
																									 segment, Strand.toStrand(segCigarStrand[2]),
																									 PROBE_TAG.valueOf(segCigarStrand[3]),
																									 Double.parseDouble(segCigarStrand[4])));
					}
				}
			}
		}
		blastAlignmentHistogram = new MarkerBlastHistogramAnnotation(MarkerBlastHistogramAnnotation.DEFAULT_NAME,
																																 MarkerBlastHistogramAnnotation.DEFAULT_DESCRIPTION,
																																 null);
		blastAlignmentHistogram.parseAnnotation(vc, log);
		markerSeqAnnotation = new MarkerSeqAnnotation();
		markerSeqAnnotation.parseAnnotation(vc, log);
		markerEvalueHistogramAnnotation = new MarkerEvalueHistogramAnnotation(MarkerEvalueHistogramAnnotation.DEFAULT_NAME,
																																					MarkerEvalueHistogramAnnotation.DEFAULT_DESCRIPTION);
		markerEvalueHistogramAnnotation.parseAnnotation(vc, log);
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

	public static Map<String, MarkerBlastAnnotation> initForMarkers(final String[] markers) {
		Map<String, MarkerBlastAnnotation> blastAnnotations = Maps.newHashMapWithExpectedSize(markers.length);
		for (String marker : markers) {
			blastAnnotations.put(marker, new MarkerBlastAnnotation(marker));
		}
		return blastAnnotations;
	}

}
