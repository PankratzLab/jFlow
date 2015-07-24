package cnv.annotation;

import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import common.Logger;
import common.ext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;
import filesys.Segment;

public class MarkerGCAnnotation extends LocusAnnotation implements AnnotationParser {

	public enum GC_TYPE {
		MARKER_GC_CONTENT;
	}

	private static final String DEFAULT_VALUE = ".";

	public MarkerGCAnnotation(Builder builder, String locusName, Segment seg) {
		super(builder, locusName, seg);

	}

	public static AnnotationData getGCAnnotationDatas() {
		return new AnnotationData(VCFHeaderLineType.Float, null, 1, GC_TYPE.MARKER_GC_CONTENT.toString(), "The gc content of the marker, not including the interrogation position", DEFAULT_VALUE, DEFAULT_VALUE);
	}

	@Override
	public void parseAnnotation(VariantContext vc, Logger log) {
		for (int i = 0; i < getAnnotations().length; i++) {
			if (vc.hasAttribute(getAnnotations()[i].getName())) {
				getAnnotations()[i].setData(vc.getAttributeAsString(getAnnotations()[i].getName(), getAnnotations()[i].getDefaultValue()));
			}
		}
	}

	@Override
	public boolean shouldAnnotateWith(VariantContext vc, Logger log) {
		return getLocusName().equals(vc.getID());
	}

	public static MarkerGCAnnotation[] initForMarkers(Project proj, String[] markers, MarkerSet markerSet, int[] markerIndicesInProject) {
		if (markerSet == null) {
			markerSet = proj.getMarkerSet();
		}
		if (markerIndicesInProject == null) {
			markerIndicesInProject = ext.indexLargeFactors(markers, proj.getMarkerNames(), true, proj.getLog(), true, false);
		}

		MarkerGCAnnotation[] markerGCAnnotations = new MarkerGCAnnotation[markers.length];
		for (int i = 0; i < markerGCAnnotations.length; i++) {
			Builder builder = new Builder();
			builder.annotations(new AnnotationData[] { getGCAnnotationDatas() });
			markerGCAnnotations[i] = new MarkerGCAnnotation(builder, markers[i], new Segment(markerSet.getChrs()[markerIndicesInProject[i]], markerSet.getPositions()[markerIndicesInProject[i]], markerSet.getPositions()[markerIndicesInProject[i]]));
		}
		return markerGCAnnotations;
	}

}
