package cnv.annotation;

import htsjdk.variant.vcf.VCFHeaderLineType;
import filesys.Segment;

public class MarkerGCAnnotation extends LocusAnnotation {

	private static final String DEFAULT_VALUE = ".";

	public MarkerGCAnnotation(Builder builder, String locusName, Segment seg) {
		super(builder, locusName, seg);

	}

	public static AnnotationData getGCAnnotationDatas() {
		return new AnnotationData(VCFHeaderLineType.Float, null, 1, "MARKER_GC_CONTENT", "The gc content of the marker, not including the interrogation position", DEFAULT_VALUE, DEFAULT_VALUE);
	}

}
