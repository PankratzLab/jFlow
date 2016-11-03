package org.genvisis.cnv.annotation.markers;

import java.util.Hashtable;
import java.util.Map;

import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Logger;
import org.genvisis.filesys.Segment;

import com.google.common.collect.Maps;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;

public class MarkerGCAnnotation extends LocusAnnotation implements AnnotationParser {

	public enum GC_TYPE {
												MARKER_GC_CONTENT;
	}

	private boolean found;
	private static final String DEFAULT_VALUE = ".";

	public MarkerGCAnnotation(Builder builder, String locusName, Segment seg) {
		super(builder, locusName, seg);

	}

	public static AnnotationData getGCAnnotationDatas() {
		return new AnnotationData(VCFHeaderLineType.Float, null, 1,
															GC_TYPE.MARKER_GC_CONTENT.toString(),
															"The gc content of the marker, not including the interrogation position",
															DEFAULT_VALUE, DEFAULT_VALUE);
	}

	@Override
	public void parseAnnotation(VariantContext vc, Logger log) {
		for (int i = 0; i < getAnnotations().length; i++) {
			if (vc.hasAttribute(getAnnotations()[i].getName())) {
				getAnnotations()[i].setData(vc.getAttributeAsString(getAnnotations()[i].getName(),
																														getAnnotations()[i].getDefaultValue()));
			}
		}
	}

	@Override
	public boolean shouldAnnotateWith(VariantContext vc, Logger log) {
		return getLocusName().equals(vc.getID());
	}

	public static Map<String, MarkerGCAnnotation> initForMarkers(Project proj, String[] markers,
																										MarkerSet markerSet,
																										Map<String, Integer> indices) {
		if (markerSet == null) {
			markerSet = proj.getMarkerSet();
		}
		if (indices == null) {
			indices = proj.getMarkerIndices();
		}

		Map<String, MarkerGCAnnotation> markerGCAnnotations = Maps.newHashMapWithExpectedSize(markers.length);
		for (String marker : markers) {
			Builder builder = new Builder();
			builder.annotations(new AnnotationData[] {getGCAnnotationDatas()});
			int index = indices.get(marker);
			markerGCAnnotations.put(marker, new MarkerGCAnnotation(builder, marker,
																											new Segment(markerSet.getChrs()[index],
																																	markerSet.getPositions()[index],
																																	markerSet.getPositions()[index])));
		}
		return markerGCAnnotations;
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
