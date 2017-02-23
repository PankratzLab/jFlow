package org.genvisis.cnv.annotation.markers;

import org.genvisis.cnv.annotation.markers.AnnotationFileLoader.AnnotationQuery;
import org.genvisis.common.Logger;

import htsjdk.variant.variantcontext.VariantContext;

/**
 * @author lane0212
 *
 */
public interface AnnotationParser {

	/**
	 * So that promiscuous methods can parse from an {@link AnnotationQuery} and
	 * {@link VariantContext}
	 *
	 */
	public void parseAnnotation(VariantContext vc, Logger log);

	/**
	 * @param vc
	 * @return true if to use this {@link VariantContext } for parsing
	 */
	public boolean shouldAnnotateWith(VariantContext vc, Logger log);

	/**
	 * @param found can be used to store whether this annotation was found in an annotation file
	 */
	public void setFound(boolean found);

	/**
	 * @return whether the annotation was found or not
	 */
	public boolean isFound();
}
