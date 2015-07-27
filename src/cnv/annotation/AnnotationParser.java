package cnv.annotation;

import cnv.annotation.AnnotationFileLoader.AnnotationQuery;
import common.Logger;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * @author lane0212
 *
 */
public interface AnnotationParser {

	/**
	 * So that promiscuous methods can parse from an {@link AnnotationQuery} and {@link VariantContext}
	 * 
	 */
	public void parseAnnotation(VariantContext vc, Logger log);

	/**
	 * @param vc
	 * @return true if to use this {@link VariantContext } for parsing
	 */
	public boolean shouldAnnotateWith(VariantContext vc, Logger log);

	public void setFound(boolean found);

	public boolean isFound();
}
