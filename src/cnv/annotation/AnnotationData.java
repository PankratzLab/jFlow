package cnv.annotation;

import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeaderLineType;

/**
 * @author lane0212 Stores the annotation value pair
 */
public class AnnotationData extends Annotation {
	private String data;

	public AnnotationData(VCFHeaderLineType type, String name, String description, String data) {
		super(type, name, description);
		this.data = data;
	}

	public String getData() {
		return data;
	}

	public void setData(String data) {
		this.data = data;
	}

	public void addAnnotation(VariantContextBuilder vBuilder) {
		vBuilder.attribute(name, data);
	}
}
