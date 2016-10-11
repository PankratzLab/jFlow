package org.genvisis.cnv.annotation;

import java.util.Arrays;
import java.util.List;

import org.genvisis.common.Logger;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;

/**
 * @author lane0212 Stores the annotation value pair
 */
public class AnnotationData extends Annotation implements AnnotationParser {
	protected static final String DEFUALT_VALUE = ".";
	protected static final String DEFUALT_DELIMITER = ",";
	private String data;
	private boolean found;

	public AnnotationData(VCFHeaderLineType type, VCFHeaderLineCount count, int number, String name,
												String description, String data, String defaultValue) {
		super(type, count, number, name, description, defaultValue);
		this.data = data;
		found = false;
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

	public List<String> getDataAsList() {
		return Arrays.asList(data	.replaceAll("\\[", "").replaceAll("\\]", "")
															.split("\\s*" + DEFUALT_DELIMITER + "\\s*"));
	}

	@Override
	public void parseAnnotation(VariantContext vc, Logger log) {

	}

	@Override
	public boolean shouldAnnotateWith(VariantContext vc, Logger log) {
		return vc.getID().equals(name);
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
