package cnv.annotation;

import htsjdk.variant.vcf.VCFHeaderLineType;

/**
 * @author lane0212 Used for forming the annotation header of an {@link AnnotationFile}
 */
public abstract class Annotation {

	private VCFHeaderLineType type;
	protected String name;
	private String description;
	private String defaultValue;

	public Annotation(VCFHeaderLineType type, String name, String description, String defaultValue) {
		super();
		this.type = type;
		this.name = name;
		this.description = description;
		this.defaultValue = defaultValue;
	}

	public VCFHeaderLineType getType() {
		return type;
	}

	public String getName() {
		return name;
	}

	public String getDescription() {
		return description;
	}

	public String getDefaultValue() {
		return defaultValue;
	}
	

}
