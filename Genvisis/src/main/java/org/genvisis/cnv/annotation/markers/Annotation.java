package org.genvisis.cnv.annotation.markers;

import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;

/**
 * @author lane0212 Used for forming the annotation header of an {@link AnnotationFile}
 */
/**
 * The Number entry is an Integer that describes the number of values that can be included with the
 * INFO field.<br>
 * For example, if the INFO field contains a single number, then this value should be 1.<br>
 * However, if the INFO field describes a pair of numbers, then this value should be 2 and so on.
 * <br>
 * If the number of possible values varies, is unknown, or is unbounded, then this value should be
 * '.'.<br>
 * Possible Types are: Integer, Float, Character, String and Flag. The 'Flag' type indicates that
 * the INFO field does not contain a Value entry, and hence the Number should be 0 in this case.<br>
 * The Description value must be surrounded by double-quotes.
 *
 */
public abstract class Annotation {

	private final VCFHeaderLineType type;
	protected String name;
	private final String description;
	private final String defaultValue;
	private final int number;
	private final VCFHeaderLineCount count;

	/**
	 * @param type {@link VCFHeaderLineType} , currently not used except by htsjdk error checks
	 * @param count {@link VCFHeaderLineCount}
	 * @param number number of entries, unbounded number can be specified with count
	 * @param name name of the entry
	 * @param description its description
	 * @param defaultValue
	 */
	public Annotation(VCFHeaderLineType type, VCFHeaderLineCount count, int number, String name,
										String description, String defaultValue) {
		super();
		this.type = type;
		this.count = count;
		this.number = number;
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

	public int getNumber() {
		return number;
	}

	public VCFHeaderLineCount getCount() {
		return count;
	}

}
