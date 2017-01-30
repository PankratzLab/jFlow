/**
 * 
 */
package org.genvisis.cnv.annotation.segments;

/**
 * @author Kitty
 *
 */
public enum SegmentAnnotationKeys {

																		/**
																		 * Genes overlapped
																		 */
																		GENE("NA", false),
																		/**
																		 * GDI of genes overlapped
																		 */
																		GDI("NaN", true),
																		/**
																		 * Beast score of this region
																		 */
																		BEAST("NaN", true),
																		/**
																		 * Average mappability of this region
																		 */
																		MAPPABILITY("NaN", true);

	private String missingValue;
	private boolean isNumeric;

	private SegmentAnnotationKeys(String missingValue, boolean isNumeric) {
		this.missingValue = missingValue;
		this.isNumeric = isNumeric;
	}


	public boolean isNumeric() {
		return isNumeric;
	}

	public String getMissingValue() {
		return missingValue;
	}



}
