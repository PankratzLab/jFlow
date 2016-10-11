package org.genvisis.cnv.annotation;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

/**
 * @author lane0212 Used to read and write custom genvisis analysis parameters (like blast results)
 *         in {@link AnnotationFileWriter } and {@link AnnotationFileLoader}
 */
public interface AnalysisParams {

	/**
	 * @return a {@link VCFHeaderLine } to be added to a {@link VCFHeader}
	 */
	public VCFHeaderLine developHeaderLine();

	/**
	 * @return the key for the header
	 */
	public String getKey();

	/**
	 * @param vcfHeaderLine parse this header, ensure proper key match
	 */
	public void parseHeaderLine(VCFHeaderLine vcfHeaderLine);
}
