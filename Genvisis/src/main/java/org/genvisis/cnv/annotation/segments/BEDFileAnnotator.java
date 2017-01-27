/**
 * 
 */
package org.genvisis.cnv.annotation.segments;

import org.genvisis.seq.manage.BEDFileReader;

/**
 * Abstract class for annotating using content loaded from a bed file
 *
 *
 */

public abstract class BEDFileAnnotator extends BEDFileReader implements SegmentAnnotator {

	/**
	 * @param file the bed file to annotate with
	 * @param requireIndex
	 */
	public BEDFileAnnotator(String file) {
		super(file, true);
	}



}
