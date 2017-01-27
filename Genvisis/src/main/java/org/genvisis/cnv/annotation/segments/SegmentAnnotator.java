/**
 * 
 */
package org.genvisis.cnv.annotation.segments;

import org.genvisis.filesys.Segment;

/**
 * @author Kitty
 *
 */
public interface SegmentAnnotator {


	/**
	 * @param segment the segment to be annotated
	 * @param segmentAnotation this method should add to this {@link SegmentAnotation}
	 * @return
	 */
	public void annotate(Segment segment, SegmentAnotation segmentAnotation);
		

}
