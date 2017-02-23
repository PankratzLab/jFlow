/**
 * 
 */
package org.genvisis.cnv.annotation.segments;

import java.io.Serializable;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Kitty
 *
 */
public class SegmentAnotation implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private Map<String, List<String>> attributes;



	/**
	 * Currently designed to store keys associated with {@link SegmentAnnotationKeys}
	 */
	public SegmentAnotation() {
		super();
		this.attributes = new HashMap<String, List<String>>();
	}



	/**
	 * @return the attribute map
	 */
	public Map<String, List<String>> getAttributes() {
		return attributes;
	}



}
