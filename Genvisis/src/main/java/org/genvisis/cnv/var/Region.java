/**
 * 
 */
package org.genvisis.cnv.var;

/**
 * @author Michael Vieths
 * 
 *         Store the region, which can be parsed by Positions.parseUCSClocation(), and the label, which may or may not be set
 * 
 */
public class Region {
	private String label = "";
	private String region = "";

	/**
	 * Set the region based on a tab-delimited string
	 * 
	 * @param line
	 */
	public Region(String line) {
		String[] tokens = line.split("\t");
		if (tokens.length > 0) {
			if (tokens[0] != "") {
				region = tokens[0];
			}
		}
		if (tokens.length > 1) {
			label = tokens[1];
		}
	}

	/**
	 * Set the region based on a UCSC-compatible location array
	 * 
	 * @param location
	 */
	public Region(int[] location) {
		region = "chr" + location[0] + ":" + location[1] + "-" + location[2];
	}

	/**
	 * @return A parseable UCSC location
	 */
	public String getRegion() {
		return region;
	}

	/**
	 * @return A label for this region, defined in the region file
	 */
	public String getLabel() {
		return label;
	}
}
