/**
 *
 */
package org.genvisis.cnv.var;

import org.genvisis.common.Positions;
import org.genvisis.filesys.GeneTrack;

/**
 * @author Michael Vieths
 *
 *         Store the region, which can be parsed by Positions.parseUCSClocation(), and the label,
 *         which may or may not be set
 *
 */
public class Region {
	private String label = "";
	private int chr;
	private int start;
	private int stop;

	/**
	 * Set the region based on a tab-delimited string. The first column is the UCSC location, while
	 * the second (optional) token is a label.
	 *
	 * @param line
	 */
	public Region(String line) {
		this(line, null);
	}

	public Region(String line, GeneTrack track) {
		if (line == null) {
			throw new IllegalArgumentException("Region string is null");
		}

		String[] tokens = line.split("\t");
		if (tokens.length > 0 && !tokens[0].isEmpty()) {
			parseLoc(tokens[0], track);
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
		setLoc(location);
	}

	/**
	 * @return A parseable UCSC location
	 */
	public String getRegion() {
		return "chr" + chr + ":" + start + "-" + stop;
	}

	public int getChr() {
		return chr;
	}

	public int getStart() {
		return start;
	}

	public int getStop() {
		return stop;
	}

	/**
	 * @return A label for this region, defined in the region file
	 */
	public String getLabel() {
		return label;
	}

	private void setLoc(int[] location) {
		chr = location[0];
		start = location[1];
		stop = location[2];
	}

	private void parseLoc(String rawLoc, GeneTrack track) {
		String location = rawLoc.toLowerCase();
		int[] loc;

		if (!location.startsWith("chr")) {
			if (track == null) {
				throw new IllegalStateException("Cannot parse gene name: '"	+ location
																				+ "' since the gene track has either not been installed or did not load properly.");
			}
			loc = track.lookupPosition(location);
		} else {
			loc = Positions.parseUCSClocation(location);
		}
		if (loc[0] == -1) {
			throw new IllegalArgumentException("'"	+ location
																					+ "' is not a valid gene name or is not a valid UCSC location.");
		}
		setLoc(loc);
	}
}
