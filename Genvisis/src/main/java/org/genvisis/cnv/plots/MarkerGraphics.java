package org.genvisis.cnv.plots;

import java.awt.Color;
import java.awt.Graphics;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.genvisis.cnv.analysis.MarkerStats;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.gui.ColorSequence;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Numbers;
import org.genvisis.common.Positions;
import org.genvisis.stats.BinnedMovingAverage;


/**
 * A collection of {@link MarkerStats} data for use in visualization.
 */
public class MarkerGraphics {
	/**
	 * Smoothing strategy for moving averages - dictates the unit of the window distance.
	 */
	public enum Smoothing {
													KBASEPAIRS, MARKERS, PIXELS
	}

	public static final int DEFAULT_MOVING_AVG = 0;

	// Info about the current chromosome with loaded markers
	private int centromereStart = -1;
	private int centromereEnd = -1;
	private int loadedChr = -1;

	// Y-axis buffer to prevent graph from colliding with other components.
	private final int heightBuffer;


	// Sorted list of markers for this chromosome
	private List<MarkerCols> markers = new ArrayList<MarkerCols>();
	private String markerStatsFile;

	// List of components we can render (columns not in MarkerStats.ID_COLUMNS)
	private List<String> markerComponents = new ArrayList<String>();
	private Map<String, Color> markerColors = new HashMap<String, Color>();
	private Map<String, Integer> columnMap = new HashMap<String, Integer>();
	private int[][] centromereBoundaries;
	private final Logger log;

	/**
	 * Create a MarkerGraphics object that can be used to draw statistics from the specified file.
	 */
	public MarkerGraphics(String markerStatsFile, Project proj) {
		this(markerStatsFile, proj, 16);
	}

	public MarkerGraphics(String markerStatsFile, Project proj, int heightBuffer) {
		this.heightBuffer = heightBuffer;
		this.markerStatsFile = markerStatsFile;
		log = proj.getLog();
		String markerSetFilename = proj.PLINK_DIR_FILEROOTS.getValue() + "plink.bim";

		int genomeBuild = proj.GENOME_BUILD_VERSION.getValue().getBuildInt();
		centromereBoundaries = Positions.determineCentromereBoundariesFromMarkerSet(markerSetFilename,
		                                                                            genomeBuild, log);
		int count = 0;
		for (String col : Files.getHeaderOfFile(markerStatsFile, log)) {
			columnMap.put(col, count++);
			if (!MarkerStats.ID_COLUMNS.contains(col)) {
				markerComponents.add(col);
				markerColors.put(col, ColorSequence.get(markerColors.size()));
			}
		}
	}

	/**
	 * @return The chromosome index currently loaded into memory.
	 */
	public int getChr() {
		return loadedChr;
	}

	/**
	 * @return How many different lines can be drawn
	 */
	public int getNumComponents() {
		return markerComponents.size();
	}

	/**
	 * @return A list of the component names available for this MarkerGraphics
	 */
	public List<String> getComponents() {
		return Collections.unmodifiableList(markerComponents);
	}

	// TODO extract API for just looking at regions for genome start/end?
	/**
	 * Draw all components with no moving average
	 *
	 * @see {@link #draw(Graphics, int, int, int, int, int, int, Smoothing, List)}
	 */
	public void draw(Graphics g, int chr, int genomeStart, int genomeEnd, int pixelWidth,
	                 int height) {
		draw(g, chr, genomeStart, genomeEnd, pixelWidth, height, markerComponents);
	}

	/**
	 * Draw with specified moving average
	 *
	 * @see {@link #draw(Graphics, int, int, int, int, int, int, Smoothing, List)}
	 */
	public void draw(Graphics g, int chr, int genomeStart, int genomeEnd, int pixelWidth, int height,
	                 String... components) {
		draw(g, chr, genomeStart, genomeEnd, pixelWidth, height, DEFAULT_MOVING_AVG, components);
	}

	/**
	 * Convenience method for converting varargs to {@code List<String>}
	 *
	 * @see {@link #draw(Graphics, int, int, int, int, int, int, Smoothing, List)}
	 */
	public void draw(Graphics g, int chr, int genomeStart, int genomeEnd, int pixelWidth, int height,
	                 int movingAverage, String... components) {
		draw(g, chr, genomeStart, genomeEnd, pixelWidth, height, movingAverage,
		     Arrays.asList(components));
	}

	/**
	 * Use default moving average
	 *
	 * @see {@link #draw(Graphics, int, int, int, int, int, int, Smoothing, List)}
	 */
	public void draw(Graphics g, int chr, int genomeStart, int genomeEnd, int pixelWidth, int height,
	                 List<String> components) {
		draw(g, chr, genomeStart, genomeEnd, pixelWidth, height, DEFAULT_MOVING_AVG, components);
	}

	/**
	 * Use {@link Smoothing#MARKERS} strategy
	 *
	 * @see {@link #draw(Graphics, int, int, int, int, int, int, Smoothing, List)}
	 */
	public void draw(Graphics g, int chr, int genomeStart, int genomeEnd, int pixelWidth, int height,
	                 int movingAverage, List<String> components) {

		draw(g, chr, genomeStart, genomeEnd, pixelWidth, height, movingAverage, Smoothing.MARKERS,
		     components);
	}

	/**
	 * @param g AWT Graphics object to draw in
	 * @param chr Chromosome of interest
	 * @param genomeStart First genome position, in base pairs, that is in view
	 * @param genomeEnd Last genome position, in base pairs, that is in view
	 * @param pixelWidth Width of the UI component to fill
	 * @param height Height of the UI component to fill
	 * @param movingAverage How many positions to use for smoothing (half this value in both
	 *        directions)
	 * @param smoothing - Smoothing strategy to use when creating moving averages
	 * @param components List of one or more {@link MarkerStats} components to draw
	 */
	public void draw(Graphics g, int chr, int genomeStart, int genomeEnd, int pixelWidth, int height,
	                 int movingAverage, Smoothing smoothing, List<String> components) {
		if (this.loadedChr != chr) {
			load(chr);
		}
		int genomeWidth = genomeEnd - genomeStart;

		Color c = g.getColor();

		// TODO binary search for start point and stop loops when passed end
		for (String comp : components) {
			// Draw a single point for each stat, for each marker
			g.setColor(markerColors.get(comp));

			if (movingAverage <= 0) {
				for (int i = 0; i < markers.size(); i++) {
					MarkerCols marker = markers.get(i);
					if (marker.pos >= genomeStart && marker.pos <= genomeEnd) {
						int y = getScaledY(marker.get(comp), height, heightBuffer);
						int x = getX(marker.pos, pixelWidth, genomeStart, genomeWidth);
						g.fillOval(x, y, 4, 4);
					}
				}
			} else {
				// Smooth using a moving average
				BinnedMovingAverage<Double> bma = new BinnedMovingAverage<Double>(movingAverage);

				// If pixel smoothing, each bin covers 1 X position
				// If marker smoothing, each bin covers 1 marker
				// If KBP smoothing, each bin covers 1KB.

				for (int i = 0; i < markers.size(); i++) {
					MarkerCols marker = markers.get(i);
					int pos = marker.pos;

					if (pos >= centromereStart && pos <= centromereEnd) {
						// traversing the centromere
						bma.clear();
					} else {
						int bin = 0;
						// Identify which bin the current marker belongs in
						switch (smoothing) {
							case KBASEPAIRS:
								bin = pos / 1000;
								break;
							case MARKERS:
								bin = i;
								break;
							case PIXELS:
								bin = getX(pos, pixelWidth, genomeStart, genomeWidth);
								break;
						}
						// if adding the current value will add a new bin, we want to draw the current stats
						if (!bma.lastBin(bin)) {
							drawPoint(g, genomeStart, pixelWidth, height, smoothing, genomeWidth, bma, false);
						}
						Double val = marker.get(comp);
						// prune NaNs
						if (!val.isNaN()) {
							bma.add(val, bin);
						}
					}
				}
				// Draw the last point if needed
				drawPoint(g, genomeStart, pixelWidth, height, smoothing, genomeWidth, bma, true);
			}
		}
		// Restore original grahpics color
		g.setColor(c);

	}

	private void drawPoint(Graphics g, int genomeStart, int pixelWidth, int height,
	                      Smoothing smoothing, int genomeWidth, BinnedMovingAverage<Double> bma, boolean forceCurrentBin) {
		double v = forceCurrentBin ? bma.getValWithCurrentBin() : bma.getValue();
		if (v > -1) {
			int x = bma.mid();
			switch (smoothing) {
				case MARKERS:
					// X is a marker index.. we want to convert it to a position
					x = markers.get(x).pos;
					x = getX(x, pixelWidth, genomeStart, genomeWidth);
					break;
				case KBASEPAIRS:
					// X already was a position in KBP. Convert back to BP
					x *= 1000;
					x = getX(x, pixelWidth, genomeStart, genomeWidth);
					break;
				case PIXELS:
					// X is in pixels and thus already set
					break;
			}
			int y = getScaledY(v, height, heightBuffer);
			// TODO draw lines instead - just need to remember two points at a time and draw
			// between them. Maybe also with thickness ~ SD?
			g.fillOval(x, y, 4, 4);
		}
	}

	/**
	 * Load marker data for the given chromosome into memory
	 */
	private void load(int targetChr) {
		try {
			BufferedReader reader = Files.getAppropriateReader(markerStatsFile);
			// Read all markers from the target chromosome into memory

			// Skip header
			reader.readLine();

			while (reader.ready()) {
				String[] line = reader.readLine().split(MarkerStats.DELIM);
				int markerChr = Integer.parseInt(line[columnMap.get(MarkerStats.ID_CHR)]);
				// If this marker matches the target chromosome, cache it
				if (markerChr == targetChr) {
					MarkerCols m = new MarkerCols(markerComponents, columnMap, line);
					markers.add(m);
				}
			}

			Collections.sort(markers);
			loadedChr = targetChr;

			centromereStart = centromereBoundaries[loadedChr][0];
			centromereEnd = centromereBoundaries[loadedChr][1];

			reader.close();
		} catch (FileNotFoundException e) {
			log.reportException(e);
		} catch (IOException e) {
			log.reportException(e);
		}
	}

	/**
	 * Converts an X-axis position in base pairs to a pixel position
	 */
	private static int getX(int pos, int pixelWidth, int genomeStart, int genomeWidth) {
		return (int) (((double) (pos - genomeStart) / genomeWidth) * pixelWidth) + Trailer.WIDTH_BUFFER;
	}

		/**
	 * Scales y for rendering in the given height
	 * 
	 * @return The column value for the given component, scaled to the specified height. NB: smaller
	 *         component values will return a "larger" height, since it is a distance from the top (0)
	 *         position.
	 */
	private static int getScaledY(Double y, int height, int buffer) {
		int adjustedHeight = height - buffer;
		return buffer + (adjustedHeight - (int) (y * adjustedHeight));
	}

	/**
	 * Helper class for sorting marker stats by physical position. NB: component values are assumed to
	 * range from 0-1
	 */
	private static class MarkerCols implements Comparable<MarkerCols> {
		private String name;
		private int pos;
		private Map<String, Double> componentValues = new HashMap<String, Double>();

		public MarkerCols(List<String> components, Map<String, Integer> columnMap, String[] line) {
			name = line[columnMap.get(MarkerStats.ID_MARKER_NAME)];
			pos = Integer.parseInt(line[columnMap.get(MarkerStats.ID_POS)]);

			for (String comp : components) {
				String v = line[columnMap.get(comp)];
				double doubleV = Double.parseDouble(v);
				// All values should be in a range from 0-1. If a number is > 1, assume it's scaled 0-100
				// and adjust.
				if (doubleV > 1) {
					doubleV /= 100;
				}
				componentValues.put(comp, doubleV);
			}
		}

		/**
		 * @return Value for this marker of the specified component
		 */
		public Double get(String component) {
			return componentValues.get(component);
		}

		@Override
		public int compareTo(MarkerCols o) {
			int c = Numbers.compare(pos, o.pos);
			if (c == 0) {
				c = name.compareTo(o.name);
			}
			return c;
		}
	}
}
