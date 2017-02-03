package org.genvisis.cnv.plots;

import java.awt.Color;
import java.awt.Graphics;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
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
import org.genvisis.stats.BinnedMovingStatistic;
import org.genvisis.stats.BinnedMovingStatistic.MovingStat;


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
	public List<RenderParams> getParams() {
		List<RenderParams> params = new ArrayList<RenderParams>();
		for (String comp : markerComponents) {
			params.add(new RenderParams(comp));
			
		}
		return params;
	}

	// TODO extract API for just looking at regions for genome start/end?

	/**
	 * @param g AWT Graphics object to draw in
	 * @param chr Chromosome of interest
	 * @param genomeStart First genome position, in base pairs, that is in view
	 * @param genomeEnd Last genome position, in base pairs, that is in view
	 * @param pixelWidth Width of the UI component to fill
	 * @param height Height of the UI component to fill
	 * @param params A list of {@link RenderParams} indicating component, moving average, smoothing, etc... to be drawn
	 */
	public void draw(Graphics g, int chr, int genomeStart, int genomeEnd, int pixelWidth, int height, List<RenderParams> params) {
		if (this.loadedChr != chr) {
			load(chr);
		}
		int genomeWidth = genomeEnd - genomeStart;

		Color c = g.getColor();

		// TODO binary search for start point and stop loops when passed end
		for (RenderParams p : params) {
			// Draw a single point for each stat, for each marker
			g.setColor(markerColors.get(p.getComponent()));

			if (p.getMovingWindow() <= 0) {
				for (int i = 0; i < markers.size(); i++) {
					MarkerCols marker = markers.get(i);
					if (marker.pos >= genomeStart && marker.pos <= genomeEnd) {
						int y = getScaledY(marker.get(p.getComponent()), height, heightBuffer);
						int x = getX(marker.pos, pixelWidth, genomeStart, genomeWidth);
						g.fillOval(x, y, 4, 4);
					}
				}
			} else {
				// Smooth using a moving average
				BinnedMovingStatistic<Double> bma = new BinnedMovingStatistic<Double>(p.getMovingWindow(), p.getMovingStat());

				// Whether the last marker was drawn
				boolean drewPoint = false;

				// If pixel smoothing, each bin covers 1 X position
				// If marker smoothing, each bin covers 1 marker
				// If KBP smoothing, each bin covers 1KB.
				for (int i = 0; i < markers.size(); i++) {
					MarkerCols marker = markers.get(i);
					int pos = marker.pos;
					drewPoint = false;

					if (pos >= centromereStart && pos <= centromereEnd) {
						// traversing the centromere
						bma.clear();
					} else {
						int bin = 0;
						// Identify which bin the current marker belongs in
						switch (p.getSmoothing()) {
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
						if (!bma.inBin(bin)) {
							drawPoint(g, genomeStart, pixelWidth, height, p.getSmoothing(), genomeWidth, bma);
							drewPoint = true;
						}
						Double val = marker.get(p.getComponent());
						// prune NaNs
						if (!val.isNaN()) {
							bma.add(val, bin);
						}
					}
				}
				// Draw the last point if needed
				if (!drewPoint) {
					drawPoint(g, genomeStart, pixelWidth, height, p.getSmoothing(), genomeWidth, bma);
				}
			}
		}
		// Restore original grahpics color
		g.setColor(c);

	}

	private void drawPoint(Graphics g, int genomeStart, int pixelWidth, int height,
	                      Smoothing smoothing, int genomeWidth, BinnedMovingStatistic<Double> bma) {
		double v = bma.getValue();
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

	/**
	 * Container class for parameters for {@link MarkerGraphics#draw(Graphics, int, int, int, int, int, List)}
	 */
	public static class RenderParams {
		private int movingWindow;
		private Smoothing smoothing;
		private MovingStat movingType;
		private String component;

		public RenderParams(String component) {
			this(0, Smoothing.MARKERS, MovingStat.MEAN, component);
		}

		public RenderParams(int movingWindow, Smoothing smoothing, MovingStat movingType, String component) {
			this.movingWindow = movingWindow;
			this.smoothing = smoothing;
			this.movingType = movingType;
			this.component = component;
		}

		public int getMovingWindow() {
			return movingWindow;
		}

		public void setMovingWindow(int movingWindow) {
			this.movingWindow = movingWindow;
		}

		public Smoothing getSmoothing() {
			return smoothing;
		}

		public void setSmoothing(Smoothing smoothing) {
			this.smoothing = smoothing;
		}

		public MovingStat getMovingStat() {
			return movingType;
		}

		public void setMovingType(MovingStat movingType) {
			this.movingType = movingType;
		}

		public String getComponent() {
			return component;
		}

		public void setComponent(String component) {
			this.component = component;
		}
	}
}
