package cnv.plots;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;

import stats.Maths;

import common.*;

/**
 * Forest Panel
 *
 * @author Rohit Sinha
 */
public class ForestPanel extends AbstractPanel {

	public static final Color[] DEFAULT_COLORS = { new Color(33, 31, 53), // dark dark
			new Color(23, 58, 172), // dark blue
			new Color(201, 30, 10), // deep red
			new Color(140, 20, 180), // deep purple
			new Color(33, 87, 0), // dark green
			new Color(55, 129, 252), // light blue
			new Color(94, 88, 214), // light purple
			new Color(189, 243, 61), // light green
			new Color(217, 109, 194), // pink
			new Color(0, 0, 128), // ALL KINDS OF BLUES
			new Color(100, 149, 237), new Color(72, 61, 139), new Color(106, 90, 205), new Color(123, 104, 238), new Color(132, 112, 255), new Color(0, 0, 205), new Color(65, 105, 225), new Color(0, 0, 255), new Color(30, 144, 255), new Color(0, 191, 255), new Color(135, 206, 250), new Color(135, 206, 250), new Color(70, 130, 180), new Color(176, 196, 222), new Color(173, 216, 230), new Color(176, 224, 230), new Color(175, 238, 238), new Color(0, 206, 209), new Color(72, 209, 204), new Color(64, 224, 208), new Color(0, 255, 255), new Color(224, 255, 255),

	};
	protected ForestPlot forestPlot;
	private Logger log;
	private boolean swapAxes;
	private boolean rectangleGeneratable;
	DecimalFormat precision2Decimal;

	public ForestPanel(ForestPlot forestPlot, Logger log) {
		super();
		rectangleGeneratable = true;

		this.forestPlot = forestPlot;
		this.log = log;
		setColorScheme(DEFAULT_COLORS);
		precision2Decimal = new DecimalFormat("##.00");

		setZoomable(true, true);
		setColorScheme(DEFAULT_COLORS);
	}

	public boolean isRectangleGeneratable() {
		return rectangleGeneratable;
	}

	public void setRectangleGeneratable(boolean rectangleGeneratable) {
		this.rectangleGeneratable = rectangleGeneratable;
	}

	@Override
	public void generatePoints() {
		ArrayList<ForestTree> currentData = forestPlot.curTrees;
		PlotPoint[] tempPoints = new PlotPoint[currentData.size()];
		ArrayList<GenericLine> linesData = new ArrayList<GenericLine>();

		float xAxisValue, yAxisValue;
		lines = new GenericLine[0];

		for (int i = 0; i < currentData.size(); i++) {
			if(currentData.get(i).getBeta() != 0 && currentData.get(i).getStderr() != 0){
				xAxisValue = currentData.get(i).getConfInterval()[0];
				yAxisValue = (float) i + 1;
				PlotPoint leftEnd = new PlotPoint(currentData.get(i).getLabel(), currentData.get(i).getShape(), xAxisValue, yAxisValue, (byte) 5, (byte) 0, (byte) 0);

				xAxisValue = currentData.get(i).getConfInterval()[1];
				yAxisValue = (float) i + 1;
				PlotPoint rightEnd = new PlotPoint(currentData.get(i).getLabel(), currentData.get(i).getShape(), xAxisValue, yAxisValue, (byte) 5, (byte) 0, (byte) 0);

				linesData.add(new GenericLine(leftEnd, rightEnd, (byte) 1, (byte) 0, (byte) 0, false));
			}
			xAxisValue = currentData.get(i).getBeta();
			yAxisValue = (float) i + 1;
			tempPoints[i] = new PlotPoint(currentData.get(i).getLabel() + "|" + prepareRightMarkers(currentData.get(i)), currentData.get(i).getShape(), xAxisValue, yAxisValue, (byte) 3, (byte) 0, (byte) 0);
			tempPoints[i].setVisible(false);
		}
		points = tempPoints;
		lines = Array.concatAll(lines, linesData.toArray(new GenericLine[linesData.size()]));
	}

	private String prepareRightMarkers(ForestTree forestTree) {
		return precision2Decimal.format(forestTree.getBeta()) + " (" + precision2Decimal.format(forestTree.confInterval[0]) + " , " + precision2Decimal.format(forestTree.confInterval[1]) + " )";
	}

	private void generateRectangles() {
		ArrayList<ForestTree> currentData = forestPlot.curTrees;
		ArrayList<GenericRectangle> rectData = new ArrayList<GenericRectangle>();
		rectangles = new GenericRectangle[0];

		float xAxisValue, yAxisValue;

		float xAxisStep = (float) calcStepStep(plotXmax - plotXmin);
		float yAxisStep = (float) calcStepStep(plotYmax - plotYmin);

		for (int i = 0; i < currentData.size(); i++) {
			if(currentData.get(i).getBeta() != 0 && currentData.get(i).getStderr() != 0){
				xAxisValue = currentData.get(i).getBeta();
				yAxisValue = (float) i + 1;
				float scale = currentData.get(i).getzScore() / forestPlot.getMaxZScore() / 4;
				rectData.add(new GenericRectangle(xAxisValue - xAxisStep * scale, yAxisValue - yAxisStep * scale, xAxisValue + xAxisStep * scale, yAxisValue + yAxisStep * scale, (byte) 5, true, false, (byte) 0, (byte) 0));
			}
		}
		rectangles = Array.concatAll(rectangles, rectData.toArray(new GenericRectangle[rectData.size()]));
	}

	@Override
	public void highlightPoints() {
		byte defaultSize;
		defaultSize = 8;
		for (int i = 0; i < points.length; i++) {
			if (points[i].isHighlighted()) {
				points[i].setSize((byte) (defaultSize * 1.5));
			} else {
				points[i].setSize((byte) (defaultSize));
			}

		}
	}

	@Override
	public void assignAxisLabels() {
		displayXaxis = displayYaxis = true;
		xAxisLabel = forestPlot.plotLabel;
		yAxisLabel = " ";
	}

	public void setSwapAxes(boolean swapAxes) {
		this.swapAxes = swapAxes;
	}

	@Override
	public void drawAll(Graphics g, boolean base) {
		float minimumObservedRawX, maximumObservedRawX, minimumObservedRawY, maximumObservedRawY;
		double[] plotMinMaxStep; // needs to be double, else x <= plotXmax can be inexact and leave off the last tick mark
		int sigFigs;
		String str, pos;
		int xLook, yLook;
		BufferedImage yLabel;
		FontMetrics fontMetrics;
		Graphics gfx;
		Hashtable<String, Vector<PlotPoint>> layers;
		Vector<PlotPoint> layer;
		String trav;
		String[] keys;
		int[] order;
		int step;
		long time;
		ProgressBarDialog prog;// zx
		int rectangleXPixel, rectangleYPixel, rectangleWidthPixel, rectangleHeightPixel;

		setFinalImage(false);
		if (isRandomTest()) {
			points = new PlotPoint[1000000];
			for (int i = 0; i < points.length; i++) {
				points[i] = new PlotPoint("", (byte) 1, (float) Math.random(), (float) Math.random(), (byte) 5, (byte) 0, (byte) 0);
			}
		} else if (isPointsGeneratable()) {
			generatePoints();
		}
		highlightPoints();

		if (points.length == 0) {
			g.setColor(Color.WHITE);
			g.fillRect(0, 0, getWidth(), getHeight());
			if (getNullMessage() != null) {
				g.setColor(Color.BLACK);
				g.drawString(getNullMessage(), getWidth() / 2 - g.getFontMetrics(g.getFont()).stringWidth(getNullMessage()) / 2, getHeight() / 2);
			}
			setFinalImage(true);
			return;
		}

		setLookupResolution(DEFAULT_LOOKUP_RESOLUTION);
		assignAxisLabels();

		// Scan for rawX, rawY range of the data points
		minimumObservedRawX = Float.MAX_VALUE;
		maximumObservedRawX = Float.MIN_VALUE;
		minimumObservedRawY = Float.MAX_VALUE;
		maximumObservedRawY = Float.MIN_VALUE;
		for (int i = 0; i < points.length && isFlow(); i++) {
			if (points[i] != null) {
				minimumObservedRawX = Maths.min(minimumObservedRawX, points[i].getRawX());
				maximumObservedRawX = Maths.max(maximumObservedRawX, points[i].getRawX());
				minimumObservedRawY = Maths.min(minimumObservedRawY, points[i].getRawY());
				maximumObservedRawY = Maths.max(maximumObservedRawY, points[i].getRawY());
			}
		}

		for (int i = 0; lines != null && i < lines.length && isFlow(); i++) {
			if (lines[i] != null) {
				minimumObservedRawX = Maths.min(minimumObservedRawX, lines[i].getStartX());
				maximumObservedRawX = Maths.max(maximumObservedRawX, lines[i].getStartX());
				minimumObservedRawY = Maths.min(minimumObservedRawY, lines[i].getStartY());
				maximumObservedRawY = Maths.max(maximumObservedRawY, lines[i].getStartY());

				minimumObservedRawX = Maths.min(minimumObservedRawX, lines[i].getStopX());
				maximumObservedRawX = Maths.max(maximumObservedRawX, lines[i].getStopX());
				minimumObservedRawY = Maths.min(minimumObservedRawY, lines[i].getStopY());
				maximumObservedRawY = Maths.max(maximumObservedRawY, lines[i].getStopY());
			}
		}

		minimumObservedRawX = minimumObservedRawX == Float.MAX_VALUE ? 0 : minimumObservedRawX;
		maximumObservedRawX = maximumObservedRawX == Float.MIN_VALUE ? 1 : maximumObservedRawX;
		minimumObservedRawY = minimumObservedRawY == Float.MAX_VALUE ? 0 : minimumObservedRawY;
		maximumObservedRawY = maximumObservedRawY == Float.MIN_VALUE ? 1 : maximumObservedRawY;

		// otherwise step is off
		minimumObservedRawX = minimumObservedRawX > 0 ? 0 : minimumObservedRawX;
		minimumObservedRawY = minimumObservedRawY > 0 ? 0 : minimumObservedRawY;

		minimumObservedRawX = Float.isNaN(forcePlotXmin) ? minimumObservedRawX : forcePlotXmin;
		maximumObservedRawX = Float.isNaN(forcePlotXmax) ? (maximumObservedRawX + (maximumObservedRawX - minimumObservedRawX) * (float) 0.01) : forcePlotXmax;
		minimumObservedRawY = Float.isNaN(forcePlotYmin) ? minimumObservedRawY : forcePlotYmin;
		maximumObservedRawY = Float.isNaN(forcePlotYmax) ? (maximumObservedRawY + (maximumObservedRawY - minimumObservedRawY) * (float) 0.01) : forcePlotYmax;

		if (makeSymmetric) {
			maximumObservedRawX = Math.max(maximumObservedRawX, maximumObservedRawY);
			maximumObservedRawY = maximumObservedRawX;
			minimumObservedRawX = Math.min(minimumObservedRawX, minimumObservedRawY);
			minimumObservedRawY = minimumObservedRawX;
		}

		setNumberOfNaNSamples(0);
		if (base) {
			g.fillRect(0, 0, getWidth(), getHeight());
			g.setFont(new Font("Arial", 0, AXIS_FONT_SIZE));

			fontMetrics = g.getFontMetrics(g.getFont());
			missingWidth = fontMetrics.stringWidth("X");
			missingWidth = fontMetrics.stringWidth("X");

			// Calculate the plot area's range (X-axis, Y-axis)
			plotMinMaxStep = null;
			if (displayXaxis) {
				canvasSectionMinimumX = WIDTH_Y_AXIS;
				canvasSectionMaximumX = getWidth() - WIDTH_BUFFER;
				canvasSectionMinimumY = 0;
				canvasSectionMaximumY = HEIGHT_X_AXIS;
				plotMinMaxStep = getPlotMinMaxStep(minimumObservedRawX, maximumObservedRawX, g, true);
				plotXmin = plotMinMaxStep[0];
				plotXmax = plotMinMaxStep[1];

				sigFigs = getNumSigFig(plotMinMaxStep[2]);
				for (double x = plotMinMaxStep[3]; x <= plotXmax; x += plotMinMaxStep[2]) {
					if (x >= plotXmin || !truncate) {
						Grafik.drawThickLine(g, getXPixel(x), getHeight() - canvasSectionMaximumY, getXPixel(x), getHeight() - (canvasSectionMaximumY - TICK_LENGTH), TICK_THICKNESS, Color.BLACK);
						str = ext.formDeci(Math.abs(x) < DOUBLE_INACCURACY_HEDGE ? 0 : x, sigFigs, true);
						g.drawString(str, getXPixel(x) - str.length() * 8, getHeight() - (canvasSectionMaximumY - TICK_LENGTH - 30));
					}
				}
				Grafik.drawThickLine(g, canvasSectionMinimumX - (int) Math.ceil((double) AXIS_THICKNESS / 2.0), getHeight() - canvasSectionMaximumY, canvasSectionMaximumX + (int) Math.ceil((double) AXIS_THICKNESS / 2.0), getHeight() - canvasSectionMaximumY, AXIS_THICKNESS, Color.BLACK);
				g.drawString(xAxisLabel, (getWidth() - WIDTH_Y_AXIS) / 2 - fontMetrics.stringWidth(xAxisLabel) / 2 + WIDTH_Y_AXIS, getHeight() - 20);
			}

			if (displayYaxis) {
				canvasSectionMinimumX = 0;
				canvasSectionMaximumX = WIDTH_Y_AXIS;
				canvasSectionMinimumY = HEIGHT_X_AXIS;
				canvasSectionMaximumY = getHeight() - HEAD_BUFFER;
				if (!makeSymmetric || plotMinMaxStep == null) {
					plotMinMaxStep = getPlotMinMaxStep(minimumObservedRawY, maximumObservedRawY, g, false);
				}
				plotYmin = plotMinMaxStep[0];
				plotYmax = plotMinMaxStep[1];
				sigFigs = getNumSigFig(plotMinMaxStep[2]);
				for (double y = plotMinMaxStep[3]; y <= plotYmax; y += plotMinMaxStep[2]) {
					if ((y >= plotYmin && y == (int) y && y <= points.length && y > 0) || !truncate) {
						//Grafik.drawThickLine(g, canvasSectionMaximumX - TICK_LENGTH, getYPixel(y), canvasSectionMaximumX, getYPixel(y), TICK_THICKNESS, Color.BLACK);
						str = ext.formDeci(Math.abs(y) < DOUBLE_INACCURACY_HEDGE ? 0 : y, sigFigs, true);
						str = str.split("\\.")[0];
						String left = points[Integer.parseInt(str) - 1].getId().split("\\|")[0];
						g.drawString(left, canvasSectionMaximumX - TICK_LENGTH - left.length() * 15 - 5, getYPixel(y) + 9);
						String right = points[Integer.parseInt(str) - 1].getId().split("\\|")[1];
						g.drawString(right, getWidth() - TICK_LENGTH - right.length() * 15 - 5, getYPixel(y) + 9);
					}
				}
				//Grafik.drawThickLine(g, canvasSectionMaximumX, getYPixel(plotYmin), canvasSectionMaximumX, getYPixel(plotYmax) - (int) Math.ceil((double) TICK_THICKNESS / 2.0), AXIS_THICKNESS, Color.BLACK);

				yLabel = new BufferedImage(fontMetrics.stringWidth(yAxisLabel), 36, BufferedImage.TYPE_INT_RGB);
				gfx = yLabel.createGraphics();
				gfx.setFont(new Font("Arial", 0, 28));
				gfx.setColor(Color.WHITE);
				gfx.fillRect(0, 0, getWidth(), getHeight());
				gfx.setColor(Color.BLACK);
				gfx.drawString(yAxisLabel, 0, yLabel.getHeight() - 6);

				g.drawImage(Grafik.rotateImage(yLabel, true), 10, (getHeight() - HEIGHT_X_AXIS) / 2 - fontMetrics.stringWidth(yAxisLabel) / 2, this);
			}

			if (errorMessage != null) {
				g.drawString(errorMessage, (getWidth() - WIDTH_Y_AXIS) / 2 - fontMetrics.stringWidth(errorMessage) / 2 + WIDTH_Y_AXIS, (getHeight() - HEAD_BUFFER - HEIGHT_X_AXIS) / 2 - 20 + HEAD_BUFFER);
			}

		}

		// TODO outercoordinates
		canvasSectionMinimumX = WIDTH_Y_AXIS;
		canvasSectionMaximumX = getWidth() - WIDTH_BUFFER;
		canvasSectionMinimumY = HEIGHT_X_AXIS;
		canvasSectionMaximumY = getHeight() - HEAD_BUFFER;

		// Draw the lines
		for (int i = 0; lines != null && i < lines.length && isFlow(); i++) {
			if ((base && (getLayersInBase() == null || Array.indexOfByte(getLayersInBase(), lines[i].getLayer()) >= 0)) || (!base && Array.indexOfByte(getExtraLayersVisible(), lines[i].getLayer()) >= 0)) {
				Grafik.drawThickLine(g, getXPixel(lines[i].getStartX()), getYPixel(lines[i].getStartY()), getXPixel(lines[i].getStopX()), getYPixel(lines[i].getStopY()), (int) lines[i].getThickness(), colorScheme[lines[i].getColor()]);
			}
		}

		if (isRectangleGeneratable()) {
			generateRectangles();
		}

		// Draw the rectangles for clusterFilters
		for (int i = 0; rectangles != null && i < rectangles.length && isFlow(); i++) {
			if ((base && (getLayersInBase() == null || Array.indexOfByte(getLayersInBase(), rectangles[i].getLayer()) >= 0)) || (!base && Array.indexOfByte(getExtraLayersVisible(), rectangles[i].getLayer()) >= 0)) {
				rectangleXPixel = Math.min(getXPixel(rectangles[i].getStartXValue()), getXPixel(rectangles[i].getStopXValue()));
				rectangleYPixel = Math.min(getYPixel(rectangles[i].getStartYValue()), getYPixel(rectangles[i].getStopYValue()));
				rectangleWidthPixel = Math.abs(getXPixel(rectangles[i].getStartXValue()) - getXPixel(rectangles[i].getStopXValue()));
				rectangleHeightPixel = (Math.abs(getYPixel(rectangles[i].getStartYValue()) - getYPixel(rectangles[i].getStopYValue())));
				g.setColor(colorScheme[rectangles[i].getColor()]);

				if (rectangles[i].getFill()) {
					if (rectangles[i].getRoundedCorners()) {
						g.fillRoundRect(rectangleXPixel, rectangleYPixel, rectangleWidthPixel, rectangleHeightPixel, 2, 2);
					} else {
						g.fillRect(rectangleXPixel, rectangleYPixel, rectangleWidthPixel, rectangleHeightPixel);
					}
				} else {
					if (rectangles[i].getRoundedCorners()) {
						g.drawRoundRect(rectangleXPixel, rectangleYPixel, rectangleWidthPixel, rectangleHeightPixel, 2, 2);
					} else {
						drawRectThick(g, rectangleXPixel, rectangleYPixel, rectangleWidthPixel, rectangleHeightPixel, rectangles[i].getThickness());
					}
				}
			}
		}

		// Draw the rectangle outlined by dragging the mouse
		if (highlightRectangle != null) {
			rectangleXPixel = Math.min(getXPixel(highlightRectangle.getStartXValue()), getXPixel(highlightRectangle.getStopXValue()));
			rectangleYPixel = Math.min(getYPixel(highlightRectangle.getStartYValue()), getYPixel(highlightRectangle.getStopYValue()));
			rectangleWidthPixel = Math.abs(getXPixel(highlightRectangle.getStartXValue()) - getXPixel(highlightRectangle.getStopXValue()));
			rectangleHeightPixel = (Math.abs(getYPixel(highlightRectangle.getStartYValue()) - getYPixel(highlightRectangle.getStopYValue())));
			g.setColor(colorScheme[0]);
			drawRectThick(g, rectangleXPixel, rectangleYPixel, rectangleWidthPixel, rectangleHeightPixel, (byte) 1);
		}

		// Draw data points, also build the lookup matrix for nearby points.
		locLookup.clear();
		prog = null;
		time = new Date().getTime();
		step = Math.max((points.length) / 100, 1);
		layers = new Hashtable<String, Vector<PlotPoint>>();

		if (chartType == HEAT_MAP_TYPE) {
			drawHeatMap(g, null);
		} else if (chartType == SCATTER_PLOT_TYPE) {
			for (int i = 0; i < points.length && isFlow(); i++) {
				if (base && i % step == 0) {
					if (new Date().getTime() - time > 1000) {
						if (prog == null) {
							prog = new ProgressBarDialog("Generating image...", 0, points.length, getWidth(), getHeight(), 5000);// zx
						}
						prog.setProgress(i);// zx
					}
				}
				if (points[i] == null || points[i].getColor() == -1 || !points[i].isVisble()) {

				} else if (truncate && (points[i].getRawX() < plotXmin || points[i].getRawX() - plotXmax > plotXmax / 1000.0 || points[i].getRawY() < plotYmin || points[i].getRawY() > plotYmax)) {
				} else {
					trav = points[i].getLayer() + "";
					if (points[i].isHighlighted() || (base && (getLayersInBase() == null || Array.indexOfByte(getLayersInBase(), points[i].getLayer()) >= 0)) || (!base && Array.indexOfByte(getExtraLayersVisible(), points[i].getLayer()) >= 0)) {
						if (trav.equals("0")) {
							if (points[i].getType() != PlotPoint.NOT_A_NUMBER) {
								drawPoint(g, points[i]);
							} else if (base) {
								setNumberOfNaNSamples(getNumberOfNaNSamples() + 1);
							}
						} else {
							if (layers.containsKey(trav)) {
								layer = layers.get(trav);
							} else {
								layers.put(trav, layer = new Vector<PlotPoint>());
							}
							layer.add(points[i]);
						}
					}
					if (createLookup && points[i] != null) {
						xLook = (int) Math.floor(getXPixel(points[i].getRawX()) / getLookupResolution());
						yLook = (int) Math.floor(getYPixel(points[i].getRawY()) / getLookupResolution());
						for (int j = xLook - 1; j <= xLook + 1; j++) {
							for (int k = yLook - 1; k <= yLook + 1; k++) {
								pos = j + "x" + k;
								if (locLookup.containsKey(pos)) {
									locLookup.get(pos).add(i);
								} else {
									locLookup.put(pos, new IntVector(new int[] { i }));
								}
							}
						}
					}
				}
			}

			// Draw those points with layer>0.
			keys = HashVec.getKeys(layers);
			order = Sort.quicksort(Array.toIntArray(keys));
			for (int i = 0; i < keys.length && isFlow(); i++) {
				layer = layers.get(keys[order[i]]);
				for (int j = 0; j < layer.size(); j++) {
					if (layer.elementAt(j).getType() != PlotPoint.NOT_A_NUMBER) {
						drawPoint(g, layer.elementAt(j));
					} else {
						setNumberOfNaNSamples(getNumberOfNaNSamples() + 1);
					}
				}
			}
		} else {
			System.err.println("Error - invalid chart type: " + chartType);
		}

		if (getNumberOfNaNSamples() > 0) {
			g.drawString(PlotPoint.NAN_STR + " (n=" + getNumberOfNaNSamples() + ")", getXPixel(0) - nanWidth / 2, getYPixel(0) + 60 + points[0].getSize() / 2);
		}

		if (base && displayGrid) {
			for (double d = 0; d < 1.0; d += 0.1) {
				g.drawLine(getXPixel(d), getYPixel(0), getXPixel(d), getYPixel(canvasSectionMaximumY));
			}
			for (double d = -0.5; d < 0.5; d += 0.1) {
				g.drawLine(getXPixel(0), getYPixel(d), getXPixel(canvasSectionMaximumX), getYPixel(d));
			}
		}
		setFinalImage(true);

		if (base && prog != null) {
			prog.close();// zxu
		}
		refreshOtherComponents();
	}
}
