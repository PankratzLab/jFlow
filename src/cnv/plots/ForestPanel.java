package cnv.plots;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;

import stats.Maths;
import common.Array;
import common.Grafik;
import common.Logger;
import common.ext;

/**
 * Forest Panel
 */
public class ForestPanel extends AbstractPanel {
	static final long serialVersionUID = 1L;

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

	public static final String META_LABEL = "Overall";
	public static final Color META_COLOR = Color.BLACK;
	
	private ForestPlot forestPlot;
	private Logger log;
	private boolean rectangleGeneratable;
	private boolean sortedDisplay = false;
	private DecimalFormat precision2Decimal;

	public ForestPanel(ForestPlot forestPlot, Logger log) {
		super();
		rectangleGeneratable = true;
		
		this.forestPlot = forestPlot;
		this.log = log;
		setColorScheme(DEFAULT_COLORS);
		precision2Decimal = new DecimalFormat("#0.00");

		setZoomable(false, true);
//		setZoomable(true, true);
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
		if (forestPlot.getCurrentMetaStudy() == null) {
			ForestInput input;
			String errorMessage;
			
			if (forestPlot.getDataIndices().size() == 0) {
				errorMessage = "No data file selected, or no data found in input file.  Please select a data file.";
			} else {
				input = forestPlot.getDataIndices().get(forestPlot.getCurrentDataIndex());
				errorMessage = "Cannot generate points for marker "+input.marker+" because the data did not load; check to see if file \""+input.file+"\" actually exists and if the beta/stderr columns are named as expected";
			}
			setNullMessage(errorMessage);
			log.reportError("Error - "+errorMessage);
			
			points = new PlotPoint[0];
			setPointsGeneratable(false);
			setRectangleGeneratable(false);
			
			return;
		}
		ArrayList<StudyData> currentData = forestPlot.getCurrentMetaStudy().getStudies(sortedDisplay);
		PlotPoint[] tempPoints = new PlotPoint[currentData.size()];
		ArrayList<GenericLine> linesData = new ArrayList<GenericLine>();

		float xAxisValue, yAxisValue;
		lines = new GenericLine[0];
		
		for (int i = 0; i < currentData.size(); i++) {
//			if(currentData.get(i).getBeta() != 0.0 && currentData.get(i).getStderr() != 0.0){
				xAxisValue = currentData.get(i).getConfInterval()[0];
				yAxisValue = (float) i + 1;
				PlotPoint leftEnd = new PlotPoint(currentData.get(i).getLabel(), currentData.get(i).getShape(), xAxisValue, yAxisValue, (byte) 5, (byte) 0, (byte) 0);

				xAxisValue = currentData.get(i).getConfInterval()[1];
				yAxisValue = (float) i + 1;
				PlotPoint rightEnd = new PlotPoint(currentData.get(i).getLabel(), currentData.get(i).getShape(), xAxisValue, yAxisValue, (byte) 5, (byte) 0, (byte) 0);

				linesData.add(new GenericLine(leftEnd, rightEnd, (byte) 1, (byte) 0, (byte) 0, false));
				
				xAxisValue = currentData.get(i).getBeta();
				yAxisValue = (float) i + 1;
				tempPoints[i] = new PlotPoint(currentData.get(i).getLabel() + "|" + prepareRightMarkers(currentData.get(i)), currentData.get(i).getShape(), xAxisValue, yAxisValue, (byte) 3, (byte) 0, (byte) 0);
				tempPoints[i].setVisible(false);
//			}
		}
		points = tempPoints;
		lines = Array.concatAll(lines, linesData.toArray(new GenericLine[linesData.size()]));
	}

	private String prepareRightMarkers(float beta, float conf0, float conf1) {
		if (beta == 0.0f && conf0 == 0.0f && conf1 == 0.0f) {
			return "  ";
		}
		return String.format("%1$4s (%2$4s, %3$4s)", 
								precision2Decimal.format(beta), 
								precision2Decimal.format(conf0), 
								precision2Decimal.format(conf1));
	}
	
	private String prepareRightMarkers(StudyData forestTree) {
		return prepareRightMarkers(forestTree.getBeta(), 
									forestTree.getConfInterval()[0], 
									forestTree.getConfInterval()[1]);
	}

	private void generateRectangles(Graphics g) {
		ArrayList<StudyData> currentData = forestPlot.getCurrentMetaStudy().getStudies(sortedDisplay);
		ArrayList<GenericRectangle> rectData = new ArrayList<GenericRectangle>();
		rectangles = new GenericRectangle[0];

		float xAxisValue, yAxisValue;

		float xAxisStep = (float) calcStepStep(plotXmax - plotXmin);
		xAxisStep = Math.max(xAxisStep, 0.1f);
		float yAxisStep = (float) calcStepStep(plotYmax - plotYmin);
		yAxisStep = Math.max(yAxisStep, 0.1f);
		
		for (int i = 0; i < currentData.size(); i++) {
			if(currentData.get(i).getBeta() != 0 && currentData.get(i).getStderr() != 0){
				xAxisValue = currentData.get(i).getBeta();
				yAxisValue = (float) i + 1;
				// max scale = .25
//				float scale = currentData.get(i).getZScore() / forestPlot.getMaxZScore() / 4;
				
				float scale = (float) (Math.sqrt(currentData.get(i).getZScore()) / Math.sqrt(forestPlot.getMaxZScore()) / 4);
				float xDelta = xAxisStep * scale;
				float yDelta = yAxisStep * scale;
				xDelta = Math.max(xDelta, 0.05f);
				yDelta = Math.max(yDelta, 0.05f);
				rectData.add(new GenericRectangle(xAxisValue - xDelta, 
													yAxisValue - yDelta, 
													xAxisValue + xDelta, 
													yAxisValue + yDelta, 
													(byte) 5, true, false, (byte) 0, (byte) 0));
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
		xAxisLabel = forestPlot.getPlotLabel();
		yAxisLabel = " ";
	}
	
	@Override
	public double[] getPlotMinMaxStep(double min, double max, Graphics g, boolean xAxis) {
		double range, plotStep, stepStep, plotMin, plotMax;
		double zoomMin, zoomMax, dist;
		int numHashes, wid;
		FontMetrics fontMetrics = g.getFontMetrics(g.getFont());;
		int sf;

		range = max-min;

		plotStep = stepStep = calcStepStep(range);
		sf = ext.getNumSigFig(stepStep);

		if (xAxis) {
			wid = Math.max(fontMetrics.stringWidth(ext.formDeci(min, sf)), fontMetrics.stringWidth(ext.formDeci(max, sf)));
			numHashes = 12;
			while ((wid+30)*numHashes>canvasSectionMaximumX-canvasSectionMinimumX) {
				numHashes -= 1;
			}
		} else {
			numHashes = 10;
		}
		numHashes = Math.max(numHashes, 1);

		while (range / plotStep >= numHashes) {
			plotStep += stepStep;
		}
		plotMin = plotMax = 0;
		while (max - plotMax > DOUBLE_INACCURACY_HEDGE) {
			plotMax += plotStep;
		}
		while (min - plotMin < -1*DOUBLE_INACCURACY_HEDGE) { // double check this, untested
			plotMin -= plotStep;
		}

		dist = plotMax - plotMin;
		zoomMin = plotMin+zoomSubsets[xAxis?0:1][0]*dist;
		zoomMax = plotMax-(1-zoomSubsets[xAxis?0:1][1])*dist;
		
		range = zoomMax-zoomMin;
		plotStep = stepStep = calcStepStep(range);
		sf = ext.getNumSigFig(stepStep);

		if (xAxis) {
			fontMetrics = g.getFontMetrics(g.getFont());
			wid = Math.max(fontMetrics.stringWidth(ext.formDeci(min, sf)), fontMetrics.stringWidth(ext.formDeci(max, sf)));
			numHashes = 12;
			while ((wid+30)*numHashes>canvasSectionMaximumX-canvasSectionMinimumX) {
				numHashes -= 1;
			}
		} else {
			numHashes = (int) ((canvasSectionMaximumY - canvasSectionMinimumY) / (fontMetrics.getHeight()));
		}
		numHashes = Math.max(numHashes, 1);
		
		while (range / plotStep >= numHashes) {
			plotStep += stepStep;
		}

		return new double[] {zoomMin, zoomMax, Double.parseDouble(ext.formDeci(plotStep, sf)), Double.parseDouble(ext.formDeci(plotMin, sf))};

	}
	
	public static double calcStepStep(double range) {
		String[] line;
		
		try {
			line = new DecimalFormat("0.0E0").format(range).split("E");
			return Math.pow(10, Integer.parseInt(line[1])-1)*1;//(Double.parseDouble(line[0]));//>2.0?5:1);
		} catch (Exception e) {
			System.err.println("Error - could not parse stepStep from range '"+range+"'");
			return Double.NaN;
		}
	}
	
	@Override
	public void drawAll(Graphics g, boolean base) {
		float minimumObservedRawX, maximumObservedRawX, minimumObservedRawY, maximumObservedRawY;
		double[] plotMinMaxStep; // needs to be double, else x <= plotXmax can be inexact and leave off the last tick mark
		int sigFigs;
		String str;
		FontMetrics fontMetrics = g.getFontMetrics();
		int[] order = null;
		int rectangleXPixel, rectangleYPixel, rectangleWidthPixel, rectangleHeightPixel;
		
		setFinalImage(false);
//		if (isRandomTest()) {
//			points = new PlotPoint[1000000];
//			for (int i = 0; i < points.length; i++) {
//				points[i] = new PlotPoint("", (byte) 1, (float) Math.random(), (float) Math.random(), (byte) 5, (byte) 0, (byte) 0);
//			}
//		} else 
		if (isPointsGeneratable()) {
			generatePoints();
			if (isPointsGeneratable()) {
				TreeMap<Float, Integer> scoreIndexMap = new TreeMap<Float, Integer>();
				TreeSet<Integer> zeros = new TreeSet<Integer>();
				for (int i = 0; i < forestPlot.getCurrentMetaStudy().getStudies().size(); i++) {
					if (forestPlot.getCurrentMetaStudy().getStudies().get(i).getZScore() != 0.0) {
						scoreIndexMap.put((forestPlot.getCurrentMetaStudy().getStudies().get(i).getZScore() / forestPlot.getSumZScore()) * 100, i);
					} else {
						zeros.add(i);
					}
				}
				order = new int[forestPlot.getCurrentMetaStudy().getStudies().size()];
				int ind = scoreIndexMap.size() - 1;
				for (java.util.Map.Entry<Float, Integer> entry : scoreIndexMap.entrySet()) {
					order[ind] = entry.getValue().intValue();
					ind--;
				}
				Iterator<Integer> iter = zeros.iterator();
				for (int i = scoreIndexMap.size(); i < order.length; i++) {
					order[i] = iter.next().intValue();
				}
			}
		}
//		highlightPoints();

		if (!isPointsGeneratable() || points.length == 0) {
			g.setColor(Color.WHITE);
			g.fillRect(0, 0, getWidth(), getHeight());
			if (getNullMessage() != null) {
				g.setColor(Color.BLACK);
				String msg = getNullMessage();
				boolean skip = false;
				int width = g.getFontMetrics().stringWidth(msg);
				int index = msg.length() / 2;
				if (width >= getWidth() - 20) {
					// search for a space
					while (msg.charAt(index) != ' ' && index < msg.length()) {
						index++;
					}
					if (index == msg.length()) {
						index = msg.length() / 2;
						// search for a space other direction
						while (msg.charAt(index) != ' ' && index >= 0) {
							index--;
						}
					}
					if (index != 0) {
						skip = true;
					}
				}
				
				int x = getWidth() / 2 - g.getFontMetrics(g.getFont()).stringWidth(msg) / 2;
				int y = getHeight() / 2;
				
				if (!skip) {
					g.drawString(msg, x, y);
				} else {
					String s1, s2;
					int x1, x2, y1, y2;
					s1 = msg.substring(0, index);
					s2 = msg.substring(index);
					x1 = getWidth() / 2 - (g.getFontMetrics().stringWidth(s1) / 2);
					x2 = getWidth() / 2 - (g.getFontMetrics().stringWidth(s2) / 2);
					y1 = getHeight() / 2;
					y2 = getHeight() / 2 + g.getFontMetrics().getHeight() + 3;
					g.drawString(s1, x1, y1);
					g.drawString(s2, x2, y2);
				}
			}
			
			setFinalImage(true);
			return;
		}

		setLookupResolution(DEFAULT_LOOKUP_RESOLUTION);
		assignAxisLabels();

		// TODO simplify and condense:
		// TODO THESE NEVER CHANGE - DON'T RE-CALC EVERY TIME!
		
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
		int leftsize = 0;
		int rightsize = 0;
		if (base) {
			g.fillRect(0, 0, getWidth(), getHeight());
			g.setFont(new Font("Arial", 0, axisFontSize));

			fontMetrics = g.getFontMetrics(g.getFont());
			missingWidth = fontMetrics.stringWidth("X");
			missingWidth = fontMetrics.stringWidth("X");

			// Calculate the plot area's range (X-axis, Y-axis)
			plotMinMaxStep = null;
			
			float xRange = getXPixel(maximumObservedRawX) - getXPixel(minimumObservedRawX);
			float yRange = getYPixel(minimumObservedRawY) - getYPixel(maximumObservedRawY);
			
			leftsize = determineLongestLeft(g, getMarkerFontSize(g));
			rightsize = determineRightBorder(g, getMarkerFontSize(g));
//			System.out.println("size: " + getMarkerFontSize());

			if (displayXaxis) {
				canvasSectionMinimumX = WIDTH_BUFFER + leftsize;
				canvasSectionMaximumX = getWidth() - rightsize;
				canvasSectionMinimumY = 0;
				canvasSectionMaximumY = HEIGHT_X_AXIS;
				plotMinMaxStep = getPlotMinMaxStep(minimumObservedRawX, maximumObservedRawX, g, true);
				plotXmin = plotMinMaxStep[0];
				plotXmax = plotMinMaxStep[1];

				sigFigs = ext.getNumSigFig(plotMinMaxStep[2]);
				for (double x = plotMinMaxStep[3]; x <= plotXmax; x += plotMinMaxStep[2]) {
					if (x >= plotXmin || !truncate) {
						Grafik.drawThickLine(g, 
											 getXPixel(x), 
											 getHeight() - canvasSectionMaximumY, 
											 getXPixel(x), 
											 getHeight() - (canvasSectionMaximumY - TICK_LENGTH), 
											 TICK_THICKNESS, 
											 Color.BLACK);
						str = ext.formDeci(Math.abs(x) < DOUBLE_INACCURACY_HEDGE ? 0 : x, sigFigs, true);
						g.drawString(str, getXPixel(x) - str.length() * 8, getHeight() - (canvasSectionMaximumY - TICK_LENGTH - 30));
					}
				}
				Grafik.drawThickLine(g, canvasSectionMinimumX - (int) Math.ceil((double) AXIS_THICKNESS / 2.0), getHeight() - canvasSectionMaximumY, canvasSectionMaximumX + (int) Math.ceil((double) AXIS_THICKNESS / 2.0), getHeight() - canvasSectionMaximumY, AXIS_THICKNESS, Color.BLACK);
				g.drawString(xAxisLabel, 
						// graph-centered plot name:
							 canvasSectionMinimumX + ((canvasSectionMaximumX - canvasSectionMinimumX) / 2 - (fontMetrics.stringWidth(xAxisLabel) / 2)),
						// window-centered plot name:
//							 (getWidth() - WIDTH_Y_AXIS) / 2 - fontMetrics.stringWidth(xAxisLabel) / 2 + WIDTH_Y_AXIS, 
							 getHeight() - 20);
			}

			if (displayYaxis) {
				g.setFont(new Font("Arial", 0, (int)getMarkerFontSize(g)));
				fontMetrics = g.getFontMetrics();
				canvasSectionMinimumX = 0;
				canvasSectionMaximumX = WIDTH_BUFFER + rightsize;
				canvasSectionMinimumY = HEIGHT_X_AXIS + (fontMetrics.getHeight() * 2);
				canvasSectionMaximumY = getHeight() - (4 * HEAD_BUFFER);
				
//				System.out.println("points: " + (points.length + 1) + "y axis pixel:" + (canvasSectionMaximumY - canvasSectionMinimumY)/(points.length+1)*0.70);
				if (!makeSymmetric || plotMinMaxStep == null) {
					plotMinMaxStep = getPlotMinMaxStep(minimumObservedRawY, maximumObservedRawY, g, false);
//					System.out.println(Array.toStr(plotMinMaxStep));
				}
				plotYmin = plotMinMaxStep[0];
				plotYmax = plotMinMaxStep[1];
				
				sigFigs = ext.getNumSigFig(plotMinMaxStep[2]);
				double step = 1;//Math.max(1, Math.round(plotMinMaxStep[2] * 2) / 2.0f);
				for (double y = plotMinMaxStep[3]; y <= plotYmax; y += step) {
					if ((y >= plotYmin && y == (int) y && y <= points.length && y > 0) || !truncate) {
//						Grafik.drawThickLine(g, canvasSectionMaximumX - TICK_LENGTH, getYPixel(y), canvasSectionMaximumX, getYPixel(y), TICK_THICKNESS, Color.BLACK);
						str = ext.formDeci(Math.abs(y) < DOUBLE_INACCURACY_HEDGE ? 0 : y, sigFigs, true);
						str = str.split("\\.")[0];
						int index = Integer.parseInt(str) - 1;
						String left = points[index].getId().split("\\|")[0];
						g.drawString(left, WIDTH_BUFFER + leftsize - fontMetrics.stringWidth(left) - 15, getYPixel(y) + 7);
						String right = points[index].getId().split("\\|")[1];
						g.drawString(right, getWidth() - rightsize + 15, getYPixel(y) + 7);
					}
				}
				g.drawString(META_LABEL, WIDTH_BUFFER + leftsize - fontMetrics.stringWidth(META_LABEL) - 15, getHeight() - HEIGHT_X_AXIS - fontMetrics.getHeight() - 10);
				
				g.drawString(prepareRightMarkers(forestPlot.getCurrentMetaStudy().getMetaBeta(), forestPlot.getCurrentMetaStudy().getMetaConf()[0], forestPlot.getCurrentMetaStudy().getMetaConf()[1]), getWidth() - rightsize + 15, getHeight() - HEIGHT_X_AXIS - fontMetrics.getHeight() - 10);
//				Grafik.drawThickLine(g, canvasSectionMaximumX, getYPixel(plotYmin), canvasSectionMaximumX, getYPixel(plotYmax) - (int) Math.ceil((double) TICK_THICKNESS / 2.0), AXIS_THICKNESS, Color.BLACK);
//				g.setFont(new Font("Arial", 0, AXIS_FONT_SIZE));
//				yLabel = new BufferedImage(fontMetrics.stringWidth(yAxisLabel), 36, BufferedImage.TYPE_INT_RGB);
//				gfx = yLabel.createGraphics();
//				gfx.setFont(new Font("Arial", 0, 28));
//				gfx.setColor(Color.WHITE);
//				gfx.fillRect(0, 0, getWidth(), getHeight());
//				gfx.setColor(Color.BLACK);
//				gfx.drawString(yAxisLabel, 0, yLabel.getHeight() - 6);
//
//				g.drawImage(Grafik.rotateImage(yLabel, true), 10, (getHeight() - HEIGHT_X_AXIS) / 2 - fontMetrics.stringWidth(yAxisLabel) / 2, this);
			}
//
//			if (errorMessage != null) {
//				g.drawString(errorMessage, (getWidth() - WIDTH_Y_AXIS) / 2 - fontMetrics.stringWidth(errorMessage) / 2 + WIDTH_Y_AXIS, (getHeight() - HEAD_BUFFER - HEIGHT_X_AXIS) / 2 - 20 + HEAD_BUFFER);
//			}

		}

		// TODO outercoordinates		
		canvasSectionMinimumX = WIDTH_BUFFER + leftsize;
		canvasSectionMaximumX = getWidth() - rightsize;
		canvasSectionMinimumY = HEIGHT_X_AXIS + (fontMetrics.getHeight() * 2);
		canvasSectionMaximumY = getHeight() - (4 * HEAD_BUFFER);
		
		
		// Draw the lines
		for (int i = 0; lines != null && i < lines.length && isFlow(); i++) {
			if ((base && (getLayersInBase() == null || Array.indexOfByte(getLayersInBase(), lines[i].getLayer()) >= 0)) || (!base && Array.indexOfByte(getExtraLayersVisible(), lines[i].getLayer()) >= 0)) {
				Grafik.drawThickLine(g, 
										getXPixel(lines[i].getStartX()), 
										getYPixel(lines[i].getStartY()), 
										getXPixel(lines[i].getStopX()), 
										getYPixel(lines[i].getStopY()), 
										(int) lines[i].getThickness(), 
										colorScheme[lines[i].getColor()]);
			}
		}

		if (isRectangleGeneratable()) {
			generateRectangles(g);
		}

		// Draw the rectangles for clusterFilters
		int actWidth, actHeight, wDiff, hDiff;
		for (int i = 0; rectangles != null && i < rectangles.length && isFlow(); i++) {
			if ((base && (getLayersInBase() == null || Array.indexOfByte(getLayersInBase(), rectangles[i].getLayer()) >= 0)) || (!base && Array.indexOfByte(getExtraLayersVisible(), rectangles[i].getLayer()) >= 0)) {
				rectangleXPixel = Math.min(getXPixel(rectangles[i].getStartXValue()), getXPixel(rectangles[i].getStopXValue()));
				rectangleYPixel = Math.min(getYPixel(rectangles[i].getStartYValue()), getYPixel(rectangles[i].getStopYValue()));
				rectangleWidthPixel = Math.abs(getXPixel(rectangles[i].getStartXValue()) - getXPixel(rectangles[i].getStopXValue()));
				rectangleHeightPixel = Math.abs(getYPixel(rectangles[i].getStartYValue()) - getYPixel(rectangles[i].getStopYValue()));
				
				actWidth = Math.min(rectangleWidthPixel, rectangleHeightPixel);
				actHeight = Math.min(rectangleWidthPixel, rectangleHeightPixel);
				
				wDiff = rectangleWidthPixel - actWidth;
				hDiff = rectangleHeightPixel - actHeight;
				
				rectangleXPixel += wDiff / 2;
				rectangleYPixel += hDiff / 2;

				g.setColor(colorScheme[rectangles[i].getColor()]);

				if (rectangles[i].getFill()) {
					if (rectangles[i].getRoundedCorners()) {
						g.fillRoundRect(rectangleXPixel, rectangleYPixel, actWidth, actHeight, 2, 2);
					} else {
						g.fillRect(rectangleXPixel, rectangleYPixel, actWidth, actHeight);
					}
				} else {
					if (rectangles[i].getRoundedCorners()) {
						g.drawRoundRect(rectangleXPixel, rectangleYPixel, actWidth, actHeight, 2, 2);
					} else {
						drawRectThick(g, rectangleXPixel, rectangleYPixel, actWidth, actHeight, rectangles[i].getThickness());
					}
				}
			}
		}
		
		if(base) {
			g.setColor(Color.BLACK);
			int xL = getXPixel(forestPlot.getCurrentMetaStudy().getMetaBeta() - 1.96 * forestPlot.getCurrentMetaStudy().getMetaStderr());
			int xM = getXPixel(forestPlot.getCurrentMetaStudy().getMetaBeta());
			int xR = getXPixel(forestPlot.getCurrentMetaStudy().getMetaBeta() + 1.96 * forestPlot.getCurrentMetaStudy().getMetaStderr());
			
			int yM = getHeight() - HEIGHT_X_AXIS - fontMetrics.getHeight() - 15;
			int yU = yM - (fontMetrics.getHeight() / 2) - 1;
			int yD = yM + (fontMetrics.getHeight() / 2) + 1;
			
			Grafik.drawThickLine(g, xL, yM, xM, yU, 2, META_COLOR);
			Grafik.drawThickLine(g, xM, yU, xR, yM, 2, META_COLOR);
			Grafik.drawThickLine(g, xL, yM, xM, yD, 2, META_COLOR);
			Grafik.drawThickLine(g, xM, yD, xR, yM, 2, META_COLOR);
			
			int yMin = (4 * HEAD_BUFFER) - 5;
			int yMax = getHeight() - HEIGHT_X_AXIS;
			
			Grafik.drawThickLine(g, getXPixel(0.0), yMin, getXPixel(0.0), yMax, 3, Color.BLACK);
			
			int dashSize = 10;
			int dashSpacing = 5;
			int yStart;
			for (int i = 0; i < ((yMax - yMin) / (dashSize + dashSpacing)) + 1; i++) {
				yStart = yMin + (dashSize * i) + (dashSpacing * i);
				Grafik.drawThickLine(g, xM+1, yStart, xM+1, yStart + dashSize, 2, Color.GRAY);
			}

			g.setColor(Color.BLACK);
			g.setFont(new Font("Arial", Font.ITALIC, 16));
			if (forestPlot.getDataIndices().size() > 0) {
				String comm = forestPlot.getDataIndices().get(forestPlot.getCurrentDataIndex()).comment;
				if (!"".equals(comm)) {
					int w = fontMetrics.stringWidth(comm) / 2;
					g.drawString(comm, getWidth() / 2 - w, 3 * HEAD_BUFFER + 14);
				} 
			}
		}
		
		
		
//		// Draw the rectangle outlined by dragging the mouse
//		if (highlightRectangle != null) {
//			rectangleXPixel = Math.min(getXPixel(highlightRectangle.getStartXValue()), getXPixel(highlightRectangle.getStopXValue()));
//			rectangleYPixel = Math.min(getYPixel(highlightRectangle.getStartYValue()), getYPixel(highlightRectangle.getStopYValue()));
//			rectangleWidthPixel = Math.abs(getXPixel(highlightRectangle.getStartXValue()) - getXPixel(highlightRectangle.getStopXValue()));
//			rectangleHeightPixel = (Math.abs(getYPixel(highlightRectangle.getStartYValue()) - getYPixel(highlightRectangle.getStopYValue())));
//			g.setColor(colorScheme[0]);
//			drawRectThick(g, rectangleXPixel, rectangleYPixel, rectangleWidthPixel, rectangleHeightPixel, (byte) 1);
//		}

//		// Draw data points, also build the lookup matrix for nearby points.
//		locLookup.clear();
//		prog = null;
//		time = new Date().getTime();
//		step = Math.max((points.length) / 100, 1);
//		layers = new Hashtable<String, Vector<PlotPoint>>();
//
//		if (chartType == HEAT_MAP_TYPE) {
//			drawHeatMap(g, null);
//		} else if (chartType == SCATTER_PLOT_TYPE) {
//			for (int i = 0; i < points.length && isFlow(); i++) {
//				if (base && i % step == 0) {
//					if (new Date().getTime() - time > 1000) {
//						if (prog == null) {
//							prog = new ProgressBarDialog("Generating image...", 0, points.length, getWidth(), getHeight(), 5000);// zx
//						}
//						prog.setProgress(i);// zx
//					}
//				}
//				if (points[i] == null || points[i].getColor() == -1 || !points[i].isVisble()) {
//
//				} else if (truncate && (points[i].getRawX() < plotXmin || points[i].getRawX() - plotXmax > plotXmax / 1000.0 || points[i].getRawY() < plotYmin || points[i].getRawY() > plotYmax)) {
//				} else {
//					trav = points[i].getLayer() + "";
//					if (points[i].isHighlighted() || (base && (getLayersInBase() == null || Array.indexOfByte(getLayersInBase(), points[i].getLayer()) >= 0)) || (!base && Array.indexOfByte(getExtraLayersVisible(), points[i].getLayer()) >= 0)) {
//						if (trav.equals("0")) {
//							if (points[i].getType() != PlotPoint.NOT_A_NUMBER) {
//								drawPoint(g, points[i]);
//							} else if (base) {
//								setNumberOfNaNSamples(getNumberOfNaNSamples() + 1);
//							}
//						} else {
//							if (layers.containsKey(trav)) {
//								layer = layers.get(trav);
//							} else {
//								layers.put(trav, layer = new Vector<PlotPoint>());
//							}
//							layer.add(points[i]);
//						}
//					}
//					if (createLookup && points[i] != null) {
//						xLook = (int) Math.floor(getXPixel(points[i].getRawX()) / getLookupResolution());
//						yLook = (int) Math.floor(getYPixel(points[i].getRawY()) / getLookupResolution());
//						for (int j = xLook - 1; j <= xLook + 1; j++) {
//							for (int k = yLook - 1; k <= yLook + 1; k++) {
//								pos = j + "x" + k;
//								if (locLookup.containsKey(pos)) {
//									locLookup.get(pos).add(i);
//								} else {
//									locLookup.put(pos, new IntVector(new int[] { i }));
//								}
//							}
//						}
//					}
//				}
//			}
//
//			// Draw those points with layer>0.
//			keys = HashVec.getKeys(layers);
//			order = Sort.quicksort(Array.toIntArray(keys));
//			for (int i = 0; i < keys.length && isFlow(); i++) {
//				layer = layers.get(keys[order[i]]);
//				for (int j = 0; j < layer.size(); j++) {
//					if (layer.elementAt(j).getType() != PlotPoint.NOT_A_NUMBER) {
//						drawPoint(g, layer.elementAt(j));
//					} else {
//						setNumberOfNaNSamples(getNumberOfNaNSamples() + 1);
//					}
//				}
//			}
//		} else {
//			log.reportError("Error - invalid chart type: " + chartType);
//		}

//		if (getNumberOfNaNSamples() > 0) {
//			g.drawString(PlotPoint.NAN_STR + " (n=" + getNumberOfNaNSamples() + ")", getXPixel(0) - nanWidth / 2, getYPixel(0) + 60 + points[0].getSize() / 2);
//		}

//		if (base && displayGrid) {
//			for (double d = 0; d < 1.0; d += 0.1) {
//				g.drawLine(getXPixel(d), getYPixel(0), getXPixel(d), getYPixel(canvasSectionMaximumY));
//			}
//			for (double d = -0.5; d < 0.5; d += 0.1) {
//				g.drawLine(getXPixel(0), getYPixel(d), getXPixel(canvasSectionMaximumX), getYPixel(d));
//			}
//		}
		setFinalImage(true);

//		if (base && prog != null) {
//			prog.close();// zxu
//		}
		refreshOtherComponents();
	}

	private int determineRightBorder(Graphics g, double markerFontSize) {
		Font        defaultFont = new Font("Arial", 0, (int) markerFontSize);
		FontMetrics fmt =  g.getFontMetrics(defaultFont);
		int width = fmt.stringWidth("2222   2222  2222 2222");
		return (int) Math.floor(width);
	}

	private int determineLongestLeft(Graphics g, double markerFontSize) {
		Font        defaultFont = new Font("Arial", 0, (int) markerFontSize);
		FontMetrics fmt =  g.getFontMetrics(defaultFont);
		int width = fmt.stringWidth(forestPlot.getLongestStudyName());
		return (int) Math.floor(width);
	}

	private double getMarkerFontSize(Graphics g) {
//		double scale = 0.60;
//		double mkSz = ((getWidth() - WIDTH_BUFFER) - WIDTH_Y_AXIS) / (points.length * 2) * scale;
//		double mkSz = 20 - Math.atan(points.length - 15);
//		return mkSz;
		String largestLbl = forestPlot.getLongestStudyName();
		double yRange = getHeight() - (2 * HEAD_BUFFER) - (2 * HEIGHT_X_AXIS);
		double studies = points.length;
		double xMax = getWidth() / 4;
		double yMax = yRange / studies;
		
		int fontSz = 20;
		FontMetrics fm = g.getFontMetrics(new Font("Arial", 0, fontSz));
		while (fm.getHeight() > yMax || fm.stringWidth(largestLbl) > xMax) {
			fontSz--;
			fm = g.getFontMetrics(new Font("Arial", 0, fontSz));
		}
		
		return fontSz;
	}

	public void setSortedDisplay(boolean sorted) {
		this.sortedDisplay = sorted;
	}
}
