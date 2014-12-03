package cnv.plots;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;

import javax.imageio.ImageIO;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.Timer;

import mining.Distance;

import common.Array;
import common.Grafik;
import common.HashVec;
import common.IntVector;
import common.ProgressBarDialog;
import common.Sort;
import common.ext;
import stats.Maths;

public abstract class AbstractPanel extends JPanel implements MouseListener, MouseMotionListener, MouseWheelListener, ComponentListener, ActionListener {
	public static final long serialVersionUID = 1L;

	public static final int HEAD_BUFFER = 25;
	public static final int HEIGHT_X_AXIS = 105;
	public static final int WIDTH_Y_AXIS = 140;
	public static final int WIDTH_BUFFER = 50;
	public static final int AXIS_THICKNESS = 4;
	public static final int TICK_THICKNESS = 3;
	public static final int TICK_LENGTH = 15;
	public static final int DEFAULT_LOOKUP_RESOLUTION = 20;
	public static final int AXIS_FONT_SIZE = 28;
	public static final double DOUBLE_INACCURACY_HEDGE = 0.00001;
	public static final double MINIMUM_ZOOM_PROPORTION_WINDOW = 0.0001;
	public static final float DEFAULT_MOUSE_WHEEL_MULTIPLIER = 0.5f;
	public static final int DEFAULT_PLOTPOINTSET_SIZE = 1000000;
	public static final int SIZE = 12;
	public static final double HIGHLIGHT_DISTANCE = 20;//= Math.sqrt(SIZE*SIZE/2);
	public static final int DELAY = 0;	//A control variable to reduce the repaint() operations during component resizing;
	public static final int SCATTER_PLOT_TYPE = 1; 
	public static final int HEAT_MAP_TYPE = 2;

	public static final int DEFAULT_TYPE = SCATTER_PLOT_TYPE;
	
	protected Color[] colorScheme;
	protected int canvasSectionMinimumX;
	protected int canvasSectionMaximumX;
	protected int canvasSectionMinimumY;
	protected int canvasSectionMaximumY;
	protected double plotXmax, plotYmax;
	protected double plotXmin, plotYmin;
	protected BufferedImage image;
	protected String prevPos = "";
	protected Hashtable<String,IntVector> locLookup;
	protected PlotPoint[] points;						// make private when worked out
	protected GenericLine[] lines;
	protected GenericRectangle[] rectangles;
	protected GenericRectangle highlightRectangle;
	protected String xAxisLabel;
	protected String yAxisLabel;
	protected boolean displayXaxis;
	protected boolean displayYaxis;
	protected boolean displayGrid;
	protected int missingWidth;
	protected int nanWidth;
	protected int axisFontSize = AXIS_FONT_SIZE;
	protected float forcePlotXmax, forcePlotYmax;
	protected float forcePlotXmin, forcePlotYmin;
	protected boolean createLookup;
	protected boolean invertX;
	protected boolean invertY;
	protected boolean makeSymmetric;
	protected String errorMessage;
	protected float mouseWheelMultiplier;
	protected boolean zoomable;
	protected boolean swapable;		//4/27/2012
	protected boolean invertable;	//4/27/2012
	protected boolean truncate;
	protected float[][] zoomSubsets;
	protected IntVector prox;
	protected int chartType;

	private boolean inDrag;
	private int startX, startY;
	private int plotPointSetSize;
	private int totalNumPlotPointSets;
	private int currentPlotPointSet;
	private int lastIndexInPlotPointSet;
	private int currentIndexInPlotPointSet;
	private int lookupResolution;
	private boolean flow;			//A control variable. If resizing is not yet done, don't start generatePoints() or drawAll();
	private boolean finalImage;		//A control variable. If drawAll() is not yet done, don't start paintComponent();
	private byte[] layersInBase;
	private byte[] extraLayersVisible;
	private boolean pointsGeneratable;
	private Timer waitingTimer;		//A control variable to reduce the repaint() operations during component resizing;
	private String nullMessage;
	private boolean randomTest;
	private int numberOfNaNSamples;
	
	public AbstractPanel() {
		canvasSectionMinimumX = 0;
		canvasSectionMaximumX = getWidth();
		canvasSectionMinimumY = 0;
		canvasSectionMaximumY = getWidth();
		displayXaxis = true;
		displayYaxis = true;
		displayGrid = false;
		createLookup = true;
		missingWidth = -1;
		nanWidth = -1;
		forcePlotXmax = forcePlotYmax = forcePlotXmin = forcePlotYmin = Float.NaN;
		mouseWheelMultiplier = DEFAULT_MOUSE_WHEEL_MULTIPLIER;
		zoomable = false;
		swapable = true;
		resetZoomProportions();
		plotPointSetSize = DEFAULT_PLOTPOINTSET_SIZE;
		points = new PlotPoint[plotPointSetSize];
		totalNumPlotPointSets = 1;
		currentPlotPointSet = 0;
		lastIndexInPlotPointSet = -1;
		currentIndexInPlotPointSet = -1;
		randomTest = false;
		setChartType(DEFAULT_TYPE);
		
		layersInBase = null;
		extraLayersVisible = null;

		image = null;
		locLookup = new Hashtable<String,IntVector>();
		finalImage=false; // if true then double renders
		flow=true;
		pointsGeneratable = true;
		
		colorScheme = new Color[] {Color.BLACK, Color.GRAY};
		addMouseListener(this);
		addMouseMotionListener(this);
		addMouseWheelListener(this);
		addComponentListener(this);
	}
	
	public void resetCurrentIndexInPlotPointSet() {
		currentIndexInPlotPointSet = -1;
	}
	
	public boolean morePlotPointsExist() {
		return currentPlotPointSet<totalNumPlotPointSets && currentIndexInPlotPointSet<=lastIndexInPlotPointSet;
	}

	public void setNullMessage(String str) {
		nullMessage = str;
	}

	public void setColorScheme(Color[] scheme) {
		colorScheme = scheme;
	}
	
	public void setZoomable(boolean zoomable, boolean truncate) {
		this.zoomable = zoomable;
		this.truncate = truncate;
	}

	public void paintAgain() {
		image = null;
		repaint();
	}

	public boolean isRandomTest() {
		return randomTest;
	}

	public String getNullMessage() {
		return nullMessage;
	}

	public boolean isFlow() {
		return flow;
	}

	public int getNumberOfNaNSamples() {
		return numberOfNaNSamples;
	}

	public void setNumberOfNaNSamples(int numberOfNaNSamples) {
		this.numberOfNaNSamples = numberOfNaNSamples;
	}

	public int getLookupResolution() {
		return lookupResolution;
	}

	public byte[] getLayersInBase() {
		return layersInBase;
	}

	public byte[] getExtraLayersVisible() {
		return extraLayersVisible;
	}

	public void setLayersInBase(byte[] layers) {
		layersInBase = layers;
	}

	public void setExtraLayersVisible(byte[] layers) {
		extraLayersVisible = layers;
	}

	public void paintComponent(Graphics g) {
		if (getFinalImage()&&image==null) {
			createImage();
			// TODO extra paint appears to be unnecessary
//			repaint();
		}
		g.drawImage(image, 0, 0, this);
		if (extraLayersVisible != null && extraLayersVisible.length > 0) {
			drawAll(g, false);
		}
//		g.setColor(Color.WHITE);
//		g.fillRect(0, 0, getWidth(), getHeight());
//		if (finalImage&&image!=null) {
//			g.drawImage(image, 0, 0, this);
//		}

	}

	public void screenCapture(String filename) {
		try {
			ImageIO.write(image, "png", new File(filename));
		} catch (IOException ie) {
			JOptionPane.showMessageDialog(null, "Error while trying to save the plot");
		}
	}
	
	public void setForcePlotXmin(float forcePlotXmin) {
		this.forcePlotXmin = forcePlotXmin;
	}

	public void setForcePlotXmax(float forcePlotXmax) {
		this.forcePlotXmax = forcePlotXmax;
	}

	public void setForcePlotYmin(float forcePlotYmin) {
		this.forcePlotYmin = forcePlotYmin;
	}

	public void setForcePlotYmax(float forcePlotYmax) {
		this.forcePlotYmax = forcePlotYmax;
	}

	public abstract void generatePoints();

	public abstract void highlightPoints();

	public abstract void assignAxisLabels();

	public void setAxisFontSize(int sz) {
		if (sz > 0) {
			this.axisFontSize = sz;
		}
	}
	
	public void setChartType(int chartType) {
		this.chartType = chartType;
	}
	
	public void setXinversion(boolean b) {
		invertX = b;
	}

	public void setYinversion(boolean b) {
		invertY = b;
	}
	
	public void toggleXinversion() {
		invertX = !invertX;
	}

	public void toggleYinversion() {
		invertY = !invertY;
	}	

	public void setSymmetricAxes(boolean b) {
		makeSymmetric = b;
	}

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
		ProgressBarDialog prog;//zx
    	int rectangleXPixel, rectangleYPixel, rectangleWidthPixel, rectangleHeightPixel;
//    	int index;
    	
		// Set control variables; Generate data for the plot;  set Lookup Resolution; Prepare AxisLabels.
		setFinalImage(false);
		if (randomTest) {
			points = new PlotPoint[1000000];
			for (int i = 0; i < points.length; i++) {
				points[i] = new PlotPoint("", (byte)1, (float)Math.random(), (float)Math.random(), (byte)5, (byte)0, (byte)0);
			}
		} else if (pointsGeneratable) {
			generatePoints();
//			pointsGeneratable = false;
		}
		highlightPoints();
		
//		System.out.println("#points= "+points.length+" flow="+flow+"; base="+base);
		
		if (points.length==0) {
			g.setColor(Color.WHITE);
			g.fillRect(0, 0, getWidth(), getHeight());
			if (nullMessage != null) {
				g.setColor(Color.BLACK);
				g.drawString(nullMessage, getWidth()/2-g.getFontMetrics(g.getFont()).stringWidth(nullMessage)/2, getHeight()/2);
			}
//			System.err.println("Error: no data. The cnv.plots.AbstractPanel.points is null.");
			setFinalImage(true);
			return;
		}

		setLookupResolution(DEFAULT_LOOKUP_RESOLUTION);
		assignAxisLabels();

		//zx 4/30/2012 swap X Y
		if (invertable) {
			for (int i=0; i<points.length; i++) {
//				float temp = points[i].getRawX();
//				points[i].setRawX() = points[i].getRawY();
//				points[i].getRawY() = temp;
			}
		}

		// Scan for rawX, rawY range of the data points
		minimumObservedRawX = Float.MAX_VALUE;
		maximumObservedRawX = Float.MIN_VALUE;
		minimumObservedRawY = Float.MAX_VALUE;
		maximumObservedRawY = Float.MIN_VALUE;
		for (int i = 0; i<points.length&&flow; i++) {
//		for (int i = 0; i<points.length; i++) {
			if (points[i] != null) {
				minimumObservedRawX = Maths.min(minimumObservedRawX, points[i].getRawX());
				maximumObservedRawX = Maths.max(maximumObservedRawX, points[i].getRawX());
				minimumObservedRawY = Maths.min(minimumObservedRawY, points[i].getRawY());
				maximumObservedRawY = Maths.max(maximumObservedRawY, points[i].getRawY());
			}
        }

		for (int i = 0; lines != null && i<lines.length&&flow; i++) {
//		for (int i = 0; lines != null && i<lines.length; i++) {
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
		
		minimumObservedRawX = minimumObservedRawX==Float.MAX_VALUE?0:minimumObservedRawX;
		maximumObservedRawX = maximumObservedRawX==Float.MIN_VALUE?1:maximumObservedRawX;
		minimumObservedRawY = minimumObservedRawY==Float.MAX_VALUE?0:minimumObservedRawY;
		maximumObservedRawY = maximumObservedRawY==Float.MIN_VALUE?1:maximumObservedRawY;

//		otherwise step is off
		minimumObservedRawX = minimumObservedRawX>0?0:minimumObservedRawX;
		minimumObservedRawY = minimumObservedRawY>0?0:minimumObservedRawY;
		
		minimumObservedRawX = Float.isNaN(forcePlotXmin)?minimumObservedRawX:forcePlotXmin;
//		minimumObservedRawX = Float.isNaN(forcePlotXmin)?minimumObservedRawX-(maximumObservedRawX-minimumObservedRawX)*(float)0.01:forcePlotXmin;
//		maximumObservedRawX = Float.isNaN(forcePlotXmax)?maximumObservedRawX:forcePlotXmax;
		maximumObservedRawX = Float.isNaN(forcePlotXmax)?(maximumObservedRawX+(maximumObservedRawX-minimumObservedRawX)*(float)0.01):forcePlotXmax;
		minimumObservedRawY = Float.isNaN(forcePlotYmin)?minimumObservedRawY:forcePlotYmin;
//		minimumObservedRawY = Float.isNaN(forcePlotYmin)?(minimumObservedRawY-(maximumObservedRawY-minimumObservedRawY)*(float)0.01):forcePlotYmin;
//		maximumObservedRawY = Float.isNaN(forcePlotYmax)?maximumObservedRawY:forcePlotYmax;
		maximumObservedRawY =  Float.isNaN(forcePlotYmax)?(maximumObservedRawY+(maximumObservedRawY-minimumObservedRawY)*(float)0.01):forcePlotYmax;
		
		if (makeSymmetric) {
			maximumObservedRawX = Math.max(maximumObservedRawX, maximumObservedRawY);
			maximumObservedRawY = maximumObservedRawX;
			minimumObservedRawX = Math.min(minimumObservedRawX, minimumObservedRawY);
			minimumObservedRawY = minimumObservedRawX;
		}
		
		numberOfNaNSamples = 0;
		if (base) {
			//g.setColor(Color.WHITE);
			g.fillRect(0, 0, getWidth(), getHeight());
			g.setFont(new Font("Arial", 0, axisFontSize));
	//		System.out.println("getWidth: "+getWidth()+"\t getHeight: "+getHeight());
			
			fontMetrics = g.getFontMetrics(g.getFont());
			missingWidth = fontMetrics.stringWidth("X");
			missingWidth = fontMetrics.stringWidth("X");

			// Calculate the plot area's range (X-axis, Y-axis)
			plotMinMaxStep = null;
			if (displayXaxis) {
				canvasSectionMinimumX = WIDTH_Y_AXIS;
				canvasSectionMaximumX = getWidth()-WIDTH_BUFFER;
				canvasSectionMinimumY = 0;
				canvasSectionMaximumY = HEIGHT_X_AXIS;
				plotMinMaxStep = getPlotMinMaxStep(minimumObservedRawX, maximumObservedRawX, g, true);
				plotXmin = plotMinMaxStep[0];
				plotXmax = plotMinMaxStep[1];

				sigFigs = getNumSigFig(plotMinMaxStep[2]);
				for (double x = plotMinMaxStep[3]; x<=plotXmax; x += plotMinMaxStep[2]) {
					if (x >= plotXmin || !truncate) {
						Grafik.drawThickLine(g, getXPixel(x), getHeight()-canvasSectionMaximumY, getXPixel(x), getHeight()-(canvasSectionMaximumY-TICK_LENGTH), TICK_THICKNESS, Color.BLACK);
						str = ext.formDeci(Math.abs(x) < DOUBLE_INACCURACY_HEDGE ? 0 : x, sigFigs, true);
						g.drawString(str, getXPixel(x)-str.length()*8, getHeight()-(canvasSectionMaximumY-TICK_LENGTH-30));
					}
				}
				Grafik.drawThickLine(g, canvasSectionMinimumX-(int)Math.ceil((double)AXIS_THICKNESS/2.0), getHeight()-canvasSectionMaximumY, canvasSectionMaximumX+(int)Math.ceil((double)AXIS_THICKNESS/2.0), getHeight()-canvasSectionMaximumY, AXIS_THICKNESS, Color.BLACK);
				g.drawString(xAxisLabel, (getWidth()-WIDTH_Y_AXIS)/2-fontMetrics.stringWidth(xAxisLabel)/2+WIDTH_Y_AXIS, getHeight()-20);
			}
			
			if (displayYaxis) {
				canvasSectionMinimumX = 0;
				canvasSectionMaximumX = WIDTH_Y_AXIS;
				canvasSectionMinimumY = HEIGHT_X_AXIS;
				canvasSectionMaximumY = getHeight()-HEAD_BUFFER;
				if (!makeSymmetric || plotMinMaxStep == null) {
					plotMinMaxStep = getPlotMinMaxStep(minimumObservedRawY, maximumObservedRawY, g, false);
				}
				plotYmin = plotMinMaxStep[0];
				plotYmax = plotMinMaxStep[1];
				sigFigs = getNumSigFig(plotMinMaxStep[2]);
				for (double y = plotMinMaxStep[3]; y<=plotYmax; y += plotMinMaxStep[2]) {
					if (y >= plotYmin || !truncate) {
						Grafik.drawThickLine(g, canvasSectionMaximumX-TICK_LENGTH, getYPixel(y), canvasSectionMaximumX, getYPixel(y), TICK_THICKNESS, Color.BLACK);
						str = ext.formDeci(Math.abs(y) < DOUBLE_INACCURACY_HEDGE ? 0 : y, sigFigs, true);
						Font prevFont = g.getFont();
						if (str.length() == 5) {
							g.setFont(prevFont.deriveFont(prevFont.getSize2D() - 5));
							g.drawString(str, canvasSectionMaximumX - TICK_LENGTH - str.length() * 15 + 17, getYPixel(y) + 8);
							g.setFont(prevFont);
						} else if (str.length() >= 6) {
							g.setFont(prevFont.deriveFont(prevFont.getSize2D() - 10));
							g.drawString(str, canvasSectionMaximumX - TICK_LENGTH - str.length() * 15 + 35, getYPixel(y) + 6);
							g.setFont(prevFont);
						} else {
							g.drawString(str, canvasSectionMaximumX - TICK_LENGTH - str.length() * 15 + 5, getYPixel(y) + 9);
						}
					}
				}
				Grafik.drawThickLine(g, canvasSectionMaximumX, getYPixel(plotYmin), canvasSectionMaximumX, getYPixel(plotYmax)-(int)Math.ceil((double)TICK_THICKNESS/2.0), AXIS_THICKNESS, Color.BLACK);

				yLabel = new BufferedImage(fontMetrics.stringWidth(yAxisLabel), 36, BufferedImage.TYPE_INT_RGB);
				gfx = yLabel.createGraphics();
				gfx.setFont(new Font("Arial", 0, axisFontSize));
				gfx.setColor(Color.WHITE);
				gfx.fillRect(0, 0, getWidth(), getHeight());
				gfx.setColor(Color.BLACK);
				gfx.drawString(yAxisLabel, 0, yLabel.getHeight()-6);

				g.drawImage(Grafik.rotateImage(yLabel, true), 10, (getHeight()-HEIGHT_X_AXIS)/2-fontMetrics.stringWidth(yAxisLabel)/2, this);
			}

			/*
			if (swapable) {
				// these images are here to mask the flicker that occurs due to the slight delay in when the plot is painted and when the buttons are painted
				BufferedImage img;
				File file = new File("/workspace/Genvisis/images/flip_and_invert/flip_10p.jpg");
				new javax.swing.JLabel();
				try {
					img = ImageIO.read(file);
					g.drawImage(img, 70, getHeight()-75, null);
				} catch (IOException e) {
					e.printStackTrace();
				}
				file = new File("/workspace/Genvisis/images/flip_and_invert/right_10.gif");
				try {
					img = ImageIO.read(file);
					g.drawImage(img, 70, getHeight()-35, null);
				} catch (IOException e) {
					e.printStackTrace();
				}
				file = new File("/workspace/Genvisis/images/flip_and_invert/up_10.gif");
				try {
					img = ImageIO.read(file);
					g.drawImage(img, 55, getHeight()-75, null);
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			*/


			if (errorMessage != null) {
				g.drawString(errorMessage, (getWidth()-WIDTH_Y_AXIS)/2-fontMetrics.stringWidth(errorMessage)/2+WIDTH_Y_AXIS, (getHeight()-HEAD_BUFFER-HEIGHT_X_AXIS)/2-20+HEAD_BUFFER);
			}
		
		}

		
		//TODO outercoordinates
		canvasSectionMinimumX = WIDTH_Y_AXIS;
		canvasSectionMaximumX = getWidth()-WIDTH_BUFFER;
		canvasSectionMinimumY = HEIGHT_X_AXIS;
		canvasSectionMaximumY = getHeight()-HEAD_BUFFER;
//		System.err.println("("+canvasSectionMinimumX+"-"+canvasSectionMaximumX+","+canvasSectionMinimumY+"-"+canvasSectionMaximumY+")");

		// Draw the lines
		for (int i = 0; lines!=null && i<lines.length && flow; i++) {
			if ((base && (layersInBase == null || Array.indexOfByte(layersInBase, lines[i].getLayer()) >= 0)) || (!base && Array.indexOfByte(extraLayersVisible, lines[i].getLayer()) >= 0)) {
				Grafik.drawThickLine(g, getXPixel(lines[i].getStartX()), getYPixel(lines[i].getStartY()), getXPixel(lines[i].getStopX()), getYPixel(lines[i].getStopY()), (int)lines[i].getThickness(), colorScheme[lines[i].getColor()]);
			}
        }

		// Draw the rectangles for clusterFilters
		for (int i = 0; rectangles!=null && i<rectangles.length && flow; i++) {
			if ((base && (layersInBase == null || Array.indexOfByte(layersInBase, rectangles[i].getLayer()) >= 0)) || (!base && Array.indexOfByte(extraLayersVisible, rectangles[i].getLayer()) >= 0)) {
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
//						g.drawRect(rectangleXPixel, rectangleYPixel, rectangleWidthPixel, rectangleHeightPixel);
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
//	    	g.drawRect(rectangleXPixel, rectangleYPixel, rectangleWidthPixel, rectangleHeightPixel);
	    	drawRectThick(g, rectangleXPixel, rectangleYPixel, rectangleWidthPixel, rectangleHeightPixel, (byte) 1);
		}

//		// Draw progress bar
//		if (base) {
//			time = new Date().getTime();
//			prog = new ProgressBarDialog("Generating image...", 0, points.length, getWidth(), getHeight(), 5000);//zx
////			System.out.println("points.length: "+(points.length)+"\3*points.length: "+3*(points.length));
//		} else {
//			prog = null;
//		}

		// Draw data points, also build the lookup matrix for nearby points.
		locLookup.clear();
		prog = null;
		time = new Date().getTime();
		step = Math.max((points.length)/100, 1);
		layers = new Hashtable<String,Vector<PlotPoint>>();

		if (chartType == HEAT_MAP_TYPE) {
			drawHeatMap(g, null);
		} else if (chartType == SCATTER_PLOT_TYPE) {
			for (int i = 0; i<points.length && flow; i++) {
				if (base && i%step==0){
					if (new Date().getTime() - time > 1000) {
						if (prog == null) {
							prog = new ProgressBarDialog("Generating image...", 0, points.length, getWidth(), getHeight(), 5000);//zx
						}
						prog.setProgress(i);
					}
				}
				if (points[i] == null || points[i].getColor() == -1 || !points[i].isVisble()) {
					
				} else if (truncate && (points[i].getRawX() < plotXmin || points[i].getRawX() - plotXmax > plotXmax/1000.0 || points[i].getRawY() < plotYmin || points[i].getRawY() > plotYmax)) {
//					System.err.println("error: data point ("+points[i].getRawX()+","+points[i].getRawY()+") is outside of plot range.");
				} else {
					trav = points[i].getLayer()+"";
					if (points[i].isHighlighted() || (base && (layersInBase == null || Array.indexOfByte(layersInBase, points[i].getLayer()) >= 0)) || (!base && Array.indexOfByte(extraLayersVisible, points[i].getLayer()) >= 0)) {
						if (trav.equals("0")) {
							if (points[i].getType()!=PlotPoint.NOT_A_NUMBER) {
								drawPoint(g, points[i]);
							} else if (base) {
								numberOfNaNSamples++;
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
						xLook = (int)Math.floor(getXPixel(points[i].getRawX())/lookupResolution);
						yLook = (int)Math.floor(getYPixel(points[i].getRawY())/lookupResolution);
						for (int j = xLook-1; j<=xLook+1; j++) {
							for (int k = yLook-1; k<=yLook+1; k++) {
								pos = j+"x"+k;
								if (locLookup.containsKey(pos)) {
									locLookup.get(pos).add(i);
								} else {
									locLookup.put(pos, new IntVector(new int[] {i}));
								}
							}
						}
					}
				}
			}
			
			// Draw those points with layer>0.
			keys = HashVec.getKeys(layers);
			order = Sort.quicksort(Array.toIntArray(keys));
			for (int i = 0; i<keys.length&&flow; i++) {
				layer = layers.get(keys[order[i]]);
				for (int j = 0; j<layer.size(); j++) {
					if (layer.elementAt(j).getType()!=PlotPoint.NOT_A_NUMBER) {
						drawPoint(g, layer.elementAt(j));
					} else {
						numberOfNaNSamples++;
					}
	            }
	        }
		} else {
			System.err.println("Error - invalid chart type: "+chartType);
		}


		if (numberOfNaNSamples > 0) {
			g.drawString(PlotPoint.NAN_STR+" (n="+numberOfNaNSamples+")", getXPixel(0)-nanWidth/2, getYPixel(0)+60+points[0].getSize()/2);
		}
		
		if (base && displayGrid) {
			for (double d = 0; d < 1.0; d+=0.1) {
				g.drawLine(getXPixel(d), getYPixel(0), getXPixel(d), getYPixel(canvasSectionMaximumY));
			}
			for (double d = -0.5; d < 0.5; d+=0.1) {
				g.drawLine(getXPixel(0), getYPixel(d), getXPixel(canvasSectionMaximumX), getYPixel(d));
			}
		}
		setFinalImage(true);
		
		if (base && prog != null) {
			prog.close();
		}
		//System.out.println("Paint time: "+ext.getTimeElapsed(time));
		
		//test out
		/*
		BufferedImage temp = new BufferedImage(250, 250, BufferedImage.TYPE_INT_RGB);
		Graphics temp1 = temp.createGraphics();
//		temp1.setFont(new Font("Arial", 0, 28));
		temp1.setColor(Color.WHITE);
		temp1.fillRect(0, 0, getWidth(), getHeight());
		temp1.drawArc(0, 0, 135, 235, 25, 15);
		temp1.drawLine(40, 40, 160, 160);
		temp1.drawRect(100, 100, 50, 50);
//		temp1.setColor(Color.BLACK);
		temp1.drawString("This is temp", 0, temp.getHeight()-6);
//		g.drawImage(Grafik.rotateImage(yLabel, true), 10, (getHeight()-HEIGHT_X_AXIS)/2-fontMetrics.stringWidth(yAxisLabel)/2, this);
		temp1.dispose();
		g.drawImage(temp, 150, 100, null);
		*/
		
		refreshOtherComponents();
	}
	
	public void refreshOtherComponents() {
	}

	public void drawHeatMap(Graphics g, byte[] clusters) {
		int nRows, nColumns;
		int[][] gridIntensities;
		int[][][] gridColors;

		nRows = getHeight();
		nColumns = getWidth();
		gridIntensities = getGridIntensityForHeapMapGrid(nRows, nColumns, 1, 1, 3);
		gridColors = getColorFromIntensityForHeapMapGrid(gridIntensities);

		for (int i = 0; i < nColumns; i++) {
			for (int j = 0; j < nRows; j++) {
				if (gridIntensities[i][j] != 0) {
					g.setColor(new Color(gridColors[i][j][0], gridColors[i][j][1], gridColors[i][j][2]));
//					g.fillRect(i + canvasSectionMinimumX, j - canvasSectionMinimumY, 1, 1);
					g.fillOval(i + canvasSectionMinimumX, j - canvasSectionMinimumY, 2, 2);
				}
			}
		}
	}


	public int[][] getGridIntensityForHeapMapGrid(int nRows, int nColumns, int cellWidth, int cellHeight, int neighbor) {
		int xPixel, yPixel;
		int[][] intensities;
		boolean zoomedIn;
		int[] origin;

		origin = new int[] {0,0};

		zoomedIn = (Math.abs(getXValueFromXPixel(5) - getXValueFromXPixel(0)) < 0.002) || (Math.abs(getYValueFromYPixel(5) - getYValueFromYPixel(0)) < 0.002);
		
		intensities = new int[nColumns][nRows];
		for (int i = 0; i < points.length; i++) {
			if (points[i] != null) {
				xPixel = getXPixel(points[i].getRawX()) - canvasSectionMinimumX;
				yPixel = getYPixel(points[i].getRawY()) + canvasSectionMinimumY;

				for (int j = -neighbor; j <= neighbor; j++) {
					for (int k = -neighbor; k <= neighbor; k++) {
						if ((xPixel + j) >= 0 && (xPixel + j) < nColumns && (yPixel + k) >= 0 && (yPixel + k) < nRows) { //  && Distance.euclidean(new int[] {Math.abs(j), Math.abs(k)}, origin) < neighbor*(neighbor-1)
							if (zoomedIn) {
								intensities[xPixel + j][yPixel + k]++;
							} else {
								intensities[xPixel + j][yPixel + k] += neighbor*neighbor - Distance.euclidean(new int[] {Math.abs(j), Math.abs(k)}, origin);
							}
						}
					}
				}
			}
		}

		return intensities;
	}

	public int[][][] getColorFromIntensityForHeapMapGrid(int[][] intensities) {
		int[][][] color;
		int max;

		max = 0;
		for (int i = 0; i < intensities.length; i++) {
			for (int j = 0; j < intensities[i].length; j++) {
				if (max < intensities[i][j]) {
					max = intensities[i][j];
				}
			}
		}

		color = new int[intensities.length][intensities[0].length][3];
		for (int i = 0; i < intensities.length; i++) {
			for (int j = 0; j < intensities[i].length; j++) {
				if (intensities[i][j] != 0) {
					color[i][j] = Grafik.getHeatmapColor((double)intensities[i][j] / (double)max);
				}
			}
		}

		return color;
	}

	public void mouseClicked(MouseEvent e) {}

	public void mouseEntered(MouseEvent e) {}

	public void mouseExited(MouseEvent e) {}

//	public void mouseMoved(MouseEvent e) {}
	
	public void mouseMoved(MouseEvent event) {
		Graphics g = getGraphics();
		IntVector indicesOfNearbyPoints;
		String pos;
		int x, y, dataPointIndex;
		byte size;

		if (getFinalImage()) {
			x = event.getX();
			y = event.getY();

			canvasSectionMinimumX = WIDTH_Y_AXIS;
			canvasSectionMaximumX = getWidth() - WIDTH_BUFFER;
			canvasSectionMinimumY = HEIGHT_X_AXIS;
			canvasSectionMaximumY = getHeight() - HEAD_BUFFER;
			pos = (int)Math.floor(x / DEFAULT_LOOKUP_RESOLUTION) + "x" + (int)Math.floor(y / DEFAULT_LOOKUP_RESOLUTION);
			if (!pos.equals(prevPos)) {
				repaint();
			}

			indicesOfNearbyPoints = lookupNearbyPoints(x, y, pos);
			prox = new IntVector();

			size = SIZE * 2;
			g.setColor(Color.GRAY);
			for (int i = 0; indicesOfNearbyPoints!=null && i<indicesOfNearbyPoints.size(); i++) {
				dataPointIndex = indicesOfNearbyPoints.elementAt(i);
				if (Distance.euclidean(new int[] {x, y}, new int[] {getXPixel(points[dataPointIndex].getRawX()), getYPixel(points[dataPointIndex].getRawY())}) < HIGHLIGHT_DISTANCE) {
					prox.add(dataPointIndex);
//					g.setColor(Color.YELLOW);
//					g.fillOval(getX(points[dataPointIndex].getRawX()) - size/2, getY(points[dataPointIndex].getRawY()) - size/2, size, size);
					g.drawOval(getXPixel(points[dataPointIndex].getRawX()) - size/2, getYPixel(points[dataPointIndex].getRawY()) - size/2, size, size);

				}
			}

			prevPos = pos;
			
			//TODO 		indicesOfNearbySamples = lookupNearbyPoints(x, y, pos);
		}
	}

	public void mousePressed(MouseEvent e) {
		startX = e.getPoint().x;
		startY = e.getPoint().y;
		inDrag = true;
	}

	public void mouseReleased(MouseEvent e) {
		inDrag = false;
	}

	public void mouseDragged(MouseEvent e) {
		int curX, curY;
		double distance;
		
		System.err.println("AbstractPanel mouseDragged has been called");

		curX = e.getPoint().x;
		curY = e.getPoint().y;
		
		for (int i = 0; i < 2; i++) {
			if (i==0) {
				distance = (startX-curX) * (zoomSubsets[0][1]-zoomSubsets[0][0])/(getWidth()-WIDTH_BUFFER-WIDTH_Y_AXIS);
			} else {
				distance = (curY-startY) * (zoomSubsets[1][1]-zoomSubsets[1][0])/(getHeight()-HEAD_BUFFER-HEIGHT_X_AXIS);
			}
			
			if (distance<0) {
				distance = Math.max(distance, -1*zoomSubsets[i][0]);
			} else {
				distance = Math.min(distance, 1-zoomSubsets[i][1]);
			}

			if ((zoomSubsets[i][0]<=0&&distance<0)||(zoomSubsets[i][1]>=1&&distance>0)) {

			} else {
				zoomSubsets[i][0] += distance;
				zoomSubsets[i][1] += distance;
			}
		}
		
		if (inDrag) {
			paintAgain();
			startX = curX;
			startY = curY;
		}
	}

	public void mouseWheelMoved(MouseWheelEvent e) {
		if (zoomable) {
			if (e.getWheelRotation()<0 && zoomSubsets[0][1]-zoomSubsets[0][0] < MINIMUM_ZOOM_PROPORTION_WINDOW) {
				return;
			}			
			zoomProportionally(e.getWheelRotation()>0, e.getPoint(), false);
		}
	}

	public void componentHidden(ComponentEvent e) {}

	public void componentMoved(ComponentEvent e) {}

	public void componentResized(ComponentEvent e) {
		if (this.waitingTimer==null) {
			/* Start waiting for DELAY to elapse. */
			this.waitingTimer = new Timer(DELAY,this);
			this.waitingTimer.start();
		}
		else {
			/* Event came too soon, swallow it by resetting the timer.. */
			this.waitingTimer.restart();
		}
	}

	public void componentShown(ComponentEvent e) {}

	public void actionPerformed(ActionEvent ae) {
		/* Timer finished? */
		if (ae.getSource()==this.waitingTimer) {
			/* Stop timer */
			this.waitingTimer.stop();
			this.waitingTimer = null;
			/* Resize */
			setFlow(true);
//			System.err.println("Action performed");
			createImage();
			repaint();
		}
		
	}

	public void zoomProportionally(boolean outNotIn, Point p, boolean center) {
		float x, y, multiplier, dist;
		float[][] proportions;
		int width, height;
		boolean changed;
		
		proportions = new float[2][2];
		
		width = getWidth()-WIDTH_Y_AXIS-WIDTH_BUFFER;
		x = (float)p.getX()-WIDTH_Y_AXIS;
		proportions[0][0] = x/width;
		proportions[0][1] = (width-x)/width;

		height = getHeight()-HEIGHT_X_AXIS-HEAD_BUFFER;
		y = (float)p.getY()-HEAD_BUFFER; // could be HEAD_BUFFER
		proportions[1][0] = (height-y)/height; // reversed because of top down
		proportions[1][1] = y/height;

		multiplier = mouseWheelMultiplier/(outNotIn?1:-2);
		
		if (!outNotIn&&center) {
			for (int i = 0; i < proportions.length; i++) {
				for (int j = 0; j < proportions[i].length; j++) {
					proportions[i][j] = 0.25f-proportions[i][j];
				}
			}
			multiplier = mouseWheelMultiplier;
		}
		
		for (int i = 0; i < proportions.length; i++) {
			dist = zoomSubsets[i][1]-zoomSubsets[i][0];
			for (int j = 0; j < proportions[i].length; j++) {
				zoomSubsets[i][j] = zoomSubsets[i][j] + (j==0?-1:1)*(proportions[i][j]*multiplier*dist);
			}
		}
		
		if (zoomSubsets[0][1] - zoomSubsets[0][0] > 1) {
			resetZoomProportions();
		}

		// necessary to prevent distortion in X versus Y when zooming in/out near an edge  
		do {
			changed = false;
			for (int i = 0; i < 2; i++) {
				if (zoomSubsets[i][0] < -0.0001) {
					zoomSubsets[i][1] -= zoomSubsets[i][0];
					zoomSubsets[i][0] = 0;
					changed = true;
				}
				if (zoomSubsets[i][1] > 1.0001) {
					zoomSubsets[i][0] -= zoomSubsets[i][1]-1;
					zoomSubsets[i][1] = 1;
					changed = true;
				}
			}
		} while (changed);
		
		paintAgain();
	}

	public void resetZoomProportions() {
		zoomSubsets = new float[2][2];
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				zoomSubsets[j][i] = i;
			}
		}
	}
	
	public void drawPoint(Graphics g, PlotPoint point) {
		g.setColor(colorScheme[point.getColor()]);
		
		switch (point.getType()) {
		case PlotPoint.FILLED_CIRCLE:
			g.fillOval(getXPixel(point.getRawX())-point.getSize()/2, getYPixel(point.getRawY())-point.getSize()/2, point.getSize(), point.getSize());
			break;
		case PlotPoint.OPEN_CIRCLE:
			g.drawOval(getXPixel(point.getRawX())-point.getSize()/2, getYPixel(point.getRawY())-point.getSize()/2, point.getSize(), point.getSize());
			break;
		case PlotPoint.FILLED_SQUARE:
			g.fillPolygon(new int[] {getXPixel(point.getRawX())-point.getSize()/2, getXPixel(point.getRawX())-point.getSize()/2, getXPixel(point.getRawX())+point.getSize()/2, getXPixel(point.getRawX())+point.getSize()/2},
						  new int[] {getYPixel(point.getRawY())-point.getSize()/2, getYPixel(point.getRawY())+point.getSize()/2, getYPixel(point.getRawY())+point.getSize()/2, getYPixel(point.getRawY())-point.getSize()/2},
						  4);
			break;
		case PlotPoint.OPEN_SQUARE:
			g.drawPolygon(new int[] {getXPixel(point.getRawX())-point.getSize()/2, getXPixel(point.getRawX())-point.getSize()/2, getXPixel(point.getRawX())+point.getSize()/2, getXPixel(point.getRawX())+point.getSize()/2},
						  new int[] {getYPixel(point.getRawY())-point.getSize()/2, getYPixel(point.getRawY())+point.getSize()/2, getYPixel(point.getRawY())+point.getSize()/2, getYPixel(point.getRawY())-point.getSize()/2},
						  4);
			break;
		case PlotPoint.FILLED_TRIANGLE:
			Grafik.drawTriangle(g, getXPixel(point.getRawX()), getYPixel(point.getRawY()), point.getSize(), true);
//			g.drawPolygon(new int[] {getXPixel(point.getRawX())-point.getSize()/2, getXPixel(point.getRawX()), getXPixel(point.getRawX())+point.getSize()/2},
//						  new int[] {getYPixel(point.getRawY())+point.getSize()/2, getYPixel(point.getRawY())-point.getSize()/2, getYPixel(point.getRawY())+point.getSize()/2},
//						  3);
			break;
		case PlotPoint.MISSING:
			//g.drawString(PlotPoint.MISSING_STR, getX(point.getRawX())-missingWidth/2, getY(point.getRawY())+point.getSize()/2);
			if (PlotPoint.MISSING_STR.equals("X") || PlotPoint.MISSING_STR.equals("x")){
				g.drawLine(getXPixel(point.getRawX())-point.getSize()/4, getYPixel(point.getRawY())-point.getSize()/4, getXPixel(point.getRawX())+point.getSize()/4, getYPixel(point.getRawY())+point.getSize()/4);//zx
				g.drawLine(getXPixel(point.getRawX())-point.getSize()/4, getYPixel(point.getRawY())+point.getSize()/4, getXPixel(point.getRawX())+point.getSize()/4, getYPixel(point.getRawY())-point.getSize()/4);//zx
			} else {
				g.setFont(new Font("Arial", 0, point.getSize()));//zx
				g.drawString(PlotPoint.MISSING_STR, getXPixel(point.getRawX())-point.getSize()/2, getYPixel(point.getRawY())+point.getSize()/2);//zx
				g.setFont(new Font("Arial", 0, axisFontSize));//zx
			}
			break;
		case PlotPoint.NOT_A_NUMBER:
//			g.drawString(PlotPoint.NAN_STR, getX(point.getRawX())-nanWidth/2, getY(point.getRawY())-30+point.getSize()/2);
			break;
		default:
			System.err.println("Error - invalid PlotPoint type");
		}
	}
	
	public static double calcStepStep(double range) {
		String[] line;
		
		try {
			line = new DecimalFormat("0.0E0").format(range).split("E");
			return Math.pow(10, Integer.parseInt(line[1])-1)*(Double.parseDouble(line[0])>2.0?5:1);
		} catch (Exception e) {
			System.err.println("Error - could not parse stepStep from range '"+range+"'");
			return Double.NaN;
		}
	}

	public static int getNumSigFig(double num) {
		String str = ext.formDeci(num, 5);

		if (str.contains(".")) {
			return str.length()-str.indexOf(".")-1;
		} else {
			return 0;
		}
	}

	public double[] getPlotMinMaxStep(double min, double max, Graphics g, boolean xAxis) {
		double range, plotStep, stepStep, plotMin, plotMax;
		double zoomMin, zoomMax, dist;
		int numHashes, wid;
		FontMetrics fontMetrics;
		int sf;

		range = max-min;
//		System.out.println(min+"\t"+max+"\t"+range);
		plotStep = stepStep = calcStepStep(range);
		sf = getNumSigFig(stepStep);

		if (xAxis) {
			fontMetrics = g.getFontMetrics(g.getFont());
			wid = Math.max(fontMetrics.stringWidth(ext.formDeci(min, sf)), fontMetrics.stringWidth(ext.formDeci(max, sf)));
			numHashes = 12;
			while ((wid+30)*numHashes>canvasSectionMaximumX-canvasSectionMinimumX) {
				numHashes -= 2;
			}
		} else {
			numHashes = 10;
		}

		while (range/plotStep>numHashes) {
			plotStep += stepStep;
		}
		plotMin = plotMax = 0;
		while (max - plotMax > DOUBLE_INACCURACY_HEDGE) {
			plotMax += plotStep;
		}
		while (min - plotMin < -1*DOUBLE_INACCURACY_HEDGE) { // double check this, untested
			plotMin -= plotStep;
		}
//		if (min >= 0 && plotMin < 0) {
//			plotMin = 0;
//		}
		
//		System.out.println(Float.parseFloat(ext.formDeci(plotMin, sf))+"\t"+Float.parseFloat(ext.formDeci(plotMax, sf))+"\t"+Float.parseFloat(ext.formDeci(plotStep, sf)));
		
		if (zoomable) {
			dist = plotMax - plotMin;
			zoomMin = plotMin+zoomSubsets[xAxis?0:1][0]*dist;
			zoomMax = plotMax-(1-zoomSubsets[xAxis?0:1][1])*dist;
			
			range = zoomMax-zoomMin;
			plotStep = stepStep = calcStepStep(range);
			sf = getNumSigFig(stepStep);

			if (xAxis) {
				fontMetrics = g.getFontMetrics(g.getFont());
				wid = Math.max(fontMetrics.stringWidth(ext.formDeci(min, sf)), fontMetrics.stringWidth(ext.formDeci(max, sf)));
				numHashes = 12;
				while ((wid+30)*numHashes>canvasSectionMaximumX-canvasSectionMinimumX) {
					numHashes -= 2;
				}
			} else {
				numHashes = 10;
			}

			while (range/plotStep>numHashes) {
				plotStep += stepStep;
			}
			
			return new double[] {zoomMin, zoomMax, Double.parseDouble(ext.formDeci(plotStep, sf)), Double.parseDouble(ext.formDeci(plotMin, sf))};
		} else {
			return new double[] {Double.parseDouble(ext.formDeci(plotMin, sf)), Double.parseDouble(ext.formDeci(plotMax, sf)), Double.parseDouble(ext.formDeci(plotStep, sf)), Double.parseDouble(ext.formDeci(plotMin, sf))};
		}
	}

	public int getXPixel(double x) {
		if (invertX) {
			return (int)((plotXmax-x)/(plotXmax-plotXmin)*(double)(canvasSectionMaximumX-canvasSectionMinimumX))+canvasSectionMinimumX;
		} else {
			return (int)((x-plotXmin)/(plotXmax-plotXmin)*(double)(canvasSectionMaximumX-canvasSectionMinimumX))+canvasSectionMinimumX;
		}
	}

	public int getYPixel(double y) {
		if (invertY) {
			return getHeight()-(int)((plotYmax-y)/(plotYmax-plotYmin)*(double)(canvasSectionMaximumY-canvasSectionMinimumY)+canvasSectionMinimumY);
		} else {
			return getHeight()-(int)((y-plotYmin)/(plotYmax-plotYmin)*(double)(canvasSectionMaximumY-canvasSectionMinimumY)+canvasSectionMinimumY);
		}
	}

	/**
	 * Converts mouse location int X,Y into data points' value double rawX,rawY.
	 * The control variable invertX is assigned elsewhere in the class.
	 * @param mouseX the mouse location X
	 * @return the rawX value of the corresponding data point.
	 */
	public double getXValueFromXPixel(int mouseX) {
		if (invertX) {
			return plotXmax - ((double)(mouseX-canvasSectionMinimumX)/(double)(canvasSectionMaximumX-canvasSectionMinimumX)*(double)(plotXmax-plotXmin));
		} else {
			return plotXmin + ((double)(mouseX-canvasSectionMinimumX)/(double)(canvasSectionMaximumX-canvasSectionMinimumX)*(double)(plotXmax-plotXmin));
		}
	}

	/**
	 * Converts mouse location int X,Y into data points' value double rawX,rawY
	 * The control variable invertY is assigned elsewhere in the class.
	 * @param mouseY the mouse location Y
	 * @return the rawY value of the corresponding data point.
	 */
	public double getYValueFromYPixel(int mouseY) {
		if (invertY) {
			return plotYmax + ((double)(mouseY+canvasSectionMinimumY-getHeight())/(double)(canvasSectionMaximumY-canvasSectionMinimumY)*(double)(plotYmax-plotYmin));
		} else {
			return plotYmin - ((double)(mouseY+canvasSectionMinimumY-getHeight())/(double)(canvasSectionMaximumY-canvasSectionMinimumY)*(double)(plotYmax-plotYmin));
		}
	}
	
//	/**
//	 * Screens the data points to find out those that fall into the range of (double rawXmin, double rawXmax, double rawYmin, double rawYmax)
//	 * @param rawXmin
//	 * @param rawXmax
//	 * @param rawYmin
//	 * @param rawYmax
//	 */
//	public void highlightPoints(double rawXmin, double rawYmin, double rawXmax, double rawYmax) {
//		if (points.length==100){
//		for (int i=0; i<points.length; i++) {
//			if (points[i].getRawX()>=rawXmin
//					&& points[i].getRawY()>=rawYmin
//					&& points[i].getRawX()<=rawXmax
//					&& points[i].getRawY()<=rawYmax) {
//				points[i].setHighlighted(true);
//				System.out.println("Highlighting: "+points[i].getRawX()+","+points[i].getRawY());
//			} else {
//				points[i].setHighlighted(false);
//			}
//			
//		}
//		}
//	}
//	
	/**
	 * Highlights those points that need to be highlighted
	 * @param array
	 */
	public void highlightPoints(boolean[] array) {
		if (points.length != array.length) {
			System.err.println("Error - mismatched array size when highlighting");
		} else {
			for (int i=0; i<points.length; i++) {
				if (points[i]!=null && array[i]) {
					points[i].setHighlighted(true);
				} else if (points[i]!=null) {
					points[i].setHighlighted(false);
				}
				
			}
		}
	}
	
	public void createImage() {
//		System.err.println("Creating Image");
		if (getWidth() > 350 && getHeight() > 0 ) {
			image = new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_INT_RGB);
			flow = true;
			drawAll(image.createGraphics(), true);
//			repaint();
//			image = new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_INT_RGB);
		}
	}
	
	public void setLookupResolution(int lookupResolution) {
		this.lookupResolution = lookupResolution;
	}
	
	public IntVector lookupNearbyPoints(int x, int y, String pos) {
		IntVector iv = locLookup.get(pos);
		IntVector indicesOfDataPoints = new IntVector();

		for (int i = 0; iv!=null && i<iv.size(); i++) {
			if (Distance.euclidean(new int[] {x, y}, new int[] {getXPixel(points[iv.elementAt(i)].getRawX()), getYPixel(points[iv.elementAt(i)].getRawY())})<HIGHLIGHT_DISTANCE) {
				indicesOfDataPoints.add(iv.elementAt(i));
			}
		}

		return indicesOfDataPoints;
	}

	public byte lookupNearbyRectangles(int xPixel, int yPixel) {
		int rectangleStartXPixel, rectangleStartYPixel, rectangleStopXPixel, rectangleStopYPixel;
		byte indicesRectangle;
		int distance;
		int minDistance;
		int[] rectangleEdges;

		indicesRectangle = -1;
		minDistance = Integer.MAX_VALUE;
		for (byte i=0; rectangles != null && i < rectangles.length; i++) {
			rectangleStartXPixel = getXPixel(rectangles[i].getStartXValue());
			rectangleStartYPixel = getYPixel(rectangles[i].getStartYValue());
			rectangleStopXPixel = getXPixel(rectangles[i].getStopXValue());
			rectangleStopYPixel = getYPixel(rectangles[i].getStopYValue());

	    	rectangleEdges = new int[] {Math.min(rectangleStartXPixel, rectangleStopXPixel), Math.max(rectangleStartXPixel, rectangleStopXPixel), Math.min(rectangleStartYPixel, rectangleStopYPixel), Math.max(rectangleStartYPixel, rectangleStopYPixel)};
	    	for (int j = 0; j < 2; j ++) {
	    		distance = Math.abs(xPixel - rectangleEdges[j]);
	    		if (distance < HIGHLIGHT_DISTANCE && distance < minDistance && yPixel > rectangleEdges[2] - HIGHLIGHT_DISTANCE && yPixel < rectangleEdges[3] + HIGHLIGHT_DISTANCE) {
		    		minDistance = distance;
		    		indicesRectangle = i;
	    		}
			}
	    	for (int j = 2; j < 4; j ++) {
	    		distance = Math.abs(yPixel - rectangleEdges[j]);
	    		if (distance < HIGHLIGHT_DISTANCE && distance < minDistance && xPixel > rectangleEdges[0] - HIGHLIGHT_DISTANCE && xPixel < rectangleEdges[1] + HIGHLIGHT_DISTANCE) {
		    		minDistance = distance;
		    		indicesRectangle = i;
	    		}
	    	}
		}

		return indicesRectangle;
	}

	public byte lookupResidingRectangles(int xPixel, int yPixel) {
		int rectangleStartXPixel, rectangleStartYPixel, rectangleStopXPixel, rectangleStopYPixel;
		byte indicesRectangle;
		int distance;
		int minDistance;
		int[] rectangleEdges;
		int[] rectangleSearchRange;

		indicesRectangle = -1;
		minDistance = Integer.MAX_VALUE;
		for (byte i=0; rectangles != null && i < rectangles.length; i++) {
			rectangleStartXPixel = getXPixel(rectangles[i].getStartXValue());
			rectangleStartYPixel = getYPixel(rectangles[i].getStartYValue());
			rectangleStopXPixel = getXPixel(rectangles[i].getStopXValue());
			rectangleStopYPixel = getYPixel(rectangles[i].getStopYValue());

	    	rectangleEdges = new int[] {(int) Math.min(rectangleStartXPixel, rectangleStopXPixel),
						    			(int) Math.max(rectangleStartXPixel, rectangleStopXPixel),
						    			(int) Math.min(rectangleStartYPixel, rectangleStopYPixel),
						    			(int) Math.max(rectangleStartYPixel, rectangleStopYPixel)};
	    	rectangleSearchRange = new int[] {(int) (rectangleEdges[0] - HIGHLIGHT_DISTANCE),
	    									  (int) (rectangleEdges[1] + HIGHLIGHT_DISTANCE),
	    									  (int) (rectangleEdges[2] - HIGHLIGHT_DISTANCE),
	    									  (int) (rectangleEdges[3] + HIGHLIGHT_DISTANCE)};
	    	if (xPixel >= rectangleSearchRange[0] && xPixel <= rectangleSearchRange[1] && yPixel >= rectangleSearchRange[2] && yPixel <= rectangleSearchRange[3]) {
		    	for (int j = 0; j < 2; j ++) {
		    		distance = Math.abs(xPixel - rectangleEdges[j]);
		    		if (minDistance > distance) {
			    		minDistance = distance;
			    		indicesRectangle = i;
		    		}
				}
		    	for (int j = 2; j < 4; j ++) {
		    		distance = Math.abs(yPixel - rectangleEdges[j]);
		    		if (minDistance > distance) {
			    		minDistance = distance;
			    		indicesRectangle = i;
		    		}
				}
	    	}
		}

		return indicesRectangle;
	}

	public void setFinalImage(boolean finalImage) {
		this.finalImage = finalImage;
	}

	public boolean getFinalImage() {
		return finalImage;
	}

	public void setFlow(boolean flow) {
		this.flow = flow;
	}

	public boolean getFlow() {
		return flow;
	}
	
	public void setPointsGeneratable(boolean pointsGeneratable) {
		this.pointsGeneratable = pointsGeneratable;
	}

	public boolean isPointsGeneratable() {
		return pointsGeneratable;
	}
	
	public void setSwapable(boolean swapable) {
		this.swapable = swapable;
	}

	public boolean isSwapable() {
		return swapable;
	}

	public void drawRectThick (Graphics g, int x, int y, int width, int height, byte thickness) {
    	for (byte i=0; i<thickness; i++) {
    		g.drawRect(x+i, y+i, width-2*i, height-2*i);
    	}
	}
}
