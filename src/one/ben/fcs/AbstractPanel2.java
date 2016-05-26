package one.ben.fcs;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.RenderingHints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Date;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

import javax.imageio.ImageIO;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.Timer;

import mining.Distance;
import stats.Maths;
import cnv.plots.GenericLine;
import cnv.plots.GenericPath;
import cnv.plots.GenericRectangle;
import cnv.plots.PlotPoint;
import common.Array;
import common.Grafik;
import common.HashVec;
import common.IntVector;
import common.ProgressBarDialog;
import common.Sort;
import common.ext;
import edu.stanford.facs.logicle.FastLogicle;
import edu.stanford.facs.logicle.Logicle;

public abstract class AbstractPanel2 extends JPanel implements MouseListener, MouseMotionListener, MouseWheelListener, ComponentListener, ActionListener {
	public static final long serialVersionUID = 1L;

	public static final boolean DEBUGGING = false;

	public static final int HEAD_BUFFER = 25;
	public static final int HEIGHT_X_AXIS = 105;
	public static final int WIDTH_Y_AXIS = 140;
	public static final int WIDTH_BUFFER = 50;
	public static final int AXIS_THICKNESS = 4;
	public static final int TICK_THICKNESS = 2;
	public static final int TICK_LENGTH = 15;
	public static final int DEFAULT_LOOKUP_RESOLUTION = 20;
	public static final int AXIS_FONT_SIZE = 28;
	public static final double DOUBLE_INACCURACY_HEDGE = 0.00001;
	public static final double MINIMUM_ZOOM_PROPORTION_WINDOW = 0.0001;
	public static final float DEFAULT_MOUSE_WHEEL_MULTIPLIER = 0.5f;
	public static final int DEFAULT_PLOTPOINTSET_SIZE = 10;
	public static final int SIZE = 12;
	public static final double HIGHLIGHT_DISTANCE = 20;//= Math.sqrt(SIZE*SIZE/2);
	public static final int DELAY = 2;	//A control variable to reduce the repaint() operations during component resizing;

    public static final String X_MIN = "xMin";
    public static final String X_MAX = "xMax";
    public static final String Y_MIN = "yMin";
    public static final String Y_MAX = "yMax";
	
	public static enum PLOT_TYPE {
	    DOT_PLOT,
	    HEATMAP,
        HISTOGRAM,
        CONTOUR;
	}
	
    public static enum AXIS_SCALE {
        LIN("Linear"),
        LOG("Logarithmic")
        ,
        BIEX("Biexponential")
        ;
        
        private AXIS_SCALE(String longName) {
            this.longName = longName;
        }
        public String longName;
        @Override
        public String toString() {
            return longName;
        }
    }
//    http://superliminal.com/sources/Axis.java.html
	
//	public static final int SCATTER_PLOT_TYPE = 1; 
//	public static final int HEAT_MAP_TYPE = 2;

	public static final PLOT_TYPE DEFAULT_TYPE = PLOT_TYPE.DOT_PLOT;

	public static final int IMAGE_NULL = 0;
	public static final int IMAGE_STARTED = 1;
	public static final int IMAGE_COMPLETE = 2;
	
	protected Color[] colorScheme;
	protected int canvasSectionMinimumX;
	protected int canvasSectionMaximumX;
	protected int canvasSectionMinimumY;
	protected int canvasSectionMaximumY;
	protected int axisYWidth = WIDTH_Y_AXIS;
	protected int axisXHeight = HEIGHT_X_AXIS;
	protected AXIS_SCALE xAxis = AXIS_SCALE.LOG;
	protected AXIS_SCALE yAxis = AXIS_SCALE.LOG;
	protected double plotXmax, plotYmax;
	protected double plotXmin, plotYmin;
	protected volatile BufferedImage image;
	protected String prevPos = "";
	protected Hashtable<String,IntVector> locLookup;
	protected PlotPoint[] points;						// make private when worked out
	protected GenericLine[] lines;
	protected GenericRectangle[] rectangles;
	protected GenericPath[] polygons;
	protected GenericLine highlightLine;
	protected GenericRectangle highlightRectangle;
	protected GenericPath highlightPoly;
	protected String xAxisLabel;
	protected String yAxisLabel;
	protected String title;
	protected boolean displayXaxis;
	protected boolean displayYaxis;
	protected boolean displayGrid;
	protected boolean displayTitle;
	protected boolean xAxisWholeNumbers;
	protected int missingWidth;
	protected int nanWidth;
	protected int axisFontSize = AXIS_FONT_SIZE;
	protected float forcePlotXmax, forcePlotYmax;
	protected float forcePlotXmin, forcePlotYmin;
	protected boolean createLookup;
	protected volatile boolean invertX;
	protected volatile boolean invertY;
	protected boolean makeSymmetric;
	protected String errorMessage;
	protected float mouseWheelMultiplier;
	protected boolean zoomable;
	protected boolean swapable;		//4/27/2012
	protected boolean invertable;	//4/27/2012
	protected boolean truncate;
	protected float[][] zoomSubsets;
	protected IntVector prox;
	protected PLOT_TYPE chartType;

	private boolean inDrag;
	protected volatile int startX;
    protected volatile int startY;
//	private int titleX, titleY;
	private int titleLocation;
	private int plotPointSetSize;
	private int totalNumPlotPointSets;
	private int currentPlotPointSet;
	private int lastIndexInPlotPointSet;
	private int currentIndexInPlotPointSet;
	private int lookupResolution;
	private boolean flow;			//A control variable. If resizing is not yet done, don't start generatePoints() or drawAll();
	private volatile int imageStatus;		//A control variable. If drawAll() is not yet done, don't start paintComponent();
	private byte[] layersInBase;
	private byte[] extraLayersVisible;
//	private boolean pointsGeneratable;
	protected Timer waitingTimer;		//A control variable to reduce the repaint() operations during component resizing;
	private String nullMessage;
	private boolean randomTest;
	private int numberOfNaNSamples;
	private boolean antiAlias = true;
	private volatile boolean beEfficient;
	private HashSet<String> pointsPlotted;
	
	public AbstractPanel2() {
		canvasSectionMinimumX = 0;
		canvasSectionMaximumX = getWidth();
		canvasSectionMinimumY = 0;
		canvasSectionMaximumY = getWidth();
		displayXaxis = true;
		displayYaxis = true;
		displayGrid = false;
		createLookup = false;
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
		imageStatus = IMAGE_NULL;
		flow=true;
//		pointsGeneratable = true;
		beEfficient = true;
		
		colorScheme = new Color[] {Color.BLACK, Color.GRAY};
		addMouseListener(this);
		addMouseMotionListener(this);
		addMouseWheelListener(this);
		addComponentListener(this);
	}
	
	
	public void createLookup(boolean value) {
		this.createLookup = value;
	}
	
	public void setImageStatus(int status) {
	    if (DEBUGGING) {
	        System.out.println("Set image status to " + status);
	    }
		imageStatus = status;
	}

	protected boolean imageIsFinal() {
		return imageStatus == IMAGE_COMPLETE;
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

	public void beEfficient(boolean value) {
		this.beEfficient = value;
	}

	public void paintAgain() {
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
        		image = null;
        		setImageStatus(IMAGE_NULL);
        		repaint();
            }
        });
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

	public void paintComponent(final Graphics g) {
	    if (waitingTimer != null) return;
	    if (imageStatus == IMAGE_STARTED) return;
		// this either doesn't affect anything or gets caught in an infinite loop while the lookup is being created
//		while (imageStatus == IMAGE_STARTED) {
//			if (DEBUGGING) {
//				System.out.println("Additional call to paint before the first was completed; sleeping 100ms");
//			}
//			try {
//				Thread.sleep(100);
//			} catch (InterruptedException ie) {
//			}
//		}
		if (image == null) {
			if (DEBUGGING) {
				System.out.println("createImage() being called from paintComponent()");
			}    
            createImage(); // if you remove this, you get a blank screen and at least QQPlot ends up with a double title panel
		} else if (DEBUGGING) {
	        System.out.println("Skipping image creation");
	    }
		g.drawImage(image, 0, 0, AbstractPanel2.this);
		if (extraLayersVisible != null && extraLayersVisible.length > 0) {
			drawAll(g, false);
		}
	}
	
	public void screenCapture(String filename) {
		try {
		    File imgFile = new File(filename);
		    boolean mkdirs = imgFile.mkdirs();
		    if (mkdirs) {
		        ImageIO.write(image, "png", imgFile);
		    } else {
		        JOptionPane.showMessageDialog(null, "Error while trying to save the plot");
		    }
		} catch (IOException ie) {
			JOptionPane.showMessageDialog(null, "Error while trying to save the plot");
		}
	}
	
	public void setPlotXMin(float plotXmin) {
	    double old = this.plotXmin;
	    this.plotXmin = plotXmin;
	    firePropertyChange(X_MIN, old, plotXmin);
	}
	public void setPlotXMax(float plotXmax) {
	    double old = this.plotXmax;
	    this.plotXmax = plotXmax;
	    firePropertyChange(X_MAX, old, plotXmax);
	}
	public void setPlotYMin(float plotYmin) {
	    double old = this.plotYmin;
	    this.plotYmin = plotYmin;
	    firePropertyChange(Y_MIN, old, plotYmin);
	}
	public void setPlotYMax(float plotYmax) {
        double old = this.plotYmax;
	    this.plotYmax = plotYmax;
        firePropertyChange(Y_MAX, old, plotYmax);
	}
	
	public void setForcePlotXMin(float forcePlotXmin) {
		this.forcePlotXmin = forcePlotXmin;
	}

	public void setForcePlotXMax(float forcePlotXmax) {
		this.forcePlotXmax = forcePlotXmax;
	}

	public void setForcePlotYMin(float forcePlotYmin) {
		this.forcePlotYmin = forcePlotYmin;
	}

	public void setForcePlotYMax(float forcePlotYmax) {
		this.forcePlotYmax = forcePlotYmax;
	}
	
	public void setForceXAxisWholeNumbers(boolean whole) {
	    this.xAxisWholeNumbers = whole;
	}

	public AXIS_SCALE getXAxis() {
		return xAxis;
	}
	public void setXAxis(AXIS_SCALE xAxis) {
		this.xAxis = xAxis;
	}
	public AXIS_SCALE getYAxis() {
		return yAxis;
	}
	public void setYAxis(AXIS_SCALE yAxis) {
		this.yAxis = yAxis;
	}

	public abstract void generatePointsRectanglesAndLines();

	public abstract void highlightPoints();

	public abstract void assignAxisLabels();

	public void setAxisFontSize(int sz) {
		if (sz > 0) {
			this.axisFontSize = sz;
		}
	}
	
	public void setChartType(PLOT_TYPE chartType) {
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
	
	public Path2D getPixelPath(Path2D dataPath) {
        PathIterator pathIter = dataPath.getPathIterator(null);
        Path2D path = new Path2D.Double();
        while (!pathIter.isDone()) {
            double[] pts = new double[6];
            int type = pathIter.currentSegment(pts);
            switch (type) {
            case PathIterator.SEG_MOVETO:
                path.moveTo(getXPixel(pts[0]), getYPixel(pts[1]));
                break;
            case PathIterator.SEG_LINETO:
                path.lineTo(getXPixel(pts[0]), getYPixel(pts[1]));
            case PathIterator.SEG_QUADTO:
                path.quadTo(getXPixel(pts[0]), getYPixel(pts[1]), getXPixel(pts[2]), getYPixel(pts[3]));
                break;
            case PathIterator.SEG_CLOSE:
                path.closePath();
            }
            pathIter.next();
        }
        return path;
	}

	public void drawAll(Graphics g, boolean base) {
		float minimumObservedRawX, maximumObservedRawX, minimumObservedRawY, maximumObservedRawY;
		double[] plotMinMaxStep; // needs to be double, else x <= plotXmax can be inexact and leave off the last tick mark 
		int sigFigs;
		String str, pos;
		int xLook, yLook;
		BufferedImage yLabel;
		FontMetrics fontMetrics = null;
		Graphics gfx;
		Hashtable<String, Vector<PlotPoint>> layers;
		Vector<PlotPoint> layer;
		String trav;
		String[] keys;
		int[] order;
		int step;
		long time;
		ProgressBarDialog prog;
    	int xPixel, yPixel, widthPixel, heightPixel;

		setImageStatus(IMAGE_STARTED);
    	
    	long fullTime = System.currentTimeMillis();
    	
    	pointsPlotted = new HashSet<String>();
    	
    	if (g instanceof Graphics2D) {
    		((Graphics2D)g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, antiAlias ? RenderingHints.VALUE_ANTIALIAS_ON : RenderingHints.VALUE_ANTIALIAS_OFF);
    		((Graphics2D)g).setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, antiAlias ? RenderingHints.VALUE_TEXT_ANTIALIAS_ON : RenderingHints.VALUE_TEXT_ANTIALIAS_OFF);
    	}
    	
        generatePointsRectanglesAndLines();
		highlightPoints();
		
		if (points.length == 0 && (rectangles == null || rectangles.length == 0) && (polygons == null || polygons.length == 0) && (lines == null || lines.length == 0)) {
		    locLookup.clear();
			g.setColor(Color.WHITE);
			g.fillRect(0, 0, getWidth(), getHeight());
			if (nullMessage != null) {
				g.setColor(Color.BLACK);
				g.drawString(nullMessage, getWidth()/2-g.getFontMetrics(g.getFont()).stringWidth(nullMessage)/2, getHeight()/2);
			}
			setImageStatus(IMAGE_COMPLETE);
			return;
		}

		setLookupResolution(DEFAULT_LOOKUP_RESOLUTION);
		assignAxisLabels();

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

		for (int i = 0; lines != null && i<lines.length && flow; i++) {
			if (lines[i] != null && lines[i].getLayer() == 0) {
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
		for (int i = 0; rectangles != null && i < rectangles.length && flow; i++) {
		    if (rectangles[i] != null && rectangles[i].getLayer() == 0) {
		        minimumObservedRawX = Maths.min(minimumObservedRawX, rectangles[i].getStartXValue());
		        maximumObservedRawX = Maths.max(maximumObservedRawX, rectangles[i].getStartXValue());
		        minimumObservedRawY = Maths.min(minimumObservedRawY, rectangles[i].getStartYValue());
		        maximumObservedRawY = Maths.max(maximumObservedRawY, rectangles[i].getStartYValue());
		        
		        minimumObservedRawX = Maths.min(minimumObservedRawX, rectangles[i].getStopXValue());
		        maximumObservedRawX = Maths.max(maximumObservedRawX, rectangles[i].getStopXValue());
		        minimumObservedRawY = Maths.min(minimumObservedRawY, rectangles[i].getStopYValue());
		        maximumObservedRawY = Maths.max(maximumObservedRawY, rectangles[i].getStopYValue());
		    }
		}
		
		minimumObservedRawX = minimumObservedRawX == Float.MAX_VALUE ? 0 : minimumObservedRawX;
		maximumObservedRawX = maximumObservedRawX == Float.MIN_VALUE ? 1 : maximumObservedRawX;
		minimumObservedRawY = minimumObservedRawY == Float.MAX_VALUE ? 0 : minimumObservedRawY;
		maximumObservedRawY = maximumObservedRawY == Float.MIN_VALUE ? 1 : maximumObservedRawY;

//		otherwise step is off
		minimumObservedRawX = minimumObservedRawX > 0 ? 0 : minimumObservedRawX;
		minimumObservedRawY = minimumObservedRawY > 0 ? 0 : minimumObservedRawY;
		
		if (Float.isNaN(forcePlotXmin)) {
//			minimumObservedRawX = minimumObservedRawX;
		} else {
			if (forcePlotXmin > minimumObservedRawX) {
			    if (DEBUGGING) {
			        System.err.println("WARNING - specified [minimum X boundary : " + forcePlotXmin + "] is higher than the data point with the [lowest X value : " + minimumObservedRawX + "]");
			    }
			}
			minimumObservedRawX = forcePlotXmin;
		}
		if (Float.isNaN(forcePlotXmax)) {
			maximumObservedRawX = (maximumObservedRawX + (maximumObservedRawX - minimumObservedRawX) * (float) 0.01);
		} else {
			if (forcePlotXmax < maximumObservedRawX) {
			    if (DEBUGGING) {
			        System.err.println("WARNING - specified [maximum X boundary : " + forcePlotXmax + "] is lower than the data point with the [highest X value : " + maximumObservedRawX + "]");
			    }
			}
			maximumObservedRawX = forcePlotXmax;
		}
		if (Float.isNaN(forcePlotYmin)) {
//			minimumObservedRawY = minimumObservedRawY;
		} else {
			if (forcePlotYmin > minimumObservedRawY) {
			    if (DEBUGGING) {
			        System.err.println("WARNING - specified [minimum Y boundary : " + forcePlotYmin + "] is higher than the data point with the [lowest Y value : " + minimumObservedRawY + "]");
			    }
			}
			minimumObservedRawY = forcePlotYmin;
		}
		if (Float.isNaN(forcePlotYmax)) {
			maximumObservedRawY = (maximumObservedRawY + (maximumObservedRawY - minimumObservedRawY) * (float) 0.01);
		} else {
			if (forcePlotYmax < maximumObservedRawY) {
			    if (DEBUGGING) {
			        System.err.println("WARNING - specified [maximum Y boundary : " + forcePlotYmax + "] is lower than the data point with the [highest Y value : " + maximumObservedRawY + "]");
			    }
			}
			maximumObservedRawY = forcePlotYmax;
		}

//        minimumObservedRawX = Float.isNaN(forcePlotXmin) ? minimumObservedRawX : forcePlotXmin;
//        maximumObservedRawX = Float.isNaN(forcePlotXmax) ? (maximumObservedRawX + (maximumObservedRawX - minimumObservedRawX) * (float) 0.01) : forcePlotXmax;
//        minimumObservedRawY = Float.isNaN(forcePlotYmin) ? minimumObservedRawY : forcePlotYmin;
//        maximumObservedRawY = Float.isNaN(forcePlotYmax) ? (maximumObservedRawY + (maximumObservedRawY - minimumObservedRawY) * (float) 0.01) : forcePlotYmax;
		
		if (makeSymmetric) {
			maximumObservedRawX = Math.max(maximumObservedRawX, maximumObservedRawY);
			maximumObservedRawY = maximumObservedRawX;
			minimumObservedRawX = Math.min(minimumObservedRawX, minimumObservedRawY);
			minimumObservedRawY = minimumObservedRawX;
		}
		
		numberOfNaNSamples = 0;
		if (base) {
		    if (DEBUGGING) {
		        System.out.println("Drawing base image.");
		    }
		    
			//g.setColor(Color.WHITE);
			g.fillRect(0, 0, getWidth(), getHeight());
			g.setFont(new Font("Arial", 0, axisFontSize));
	//		System.out.println("getWidth: "+getWidth()+"\t getHeight: "+getHeight());
			
			fontMetrics = g.getFontMetrics(g.getFont());
			missingWidth = fontMetrics.stringWidth("X");
			missingWidth = fontMetrics.stringWidth("X");

			// Calculate the plot area's range (X-axis, Y-axis)
			plotMinMaxStep = null;
			canvasSectionMinimumX = axisYWidth;//WIDTH_Y_AXIS;
			canvasSectionMaximumX = getWidth()-WIDTH_BUFFER;
			canvasSectionMinimumY = 0;
			canvasSectionMaximumY = axisXHeight;//HEIGHT_X_AXIS;
			plotMinMaxStep = getPlotMinMaxStep(minimumObservedRawX, maximumObservedRawX, g, true);
			if (plotMinMaxStep[1] > maximumObservedRawX) {
			    plotMinMaxStep[1] = maximumObservedRawX;
			}
			if (xAxisWholeNumbers) {
			    if (plotMinMaxStep[2] < 1) {
			        plotMinMaxStep[2] = 1;
			    } else {
			        plotMinMaxStep[2] = Math.round(plotMinMaxStep[2]);
			    }
			    if (plotMinMaxStep[3] >= (minimumObservedRawX - plotMinMaxStep[2])) {
			        plotMinMaxStep[3] = plotMinMaxStep[3] - plotMinMaxStep[2];
			    }
			}
			setPlotXMin((float) (Float.isNaN(forcePlotXmin) ? plotMinMaxStep[0] : forcePlotXmin));
			setPlotXMax((float) (Float.isNaN(forcePlotXmax) ? plotMinMaxStep[1] : forcePlotXmax));
			
	        if (displayXaxis) {
	            drawXAxis(g, plotMinMaxStep, fontMetrics);
	        } 
	        if (xAxisLabel != null) {
	            Grafik.drawThickLine(g, canvasSectionMinimumX-(int)Math.ceil((double)AXIS_THICKNESS/2.0), getHeight()-canvasSectionMaximumY, canvasSectionMaximumX+(int)Math.ceil((double)AXIS_THICKNESS/2.0), getHeight()-canvasSectionMaximumY, AXIS_THICKNESS, Color.BLACK);
	            g.drawString(xAxisLabel, (getWidth()-axisYWidth)/2 - fontMetrics.stringWidth(xAxisLabel)/2 + 2*(axisYWidth/3), getHeight()-20);
	        }
		
	        if (DEBUGGING) {
	            System.out.println("Finished drawing x-axis");
	        }
	        
			canvasSectionMinimumX = 0;
			canvasSectionMaximumX = axisYWidth;//WIDTH_Y_AXIS;
			canvasSectionMinimumY = axisXHeight;//HEIGHT_X_AXIS;
			canvasSectionMaximumY = getHeight()-HEAD_BUFFER;
			if (!makeSymmetric || plotMinMaxStep == null) {
				plotMinMaxStep = getPlotMinMaxStep(minimumObservedRawY, maximumObservedRawY, g, false);
	            if (plotMinMaxStep[1] > maximumObservedRawY) {
	                plotMinMaxStep[1] = maximumObservedRawY;
	            }
			}
			setPlotYMin((float) plotMinMaxStep[0]);
			setPlotYMax((float) plotMinMaxStep[1]);
	        if (displayYaxis && yAxisLabel != null) {
	            drawYAxis(g, plotMinMaxStep);
	        }

			if (errorMessage != null) {
				g.drawString(errorMessage, (getWidth()-axisYWidth/*WIDTH_Y_AXIS*/)/2-fontMetrics.stringWidth(errorMessage)/2+axisYWidth/*WIDTH_Y_AXIS*/, (getHeight()-HEAD_BUFFER-axisXHeight/*HEIGHT_X_AXIS*/)/2-20+HEAD_BUFFER);
			}
			
		}
		
		if (DEBUGGING) {
		    System.out.println("Finished drawing axes");
		}
		
		//TODO outercoordinates
		canvasSectionMinimumX = axisYWidth;//WIDTH_Y_AXIS;
		canvasSectionMaximumX = getWidth()-WIDTH_BUFFER;
		canvasSectionMinimumY = axisXHeight;//HEIGHT_X_AXIS;
		canvasSectionMaximumY = getHeight()-HEAD_BUFFER;


        // Draw data points, also build the lookup matrix for nearby points.
        locLookup.clear();
        time = new Date().getTime();
        step = Math.max((points.length)/100, 1);
        layers = new Hashtable<String,Vector<PlotPoint>>();

//        long t1 = System.currentTimeMillis();

        if (chartType == PLOT_TYPE.HEATMAP) {
            drawHeatMap(g);
        } else if (chartType == PLOT_TYPE.DOT_PLOT) {
            for (int i = 0; i<points.length && flow; i++) {
                if (points[i] == null || points[i].getColor() == -1 || !points[i].isVisible()) {
                    
                } else if (truncate && (points[i].getRawX() < plotXmin || points[i].getRawX() - plotXmax > plotXmax/1000.0 || points[i].getRawY() < plotYmin || points[i].getRawY() > plotYmax)) {
//                  System.err.println("error: data point ("+points[i].getRawX()+","+points[i].getRawY()+") is outside of plot range.");
                } else {
                    
                    if (points[i].isHighlighted() || (base && (layersInBase == null || Array.indexOfByte(layersInBase, points[i].getLayer()) >= 0)) || (!base && Array.indexOfByte(extraLayersVisible, points[i].getLayer()) >= 0)) {
                        if (points[i].getLayer() == 0) {
                            if (points[i].getType()!=PlotPoint.NOT_A_NUMBER) {
                                drawPoint(g, points[i]);
                            } else if (base) {
                                numberOfNaNSamples++;
                            }
                        } else {
                            trav = points[i].getLayer()+"";
                            if (layers.containsKey(trav)) {
                                layer = layers.get(trav);
                            } else {
                                layers.put(trav, layer = new Vector<PlotPoint>(points.length));
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
        }/* else {
            System.err.println("Error - invalid chart type: "+chartType);
        }*/
        
//        System.out.println("pts: " + ext.getTimeElapsed(t1));
        
        g.setClip((int)canvasSectionMinimumX, HEAD_BUFFER, (int)(canvasSectionMaximumX - canvasSectionMinimumX) + 1, (int)(getHeight() - axisXHeight - 24));
        
		// Draw the lines
        for (int i = 0; lines!=null && i<lines.length && flow; i++) {
            if ((base && (layersInBase == null || Array.indexOfByte(layersInBase, lines[i].getLayer()) >= 0)) || (!base && Array.indexOfByte(extraLayersVisible, lines[i].getLayer()) >= 0)) {
                double x1raw = lines[i].getStartX();
                double y1raw = lines[i].getStartY();
                double x2raw = lines[i].getStopX();
                double y2raw = lines[i].getStopY();

                double x1final = x1raw, y1final = y1raw, x2final = x2raw, y2final = y2raw;
                
                int x1 = getXPixel(x1final);
                int y1 = getYPixel(y1final);
                int x2 = getXPixel(x2final);
                int y2 = getYPixel(y2final);
                if (lines[i].getThickness() == 1) {
                    Color c = g.getColor();
                    g.setColor(colorScheme[lines[i].getColor()]);
                    g.drawLine(x1, y1, x2, y2);
                    g.setColor(c);
                } else {
                    Grafik.drawThickLine(g, x1, y1, 
                            x2, y2, 
                            (int)lines[i].getThickness(), 
                            colorScheme[lines[i].getColor()], 
                            lines[i].getDirection(), lines[i].getShouldScale());
                }
            }
        }

		for (int i = 0; rectangles!=null && i<rectangles.length && flow; i++) {
			if ((base && (layersInBase == null || Array.indexOfByte(layersInBase, rectangles[i].getLayer()) >= 0)) || (!base && Array.indexOfByte(extraLayersVisible, rectangles[i].getLayer()) >= 0)) {
				xPixel = Math.min(getXPixel(rectangles[i].getStartXValue()), getXPixel(rectangles[i].getStopXValue()));
				yPixel = Math.min(getYPixel(rectangles[i].getStartYValue()), getYPixel(rectangles[i].getStopYValue()));
		    	widthPixel = Math.abs(getXPixel(rectangles[i].getStartXValue()) - getXPixel(rectangles[i].getStopXValue()));
		    	heightPixel = (Math.abs(getYPixel(rectangles[i].getStartYValue()) - getYPixel(rectangles[i].getStopYValue())));
		    	g.setColor(colorScheme[rectangles[i].getColor()]);
				
				if (rectangles[i].getFill()) {
					if (rectangles[i].getColor() != rectangles[i].getFillColor()) {
						if (rectangles[i].getRoundedCorners()) {
							g.drawRoundRect(xPixel, yPixel, widthPixel, heightPixel, 2, 2);
							g.setColor(colorScheme[rectangles[i].getFillColor()]);
							g.fillRoundRect(xPixel + rectangles[i].getThickness(), yPixel + rectangles[i].getThickness(), widthPixel - (rectangles[i].getThickness() * 2) + 1, heightPixel - (rectangles[i].getThickness() * 2) + 1, 2, 2);
							g.setColor(colorScheme[rectangles[i].getColor()]);
						} else {
					    	drawRectThick(g, xPixel, yPixel, widthPixel, heightPixel, rectangles[i].getThickness());
					    	g.setColor(colorScheme[rectangles[i].getFillColor()]);
							g.fillRect(xPixel + rectangles[i].getThickness(), yPixel + rectangles[i].getThickness(), widthPixel - (rectangles[i].getThickness() * 2) + 1, heightPixel - (rectangles[i].getThickness() * 2) + 1);
							g.setColor(colorScheme[rectangles[i].getColor()]);
						}
					} else {
						if (rectangles[i].getRoundedCorners()) {
							g.fillRoundRect(xPixel, yPixel, widthPixel, heightPixel, 2, 2);
						} else {
							g.fillRect(xPixel, yPixel, widthPixel, heightPixel);
						}
					}
				} else {
					if (rectangles[i].getRoundedCorners()) {
						g.drawRoundRect(xPixel, yPixel, widthPixel, heightPixel, 2, 2);
					} else {
				    	drawRectThick(g, xPixel, yPixel, widthPixel, heightPixel, rectangles[i].getThickness());
					}
				}
				if (rectangles[i].getEditable()) { // draw corner/control points
				    g.setColor(Color.black);
				    g.fillRect(xPixel - 2, yPixel - 2, 4, 4);
				    g.fillRect(xPixel + widthPixel - 2, yPixel - 2, 4, 4);
				    g.fillRect(xPixel + widthPixel - 2, yPixel + heightPixel - 2, 4, 4);
				    g.fillRect(xPixel - 2, yPixel + heightPixel - 2, 4, 4);
				}
			}
        }
        
		class PolyTransformer {
		    
		    Path2D transform(Path2D poly) {
		        Path2D pixPoly = new Path2D.Double(Path2D.WIND_EVEN_ODD);
		        double[] coords = new double[6];
		        PathIterator path = poly.getPathIterator(null);
		        while(!path.isDone()) {
		            int code = path.currentSegment(coords);
		            switch(code) {
		                case PathIterator.SEG_MOVETO:
		                    pixPoly.moveTo(getXPixel(coords[0]), getYPixel(coords[1]));
		                    break;
		                case PathIterator.SEG_LINETO:
		                    pixPoly.lineTo(getXPixel(coords[0]), getYPixel(coords[1]));
		                    break;
		                case PathIterator.SEG_CLOSE:
		                    pixPoly.closePath();
		            }
	                path.next();
		        }
		        return pixPoly;
		    }
		    
		}
		
		PolyTransformer pt = new PolyTransformer();
		
		for (int i = 0; polygons !=null && i < polygons.length && flow; i++) {
		    Graphics2D g2d = (Graphics2D) g;
		    if ((base && (layersInBase == null || Array.indexOfByte(layersInBase, polygons[i].getLayer()) >= 0)) || (!base && Array.indexOfByte(extraLayersVisible, polygons[i].getLayer()) >= 0)) {
		        g.setColor(colorScheme[polygons[i].getColor()]);
		        Path2D drawPoly = pt.transform(polygons[i].myPath);
		        if (polygons[i].getFill()) {
		            if (polygons[i].getColor() != polygons[i].getFillColor()) {
                        g.setColor(colorScheme[polygons[i].getFillColor()]);
                        g2d.fill(drawPoly);
                        g.setColor(colorScheme[polygons[i].getColor()]);
                        g2d.draw(drawPoly);
		            } else {
	                    g2d.fill(drawPoly);
		            }
		        } else {
		            g2d.draw(drawPoly);
		        }
		        if (polygons[i].getEditable()) {
    	            PathIterator path = polygons[i].myPath.getPathIterator(null);
    	            while (!path.isDone()) {
    	                double[] coords = new double[6];
    	                path.currentSegment(coords);
    	                g.fillRect(getXPixel(coords[0]) - 2, getYPixel(coords[1]) - 2, 4, 4);
    	                path.next();
    	            }
		        }
		    }
		}
		
		// Draw the rectangle outlined by dragging the mouse
		if (highlightRectangle != null) {
			xPixel = Math.min(getXPixel(highlightRectangle.getStartXValue()), getXPixel(highlightRectangle.getStopXValue()));
			yPixel = Math.min(getYPixel(highlightRectangle.getStartYValue()), getYPixel(highlightRectangle.getStopYValue()));
	    	widthPixel = Math.abs(getXPixel(highlightRectangle.getStartXValue()) - getXPixel(highlightRectangle.getStopXValue()));
	    	heightPixel = (Math.abs(getYPixel(highlightRectangle.getStartYValue()) - getYPixel(highlightRectangle.getStopYValue())));
	    	g.setColor(colorScheme[0]);
	    	drawRectThick(g, xPixel, yPixel, widthPixel, heightPixel, (byte) 1);
		}

//		if (highlightLine != null) {
//		    xPixel = Math.min(getXPixel(highlightLine.getStartX()), getXPixel(highlightLine.getStopX()));
//		    yPixel = Math.min(getYPixel(highlightLine.getStartY()), getYPixel(highlightLine.getStopY()));
//		    widthPixel = Math.abs(getXPixel(highlightLine.getStartX()) - getXPixel(highlightLine.getStopX()));
//		    heightPixel = (Math.abs(getYPixel(highlightLine.getStartY()) - getYPixel(highlightLine.getStopY())));
//		    g.setColor(colorScheme[0]);
//		    drawRectThick(g, xPixel, yPixel, widthPixel, heightPixel, (byte) 1);
//		}
//
		if (highlightPoly != null) {
		    g.setColor(colorScheme[0]);
		    ((Graphics2D) g).draw(pt.transform(highlightPoly.myPath));

            PathIterator path = highlightPoly.myPath.getPathIterator(null);
            while (!path.isDone()) {
                double[] coords = new double[6];
                path.currentSegment(coords);
                g.fillRect(getXPixel(coords[0]) - 2, getYPixel(coords[1]) - 2, 4, 4);
                path.next();
            }
		}

		g.setClip(null);
		
		
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
		
		if (base && displayTitle && title != null && !title.equals("")) {
		    int titleX, titleY;
		    int PAD = 5;
		    int fontHeight = (fontMetrics == null ? 25 : fontMetrics.getHeight());
		    int titleWidth = (fontMetrics == null ? title.length() * PAD : fontMetrics.stringWidth(title));
		    
	        switch (titleLocation) {
	            default: // DEFAULT TO TOP
	            case SwingConstants.NORTH:
	                titleX = (getWidth()-axisYWidth) / 2 + axisYWidth - titleWidth / 2;
	                titleY = PAD + fontHeight;
	                break;
	            case SwingConstants.NORTH_EAST:
	                titleX = getWidth() - 2 * PAD - titleWidth;
	                titleY = PAD + fontHeight;
	                break;
	            case SwingConstants.NORTH_WEST:
	                titleX = axisYWidth + 2 * PAD;
	                titleY = PAD + fontHeight;
	                break;
	        }
	        Color currColor = g.getColor();
	        g.setColor(Color.black);
            g.drawString(title, titleX, titleY);
            g.setColor(currColor);
		}
		
		setImageStatus(IMAGE_COMPLETE);

		refreshOtherComponents();

		if (DEBUGGING) {
			System.out.println("Took " + ext.getTimeElapsed(fullTime)+" to draw "+(createLookup?"(and create lookup for) ":"")+points.length+" points");
		}
	}


    private void drawYAxis(Graphics g, double[] plotMinMaxStep) {
        int sigFigs;
        String str;
        BufferedImage yLabel;
        Graphics gfx;
        FontMetrics fontMetrics;

        sigFigs = ext.getNumSigFig(plotMinMaxStep[2]);
    	float minSize = 14;// g.getFont().getSize2D();
    	Font prevFont = g.getFont();
    	Font minFont = prevFont.deriveFont(minSize);
    	g.setFont(minFont);
    	fontMetrics = g.getFontMetrics();

    	if (getYAxis() == AXIS_SCALE.LIN) {
	    	for (double y = plotMinMaxStep[3]; y <= plotYmax; y += plotMinMaxStep[2]) {
	    	    if (y >= plotYmin || !truncate) {
	    	        Grafik.drawThickLine(g, canvasSectionMaximumX-TICK_LENGTH, getYPixel(y), canvasSectionMaximumX, getYPixel(y), TICK_THICKNESS, Color.BLACK);
	    	        str = ext.formDeci(Math.abs(y) < DOUBLE_INACCURACY_HEDGE ? 0 : y, sigFigs, true);
	    	        g.drawString(str, canvasSectionMaximumX - TICK_LENGTH - 5 - fontMetrics.stringWidth(str), getYPixel(y) + fontMetrics.getHeight() / 2);
	    	    }
	    	}
    	} else if (getYAxis() == AXIS_SCALE.LOG) {
            double y;
            str = "-10";
            int strWid = fontMetrics.stringWidth(str);
            if (plotYmin < 0) {
                for (double exp = 1; (y = -1 * Math.pow(10, exp)) >= plotYmin; exp++) {
                    if (y <= plotYmax) {
                        Grafik.drawThickLine(g, canvasSectionMaximumX - TICK_LENGTH, getYPixel(y), canvasSectionMaximumX, getYPixel(y), TICK_THICKNESS, Color.BLACK);
                        g.drawString(str, canvasSectionMaximumX - TICK_LENGTH - 13 - strWid, getYPixel(y) + fontMetrics.getHeight() / 2);
                        g.drawString("" + ((int) exp), canvasSectionMaximumX - TICK_LENGTH - strWid + 9, getYPixel(y) + fontMetrics.getHeight() / 2 - 5);
                        for (double yPow = y; yPow > -1 * Math.pow(10, exp+1) && yPow >= plotYmin; yPow += y) {
                            Grafik.drawThickLine(g, canvasSectionMaximumX - (TICK_LENGTH/3 *2), getYPixel(yPow), canvasSectionMaximumX, getYPixel(yPow), TICK_THICKNESS, Color.BLACK);
                        }
                    }
                }
            }
            if (plotYmin <= 0 && plotYmax >= 0) {
                y = 0;
                Grafik.drawThickLine(g, canvasSectionMaximumX - TICK_LENGTH, getYPixel(y), canvasSectionMaximumX, getYPixel(y), TICK_THICKNESS, Color.BLACK);
                g.drawString("0", canvasSectionMaximumX - TICK_LENGTH - 5 - fontMetrics.stringWidth("0"), getYPixel(y) + fontMetrics.getHeight() / 2);
            }
            str = "10";
            strWid = fontMetrics.stringWidth(str);
            if (plotXmax > 0) {
                for (double exp = 1; (y = Math.pow(10, exp)) <= plotYmax; exp++) {
                    int yPix = getYPixel(y);
                    if (y >= plotYmin) { 
                        Grafik.drawThickLine(g, canvasSectionMaximumX - TICK_LENGTH, yPix, canvasSectionMaximumX, yPix, TICK_THICKNESS, Color.BLACK);
                        g.drawString(str, canvasSectionMaximumX - TICK_LENGTH - 13 - strWid, yPix + fontMetrics.getHeight() / 2);
                        g.drawString("" + ((int) exp), canvasSectionMaximumX - TICK_LENGTH - strWid + 4, yPix + fontMetrics.getHeight() / 2 - 5);
                        for (double yPow = y; yPow < Math.pow(10, exp+1) && yPow <= plotYmax; yPow += y) {
                            Grafik.drawThickLine(g, canvasSectionMaximumX - (TICK_LENGTH/3 *2), getYPixel(yPow), canvasSectionMaximumX, getYPixel(yPow), TICK_THICKNESS, Color.BLACK);
                        }
                    }
                }
            }
    	} else {
            Logicle fl = getBiexScale(false);
            double[] lbls = fl.axisLabels();
            int lin = (int) Math.ceil(fl.W) + 1;
            int maxLin = (int) Math.pow(10, lin);
            
            for (int i = 0; i < lbls.length; i++) {
                int yPix = getYPixel(lbls[i]);
                int exp = getExp((int)Math.abs(lbls[i]));
                if (lbls[i] < 0) {
                    str = "-10";
                    int strWid = fontMetrics.stringWidth(str);
                    if (yPix < canvasSectionMaximumY - axisXHeight + HEAD_BUFFER) {
                        Grafik.drawThickLine(g, canvasSectionMaximumX - TICK_LENGTH, yPix, canvasSectionMaximumX, yPix, TICK_THICKNESS, Color.BLACK);
                        g.drawString(str, canvasSectionMaximumX - TICK_LENGTH - 13 - strWid, yPix + fontMetrics.getHeight() / 2);
                        g.drawString("" + ((int) exp), canvasSectionMaximumX - TICK_LENGTH - strWid + 9, yPix + fontMetrics.getHeight() / 2 - 5);
                    }
                    for (double yPow = -Math.pow(10, exp); yPow > -1 * Math.pow(10, exp+1) && yPow >= plotYmin; yPow -= Math.pow(10, exp)) {
                        if (getYPixel(yPow) < canvasSectionMaximumY - axisXHeight + HEAD_BUFFER) {
                            Grafik.drawThickLine(g, canvasSectionMaximumX - (TICK_LENGTH/3 *2), getYPixel(yPow), canvasSectionMaximumX, getYPixel(yPow), TICK_THICKNESS, Color.BLACK);
                        }
                    }
                } else if (lbls[i] == 0) {
                    str = "0";
                    Grafik.drawThickLine(g, canvasSectionMaximumX - TICK_LENGTH, yPix, canvasSectionMaximumX, yPix, TICK_THICKNESS, Color.BLACK);
                    g.drawString("0", canvasSectionMaximumX - TICK_LENGTH - 5 - fontMetrics.stringWidth("0"), yPix + fontMetrics.getHeight() / 2);
                    if (i > 0) {
                        int next = getExp((int) lbls[i + 1]);
                        for (int e = 0; e < next; e++) {
                            for (double yPow = -Math.pow(10, e); yPow > -Math.pow(10, e+1) && yPow <= plotXmax; yPow -= Math.pow(10, e)) {
                                if (getYPixel(yPow) < canvasSectionMaximumY - axisXHeight + HEAD_BUFFER) {
                                    Grafik.drawThickLine(g, canvasSectionMaximumX - TICK_LENGTH, getYPixel(yPow), canvasSectionMaximumX, getYPixel(yPow), TICK_THICKNESS, Color.BLACK);
                                }
                            }
                        }
                    }
                    if (i < lbls.length - 1) {
                        int next = getExp((int) lbls[i + 1]);
                        for (int e = 0; e < next; e++) {
                            for (double yPow = Math.pow(10, e); yPow < Math.pow(10, e+1) && yPow <= plotXmax; yPow += Math.pow(10, e)) {
                                Grafik.drawThickLine(g, canvasSectionMaximumX - TICK_LENGTH, getYPixel(yPow), canvasSectionMaximumX, getYPixel(yPow), TICK_THICKNESS, Color.BLACK);
                            }
                        }
                    }
                } else {
                    str = "10";
                    int strWid = fontMetrics.stringWidth(str);
                    Grafik.drawThickLine(g, canvasSectionMaximumX - TICK_LENGTH, yPix, canvasSectionMaximumX, yPix, TICK_THICKNESS, Color.BLACK);
                    g.drawString(str, canvasSectionMaximumX - TICK_LENGTH - 13 - strWid, yPix + fontMetrics.getHeight() / 2);
                    g.drawString("" + ((int) exp), canvasSectionMaximumX - TICK_LENGTH - strWid + 4, yPix + fontMetrics.getHeight() / 2 - 5);
                    for (double yPow = Math.pow(10, exp); yPow < Math.pow(10, exp+1) && yPow <= plotYmax; yPow += Math.pow(10, exp)) {
                        Grafik.drawThickLine(g, canvasSectionMaximumX - (TICK_LENGTH/3 *2), getYPixel(yPow), canvasSectionMaximumX, getYPixel(yPow), TICK_THICKNESS, Color.BLACK);
                    }
                }
            }
            if (lbls[0] < -maxLin) {
                str = "-10";
                int strWid = fontMetrics.stringWidth(str);
                int yPix = getYPixel(-maxLin);
                Grafik.drawThickLine(g, canvasSectionMaximumX - TICK_LENGTH, yPix, canvasSectionMaximumX, yPix, TICK_THICKNESS, Color.BLACK);
                g.drawString(str, canvasSectionMaximumX - TICK_LENGTH - 13 - strWid, yPix + fontMetrics.getHeight() / 2);
                g.drawString("" + (int) lin, canvasSectionMaximumX - TICK_LENGTH - strWid + 9, yPix + fontMetrics.getHeight() / 2 - 5);
            }
            if (plotXmax > maxLin) {
                str = "10";
                int strWid = fontMetrics.stringWidth(str);
                int yPix = getYPixel(maxLin);
                Grafik.drawThickLine(g, canvasSectionMaximumX - TICK_LENGTH, yPix, canvasSectionMaximumX, yPix, TICK_THICKNESS, Color.BLACK);
                g.drawString(str, canvasSectionMaximumX - TICK_LENGTH - 13 - strWid, yPix + fontMetrics.getHeight() / 2);
                g.drawString("" + (int) lin, canvasSectionMaximumX - TICK_LENGTH - strWid + 4, yPix + fontMetrics.getHeight() / 2 - 5);
            }
    	}
    	
    	g.setFont(prevFont);
    	fontMetrics = g.getFontMetrics();
    	Grafik.drawThickLine(g, canvasSectionMaximumX, getHeight() - canvasSectionMaximumY, canvasSectionMaximumX, canvasSectionMaximumY - axisXHeight + HEAD_BUFFER, AXIS_THICKNESS, Color.BLACK);
    	int strWidth = fontMetrics.stringWidth(yAxisLabel);
    	if (strWidth > 0) {
    		yLabel = new BufferedImage(strWidth, 36, BufferedImage.TYPE_INT_RGB);
    		gfx = yLabel.createGraphics();
    		gfx.setFont(new Font("Arial", 0, axisFontSize));
    		gfx.setColor(Color.WHITE);
    		gfx.fillRect(0, 0, getWidth(), getHeight());
    		gfx.setColor(Color.BLACK);
    		gfx.drawString(yAxisLabel, 0, yLabel.getHeight()-6);
    		
    		g.drawImage(Grafik.rotateImage(yLabel, true), 10, (getHeight()-axisXHeight/*HEIGHT_X_AXIS*/)/2-fontMetrics.stringWidth(yAxisLabel)/2, this);
        }
    }


    private void drawXAxis(Graphics g, double[] plotMinMaxStep, FontMetrics fontMetrics) {
        int sigFigs;
        String str;
        
    	float minSize = 14;
    	Font prevFont = g.getFont();
        Font minFont = prevFont.deriveFont(minSize);
        g.setFont(minFont);
        fontMetrics = g.getFontMetrics();
    	
        int fontHgt = fontMetrics.getHeight();
        
        if (getXAxis() == AXIS_SCALE.LIN) {
	    	sigFigs = ext.getNumSigFig(plotMinMaxStep[2]);
	    	for (double x = plotMinMaxStep[3]; x<=plotXmax; x += plotMinMaxStep[2]) {
	    		if (x >= plotXmin || !truncate) {
	    			Grafik.drawThickLine(g, getXPixel(x), getHeight()-canvasSectionMaximumY, getXPixel(x), getHeight()-(canvasSectionMaximumY-TICK_LENGTH), TICK_THICKNESS, Color.BLACK);
	    			str = ext.formDeci(Math.abs(x) < DOUBLE_INACCURACY_HEDGE ? 0 : x, sigFigs, true);
	    			g.drawString(str, getXPixel(x)-str.length()*4, getHeight()-(canvasSectionMaximumY-TICK_LENGTH-fontHgt));
	    		}
	    	}
        } else if (getXAxis() == AXIS_SCALE.LOG){
            double x;
            str = "-10";
            int strWid = fontMetrics.stringWidth(str);
            if (plotXmin < 0) {
                for (double exp = 1; (x = -1 * Math.pow(10, exp)) >= plotXmin; exp++) {
                    if (x <= plotXmax) {
                        Grafik.drawThickLine(g, getXPixel(x), getHeight()-canvasSectionMaximumY, getXPixel(x), getHeight()-(canvasSectionMaximumY-TICK_LENGTH), TICK_THICKNESS, Color.BLACK);
                        g.drawString(str, getXPixel(x) - strWid / 2, getHeight()-(canvasSectionMaximumY-TICK_LENGTH-fontHgt));
                        g.drawString("" + ((int) exp), getXPixel(x) + (strWid / 2) + 1, getHeight()-(canvasSectionMaximumY-TICK_LENGTH-fontHgt + 5));
                        for (double xPow = x; xPow > -1 * Math.pow(10, exp+1) && xPow >= plotXmin; xPow += x) {
                            Grafik.drawThickLine(g, getXPixel(xPow), getHeight()-canvasSectionMaximumY, getXPixel(xPow), getHeight()-(canvasSectionMaximumY-(TICK_LENGTH/3 *2)), TICK_THICKNESS, Color.BLACK);
                        }
                    }
                }
            }
            if (plotXmin <= 0 && plotXmax >= 0) {
                x = 0;
                Grafik.drawThickLine(g, getXPixel(x), getHeight()-canvasSectionMaximumY, getXPixel(x), getHeight()-(canvasSectionMaximumY-TICK_LENGTH), TICK_THICKNESS, Color.BLACK);
                g.drawString("0", getXPixel(x)-fontMetrics.stringWidth("0")/2, getHeight()-(canvasSectionMaximumY-TICK_LENGTH-fontHgt));
            }
            str = "10";
            strWid = fontMetrics.stringWidth(str);
            if (plotXmax > 0) {
                for (double exp = 1; (x = Math.pow(10, exp)) <= plotXmax; exp++) {
                    int xPix = getXPixel(x);
                    if (x >= plotXmin) { 
                        Grafik.drawThickLine(g, xPix, getHeight()-canvasSectionMaximumY, xPix, getHeight()-(canvasSectionMaximumY-TICK_LENGTH), TICK_THICKNESS, Color.BLACK);
                        g.drawString(str, xPix - strWid / 2, getHeight()-(canvasSectionMaximumY-TICK_LENGTH-fontHgt));
                        g.drawString("" + (int) exp, xPix + (strWid / 2) + 1, getHeight()-(canvasSectionMaximumY-TICK_LENGTH-fontHgt + 5));
                        for (double xPow = x; xPow < Math.pow(10, exp+1) && xPow <= plotXmax; xPow += x) {
                            Grafik.drawThickLine(g, getXPixel(xPow), getHeight()-canvasSectionMaximumY, getXPixel(xPow), getHeight()-(canvasSectionMaximumY-(TICK_LENGTH/3 *2)), TICK_THICKNESS, Color.BLACK);
                        }
                    }
                }
            }
        } else {
            Logicle fl = getBiexScale(true);
            double[] lbls = fl.axisLabels();
            int lin = (int) Math.ceil(fl.W) + 1;
            int maxLin = (int) Math.pow(10, lin);
            
            for (int i = 0; i < lbls.length; i++) {
                int xPix = getXPixel(lbls[i]);
                int exp = getExp((int)Math.abs(lbls[i]));
                if (lbls[i] < 0) {
                    str = "-10";
                    int strWid = fontMetrics.stringWidth(str);
                    Grafik.drawThickLine(g, xPix, getHeight()-canvasSectionMaximumY, xPix, getHeight()-(canvasSectionMaximumY-TICK_LENGTH), TICK_THICKNESS, Color.BLACK);
                    g.drawString(str, xPix - strWid / 2, getHeight()-(canvasSectionMaximumY-TICK_LENGTH-fontHgt));
                    g.drawString("" + (int) exp, xPix + (strWid / 2) + 1, getHeight()-(canvasSectionMaximumY-TICK_LENGTH-fontHgt + 5));
                    for (double xPow = -Math.pow(10, exp); xPow > -1 * Math.pow(10, exp+1) && xPow >= plotXmin; xPow -= Math.pow(10, exp)) {
                        Grafik.drawThickLine(g, getXPixel(xPow), getHeight()-canvasSectionMaximumY, getXPixel(xPow), getHeight()-(canvasSectionMaximumY-(TICK_LENGTH/3 *2)), TICK_THICKNESS, Color.BLACK);
                    }
                } else if (lbls[i] == 0) {
                    str = "0";
                    int strWid = fontMetrics.stringWidth(str);
                    Grafik.drawThickLine(g, xPix, getHeight()-canvasSectionMaximumY, xPix, getHeight()-(canvasSectionMaximumY-TICK_LENGTH), TICK_THICKNESS, Color.BLACK);
                    g.drawString(str, xPix - strWid / 2, getHeight()-(canvasSectionMaximumY-TICK_LENGTH-fontHgt));
                    if (i > 0) {
                        int next = getExp((int) lbls[i + 1]);
                        for (int e = 0; e < next; e++) {
                            for (double xPow = -Math.pow(10, e); xPow > -Math.pow(10, e+1) && xPow <= plotXmax; xPow -= Math.pow(10, e)) {
                                Grafik.drawThickLine(g, getXPixel(xPow), getHeight()-canvasSectionMaximumY, getXPixel(xPow), getHeight()-(canvasSectionMaximumY-(TICK_LENGTH/3 *2)), TICK_THICKNESS, Color.BLACK);
                            }
                        }
                    }
                    if (i < lbls.length - 1) {
                        int next = getExp((int) lbls[i + 1]);
                        for (int e = 0; e < next; e++) {
                            for (double xPow = Math.pow(10, e); xPow < Math.pow(10, e+1) && xPow <= plotXmax; xPow += Math.pow(10, e)) {
                                Grafik.drawThickLine(g, getXPixel(xPow), getHeight()-canvasSectionMaximumY, getXPixel(xPow), getHeight()-(canvasSectionMaximumY-(TICK_LENGTH/3 *2)), TICK_THICKNESS, Color.BLACK);
                            }
                        }
                    }
                } else {
                    str = "10";
                    int strWid = fontMetrics.stringWidth(str);
                    Grafik.drawThickLine(g, xPix, getHeight()-canvasSectionMaximumY, xPix, getHeight()-(canvasSectionMaximumY-TICK_LENGTH), TICK_THICKNESS, Color.BLACK);
                    g.drawString(str, xPix - strWid / 2, getHeight()-(canvasSectionMaximumY-TICK_LENGTH-fontHgt));
                    g.drawString("" + (int) exp, xPix + (strWid / 2) + 1, getHeight()-(canvasSectionMaximumY-TICK_LENGTH-fontHgt + 5));
                    for (double xPow = Math.pow(10, exp); xPow < Math.pow(10, exp+1) && xPow <= plotXmax; xPow += Math.pow(10, exp)) {
                        Grafik.drawThickLine(g, getXPixel(xPow), getHeight()-canvasSectionMaximumY, getXPixel(xPow), getHeight()-(canvasSectionMaximumY-(TICK_LENGTH/3 *2)), TICK_THICKNESS, Color.BLACK);
                    }
                }
            }
            if (lbls[0] < -maxLin) {
                str = "-10";
                int strWid = fontMetrics.stringWidth(str);
                int xPix = getXPixel(-maxLin);
                Grafik.drawThickLine(g, xPix, getHeight()-canvasSectionMaximumY, xPix, getHeight()-(canvasSectionMaximumY-TICK_LENGTH), TICK_THICKNESS, Color.BLACK);
                g.drawString(str, xPix - strWid / 2, getHeight()-(canvasSectionMaximumY-TICK_LENGTH-fontHgt));
                g.drawString("" + (int) lin, xPix + (strWid / 2) + 1, getHeight()-(canvasSectionMaximumY-TICK_LENGTH-fontHgt + 5));
            }
            if (plotXmax > maxLin) {
                str = "10";
                int strWid = fontMetrics.stringWidth(str);
                int xPix = getXPixel(maxLin);
                Grafik.drawThickLine(g, xPix, getHeight()-canvasSectionMaximumY, xPix, getHeight()-(canvasSectionMaximumY-TICK_LENGTH), TICK_THICKNESS, Color.BLACK);
                g.drawString(str, xPix - strWid / 2, getHeight()-(canvasSectionMaximumY-TICK_LENGTH-fontHgt));
                g.drawString("" + (int) lin, xPix + (strWid / 2) + 1, getHeight()-(canvasSectionMaximumY-TICK_LENGTH-fontHgt + 5));
            }
        }
        g.setFont(prevFont);
        fontMetrics = g.getFontMetrics();
        
    }
    
    private int getExp(int l) {
        int e = 0;
        double pow = Math.pow(10, e);
        while ((l / pow) > 1) {
            e++;
            pow = Math.pow(10, e);
        }
        return e;
    }
	
	/**
	 * Currently works with:<br />
	 * <ul>
	 * <li>SwingConstants.NORTH</li>
	 * <li>SwingConstants.NORTH_EAST</li>
	 * <li>SwingConstants.NORTH_WEST</li>
	 * </ul>
	 * @param loc
	 */
	public void setTitleLocation(int loc) {
	    this.titleLocation = loc;
	}

	public void setTitle(String title) {
	    this.title = title;
	}
	
	public void setDisplayTitle(boolean show) {
	    this.displayTitle = show;
	}
	
	public boolean getDisplayTitle() {
	    return this.displayTitle;
	}
	
	public String getTitle() { return this.title; }
	
	public void refreshOtherComponents() {}

	public void drawHeatMap(Graphics g) {
		int nRows, nColumns;
		int[][] gridIntensities;
		int[][][] gridColors;
		
		int width = 1, height = 1, radius = 2;
		
		nRows = getHeight();
		nColumns = getWidth();
		gridIntensities = getGridIntensityForHeapMapGrid(nRows, nColumns, width, height, radius);
		gridColors = getColorFromIntensityForHeapMapGrid(gridIntensities);
		for (int i = 0; i < nColumns; i++) {
			for (int j = 0; j < nRows; j++) {
				if (gridIntensities[i][j] != 0) {
					g.setColor(new Color(gridColors[i][j][0], gridColors[i][j][1], gridColors[i][j][2]));
					g.drawRect(i + canvasSectionMinimumX, j - canvasSectionMinimumY, 1, 1);
//					g.fillRect(i + canvasSectionMinimumX, j - canvasSectionMinimumY, 2, 2);
				}
			}
		}
	}


	public int[][] getGridIntensityForHeapMapGrid(int nRows, int nColumns, int cellWidth, int cellHeight, int neighbor) {
		int xPixel, yPixel;
		int[][] intensities;
//		boolean zoomedIn;
//		int[] origin;

//		origin = new int[] {0,0};

//		zoomedIn = (Math.abs(getXValueFromXPixel(5) - getXValueFromXPixel(0)) < 0.002) || (Math.abs(getYValueFromYPixel(5) - getYValueFromYPixel(0)) < 0.002);
		intensities = new int[nColumns][nRows];
		for (int i = 0; i < points.length; i++) {
			if (points[i] != null) {
				xPixel = getXPixel(points[i].getRawX()) - canvasSectionMinimumX;
				yPixel = getYPixel(points[i].getRawY()) + canvasSectionMinimumY;
				for (int j = -neighbor; j <= neighbor; j++) {
					for (int k = -neighbor; k <= neighbor; k++) {
						if ((xPixel + j) >= 0 && (xPixel + j) < nColumns && (yPixel + k) >= 0 && (yPixel + k) < nRows) { //  && Distance.euclidean(new int[] {Math.abs(j), Math.abs(k)}, origin) < neighbor*(neighbor-1)
//							if (zoomedIn) {
								intensities[xPixel + j][yPixel + k]++;
//							} else {
//								intensities[xPixel + j][yPixel + k] += neighbor*neighbor;// - Distance.euclidean(new int[] {Math.abs(j), Math.abs(k)}, origin);
//							}
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
					color[i][j] = Grafik.getHeatmapColor(Math.min(1, ((double)intensities[i][j] / (double)max) * 1.4));
				}
			}
		}

		return color;
	}

	public void mouseClicked(MouseEvent e) {}

	public void mouseEntered(MouseEvent e) {}

	public void mouseExited(MouseEvent e) {}
	
	public void mouseMoved(MouseEvent event) {
		Graphics g = getGraphics();
		IntVector indicesOfNearbyPoints;
		String pos;
		int x, y, dataPointIndex;
		byte size;

//		if (imageIsFinal() && chartType != PLOT_TYPE.HEATMAP) {
//			x = event.getX();
//			y = event.getY();
//
//			canvasSectionMinimumX = axisYWidth;//WIDTH_Y_AXIS;
//			canvasSectionMaximumX = getWidth() - WIDTH_BUFFER;
//			canvasSectionMinimumY = axisXHeight;//HEIGHT_X_AXIS;
//			canvasSectionMaximumY = getHeight() - HEAD_BUFFER;
//			pos = (int)Math.floor(x / DEFAULT_LOOKUP_RESOLUTION) + "x" + (int)Math.floor(y / DEFAULT_LOOKUP_RESOLUTION);
////			if (!pos.equals(prevPos)) {
////			    System.out.println("Repainting");
////				repaint();
////			}
//
//			indicesOfNearbyPoints = lookupNearbyPoints(x, y, pos);
//			
//			if (!pos.equals(prevPos)) {
//			    if (indicesOfNearbyPoints.size() != 0 || (prox != null && prox.size() != 0)) {
//                    repaint();
//                }
//			}
//
//			prox = new IntVector();
//
//			size = SIZE * 2;
//			g.setColor(Color.GRAY);
//			for (int i = 0; indicesOfNearbyPoints!=null && i<indicesOfNearbyPoints.size(); i++) {
//				dataPointIndex = indicesOfNearbyPoints.elementAt(i);
//				if (Distance.euclidean(new int[] {x, y}, new int[] {getXPixel(points[dataPointIndex].getRawX()), getYPixel(points[dataPointIndex].getRawY())}) < HIGHLIGHT_DISTANCE) {
//					prox.add(dataPointIndex);
////					g.setColor(Color.YELLOW);
////					g.fillOval(getX(points[dataPointIndex].getRawX()) - size/2, getY(points[dataPointIndex].getRawY()) - size/2, size, size);
//					g.drawOval(getXPixel(points[dataPointIndex].getRawX()) - size/2, getYPixel(points[dataPointIndex].getRawY()) - size/2, size, size);
//
//				}
//			}
//
//			prevPos = pos;
//			
//			//TODO 		indicesOfNearbySamples = lookupNearbyPoints(x, y, pos);
//		}
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
		double distance = -1;
		
		
		check : {
		    // if we're dragging when bounds are forced, set current bounds to forced bounds and remove forced bounds
		    if (!Float.isNaN(forcePlotXmax)) {
		        setPlotXMax(forcePlotXmax);
		    }
		    if (!Float.isNaN(forcePlotXmin)) {
		        setPlotXMin(forcePlotXmin);
		    }
		    if (!Float.isNaN(forcePlotYmax)) {
		        setPlotYMax(forcePlotYmax);
		    }
		    if (!Float.isNaN(forcePlotYmin)) {
		        setPlotYMin(forcePlotYmin);
		    }
            setForcePlotXMin(Float.NaN);
            setForcePlotXMax(Float.NaN);
            setForcePlotYMin(Float.NaN);
            setForcePlotYMax(Float.NaN);
		}
		
//		System.err.println("AbstractPanel mouseDragged has been called");

		curX = e.getPoint().x;
		curY = e.getPoint().y;
//		if (invertX) {
//			curX = (2 * startX) - curX;
//		}
//		if (invertY) {
//			curY = (2 * startY) - curY;
//		}
		
		for (int i = 0; i < 2; i++) {
			if (i==0) {
				distance = (startX-curX) * (zoomSubsets[0][1]-zoomSubsets[0][0])/(getWidth()-WIDTH_BUFFER-axisYWidth/*WIDTH_Y_AXIS*/);
			} else {
				distance = (curY-startY) * (zoomSubsets[1][1]-zoomSubsets[1][0])/(getHeight()-HEAD_BUFFER-axisXHeight/*HEIGHT_X_AXIS*/);
			}
			
			if (distance < 0) {
				distance = Math.max(distance, -1*zoomSubsets[i][0]);
			} else {
				distance = Math.min(distance, 1-zoomSubsets[i][1]);
			}
//			if ((zoomSubsets[i][0] <= 0 && distance < 0) || (zoomSubsets[i][1] >= 1 && distance > 0)) {
//
//			} else {
//				if ((invertX && i == 0) || (invertY && i == 1)) {
//					zoomSubsets[i][0] -= distance;
//					zoomSubsets[i][1] -= distance;
//				} else {
//					zoomSubsets[i][0] += distance;
//					zoomSubsets[i][1] += distance;
//				}
				zoomSubsets[i][0] += distance;
				zoomSubsets[i][1] += distance;
//			}
		}
		
		if (inDrag) {
//		    if (distance > 0.00001 || distance < -0.00001) {
		        paintAgain();
//		    }
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
        int w = (int) e.getComponent().getSize().getWidth();
        int h = (int) e.getComponent().getSize().getHeight();
        if (trackedW == -1 || trackedH == -1) {
            setTrackedSize(w, h);
            return;
        }
        if (w != trackedW || h != trackedH) {
            setTrackedSize(w, h);
        }
	    
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
	
    private volatile int trackedW = -1;
    private volatile int trackedH = -1; 
    public void setTrackedSize(int w, int h) {
        this.trackedW = w;
        this.trackedH = h;
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
			if (DEBUGGING) {
				System.err.println("Action performed in AbstractPanel");
			}
			paintAgain();
		}
		
	}

	public void zoomProportionally(boolean outNotIn, Point p, boolean center) {
		float x, y, multiplier, dist;
		float[][] proportions;
		int width, height;
		boolean changed;
		
		proportions = new float[2][2];
		
		width = getWidth()-/*WIDTH_Y_AXIS*/axisYWidth-WIDTH_BUFFER;
		x = (float)p.getX()-axisYWidth;//WIDTH_Y_AXIS;
		proportions[0][0] = x/width;
		proportions[0][1] = (width-x)/width;

		height = getHeight()-/*HEIGHT_X_AXIS*/axisXHeight-HEAD_BUFFER;
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
	
	public String getNumPointsPlottedEfficiently() {
		if (beEfficient) {
			return pointsPlotted.size()+"";
		} else {
			return "[not set to be efficient]";
		}
	}
	
	public void drawPoint(Graphics g, PlotPoint point) {
		int x, y, size, color;
		
		x = getXPixel(point.getRawX());
		y = getYPixel(point.getRawY());
		size = point.getSize();
		color = point.getColor();
		

		if (size == 0) {
			return;
		}

		if (beEfficient) {
			if (pointsPlotted.contains(x+":"+y+":"+size+":"+color)) {
				return;
			} else {
				pointsPlotted.add(x+":"+y+":"+size+":"+color);
			}
//	    	Files.appendStringToFile("listOfPoints.out", x+":"+y+":"+size+":"+color);
		}
		g.setColor(colorScheme[color]);
		
		switch (point.getType()) {
			case PlotPoint.FILLED_CIRCLE:
				g.fillOval(x-size/2, y-size/2, size, size);
				break;  
			case PlotPoint.OPEN_CIRCLE:
				g.drawOval(x-size/2, y-size/2, size, size);
				break;
			case PlotPoint.FILLED_SQUARE:
				g.fillPolygon(new int[] {x-size/2, x-size/2, x+size/2, x+size/2},
							  new int[] {y-size/2, y+size/2, y+size/2, y-size/2},
							  4);
				break;
			case PlotPoint.OPEN_SQUARE:
				g.drawPolygon(new int[] {x-size/2, x-size/2, x+size/2, x+size/2},
							  new int[] {y-size/2, y+size/2, y+size/2, y-size/2},
							  4);
				break;
			case PlotPoint.FILLED_TRIANGLE:
				Grafik.drawTriangle(g, x, y, size, true);
	//			g.drawPolygon(new int[] {x-size/2, x, +size/2},
	//						  new int[] {y+size/2, y-size/2, y+size/2},
	//						  3);
				break;
			case PlotPoint.MISSING:
				//g.drawString(PlotPoint.MISSING_STR, getX(point.getRawX())-missingWidth/2, getY(point.getRawY())+size/2);
				if (PlotPoint.MISSING_STR.equals("X") || PlotPoint.MISSING_STR.equals("x")){
					g.drawLine(x-size/4, y-size/4, x+size/4, y+size/4);
					g.drawLine(x-size/4, y+size/4, x+size/4, y-size/4);
				} else {
					g.setFont(new Font("Arial", 0, size));
					g.drawString(PlotPoint.MISSING_STR, -size/2, y+size/2);
					g.setFont(new Font("Arial", 0, axisFontSize));
				}
				break;
			case PlotPoint.NOT_A_NUMBER:
	//			g.drawString(PlotPoint.NAN_STR, getX(point.getRawX())-nanWidth/2, getY(point.getRawY())-30+size/2);
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

	public double[] getPlotMinMaxStep(double min, double max, Graphics g, boolean xAxis) {
		double range, plotStep, stepStep, plotMin, plotMax;
		double zoomMin, zoomMax, dist, tempD;
		int numHashes, wid, hgt, canvasRange;
		FontMetrics fontMetrics;
		int sf, temp;

		range = max-min;

		plotStep = stepStep = calcStepStep(range);
		sf = ext.getNumSigFig(stepStep);

		fontMetrics = g.getFontMetrics(g.getFont());
		if (xAxis) {
			wid = Math.max(fontMetrics.stringWidth(ext.formDeci(min, sf)), fontMetrics.stringWidth(ext.formDeci(max, sf)));
			numHashes = 12;
            canvasRange = canvasSectionMaximumX-canvasSectionMinimumX;
            temp = (wid + 30) * numHashes;
            while (temp > canvasRange) {
                numHashes -= 1;
                temp = (wid + 30) * numHashes;
            }
		} else {
		    hgt = fontMetrics.getHeight();
		    numHashes = 10;
		    canvasRange = canvasSectionMaximumY - canvasSectionMinimumY;
		    temp = (hgt + 10) * numHashes;
		    while (temp > canvasRange) {
		        numHashes -= 1;
		        temp = (hgt + 10) * numHashes;
		    }
		}
		numHashes = Math.max(numHashes, 1);
		while (range/plotStep>numHashes) {
			plotStep += stepStep;
		}
		
        plotMin = plotMax = 0;
        while (max - plotMax > DOUBLE_INACCURACY_HEDGE) {
            plotMax += plotStep;
        }
        
        if (min > plotMin) {
            tempD = min - (plotMin + 2 * plotStep);
            while (tempD > DOUBLE_INACCURACY_HEDGE) {
                plotMin += plotStep;
                tempD = min - (plotMin + 2 * plotStep);
            }
        } else {
            tempD = min - plotMin;
            while (tempD < -DOUBLE_INACCURACY_HEDGE) {
                plotMin -= plotStep;
                tempD = min - plotMin;
            }
        }
		
//		plotMin = plotMax = 0;
//		while (max - plotMax > DOUBLE_INACCURACY_HEDGE) {
//			plotMax += plotStep;
//		}
//		while (min - plotMin < -1*DOUBLE_INACCURACY_HEDGE) { // double check this, untested
//			plotMin -= plotStep;
//		}
////		if (min >= 0 && plotMin < 0) {
////			plotMin = 0;
////		}
//		
////		System.out.println(Float.parseFloat(ext.formDeci(plotMin, sf))+"\t"+Float.parseFloat(ext.formDeci(plotMax, sf))+"\t"+Float.parseFloat(ext.formDeci(plotStep, sf)));
//		
		if (zoomable) {
			dist = plotMax - plotMin;
			zoomMin = plotMin + zoomSubsets[xAxis?0:1][0] * dist;
			zoomMax = plotMax - (1-zoomSubsets[xAxis?0:1][1]) * dist;
			
			range = zoomMax-zoomMin;
			plotStep = stepStep = calcStepStep(range);
			sf = ext.getNumSigFig(stepStep);

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
			
			numHashes = Math.max(numHashes, 1);
			while (range/plotStep>numHashes) {
				plotStep += stepStep;
			}
//			System.out.println("2");
			return new double[] {zoomMin, zoomMax, Double.parseDouble(ext.formDeci(plotStep, sf)), Double.parseDouble(ext.formDeci(plotMin, sf))};
		} else {
			return new double[] {Double.parseDouble(ext.formDeci(plotMin, sf)), Double.parseDouble(ext.formDeci(plotMax, sf)), Double.parseDouble(ext.formDeci(plotStep, sf)), Double.parseDouble(ext.formDeci(plotMin, sf))};
		}
	}

	private static final double LOG10SCALE = 1/Math.log(10);
	// handy static methods
	private static double log10(double val) { return Math.log(val) * LOG10SCALE; }
	private static float flog10(double val) { return (float)log10(val); }


	private int getLogPixel(double val, boolean xAxis) {
	    int screenMax, screenMin;
	    double plotMax, plotMin;
	    
	    screenMin = xAxis ? canvasSectionMinimumX : canvasSectionMinimumY;
	    screenMax = xAxis ? canvasSectionMaximumX : canvasSectionMaximumY;
	    plotMin = xAxis ? plotXmin : plotYmin;
	    plotMax = xAxis ? plotXmax : plotYmax;
	    
	    float log_low = plotMin == 0 ? 0 : (plotMin < 0 ? -1 * flog10(plotMin * -1) : flog10(plotMin));
	    float log_high = plotMax == 0 ? 0 : (plotMax < 0 ? -1 * flog10(plotMax * -1) : flog10(plotMax));
	    float log_val = val == 0 ? 0 : (val < 0 ? -1 * flog10(val * -1) : flog10(val));
	    float pixels_per_log_unit = (screenMax - screenMin) / (log_high - log_low);
	    int retVal = (int)((log_val - log_low) * pixels_per_log_unit) + screenMin;
	    if (xAxis) return retVal;
	    else return getHeight() - retVal;
	}

	private int getLinPixel(double val, boolean xAxis) {
	    if (xAxis) {    
    	    if (invertX) {
                return (int)((plotXmax-val)/(plotXmax-plotXmin)*(double)(canvasSectionMaximumX-canvasSectionMinimumX))+canvasSectionMinimumX;
            } else {
                return (int)((val-plotXmin)/(plotXmax-plotXmin)*(double)(canvasSectionMaximumX-canvasSectionMinimumX))+canvasSectionMinimumX;
            }
	    } else {
            if (invertY) {
                return getHeight()-(int)((plotYmax-val)/(plotYmax-plotYmin)*(double)(canvasSectionMaximumY-canvasSectionMinimumY)+canvasSectionMinimumY);
            } else {
                return getHeight()-(int)((val-plotYmin)/(plotYmax-plotYmin)*(double)(canvasSectionMaximumY-canvasSectionMinimumY)+canvasSectionMinimumY);
            }
	    }
	}
	
	protected int getXPixel(double x) {
		if (getXAxis() == AXIS_SCALE.LIN) {
			return getLinPixel(x, true);
		} else if (getXAxis() == AXIS_SCALE.LOG){
			return getLogPixel(x, true);
		} else {
		    Logicle l = getBiexScale(true);
		    double bi = l.scale(x);
		    int dis = (int) (canvasSectionMaximumX * bi) - TICK_LENGTH - AXIS_THICKNESS;
		    return dis; // TODO
		}
	}
	
	protected Logicle biexScaleX = null;
	protected Logicle biexScaleY = null;
	private Logicle getBiexScale(boolean xAxis) {
	    double W = 2; // linear decades
	    double M = 5.5; // log decades (/total decades)
	    double A = -1.35; // negative decades
	    if (xAxis) {
    	    if (biexScaleX != null) return biexScaleX;
//    	    Logicle fl = new FastLogicle(plotXmax, W, M, A, 1); // FastLogicle causes issues (possible that parameters are incorrectly set)
            Logicle fl = new Logicle(plotXmax, W, M, A);
    	    biexScaleX = fl;
    	    return fl;
	    } else {
    	    if (biexScaleY != null) return biexScaleY;
//            Logicle fl = new FastLogicle(plotYmax, W, M, A, 1);
            Logicle fl = new Logicle(plotYmax, W, M, A);
            biexScaleY = fl;
            return fl;
	    }
	}
	
	protected int getYPixel(double y) {
		if (getYAxis() == AXIS_SCALE.LIN) {
		    return getLinPixel(y, false);
		} else if (getYAxis() == AXIS_SCALE.LOG) {
			return getLogPixel(y, false);
        } else {
            Logicle l = getBiexScale(false);
            double bi = l.scale(y);
            int dis = getHeight() - (int) (canvasSectionMaximumY * bi) - TICK_LENGTH;
            return dis; // TODO biex
		}
	}
    /**
     * Converts mouse location int X,Y into data points' value double rawX,rawY.
     * The control variable invertX is assigned elsewhere in the class.
     * @param mouseX the mouse location X
     * @return the rawX value of the corresponding data point.
     */
    public double getXValueFromXPixel(int mouseX) {
        if (getXAxis() == AXIS_SCALE.LIN) {
            return getLinValueFromPixel(mouseX);
        } else if (getXAxis() == AXIS_SCALE.LOG) {
            return getLogValueFromPixel(mouseX, true);
        } else {
            double x = (mouseX + TICK_LENGTH + AXIS_THICKNESS) / (double) canvasSectionMaximumX;
            Logicle l = getBiexScale(true);
            double xVal;
            xVal = l.inverse(x);
            return xVal;
        }
    }
    
    private double getLinValueFromPixel(int mouseX) {
        if (invertX) {
            return plotXmax - ((double)(mouseX-canvasSectionMinimumX)/(double)(canvasSectionMaximumX-canvasSectionMinimumX)*(double)(plotXmax-plotXmin));
        } else {
            return plotXmin + ((double)(mouseX-canvasSectionMinimumX)/(double)(canvasSectionMaximumX-canvasSectionMinimumX)*(double)(plotXmax-plotXmin));
        }
    }


    private double getLogValueFromPixel(int value, boolean xAxis) {
        int screenMax, screenMin;
        double plotMax, plotMin;
        
        screenMin = xAxis ? canvasSectionMinimumX : canvasSectionMinimumY;
        screenMax = xAxis ? canvasSectionMaximumX : canvasSectionMaximumY;
        plotMin = xAxis ? plotXmin : plotYmin;
        plotMax = xAxis ? plotXmax : plotYmax;
        
        float log_low = plotMin == 0 ? 0 : (plotMin < 0 ? -1 * flog10(plotMin * -1) : flog10(plotMin));
        float log_high = plotMax == 0 ? 0 : (plotMax < 0 ? -1 * flog10(plotMax * -1) : flog10(plotMax));
        
        double ret = xAxis ? (value - screenMin) : (getHeight() - value - screenMin);
        ret = log_low + (ret * (log_high - log_low)) / (screenMax - screenMin);
        ret = Math.pow(10, ret);
        
        return ret;
    }
    
    /**
     * Converts mouse location int X,Y into data points' value double rawX,rawY
     * The control variable invertY is assigned elsewhere in the class.
     * @param mouseY the mouse location Y
     * @return the rawY value of the corresponding data point.
     */
    public double getYValueFromYPixel(int mouseY) {
        if (getYAxis() == AXIS_SCALE.LIN) {
            if (invertY) {
                return plotYmax + ((double)(mouseY+canvasSectionMinimumY-getHeight())/(double)(canvasSectionMaximumY-canvasSectionMinimumY)*(double)(plotYmax-plotYmin));
            } else {
                return plotYmin - ((double)(mouseY+canvasSectionMinimumY-getHeight())/(double)(canvasSectionMaximumY-canvasSectionMinimumY)*(double)(plotYmax-plotYmin));
            }
        } else if (getYAxis() == AXIS_SCALE.LOG){
            return getLogValueFromPixel(mouseY, false);
        } else {
            Logicle l = getBiexScale(false);
            double y = (getHeight() - mouseY - TICK_LENGTH) / (double) canvasSectionMaximumY;
            double yVal = l.inverse(y);
            return yVal; 
        }
    }

    /**
	 * Highlights those points that need to be highlighted
	 * @param array
	 */
	public void highlightPoints(boolean[] array) {
		if (points.length != array.length) {
		    if (DEBUGGING) {
		        System.err.println("Error - mismatched array size when highlighting");
		    }
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
	    setImageStatus(IMAGE_STARTED);
        image = new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_INT_RGB);
        flow = true;
        drawAll(image.createGraphics(), true);
	}
	
	public void setLookupResolution(int lookupResolution) {
		this.lookupResolution = lookupResolution;
	}
	
	public IntVector lookupNearbyPoints(int x, int y, String pos) {
		IntVector iv = locLookup.get(pos);
		IntVector indicesOfDataPoints = new IntVector();

		for (int i = 0; iv != null && i < iv.size(); i++) {
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

	public void setFlow(boolean flow) {
		this.flow = flow;
	}

	public boolean getFlow() {
		return flow;
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
