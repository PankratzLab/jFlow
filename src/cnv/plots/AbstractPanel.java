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
import java.net.URL;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
//import java.util.Date;
import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;

import javax.imageio.ImageIO;
import javax.swing.JButton;
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
	public final int DELAY = 0;	//A control variable to reduce the repaint() operations during component resizing;
	
	protected Color[] colorScheme;
	protected int canvasSectionMinimumX = 0;
	protected int canvasSectionMaximumX = getWidth();
	protected int canvasSectionMinimumY = 0;
	protected int canvasSectionMaximumY = getWidth();
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
	protected IntVector indicesOfNaNSamples;	//zx
	private boolean inDrag;
	private int startX, startY;
	private int plotPointSetSize;
	private int totalNumPlotPointSets;
	private int currentPlotPointSet;
	private int lastIndexInPlotPointSet;
	private int currentIndexInPlotPointSet;
	private String tempDirectory;
	private int lookupResolution;
	private boolean flow;			//A control variable. If resizing is not yet done, don't start generatePoints() or drawAll();
	private boolean finalImage;		//A control variable. If drawAll() is not yet done, don't start paintComponent();
	private byte[] layersInBase;
	private byte[] extraLayersVisible;
	private boolean pointsGeneratable;
	private Timer waitingTimer;		//A control variable to reduce the repaint() operations during component resizing;
	private String nullMessage;
	private boolean randomTest;
	
	public AbstractPanel() {
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
		tempDirectory = "";
		randomTest = false;
		
		layersInBase = null;
		extraLayersVisible = null;

		image = null;
		locLookup = new Hashtable<String,IntVector>();
		finalImage=true;
		flow=true;
		pointsGeneratable = true;
		
		colorScheme = new Color[] {Color.BLACK, Color.GRAY};
		addMouseListener(this);
		addMouseMotionListener(this);
		addMouseWheelListener(this);
	}
	
	public void addPlotPoint(PlotPoint point) {
		lastIndexInPlotPointSet++;
		points[lastIndexInPlotPointSet] = point;
		if (lastIndexInPlotPointSet == plotPointSetSize) {
			new PlotPointSet(points).serialize(tempDirectory+"PlotPointSet."+currentPlotPointSet+".ser");
		}
		points = new PlotPoint[plotPointSetSize];
		totalNumPlotPointSets++;
		currentPlotPointSet++;
		lastIndexInPlotPointSet = 0;
	}
	
	public void resetCurrentIndexInPlotPointSet() {
		currentIndexInPlotPointSet = -1;
	}
	
	public boolean morePlotPointsExist() {
		return currentPlotPointSet<totalNumPlotPointSets && currentIndexInPlotPointSet<=lastIndexInPlotPointSet;
	}

	public PlotPoint nextPlotPoint() {
		PlotPoint point;
		
		currentIndexInPlotPointSet++;
		if (currentIndexInPlotPointSet == plotPointSetSize) {
			
		}
		
//		currentPlotPointSet<totalNumPlotPointSets && currentIndexInPlotPointSet<=lastIndexInPlotPointSet
//		boolean hi = true;
		point = null;
		
		return point;
	}
	
	public void setNullMessage(String str) {
		nullMessage = str;
	}

	public void setTempDirectory(String dir) {
		tempDirectory = dir;
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
	
	public void setLayersInBase(byte[] layers) {
		layersInBase = layers;
	}

	public void setExtraLayersVisible(byte[] layers) {
		extraLayersVisible = layers;
	}

	public void paintComponent(Graphics g) {
		if (getFinalImage()&&image==null) {
			createImage();
			repaint();
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

	abstract void generatePoints();

	abstract void highlightPoints();

	abstract void assignAxisLabels();

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
    	int recX, recY, recWidth, recHeight;
    	
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
		//System.out.println(minimumObservedRawX+"\t"+maximumObservedRawX+"\t"+minimumObservedRawY+"\t"+maximumObservedRawY);//zx

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
		
		if (base) {
			indicesOfNaNSamples = new IntVector();
		
			//g.setColor(Color.WHITE);
			g.fillRect(0, 0, getWidth(), getHeight());
			g.setFont(new Font("Arial", 0, AXIS_FONT_SIZE));
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
						Grafik.drawThickLine(g, getX(x), getHeight()-canvasSectionMaximumY, getX(x), getHeight()-(canvasSectionMaximumY-TICK_LENGTH), TICK_THICKNESS, Color.BLACK);
						str = ext.formDeci(Math.abs(x)<DOUBLE_INACCURACY_HEDGE?0:x, sigFigs, true);
						g.drawString(str, getX(x)-str.length()*8, getHeight()-(canvasSectionMaximumY-TICK_LENGTH-30));
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
						Grafik.drawThickLine(g, canvasSectionMaximumX-TICK_LENGTH, getY(y), canvasSectionMaximumX, getY(y), TICK_THICKNESS, Color.BLACK);
						str = ext.formDeci(Math.abs(y)<DOUBLE_INACCURACY_HEDGE?0:y, sigFigs, true);
						g.drawString(str, canvasSectionMaximumX-TICK_LENGTH-str.length()*15-5, getY(y)+9);
					}
				}
				Grafik.drawThickLine(g, canvasSectionMaximumX, getY(plotYmin), canvasSectionMaximumX, getY(plotYmax)-(int)Math.ceil((double)TICK_THICKNESS/2.0), AXIS_THICKNESS, Color.BLACK);
	
				yLabel = new BufferedImage(fontMetrics.stringWidth(yAxisLabel), 36, BufferedImage.TYPE_INT_RGB);
				gfx = yLabel.createGraphics();
				gfx.setFont(new Font("Arial", 0, 28));
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

		canvasSectionMinimumX = WIDTH_Y_AXIS;
		canvasSectionMaximumX = getWidth()-WIDTH_BUFFER;
		canvasSectionMinimumY = HEIGHT_X_AXIS;
		canvasSectionMaximumY = getHeight()-HEAD_BUFFER;
//		System.err.println("("+canvasSectionMinimumX+"-"+canvasSectionMaximumX+","+canvasSectionMinimumY+"-"+canvasSectionMaximumY+")");

		// Draw the lines
		for (int i = 0; lines!=null && i<lines.length && flow; i++) {
//		for (int i = 0; lines!=null&&i<lines.length; i++) {
			if ((base && (layersInBase == null || Array.indexOfByte(layersInBase, lines[i].getLayer()) >= 0)) || (!base && Array.indexOfByte(extraLayersVisible, lines[i].getLayer()) >= 0)) {
				Grafik.drawThickLine(g, getX(lines[i].getStartX()), getY(lines[i].getStartY()), getX(lines[i].getStopX()), getY(lines[i].getStopY()), (int)lines[i].getThickness(), colorScheme[lines[i].getColor()]);
			}
        }

		// Draw the rectangles for clusterFilters
		for (int i = 0; rectangles!=null && i<rectangles.length && flow; i++) {
//		for (int i = 0; rectangles!=null&&i<rectangles.length; i++) {
			if ((base && (layersInBase == null || Array.indexOfByte(layersInBase, rectangles[i].getLayer()) >= 0)) || (!base && Array.indexOfByte(extraLayersVisible, rectangles[i].getLayer()) >= 0)) {
//				recX = getX(Math.min((int)rectangles[i].getStartX(), (int)rectangles[i].getStopX()));
//				recY = getY(Math.min((int)rectangles[i].getStartY(), (int)rectangles[i].getStopY()));
//		    	recWidth = getX(Math.abs((int)rectangles[i].getStartX() - (int)rectangles[i].getStopX()));
//		    	recHeight = getY((Math.abs((int)rectangles[i].getStartY() - (int)rectangles[i].getStopY())));
				recX = Math.min(getX(rectangles[i].getStartX()), getX(rectangles[i].getStopX()));
				recY = Math.min(getY(rectangles[i].getStartY()), getY(rectangles[i].getStopY()));
		    	recWidth = Math.abs(getX(rectangles[i].getStartX()) - getX(rectangles[i].getStopX()));
		    	recHeight = (Math.abs(getY(rectangles[i].getStartY()) - getY(rectangles[i].getStopY())));
//		    	g.setColor(Color.GRAY);
		    	g.setColor(colorScheme[rectangles[i].getColor()]);
				
				if (rectangles[i].getFill()) {
					if (rectangles[i].getRoundedCorners()) {
						g.fillRoundRect(recX, recY, recWidth, recHeight, 2, 2);
					} else {
						g.fillRect(recX, recY, recWidth, recHeight);
					}
				} else {
					if (rectangles[i].getRoundedCorners()) {
						g.drawRoundRect(recX, recY, recWidth, recHeight, 2, 2);
					} else {
//						g.drawRect(recX, recY, recWidth, recHeight);
				    	drawRectThick(g, recX, recY, recWidth, recHeight, rectangles[i].getThickness());
					}
				}
			}
        }

		// Draw the rectangle outlined by dragging the mouse
		if (highlightRectangle != null) {
			recX = Math.min(getX(highlightRectangle.getStartX()), getX(highlightRectangle.getStopX()));
			recY = Math.min(getY(highlightRectangle.getStartY()), getY(highlightRectangle.getStopY()));
	    	recWidth = Math.abs(getX(highlightRectangle.getStartX()) - getX(highlightRectangle.getStopX()));
	    	recHeight = (Math.abs(getY(highlightRectangle.getStartY()) - getY(highlightRectangle.getStopY())));
//	    	g.setColor(Color.BLACK);
	    	g.setColor(colorScheme[0]);
//	    	g.drawRect(recX, recY, recWidth, recHeight);
	    	drawRectThick(g, recX, recY, recWidth, recHeight, (byte) 1);
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
		locLookup.clear();	// -- This was here in Nathan's original code.
		prog = null;
		time = new Date().getTime();
		step = Math.max((points.length)/100, 1);
		layers = new Hashtable<String,Vector<PlotPoint>>();
		
		for (int i = 0; i<points.length && flow; i++) {
//			System.out.println("loop");
//		for (int i = 0; i<points.length; i++) {
			if (base && i%step==0){
				if (new Date().getTime() - time > 1000) {
					if (prog == null) {
						prog = new ProgressBarDialog("Generating image...", 0, points.length, getWidth(), getHeight(), 5000);//zx
					}
					prog.setProgress(i);//zx
				}
			}
			if (points[i] == null || points[i].getColor() == -1) {
//			if (points[i] == null || points[i].getColor() == -1 || (points[i].getRawX() < 1 && points[i].getRawY() < 1)) {
				
			} else if (truncate && (points[i].getRawX() < plotXmin || points[i].getRawX() > plotXmax || points[i].getRawY() < plotYmin || points[i].getRawY() > plotYmax)) {
//				System.err.println("error: data point ("+points[i].getRawX()+","+points[i].getRawY()+") is outside of plot range.");
			} else {
				trav = points[i].getLayer()+"";
				if (points[i].isHighlighted() || (base && (layersInBase == null || Array.indexOfByte(layersInBase, points[i].getLayer()) >= 0)) || (!base && Array.indexOfByte(extraLayersVisible, points[i].getLayer()) >= 0)) {
					if (trav.equals("0")) {
						if (points[i].getType()!=PlotPoint.NOT_A_NUMBER) {
							drawPoint(g, points[i]);
						} else if (base) {
							indicesOfNaNSamples.add(i);
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
				if (createLookup) {
					xLook = (int)Math.floor(getX(points[i].getRawX())/lookupResolution);
					yLook = (int)Math.floor(getY(points[i].getRawY())/lookupResolution);
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
		
		//buildLookupTableOfNearbyPoints();//zx Looks like it gets duplicated with the above
		//for (int i=0; i<locLookup.size(); i++) {
		//	System.out.println(locLookup.get(i));
		//}

		// Draw those points with layer>0.
		keys = HashVec.getKeys(layers);
		order = Sort.quicksort(Array.toIntArray(keys));
		for (int i = 0; i<keys.length&&flow; i++) {
//		for (int i = 0; i<keys.length; i++) {
			layer = layers.get(keys[order[i]]);
			for (int j = 0; j<layer.size(); j++) {
				if (points[i].getType()!=PlotPoint.NOT_A_NUMBER) {
					drawPoint(g, layer.elementAt(j));
				} else {
					indicesOfNaNSamples.add(j);
					//TODO This is a problem
				}
            }
        }

		if (indicesOfNaNSamples!=null && indicesOfNaNSamples.size()>0) {
			g.drawString(PlotPoint.NAN_STR+" (n="+indicesOfNaNSamples.size()+")", getX(0)-nanWidth/2, getY(0)+60+points[0].getSize()/2);
		}
		
		if (base && displayGrid) {
			for (double d = 0; d < 1.0; d+=0.1) {
				g.drawLine(getX(d), getY(0), getX(d), getY(canvasSectionMaximumY));
			}
			for (double d = -0.5; d < 0.5; d+=0.1) {
				g.drawLine(getX(0), getY(d), getX(canvasSectionMaximumX), getY(d));
			}
		}
		setFinalImage(true);
		
		if (base && prog != null) {
			prog.close();//zxu
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
	
	public void mouseClicked(MouseEvent e) {}

	public void mouseEntered(MouseEvent e) {}

	public void mouseExited(MouseEvent e) {}

	public void mouseMoved(MouseEvent e) {}
	
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
//		setFlow(false);//zx
		if (this.waitingTimer==null) {
			/* Start waiting for DELAY to elapse. */
			this.waitingTimer = new Timer(DELAY,this);
			this.waitingTimer.start();
		}
		else {
			/* Event came too soon, swallow it by resetting the timer.. */
			this.waitingTimer.restart();
		}

//		createImage();

//		if (patience!=null) {
//			patience.cancel();
//			patience = null;
//		}
//		if (prevWidth==-1) {
//			prevWidth = getWidth();
//			prevHeight = getHeight();
//		} else if (prevWidth!=getWidth()||prevHeight!=getHeight()) {
////			image = null;
////			repaint();
////			new Thread(patience = new Repress(this, 100)).start();
//
//			prevWidth = getWidth();
//			prevHeight = getHeight();
//		}
//		setFlow(true);//zx
	}

	public void componentShown(ComponentEvent e) {}

	public void actionPerformed(ActionEvent ae)
	{
	  /* Timer finished? */
	  if (ae.getSource()==this.waitingTimer)
	  {
	    /* Stop timer */
	    this.waitingTimer.stop();
	    this.waitingTimer = null;
	    /* Resize */
		setFlow(true);//zx
	    createImage();
	    repaint();
//	    setFlow(false);//zx
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
		
//		getGraphics().fillRect(mouseStartX, mouseStartY, mouseEndX-mouseStartX, mouseEndY-mouseStartY);
		
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
			g.fillOval(getX(point.getRawX())-point.getSize()/2, getY(point.getRawY())-point.getSize()/2, point.getSize(), point.getSize());
			break;
		case PlotPoint.OPEN_CIRCLE:
			g.drawOval(getX(point.getRawX())-point.getSize()/2, getY(point.getRawY())-point.getSize()/2, point.getSize(), point.getSize());
			break;
		case PlotPoint.MISSING:
			//g.drawString(PlotPoint.MISSING_STR, getX(point.getRawX())-missingWidth/2, getY(point.getRawY())+point.getSize()/2);
			if (PlotPoint.MISSING_STR=="X" || PlotPoint.MISSING_STR=="x"){
				g.drawLine(getX(point.getRawX())-point.getSize()/4, getY(point.getRawY())-point.getSize()/4, getX(point.getRawX())+point.getSize()/4, getY(point.getRawY())+point.getSize()/4);//zx
				g.drawLine(getX(point.getRawX())-point.getSize()/4, getY(point.getRawY())+point.getSize()/4, getX(point.getRawX())+point.getSize()/4, getY(point.getRawY())-point.getSize()/4);//zx
			}
			else {
				g.setFont(new Font("Arial", 0, point.getSize()));//zx
				g.drawString(PlotPoint.MISSING_STR, getX(point.getRawX())-point.getSize()/2, getY(point.getRawY())+point.getSize()/2);//zx
				g.setFont(new Font("Arial", 0, AXIS_FONT_SIZE));//zx
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

	public int getX(double x) {
		if (invertX) {
			return (int)((plotXmax-x)/(plotXmax-plotXmin)*(double)(canvasSectionMaximumX-canvasSectionMinimumX))+canvasSectionMinimumX;
		} else {
			return (int)((x-plotXmin)/(plotXmax-plotXmin)*(double)(canvasSectionMaximumX-canvasSectionMinimumX))+canvasSectionMinimumX;
		}
	}

	public int getY(double y) {
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
	public double getRawX(int mouseX) {
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
	public double getRawY(int mouseY) {
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
				if (array[i]) {
					points[i].setHighlighted(true);
				} else {
					points[i].setHighlighted(false);
				}
				
			}
		}
	}
	
	public void createImage() {
		if (getWidth() > 350 && getHeight() > 0 ) {
			image = new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_INT_RGB);
			flow = true;
			drawAll(image.createGraphics(), true);
//			repaint();
//			image = new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_INT_RGB);
		}
	}
	
	//zx
	public void setLookupResolution(int lookupResolution) {
		this.lookupResolution = lookupResolution;
	}
	
	//zx
	/*
	public void buildLookupTableOfNearbyPoints() {
		int x, y;
		String pos;
		locLookup.clear();
		//System.out.println("---Beginning of locLookup---\npos\tLocation\tSample Index");//test point
		for (int i = 0; i<points.length; i++) {
			//x = (int)Math.floor(getX(data[i][0])/LOOKUP_RESOLUTION);
			//y = (int)Math.floor(getY(data[i][1])/LOOKUP_RESOLUTION);
			x = (int)Math.floor(getX(points[i].getRawX())/lookupResolution);
			y = (int)Math.floor(getY(points[i].getRawY())/lookupResolution);
			for (int j = x-1; j<=x+1; j++) {
				for (int k = y-1; k<=y+1; k++) {
					pos = j+"x"+k;
					if (locLookup.containsKey(pos)) {
						locLookup.get(pos).add(i);
					} else {
						locLookup.put(pos, new IntVector(new int[] {i}));
					}
					//System.out.println(pos+"\t("+getX(points[i].getRawX())+", "+getY(points[i].getRawY())+")\t"+i);//test point
				}
			}
		}
		//System.out.println("---End of locLookup---\n\n\n");//test point
	}
	*/

	//zx
	public IntVector lookupNearbyPoints(int x, int y, String pos) {
		IntVector iv = locLookup.get(pos);
		//System.out.print("\t iv size: "+iv.size());
		IntVector indeciesOfDataPoints = new IntVector();
		//System.out.println("---Log of lookup process---\nMouse pos\tMouse location\tNumber of nearby samples\tNearby sample index\tNearby sample location\tDistance");//test point
		for (int i = 0; iv!=null && i<iv.size(); i++) {
			//System.out.println(pos+"\t("+x+", "+y+")\t"+(iv==null?"null":iv.size())+"\t"+iv.elementAt(i)+"\t("+getX(points[iv.elementAt(i)].getRawX())+", "+getY(points[iv.elementAt(i)].getRawY())+")\t"+Distance.euclidean(new int[] {x, y}, new int[] {getX(points[iv.elementAt(i)].getRawX()), getY(points[iv.elementAt(i)].getRawY())}));//test point
			if (Distance.euclidean(new int[] {x, y}, new int[] {getX(points[iv.elementAt(i)].getRawX()), getY(points[iv.elementAt(i)].getRawY())})<HIGHLIGHT_DISTANCE) {
				indeciesOfDataPoints.add(iv.elementAt(i));
			}
		}
		//System.out.println("\t indeciesOfDataPoints: "+indeciesOfDataPoints.size());
		return indeciesOfDataPoints;
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
	
//	public void interruptFlow1() {
//		flow = false;
//	}
//	
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
