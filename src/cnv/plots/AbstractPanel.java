package cnv.plots;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
//import java.util.Date;
import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;

import javax.imageio.ImageIO;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import mining.Distance;

import common.Array;
import common.Grafik;
import common.HashVec;
import common.IntVector;
import common.ProgressBarDialog;
import common.Sort;
import common.ext;
import stats.Maths;

public abstract class AbstractPanel extends JPanel implements MouseListener, MouseMotionListener, MouseWheelListener {
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
	public static final int SIZE = 12;//zx
	public static final double HIGHLIGHT_DISTANCE = 20;//= Math.sqrt(SIZE*SIZE/2);//zx

	
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
	protected AbstractLine[] lines;
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
	protected boolean truncate;
	protected float[][] zoomSubsets;
	private boolean inDrag;
	private int startX, startY;
	private int plotPointSetSize;
	private int totalNumPlotPointSets;
	private int currentPlotPointSet;
	private int lastIndexInPlotPointSet;
	private int currentIndexInPlotPointSet;
	private String tempDirectory;
	private int lookupResolution;	   //zx
	private boolean flow;			   //zx: If resizing is not yet done, don't start generatePoints() or drawAll();
	private boolean finalImage;		   //zx: If drawAll() is not yet done, don't start paintComponent();

	
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
		resetZoomProportions();
		plotPointSetSize = DEFAULT_PLOTPOINTSET_SIZE;
		points = new PlotPoint[plotPointSetSize];
		totalNumPlotPointSets = 1;
		currentPlotPointSet = 0;
		lastIndexInPlotPointSet = -1;
		currentIndexInPlotPointSet = -1;
		tempDirectory = "";

		image = null;
		locLookup = new Hashtable<String,IntVector>();
		finalImage=true;//zx
		flow=true;//zx
		
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

	public void paintComponent(Graphics g) {
		if (getFinalImage()&&image==null) {
			createImage();
			repaint();
		}
		g.drawImage(image, 0, 0, this);
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

	public void drawAll(Graphics g) {
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
		
		//zx
		if (points.length==0) {
			System.out.println("Error: no data. The cnv.plots.AbstractPanel.points is null.");
			return;
		}
		
		setFinalImage(false);
		generatePoints();
		setLookupResolution(DEFAULT_LOOKUP_RESOLUTION);//zx
		assignAxisLabels();

		minimumObservedRawX = Float.MAX_VALUE;
		maximumObservedRawX = Float.MIN_VALUE;
		minimumObservedRawY = Float.MAX_VALUE;
		maximumObservedRawY = Float.MIN_VALUE;
		for (int i = 0; i<points.length&&flow; i++) {
			if (points[i] != null) {
				minimumObservedRawX = Maths.min(minimumObservedRawX, points[i].getRawX());
				maximumObservedRawX = Maths.max(maximumObservedRawX, points[i].getRawX());
				minimumObservedRawY = Maths.min(minimumObservedRawY, points[i].getRawY());
				maximumObservedRawY = Maths.max(maximumObservedRawY, points[i].getRawY());
				//System.out.println(minimumObservedRawX+"\t"+maximumObservedRawX+"\t"+minimumObservedRawY+"\t"+maximumObservedRawY);//zx
				//System.out.println(points[i].getRawX()+"\t"+points[i].getRawY());//zx
			}
        }
		for (int i = 0; lines != null && i<lines.length&&flow; i++) {
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
		maximumObservedRawX = Float.isNaN(forcePlotXmax)?maximumObservedRawX:forcePlotXmax;
		minimumObservedRawY = Float.isNaN(forcePlotYmin)?minimumObservedRawY:forcePlotYmin;
		maximumObservedRawY = Float.isNaN(forcePlotYmax)?maximumObservedRawY:forcePlotYmax;
		
		if (makeSymmetric) {
			maximumObservedRawX = Math.max(maximumObservedRawX, maximumObservedRawY);
			maximumObservedRawY = maximumObservedRawX;
			minimumObservedRawX = Math.min(minimumObservedRawX, minimumObservedRawY);
			minimumObservedRawY = minimumObservedRawX;
		}
		
		//g.setColor(Color.WHITE);
//		g.setColor(Color.BLACK);//zx tmp
		g.fillRect(0, 0, getWidth(), getHeight());
		g.setFont(new Font("Arial", 0, AXIS_FONT_SIZE));
//		System.out.println("getWidth: "+getWidth()+"\t getHeight: "+getHeight());
		
		fontMetrics = g.getFontMetrics(g.getFont());
		missingWidth = fontMetrics.stringWidth("X");
		missingWidth = fontMetrics.stringWidth("X");

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
		
		
		canvasSectionMinimumX = WIDTH_Y_AXIS;
		canvasSectionMaximumX = getWidth()-WIDTH_BUFFER;
		canvasSectionMinimumY = HEIGHT_X_AXIS;
		canvasSectionMaximumY = getHeight()-HEAD_BUFFER;
		
		if (errorMessage != null) {
			g.drawString(errorMessage, (getWidth()-WIDTH_Y_AXIS)/2-fontMetrics.stringWidth(errorMessage)/2+WIDTH_Y_AXIS, (getHeight()-HEAD_BUFFER-HEIGHT_X_AXIS)/2-20+HEAD_BUFFER);
		}
		
		for (int i = 0; lines!=null&&i<lines.length&&flow; i++) {
			Grafik.drawThickLine(g, getX(lines[i].getStartX()), getY(lines[i].getStartY()), getX(lines[i].getStopX()), getY(lines[i].getStopY()), (int)lines[i].getThickness(), colorScheme[lines[i].getColor()]);
        }
		time = new Date().getTime();
		prog = new ProgressBarDialog("Generating image...", 0, points.length, getWidth(), getHeight(), 500);//zx
//		System.out.println("points.length: "+(points.length)+"\3*points.length: "+3*(points.length));
		step = (points.length)/100;

		locLookup.clear();	// -- This was here since Nathan's original code.
		layers = new Hashtable<String,Vector<PlotPoint>>();
		for (int i = 0; i<points.length&&flow; i++) {
			if (i%step==0){
				prog.setProgress(i);//zx
			}
			if (points[i] == null || points[i].getColor() == -1) {
//			if (points[i] == null || points[i].getColor() == -1 || (points[i].getRawX() < 1 && points[i].getRawY() < 1)) {
				
			} else if (truncate && (points[i].getRawX() <= plotXmin || points[i].getRawX() > plotXmax || points[i].getRawY() <= plotYmin || points[i].getRawY() > plotYmax)) {
				
			} else {
				trav = points[i].getLayer()+"";
				if (trav.equals("0")) {
					drawPoint(g, points[i]);
				} else {
					if (layers.containsKey(trav)) {
						layer = layers.get(trav);
					} else {
						layers.put(trav, layer = new Vector<PlotPoint>());
					}
					layer.add(points[i]);
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
							//System.out.println("pos: "+pos+"\t sample id: "+i);//zx
						}
					}
				}
			}
		}
		
		//buildLookupTableOfNearbyPoints();//zx Looks like it gets duplicated with the above
		//for (int i=0; i<locLookup.size(); i++) {
		//	System.out.println(locLookup.get(i));
		//}

		keys = HashVec.getKeys(layers);
		order = Sort.quicksort(Array.toIntArray(keys));
		for (int i = 0; i<keys.length&&flow; i++) {
			layer = layers.get(keys[order[i]]);
			for (int j = 0; j<layer.size(); j++) {
				drawPoint(g, layer.elementAt(j));
            }
        }		

		if (displayGrid) {
			for (double d = 0; d < 1.0; d+=0.1) {
				g.drawLine(getX(d), getY(0), getX(d), getY(canvasSectionMaximumY));
			}
			for (double d = -0.5; d < 0.5; d+=0.1) {
				g.drawLine(getX(0), getY(d), getX(canvasSectionMaximumX), getY(d));
			}
		}
		setFinalImage(true);
		prog.close();//zxu
		//System.out.println("Paint time: "+ext.getTimeElapsed(time));		
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
			g.drawString(PlotPoint.NAN_STR, getX(point.getRawX())-nanWidth/2, getY(point.getRawY())-30+point.getSize()/2);
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

	public void createImage() {
		image = new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_INT_RGB);
		drawAll(image.createGraphics());
		flow = true;
//		repaint();
//		image = new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_INT_RGB);
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
		for (int i = 0; iv!=null&&i<iv.size(); i++) {
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
	
	public void interruptFlow() {
		flow = false;
	}

}
