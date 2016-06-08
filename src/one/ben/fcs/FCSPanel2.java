package one.ben.fcs;

import java.awt.Color;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.AffineTransform;
import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;

import javax.swing.SwingUtilities;

import one.ben.fcs.gating.Gate;
import one.ben.fcs.gating.Gate.PolygonGate;
import one.ben.fcs.gating.Gate.RectangleGate;
import one.ben.fcs.gating.GateDimension.RectangleGateDimension;
import one.ben.fcs.gating.GateDimension;
import stats.Histogram;
import cnv.plots.GenericLine;
import cnv.plots.GenericPath;
//import cnv.filesys.MarkerLookup;
import cnv.plots.GenericRectangle;
import cnv.plots.PlotPoint;
import common.Array;

public class FCSPanel2 extends AbstractPanel2 implements MouseListener, MouseMotionListener {
	public static final long serialVersionUID = 3L;
	public static final int LOOKUP_RESOLUTION = 20;
	public static final Color[] DEFAULT_COLORS = {
                                                new Color(33, 31, 53), // dark dark
			   									new Color(201, 30, 10), // deep red
			   									new Color(94, 88, 214), // light purple
			   									new Color(189, 243, 61), // light green
			   									new Color(217, 109, 194), // pink
			   									new Color(33, 87, 0), // dark green
			   									new Color(23, 58, 172), // dark blue
			   									new Color(140, 20, 180), // deep purple
			   									new Color(182, 182, 182), // light grey
			   									new Color(220, 220, 220), // very light grey
			   									new Color(0, 0, 128), // ALL KINDS OF BLUES
			   									new Color(55, 129, 252), // light blue
			   									new Color(100, 149, 237),
			   									new Color(72, 61, 139),
			   									new Color(106, 90, 205),
			   									new Color(123, 104, 238),
			   									new Color(132, 112, 255),
			   									new Color(0, 0, 205),
			   									new Color(65, 105, 225),
			   									new Color(0, 0, 255),
			   									new Color(30, 144, 255),
			   									new Color(0, 191, 255),
			   									new Color(135, 206, 250),
			   									new Color(135, 206, 250),
			   									new Color(70, 130, 180),
			   									new Color(176, 196, 222),
			   									new Color(173, 216, 230),
			   									new Color(176, 224, 230),
			   									new Color(175, 238, 238),
			   									new Color(0, 206, 209),
			   									new Color(72, 209, 204),
			   									new Color(64, 224, 208),
			   									new Color(0, 255, 255),
			   									new Color(224, 255, 255),

	};
    private static final byte POINT_SIZE = 5;
	
	protected FCSPlot fcp;
	
	public FCSPanel2(FCSPlot fcsPlot) {
		super();
		setDoubleBuffered(false);
		
		this.fcp = fcsPlot;
		this.setAxisFontSize(24);
		this.setSymmetricAxes(false);
		setZoomable(true, true);
//		setZoomable(false, true);

		setColorScheme(DEFAULT_COLORS);

		setNullMessage("Select two variables to plot");
		
	}

    int dataCount = -1;
	String xCol = null, yCol = null;
	double xMed = Double.NaN, xMin = Double.NaN, xMax = Double.NaN, yMed = Double.NaN, yMin = Double.NaN, yMax = Double.NaN, xSD = Double.NaN, ySD = Double.NaN;
	PLOT_TYPE prevPlotType;
	boolean[] showMedSD = {false, false, false, false}; // xMed, xSD, yMed, ySD
	float[] xData;
	float[] yData;
	ArrayList<GenericRectangle> rects = new ArrayList<GenericRectangle>();
	ArrayList<GenericPath> polys = new ArrayList<GenericPath>();
	
	ArrayList<GenericLine> histLines = new ArrayList<GenericLine>();
	ArrayList<double[]> tempPoly = new ArrayList<double[]>();
    volatile boolean forceGatesChanged = false;
	
	public static final int RECT_TOOL = 0;
	public static final int POLY_TOOL = 1;
	
	private volatile int currentTool = POLY_TOOL;
//	private volatile int currentTool = RECT_TOOL;
	
	private int dragInd = -1;
    private int polyDragVertInd = -1;

    private boolean isHistogram() {
        return yCol != null && yCol.equals(FCSPlot.HISTOGRAM_COL);
    }

    private void setLines() {
        ArrayList<GenericLine> allLines = new ArrayList<GenericLine>();
        allLines.addAll(histLines);
        if (highlightRectangle != null ) {
            allLines.add(new GenericLine(
                    highlightRectangle.getStartXValue(),
                    (float)plotYmin,
                    highlightRectangle.getStartXValue(),
                    (float)plotYmax,
                    (byte)1,
                    (byte)1,
                    (byte)99
                    ));
            allLines.add(new GenericLine(
                    highlightRectangle.getStopXValue(),
                    (float)plotYmin,
                    highlightRectangle.getStopXValue(),
                    (float)plotYmax,
                    (byte)1,
                    (byte)1,
                    (byte)99
                    ));
        }
        for (GenericRectangle rect : rects) {
            allLines.add(new GenericLine(
                    rect.getStartXValue(),
                    (float)plotYmin,
                    rect.getStartXValue(),
                    (float)plotYmax,
                    (byte)1,
                    (byte)1,
                    (byte)99
                    ));
            allLines.add(new GenericLine(
                    rect.getStopXValue(),
                    (float)plotYmin,
                    rect.getStopXValue(),
                    (float)plotYmax,
                    (byte)1,
                    (byte)1,
                    (byte)99
                    ));
        }
        
        lines = allLines.toArray(new GenericLine[allLines.size()]);
    }
    
    public void generatePointsRectanglesAndLines() {
		byte type;
//		String[] line;
		float xAxisValue, yAxisValue;
		byte size = POINT_SIZE;
		
		if (fcp.dataLoader == null && !fcp.isLoading) {
		    setNullMessage("Please load an FCS file..");
		    points = new PlotPoint[0];
            rectangles = new GenericRectangle[0];
            lines = new GenericLine[0];
            polygons = new GenericPath[0];
		    return;
		}
		if (!fcp.isCurrentDataDisplayable()) {
		    setNullMessage("Please wait, data is loading...");
		    points = new PlotPoint[0];
            rectangles = new GenericRectangle[0];
            lines = new GenericLine[0];
            polygons = new GenericPath[0];
		    return;
		}
		
		boolean columnsChangedX = false;
		boolean columnsChangedY = false;
		boolean dataChanged = false;
		boolean optionsChanged = false;
		boolean gatesChanged = false;
		boolean typeChanged = false;
		{
    		String newX = fcp.getXDataName();
    		String newY = fcp.getYDataName();
    		if (xCol == null || !fcp.getXDataName().equals(xCol)) {
    		    columnsChangedX = true;
    		}
    		if (yCol == null || !fcp.getYDataName().equals(yCol)) {
    		    columnsChangedY = true;
    		}
            xCol = newX;
            yCol = newY;
		}
		{
    		int count = fcp.getDataCount();
    	    if (count != dataCount) {
    	        dataChanged = true;
    	    }
    		dataCount = count;
		}
		{
    		PLOT_TYPE plotType = fcp.getPlotType();
    		if (plotType != prevPlotType) {
    		    typeChanged = true;
    		}
    		prevPlotType = plotType;
		}
		{
            boolean mX = fcp.showMedian(false), mY = fcp.showMedian(true), sdX = fcp.showSD(false), sdY = fcp.showSD(true); 
            if (mX != showMedSD[0] || sdX != showMedSD[1] || mY != showMedSD[2] || sdY != showMedSD[3]) {
                optionsChanged = true;
            }
            showMedSD[0] = mX;
            showMedSD[1] = sdX;
            showMedSD[2] = mY;
            showMedSD[3] = sdY;
		}
		{
		    gatesChanged = rectangles != null && rectangles.length != rects.size() || polygons != null && polygons.length != polys.size() || forceGatesChanged; 
		    forceGatesChanged = false;
		}
		
		// TODO replace with loading gates from FCSPlot (and creating gates from mouse actions)
        rectangles = rects.toArray(new GenericRectangle[rects.size()]);
        polygons = polys.toArray(new GenericPath[polys.size()]);
        
        if (isHistogram()) {
            setLines();
        }
        
//        setForcePlotXMin(0);
//        setForcePlotYMin(0);
        
		boolean skip = !columnsChangedX && !columnsChangedY && !dataChanged && !optionsChanged && !gatesChanged/* && !typeChanged /* don't need to regen if only type has changed, for now */;
		if (skip) return;
		
		xData = columnsChangedX || dataChanged || xData == null ? fcp.getAxisData(false, true) : xData;
		yData = columnsChangedY || dataChanged || yData == null ? isHistogram() ? null : fcp.getAxisData(false, false) : yData;

        ArrayList<boolean[]> gates = new ArrayList<boolean[]>();
        ArrayList<Gate> gating;
        if (fcp.gating == null) {
             gating = new ArrayList<Gate>();
        } else {
            
            // TODO only retrieve gating for a particular branch
            // TODO  --  determine which file is on which branch based on gating file
            
            if (isHistogram()) {
                gating = fcp.gating.getGatesForParamOnly(xCol);
            } else {
                gating = fcp.gating.getGatesForParams(xCol, yCol);
            }
        }
//        rects.clear();
//        polys.clear();
        for (Gate g : gating) {
            gates.add(g.gate(fcp.dataLoader));
            if (g instanceof RectangleGate) {
                RectangleGate rg = (RectangleGate) g;
                RectangleGateDimension rgdX = rg.getDimension(xCol);
                RectangleGateDimension rgdY = isHistogram() ? null : rg.getDimension(yCol);
                rects.add(new GenericRectangle(
                      rgdX.getMin(), 
                      (float)(isHistogram() ? plotYmin + (plotYmax - plotYmin) / 2 : rgdY.getMin()), 
                      rgdX.getMax(),
                      (float)(isHistogram() ? plotYmin + (plotYmax - plotYmin) / 2 : rgdY.getMax()), 
                      (byte)1, false, false, (byte)0, (byte)99, true));
            } else if (g instanceof PolygonGate) {
                polys.add(new GenericPath(((PolygonGate)g).getPath(), (byte)0, (byte)0, (byte)99, false, true));
            }
        }
        
		if (isHistogram()) {
		    points = new PlotPoint[0];
		    if (!columnsChangedX && !dataChanged && histLines != null && histLines.size() > 0) return;
		    
		    float[] minMax = Array.minMax(xData);
		    int range = (int)Math.ceil(minMax[1]) - (int)Math.floor(minMax[0]) + 1;
//		    int step1 = (int) Math.ceil(range / (2 * Array.iqr(xData) / Math.pow(xData.length, 1/3))); // bin size too small with this formula
		    int step = Math.max(2, (int) (range / Math.sqrt(xData.length)));
		    float[] histData = new float[range / step + 1];
		    for (float x : xData) {
		        histData[(int) (x - minMax[0]) / step]++;
		    }
		    
		    histLines = new ArrayList<GenericLine>();
		    for (int i = 0; i < histData.length - 1; i++) {
		        histLines.add(new GenericLine((i * step) + minMax[0], histData[i], ((i + 1) * step) + minMax[0], histData[i + 1], (byte)2, (byte) 0, (byte)0));
		    }
		    setLines();
		    
		    return;
		}
		
		ArrayList<GenericLine> lineList = new ArrayList<GenericLine>();
		if (showMedSD[0] || showMedSD[1]) {
		    xMed = columnsChangedX || dataChanged || Double.isNaN(xMed) ? Array.median(xData) : xMed;
		    xMin = columnsChangedX || dataChanged || Double.isNaN(xMin) ? Math.min(Math.min(0, plotXmin) - xMed, Array.min(xData) - xMed) : xMin;
		    xMax = columnsChangedX || dataChanged || Double.isNaN(xMax) ? Math.max(this.plotXmax + xMed, Array.max(xData) + xMed) : xMax;
		    if (showMedSD[0]) {
		        lineList.add(new GenericLine((float)xMed, (float)xMin, (float)xMed, (float)xMax, (byte)1, (byte) 8, (byte)1, 0, false));  
		    }
		    if (showMedSD[1]) {
    		    xSD = columnsChangedX || dataChanged || Double.isNaN(xSD) ? Array.stdev(xData, false) : xSD;
    		    lineList.add(new GenericLine((float)(xMed - xSD), (float)xMin, (float)(xMed - xSD), (float)xMax, (byte)1, (byte) 9, (byte)1, 0, false));  
    		    lineList.add(new GenericLine((float)(xMed + xSD), (float)xMin, (float)(xMed + xSD), (float)xMax, (byte)1, (byte) 9, (byte)1, 0, false));  
		    }
		}
		if (showMedSD[2] || showMedSD[3]) {
		    yMed = columnsChangedY || dataChanged || Double.isNaN(yMed) ? Array.median(yData) : yMed;
		    yMin = columnsChangedY || dataChanged || Double.isNaN(yMax) ? Math.min(Math.min(0, plotYmin) - yMed, Array.min(yData)) : yMin;
		    yMax = columnsChangedY || dataChanged || Double.isNaN(yMin) ? Math.max(this.plotYmax + yMed, Array.max(yData)) : yMax;
		    if (showMedSD[2]) {
		        lineList.add(new GenericLine((float)yMin, (float)yMed, (float)yMax, (float)yMed, (byte)1, (byte) 8, (byte)1, 0, false));  
		    }
		    if (showMedSD[3]) {
		        ySD = columnsChangedY || dataChanged || Double.isNaN(ySD) ? Array.stdev(yData, false) : ySD;
		        lineList.add(new GenericLine((float)yMin, (float)(yMed - ySD), (float)yMax, (float)(yMed - ySD), (byte)1, (byte) 9, (byte)1, 0, false));  
		        lineList.add(new GenericLine((float)yMin, (float)(yMed + ySD), (float)yMax, (float)(yMed + ySD), (byte)1, (byte) 9, (byte)1, 0, false));  
		    }
		    
		}
		
        lines = lineList.toArray(new GenericLine[lineList.size()]);
        lineList = null;
        
//        ArrayList<boolean[]> gates = new ArrayList<boolean[]>();
//        for (int i = 0; i < polygons.length; i++) {
//            PolygonGate pg = new PolygonGate();
//            pg.addDimension(new GateDimension(xCol));
//            pg.addDimension(new GateDimension(yCol));
//            pg.setPath(polygons[i].myPath);
//            gates.add(pg.gate(fcp.dataLoader));
//        }
//        for (int i = 0; i < rectangles.length; i++) {
//            RectangleGate rg = new RectangleGate();
//            rg.addDimension(new GateDimension.RectangleGateDimension(xCol, Math.min(rectangles[i].getStartXValue(), rectangles[i].getStopXValue()), Math.max(rectangles[i].getStartXValue(), rectangles[i].getStopXValue())));
//            rg.addDimension(new GateDimension.RectangleGateDimension(yCol, Math.min(rectangles[i].getStartYValue(), rectangles[i].getStopYValue()), Math.max(rectangles[i].getStartYValue(), rectangles[i].getStopYValue())));
//            gates.add(rg.gate(fcp.dataLoader));
//        }
        
        byte color = 0;
        if (columnsChangedX || columnsChangedY || dataChanged || gatesChanged) {
    		points = new PlotPoint[dataCount];
    		for (int i = 0; i < points.length; i++) {
    			xAxisValue = (float) xData[i];
    			yAxisValue = (float) yData[i];
    			if (Float.isNaN(xAxisValue) || Float.isNaN(yAxisValue)) {
    				type = PlotPoint.NOT_A_NUMBER;
    			} else {
    				type = PlotPoint.FILLED_CIRCLE;
    			}
    			
    			color = (byte) 0; // TODO apply gating for colors
    			for (int g = 0; g < gates.size(); g++) {
    			    if (gates.get(g)[i]) {
    			        color = (byte) 3;
    			    }
    			}
    			points[i] = new PlotPoint(i + "", type, xAxisValue, yAxisValue, size, color, (byte)0);
    		}
//            rects.clear();
//            ellipses.clear();
            
            rectangles = rects.toArray(new GenericRectangle[rects.size()]);
            polygons = polys.toArray(new GenericPath[polys.size()]);
        }
	}
	
    public void mouseMoved(MouseEvent e) {
        super.mouseMoved(e);
        int tempX = e.getX();
        int tempY = e.getY();
        for (int i = 0; i < rects.size(); i++) {
            GenericRectangle rect = rects.get(i);
            int xPix1, xPix2, yPix1, yPix2;
            xPix1 = getXPixel(rect.getStartXValue());
            xPix2 = getXPixel(rect.getStopXValue());
            yPix1 = getYPixel(rect.getStartYValue());
            yPix2 = getYPixel(rect.getStopYValue());
            Rectangle myRect = new Rectangle(Math.min(xPix1, xPix2) - 10, 
                                            Math.min(yPix1, yPix2) - 10, 
                                            Math.max(xPix1, xPix2) 
                                                - Math.min(xPix1, xPix2) + 10, 
                                            Math.max(yPix1, yPix2)
                                                - Math.min(yPix1, yPix2) + 10);
            rect.setEditable(myRect.contains(tempX, tempY));
        }
        for (int i = 0; i < polys.size(); i++) {
            GenericPath poly = polys.get(i);
            Path2D path = ((Path2D)poly.myPath.clone());
            poly.setEditable(path.contains(getXValueFromXPixel(tempX), getYValueFromYPixel(tempY)));
        }
        paintAgain();
        
    }
    
	private void leftMousePressedRect(MouseEvent e) {
        int tempX = e.getX();
        int tempY = e.getY();
        int toRemove = -1;

        for (int i = 0; i < rects.size(); i++) {
            GenericRectangle rect = rects.get(i);
            boolean closeToStartX = Math.abs(getXPixel(rect.getStartXValue()) - tempX) < 4; 
            boolean closeToStartY = Math.abs(getYPixel(rect.getStartYValue()) - tempY) < 4; 
            boolean closeToStopX = Math.abs(getXPixel(rect.getStopXValue()) - tempX) < 4; 
            boolean closeToStopY = Math.abs(getYPixel(rect.getStopYValue()) - tempY) < 4; 
            
            double tempStartX = Double.NaN, tempStartY = Double.NaN, tempStopX = Double.NaN, tempStopY = Double.NaN;
            boolean drag = true;
            if (closeToStartX && closeToStartY) {
                tempStartX = rect.getStopXValue();
                tempStartY = rect.getStopYValue();
                tempStopX = rect.getStartXValue();
                tempStopY = rect.getStartYValue();
                drag = false;
            } else if (closeToStartX && closeToStopY) {
                tempStartX = rect.getStopXValue();
                tempStartY = rect.getStartYValue();
                tempStopX = rect.getStartXValue();
                tempStopY = rect.getStopYValue();
                drag = false;
            } else if (closeToStopX && closeToStartY) {
                tempStartX = rect.getStartXValue();
                tempStartY = rect.getStopYValue();
                tempStopX = rect.getStopXValue();
                tempStopY = rect.getStartYValue();
                drag = false;
            } else if (closeToStopX && closeToStopY) {
                tempStartX = rect.getStartXValue();
                tempStartY = rect.getStartYValue();
                tempStopX = rect.getStopXValue();
                tempStopY = rect.getStopYValue();
                drag = false;
            }
            if (!Double.isNaN(tempStartX)) {
                startX = getXPixel(tempStartX);
                startY = getYPixel(tempStartY);
                highlightRectangle = new GenericRectangle((float)tempStartX, 
                                                            (float)tempStartY,
                                                            (float)tempStopX, 
                                                            (float)tempStopY,
                                                            (byte)1, false, false, (byte)0, (byte)99, true);
                toRemove = i;
                break;
            } else {
                Rectangle2D myRect = new Rectangle2D.Float(Math.min(rect.getStartXValue(), rect.getStopXValue()), Math.min(rect.getStartYValue(), rect.getStopYValue()), Math.max(rect.getStartXValue(), rect.getStopXValue()) - Math.min(rect.getStartXValue(), rect.getStopXValue()), Math.max(rect.getStartYValue(), rect.getStopYValue()) - Math.min(rect.getStartYValue(), rect.getStopYValue()));
                if (myRect.contains(getXValueFromXPixel(tempX), getYValueFromYPixel(tempY)) && drag) {
                    dragInd = i;
                    startX = getXPixel(tempStartX);
                    startY = getYPixel(tempStartY);
                }
            }
        }
        if (toRemove >= 0) {
            rects.remove(toRemove);
        } else {
            startX = e.getX();
            startY = e.getY();
        }
        paintAgain();
	}
	
	private void leftMousePressedPoly(MouseEvent e) {
	    int tempX = e.getX();
	    int tempY = e.getY();
	    
	    double[] coords = new double[6];
	    int polyInd = -1;
	    int vertInd = -1;
	    for (int i = 0; i < polys.size(); i++) {
	    	Path2D path = polys.get(i).myPath;
	    	if (path.contains(getXValueFromXPixel(tempX), getYValueFromYPixel(tempY))) {
	    	    polyInd = i;
	    	}
	    	PathIterator pi = path.getPathIterator(null);
	    	int vInd = 0;
	    	while (!pi.isDone()) {
	    	    pi.currentSegment(coords);
	    	    if (Math.abs(tempX - getXPixel(coords[0])) < 4 && Math.abs(tempY - getYPixel(coords[1])) < 4) {
	    	        polyInd = i;
	    	        vertInd = vInd;
	    	        break;
	    	    }
	    	    vInd++;
	    	    pi.next();
	    	}
	    }
	    if (polyInd >= 0) {
	        dragInd = polyInd;
	    }
	    if (vertInd >= 0) {
	        polyDragVertInd = vertInd;
	    }
        startX = e.getX();
        startY = e.getY();
	    paintAgain();
	}
	
	public void mousePressed(MouseEvent e) {
        if (SwingUtilities.isLeftMouseButton(e) && !e.isControlDown()) {
            
            if (currentTool == RECT_TOOL) {
                leftMousePressedRect(e);
            } else if (currentTool == POLY_TOOL) {
            	leftMousePressedPoly(e);
            }
            
            
        } else {
            super.mousePressed(e);
        }
    }

	private void rightMouseClickedRect(MouseEvent e) {
        int tempX = e.getX();
        int tempY = e.getY();
        int toRemove = -1;

        for (int i = rects.size() - 1; i >= 0; i--) {
            GenericRectangle rect = rects.get(i);
            double xLow, xHigh, yLow, yHigh;
            xLow = getXPixel(Math.min(rect.getStartXValue(), rect.getStopXValue()));
            xHigh = getXPixel(Math.max(rect.getStartXValue(), rect.getStopXValue()));
            yHigh = getYPixel(Math.min(rect.getStartYValue(), rect.getStopYValue()));
            yLow = getYPixel(Math.max(rect.getStartYValue(), rect.getStopYValue()));
            if (isHistogram()) {
                if (xLow <= tempX && xHigh >= tempX) {
                    toRemove = i;
                    break;
                }
            } else {
                if (xLow <= tempX && xHigh >= tempX && yLow <= tempY && yHigh >= tempY) {
                    toRemove = i;
                    break;
                }
            }
        }
        if (toRemove != -1) {
            rects.remove(toRemove);
        }
        // TODO remove gate in gating data struct
        // TODO confirm gate removal
        super.mouseClicked(e);
        paintAgain();
	}
	
	private void rightMouseClickedPoly(MouseEvent e) {
		int tempX = e.getX();
		int tempY = e.getY();
		int toRemove = -1;
		
		if (tempPoly.size() > 0) {
		    tempPoly.clear();
		    highlightPoly = null;
		    return;
		}
		
		double tempValX = getXValueFromXPixel(tempX);
		double tempValY = getYValueFromYPixel(tempY);
	    for (int i = polys.size() - 1; i >= 0; i--) {
	    	if (polys.get(i).myPath.contains(tempValX, tempValY)) {
	    		toRemove = i;
	    		break;
	    	}
	    }
		if (toRemove != -1) {
			polys.remove(toRemove);
		}

        // TODO remove gate in gating data struct
        // TODO confirm gate removal
		
		super.mouseClicked(e);
		paintAgain();
	}
	
	public void mouseClicked(MouseEvent e) {   
		if (SwingUtilities.isLeftMouseButton(e) && !e.isControlDown()) {
		    if (currentTool == POLY_TOOL) {
		        
		    }
			// 
		} else {
		    
		    if (currentTool == RECT_TOOL) {
		        rightMouseClickedRect(e);
		    } else if (currentTool == POLY_TOOL) {
		        rightMouseClickedPoly(e);
		    }
		    
		}
		paintAgain();
	}
    
	private void setForceGatesChanged() {
	    this.forceGatesChanged = true;
	}
	
	private void leftMouseReleasedPoly(MouseEvent e) {
	    if (dragInd >= 0) {
	        dragInd = -1;
	        polyDragVertInd = -1;
	        setForceGatesChanged();
	        paintAgain();
	        return;
	    }
	    
        int tempPxX = e.getX();
        int tempPxY = e.getY();
        double tempValX = getXValueFromXPixel(tempPxX);
        double tempValY = getYValueFromYPixel(tempPxY);
        
        if (!tempPoly.isEmpty()) {
            int initPtX = getXPixel(tempPoly.get(0)[0]);
            int initPtY = getYPixel(tempPoly.get(0)[1]);
            if (Math.abs(initPtX - tempPxX) < 5 && Math.abs(initPtY - tempPxY) < 5) {
                Path2D path = new Path2D.Double(Path2D.WIND_EVEN_ODD);
                path.moveTo(tempPoly.get(0)[0], tempPoly.get(0)[1]);
                for(int i = 1; i < tempPoly.size(); ++i) {
                   path.lineTo(tempPoly.get(i)[0], tempPoly.get(i)[1]);
                }
                path.closePath();
                PolygonGate pg = new PolygonGate();
                pg.addDimension(new GateDimension(xCol));
                pg.addDimension(new GateDimension(yCol));
                pg.setPath(path);
//                fcp.addGate(pg);
//                 TODO add new PolygonGate to fcp gating strategy, instead of adding to polys
                polys.add(new GenericPath(path, (byte)0, (byte)0, (byte)99, false, false));
                tempPoly.clear();
                highlightPoly = null;
                paintAgain();
                return;
            }
        } 
        tempPoly.add(new double[]{tempValX, tempValY});
        Path2D path = new Path2D.Double(Path2D.WIND_EVEN_ODD);
        path.moveTo(tempPoly.get(0)[0], tempPoly.get(0)[1]);
        for(int i = 1; i < tempPoly.size(); ++i) {
           path.lineTo(tempPoly.get(i)[0], tempPoly.get(i)[1]);
        }
        highlightPoly = new GenericPath(path, (byte)0, (byte)0, (byte)99, false, true);
        setForceGatesChanged();
        paintAgain();
	}
	
	private void leftMouseReleasedRect(MouseEvent e) {
        if (dragInd >= 0) {
            dragInd = -1;
            setForceGatesChanged();
            paintAgain();
            return;
        }
        int mouseEndX;
        int mouseEndY;
        mouseEndX = e.getX();
        mouseEndY = e.getY();
        highlightRectangle = null;
        // TODO add new RectangleGate to fcp gating strategy, instead of adding to rects
        RectangleGate rg = new RectangleGate();
        rg.addDimension(new GateDimension.RectangleGateDimension(xCol, (float)getXValueFromXPixel(startX), (float)getXValueFromXPixel(mouseEndX)));
        if (!isHistogram()) {
            rg.addDimension(new GateDimension.RectangleGateDimension(yCol, (float)getYValueFromYPixel(startY), (float)getYValueFromYPixel(mouseEndY)));
        }
        // for testing:
        rects.add(new GenericRectangle((float)getXValueFromXPixel(startX), 
                (float)getYValueFromYPixel(startY), (float)getXValueFromXPixel(mouseEndX), (float)getYValueFromYPixel(mouseEndY), (byte)1, false, false, (byte)0, (byte)0, false));
        // TODO fcp.addGate(rg);
        paintAgain();
	}

	public void mouseReleased(MouseEvent e) {
        if (SwingUtilities.isLeftMouseButton(e) && !e.isControlDown()) {
            
            if (currentTool == RECT_TOOL) {
                leftMouseReleasedRect(e);
            } else if (currentTool == POLY_TOOL) {
                leftMouseReleasedPoly(e);
            }
            
        } else {
            super.mouseReleased(e);
        }
    }
    
	private void leftMouseDraggedRect(MouseEvent e) {
        int mouseEndX;
        int mouseEndY;
        mouseEndX = e.getX();
        mouseEndY = e.getY();
        if (dragInd >= 0) {
            float dx, dy;
            dx = (float) (getXValueFromXPixel(mouseEndX) - getXValueFromXPixel(startX));
            dy = (float) (getYValueFromYPixel(mouseEndY) - getYValueFromYPixel(startY));
            GenericRectangle gr = rects.get(dragInd);
            rects.set(dragInd, new GenericRectangle(gr.getStartXValue() + dx, gr.getStartYValue() + dy, gr.getStopXValue() + dx, gr.getStopYValue() + dy, gr.getThickness(), gr.getFill(), gr.getRoundedCorners(), gr.getColor(), gr.getFillColor(), gr.getLayer(), gr.getEditable()));
            startX = mouseEndX;
            startY = mouseEndY;
            setForceGatesChanged();
        } else {
            highlightRectangle = new GenericRectangle(
                    (float)getXValueFromXPixel(startX), 
                    (float)getYValueFromYPixel(startY), 
                    (float)getXValueFromXPixel(mouseEndX), 
                    isHistogram() ? (float)getYValueFromYPixel(startY) : (float)getYValueFromYPixel(mouseEndY), 
                    (byte)1, false, false, (byte)0, (byte)99, true);
        }
        paintAgain();
    }
    
    private void leftMouseDraggedPoly(MouseEvent e) {
        if (dragInd >= 0) {
            int mouseEndX;
            int mouseEndY;
            mouseEndX = e.getX();
            mouseEndY = e.getY();
            
            // TODO don't rely on polys here, instead use fcp.gating object
            GenericPath gp = polys.get(dragInd);
            if (polyDragVertInd >= 0) {
                Path2D newPath = new Path2D.Double();
                PathIterator pi = gp.myPath.getPathIterator(null);
                double[] coords = new double[6];
                int v = 0;
                while(!pi.isDone()) {
                    int code = pi.currentSegment(coords);
                    if (v == polyDragVertInd) {
                        coords[0] = getXValueFromXPixel(mouseEndX);
                        coords[1] = getYValueFromYPixel(mouseEndY);
                    }
                    switch(code) {
                        case PathIterator.SEG_MOVETO:
                            newPath.moveTo(coords[0], coords[1]);
                            break;
                        case PathIterator.SEG_LINETO:
                            newPath.lineTo(coords[0], coords[1]);
                            break;
                        case PathIterator.SEG_CLOSE:
                            newPath.closePath();
                    }
                    v++;
                    pi.next();
                }
                polys.get(dragInd).myPath = newPath;
            } else {
                double dx, dy;
                dx = getXValueFromXPixel(mouseEndX) - getXValueFromXPixel(startX);
                dy = getYValueFromYPixel(mouseEndY) - getYValueFromYPixel(startY);
                AffineTransform at = AffineTransform.getTranslateInstance(dx, dy);
                gp.myPath.transform(at);
                startX = mouseEndX;
                startY = mouseEndY;
                // TODO update gating??
            }
            setForceGatesChanged();
            paintAgain();
        }
    }
    
    public void mouseDragged(MouseEvent e) {
        if (SwingUtilities.isLeftMouseButton(e) && !e.isControlDown()) {
            if (currentTool == RECT_TOOL) {
                leftMouseDraggedRect(e);
            } else if (currentTool == POLY_TOOL) {
                leftMouseDraggedPoly(e);
            }
        } else {
            super.mouseDragged(e);
        }
    }
    

    @Override
    public void highlightPoints() {
        return;
    }

    @Override
    public void assignAxisLabels() {
        xAxisLabel = fcp.getXDataName();
        yAxisLabel = fcp.getYDataName();
    }
}
