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
import java.util.HashSet;

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

public class FCSPanel extends AbstractPanel2 implements MouseListener, MouseMotionListener {
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
	
	public FCSPanel(FCSPlot fcsPlot) {
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
	
    public static enum GATING_TOOL {
        RECT_TOOL,
        POLY_TOOL;
    };
	
//	private volatile GATING_TOOL currentTool = GATING_TOOL.POLY_TOOL;
	private volatile GATING_TOOL currentTool = GATING_TOOL.RECT_TOOL;
	
//	private int dragInd = -1;
//    private int polyDragVertInd = -1;

    public void setGatingTool(GATING_TOOL tool) {
        // TODO cancel drawing any gates
        currentTool = tool;
    }

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
		
		if (fcp.dataLoader == null) {
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
            
    		ArrayList<GenericRectangle> temp = new ArrayList<GenericRectangle>();
    		temp.addAll(rects);
    		temp.addAll(draggingVertexRects);
            rectangles = temp.toArray(new GenericRectangle[temp.size()]);
            polygons = polys.toArray(new GenericPath[polys.size()]);
        }
	}
	
    
    private ArrayList<GenericPath> getPolysWithVerticesNearClick(MouseEvent e) {
        int tempX = e.getX();
        int tempY = e.getY();
        
        ArrayList<GenericPath> retPolys = new ArrayList<GenericPath>();
        
        double[] coords = new double[6];
        for (int i = 0; i < polys.size(); i++) {
            Path2D path = polys.get(i).myPath;
            PathIterator pi = path.getPathIterator(null);
            while (!pi.isDone()) {
                pi.currentSegment(coords);
                if (Math.abs(tempX - getXPixel(coords[0])) < DEFAULT_NEARBY_DIST  && Math.abs(tempY - getYPixel(coords[1])) < DEFAULT_NEARBY_DIST) {
                    retPolys.add(polys.get(i));
                    break;
                }
                pi.next();
            }
        }
        
        return retPolys;
    }
    
    private ArrayList<GenericPath> getPolysContainingClick(MouseEvent e) {
        int tempX = e.getX();
        int tempY = e.getY();
        
        ArrayList<GenericPath> retPolys = new ArrayList<GenericPath>();
        
        for (int i = 0; i < polys.size(); i++) {
            Path2D path = polys.get(i).myPath;
            if (path.contains(getXValueFromXPixel(tempX), getYValueFromYPixel(tempY))) {
                retPolys.add(polys.get(i));
            }
        }
        
        return retPolys;
    }
    
    private ArrayList<GenericRectangle> getRectsWithVerticesNearClick(MouseEvent e) {
        int tempX = e.getX();
        int tempY = e.getY();
    
        ArrayList<GenericRectangle> retRects = new ArrayList<GenericRectangle>();
        for (int i = 0; i < rects.size(); i++) {
            GenericRectangle rect = rects.get(i);
            boolean aX = Math.abs(getXPixel(rect.getStartXValue()) - tempX) < DEFAULT_NEARBY_DIST; 
            boolean aY = Math.abs(getYPixel(rect.getStartYValue()) - tempY) < DEFAULT_NEARBY_DIST; 
            boolean bX = Math.abs(getXPixel(rect.getStopXValue()) - tempX) < DEFAULT_NEARBY_DIST; 
            boolean bY = Math.abs(getYPixel(rect.getStopYValue()) - tempY) < DEFAULT_NEARBY_DIST; 
            if (aX && aY || aX && bY || bX && aY || bX && bY) {
                retRects.add(rect);
            }
        }
        return retRects;
    }
    
    private ArrayList<GenericRectangle> getRectsContainingClick(MouseEvent e) {
        int tempX = e.getX();
        int tempY = e.getY();
        
        ArrayList<GenericRectangle> retRects = new ArrayList<GenericRectangle>();
        
        for (int i = rects.size() - 1; i >= 0; i--) {
            GenericRectangle rect = rects.get(i);
            double xLow, xHigh, yLow, yHigh;
            xLow = getXPixel(Math.min(rect.getStartXValue(), rect.getStopXValue()));
            xHigh = getXPixel(Math.max(rect.getStartXValue(), rect.getStopXValue()));
            if (xLow <= tempX && xHigh >= tempX) {
                if (!isHistogram()) {
                    yHigh = getYPixel(Math.min(rect.getStartYValue(), rect.getStopYValue()));
                    yLow = getYPixel(Math.max(rect.getStartYValue(), rect.getStopYValue()));
                    if (yLow <= tempY && yHigh >= tempY) {
                        retRects.add(rect);
                    }
                } else {
                    retRects.add(rect);
                }
            }
        }
        
        return retRects;
    }
    

    
    
    
    /*

    - starting a new gate clears selections
    - clicking a selected gate deselects it
    - press and drag on a selected gate drags all selected gates
    - event order: press, drag, click(?), release(?) (Or release, click?)
    - drag flag? (InDrag)
    - if release and not in drag, and not creating new gate, either select (actual logic in mousePressed, not mouseReleased) 
                    or deselect (which has to be in mouseReleased cuz drag) single gate if inside, or deselect all if outside
    - if release and inDrag, finalize drag and clear flag (either create rect gate


    //  mousePressed():
    // if drawing new poly, click == new point (i.e., set startX, startY, return) 
     * else
    // check mouse location vs all polys
    //     if not within any shapes, deselect all
    //     -   and start drawing next shape (i.e. set startX, startY) [more important in mouseReleased]
    //     else if within shape, select shape (begin drag)
    //     -   followed by either mouseDragged or mouseReleased
    //     -   mouseDragged == drag all selected polys/rects and mouse polys/rects
    //     -   mouseReleased == add poly to selected polys if still within poly and not dragging
    //     if within N pixels of vertex, start vertex drag (don't deselect everything; can drag/edit while preserving multiselect)
    //      
     */
    
    private static final int DEFAULT_NEARBY_DIST = 4;
    private HashSet<GenericPath> selectedPolys = new HashSet<GenericPath>();
    private HashSet<GenericRectangle> selectedRects = new HashSet<GenericRectangle>();

    // in mouseReleased, if !drag && tempPoly.isEmpty (and not creating a rect gate), XAND mousePoly/Rect from selectedPoly/Rect (remove if in selected, add if not) 
    private HashSet<GenericPath> mousePolys = new HashSet<GenericPath>();
    private HashSet<GenericRectangle> mouseRects = new HashSet<GenericRectangle>();
    
    private ArrayList<GenericRectangle> draggingVertexRects = new ArrayList<GenericRectangle>();
    private ArrayList<GenericPath> draggingPolys = new ArrayList<GenericPath>();
    private ArrayList<Integer> draggingPolyInds = new ArrayList<Integer>();
    
    public void mouseReleased(MouseEvent e) {
        if (SwingUtilities.isLeftMouseButton(e) && !e.isControlDown()) {
            int mouseEndX;
            int mouseEndY;
            mouseEndX = e.getX();
            mouseEndY = e.getY();
            
            boolean didSelect = mouseRects.size() > 0 || mousePolys.size() > 0;
            boolean didClear = !didSelect && (selectedRects.size() > 0 || selectedPolys.size() > 0); 
            boolean wasDrag = drag;
            if (!drag/* && (currentTool == GATING_TOOL.RECT_TOOL || tempPoly.isEmpty())*/) {
                if (!didSelect) {
                    for (GenericRectangle rect : selectedRects) {
                        rect.setEditable(false);
                    }
                    for (GenericPath poly : selectedPolys) {
                        poly.setEditable(false);
                    }
                    selectedRects.clear();
                    selectedPolys.clear();
//                    for (GenericRectangle rect : mouseRects) {
//                        rect.setEditable(false);
//                    }
//                    for (GenericPath poly : mousePolys) {
//                        poly.setEditable(false);
//                    }
//                    mouseRects.clear();
//                    mousePolys.clear();
                } else {
                    for (GenericRectangle rect : mouseRects) {
                        if (selectedRects.contains(rect)) {
                            rect.setEditable(false);
                            selectedRects.remove(rect);
                        } else {
                            rect.setEditable(true);
                            selectedRects.add(rect);
                        }
                    }
                    mouseRects.clear();
                    for (GenericPath path : mousePolys) {
                        if (selectedPolys.contains(path)) {
                            path.setEditable(false);
                            selectedPolys.remove(path);
                        } else {
                            path.setEditable(true);
                            selectedPolys.add(path);
                        }
                    }
                    mousePolys.clear();
                }
            } else {
                drag = false;
                for (GenericRectangle rect : mouseRects) {
                    if (selectedRects.contains(rect)) {
                    } else {
                        rect.setEditable(true);
                        selectedRects.add(rect);
                    }
                }
                mouseRects.clear();
                for (GenericPath path : mousePolys) {
                    if (selectedPolys.contains(path)) {
                    } else {
                        path.setEditable(true);
                        selectedPolys.add(path);
                    }
                }
                mousePolys.clear();
                if (!draggingPolys.isEmpty()) {
                    draggingPolys.clear();  
                    draggingPolyInds.clear();
                }
                if (draggingVertexRects.isEmpty()) {
                    if (currentTool == GATING_TOOL.RECT_TOOL) {
                        highlightRectangle = null;
        //                    RectangleGate rg = new RectangleGate();
        //                    rg.addDimension(new GateDimension.RectangleGateDimension(xCol, (float)getXValueFromXPixel(startX), (float)getXValueFromXPixel(mouseEndX)));
        //                    if (!isHistogram()) {
        //                        rg.addDimension(new GateDimension.RectangleGateDimension(yCol, (float)getYValueFromYPixel(startY), (float)getYValueFromYPixel(mouseEndY)));
        //                    }
        //                    fcp.addGate(rg);
                        
                        // TODO for testing:
                        if (Math.abs(mouseEndX - startX) > DEFAULT_NEARBY_DIST) {
                            if (isHistogram() || (Math.abs(mouseEndY - startY) > DEFAULT_NEARBY_DIST)) {
                                GenericRectangle rect = new GenericRectangle((float)getXValueFromXPixel(startX), 
                                                            (float)getYValueFromYPixel(startY), 
                                                            (float)getXValueFromXPixel(mouseEndX), 
                                                            (float)getYValueFromYPixel(mouseEndY), 
                                                            (byte)1, false, false, (byte)0, (byte)0, true);
                                rects.add(rect);
                                selectedRects.add(rect);
                            }
                        }
                    }
                } else {
                    rects.addAll(draggingVertexRects);
                    selectedRects.addAll(draggingVertexRects);
                    draggingVertexRects.clear();
                }
            }
            if (!didSelect && !didClear && !wasDrag) {
                if (currentTool == GATING_TOOL.POLY_TOOL) {
                    double tempValX = getXValueFromXPixel(mouseEndX);
                    double tempValY = getYValueFromYPixel(mouseEndY);
                    if (!tempPoly.isEmpty()) {
                        int initPtX = getXPixel(tempPoly.get(0)[0]);
                        int initPtY = getYPixel(tempPoly.get(0)[1]);
                        if (Math.abs(initPtX - mouseEndX) < 5 && Math.abs(initPtY - mouseEndY) < 5) {
                            if (tempPoly.size() > 2) {
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
        //                        fcp.addGate(pg);
        //                         TODO add new PolygonGate to fcp gating strategy, instead of adding to polys
                                GenericPath newPoly = new GenericPath(path, (byte)0, (byte)0, (byte)99, false, true);
                                polys.add(newPoly);
                                selectedPolys.add(newPoly);
                            }
                            tempPoly.clear();
                            highlightPoly = null;
                            setForceGatesChanged();
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
                }
                
            }
            setForceGatesChanged();
            paintAgain();
            
            
        } else {
            super.mouseReleased(e);
        }
    }
    
    public void mousePressed(MouseEvent e) {
        if (SwingUtilities.isLeftMouseButton(e) && !e.isControlDown()) {
            int tempX = e.getX();
            int tempY = e.getY();
            if (tempPoly.isEmpty()) {
                ArrayList<GenericPath> insidePolys = getPolysContainingClick(e);
                ArrayList<GenericPath> closePolys = getPolysWithVerticesNearClick(e);
                ArrayList<GenericRectangle> insideRects = getRectsContainingClick(e);
                ArrayList<GenericRectangle> closeRects = getRectsWithVerticesNearClick(e);
                
//                if (insidePolys.isEmpty() && insideRects.isEmpty() && closePolys.isEmpty() && closeRects.isEmpty()) {
//                    selectedPolys.clear();
//                    selectedRects.clear();
//                }
                if (closeRects.isEmpty()) {
                    if (!insideRects.isEmpty()) {
                        for (GenericRectangle gr : insideRects) {
                            gr.setEditable(true);
                            if (!selectedRects.contains(gr)) {
                                mouseRects.add(gr);
                            }
                        }
                    }
                } else {
                    for (GenericRectangle rect : closeRects) {
                        rects.remove(rect);
                        selectedRects.remove(rect);
                        mouseRects.remove(rect);

                        boolean closeToStartX = Math.abs(getXPixel(rect.getStartXValue()) - tempX) < 4; 
                        boolean closeToStartY = Math.abs(getYPixel(rect.getStartYValue()) - tempY) < 4; 
                        boolean closeToStopX = Math.abs(getXPixel(rect.getStopXValue()) - tempX) < 4; 
                        boolean closeToStopY = Math.abs(getYPixel(rect.getStopYValue()) - tempY) < 4; 
                        
                        double tempStartX = Double.NaN, tempStartY = Double.NaN, tempStopX = Double.NaN, tempStopY = Double.NaN;
                        if (closeToStartX && closeToStartY) {
                            tempStartX = rect.getStopXValue();
                            tempStartY = rect.getStopYValue();
                            tempStopX = rect.getStartXValue();
                            tempStopY = rect.getStartYValue();
                        } else if (closeToStartX && closeToStopY) {
                            tempStartX = rect.getStopXValue();
                            tempStartY = rect.getStartYValue();
                            tempStopX = rect.getStartXValue();
                            tempStopY = rect.getStopYValue();
                        } else if (closeToStopX && closeToStartY) {
                            tempStartX = rect.getStartXValue();
                            tempStartY = rect.getStopYValue();
                            tempStopX = rect.getStopXValue();
                            tempStopY = rect.getStartYValue();
                        } else if (closeToStopX && closeToStopY) {
                            tempStartX = rect.getStartXValue();
                            tempStartY = rect.getStartYValue();
                            tempStopX = rect.getStopXValue();
                            tempStopY = rect.getStopYValue();
                        }
                        if (!Double.isNaN(tempStartX)) {
                            GenericRectangle dragRect = new GenericRectangle((float)tempStartX, 
                                                                        (float)tempStartY,
                                                                        (float)tempStopX, 
                                                                        (float)tempStopY,
                                                                        (byte)1, false, false, (byte)0, (byte)99, true);
                            draggingVertexRects.add(dragRect);
                        }
                    }
                }
                if (closePolys.isEmpty()) {
                    // TODO should this affect ALL of these shapes?
                    if (!insidePolys.isEmpty()) {
                        for (GenericPath gp : insidePolys) {
                            gp.setEditable(true);
                            mousePolys.add(gp);
                        }
                    }
                } else {
//                    polys.removeAll(closePolys);
//                    selectedPolys.removeAll(closePolys);
                    selectedPolys.addAll(closePolys);
                    mousePolys.removeAll(closePolys);
                    
                    draggingPolys = closePolys;
                    draggingPolyInds.clear();
                    
                    double[] coords = new double[6];
                    polyLoop : for (int i = 0; i < closePolys.size(); i++) {
                        Path2D path = closePolys.get(i).myPath;
                        PathIterator pi = path.getPathIterator(null);
                        int vInd = 0;
                        while (!pi.isDone()) {
                            pi.currentSegment(coords);
                            if (Math.abs(tempX - getXPixel(coords[0])) < 4 && Math.abs(tempY - getYPixel(coords[1])) < 4) {
                                draggingPolyInds.add(vInd);
                                continue polyLoop;
                            }
                            vInd++;
                            pi.next();
                        }
                        draggingPolyInds.add(-1); // this shouldn't technically happen, since closePolys only has polys with at least one vertex close to the mouse
                    }
                    
                }
            } else {
                // do nothing, wait for mouseReleased to set new point or create poly
            }
            startX = e.getX();
            startY = e.getY();
            paintAgain();
        } else {
            super.mousePressed(e);
        }
    }
	
	

//	public void mouseMoved(MouseEvent e) {
//        super.mouseMoved(e);
//        int tempX = e.getX();
//        int tempY = e.getY();
//        for (int i = 0; i < rects.size(); i++) {
//            GenericRectangle rect = rects.get(i);
//            int xPix1, xPix2, yPix1, yPix2;
//            xPix1 = getXPixel(rect.getStartXValue());
//            xPix2 = getXPixel(rect.getStopXValue());
//            yPix1 = getYPixel(rect.getStartYValue());
//            yPix2 = getYPixel(rect.getStopYValue());
//            Rectangle myRect = new Rectangle(Math.min(xPix1, xPix2) - 10, 
//                                            Math.min(yPix1, yPix2) - 10, 
//                                            Math.max(xPix1, xPix2) 
//                                                - Math.min(xPix1, xPix2) + 10, 
//                                            Math.max(yPix1, yPix2)
//                                                - Math.min(yPix1, yPix2) + 10);
//            rect.setEditable(myRect.contains(tempX, tempY));
//        }
//        for (int i = 0; i < polys.size(); i++) {
//            GenericPath poly = polys.get(i);
//            Path2D path = ((Path2D)poly.myPath.clone());
//            poly.setEditable(path.contains(getXValueFromXPixel(tempX), getYValueFromYPixel(tempY)));
//        }
//        paintAgain();
//    }

    public void mouseClicked(MouseEvent e) {   
    	if (SwingUtilities.isLeftMouseButton(e) && !e.isControlDown()) {
    	    if (currentTool == GATING_TOOL.POLY_TOOL) {
    	        
    	    }
    		// 
    	} else {
    
            // check mouse location vs all shapes
            //     else if within shape, delete shape (and gate attached - CONFIRM DELETE)
    	    
    	    if (currentTool == GATING_TOOL.RECT_TOOL) {
    	        rightMouseClickedRect(e);
    	    } else if (currentTool == GATING_TOOL.POLY_TOOL) {
    	        rightMouseClickedPoly(e);
    	    }
    	    
    	}
    	paintAgain();
    }

    public void mouseDragged(MouseEvent e) {
        if (SwingUtilities.isLeftMouseButton(e) && !e.isControlDown()) {
            drag = true;
            
            boolean updateStartCoords = false;
            boolean updateForceGating = false;
            
            int mouseEndX;
            int mouseEndY;
            mouseEndX = e.getX();
            mouseEndY = e.getY();
            
            if (draggingVertexRects.isEmpty() && draggingPolys.isEmpty()) {
                if (selectedRects.isEmpty() && mouseRects.isEmpty() && selectedPolys.isEmpty() && mousePolys.isEmpty()) {
                    if (currentTool == GATING_TOOL.RECT_TOOL) {
                        highlightRectangle = new GenericRectangle((float) getXValueFromXPixel(startX), 
                                                                    (float) getYValueFromYPixel(startY), 
                                                                    (float) getXValueFromXPixel(mouseEndX), 
                                                                    isHistogram() ? (float) getYValueFromYPixel(startY) : (float) getYValueFromYPixel(mouseEndY), 
                                                                    (byte) 1, false, false, (byte) 0, (byte) 99, true);
                    }
                } else {
                    
                    float dx, dy;
                    dx = (float) (getXValueFromXPixel(mouseEndX) - getXValueFromXPixel(startX));
                    dy = (float) (getYValueFromYPixel(mouseEndY) - getYValueFromYPixel(startY));
                    HashSet<GenericRectangle> newSel = new HashSet<GenericRectangle>();
                    HashSet<GenericRectangle> newMou = new HashSet<GenericRectangle>();
                    HashSet<GenericRectangle> newRec = new HashSet<GenericRectangle>();
                    for (GenericRectangle gr : selectedRects) {
                        rects.remove(gr);
                        GenericRectangle newGr = new GenericRectangle(gr.getStartXValue() + dx, gr.getStartYValue() + dy, gr.getStopXValue() + dx, gr.getStopYValue() + dy, gr.getThickness(), gr.getFill(), gr.getRoundedCorners(), gr.getColor(), gr.getFillColor(), gr.getLayer(), gr.getEditable());
                        if (mouseRects.contains(gr)) {
                            newMou.add(newGr);
                            mouseRects.remove(gr);
                        }
                        newSel.add(newGr);
                        newRec.add(newGr);
                    }
                    selectedRects.clear();
                    selectedRects.addAll(newSel);
                    for (GenericRectangle gr : mouseRects) {
                        rects.remove(gr);
                        GenericRectangle newGr = new GenericRectangle(gr.getStartXValue() + dx, gr.getStartYValue() + dy, gr.getStopXValue() + dx, gr.getStopYValue() + dy, gr.getThickness(), gr.getFill(), gr.getRoundedCorners(), gr.getColor(), gr.getFillColor(), gr.getLayer(), gr.getEditable());
                        newMou.add(newGr);
                        newRec.add(newGr);
                    }
                    mouseRects.clear();
                    mouseRects.addAll(newMou);
                    rects.addAll(newRec);
                    AffineTransform at = AffineTransform.getTranslateInstance(dx, dy);
                    for (GenericPath gp : selectedPolys) {
                        gp.myPath.transform(at);
                    }
                    for (GenericPath gp : mousePolys) {
                        if (!selectedPolys.contains(gp)) {
                            gp.myPath.transform(at);
                        }
                    }
                    updateStartCoords = true;
                    updateForceGating = true;
                }
            }
            if (!draggingVertexRects.isEmpty()) {
                for (int i = 0; i < draggingVertexRects.size(); i++) {
                    float tempStartX, tempStartY;
                    tempStartX = draggingVertexRects.get(i).getStartXValue();
                    tempStartY = draggingVertexRects.get(i).getStartYValue();
                    GenericRectangle gr = new GenericRectangle(tempStartX, 
                            tempStartY, 
                            (float) getXValueFromXPixel(mouseEndX), 
                            isHistogram() ? (float) getYValueFromYPixel(startY) : (float) getYValueFromYPixel(mouseEndY), 
                            (byte) 1, false, false, (byte) 0, (byte) 99, true);
                    draggingVertexRects.set(i, gr);
                }
            }
            if (!draggingPolys.isEmpty()) {
                for (int i = 0; i < draggingPolys.size(); i++) {
                    GenericPath gp = draggingPolys.get(i);
                    int vertexInd = draggingPolyInds.get(i);
                    Path2D newPath = new Path2D.Double();
                    PathIterator pi = gp.myPath.getPathIterator(null);
                    double[] coords = new double[6];
                    int v = 0;
                    while (!pi.isDone()) {
                        int code = pi.currentSegment(coords);
                        if (v == vertexInd) {
                            coords[0] = getXValueFromXPixel(mouseEndX);
                            coords[1] = getYValueFromYPixel(mouseEndY);
                        }
                        switch (code) {
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
                    draggingPolys.get(i).myPath = newPath;
                }
                updateStartCoords = true;
                updateForceGating = true;
            }
            
            if (updateStartCoords) {
                startX = e.getX();
                startY = e.getY();
            }
            if (updateForceGating) {
                setForceGatesChanged();
            }
            paintAgain();
        } else {
            super.mouseDragged(e);
        }
    }

    private void rightMouseClickedRect(MouseEvent e) {
        int tempX = e.getX();
        int tempY = e.getY();
        GenericRectangle toRemove = null;
        
        for (int i = rects.size() - 1; i >= 0; i--) {
            GenericRectangle rect = rects.get(i);
            double xLow, xHigh, yLow, yHigh;
            xLow = getXPixel(Math.min(rect.getStartXValue(), rect.getStopXValue()));
            xHigh = getXPixel(Math.max(rect.getStartXValue(), rect.getStopXValue()));
            yHigh = getYPixel(Math.min(rect.getStartYValue(), rect.getStopYValue()));
            yLow = getYPixel(Math.max(rect.getStartYValue(), rect.getStopYValue()));
            if (isHistogram()) {
                if (xLow <= tempX && xHigh >= tempX) {
                    toRemove = rect;
                    break;
                }
            } else {
                if (xLow <= tempX && xHigh >= tempX && yLow <= tempY && yHigh >= tempY) {
                    toRemove = rect;
                    break;
                }
            }
        }
        if (toRemove == null) {
            for (GenericRectangle rect : selectedRects) {
                double xLow, xHigh, yLow, yHigh;
                xLow = getXPixel(Math.min(rect.getStartXValue(), rect.getStopXValue()));
                xHigh = getXPixel(Math.max(rect.getStartXValue(), rect.getStopXValue()));
                yHigh = getYPixel(Math.min(rect.getStartYValue(), rect.getStopYValue()));
                yLow = getYPixel(Math.max(rect.getStartYValue(), rect.getStopYValue()));
                if (isHistogram()) {
                    if (xLow <= tempX && xHigh >= tempX) {
                        toRemove = rect;
                        break;
                    }
                } else {
                    if (xLow <= tempX && xHigh >= tempX && yLow <= tempY && yHigh >= tempY) {
                        toRemove = rect;
                        break;
                    }
                }
            }
            
        }
        if (toRemove != null) {
            rects.remove(toRemove);
            selectedRects.remove(toRemove);
        }
        // TODO remove gate in gating data struct
        // TODO confirm gate removal
        super.mouseClicked(e);
        paintAgain();
	}
	
	private void rightMouseClickedPoly(MouseEvent e) {
		int tempX = e.getX();
		int tempY = e.getY();
		GenericPath toRemove = null;
		
		if (tempPoly.size() > 0) {
		    tempPoly.clear();
		    highlightPoly = null;
		    return;
		}
		
		double tempValX = getXValueFromXPixel(tempX);
		double tempValY = getYValueFromYPixel(tempY);
	    for (int i = polys.size() - 1; i >= 0; i--) {
	    	if (polys.get(i).myPath.contains(tempValX, tempValY)) {
	    		toRemove = polys.get(i);
	    		break;
	    	}
	    }
		if (toRemove != null) {
			polys.remove(toRemove);
			selectedPolys.remove(toRemove);
		}

        // TODO remove gate in gating data struct
        // TODO confirm gate removal
		
		super.mouseClicked(e);
		paintAgain();
	}
	
    boolean drag = false;
    
    private void setForceGatesChanged() {
        this.forceGatesChanged = true;
    }

    @Override
    public void highlightPoints() {
        return;
    }

    @Override
    public void assignAxisLabels() {
        setXAxisLabel(fcp.getXDataName());
        setYAxisLabel(fcp.getYDataName());
    }
}
