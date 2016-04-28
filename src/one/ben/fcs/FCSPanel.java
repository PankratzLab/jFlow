package one.ben.fcs;

import java.awt.Color;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Hashtable;

import javax.swing.SwingUtilities;

import one.ben.fcs.gating.Gate;
import one.ben.fcs.gating.Gate.*;
import one.ben.fcs.gating.GateDimension.RectangleGateDimension;
import stats.Histogram;
import cnv.filesys.ClusterFilter;
import cnv.plots.GenericLine;
//import cnv.filesys.MarkerLookup;
import cnv.plots.GenericRectangle;
import cnv.plots.PlotPoint;
import common.Array;
import common.CountVector;
import common.IntVector;
import common.ext;

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
	
	public void generatePointsRectanglesAndLines() {
		byte type;
//		String[] line;
		float xAxisValue, yAxisValue;
		byte size = POINT_SIZE;
		
		if (fcp.dataLoader == null && !fcp.isLoading) {
		    setNullMessage("Please load an FCS file..");
		    points = new PlotPoint[0];
            rectangles = new GenericRectangle[0];
		    return;
		    
		}
		if (!fcp.isCurrentDataDisplayable()) {
		    setNullMessage("Please wait, data is loading...");
		    points = new PlotPoint[0];
            rectangles = new GenericRectangle[0];
		    return;
		}
		
		boolean columnsChangedX = false;
		boolean columnsChangedY = false;
		boolean dataChanged = false;
		boolean optionsChanged = false;
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

        rectangles = rects.toArray(new GenericRectangle[rects.size()]);
		
		boolean skip = !columnsChangedX && !columnsChangedY && !dataChanged && !optionsChanged/* && !typeChanged /* don't need to regen if only type has changed, for now */;
		if (skip) return;
		
		xData = columnsChangedX || dataChanged || xData == null ? fcp.getAxisData(false, true) : xData;
		yData = columnsChangedY || dataChanged || yData == null ? fcp.getAxisData(false, false) : yData;
		
//		zoomable = true;
//		setForcePlotXMin(0);// TODO fix this to be more intelligent
//		setForcePlotYMin(0);//
		
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

        byte color = 0;
        if (columnsChangedX || columnsChangedY || dataChanged) {
    		points = new PlotPoint[dataCount];
    		for (int i = 0; i < points.length; i++) {
    			xAxisValue = (float) xData[i];
    			yAxisValue = (float) yData[i];
    			if (Float.isNaN(xAxisValue) || Float.isNaN(yAxisValue)) {
    				type = PlotPoint.NOT_A_NUMBER;
    			} else {
    				type = PlotPoint.FILLED_CIRCLE;
    			}
    			color = 0; // TODO apply gating for colors
    			points[i] = new PlotPoint(i + "", type, xAxisValue, yAxisValue, size, color, (byte)0);
    		}
            rects.clear();
            
            rectangles = rects.toArray(new GenericRectangle[rects.size()]);
        }
	}
	

	public void mousePressed(MouseEvent e) {
        int tempX = e.getX();
        int tempY = e.getY();
        int toRemove = -1;
        if (SwingUtilities.isLeftMouseButton(e) && !e.isControlDown()) {
            for (int i = 0; i < rects.size(); i++) {
                GenericRectangle rect = rects.get(i);
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
                    startX = getXPixel(tempStartX);
                    startY = getYPixel(tempStartY);
                    highlightRectangle = new GenericRectangle((float)tempStartX, 
                                                                (float)tempStartY,
                                                                (float)tempStopX, 
                                                                (float)tempStopY,
                                                                (byte)1, false, false, (byte)0, (byte)99);
                    toRemove = i;
                    break;
                }
            }
            if (toRemove >= 0) {
                rects.remove(toRemove);
            } else {
                startX = e.getX();
                startY = e.getY();
            }
            paintAgain();
        } else {
//        	double tempValX = getXValueFromXPixel(tempX);
//        	double tempValY = getYValueFromYPixel(tempY);
//            for (int i = 0; i < rects.size(); i++) {
//            	GenericRectangle rect = rects.get(i);
//            	if (rect.getStartXValue() <= tempValX && rect.getStopXValue() >= tempValX && rect.getStartYValue() <= tempValY && rect.getStopYValue() >= tempValY) {
//            		toRemove = i;
//            		break;
//            	}
//            }
//            if (toRemove != -1) {
//            	rects.remove(toRemove);
//            }
            super.mousePressed(e);
        }
    }

	public void mouseClicked(MouseEvent e) {
		int tempX = e.getX();
		int tempY = e.getY();
		int toRemove = -1;
		if (SwingUtilities.isLeftMouseButton(e) && !e.isControlDown()) {
			// 
		} else {
			double tempValX = getXValueFromXPixel(tempX);
			double tempValY = getYValueFromYPixel(tempY);
			for (int i = rects.size() - 1; i >= 0; i--) {
				GenericRectangle rect = rects.get(i);
				double xLow, xHigh, yLow, yHigh;
				xLow = Math.min(rect.getStartXValue(), rect.getStopXValue());
				xHigh = Math.max(rect.getStartXValue(), rect.getStopXValue());
				yLow = Math.min(rect.getStartYValue(), rect.getStopYValue());
				yHigh = Math.max(rect.getStartYValue(), rect.getStopYValue());
				if (xLow <= tempValX && xHigh >= tempValX && yLow <= tempValY && yHigh >= tempValY) {
					toRemove = i;
					break;
				}
			}
			if (toRemove != -1) {
				rects.remove(toRemove);
			}
			super.mouseClicked(e);
			paintAgain();
		}
	}
    

    public void mouseReleased(MouseEvent e) {
        int mouseEndX;
        int mouseEndY;
        mouseEndX = e.getX();
        mouseEndY = e.getY();
        highlightRectangle = null;
        if (SwingUtilities.isLeftMouseButton(e) && !e.isControlDown()) {
            rects.add(new GenericRectangle((float)getXValueFromXPixel(startX), 
                    (float)getYValueFromYPixel(startY), (float)getXValueFromXPixel(mouseEndX), 
                    (float)getYValueFromYPixel(mouseEndY), (byte)1, false, false, (byte)0, (byte)99));
            paintAgain();
        } else {
            super.mouseReleased(e);
        }
    }
    
    public void mouseDragged(MouseEvent e) {
        int mouseEndX;
        int mouseEndY;
        if (SwingUtilities.isLeftMouseButton(e) && !e.isControlDown()) {
            mouseEndX = e.getX();
            mouseEndY = e.getY();
            highlightRectangle = new GenericRectangle((float)getXValueFromXPixel(startX), 
                    (float)getYValueFromYPixel(startY), (float)getXValueFromXPixel(mouseEndX), 
                    (float)getYValueFromYPixel(mouseEndY), (byte)1, false, false, (byte)0, (byte)99);
            paintAgain();
        } else {
            super.mouseDragged(e);
        }
    }
    
//    public void mouseClicked(MouseEvent event) {}

    public BufferedImage getImage() {
        return image;
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
