package one.ben.fcs;

import java.awt.Color;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Hashtable;

import stats.Histogram;
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
	                                            new Color(182, 182, 182), // light grey
	                                            new Color(220, 220, 220), // very light grey
                                                new Color(33, 31, 53), // dark dark
			   									new Color(201, 30, 10), // deep red
			   									new Color(94, 88, 214), // light purple
			   									new Color(189, 243, 61), // light green
			   									new Color(217, 109, 194), // pink
			   									new Color(33, 87, 0), // dark green
			   									new Color(23, 58, 172), // dark blue
			   									new Color(140, 20, 180), // deep purple
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
	
	protected FCSPlot fcp;
	
	public FCSPanel(FCSPlot fcsPlot) {
		super();
		setDoubleBuffered(false);
		
		this.fcp = fcsPlot;
		this.setAxisFontSize(24);
		this.setSymmetricAxes(false);
		setZoomable(true, true);

		setColorScheme(DEFAULT_COLORS);

		setNullMessage("Select two variables to plot");
	}

	String xCol = null, yCol = null;
	float xMed, yMed, xSD, ySD;
	
	public void generatePoints() {
	    ArrayList<String[]> currentData;
		CountVector uniqueValueCounts;
		byte type;
//		String[] line;
		float xAxisValue, yAxisValue;
		byte size = 5;
		byte color = 0;
		
		if (fcp.dataLoader == null && !fcp.isLoading) {
		    setNullMessage("Please load an FCS file..");
		    points = new PlotPoint[0];
		    return;
		    
		}
		if (!fcp.isCurrentDataDisplayable()) {
		    setNullMessage("Please wait, data is loading...");
		    points = new PlotPoint[0];
		    return;
		}
		
		boolean generatePoints = true;
		
		int count = -1;
		if (xCol != null && yCol != null && fcp.getXDataName().equals(xCol) && fcp.getYDataName().equals(yCol) && (count = fcp.getDataCount()) > 0) {
		    generatePoints = false; // No need to re-create points[] if nothing has changed
		}
		if (generatePoints) {
    		if (count == -1) {
    		    count = fcp.getDataCount();
    		}
    		if (count <= 0) {
    		    generatePoints = false;
    		}
		} else if (count >= points.length) {
		    generatePoints = true;
		}
		
		PLOT_TYPE plotType = fcp.getPlotType();
		
		double[] xData = fcp.getAxisData(false, true);
		double[] yData = fcp.getAxisData(false, false);
		xCol = fcp.getXDataName();
		yCol = fcp.getYDataName();
		
		zoomable = true;
		rectangles = new GenericRectangle[0];
		setForcePlotXMin(0);// TODO fix this to be more intelligent
		setForcePlotYMin(0);//
		
		ArrayList<GenericLine> lineList = new ArrayList<GenericLine>();
		boolean mX = fcp.showMedian(false), mY = fcp.showMedian(true), sdX = fcp.showSD(false), sdY = fcp.showSD(true); 
		if (mX || sdX) {
		    double med = Array.median(xData);
		    double min = Math.min(Math.min(0, plotXmin) - med, Array.min(xData));
		    double max = Math.max(this.plotXmax + med, Array.max(xData));
		    if (mX) {
		        lineList.add(new GenericLine((float)med, (float)min, (float)med, (float)max, (byte)1, (byte) 0, (byte)1, 0, false));  
		    }
		    if (sdX) {
    		    double sd = Array.stdev(xData);
    		    lineList.add(new GenericLine((float)(med - sd), (float)min, (float)(med - sd), (float)max, (byte)1, (byte) 1, (byte)1, 0, false));  
    		    lineList.add(new GenericLine((float)(med + sd), (float)min, (float)(med + sd), (float)max, (byte)1, (byte) 1, (byte)1, 0, false));  
		    }
		}
		if (mY || sdY) {
		    double med = Array.median(yData);
		    double min = Math.min(Math.min(0, plotYmin) - med, Array.min(yData));
		    double max = Math.max(this.plotYmax + med, Array.max(yData));
		    if (mY) {
		        lineList.add(new GenericLine((float)min, (float)med, (float)max, (float)med, (byte)1, (byte) 0, (byte)1, 0, false));  
		    }
		    if (sdY) {
		        double sd = Array.stdev(yData);
		        lineList.add(new GenericLine((float)min, (float)(med - sd), (float)max, (float)(med - sd), (byte)1, (byte) 1, (byte)1, 0, false));  
		        lineList.add(new GenericLine((float)min, (float)(med + sd), (float)max, (float)(med + sd), (byte)1, (byte) 1, (byte)1, 0, false));  
		    }
		    
		}
        lines = lineList.toArray(new GenericLine[lineList.size()]);
        lineList = null;
		
        if (generatePoints) {
    		points = new PlotPoint[count];
    		for (int i = 0; i < points.length; i++) {
    			xAxisValue = (float) xData[i];
    			yAxisValue = (float) yData[i];
    			if (Float.isNaN(xAxisValue) || Float.isNaN(yAxisValue)) {
    				type = PlotPoint.NOT_A_NUMBER;
    			} else {
    				type = PlotPoint.FILLED_CIRCLE;
    			}
    			points[i] = new PlotPoint(i + "", type, xAxisValue, yAxisValue, size, color, (byte)0);
    		}
        }

	}

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
