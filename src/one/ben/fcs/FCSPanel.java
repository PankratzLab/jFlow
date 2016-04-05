package one.ben.fcs;

import java.awt.Color;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Hashtable;

import stats.Histogram;
//import cnv.filesys.MarkerLookup;
import cnv.plots.GenericRectangle;
import cnv.plots.PlotPoint;

import common.CountVector;
import common.IntVector;
import common.ext;

public class FCSPanel extends AbstractPanel2 implements MouseListener, MouseMotionListener {
	public static final long serialVersionUID = 3L;
	public static final int LOOKUP_RESOLUTION = 20;
	public static final Color[] DEFAULT_COLORS = {new Color(33, 31, 53), // dark dark
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
		
		this.fcp = fcsPlot;
		this.setAxisFontSize(24);
		setZoomable(true, true);

		setColorScheme(DEFAULT_COLORS);

		setNullMessage("Select two variables to plot");
	}

	String xCol, yCol;
	
	public void generatePoints() {
	    ArrayList<String[]> currentData;
		CountVector uniqueValueCounts;
		byte type;
//		String[] line;
		float xAxisValue, yAxisValue;
		byte size = 5;
		byte color = 0;
		
		if (fcp.isLoading) {
		    setNullMessage("Please wait, data is loading...");
		    points = new PlotPoint[0];
//		    repaint();
		    return;
		}
		
		if (xCol != null && yCol != null && fcp.getXDataName().equals(xCol) && fcp.getYDataName().equals(yCol) && points != null) {
		    return; // No need to re-create points[] if nothing has changed
		}
		xCol = fcp.getXDataName();
		yCol = fcp.getYDataName();
		
		double[] xData = fcp.getAxisData(false, true);
		double[] yData = fcp.getAxisData(false, false);
		int count = fcp.getDataCount();
//		currentData = new ArrayList<String[]>();// tdp.getDataSelected();
//		uniqueValueCounts = new CountVector();
		
		zoomable = true;
		rectangles = new GenericRectangle[0];
        forcePlotXmax = Float.NaN;
        forcePlotXmin = Float.NaN;
		
		points = new PlotPoint[count];   
		for (int i = 0; i < points.length; i++) {
//			line = currentData.get(i);
			xAxisValue = (float) xData[i];
			yAxisValue = (float) yData[i];
			if (Float.isNaN(xAxisValue) || Float.isNaN(yAxisValue)) {
				type = PlotPoint.NOT_A_NUMBER;
//				uniqueValueCounts.add("0");
			} else {
				type = PlotPoint.FILLED_CIRCLE;
//				uniqueValueCounts.add(line[3]);
			}
			points[i] = new PlotPoint(i + "", type, xAxisValue, yAxisValue, size, color, (byte)0);
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
