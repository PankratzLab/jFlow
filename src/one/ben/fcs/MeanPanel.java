package one.ben.fcs;

import java.awt.Color;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.awt.image.BufferedImage;
import java.util.ArrayList;

import javax.swing.SwingUtilities;

import one.ben.fcs.gating.Gate.PolygonGate;
import one.ben.fcs.gating.Gate.RectangleGate;
import one.ben.fcs.gating.GateDimension;
import stats.Histogram;
import cnv.plots.GenericLine;
import cnv.plots.GenericPath;
//import cnv.filesys.MarkerLookup;
import cnv.plots.GenericRectangle;
import cnv.plots.PlotPoint;
import common.Array;

public class MeanPanel extends AbstractPanel2  {
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
	
	
	public MeanPanel() {
		super();
		setDoubleBuffered(false);
		
		this.setAxisFontSize(24);
		this.setSymmetricAxes(false);
//		setZoomable(true, true);
		setZoomable(false, true);

		setColorScheme(DEFAULT_COLORS);

//		setNullMessage("Select two variables to plot");
	}

	float[] xDataBase; // Files - in date order
	float[] xDataComp; // Files - in date order
	float[] yDataBase; // Means
	float[] yDataComp; // Means
	
	public void setData(float[] xDataBase, float[] xDataComp, float[] yDataBase, float[] yDataComp) {
	    this.xDataBase = xDataBase;
	    this.xDataComp = xDataComp;
	    this.yDataBase = yDataBase;
	    this.yDataComp = yDataComp;
	}

    public void generatePointsRectanglesAndLines() {
		byte type;
//		String[] line;
		float xAxisValue, yAxisValue;
		byte size = POINT_SIZE;
		
		if (xDataBase == null || yDataBase == null) {
		    setNullMessage("Error, no data set..");
		    points = new PlotPoint[0];
            rectangles = new GenericRectangle[0];
            lines = new GenericLine[0];
            polygons = new GenericPath[0];
		    return;
		    
		}
		
		points = new PlotPoint[xDataBase.length + xDataComp.length];
		lines = new GenericLine[xDataBase.length + xDataComp.length + 3];
		byte color;
		for (int i = 0; i < xDataBase.length; i++) {
			xAxisValue = (float) xDataBase[i];
			yAxisValue = (float) yDataBase[i];
			if (Float.isNaN(xAxisValue) || Float.isNaN(yAxisValue)) {
				type = PlotPoint.NOT_A_NUMBER;
			} else {
				type = PlotPoint.FILLED_CIRCLE;
			}
			
			color = (byte) 0; // TODO apply gating for colors
			points[i] = new PlotPoint(i + "", type, xAxisValue, yAxisValue, size, color, (byte)0);
			if (i < xDataBase.length - 1) {
			    lines[i] = new GenericLine(xAxisValue, yAxisValue, (float) xDataBase[i+1], (float) yDataBase[i+1], (byte)1, (byte)0, (byte)0);
			}
		}
		for (int i = 0; i < xDataComp.length; i++) {
		    xAxisValue = (float) xDataComp[i];
		    yAxisValue = (float) yDataComp[i];
		    if (Float.isNaN(xAxisValue) || Float.isNaN(yAxisValue)) {
		        type = PlotPoint.NOT_A_NUMBER;
		    } else {
		        type = PlotPoint.FILLED_CIRCLE;
		    }
		    color = (byte) 1; // TODO apply gating for colors
		    points[i] = new PlotPoint(i + "", type, xAxisValue, yAxisValue, size, color, (byte)0);
		    if (i < xDataBase.length - 1) {
		        lines[i] = new GenericLine(xAxisValue, yAxisValue, (float) xDataComp[i+1], (float) yDataComp[i+1], (byte)1, (byte)0, (byte)0);
		    }
		}
		float mean = Array.mean(yDataBase);
		float sd = Array.stdev(yDataBase, true);
		lines[xDataBase.length + xDataComp.length] = new GenericLine(Float.NEGATIVE_INFINITY, mean, Float.POSITIVE_INFINITY, mean, (byte)1, (byte)0, (byte)99);
		lines[xDataBase.length + xDataComp.length +1] = new GenericLine(Float.NEGATIVE_INFINITY, mean - sd, Float.POSITIVE_INFINITY, mean - sd, (byte)1, (byte)0, (byte)99);
		lines[xDataBase.length + xDataComp.length +2] = new GenericLine(Float.NEGATIVE_INFINITY, mean + sd, Float.POSITIVE_INFINITY, mean + sd, (byte)1, (byte)0, (byte)99);
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
        return;
    }
    
}
