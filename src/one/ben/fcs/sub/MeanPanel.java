package one.ben.fcs.sub;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;

import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;

import one.ben.fcs.AbstractPanel2;
import cnv.plots.GenericLine;
import cnv.plots.GenericPath;
import cnv.plots.GenericRectangle;
import cnv.plots.PlotPoint;
import common.Array;
import common.ext;

public class MeanPanel extends AbstractPanel2  {
	public static final long serialVersionUID = 3L;
	public static final int LOOKUP_RESOLUTION = 20;
	public static final Color[] DEFAULT_COLORS = {
                                                new Color(33, 31, 53), // dark dark
			   									new Color(201, 30, 10), // deep red
			   									new Color(182, 182, 182), // light grey
			   									new Color(94, 88, 214), // light purple
			   									new Color(189, 243, 61), // light green
			   									new Color(217, 109, 194), // pink
			   									new Color(33, 87, 0), // dark green
			   									new Color(23, 58, 172), // dark blue
			   									new Color(140, 20, 180), // deep purple
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
		setZoomable(false, true);

		setColorScheme(DEFAULT_COLORS);
		createLookup(true);
	}

	float[] xDataBase; // Files - in date order
	float[] xDataComp; // Files - in date order
	float[] yDataBase; // Means
	float[] yDataComp; // Means
	String[] compLbls;
	String[] baseLbls;
	String col;
	
	public void setData(String col, String[] baseLbls, float[] xDataBase, float[] xDataComp, String[] compLbls, float[] yDataBase, float[] yDataComp) {
	    this.col = col;
	    this.compLbls = compLbls;
	    this.baseLbls = baseLbls;
	    this.xDataBase = xDataBase;
	    this.xDataComp = xDataComp;
	    this.yDataBase = yDataBase;
	    this.yDataComp = yDataComp;
	}

    public void generatePointsRectanglesAndLines() {
		byte type;
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
			points[i] = new PlotPoint(baseLbls[i], type, xAxisValue, yAxisValue, size, color, (byte)0);
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
		    points[i + xDataBase.length] = new PlotPoint(compLbls[i], type, xAxisValue, yAxisValue, size, color, (byte)0);
		    if (i < xDataComp.length - 1) {
		        lines[i - 1 + xDataBase.length] = new GenericLine(xAxisValue, yAxisValue, (float) xDataComp[i+1], (float) yDataComp[i+1], (byte)1, (byte)1, (byte)0);
		    }
		}
		float mean = Array.mean(yDataBase);
		float sd = Array.stdev(yDataBase, true);
		lines[xDataBase.length + xDataComp.length - 2] = new GenericLine(-1, mean, xDataBase.length + xDataComp.length + 1, mean, (byte)1, (byte)0, (byte)99);
		lines[xDataBase.length + xDataComp.length - 1] = new GenericLine(-1, mean - sd, xDataBase.length + xDataComp.length + 1, mean - sd, (byte)1, (byte)2, (byte)99);
		lines[xDataBase.length + xDataComp.length] = new GenericLine(-1, mean + sd, xDataBase.length + xDataComp.length + 1, mean + sd, (byte)1, (byte)2, (byte)99);
		lines[xDataBase.length + xDataComp.length + 1] = new GenericLine(-1, mean - 2*sd, xDataBase.length + xDataComp.length + 1, mean - 2*sd, (byte)1, (byte)2, (byte)99);
		lines[xDataBase.length + xDataComp.length + 2] = new GenericLine(-1, mean + 2*sd, xDataBase.length + xDataComp.length + 1, mean + 2*sd, (byte)1, (byte)2, (byte)99);
		
		setForcePlotYMin(Math.min(Array.min(yDataBase), Array.min(yDataComp)));
		setForcePlotYMax(Math.max(Array.max(yDataBase), Array.max(yDataComp)));
		
		setYAxis(AXIS_SCALE.LIN);
		setXAxis(AXIS_SCALE.LIN);
	}
	

	public BufferedImage getImage() {
        return image;
    }
    
    public void mouseClicked(MouseEvent e) {
        JPopupMenu menu;
        byte maxNumPoints;

        if (prox != null && prox.size() > 0) {
            menu = new JPopupMenu();
            maxNumPoints = (byte) Math.min(20, prox.size());
            for (int i = 0; i < maxNumPoints; i++) {
                final int ind = prox.elementAt(i);
                JMenuItem jmi = new JMenuItem(points[ind].getId() + " -- " + points[ind].getRawY());
                jmi.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        ext.setClipboard(points[ind].getId());                        
                    }
                });
                menu.add(jmi);
            }
            menu.show(this, e.getX(), e.getY());
        }
    }
	
    @Override
    public void highlightPoints() {
        byte defaultSize;
        
        defaultSize = POINT_SIZE;
        for (int i = 0; i < points.length; i++) {
            if (points[i].isHighlighted()) {
                points[i].setSize((byte)(defaultSize * 1.5));
            } else {
                points[i].setSize((byte)(defaultSize));
            }
            
        }
    }


    @Override
    public void assignAxisLabels() {
        xAxisLabel = "File by Date";
        yAxisLabel = "Mean - " + col;
        return;
    }
    
}
