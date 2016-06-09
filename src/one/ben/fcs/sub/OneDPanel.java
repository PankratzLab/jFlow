package one.ben.fcs.sub;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.util.ArrayList;

import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;

import one.ben.fcs.AbstractPanel2;
import one.ben.fcs.AbstractPanel2.AXIS_SCALE;
import cnv.plots.GenericLine;
import cnv.plots.GenericPath;
import cnv.plots.GenericRectangle;
import cnv.plots.PlotPoint;
import common.Array;
import common.ext;

public class OneDPanel extends AbstractPanel2  {
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
    
    public static enum PLOT_TYPE {
        BOX_PLOT,
        DOT_LINE_PLOT,
    }
    
    private PLOT_TYPE currentPlot;
    
    public void setPlotType(PLOT_TYPE type) {
        this.currentPlot = type;
        switch(type) {
            case BOX_PLOT:
                this.setAxisFontSize(12);

                setForcePlotXMax(20);
                setDisplayXAxis(false);
                break;
            case DOT_LINE_PLOT:
                
                setAxisFontSize(24);
                
                setForcePlotXMax(Float.NaN);
                setDisplayXAxis(true);

                break;
        }
        
        // for all:
        setDoubleBuffered(false);
        setSymmetricAxes(false);
        setYAxis(AXIS_SCALE.LIN);
        setXAxis(AXIS_SCALE.LIN);
        setZoomable(false, false);
        setColorScheme(DEFAULT_COLORS);
        createLookup(true);
        
    }
    
	public OneDPanel() {
		super();
		setPlotType(PLOT_TYPE.BOX_PLOT);
	}

	double[][] data;// = {11.8, 0.93, 1.76, 14, 16.5, 17.1, 32.5, 33.4, 16.8, 21.5, 13.1, 22.2, 22.2, 16, 16.2};
	String[][] dataLabels; 
	String plotLabel;
	
	public void setData(String dataName, String[] files, double[] data) {
	    this.plotLabel = dataName;
	    this.dataLabels = new String[][]{files};
	    this.data = new double[][]{data};
	}
	
	public void setData(String dataName, String[][] files, double[][] data) {
	    this.plotLabel = dataName;
	    this.dataLabels = files;
	    this.data = data;
	}
	
    public void generatePointsRectanglesAndLines() {
		if (data == null) {
		    setNullMessage("Error, no data set..");
		    points = new PlotPoint[0];
            rectangles = new GenericRectangle[0];
            lines = new GenericLine[0];
            polygons = new GenericPath[0];
		    return;
		}
		
		switch(currentPlot) {
		    case BOX_PLOT:
		        generateBoxPlot();
		        break;
		    case DOT_LINE_PLOT:
		        generateDotLinePlot();
		        break;
		}
		
    }
    
    private void generateDotLinePlot() {
        byte type;
        float xAxisValue, yAxisValue;
        byte size = POINT_SIZE;
        
        int numPoints = 0;
        for (double[] dataArr : data) {
            numPoints += dataArr.length;
        }
        points = new PlotPoint[numPoints];
        lines = new GenericLine[numPoints + 3];
        
        int ind = 0;
        int lInd = 0;
        byte color;
        for (int d = 0; d < data.length; d++) {
            for (int i = 0; i < data[d].length; i++) {
                xAxisValue = (float) ind;
                yAxisValue = (float) data[d][i];
                if (Float.isNaN(xAxisValue) || Float.isNaN(yAxisValue)) {
                    type = PlotPoint.NOT_A_NUMBER;
                } else {
                    type = PlotPoint.FILLED_CIRCLE;
                }
                
                color = (byte) 0; // TODO apply gating for colors
                points[ind] = new PlotPoint(dataLabels[d][i], type, xAxisValue, yAxisValue, size, color, (byte)0);
                if (i < data[d].length - 1) {
                    lines[lInd++] = new GenericLine(xAxisValue, yAxisValue, (float) ind+1, (float) data[d][i+1], (byte)1, (byte)d, (byte)0);
                }
                ind++;
            }
        }
        
        
        float mean = (float)Array.mean(data[0]);
        float sd = (float)Array.stdev(data[0], true);
        lines[lInd++] = new GenericLine(-1, mean, numPoints + 1, mean, (byte)1, (byte)0, (byte)99);
        lines[lInd++] = new GenericLine(-1, mean - sd, numPoints + 1, mean - sd, (byte)1, (byte)2, (byte)99);
        lines[lInd++] = new GenericLine(-1, mean + sd, numPoints + 1, mean + sd, (byte)1, (byte)2, (byte)99);
        lines[lInd++] = new GenericLine(-1, mean - 2*sd, numPoints + 1, mean - 2*sd, (byte)1, (byte)2, (byte)99);
        lines[lInd++] = new GenericLine(-1, mean + 2*sd, numPoints + 1, mean + 2*sd, (byte)1, (byte)2, (byte)99);
        
        double dataMin = Array.min(data[0]), dataMax = Array.max(data[0]);
        for (int i = 1; i < data.length; i++) {
            dataMin = Math.min(dataMin, Array.min(data[i]));
            dataMax = Math.max(dataMax, Array.max(data[i]));
        }
        setForcePlotYMin((float)dataMin);
        setForcePlotYMax((float)dataMax);
        
        setYAxis(AXIS_SCALE.LIN);
        setXAxis(AXIS_SCALE.LIN);
    }
    
    private void generateBoxPlot() {
        //  points for any data above/below wiskLow/wiskHigh
		ArrayList<GenericLine> lns = new ArrayList<GenericLine>();
        ArrayList<PlotPoint> pts = new ArrayList<PlotPoint>();
        
        for (int i = 0; i < data.length; i++) {
            float xLow = 20 * i + 2;
            float xHigh = xLow + 18; 
            float xMed = xLow + (xHigh-xLow)/2;
            
            double med = Array.median(data[i]);
            double qr25 = Array.quantExclusive(data[i], 0.25);
            double qr75 = Array.quantExclusive(data[i], 0.75);
            double iqr = Array.iqrExclusive(data[i]);
            double wiskLow = qr25 - 1.5 * iqr; 
            double wiskHigh = qr75 + 1.5 * iqr;
            double min = Math.min(wiskLow, Array.min(data[i]));
            double max = Math.max(wiskHigh, Array.max(data[i]));
            setForcePlotYMin((float) (min));
            setForcePlotYMax((float) (max));
            setPlotYMin((float) (min));
            setPlotYMax((float) (max));
            
            //  line @ med
            lns.add(new GenericLine(xLow, (float)med, xHigh, (float)med, (byte)4, (byte)0, (byte)0));
            //  line @ qr25
            lns.add(new GenericLine(xLow, (float)qr25, xHigh, (float)qr25, (byte)2, (byte)0, (byte)0));
            //  line @ qr75
            lns.add(new GenericLine(xLow, (float)qr75, xHigh, (float)qr75, (byte)2, (byte)0, (byte)0));
            //  small line at wiskLow
            lns.add(new GenericLine(xMed - (xMed - xLow) / 2, (float)wiskLow, xMed + (xMed - xLow) / 2, (float)wiskLow, (byte)2, (byte)0, (byte)0));
            //  small line at wiskHigh
            lns.add(new GenericLine(xMed - (xMed - xLow) / 2, (float)wiskHigh, xMed + (xMed - xLow) / 2, (float)wiskHigh, (byte)2, (byte)0, (byte)0));
            //  line from qr25 -> wiskLow
            lns.add(new GenericLine(xMed, (float)qr25, xMed, (float)wiskLow, (byte)1, (byte)0, (byte)0));
            //  line from qr75 -> wiskHigh
            lns.add(new GenericLine(xMed, (float)qr75, xMed, (float)wiskHigh, (byte)1, (byte)0, (byte)0));
            //  two lines vert, from qr25 to qr75
            lns.add(new GenericLine(xLow, (float)qr25, xLow, (float)qr75, (byte)1, (byte)0, (byte)0));
            lns.add(new GenericLine(xHigh, (float)qr25, xHigh, (float)qr75, (byte)1, (byte)0, (byte)0));
            
            for (int j = 0; j < data[i].length; j++) {
                if (data[i][j] < wiskLow || data[i][j] > wiskHigh) {
                    pts.add(new PlotPoint(dataLabels[i][j], PlotPoint.FILLED_CIRCLE, xMed, (float)data[i][j], (byte)POINT_SIZE, (byte)0, (byte)0));
                }
            }
        }
        
		
		lines = lns.toArray(new GenericLine[lns.size()]);
		points = pts.toArray(new PlotPoint[pts.size()]);
		
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
            points[i].setSize((byte)(defaultSize * (points[i].isHighlighted() ? 1.5 : 1)));
        }
    }
    
    @Override
    public void assignAxisLabels() {
//        String[] pts = plotLabel.split("\\|");
//        setXAxisLabel("");//pts[0].trim().replaceAll("/", " /\n");
//        setYAxisLabel(pts[1].trim());

//        setXAxisLabel("File by Date");
//        setYAxisLabel("Mean - " + plotLabel);
    }
    
}
