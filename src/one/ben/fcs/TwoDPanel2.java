package one.ben.fcs;

import java.awt.Color;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Hashtable;

import stats.Histogram;
//import cnv.filesys.MarkerLookup;
import cnv.plots.AbstractPanel;
import cnv.plots.GenericRectangle;
import cnv.plots.PlotPoint;

import common.CountVector;
import common.IntVector;
import common.ext;

public class TwoDPanel2 extends AbstractPanel implements MouseListener, MouseMotionListener {
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
	private static final double DEFAULT_HALF_BIN_SIZE = 0.001;
	
	
	protected TwoDPlot2 tdp;
	IntVector indicesOfNearbySamples;
	private boolean updateQcPanel;
	private boolean swapAxes;
	private String[][] linkerData;
	private int dataHash = -1;
	private stats.Histogram currentHistogram;
	private boolean histogramOverride = false;
    private boolean overrideAxisLabels = false;
	
	public TwoDPanel2(TwoDPlot2 twoDPlot) {
		super();
		
		this.tdp = twoDPlot;
		locLookup = new Hashtable<String,IntVector>();//??? zx
		linkerData = new String[0][0];
		this.setAxisFontSize(24);
		setZoomable(true, true);


		setColorScheme(DEFAULT_COLORS);

		setNullMessage("Select two variables to plot");
	}
	
	public void refreshOtherComponents() {
		tdp.refreshOtherButtons();
	}

	public void overrideAxisLabels(String xAxis, String yAxis) {
	    this.overrideAxisLabels  = xAxis != null && yAxis != null;
	    if (this.overrideAxisLabels) {
	        xAxisLabel = xAxis;
	        yAxisLabel = yAxis;
	    }
	}
	
	public void assignAxisLabels() {
		displayXaxis = displayYaxis = true;
		if (!overrideAxisLabels) {
//		    xAxisLabel = "Bins";
//		    yAxisLabel = "";
//		} else {
		    String[] data = tdp.getNamesSelected();
		    xAxisLabel = tdp.isHistogram() ? data[0] : swapAxes ? data[1] : data[0];
		    yAxisLabel = tdp.isHistogram() ?      "" : swapAxes ? data[0] : data[1];
		}
	}
	
	public boolean invertX() {
		return false;
	}

	public boolean invertY() {
		return false;
	}
	
//	public void toggleMasking() {
//		
//	}	
	
	public void highlightPoints() {
		byte defaultSize;
		
		defaultSize = tdp.getPointSize();
		for (int i = 0; i < points.length; i++) {
			if (points[i].isHighlighted()) {
				points[i].setSize((byte)(defaultSize * 1.5));
			} else {
				points[i].setSize((byte)(defaultSize));
			}
			
		}
	}

	public void generatePoints() {
	    ArrayList<String[]> currentData;
		CountVector uniqueValueCounts;
		boolean includeColorKeyValue;
		byte type;
		String[] line;
		float xAxisValue, yAxisValue;
		byte index;

		includeColorKeyValue = true;
		currentData = tdp.getDataSelected(includeColorKeyValue);
		uniqueValueCounts = new CountVector();
		
		if (tdp.isHistogram()) {
			zoomable = false;
//			tdp.setPointSize((byte) 0);
			points = new PlotPoint[0];//currentData.size()];
			for (int i = 0; i < points.length; i++) {
				points[i] = new PlotPoint("" + Float.parseFloat(currentData.get(i)[1]), PlotPoint.FILLED_SQUARE, Float.parseFloat(currentData.get(i)[1]), Float.parseFloat(currentData.get(i)[2]), (byte) 0, (byte) 0, (byte) 0);
			}
			generateRectangles();
			return;
		} else {
			zoomable = true;
//			tdp.setPointSize(tdp.getPointSize());
			rectangles = new GenericRectangle[0];
            forcePlotXmax = Float.NaN;
            forcePlotXmin = Float.NaN;
		}
		
		points = new PlotPoint[currentData.size()];
		index = (byte) (includeColorKeyValue? 4 : 3);
		if (currentData.size()>0) {
			linkerData = new String[currentData.size()][currentData.get(0).length - index];
		}
		for (int i = 0; i < points.length; i++) {
			line = currentData.get(i);
			boolean missing = false;
			for (String miss : TwoDPlot2.MISSING_VALUES) {
				if (miss.equals(line[1]) || miss.equals(line[2])) {
					missing = true;
					break;
				}	
			}
			if (missing) {
				xAxisValue = Float.NaN;
				yAxisValue = Float.NaN;
				type = PlotPoint.MISSING;
				uniqueValueCounts.add("0");
			} else {
				xAxisValue = Float.parseFloat(line[1]);
				yAxisValue = Float.parseFloat(line[2]);
				if (Float.isNaN(xAxisValue) || Float.isNaN(xAxisValue)) {
					type = PlotPoint.NOT_A_NUMBER;
					uniqueValueCounts.add("0");
	//			} else if (alleleCounts[i]==-1) {
	//				type = PlotPoint.MISSING;
	//				uniqueValueCounts.add("0");
				} else {
					type = PlotPoint.FILLED_CIRCLE;
					uniqueValueCounts.add(line[3]);
				}
			}
			if (swapAxes) {
				points[i] = new PlotPoint(line[0], type, yAxisValue, xAxisValue, (byte)5, Byte.parseByte(line[3]), (byte)0);
			} else {
				points[i] = new PlotPoint(line[0], type, xAxisValue, yAxisValue, (byte)5, Byte.parseByte(line[3]), (byte)0);
			}

			for (int j = 0; j < linkerData[i].length; j ++) {
				linkerData[i][j] = line[j + index];
			}
		}

		tdp.updateColorKey(uniqueValueCounts.convertToHash());
	}

	private void generateRectangles() {
	    stats.Histogram hist;
	    if ((!isHistogramOverride() && dataHash != tdp.getSelectedDataHash()) || currentHistogram == null) {
    	    ArrayList<String[]> currentData;
    		boolean includeColorKeyValue;
    		int index;
    		
    		includeColorKeyValue = true;
    		currentData = tdp.getDataSelected(includeColorKeyValue);
    		index = includeColorKeyValue? 4 : 3;
    		if (currentData.size() > 0) {
    			linkerData = new String[currentData.size()][currentData.get(0).length - index];
    		} else {
    		    rectangles = new GenericRectangle[0];
    		    return;
    		}
    		
    		double[] dataArray = new double[currentData.size()];
    		for (int i = 0; i < dataArray.length; i++) {
    		    dataArray[i] = ext.isMissingValue(currentData.get(i)[1]) ? Double.NaN : Double.parseDouble(currentData.get(i)[1]);
    		}
    		hist = new stats.Histogram(dataArray);
    		setHistogram(hist);
	    } else {
            hist = getHistogram();
	    }

		double min = hist.getMin();
		double max = hist.getMax();
		double minDiff = hist.determineStep();
		int sig = hist.getSigfigs();

		float binHalf = (float) (hist.getBins().length > 1 ? ext.roundToSignificantFigures(minDiff / 2.0, sig + 1) : DEFAULT_HALF_BIN_SIZE);
		
		forcePlotXmax = (float) (max + minDiff);
		forcePlotXmin = (float) (min - minDiff);
		
		rectangles = new GenericRectangle[hist.getBins().length];
		for (int i = 0; i < rectangles.length; i++) {
		
			float startX, startY, stopX, stopY;
			if (swapAxes) {
				startX = 0f;
				startY = sig > 0 ? (float) ext.roundToSignificantFigures((float)(hist.getBins()[i] - binHalf), sig) : (float)(hist.getBins()[i] - binHalf);
				stopX = hist.getCounts()[i];
				stopY = sig > 0 ? (float) ext.roundToSignificantFigures((float)(hist.getBins()[i] + binHalf), sig) : (float)(hist.getBins()[i] + binHalf);
			} else {
				startX = sig > 0 ? (float) ext.roundToSignificantFigures((float)(hist.getBins()[i] - binHalf), sig) : (float)(hist.getBins()[i] - binHalf);
				startY = 0f; 
				stopX = sig > 0 ? (float) ext.roundToSignificantFigures((float)(hist.getBins()[i] + binHalf), sig) : (float)(hist.getBins()[i] + binHalf);
				stopY = hist.getCounts()[i];
			}
			rectangles[i] = new GenericRectangle(startX, startY, stopX, stopY, (byte) 5, true, false, (byte) 0, (byte) 2, (byte) 0);

		}
		
	}
	
	public void setUpdateQcPanel(boolean updateQcPanel) {
		this.updateQcPanel = updateQcPanel;
	}

	public boolean getUpdateQcPanel() {
		return updateQcPanel;
	}
	
    public boolean isHistogramOverride() {
        return histogramOverride;
    }

    public void setHistogramOverride(boolean histogramOverride) {
        this.histogramOverride = histogramOverride;
    }
    
    private Histogram getHistogram() {
        return currentHistogram;
    }

    public void setHistogram(stats.Histogram hist) {
        this.currentHistogram = hist;
    }
    
    public void setSwapAxes(boolean swapAxes) {
		this.swapAxes = swapAxes;
	}

	public boolean isSwapAxes() {
		return swapAxes;
	}

    public BufferedImage getImage() {
        return image;
    }
}
