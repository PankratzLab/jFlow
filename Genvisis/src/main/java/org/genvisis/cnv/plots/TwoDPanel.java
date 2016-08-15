package org.genvisis.cnv.plots;

import java.awt.Color;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Hashtable;

import javax.swing.JPopupMenu;

import org.genvisis.cnv.filesys.MarkerLookup;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.gui.LaunchAction;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.CountVector;
import org.genvisis.common.Files;
import org.genvisis.common.IntVector;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.stats.Histogram;

public class TwoDPanel extends AbstractPanel implements MouseListener, MouseMotionListener {
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
	
	
	protected TwoDPlot tdp;
	IntVector indicesOfNearbySamples;
	private boolean updateQcPanel;
	private boolean swapAxes;
	private MarkerLookup markerLookup;
	private SampleData sampleData;
	private Project proj;
	private String[][] linkerData;
	private int dataHash = -1;
	private org.genvisis.stats.Histogram currentHistogram;
	private boolean histogramOverride = false;
    private boolean overrideAxisLabels = false;
	
	public TwoDPanel(TwoDPlot twoDPlot) {
		super();
		
		this.tdp = twoDPlot;
		this.proj = tdp.getProject();
//		this.samples = twoDPlot.getSamples();
//		this.markerData = twoDPlot.getMarkerData();
//		this.sampleData = twoDPlot.getSampleData();
		locLookup = new Hashtable<String,IntVector>();//??? zx
		linkerData = new String[0][0];
//		this.updateQcPanel = true;//zx
		this.setAxisFontSize(24);
//		setColorScheme(DEFAULT_COLORS);

		// taken care of in AbstractPanel constructor
//		addMouseListener(this);
//		addMouseMotionListener(this);
//		addComponentListener(this);
		setZoomable(true, true);

//		this.names = names;
//		this.hash = hash;
//		this.finalImage = false;
//		this.classCounts = new CountVector();

		setColorScheme(DEFAULT_COLORS);

//		sampleList = HashVec.getKeys(hash);
		
		setNullMessage("Select two variables to plot");
		
		sampleData = null;
		if (proj == null) {
			markerLookup = new MarkerLookup(new Hashtable<String, String>());
		} else {
			if (Files.exists(proj.MARKERLOOKUP_FILENAME.getValue(false, false), proj.JAR_STATUS.getValue())) {
				markerLookup = proj.getMarkerLookup();
				proj.getLog().report("Marker data is available for this project");
			} else {
				markerLookup = new MarkerLookup(new Hashtable<String, String>());
			}
			if (Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false), proj.JAR_STATUS.getValue())) {
				sampleData = proj.getSampleData(1, false);
				proj.getLog().report("Sample lookup is available for this project");
			}
		}
		
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

//	public void generatePoints() {
//		float[][] currentData;
//		CountVector uniqueValueCounts;
//		byte type;
////		String[] twoDPlot.get
//
//		currentData= tdp.getDataSelected(true);
//		uniqueValueCounts = new CountVector();
////		sampleData.getClass();
//
//		points = new PlotPoint[currentData.length];
//		for (int i = 0; i < points.length; i++) {
//			if (Float.isNaN(currentData[i][0]) || Float.isNaN(currentData[i][1])) {
//				type = PlotPoint.NOT_A_NUMBER;
//				uniqueValueCounts.add("0");
////			} else if (alleleCounts[i]==-1) {
////				type = PlotPoint.MISSING;
////				uniqueValueCounts.add("0");
//			} else {
//				type = PlotPoint.FILLED_CIRCLE;
//				uniqueValueCounts.add((byte)currentData[i][2] + "");
//			}
//
//			if (swapAxes) {
//				points[i] = new PlotPoint(i+"", type, currentData[i][1], currentData[i][0], (byte)5, (byte)currentData[i][2], (byte)0);
//			} else {
//				points[i] = new PlotPoint(i+"", type, currentData[i][0], currentData[i][1], (byte)5, (byte)currentData[i][2], (byte)0);
//			}
//
//		}
//		
//		tdp.updateColorKey(uniqueValueCounts.convertToHash());
//	}

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
			setZoomable(false, false);
//			tdp.setPointSize((byte) 0);
			points = new PlotPoint[0];//currentData.size()];
			for (int i = 0; i < points.length; i++) {
				points[i] = new PlotPoint("" + Float.parseFloat(currentData.get(i)[1]), PlotPoint.FILLED_SQUARE, Float.parseFloat(currentData.get(i)[1]), Float.parseFloat(currentData.get(i)[2]), (byte) 0, (byte) 0, (byte) 0);
			}
			generateRectangles();
			return;
		} else {
			setZoomable(true, true);
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
			for (String miss : TwoDPlot.MISSING_VALUES) {
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
				if (Float.isNaN(xAxisValue) || Float.isNaN(yAxisValue)) {
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
	    org.genvisis.stats.Histogram hist;
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
    		hist = new org.genvisis.stats.Histogram(dataArray);
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
			rectangles[i] = new GenericRectangle(startX, startY, stopX, stopY, (byte) 5, true, false, (byte) 0, (byte) 2, (byte) 0, false);

//			for (int j = 0; j < linkerData[i].length; j ++) {
//				linkerData[i][j] = line[j + index];
//			}
			
		}
		
	}
	


//	private void generateRectangles() {
//	    Vector<String[]> currentData;
//	    boolean includeColorKeyValue;
//	    int index;
//	    String[] line;
//	    float xAxisValue, yAxisValue;
//	    CountVector uniqueValueCounts;
//	    
//	    includeColorKeyValue = true;
//	    currentData = tdp.getDataSelected(includeColorKeyValue);
//	    uniqueValueCounts = new CountVector();
//	    index = includeColorKeyValue? 4 : 3;
//	    if (currentData.size() > 0) {
//	        linkerData = new String[currentData.size()][currentData.elementAt(0).length - index];
//	    } else {
//	        rectangles = new GenericRectangle[0];
//	        return;
//	    }
//	    
//	    double minDiff = 9999.0;
//	    double val1 = Double.parseDouble(currentData.get(0)[1]); // TODO could be missing
//	    double min = val1;
//	    double max = val1;
//	    for (int i = 1; i < currentData.size(); i++) {
//	        double val2 = Double.parseDouble(currentData.get(i)[1]);
//	        if (Math.abs(val1 - val2) < minDiff) {
//	            minDiff = Math.abs(val1 - val2);
//	        }
//	        min = Math.min(min, Math.min(val1, val2));
//	        max = Math.max(max, Math.max(val1, val2));
//	    }
//	    
//	    int sig = ext.getNumSigFig(minDiff);
//	    if (sig > 0) { 
//	        minDiff = ext.roundToSignificantFigures(minDiff, sig);
//	    }
//	    float binHalf = (float) (currentData.size() > 1 ? ext.roundToSignificantFigures(minDiff / 2.0, sig + 1) : DEFAULT_HALF_BIN_SIZE);
//	    
////		forcePlotXmax = (float) (max + (minDiff / 2.0));
////		forcePlotXmin = (float) (min - (minDiff / 2.0));
//	    forcePlotXmax = (float) (max + minDiff);
//	    forcePlotXmin = (float) (min - minDiff);
//	    
//	    rectangles = new GenericRectangle[currentData.size()];
//	    for (int i = 0; i < rectangles.length; i++) {
//	        line = currentData.get(i);
//	        boolean missing = false;
//	        
//	        for (String miss : TwoDPlot.MISSING_VALUES) {
//	            if (miss.equals(line[1]) || miss.equals(line[2])) {
//	                missing = true;
//	                break;
//	            }	
//	        }
//	        if (missing) {
//	            xAxisValue = Float.NaN;
//	            yAxisValue = Float.NaN;
//	            uniqueValueCounts.add("0");
//	        } else {
//	            xAxisValue = Float.parseFloat(line[1]);
//	            yAxisValue = Float.parseFloat(line[2]);
//	            if (Float.isNaN(xAxisValue) || Float.isNaN(xAxisValue)) {
//	                uniqueValueCounts.add("0");
//	                //			} else if (alleleCounts[i]==-1) {
//	                //				type = PlotPoint.MISSING;
//	                //				uniqueValueCounts.add("0");
//	            } else {
//	                uniqueValueCounts.add(line[3]);
//	            }
//	        }
//	        float startX, startY, stopX, stopY;
//	        if (swapAxes) {
//	            startX = 0f;
//	            startY = sig > 0 ? (float) ext.roundToSignificantFigures((float)(xAxisValue - binHalf), sig) : (float)(xAxisValue - binHalf);
//	            stopX = yAxisValue;
//	            stopY = sig > 0 ? (float) ext.roundToSignificantFigures((float)(xAxisValue + binHalf), sig) : (float)(xAxisValue + binHalf);
//	        } else {
//	            startX = sig > 0 ? (float) ext.roundToSignificantFigures((float)(xAxisValue - binHalf), sig) : (float)(xAxisValue - binHalf);
//	            startY = 0f; 
//	            stopX = sig > 0 ? (float) ext.roundToSignificantFigures((float)(xAxisValue + binHalf), sig) : (float)(xAxisValue + binHalf);
//	            stopY = yAxisValue;
//	        }
//	        rectangles[i] = new GenericRectangle(startX, startY, stopX, stopY, (byte) 5, true, false, Byte.parseByte(line[3]), (byte) 2, (byte) 0);
//	        
//	        for (int j = 0; j < linkerData[i].length; j ++) {
//	            linkerData[i][j] = line[j + index];
//	        }
//	        
//	    }
//	    
//	    
//	}



//	public byte determineCodeFromClass(int currentClass, byte alleleCount, IndiPheno indi, byte chr, int position) {
//		int[] classes, indices;
//		CNVariant[] segs;
//		int index;
//		
//		indices = sampleData.getClassCategoryAndIndex(currentClass);
//		switch (indices[0]) {
//        case 0:
//			if (SampleData.BASIC_CLASSES[indices[1]].equals("All")) {
//				return 0;
//			} else if (SampleData.BASIC_CLASSES[indices[1]].equals("Genotype")) {
//				return (byte)(alleleCount+1);
//			} else {
//				return 0;
//			}
//        case 1:
//    		classes = indi.getClasses();
//			if (classes[indices[1]] == Integer.MIN_VALUE) {
//				return -1;
//			} else {
//				return (byte)classes[indices[1]];
//			}
//        case 2:
//			segs = indi.getCNVs(indices[1], chr);
//			if (segs == null) {
//				return 0;
//			} else {
//				index = Segment.binarySearchForOverlap(new Segment((byte)-1, position, position), segs); 
//				if (index == -1) {
//					return 0;
//				} else {
//					return (byte)(segs[index].getChr()+1);
//				}
//			}
//        default:
//        	System.err.println("Error - invalid class index");
//        	return 0;
//        }
//	}	


//	public void mouseMoved(MouseEvent event) {
//		Graphics g = getGraphics();
//		IntVector iv;
//		String pos;
//		int x, y, dataPointIndex;
//		byte size, xFontSize;
//
//		if (getFinalImage()) {
//			x = event.getX();
//			y = event.getY();
//
//			canvasSectionMinimumX = WIDTH_Y_AXIS;
//			canvasSectionMaximumX = getWidth() - WIDTH_BUFFER;
//			canvasSectionMinimumY = HEIGHT_X_AXIS;
//			canvasSectionMaximumY = getHeight() - HEAD_BUFFER;
//			pos = (int)Math.floor(x / LOOKUP_RESOLUTION) + "x" + (int)Math.floor(y / LOOKUP_RESOLUTION);
//			if (!pos.equals(prevPos)) {
//				repaint();
//			}
//
//			iv = lookupNearbyPoints(x, y, pos);
//			prox = new IntVector();
//
//			size = SIZE * 2;
//			xFontSize = (byte)(size*2);
//			g.setColor(Color.RED);
//			for (int i = 0; iv!=null && i<iv.size(); i++) {
//				dataPointIndex = iv.elementAt(i);
//				if (Distance.euclidean(new int[] {x, y}, new int[] {getX(points[dataPointIndex].getRawX()), getY(points[dataPointIndex].getRawY())}) < HIGHLIGHT_DISTANCE) {
//					g.setColor(Color.YELLOW);
//					prox.add(dataPointIndex);
//					g.fillOval(getX(points[dataPointIndex].getRawX()) - size/2, getY(points[dataPointIndex].getRawY()) - size/2, size, size);
//
//					// } else {
//					// g.setColor(Color.BLACK);
//					// g.fillOval(getX(data[iv.elementAt(i)][0])-SIZE/2,
//					// getY(data[iv.elementAt(i)][1])-SIZE/2, SIZE, SIZE);
//				}
//			}
//
//			prevPos = pos;
//		}
//	}

	// Begin of original section
	/*
	public void mouseMoved(MouseEvent event) {
		Graphics g = getGraphics();
		IntVector iv;
		String pos;
		int x, y;

		float[][] datapoints;
		IndiPheno indi;
//		Hashtable<String,IndiPheno> sampleHash = sampleData.getSampleHash();
		float[] gcScores;
		byte[] alleleCounts;
		float gcThreshold;
		int xWidth;
		int plotType, currentClass;
		int i;
		byte chr;
		int position;
		int markerIndex;
		byte size, xFontSize;

		x = event.getX();
		y = event.getY();

		canvasSectionMinimumX = WIDTH_Y_AXIS;
		canvasSectionMaximumX = getWidth()-WIDTH_BUFFER;
		canvasSectionMinimumY = HEIGHT_X_AXIS;
		canvasSectionMaximumY = getHeight()-HEAD_BUFFER;
		pos = (int)Math.floor(x/DEFAULT_LOOKUP_RESOLUTION)+"x"+(int)Math.floor(y/DEFAULT_LOOKUP_RESOLUTION);
		if (!pos.equals(prevPos)) {
			repaint();
		}
		iv = locLookup.get(pos);
		prox = new IntVector();

		plotType = sp.getPlotType();
		currentClass = sp.getCurrentClass();
		markerIndex = sp.getMarkerIndex();
		datapoints = markerData[markerIndex].getDatapoints(plotType);
		gcScores = markerData[markerIndex].getGCs();
		alleleCounts = markerData[markerIndex].getAB_Genotypes();
		chr = markerData[markerIndex].getChr();
		position = markerData[markerIndex].getPosition();
		gcThreshold = sp.getGCthreshold();

		size = sp.getPointSize();
		xFontSize = (byte)(size*2);

		g.setFont(new Font("Arial", 0, (int)(xFontSize*1.5)));
		xWidth = g.getFontMetrics(g.getFont()).stringWidth("X");

		for (int l = 0; iv!=null&&l<iv.size(); l++) {
			i = iv.elementAt(l);
			if (Distance.euclidean(new int[] {x, y}, new int[] {getX(datapoints[0][i]), getY(datapoints[1][i])})<Math.sqrt(size*size/2)) {
				indi = sampleHash.get(samples[i]);
				if (indi!=null) {
					g.setColor(colorScheme[determineCodeFromClass(currentClass, alleleCounts[i], indi, chr, position)]);
					if (gcScores[i]<gcThreshold) {
						g.drawString("X", getX(datapoints[0][i])-xWidth/2, getY(datapoints[1][i])+(int)(xFontSize*1.5/2.0)-1);
					} else {
						g.fillOval(getX(datapoints[0][i])-(int)(size*1.5)/2, getY(datapoints[1][i])-(int)(size*1.5)/2, (int)(size*1.5), (int)(size*1.5));
					}
				}

				prox.add(i);
			}
		}

		prevPos = pos;
	}
	// End of original section
	*/

//    public void mousePressed(MouseEvent e) {}
//
//    public void mouseReleased(MouseEvent e) { }
//
//    public void mouseDragged(MouseEvent e) { }

	public void mouseClicked(MouseEvent e) {
		JPopupMenu menu;
		int[] linkKeyIndicies;
//		Vector<String[]> linkKeyValues;
//		boolean scatter, trailer;
		String[] ids = null;
		String markerName;
		String sample, region, region2;
		int[] positions;
		byte maxNumPoints;

//		if (e.getButton()==MouseEvent.BUTTON1) { // left click
//		} else if (e.getButton()==MouseEvent.BUTTON3) { // right click
//		}

		linkKeyIndicies = tdp.getCurrentLinkKeyColumnIndices();
//		linkKeyValues = tdp.getCurrentLinkKeyValues();
//		if (linkKeyValues == null) {
//			return;
//		}
		
		if (prox != null && prox.size() > 0) {
			menu = new JPopupMenu();
			maxNumPoints = (byte) Math.min(20, prox.size());
			for (int i = 0; i < maxNumPoints; i++) {
				String[] linkerDataElem = linkerData[prox.elementAt(i)];
				
				if (linkKeyIndicies[TwoDPlot.MARKER_INDEX_IN_LINKERS] >= 0) {
					markerName = linkerDataElem[TwoDPlot.MARKER_INDEX_IN_LINKERS];
					if (markerLookup.get(markerName) != null) {
						menu.add(new LaunchAction(proj, markerName, Color.CYAN));
					}
				}
				
				sample = null;
				// TODO this check will ALWAYS fail!
				if (linkKeyIndicies[TwoDPlot.DNA_INDEX_IN_LINKERS] >= 0 && Files.exists(proj.SAMPLE_DIRECTORY.getValue(false, false) + sample + Sample.SAMPLE_FILE_EXTENSION, proj.JAR_STATUS.getValue())) {
					sample = linkerDataElem[TwoDPlot.DNA_INDEX_IN_LINKERS];
				}
				if (sample == null && sampleData != null) { // if Sample not already identified and if a sample lookup exists
					ids = null;
					if (linkKeyIndicies[TwoDPlot.FID_INDEX_IN_LINKERS] >= 0) { // if FID present
						ids = sampleData.lookup(linkerDataElem[TwoDPlot.FID_INDEX_IN_LINKERS] + "\t" + linkerDataElem[TwoDPlot.IID_INDEX_IN_LINKERS]);
					}
					if (ids == null) {
						ids = sampleData.lookup(linkerDataElem[TwoDPlot.IID_INDEX_IN_LINKERS]);
					}
					if (ids != null && Files.exists(proj.SAMPLE_DIRECTORY.getValue(false, false) + ids[0] + Sample.SAMPLE_FILE_EXTENSION, proj.JAR_STATUS.getValue())) {
						sample = ids[0];
					}
				}

				positions = new int[] {-1,-1,-1};
				region = null;
				region2 = null;
				String[][] metaData = tdp.getCurrentColumnMetaData();
				if (metaData != null && metaData.length > 0 && !(metaData[0] == null && (metaData.length > 1 ? metaData[1] == null : true))) {
				    if (metaData[0] != null) {
				        int[] tempPositions = new int[3];
				        tempPositions[0] = Positions.chromosomeNumber(metaData[0][0]);
				        tempPositions[1] = Integer.parseInt(metaData[0][2]);
				        tempPositions[2] = Integer.parseInt(metaData[0][3]);
				        region = Positions.getUCSCformat(tempPositions);
//							if (sample != null && tempRegion != null) {
//								menu.add(new LaunchAction(proj, sample, tempRegion, Color.GRAY));
//							}
				    }
				    if (metaData.length > 1 && metaData[1] != null) {
				        positions[0] = Positions.chromosomeNumber(metaData[1][0]);
				        positions[1] = Integer.parseInt(metaData[1][2]);
				        positions[2] = Integer.parseInt(metaData[1][3]);
				        region2 = Positions.getUCSCformat(positions);
				    }
				} else {
    				if (linkKeyIndicies[TwoDPlot.REGION_INDEX_IN_LINKERS] >= 0) {
    					region = linkerDataElem[TwoDPlot.REGION_INDEX_IN_LINKERS];
    				} else if (linkKeyIndicies[TwoDPlot.CHR_INDEX_IN_LINKERS] >= 0) {
    					positions[0] = Positions.chromosomeNumber(linkerDataElem[TwoDPlot.CHR_INDEX_IN_LINKERS]);
    					if (positions[0] != -1) {
    						if (linkKeyIndicies[TwoDPlot.POS_INDEX_IN_LINKERS] >= 0) {
    							try {
    								positions[1] = Integer.parseInt(linkerDataElem[TwoDPlot.POS_INDEX_IN_LINKERS]);
    							} catch (NumberFormatException nfe) {}
    						}
    						if (linkKeyIndicies[TwoDPlot.STOP_POS_INDEX_IN_LINKERS] >= 0) {
    							try {
    								positions[2] = Integer.parseInt(linkerDataElem[TwoDPlot.STOP_POS_INDEX_IN_LINKERS]);
    							} catch (NumberFormatException nfe) {}
    						}
    						region = Positions.getUCSCformat(positions);
    					}
    				}
				}
				
				Color sampColor = null;
				if (sampleData != null && ids != null) {
					sampColor = this.colorScheme[sampleData.determineCodeFromClass(tdp.colorKeyPanel.getCurrentClass(), (byte) 0, sampleData.getIndiFromSampleHash(ids[0]), (byte) 0, 0)];
				}
				if (sample != null) {
					if (region != null) {
						if(region2 != null && !region.equals(region2)) {
							menu.add(new LaunchAction(proj, sample, new String[]{region, region2}, sampColor == null ? Color.GRAY : sampColor));
						} else {
							menu.add(new LaunchAction(proj, sample, region, sampColor == null ? Color.GRAY : sampColor));
						}
					} else {
						menu.add(new LaunchAction(proj, sample, Trailer.DEFAULT_LOCATION, sampColor));
					}
				}
			}
			menu.show(this, e.getX(), e.getY());
		}
	}

//	public void mouseEntered(MouseEvent e) {}
//
//	public void mouseExited(MouseEvent e) {}
//
//	public void mousePressed(MouseEvent e) {}
//
//	public void mouseReleased(MouseEvent e) {}

//	public void componentHidden(ComponentEvent e) {}
//
//	public void componentMoved(ComponentEvent e) {}
//
//	public void componentResized(ComponentEvent e) {
//		paintAgain();
//	}
//
//	public void componentShown(ComponentEvent e) {}

//	public static void main(String[] args) {
//		ScatterPlot.main(new String[] {"-notJar"});
//	}

	public SampleData getSampleData(){
		return sampleData;
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

    public void setHistogram(org.genvisis.stats.Histogram hist) {
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
