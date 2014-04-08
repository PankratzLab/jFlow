package cnv.plots;

import java.awt.Color;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JPopupMenu;

import cnv.filesys.MarkerLookup;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.gui.LaunchAction;
import cnv.var.SampleData;
import common.CountVector;
import common.Files;
import common.IntVector;
import common.Logger;
import common.Positions;

public class TwoDPanel extends AbstractPanel implements MouseListener, MouseMotionListener {
	public static final long serialVersionUID = 3L;
	public static final int LOOKUP_RESOLUTION = 20;
	public static final Color[] DEFAULT_COLORS = {new Color(33, 31, 53), // dark dark
			   									  new Color(23, 58, 172), // dark blue
			   									  new Color(201, 30, 10), // deep red
			   									  new Color(140, 20, 180), // deep purple
			   									  new Color(33, 87, 0), // dark green
			   									  new Color(55, 129, 252), // light blue
			   									  new Color(94, 88, 214), // light purple
			   									  new Color(189, 243, 61), // light green
			   									  new Color(217, 109, 194), // pink
			   									new Color(0, 0, 128), // ALL KINDS OF BLUES
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
	
	protected TwoDPlot tdp;
	IntVector indicesOfNearbySamples;
	private boolean updateQcPanel;
	private boolean swapAxes;
	private Logger log;
	private MarkerLookup markerLookup;
	private SampleData sampleData;
	private Project proj;
	private String[][] setOfKeys;

	public TwoDPanel(TwoDPlot twoDPlot, Logger log) {
		super();
		
		this.tdp = twoDPlot;
		this.log = log;
		this.proj = tdp.getProject();
//		this.samples = twoDPlot.getSamples();
//		this.markerData = twoDPlot.getMarkerData();
//		this.sampleData = twoDPlot.getSampleData();
		locLookup = new Hashtable<String,IntVector>();//??? zx
		setOfKeys = new String[0][0];
//		this.updateQcPanel = true;//zx
		
		setColorScheme(DEFAULT_COLORS);

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
			if (Files.exists(proj.getFilename(Project.MARKERLOOKUP_FILENAME, false, false), proj.getJarStatus())) {
				markerLookup = proj.getMarkerLookup();
				System.out.println("Marker data is available for this project");
			} else {
				markerLookup = new MarkerLookup(new Hashtable<String, String>());
			}
			if (Files.exists(proj.getFilename(Project.SAMPLE_DATA_FILENAME, false, false), proj.getJarStatus())) {
				sampleData = proj.getSampleData(1, false);
				System.out.println("Sample lookup is available for this project");
			}
		}
		
	}
	
	public void refreshOtherComponents() {
		tdp.refreshOtherButtons();
	}

	public void assignAxisLabels() {
		displayXaxis = displayYaxis = true;
		xAxisLabel = swapAxes? tdp.getNamesSelected()[1]:tdp.getNamesSelected()[0];
		yAxisLabel = swapAxes? tdp.getNamesSelected()[0]:tdp.getNamesSelected()[1];
	}
	
	public boolean invertX() {
		return false;
	}

	public boolean invertY() {
		return false;
	}
	
	public void toggleMasking() {
		
	}	
	
	public void highlightPoints() {
		byte defaultSize;
		
		defaultSize = tdp.getPointSize();
		for (int i = 0; i < points.length; i++) {
			if (points[i].isHighlighted()) {
				points[i].setSize((byte)(defaultSize*1.5));
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
		Vector<String[]> currentData;
		CountVector uniqueValueCounts;
		boolean includeColorKeyValue;
		byte type;
		String[] line;
		float xAxisValue, yAxisValue;
		byte index;

		includeColorKeyValue = true;
		currentData= tdp.getDataSelected(includeColorKeyValue);
		uniqueValueCounts = new CountVector();
		
		points = new PlotPoint[currentData.size()];
		index = (byte) (includeColorKeyValue? 4 : 3);
		if (currentData.size()>0) {
			setOfKeys = new String[currentData.size()][currentData.elementAt(0).length - index];
		}
		for (int i = 0; i < points.length; i++) {
			line = currentData.elementAt(i);
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

			if (swapAxes) {
				points[i] = new PlotPoint(line[0], type, yAxisValue, xAxisValue, (byte)5, Byte.parseByte(line[3]), (byte)0);
			} else {
				points[i] = new PlotPoint(line[0], type, xAxisValue, yAxisValue, (byte)5, Byte.parseByte(line[3]), (byte)0);
			}

			for (int j = 0; j < setOfKeys[i].length; j ++) {
				setOfKeys[i][j] = line[j + index];
			}
		}
		
		tdp.updateColorKey(uniqueValueCounts.convertToHash());
	}





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

    public void mousePressed(MouseEvent e) {}

    public void mouseReleased(MouseEvent e) { }

    public void mouseDragged(MouseEvent e) { }

	public void mouseClicked(MouseEvent e) {
		JPopupMenu menu;
		int[] linkKeyIndicies;
//		Vector<String[]> linkKeyValues;
//		boolean scatter, trailer;
		String[] ids;
		String markerName;
		String sample, region;
		int[] positions;
		byte maxNumPoints;

//		if (e.getButton()==MouseEvent.BUTTON1) { // left click
//		} else if (e.getButton()==MouseEvent.BUTTON3) { // right click
//		}

		linkKeyIndicies = tdp.getCurrentLinkKeyColumnLabels();
//		linkKeyValues = tdp.getCurrentLinkKeyValues();
//		if (linkKeyValues == null) {
//			return;
//		}
		
		if (prox != null && prox.size() > 0) {
			menu = new JPopupMenu();
			maxNumPoints = (byte) Math.min(20, prox.size());
			for (int i = 0; i < maxNumPoints; i++) {
//				menu.add(new LaunchAction(linkKeyValues[prox.elementAt(i)][0] + "\t" + points[prox.elementAt(i)].getRawX() + "\t" + points[prox.elementAt(i)].getRawY(), true));
				menu.add(new LaunchAction(points[prox.elementAt(i)].getId() + "\t" + points[prox.elementAt(i)].getRawX() + "\t" + points[prox.elementAt(i)].getRawY(), true));

				if (linkKeyIndicies[3] >= 0) {
					markerName = setOfKeys[prox.elementAt(i)][3];
					if (markerLookup.get(markerName) != null) {
						menu.add(new LaunchAction(proj, markerName, Color.CYAN));
					}
				}
				
				sample = null;
				if (linkKeyIndicies[2] >= 0 && Files.exists(proj.getDir(Project.SAMPLE_DIRECTORY, false, log, false) + sample + Sample.SAMPLE_DATA_FILE_EXTENSION, proj.getJarStatus())) {
					sample = setOfKeys[prox.elementAt(i)][2];
				}
				if (sample == null && sampleData != null) { // if Sample not already identified and if a sample lookup exists
					ids = null;
					if (linkKeyIndicies[1] >= 0) { // if FID present
						ids = sampleData.lookup(setOfKeys[prox.elementAt(i)][1] + "\t" + setOfKeys[prox.elementAt(i)][0]);
					}
					if (ids == null) {
						ids = sampleData.lookup(setOfKeys[prox.elementAt(i)][0]);
					}
					if (ids != null && Files.exists(proj.getDir(Project.SAMPLE_DIRECTORY, false, log, false) + ids[0] + Sample.SAMPLE_DATA_FILE_EXTENSION, proj.getJarStatus())) {
						sample = ids[0];
					}
				}

				positions = new int[] {-1,-1,-1};
				region = null;
				if (linkKeyIndicies[4] >= 0) {
					region = setOfKeys[prox.elementAt(i)][4];
				} else if (linkKeyIndicies[5] >= 0) {
					positions[0] = Positions.chromosomeNumber(setOfKeys[prox.elementAt(i)][5]);
					if (positions[0] != -1) {
						if (linkKeyIndicies[6] >= 0) {
							try {
								positions[1] = Integer.parseInt(setOfKeys[prox.elementAt(i)][6]);
							} catch (NumberFormatException nfe) {
							}
						}
						if (linkKeyIndicies[7] >= 0) {
							try {
								positions[2] = Integer.parseInt(setOfKeys[prox.elementAt(i)][7]);
							} catch (NumberFormatException nfe) {
							}
						}
						region = Positions.getUCSCformat(positions);
					}
				}
				
				if (sample != null && region != null) {
					menu.add(new LaunchAction(proj, sample, region, Color.GRAY));
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
	
	public void setSwapAxes(boolean swapAxes) {
		this.swapAxes = swapAxes;
	}

	public boolean isSwapAxes() {
		return swapAxes;
	}
}
