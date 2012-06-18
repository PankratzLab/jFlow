/*
 * Revised in Feb, 2012. Add mouse dragged feature.
 */
package cnv.plots;

//import java.util.Date;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.Hashtable;

import javax.swing.JOptionPane;
import javax.swing.JPopupMenu;

import cnv.filesys.MarkerData;
import cnv.filesys.Project;
import cnv.gui.LaunchAction;
import cnv.var.CNVariant;
import cnv.var.SampleData;
import cnv.var.IndiPheno;
import filesys.Segment;
//import common.CountVector;
import common.CountVector;
import common.HashVec;
import common.IntVector;
import common.Sort;
//import common.ext;
//import mining.Distance;
import common.Array;
//import mining.Distance;

//public class ScatterPanel extends AbstractPanel implements MouseListener, MouseMotionListener, ComponentListener {
public class TwoDPanel extends AbstractPanel2 implements MouseListener, MouseMotionListener {
	public static final long serialVersionUID = 3L;
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
	
//	protected MarkerData[] markerData;
//	byte[] alleleCounts;				//zx
	protected TwoDPlot tdp;
//	protected String[] samples;
	protected IntVector prox;
	protected SampleData sampleData;
	IntVector indeciesOfNearbySamples;	//zx
	private boolean updateQcPanel;		//zx: A control variable. Do not update QcPanel when resizing, or etc.
	private int mouseStartX ;
	private int mouseStartY ;
	private int mouseEndX ;
	private int mouseEndY ;
	private boolean swapAxes;

	public TwoDPanel(TwoDPlot twoDPlot) {
		super();
		
		this.tdp = twoDPlot;
//		this.samples = twoDPlot.getSamples();
//		this.markerData = twoDPlot.getMarkerData();
		this.sampleData = twoDPlot.getSampleData();
		locLookup = new Hashtable<String,IntVector>();//??? zx
//		this.updateQcPanel = true;//zx
		
		setColorScheme(DEFAULT_COLORS);

//		addMouseListener(this);
//		addMouseMotionListener(this);
		addComponentListener(this);
		setZoomable(true, true);

//		this.names = names;
//		this.hash = hash;
//		this.finalImage = false;
//		this.classCounts = new CountVector();

		setColorScheme(DEFAULT_COLORS);

//		sampleList = HashVec.getKeys(hash);
		
		setNullMessage("Select two variables to plot");
		
		addMouseListener(this);
		addMouseMotionListener(this);
		addComponentListener(this);
		setZoomable(true, true);
	}
	
	public void refreshOtherComponents() {
		tdp.refreshOtherButtons();
	}

	public void assignAxisLabels() {
		displayXaxis = displayYaxis = true;
		xAxisLabel = swapAxes? tdp.getNamesSelected()[1]:tdp.getNamesSelected()[0];
		yAxisLabel = swapAxes? tdp.getNamesSelected()[0]:tdp.getNamesSelected()[1];
		
//		xAxisLabel = MarkerData.TYPES[sp.getPlotType()][swapAxes?0:1];
//		yAxisLabel = MarkerData.TYPES[sp.getPlotType()][swapAxes?1:0];

//		xAxisLabel = swapAxes?
//						tdp.getNamesHash().keySet().toArray()[tdp.getCurrentPair()[0][0]][tdp.getCurrentPair()[0][1]]:
//						tdp.getNamesHash().keySet().toArray()[tdp.getCurrentPair()[1][0]][tdp.getCurrentPair()[1][1]];
//		yAxisLabel = swapAxes?
//						tdp.getNamesHash().keys().toString()[tdp.getCurrentPair()[1][0]][tdp.getCurrentPair()[1][1]]:
//						tdp.getNamesHash().keys().toString()[tdp.getCurrentPair()[0][0]][tdp.getCurrentPair()[0][1]];

		/*
		int[][] currentPair;
		
		currentPair = sp.getCurrentPair();
		if (names.length == 0) {
			sp.setDescriptor("Error - no files with .mds extension were present in the project directory");
			displayXaxis = displayYaxis = false;
		} else if (currentPair[0][0] == -1 || currentPair[1][0] == -1) {
			sp.setDescriptor("Double click on one of the files within the box on the left-hand side and select two of its components");
			if (currentPair[0][0] == -1) {
				displayXaxis = false;
			} else {
				xAxisLabel = names[currentPair[0][0]][0]+"_"+names[currentPair[0][0]][currentPair[0][1]+1];
				displayXaxis = true;
			}
			if (currentPair[1][0] == -1) {
				displayYaxis = false;
			} else {
				yAxisLabel = names[currentPair[1][0]][0]+"_"+names[currentPair[1][0]][currentPair[1][1]+1];
				displayYaxis = true;
			}
		} else {
			xAxisLabel = names[currentPair[0][0]][0]+"_"+names[currentPair[0][0]][currentPair[0][1]+1];
			yAxisLabel = names[currentPair[1][0]][0]+"_"+names[currentPair[1][0]][currentPair[1][1]+1];
			sp.setDescriptor(xAxisLabel+" x "+yAxisLabel);
			displayXaxis = displayYaxis = true;
		}
		*/

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

	public void generatePoints() {
		float[][] currentData = tdp.getDataSelected();
		points = new PlotPoint[currentData.length];
		for (int i = 0; i < points.length; i++) {
			if (swapAxes) {
				points[i] = new PlotPoint(i+"", PlotPoint.FILLED_CIRCLE, currentData[i][1], currentData[i][0], (byte)5, (byte)0, (byte)0);
			} else {
				points[i] = new PlotPoint(i+"", PlotPoint.FILLED_CIRCLE, currentData[i][0], currentData[i][1], (byte)5, (byte)0, (byte)0);
			}
		}
		if (points.length > 0) {
			System.out.println(points.length);
		}

//		CountVector classCounts;
//
//		if (sampleData.getActualClassName(currentClass).equals("Site")) {
//			setColorScheme(BLUES);
//		} else {
//			setColorScheme(DEFAULT_COLORS);
//		}
//		
//		points = new PlotPoint[sampleList.length]; 
//		classCounts.clear();
//		for (int i = 0; i<sampleList.length; i++) {
//			data = hash.get(sampleList[i]);
//			if (data[currentPair[0][0]] != null && data[currentPair[1][0]] != null && !Float.isNaN(data[currentPair[0][0]][currentPair[0][1]])  && !Float.isNaN(data[currentPair[1][0]][currentPair[1][1]])) {
//				trav = sampleData.lookup(sampleList[i]);
//				if (trav == null) {
//					System.err.println("Error - could not look up "+sampleList[i]); // looks up any individual present in any .mds file that was loaded, even those not in the current file
//					tagalong = true;
//					colorCode = 0;
//				} else {
//					tagalong = false;
//					colorCode = sampleData.getClassForInd(trav, currentClass);
//					if (colorCode == -1 && !sp.maskMissing()) {
//						colorCode = 0;
//					}
//				}
//
//				points[i] = new PlotPoint(sampleList[i], PlotPoint.FILLED_CIRCLE, data[currentPair[0][0]][currentPair[0][1]], data[currentPair[1][0]][currentPair[1][1]], tagalong?SIZE_TAGALONGS:SIZE, colorCode, (byte)(colorCode>0?1:0));
//				classCounts.add(colorCode+"");
//			}
//		}
//		sp.updateColorKey(getClassCounts().convertToHash());

	}

	public byte determineCodeFromClass(int currentClass, byte alleleCount, IndiPheno indi, byte chr, int position) {
		int[] classes, indices;
		CNVariant[] segs;
		int index;
		
		indices = sampleData.getClassCategoryAndIndex(currentClass);
		switch (indices[0]) {
        case 0:
			if (SampleData.BASIC_CLASSES[indices[1]].equals("All")) {
				return 0;
			} else if (SampleData.BASIC_CLASSES[indices[1]].equals("Genotype")) {
				return (byte)(alleleCount+1);
			} else {
				return 0;
			}
        case 1:
    		classes = indi.getClasses();
			if (classes[indices[1]] == Integer.MIN_VALUE) {
				return -1;
			} else {
				return (byte)classes[indices[1]];
			}
        case 2:
			segs = indi.getCNVs(indices[1], chr);
			if (segs == null) {
				return 0;
			} else {
				index = Segment.binarySearchForOverlap(new Segment((byte)-1, position, position), segs); 
				if (index == -1) {
					return 0;
				} else {
					return (byte)(segs[index].getChr()+1);
				}
			}
        default:
        	System.err.println("Error - invalid class index");
        	return 0;
        }
	}	


	public void mouseMoved(MouseEvent event) {
		Graphics g = getGraphics();
		String pos;
		int x, y;

		float[][] datapoints;
		IndiPheno indi;
		float[] gcScores;
//		byte[] alleleCounts;
//		float gcThreshold;
		int xWidth;
		int plotType, currentClass;
		int i;
		byte chr;
		int position;
		int markerIndex;
		byte size, xFontSize;
		
//		System.out.println("mouse moved");
		
//		
//		//IntVector indeciesOfDataPoint;//zx
//
//		x = event.getX();
//		y = event.getY();
//
//		canvasSectionMinimumX = WIDTH_Y_AXIS;
//		canvasSectionMaximumX = getWidth()-WIDTH_BUFFER;
//		canvasSectionMinimumY = HEIGHT_X_AXIS;
//		canvasSectionMaximumY = getHeight()-HEAD_BUFFER;
//		pos = (int)Math.floor(x/DEFAULT_LOOKUP_RESOLUTION)+"x"+(int)Math.floor(y/DEFAULT_LOOKUP_RESOLUTION);
//		if (!pos.equals(prevPos)) {
//			repaint();
//		}
//		//iv = locLookup.get(pos);
//		//indeciesOfDataPoint = lookupNearbyPoints(x, y, pos);
//		indeciesOfNearbySamples = lookupNearbyPoints(x, y, pos);
//		//System.out.println("Number of nearby samples: "+(indeciesOfNearbySamples==null?0:indeciesOfNearbySamples.size()));//zx test point
//		//prox = new IntVector();
//
//		plotType = tdp.getPlotType();
//		currentClass = tdp.getCurrentClass();
//		markerIndex = tdp.getMarkerIndex();
//		datapoints = markerData[markerIndex].getDatapoints(plotType);
//		gcScores = markerData[markerIndex].getGCs();
////		alleleCounts = markerData[markerIndex].getAB_Genotypes();
//		chr = markerData[markerIndex].getChr();
//		position = markerData[markerIndex].getPosition();
////		gcThreshold = sp.getGCthreshold();
//
//		size = tdp.getPointSize();
//		xFontSize = (byte)(size*2);
//
//		g.setFont(new Font("Arial", 0, (int)(xFontSize*1.5)));
//		xWidth = g.getFontMetrics(g.getFont()).stringWidth("X");
//
//		//System.out.println("pos: "+pos+"\t iv.size():"+(indeciesOfNearbySamples==null?"null":indeciesOfNearbySamples.size()));//zx test point
//		for (int l = 0; indeciesOfNearbySamples!=null&&l<indeciesOfNearbySamples.size(); l++) {
//			i = indeciesOfNearbySamples.elementAt(l);
//			indi = sampleData.getIndiFromSampleHash(samples[i]);
//			g.setColor(colorScheme[determineCodeFromClass(currentClass, alleleCounts[i], indi, chr, position)]);
//			//g.setColor(Color.YELLOW);
////			if (gcScores[i]<gcThreshold) {
//			if (alleleCounts[i]==-1) {
//				g.drawString("X", getX(datapoints[0][i])-xWidth/2, getY(datapoints[1][i])+(int)(xFontSize/2.0));
//			} else {
//				g.fillOval(getX(datapoints[0][i])-(int)(size*2)/2, getY(datapoints[1][i])-(int)(size*2)/2, (int)(size*2), (int)(size*2));
//			}
//		}
//		prevPos = pos;
	}

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

    public void mousePressed(MouseEvent e) {
    	mouseStartX = e.getX();
    	mouseStartY = e.getY();
    }

    public void mouseReleased(MouseEvent e) { }

    public void mouseDragged(MouseEvent e) { }

//	public void paintComponent(Graphics g) {
//		System.out.println("inner");
//		this.paintComponent(g);
////		super(g);
//	}
    
	public void mouseClicked(MouseEvent event) {
//		System.out.println("mouse location: ("+event.getPoint().x+", "+event.getPoint().y+")");
		if(event.getPoint().x>=70 && event.getPoint().x<=90 && event.getPoint().y>=(getHeight()-55) && event.getPoint().y<=(getHeight()-75)) {
			JOptionPane.showMessageDialog(null, "You've just clicked the invert button", "Message", JOptionPane.PLAIN_MESSAGE);
		}
		
		JPopupMenu menu;
		MarkerData mData;
		String markerPosition;
		int window;
		// float[][] datapoints = markerData[sp.getMarkerIndex()].getDatapoints(sp.getPlotType());

		/*
		if (event.getButton()==MouseEvent.BUTTON1) { // left click
		} else if (event.getButton()==MouseEvent.BUTTON3) { // right click
		}
		*/

		System.out.println("mouse clicked");
		
//		window = Integer.parseInt(tdp.getProject().getProperty(Project.NUM_MARKERS_PER_FILE));
//		mData = markerData[tdp.getMarkerIndex()];
//		markerPosition = "chr"+mData.getChr()+":"+(mData.getPosition()-window)+"-"+(mData.getPosition()+window);
//		if (indeciesOfNearbySamples!=null&&indeciesOfNearbySamples.size()>0) {
//			menu = new JPopupMenu();
//			for (int i = 0; i<indeciesOfNearbySamples.size(); i++) {
//				// menu.add(samples[prox.elementAt(i)] +"
//				// ("+datapoints[0][prox.elementAt(i)]+",
//				// "+datapoints[1][prox.elementAt(i)]+")");
//
//				menu.add(new LaunchAction(tdp.getProject(), samples[indeciesOfNearbySamples.elementAt(i)], markerPosition, Color.BLACK));
//			}
//			menu.show(this, event.getX(), event.getY());
//		}
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
