/*
 * Revised in Feb, 2012. Add mouse dragged feature.
 */
package cnv.plots;

//import java.util.Date;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;
import java.util.Hashtable;
import javax.swing.JPopupMenu;

import cnv.filesys.ClusterFilter;
import cnv.filesys.ClusterFilterCollection;
import cnv.filesys.MarkerData;
import cnv.filesys.Project;
import cnv.gui.LaunchAction;
import cnv.var.CNVariant;
import cnv.var.SampleData;
import cnv.var.IndiPheno;
import filesys.Segment;
import filesys.SerialHash;
//import common.CountVector;
import common.CountVector;
import common.HashVec;
import common.IntVector;
import common.Sort;
import common.ext;
//import common.ext;
//import mining.Distance;
import common.Array;
//import mining.Distance;

public class ScatterPanel extends AbstractPanel implements MouseListener, MouseMotionListener, ComponentListener {
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
	
	protected MarkerData[] markerData;
	byte[] alleleCounts;				//zx
	protected ScatterPlot sp;
	protected String[] samples;
	protected IntVector prox;
	protected SampleData sampleData;
	IntVector indeciesOfNearbySamples;	//zx
	private boolean updateQcPanel;		//zx: A control variable. Do not update QcPanel when resizing, or etc.
	private int mouseStartX ;
	private int mouseStartY ;
	private int mouseEndX ;
	private int mouseEndY ;


	public ScatterPanel(ScatterPlot sp) {
		super();
		
		this.sp = sp;
		this.samples = sp.getSamples();
		this.markerData = sp.getMarkerData();
		this.sampleData = sp.getSampleData();
//		locLookup = new Hashtable<String,IntVector>();//zx
		this.updateQcPanel = true;//zx
		
		setColorScheme(DEFAULT_COLORS);

//		addMouseListener(this);
//		addMouseMotionListener(this);
		addComponentListener(this);
		setZoomable(true, true);
	}

	public void assignAxisLabels() {
		xAxisLabel = MarkerData.TYPES[sp.getPlotType()][0];
		yAxisLabel = MarkerData.TYPES[sp.getPlotType()][1];
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
		
		defaultSize = sp.getPointSize();
		for (int i = 0; i < points.length; i++) {
			if (points[i].isHighlighted()) {
				points[i].setSize((byte)(defaultSize*1.5));
			} else {
				points[i].setSize((byte)(defaultSize));
			}
			
		}
	}

	public void generatePoints() {
		int position, markerIndex, plotType, currentClass;
		Hashtable<String,IndiPheno> sampleHash;
		byte chr, genotypeCode, sexCode, classCode, type;
		float[][] datapoints;
//		byte[] alleleCounts;
//		float gcThreshold;
		byte layer;
		float[] gcScores;
		IndiPheno indi;
//		long time;
		byte size, xFontSize;
		boolean[] displayCents;
		float[][][][] cents;
		int numCents, count;
		byte centSize;
		float x, y;
		//int[][] dataForQc;//zx
		int[] genotype;
		String[] sex;
		String[] otherClass;
		CountVector classCounts;//zx
//		ClusterFilterCollection clusterFilterCollection;//zx

//		time = new Date().getTime();

		plotType = sp.getPlotType();
		currentClass = sp.getCurrentClass();
		markerIndex = sp.getMarkerIndex();
//		gcThreshold = sp.getGCthreshold();
		datapoints = markerData[markerIndex].getDatapoints(plotType);
		gcScores = markerData[markerIndex].getGCs();
//		alleleCounts = markerData[markerIndex].getAB_Genotypes();//zx
//		alleleCounts = sp.getClusterFilterCollection().filterMarker(markerData[markerIndex], sp.getGCthreshold());
		alleleCounts = markerData[markerIndex].getAB_GenotypesAfterFilters(sp.getClusterFilterCollection(), sp.getMarkerName(), sp.getGCthreshold());//zx
		chr = markerData[markerIndex].getChr();
		position = markerData[markerIndex].getPosition();
		size = sp.getPointSize();
		xFontSize = (byte)(size*2);
		sampleHash = sampleData.getSampleHash();
		displayCents = sp.getDisplayCents();
		cents = sp.getCents();
		centSize = 20;
		
		if (datapoints[0] == null || datapoints[1] == null) {
			errorMessage = "Data not available:";
			points = new PlotPoint[0]; 
			return;
		} else {
			errorMessage = null;
		}
		
		if (plotType == 1 || plotType == 2) {
			numCents = Array.booleanArraySum(displayCents);
			points = new PlotPoint[samples.length+numCents*3];

			count = 0;
			for (int i = 0; i<displayCents.length; i++) {
				if (displayCents[i]) {
					for (int j = 0; j<3; j++) {
						if (cents[i][markerIndex][j] == null) {
							x = 0;
							y = 0;
						} else if (plotType == 1) {
							x = (float)(cents[i][markerIndex][j][1] /(1+ Math.sin(cents[i][markerIndex][j][0]*Math.PI/2)/Math.cos(cents[i][markerIndex][j][0]*Math.PI/2)));
							y = (float)(cents[i][markerIndex][j][1] /(1+ Math.cos(cents[i][markerIndex][j][0]*Math.PI/2)/Math.sin(cents[i][markerIndex][j][0]*Math.PI/2)));
						} else {
							x = cents[i][markerIndex][j][0];
							y = cents[i][markerIndex][j][1];
						}
						points[count*3+j] = new PlotPoint("Centoids", PlotPoint.FILLED_CIRCLE, x, y, centSize, (byte)(5+i), (byte)10);

	                }
					count++;
				}
	        }
		} else {
			points = new PlotPoint[samples.length];
			numCents = 0;
		}
		
		if (plotType < 2) {
			forcePlotXmax = Float.NaN;
		} else {
			forcePlotXmax = 1;
		}

		//dataForQc = new int[3][samples.length];//zx
		//genotype = markerData[markerIndex].getAB_GenotypesAfterFilters(null, sp.getGCthreshold());
		genotype = new int[samples.length];
		sex = new String[samples.length];
		otherClass = new String[samples.length];
		classCounts = new CountVector();
		for (int i = 0; i<samples.length; i++) {
			indi = sampleHash.get(samples[i]);
			if (indi!=null) {
				genotypeCode = (byte)(alleleCounts[i]+1);
//				genotypeCode = determineCodeFromClass(1, alleleCounts[i], indi, chr, position);
//				if (gcScores[i]<gcThreshold) {
//					genotypeCode = 0;
//				}
//				clusterFilterCollection = sp.getClusterFilterCollection();
//				clusterFilterCollection.filterMarker(markerData[markerIndex]);
				
				
				// additional genotypeFilters
				if (currentClass == 1) {
					classCode = genotypeCode;
				} else {
					classCode = determineCodeFromClass(currentClass, alleleCounts[i], indi, chr, position);
				}
				if (sampleData.getSexClassIndex() == -1) {
					sexCode = 0;
				} else {
					sexCode = determineCodeFromClass(sampleData.getSexClassIndex()+SampleData.BASIC_CLASSES.length, alleleCounts[i], indi, chr, position);
				}
				
				//System.out.println("Current loop:\t"+i+"\t code:\t"+code+"\t alleleCounts: "+alleleCounts[i]+"\t gcScores: "+gcScores[i]);//zx
				if (classCode == -1 && !sp.maskMissing()) {
					classCode = 0;
				}
				if (Float.isNaN(datapoints[0][i]) || Float.isNaN(datapoints[1][i])) {
					type = PlotPoint.NOT_A_NUMBER;
				} else if (alleleCounts[i]==-1) {
					type = PlotPoint.MISSING;
				} else {
					type = PlotPoint.FILLED_CIRCLE;
				}
				layer = (byte)((sampleData.getClassCategoryAndIndex(currentClass)[0]==2 && classCode > 0)?1:0);
				if (type == PlotPoint.MISSING || type == PlotPoint.NOT_A_NUMBER) {
					classCounts.add(0+"");
					genotype[i]=0;//zx
				} else {
					classCounts.add(classCode+"");
				}
				if (classCode < 0) {
					System.err.println("Error - classCOde is less than 0 ("+classCode+")");
				}
				points[numCents*3+i] = new PlotPoint(samples[i], type, datapoints[0][i], datapoints[1][i], type==PlotPoint.FILLED_CIRCLE?size:xFontSize, classCode, layer);
				genotype[i]=genotypeCode;
				//sex[i]=(sexCode==1?"Female":(sexCode==2?"Male":"Missing"));
				//for (int j=0; j<sampleData.getActualClassColorKey(0).length; j++) {
				//	if (sampleData.getActualClassColorKey(0)[j][0].equals(determineCodeFromClass(2, alleleCounts[i], indi, chr, position)+"")){
				//		sex[i]=sampleData.getActualClassColorKey(0)[j][1];
				//		break;
				//	}
				//	sex[i]="Missing";
				//}
				sex[i] = determineCodeFromClass(2, alleleCounts[i], indi, chr, position)+"";
				
				//for (int j=0; j<sampleData.getActualClassColorKey(1).length; j++) {
				//	if (sampleData.getActualClassColorKey(1)[j][0].equals(classCode+"")){
				//		otherClass[i]=sampleData.getActualClassColorKey(1)[j][1];
				//		break;
				//	}
				//	otherClass[i]="Missing";
				//}
				//classCounts.add(code+"");//np
				//if (type == PlotPoint.MISSING || type == PlotPoint.NOT_A_NUMBER) callRate++;//zx
				otherClass[i] = determineCodeFromClass(3, alleleCounts[i], indi, chr, position)+"";
			} else {
				System.err.println("Error - no data for "+samples[i]);
			}
			
			// create grid
		}
		//callRate=(samples.length-callRate)*100/samples.length;//zx
		
		if (getUpdateQcPanel()) {
			sp.updateQcPanel(genotype, sex, otherClass);//zx
			setUpdateQcPanel(false);
		}
//		sp.updateColorKey(classCounts.convertToHash());
		
		Hashtable<String, String> hash = new Hashtable<String, String>();
		for (int i = 0; i < points.length; i++) {
			hash.put(points[i].getLayer()+"", "");
		}
		setLayersInBase(Array.toByteArray(Sort.putInOrder(Array.toIntArray(HashVec.getKeys(hash)))));
//    	rectangles = sp.getClusterFilterCollection().getRectangles(sp.getMarkerName(), sp.getCurrentClusterFilter(), (byte) plotType, (byte)1, false, false, (byte)0, (byte)99);
    	generateRectangles();
//    	System.out.println("Rectangles length: "+rectangles.length+" size: "+sp.getClusterFilterCollection().getSize(sp.getMarkerName())+" currentClusterFilter: "+sp.getCurrentClusterFilter());
//		if (sp.getClusterFilterCollection().getSize(sp.getMarkerName())>0) {
//			rectangles[sp.getCurrentClusterFilter()].setColor((byte)0);
//		}
		if (sp.getCurrentClusterFilter()>=0) {
			rectangles[sp.getCurrentClusterFilter()].setColor((byte)0);
		}
		sp.setCurrentClusterFilter(sp.getCurrentClusterFilter());
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
				index = Segment.binarySearch(new Segment((byte)-1, position, position), segs); 
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
		Hashtable<String,IndiPheno> sampleHash = sampleData.getSampleHash();
		float[] gcScores;
//		byte[] alleleCounts;
		float gcThreshold;
		int xWidth;
		int plotType, currentClass;
		int i;
		byte chr;
		int position;
		int markerIndex;
		byte size, xFontSize;
		
		//IntVector indeciesOfDataPoint;//zx

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
		//iv = locLookup.get(pos);
		//indeciesOfDataPoint = lookupNearbyPoints(x, y, pos);
		indeciesOfNearbySamples = lookupNearbyPoints(x, y, pos);
		//System.out.println("Number of nearby samples: "+(indeciesOfNearbySamples==null?0:indeciesOfNearbySamples.size()));//zx test point
		//prox = new IntVector();

		plotType = sp.getPlotType();
		currentClass = sp.getCurrentClass();
		markerIndex = sp.getMarkerIndex();
		datapoints = markerData[markerIndex].getDatapoints(plotType);
		gcScores = markerData[markerIndex].getGCs();
//		alleleCounts = markerData[markerIndex].getAB_Genotypes();
		chr = markerData[markerIndex].getChr();
		position = markerData[markerIndex].getPosition();
		gcThreshold = sp.getGCthreshold();

		size = sp.getPointSize();
		xFontSize = (byte)(size*2);

		g.setFont(new Font("Arial", 0, (int)(xFontSize*1.5)));
		xWidth = g.getFontMetrics(g.getFont()).stringWidth("X");

		//System.out.println("pos: "+pos+"\t iv.size():"+(indeciesOfNearbySamples==null?"null":indeciesOfNearbySamples.size()));//zx test point
		for (int l = 0; indeciesOfNearbySamples!=null&&l<indeciesOfNearbySamples.size(); l++) {
			i = indeciesOfNearbySamples.elementAt(l);
			indi = sampleHash.get(samples[i]);
			g.setColor(colorScheme[determineCodeFromClass(currentClass, alleleCounts[i], indi, chr, position)]);
			//g.setColor(Color.YELLOW);
//			if (gcScores[i]<gcThreshold) {
			if (alleleCounts[i]==-1) {
				g.drawString("X", getX(datapoints[0][i])-xWidth/2, getY(datapoints[1][i])+(int)(xFontSize/2.0));
			} else {
				g.fillOval(getX(datapoints[0][i])-(int)(size*2)/2, getY(datapoints[1][i])-(int)(size*2)/2, (int)(size*2), (int)(size*2));
			}
		}
		prevPos = pos;
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
		Hashtable<String,IndiPheno> sampleHash = sampleData.getSampleHash();
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

    public void mouseReleased(MouseEvent e) {
    	mouseEndX = e.getX();
    	mouseEndY = e.getY();
    	highlightRectangle = null;
    	
    	if (Math.abs(mouseEndX-mouseStartX)>(sp.getPointSize()/2) || Math.abs(mouseEndX-mouseStartX)>(sp.getPointSize()/2)) {
	    	// Automatically predict the new genotype and assigns to the last filter.
	    	sp.getClusterFilterCollection().addClusterFilter(sp.getMarkerName(),
	    											  new ClusterFilter((byte)sp.getPlotType(),
																		(float)Math.min(getRawX(mouseStartX), getRawX(mouseEndX)),
																		(float)Math.min(getRawY(mouseStartY), getRawY(mouseEndY)),
																		(float)Math.max(getRawX(mouseStartX), getRawX(mouseEndX)),
																		(float)Math.max(getRawY(mouseStartY), getRawY(mouseEndY)),
																		markerData[sp.getMarkerIndex()]));
	    	sp.saveClusterFilterCollection();
	    	sp.clusterFilterCollectionUpdated(true);
			setPointsGenerated(false);
			setUpdateQcPanel(true);
	    	generateRectangles();
			sp.setCurrentClusterFilter((byte)(sp.getClusterFilterCollection().getSize(sp.getMarkerName())-1));
			sp.displayClusterFilterIndex();
			paintAgain();
	    }

    	// Save the filters into files on hard drives;

    }

    public void mouseDragged(MouseEvent e) {
    	ClusterFilter clusterFilter;
    	mouseEndX = e.getX();
    	mouseEndY = e.getY();
    	highlightRectangle = new GenericRectangle((float)getRawX(mouseStartX), (float)getRawY(mouseStartY), (float)getRawX(mouseEndX), (float)getRawY(mouseEndY), (byte)1, false, false, (byte)0, (byte)99);

    	clusterFilter = new ClusterFilter((byte)sp.getPlotType(),
    			(float)Math.min(getRawX(mouseStartX), getRawX(mouseEndX)),
    			(float)Math.min(getRawY(mouseStartY), getRawY(mouseEndY)),
    			(float)Math.max(getRawX(mouseStartX), getRawX(mouseEndX)),
    			(float)Math.max(getRawY(mouseStartY), getRawY(mouseEndY)),
    			(byte)0);
    	highlightPoints(markerData[sp.getMarkerIndex()].getHighlightStatus(clusterFilter));
    	setExtraLayersVisible(new byte[] {99});
    	repaint();
    }

//	public void paintComponent(Graphics g) {
//		System.out.println("inner");
//		this.paintComponent(g);
////		super(g);
//	}
    
	public void mouseClicked(MouseEvent event) {
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

		
		window = Integer.parseInt(sp.getProject().getProperty(Project.NUM_MARKERS_PER_FILE));
		mData = markerData[sp.getMarkerIndex()];
		markerPosition = "chr"+mData.getChr()+":"+(mData.getPosition()-window)+"-"+(mData.getPosition()+window);
		if (indeciesOfNearbySamples!=null&&indeciesOfNearbySamples.size()>0) {
			menu = new JPopupMenu();
			for (int i = 0; i<indeciesOfNearbySamples.size(); i++) {
				// menu.add(samples[prox.elementAt(i)] +"
				// ("+datapoints[0][prox.elementAt(i)]+",
				// "+datapoints[1][prox.elementAt(i)]+")");

				menu.add(new LaunchAction(sp.getProject(), samples[indeciesOfNearbySamples.elementAt(i)], markerPosition, Color.BLACK));
			}
			menu.show(this, event.getX(), event.getY());
		}
	}

//	public void mouseEntered(MouseEvent e) {}
//
//	public void mouseExited(MouseEvent e) {}
//
//	public void mousePressed(MouseEvent e) {}
//
//	public void mouseReleased(MouseEvent e) {}

	public void componentHidden(ComponentEvent e) {}

	public void componentMoved(ComponentEvent e) {}

	public void componentResized(ComponentEvent e) {
		paintAgain();
	}

	public void componentShown(ComponentEvent e) {}

//	public static void main(String[] args) {
//		ScatterPlot.main(new String[] {"-notJar"});
//	}

	public void setUpdateQcPanel(boolean updateQcPanel) {
		this.updateQcPanel = updateQcPanel;
	}

	public boolean getUpdateQcPanel() {
		return updateQcPanel;
	}

	public void generateRectangles() {
    	rectangles = sp.getClusterFilterCollection().getRectangles(sp.getMarkerName(), (byte) sp.getPlotType(), (byte)1, false, false, (byte)7, (byte)99);
	}

	public GenericRectangle[] getRectangles() {
    	return rectangles;
	}
}
