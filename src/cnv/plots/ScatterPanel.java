/*
 * Revised in Feb, 2012. Add mouse dragged feature.
 */
package cnv.plots;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.Hashtable;

import javax.swing.JPopupMenu;

import cnv.filesys.ClusterFilter;
import cnv.filesys.ClusterFilterCollection;
import cnv.filesys.MarkerData;
import cnv.filesys.Project;
import cnv.gui.LaunchAction;
import cnv.var.SampleData;
import cnv.var.IndiPheno;
import common.CountVector;
import common.HashVec;
import common.IntVector;
import common.Sort;
import common.Array;

// TODO Needs some cleanup, especially MouseMoved, MouseClicked, and generatePoints
public class ScatterPanel extends AbstractPanel implements MouseListener, MouseMotionListener {
	public static final long serialVersionUID = 3L;
	public static Color[] generateDefaultColors() { // new after merging with npBranch, feel free to delete if better solution is found
		return	new Color[] {new Color(33, 31, 53), // dark dark
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
			   	new Color(255, 192, 0),	// yellowy orange
			   	new Color(227, 108, 9), // halloween orange
			};
	}
	public static String[] comments = new String[] {"Ugly", "Monomorphic", "Off Y-axis"};
	
	protected MarkerData[] markerData1;
	byte[] alleleCounts;				//zx
	protected ScatterPlot sp;
	protected String[] samples;
	protected IntVector prox;
	protected SampleData sampleData;
	protected Color[] colorScheme; // new after merging with npBranch, feel free to delete if better solution is found
	protected IntVector indicesOfNearbySamples;	//zx
//	protected IntVector indeciesOfNaNSamples;	//zx
	private boolean updateQcPanel;		//zx: A control variable. Do not update QcPanel when resizing, or etc.
	private int mouseStartX;
	private int mouseStartY;
	private int mouseEndX;
	private int mouseEndY;
	
	public ScatterPanel(ScatterPlot sp) {
		super();
		
		this.sp = sp;
		this.samples = sp.getSamples();
//		this.markerData = sp.getMarkerData();
		this.sampleData = sp.getSampleData();
//		locLookup = new Hashtable<String,IntVector>();//zx
		this.updateQcPanel = true;//zx
//		this.indeciesOfNaNSamples = new IntVector();
		
//		setColorScheme(ColorKeyPanel.DEFAULT_COLORS);
		colorScheme = generateDefaultColors();   // new after merging with npBranch, feel free to delete if better solution is found
		setColorScheme(colorScheme);
		

//		addMouseListener(this);
//		addMouseMotionListener(this);
		addComponentListener(this);
		setZoomable(true, true);
	}

	public Color[] getColorScheme() {
		return colorScheme;
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
			if (points[i] != null) {
				if (points[i].isHighlighted()) {
					points[i].setSize((byte)(defaultSize*1.5));
				} else {
					points[i].setSize((byte)(defaultSize));
				}
			}
		}
	}

	public void generatePoints() {
		int position, markerIndex, plotType;
		byte chr, genotypeCode, classCode, type;
		float[][] datapoints;
//		byte[] alleleCounts;
//		float gcThreshold;
		byte layer;
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
		CountVector uniqueValueCounts;//zx
//		ClusterFilterCollection clusterFilterCollection;//zx
		MarkerData markerData;
		int currentClass;
		Hashtable<String, String> disabledClassValues;
		
		disabledClassValues = sp.getDisabledClassValues();

//		time = new Date().getTime();
		
		if (!sp.markerDataIsActive()) {
			return;
		}

		plotType = sp.getPlotType();
		currentClass = sp.getCurrentClass();
		if (currentClass < SampleData.BASIC_CLASSES.length && SampleData.BASIC_CLASSES[currentClass].equals(SampleData.HEATMAP)) {
	    	chartType = HEAT_MAP_TYPE;
		} else {
	    	chartType = SCATTER_PLOT_TYPE;
		}
		markerIndex = sp.getMarkerIndex();
//		gcThreshold = sp.getGCthreshold();
		markerData = sp.getCurrentMarkerData();
		datapoints = markerData.getDatapoints(plotType);
//		alleleCounts = markerData[markerIndex].getAB_Genotypes();//zx
//		alleleCounts = sp.getClusterFilterCollection().filterMarker(markerData[markerIndex], sp.getGCthreshold());
		alleleCounts = markerData.getAbGenotypesAfterFilters(sp.getClusterFilterCollection(), sp.getMarkerName(), sp.getGCthreshold());//zx
//		sp.setCurrentClusterFilter(sp.getCurrentClusterFilter()); // what did this patch? this causes a continuous loop
		sp.displayClusterFilterIndex();
		chr = markerData.getChr();
		position = markerData.getPosition();
		size = sp.getPointSize();
		xFontSize = (byte)(size*2);
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
		uniqueValueCounts = new CountVector();
		for (int i = 0; i<samples.length; i++) {
			indi = sampleData.getIndiFromSampleHash(samples[i]);

			if (indi != null && sampleData.individualShouldBeExcluded(samples[i])) {
				// if sample should be excluded then do nothing
				genotype[i]=-3;
				sex[i] = "e";
				otherClass[i] = "e";
				
			} else if (indi!=null) {
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
					classCode = sampleData.determineCodeFromClass(currentClass, alleleCounts[i], indi, chr, position);
				}
				if (sampleData.getSexClassIndex() == -1) {
				} else {
				}
				
				//System.out.println("Current loop:\t"+i+"\t code:\t"+code+"\t alleleCounts: "+alleleCounts[i]+"\t gcScores: "+gcScores[i]);//zx
				if (classCode <= -1 && !sp.maskMissing()) {
					classCode = 0;
				}
				if (Float.isNaN(datapoints[0][i]) || Float.isNaN(datapoints[1][i])) {
					type = PlotPoint.NOT_A_NUMBER;
//				} else if (currentClass==1 && alleleCounts[i]==-1) {
				} else if (sp.getGCthreshold() > 0 && alleleCounts[i]==-1) {
					type = PlotPoint.MISSING;
				} else {
					type = PlotPoint.FILLED_CIRCLE;
				}
				if (classCode == 100) {
					classCode = genotypeCode;
					if (classCode == 0) {
						type = PlotPoint.FILLED_TRIANGLE;
						classCode = (byte)(colorScheme.length-1);
					}
				}
				
				layer = (byte)((sampleData.getClassCategoryAndIndex(currentClass)[0]==2 && classCode > 0)?1:0);
				layer = classCode; // TODO temporary fix, since was always zero otherwise
				
				if (type == PlotPoint.NOT_A_NUMBER || type == PlotPoint.MISSING) {
					uniqueValueCounts.add(0+"");
					genotype[i]=0;//zx
				} else {
					uniqueValueCounts.add(classCode+"");
				}
				if (classCode < 0) {
					System.err.println("Error - classCode is less than 0 ("+classCode+")");
				}
				points[numCents*3+i] = new PlotPoint(samples[i], type, datapoints[0][i], datapoints[1][i], type==PlotPoint.FILLED_CIRCLE?size:(type==PlotPoint.FILLED_TRIANGLE?(byte)(size+5):xFontSize), classCode, layer);
				genotype[i]=genotypeCode;
				//sex[i]=(sexCode==1?"Female":(sexCode==2?"Male":"Missing"));
				//for (int j=0; j<sampleData.getActualClassColorKey(0).length; j++) {
				//	if (sampleData.getActualClassColorKey(0)[j][0].equals(determineCodeFromClass(2, alleleCounts[i], indi, chr, position)+"")){
				//		sex[i]=sampleData.getActualClassColorKey(0)[j][1];
				//		break;
				//	}
				//	sex[i]="Missing";
				//}
				sex[i] = sampleData.determineCodeFromClass(2, alleleCounts[i], indi, chr, position)+"";
				
				//for (int j=0; j<sampleData.getActualClassColorKey(1).length; j++) {
				//	if (sampleData.getActualClassColorKey(1)[j][0].equals(classCode+"")){
				//		otherClass[i]=sampleData.getActualClassColorKey(1)[j][1];
				//		break;
				//	}
				//	otherClass[i]="Missing";
				//}
				//classCounts.add(code+"");
				//if (type == PlotPoint.MISSING || type == PlotPoint.NOT_A_NUMBER) callRate++;
				otherClass[i] = sampleData.determineCodeFromClass(currentClass, alleleCounts[i], indi, chr, position) + "";
			} else {
				System.err.println("Error - no data pts for "+samples[i]);
				sex[i] = "missing";
				points[numCents*3+i] = new PlotPoint(samples[i], PlotPoint.MISSING, datapoints[0][i], datapoints[1][i], (byte)(xFontSize*2), (byte)0, (byte)99);
			}
			
			// create grid
		}
		//callRate=(samples.length-callRate)*100/samples.length;
		
		if (getUpdateQcPanel()) {
			sp.updateQcPanel(chr, genotype, sex, otherClass);
			setQcPanelUpdatable(false);
		}
		sp.updateColorKey(uniqueValueCounts.convertToHash());
		
		Hashtable<String, String> hash = new Hashtable<String, String>();
		for (int i = 0; i < points.length; i++) {
			if (points[i] != null && disabledClassValues.containsKey(currentClass+"\t"+points[i].getColor())) {
				points[i].setVisible(false);
			}
			if (points[i] != null) {
				hash.put(points[i].getLayer()+"", "");
			}
		}
		setLayersInBase(Array.toByteArray(Sort.putInOrder(Array.toIntArray(HashVec.getKeys(hash)))));
    	generateRectangles();
    	setSwapable(false);
		if (sp.getCurrentClusterFilter()>=0) {
			rectangles[sp.getCurrentClusterFilter()].setColor((byte)0);
		}
//		sp.setCurrentClusterFilter(sp.getCurrentClusterFilter()); // what did this patch? this causes a continuous loop
	}

	public void mouseMoved(MouseEvent event) {
		Graphics g = getGraphics();
		String pos;
		int x, y;

		float[][] datapoints;
		IndiPheno indi;
//		float[] gcScores;
//		byte[] alleleCounts;
//		float gcThreshold;
		int xWidth;
		int plotType, currentClass;
		int i;
		byte chr;
		int position;
//		int markerIndex;
		byte size, xFontSize;
		MarkerData mData;
		
		//IntVector indeciesOfDataPoint;//zx
		
		
		plotType = sp.getPlotType();
		currentClass = sp.getCurrentClass();
//		markerIndex = sp.getMarkerIndex();

//		if (markerData == null || markerData[markerIndex] == null) {
		if (!sp.markerDataIsActive()) {
			return;
		}

		x = event.getX();
		y = event.getY();

//		canvasSectionMinimumX = WIDTH_Y_AXIS;
//		canvasSectionMaximumX = getWidth()-WIDTH_BUFFER;
//		canvasSectionMinimumY = HEIGHT_X_AXIS;
//		canvasSectionMaximumY = getHeight()-HEAD_BUFFER;
		pos = (int) Math.floor(x / DEFAULT_LOOKUP_RESOLUTION) + "x" + (int) Math.floor(y / DEFAULT_LOOKUP_RESOLUTION);
		if (!pos.equals(prevPos)) {
			repaint();
		}
		//iv = locLookup.get(pos);
		//indeciesOfDataPoint = lookupNearbyPoints(x, y, pos);
		indicesOfNearbySamples = lookupNearbyPoints(x, y, pos);
		//System.out.println("Number of nearby samples: "+(indeciesOfNearbySamples==null?0:indeciesOfNearbySamples.size()));//zx test point
		//prox = new IntVector();

		mData = sp.getCurrentMarkerData();
		datapoints = mData.getDatapoints(plotType);
//		gcScores = mData.getGCs();
//		alleleCounts = markerData[markerIndex].getAB_Genotypes();
		chr = mData.getChr();
		position = mData.getPosition();
//		gcThreshold = sp.getGCthreshold();

		size = sp.getPointSize();
		xFontSize = (byte)(size*2);

		g.setFont(new Font("Arial", 0, (int)(xFontSize*1.5)));
		xWidth = g.getFontMetrics(g.getFont()).stringWidth("X");

		//System.out.println("pos: "+pos+"\t iv.size():"+(indeciesOfNearbySamples==null?"null":indeciesOfNearbySamples.size()));//zx test point
		for (int l = 0; indicesOfNearbySamples!=null && l<indicesOfNearbySamples.size(); l++) {
			i = indicesOfNearbySamples.elementAt(l);
			if (i < samples.length) { // can also be centroids or other points
				indi = sampleData.getIndiFromSampleHash(samples[i]);
				g.setColor(colorScheme[Math.min(sampleData.determineCodeFromClass(currentClass, alleleCounts[i], indi, chr, position), colorScheme.length-1)]);
				//g.setColor(Color.YELLOW);
	//			if (gcScores[i]<gcThreshold) {
	//			if (currentClass==1 && alleleCounts[i]==-1) {
				if (sp.getGCthreshold() > 0 && alleleCounts[i]==-1) {
					g.drawString("X", getXPixel(datapoints[0][i])-xWidth/2, getYPixel(datapoints[1][i])+(int)(xFontSize/2.0));
				} else {
					g.fillOval(getXPixel(datapoints[0][i])-(int)(size*2)/2, getYPixel(datapoints[1][i])-(int)(size*2)/2, (int)(size*2), (int)(size*2));
				}
			}
		}
		prevPos = pos;
	}

    public void mousePressed(MouseEvent e) {
    	mouseStartX = e.getX();
    	mouseStartY = e.getY();
    }

    public void mouseReleased(MouseEvent e) {
    	mouseEndX = e.getX();
    	mouseEndY = e.getY();
    	highlightRectangle = null;
    	
    	if (Math.abs(mouseEndX - mouseStartX) > (sp.getPointSize() / 2) || Math.abs(mouseEndX - mouseStartX) > (sp.getPointSize() / 2)) {
	    	// Automatically predict the new genotype and assigns to the last filter.
	    	sp.getClusterFilterCollection().addClusterFilter(sp.getMarkerName(),
	    											  new ClusterFilter((byte)sp.getPlotType(),
																		(float)Math.min(getXValueFromXPixel(mouseStartX), getXValueFromXPixel(mouseEndX)),
																		(float)Math.min(getYValueFromYPixel(mouseStartY), getYValueFromYPixel(mouseEndY)),
																		(float)Math.max(getXValueFromXPixel(mouseStartX), getXValueFromXPixel(mouseEndX)),
																		(float)Math.max(getYValueFromYPixel(mouseStartY), getYValueFromYPixel(mouseEndY)),
																		sp.getCurrentMarkerData()));
//	    	sp.startAutoSaveToTempFile();
	    	sp.updateClusterFilterStatus(true);
			setPointsGeneratable(true);
			setQcPanelUpdatable(true);
	    	generateRectangles();
			sp.setCurrentClusterFilter((byte)(sp.getClusterFilterCollection().getSize(sp.getMarkerName())-1));
			sp.displayClusterFilterIndex();
			paintAgain();
	    }

    	//TODO Save the filters into files on hard drives;
    }

    public void mouseDragged(MouseEvent e) {
    	ClusterFilter clusterFilter;
    	mouseEndX = e.getX();
    	mouseEndY = e.getY();
    	highlightRectangle = new GenericRectangle((float)getXValueFromXPixel(mouseStartX), (float)getYValueFromYPixel(mouseStartY), (float)getXValueFromXPixel(mouseEndX), (float)getYValueFromYPixel(mouseEndY), (byte)1, false, false, (byte)0, (byte)99);

    	clusterFilter = new ClusterFilter((byte)sp.getPlotType(),
    			(float)Math.min(getXValueFromXPixel(mouseStartX), getXValueFromXPixel(mouseEndX)),
    			(float)Math.min(getYValueFromYPixel(mouseStartY), getYValueFromYPixel(mouseEndY)),
    			(float)Math.max(getXValueFromXPixel(mouseStartX), getXValueFromXPixel(mouseEndX)),
    			(float)Math.max(getYValueFromYPixel(mouseStartY), getYValueFromYPixel(mouseEndY)),
    			(byte)0);
    	highlightPoints(sp.getCurrentMarkerData().getHighlightStatus(clusterFilter));
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
		int numberToInclude;
		byte newClusterFilter;
//		byte currentGenotype;
//		ClusterFilter currentClusterFilter;

		if (event.getButton() == MouseEvent.BUTTON1) {
			window = Integer.parseInt(sp.getProject().getProperty(Project.WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER));
			mData = sp.getCurrentMarkerData();
			markerPosition = "chr"+mData.getChr()+":"+(mData.getPosition()-window)+"-"+(mData.getPosition()+window);
			if (indicesOfNearbySamples!=null && indicesOfNearbySamples.size()>0) {
				menu = new JPopupMenu();
				numberToInclude = Math.min(50, indicesOfNearbySamples.size());
				for (int i = 0; i<numberToInclude; i++) {
					// menu.add(samples[prox.elementAt(i)] +"
					// ("+datapoints[0][prox.elementAt(i)]+",
					// "+datapoints[1][prox.elementAt(i)]+")");
	
					menu.add(new LaunchAction(sp.getProject(), samples[indicesOfNearbySamples.elementAt(i)], markerPosition, Color.BLACK));
				}
				if (indicesOfNearbySamples.size() > 50) {
					menu.add(new LaunchAction("Plus "+(indicesOfNearbySamples.size() - 50)+" additional samples"));
				}
	
				menu.show(this, event.getX(), event.getY());
			}

			newClusterFilter = lookupNearbyRectangles(event.getX(), event.getY());
			if (newClusterFilter >= 0) {
				sp.setCurrentClusterFilter(newClusterFilter);
			}

		} else if (event.getButton() == MouseEvent.BUTTON2) {
			newClusterFilter = lookupNearbyRectangles(event.getX(), event.getY());
			if (newClusterFilter >= 0) {
				ClusterFilterCollection clusterFilterCollection;
				clusterFilterCollection = sp.getClusterFilterCollection();
				clusterFilterCollection.deleteClusterFilter(sp.getMarkerName(), newClusterFilter);
				sp.setCurrentClusterFilter((byte) Math.min(newClusterFilter, clusterFilterCollection.getSize(sp.getMarkerName())-1));
//				sp.startAutoSaveToTempFile();
				sp.updateClusterFilterStatus(true);
				//clusterFilterNavigation.setText((clusterFilterCollection.getSize(getMarkerName())==0?0:(currentClusterFilter+1))+" of "+clusterFilterCollection.getSize(getMarkerName()));
				sp.displayClusterFilterIndex();
//					annotationPanelLowerPart();
				setPointsGeneratable(true);
				setQcPanelUpdatable(true);
				//repaint();
				//paintAgain();
				generateRectangles();
				sp.updateGUI();
			}

		} else if (event.getButton() == MouseEvent.BUTTON3) {
			newClusterFilter = lookupResidingRectangles(event.getX(), event.getY());
			if (newClusterFilter >= 0) {
				sp.setCurrentClusterFilter(newClusterFilter);
				sp.updateCurrentClusterFiltersGenotype();
			}
		}
	}

//	public static void main(String[] args) {
//		ScatterPlot.main(new String[] {"-notJar"});
//	}

	public void setQcPanelUpdatable(boolean updateQcPanel) {
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

//	public void setCurrentClass (byte newCurrentClass) {
//		currentClass = newCurrentClass;
//	}
//


//	public boolean isCNV () {
//		DoubleVector x = new DoubleVector();
//		DoubleVector y = new DoubleVector();
//		float[][] datapoints;
//		
//		if (sp.getPlotType()==1) {
//			datapoints = markerData[sp.getMarkerIndex()].getDatapoints(sp.getPlotType());
//			for (int i=0; i<alleleCounts.length; i++) {
//				if (alleleCounts[i]==0 || alleleCounts[i]==1) {
//					x.add(datapoints[0][i]);
//				} else if (alleleCounts[i]==2 || alleleCounts[i]==1) {
//					y.add(datapoints[1][i]);
//				}
//			}
//			if (Array.isBimodal(x.toArray()) || Array.isBimodal(y.toArray())) {
//				return true;
//			} else {
//				return false;
//			}
//		}
//		return (Boolean) null;
//	}

}
