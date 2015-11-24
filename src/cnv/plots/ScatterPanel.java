package cnv.plots;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import javax.swing.JPopupMenu;
import javax.swing.SwingUtilities;

import cnv.filesys.ClusterFilter;
import cnv.filesys.ClusterFilterCollection;
import cnv.filesys.MarkerData;
import cnv.filesys.Pedigree;
import cnv.gui.LaunchAction;
import cnv.qc.GcAdjustorParameter.GcAdjustorParameters;
import cnv.qc.MendelErrors.MendelErrorCheck;
import cnv.var.IndiPheno;
import cnv.var.SampleData;
import common.Array;
import common.CountVector;
import common.HashVec;
import common.IntVector;
import common.Logger;
import common.Sort;

// TODO Needs some cleanup, especially MouseMoved, MouseClicked, and generatePoints
public class ScatterPanel extends AbstractPanel implements MouseListener, MouseMotionListener {
	public static final long serialVersionUID = 3L;
	public static Color[] DEFAULT_COLORS = new Color[] {
			new Color(33, 31, 53), // dark dark
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
	
	
	byte[] alleleCounts;
	protected ScatterPlot sp;
	protected String[] samples;
	protected IntVector prox;
	protected SampleData sampleData;
	protected IntVector indicesOfNearbySamples;
	private boolean updateQcPanel;		//A control variable. Do not update QcPanel when resizing, or etc.
	private int mouseStartX;
	private int mouseStartY;
	private int mouseEndX;
	private int mouseEndY;
	private int panelIndex;
	CountVector uniqueValueCounts;
	protected boolean shrunk = false;
	
	public ScatterPanel(ScatterPlot sp, int index) {
		super();
		
		this.sp = sp;
		this.panelIndex = index;
		this.samples = sp.getSamples();
		this.sampleData = sp.getSampleData();
		this.updateQcPanel = true;
		
		setColorScheme(DEFAULT_COLORS);

		// taken care of in AbstractPanel constructor
//		addComponentListener(this);
		setZoomable(true, true);
	}
	
	public Color[] getColorScheme() {
		return colorScheme;
	}

	public void assignAxisLabels() {
		xAxisLabel = shrunk ? " " : ScatterPlot.TYPES[sp.getPlotType(panelIndex)][0];
		yAxisLabel = shrunk ? " " : ScatterPlot.TYPES[sp.getPlotType(panelIndex)][1];
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
		byte layer;
		IndiPheno indi;
		byte size, xFontSize;
		boolean[] displayCents;
		float[][][][] cents;
		int numCents, count;
		byte centSize;
		float x, y;
		int[] genotype;
		String[] sex;
		String[] otherClass;
		MarkerData markerData;
		int currentClass;
		Hashtable<String, String> disabledClassValues;
		boolean shiftColorOfSexChromosomes;
		String newGenotypingFilename;
		boolean[] isNewGenotypingDifferent = null;
		Logger log;
		
		log = sp.getProject().getLog();
		disabledClassValues = sp.getDisabledClassValues();//panelIndex);
		
		shiftColorOfSexChromosomes = sp.getProject().SHIFT_SEX_CHR_COLORS_YESNO.getValue();

		if (!sp.markerDataIsActive()) {
			return;
		}

		plotType = sp.getPlotType(panelIndex);
		currentClass = sp.getCurrentClass(panelIndex);
		if (currentClass < SampleData.BASIC_CLASSES.length && SampleData.BASIC_CLASSES[currentClass].equals(SampleData.HEATMAP)) {
	    	chartType = HEAT_MAP_TYPE;
		} else {
	    	chartType = SCATTER_PLOT_TYPE;
		}
		markerIndex = sp.getMarkerIndex();//this is the index from the markers loaded perspective
//		gcThreshold = sp.getGCthreshold();
		markerData = sp.getCurrentMarkerData();
		int markerProjectIndex = sp.getMarkerProjectIndices().get(markerData.getMarkerName()); //index of the marker in the project

		boolean[] toInclude = sp.hideExcludedSamples(panelIndex) ? sp.getProject().getSamplesToInclude(null, false) : Array.booleanArray(samples.length, true);
		out: if (plotType == 2 && sp.getDisplaygcAdjustor() != null && Array.booleanArraySum(sp.getDisplaygcAdjustor()) == 1) {// plot types should be changed to enums sometime
			for (int i = 0; i < sp.getDisplaygcAdjustor().length; i++) {
				if (sp.getDisplaygcAdjustor()[i]) {
					if (sp.getGcAdjustorParameters()[i] == null) {
						sp.getGcAdjustorParameters()[i] = GcAdjustorParameters.readSerial(sp.getGcCList()[1][i], log);
						log.reportTimeInfo("Lazy loading " + sp.getGcCList()[1][i]);
					}
					datapoints = markerData.getGCCorrectedLRR(sp.getGcAdjustorParameters()[i], markerProjectIndex, log);
					break out;
				}
			}
			datapoints = markerData.getDatapoints(plotType, null, toInclude, false, 1, sp.getGCthreshold(), sp.getClusterFilterCollection(), true, sp.getPcResids(), sp.getNumComponents(), 5, sp.getstdevFilter(), sp.getCorrectionRatio(), sp.getProject().getProperty(sp.getProject().NUM_THREADS), sp.getCorrection(panelIndex), sp.getProject().getLog());

		} else {

			datapoints = markerData.getDatapoints(plotType, null, toInclude, false, 1, sp.getGCthreshold(), sp.getClusterFilterCollection(), true, sp.getPcResids(), sp.getNumComponents(), 5, sp.getstdevFilter(), sp.getCorrectionRatio(), sp.getProject().getProperty(sp.getProject().NUM_THREADS), sp.getCorrection(panelIndex), sp.getProject().getLog());
		}
		
		// alleleCounts = markerData[markerIndex].getAB_Genotypes();
		//		alleleCounts = sp.getClusterFilterCollection().filterMarker(markerData[markerIndex], sp.getGCthreshold());
		alleleCounts = markerData.getAbGenotypesAfterFilters(sp.getClusterFilterCollection(), sp.getMarkerName(), sp.getGCthreshold(), log);
//		Project r = sp.getProject();
		newGenotypingFilename = sp.getProject().DATA_DIRECTORY.getValue(false, true) + sp.getMarkerName() + "_newGenotyping.xln";
		if(new File(newGenotypingFilename).exists()) {
			isNewGenotypingDifferent = loadNewGenotyping(newGenotypingFilename);
		}
//		sp.setCurrentClusterFilter(sp.getCurrentClusterFilter()); // what did this patch? this causes a continuous loop
		sp.displayClusterFilterIndex();
		chr = markerData.getChr();
		position = markerData.getPosition();
		size = sp.getPointSize();
		xFontSize = (byte)(size*2);
		displayCents = sp.getDisplayCentroids();
		cents = sp.getCentroids();
		centSize = 20;
		
		if (datapoints[0] == null || datapoints[1] == null) {
			errorMessage = "Data not available:";
			points = new PlotPoint[0]; 
			return;
		} else {
			errorMessage = null;
		}
		
		if (plotType == 0 || plotType == 1 || plotType >= 4) {
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
						points[count*3+j] = new PlotPoint("Centroids", PlotPoint.FILLED_CIRCLE, x, y, centSize, (byte)(5+i), (byte)10);

	                }
					count++;
				}
	        }
		} else {
			points = new PlotPoint[samples.length];
			numCents = 0;
		}
		
		
		if (plotType < 1 || plotType >= 4) {
			forcePlotXmax = Float.NaN;
		} else {
			forcePlotXmax = 1;
		}

		//dataForQc = new int[3][samples.length];
		//genotype = markerData[markerIndex].getAB_GenotypesAfterFilters(null, sp.getGCthreshold());
		genotype = new int[samples.length];
		sex = new String[samples.length];
		otherClass = new String[samples.length];
		uniqueValueCounts = new CountVector();
		
		for (int i = 0; i<samples.length; i++) {
			indi = sampleData.getIndiFromSampleHash(samples[i]);

			if (indi != null && (sp.hideExcludedSamples(panelIndex) && sampleData.individualShouldBeExcluded(samples[i]))) {
				// if sample should be excluded then do nothing
				genotype[i]=-3;
				sex[i] = "e";
				otherClass[i] = "e";
				
			} else if (indi != null) {
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
				} else if (sampleData.getClassName(currentClass).startsWith(SampleData.PLINK_CLASS_PREFIX)) {
				    byte indiCode = sp.getPlinkGenotypeForIndi(samples[i], currentClass);// chr1:159,937,288-159,945,728
					classCode = (byte)((int)indiCode + 1);
				} else {
					classCode = sampleData.determineCodeFromClass(currentClass, alleleCounts[i], indi, chr, position);
				}
				
				if (classCode <= -1 && !sp.maskMissing(panelIndex)) {
					classCode = 0;
				}
				if (Float.isNaN(datapoints[0][i]) || Float.isNaN(datapoints[1][i])) {
					type = PlotPoint.NOT_A_NUMBER;
//				} else if (currentClass==1 && alleleCounts[i]==-1) {
				} else if (sp.getGCthreshold() > 0 && alleleCounts[i]==-1) {
					type = PlotPoint.MISSING;
				} else if (isNewGenotypingDifferent != null && isNewGenotypingDifferent[i]) {
					type = PlotPoint.OPEN_SQUARE;
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
					genotype[i]=0;
				} else {
					uniqueValueCounts.add(classCode+"");
				}
				if (classCode < 0) {
					log.reportError("Error - classCode is less than 0 ("+classCode+")");
				}
				if (currentClass < SampleData.BASIC_CLASSES.length && SampleData.BASIC_CLASSES[currentClass].equals(SampleData.GENOTYPE) && chr > 22 && shiftColorOfSexChromosomes) {
					points[numCents*3+i] = new PlotPoint(samples[i], type, datapoints[0][i], datapoints[1][i], type==PlotPoint.FILLED_CIRCLE?size:(type==PlotPoint.FILLED_TRIANGLE?(byte)(size+5):xFontSize), classCode==0? 0 : (byte) (classCode +  3), layer);
				} else {
					points[numCents*3+i] = new PlotPoint(samples[i], type, datapoints[0][i], datapoints[1][i], type==PlotPoint.FILLED_CIRCLE?size:(type==PlotPoint.FILLED_TRIANGLE?(byte)(size+5):xFontSize), classCode, layer);
				}
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
				log.reportError("Error - no data pts for "+samples[i]);
				sex[i] = "missing";
				points[numCents*3+i] = new PlotPoint(samples[i], PlotPoint.MISSING, datapoints[0][i], datapoints[1][i], (byte)(xFontSize*2), (byte)0, (byte)99);
			}
			
			// create grid
		}
		//callRate=(samples.length-callRate)*100/samples.length;
		if (getUpdateQcPanel()) {
			sp.updateQcPanel(chr, genotype, sex, otherClass, panelIndex);
			setUpdateQCPanel(false);
		}
		sp.updateColorKey(uniqueValueCounts.convertToHash(), panelIndex);
		
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
		if (sp.displayMendelianErrors(panelIndex)) {
		    generateLines(sex);
		} else {
		    lines = new GenericLine[0];
		}
//		sp.setCurrentClusterFilter(sp.getCurrentClusterFilter()); // what did this patch? this causes a continuous loop
	}

	private void generateLines(String[] sex) {
		ArrayList<GenericLine> linesList = new ArrayList<GenericLine>();
		int plotType = sp.getPlotType(panelIndex);
		int centroidOffset = (plotType == 0 || plotType == 1 || plotType >= 4 ? Array.booleanArraySum(sp.getDisplayCentroids()) : 0);
		centroidOffset *= 3;
		byte size = (byte) 3;
		byte momColor = (byte) 6;
		byte dadColor = (byte) 7;
		byte layer = (byte) 1;
		boolean swapAxes = false;
		
		if (sp.getPedigree() != null) {
			MendelErrorCheck[] mendelErrorChecks = Pedigree.PedigreeUtils.checkMendelErrors(sp.getPedigree(), sp.getCurrentMarkerData(), sp.hideExcludedSamples(panelIndex) ? sp.getProject().getSamplesToInclude(null,false) : null, sex, sp.getClusterFilterCollection(), sp.getGCthreshold(), sp.getProject().getLog());
			if (mendelErrorChecks != null) {
    			for (int i = 0; i < samples.length; i++) {
    				PlotPoint indiPoint = points[centroidOffset + i];
    				if (mendelErrorChecks[i].hasMoMendelError()) {
    					PlotPoint momPoint = points[centroidOffset + sp.getPedigree().getMoDNAIndex(i)];
    					GenericLine gl = new GenericLine(momPoint, indiPoint, size, momColor, layer, swapAxes, 1);
    					linesList.add(gl);
    				}
    				if (mendelErrorChecks[i].hasFaMendelError()) {
    					PlotPoint dadPoint = points[centroidOffset + sp.getPedigree().getFaDNAIndex(i)];
    					GenericLine gl = new GenericLine(dadPoint, indiPoint, size, dadColor, layer, swapAxes, 1);
    					linesList.add(gl);
    				}
    
    		    }
			}
		}
	    
	    lines = linesList.toArray(new GenericLine[linesList.size()]);
    }

    private boolean[] loadNewGenotyping(String alternativeGenotypingFilename) {
		Hashtable<String, Byte> alleleCountsNew;
		boolean[] isNewGenotypingDifferent;
		BufferedReader reader;
		String[] line;
		
		try {
			alleleCountsNew = new Hashtable<String, Byte> (alleleCounts.length);
			isNewGenotypingDifferent = new boolean[alleleCounts.length];
			reader = new BufferedReader(new FileReader(alternativeGenotypingFilename));
			reader.readLine();
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				alleleCountsNew.put(line[0], Byte.parseByte(line[9]));
			}
			reader.close();
			
			for (int i = 0; i < samples.length; i++) {
				if (alleleCountsNew.containsKey(samples[i]) && alleleCounts[i] != alleleCountsNew.get(samples[i])) {
					alleleCounts[i] = alleleCountsNew.get(samples[i]);
					isNewGenotypingDifferent[i] = true;
				}
			}
			
		} catch (FileNotFoundException fnfe) {
//			log.reportError("Error: file \"" + name + "\" not found in current directory");
			return null;
		} catch (IOException ioe) {
//			log.reportError("Error reading file \"" + name + "\"");
			return null;
		}
		return isNewGenotypingDifferent;
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
		
		//IntVector indicesOfDataPoint;
		
		
		plotType = sp.getPlotType(panelIndex);
		currentClass = sp.getCurrentClass(panelIndex);
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
		//indicesOfDataPoint = lookupNearbyPoints(x, y, pos);
		indicesOfNearbySamples = lookupNearbyPoints(x, y, pos);
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

		for (int l = 0; indicesOfNearbySamples!=null && l<indicesOfNearbySamples.size(); l++) {
			i = indicesOfNearbySamples.elementAt(l);
			if (i < samples.length) { // can also be centroids or other points
				indi = sampleData.getIndiFromSampleHash(samples[i]);
				byte classCode;
				if (sampleData.getClassName(currentClass).startsWith(SampleData.PLINK_CLASS_PREFIX))  {
				    classCode = (byte)((int)sp.getPlinkGenotypeForIndi(samples[i], currentClass) + 1);
				} else {
				    classCode = sampleData.determineCodeFromClass(currentClass, alleleCounts[i], indi, chr, position);
				}
				g.setColor(colorScheme[Math.max(0, Math.min(classCode, colorScheme.length-1))]);
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
        if (SwingUtilities.isLeftMouseButton(e) && !e.isControlDown()) {
        	mouseStartX = e.getX();
        	mouseStartY = e.getY();
        } else {
            super.mousePressed(e);
        }
    }

    public void mouseReleased(MouseEvent e) {
    	mouseEndX = e.getX();
    	mouseEndY = e.getY();
    	highlightRectangle = null;
    	if (SwingUtilities.isLeftMouseButton(e) && !e.isControlDown()) {
        	if (Math.abs(mouseEndX - mouseStartX) > (sp.getPointSize() / 2) || Math.abs(mouseEndX - mouseStartX) > (sp.getPointSize() / 2)) {
    	    	// Automatically predict the new genotype and assigns to the last filter.
    	    	sp.getClusterFilterCollection().addClusterFilter(sp.getMarkerName(),
    	    											  new ClusterFilter((byte)sp.getPlotType(panelIndex),
    																		(float)Math.max(plotXmin, Math.min(getXValueFromXPixel(mouseStartX), getXValueFromXPixel(mouseEndX))),
    																		(float)Math.max(plotYmin, Math.min(getYValueFromYPixel(mouseStartY), getYValueFromYPixel(mouseEndY))),
    																		(float)Math.min(plotXmax, Math.max(getXValueFromXPixel(mouseStartX), getXValueFromXPixel(mouseEndX))),
    																		(float)Math.min(plotYmax, Math.max(getYValueFromYPixel(mouseStartY), getYValueFromYPixel(mouseEndY))),
    																		sp.getCurrentMarkerData()));
    //	    	sp.startAutoSaveToTempFile();
    	    	sp.setClusterFilterUpdated(true);
    			setPointsGeneratable(true);
    			setUpdateQCPanel(true);
    	    	generateRectangles();
    			sp.setCurrentClusterFilter((byte)(sp.getClusterFilterCollection().getSize(sp.getMarkerName())-1));
    			sp.displayClusterFilterIndex();
    			paintAgain();
    	    }
    	} else {
    	    super.mouseReleased(e);
    	}

    	//TODO Save the filters into files on hard drives;
    }

    public void mouseDragged(MouseEvent e) {
        if (SwingUtilities.isLeftMouseButton(e) && !e.isControlDown()) {
        	ClusterFilter clusterFilter;
        	mouseEndX = e.getX();
        	mouseEndY = e.getY();
        	highlightRectangle = new GenericRectangle((float)getXValueFromXPixel(mouseStartX), (float)getYValueFromYPixel(mouseStartY), (float)getXValueFromXPixel(mouseEndX), (float)getYValueFromYPixel(mouseEndY), (byte)1, false, false, (byte)0, (byte)99);
        	
        	clusterFilter = new ClusterFilter((byte)sp.getPlotType(panelIndex),
        			(float)Math.min(getXValueFromXPixel(mouseStartX), getXValueFromXPixel(mouseEndX)),
        			(float)Math.min(getYValueFromYPixel(mouseStartY), getYValueFromYPixel(mouseEndY)),
        			(float)Math.max(getXValueFromXPixel(mouseStartX), getXValueFromXPixel(mouseEndX)),
        			(float)Math.max(getYValueFromYPixel(mouseStartY), getYValueFromYPixel(mouseEndY)),
        			(byte)0);
        	highlightPoints(sp.getCurrentMarkerData().getHighlightStatus(clusterFilter));
        	setExtraLayersVisible(new byte[] {99});
            repaint();
        } else {
            super.mouseDragged(e);
        }
    }
    
	public void mouseClicked(MouseEvent event) {
		JPopupMenu menu;
		MarkerData mData;
		String markerPosition;
		int window, position;
		int numberToInclude, currentClass;
		byte newClusterFilter, chr;

		if (SwingUtilities.isRightMouseButton(event)) {
			window = sp.getProject().WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER.getValue();
			mData = sp.getCurrentMarkerData();
			markerPosition = "chr"+mData.getChr()+":"+(mData.getPosition()-window)+"-"+(mData.getPosition()+window);
	        currentClass = sp.getCurrentClass(panelIndex);
	        chr = mData.getChr();
	        position = mData.getPosition();
			if (indicesOfNearbySamples!=null && indicesOfNearbySamples.size()>0) {
				menu = new JPopupMenu();
				numberToInclude = Math.min(50, indicesOfNearbySamples.size());
				for (int i = 0; i<numberToInclude; i++) {
				    int sampleIndex = indicesOfNearbySamples.elementAt(i);
				    IndiPheno indi = sampleData.getIndiFromSampleHash(samples[sampleIndex]);
	                byte classCode;
	                if (sampleData.getClassName(currentClass).startsWith(SampleData.PLINK_CLASS_PREFIX))  {
	                    classCode = (byte)((int)sp.getPlinkGenotypeForIndi(samples[sampleIndex], currentClass) + 1);
	                } else {
	                    classCode = sampleData.determineCodeFromClass(currentClass, alleleCounts[sampleIndex], indi, chr, position);
	                }
					menu.add(new LaunchAction(sp.getProject(), samples[sampleIndex], markerPosition, colorScheme[Math.max(0, Math.min(classCode, colorScheme.length-1))]));
				}
				if (indicesOfNearbySamples.size() > 50) {
					menu.add(new LaunchAction("Plus "+(indicesOfNearbySamples.size() - 50)+" additional samples"));
				}
	
				menu.show(this, event.getX(), event.getY());
			}

		} else if (SwingUtilities.isLeftMouseButton(event)) {
			newClusterFilter = lookupNearbyRectangles(event.getX(), event.getY());
			if (newClusterFilter >= 0) {
				ClusterFilterCollection clusterFilterCollection;
				clusterFilterCollection = sp.getClusterFilterCollection();
				clusterFilterCollection.deleteClusterFilter(sp.getMarkerName(), newClusterFilter);
				sp.setCurrentClusterFilter((byte) Math.min(newClusterFilter, clusterFilterCollection.getSize(sp.getMarkerName())-1));
		    	sp.setClusterFilterUpdated(true);
				sp.displayClusterFilterIndex();
				setPointsGeneratable(true);
				setUpdateQCPanel(true);
				generateRectangles();
				sp.updateGUI();
			}

		}
	}

	public void setUpdateQCPanel(boolean updateQcPanel) {
		this.updateQcPanel = updateQcPanel;
	}

	public boolean getUpdateQcPanel() {
		return updateQcPanel;
	}

	public void generateRectangles() {
	    rectangles = sp.getClusterFilterCollection().getRectangles(sp.getMarkerName(), (byte) sp.getPlotType(panelIndex), (byte)1, false, false, (byte)7, (byte)99);
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
