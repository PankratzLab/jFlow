package cnv.plots;

import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;

import javax.swing.*;
import javax.swing.border.EmptyBorder;

import common.*;
import cnv.filesys.*;
import cnv.gui.SingleClick;
import cnv.gui.ClickListener;
import cnv.manage.Transforms;
import cnv.qc.GcAdjustor;
import cnv.qc.GcAdjustor.GcModel;
import cnv.var.CNVariant;
import cnv.var.IndiPheno;
import cnv.var.SampleData;
import filesys.*;

public class Trailer extends JFrame implements ActionListener, ClickListener, MouseListener, MouseMotionListener, MouseWheelListener {
	public static final long serialVersionUID = 1L;

//	public static final String DEFAULT_LOCATION = "chr6:161,739,001-163,119,245"; // parkin
//	public static final String DEFAULT_LOCATION = "chr5:151497267-151499003"; // hit
//	public static final String DEFAULT_LOCATION = "chr5:151,175,382-151,203,885"; // 1 gene
	public static final String DEFAULT_LOCATION = "chr17:55,609,472-55,824,368"; // USP32
//	public static final String DEFAULT_LOCATION = "chr12"; // USP32

	public static final String DEFAULT_SAMPLE = null;
//	public static final String DEFAULT_SAMPLE = "ND12014";
	
	// public static final String DEFAULT_LOCATION = "chr10:96,040,187-97,106,178";
	// public static final String DEFAULT_LOCATION = "chr23";

	// public static final String[] DEFAULT_CNV_FILES = {"data/penncnv.cnv", "data/quantisnp_windows.cnv", "data/quantisnp_linux.cnv", "data/dnacopy10.cnv"};
	// public static final String[] DEFAULT_CNV_FILES = {"data/conf.cnv", "data/allMarkers.cnv"};
	// public static final String[] DEFAULT_CNV_FILES = {"data/conf_100kb_5SNP_10.0.cnv", "data/allMarkers_100kb_5SNP_10.0.cnv"};
	// public static final String[] DEFAULT_CNV_FILES = {"data/conf_100bp_5SNP_10.0.cnv"};
	// public static final String[] DEFAULT_CNV_FILES = {"conf_100bp_5SNP_10.0.cnv"};
//	public static final String[] DEFAULT_CNV_FILES = {"data/penncnv_1SNP.cnv", "data/quantisnp_1SNP.cnv"};
//	public static final String DEFAULT_GENE_TRACK = GeneSet.DIRECTORY+GeneSet.REFSEQ_TRACK;
//	public static final String DEFAULT_GENE_TRACK = GeneSet.REFSEQ_TRACK;
	public static final boolean SHOW_MIDLINE = true;

	public static final double MOUSE_WHEEL_MULTIPLIER = 0.5;
	public static final int WIDTH_BUFFER = 25;
	public static final int HEIGHT_BUFFER = 10;
//	public static final int DYNAMIC_HEIGHT_LIMIT = 2000;
	public static final int DYNAMIC_HEIGHT_LIMIT = 0;
	public static final int DOUBLE_CLICK_INTERVAL = 500;
	public static final int SIZE = 4;
	public static final int MIN_BUFFER = 10000;
	// public static final String DEFAULT_SAMPLE_DATA = "data/SampleData.dat";
	public static final int DEFAULT_STARTX = 20;
	public static final int DEFAULT_STARTY = 20;
	public static final int DEFAULT_WIDTH = 1100;
	public static final int DEFAULT_HEIGHT = 720;

	// private static final String ALT_UP = "ALT UP";
	// private static final String ALT_DOWN = "ALT DOWN";
	// private static final String ALT_LEFT = "ALT LEFT";
	// private static final String ALT_RIGHT = "ALT RIGHT";

	private static final String FIRST_CHR = "First chr";
	private static final String PREVIOUS_CHR = "Previous chr";
	private static final String NEXT_CHR = "Next chr";
	private static final String LAST_CHR = "Last chr";
	private static final String FIRST_REGION = "First region";
	private static final String PREVIOUS_REGION = "Previous region";
	private static final String NEXT_REGION = "Next region";
	private static final String LAST_REGION = "Last region";
	private static final String TO_SCATTER_PLOT = "To Scatter Plot";
	private static final String REGION_LIST_NEW_FILE = "Load Region File";
//	private static final String REGION_LIST_NEW_FILE = "Add New File(s)...";
	private static final String REGION_LIST_USE_CNVS = "Use CNVs as Regions...";
	private static final String REGION_LIST_PLACEHOLDER = "Select Region File...";

	private JComboBox<String> sampleList;
	private String[] samplesPresent;
	private JTextField navigationField;
//	private JTextField regionsField;
	private JButton firstChr, previousChr, nextChr, lastChr, previousRegion, nextRegion, firstRegion, lastRegion;
	private Project proj;
	private String sample;
	private IndiPheno indiPheno;
	private boolean jar;
	private String[] markerNames;
	private MarkerSet markerSet;
	private long fingerprint;
	private int[] positions;
	private boolean[] dropped;
	private int[][] chrBoundaries;
	private float[] lrrs, lrrValues;
	private float[] bafs, originalBAFs;
	private byte[] genotypes;
	private byte chr;
	private int start, startMarker;
	private int stop, stopMarker;
//	private float lrrMin, lrrMax;
	private boolean inDrag;
	private int startX;
	private String[] cnvFilenames;
	private String[] cnvLabels;
	private CNVariant[][] cnvs;
//	private String[] regionsList;
//	private int regionsListIndex;
	private String[][] regions;
	private int regionIndex;
	private JPanel lrrPanel;
	private SingleClick leftClick;
	private SingleClick rightClick;
	private GeneTrack track;
	private SampleData sampleData;
	private JLabel commentLabel;
	private int transformation_type;
	private boolean transformSeparatelyByChromosome;
	private GcModel gcModel;
	private JCheckBoxMenuItem gcCorrectButton;
	private Hashtable<String, String> namePathMap;
//	private JComboBox<String> centroidsSelection;
	private Logger log;
	private boolean fail;
	boolean isSettingCentroid = false;
	private float[][][] centroids;
	private String currentCentroid;
	private static final String SEX_CENT = "Sex-Specific";
	private JMenu loadRecentFileMenu;
	private ButtonGroup regionButtonGroup;
	private AbstractAction markerFileSelectAction;
	private JMenuItem launchScatter;
	
	// private Color[] colorScheme = {new Color(33, 31, 53), // dark dark
	// new Color(23, 58, 172), // dark blue
	// new Color(201, 30, 10), // deep red
	// new Color(140, 20, 180), // deep purple
	// new Color(33, 87, 0), // dark green
	// new Color(55, 129, 252), // light blue
	// new Color(217, 109, 194), // pink
	// new Color(94, 88, 214), // light purple
	// new Color(189, 243, 61), // light green
	// };

	private Color[][] colorScheme = {
			{new Color(23, 58, 172), new Color(55, 129, 252)}, // dark/light blue
			{new Color(140, 20, 180), new Color(94, 88, 214)}, // deep/light purple
			{new Color(33, 87, 0), new Color(189, 243, 61)}, // dark green
			{new Color(201, 30, 10), new Color(217, 109, 194)}, // deep red/pink
			{new Color(33, 31, 53), new Color(255, 255, 255)} // dark dark/ light light
	};

	final ArrayList<int[]> activeCNVs = new ArrayList<int[]>();
	volatile int[] selectedCNV = null;
	private HashMap<String, JCheckBoxMenuItem> regionFileNameBtn = new HashMap<String, JCheckBoxMenuItem>();
	private HashMap<String, String> regionFileNameLoc = new HashMap<String, String>();
//	private JComboBox<String> regionFileList;
	private String regionFileName;
	private volatile boolean loadingFile = false;

	private JTextField regionField;

	private HashMap<String, JCheckBoxMenuItem> centButtonMap;

	private JCheckBoxMenuItem autoSwitch;
	
	public Trailer(Project proj, String selectedSample, String[] filenames, String location) {
		this(proj, selectedSample, filenames, location, DEFAULT_STARTX, DEFAULT_STARTX, DEFAULT_WIDTH, DEFAULT_HEIGHT);
	}

	// TODO Trailer should have a createAndShowGUI, same as all the other plots, as opposed to being its own frame 
	public Trailer(Project proj, String selectedSample, String[] filenames, String location, int startX, int startY, int width, int height) {
		super("Genvisis - Trailer - " + proj.getNameOfProject());
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		addWindowListener(new WindowAdapter() {
			@Override
			public void windowClosing(WindowEvent e) {
				if (Trailer.this.proj != null) {
					ArrayList<String> files = new ArrayList<String>(regionFileNameLoc.values());
					String[] currSet = Trailer.this.proj.INDIVIDUAL_CNV_LIST_FILENAMES.getValue();
					
					ArrayList<String> newSet = new ArrayList<String>();
					outer: for (String s : files) {
						for (int i = 0; i < currSet.length; i++) {
							if (currSet[i].equals(s)) {
								continue outer;
							}
						}
						newSet.add(s);
					}
					
					if (newSet.size() > 0) {
					    String[] newList = files.toArray(new String[]{});

				        String message = newSet.size() + " files have been added.  ";
				        int choice = JOptionPane.showOptionDialog(null, message+" Would you like to keep this configuration for the next time Trailer is loaded?", "Preserve Trailer workspace?", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, null, null);
				        if (choice == 0) {
				            Trailer.this.proj.INDIVIDUAL_CNV_LIST_FILENAMES.setValue(newList);
				            Trailer.this.proj.saveProperties();
				        }
					}
				}
				super.windowClosing(e);
			}
		});
		
		System.out.println("startX: "+startX+"\t startY: "+startY+"\t width: "+width+"\t height: "+height);
		System.out.println(Array.toStr(filenames, "; "));

		long time;
		String trackFilename;

		this.proj = proj;
		this.log = proj.getLog();
		jar = proj.JAR_STATUS.getValue();
		cnvFilenames = filenames;
		fail = false;
		
		time = new Date().getTime();

		chr = (byte)Positions.parseUCSClocation(location)[0]; 
//		parseLocation(location, false);

		fail = !loadMarkers();
		if (fail) {
			return;
		}
		if (Files.exists(proj.GC_MODEL_FILENAME.getValue(false, false))) {
			gcModel = GcAdjustor.GcModel.populateFromFile(proj.GC_MODEL_FILENAME.getValue(false, false), true, proj.getLog());
		} else {
			gcModel = null;
		}
		generateComponents();
		this.setJMenuBar(createMenuBar());
			
		sample = selectedSample==null?samplesPresent[0]:selectedSample;
		sampleData = proj.getSampleData(2, cnvFilenames);
		if (sampleData.failedToLoad()) {
			proj.getLog().reportError("Without a SampleData file, Trailer will not start");
			return;
		}
		
		cnvLabels = sampleData.getCnvClasses();
//		System.err.println("Error - ");
//		System.err.println("Error - there are "+cnvLabels.length+" cnvLabels: "+Array.toStr(cnvFilenames));
//		System.err.println("Error - ");
//		System.exit(1);
//		loadCNVfiles(proj, cnvFilenames);

//		regionsList = proj.getIndividualRegionLists();
//		if (regionsList.length > 1) {
//			JOptionPane.showMessageDialog(null, "Warning - only one list file is currently supported within Trailer", "Warning", JOptionPane.ERROR_MESSAGE);
//		}
//		regionsListIndex = 0;
//		if (regionsList.length > 0) {
//			if (Files.exists(regionsList[regionsListIndex], jar)) {
//				loadRegions();
//			} else {
//				System.err.println("Error - couldn't find '"+regionsList[regionsListIndex]+"' in data directory; populating with CNVs of current subject");
//			}
//		} else {
//			System.err.println("Warning - no region list was provided; populating regions with CNVs of current subject, if any exist");
//		}
//		regionIndex = -1;

        time = new Date().getTime();
//		track = GeneTrack.load(proj.getDir(Project.DATA_DIRECTORY)+GeneSet.REFSEQ_TRACK, jar);
        if (new File(proj.GENETRACK_FILENAME.getValue(false, false)).isFile()) {
//        	trackFilename = proj.getFilename(proj.GENETRACK_FILENAME);
        	trackFilename = proj.GENETRACK_FILENAME.getValue();
        } else if (new File(GeneSet.DIRECTORY+GeneSet.REFSEQ_TRACK).exists()) {
        	trackFilename = GeneSet.DIRECTORY+GeneSet.REFSEQ_TRACK;
        } else if (new File(GeneSet.REFSEQ_TRACK).exists()) {
        	trackFilename = GeneSet.REFSEQ_TRACK;
        } else {
//			JOptionPane.showMessageDialog(this, "Gene track is not installed. Gene boundaries will not be displayed.", "FYI", JOptionPane.INFORMATION_MESSAGE);
        	trackFilename = null;
        }
        if (trackFilename != null) {
        	log.report("Loading track from "+trackFilename);	
        	track = GeneTrack.load(trackFilename, jar);
        }
        log.report("Loaded track in "+ext.getTimeElapsed(time));
		
		
		updateSample(sample);
		
		System.out.println("All in "+ext.getTimeElapsed(time));

		parseLocation(location);
//		setBounds(20, 20, 1000, 720);
		setBounds(startX, startY, width, height);
		setVisible(true);
		
		Trailer.this.regionFileName = REGION_LIST_USE_CNVS;
		loadCNVsAsRegions();
		procCNVs(chr);
		updateGUI();
		regionIndex = -1;
		if (location.equals(DEFAULT_LOCATION)) {
			regionIndex = 0;
			showRegion();
		}
		
	}
	
	private void paintLRRPanel(Graphics g) {
		// TODO moving paintComponent code here breaks drawing, and I haven't figured out why.. (cole - 3/6/15)
	}
	
	private void paintCNVGeneTrackPanel(Graphics g) {

		GeneData[] genes;
		int[][] exons;
		Vector<Segment> v = new Vector<Segment>();
		Segment[] segs;
		int width, begin, end, source;
		Segment currentView;
		String text;
		
//		g.drawRect(0, 0, this.getWidth()-1, this.getHeight()-1);
		
		
		if (track == null) {
			text = "Gene track is not installed";
			width = g.getFontMetrics(g.getFont()).stringWidth(text);
			g.drawString(text, this.getWidth()/2-width/2, 10);
		} else {
			if (stop-start > 10000000) {
				g.drawString("Zoom in to see genes", 10, 10);
			} else {
				genes = track.getBetween(chr, start, stop, 30);
//				System.out.println(ext.getUCSCformat(new int[] {chr, start, stop}));
				g.setColor(Color.BLACK);
				for (int i = 0; i<genes.length; i++) {
					begin = getX(genes[i].getStart());
					end = getX(genes[i].getStop());
					g.drawRoundRect(begin, 0*15, end-begin, 10, 2, 2);
					v.add(new Segment(begin, end));
					exons = genes[i].getExonBoundaries();
					for (int j = 0; j<exons.length; j++) {
						begin = getX(exons[j][0]);
						end = getX(exons[j][1]);
						if (j==0 || j==exons.length-1) {
							g.fillRoundRect(begin, 0*15, end-begin+1, 10, 2, 2);
						} else {
							g.fillRect(begin, 0*15, end-begin+1, 10);
						}
						
					}
//					System.out.println(genes[i].getGeneName()+"\t"+genes[i].getStart()+"\t"+genes[i].getStop());
                }
//				System.out.println();
				Segment.mergeOverlapsAndSort(v);
				segs = Segment.toArray(v);
				g.setFont(new Font("Arial", 0, 14));
				
				for (int i = 0; i<genes.length; i++) {
					begin = getX(genes[i].getStart());
					width = g.getFontMetrics(g.getFont()).stringWidth(genes[i].getGeneName());
					if (!Segment.overlapsAny(new Segment(begin-width-5, begin-1), segs)) {
						g.drawString(genes[i].getGeneName(), begin-width-3, 0*15+10);
					}
				}
			}
		}
		currentView = new Segment(chr, start, stop);
		activeCNVs.clear();
		int firstBegin;
		for (int i = 0; cnvs != null && i<cnvs.length; i++) {
			source = i;
			firstBegin = Integer.MAX_VALUE;
			for (int j = 0; j<cnvs[i].length; j++) {
				if (cnvs[i][j].overlaps(currentView)) {
					activeCNVs.add(new int[]{i, j});
					begin = getX(cnvs[i][j].getStart());
					if (begin < firstBegin) {
						firstBegin = begin;
					}
					end = getX(cnvs[i][j].getStop());
					g.setColor(colorScheme[source][cnvs[i][j].getCN()<2?0:1]);
					g.fillRoundRect(begin, (source+2)*15, end-begin+1, 10, 2, 2);
					g.setColor(Color.BLACK);
					if (selectedCNV != null && selectedCNV[0] == i && selectedCNV[1] == j) {
						g.drawRoundRect(begin-1, (source+2)*15-1, end-begin+2, 11, 5, 5);
					}
//					g.drawString(ext.rootOf(cnvLabels[source]), begin+2, (source+1)*15+10);
				}
			}
			if (firstBegin != Integer.MAX_VALUE) {
				g.drawString(ext.rootOf(cnvLabels[source]), firstBegin-Grafik.getTextWidth(cnvLabels[source], g)-3, (source+2)*15+10);
			}
		}
	}
	
	private int getX(int pos) {
		return (int)((double)(pos-start)/(double)(stop-start)*(double)(getWidth()-2*WIDTH_BUFFER))+WIDTH_BUFFER;
	}
	
	private void chooseNewFiles() {
		JFileChooser jfc = new JFileChooser((proj != null || regionFileName == null ? proj.PROJECT_DIRECTORY.getValue() : ext.parseDirectoryOfFile(regionFileName)));
		jfc.setMultiSelectionEnabled(true);
		if (jfc.showOpenDialog(Trailer.this) == JFileChooser.APPROVE_OPTION) {
			File[] files = jfc.getSelectedFiles();
			if (files.length > 0) {
				boolean[] keep = Array.booleanArray(files.length, true);
				for (int i = 0; i < files.length; i++) {
					for (String fileName : regionFileNameLoc.keySet()) {
						if (ext.rootOf(files[i].toString()).equals(fileName)) {
							keep[i] = false;
						}
					}
				}
				File[] keptFiles = Array.subArray(files, keep);
				File[] discards = Array.subArray(files, Array.booleanNegative(keep));
				
				if (discards.length > 0) {
					StringBuilder msg = new StringBuilder("The following data file(s) are already present:");
					for (File disc : discards) {
						msg.append("\n").append(disc.getName());
					}
					JOptionPane.showMessageDialog(Trailer.this, msg.toString()); 
				}
				
				for (File kept : keptFiles) {
					addFileToList(kept.getAbsolutePath());
				}
			} else {
				File file = jfc.getSelectedFile();
				boolean keep = true;
				for (String fileName : regionFileNameLoc.keySet()) {
					if (ext.rootOf(file.toString()).equals(fileName)) {
						keep = false;
					}
				}
				
				if (!keep) {
					StringBuilder msg = new StringBuilder("The following data file is already present:\n").append(file.getName());
					JOptionPane.showMessageDialog(Trailer.this, msg.toString()); 
				} else {
					addFileToList(file.getAbsolutePath());
				}
				
			}
			
		}
	}
	
	private void addFileToList(String rawfile) {
		String file = ext.verifyDirFormat(rawfile);
		file = file.substring(0, file.length() - 1);
		String name = ext.rootOf(file);
		regionFileNameLoc.put(name, file);
		
		JCheckBoxMenuItem item = new JCheckBoxMenuItem();
		item.setAction(markerFileSelectAction);
		item.setText(name);

		regionFileNameBtn.put(name, item);
		regionButtonGroup.add(item);
		loadRecentFileMenu.add(item);
		
	}
	
	public void generateComponents() {
		JPanel dataPanel = new JPanel();
		dataPanel.setLayout(new GridLayout(3, 1, 5, 5));

		lrrPanel = new JPanel() {
			public static final long serialVersionUID = 2L;
			public void paintComponent(Graphics g) {
				paintLRRPanel(g);
				

				float min, max;

				if (lrrValues != null) {
					min = Array.min(lrrValues);
					max = Array.max(lrrValues);
					
					// System.out.println("Displaying "+(stopMarker-startMarker)+"
					// markers");
					if (stopMarker-startMarker<DYNAMIC_HEIGHT_LIMIT) {
						min = Float.POSITIVE_INFINITY;
						max = Float.NEGATIVE_INFINITY;
						for (int i = startMarker; i<=stopMarker; i++) {
							if (lrrValues[i]<min) {
								min = lrrValues[i];
							}
							if (lrrValues[i]>max) {
								max = lrrValues[i];
							}
						}
						min = (float)Math.floor(min);
						max = (float)Math.ceil(max);
					}

					// min = -1.5;
					// max = 1.5;

					switch (transformation_type) {
					case -1:
					case 0:
						// Illumina
						min = -3;
						max = 3;
						// Agilent
//						min = -1.5f;
//						max = 1.5f;
						break;
					case 1:
						min = 0;
						max = 1;
						break;
					case 2:
						min = -8;
						max = 8;
						break;
					case 3:
						min = -12;
						max = 12;
						break;
						
					}
					
					g.setFont(new Font("Arial", 0, 20));
					g.drawString("Log R Ratio", WIDTH_BUFFER, 20);
					
					if (SHOW_MIDLINE) {
						g.setColor(Color.LIGHT_GRAY);
						g.drawLine(WIDTH_BUFFER, getHeight()-(int)((double)(0-min)/(double)(max-min)*(double)(getHeight()-2*HEIGHT_BUFFER))-HEIGHT_BUFFER, getWidth()-WIDTH_BUFFER, getHeight()-(int)((double)(0-min)/(double)(max-min)*(double)(getHeight()-2*HEIGHT_BUFFER))-HEIGHT_BUFFER);
					}

					g.setFont(new Font("Arial", 0, 12));
					for (int i = startMarker; i<=stopMarker; i++) {
						// if (genotypes[i] == 1) {
						if (bafs != null && bafs.length > i && bafs[i] > 0.2 && bafs[i] < 0.8) {
							g.setColor(Color.RED);
							// colorScheme[2]
						} else {
							g.setColor(Color.BLACK);
						}
						if (!Float.isNaN(lrrValues[i])) {
							if (dropped[i]) {
//								g.drawString("X", getX(positions[i]), getHeight()-(int)((double)(lrrValues[i]-min)/(double)(max-min)*(double)(getHeight()-2*HEIGHT_BUFFER))-HEIGHT_BUFFER);
							} else if (lrrValues[i] < min){
								g.drawString("v", Trailer.this.getX(positions[i]), getHeight()-(int)((double)(min-min)/(double)(max-min)*(double)(getHeight()-2*HEIGHT_BUFFER))-HEIGHT_BUFFER);
							} else if (lrrValues[i] > max){
								g.drawString("^", Trailer.this.getX(positions[i]), getHeight()-(int)((double)(max-min)/(double)(max-min)*(double)(getHeight()-2*HEIGHT_BUFFER))-HEIGHT_BUFFER);
							} else {
								g.fillOval(Trailer.this.getX(positions[i]), getHeight()-(int)((double)(lrrValues[i]-min)/(double)(max-min)*(double)(getHeight()-2*HEIGHT_BUFFER))-HEIGHT_BUFFER, SIZE, SIZE);
							}
						}
					}
				}
			}
		};
		lrrPanel.addMouseListener(this);
		lrrPanel.addMouseMotionListener(this);
		lrrPanel.addMouseWheelListener(this);
		dataPanel.add(lrrPanel);

		final JPanel cnvPanel = new JPanel() {
			public static final long serialVersionUID = 8L;
			public void paintComponent(Graphics g) {
				paintCNVGeneTrackPanel(g);
			}
		};
		cnvPanel.setMaximumSize(new Dimension(getWidth(), 20));
		MouseAdapter cnvAdapter = new MouseAdapter() {
			int defaultInitial = ToolTipManager.sharedInstance().getInitialDelay();
			int defaultReshow = ToolTipManager.sharedInstance().getReshowDelay();
//			Object defaultBG = UIManager.get("ToolTip.background");
			@Override
			public void mouseEntered(MouseEvent e) {
				super.mouseEntered(e);
				ToolTipManager.sharedInstance().setReshowDelay(3);
				ToolTipManager.sharedInstance().setInitialDelay(3);
			}
			@Override
			public void mouseExited(MouseEvent e) {
				super.mouseExited(e);
				ToolTipManager.sharedInstance().setReshowDelay(defaultReshow);
				ToolTipManager.sharedInstance().setInitialDelay(defaultInitial);
			}
			@Override
			public void mouseMoved(MouseEvent e) {
				super.mouseMoved(e);
				if (selectedCNV == null) return;
				int x = e.getX();
				int y = e.getY();
				
				if (cnvs.length <= selectedCNV[0] || cnvs[selectedCNV[0]].length <= selectedCNV[1]) {
//					UIManager.put("ToolTip.background", defaultBG);
					cnvPanel.setToolTipText(null);
					return;
				}
				
				CNVariant cnv = cnvs[selectedCNV[0]][selectedCNV[1]];
				int cnvX1 = getX(cnv.getStart());
				int cnvX2 = getX(cnv.getStop());
				int cnvY1 = (selectedCNV[0] + 2) * 15;
				int cnvY2 = cnvY1 + 10;
				if (x >= cnvX1 && x <= cnvX2 && y > cnvY1 && y < cnvY2) {
					if (cnvPanel.getToolTipText() != null) {
//						Point locationOnScreen = new Point(e.getXOnScreen(), e.getYOnScreen());
//		                Point locationOnComponent = new Point(e.getX(), e.getY());
//						final MouseEvent me = new MouseEvent(cnvPanel, -1, System.currentTimeMillis(), 0, locationOnComponent.x, locationOnComponent.y, locationOnScreen.x, locationOnScreen.y, 0, false, 0);
//						SwingUtilities.invokeLater(new Runnable() {
//							@Override
//							public void run() {
//								ToolTipManager.sharedInstance().mouseMoved(me);
//							}
//						});
					} else {
						StringBuilder txtBld = new StringBuilder();
						txtBld.append("<html>Start: ").append(ext.addCommas(cnv.getStart())).append("<br/>");
						txtBld.append(" Stop: ").append(ext.addCommas(cnv.getStop())).append("<br/>");
						txtBld.append("# Mkrs: ").append(cnv.getNumMarkers()).append("<br/>");
						txtBld.append("Score: ").append(cnv.getScore()).append("</html>");
//						UIManager.put("ToolTip.background", new javax.swing.plaf.ColorUIResource(colorScheme[selectedCNV[0]][cnv.getCN() < 2 ? 0 : 1]));
	//					cnvPanel.createToolTip();
						cnvPanel.setToolTipText(txtBld.toString());
					}
				} else {
//					UIManager.put("ToolTip.background", defaultBG);
					cnvPanel.setToolTipText(null);
				}
			}
			@Override
			public void mouseClicked(MouseEvent e) {
				super.mouseClicked(e);
				int x = e.getX();
				int y = e.getY();
				
				for (int[] cnvInd : activeCNVs) {
					int yMin = (cnvInd[0] + 2) * 15;
					int yMax = yMin + 10;
					int xBegin = getX(cnvs[cnvInd[0]][cnvInd[1]].getStart());
					int xEnd = getX(cnvs[cnvInd[0]][cnvInd[1]].getStop());
					
					if (x >= xBegin && x <= xEnd && y > yMin && y < yMax) {
						selectedCNV = cnvInd;
						Trailer.this.repaint();
						return;
					}
					
				}
				
				selectedCNV = null;
				
//				cnvPanel.revalidate();
//				cnvPanel.repaint();
//				Trailer.this.revalidate();
				Trailer.this.repaint();
			}
		};
		cnvPanel.addMouseListener(cnvAdapter);
		cnvPanel.addMouseMotionListener(cnvAdapter);
		dataPanel.add(cnvPanel);

		JPanel panel = new JPanel() {
			public static final long serialVersionUID = 7L;

			public void paintComponent(Graphics g) {
				g.setFont(new Font("Arial", 0, 20));
				g.drawString("B Allele Frequency", WIDTH_BUFFER, 20);

				g.setFont(new Font("Arial", 0, 10));
				if (lrrs != null) {
					for (int i = startMarker; i<=stopMarker; i++) {
						if (!Float.isNaN(lrrs[i])) {
							if (dropped[i]) {
								g.drawString("X", Trailer.this.getX(positions[i]), getHeight()-(int)(bafs[i]*(double)(getHeight()-2*HEIGHT_BUFFER))-HEIGHT_BUFFER+5);
							} else if (genotypes != null && genotypes[i]==-1) {
								g.drawString("+", Trailer.this.getX(positions[i]), getHeight()-(int)(bafs[i]*(double)(getHeight()-4*HEIGHT_BUFFER))-HEIGHT_BUFFER/2);
							} else if (bafs != null && bafs.length > i) {
								g.fillOval(Trailer.this.getX(positions[i]), getHeight()-(int)(bafs[i]*(double)(getHeight()-4*HEIGHT_BUFFER))-HEIGHT_BUFFER, SIZE, SIZE);
							} else {
								g.setFont(new Font("Arial", 0, 20));
								g.drawString("Error - no BAF values present for sample.", (this.getWidth() / 2) - 175, (this.getHeight() / 2) - 10);
							}
						}
					}
				}
			}
		};
		dataPanel.add(panel);

		// getContentPane().setBackground(Color.WHITE);
		getContentPane().add(dataPanel, BorderLayout.CENTER);
		
		// JLabel sampleName = new JLabel(sample, JLabel.CENTER);
		// sampleName.setFont(new Font("Arial", 0, 20));
		// descrPanel.add(sampleName);

		JPanel sampPanel = new JPanel();
		previousRegion = new JButton(Grafik.getImageIcon("images/firstLast/Left.gif", true));
		previousRegion.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLeft.gif", true));
		previousRegion.addActionListener(this);
		previousRegion.setActionCommand(PREVIOUS_REGION);
		previousRegion.setPreferredSize(new Dimension(25, 25));
//		sampPanel.add(previousRegion);

		sampleList = new JComboBox<String>();
		sampleList.setFont(new Font("Arial", 0, 20));
		createSampleList();
		
		DefaultListCellRenderer dlcr = new DefaultListCellRenderer();
	    dlcr.setHorizontalAlignment(DefaultListCellRenderer.CENTER);
	    sampleList.setRenderer(dlcr);
		sampleList.setBorder(BorderFactory.createEtchedBorder());
		sampleList.setEditable(false);
		sampleList.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				@SuppressWarnings("unchecked")
				JComboBox<String> jcb = (JComboBox<String>)e.getSource();
				int index = jcb.getSelectedIndex();
//				System.out.println("Selected Index = "+index);
				if (index == samplesPresent.length-1) {
					createSampleList();
				} else if (sample != samplesPresent[index]) {
					updateSample(samplesPresent[index]);
				}
			}
		});
		sampPanel.add(sampleList);
		
		JPanel descrPanel = new JPanel();
		descrPanel.setLayout(new BoxLayout(descrPanel, BoxLayout.Y_AXIS));
		
		nextRegion = new JButton(Grafik.getImageIcon("images/firstLast/Right.gif", true));
		nextRegion.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dRight.gif", true));
		nextRegion.addActionListener(this);
		nextRegion.setActionCommand(NEXT_REGION);
		nextRegion.setPreferredSize(new Dimension(25, 25));
//		sampPanel.add(nextRegion);
		sampPanel.setPreferredSize(new Dimension(sampPanel.getPreferredSize().width, sampleList.getPreferredSize().height + 5));
		descrPanel.add(sampPanel);
		
//		JPanel ctrlPanel = new JPanel(new GridLayout(1, 2, 5, 5));
//		ctrlPanel.setMinimumSize(new Dimension(200, 10));
//		ctrlPanel.setPreferredSize(new Dimension(355, 20));
//		ctrlPanel.setMaximumSize(new Dimension(355, 40));
		
//		JButton clipBtn = new JButton();
//		clipBtn.setFont(clipBtn.getFont().deriveFont(8f));
//		clipBtn.setText("Copy to Clipboard");
//		clipBtn.setMinimumSize(new Dimension(150, 10));
//		clipBtn.setPreferredSize(new Dimension(150, 10));
//		clipBtn.setMaximumSize(new Dimension(150, 10));
//		JButton clipPlusBtn = new JButton();
//		clipPlusBtn.setFont(clipPlusBtn.getFont().deriveFont(8f));
//		clipPlusBtn.setText("Copy to Clipboard w/ Position");
//		clipPlusBtn.setMinimumSize(new Dimension(150, 10));
//		clipPlusBtn.setPreferredSize(new Dimension(150, 10));
//		clipPlusBtn.setMaximumSize(new Dimension(150, 10));
//
//		ctrlPanel.add(clipBtn);
//		ctrlPanel.add(clipPlusBtn);
//		descrPanel.add(ctrlPanel);
		descrPanel.add(Box.createVerticalGlue());
		
//		Vector<String> items = new Vector<String>();
//		items.add(REGION_LIST_PLACEHOLDER);
//		if (proj != null) {
//			String[] files = proj.getIndividualRegionLists();
//			String name;
//			for (String file : files) {
//				name = ext.rootOf(file);
//				regionFileNameLoc.put(name, file);
//				items.add(name);
//			}
//		}
//		items.add(REGION_LIST_USE_CNVS);
//		items.add(REGION_LIST_NEW_FILE);
//		regionFileList = new JComboBox<String>(items);
//		regionFileList.setFont(new Font("Arial", 0, 12));
//		regionFileList.setMinimumSize(new Dimension(200, 20));
//		regionFileList.setPreferredSize(new Dimension(200, 20));
//		regionFileList.setMaximumSize(new Dimension(200, 20));
//		AbstractAction markerFileSelectAction = new AbstractAction() {
//			private static final long serialVersionUID = 1L;
//	
//			@SuppressWarnings("unchecked")
//			@Override
//			public void actionPerformed(ActionEvent e) {
//				String shortName = (String) ((JComboBox<String>)e.getSource()).getSelectedItem();
//				if (!loadingFile 
//						&& !REGION_LIST_NEW_FILE.equals(shortName) 
//						&& !REGION_LIST_PLACEHOLDER.equals(shortName)
//						&& !REGION_LIST_USE_CNVS.equals(shortName)) {
//					String file = regionFileNameLoc.get(shortName);
//					if (file != null && file.equals(Trailer.this.regionFileName)) {
//						return;
//					}
//					Trailer.this.regionFileName = file;
//					loadRegions();
//				} /*else if (loadingFile && REGION_LIST_PLACEHOLDER.equals(shortName)) {
//					// do nothing
//				} */else if (loadingFile || REGION_LIST_PLACEHOLDER.equals(shortName)) {
//					// leave as currently selected marker
//					if (Trailer.this.regionFileName != "" && Trailer.this.regionFileName != null) {
//						String file = Trailer.this.regionFileName;
//						if (!REGION_LIST_USE_CNVS.equals(Trailer.this.regionFileName)) {
//							file = ext.rootOf(Trailer.this.regionFileName);
//						}
//						((JComboBox<String>)e.getSource()).setSelectedItem(file);
//					}
//					return;
//				} else if (REGION_LIST_USE_CNVS.equals(shortName)) {
//					if (!shortName.equals(Trailer.this.regionFileName)) {
//						Trailer.this.regionFileName = REGION_LIST_USE_CNVS;
//						loadCNVsAsRegions();
//						procCNVs(chr);
//						updateGUI();
//					}
//				} else if (REGION_LIST_NEW_FILE.equals(shortName)) {
//					chooseNewFiles();
//					if (Trailer.this.regionFileName != null && !"".equals(Trailer.this.regionFileName)) {
//						if (REGION_LIST_USE_CNVS.equals(Trailer.this.regionFileName)) {
//							((JComboBox<String>)e.getSource()).setSelectedItem(REGION_LIST_USE_CNVS);
//						} else {
//							((JComboBox<String>)e.getSource()).setSelectedItem(ext.rootOf(Trailer.this.regionFileName));
//						}
//					} else {
//						((JComboBox<String>)e.getSource()).setSelectedItem(REGION_LIST_PLACEHOLDER);
//					}
//				}
//			}
//		};
//		regionFileList.addActionListener(markerFileSelectAction);
		
		JPanel compPanel = new JPanel(new GridLayout(2, 1, 5, 5));
//		compPanel.setBorder(new LineBorder(Color.BLACK, 1));

		JPanel regionPanel = new JPanel();
		regionField = new JTextField("", 8);
		regionField.setHorizontalAlignment(JTextField.CENTER);
		regionField.setFont(new Font("Arial", 0, 14));
		//navigationField.setEditable(false);
//		regionField.setBackground(BACKGROUND_COLOR);
		regionField.addActionListener(new ActionListener() {
	        public void actionPerformed(ActionEvent e) {
				try {
					int trav = Integer.valueOf(((JTextField)e.getSource()).getText().split("[\\s]+")[0]).intValue()-1;
					if (trav >=0 && trav < regions.length) {
						regionIndex = trav;
						showRegion();
					}
				} catch (NumberFormatException nfe) {}
				updateGUI();
	        }
		});
		regionPanel.add(previousRegion);
		regionPanel.add(regionField);
		regionField.setPreferredSize(new Dimension(regionField.getPreferredSize().width, 26));
		regionPanel.add(nextRegion);
		compPanel.add(regionPanel);
		
//		launchScatterButton = new JButton();
//		launchScatterButton.setFont(launchScatterButton.getFont().deriveFont(10f));
//		launchScatterButton.setText(TO_SCATTER_PLOT);
//		launchScatterButton.setActionCommand(TO_SCATTER_PLOT);
//		launchScatterButton.addActionListener(this);
//		compPanel.setMaximumSize(new Dimension(200, 20));
//		compPanel.add(launchScatterButton);
		
		commentLabel = new JLabel(" ", JLabel.CENTER);
		commentLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
		commentLabel.setFont(new Font("Arial", 0, 14));
//		commentLabel.setBorder(new LineBorder(Color.RED, 1));
		compPanel.add(commentLabel);
		
		descrPanel.add(compPanel);
		compPanel.setPreferredSize(new Dimension(compPanel.getPreferredSize().width, 75));


		JPanel navigateChrPanel = new JPanel();
		firstChr = new JButton(Grafik.getImageIcon("images/firstLast/First.gif", true));
		firstChr.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dFirst.gif", true));
		firstChr.addActionListener(this);
		firstChr.setActionCommand(FIRST_CHR);
		firstChr.setPreferredSize(new Dimension(20, 20));
		previousChr = new JButton(Grafik.getImageIcon("images/firstLast/Left.gif", true));
		previousChr.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLeft.gif", true));
		previousChr.addActionListener(this);
		previousChr.setActionCommand(PREVIOUS_CHR);
		previousChr.setPreferredSize(new Dimension(20, 20));
		navigationField = new JTextField("", 20);
		navigationField.setHorizontalAlignment(JTextField.CENTER);
		navigationField.setFont(new Font("Arial", 0, 14));
		navigationField.addActionListener(new ActionListener() {
	        public void actionPerformed(ActionEvent e) {
	        	JTextField navField;
	        	
	        	navField = (JTextField)e.getSource();
	        	navField.setText(navField.getText().trim());
				parseLocation(navField.getText());
	        }
		});

//		navigationField.addFocusListener(new FocusListener() {
//			public void focusGained(FocusEvent focusevent) {}
//
//			public void focusLost(FocusEvent fe) {
//				parseLocation(((JTextField)fe.getSource()).getText());
//			}
//		});

		nextChr = new JButton(Grafik.getImageIcon("images/firstLast/Right.gif", true));
		nextChr.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dRight.gif", true));
		nextChr.addActionListener(this);
		nextChr.setActionCommand(NEXT_CHR);
		nextChr.setPreferredSize(new Dimension(20, 20));
		lastChr = new JButton(Grafik.getImageIcon("images/firstLast/Last.gif", true));
		lastChr.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLast.gif", true));
		lastChr.addActionListener(this);
		lastChr.setActionCommand(LAST_CHR);
		lastChr.setPreferredSize(new Dimension(20, 20));
		navigateChrPanel.add(firstChr);
		navigateChrPanel.add(previousChr);
		navigateChrPanel.add(navigationField);
		navigateChrPanel.add(nextChr);
		navigateChrPanel.add(lastChr);
		descrPanel.add(navigateChrPanel);

//		JPanel dataOptionPanel = new JPanel(new GridLayout(2, 2, 5, 5));
//		dataOptionPanel.setBorder(new EmptyBorder(0, 25, 0, 10));
//		JPanel transformationPanel = new JPanel();
//		transformationPanel.setLayout(new BoxLayout(transformationPanel, BoxLayout.PAGE_AXIS));
//		JLabel label = new JLabel("Log R Ratio transformation: ");
//		label.setFont(new Font("Arial", 0, 14));
//		label.setAlignmentX(Component.LEFT_ALIGNMENT);
//		label.setBorder(new EmptyBorder(0, 0, 5, 0));
//		transformationPanel.add(label);
//		transformationPanel.add(Box.createHorizontalGlue());
//		
//		ButtonGroup typeRadio = new ButtonGroup();
//		final JRadioButton[] transformationRadioButtons = new JRadioButton[Transforms.TRANFORMATIONS.length];
//		ItemListener typeListener = new ItemListener() {
//			public void itemStateChanged(ItemEvent ie) {
//				JRadioButton jrb = (JRadioButton) ie.getItem();
//				if (jrb.isSelected()) {
//					centroids = null;
//					bafs = originalBAFs;
//					for (int i = 0; i < Transforms.TRANFORMATIONS.length; i++) {
//						if (jrb.getText().equals(Transforms.TRANFORMATIONS[i])) {
//							transformation_type = i;
//							lrrValues = getNewLRRs(proj, lrrs, transformation_type, transformSeparatelyByChromosome, markerSet, gcModel, gcCorrectButton.isSelected(), true, log);
//							updateGUI();
//						}
//					}
//				}
//			}
//		};
//		ItemListener centListener = new ItemListener() {
//			@Override
//			public void itemStateChanged(ItemEvent ie) {
//				JRadioButton jrb = (JRadioButton) ie.getItem();
//				if (jrb.isSelected()) {
//					
//					if (namePathMap == null || namePathMap.isEmpty()) {
//						jrb.setSelected(false);
//						if (transformationRadioButtons != null && transformationRadioButtons[0] != null) {
//							transformationRadioButtons[0].setSelected(true);
//						}
//						return;
//					}
//					
//					transformation_type = -1;
//					setCentroid();
//					loadValues();
//					updateGUI();
//					repaint();
//				} else {
//					currentCentroid = null;
//					centroids = null;
//					bafs = originalBAFs;
//				}
//			}
//		};
//		
//		
//		for (int i = 0; i<Transforms.TRANFORMATIONS.length; i++) {
//			transformationRadioButtons[i] = new JRadioButton(Transforms.TRANFORMATIONS[i], false);
//			transformationRadioButtons[i].setFont(new Font("Arial", 0, 14));
//			typeRadio.add(transformationRadioButtons[i]);
//			transformationRadioButtons[i].addItemListener(typeListener);
//			transformationRadioButtons[i].setMargin(new Insets(0,0,0,0));
////			transformationRadioButtons[i].setBackground(BACKGROUND_COLOR);
//
//			JPanel transOptPanel = new JPanel();
//			transOptPanel.setLayout(new BoxLayout(transOptPanel, BoxLayout.LINE_AXIS));
////			transformationPanel.add(transformationRadioButtons[i]);
//			transOptPanel.add(Box.createRigidArea(new Dimension(25, 1)));
//			transOptPanel.add(transformationRadioButtons[i]);
//			transOptPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
//			transformationPanel.add(transOptPanel);
//		}
//		
//		ItemListener gcListener = new ItemListener() {
//			public void itemStateChanged(ItemEvent ie) {
//				JCheckBox jrb = (JCheckBox) ie.getItem();
//				lrrValues = getNewLRRs(proj, lrrs, transformation_type, transformSeparatelyByChromosome, markerSet, gcModel, jrb.isSelected(), jrb.isSelected(), log);
////				if (jrb.isSelected()) {
////					// gc correct, and apply any transformation
////					lrrValues = getNewLRRs(proj, lrrs, transformation_type, transformSeparatelyByChromosome, markerSet, gcModel, true, true, log);
////				} else {
////					// reset any transformation
////					lrrValues = getNewLRRs(proj, lrrs, transformation_type, transformSeparatelyByChromosome, markerSet, gcModel, false, false, log);
////				}
//				updateGUI();
//			}
//		};

//		gcCorrectButton = new JCheckBox(GcAdjustor.GC_ADJUSTOR_TITLE[0], false);// stays hidden if gcModel is not detected
//		if (gcModel != null) {
//			gcCorrectButton.setFont(new Font("Arial", 0, 14));
//			gcCorrectButton.setToolTipText("GC correction will be applied prior to any transformation");
//			gcCorrectButton.addItemListener(gcListener);
////			transformationPanel.add(gcCorrectButton);
//		}

//		transformationRadioButtons[0].setSelected(true);
//		dataOptionPanel.add(transformationPanel);
		//descrPanel.add(transformationPanel);
		

//		JPanel centroidPanel = new JPanel();
//		centroidPanel.setLayout(new BoxLayout(centroidPanel, BoxLayout.PAGE_AXIS));
//		label = new JLabel("Derive from Centroids:");
//		label.setFont(new Font("Arial", 0, 14));
//		label.setAlignmentX(Component.LEFT_ALIGNMENT);
//		label.setBorder(new EmptyBorder(0, 0, 5, 0));
//		centroidPanel.add(label);
//		JPanel centSubPanel = new JPanel();
//		centSubPanel.setLayout(new BoxLayout(centSubPanel, BoxLayout.LINE_AXIS));
//		centSubPanel.add(Box.createRigidArea(new Dimension(25, 1)));
//		centSubPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
//		JRadioButton centTransformButton = new JRadioButton("", false);
//		centTransformButton.setAlignmentX(Component.LEFT_ALIGNMENT);
//		centTransformButton.setMargin(new Insets(0,0,0,0));
//		centTransformButton.addItemListener(centListener);
//		typeRadio.add(centTransformButton);
//		centSubPanel.add(centTransformButton);
		
//		namePathMap = new Hashtable<String, String>();
//		Vector<String> centFiles = new Vector<String>();
//		centFiles.add(proj.getFilename(Project.ORIGINAL_CENTROIDS_FILENAME));
//		centFiles.add(proj.getFilename(Project.GENOTYPE_CENTROIDS_FILENAME));
//		centFiles.add(proj.getFilename(Project.CUSTOM_CENTROIDS_FILENAME));
//		centFiles.add(proj.getFilename(Project.CHIMERA_CENTROIDS_FILENAME));
//		
//		String[] tempFiles = proj.getFilenames(Project.SEX_CENTROIDS_FILENAMES);
//		if (tempFiles != null && tempFiles.length > 0) {
//			centFiles.add(SEX_CENT);
//		}
//		
//		for (String file : centFiles) {
//			if (SEX_CENT.equals(file)) {
//				namePathMap.put(SEX_CENT + " - Male", tempFiles[0]);
//				namePathMap.put(SEX_CENT + " - Female", tempFiles[1]);
//			} else if (Files.exists(file)) {
//				String name = file.substring(file.lastIndexOf("/") + 1);
//				name = name.substring(0, name.lastIndexOf("."));
//				namePathMap.put(name, file);
//			}
//		}
//		
//		centroidsSelection = new JComboBox<String>((String[]) namePathMap.keySet().toArray(new String[]{}));
//		centroidsSelection.setMaximumSize(new Dimension(160, 25));
//		centroidsSelection.addActionListener(new ActionListener() {
//			@Override
//			public void actionPerformed(ActionEvent e) {
//				if (!isSettingCentroid) {
//					isSettingCentroid = true;
//					if (transformation_type == -1) {
//						setCentroid();
//						loadValues();
//						updateGUI();
//						repaint();
//					}
//					isSettingCentroid = false;
//				} else {
//					setCentroid();
//				}
//			}
//		});
//		centSubPanel.add(centroidsSelection);
//		centroidPanel.add(centSubPanel);
//		dataOptionPanel.add(centroidPanel);
		
		
//		JPanel scopePanel = new JPanel();
//		scopePanel.setLayout(new BoxLayout(scopePanel, BoxLayout.PAGE_AXIS));
//		label = new JLabel("Transform by: ");
//		label.setFont(new Font("Arial", 0, 14));
//		label.setAlignmentX(Component.LEFT_ALIGNMENT);
//		label.setBorder(new EmptyBorder(0, 0, 5, 0));
//		scopePanel.add(label);
//		ButtonGroup scopeRadio = new ButtonGroup();
//		JRadioButton[] scopeRadioButtons = new JRadioButton[Transforms.SCOPES.length];
//		ItemListener scopeListener = new ItemListener() {
//			public void itemStateChanged(ItemEvent ie) {
//				JRadioButton jrb = (JRadioButton)ie.getItem();
//				if (jrb.isSelected()) {
//					transformSeparatelyByChromosome = jrb.getText().equals(Transforms.SCOPES[1]); 
//					lrrValues = getNewLRRs(proj, lrrs, transformation_type, transformSeparatelyByChromosome, markerSet, gcModel, gcCorrectButton.isSelected(), true, log);
//					updateGUI();
//				}
//			}
//		};
//		for (int i = 0; i<Transforms.SCOPES.length; i++) {
//			scopeRadioButtons[i] = new JRadioButton(Transforms.SCOPES[i], false);
//			scopeRadioButtons[i].setFont(new Font("Arial", 0, 14));
//			scopeRadio.add(scopeRadioButtons[i]);
//			scopeRadioButtons[i].addItemListener(scopeListener);
//			scopeRadioButtons[i].setAlignmentX(Component.LEFT_ALIGNMENT);
////			scopeRadioButtons[i].setBackground(BACKGROUND_COLOR);
//			scopeRadioButtons[i].setMargin(new Insets(0,0,0,0));
//
//			JPanel scopeOptPanel = new JPanel();
//			scopeOptPanel.setLayout(new BoxLayout(scopeOptPanel, BoxLayout.LINE_AXIS));
//			scopeOptPanel.setBorder(null);
////			transformationPanel.add(transformationRadioButtons[i]);
//			scopeOptPanel.add(Box.createRigidArea(new Dimension(25, 1)));
//			scopeOptPanel.add(scopeRadioButtons[i]);
//			scopeOptPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
//			scopePanel.add(scopeOptPanel);
//			
//			
////			scopePanel.add(scopeRadioButtons[i]);
//		}
//		scopeRadioButtons[0].setSelected(true);
//		//descrPanel.add(scopePanel);
//		dataOptionPanel.add(scopePanel);
		
		
		JPanel overPanel = new JPanel();
//		overPanel.setLayout(new FlowLayout(FlowLayout.LEFT, 15, 2));
		overPanel.setLayout(new BoxLayout(overPanel, BoxLayout.LINE_AXIS));
//		overPanel.setLayout(new GridLayout(3, 1, 5, 5));
		
		JSeparator sep = new JSeparator(SwingConstants.VERTICAL);
		sep.setMaximumSize(new Dimension(1, 150));
		overPanel.setBorder(new EmptyBorder(5, 0, 5, 0));
		
		overPanel.add(Box.createHorizontalGlue());
		overPanel.add(descrPanel);
//		overPanel.add(sep);
//		overPanel.add(dataOptionPanel);
		overPanel.add(Box.createHorizontalGlue());
		
//		getContentPane().add(descrPanel, BorderLayout.NORTH);
		getContentPane().add(overPanel, BorderLayout.NORTH);

		InputMap inputMap = panel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
		// inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_UP,
		// InputEvent.ALT_MASK), ALT_UP);
		// inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_DOWN,
		// InputEvent.ALT_MASK), ALT_DOWN);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, InputEvent.ALT_MASK), PREVIOUS_REGION);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, InputEvent.ALT_MASK), NEXT_REGION);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_UP, InputEvent.CTRL_MASK), PREVIOUS_CHR);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_PAGE_DOWN, InputEvent.CTRL_MASK), NEXT_CHR);

		ActionMap actionMap = panel.getActionMap();
		actionMap.put(FIRST_CHR, new AbstractAction() {
			public static final long serialVersionUID = 5L;
			public void actionPerformed(ActionEvent e) {
				parseLocation("chr1");
			}
		});
		actionMap.put(PREVIOUS_CHR, new AbstractAction() {
			public static final long serialVersionUID = 6L;
			public void actionPerformed(ActionEvent e) {
				parseLocation("chr"+Math.max(chr-1, 1));
			}
		});
		actionMap.put(NEXT_CHR, new AbstractAction() {
			public static final long serialVersionUID = 7L;
			public void actionPerformed(ActionEvent e) {
				parseLocation("chr"+Math.min(chr+1, 26));
			}
		});
		actionMap.put(LAST_CHR, new AbstractAction() {
			public static final long serialVersionUID = 8L;
			public void actionPerformed(ActionEvent e) {
				parseLocation("chr26");
			}
		});
		actionMap.put(FIRST_REGION, new AbstractAction() {
			public static final long serialVersionUID = 9L;
			public void actionPerformed(ActionEvent e) {
				regionIndex = 0;
				showRegion();
			}
		});
		actionMap.put(PREVIOUS_REGION, new AbstractAction() {
			public static final long serialVersionUID = 10L;
			public void actionPerformed(ActionEvent e) {
				regionIndex = Math.max(regionIndex-1, 0);
				showRegion();
			}
		});
		actionMap.put(NEXT_REGION, new AbstractAction() {
			public static final long serialVersionUID = 11L;
			public void actionPerformed(ActionEvent e) {
				regionIndex = Math.min(regionIndex+1, regions.length-1);
				showRegion();
			}
		});
		actionMap.put(LAST_REGION, new AbstractAction() {
			public static final long serialVersionUID = 12L;
			public void actionPerformed(ActionEvent e) {
				regionIndex = regions.length-1;
				showRegion();
			}
		});
		panel.setActionMap(actionMap);

		// doesn't seem to get captured properly...
		nextChr.getInputMap().put(KeyStroke.getKeyStroke("space"), NEXT_REGION);
		nextChr.setActionMap(actionMap);
		previousChr.setActionMap(actionMap);
	}
	
	private JMenuBar createMenuBar() {
		JMenuBar menuBar = new JMenuBar();
		
		JMenu fileMenu = new JMenu("File");
		fileMenu.setMnemonic(KeyEvent.VK_F);
		AbstractAction loadNewFileAction = new AbstractAction() {
			private static final long serialVersionUID = 1L;
			@Override
			public void actionPerformed(ActionEvent e) {
				chooseNewFiles();
				if (Trailer.this.regionFileName != null && !"".equals(Trailer.this.regionFileName)) {
					if (REGION_LIST_USE_CNVS.equals(Trailer.this.regionFileName)) {
						regionFileNameBtn.get(REGION_LIST_USE_CNVS).setSelected(true);
//						((JComboBox<String>)e.getSource()).setSelectedItem(REGION_LIST_USE_CNVS);
					} else {
						regionFileNameBtn.get(ext.rootOf(Trailer.this.regionFileName)).setSelected(true);
//						((JComboBox<String>)e.getSource()).setSelectedItem(ext.rootOf(Trailer.this.regionFileName));
					}
				}/* else {
					((JComboBox<String>)e.getSource()).setSelectedItem(REGION_LIST_PLACEHOLDER);
				}*/	
			}
		};
		markerFileSelectAction = new AbstractAction() {
			private static final long serialVersionUID = 1L;
			@Override
			public void actionPerformed(ActionEvent e) {
				String shortName = ((JCheckBoxMenuItem)e.getSource()).getText();
				if (!loadingFile 
						&& !REGION_LIST_NEW_FILE.equals(shortName) 
						&& !REGION_LIST_PLACEHOLDER.equals(shortName)
						&& !REGION_LIST_USE_CNVS.equals(shortName)) {
					String file = regionFileNameLoc.get(shortName);
					if (file != null && file.equals(Trailer.this.regionFileName)) {
						return;
					}
					String tempFile = file.startsWith("./") ? proj.PROJECT_DIRECTORY.getValue() + file : file;
					if (!Files.exists(tempFile)) {
						proj.message("Error - region file '" + shortName + "' doesn't exist.");
						regionFileNameBtn.get(regionFileName).setSelected(true);
					} else {
						Trailer.this.regionFileName = file;
						loadRegions();
						regionIndex = 0;
						showRegion();
					}
				} /*else if (loadingFile && REGION_LIST_PLACEHOLDER.equals(shortName)) {
					// do nothing
				} */else if (loadingFile || REGION_LIST_PLACEHOLDER.equals(shortName)) {
					// leave as currently selected marker
					if (Trailer.this.regionFileName != "" && Trailer.this.regionFileName != null) {
						String file = Trailer.this.regionFileName;
						if (!REGION_LIST_USE_CNVS.equals(Trailer.this.regionFileName)) {
							file = ext.rootOf(Trailer.this.regionFileName);
						}
//						((JComboBox<String>)e.getSource()).setSelectedItem(file);
						regionFileNameBtn.get(file).setSelected(true);
					}
					return;
				} else if (REGION_LIST_USE_CNVS.equals(shortName)) {
					if (!shortName.equals(Trailer.this.regionFileName)) {
						Trailer.this.regionFileName = REGION_LIST_USE_CNVS;
						loadCNVsAsRegions();
						procCNVs(chr);
						updateGUI();
						regionIndex = 0;
						showRegion();
					}
				} /*else if (REGION_LIST_NEW_FILE.equals(shortName)) {
					chooseNewFiles();
					if (Trailer.this.regionFileName != null && !"".equals(Trailer.this.regionFileName)) {
						if (REGION_LIST_USE_CNVS.equals(Trailer.this.regionFileName)) {
							((JComboBox<String>)e.getSource()).setSelectedItem(REGION_LIST_USE_CNVS);
						} else {
							((JComboBox<String>)e.getSource()).setSelectedItem(ext.rootOf(Trailer.this.regionFileName));
						}
					} else {
						((JComboBox<String>)e.getSource()).setSelectedItem(REGION_LIST_PLACEHOLDER);
					}
				}*/
			}
		};
		JMenuItem loadRegionFile = new JMenuItem();
		loadRegionFile.setAction(loadNewFileAction);
		loadRegionFile.setText("Load Region File");
		loadRegionFile.setMnemonic(KeyEvent.VK_L);
		fileMenu.add(loadRegionFile);
		loadRecentFileMenu = new JMenu("Load Recent...");
		loadRecentFileMenu.setMnemonic(KeyEvent.VK_R);
		fileMenu.add(loadRecentFileMenu);

		menuBar.add(fileMenu);
		
		regionButtonGroup = new ButtonGroup();
		if (proj != null) {
			String[] files = proj.getIndividualRegionLists();
			String name;
			for (String file : files) {
				name = ext.rootOf(file);
				regionFileNameLoc.put(name, file);
				JCheckBoxMenuItem menuItem = new JCheckBoxMenuItem();
				menuItem.setAction(markerFileSelectAction);
				menuItem.setText(name);
				regionFileNameBtn.put(name, menuItem);
				regionButtonGroup.add(menuItem);
				loadRecentFileMenu.add(menuItem);
			}
		}

		JCheckBoxMenuItem item1 = new JCheckBoxMenuItem();
		item1.setAction(markerFileSelectAction);
		item1.setText(REGION_LIST_USE_CNVS);
		regionFileNameBtn.put(REGION_LIST_USE_CNVS, item1);
		regionButtonGroup.add(item1);
		loadRecentFileMenu.add(item1);
		item1.setSelected(true);
		
		gcCorrectButton = new JCheckBoxMenuItem(GcAdjustor.GC_ADJUSTOR_TITLE[0], false);// stays hidden if gcModel is not detected
		
		JMenu act = new JMenu("Actions");
        act.setMnemonic(KeyEvent.VK_A);
		menuBar.add(act);
		
		launchScatter = new JMenuItem();
		launchScatter.setText(TO_SCATTER_PLOT);
		launchScatter.setMnemonic(KeyEvent.VK_S);
		launchScatter.setFont(new Font("Arial", 0, 12));
		launchScatter.addActionListener(this);
		act.add(launchScatter);
		act.addSeparator();
		
//		launchScatterButton = new JButton();
//		launchScatterButton.setFont(launchScatterButton.getFont().deriveFont(10f));
//		launchScatterButton.setText(TO_SCATTER_PLOT);
//		launchScatterButton.setActionCommand(TO_SCATTER_PLOT);
//		launchScatterButton.addActionListener(this);
//		compPanel.setMaximumSize(new Dimension(200, 20));
//		compPanel.add(launchScatterButton);
		
		ButtonGroup lrrBtnGrp = new ButtonGroup();
		ButtonGroup transBtnGrp = new ButtonGroup();
		
		JMenuItem lbl1 = act.add("LRR Transforms:");
		lbl1.setEnabled(false);
		lbl1.setFont(new Font("Arial", 0, 12));
		
		final JRadioButtonMenuItem[] transformBtns = new JRadioButtonMenuItem[Transforms.TRANFORMATIONS.length];
		ItemListener typeListener = new ItemListener() {
			public void itemStateChanged(ItemEvent ie) {
				JRadioButtonMenuItem jrb = (JRadioButtonMenuItem) ie.getItem();
				if (jrb.isSelected()) {
					centroids = null;
					bafs = originalBAFs;
					for (int i = 0; i < Transforms.TRANFORMATIONS.length; i++) {
						if (jrb.getText().equals(Transforms.TRANFORMATIONS[i])) {
							transformation_type = i;
							lrrValues = getNewLRRs(proj, lrrs, transformation_type, transformSeparatelyByChromosome, markerSet, gcModel, gcCorrectButton.isSelected(), true, log);
							updateGUI();
						}
					}
				}
			}
		};
		HashSet<String> mnems = new HashSet<String>();
		for (int i = 0; i < Transforms.TRANFORMATIONS.length; i++) {
			transformBtns[i] = new JRadioButtonMenuItem(Transforms.TRANFORMATIONS[i]);
			transformBtns[i].addItemListener(typeListener);
			transformBtns[i].setFont(new Font("Arial", 0, 12));
			String[] wds = Transforms.TRANFORMATIONS[i].split("[\\s]+");
			int ind = 0;
			String mnem = wds[ind].substring(0, 1);
			while(mnems.contains(mnem)) {
				mnem = wds[++ind].substring(0, 1);
			}
			mnems.add(mnem);
			transformBtns[i].setMnemonic(mnem.charAt(0));
			lrrBtnGrp.add(transformBtns[i]);
			act.add(transformBtns[i]);
		}
		transformBtns[0].setSelected(true);
		
		act.addSeparator();
		JMenuItem lbl2 = act.add("Transform By:");
		lbl2.setEnabled(false);
		lbl2.setFont(new Font("Arial", 0, 12));
		
		JRadioButtonMenuItem[] scopeBtns = new JRadioButtonMenuItem[Transforms.SCOPES.length];
		ItemListener scopeListener = new ItemListener() {
			public void itemStateChanged(ItemEvent ie) {
				JRadioButtonMenuItem jrb = (JRadioButtonMenuItem)ie.getItem();
				if (jrb.isSelected()) {
					transformSeparatelyByChromosome = jrb.getText().equals(Transforms.SCOPES[1]); 
					lrrValues = getNewLRRs(proj, lrrs, transformation_type, transformSeparatelyByChromosome, markerSet, gcModel, gcCorrectButton.isSelected(), true, log);
					updateGUI();
				}
			}
		};
		mnems = new HashSet<String>();
		for (int i = 0; i<Transforms.SCOPES.length; i++) {
			scopeBtns[i] = new JRadioButtonMenuItem(Transforms.SCOPES[i]);
			scopeBtns[i].addItemListener(scopeListener);
			scopeBtns[i].setFont(new Font("Arial", 0, 12));
			
			String[] wds = Transforms.SCOPES[i].split("[\\s]+");
			int ind = 0;
			String mnem = wds[ind].substring(0, 1);
			while(mnems.contains(mnem)) {
				mnem = wds[++ind].substring(0, 1);
			}
			mnems.add(mnem);
			scopeBtns[i].setMnemonic(mnem.charAt(0));
			transBtnGrp.add(scopeBtns[i]);
			act.add(scopeBtns[i]);
		}
		scopeBtns[0].setSelected(true);
		
		ItemListener gcListener = new ItemListener() {
			public void itemStateChanged(ItemEvent ie) {
				JCheckBoxMenuItem jrb = (JCheckBoxMenuItem) ie.getItem();
				lrrValues = getNewLRRs(proj, lrrs, transformation_type, transformSeparatelyByChromosome, markerSet, gcModel, jrb.isSelected(), jrb.isSelected(), log);
				updateGUI();
			}
		};
		
		gcCorrectButton.setToolTipText("GC correction will be applied prior to any transformation");
		gcCorrectButton.addItemListener(gcListener);
		gcCorrectButton.setFont(new Font("Arial", 0, 12));
		act.addSeparator();
		act.add(gcCorrectButton).setEnabled(gcModel != null);
		
		act.addSeparator();

		autoSwitch = new JCheckBoxMenuItem();
		autoSwitch.setText("Auto-Select Sex Centroid by Sample Sex");
		autoSwitch.setMnemonic(KeyEvent.VK_A);
		act.add(autoSwitch);
		
		ItemListener centListener = new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent ie) {
				JCheckBoxMenuItem jrb = (JCheckBoxMenuItem) ie.getItem();
				if (jrb.isSelected() && !"None".equals(jrb.getText())) {
					if (namePathMap == null || namePathMap.isEmpty()) {
						jrb.setSelected(false);
						if (transformBtns != null && transformBtns[0] != null) {
							transformBtns[0].setSelected(true);
						}
						return;
					}
					if (!isSettingCentroid) {
						isSettingCentroid = true;
						transformation_type = -1;
						setCentroid(jrb.getText());
						loadValues();
						updateGUI();
						repaint();
						isSettingCentroid = false;
					} else {
						setCentroid(jrb.getText());
					}
				} else {
					currentCentroid = null;
					centroids = null;
					bafs = originalBAFs;
					loadValues();
					updateGUI();
					repaint();
				}
			}
		};
		namePathMap = new Hashtable<String, String>();
		Vector<String> centFiles = new Vector<String>();
		centFiles.add(proj.ORIGINAL_CENTROIDS_FILENAME.getValue());
		centFiles.add(proj.GENOTYPE_CENTROIDS_FILENAME.getValue());
		centFiles.add(proj.CUSTOM_CENTROIDS_FILENAME.getValue());
		centFiles.add(proj.CHIMERA_CENTROIDS_FILENAME.getValue());
		
		String tempFileMale = proj.SEX_CENTROIDS_MALE_FILENAME.getValue();
		String tempFileFemale = proj.SEX_CENTROIDS_FEMALE_FILENAME.getValue();
		if (tempFileMale != null && tempFileFemale != null && !"".equals(tempFileMale) && !"".equals(tempFileFemale)) {
			centFiles.add(SEX_CENT);
		}
		
		for (String file : centFiles) {
			if (SEX_CENT.equals(file)) {
				namePathMap.put(SEX_CENT + " - Male", tempFileMale);
				namePathMap.put(SEX_CENT + " - Female", tempFileFemale);
			} else if (Files.exists(file)) {
				String name = file.substring(file.lastIndexOf("/") + 1);
				name = name.substring(0, name.lastIndexOf("."));
				namePathMap.put(name, file);
			}
		}
		
		centButtonMap = new HashMap<String, JCheckBoxMenuItem>();
		ButtonGroup centButtons = new ButtonGroup();
		JCheckBoxMenuItem centBox = new JCheckBoxMenuItem("None");
		centBox.addItemListener(centListener);
		centBox.setFont(new Font("Arial", 0, 12));
		centBox.setSelected(true);
		centButtons.add(centBox);
		
		JMenu lbl3 = new JMenu("Derive from Centroids");
		lbl3.setMnemonic(KeyEvent.VK_D);
		act.add(lbl3);
		lbl3.add(centBox);
		String[] centKeys = (String[]) namePathMap.keySet().toArray(new String[]{});
		for (String key : centKeys) {
			centBox = new JCheckBoxMenuItem(key);
			centBox.addItemListener(centListener);
			centBox.setFont(new Font("Arial", 0, 12));
			centButtons.add(centBox);
			lbl3.add(centBox);
			centButtonMap.put(key, centBox);
		}
		
		return menuBar;
	}
	
	public void actionPerformed(ActionEvent ae) {
		String command = ae.getActionCommand();
//		String[] filenames;

		if (command.equals(FIRST_CHR)) {
			parseLocation("chr1");
		} else if (command.equals(PREVIOUS_CHR)) {
			parseLocation("chr"+Math.max(chr-1, 1));
		} else if (command.equals(NEXT_CHR)) {
			parseLocation("chr"+Math.min(chr+1, 26));
		} else if (command.equals(LAST_CHR)) {
			parseLocation("chr26");
		} else if (command.equals(FIRST_REGION)) {
			if (regions == null || regions.length == 0) {
				JOptionPane.showMessageDialog(null, "Error - No regions have been loaded", "Error", JOptionPane.ERROR_MESSAGE);
				return;
			}
			regionIndex = 0;
			showRegion();
		} else if (command.equals(PREVIOUS_REGION)) {
			if (regions == null || regions.length == 0) {
				JOptionPane.showMessageDialog(null, "Error - No regions have been loaded", "Error", JOptionPane.ERROR_MESSAGE);
				return;
			}
			regionIndex = Math.max(regionIndex-1, 0);
			showRegion();
		} else if (command.equals(NEXT_REGION)) {
			System.out.println("next");
			if (regions == null || regions.length == 0) {
//				filenames = proj.getIndividualRegionLists();
//				if (filenames.length == 0) {
//					JOptionPane.showMessageDialog(null, "Error - No regions have been loaded, since there no individual CNV region files defined in the properties file", "Error", JOptionPane.ERROR_MESSAGE);
					JOptionPane.showMessageDialog(null, "Error - No regions have been loaded", "Error", JOptionPane.ERROR_MESSAGE);
//				} else {
//					JOptionPane.showMessageDialog(null, "Error - No regions have been loaded; files include: "+Array.toStr(filenames, ", "), "Error", JOptionPane.ERROR_MESSAGE);
//				}
				return;
			}
			regionIndex = Math.min(regionIndex+1, regions.length-1);
			showRegion();
		} else if (command.equals(LAST_REGION)) {
			if (regions == null || regions.length == 0) {
				JOptionPane.showMessageDialog(null, "Error - No regions have been loaded", "Error", JOptionPane.ERROR_MESSAGE);
				return;
			}
			regionIndex = regions.length-1;
			showRegion();
		} else if (command.equals(TO_SCATTER_PLOT)) {
			String[] listOfMarkers;
			
			listOfMarkers = new String[stopMarker-startMarker];
			if (listOfMarkers.length == 0) {
				JOptionPane.showMessageDialog(null, "There are no markers within this region; ScatterPlot will not bother launching", "Error", JOptionPane.ERROR_MESSAGE);
			} else {
				for (int i = startMarker; i < stopMarker; i++) {
					listOfMarkers[i-startMarker] = markerNames[i];
				}
				ScatterPlot.createAndShowGUI(proj, listOfMarkers, null, false);
			}
		} else {
			System.err.println("Error - unknown command '"+command+"'");
		}
	}

	public void mouseClicked(MouseEvent e) {
		if (e.getButton()==MouseEvent.BUTTON1) {
			if (leftClick!=null && leftClick.isAlive()) {
				leftClick.cancel();
				leftClick = null;
				zoomProportionally(false, e.getPoint(), true);
			} else {
				new Thread(leftClick = new SingleClick(this, MouseEvent.BUTTON1, DOUBLE_CLICK_INTERVAL)).start();
			}

		} else if (e.getButton()==MouseEvent.BUTTON3) {
			if (rightClick!=null && rightClick.isAlive()) {
				rightClick.cancel();
				rightClick = null;
				zoom(1, 1);
			} else {
				new Thread(rightClick = new SingleClick(this, MouseEvent.BUTTON3, DOUBLE_CLICK_INTERVAL)).start();
			}
		}
	}

	public void singleLeftClick() {
	// System.out.println("Single left click");
	}

	public void singleRightClick() {
	// System.out.println("Single right click");
	}

	public void mouseEntered(MouseEvent e) {}

	public void mouseExited(MouseEvent e) {}

	public void mousePressed(MouseEvent e) {
		startX = e.getPoint().x;
		inDrag = true;
	}

	public void mouseReleased(MouseEvent e) {
		inDrag = false;
	}

	public void mouseDragged(MouseEvent e) {
		int curX = e.getPoint().x;
		int distance = startX-curX;

		distance *= (stop-start)/(getWidth()-2*WIDTH_BUFFER);

		if (distance<0) {
			distance = Math.max(distance, 1-start);
		} else {
			distance = Math.min(distance, positions[chrBoundaries[chr][1]]-stop);
		}

		if ((start<=1&&distance<0)||(stop>=positions[chrBoundaries[chr][1]]&&distance>0)) {

		} else {
			start += distance;
			stop += distance;
		}

		if (inDrag) {
			updateGUI();
			startX = curX;
		}
	}

	public void mouseMoved(MouseEvent e) {}

	public void mouseWheelMoved(MouseWheelEvent e) {
		zoomProportionally(e.getWheelRotation()>0, e.getPoint(), false);
	}

	public void zoomProportionally(boolean outNotIn, Point p, boolean center) {
		int width = lrrPanel.getWidth()-2*WIDTH_BUFFER;
		double x = p.getX()-WIDTH_BUFFER;
		double xHat = width-x;
		double left = x/width;
		double right = xHat/width;
		double multiplier = MOUSE_WHEEL_MULTIPLIER/(outNotIn?1:-2);

		if (!outNotIn&&center) {
			right = 0.25-right;
			left = 0.25-left;
			multiplier = MOUSE_WHEEL_MULTIPLIER;
		}

		zoom(left*multiplier, right*multiplier);

	}

	public void zoom(double leftProportion, double rightProportion) {
		int dist = stop-start;
		start = start-(int)(leftProportion*dist);
		stop = stop+(int)(rightProportion*dist);
		updateGUI();
	}
	
	class MessageOfEncouragment implements Runnable {
		private String message;
		private Project proj;
		private boolean noLongerNecessary;

		public MessageOfEncouragment(String message, Project proj) {
			this.message = message;
			this.proj = proj;
			this.noLongerNecessary = false;
		}
		
		@Override
		public void run() {
			int count;
			
			count = 0;
			while (!noLongerNecessary && count < 6) {
				try {
					Thread.sleep(500);
				} catch (InterruptedException ie) {}
				count++;
			}
			
			if (!noLongerNecessary) {
				proj.message(message, "Patience...", JOptionPane.INFORMATION_MESSAGE);
			}
		}
		
		public void disregard() {
			noLongerNecessary = true;
		}
		
	}

	public void createSampleList() {
		long time;
		String[] filesPresent;
		FontMetrics fontMetrics;
		String refresh;
		int maxWidth;
		MessageOfEncouragment mess;

		time = new Date().getTime();
		log.report("  Getting a list of all files with extension "+Sample.SAMPLE_DATA_FILE_EXTENSION+" (if the process hangs here the first time after reverse transposing, please be patient, the operating system is busy indexing the new files) ...");
		mess = new MessageOfEncouragment("Getting a list of all sample files is taking longer than usual and probably means that your recently created files are still being indexed on the hard drive. Please be patient...", proj);
		new Thread(mess).start();
		filesPresent = Files.list(proj.SAMPLE_DIRECTORY.getValue(false, true), Sample.SAMPLE_DATA_FILE_EXTENSION, jar);
		log.report("Getting list of files took "+ext.getTimeElapsed(time));
		time = new Date().getTime();
		fontMetrics = sampleList.getFontMetrics(sampleList.getFont());
		refresh = "refresh list";
		maxWidth = fontMetrics.stringWidth(refresh);
		System.out.println("  computing font metrics took "+ext.getTimeElapsed(time));
		System.out.println("Determined sample list in "+ext.getTimeElapsed(time));
		mess.disregard();

		if (filesPresent == null || filesPresent.length == 0) {
//			samplesPresent = new String[] {proj.get(Project.SAMPLE_DIRECTORY)+" directory is empty", refresh};
			samplesPresent = new String[] {proj.SAMPLE_DIRECTORY.getValue(false, true)+" directory is empty", refresh};
			maxWidth = Math.max(maxWidth, fontMetrics.stringWidth(samplesPresent[0]));
		} else {
			samplesPresent = new String[filesPresent.length+1];
			for (int i = 0; i<filesPresent.length; i++) {
				samplesPresent[i] = filesPresent[i].substring(0, filesPresent[i].lastIndexOf("."));
				maxWidth = Math.max(maxWidth, fontMetrics.stringWidth(samplesPresent[i]));
            }
			samplesPresent[filesPresent.length] = refresh;
		}
		
		sampleList.setModel(new DefaultComboBoxModel<String>(samplesPresent));
		sampleList.setPreferredSize(new Dimension(maxWidth+50, 30));

	}

	public boolean loadMarkers() {
		Hashtable<String,String> hash;
		byte[] chrs;
		long time;
		int chr;

		time = new Date().getTime();

		hash = proj.getFilteredHash();
		markerSet = proj.getMarkerSet();
		if (markerSet == null) {
			JOptionPane.showMessageDialog(null, "Error - Failed to load the MarkerSet file; make sure the raw data is parsed", "Error", JOptionPane.ERROR_MESSAGE);
			log.reportError("Error - failed to load MarkerSet for project "+proj.getNameOfProject()+"; make sure the raw data is parsed");
			return false;
		}
		fingerprint = markerSet.getFingerprint();
		markerNames = markerSet.getMarkerNames();
		chrs = markerSet.getChrs();
		positions = markerSet.getPositions();

		dropped = new boolean[markerNames.length];
		chrBoundaries = new int[27][2];
		for (int i = 0; i<chrBoundaries.length; i++) {
//			chrBoundaries[i][0] = chrBoundaries[i][1] = -1;
			chrBoundaries[i][0] = chrBoundaries[i][1] = 0;
		}
		chr = 0;
		for (int i = 0; i<markerNames.length; i++) {
			dropped[i] = hash.containsKey(markerNames[i]);
			if (chrs[i]>chr) {
				if (chr!=0) {
					chrBoundaries[chr][1] = i-1;
				}
				chr = chrs[i];
				chrBoundaries[chr][0] = i;
			}
		}
		chrBoundaries[chr][1] = markerNames.length-1;
		chrBoundaries[0][0] = 0;
		chrBoundaries[0][1] = markerNames.length-1;

		System.out.println("Read in data for "+markerNames.length+" markers in "+ext.getTimeElapsed(time));
		return true;
	}

	public void loadValues() {
		Sample samp;
		
		samp = proj.getPartialSampleFromRandomAccessFile(sample, false, true, true, true, false);
		if (samp == null) {
			System.err.println("Error - sample '"+sample+"' not found in "+proj.SAMPLE_DIRECTORY.getValue(false, true));
		} else if ( samp.getFingerprint()!=fingerprint) {
			System.err.println("Error - Sample "+proj.SAMPLE_DIRECTORY.getValue(false, true)+sample+Sample.SAMPLE_DATA_FILE_EXTENSION+" has a different fingerprint ("+samp.getFingerprint()+") than the MarkerSet ("+fingerprint+")");
		} else {
			
			if (currentCentroid != null && currentCentroid.startsWith(SEX_CENT) && autoSwitch.isSelected()) {
				SampleData sampleData = proj.getSampleData(0, false);
				int sex = sampleData.getSexForIndividual(samp.getSampleName());
				if (sex == 1) {
					if (currentCentroid.endsWith("Female") && !isSettingCentroid) {
						centButtonMap.get(SEX_CENT + " - Male").setSelected(true);
						log.report("Switching to specified male centroid file");
					}
				} else if (sex == 2) {
					if (currentCentroid.endsWith("Male") && !isSettingCentroid) {
						centButtonMap.get(SEX_CENT + " - Female").setSelected(true);
						log.report("Switching to specified female centroid file");
					}
				} else {
					log.report("Warning - no sex specified for sample " + samp.getSampleName() + "; using currently selected centroid file");
				}
			}
			
			lrrs = samp.getLRRs();
			bafs = samp.getBAFs();
			genotypes = samp.getAB_Genotypes();
			
			if (transformation_type > 0) {
				lrrValues = Transforms.transform(lrrs, transformation_type, transformSeparatelyByChromosome, markerSet);
			} else if (transformation_type < 0 && centroids != null) {
				lrrValues = samp.getLRRs(centroids);
				originalBAFs = bafs;
				bafs = samp.getBAFs(centroids);
			} else {
				lrrValues = lrrs;
				originalBAFs = bafs;
			}
			
		}
		// lrrMin = Math.floor(lrrMin);
		// lrrMax = Math.ceil(lrrMax);
	}

	public void parseLocation(String location) {
		byte oldChr;
		int[] loc;

		oldChr = chr;

		if (location == null) {
			System.err.println("Error - null location");
			return;
		}
		
		if (!location.startsWith("chr")) {
			if (track == null) {
				JOptionPane.showMessageDialog(this, "Cannot parse '"+location+"' since the gene track has either not been installed or did not load properly.", "Error", JOptionPane.ERROR_MESSAGE);
				return;
			} else {
				loc = track.lookupPosition(location);
				if (loc[0] == -1) {
					JOptionPane.showMessageDialog(this, "'"+location+"' is not a valid gene name and is not a valid UCSC location.", "Error", JOptionPane.ERROR_MESSAGE);
					return;
				}
			}
		} else {
			loc = Positions.parseUCSClocation(location);
		}
		
		if (loc == null) {
			JOptionPane.showMessageDialog(this, "'"+location+"' is not a valid UCSC location.", "Error", JOptionPane.ERROR_MESSAGE);
			return;
		}

		chr = (byte)loc[0];
		if (chr==-1) {
			chr = oldChr;
			return;
		}
		if (chr!=oldChr) {
			start = stop = -1;
			procCNVs(chr);
		}
		start = loc[1];
		stop = loc[2];
		if (start==-1||start<0) {
			start = 1;
		}
		if (stop==-1||stop>positions[chrBoundaries[chr][1]]) {
			stop = positions[chrBoundaries[chr][1]];
		}
		
		updateGUI();
	}

	public void updateSample(String newSample) {
		boolean found;
		long time;

		found = sampleList.getSelectedItem().equals(newSample);
		if (found) {
			time = new Date().getTime();
//			System.out.println("Loading CNVs...");
//			System.out.println("took "+ext.getTimeElapsed(time));
			sample = newSample;
			time = new Date().getTime();
			System.out.print("Found "+sample+"...");
			indiPheno = sampleData.getIndiPheno(sample.toLowerCase());
			if (indiPheno == null) {
//				if (!sample.equals(proj.get(Project.SAMPLE_DIRECTORY)+" directory is empty")) {
				if (!sample.equals(proj.SAMPLE_DIRECTORY.getValue(false, true)+" directory is empty")) {
					JOptionPane.showMessageDialog(this, "Sample '"+sample+"' was not present in the SampleData file", "Error", JOptionPane.ERROR_MESSAGE);
				}
				return;
			}
			loadValues();
//			if (regionsList.length == 0 || !Files.exists(regionsList[regionsListIndex], jar)) {
//				loadCNVsAsRegions();
//			}
			if (REGION_LIST_USE_CNVS.equals(Trailer.this.regionFileName)) {
				loadCNVsAsRegions();
				procCNVs(chr);
				updateGUI();
				regionIndex = 0;
//				showRegion();
			} else {
				procCNVs(chr);
				updateGUI();
			}
			System.out.println("updated in "+ext.getTimeElapsed(time));
		} else {
			for (int i = 0; i<samplesPresent.length && !found; i++) {
				if (samplesPresent[i].equals(newSample)) {
					sampleList.setSelectedIndex(i);
					found = true;
				}
	        }
			if (!found) {
				if (Files.exists(proj.SAMPLE_DIRECTORY.getValue(false, true)+newSample+Sample.SAMPLE_DATA_FILE_EXTENSION, jar)) {
					createSampleList();
					updateSample(newSample);
				} else {
					JOptionPane.showMessageDialog(this, "Data was not found for the next sample in the list ("+newSample+").", "Error", JOptionPane.ERROR_MESSAGE);
				}
			}
		}
		
	}
	
	public void updateGUI() {
		if (start<=1) {
			startMarker = chrBoundaries[chr][0];
			start = 1;
		} else {
			startMarker = Array.binarySearch(positions, start, chrBoundaries[chr][0], chrBoundaries[chr][1], false);
		}

		if (stop>=positions[chrBoundaries[chr][1]]) {
			stop = positions[chrBoundaries[chr][1]];
			stopMarker = chrBoundaries[chr][1];
		} else {
			stopMarker = Array.binarySearch(positions, stop, chrBoundaries[chr][0], chrBoundaries[chr][1], false);
		}

		if (startMarker==-1) {
			System.err.println("Error - failed to find startMarker");
			startMarker = chrBoundaries[chr][0];
		}

		if (stopMarker==-1) {
			System.err.println("Error - failed to find stopMarker");
			stopMarker = chrBoundaries[chr][1];
		}
		
		displayIndex();
		repaint();
	}

	public void loadRegions() {
		BufferedReader reader;
        Vector<String[]> v;
        String[] line;
        int ignoredLines, countMissingRegions, invalidSamples;
		
		try {
//			reader = Files.getReader(regionsList[regionsListIndex], jar, false, false);
			String file = regionFileName.startsWith("./") ? proj.PROJECT_DIRECTORY.getValue() + regionFileName : regionFileName;
			reader = Files.getReader(file, jar, false, false);
//			System.out.print("Loading regions from "+regionsList[regionsListIndex]+"...");
			System.out.print("Loading regions from " + regionFileName + "...");
	        v = new Vector<String[]>();
	        ignoredLines = countMissingRegions = invalidSamples = 0;
            while (reader.ready()) {
            	line = reader.readLine().trim().split("\t");
            	if (sampleData.lookup(line[0]) == null) {
            		log.reportError("Error - '"+line[0]+"' is not a valid sample id");
            		invalidSamples++;
            	} else if (line.length == 1) {
            		v.add(new String[] {line[0], "chr1"});
            		countMissingRegions++;
            	} else if (line.length > 1 && line[1].startsWith("chr")) {
            		v.add(line);
            	} else {
            		ignoredLines++;
            	}
            }
            System.out.println(" loaded "+v.size()+" regions");
            regions = Matrix.toStringArrays(v);
            if (invalidSamples > 0) {
            	JOptionPane.showMessageDialog(null, "Error - there were "+invalidSamples+" invalid samples in '"+regionFileName+"' that were ignored because they could not be found", "Error", JOptionPane.ERROR_MESSAGE);
            }
            if (countMissingRegions > 0) {
            	JOptionPane.showMessageDialog(null, "Warning - there were "+countMissingRegions+" lines in '"+regionFileName+"' without a chromosomal region listed; using \"chr1\" for all missing values", "Warning", JOptionPane.ERROR_MESSAGE);
            }
            if (ignoredLines > 1) {
            	JOptionPane.showMessageDialog(null, "Error - there were "+ignoredLines+" regions in '"+regionFileName+"' that were ignored due to improper formatting", "Error", JOptionPane.ERROR_MESSAGE);
            }
            reader.close();
            
            
        } catch (FileNotFoundException fnfe) {
            System.err.println("Error: file \""+regionFileName+"\" not found in data directory");
            System.exit(1);
        } catch (IOException ioe) {
            System.err.println("Error reading file \""+regionFileName+"\"");
            System.exit(2);
        }
	}
	
	public void loadCNVsAsRegions() {
		int buffer;
		Segment[][] segs;
		int count;
		
		count = 0;
		segs = new Segment[26][];
		for (int i = 0; i<segs.length; i++) {
			procCNVs((byte)i);
			segs[i] = findUniqueRegions(cnvs);
			count += segs[i].length;
        }

		regions = new String[count][2];
		count = 0;
		for (int i = 0; i<segs.length; i++) {
			for (int j = 0; j<segs[i].length; j++) {
				regions[count][0] = sample;
				buffer = Math.max(segs[i][j].getSize() / 2, MIN_BUFFER);
				regions[count][1] = Positions.getUCSCformat(new int[] {segs[i][j].getChr(), segs[i][j].getStart()-buffer, segs[i][j].getStop()+buffer});
				count++;
            }
		}
	}

	public static Segment[] findUniqueRegions(CNVariant[][] cnvs) { // haven't actually coded the collapsing of the segments yet
		Segment[] segs;
		int count;
		
		count = 0;
		for (int i = 0; i<cnvs.length; i++) {
			count += cnvs[i].length;
        }
		
		segs = new Segment[count];
		count = 0;
		for (int i = 0; i<cnvs.length; i++) {
			for (int j = 0; j<cnvs[i].length; j++) {
				segs[count] = cnvs[i][j];
				count++;
            }
		}

		return segs;
	}

	public void showRegion() {
		if (regions == null || regions.length == 0) {
			regionField.setText("");
			commentLabel.setText("");
			return;
		}
//		System.out.println("regionIndex="+regionIndex+"\t"+"regions.length="+regions.length);
		parseLocation(regions[regionIndex][1]);
		if (regions[regionIndex].length > 2) {
			commentLabel.setText("region #"+(regionIndex+1)+":  "+ regions[regionIndex][2]);
		} else {
			commentLabel.setText(" -- no comment -- ");
		}
		
		if (!regions[regionIndex][0].equals(sample)) {
			updateSample(regions[regionIndex][0]);
		}

		regionField.setText((regionIndex + 1) + " of " + regions.length);
	}
	
//	public static Vector<String[]> loadFileToVec(String filename, boolean demo, boolean ignoreFirstLine, int[] cols, boolean onlyIfAbsent) {
//		BufferedReader reader = null;
//		Vector<String[]> v = new Vector<String[]>();
//		String[] trav;
//		String[] line;
//
//		try {
//			if (demo) {
//				reader = new BufferedReader(new InputStreamReader(ClassLoader.getSystemResourceAsStream(filename)));
//			} else {
//				reader = new BufferedReader(new FileReader(filename));
//			}
//			if (ignoreFirstLine) {
//				reader.readLine();
//			}
//			while (reader.ready()) {
//				trav = line = reader.readLine().trim().split("\t", -1);
//				if (cols!=null) {
//					trav = new String[cols.length];
//					for (int i = 0; i<cols.length; i++) {
//						trav[i] = line[cols[i]];
//					}
//				}
//				if (onlyIfAbsent) {
//					if (!v.contains(trav)) {
//						v.add(trav);
//					}
//				} else {
//					v.add(trav);
//				}
//			}
//			reader.close();
//		} catch (FileNotFoundException fnfe) {
//			System.err.println("Error: file \""+filename+"\" not found in current directory");
//			System.exit(1);
//		} catch (Exception e) {
//			System.err.println("Error reading file \""+filename+"\"");
//			e.printStackTrace();
//			System.exit(2);
//		}
//
//		return v;
//	}

	public void displayIndex() {
		if (start<0) {
			start = 1;
		}
		if (stop>positions[chrBoundaries[chr][1]]) {
			stop = positions[chrBoundaries[chr][1]];
		}

		navigationField.setText(chr==0?"all":"chr"+chr+":"+ext.addCommas(start)+"-"+ext.addCommas(stop));
	}

//	public static CNVariant[] getCNVlist(Hashtable<String,Hashtable<String,CNVariant[]>> cnvHashes, String[] filenames, String sample) {
//		Hashtable<String,CNVariant[]> hash;
//		Vector<CNVariant> v;
//		CNVariant[] set;
//		
//		v = new Vector<CNVariant>();
//		
//		for (int i = 0; i<filenames.length; i++) {
//			hash = cnvHashes.get(filenames[i]);
//			if (hash != null) {
//				set = hash.get(sample);
//				if (set != null) {
//					for (int j = 0; j<set.length; j++) {
//						set[j].setSource(i);
//						v.add(set[j]);
//                    }
//				}
//			}
//        }
//		
//		return CNVariant.toArray(v);
//	}
	
	public void procCNVs(byte chr) {
		cnvs = new CNVariant[cnvLabels.length][];
		for (int i = 0; i<cnvLabels.length; i++) {
			cnvs[i] = indiPheno.getCNVs(i, chr);
			if (cnvs[i] == null) {
				cnvs[i] = new CNVariant[0];
			}
			System.out.println("Proccessed "+cnvs[i].length+" "+cnvLabels[i]+" CNVs for chr"+chr);
        }
	}
	
	/**
	 * @param proj
	 * @param lrrsToTransform
	 *            the log R ratio input array
	 * @param transformation_type
	 *            transformation type for {@link Transforms#transform(float[], int, boolean, MarkerSet)}
	 * @param transformSeparatelyByChromosome
	 *            transform log R ratios separately
	 * @param markerSet
	 * @param gcModel
	 *            a gc model, if the gcmodel is null and correctGC is true, we report an error
	 * @param correctGC
	 *            whether to perform gc correction, must have a valid {@link GcModel} to correct
	 * @param correctGCFirst
	 *            perform the gc correction first, and then any transformations
	 * @param log
	 * @return
	 */
	private static float[] getNewLRRs(Project proj, float[] lrrsToTransform, int transformation_type, boolean transformSeparatelyByChromosome, MarkerSet markerSet, GcModel gcModel, boolean correctGC, boolean correctGCFirst, Logger log) {
		float[] tmpLrrs = lrrsToTransform; // make sure not to modify
		if (gcModel == null && correctGC) {
			log.reportError("Error - gc Correction was flagged and the model was null, this should not happen...skipping gc correction");
			correctGC = false;
		}
		if (correctGC && correctGCFirst) {
			tmpLrrs = Array.toFloatArray(GcAdjustor.getComputedAdjustor(proj, tmpLrrs, gcModel, false, false, false, true).getCorrectedIntensities());
		}

		if (transformation_type > 0) {
			tmpLrrs = Transforms.transform(tmpLrrs, transformation_type, transformSeparatelyByChromosome, markerSet);
		}
		if (correctGC && !correctGCFirst) {
			tmpLrrs = Array.toFloatArray(GcAdjustor.getComputedAdjustor(proj, tmpLrrs, gcModel, false, false, false, true).getCorrectedIntensities());
		}
		return tmpLrrs;
	}
	
	private void setCentroid(String name) {
//		String name = centroidsSelection.getSelectedItem().toString();
		if (!name.equals(currentCentroid)) {
			String path = namePathMap.get(name);
			Centroids cent = Centroids.load(path, false);
			centroids = cent.getCentroids();
			currentCentroid = name;
		}
	}
	
//	public static CNVariant[] loadCNVfiles(Project proj, String[] filenames) {
//		BufferedReader reader;
//		Vector<CNVariant> v = null;
//		String[] line;
//		String id = null;
//		boolean jar;
//
//		jar = proj.getJarStatus();
//		
//		v = new Vector<CNVariant>();
//		for (int i = 0; i<filenames.length; i++) {
//			try {
//				reader = Files.getReader(filenames[i], jar, true, false);
//				if (reader!=null) {
//					reader.mark(1000);
//					line = reader.readLine().trim().split("[\\s]+");
//					if (!line[2].toLowerCase().equals("chr")&&ext.chromosomeNumber(line[2])!=-1) {
//						reader.reset();
//					}
//					while (reader.ready()) {
//						line = reader.readLine().trim().split("[\\s]+");
//						if (line[1].equals(id)) {
//							v.add(new CNVariant(line, i));
//						}
//					}
//					reader.close();
//				}
//			} catch (IOException ioe) {
//				System.err.println("Error reading file \""+filenames[i]+"\"");
//				ioe.printStackTrace();
//			}
//
//		}
//
//		return CNVariant.sortCNVs(CNVariant.toArray(v));
//	}
	
	
	public static void main(String[] args) {
		Project proj;
		
		proj = new Project(cnv.Launch.getDefaultDebugProjectFile(true), false);
		new Trailer(proj, DEFAULT_SAMPLE, proj.CNV_FILENAMES.getValue(), DEFAULT_LOCATION);
	}
}


