package cnv.plots;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Map.Entry;

import javax.imageio.ImageIO;
import javax.swing.AbstractAction;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.DefaultComboBoxModel;
import javax.swing.DefaultListCellRenderer;
import javax.swing.JButton;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.ToolTipManager;
import javax.swing.border.EmptyBorder;

import net.miginfocom.swing.MigLayout;
import seq.manage.VCFOps;
import seq.manage.VCOps;
import cnv.filesys.MarkerSet.PreparedMarkerSet;
import cnv.filesys.Project;
import cnv.gui.NewRegionListDialog;
import cnv.var.Region;
import common.Aliases;
import common.Array;
import common.Files;
import common.Fonts;
import common.Grafik;
import common.Logger;
import common.Positions;
import common.TransferableImage;
import common.ext;
import filesys.GeneData;
import filesys.GeneTrack;
import filesys.Segment;

public class VariantViewer extends JFrame implements ActionListener, MouseListener, MouseMotionListener, MouseWheelListener {
	private static final String COLLAPSE_ISOFORMS_KEY = "Collapse Isoforms";

    public static final long serialVersionUID = 1L;

	public static final String DEFAULT_LOCATION = "chr17:55,609,472-55,824,368"; // USP32

	public static final String DEFAULT_SAMPLE = null;
	public static final boolean SHOW_MIDLINE = true;

	public static final double MOUSE_WHEEL_MULTIPLIER = 0.5;
	public static final int WIDTH_BUFFER = 25;
	public static final int HEIGHT_BUFFER = 10;
	public static final int DYNAMIC_HEIGHT_LIMIT = 0;
	public static final int DOUBLE_CLICK_INTERVAL = 500;
	public static final int MIN_BUFFER = 1500;
	public static final int DEFAULT_STARTX = 20;
	public static final int DEFAULT_STARTY = 20;
	public static final int DEFAULT_WIDTH = 1200;
	public static final int DEFAULT_HEIGHT = 600;

	private static final String FIRST_REGION = "First region";
	private static final String PREVIOUS_REGION = "Previous region";
	private static final String NEXT_REGION = "Next region";
	private static final String LAST_REGION = "Last region";
	private static final String TO_TRAILER_PLOT = "To Trailer";
	private static final String TO_SCATTER_PLOT = "To Scatter Plot";
	private static final String TO_COMP_PLOT = "To Comp Plot";
	private static final String REGION_LIST_NEW_FILE = "Load Region File";
	private static final String REGION_LIST_PLACEHOLDER = "Select Region File...";

    protected static final int COLLAPSE_ISOFORMS = -9;
    protected static final int EXON_LBL_LEFT = 1;
    protected static final int EXON_LBL_CENTER = 2;
    protected static final int EXON_LBL_RIGHT = 3;

    private static final int GENE_HEIGHT = 30;
    private static final int EQUALIZED_EXON_BP_LENGTH = 300;
    private static final int INTRON_PLACEHOLDER_BP_LENGTH = 25;
    private static final String EXON_PREFIX = "";
    private static final int Y_START = 2*15;
    private static final Color FILLED_EXON_COLOR = Color.GRAY;
    private static final int DATA_PNT_SIZE = 8;

    private volatile int intronBPWidth = INTRON_PLACEHOLDER_BP_LENGTH;
    private volatile boolean equalizeExonLength = false;
    private volatile boolean paintExonNumbers = false;
    private volatile boolean paintExonBoundaries = true;
    private volatile boolean paintIntrons = false;
    private volatile boolean fillExons = true;
    private volatile boolean paintExonBoxes = true;
    private volatile boolean paintInternalLine = false;
    private volatile int exonLabelLocation = EXON_LBL_LEFT;
    private volatile boolean showExcludes = false;
    private volatile boolean showLegend = true;
    private volatile boolean drawMAF = false;
    private volatile boolean drawMAC = false;
    
    private volatile int dataPntSize = DATA_PNT_SIZE;
    private volatile int yStart = Y_START;
    
	private JComboBox<String> isoformList;
	private String[] isoformsPresent;
	private JButton previousGene, nextGene;
	private Project proj;
	private boolean jar;
	private int[] positions;
	private boolean[] dropped;
	private int[][] chrBoundaries;
	private byte chr;
	private int start, startMarker;
	private int stop, stopMarker;
	private boolean inDrag;
	private int startX;
	private int geneIndex;
//	private int isoformIndex;
//	private JPanel lrrPanel;
//	private JPanel bafPanel;
	private JPanel genePanel;
	private GeneTrack track;
	HashMap<String, String> geneToCommentMap;
	HashMap<String, HashMap<String, String>> geneToRegionMap;
	HashMap<String, HashMap<String, Segment[]>> geneToExonSegmentMap;
	HashMap<String, HashMap<String, GeneData>> geneToIsoformMap;
	HashMap<String, HashMap<String, ArrayList<ArrayList<VariantContextWithFile>>>> loadedVCFData;
	VCFHeader vcfHeader;
	HashSet<String> popSet;
	HashSet<String> excluded;
	HashMap<String, String> popMap;
    HashMap<String, String> superPopMap;
//    HashMap<String, String> popColorCodeMap;
//    HashMap<String, Color> colorCodeMap;
    HashMap<String, Color> popColorMap;
    HashMap<String, HashSet<String>> popIndiMap;
    HashMap<String, HashSet<String>> popIndiMapWithExcludes;
	HashMap<String, VCFHeader> headerMap;
	private JLabel commentLabel;
	private JTextField commentField;
	private String isoform;
	private ArrayList<String> geneList;
	
	private String[] vcfFiles;
	private String popFile;
	
	private Hashtable<String, String> namePathMap;
	private Logger log;
	private boolean fail;
	private JMenu loadRecentFileMenu;
	private ButtonGroup regionButtonGroup;
    private JCheckBoxMenuItem chkbxDisplayExcludes;
	
	private static final int DRAW_AS_CIRCLES = 1;
	private static final int DRAW_AS_BLOCKS = 2;
	
	private int drawType = DRAW_AS_BLOCKS;
	
	private enum DrawType { 
	    FILLED_CIRCLE,
	    EMPTY_CIRCLE,
	    X;
	    
	    public static DrawType getDrawType(VariantContext vc) {
	        String impAttr = vc.getAttribute("SNPEFF_IMPACT").toString();

	        if ("LOW".equals(impAttr)) {
	            return DrawType.EMPTY_CIRCLE;
	        } else if ("MODERATE".equals(impAttr)) {
	            return DrawType.FILLED_CIRCLE;
	        } else if ("HIGH".equals(impAttr)) {
	            return DrawType.X;
	        } else {
	            return null;
	        }
	    }
	}
	private static class BlockDraw {
	    public BlockDraw(int basePairLoc, int xPixel, int numGenotypes, int totalAffected, HashMap<String, Integer> popCnts, DrawType drawType, VariantContextWithFile vc) {
	        this.bpX = basePairLoc;
	        this.x = xPixel;
	        this.g = numGenotypes;
	        this.aff = totalAffected;
	        this.gen = popCnts;
	        this.dt = drawType;
	        this.vcRecord = vc;
        }
	    int bpX;
	    int x;
	    int g;
	    int aff;
	    HashMap<String, Integer> gen;
	    DrawType dt;
	    VariantContextWithFile vcRecord;
	}
	private static class DrawPoint {
	    public DrawPoint(int x2, int y2, int height, int width, DrawType drawType, Color color, String sampleID, VariantContextWithFile vc) {
            this.x = x2;
            this.y = y2;
            this.h = height;
            this.w = width;
            this.type = drawType;
            this.c = color;
            this.vcRecord = vc;
            this.sampleID = sampleID;
        }
        int x;
	    int y;
	    int h;
	    int w;
	    DrawType type;
	    Color c;
	    VariantContextWithFile vcRecord;
	    String sampleID;
	}
	
    private static class VCFLocation {
        public VCFLocation(int x, VariantContext vc) {
            this.x = x;
            this.vc = vc;
        }
        int x;
        VariantContext vc;
        HashMap<String, Double> mafMap = new HashMap<String, Double>();
        HashMap<String, Double> macMap = new HashMap<String, Double>();
    }
    
	private AbstractAction geneFileSelectAction = new AbstractAction() {
        private static final long serialVersionUID = 1L;
        @Override
        public void actionPerformed(ActionEvent e) {
            String shortName = ((JCheckBoxMenuItem)e.getSource()).getText();
            if (!loadingFile 
                    && !REGION_LIST_NEW_FILE.equals(shortName) 
                    && !REGION_LIST_PLACEHOLDER.equals(shortName)) {
                String file = regionFileNameLoc.get(shortName);
                if (file != null && file.equals(VariantViewer.this.geneFileName)) {
                    return;
                }
                String tempFile = file.startsWith("./") ? proj.PROJECT_DIRECTORY.getValue() + file : file;
                if (!Files.exists(tempFile)) {
                    proj.message("Error - region file '" + shortName + "' doesn't exist.");
                    regionFileNameBtn.get(shortName).setSelected(true);
                } else {
                    VariantViewer.this.geneFileName = file;
                    loadGenes(VariantViewer.this.geneFileName);
                    showGene(0);
                }
            } /*else if (loadingFile && REGION_LIST_PLACEHOLDER.equals(shortName)) {
                // do nothing
            } */else if (loadingFile || REGION_LIST_PLACEHOLDER.equals(shortName)) {
                // leave as currently selected marker
                if (VariantViewer.this.geneFileName != "" && VariantViewer.this.geneFileName != null) {
                    String file = VariantViewer.this.geneFileName;
                    file = ext.rootOf(VariantViewer.this.geneFileName);
                    regionFileNameBtn.get(file).setSelected(true);
                }
                return;
            } 
        }
    };
	
    AbstractAction loadGeneListFileAction = new AbstractAction() {
        private static final long serialVersionUID = 1L;
        @Override
        public void actionPerformed(ActionEvent e) {
            String newFile = chooseNewFiles();
            if (newFile == null) {
                if (VariantViewer.this.geneFileName != null && !"".equals(VariantViewer.this.geneFileName)) {
                    regionFileNameBtn.get(ext.rootOf(VariantViewer.this.geneFileName)).setSelected(true);
                }
            } else {
                String file = ext.verifyDirFormat(newFile);
                file = file.substring(0, file.length() - 1);
                String name = ext.rootOf(file);
                regionFileNameBtn.get(name).setSelected(true);
                regionFileNameBtn.get(name).doClick();
            }
        }
    };
	
    AbstractAction screencapAction = new AbstractAction() {
        private static final long serialVersionUID = 1L;
        @Override
        public void actionPerformed(ActionEvent e) {
            JFileChooser jfc = new JFileChooser(proj != null ? proj.PROJECT_DIRECTORY.getValue() : ".");
            jfc.setMultiSelectionEnabled(false);
            jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
            jfc.setDialogTitle("Save Screen Capture...");
            jfc.setDialogType(JFileChooser.SAVE_DIALOG);
            int code = jfc.showSaveDialog(VariantViewer.this);
            if (code == JFileChooser.APPROVE_OPTION) {
                String filename = jfc.getSelectedFile().getAbsolutePath();
                doScreenCapture(filename);
            }
        }
    };
    
    AbstractAction screencapClipboardAction = new AbstractAction() {
        private static final long serialVersionUID = 1L;
        @Override
        public void actionPerformed(ActionEvent e) {
            doScreenCapture(null);
        }
    };

	private ArrayList<Color[]> colorScheme = Trailer.getColor();

	private HashMap<String, JCheckBoxMenuItem> regionFileNameBtn = new HashMap<String, JCheckBoxMenuItem>();
	private HashMap<String, String> regionFileNameLoc = new HashMap<String, String>();
	private String geneFileName;
	private volatile boolean loadingFile = false;
	private JComboBox<String> geneListCmb;

    private PreparedMarkerSet markerSet;
    private String[] markerNames;
	
	public VariantViewer(Project proj, String[] vcfFiles, String popFile) {
		this(proj, proj.GENE_LIST_FILENAMES.getValue()[0], vcfFiles, popFile);
	}

	// TODO TrailerClone should have a createAndShowGUI, same as all the other plots, as opposed to being its own frame 
	public VariantViewer(Project proj, String geneListFile, String[] vcfFiles, String popFile) {
		super("Genvisis - Trailer - " + proj.PROJECT_NAME.getValue());
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		addWindowListener(new WindowAdapter() {
			@Override
			public void windowClosing(WindowEvent e) {
				if (VariantViewer.this.proj != null) {
					ArrayList<String> files = new ArrayList<String>(regionFileNameLoc.values());
					String[] currSet = VariantViewer.this.proj.GENE_LIST_FILENAMES.getValue();
					
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
				        int choice = JOptionPane.showOptionDialog(null, message+" Would you like to keep this configuration for the next time TrailerClose is loaded?", "Preserve TrailerClone workspace?", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, null, null);
				        if (choice == 0) {
				            VariantViewer.this.proj.GENE_LIST_FILENAMES.setValue(newList);
				            VariantViewer.this.proj.saveProperties();
				        }
					}
				}
				super.windowClosing(e);
			}
		});

		long time;
		String trackFilename;

		this.proj = proj;
		this.log = proj.getLog();
		jar = proj.JAR_STATUS.getValue();
		fail = false;
		this.vcfFiles = vcfFiles;
		this.popFile = popFile;
		
		time = new Date().getTime();

		fail = !loadMarkers();
		if (fail) {
			return;
		}        
//		trackFilename = proj.getGeneTrackFilename(false);
		trackFilename = "N:/statgen/NCBI/fullIsoforms/RefSeq.gtrack";
        if (trackFilename != null) {
            log.report("Loading track from "+trackFilename);    
            track = GeneTrack.load(trackFilename, jar);
            log.report("Loaded track in "+ext.getTimeElapsed(time));
        } else {
            log.report("Cannot create display without GeneTrack");
            return;
        }
        
		generateComponents();
		this.setJMenuBar(createMenuBar());
		
		time = new Date().getTime();

		setBounds(startX, DEFAULT_STARTX, DEFAULT_WIDTH, DEFAULT_HEIGHT);
		setVisible(true);
		
		updateGUI();
		loadGenes(geneListFile);
        createIsoformList();
		geneIndex = 0;
        try {
            loadPopulationFile();
        } catch (IOException e) {
            System.err.println("Error - problem occurred when loading population file: " + e.getMessage());
            e.printStackTrace();
            return;
        }
		showGene(0);
	}
	
	private String[][] readFile(String file) {
	    BufferedReader reader;
	    String line;
	    ArrayList<String[]> data;
	    
	    data = new ArrayList<String[]>();
	    try {
	        reader = Files.getAppropriateReader(file);
            while ((line = reader.readLine()) != null) {
                data.add(line.trim().split("\t"));
            }
        } catch (FileNotFoundException fnfe) {
            System.err.println("Error: file \"" + file + "\" not found in data directory");
        } catch (IOException ioe) {
            System.err.println("Error reading file \"" + file + "\"");
        }
	    
	    return data.toArray(new String[data.size()][]);
	}

	private void loadGenes(String file) {
		if (!Files.exists(file)) {// TODO, can put somewhere else
			proj.getLog().reportTimeWarning("Generating " + file + " using all genes in " + vcfFiles[0]);
			VCFOps.dumpSnpEffGenes(vcfFiles[0], file, proj.getLog());
		}
		proj.getLog().reportTimeWarning("Loading " + file );

		String[][] geneFile = readFile(file);
	    geneList = new ArrayList<String>();
	    geneToIsoformMap = new HashMap<String, HashMap<String, GeneData>>();
	    geneToRegionMap = new HashMap<String, HashMap<String, String>>();
	    geneToExonSegmentMap = new HashMap<String, HashMap<String, Segment[]>>();
	    geneToCommentMap = new HashMap<String, String>();
	    loadedVCFData = new HashMap<String, HashMap<String,ArrayList<ArrayList<VariantContextWithFile>>>>();
	    headerMap = new HashMap<String, VCFHeader>();
	    String[] genes = Array.extract(geneFile, 0);
	    GeneData[][] geneData = track.lookupAllGeneData(genes);
	    
	    for (int i = 0; i < geneFile.length; i++) {
	        // first, check if genes are on multiple chrs:
	        boolean multiChr = false;
	        int chr = -1;
	        for (int g = 0; g < geneData[i].length; g++) {
	            if (chr == -1) {
	                chr = geneData[i][g].getChr();
	            }
	            if (chr != geneData[i][g].getChr()) {
	                multiChr = true;
	                break;
	            }
	        }
	        if (multiChr) {
	            HashMap<Integer, HashMap<String, GeneData>> chrMap = new HashMap<Integer, HashMap<String,GeneData>>();
	            for (int g = 0; g < geneData[i].length; g++) {
                    HashMap<String, GeneData> isoMap = chrMap.get(geneData[i][g].getChr());
                    if (isoMap == null) {
                        isoMap = new HashMap<String, GeneData>();
                        chrMap.put((int) geneData[i][g].getChr(), isoMap);
                    }
                    isoMap.put(geneData[i][g].isCollapsedIsoforms() ? COLLAPSE_ISOFORMS_KEY : geneData[i][g].getNcbiAssessionNumbers()[0], geneData[i][g]);
                }
	            for (Entry<Integer, HashMap<String, GeneData>> chrEntry : chrMap.entrySet()) {
	                geneToIsoformMap.put(genes[i].toUpperCase() + " - chr" + chrEntry.getKey().intValue(), chrEntry.getValue());
	                geneList.add(genes[i].toUpperCase() + " - chr" + chrEntry.getKey().intValue());
	            }
	            
	        } else {
	            HashMap<String, GeneData> isoformMap = new HashMap<String, GeneData>();
	            geneToIsoformMap.put(genes[i].toUpperCase(), isoformMap);
	            geneList.add(genes[i].toUpperCase());
    	        for (int g = 0; g < geneData[i].length; g++) {
    	            if (geneData[i][g].isCollapsedIsoforms()) {
    	                isoformMap.put(COLLAPSE_ISOFORMS_KEY, geneData[i][g]);
    	            } else {
    	                isoformMap.put(geneData[i][g].getNcbiAssessionNumbers()[0], geneData[i][g]);
    	            }
    	        }
	        }
	        if (geneFile[i].length >= 2) {
	            geneToCommentMap.put(geneFile[i][0], geneFile[i][1]);
	        }
	    }
	    geneToRegionMap = new HashMap<String, HashMap<String, String>>();
	    geneToExonSegmentMap = new HashMap<String, HashMap<String,Segment[]>>();
	    for (String gene : geneList) {
	        HashMap<String, GeneData> isoMap = geneToIsoformMap.get(gene);
	        HashMap<String, String> isoPosMap = new HashMap<String, String>();
	        HashMap<String, Segment[]> isoSegMap = new HashMap<String, Segment[]>();
	        geneToRegionMap.put(gene, isoPosMap);
	        geneToExonSegmentMap.put(gene, isoSegMap);
	        for (Entry<String, GeneData> isoEntry : isoMap.entrySet()) {
	            GeneData value = isoEntry.getValue();
                isoPosMap.put(isoEntry.getKey(), "chr" + value.getChr() + ":" + value.getStart() + "-" + value.getStop());                
                int[][] exons = value.getExonBoundaries();
	            ArrayList<Segment> segList = new ArrayList<Segment>();
	            for (int[] i : exons) {
	                segList.add(new Segment(value.getChr(), i[0], i[1]));
	            }
	            isoSegMap.put(isoEntry.getKey(), segList.toArray(new Segment[segList.size()]));
	        }
	    }
        VariantViewer.this.geneFileName = file;
        updateGeneList();
	}
	
	private void updateGeneList() {
	    FontMetrics fontMetrics;
        int maxWidth;

        fontMetrics = geneListCmb.getFontMetrics(geneListCmb.getFont());
        maxWidth = fontMetrics.stringWidth("----------");
        
        String[] geneNames = new String[geneList.size()];
        for (int i = 0; i < geneList.size(); i++) {
            geneNames[i] = geneList.get(i);
            maxWidth = Math.max(maxWidth, fontMetrics.stringWidth(geneNames[i]));
        }

        geneListCmb.setModel(new DefaultComboBoxModel<String>(geneNames));
        geneListCmb.setPreferredSize(new Dimension(maxWidth + 50, 30));
	}
	
	private ArrayList<VariantContextWithFile> filter(Segment exon, ArrayList<VariantContextWithFile> all) {
	    ArrayList<VariantContextWithFile> retArr = new ArrayList<VariantContextWithFile>();
	    for (VariantContextWithFile vc : all) {
	        if (exon.overlaps(new Segment(vc.vc.getContig(), vc.vc.getStart(), vc.vc.getEnd()))) {
	            retArr.add(vc);
	        }
	    }
	    return retArr;
	}
	
	private GeneData getCurrentGeneData() {
	    if (geneToIsoformMap == null || geneList == null || geneList.isEmpty() || isoformList == null) return null;
	    return geneToIsoformMap.get(geneList.get(geneIndex)).get(isoformList.getSelectedItem());
	}
	
	private void paintPanel(Graphics g) {
		GeneData gene;
		ArrayList<VariantContextWithFile> vcfInSeg;
		int[][] exons;
		int width, begin, tempX, tempPx, len, lenPx, height;
		
		ArrayList<VCFLocation> freqLocs = new ArrayList<VariantViewer.VCFLocation>();
		
		boolean antiAlias = true;
        if (g instanceof Graphics2D) {
            ((Graphics2D)g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, antiAlias ? RenderingHints.VALUE_ANTIALIAS_ON : RenderingHints.VALUE_ANTIALIAS_OFF);
            ((Graphics2D)g).setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, antiAlias ? RenderingHints.VALUE_TEXT_ANTIALIAS_ON : RenderingHints.VALUE_TEXT_ANTIALIAS_OFF);
        }
        
		if (stop-start > 10000000) {
			g.drawString("Zoom in to see genes", 10, 10);
		} else {
		    height = genePanel.getHeight() - yStart - GENE_HEIGHT;
			gene = getCurrentGeneData();
			if (gene == null) {
			    return;
			}
			g.setColor(Color.BLACK);
			
            g.setFont(new Font("Arial", 0, 12));
	        tempX = begin = gene.getStart();
	        exons = gene.getExonBoundaries();
            activeBlocks.clear();
            activePoints.clear();
            activeRects.clear();
            drawnFreqs.clear();
            for (int j = 0; j < exons.length; j++) {
                tempPx = getX(tempX);
                len = equalizeExonLength ? EQUALIZED_EXON_BP_LENGTH : exons[j][1] - exons[j][0];
                lenPx = getX(tempX + len) - tempPx;
                // vertical line: 
                if (paintExonBoundaries) {
                    g.fillRect(tempPx, height, 1, GENE_HEIGHT);
                }
                if (fillExons) {
                    g.setColor(FILLED_EXON_COLOR);
                    g.fillRect(tempPx, height, lenPx, GENE_HEIGHT);
                    g.setColor(Color.BLACK);
                }
                if (paintExonBoxes) {
                    // horiz. line: for individual lines
                    g.drawRect(tempPx, height, lenPx, GENE_HEIGHT);
                }
                if (paintExonNumbers) {
                    int exonNumber = determineExonNumber(geneToIsoformMap.get(geneList.get(geneIndex)).get(COLLAPSE_ISOFORMS_KEY), exons[j]);
                    width = g.getFontMetrics().stringWidth(EXON_PREFIX + exonNumber);
                    if (width < lenPx - 2) {
                        if (fillExons) {
                            g.setColor(Color.WHITE);
                        }
                        int exonX = tempPx + 2;
                        if (exonLabelLocation == EXON_LBL_CENTER) {
                            exonX = (int)(tempPx + (lenPx / 2d)) - (int)(width / 2d);
                        } else if (exonLabelLocation == EXON_LBL_RIGHT) {
                            exonX = getX(tempX + len) - width - 2; 
                        }
                        g.drawString(EXON_PREFIX + exonNumber, exonX, height + GENE_HEIGHT / 2 + g.getFontMetrics().getHeight() - 2);
                        if (fillExons) {
                            g.setColor(Color.BLACK);
                        }
                    }
                }
                vcfInSeg = getExonVCFRecords(j);
                if (vcfInSeg.size() > 0) {
                    if (lenPx <= dataPntSize + 2) {
                        g.setColor(Color.RED);
                        g.drawLine(tempPx + (lenPx / 2), height - 30, tempPx + (lenPx / 2), height - 5);
                        g.setColor(Color.BLACK);
                    } else {
                        String selIso = isoformList.getSelectedItem().toString();
                        if (drawType == DRAW_AS_BLOCKS) {
                            // draw populations
                            ArrayList<BlockDraw> toDraw = new ArrayList<VariantViewer.BlockDraw>();
                            for (VariantContextWithFile vc : vcfInSeg) {
                                String isoAttr = vc.vc.getAttribute("SNPEFF_TRANSCRIPT_ID").toString();
                                // only draw if collapsed, or showing affected isoform
                                if (!selIso.equals(COLLAPSE_ISOFORMS_KEY) && !(selIso.equals(isoAttr) || selIso.equals(isoAttr.split("\\.")[0]))) {
                                    continue;
                                }
                                
                                GenotypesContext gctx = vc.vc.getGenotypes();
                                HashMap<String, Integer> popGenoCnt = new HashMap<String, Integer>();
                                int totAff = 0;
                                for (int i = 0; i < gctx.size(); i++) {
                                    Genotype geno = gctx.get(i);
                                    if (excluded.contains(geno.getSampleName()) && !showExcludes) {
                                        continue;
                                    }
                                    if (geno.getType() != GenotypeType.HOM_REF && geno.getType() != GenotypeType.NO_CALL) { // anything besides Homozygous Reference
                                        totAff++;
                                        String pop = popMap.get(geno.getSampleName());
                                        if (pop == null) {
                                            System.err.println("Error - no population entry found for ID: " + geno.getSampleName());
                                        }
                                        Integer cnt = popGenoCnt.get(pop);
                                        if (cnt == null) {
                                            cnt = Integer.valueOf(0);
                                        }
                                        popGenoCnt.put(pop, cnt + 1);
                                    }
                                }
                                
                                DrawType vcType = DrawType.getDrawType(vc.vc);
                                
                                int diffA = vc.vc.getStart() - exons[j][0];
                                double prop = diffA / (double) len;
                                int diffPx = (int) (lenPx * prop);
                                if (diffPx + dataPntSize + 2 > lenPx) { // don't let the markings go past the exon
                                    diffPx = lenPx - dataPntSize - 2;
                                }

                                freqLocs.add(new VCFLocation(tempPx + diffPx, vc.vc));
                                BlockDraw bd = new BlockDraw(vc.vc.getStart(), tempPx + diffPx, gctx.size(), totAff, popGenoCnt, vcType, vc);
                                toDraw.add(bd);
                            }
                            for (BlockDraw bd : toDraw) {
                                drawBlock(g, bd);
                            }
                        } else if (drawType == DRAW_AS_CIRCLES) {
                            // draw all in relative position, pushing up if overlapping
                            ArrayList<DrawPoint> drawn = new ArrayList<DrawPoint>();
                            ArrayList<Rectangle> plotted = new ArrayList<Rectangle>();
    //                        HashMap<String, Integer> plottedCnt = new HashMap<String, Integer>();
                            for (VariantContextWithFile vc : vcfInSeg) {
                                String isoAttr = vc.vc.getAttribute("SNPEFF_TRANSCRIPT_ID").toString();
                                // only draw if collapsed, or showing affected isoform
                                if (!selIso.equals(COLLAPSE_ISOFORMS_KEY) && !(selIso.equals(isoAttr) || selIso.equals(isoAttr.split("\\.")[0]))) {
                                    continue;
                                }
                                GenotypesContext gctx = vc.vc.getGenotypes();
                                int xOffset = 0;
                                int diffA = vc.vc.getStart() - exons[j][0];
                                double prop = diffA / (double)len;
                                int diffPx = (int) (lenPx * prop);
                                if (diffPx + dataPntSize + 2 > lenPx) { // don't let the markings go past the exon
                                    diffPx = lenPx - dataPntSize - 2;
                                }
                                Rectangle vcRect = new Rectangle(tempPx + diffPx + xOffset, 0, dataPntSize + 2, dataPntSize + 2);
                                int genoCnt = 0;
                                for (int i = 0; i < gctx.size(); i++) {
                                    Genotype geno = gctx.get(i);
                                    if (excluded.contains(geno.getSampleName()) && !showExcludes) {
                                        continue;
                                    }
                                    if (geno.getType() != GenotypeType.HOM_REF && geno.getType() != GenotypeType.NO_CALL) { // anything besides Homozygous Reference
                                        genoCnt++;
                                        boolean overlap = false;
                                        do {
                                            overlap = false;
                                            for (Rectangle rect : plotted) {
                                                if (rect.intersects(vcRect)) {
                                                    overlap = true;
                                                    break;
                                                }
                                            }
                                            if (overlap) {
                                                vcRect = new Rectangle(vcRect.x, vcRect.y + vcRect.height, dataPntSize + 2, dataPntSize + 2);
    //                                            if (vcRect.y + vcRect.height + 110 >= getHeight()) {
    //                                                vcRect = new Rectangle(tempPx + diffPx + (xOffset += dataPntSize), 0, dataPntSize + 2, dataPntSize + 2);
    //                                            }
                                            }
                                        } while (overlap && vcRect.y + vcRect.height + 110 < getHeight()/* && xOffset <= 20*/); // running off top of screen
                                        plotted.add(vcRect);
                                        DrawPoint dp = new DrawPoint(vcRect.x, vcRect.y, vcRect.height, vcRect.width, DrawType.getDrawType(vc.vc), popColorMap.get(popMap.get(geno.getSampleName())), geno.getSampleName(), vc);
                                        activePoints.add(dp);
                                        activeRects.add(vcRect);
                                        if (selectedDrawPoint != null && dp.sampleID.equals(selectedDrawPoint.sampleID)) {
                                            selectedRect = vcRect;
                                        }
                                        drawn.add(dp);
                                        if (vcRect.y + vcRect.height + 110 >= getHeight()) {
                                            break; // stop doing things if we're off the screen
                                        }
                                    }
                                }
                                if (genoCnt > 0) {
                                    freqLocs.add(new VCFLocation(vcRect.x, vc.vc));
                                }
                            }
                            for (int i = 0; i < drawn.size(); i++) {
                                drawVcfEntry(g, drawn.get(i));
                            }
                        }
                        g.setColor(Color.BLACK);
                    }
                }
                // move tempX to other side of exon:
                tempX += len;
                if (!fillExons && !paintExonBoxes) {
                    if (paintExonBoundaries) {
                        g.fillRect(getX(tempX)-1, height, 1, GENE_HEIGHT);
                    }
                }
                if (intronBPWidth > 0 || paintIntrons) {
                    if (j < exons.length - 1) {
                        int intronLen = paintIntrons ? exons[j+1][0] - exons[j][1] : intronBPWidth + 1;
                        if (!paintInternalLine) {
                            g.fillRect(getX(tempX), height + GENE_HEIGHT / 2 - 1, getX(tempX + intronLen) - getX(tempX), 2);
                        }
                        tempX += intronLen;
                    }
                }
            }
            if (paintInternalLine) {
                // horiz. line connecting all exons 
                g.fillRect(getX(begin), height + GENE_HEIGHT / 2 - 1, getX(tempX) - getX(begin), 2);
            }

			g.setFont(new Font("Arial", 0, 14));
            begin = getX(gene.getStart());
            width = g.getFontMetrics().stringWidth(gene.getGeneName());
			g.drawString(gene.getGeneName(), begin-width-3, height + GENE_HEIGHT / 2 + g.getFontMetrics().getHeight() / 2 - 3);
		    drawFreqs(g, freqLocs, begin);
		}
		if (showLegend) {
		    drawLegend(g);
		}
	}
	
	private void drawBlock(Graphics g, BlockDraw bd) {
        int height = genePanel.getHeight() - yStart - GENE_HEIGHT - 4;
        int width = dataPntSize;
        int x = bd.x;
        int tempY = 2;
        int scav = 0;
        HashMap<String, Integer> drawPop = new HashMap<String, Integer>();
        int drawMax = 0;
        activeRects.add(new Rectangle(x, 2, width, height));
        activeBlocks.add(bd);
        if (selectedBlockDraw != null && bd.bpX == selectedBlockDraw.bpX) {
            selectedRect = activeRects.get(activeRects.size() - 1);
        }
        int popCnt = 0;
        for (String pop : popSet) {
            Integer propInt = bd.gen.get(pop);
            if (propInt == null) {
                continue;
            }
            popCnt++;
        }
        int index = 0;
        int drawCumu = 0;
        for (String pop : popSet) {
            Integer propInt = bd.gen.get(pop);
            if (propInt == null) {
                continue;
            }            
            index++;
            double prop = propInt.doubleValue() / bd.aff;
            int minDraw = 2;
            int drawLen = (int) (height * prop);
            if (drawLen == 0) {
                scav += minDraw;
                drawLen = minDraw;
            } else if (bd.dt == DrawType.EMPTY_CIRCLE || bd.dt == DrawType.X) {
                drawLen -= 1;
            }
            drawCumu += drawLen;
            if (popCnt > 1 && index == popCnt) {
                if (height - drawCumu > 0 || height - drawCumu < 0) {
                    drawLen += height - drawCumu;
                    if (bd.dt == DrawType.EMPTY_CIRCLE || bd.dt == DrawType.X) {
                        drawLen -= 1;
                    }
                }
            }
            drawPop.put(pop, drawLen);
            drawMax = Math.max(drawMax, drawLen);
        }
        
        for (String pop : popSet) {
            if (!drawPop.containsKey(pop)) {
                continue;
            }
            int drawLen = drawPop.get(pop);
            if (drawLen == drawMax && scav > 0) {
                drawLen -= scav;
            }
            
            g.setColor(popColorMap.get(pop));

            if (bd.dt == DrawType.FILLED_CIRCLE) {
                g.fillRect(x, tempY, width + 1, drawLen);
            } else if (bd.dt == DrawType.EMPTY_CIRCLE) {
                g.drawRect(x, tempY, width, drawLen);
                g.drawRect(x+1, tempY + 1, width - 2, drawLen - 2);
                g.drawRect(x+2, tempY + 2, width - 4, drawLen - 4);
            } else if (bd.dt == DrawType.X) {
                g.drawRect(x, tempY, width, drawLen);
                int tickLen = 7;
                int ticks = drawLen / tickLen;
                for (int i = 0; i < ticks; i++) {
                    g.drawLine(x + 1, tempY + (i * tickLen), x + width - 1, tempY + ((i + 1) * tickLen));
                }
            } else {
                // TODO Error
            }
            tempY += drawLen;
            
        }
        g.setColor(Color.BLACK);
        if (selectedBlockDraw != null && bd.bpX == selectedBlockDraw.bpX) {
            g.drawRect(x - 1, 2, width + 2, tempY - 2);
        }
	}
	
    private void drawVcfEntry(Graphics g, DrawPoint dp) {
        int height = genePanel.getHeight() - yStart - GENE_HEIGHT;
        g.setColor(dp.c);
        
        int x = dp.x;
        int y = height - dataPntSize - 2 - dp.y;
        if (dp.type == DrawType.FILLED_CIRCLE) {
            g.fillOval(x, y, dataPntSize, dataPntSize);
        } else if (dp.type == DrawType.EMPTY_CIRCLE) {
            g.drawOval(x, y, dataPntSize, dataPntSize);
            g.drawOval(x + 1, y + 1, dataPntSize - 2, dataPntSize - 2);
        } else if (dp.type == DrawType.X) {
            Grafik.drawThickLine(g, x, y, x + dataPntSize, y + dataPntSize, 2, dp.c);
            Grafik.drawThickLine(g, x, y + dataPntSize, x + dataPntSize, y, 2, dp.c);
        } else {
            // TODO Error
        }
        
        int offset = 1;
        if (dp.type == DrawType.EMPTY_CIRCLE) {
            offset = 2;
        }
        if (selectedDrawPoint != null) {
            if (dp.sampleID.equals(selectedDrawPoint.sampleID) && dp.vcRecord.equals(selectedDrawPoint.vcRecord)) {
                g.setColor(Color.black);
                g.drawOval(x-1, y-1, dataPntSize + offset, dataPntSize + offset);
            }
        }
        
        g.setColor(Color.BLACK);
    }
    
    private void drawLegend(Graphics g) {
        FontMetrics fm = g.getFontMetrics();
        int lblMaxWidth = 0;
        int maxWidth;
        HashMap<String, String> lblMap = new HashMap<String, String>();

        lblMaxWidth = fm.stringWidth("Moderate Impact");
        if (popColorMap != null) {
            for (Entry<String, Color> colEntry : popColorMap.entrySet()) {
                HashSet<String> popSet = (showExcludes ? popIndiMapWithExcludes : popIndiMap).get(colEntry.getKey());
                String lbl = colEntry.getKey() + " (n=" + (popSet == null ? 0 : popSet.size()) + ")";
                lblMap.put(colEntry.getKey(), lbl);
                lblMaxWidth = Math.max(lblMaxWidth, fm.stringWidth(lbl));
            }
        }
        maxWidth = lblMaxWidth + 25 + dataPntSize;
        g.setColor(genePanel.getBackground());
        g.fillRect(10, 10, maxWidth, 200);
        g.setColor(Color.black);
        
        int x = 15;
        int y = 15 + fm.getHeight();
        if (drawType == DRAW_AS_CIRCLES) {
            g.drawOval(x, y, dataPntSize, dataPntSize);
            g.drawOval(x + 1, y + 1, dataPntSize - 2, dataPntSize - 2);
        } else {
            g.drawRect(x, y, dataPntSize, dataPntSize);
            g.drawRect(x+1, y + 1, dataPntSize - 2, dataPntSize - 2);
            g.drawRect(x+2, y + 2, dataPntSize - 4, dataPntSize - 4);
        }
        g.drawString("Low Impact", x + dataPntSize + 5, y + fm.getHeight() / 2 + 2);
        y += 15;
        if (drawType == DRAW_AS_CIRCLES) {
            g.fillOval(x, y, dataPntSize, dataPntSize);
        } else {
            g.fillRect(x, y, dataPntSize, dataPntSize);
        }
        g.drawString("Moderate Impact", x + dataPntSize + 5, y + fm.getHeight() / 2 + 2);
        y += 15;
        if (drawType == DRAW_AS_CIRCLES) {
            Grafik.drawThickLine(g, x, y, x + dataPntSize, y + dataPntSize, 2, Color.BLACK);
            Grafik.drawThickLine(g, x, y + dataPntSize, x + dataPntSize, y, 2, Color.BLACK);
        } else {
            g.drawRect(x, y, dataPntSize, dataPntSize);
            int tickLen = 4;
            int ticks = dataPntSize / tickLen;
            for (int i = 0; i < ticks; i++) {
                g.drawLine(x + 1, y + (i * tickLen), x + dataPntSize - 1, y + ((i + 1) * tickLen));
            }
        }
        g.drawString("High Impact", x + dataPntSize + 5, y + fm.getHeight() / 2 + 2);
        y += 15;
//        g.drawLine(10, y, 110, y);
        y += 10;
        if (popColorMap != null) {
            for (Entry<String, Color> colEntry : popColorMap.entrySet()) {
                String lbl = lblMap.containsKey(colEntry.getKey()) ? lblMap.get(colEntry.getKey()) : colEntry.getKey() + " (n=" + (showExcludes ? popIndiMapWithExcludes : popIndiMap).get(colEntry.getKey()).size() + ")";
                g.setColor(Color.BLACK);
                g.drawRect(x, y, dataPntSize, dataPntSize);
                g.drawString(lbl, x + dataPntSize + 5, y + fm.getHeight() / 2 + 2);
                g.setColor(colEntry.getValue());
                g.fillRect(x+1, y+1, dataPntSize - 1, dataPntSize - 1);
                y += dataPntSize * 2 - 1;
            }
        }
        
        
        g.setColor(Color.BLACK);
        g.drawRect(10, 10, maxWidth, 200);
        g.drawString("Key:", 15, 7 + fm.getHeight());
        g.drawLine(10, 10 + fm.getHeight(), 10 + maxWidth, 10 + fm.getHeight());
    }
    
    private void drawFreqs(Graphics g, ArrayList<VCFLocation> vcfLocs, int lblX) {
        if (!drawMAF && !drawMAC) {
            return;
        }
        Font prevFont = g.getFont();
        int spanAll = getStop() - getStart(true);
        int spanCur = stop - start;
        double prop = spanCur / (double) spanAll;
        float fsz = (float) Math.min(14d, 11 / prop);
//        float fsz = 11f;
        Font newFont = (Fonts.SOURCE_CODE_PRO_REGULAR == null ? Font.decode(Font.MONOSPACED) : Fonts.SOURCE_CODE_PRO_REGULAR).deriveFont(fsz);
        g.setFont(newFont);
        FontMetrics fm = g.getFontMetrics();
        
        if ((showExcludes ? popIndiMapWithExcludes : popIndiMap) == null) return;
        
        ArrayList<String> pops = new ArrayList<String>();

        for (Entry<String, HashSet<String>> popSetEntry : (showExcludes ? popIndiMapWithExcludes.entrySet() : popIndiMap.entrySet())) {
            for (VCFLocation vcfLoc : vcfLocs) {
                double maf = VCOps.getMAF(vcfLoc.vc, popSetEntry.getValue());
                double mac = VCOps.getMAC(vcfLoc.vc, popSetEntry.getValue());
                if ((Double.isNaN(maf) || maf == 0d) && (mac == 0d || Double.isNaN(mac))) continue;
                if (!pops.contains(popSetEntry.getKey())) {
                    pops.add(popSetEntry.getKey());
                }
                vcfLoc.mafMap.put(popSetEntry.getKey(), maf);
                vcfLoc.macMap.put(popSetEntry.getKey(), mac);
            }
        }
        
        int y = genePanel.getHeight() - yStart + 3;
        
        if (drawMAF) {
            for (int i = 0; i < pops.size(); i++) {
                for (int v = 0; v < vcfLocs.size(); v++) {
                    String draw = vcfLocs.get(v).mafMap.containsKey(pops.get(i)) ? ext.formDeci(vcfLocs.get(v).mafMap.get(pops.get(i)), 4) : "--";
                    int width = fm.stringWidth(draw);
                    if (v < vcfLocs.size() - 1) {
                        if (vcfLocs.get(v).x + width > vcfLocs.get(v + 1).x) {
                            draw = "*";
                        }
                    } else if (v > 0) {
                        String prevDraw = vcfLocs.get(v - 1).mafMap.containsKey(pops.get(i)) ? ext.formDeci(vcfLocs.get(v - 1).mafMap.get(pops.get(i)), 4) : "--";
                        if (vcfLocs.get(v).x < vcfLocs.get(v - 1).x + fm.stringWidth(prevDraw)) {
                            draw = "*";
                        }
                    }
                    
                    int yTemp = y + i * (dataPntSize * 2 - 1);
                    g.drawString(draw, vcfLocs.get(v).x, yTemp + fm.getHeight() / 2 + 1);
                }
            }
            for (int i = 0; i < pops.size(); i++) {
                int yTemp = y + i * (dataPntSize * 2 - 1) + 1;
                g.setColor(Color.BLACK);
                g.drawRect(lblX - dataPntSize - 5, yTemp, dataPntSize, dataPntSize);
                g.setColor(popColorMap.get(pops.get(i)));
                g.fillRect(lblX - dataPntSize - 4, yTemp+1, dataPntSize - 1, dataPntSize - 1);
            }
            g.setColor(Color.BLACK);
            g.setFont(prevFont);
            if (pops.size() > 0) {
                String lbl = "MAF:";
                fm = g.getFontMetrics();
                int yHalf = y + (pops.size() * (dataPntSize * 2 - 1)) / 2 + fm.getHeight() / 3 - 2;
                int tempX = lblX - dataPntSize - 10;
                g.drawString(lbl, tempX - fm.stringWidth(lbl), yHalf);
            }
            y += dataPntSize * pops.size();
            y += fm.getHeight() + 4;
        }
        if (drawMAC) {
            g.setFont(newFont);
            fm = g.getFontMetrics();
            for (int i = 0; i < pops.size(); i++) {
                for (int v = 0; v < vcfLocs.size(); v++) {
                    String draw = vcfLocs.get(v).macMap.containsKey(pops.get(i)) ? ext.formDeci(vcfLocs.get(v).macMap.get(pops.get(i)), 4) : "--";
                    int width = fm.stringWidth(draw);
                    if (v < vcfLocs.size() - 1) {
                        if (vcfLocs.get(v).x + width > vcfLocs.get(v + 1).x) {
                            draw = "*";
                        }
                    } else if (v > 0) {
                        String prevDraw = vcfLocs.get(v - 1).macMap.containsKey(pops.get(i)) ? ext.formDeci(vcfLocs.get(v - 1).macMap.get(pops.get(i)), 4) : "--";
                        if (vcfLocs.get(v).x < vcfLocs.get(v - 1).x + fm.stringWidth(prevDraw)) {
                            draw = "*";
                        }
                    }
                    int yTemp = y + i * (dataPntSize * 2 - 1);
                    g.drawString(draw, vcfLocs.get(v).x, yTemp + fm.getHeight() / 2 + 1);
                }
            }
            for (int i = 0; i < pops.size(); i++) {
                int yTemp = y + i * (dataPntSize * 2 - 1) + 1;
                g.setColor(Color.BLACK);
                g.drawRect(lblX - dataPntSize - 5, yTemp, dataPntSize, dataPntSize);
                g.setColor(popColorMap.get(pops.get(i)));
                g.fillRect(lblX - dataPntSize - 4, yTemp+1, dataPntSize - 1, dataPntSize - 1);
            }
            g.setColor(Color.BLACK);
            g.setFont(prevFont);
            if (pops.size() > 0) {
                String lbl = "MAC:";
                fm = g.getFontMetrics();
                int yHalf = y + (pops.size() * (dataPntSize * 2 - 1)) / 2 + fm.getHeight() / 3 - 2;
                int tempX = lblX - dataPntSize - 10;
                g.drawString(lbl, tempX - fm.stringWidth(lbl), yHalf);
            }
        }
        
    }
    
    private int determineExonNumber(GeneData geneData, int[] exonBnds) {
        Segment seg = new Segment(exonBnds[0], exonBnds[1]);
        for (int i = 0; i < geneData.getExonBoundaries().length; i++) {
            Segment exSeg = new Segment(geneData.getExonBoundaries()[i][0], geneData.getExonBoundaries()[i][1]);
            if(exSeg.significantOverlap(seg)) {
                return i + 1;
            }
        }
        return -1;
    }

    private Color[] getAColor(int index) {
		while (index >= colorScheme.size()) {
			colorScheme.add(ColorExt.generatePastelShades());
		}
		return colorScheme.get(index);
	}

	private int getX(int pos) {
		return (int)((double)(pos-start)/(double)(stop-start)*(double)(getWidth()-2*WIDTH_BUFFER))+WIDTH_BUFFER;
	}
	
	private String chooseNewFiles() {
		JFileChooser jfc = new JFileChooser((proj != null || geneFileName == null ? proj.PROJECT_DIRECTORY.getValue() : ext.parseDirectoryOfFile(geneFileName)));
		jfc.setMultiSelectionEnabled(true);
		if (jfc.showOpenDialog(VariantViewer.this) == JFileChooser.APPROVE_OPTION) {
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
					JOptionPane.showMessageDialog(VariantViewer.this, msg.toString()); 
				}
				
				for (File kept : keptFiles) {
					addFileToList(kept.getAbsolutePath());
				}
				return keptFiles[0].getAbsolutePath();
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
					JOptionPane.showMessageDialog(VariantViewer.this, msg.toString()); 
					return null;
				} else {
					addFileToList(file.getAbsolutePath());
					return file.getAbsolutePath();
				}
				
			}
			
		}
		return null;
	}
	
	private void addFileToList(String rawfile) {
		String file = ext.verifyDirFormat(rawfile);
		file = file.substring(0, file.length() - 1);
		String name = ext.rootOf(file);
		regionFileNameLoc.put(name, file);
		
		JCheckBoxMenuItem item = new JCheckBoxMenuItem();
		item.setAction(geneFileSelectAction);
		item.setText(name);

		regionFileNameBtn.put(name, item);
		regionButtonGroup.add(item);
		loadRecentFileMenu.add(item);
	}
	
	public void generateComponents() {
		JPanel dataPanel = new JPanel();
		dataPanel.setLayout(new MigLayout("", "[grow, fill]", "[grow, fill]"));

//		lrrPanel = new JPanel() {
//			public static final long serialVersionUID = 2L;
//			public void paintComponent(Graphics g) {
//				paintLRRPanel(g);
//				
////				float min, max;
////				g.setFont(new Font("Arial", 0, 20));
////				g.drawString("Log R Ratio", WIDTH_BUFFER, 20);
////				if (SHOW_MIDLINE) {
////					g.setColor(Color.LIGHT_GRAY);
////					g.drawLine(WIDTH_BUFFER, getHeight()-(int)((double)(0-min)/(double)(max-min)*(double)(getHeight()-2*HEIGHT_BUFFER))-HEIGHT_BUFFER, getWidth()-WIDTH_BUFFER, getHeight()-(int)((double)(0-min)/(double)(max-min)*(double)(getHeight()-2*HEIGHT_BUFFER))-HEIGHT_BUFFER);
////				}
////				g.setFont(new Font("Arial", 0, 12));
//			}
//		};
//		lrrPanel.addMouseListener(this);
//		lrrPanel.addMouseMotionListener(this);
//		lrrPanel.addMouseWheelListener(this);
//		dataPanel.add(lrrPanel, "cell 0 0");

		genePanel = new JPanel() {
			public static final long serialVersionUID = 8L;
			public void paintComponent(Graphics g) {
			    paintPanel(g);
		    }
		};
		genePanel.addMouseListener(this);
		genePanel.addMouseMotionListener(this);
		genePanel.addMouseWheelListener(this);
		dataPanel.add(genePanel, "cell 0 0");

		getContentPane().add(dataPanel, BorderLayout.CENTER);

		JPanel sampPanel = new JPanel();
		((FlowLayout)sampPanel.getLayout()).setVgap(0);
		previousGene = new JButton(Grafik.getImageIcon("images/firstLast/Left.gif", true));
		previousGene.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLeft.gif", true));
		previousGene.addActionListener(this);
		previousGene.setActionCommand(PREVIOUS_REGION);
		previousGene.setPreferredSize(new Dimension(25, 25));
		JPanel compPanel = new JPanel(new MigLayout("align center, fill, gap 0", "[grow, center]", "[]5[21:21:21]5[]"));
		
		JPanel regionPanel = new JPanel();
		((FlowLayout)regionPanel.getLayout()).setVgap(0);
		geneListCmb = new JComboBox<String>();
        DefaultListCellRenderer dlcr = new DefaultListCellRenderer();
        dlcr.setHorizontalAlignment(DefaultListCellRenderer.CENTER);
        geneListCmb.setRenderer(dlcr);
        geneListCmb.setBorder(BorderFactory.createEtchedBorder());
        geneListCmb.setEditable(false);
		Font font = new Font("Arial", 0, 14);
		geneListCmb.setFont(font);
		geneListCmb.setAction(new AbstractAction() {
		    private static final long serialVersionUID = 1L;
		    public void actionPerformed(ActionEvent e) {
		        showGene(geneListCmb.getSelectedIndex());
		        updateGUI();
		    }
		});
		geneListCmb.setPreferredSize(new Dimension(geneListCmb.getPreferredSize().width, 26));
		nextGene = new JButton(Grafik.getImageIcon("images/firstLast/Right.gif", true));
        nextGene.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dRight.gif", true));
        nextGene.addActionListener(this);
        nextGene.setActionCommand(NEXT_REGION);
        nextGene.setPreferredSize(new Dimension(25, 25));
        
        regionPanel.add(previousGene);
        regionPanel.add(geneListCmb);
		regionPanel.add(nextGene);

		compPanel.add(regionPanel, "cell 0 0");
		compPanel.setPreferredSize(new Dimension(compPanel.getPreferredSize().width, 105));
		
        commentLabel = new JLabel(" ", JLabel.CENTER);
        commentLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
        commentLabel.setFont(font);
        commentLabel.setToolTipText("Click to Edit Comment");
        commentLabel.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                super.mouseClicked(e);
                commentLabel.setVisible(false);
                commentField.setVisible(true);
                commentField.requestFocusInWindow();
            }
        });
        compPanel.add(commentLabel, "cell 0 1, hidemode 3");
        
        commentField = new JTextField(20);
        commentField.setAlignmentX(Component.CENTER_ALIGNMENT);
        commentField.setFont(font);
        commentField.addFocusListener(new FocusAdapter() {
            @Override
            public void focusLost(FocusEvent e) {
                super.focusLost(e);
                String newComment = commentField.getText();
                if (newComment.length() == 0) {
                    newComment = null;
                    geneToCommentMap.remove(geneList.get(geneIndex));
                } else {
                    geneToCommentMap.put(geneList.get(geneIndex), newComment);
                }
                if (newComment != null) {
                    commentLabel.setText("gene #"+(geneIndex+1)+":  "+ geneToCommentMap.get(geneList.get(geneIndex)));
                } else {
                    commentLabel.setText(" -- no comment -- ");
                }
                commentLabel.setVisible(true);
                commentField.setVisible(false);
            }
        });
        commentField.setVisible(false);
        compPanel.add(commentField, "cell 0 1, hidemode 3");
		
		isoformList = new JComboBox<String>();
		isoformList.setFont(font);
		
		dlcr = new DefaultListCellRenderer();
	    dlcr.setHorizontalAlignment(DefaultListCellRenderer.CENTER);
	    isoformList.setRenderer(dlcr);
	    isoformList.setBorder(BorderFactory.createEtchedBorder());
	    isoformList.setEditable(false);
	    isoformList.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                @SuppressWarnings("unchecked")
                JComboBox<String> jcb = (JComboBox<String>)e.getSource();
                int index = jcb.getSelectedIndex();
                if (index == isoformsPresent.length-1) {
//                    createIsoformList();
//                    isoformIndex = COLLAPSE_ISOFORMS;
                } else if (!isoformsPresent[index].equals(isoform)) {
//                    isoformIndex = index;
//                    updateSample(isoformsPresent[index]);
                }
                selectedBlockDraw = null;
                selectedDrawPoint = null;
                selectedRect = null;
                parseLocation(geneToRegionMap.get(geneList.get(geneIndex)).get(isoformList.getSelectedItem()));
                
                updateGUI();
            }
        });
        compPanel.add(isoformList, "cell 0 2");

		JPanel overPanel = new JPanel();
		overPanel.setLayout(new BoxLayout(overPanel, BoxLayout.LINE_AXIS));
		
		JSeparator sep = new JSeparator(SwingConstants.VERTICAL);
		sep.setMaximumSize(new Dimension(1, 150));
		overPanel.setBorder(new EmptyBorder(5, 0, 5, 0));
		
		overPanel.add(Box.createHorizontalGlue());
		overPanel.add(compPanel);
		overPanel.add(Box.createHorizontalGlue());
		
		getContentPane().add(overPanel, BorderLayout.NORTH);

//		InputMap inputMap = bafPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
//		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, InputEvent.ALT_MASK), PREVIOUS_REGION);
//		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, InputEvent.ALT_MASK), NEXT_REGION);
//
//		ActionMap actionMap = bafPanel.getActionMap();
//		actionMap.put(FIRST_REGION, new AbstractAction() {
//			public static final long serialVersionUID = 9L;
//			public void actionPerformed(ActionEvent e) {
//				showGene(0);
//			}
//		});
//		actionMap.put(PREVIOUS_REGION, new AbstractAction() {
//			public static final long serialVersionUID = 10L;
//			public void actionPerformed(ActionEvent e) {
//				showGene(Math.max(geneIndex-1, 0));
//			}
//		});
//		actionMap.put(NEXT_REGION, new AbstractAction() {
//			public static final long serialVersionUID = 11L;
//			public void actionPerformed(ActionEvent e) {
//				showGene(Math.min(geneIndex+1, geneRegions.length-1));
//			}
//		});
//		actionMap.put(LAST_REGION, new AbstractAction() {
//			public static final long serialVersionUID = 12L;
//			public void actionPerformed(ActionEvent e) {
//				showGene(geneRegions.length - 1);
//			}
//		});
//		bafPanel.setActionMap(actionMap);
	}
	
	private void doScreenCapture(String filename) {
	    BufferedImage cap = generateScreenshot();
	    
	    if (filename == null) {
	        TransferableImage ti = new TransferableImage(cap);
	        Clipboard c = Toolkit.getDefaultToolkit().getSystemClipboard();
            c.setContents( ti, null );
	    } else {
	        try {
	            ImageIO.write(cap, "png", new File(filename));
	        } catch (IOException e) {
	            if (proj != null) {
	                proj.getLog().reportException(e);
	                proj.message("Error occured while writing screen capture to file.  Please check log for more details.");
	            }
	        }
	    }
	}
	
	
	private BufferedImage generateScreenshot() {
//	    int lW = lrrPanel.getWidth();
//	    int bW = bafPanel.getWidth();
	    int cW = genePanel.getWidth();
//	    int lH = lrrPanel.getHeight();
//	    int bH = bafPanel.getHeight();
	    int cH = genePanel.getHeight();
//	    BufferedImage imageLrr = new BufferedImage(lW, lH, BufferedImage.TYPE_INT_RGB);
//	    BufferedImage imageBaf = new BufferedImage(bW, bH, BufferedImage.TYPE_INT_RGB);
	    BufferedImage imageCnv = new BufferedImage(cW, cH, BufferedImage.TYPE_INT_RGB);
	    
//	    Graphics g = imageLrr.getGraphics();
//        g.setColor(lrrPanel.getBackground());
//        g.fillRect(0, 0, imageLrr.getWidth(), imageLrr.getHeight());
//	    lrrPanel.paintAll(g);
	    
	    Graphics g = imageCnv.getGraphics();
        g.setColor(genePanel.getBackground());
        g.fillRect(0, 0, imageCnv.getWidth(), imageCnv.getHeight());
        genePanel.paintAll(g);
	    
//	    g = imageBaf.getGraphics();
//        g.setColor(bafPanel.getBackground());
//        g.fillRect(0, 0, imageBaf.getWidth(), imageBaf.getHeight());
//	    bafPanel.paintAll(g);
	    
//        int w = Math.max(lW, Math.max(bW, cW));
	    int w = cW;//Math.max(lW, cW);
	    int h = /*lH + bH + */cH + HEIGHT_BUFFER;
	    BufferedImage finalImage = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
	    
	    g = finalImage.getGraphics();
	    g.setColor(genePanel.getBackground());
	    g.fillRect(0, 0, finalImage.getWidth(), finalImage.getHeight());
	    g.drawImage(imageCnv, 0, 0, null);
//	    g.drawImage(imageCnv, 0, lH + 5, null);
//	    g.drawImage(imageBaf, 0, lH + cH + 10, null);
	    g.dispose();
	    
	    return finalImage;
	}
	
	private JMenuBar createMenuBar() {
		JMenuBar menuBar = new JMenuBar();
		
		JMenu fileMenu = new JMenu("File");
		fileMenu.setMnemonic(KeyEvent.VK_F);
		
		JMenuItem newRegionFile = new JMenuItem();
		newRegionFile.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
	            NewRegionListDialog newRgnList = new NewRegionListDialog(proj == null ? null : proj.getSamples(), proj == null ? null : proj.PROJECT_DIRECTORY.getValue(), true);
	            newRgnList.setModal(true);
	            newRgnList.setVisible(true);
	            if (newRgnList.getReturnCode() == JOptionPane.YES_OPTION) {
	                String rgnFile = newRgnList.getFileName();
	                addFileToList(rgnFile);
	                String file = ext.verifyDirFormat(rgnFile);
	                file = file.substring(0, file.length() - 1);
	                String name = ext.rootOf(file);
	                regionFileNameBtn.get(name).setSelected(true);
	                regionFileNameBtn.get(name).doClick();
	            }
            }
        });
		newRegionFile.setText("New Gene List File");
		newRegionFile.setMnemonic(KeyEvent.VK_N);
		Font font = new Font("Arial", 0, 12);
        newRegionFile.setFont(font);
		fileMenu.add(newRegionFile);
		
		JMenuItem loadRegionFile = new JMenuItem();
		loadRegionFile.setAction(loadGeneListFileAction);
		loadRegionFile.setText(REGION_LIST_NEW_FILE);
		loadRegionFile.setMnemonic(KeyEvent.VK_L);
		loadRegionFile.setFont(font);
		fileMenu.add(loadRegionFile);
		loadRecentFileMenu = new JMenu("Load Recent Gene List...");
		loadRecentFileMenu.setMnemonic(KeyEvent.VK_R);
		loadRecentFileMenu.setFont(font);
		fileMenu.add(loadRecentFileMenu);
		
		JMenuItem screencap1 = new JMenuItem();
		screencap1.setAction(screencapAction);
		screencap1.setMnemonic(KeyEvent.VK_S);
		screencap1.setText("Screen Capture");
		screencap1.setFont(font);
		fileMenu.add(screencap1);
		
		JMenuItem screencap2 = new JMenuItem();
		screencap2.setAction(screencapClipboardAction);
		screencap2.setMnemonic(KeyEvent.VK_C);
		screencap2.setText("Screen Capture to Clipboard");
		screencap2.setFont(font);
		fileMenu.add(screencap2);
		
		menuBar.add(fileMenu);
		
		regionButtonGroup = new ButtonGroup();
		if (proj != null) {
			String[] files = proj.GENE_LIST_FILENAMES.getValue();
			String name;
			for (String file : files) {
				name = ext.rootOf(file);
				regionFileNameLoc.put(name, file);
				JCheckBoxMenuItem menuItem = new JCheckBoxMenuItem();
				menuItem.setAction(geneFileSelectAction);
				menuItem.setFont(font);
				boolean found = Files.exists(file);
				menuItem.setText(name + (found ? "" : " -- [file not found]"));
				menuItem.setEnabled(found);
				regionFileNameBtn.put(name, menuItem);
				regionButtonGroup.add(menuItem);
				loadRecentFileMenu.add(menuItem);
			}
		}
		
		JMenu disp = new JMenu("Display");
		disp.setMnemonic(KeyEvent.VK_D);
		menuBar.add(disp);
		final JCheckBoxMenuItem paintExonBoundariesChk = new JCheckBoxMenuItem();
		disp.add(paintExonBoundariesChk);
		final JCheckBoxMenuItem paintIntronsChk = new JCheckBoxMenuItem();
		disp.add(paintIntronsChk);
		final JCheckBoxMenuItem paintExonBoxesChk = new JCheckBoxMenuItem();
		disp.add(paintExonBoxesChk);
		final JCheckBoxMenuItem fillExonsChk = new JCheckBoxMenuItem();
		disp.add(fillExonsChk);
		final JCheckBoxMenuItem paintExonLinesChk = new JCheckBoxMenuItem();
		disp.add(paintExonLinesChk);
		final JCheckBoxMenuItem equalizeExonsChk = new JCheckBoxMenuItem();
        disp.add(equalizeExonsChk);
        disp.add(new JSeparator());
        final JCheckBoxMenuItem paintExonNumbersChk = new JCheckBoxMenuItem();
        disp.add(paintExonNumbersChk);
        ButtonGroup exonLblGroup = new ButtonGroup();
        final JRadioButtonMenuItem exonLblLeft = new JRadioButtonMenuItem();
        exonLblGroup.add(exonLblLeft);
        disp.add(exonLblLeft);
        final JRadioButtonMenuItem exonLblCenter = new JRadioButtonMenuItem();
        exonLblGroup.add(exonLblCenter);
        disp.add(exonLblCenter);
        final JRadioButtonMenuItem exonLblRight = new JRadioButtonMenuItem();
        exonLblGroup.add(exonLblRight);
        disp.add(exonLblRight);
        disp.add(new JSeparator());
        final JCheckBoxMenuItem chkbxDisplayLegend = new JCheckBoxMenuItem();
        disp.add(chkbxDisplayLegend);
        chkbxDisplayExcludes = new JCheckBoxMenuItem();
        disp.add(chkbxDisplayExcludes);
        disp.add(new JSeparator());
        ButtonGroup displayGroup = new ButtonGroup();
        final JRadioButtonMenuItem displayPopBlocks = new JRadioButtonMenuItem();
        displayGroup.add(displayPopBlocks);
        disp.add(displayPopBlocks);
        final JRadioButtonMenuItem displayIndiCircles = new JRadioButtonMenuItem();
        displayGroup.add(displayIndiCircles);
        disp.add(displayIndiCircles);
        disp.add(new JSeparator());
        final JCheckBoxMenuItem displayMAF = new JCheckBoxMenuItem();
        disp.add(displayMAF);
        final JCheckBoxMenuItem displayMAC = new JCheckBoxMenuItem();
        disp.add(displayMAC);
        
        paintExonBoundariesChk.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent arg0) {
                paintExonBoundaries = paintExonBoundariesChk.isSelected();
                updateGUI();
            }
        });
        paintExonBoundariesChk.setText("Paint Exon Boundaries");
        paintExonBoundariesChk.setSelected(paintExonBoundaries);
        
		paintIntronsChk.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                paintIntrons = paintIntronsChk.isSelected();
//                parseLocation(geneRegions[geneIndex][1]);
                parseLocation(geneToRegionMap.get(geneList.get(geneIndex)).get(isoformList.getSelectedItem()));
                updateGUI();
            }
        });
        paintIntronsChk.setText("Paint Introns");
        paintIntronsChk.setSelected(paintIntrons);
        
        paintExonBoxesChk.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                paintExonBoxes = paintExonBoxesChk.isSelected();
                updateGUI();
            }
        });
        paintExonBoxesChk.setText("Paint Exon Boxes");
        paintExonBoxesChk.setSelected(paintExonBoxes);
        
        fillExonsChk.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                fillExons = fillExonsChk.isSelected();
                updateGUI();
            }
        });
        fillExonsChk.setText("Fill Exon Boxes");
        fillExonsChk.setSelected(fillExons);
        
        paintExonLinesChk.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                paintInternalLine = paintExonLinesChk.isSelected();
                updateGUI();
            }
        });
        paintExonLinesChk.setText("Paint Internal Exon Lines");
        paintExonLinesChk.setSelected(paintInternalLine);
        
		equalizeExonsChk.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                equalizeExonLength = equalizeExonsChk.isSelected();
//                parseLocation(geneRegions[geneIndex][1]);
                parseLocation(geneToRegionMap.get(geneList.get(geneIndex)).get(isoformList.getSelectedItem()));
                updateGUI();
            }
        });
		equalizeExonsChk.setText("Equalize Exon Length");
		equalizeExonsChk.setSelected(equalizeExonLength);
		
		paintExonNumbersChk.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                paintExonNumbers = paintExonNumbersChk.isSelected();
                updateGUI();
            }
        });
		paintExonNumbersChk.setText("Display Exon Numbers");
		paintExonNumbersChk.setSelected(paintExonNumbers);
		
		exonLblLeft.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                exonLabelLocation = EXON_LBL_LEFT;
                updateGUI();
            }
        });
		exonLblLeft.setText("Left-Aligned");
		exonLblLeft.setSelected(true);
		exonLblCenter.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                exonLabelLocation = EXON_LBL_CENTER;
                updateGUI();
            }
        });
        exonLblCenter.setText("Center-Aligned");
		exonLblRight.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                exonLabelLocation = EXON_LBL_RIGHT;
                updateGUI();
            }
        });
        exonLblRight.setText("Right-Aligned");
		
        chkbxDisplayLegend.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                showLegend = chkbxDisplayLegend.isSelected();
                updateGUI();
            }
        });
        chkbxDisplayLegend.setText("Show Legend");
        chkbxDisplayLegend.setSelected(showLegend);
        
        chkbxDisplayExcludes.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                showExcludes = chkbxDisplayExcludes.isSelected();
                updateGUI();
            }
        });
        chkbxDisplayExcludes.setText("Show Excluded");
        
        displayPopBlocks.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                drawType = DRAW_AS_BLOCKS;
                selectedRect = null;
                selectedDrawPoint = null;
                selectedBlockDraw = null;
                updateGUI();
            }
        });
        displayPopBlocks.setText("Display Populations");
        displayPopBlocks.setSelected(true);
        
        displayIndiCircles.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent e) {
                drawType = DRAW_AS_CIRCLES;
                selectedRect = null;
                selectedDrawPoint = null;
                selectedBlockDraw = null;
                updateGUI();
            }
        });
        displayIndiCircles.setText("Display Individuals");
        
        displayMAF.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent arg0) {
                drawMAF = displayMAF.isSelected();
                if (drawMAF) {
                    yStart = Y_START + dataPntSize * popSet.size();
                } else {
                    yStart = Y_START;
                }
                if (drawMAC) {
                    yStart += dataPntSize * popSet.size();
                }
                updateGUI();
            }
        });
        displayMAF.setText("Display Minor Allele Frequencies (MAF)");
		
        displayMAC.setAction(new AbstractAction() {
            private static final long serialVersionUID = 1L;
            @Override
            public void actionPerformed(ActionEvent arg0) {
                drawMAC = displayMAC.isSelected();
                if (drawMAC) {
                    yStart = Y_START + dataPntSize * popSet.size();
                } else {
                    yStart = Y_START;
                }
                if (drawMAF) {
                    yStart += dataPntSize * popSet.size();
                }
                updateGUI();
            }
        });
        displayMAC.setText("Display Minor Allele Counts (MAC)");

        JMenu act = new JMenu("Actions");
        act.setMnemonic(KeyEvent.VK_A);
        menuBar.add(act);

        JMenuItem launchTrailer = new JMenuItem();
        launchTrailer.setText(TO_TRAILER_PLOT);
        launchTrailer.setMnemonic(KeyEvent.VK_T);
        launchTrailer.setFont(font);
        launchTrailer.addActionListener(this);
        act.add(launchTrailer);
        JMenuItem launchScatter = new JMenuItem();
        launchScatter.setText(TO_SCATTER_PLOT);
        launchScatter.setMnemonic(KeyEvent.VK_S);
        launchScatter.setFont(font);
        launchScatter.addActionListener(this);
        act.add(launchScatter);
        JMenuItem launchComp = new JMenuItem();
        launchComp.setText(TO_COMP_PLOT);
        launchComp.setMnemonic(KeyEvent.VK_C);
        launchComp.setFont(font);
        launchComp.addActionListener(this);
        act.add(launchComp);
        // act.addSeparator();
		
		return menuBar;
	}
	
	public void actionPerformed(ActionEvent ae) {
		String command = ae.getActionCommand();
//		String[] filenames;

		if (command.equals(FIRST_REGION)) {
			if (geneToRegionMap == null || geneToRegionMap.size() == 0) {
				JOptionPane.showMessageDialog(null, "Error - No regions have been loaded", "Error", JOptionPane.ERROR_MESSAGE);
				return;
			}
			showGene(0);
		} else if (command.equals(PREVIOUS_REGION)) {
			if (geneToRegionMap == null || geneToRegionMap.size() == 0) {
				JOptionPane.showMessageDialog(null, "Error - No regions have been loaded", "Error", JOptionPane.ERROR_MESSAGE);
				return;
			}
			showGene(Math.max(geneIndex-1, 0));
		} else if (command.equals(NEXT_REGION)) {
			if (geneToRegionMap == null || geneToRegionMap.size() == 0) {
//				filenames = proj.getIndividualRegionLists();
//				if (filenames.length == 0) {
//					JOptionPane.showMessageDialog(null, "Error - No regions have been loaded, since there no individual CNV region files defined in the properties file", "Error", JOptionPane.ERROR_MESSAGE);
					JOptionPane.showMessageDialog(null, "Error - No regions have been loaded", "Error", JOptionPane.ERROR_MESSAGE);
//				} else {
//					JOptionPane.showMessageDialog(null, "Error - No regions have been loaded; files include: "+Array.toStr(filenames, ", "), "Error", JOptionPane.ERROR_MESSAGE);
//				}
				return;
			}
			showGene(Math.min(geneIndex+1, geneToRegionMap.size() - 1));
		} else if (command.equals(LAST_REGION)) {
			if (geneToRegionMap == null || geneToRegionMap.size() == 0) {
				JOptionPane.showMessageDialog(null, "Error - No regions have been loaded", "Error", JOptionPane.ERROR_MESSAGE);
				return;
			}
			showGene(geneToRegionMap.size() - 1);
//		if (command.equals(FIRST_REGION)) {
//			if (geneRegions == null || geneRegions.length == 0) {
//				JOptionPane.showMessageDialog(null, "Error - No regions have been loaded", "Error", JOptionPane.ERROR_MESSAGE);
//				return;
//			}
//			showGene(0);
//		} else if (command.equals(PREVIOUS_REGION)) {
//			if (geneRegions == null || geneRegions.length == 0) {
//				JOptionPane.showMessageDialog(null, "Error - No regions have been loaded", "Error", JOptionPane.ERROR_MESSAGE);
//				return;
//			}
//			showGene(Math.max(geneIndex-1, 0));
//		} else if (command.equals(NEXT_REGION)) {
//			if (geneRegions == null || geneRegions.length == 0) {
////				filenames = proj.getIndividualRegionLists();
////				if (filenames.length == 0) {
////					JOptionPane.showMessageDialog(null, "Error - No regions have been loaded, since there no individual CNV region files defined in the properties file", "Error", JOptionPane.ERROR_MESSAGE);
//					JOptionPane.showMessageDialog(null, "Error - No regions have been loaded", "Error", JOptionPane.ERROR_MESSAGE);
////				} else {
////					JOptionPane.showMessageDialog(null, "Error - No regions have been loaded; files include: "+Array.toStr(filenames, ", "), "Error", JOptionPane.ERROR_MESSAGE);
////				}
//				return;
//			}
//			showGene(Math.min(geneIndex+1, geneRegions.length-1));
//		} else if (command.equals(LAST_REGION)) {
//			if (geneRegions == null || geneRegions.length == 0) {
//				JOptionPane.showMessageDialog(null, "Error - No regions have been loaded", "Error", JOptionPane.ERROR_MESSAGE);
//				return;
//			}
//			showGene(geneRegions.length - 1);
		} else if (command.equals(TO_SCATTER_PLOT)) {
//		    if (proj == null) {
//		        JOptionPane.showConfirmDialog(this, "Error - a Project is required to open ScatterPlot", "Error - no Project", JOptionPane.CANCEL_OPTION, JOptionPane.ERROR_MESSAGE);
//		        return;
//		    } 
//			String[] listOfMarkers;
//			
//			listOfMarkers = new String[stopMarker-startMarker];
//			if (listOfMarkers.length == 0) {
//				JOptionPane.showMessageDialog(null, "There are no markers within this region; ScatterPlot will not bother launching", "Error", JOptionPane.ERROR_MESSAGE);
//			} else {
//				for (int i = startMarker; i < stopMarker; i++) {
//					listOfMarkers[i-startMarker] = markerNames[i];
//				}
//				ScatterPlot.createAndShowGUI(proj, listOfMarkers, null, false);
//			}
		} else if (command.equals(TO_COMP_PLOT)) {
		    if (proj == null) {
		        JOptionPane.showConfirmDialog(this, "Error - a Project is required to open CompPlot", "Error - no Project", JOptionPane.CANCEL_OPTION, JOptionPane.ERROR_MESSAGE);
		        return;
		    } 
		    final Region toRegion = new Region(new int[]{chr, start, stop});
		    SwingUtilities.invokeLater(new Runnable() {
                @Override
                public void run() {
                    CompPlot cp = new CompPlot(proj);
                    cp.setRegion(toRegion);
                }
            });
		} else {
			System.err.println("Error - unknown command '"+command+"'");
		}
	}
	
	Rectangle selectedRect;
	BlockDraw selectedBlockDraw;
	DrawPoint selectedDrawPoint;
	ArrayList<Rectangle> activeRects = new ArrayList<Rectangle>();
	ArrayList<BlockDraw> activeBlocks = new ArrayList<VariantViewer.BlockDraw>();
	ArrayList<DrawPoint> activePoints = new ArrayList<VariantViewer.DrawPoint>();
	HashSet<VariantContext> drawnFreqs = new HashSet<VariantContext>();

	public void mouseClicked(MouseEvent e) {
        int x = e.getX();
        int y = genePanel.getHeight() - yStart - GENE_HEIGHT - e.getY() - 2;

        for (int i = 0; i < activeRects.size(); i++) {
            if (activeRects.get(i).contains(x, y)) {
                selectedRect = activeRects.get(i);
                if (drawType == DRAW_AS_BLOCKS) {
                    selectedBlockDraw = activeBlocks.get(i);
                } else if (drawType == DRAW_AS_CIRCLES) {
                    selectedDrawPoint = activePoints.get(i);
                }
                MouseEvent phantom = new MouseEvent(e.getComponent(), MouseEvent.MOUSE_MOVED, System.currentTimeMillis(), 0, x, e.getY(), 0, false);
                ToolTipManager.sharedInstance().mouseMoved(phantom); // order of mouseMoved calls doesn't matter, but both are necessary
                this.mouseMoved(phantom);
                VariantViewer.this.repaint();
                return;
            }
        }
        
        selectedRect = null;
        selectedBlockDraw = null;
        selectedDrawPoint = null;
        VariantViewer.this.repaint();
//        if (e.getButton()==MouseEvent.BUTTON1) {
//            zoomProportionally(false, e.getPoint(), true);
//        } else if (e.getButton()==MouseEvent.BUTTON3) {
//            zoom(1, 1);
//        }
    }
	
    int defaultInitial = ToolTipManager.sharedInstance().getInitialDelay();
    int defaultReshow = ToolTipManager.sharedInstance().getReshowDelay();
    
    @Override
    public void mouseEntered(MouseEvent e) {
        ToolTipManager.sharedInstance().setReshowDelay(3);
        ToolTipManager.sharedInstance().setInitialDelay(3);
    }
    
    @Override
    public void mouseExited(MouseEvent e) {
        ToolTipManager.sharedInstance().setReshowDelay(defaultReshow);
        ToolTipManager.sharedInstance().setInitialDelay(defaultInitial);
    }

	public void mousePressed(MouseEvent e) {
		startX = e.getPoint().x;
		inDrag = true;
	}

	public void mouseReleased(MouseEvent e) {
		inDrag = false;
	}
	
	private int getStart(boolean buffered) {
//	    return geneData[geneIndex][0].getStart() - MIN_BUFFER;
	    GeneData gd = getCurrentGeneData();
	    return gd == null ? (buffered ? -1 * MIN_BUFFER : 0) : getCurrentGeneData().getStart() - (buffered ? MIN_BUFFER : 0);
	}

	private int getStop() {
//	    int stop = geneData[geneIndex][0].getStart();
        GeneData gd = getCurrentGeneData();
        if (gd == null) return MIN_BUFFER;
        int stop = gd.getStart();
        if (paintIntrons) {
//            stop = geneData[geneIndex][0].getStop();
            stop = gd.getStop();
        } else {
//            int[][] exons = geneData[geneIndex][0].getExonBoundaries();
            int[][] exons = gd.getExonBoundaries();
            for (int i = 0; i < exons.length; i++) {
                stop += equalizeExonLength ? EQUALIZED_EXON_BP_LENGTH : (exons[i][1] - exons[i][0]);
                if (i < exons.length - 1) {
                    if (paintIntrons) {
                        stop += exons[i + 1][0] - exons[i][1];
                    } else {
                        stop += intronBPWidth;
                    }
                }
            }
        }
        return stop + MIN_BUFFER;
    }
	
	public void mouseDragged(MouseEvent e) {
		int curX = e.getPoint().x;
		int distance = startX - curX;

		distance *= (stop - start) / (getWidth() - 2 * WIDTH_BUFFER);

		if (distance < 0) {
			distance = Math.max(distance, 1 - start);
		} else {
			distance = Math.min(distance, getStop() - stop);
		}

	    if ((start <= getStart(true) && distance < 0) || (stop >= getStop() && distance > 0)) {

		} else {
			start += distance;
			stop += distance;
		}

		if (inDrag) {
			updateGUI();
			startX = curX;
		}
	}

	public void mouseMoved(MouseEvent e) {
	    if (selectedRect == null) return;
	    int x = e.getX();
        int y = genePanel.getHeight() - yStart - GENE_HEIGHT - e.getY() - 2;
        
        if (selectedRect.contains(x, y)) {
            if (genePanel.getToolTipText() == null) {
                genePanel.setToolTipText(selectedBlockDraw == null ? buildToolTip(selectedDrawPoint) : buildToolTip(selectedBlockDraw));
            }
        } else {
            genePanel.setToolTipText(null);
        }
        
	}
	
	private String buildToolTip(DrawPoint dp) {
	    StringBuilder txtBld = new StringBuilder("<html><pre>");
	    txtBld.append("Sample ID: ").append(dp.sampleID).append("<br />");
        txtBld.append("</pre><hr><pre>");
        txtBld.append("Genotype: ").append(dp.vcRecord.vc.getGenotype(dp.sampleID).toBriefString());
        txtBld.append("</pre><pre>");
        txtBld.append("Genotype Quality: ").append(dp.vcRecord.vc.getGenotype(dp.sampleID).getGQ());
        txtBld.append("</pre><hr><pre>");
        txtBld.append("Position: ").append(dp.vcRecord.vc.getContig()).append(":").append(dp.vcRecord.vc.getStart());
        if (dp.vcRecord.vc.getStart() != dp.vcRecord.vc.getEnd()) {
            txtBld.append("-").append(dp.vcRecord.vc.getEnd());
        }
        txtBld.append("</pre><hr><pre>");
        txtBld.append("File: ").append(dp.vcRecord.file);
        txtBld.append("</pre></html>");
	    return txtBld.toString();
	}
	
	private String buildToolTip(BlockDraw bd) {
	    StringBuilder txtBld = new StringBuilder("<html><pre>");
	    txtBld.append("Total = ").append(bd.aff).append("<br />");
	    txtBld.append("</pre><hr><pre>");
	    int lenMax = 0;
	    for (String pop : bd.gen.keySet()) {
	        lenMax = Math.max(lenMax, pop == null ? 4 : pop.length());
	    }
	    int szMax = 0;
	    for (Integer t : bd.gen.values()) {
	        szMax = Math.max(szMax, ("" + t).length());
	    }
	    String format = "%1$" + lenMax + "s = %2$" + szMax + "d";
	    for (Entry<String, Integer> entry : bd.gen.entrySet()) {
	        String s = String.format(format, entry.getKey(), entry.getValue());
	        txtBld.append(s).append("<br />");
//	        txtBld.append(entry.getKey()).append(" = ").append(entry.getValue()).append("<br />");
	    }
        txtBld.append("</pre><hr><pre>");
        txtBld.append("File: ").append(bd.vcRecord.file);
        txtBld.append("</pre></html>");
	    return txtBld.toString();
	}

	public void mouseWheelMoved(MouseWheelEvent e) {
		zoomProportionally(e.getWheelRotation() > 0, e.getPoint(), false);
	}

	public void zoomProportionally(boolean outNotIn, Point p, boolean center) {
		int width = genePanel.getWidth()-2*WIDTH_BUFFER;
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
		int dist = stop - start;
		start = start - (int)(leftProportion*dist);
		stop = stop + (int)(rightProportion*dist);
		updateGUI();
	}
	
   public boolean loadMarkers() {
        Hashtable<String,String> hash;
        byte[] chrs;
        long time;
        int chr;

        time = new Date().getTime();

        hash = proj.getFilteredHash();
        markerSet = new PreparedMarkerSet(proj.getMarkerSet());
        if (markerSet == null) {
            JOptionPane.showMessageDialog(null, "Error - Failed to load the MarkerSet file; make sure the raw data is parsed", "Error", JOptionPane.ERROR_MESSAGE);
            log.reportError("Error - failed to load MarkerSet for project "+proj.PROJECT_NAME.getValue()+"; make sure the raw data is parsed");
            return false;
        }
        markerNames = markerSet.getMarkerNames();
        chrs = markerSet.getChrs();
        positions = markerSet.getPositions();

        dropped = new boolean[markerNames.length];
        chrBoundaries = new int[27][2];
        for (int i = 0; i<chrBoundaries.length; i++) {
//	          chrBoundaries[i][0] = chrBoundaries[i][1] = -1;
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

	public void createIsoformList() {
		FontMetrics fontMetrics;
		int maxWidth;

		fontMetrics = isoformList.getFontMetrics(isoformList.getFont());
		maxWidth = fontMetrics.stringWidth(COLLAPSE_ISOFORMS_KEY);
		
		HashMap<String, GeneData> isoMap = geneToIsoformMap.get(geneList.get(geneIndex));
		boolean missing = isoMap == null || isoMap.size() == 0;
		
		if (missing) {
			isoformsPresent = new String[] {"Gene Not Found"};
			maxWidth = Math.max(maxWidth, fontMetrics.stringWidth(isoformsPresent[0]));
		} else {
			isoformsPresent = new String[isoMap.size() == 2 ? 1 : isoMap.size()];
			int index = 0;
			for (Entry<String, GeneData> isoEntry : isoMap.entrySet()) {
			    if (isoEntry.getKey().equals(COLLAPSE_ISOFORMS_KEY)) {
			        continue;
			    }
				isoformsPresent[index] = isoEntry.getKey();
				maxWidth = Math.max(maxWidth, fontMetrics.stringWidth(isoformsPresent[index++]));
            }
			if (isoMap.size() > 2) {
			    isoformsPresent[index] = COLLAPSE_ISOFORMS_KEY;
			}
		}
		DefaultComboBoxModel<String> dcbm = new DefaultComboBoxModel<String>(isoformsPresent);
		isoformList.setModel(dcbm);
		isoformList.setPreferredSize(new Dimension(maxWidth+50, 30));
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
		}
		start = loc[1] - MIN_BUFFER;
		stop = getStop();
//		if (start==-1||start<0) {
//			start = 1;
//		}
//		if (stop == -1 || stop > positions[chrBoundaries[chr][1]]) {
//			stop = positions[chrBoundaries[chr][1]];
//		}
		
		updateGUI();
	};

	public void updateGUI() {
		if (start < getStart(true)) {
			start = getStart(true);
		}
		startMarker = Array.binarySearch(positions, start, chrBoundaries[chr][0], chrBoundaries[chr][1], false);

//		if (stop >= geneData[geneIndex][0].getStop() + MIN_BUFFER) {
//		    stop = geneData[geneIndex][0].getStop() + MIN_BUFFER;
//		}
		if (stop > getStop()) {
			stop = getStop();
		}
		stopMarker = Array.binarySearch(positions, stop, chrBoundaries[chr][0], chrBoundaries[chr][1], false);

		if (startMarker == -1) {
			System.err.println("Error - failed to find startMarker");
//			startMarker = chrBoundaries[chr][0];
		}

		if (stopMarker == -1) {
//			System.err.println("Error - failed to find stopMarker");
			stopMarker = chrBoundaries[chr][1];
		}
		
//		displayIndex();
		
		repaint();
	}
	
	public void showGene(int geneIndex) {
	    this.geneIndex = geneIndex;
		if (geneToRegionMap == null || geneToRegionMap.size() == 0) {
			geneListCmb.setSelectedIndex(geneIndex);
			commentLabel.setText(" ");
			return;
		}
		geneListCmb.setSelectedIndex(geneIndex);
		createIsoformList();
		isoformList.setSelectedItem(COLLAPSE_ISOFORMS_KEY);
		parseLocation(geneToRegionMap.get(geneList.get(geneIndex)).get(COLLAPSE_ISOFORMS_KEY));
		if (geneToCommentMap.containsKey(geneList.get(geneIndex)) && geneToCommentMap.get(geneList.get(geneIndex)) != null) {
			commentLabel.setText("gene #"+(geneIndex+1)+":  "+ geneToCommentMap.get(geneList.get(geneIndex)));
		} else {
			commentLabel.setText(" -- no comment -- ");
		}
		loadDataIfMissing();
	}
	
	private ArrayList<VariantContextWithFile> getExonVCFRecords(int exonIndex) {
	    String gene = geneList.get(geneIndex);
	    HashMap<String, ArrayList<ArrayList<VariantContextWithFile>>> isoformData = loadedVCFData.get(gene);
	    if (isoformData == null) {
	        return new ArrayList<VariantContextWithFile>();
	    }
	    String isoform = (String) isoformList.getSelectedItem();
	    ArrayList<ArrayList<VariantContextWithFile>> exonData = isoformData.get(isoform);
	    if (exonData == null || exonIndex >= exonData.size()) {
	        return new ArrayList<VariantContextWithFile>();
	    }
	    return exonData.get(exonIndex);
	}
	
	private void loadPopulationFile() throws IOException {
	    BufferedReader reader = Files.getAppropriateReader(popFile);
	    String[] header = reader.readLine().trim().split("\t", -1);
//	    colorCodeMap = new HashMap<String, Color>();
	    popColorMap = new HashMap<String, Color>();
	    popIndiMap = new HashMap<String, HashSet<String>>();
	    popIndiMapWithExcludes = new HashMap<String, HashSet<String>>();
	    
	    int iidCol = -1;
	    int exclCol = -1;
//	    int colorCol = -1;
	    int popCol = -1;
	    
	    int[] idCols = ext.indexFactors(Aliases.INDIVIDUAL_ID, header, false, false);
	    for (int i : idCols) {
	        if (i > -1) {
	            iidCol = i;
	            break;
	        }
	    }
	    System.out.println("Found ID column at position " + iidCol);
	    
	    for (int i = 0; i < header.length; i++) {
	        String hdr = header[i].toLowerCase();
	        if (hdr.startsWith("class=exclude")) {
	            exclCol = i;
	        } else if (hdr.startsWith("class=pop")) {
	            popCol = i;
	        }
	    }

	    chkbxDisplayExcludes.setEnabled(exclCol != -1);

        popSet = new HashSet<String>();
        excluded = new HashSet<String>();
        popMap = new HashMap<String, String>();
        superPopMap = new HashMap<String, String>();
	    String line = null;
	    String[] parts;
	    while((line = reader.readLine()) != null) {
	        parts = line.trim().split("\t", -1);
	        if (exclCol != -1 && parts[exclCol].equals("1")) {
	            excluded.add(parts[iidCol]);
	        }
	        if (popCol != -1) {
	            popMap.put(parts[iidCol], parts[popCol]);
	            popSet.add(parts[popCol]);
	            if (exclCol == -1 || !parts[exclCol].equals("1")) {
    	            HashSet<String> indiSet = popIndiMap.get(parts[popCol]);
    	            if (indiSet == null) {
    	                indiSet = new HashSet<String>();
    	                popIndiMap.put(parts[popCol], indiSet);
    	            }
    	            indiSet.add(parts[iidCol]);
	            }
	            HashSet<String> indiSet = popIndiMapWithExcludes.get(parts[popCol]);
	            if (indiSet == null) {
	                indiSet = new HashSet<String>();
	                popIndiMapWithExcludes.put(parts[popCol], indiSet);
	            }
	            indiSet.add(parts[iidCol]);
	            
	        }
	    }
	    reader.close();
	    
        if (popCol != -1) {
            popColorMap = parseColors(header[popCol]);
        } else {
            int index = 0;
            for (String pop : popSet) {
                popColorMap.put(pop, getAColor(index++)[0]);
            }
        }
	}
	
	private HashMap<String, Color> parseColors(String colorHeader) {
	    HashMap<String, Color> colorMap = new HashMap<String, Color>();
	    String[] pts = colorHeader.split(";");
	    int overflow = 9;
        for (int i = 1; i < pts.length; i++) {
            Color c = null;
            String code = pts[i].split("=")[0];
            String col = pts[i].split("=")[1].toLowerCase();
            if (col.contains(",")) {
                String[] colPts = col.split(",");
                c = new Color(Integer.parseInt(colPts[0]), Integer.parseInt(colPts[1]), Integer.parseInt(colPts[2]));
            } else {
                try {
                    c = (Color)Class.forName("java.awt.Color").getField(col).get(null);
                } catch (IllegalAccessException e) {
                    System.err.println("Error - couldn't parse color " + col);
                } catch (IllegalArgumentException e) {
                    System.err.println("Error - couldn't parse color " + col);
                } catch (NoSuchFieldException e) {
                    System.err.println("Error - couldn't parse color " + col);
                } catch (SecurityException e) {
                    System.err.println("Error - couldn't parse color " + col);
                } catch (ClassNotFoundException e) {
                    System.err.println("Error - couldn't parse color " + col);
                }
            }
            colorMap.put(code, c == null ? getAColor(i * overflow)[0] : c); // TODO check overflow colors
        }
        return colorMap;
	}
	
	private class VariantContextWithFile {
	    public VariantContextWithFile(VariantContext vc2, String vcfFile) {
	        this.vc = vc2;
	        this.file = vcfFile;
	    }
        VariantContext vc;
	    String file;
	}
	
	private void loadDataIfMissing() {
	    String gene = geneList.get(geneIndex);
	    if (loadedVCFData.containsKey(gene)) {
	        return;
	    } else {
	        GeneData gd = getCurrentGeneData();
	        ArrayList<VariantContextWithFile> data = new ArrayList<VariantContextWithFile>();
	        for (String vcfFile : vcfFiles) {
    	        VCFFileReader vcfReader = new VCFFileReader(vcfFile, true);
    	        VCFHeader header = vcfReader.getFileHeader();
    	        vcfHeader = header;
        	    CloseableIterator<VariantContext> vcIter = vcfReader.query(gd.getChromosomeUCSC(), gd.getStart(), gd.getStop());
				System.out.println(gd.getChromosomeUCSC() + ":" + gd.getStart() + "-" + gd.getStop());
        	    while (vcIter.hasNext()) {
        	        VariantContext vc = vcIter.next();
        	        VariantContextWithFile vcwf = new VariantContextWithFile(vc, vcfFile);
        	        data.add(vcwf);
        	    }
                vcIter.close();
                vcfReader.close();
				System.out.println(data.size() + " for " + gene + " in " + vcfFile);
	        }
            loadedVCFData.put(gene, sortData(data));
	    }
	}
	
	private HashMap<String, ArrayList<ArrayList<VariantContextWithFile>>> sortData(ArrayList<VariantContextWithFile> vcfEntries) {	    
	    HashMap<String, ArrayList<ArrayList<VariantContextWithFile>>> isoformMapToVCFList = new HashMap<String, ArrayList<ArrayList<VariantContextWithFile>>>();
	    HashMap<String, Segment[]> isoToExonSegMap = geneToExonSegmentMap.get(geneList.get(geneIndex));
	    
	    for (Entry<String, Segment[]> isoSegEntry : isoToExonSegMap.entrySet()) {
	        String isoKey = isoSegEntry.getKey();
	        Segment[] isoExonSegs = isoSegEntry.getValue();
	        
	        ArrayList<ArrayList<VariantContextWithFile>> isoVCFLists = new ArrayList<ArrayList<VariantContextWithFile>>();
	        isoformMapToVCFList.put(isoKey, isoVCFLists);
	        
	        for (int i = 0; i < isoExonSegs.length; i++) {
	            isoVCFLists.add(i, filter(isoExonSegs[i], vcfEntries));
	        }
	    }
	    
	    return isoformMapToVCFList;
	}

	public static void main(String[] args) {
		// Project proj = new Project("D:/projects/poynter.properties", false);
		Project proj = new Project("C:/workspace/Genvisis/projects/OSv2_hg19.properties", false);
		proj.GENE_LIST_FILENAMES.setValue(new String[] { "N:/statgen/VariantMapper/test2/genes.txt" });
		//
		String[] vcfFiles = new String[] { "N:/statgen/VariantMapper/test2/OSTEO_OFF_INHERIT.maf_0.01.final.vcf.gz", "N:/statgen/VariantMapper/test2/OSTEO_OFF_INHERIT_CONTROL.maf_0.01.final.vcf.gz" };
		String popFile = "N:/statgen/VariantMapper/test2/OSTEO_OFF_INHERIT_ALL.vpop";
		new VariantViewer(proj, vcfFiles, popFile);

		proj.GENE_LIST_FILENAMES.setValue(new String[] { "N:/statgen/VariantMapper/test3/genes.txt" });
		String[] vcfFiles2 = new String[] { "N:/statgen/VariantMapper/test3/CUSHINGS_TUMOR.maf_0.001.final.vcf.gz", "N:/statgen/VariantMapper/test3/CUSHINGS_TUMOR_CONTROL.maf_0.001.final.CUSHING_FREQ.vcf.gz" };
		String popFile2 = "N:/statgen/VariantMapper/test3/TN.vpop";
		new VariantViewer(proj, vcfFiles2, popFile2);

	}
}


