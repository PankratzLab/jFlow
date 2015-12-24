package cnv.plots;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Shape;
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
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;

import net.miginfocom.swing.MigLayout;
import cnv.filesys.MarkerSet.PreparedMarkerSet;
import cnv.filesys.Project;
import cnv.gui.NewRegionListDialog;
import cnv.var.Region;
import common.Array;
import common.Files;
import common.Grafik;
import common.Logger;
import common.Positions;
import common.TransferableImage;
import common.ext;
import filesys.GeneData;
import filesys.GeneTrack;
import filesys.Segment;

public class TrailerClone extends JFrame implements ActionListener, MouseListener, MouseMotionListener, MouseWheelListener {
	public static final long serialVersionUID = 1L;

	public static final String DEFAULT_LOCATION = "chr17:55,609,472-55,824,368"; // USP32

	public static final String DEFAULT_SAMPLE = null;
	public static final boolean SHOW_MIDLINE = true;

	public static final double MOUSE_WHEEL_MULTIPLIER = 0.5;
	public static final int WIDTH_BUFFER = 25;
	public static final int HEIGHT_BUFFER = 10;
	public static final int DYNAMIC_HEIGHT_LIMIT = 0;
	public static final int DOUBLE_CLICK_INTERVAL = 500;
	public static final int SIZE = 4;
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

    private static final int GENE_HEIGHT = 30;
    private static final int EQUALIZED_EXON_BP_LENGTH = 300;
    private static final int INTRON_PLACEHOLDER_BP_LENGTH = 25;
    private static final String EXON_PREFIX = "";
    private static final int Y_START = 3*15;
    private static final Color FILLED_EXON_COLOR = Color.GRAY;
    private static final int DATA_PNT_SIZE = 10;

    private volatile int intronBPWidth = INTRON_PLACEHOLDER_BP_LENGTH;
    private volatile boolean equalizeExonLength = false;
    private volatile boolean paintExonNumbers = false;
    private volatile boolean paintIntrons = false;
    private volatile boolean fillExons = true;
    private volatile boolean paintExonBoxes = true;
    private volatile boolean paintInternalLine = false;
    
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
	private int isoformIndex;
//	private JPanel lrrPanel;
//	private JPanel bafPanel;
	private JPanel genePanel;
	private GeneTrack track;
	HashMap<String, String> geneToCommentMap;
	HashMap<String, HashMap<String, String>> geneToRegionMap;
	HashMap<String, HashMap<String, Segment[]>> geneToSegmentMap;
	HashMap<String, HashMap<String, GeneData>> geneToIsoformMap;
	HashMap<String, HashMap<String, ArrayList<VariantContext>>> loadedVCFData;
	private JLabel commentLabel;
	private JTextField commentField;
	private String isoform;
	private ArrayList<String> geneList;
	
	private Hashtable<String, String> namePathMap;
	private Logger log;
	private boolean fail;
	private JMenu loadRecentFileMenu;
	private ButtonGroup regionButtonGroup;

	private AbstractAction geneFileSelectAction = new AbstractAction() {
        private static final long serialVersionUID = 1L;
        @Override
        public void actionPerformed(ActionEvent e) {
            String shortName = ((JCheckBoxMenuItem)e.getSource()).getText();
            if (!loadingFile 
                    && !REGION_LIST_NEW_FILE.equals(shortName) 
                    && !REGION_LIST_PLACEHOLDER.equals(shortName)) {
                String file = regionFileNameLoc.get(shortName);
                if (file != null && file.equals(TrailerClone.this.geneFileName)) {
                    return;
                }
                String tempFile = file.startsWith("./") ? proj.PROJECT_DIRECTORY.getValue() + file : file;
                if (!Files.exists(tempFile)) {
                    proj.message("Error - region file '" + shortName + "' doesn't exist.");
                    regionFileNameBtn.get(shortName).setSelected(true);
                } else {
                    TrailerClone.this.geneFileName = file;
                    loadGenes(TrailerClone.this.geneFileName);
                    showGene(0);
                }
            } /*else if (loadingFile && REGION_LIST_PLACEHOLDER.equals(shortName)) {
                // do nothing
            } */else if (loadingFile || REGION_LIST_PLACEHOLDER.equals(shortName)) {
                // leave as currently selected marker
                if (TrailerClone.this.geneFileName != "" && TrailerClone.this.geneFileName != null) {
                    String file = TrailerClone.this.geneFileName;
                    file = ext.rootOf(TrailerClone.this.geneFileName);
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
                if (TrailerClone.this.geneFileName != null && !"".equals(TrailerClone.this.geneFileName)) {
                    regionFileNameBtn.get(ext.rootOf(TrailerClone.this.geneFileName)).setSelected(true);
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
            int code = jfc.showSaveDialog(TrailerClone.this);
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
    private long fingerprint;
    private String[] markerNames;
	
	public TrailerClone(Project proj, String[] vcfFiles) {
		this(proj, proj.GENE_LIST_FILENAMES.getValue()[0], vcfFiles);
	}

	// TODO Trailer should have a createAndShowGUI, same as all the other plots, as opposed to being its own frame 
	public TrailerClone(Project proj, String file, String[] vcfs) {
		super("Genvisis - Trailer - " + proj.PROJECT_NAME.getValue());
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		addWindowListener(new WindowAdapter() {
			@Override
			public void windowClosing(WindowEvent e) {
				if (TrailerClone.this.proj != null) {
					ArrayList<String> files = new ArrayList<String>(regionFileNameLoc.values());
					String[] currSet = TrailerClone.this.proj.GENE_LIST_FILENAMES.getValue();
					
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
				            TrailerClone.this.proj.GENE_LIST_FILENAMES.setValue(newList);
				            TrailerClone.this.proj.saveProperties();
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
        loadGenes(file);
        
		generateComponents();
		this.setJMenuBar(createMenuBar());
		
		time = new Date().getTime();

		setBounds(startX, DEFAULT_STARTX, DEFAULT_WIDTH, DEFAULT_HEIGHT);
		setVisible(true);
		
		updateGUI();
		geneIndex = -1;
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
	    String[][] geneFile = readFile(file);
	    geneList = new ArrayList<String>();
	    geneToIsoformMap = new HashMap<String, HashMap<String, GeneData>>();
	    geneToRegionMap = new HashMap<String, HashMap<String, String>>();
	    geneToSegmentMap = new HashMap<String, HashMap<String, Segment[]>>();
	    geneToCommentMap = new HashMap<String, String>();
	    loadedVCFData = new HashMap<String, HashMap<String,ArrayList<VariantContext>>>();
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
                    isoMap.put(geneData[i][g].isCollapsedIsoforms() ? "Collapse Isoforms" : geneData[i][g].getNcbiAssessionNumbers()[0], geneData[i][g]);
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
    	                isoformMap.put("Collapse Isoforms", geneData[i][g]);
    	            }
    	            isoformMap.put(geneData[i][g].getNcbiAssessionNumbers()[0], geneData[i][g]);
    	        }
	        }
	        if (geneFile[i].length >= 2) {
	            geneToCommentMap.put(geneFile[i][0], geneFile[i][1]);
	        }
	    }
	    geneToRegionMap = new HashMap<String, HashMap<String, String>>();
	    geneToSegmentMap = new HashMap<String, HashMap<String,Segment[]>>();
	    for (String gene : geneList) {
	        HashMap<String, GeneData> isoMap = geneToIsoformMap.get(gene);
	        HashMap<String, String> isoPosMap = new HashMap<String, String>();
	        HashMap<String, Segment[]> isoSegMap = new HashMap<String, Segment[]>();
	        geneToRegionMap.put(gene, isoPosMap);
	        geneToSegmentMap.put(gene, isoSegMap);
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
        TrailerClone.this.geneFileName = file;
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
	
	private ArrayList<VariantContext> filter(Segment exon, ArrayList<VariantContext> all) {
	    ArrayList<VariantContext> retArr = new ArrayList<VariantContext>();
	    for (VariantContext vc : all) {
	        if (exon.overlaps(new Segment(vc.getContig(), vc.getStart(), vc.getEnd()))) {
	            retArr.add(vc);
	        }
	    }
	    return retArr;
	}
	
	private void paintGeneTrackPanel(Graphics g) {
		GeneData gene;
		ArrayList<VariantContext> vcfData;
		ArrayList<VariantContext> vcfInSeg;
		int[][] exons;
		int width, begin, tempX, tempPx, len, lenPx, height;
		
		if (stop-start > 10000000) {
			g.drawString("Zoom in to see genes", 10, 10);
		} else {
		    height = genePanel.getHeight() - Y_START - GENE_HEIGHT;
			gene = geneToIsoformMap.get(geneList.get(geneIndex)).get(isoformList.getSelectedItem());
			if (loadedVCFData.get(geneList.get(geneIndex)) != null) {
			    vcfData = loadedVCFData.get(geneList.get(geneIndex)).get(isoformList.getSelectedItem());
			} else {
			    vcfData = new ArrayList<VariantContext>();
			}
			g.setColor(Color.BLACK);
			
            g.setFont(new Font("Arial", 0, 12));
	        tempX = begin = gene.getStart();
	        exons = gene.getExonBoundaries();
            for (int j = 0; j < exons.length; j++) {
                Segment exonSeg = new Segment(gene.getChr(), exons[j][0], exons[j][1]);
                vcfInSeg = filter(exonSeg, vcfData);
                tempPx = getX(tempX);
                len = equalizeExonLength ? EQUALIZED_EXON_BP_LENGTH : exons[j][1] - exons[j][0];
                lenPx = getX(tempX + len) - tempPx;
                // vertical line: 
                g.fillRect(tempPx, height, 1, GENE_HEIGHT);
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
                    width = g.getFontMetrics().stringWidth(EXON_PREFIX + (j + 1));
                    if (width < lenPx - 2) {
                        if (fillExons) {
                            g.setColor(Color.WHITE);
                        }
                        g.drawString(EXON_PREFIX + (j + 1), tempPx + 2, height + GENE_HEIGHT / 2 + g.getFontMetrics().getHeight() - 2);
                        if (fillExons) {
                            g.setColor(Color.BLACK);
                        }
                    }
                }
                if (vcfInSeg.size() > 0) {
                    if (lenPx <= DATA_PNT_SIZE + 2) {
                        g.setColor(Color.RED);
                        g.drawLine(tempPx + (lenPx / 2), height - 30, tempPx + (lenPx / 2), height - 5);
                        g.setColor(Color.BLACK);
                    } else {
                        // draw all in rows, ignoring actual position: 
//                        int numRows = (int) Math.ceil((vcfInSeg.size() * (PNT_SIZE + 2)) / (double)px);
//                        int numInRow = (int) Math.floor(px / (double)PNT_SIZE);
//                        int index = 0;
//                        g.setColor(Color.RED);
//                        for (int r = 0; index < vcfInSeg.size() && r < numRows; r++) {
//                            for (int n = 0; index < vcfInSeg.size() && n < numInRow; n++) {
//                                g.fillOval(getX(tempX) + 1 + n * (PNT_SIZE + 2), height - ((PNT_SIZE + 2) * (r + 1)), PNT_SIZE, PNT_SIZE);
//                                index++;
//                            }
//                        }
//                        g.setColor(Color.BLACK);
                        
                        // draw all in relative position, pushing up if overlapping
                        g.setColor(Color.RED);
                        
                        ArrayList<Rectangle> plotted = new ArrayList<Rectangle>();
                        for (VariantContext vc : vcfInSeg) {
                            int diffA = vc.getStart() - exons[j][0];
                            double prop = diffA / (double)len;
                            int diffPx = (int) (lenPx * prop);
                            if (diffPx + DATA_PNT_SIZE + 2 > lenPx) { // don't let the markings go past the exon
                                diffPx = lenPx - DATA_PNT_SIZE - 2;
                            }
                            Rectangle vcRect = new Rectangle(tempPx + diffPx, 0, DATA_PNT_SIZE + 2, DATA_PNT_SIZE + 2);
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
                                    vcRect = new Rectangle(vcRect.x, vcRect.y + vcRect.height, DATA_PNT_SIZE + 2, DATA_PNT_SIZE + 2);
                                }
                            } while (overlap);
                            plotted.add(vcRect);
                        }
                        for (Rectangle rect : plotted) {
                            g.fillOval(rect.x, height - DATA_PNT_SIZE - 2 - rect.y, DATA_PNT_SIZE, DATA_PNT_SIZE);
                        }
                        g.setColor(Color.BLACK);
                        
                    }
                }
                // move tempX to other side of exon:
                tempX += len;
                if (!fillExons && !paintExonBoxes) {
                    g.fillRect(getX(tempX)-1, height, 1, GENE_HEIGHT);
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
		}
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
		if (jfc.showOpenDialog(TrailerClone.this) == JFileChooser.APPROVE_OPTION) {
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
					JOptionPane.showMessageDialog(TrailerClone.this, msg.toString()); 
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
					JOptionPane.showMessageDialog(TrailerClone.this, msg.toString()); 
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
			    paintGeneTrackPanel(g);
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
		        geneIndex = geneListCmb.getSelectedIndex();
		        showGene(geneIndex);
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
		updateGeneList();
		
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
		createIsoformList();
		
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
//                    createIsoformList(); // TODO collapse isoforms
                    isoformIndex = COLLAPSE_ISOFORMS;
                } else if (!isoformsPresent[index].equals(isoform)) {
                    isoformIndex = index;
//                    updateSample(isoformsPresent[index]);
                }
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
		JCheckBoxMenuItem paintIntronsChk = new JCheckBoxMenuItem();
		disp.add(paintIntronsChk);
		JCheckBoxMenuItem paintExonBoxesChk = new JCheckBoxMenuItem();
		disp.add(paintExonBoxesChk);
		JCheckBoxMenuItem fillExonsChk = new JCheckBoxMenuItem();
		disp.add(fillExonsChk);
		JCheckBoxMenuItem paintExonLinesChk = new JCheckBoxMenuItem();
		disp.add(paintExonLinesChk);
		JCheckBoxMenuItem equalizeExonsChk = new JCheckBoxMenuItem();
        disp.add(equalizeExonsChk);
        JCheckBoxMenuItem paintExonNumbersChk = new JCheckBoxMenuItem();
        disp.add(paintExonNumbersChk);
        
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
	
	public void mouseClicked(MouseEvent e) {
        if (e.getButton()==MouseEvent.BUTTON1) {
            zoomProportionally(false, e.getPoint(), true);
        } else if (e.getButton()==MouseEvent.BUTTON3) {
            zoom(1, 1);
        }
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
	
	private int getBufferedStart() {
//	    return geneData[geneIndex][0].getStart() - MIN_BUFFER;
	    return geneToIsoformMap.get(geneList.get(geneIndex)).get(isoformList.getSelectedItem()).getStart() - MIN_BUFFER;
	}

	private int getBufferedStop() {
//	    int stop = geneData[geneIndex][0].getStart();
        int stop = geneToIsoformMap.get(geneList.get(geneIndex)).get(isoformList.getSelectedItem()).getStart();
        if (paintIntrons) {
//            stop = geneData[geneIndex][0].getStop();
            stop = geneToIsoformMap.get(geneList.get(geneIndex)).get(isoformList.getSelectedItem()).getStop();
        } else {
//            int[][] exons = geneData[geneIndex][0].getExonBoundaries();
            int[][] exons = geneToIsoformMap.get(geneList.get(geneIndex)).get(isoformList.getSelectedItem()).getExonBoundaries();
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
			distance = Math.min(distance, getBufferedStop() - stop);
		}

	    if ((start <= getBufferedStart() && distance < 0) || (stop >= getBufferedStop() && distance > 0)) {

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
        fingerprint = markerSet.getFingerprint();
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
		String collapse;
		int maxWidth;

		fontMetrics = isoformList.getFontMetrics(isoformList.getFont());
		collapse = "Collapse Isoforms";
		maxWidth = fontMetrics.stringWidth(collapse);
		
		HashMap<String, GeneData> isoMap = geneToIsoformMap.get(geneList.get(geneIndex));
		boolean missing = isoMap == null || isoMap.size() == 0;
		
		if (missing) {
			isoformsPresent = new String[] {"Gene Not Found"};
			maxWidth = Math.max(maxWidth, fontMetrics.stringWidth(isoformsPresent[0]));
		} else {
			isoformsPresent = new String[isoMap.size()];
			int index = 0;
			for (Entry<String, GeneData> isoEntry : isoMap.entrySet()) {
			    if (isoEntry.getKey().equals(collapse)) {
			        continue;
			    }
				isoformsPresent[index] = isoEntry.getKey();
				maxWidth = Math.max(maxWidth, fontMetrics.stringWidth(isoformsPresent[index++]));
            }
			isoformsPresent[index] = collapse;
		}
		
		isoformList.setModel(new DefaultComboBoxModel<String>(isoformsPresent));
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
		stop = getBufferedStop();
//		if (start==-1||start<0) {
//			start = 1;
//		}
//		if (stop == -1 || stop > positions[chrBoundaries[chr][1]]) {
//			stop = positions[chrBoundaries[chr][1]];
//		}
		
		updateGUI();
	};

	public void updateGUI() {
		if (start < getBufferedStart()) {
			start = getBufferedStart();
		}
		startMarker = Array.binarySearch(positions, start, chrBoundaries[chr][0], chrBoundaries[chr][1], false);

//		if (stop >= geneData[geneIndex][0].getStop() + MIN_BUFFER) {
//		    stop = geneData[geneIndex][0].getStop() + MIN_BUFFER;
//		}
		if (stop > getBufferedStop()) {
			stop = getBufferedStop();
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
		isoformIndex = 0;
		isoformList.setSelectedIndex(isoformIndex);
		parseLocation(geneToRegionMap.get(geneList.get(geneIndex)).get(isoformList.getItemAt(0)));
		if (geneToCommentMap.containsKey(geneList.get(geneIndex)) && geneToCommentMap.get(geneList.get(geneIndex)) != null) {
			commentLabel.setText("gene #"+(geneIndex+1)+":  "+ geneToCommentMap.get(geneList.get(geneIndex)));
		} else {
			commentLabel.setText(" -- no comment -- ");
		}
		loadData();
	}
	
	private HashMap<String, ArrayList<VariantContext>> sortData(ArrayList<VariantContext> vcfEntries) {
	    ArrayList<Segment> vcfSegments = new ArrayList<Segment>();
	    for (VariantContext vc : vcfEntries) {
            Segment vcfSeg = new Segment(vc.getContig(), vc.getStart(), vc.getEnd());
            vcfSegments.add(vcfSeg);
        }
	    
	    HashMap<String, ArrayList<VariantContext>> isoformToVCFMap = new HashMap<String, ArrayList<VariantContext>>();
	    
	    HashMap<String, Segment[]> isoToSegMap = geneToSegmentMap.get(geneList.get(geneIndex));
	    for (Entry<String, Segment[]> isoSegEntry : isoToSegMap.entrySet()) {
	        String isoKey = isoSegEntry.getKey();
	        Segment[] isoSegs = isoSegEntry.getValue();
	        
	        ArrayList<VariantContext> isoVCFs = new ArrayList<VariantContext>();
	        isoformToVCFMap.put(isoKey, isoVCFs);
	        
	        for (int i = 0; i < vcfSegments.size(); i++) {
	            if (Segment.overlapsAny(vcfSegments.get(i), isoSegs)) {
	                isoVCFs.add(vcfEntries.get(i));
	            }
	        }
	    }
	    
	    return isoformToVCFMap;
	}
	
	private void loadData() {
	    String[] vcfFiles = new String[]{"N:/statgen/VariantMapper/OSTEO_OFF_INHERIT.final.vcf.gz", "N:/statgen/VariantMapper/OSTEO_OFF_INHERIT_CONTROL.final.vcf.gz"};
	    String gene = geneList.get(geneIndex);
	    if (loadedVCFData.containsKey(gene)) {
	        return;
	    } else {
	        ArrayList<VariantContext> data = new ArrayList<VariantContext>();
	        for (String vcfFile : vcfFiles) {
    	        VCFFileReader vcfReader = new VCFFileReader(vcfFile, true);
        	    CloseableIterator<VariantContext> vcIter = vcfReader.query("chr" + chr, start, stop);
        	    while (vcIter.hasNext()) {
        	        data.add(vcIter.next());
        	    }
                vcIter.close();
                vcfReader.close();
	        }
            loadedVCFData.put(gene, sortData(data));
	    }
	}
	
	public static void main(String[] args) {
		Project proj = new Project("D:/projects/poynter.properties", false);

		new TrailerClone(proj, new String[]{});
	}
}


