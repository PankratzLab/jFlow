package cnv.plots;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.GridLayout;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Vector;

import javax.imageio.ImageIO;
import javax.swing.AbstractAction;
import javax.swing.ActionMap;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.DefaultComboBoxModel;
import javax.swing.DefaultListCellRenderer;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.KeyStroke;
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
import common.HashVec;
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
	public static final int MIN_BUFFER = 10000;
	public static final int DEFAULT_STARTX = 20;
	public static final int DEFAULT_STARTY = 20;
	public static final int DEFAULT_WIDTH = 1100;
	public static final int DEFAULT_HEIGHT = 720;

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
	private String[][] geneRegions;
	private int geneIndex;
	private int isoformIndex;
	private JPanel lrrPanel;
	private JPanel bafPanel;
	private JPanel cnvPanel;
	private GeneTrack track;
	GeneData[][] geneData;
	private JLabel commentLabel;
	private String isoform;
	
	private Hashtable<String, String> namePathMap;
	private Logger log;
	private boolean fail;
	private JMenu loadRecentFileMenu;
	private ButtonGroup regionButtonGroup;

	private AbstractAction markerFileSelectAction = new AbstractAction() {
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
	private JComboBox<String> geneList;

    private PreparedMarkerSet markerSet;
    private long fingerprint;
    private String[] markerNames;
	
	public TrailerClone(Project proj, String[] vcfFiles) {
		this(proj, proj.DATA_DIRECTORY.getValue()  + "genes.txt"/*proj.GENE_LIST_FILENAMES.getValue()[0]*/, vcfFiles);
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
					String[] currSet = TrailerClone.this.proj.INDIVIDUAL_CNV_LIST_FILENAMES.getValue(); // TODO gene regions file property
					
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
//				            TrailerClone.this.proj.INDIVIDUAL_CNV_LIST_FILENAMES.setValue(newList);
//				            TrailerClone.this.proj.saveProperties();
				        }
					}
				}
				super.windowClosing(e);
			}
		});
		
//		System.out.println("startX: "+startX+"\t startY: "+startY+"\t width: "+width+"\t height: "+height);

		long time;
		String trackFilename;

		this.proj = proj;
		this.log = proj.getLog();
		jar = proj.JAR_STATUS.getValue();
		fail = false;
		
		time = new Date().getTime();

//		chr = (byte)Positions.parseUCSClocation(location)[0]; 

		fail = !loadMarkers();
		if (fail) {
			return;
		}        
		trackFilename = proj.getGeneTrackFilename(false);
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
	
	private void loadGenes(String file) {
	    String[] genes = HashVec.loadFileToStringArray(file, false, new int[]{0}, false);
	    geneData = track.lookupAllGeneData(genes); 
	    geneRegions = new String[geneData.length][];
	    for (int g = 0; g < geneData.length; g++) {
	        // TODO deal with multiple geneDatas
	        // TODO deal with comments
	        geneRegions[g] = new String[]{geneData[g][0].getGeneName(), geneData[g][0].getUCSClocation(), ""};
	    }
        TrailerClone.this.geneFileName = file;
	}
	
	private void updateGeneList() {
	    FontMetrics fontMetrics;
        int maxWidth;

        fontMetrics = geneList.getFontMetrics(geneList.getFont());
        maxWidth = fontMetrics.stringWidth("----------");
        
        String[] geneNames = new String[geneData.length];
        for (int i = 0; i < geneData.length; i++) {
            geneNames[i] = geneData[i][0].getGeneName();
            maxWidth = Math.max(maxWidth, fontMetrics.stringWidth(geneNames[i]));
        }

        geneList.setModel(new DefaultComboBoxModel<String>(geneNames));
        geneList.setPreferredSize(new Dimension(maxWidth + 50, 30));
	}
	
	private void paintLRRPanel(Graphics g) {
		// TODO moving paintComponent code here breaks drawing, and I haven't figured out why.. (cole - 3/6/15)
	}
	
	private void paintGeneTrackPanel(Graphics g) {
		GeneData[] genes;
		int[][] exons;
		Vector<Segment> v = new Vector<Segment>();
		Segment[] segs;
		int width, begin, end;
		
		if (stop-start > 10000000) {
			g.drawString("Zoom in to see genes", 10, 10);
		} else {
			genes = geneData[geneIndex];// track.getBetween(chr, start, stop, 30);
			g.setColor(Color.BLACK);
			for (int i = 0; i<genes.length; i++) {
				begin = getX(genes[i].getStart());
				end = getX(genes[i].getStop());
				g.drawRoundRect(begin, 0*15, end-begin, 10, 2, 2);
				v.add(new Segment(begin, end));
				exons = genes[i].getExonBoundaries();
				for (int j = 0; j < exons.length; j++) {
					begin = getX(exons[j][0]);
					end = getX(exons[j][1]);
					if (j==0 || j==exons.length-1) {
						g.fillRoundRect(begin, 0*15, end-begin+1, 10, 2, 2);
					} else {
						g.fillRect(begin, 0*15, end-begin+1, 10);
					}
					
				}
            }
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

				g.setFont(new Font("Arial", 0, 20));
				g.drawString("Log R Ratio", WIDTH_BUFFER, 20);
				
//				if (SHOW_MIDLINE) {
//					g.setColor(Color.LIGHT_GRAY);
//					g.drawLine(WIDTH_BUFFER, getHeight()-(int)((double)(0-min)/(double)(max-min)*(double)(getHeight()-2*HEIGHT_BUFFER))-HEIGHT_BUFFER, getWidth()-WIDTH_BUFFER, getHeight()-(int)((double)(0-min)/(double)(max-min)*(double)(getHeight()-2*HEIGHT_BUFFER))-HEIGHT_BUFFER);
//				}

				g.setFont(new Font("Arial", 0, 12));
			}
		};
		lrrPanel.addMouseListener(this);
		lrrPanel.addMouseMotionListener(this);
		lrrPanel.addMouseWheelListener(this);
		dataPanel.add(lrrPanel);

		cnvPanel = new JPanel() {
			public static final long serialVersionUID = 8L;
			public void paintComponent(Graphics g) {
			    paintGeneTrackPanel(g);
		    }
		};
		cnvPanel.setMaximumSize(new Dimension(getWidth(), 20));
		cnvPanel.addMouseListener(this);
		cnvPanel.addMouseMotionListener(this);
		cnvPanel.addMouseWheelListener(this);
		dataPanel.add(cnvPanel);

		bafPanel = new JPanel() {
			public static final long serialVersionUID = 7L;

			public void paintComponent(Graphics g) {
				g.setFont(new Font("Arial", 0, 20));
				g.drawString("B Allele Frequency", WIDTH_BUFFER, 20);

				g.setFont(new Font("Arial", 0, 10));
			}
		};
		dataPanel.add(bafPanel);

		getContentPane().add(dataPanel, BorderLayout.CENTER);

		JPanel sampPanel = new JPanel();
		((FlowLayout)sampPanel.getLayout()).setVgap(0);
		previousGene = new JButton(Grafik.getImageIcon("images/firstLast/Left.gif", true));
		previousGene.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dLeft.gif", true));
		previousGene.addActionListener(this);
		previousGene.setActionCommand(PREVIOUS_REGION);
		previousGene.setPreferredSize(new Dimension(25, 25));
		JPanel compPanel = new JPanel(new MigLayout("align center, fill, gap 0", "[grow, center]", "[][][]"));
		
		JPanel regionPanel = new JPanel();
		((FlowLayout)regionPanel.getLayout()).setVgap(0);
		geneList = new JComboBox<String>();
        DefaultListCellRenderer dlcr = new DefaultListCellRenderer();
        dlcr.setHorizontalAlignment(DefaultListCellRenderer.CENTER);
        geneList.setRenderer(dlcr);
        geneList.setBorder(BorderFactory.createEtchedBorder());
        geneList.setEditable(false);
		Font font = new Font("Arial", 0, 14);
		geneList.setFont(font);
		geneList.setAction(new AbstractAction() {
		    private static final long serialVersionUID = 1L;
		    public void actionPerformed(ActionEvent e) {
		        geneIndex = geneList.getSelectedIndex();
		        showGene(geneIndex);
		        updateGUI();
		    }
		});
		geneList.setPreferredSize(new Dimension(geneList.getPreferredSize().width, 26));
		nextGene = new JButton(Grafik.getImageIcon("images/firstLast/Right.gif", true));
        nextGene.setDisabledIcon(Grafik.getImageIcon("images/firstLast/dRight.gif", true));
        nextGene.addActionListener(this);
        nextGene.setActionCommand(NEXT_REGION);
        nextGene.setPreferredSize(new Dimension(25, 25));
        
        regionPanel.add(previousGene);
        regionPanel.add(geneList);
		regionPanel.add(nextGene);
		updateGeneList();
		
		compPanel.add(regionPanel, "cell 0 0");
		compPanel.setPreferredSize(new Dimension(compPanel.getPreferredSize().width, 95));
		
		commentLabel = new JLabel(" ", JLabel.CENTER);
		commentLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
		commentLabel.setFont(font);
		compPanel.add(commentLabel, "cell 0 1");

        JPanel descrPanel = new JPanel();
        descrPanel.setLayout(new MigLayout("gap 0", "[grow, center]", "[]0[]0[]"));
		descrPanel.add(compPanel, "cell 0 0");
		
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
                updateGUI();
            }
        });
		sampPanel.add(isoformList);
        sampPanel.setPreferredSize(new Dimension(sampPanel.getPreferredSize().width, isoformList.getPreferredSize().height + 5));
        descrPanel.add(sampPanel, "cell 0 1");

        
		JPanel overPanel = new JPanel();
		overPanel.setLayout(new BoxLayout(overPanel, BoxLayout.LINE_AXIS));
		
		JSeparator sep = new JSeparator(SwingConstants.VERTICAL);
		sep.setMaximumSize(new Dimension(1, 150));
		overPanel.setBorder(new EmptyBorder(5, 0, 5, 0));
		
		overPanel.add(Box.createHorizontalGlue());
		overPanel.add(descrPanel);
		overPanel.add(Box.createHorizontalGlue());
		
		getContentPane().add(overPanel, BorderLayout.NORTH);

		InputMap inputMap = bafPanel.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, InputEvent.ALT_MASK), PREVIOUS_REGION);
		inputMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, InputEvent.ALT_MASK), NEXT_REGION);

		ActionMap actionMap = bafPanel.getActionMap();
		actionMap.put(FIRST_REGION, new AbstractAction() {
			public static final long serialVersionUID = 9L;
			public void actionPerformed(ActionEvent e) {
				showGene(0);
			}
		});
		actionMap.put(PREVIOUS_REGION, new AbstractAction() {
			public static final long serialVersionUID = 10L;
			public void actionPerformed(ActionEvent e) {
				showGene(Math.max(geneIndex-1, 0));
			}
		});
		actionMap.put(NEXT_REGION, new AbstractAction() {
			public static final long serialVersionUID = 11L;
			public void actionPerformed(ActionEvent e) {
				showGene(Math.min(geneIndex+1, geneRegions.length-1));
			}
		});
		actionMap.put(LAST_REGION, new AbstractAction() {
			public static final long serialVersionUID = 12L;
			public void actionPerformed(ActionEvent e) {
				showGene(geneRegions.length - 1);
			}
		});
		bafPanel.setActionMap(actionMap);
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
	    int lW = lrrPanel.getWidth();
	    int bW = bafPanel.getWidth();
	    int cW = cnvPanel.getWidth();
	    int lH = lrrPanel.getHeight();
	    int bH = bafPanel.getHeight();
	    int cH = cnvPanel.getHeight();
	    BufferedImage imageLrr = new BufferedImage(lW, lH, BufferedImage.TYPE_INT_RGB);
	    BufferedImage imageBaf = new BufferedImage(bW, bH, BufferedImage.TYPE_INT_RGB);
	    BufferedImage imageCnv = new BufferedImage(cW, cH, BufferedImage.TYPE_INT_RGB);
	    
	    Graphics g = imageLrr.getGraphics();
        g.setColor(lrrPanel.getBackground());
        g.fillRect(0, 0, imageLrr.getWidth(), imageLrr.getHeight());
	    lrrPanel.paintAll(g);
	    
	    g = imageCnv.getGraphics();
        g.setColor(cnvPanel.getBackground());
        g.fillRect(0, 0, imageCnv.getWidth(), imageCnv.getHeight());
        cnvPanel.paintAll(g);
	    
	    g = imageBaf.getGraphics();
        g.setColor(bafPanel.getBackground());
        g.fillRect(0, 0, imageBaf.getWidth(), imageBaf.getHeight());
	    bafPanel.paintAll(g);
	    
	    int w = Math.max(lW, Math.max(bW, cW));
	    int h = lH + bH + cH + HEIGHT_BUFFER;
	    BufferedImage finalImage = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
	    
	    g = finalImage.getGraphics();
	    g.setColor(lrrPanel.getBackground());
	    g.fillRect(0, 0, finalImage.getWidth(), finalImage.getHeight());
	    g.drawImage(imageLrr, 0, 0, null);
	    g.drawImage(imageCnv, 0, lH + 5, null);
	    g.drawImage(imageBaf, 0, lH + cH + 10, null);
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
			String[] files = proj.INDIVIDUAL_CNV_LIST_FILENAMES.getValue();
			String name;
			for (String file : files) {
				name = ext.rootOf(file);
				regionFileNameLoc.put(name, file);
				JCheckBoxMenuItem menuItem = new JCheckBoxMenuItem();
				menuItem.setAction(markerFileSelectAction);
				menuItem.setFont(font);
				boolean found = Files.exists(file);
				menuItem.setText(name + (found ? "" : " -- [file not found]"));
				menuItem.setEnabled(found);
				regionFileNameBtn.put(name, menuItem);
				regionButtonGroup.add(menuItem);
				loadRecentFileMenu.add(menuItem);
			}
		}
		
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
			if (geneRegions == null || geneRegions.length == 0) {
				JOptionPane.showMessageDialog(null, "Error - No regions have been loaded", "Error", JOptionPane.ERROR_MESSAGE);
				return;
			}
			showGene(0);
		} else if (command.equals(PREVIOUS_REGION)) {
			if (geneRegions == null || geneRegions.length == 0) {
				JOptionPane.showMessageDialog(null, "Error - No regions have been loaded", "Error", JOptionPane.ERROR_MESSAGE);
				return;
			}
			showGene(Math.max(geneIndex-1, 0));
		} else if (command.equals(NEXT_REGION)) {
//			System.out.println("next");
			if (geneRegions == null || geneRegions.length == 0) {
//				filenames = proj.getIndividualRegionLists();
//				if (filenames.length == 0) {
//					JOptionPane.showMessageDialog(null, "Error - No regions have been loaded, since there no individual CNV region files defined in the properties file", "Error", JOptionPane.ERROR_MESSAGE);
					JOptionPane.showMessageDialog(null, "Error - No regions have been loaded", "Error", JOptionPane.ERROR_MESSAGE);
//				} else {
//					JOptionPane.showMessageDialog(null, "Error - No regions have been loaded; files include: "+Array.toStr(filenames, ", "), "Error", JOptionPane.ERROR_MESSAGE);
//				}
				return;
			}
			showGene(Math.min(geneIndex+1, geneRegions.length-1));
		} else if (command.equals(LAST_REGION)) {
			if (geneRegions == null || geneRegions.length == 0) {
				JOptionPane.showMessageDialog(null, "Error - No regions have been loaded", "Error", JOptionPane.ERROR_MESSAGE);
				return;
			}
			showGene(geneRegions.length - 1);
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

	public void mouseDragged(MouseEvent e) {
		int curX = e.getPoint().x;
		int distance = startX - curX;

		distance *= (stop - start) / (getWidth() - 2 * WIDTH_BUFFER);

		if (distance < 0) {
			distance = Math.max(distance, 1 - start);
		} else {
//			distance = Math.min(distance, positions[chrBoundaries[chr][1]]-stop);
			distance = Math.min(distance, (geneData[geneIndex][0].getStop() + MIN_BUFFER) - stop);
		}

//		if ((start<=1&&distance<0)||(stop>=positions[chrBoundaries[chr][1]]&&distance>0)) {
		if ((start <= geneData[geneIndex][0].getStart() - MIN_BUFFER && distance < 0) || (stop >= geneData[geneIndex][0].getStop() + MIN_BUFFER && distance > 0)) {

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
		
		boolean missing = geneData == null || geneData[geneIndex] == null || geneData[geneIndex].length == 0;
		String[] numbers = missing ? new String[0] : geneData[geneIndex][0].getNcbiAssessionNumbers();
		
		if (missing) {
			isoformsPresent = new String[] {"Gene Not Found"};
			maxWidth = Math.max(maxWidth, fontMetrics.stringWidth(isoformsPresent[0]));
		} else {
			isoformsPresent = new String[numbers.length + 1];
			for (int i = 0; i < numbers.length; i++) {
				isoformsPresent[i] = numbers[i];
				maxWidth = Math.max(maxWidth, fontMetrics.stringWidth(isoformsPresent[i]));
            }
			isoformsPresent[numbers.length] = collapse;
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
		start = loc[1];
		stop = loc[2];
		if (start==-1||start<0) {
			start = 1;
		}
		if (stop == -1 || stop > positions[chrBoundaries[chr][1]]) {
			stop = positions[chrBoundaries[chr][1]];
		}
		
		updateGUI();
	}

	public void updateGUI() {
		if (start <= geneData[geneIndex][0].getStart() - MIN_BUFFER) {
			start = geneData[geneIndex][0].getStart() - MIN_BUFFER;
		}
		startMarker = Array.binarySearch(positions, start, chrBoundaries[chr][0], chrBoundaries[chr][1], false);

		if (stop >= geneData[geneIndex][0].getStop() + MIN_BUFFER) {
			stop = geneData[geneIndex][0].getStop() + MIN_BUFFER;
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
	
//	public void loadRegions() {
//		BufferedReader reader;
//        Vector<String[]> v;
//        String line;
//        String[] parts;
//        int ignoredLines, countMissingRegions, invalidSamples;
//		
//		try {
//			String file = regionFileName.startsWith("./") ? proj.PROJECT_DIRECTORY.getValue() + regionFileName : regionFileName;
//			reader = Files.getAppropriateReader(file);//Files.getReader(file, jar, false, false);
//			System.out.print("Loading regions from " + regionFileName + "...");
//	        v = new Vector<String[]>();
//	        ignoredLines = countMissingRegions = invalidSamples = 0;
//	        line = null;
//            while ((line = reader.readLine()) != null) {
//            	parts = line.trim().split("\t");
//            	if (parts.length == 1) {
//            		v.add(new String[] {parts[0], "chr1"});
//            		countMissingRegions++;
//            	} else if (parts.length > 1 && parts[1].startsWith("chr")) {
//            		v.add(parts);
//            	} else {
//            		ignoredLines++;
//            	}
//            }
//            System.out.println(" loaded "+v.size()+" regions");
//            regions = Matrix.toStringArrays(v);
//            if (invalidSamples > 0) {
//            	JOptionPane.showMessageDialog(null, "Error - there were "+invalidSamples+" invalid samples in '"+regionFileName+"' that were ignored because they could not be found", "Error", JOptionPane.ERROR_MESSAGE);
//            }
//            if (countMissingRegions > 0) {
//            	JOptionPane.showMessageDialog(null, "Warning - there were "+countMissingRegions+" lines in '"+regionFileName+"' without a chromosomal region listed; using \"chr1\" for all missing values", "Warning", JOptionPane.ERROR_MESSAGE);
//            }
//            if (ignoredLines > 1) {
//            	JOptionPane.showMessageDialog(null, "Error - there were "+ignoredLines+" regions in '"+regionFileName+"' that were ignored due to improper formatting", "Error", JOptionPane.ERROR_MESSAGE);
//            }
//            reader.close();
//            
//            
//        } catch (FileNotFoundException fnfe) {
//            System.err.println("Error: file \""+regionFileName+"\" not found in data directory");
//        } catch (IOException ioe) {
//            System.err.println("Error reading file \""+regionFileName+"\"");
//        }
//	}
	
	public void showGene(int geneIndex) {
	    this.geneIndex = geneIndex;
		if (geneRegions == null || geneRegions.length == 0) {
			geneList.setSelectedIndex(geneIndex);
			commentLabel.setText(" ");
			return;
		}
		parseLocation(geneRegions[geneIndex][1]);
		if (geneRegions[geneIndex].length > 2) {
			commentLabel.setText("region #"+(geneIndex+1)+":  "+ geneRegions[geneIndex][2]);
		} else {
			commentLabel.setText(" -- no comment -- ");
		}
		geneList.setSelectedIndex(geneIndex);
		isoformIndex = 0;
		createIsoformList();
	}

	public static void main(String[] args) {
		Project proj = new Project("D:/projects/poynter.properties", false);

		new TrailerClone(proj, new String[]{});
	}
}


