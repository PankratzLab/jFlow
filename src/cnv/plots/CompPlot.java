/**
 * 
 */
package cnv.plots;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;

import javax.swing.AbstractAction;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ScrollPaneConstants;

import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.gui.ChromosomeViewer;
import cnv.gui.CompConfig;
import cnv.gui.FileNavigator;
import cnv.gui.LRRComp;
import cnv.gui.NewRegionListDialog;
import cnv.gui.RegionNavigator;
import cnv.manage.UCSCtrack;
import cnv.var.CNVRectangles;
import cnv.var.CNVariant;
import cnv.var.CNVariantHash;
import cnv.var.Region;

import common.Array;
import common.Files;
import common.Positions;
import common.ext;

import filesys.GeneTrack;

/**
 * @author Michael Vieths
 * 
 */
public class CompPlot extends JFrame {
	public static final long serialVersionUID = 1L;

	// public static final String DEFAULT_LOCATION = "chr17:55,609,472-55,824,368"; // USP32
	// public static final String DEFAULT_LOCATION = "chr17:15,609,472-40,824,368"; // USP32
	// public static final String DEFAULT_LOCATION = "chr6:161,590,461-163,364,497"; // PARK2
	public static final String DEFAULT_LOCATION = "chr6:161,624,000-163,776,000"; // PARK2 region
    private static final String REGION_LIST_NEW_FILE = "Load Region File";

	public static Color[] colorScheme = { Color.RED, Color.GREEN, Color.BLUE, Color.MAGENTA, Color.CYAN, Color.ORANGE, Color.YELLOW };

	Project proj;
	private String[] files;
	ArrayList<String> allFiles;
	ArrayList<String> filterFiles;
	GeneTrack track;

	// UI Components
	public JPanel compView;
	public RegionNavigator regionNavigator;
	public FileNavigator fileNavigator;
	public CompConfig compConfig;
	public ChromosomeViewer chromosomeViewer;
	public CompPropertyChangeListener cpcl;
	public CompPanel compPanel;

	// Variables configured via subpanels
	// From CompConfig
	int probes;
	int minSize;
	int qualityScore;
	int rectangleHeight;
	String displayMode;

	CNVRectangles cnvRects;

	// From RegionNavigator
	int[] location = new int[3];
	MarkerSet markerSet;
	int[] positions;
	String[] markerNames;
	boolean[] dropped;
	int[][] chrBoundaries;

	ArrayList<CNVariantHash> hashes;
    private JMenu delRegionFileMenu;
    private JMenu loadRecentFileMenu;
    private ButtonGroup regionButtonGroup;
    private HashMap<String, JCheckBoxMenuItem> regionFileNameBtn = new HashMap<String, JCheckBoxMenuItem>();
    private HashMap<String, String> regionFileNameLoc = new HashMap<String, String>();
    private String[] originalRegionFiles = null;
    
	public CompPlot(Project proj) {
		super("Genvisis - CompPlot - " + proj.PROJECT_NAME.getValue());
		this.proj = proj;
		this.markerSet = this.proj.getMarkerSet();
		if (markerSet != null) {
		    this.positions = this.markerSet.getPositions();
	        this.markerNames = markerSet.getMarkerNames();
	        Hashtable<String,String> hash = proj.getFilteredHash();
	        byte[] chrs = markerSet.getChrs();
		    dropped = new boolean[markerNames.length];
            chrBoundaries = new int[27][2];
            for (int i = 0; i < chrBoundaries.length; i++) {
                chrBoundaries[i][0] = chrBoundaries[i][1] = 0;
            }
            byte chr = 0;
            for (int i = 0; i < markerNames.length; i++) {
                dropped[i] = hash.containsKey(markerNames[i]);
                if (chrs[i] > chr) {
                    if (chr != 0) {
                        chrBoundaries[chr][1] = i - 1;
                    }
                    chr = chrs[i];
                    chrBoundaries[chr][0] = i;
                }
            }
            chrBoundaries[chr][1] = markerNames.length - 1;
            chrBoundaries[0][0] = 0;
            chrBoundaries[0][1] = markerNames.length-1;
		}
		init();
		
		addWindowListener(new WindowAdapter() {
            @Override
            public void windowClosing(WindowEvent e) {
                String[] curr = CompPlot.this.proj.REGION_LIST_FILENAMES.getValue();
                HashSet<String> currSet = new HashSet<String>();
                for (String s : curr) {
                    currSet.add(s);
                }
                for (String s : originalRegionFiles) {
                    currSet.remove(s);
                }
                if (currSet.size() > 0) {
                    String message = currSet.size() + " files have been added.  ";
                    int choice = JOptionPane.showOptionDialog(null, message+" Would you like to keep this configuration for the next time CompPlot is loaded?", "Preserve CompPlot workspace?", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, null, null);
                    if (choice == 0) {
                        CompPlot.this.proj.saveProperties();
                    }
                }
            }
		});
	}

	private void init() {
		// Position

		// Get a list of the .cnv files
		files = proj.CNV_FILENAMES.getValue();
		if (files.equals("")) {
			// CNV_FILENAMES is empty, throw an error and exit
			JOptionPane.showMessageDialog(null, "Error - CNV_FILENAMES property is empty");
			return;
		}
		allFiles = new ArrayList<String>();
		for (int i = 0; i < files.length; i++) {
			File file = new File(files[i]);
			String filename = file.getName();
			allFiles.add(filename);
		}
		
		originalRegionFiles = proj.REGION_LIST_FILENAMES.getValue();

		// Get the GeneTrack
		String geneTrackFile = proj.getGeneTrackFilename(false);
		if (geneTrackFile != null && !geneTrackFile.endsWith("/") && new File(geneTrackFile).exists()) {
			track = GeneTrack.load(geneTrackFile, false);
//		} else if (new File(GeneSet.REFSEQ_TRACK).exists()) {
//			track = GeneTrack.load(GeneSet.REFSEQ_TRACK, false);
		} else {
			JOptionPane.showMessageDialog(this, "Gene track is not installed. Gene boundaries will not be displayed.", "FYI", JOptionPane.INFORMATION_MESSAGE);
			track = null;
		}

		// Load the variants into memory
		hashes = new ArrayList<CNVariantHash>();
		filterFiles = new ArrayList<String>();
		for (String file : files) {
			File filename = new File(file);
			filterFiles.add(filename.getName());
			if (Files.exists(file, false)) {
				// Load the CNVs out of the files
				CNVariantHash cnvHash = CNVariantHash.load(file, CNVariantHash.CONSTRUCT_ALL, false);
				hashes.add(cnvHash);
			} else {
				JOptionPane.showMessageDialog(null, "Error - File " + file + " does not exist");
				return;
			}
		}

		setupGUI();
        this.setJMenuBar(createMenuBar());

		// Initialize the filter attributes
		probes = compConfig.getProbes();
		minSize = compConfig.getMinSize();
		qualityScore = compConfig.getQualityScore();
		rectangleHeight = compConfig.getRectangleHeight();
		setDisplayMode(compConfig.getDisplayMode());

		setRegion(regionNavigator.getRegion());
		
        JCheckBoxMenuItem jcbmi = regionFileNameBtn.get(ext.rootOf(CompPlot.this.regionNavigator.getRegionFile()));
        if (jcbmi != null) {
            jcbmi.setSelected(true);
        }
	}

	private void setupGUI() {
		// Set the default window size
		setSize(1000, 720);

		// Close this window but not the entire application on close
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		// Close the whole thing for debugging purposes
		// setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		cpcl = new CompPropertyChangeListener(this);

		// Create a new JPanel to contain everything
		compView = new JPanel();
		compView.setLayout(new BorderLayout());

		JPanel topPanel = new JPanel();
		topPanel.setLayout(new BoxLayout(topPanel, BoxLayout.Y_AXIS));

		regionNavigator = new RegionNavigator(this);
		regionNavigator.addPropertyChangeListener(cpcl);
		topPanel.add(regionNavigator);

		fileNavigator = new FileNavigator(files, colorScheme);
		fileNavigator.addPropertyChangeListener(cpcl);
		topPanel.add(fileNavigator);

		compView.add(topPanel, BorderLayout.PAGE_START);

		JPanel viewers = new JPanel();
		viewers.setLayout(new BorderLayout());

		chromosomeViewer = new ChromosomeViewer(location[0], location[1], location[2], track);

		viewers.add(chromosomeViewer, BorderLayout.NORTH);
		chromosomeViewer.setPreferredSize(new Dimension(800, 45));

		compPanel = new CompPanel(this);
		compPanel.addPropertyChangeListener(cpcl);
		compPanel.setChromosomeViewer(chromosomeViewer);

		JScrollPane jsp = new JScrollPane(compPanel);
		jsp.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
		viewers.add(jsp, BorderLayout.CENTER);

		compView.add(viewers);

		compConfig = new CompConfig(this);
		compConfig.addPropertyChangeListener(cpcl);

		compView.add(compConfig, BorderLayout.LINE_END);

		add(compView);

		// Set the panel visible
		setVisible(true);
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
                NewRegionListDialog newRgnList = new NewRegionListDialog(null, proj == null ? null : proj.getProperty(proj.PROJECT_DIRECTORY), false);
                newRgnList.setModal(true);
                newRgnList.setVisible(true);
                if (newRgnList.getReturnCode() == JOptionPane.YES_OPTION) {
                    final String rgnFile = newRgnList.getFileName();
                    addFileToList(rgnFile);
                    final JMenuItem remove = new JMenuItem();
                    remove.setAction(deleteFileAction);
                    remove.setText(rgnFile);
                    delRegionFileMenu.add(remove);
                    CompPlot.this.regionNavigator.loadRegions();
                    CompPlot.this.regionNavigator.setRegionFile(rgnFile);
                    CompPlot.this.regionNavigator.setRegion(0);
                    CompPlot.this.setRegion(CompPlot.this.regionNavigator.getRegion());
                }
            }
        });
        newRegionFile.setText("New Region List File");
        newRegionFile.setMnemonic(KeyEvent.VK_N);
        fileMenu.add(newRegionFile);

        delRegionFileMenu = new JMenu();
        delRegionFileMenu.setMnemonic(KeyEvent.VK_D);
        delRegionFileMenu.setText("Delete Region List File...");

        for (final String str : proj.REGION_LIST_FILENAMES.getValue()) {
            final JMenuItem remove = new JMenuItem();
            remove.setActionCommand(str);
            remove.setAction(deleteFileAction);
            remove.setText(str);
            delRegionFileMenu.add(remove);
        }

        fileMenu.add(delRegionFileMenu);

        JMenuItem loadRegionFile = new JMenuItem();
        loadRegionFile.setAction(loadNewFileAction);
        loadRegionFile.setText(REGION_LIST_NEW_FILE);
        loadRegionFile.setMnemonic(KeyEvent.VK_L);
        fileMenu.add(loadRegionFile);
        loadRecentFileMenu = new JMenu("Load Recent Region List...");
        loadRecentFileMenu.setMnemonic(KeyEvent.VK_R);
        fileMenu.add(loadRecentFileMenu);

        menuBar.add(fileMenu);

        regionButtonGroup = new ButtonGroup();
        if (proj != null) {
            String[] files = proj.REGION_LIST_FILENAMES.getValue();
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

        JMenu act = new JMenu("Actions");
        act.setMnemonic(KeyEvent.VK_A);
        menuBar.add(act);
        
        JMenuItem ucsc;
        JMenuItem bedUcsc;
        JMenuItem medianLRR;
        
        ucsc = new JMenuItem();
        ucsc.setAction(ucscAction);
        ucsc.setText("Open Region in UCSC");
        if (Desktop.isDesktopSupported()) {
            ucsc.setToolTipText("View this location on UCSC in a browser");
            ucsc.setEnabled(true);
        } else {
            ucsc.setToolTipText("Browser operations are not supported");
            ucsc.setEnabled(false);
        }
        act.add(ucsc);

        bedUcsc = new JMenuItem();
        bedUcsc.setAction(ucscBedAction);
        bedUcsc.setText("Upload to UCSC");
        if (Desktop.isDesktopSupported()) {
            bedUcsc.setToolTipText("Generate and upload a .BED file to UCSC");
            bedUcsc.setEnabled(true);
        } else {
            bedUcsc.setToolTipText("Browser operations are not supported");
            bedUcsc.setEnabled(false);
        }
        act.add(bedUcsc);
        
        medianLRR = new JMenuItem();
        medianLRR.setAction(lrrCompAction);
        medianLRR.setText("Median LRR");
        medianLRR.setToolTipText("Compute median Log R Ratios for a region");
        act.add(medianLRR);

        return menuBar;
    }
    
    private AbstractAction lrrCompAction = new AbstractAction() {
        private static final long serialVersionUID = 1L;
        @Override
        public void actionPerformed(ActionEvent e) {
            new Thread(new LRRComp(proj, regionNavigator.getTextField().getText())).start();
        }
    };
    
    private AbstractAction ucscBedAction = new AbstractAction() {
        private static final long serialVersionUID = 1L;
        @Override
        public void actionPerformed(ActionEvent e) {
            // Figure out which files are selected
            // Only allow upload if one file is selected (JDialog warning if multiples)
            ArrayList<String> files = getFilterFiles();
            if (files.size() != 1) {
                JOptionPane.showMessageDialog(null, "One and only one file must be selected before a .BED File can be generated", "Error", JOptionPane.ERROR_MESSAGE);
            } else {
                // Find the full path to the selected file
                String[] filePaths = proj.CNV_FILENAMES.getValue();
                String compressedFile = "";
                for (String file : filePaths) {
                    if (file.endsWith(files.get(0))) {
                        System.out.println("File path is " + file);
                        compressedFile = file + ".bed..gz"; // TODO this should be .bed.gz?
                        // Generate BED file with:
                        UCSCtrack.makeTrack(file, file, proj.getLog());
                        break;
                    }
                }

                // Direct the user to the BED upload page at UCSC Genome Browser
                Desktop desktop = Desktop.getDesktop();
                String URL = Positions.getUCSCUploadLink(Positions.parseUCSClocation(regionNavigator.getTextField().getText()), compressedFile);

                // UCSC uses chrX and chrY instead of 23 and 24
                URL = URL.replaceAll("chr23", "chrX");
                URL = URL.replaceAll("chr24", "chrY");
                try {
                    URI uri = new URI(URL);
                    System.out.println("Browsing to " + URL);
                    desktop.browse(uri);
                } catch (URISyntaxException e1) {
                    e1.printStackTrace();
                } catch (IOException e1) {
                    e1.printStackTrace();
                }
            }
        }
    };
    
    private AbstractAction ucscAction = new AbstractAction() {
        private static final long serialVersionUID = 1L;
        @Override
        public void actionPerformed(ActionEvent arg0) {
            Desktop desktop = Desktop.getDesktop();
            String URL = Positions.getUCSClink(Positions.parseUCSClocation(regionNavigator.getTextField().getText()));

            // UCSC uses chrX and chrY instead of 23 and 24
            URL = URL.replaceAll("chr23", "chrX");
            URL = URL.replaceAll("chr24", "chrY");
            try {
                URI uri = new URI(URL);
                System.out.println("Browsing to " + URL);
                desktop.browse(uri);
            } catch (URISyntaxException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }            
        }
    };
    
    private AbstractAction deleteFileAction = new AbstractAction() {
        private static final long serialVersionUID = 1L;
        @Override
        public void actionPerformed(ActionEvent e) {
            int opt = JOptionPane.showConfirmDialog(CompPlot.this, "Delete file from disk?", "Delete region file?", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null);
            switch(opt) {
                case JOptionPane.CANCEL_OPTION:
                    return;
                case JOptionPane.YES_OPTION:
                    boolean deleted = (new File(e.getActionCommand())).delete();
                    if (!deleted) {
                        JOptionPane.showMessageDialog(CompPlot.this, "Error - failed to delete file {" + e.getActionCommand() + "}", "Delete File Failed...", JOptionPane.ERROR_MESSAGE);
                    }
                    break;
                case JOptionPane.NO_OPTION:
                    break;
            }
            proj.REGION_LIST_FILENAMES.removeValue(e.getActionCommand());
//            CompPlot.this.regionNavigator.loadRegions();
            String[] val = proj.REGION_LIST_FILENAMES.getValue();
            CompPlot.this.regionNavigator.setRegionFile(val.length > 0 ? val[0] : "");
            CompPlot.this.regionNavigator.setRegion(0);
            CompPlot.this.setRegion(CompPlot.this.regionNavigator.getRegion());
            delRegionFileMenu.remove((JMenuItem) e.getSource());
            loadRecentFileMenu.remove(regionFileNameBtn.remove(ext.rootOf(e.getActionCommand())));
        }
    };
    
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
        
        final JMenuItem remove = new JMenuItem();
        remove.setActionCommand(file);
        remove.setAction(deleteFileAction);
        remove.setText(file);
        delRegionFileMenu.add(remove);
        
        proj.REGION_LIST_FILENAMES.addValue(file);
    }
    

    private AbstractAction markerFileSelectAction = new AbstractAction() {
        private static final long serialVersionUID = 1L;
        @Override
        public void actionPerformed(ActionEvent e) {
            String shortName = ((JCheckBoxMenuItem)e.getSource()).getText();
//            if (!loadingFile) {
                String file = regionFileNameLoc.get(shortName);
                if (file != null && file.equals(CompPlot.this.regionNavigator.getRegionFile())) {
                    return;
                }
                String tempFile = file.startsWith("./") ? proj.PROJECT_DIRECTORY.getValue() + file : file;
                if (!Files.exists(tempFile)) {
                    proj.message("Error - region file '" + shortName + "' doesn't exist.");
                    regionFileNameBtn.get(shortName).setSelected(true);
                } else {
//                    proj.REGION_LIST_FILENAMES.setValue(Array.insertStringAt(file, proj.REGION_LIST_FILENAMES.getValue(), 0));
//                    CompPlot.this.regionNavigator.loadRegions();
                    CompPlot.this.regionNavigator.setRegionFile(file);
                    CompPlot.this.regionNavigator.setRegion(0);
                    CompPlot.this.setRegion(CompPlot.this.regionNavigator.getRegion());
//                    regionIndex = 0;
//                    showRegion();
                }
            /*} else if (loadingFile) {
                // leave as currently selected marker
                if (CompPlot.this.regionNavigator.getRegionFile() != "" && CompPlot.this.regionNavigator.getRegionFile() != null) {
                    String file = ext.rootOf(CompPlot.this.regionNavigator.getRegionFile());
                    regionFileNameBtn.get(file).setSelected(true);
                }
                return;
            } */
        }
    };
    
    AbstractAction loadNewFileAction = new AbstractAction() {
        private static final long serialVersionUID = 1L;
        @Override
        public void actionPerformed(ActionEvent e) {
            String newFile = chooseNewFiles();
            if (newFile != null) {
                CompPlot.this.regionNavigator.loadRegions();
                CompPlot.this.regionNavigator.setRegionFile(newFile);
                CompPlot.this.regionNavigator.setRegion(0);
                regionFileNameBtn.get(ext.rootOf(CompPlot.this.regionNavigator.getRegionFile())).setSelected(true);
                CompPlot.this.setRegion(CompPlot.this.regionNavigator.getRegion());
            }
        }
    };

   
    private String chooseNewFiles() {
        JFileChooser jfc = new JFileChooser((proj != null /*|| regionFileName == null */? proj.PROJECT_DIRECTORY.getValue() : null /* ext.parseDirectoryOfFile(regionFileName)*/));
        jfc.setMultiSelectionEnabled(true);
        if (jfc.showOpenDialog(CompPlot.this) == JFileChooser.APPROVE_OPTION) {
            File[] files = jfc.getSelectedFiles();
            if (files.length > 0) {
                boolean[] keep = Array.booleanArray(files.length, true);
                for (int i = 0; i < files.length; i++) {
                    for (String fileName : proj.REGION_LIST_FILENAMES.getValue()) {
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
                    JOptionPane.showMessageDialog(CompPlot.this, msg.toString()); 
                }
                
                for (File kept : keptFiles) {
                    if (verifyValidFile(kept.getAbsolutePath())) {
                        addFileToList(kept.getAbsolutePath());
                    } else {
                        proj.getLog().reportError("Error - contents of file {" + kept.getAbsolutePath() + "} are not valid UCSC regions");
                    }
                }
                return keptFiles[0].getAbsolutePath();
            } else {
                File file = jfc.getSelectedFile();
                boolean keep = true;
                for (String fileName : proj.REGION_LIST_FILENAMES.getValue()) {
                    if (ext.rootOf(file.toString()).equals(fileName)) {
                        keep = false;
                    }
                }
                
                if (!keep) {
                    StringBuilder msg = new StringBuilder("The following data file is already present:\n").append(file.getName());
                    JOptionPane.showMessageDialog(CompPlot.this, msg.toString()); 
                } else {
                    if (verifyValidFile(file.getAbsolutePath())) {
                        addFileToList(file.getAbsolutePath());
                        return file.getAbsolutePath();
                    } else {
                        proj.getLog().reportError("Error - contents of file {" + file.getAbsolutePath() + "} are not valid UCSC regions");
                        return null;
                    }
                }
                
            }
            
        }
        return null;
    }
	    
	private boolean verifyValidFile(String file) {
	    File f = new File(file);
	    if (!f.exists()) return false;
	    try {
	        BufferedReader reader = Files.getAppropriateReader(file);
	        String line = null;
            while ((line = reader.readLine()) != null) {
                String[] parts = line.split("\t");
                int[] pos = Positions.parseUCSClocation(parts[0]);
                if (pos == null || pos[0] == -1 || pos[1] == -1 || pos[2] == -1) {
                    return false;
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
            return false;
        }
	    return true;
    }

    public void loadCNVs(int[] location) {
		// long startTime = Calendar.getInstance().getTimeInMillis();

		cnvRects = new CNVRectangles(hashes, allFiles, filterFiles, location, probes, minSize, qualityScore);
		cnvRects.setRectangleHeight(rectangleHeight);
		compPanel.setWindow(location[1], location[2]);
		cnvRects.setScalingFactor(compPanel.getScalingFactor());
		compPanel.setCNVRectangles(cnvRects);

		// long stopTime = Calendar.getInstance().getTimeInMillis();
		//
		// System.out.println("loadCNVs() took " + (stopTime - startTime) + "ms");
	}

	public String[] getFiles() {
		return files;
	}

	public void setFiles(String[] files) {
		this.files = files;
	}

	/*
	 * Methods to set values pulled from CompConfig
	 */
	public void setProbes(int p) {
		probes = p;
		loadCNVs(location);
	}

	public void setMinSize(int ms) {
		minSize = ms;
		loadCNVs(location);
	}

	public void setQualityScore(int qs) {
		qualityScore = qs;
		loadCNVs(location);
	}

	public void setRectangleHeight(int height) {
		rectangleHeight = height;
		loadCNVs(location);
	}

	public void setDisplayMode(String dm) {
		displayMode = dm;
		compPanel.setDisplayMode(displayMode);
		loadCNVs(location);
	}

	public void setSelectedCNVs(ArrayList<CNVariant> cnvs) {
		compConfig.setSelectedCNVs(cnvs);
	}

	/*
	 * Method to set values pulled from RegionNavigator
	 */
	public void setRegion(Region region) {
		location = Positions.parseUCSClocation(region.getRegion());
		byte chr;
		int start, stop;

        chr = (byte)location[0];
        if (chr==-1) {
            return;
        }
        start = location[1];
        stop = location[2];
        if (start == -1 || start < 0) {
            start = 1;
        }
        if (stop == -1 || stop > positions[chrBoundaries[chr][1]]) {
            stop = positions[chrBoundaries[chr][1]];
        }
        regionNavigator.getTextField().setText(Positions.getUCSCformat(new int[]{chr, start, stop}));
		chromosomeViewer.updateView(chr, start, stop);
		loadCNVs(location);
		chromosomeViewer.repaint();
	}

	public Region getRegion() {
		return regionNavigator.getRegion();
	}

	public void setCPLocation(int[] location) {
		this.location = location;
		regionNavigator.setLocation(location);
		chromosomeViewer.updateView(location[0], location[1], location[2]);
		loadCNVs(location);
		chromosomeViewer.repaint();
	}

	public int[] getCPLocation() {
		return location;
	}

	public Project getProject() {
		return proj;
	}

	public void setFilter(ArrayList<String> files) {
		filterFiles = files;
		loadCNVs(location);
	}

	public ArrayList<String> getFilterFiles() {
		return filterFiles;
	}
}

class CompPropertyChangeListener implements PropertyChangeListener {
	CompPlot compPlot;

	public CompPropertyChangeListener(CompPlot cp) {
		compPlot = cp;
	}

	@Override
	public void propertyChange(PropertyChangeEvent pve) {
		String propertyName;

		propertyName = pve.getPropertyName();

		if (propertyName.equals("probes")) {
			compPlot.setProbes(Integer.parseInt(pve.getNewValue().toString()));
		} else if (propertyName.equals("minSize")) {
			compPlot.setMinSize(Integer.parseInt(pve.getNewValue().toString()));
		} else if (propertyName.equals("qualityScore")) {
			compPlot.setQualityScore(Integer.parseInt(pve.getNewValue().toString()));
		} else if (propertyName.equals("rectangleHeight")) {
			compPlot.setRectangleHeight(Integer.parseInt(pve.getNewValue().toString()));
		} else if (propertyName.equals("displayMode")) {
			compPlot.setDisplayMode((String) pve.getNewValue());
			// } else if (propertyName.equals("firstRegion")) {
			// } else if (propertyName.equals("previousRegion")) {
			// } else if (propertyName.equals("nextRegion")) {
			// } else if (propertyName.equals("lastRegion")) {
		} else if (propertyName.equals("location")) {
			compPlot.setRegion((Region) pve.getNewValue());
		} else if (propertyName.equals("selectedCNV")) {
			@SuppressWarnings("unchecked")
			ArrayList<CNVariant> cnvs = (ArrayList<CNVariant>) pve.getNewValue();
			compPlot.setSelectedCNVs(cnvs);
		} else if (propertyName.equals("files")) {
			@SuppressWarnings("unchecked")
			ArrayList<String> files = (ArrayList<String>) pve.getNewValue();
			compPlot.setFilter(files);
		} else {
			// System.out.println(pve.getPropertyName() + " changed from " + pve.getOldValue() + " to " + pve.getNewValue());
		}
	}
}
