/**
 * 
 */
package cnv.plots;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;

import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ScrollPaneConstants;

import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.gui.ChromosomeViewer;
import cnv.gui.CompConfig;
import cnv.gui.FileNavigator;
import cnv.gui.RegionNavigator;
import cnv.var.CNVRectangles;
import cnv.var.CNVariant;
import cnv.var.CNVariantHash;
import cnv.var.Region;
import common.Files;
import common.Positions;
import filesys.GeneSet;
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

	public CompPlot(Project proj) {
		super("Genvisis - CompPlot - " + proj.getNameOfProject());
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

		// Get the GeneTrack
		String geneTrackFile = proj.GENETRACK_FILENAME.getValue(false, false);
		if (geneTrackFile != null && !geneTrackFile.endsWith("/") && new File(geneTrackFile).exists()) {
			track = GeneTrack.load(geneTrackFile, false);
		} else if (new File(GeneSet.REFSEQ_TRACK).exists()) {
			track = GeneTrack.load(GeneSet.REFSEQ_TRACK, false);
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

		// Initialize the filter attributes
		probes = compConfig.getProbes();
		minSize = compConfig.getMinSize();
		qualityScore = compConfig.getQualityScore();
		rectangleHeight = compConfig.getRectangleHeight();
		setDisplayMode(compConfig.getDisplayMode());

		setRegion(regionNavigator.getRegion());

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
