/**
 * 
 */
package cnv.plots;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Rectangle;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Vector;

import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import cnv.filesys.Project;
import cnv.gui.ChromosomeViewer;
import cnv.gui.CompConfig;
import cnv.gui.FileNavigator;
import cnv.gui.RegionNavigator;
import cnv.var.CNVariant;
import cnv.var.CNVariantHash;

import common.Files;
import common.Matrix;
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

	public Color[] colorScheme = { Color.RED, Color.GREEN, Color.BLUE, Color.YELLOW, Color.CYAN, Color.MAGENTA, Color.ORANGE };

	Project proj;
	private String[] files;
	GeneTrack track;
	String[] regionsList; // List of the region files
	int regionsListIndex; // Current file in the regionslist
	int regionIndex; // The current region
	String[][] regions;

	// UI Components
	public JPanel compView;
	public RegionNavigator regionNavigator;
	public FileNavigator fileNavigator;
	public CompConfig compConfig;
	public ChromosomeViewer chromosomeViewer;
	public CompPropertyChangeListener cpcl;
	public CompPanel compPanel;
	// public JPanel compPanel;

	// Variables configured via subpanels
	// From CompConfig
	int probes;
	int minSize;
	int qualityScore;
	String displayMode;
	Vector<CNVRectangle> rectangles;

	// From RegionNavigator
	private String geneLocation;
	private int[] location = new int[3];

	public CompPlot(Project proj) {
		this.proj = proj;

		init();

		// setupGUI();
	}

	private void init() {
		// Position
		// Starting position is a defined constant
		geneLocation = DEFAULT_LOCATION;

		// Get a list of the .cnv files
		files = proj.getFilenames(Project.CNV_FILENAMES);
		System.out.println("There are " + files.length + " files");
		for (int i = 0; i < files.length; i++) {
			System.out.println("  " + files[i]);
		}

		// Get a list of the regions
		regionsList = proj.getFilenames(Project.REGION_LIST_FILENAMES);
		regionsListIndex = 0;
		if (regionsList.length > 0) {
			if (Files.exists(regionsList[regionsListIndex], proj.getJarStatus())) {
				loadRegions();
			} else {
				System.err.println("Error - couldn't find '" + regionsList[regionsListIndex] + "' in data directory; populating with CNVs of current subject");
			}
		}
		regionIndex = -1;

		location = Positions.parseUCSClocation(geneLocation);

		// Get the GeneTrack
		String geneTrackFile = proj.getFilename(Project.GENETRACK_FILENAME);
		if (new File(geneTrackFile).exists()) {
			track = GeneTrack.load(geneTrackFile, false);
		} else if (new File(GeneSet.REFSEQ_TRACK).exists()) {
			track = GeneTrack.load(GeneSet.REFSEQ_TRACK, false);
		} else {
			JOptionPane.showMessageDialog(this, "Gene track is not installed. Gene boundaries will not be displayed.", "FYI", JOptionPane.INFORMATION_MESSAGE);
			track = null;
		}

		// Parse out the location chromosome/start base/end base
		int[] location = Positions.parseUCSClocation(DEFAULT_LOCATION);
		rectangles = new Vector<CNVRectangle>();

		setupGUI();

		// Initialize the filter attributes
		probes = compConfig.getProbes();
		minSize = compConfig.getMinSize();
		qualityScore = compConfig.getQualityScore();

		loadCNVs(location);
	}

	private void setupGUI() {
		// Set the default window size
		setSize(1000, 720);

		// Close this window but not the entire application on close
		// setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		// Close the whole thing for debugging purposes
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		cpcl = new CompPropertyChangeListener(this);

		// Create a new JPanel to contain everything
		compView = new JPanel();
		compView.setLayout(new BorderLayout());

		JPanel topPanel = new JPanel();
		topPanel.setLayout(new BoxLayout(topPanel, BoxLayout.Y_AXIS));

		regionNavigator = new RegionNavigator(geneLocation);
		regionNavigator.addPropertyChangeListener(cpcl);
		topPanel.add(regionNavigator);

		fileNavigator = new FileNavigator(files);
		topPanel.add(fileNavigator);

		compView.add(topPanel, BorderLayout.PAGE_START);

		JPanel viewers = new JPanel();
		viewers.setLayout(new BorderLayout());

		chromosomeViewer = new ChromosomeViewer(location[0], location[1], location[2], track);
		viewers.add(chromosomeViewer, BorderLayout.NORTH);
		chromosomeViewer.setPreferredSize(new Dimension(800, 15));

		compPanel = new CompPanel(this);
		viewers.add(compPanel, BorderLayout.CENTER);

		compView.add(viewers);

		compConfig = new CompConfig();
		compConfig.addPropertyChangeListener(cpcl);

		compView.add(compConfig, BorderLayout.LINE_END);

		add(compView);

		// Set the panel visible
		setVisible(true);
	}

	public void loadCNVs(int[] location) {
		// Read the data from the CNV files
		CNVariantHash[] cnvhs = new CNVariantHash[files.length];
		CNVariant[][] variants = new CNVariant[files.length][];

		// Clear the rectangles list
		rectangles.clear();

		for (int i = 0; i < cnvhs.length; i++) {
			// Load the CNVs out of the files
			cnvhs[i] = CNVariantHash.load(files[i], CNVariantHash.CONSTRUCT_ALL, false);
			// All CNVs are loaded into an array, which is then stored in an array by file
			// Format is:
			// file1: CNV1..CNVN
			// fileN: CNV1..CNVN
			// All CNVs in each file, when they are rendered, will be a single color
			// Any CNVs where CNVariant.getCN() > 2 will be a different shade of that color
			variants[i] = cnvhs[i].getAllInRegion((byte) location[0], location[1], location[2], probes, minSize, qualityScore);
		}

		for (int i = 0; i < cnvhs.length; i++) {
			// System.out.println("=== variant[" + i + "]");
			for (int j = 0; j < variants[i].length; j++) {
				CNVariant x = variants[i][j];
				// System.out.println("    Variant" + j + " at " + x.getUCSClocation() + " has " + x.getCN() + " copies");

				int startX = (int) x.getStart() - location[1];
				int stopX = (int) x.getStop() - location[1];

				CNVRectangle cnvRect = new CNVRectangle(startX, (j * 10) + 50, stopX, (j * 10) + 70, (byte) 2, true, true, (byte) 2, (byte) 1);
				cnvRect.setCNV(x);
				// Modulus the scheme we choose so it will wrap around in the event of too many files instead of throwing an exception
				cnvRect.setCNVColor(colorScheme[i % colorScheme.length]);
				rectangles.add(cnvRect);
			}
		}

		if (rectangles != null) {
			compPanel.setRectangles(rectangles.toArray(new CNVRectangle[0]));
			compPanel.setWindow(location[2] - location[1]);
		}
	}

	public void loadRegions() {
		BufferedReader reader;
		Vector<String[]> v;
		String[] line;

		try {
			reader = Files.getReader(regionsList[regionsListIndex], proj.getJarStatus(), false, false);
			v = new Vector<String[]>();
			while (reader.ready()) {
				line = reader.readLine().trim().split("\t");
				if (line.length > 1 && line[1].startsWith("chr")) {
					v.add(line);
				}
			}
			regions = Matrix.toStringArrays(v);
			System.out.println("Loaded " + regionsList.length + " regions");
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + regionsList[regionsListIndex] + "\" not found in data directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + regionsList[regionsListIndex] + "\"");
			System.exit(2);
		}
	}

	public String[] getFiles() {
		return files;
	}

	public void setFiles(String[] files) {
		this.files = files;
	}

	public String getGeneLocation() {
		return geneLocation;
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

	public void setDisplayMode(String dm) {
		System.out.println("Changing display mode");
		displayMode = dm;
	}

	/*
	 * Method to set values pulled from RegionNavigator
	 */
	public void setLocation(String region) {
		System.out.println("Changing location");
		location = Positions.parseUCSClocation(region);
		chromosomeViewer.updateView(location[0], location[1], location[2]);
		loadCNVs(location);
	}
}

class CompPropertyChangeListener implements PropertyChangeListener {
	CompPlot compPlot;

	public CompPropertyChangeListener(CompPlot cp) {
		compPlot = cp;
	}

	public void propertyChange(PropertyChangeEvent pve) {
		switch (pve.getPropertyName()) {
		case "probes":
			compPlot.setProbes(Integer.parseInt(pve.getNewValue().toString()));
			break;
		case "minSize":
			compPlot.setMinSize(Integer.parseInt(pve.getNewValue().toString()));
			break;
		case "qualityScore":
			compPlot.setQualityScore(Integer.parseInt(pve.getNewValue().toString()));
			break;
		case "displayMode":
			compPlot.setDisplayMode((String) pve.getNewValue());
			break;
		case "firstRegion":
			compPlot.setLocation((String) pve.getNewValue());
			break;
		case "previousRegion":
			compPlot.setLocation((String) pve.getNewValue());
			break;
		case "nextRegion":
			compPlot.setLocation((String) pve.getNewValue());
			break;
		case "lastRegion":
			compPlot.setLocation((String) pve.getNewValue());
			break;
		default:
			System.out.println(pve.getPropertyName() + " changed from " + pve.getOldValue() + " to " + pve.getNewValue());
		}
	}
}

/**
 * CNVRectangle describes a rectangle to be rendered in CompPanel
 * 
 * @author Michael Vieths
 * 
 *         Contains the CNVariant and a color based on the file from which it came
 */
class CNVRectangle extends GenericRectangle {
	private CNVariant cnv;
	private Color CNVColor;
	private Rectangle rect;

	public CNVRectangle(float startX, float startY, float stopX, float stopY, byte thickness, boolean fill, boolean roundedCorners, byte color, byte layer) {
		super(startX, startY, stopX, stopY, thickness, fill, roundedCorners, color, layer);
	}

	public CNVariant getCNV() {
		return cnv;
	}

	public void setCNV(CNVariant cnv) {
		this.cnv = cnv;
	}

	public Color getCNVColor() {
		return CNVColor;
	}

	/**
	 * Sets the color to render for the CNV. Adjusts the brightness based on number of copies.
	 * 
	 * @param CNVColor
	 */
	public void setCNVColor(Color CNVColor) {
		int copies = cnv.getCN();
		// Need to adjust the brightness
		float[] hsbVals = Color.RGBtoHSB(CNVColor.getRed(), CNVColor.getGreen(), CNVColor.getBlue(), null);
		float newBrightness = hsbVals[2];

		if (copies > 2) {
			// It's a duplication, make it brighter
			newBrightness *= (copies - 2);

		} else if (copies == 1) {
			// It's a deletion, make it darker
			newBrightness *= 0.5f;
		} else if (copies == 0) {
			// No copies, make it much darker
			newBrightness *= 0.2f;
		} else {
			// Normal number of copies, no change in brightness
		}

		this.CNVColor = Color.getHSBColor(hsbVals[0], hsbVals[1], newBrightness);
	}

	public Rectangle getRect() {
		return rect;
	}

	public void setRect(int x, int y, int width, int height) {
		rect = new Rectangle(x, y, width, height);
	}
}