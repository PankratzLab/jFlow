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
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import cnv.filesys.Project;
import cnv.gui.ChromosomeViewer;
import cnv.gui.CompConfig;
import cnv.gui.FileNavigator;
import cnv.gui.RegionNavigator;
import cnv.var.CNVariant;
import cnv.var.CNVariantHash;
import cnv.var.Region;

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
	ArrayList<CNVRectangle> rectangles;

	// From RegionNavigator
	int[] location = new int[3];

	public CompPlot(Project proj) {
		this.proj = proj;

		init();

		// setupGUI();
	}

	private void init() {
		// Position

		// Get a list of the .cnv files
		files = proj.getFilenames(Project.CNV_FILENAMES);
		// for (int i = 0; i < files.length; i++) {
		// System.out.println("  " + files[i]);
		// }

		// Get the GeneTrack
		String geneTrackFile = proj.getFilename(Project.GENETRACK_FILENAME, false, false);
		if (geneTrackFile != null && !geneTrackFile.endsWith("/") && new File(geneTrackFile).exists()) {
			track = GeneTrack.load(geneTrackFile, false);
		} else if (new File(GeneSet.REFSEQ_TRACK).exists()) {
			track = GeneTrack.load(GeneSet.REFSEQ_TRACK, false);
		} else {
			JOptionPane.showMessageDialog(this, "Gene track is not installed. Gene boundaries will not be displayed.", "FYI", JOptionPane.INFORMATION_MESSAGE);
			track = null;
		}

		// Parse out the location chromosome/start base/end base
		rectangles = new ArrayList<CNVRectangle>();

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
		// setSize(1920, 1200);

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

		regionNavigator = new RegionNavigator(proj);
		regionNavigator.addPropertyChangeListener(cpcl);
		topPanel.add(regionNavigator);

		fileNavigator = new FileNavigator(files, colorScheme);
		topPanel.add(fileNavigator);

		compView.add(topPanel, BorderLayout.PAGE_START);

		JPanel viewers = new JPanel();
		viewers.setLayout(new BorderLayout());

		chromosomeViewer = new ChromosomeViewer(location[0], location[1], location[2], track);
		viewers.add(chromosomeViewer, BorderLayout.NORTH);
		chromosomeViewer.setPreferredSize(new Dimension(800, 15));

		compPanel = new CompPanel(this);
		compPanel.addPropertyChangeListener(cpcl);
		compPanel.setRectangles(rectangles.toArray(new CNVRectangle[0]));

		JScrollPane jsp = new JScrollPane(compPanel);
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
		// Read the data from the CNV files
		CNVariantHash[] cnvhs = new CNVariantHash[files.length];
		CNVariant[][] variants = new CNVariant[files.length][];

		// Clear the rectangles list
		rectangles.clear();
		compPanel.setRectangles(rectangles.toArray(new CNVRectangle[0]));

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

		// Create a hashmap of the CNVs with a Rectangle as the key and quantity as value?
		HashMap<String, CNVRectangle> cnvMap = new HashMap<String, CNVRectangle>();
		// Outer loop represents the list of files
		for (int i = 0; i < cnvhs.length; i++) {
			// Inner loop represents the CNVs within each file
			for (int j = 0; j < variants[i].length; j++) {
				CNVariant variant = variants[i][j];

				int startX = (int) variant.getStart() - location[1];
				int stopX = (int) variant.getStop() - location[1];
				CNVRectangle cnvRect = new CNVRectangle(startX, stopX, (byte) 2, true, true, (byte) 2, (byte) 1);
				cnvRect.addCNV(variant);
				// Modulus the scheme we choose so it will wrap around in the event of too many files instead of throwing an exception
				cnvRect.setCNVColor(colorScheme[i % colorScheme.length], displayMode);

				if (displayMode.equals("Collapsed")) {
					// Collapsed display (CNVs with the same start and length are shown as one rectangle)
					String ucscLocation = variant.getUCSClocation();
					if (!cnvMap.containsKey(ucscLocation)) {
						cnvMap.put(ucscLocation, cnvRect);
					} else {
						CNVRectangle rect = cnvMap.get(ucscLocation);
						rect.addCNV(variant);
					}
				} else {
					// Full or Packed display mode (one CNV per rectangle)
					rectangles.add(cnvRect);
				}
			}
		}

		// Fill in the collapsed rectangles
		if (displayMode.equals("Collapsed")) {
			for (CNVRectangle cnvRect : cnvMap.values()) {
				rectangles.add(cnvRect);
			}
		}

		if (rectangles.size() > 0) {

			Collections.sort(rectangles);
			// Set the preferred size of the window to be large enough to encompass all of the rectangles
			compPanel.setPreferredSize(new Dimension(800, (rectangles.size() * rectangleHeight) + rectangleHeight));
			compPanel.setRectangleHeight(rectangleHeight);
			compPanel.setRectangles(rectangles.toArray(new CNVRectangle[0]));
			compPanel.setWindow(location[1], location[2]);
		}
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
		chromosomeViewer.updateView(location[0], location[1], location[2]);
		loadCNVs(location);
		chromosomeViewer.repaint();
	}

	public Region getRegion() {
		return regionNavigator.getRegion();
	}

	public int[] getCPLocation() {
		return location;
	}

	public Project getProject() {
		return proj;
	}
}

class CompPropertyChangeListener implements PropertyChangeListener {
	CompPlot compPlot;

	public CompPropertyChangeListener(CompPlot cp) {
		compPlot = cp;
	}

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
		} else {
			// System.out.println(pve.getPropertyName() + " changed from " + pve.getOldValue() + " to " + pve.getNewValue());
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
class CNVRectangle extends GenericRectangle implements Comparable<CNVRectangle> {
	private ArrayList<CNVariant> cnvs;
	private Color CNVColor;
	private Rectangle rect;
	private boolean selected;
	private int quantity; // How many CNVs are represented by this rectangle
	private boolean inUse;

	public CNVRectangle(float startX, float stopX, byte thickness, boolean fill, boolean roundedCorners, byte color, byte layer) {
		// Y coord doesn't matter, that'll get set at render time
		super(startX, 0, stopX, 0, thickness, fill, roundedCorners, color, layer);
		quantity = 1;
		selected = false;
		inUse = false;
		cnvs = new ArrayList<CNVariant>();
	}

	public ArrayList<CNVariant> getCNVs() {
		return cnvs;
	}

	public CNVariant getCNV() {
		return cnvs.get(0);
	}

	public void setCNVs(ArrayList<CNVariant> variants) {
		cnvs = variants;
	}

	public void addCNV(CNVariant variant) {
		cnvs.add(variant);
	}

	public Color getCNVColor() {
		return CNVColor;
	}

	public int getQuantity() {
		return quantity;
	}

	public void setQuantity(int newQuantity) {
		quantity = newQuantity;
	}

	/**
	 * Sets the color to render for the CNV. Adjusts the brightness based on number of copies.
	 * 
	 * @param CNVColor
	 */
	public void setCNVColor(Color CNVColor, String displayMode) {
		// Default to 2 copies
		int copies = 2;

		// Only change the brightness if we're in Full or Compressed mode.
		// There's only one CNV in this CNVRectangle for Full and Pack, could be multiple if it's Compressed
		if (cnvs.size() > 0) {
			if (displayMode.equals("Full") || displayMode.equals("Pack")) {
				copies = cnvs.get(0).getCN();
			}
		}

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

	public void setSelected(boolean sel) {
		selected = sel;
	}

	public boolean isSelected() {
		return selected;
	}

	public void setUsed(boolean used) {
		inUse = used;
	}

	public boolean isInUse() {
		return inUse;
	}

	@Override
	// Allow sorting the entire list based first on start position, then on length
	public int compareTo(CNVRectangle o) {
		int retValue = 0;
		float start1 = getStartXValue();
		float start2 = o.getStartXValue();
		float length1 = getStopXValue() - start1;
		float length2 = o.getStopXValue() - start2;

		if (start1 > start2) {
			retValue = 1;
		} else if (start1 < start2) {
			retValue = -1;
		} else {
			// They start at the same spot, but their lengths may not be the same
			if (length1 > length2) {
				retValue = 1;
			} else if (length1 < length2) {
				retValue = -1;
			} else {
				retValue = 0;
			}
		}

		return retValue;
	}
}