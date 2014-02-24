package cnv.plots;

import java.io.*;
import java.util.*;
import common.*;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;
import java.awt.image.BufferedImage;
import javax.imageio.*;

import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.gui.LaunchAction;
import mining.Distance;

public class MosaicPanel extends AbstractPanel implements MouseListener, MouseMotionListener {
	public static final long serialVersionUID = 3L;
	public static final int HEAD_BUFFER = 25;
	public static final int HEIGHT_X_AXIS = 105;
	public static final int WIDTH_Y_AXIS = 125;
	public static final int WIDTH_BUFFER = 50;
	public static final int SIZE = 12;
	public static final int SIZE_FAILED = 6;
	public static final int AXIS_THICKNESS = 4;
	public static final int TICK_THICKNESS = 3;
	public static final int TICK_LENGTH = 15;
	public static final int MARKER_SIZE = 6;
	public static final double X_STEP = 0.05;
	public static final double Y_STEP = 0.05;
	public static final int LOOKUP_RESOLUTION = 20;
	public static final double HIGHLIGHT_DISTANCE = SIZE*0.8; // used to be Math.sqrt(SIZE*SIZE/2);
	public static final int DOUBLE_CLICK_INTERVAL = 500;

	private Color[] colorScheme = {Color.BLACK,
			new Color(55, 129, 252), // dark/light blue
			new Color(140, 20, 180), // deep purple
			new Color(0, 100, 0), // green
			new Color(200, 30, 10), // red
			new Color(100, 50, 0), // brown
			new Color(255, 100, 0), // orange
			new Color(217, 109, 194), // pink
			new Color(94, 88, 214), // light purple
			new Color(189, 243, 61), // neon green
			new Color(33, 31, 53), // nearly black
			new Color(255, 255, 255), // white
	};

	private Project proj;
	private String[][] samples;
	private double[][] data;
	private double minX;
	private double maxX;
	private double minY;
	private double maxY;
	// private boolean linkSamples;
	private String prevPos = "";
	private BufferedImage image;
	private Hashtable<String,IntVector> sampLookup;
	private IntVector prox;
	private Hashtable<String,String> colorHash;
	private Hashtable<String,String> failedHash;

	public MosaicPanel(Project proj, String[][] samples, double[][] data) {
		BufferedReader reader;
		String[] line;
		int[] indices;

		this.proj = proj;
		this.samples = samples;
		this.data = data;

		colorHash = new Hashtable<String,String>();
		try {
			reader = Files.getReader(proj.getFilename(Project.MOSAIC_COLOR_CODES_FILENAME), proj.getJarStatus(), true, false);
			if (reader!=null) {
				while (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
					colorHash.put(line[0]+"\t"+line[1], line[2]);
				}
				reader.close();
			}
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+proj.getFilename(Project.MOSAIC_COLOR_CODES_FILENAME)+"\"");
			System.exit(2);
		}

		failedHash = new Hashtable<String,String>();
		try {
			reader = Files.getReader(proj.getFilename(Project.SAMPLE_DATA_FILENAME), proj.getJarStatus(), true, false);
			if (reader!=null) {
				indices = ext.indexFactors(new String[] {"Sample Name", "CLASS=Suitable for CNV"}, reader.readLine().trim().split("\t", -1), false, false);
				if (Array.min(indices)>=0) {
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						failedHash.put(line[indices[0]], line[indices[1]]);
					}
				}
				reader.close();
			}
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+proj.getFilename(Project.SAMPLE_DATA_FILENAME)+"\"");
		}

		image = null;
		locLookup = new Hashtable<String,IntVector>();
		sampLookup = new Hashtable<String,IntVector>();
		// linkSamples = true;

		for (int i = 0; i<data.length; i++) {
			if (sampLookup.containsKey(samples[i][0])) {
				sampLookup.get(samples[i][0]).add(i);
			} else {
				sampLookup.put(samples[i][0], new IntVector(new int[] {i}));
			}
		}

		for (int j = 0; j<data.length; j++) {
			minX = Math.min(data[j][0], minX);
			maxX = Math.max(data[j][0], maxX);
			minY = Math.min(data[j][1], minY);
			maxY = Math.max(data[j][1], maxY);
		}

		// taken care of in AbstractPanel constructor
//		addMouseListener(this);
//		addMouseMotionListener(this);
//		addComponentListener(this);
//		
		setFlow(true);
	}

	public void savePlotToFile() {
		try {
			ImageIO.write(image, "png", new File("tryit.png"));
		} catch (IOException ie) {
			JOptionPane.showMessageDialog(null, "Error while trying to save the plot");
		}
	}

	public void mouseMoved(MouseEvent event) {
		Graphics g = getGraphics();
		IntVector iv;
		String pos;
		int x, y;

		if (getFinalImage()) {
			x = event.getX();
			y = event.getY();

			pos = (int)Math.floor(x/LOOKUP_RESOLUTION)+"x"+(int)Math.floor(y/LOOKUP_RESOLUTION);
			if (!pos.equals(prevPos)) {
				repaint();
			}

			//iv = locLookup.get(pos);
			iv = lookupNearbyPoints(x, y, pos);
			prox = new IntVector();
			//System.out.println("pos: "+pos+"\t iv.size():"+(iv==null?"null":iv.size()));//zx test point
			g.setColor(Color.RED);
			for (int i = 0; iv!=null&&i<iv.size(); i++) {
				if (Distance.euclidean(new int[] {x, y}, new int[] {getXPixel(data[iv.elementAt(i)][0]), getYPixel(data[iv.elementAt(i)][1])})<HIGHLIGHT_DISTANCE) {
					g.setColor(Color.RED);
					prox.add(iv.elementAt(i));
					if (failedHash.containsKey(samples[iv.elementAt(i)][0])&&failedHash.get(samples[iv.elementAt(i)][0]).equals("0")) {
						g.fillOval(getXPixel(data[iv.elementAt(i)][0])-SIZE_FAILED/2, getYPixel(data[iv.elementAt(i)][1])-SIZE_FAILED/2, SIZE_FAILED, SIZE_FAILED);
					} else {
						g.fillOval(getXPixel(data[iv.elementAt(i)][0])-SIZE/2, getYPixel(data[iv.elementAt(i)][1])-SIZE/2, SIZE, SIZE);
					}

					// } else {
					// g.setColor(Color.BLACK);
					// g.fillOval(getX(data[iv.elementAt(i)][0])-SIZE/2,
					// getY(data[iv.elementAt(i)][1])-SIZE/2, SIZE, SIZE);
				}
			}
			// if (linkSamples && prox != null && prox.size() > 0) {
			if (prox!=null&&prox.size()>0) {
				for (int i = 0; i<prox.size(); i++) {
					iv = sampLookup.get(samples[prox.elementAt(i)][0]);
					g.setColor(Color.YELLOW);
					for (int j = 0; j<Math.min(iv.size(), 10); j++) {
						if (failedHash.containsKey(samples[iv.elementAt(j)][0])&&failedHash.get(samples[iv.elementAt(j)][0]).equals("0")) {
							g.fillOval(getXPixel(data[iv.elementAt(j)][0])-SIZE_FAILED/2, getYPixel(data[iv.elementAt(j)][1])-SIZE_FAILED/2, SIZE_FAILED, SIZE_FAILED);
						} else {
							g.fillOval(getXPixel(data[iv.elementAt(j)][0])-SIZE/2, getYPixel(data[iv.elementAt(j)][1])-SIZE/2, SIZE, SIZE);
						}
					}
				}
			}

			prevPos = pos;
		}
	}

	public void mouseClicked(MouseEvent event) {
		JPopupMenu menu;

		if (event.getButton()==MouseEvent.BUTTON1) { // left click
			// linkSamples = !linkSamples;
		} else if (event.getButton()==MouseEvent.BUTTON3) { // right click
		}

		//System.out.println("Click with "+prox.size()+" in proximity");
		if (prox!=null&&prox.size()>0) {
			menu = new JPopupMenu();
			for (int i = 0; i<prox.size(); i++) {
				menu.add(new LaunchAction(proj, samples[prox.elementAt(i)][0], samples[prox.elementAt(i)][1], colorHash.containsKey(samples[prox.elementAt(i)][0]+"\t"+samples[prox.elementAt(i)][1])?colorScheme[Integer.parseInt(colorHash.get(samples[prox.elementAt(i)][0]+"\t"+samples[prox.elementAt(i)][1]))]:Color.GRAY));
			}
			menu.show(this, event.getX(), event.getY());
		}
	}

	public void mouseEntered(MouseEvent e) {}

	public void mouseExited(MouseEvent e) {}

	public void mousePressed(MouseEvent e) {}

	public void mouseReleased(MouseEvent e) {}

	public void mouseDragged(MouseEvent e) {}

	public static void main(String[] args) {
		MosaicPlot.main(args);
	}

	/**
	 * Assigns x-Axis label and y-axis label to the panel.
	 */
	@Override
	public void assignAxisLabels() {
		xAxisLabel = "Standard deviation of BAF values";
		yAxisLabel = "Inter-quartile range of BAF values";
	}

	public void highlightPoints() {}

	/**
	 * Generate the array "protected PlotPoint[] points" from the following:
	 * 		samples   - the array of sample indices. One sample is represented as one data point in the graph
	 * 		type      - Missing (4), Not_A_Number (3), or Filled_Circle (1) 
	 *		data      - double[][], with data[][0] being the X, and data[][1] being the Y;
	 *		size      - the size of the points on the drawing;
	 *		classcode - All (0), GenotypeCode (1), Sex (2) or etc
	 *		layer:    - layer (0 or 1) of the drawing. For MosaicPanel, we only use 0.
	 */
	@Override
	public void generatePoints() {
		byte color;
		String[] files;
		files = new File(proj.getDir(Project.SAMPLE_DIRECTORY)).list(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(Sample.SAMPLE_DATA_FILE_EXTENSION);
			}
		});
		if (files==null) {
			files = new String[0];
		}
		
		points = new PlotPoint[data.length];
		for (int i = 0; i<data.length&&getFlow(); i++) {
//		for (int i = 0; i<data.length; i++) {
			if (colorHash.containsKey(samples[i][0]+"\t"+samples[i][1])) {
				//color = colorScheme[Integer.parseInt(colorHash.get(samples[i][0]+"\t"+samples[i][1]))];
				color = (byte) Integer.parseInt(colorHash.get(samples[i][0]+"\t"+samples[i][1]));
			} else {
				color = (byte) ((byte) ext.indexOfStr(samples[i][0]+Sample.SAMPLE_DATA_FILE_EXTENSION, files)>=0?0:1);	// What is the color code for Color.GRAY
			}
			points[i] = new PlotPoint("",
									  (byte) 1,
									  (float) data[i][0],
									  (float) data[i][1],
									  (byte) SIZE,
									  (byte) color,
									  (byte) 0
									  );
		}
		setSwapable(false);

		//SampleData sampleData = new SampleData(proj, true);
		/*
		String[] markerList = Array.toStringArray(markerNames);
		MarkerData[] markerData = MarkerSet.loadFromList(proj, markerList);
		for (int i=0; i<samples.length; i++) {
			points[i] = new PlotPoint(samples[i][0],
											 1,
											 markerData[1].getDatapoints(0)[0][i],
											 markerData[1].getDatapoints(0)[1][i],
											 SIZE,
											 0,
											 0
											);
		}
		*/
		
		//Color color;
	}
}
