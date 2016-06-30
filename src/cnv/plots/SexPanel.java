package cnv.plots;

import java.io.*;
import java.util.*;
import java.awt.*;

import javax.swing.*;

import java.awt.event.*;
import java.awt.image.BufferedImage;

import javax.imageio.*;

import common.*;
import cnv.filesys.Project;
import cnv.gui.LaunchAction;
import cnv.var.SampleData;
import mining.Distance;


//public class MosaicPanel extends JPanel implements MouseListener, MouseMotionListener, ComponentListener {
//public class SexPanel extends AbstractPanel implements MouseListener, MouseMotionListener, ComponentListener, ActionListener {
public class SexPanel extends AbstractPanel implements MouseListener, MouseMotionListener {
	public static final long serialVersionUID = 3L;
	public static final int HEAD_BUFFER = 25;
//	public static final int HEIGHT_X_AXIS = 55;
	public static final int HEIGHT_X_AXIS = 105;
//	public static final int WIDTH_Y_AXIS = 75;
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
	public static final double HIGHLIGHT_DISTANCE = Math.sqrt(SIZE*SIZE/2);
	public static final int DOUBLE_CLICK_INTERVAL = 500;
	public static final double NUM_SD_FOR_OUTLIERS = 3.0;

	public static Color[] COLOR_SCHEME = { new Color(0, 0, 0),		// black.		for 0 - missing
										  new Color(0, 0, 255),		// blue.		for 1 - normal male
										  new Color(255, 175, 175),	// pink.		for 2 - normal female
										  new Color(0, 255, 0),		// Green.	    for 3 - full Klinefelter's XXY
										  new Color(0, 150, 150),	// blue-green.	for 4 - UPD Klinefelter's XXY
										  new Color(51, 153, 255),	// light blue.	for 5 - mosaic Klinefelter's
										  new Color(255, 0, 0),		// red.			for 6 - Triple X
										  new Color(178,34,34),		// firebrick.	for 7 - mosaic Triple X
										  new Color(66, 28, 82),	// dark purple. for 8 - Turner's XO without mosaicism
										  new Color(189, 174, 198)	// light purple.for 9 - mosaic Turner's
	};

	public static String[] COLOR_SCHEME_MEANING = { "Missing",
												  "Normal Male",
												  "Normal Female",
												  "Full Klinefelter's",
												  "UPD Klinefelter's",
												  "Mosaic Klinefelter's",
												  "Triple X",
												  "Mosaic Triple X",
												  "Full Turner's",
												  "Mosaic Turner's"
	};

	private Project proj;
	private String[][] samples;
	private double[][] data;
	private byte[] sexes;
	private byte[] estimatedSexes;
	private double minX;
	private double maxX;
	private double minY;
	private double maxY;
	private String prevPos = "";
	private BufferedImage image;
	private Hashtable<String,IntVector> sampLookup;
	private IntVector prox;
	private Hashtable<String,String> colorHash;
	private SampleData sampleData;

	public SexPanel(Project proj, String[][] samples, double[][] data, byte[] sexes, byte[] estimatedSexes) {
		this(proj, samples, data);
		if (sexes!=null) {
			this.sexes = sexes;
		}
		if (estimatedSexes!=null) {
			this.estimatedSexes = estimatedSexes;
		}
	}

	public SexPanel(Project proj, String[][] samples, double[][] data) {
		BufferedReader reader;
		String[] line;

		this.proj = proj;
		this.samples = samples;
		this.data = data;

		colorHash = new Hashtable<String,String>();
		try {
//			reader = Files.getReader(proj.getFilename(proj.MOSAIC_COLOR_CODES_FILENAME), proj.getJarStatus(), true, false);
			reader = Files.getReader(proj.MOSAIC_COLOR_CODES_FILENAME.getValue(), proj.JAR_STATUS.getValue(), true, false);
			if (reader!=null) {
				while (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
					colorHash.put(line[0]+"\t"+line[1], line[2]);
				}
				reader.close();
			}
		} catch (IOException ioe) {
//			System.err.println("Error reading file \""+proj.getFilename(proj.MOSAIC_COLOR_CODES_FILENAME)+"\"");
			System.err.println("Error reading file \""+proj.MOSAIC_COLOR_CODES_FILENAME.getValue()+"\"");
			System.exit(2);
		}

		image = null;
		locLookup = new Hashtable<String,IntVector>();
		sampLookup = new Hashtable<String,IntVector>();
		
		sampleData = proj.getSampleData(0, false);
		
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

//		addMouseListener(this);
//		addMouseMotionListener(this);
		
		//setFlow(true);
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

		if (imageIsFinal()) {
			x = event.getX();
			y = event.getY();

			pos = (int)Math.floor(x/LOOKUP_RESOLUTION)+"x"+(int)Math.floor(y/LOOKUP_RESOLUTION);
			if (!pos.equals(prevPos)) {
				repaint();
			}

			//iv = locLookup.get(pos);
			iv = lookupNearbyPoints(x, y, pos);
			prox = new IntVector();
//			System.out.println("pos: "+pos+"\t iv.size():"+(iv==null?"null":iv.size()));//zx test point
			g.setColor(Color.RED);
			for (int i = 0; iv!=null&&i<iv.size(); i++) {
				if (Distance.euclidean(new int[] {x, y}, new int[] {getXPixel(data[iv.elementAt(i)][0]), getYPixel(data[iv.elementAt(i)][1])})<HIGHLIGHT_DISTANCE) {
					g.setColor(Color.RED);
					prox.add(iv.elementAt(i));
					if (sampleData.individualShouldBeExcluded(samples[iv.elementAt(i)][0])) {
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
			if (prox!=null && prox.size()>0) {
				for (int i = 0; i<prox.size(); i++) {
					iv = sampLookup.get(samples[prox.elementAt(i)][0]);
					g.setColor(Color.YELLOW);
					for (int j = 0; j<Math.min(iv.size(), 10); j++) {
						if (sampleData.individualShouldBeExcluded(samples[iv.elementAt(j)][0])) {
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
//				System.out.println(i+"\t"+samples[prox.elementAt(i)][0]);
				menu.add(new LaunchAction(proj,
										  samples[prox.elementAt(i)][0],
										  new String[] {"chr23", "chr24"},
										  colorHash.containsKey(samples[prox.elementAt(i)][0]+"\t"+samples[prox.elementAt(i)][1])?
												  COLOR_SCHEME[Integer.parseInt(colorHash.get(samples[prox.elementAt(i)][0]+"\t"+samples[prox.elementAt(i)][1]))]:
												  COLOR_SCHEME[points[prox.elementAt(i)].getColor()]));
			}
			menu.show(this, event.getX(), event.getY());
		}
	}

	public void mouseEntered(MouseEvent e) {}

	public void mouseExited(MouseEvent e) {}

	public void mouseDragged(MouseEvent e) {}

	public static void main(String[] args) {
		SexPlot.main(args);
	}

	/**
	 * Assigns x-Axis label and y-axis label to the panel.
	 */
	@Override
	public void assignAxisLabels() {
		xAxisLabel = "Median X LRR";
		yAxisLabel = "Median Y LRR";
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
		/*
		byte color;
		String[] files;
		files = new File(proj.getDir(Project.IND_DIRECTORY)).list(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith("");
			}
		});
		if (files==null) {
			files = new String[0];
		}
		*/
		
//		byte[] sexNew = estimatedSex();
		points = new PlotPoint[data.length];
//		setColorScheme(new Color[] {Color.GRAY, Color.BLUE, Color.PINK});
		setColorScheme(COLOR_SCHEME);
//		for (int i = 0; i<data.length&&getFlow(); i++) {
		for (int i = 0; i<data.length; i++) {
			/* TODO would be nice to move this functionality (grey out those samples with no sample data) to TwoDPlot
			if (colorHash.containsKey(samples[i][0]+"\t"+samples[i][1])) {
				//color = colorScheme[Integer.parseInt(colorHash.get(samples[i][0]+"\t"+samples[i][1]))];
				color = (byte) Integer.parseInt(colorHash.get(samples[i][0]+"\t"+samples[i][1]));
			} else {
				color = (byte) ((byte) ext.indexOfStr(samples[i][0]+".samp", files)>=0?0:1);	// What is the color code for Color.GRAY
			}
			*/
			if (!sampleData.individualShouldBeExcluded(samples[i][0])) {
				points[i] = new PlotPoint("",
										  (byte) 1,
										  (float) data[i][0],
										  (float) data[i][1],
										  (byte) SIZE,
										  estimatedSexes[i],	//sexes[i],
										  (byte) (estimatedSexes[i]==sexes[i]?0:2)
										  );
			}
		}
		//Color color;
		setSwapable(false);
	}
}
