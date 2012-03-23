package cnv.plots;

import java.io.*;
import java.util.*;
import common.*;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;
import java.awt.image.BufferedImage;
import javax.imageio.*;

import cnv.filesys.MarkerData;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.gui.LaunchAction;
import cnv.var.SampleData;
import mining.Distance;

//public class MosaicPanel extends JPanel implements MouseListener, MouseMotionListener, ComponentListener {
public class GenderPanel extends AbstractPanel implements MouseListener, MouseMotionListener, ComponentListener {
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

	private Color[] colorScheme = {Color.BLACK,
			new Color(55, 129, 252), // dark/light blue
			new Color(140, 20, 180), // deep purple
			new Color(189, 243, 61), // light green
			new Color(201, 30, 10), new Color(217, 109, 194), // deep red/pink
			new Color(33, 31, 53), new Color(255, 255, 255), // dark dark / light light
			new Color(94, 88, 214), // light purple
			new Color(33, 87, 0)}; // dark green

	private Project proj;
	private String[][] samples;
	private double[][] data;
	private double minX;
	private double maxX;
	private double minY;
	private double maxY;
	// private boolean linkSamples;
	private int xMin = 0;
	private int xMax = getWidth();
	private int yMin = 0;
	private int yMax = getWidth();
	private double plotXmax, plotYmax;
	private String prevPos = "";
	private BufferedImage image;
	//private Hashtable<String,IntVector> locLookup;//zx
	private Hashtable<String,IntVector> sampLookup;
	private IntVector prox;
	private Hashtable<String,String> colorHash;
	private Hashtable<String,String> failedHash;
	private int prevWidth;
	private int prevHeight;
	//private boolean flow;
	//private boolean finalImage;//zx
	private Repress patience;

	public GenderPanel(Project proj, String[][] samples, double[][] data) {
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

		prevWidth = prevHeight = -1;

		addMouseListener(this);
		addMouseMotionListener(this);
		addComponentListener(this);
		
		//setFlow(true);
	}

//	public void paintComponent(Graphics g) {
//		g.setColor(Color.WHITE);
//		g.fillRect(0, 0, getWidth(), getHeight());
//		if (finalImage&&image!=null) {
//			g.drawImage(image, 0, 0, this);
//		}
//	}

//	public void interruptFlow() {
//		flow = false;
//	}

//	public void createImage() {
//		repaint();
//		flow = true;
//		image = new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_INT_RGB);
//		drawAll(image.createGraphics());
//		repaint();
//	}

	public void savePlotToFile() {
		try {
			ImageIO.write(image, "png", new File("tryit.png"));
		} catch (IOException ie) {
			JOptionPane.showMessageDialog(null, "Error while trying to save the plot");
		}
	}

	/*
	public void drawAll(Graphics g) {
		double xStep = X_STEP;
		double yStep = Y_STEP;
		String str, pos;
		int x, y;
		long time;
		BufferedImage yLabel;
		FontMetrics fontMetrics;
		Graphics gfx;
		ProgressBarDialog prog;
		int step;
		String[] files;

		files = new File(proj.getDir(Project.IND_DIRECTORY)).list(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith("");
			}
		});
		if (files==null) {
			files = new String[0];
		}

		time = new Date().getTime();
		prog = new ProgressBarDialog("Generating image...", 0, data.length, getWidth(), getHeight());
		finalImage = false;

		g.setFont(new Font("Arial", 0, 28));
		fontMetrics = g.getFontMetrics(g.getFont());

		g.setColor(Color.WHITE);
		g.fillRect(0, 0, getWidth(), getHeight());

		while (maxX/xStep>10) {
			xStep += X_STEP;
		}
		plotXmax = Math.ceil(maxX/xStep)*xStep;

		while (maxY/yStep>10) {
			yStep += Y_STEP;
		}
		plotYmax = Math.ceil(maxY/yStep)*yStep;

		// plot
		xMin = WIDTH_Y_AXIS;
		xMax = getWidth();
		yMin = HEIGHT_X_AXIS;
		yMax = getHeight()-HEAD_BUFFER;
		locLookup.clear();
		step = data.length/100;
		if (step==0) {
			step = 1;
		}
		for (int i = 0; i<data.length&&flow; i++) {
			if (i%step==0) {
				prog.setProgress(i);
			}
			if (colorHash.containsKey(samples[i][0]+"\t"+samples[i][1])) {
				g.setColor(colorScheme[Integer.parseInt(colorHash.get(samples[i][0]+"\t"+samples[i][1]))]);
			} else {
				g.setColor(ext.indexOfStr(samples[i][0]+".samp", files)>=0?Color.BLACK:Color.GRAY);
			}
			if (failedHash.containsKey(samples[i][0])&&failedHash.get(samples[i][0]).equals("0")) {
				g.fillOval(getX(data[i][0])-SIZE_FAILED/2, getY(data[i][1])-SIZE_FAILED/2, SIZE_FAILED, SIZE_FAILED);
			} else {
				g.fillOval(getX(data[i][0])-SIZE/2, getY(data[i][1])-SIZE/2, SIZE, SIZE);
			}
			x = (int)Math.floor(getX(data[i][0])/LOOKUP_RESOLUTION);
			y = (int)Math.floor(getY(data[i][1])/LOOKUP_RESOLUTION);
			for (int j = x-1; j<=x+1; j++) {
				for (int k = y-1; k<=y+1; k++) {
					pos = j+"x"+k;
					if (locLookup.containsKey(pos)) {
						locLookup.get(pos).add(i);
					} else {
						locLookup.put(pos, new IntVector(new int[] {i}));
					}
				}
			}
		}

		// x-Axis
		xMin = WIDTH_Y_AXIS;
		xMax = getWidth()-WIDTH_BUFFER;
		yMin = 0;
		yMax = HEIGHT_X_AXIS;
		// System.out.println(xMin+"-"+xMax+"\t"+yMin+"-"+yMax);

		for (int i = 0; i<=plotXmax/xStep; i++) {
			Grafik.drawThickLine(g, getX(i*xStep), getHeight()-yMax, getX(i*xStep), getHeight()-(yMax-TICK_LENGTH), TICK_THICKNESS, Color.BLACK);
			str = ext.formDeci(i*xStep, 2);
			g.drawString(str, getX(i*xStep)-str.length()*8, getHeight()-(yMax-TICK_LENGTH-30));
		}
		System.out.println(xMin+"-"+xMax+"\t"+yMin+"-"+yMax);
		Grafik.drawThickLine(g, xMin-(int)Math.ceil((double)AXIS_THICKNESS/2.0), getHeight()-yMax, xMax+(int)Math.ceil((double)AXIS_THICKNESS/2.0), getHeight()-yMax, AXIS_THICKNESS, Color.BLACK);
		str = "Standard deviation of BAF values";
		g.drawString(str, (getWidth()-WIDTH_Y_AXIS)/2-fontMetrics.stringWidth(str)/2+WIDTH_Y_AXIS, getHeight()-25);

		// y-axis
		xMin = 0;
		xMax = WIDTH_Y_AXIS;
		yMin = HEIGHT_X_AXIS;
		yMax = getHeight()-HEAD_BUFFER;
		for (int i = 0; i<=plotYmax/yStep; i++) {
			Grafik.drawThickLine(g, xMax-TICK_LENGTH, getY(i*yStep), xMax, getY(i*yStep), TICK_THICKNESS, Color.BLACK);
			str = ext.formDeci(i*yStep, 2);
			g.drawString(str, xMax-TICK_LENGTH-str.length()*15-5, getY(i*yStep)+9);
		}
		Grafik.drawThickLine(g, xMax, getY(minY), xMax, getY(plotYmax)-(int)Math.ceil((double)TICK_THICKNESS/2.0), AXIS_THICKNESS, Color.BLACK);

		str = "Inter-quartile range of BAF values";
		yLabel = new BufferedImage(fontMetrics.stringWidth(str), 36, BufferedImage.TYPE_INT_RGB);
		gfx = yLabel.createGraphics();
		gfx.setFont(new Font("Arial", 0, 28));
		gfx.setColor(Color.WHITE);
		gfx.fillRect(0, 0, getWidth(), getHeight());
		gfx.setColor(Color.BLACK);
		gfx.drawString(str, 0, yLabel.getHeight()-6);

		g.drawImage(Grafik.rotateImage(yLabel, true), 10, (getHeight()-HEIGHT_X_AXIS)/2-fontMetrics.stringWidth(str)/2, this);

		// try {
		// ImageIO.write(yLabel, "png", new File("yLabel.png"));
		// } catch(IOException ie) {
		// JOptionPane.showMessageDialog(null, "Error while trying to save the
		// plot");
		// }
		//
		prog.close();
		finalImage = true;
		System.out.println("Took "+ext.getTimeElapsed(time)+" to paint");
	}
	*/

//	public int getX(double x) {
//		return (int)((x-minX)/(plotXmax-minX)*(double)(xMax-xMin))+xMin;
//	}
//
//	public int getY(double y) {
//		return getHeight()-(int)((y-minY)/(plotYmax-minY)*(double)(yMax-yMin)+yMin);
//	}

	public void mouseMoved(MouseEvent event) {
		Graphics g = getGraphics();
		IntVector iv;
		String pos;
		int x, y;

		if (getFinalImage()) {
			x = event.getX();
			y = event.getY();

			xMin = WIDTH_Y_AXIS;
			xMax = getWidth();
			yMin = HEIGHT_X_AXIS;
			yMax = getHeight()-HEAD_BUFFER;
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
				if (Distance.euclidean(new int[] {x, y}, new int[] {getX(data[iv.elementAt(i)][0]), getY(data[iv.elementAt(i)][1])})<HIGHLIGHT_DISTANCE) {
					g.setColor(Color.RED);
					prox.add(iv.elementAt(i));
					if (failedHash.containsKey(samples[iv.elementAt(i)][0])&&failedHash.get(samples[iv.elementAt(i)][0]).equals("0")) {
						g.fillOval(getX(data[iv.elementAt(i)][0])-SIZE_FAILED/2, getY(data[iv.elementAt(i)][1])-SIZE_FAILED/2, SIZE_FAILED, SIZE_FAILED);
					} else {
						g.fillOval(getX(data[iv.elementAt(i)][0])-SIZE/2, getY(data[iv.elementAt(i)][1])-SIZE/2, SIZE, SIZE);
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
							g.fillOval(getX(data[iv.elementAt(j)][0])-SIZE_FAILED/2, getY(data[iv.elementAt(j)][1])-SIZE_FAILED/2, SIZE_FAILED, SIZE_FAILED);
						} else {
							g.fillOval(getX(data[iv.elementAt(j)][0])-SIZE/2, getY(data[iv.elementAt(j)][1])-SIZE/2, SIZE, SIZE);
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

	public void componentHidden(ComponentEvent e) {}

	public void componentMoved(ComponentEvent e) {}

	public void componentResized(ComponentEvent e) {
		setFlow(false);//zx
		if (prevWidth==-1) {
			prevWidth = getWidth();
			prevHeight = getHeight();
		} else if (prevWidth!=getWidth()||prevHeight!=getHeight()) {
			image = null;
			repaint();

			if (patience!=null) {
				patience.cancel();
				patience = null;
			}
//			new Thread(patience = new Repress(this, 100)).start();

			prevWidth = getWidth();
			prevHeight = getHeight();
		}
		setFlow(true);//zx
	}

	public void componentShown(ComponentEvent e) {}

	public static void main(String[] args) {
		MosaicPlot.main(args);
	}

	/**
	 * Assigns x-Axis label and y-axis label to the panel.
	 */
	@Override
	void assignAxisLabels() {
		xAxisLabel = "LRR";
		yAxisLabel = "LRR";
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
	void generatePoints() {
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
		
		points = new PlotPoint[data.length];
		for (int i = 0; i<data.length&&getFlow(); i++) {
			if (colorHash.containsKey(samples[i][0]+"\t"+samples[i][1])) {
				//color = colorScheme[Integer.parseInt(colorHash.get(samples[i][0]+"\t"+samples[i][1]))];
				color = (byte) Integer.parseInt(colorHash.get(samples[i][0]+"\t"+samples[i][1]));
			} else {
				color = (byte) ((byte) ext.indexOfStr(samples[i][0]+".samp", files)>=0?0:1);	// What is the color code for Color.GRAY
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
