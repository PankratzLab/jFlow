package org.genvisis.cnv.plots;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Toolkit;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.Hashtable;
import java.util.Set;

import javax.imageio.ImageIO;
import javax.swing.JOptionPane;
import javax.swing.JPopupMenu;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.gui.LaunchAction;
import org.genvisis.cnv.plots.PlotPoint.PointType;
import org.genvisis.cnv.qc.SexChecks;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.IntVector;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.PSF.Colors.BLUES;
import org.pankratzlab.common.PSF.Colors.GREENS;
import org.pankratzlab.common.PSF.Colors.REDS;
import org.pankratzlab.common.PSF.Colors.VIOLETS;
import org.pankratzlab.common.mining.Distance;

public class SexPanel extends AbstractPanel implements MouseListener, MouseMotionListener {

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
  public static final double HIGHLIGHT_DISTANCE = Math.sqrt(SIZE * SIZE / 2);
  public static final int DOUBLE_CLICK_INTERVAL = 500;
  public static final double NUM_SD_FOR_OUTLIERS = 3.0;

  public static Color[] COLOR_SCHEME = {Color.BLACK, // black. for 0 - missing
                                        Color.BLUE, // blue. for 1 - normal male
                                        Color.PINK, // pink. for 2 - normal female
                                        Color.GREEN, // Green. for 3 - full Klinefelter's
                                        // XXY
                                        GREENS.DARK_CYAN, // blue-green. for 4 - UPD
                                        // Klinefelter's XXY
                                        BLUES.DODGER_BLUE, // light blue. for 5 - mosaic
                                        // Klinefelter's
                                        Color.RED, // red. for 6 - Triple X
                                        REDS.FIREBRICK, // firebrick. for 7 - mosaic Triple
                                        // X
                                        VIOLETS.CHRISTALLE, // dark purple. for 8 - Turner's XO
                                        // without mosaicism
                                        VIOLETS.LOLA // light purple.for 9 - mosaic
                                                     // Turner's
  };

  public static String[] COLOR_SCHEME_MEANING = {"Missing", "Normal Male", "Normal Female",
                                                 "Full Klinefelter's", "UPD Klinefelter's",
                                                 "Mosaic Klinefelter's", "Triple X",
                                                 "Mosaic Triple X", "Full Turner's",
                                                 "Mosaic Turner's"};

  private final Project proj;
  private final String[] samples;
  private final double[][] data;
  private byte[] sexes;
  private byte[] estimatedSexes;
  private final boolean[] excluded;
  private final boolean[] uncertains;
  private final String[] notes;
  private boolean showExcluded = false;

  private double minX;
  private double maxX;
  private double minY;
  private double maxY;
  private String prevPos = "";
  private final BufferedImage image;
  private final Hashtable<String, IntVector> sampLookup;
  private IntVector prox;
  private final Hashtable<String, String> colorHash;

  public SexPanel(Project proj, String[] samples, double[][] data, byte[] sexes,
                  byte[] estimatedSexes, boolean[] excluded, boolean[] uncertains, String[] notes) {
    this(proj, samples, data, excluded, uncertains, notes);
    if (sexes != null) {
      this.sexes = sexes;
    }
    if (estimatedSexes != null) {
      this.estimatedSexes = estimatedSexes;
    }
  }

  public SexPanel(Project proj, String[] samples, double[][] data, boolean[] excluded,
                  boolean[] uncertains, String[] notes) {
    BufferedReader reader;
    String[] line;

    this.proj = proj;
    this.samples = samples;
    this.data = data;
    this.excluded = excluded;
    this.uncertains = uncertains;
    this.notes = notes;

    colorHash = new Hashtable<>();
    try {
      reader = Files.getReader(proj.MOSAIC_COLOR_CODES_FILENAME.getValue(), true, false);
      if (reader != null) {
        while (reader.ready()) {
          line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
          colorHash.put(line[0] + "\t" + line[1], line[2]);
        }
        reader.close();
      }
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + proj.MOSAIC_COLOR_CODES_FILENAME.getValue()
                         + "\"");
      System.exit(2);
    }

    image = null;
    locLookup = new Hashtable<>();
    sampLookup = new Hashtable<>();

    for (int i = 0; i < data.length; i++) {
      if (sampLookup.containsKey(samples[i])) {
        sampLookup.get(samples[i]).add(i);
      } else {
        sampLookup.put(samples[i], new IntVector(new int[] {i}));
      }
    }

    for (double[] element : data) {
      minX = Math.min(element[0], minX);
      maxX = Math.max(element[0], maxX);
      minY = Math.min(element[1], minY);
      maxY = Math.max(element[1], maxY);
    }
  }

  public void setShowExcluded(boolean showExcluded) {
    this.showExcluded = showExcluded;
  }

  public void uncertainsToTrailer() {
    int numUncertains = ArrayUtils.booleanArraySum(uncertains);
    if (numUncertains == 0) {
      JOptionPane.showMessageDialog(null, "No uncertain sex calls to display");
      return;
    }
    String[][] uncertainsXRegions = new String[numUncertains][];
    String[][] uncertainsYRegions = new String[numUncertains][];
    int writeIndex = 0;
    for (int i = 0; i < samples.length; i++) {
      if (uncertains[i]) {
        uncertainsXRegions[writeIndex] = new String[] {samples[i], "chr23",
                                                       "Called as " + SexChecks.EstimatedSex.values()[estimatedSexes[i]].toString()
                                                                            + " (" + notes[i]
                                                                            + ")"};
        uncertainsYRegions[writeIndex++] = new String[] {samples[i], "chr24",
                                                         "Called as " + SexChecks.EstimatedSex.values()[estimatedSexes[i]].toString()
                                                                              + " (" + notes[i]
                                                                              + ")"};
      }
    }
    String xRegionsFile = proj.RESULTS_DIRECTORY.getValue() + "list_uncertainSexCallsXChr.txt";
    String yRegionsFile = proj.RESULTS_DIRECTORY.getValue() + "list_uncertainSexCallsYChr.txt";
    Files.writeMatrix(uncertainsXRegions, xRegionsFile, "\t");
    Files.writeMatrix(uncertainsYRegions, yRegionsFile, "\t");

    Trailer xTrailer = new Trailer(proj, uncertainsXRegions[0][0], proj.CNV_FILENAMES.getValue(),
                                   uncertainsXRegions[0][1], Trailer.DEFAULT_STARTX, 1,
                                   Toolkit.getDefaultToolkit()
                                          .getScreenSize().width - 30 - Trailer.DEFAULT_STARTX,
                                   (Toolkit.getDefaultToolkit().getScreenSize().height - 50) / 2);
    Trailer yTrailer = new Trailer(proj, uncertainsYRegions[0][0], proj.CNV_FILENAMES.getValue(),
                                   uncertainsYRegions[0][1], Trailer.DEFAULT_STARTX,
                                   1 + (Toolkit.getDefaultToolkit().getScreenSize().height - 50)
                                       / 2,
                                   Toolkit.getDefaultToolkit()
                                          .getScreenSize().width - 30 - Trailer.DEFAULT_STARTX,
                                   (Toolkit.getDefaultToolkit().getScreenSize().height - 50) / 2);

    xTrailer.loadRegionFile(xRegionsFile);
    yTrailer.loadRegionFile(yRegionsFile);

    xTrailer.setVisible(true);
    yTrailer.setVisible(true);
  }

  public void savePlotToFile() {
    try {
      ImageIO.write(image, "png", new File("tryit.png"));
    } catch (IOException ie) {
      JOptionPane.showMessageDialog(null, "Error while trying to save the plot");
    }
  }

  @Override
  public void mouseMoved(MouseEvent event) {
    Graphics g = getGraphics();
    Set<Integer> nearby;
    IntVector iv;
    String pos;
    int x, y;

    if (imageIsFinal()) {
      x = event.getX();
      y = event.getY();

      pos = (int) Math.floor(x / LOOKUP_RESOLUTION) + "x" + (int) Math.floor(y / LOOKUP_RESOLUTION);
      if (!pos.equals(prevPos)) {
        repaint();
      }

      nearby = lookupNearbyPoints(x, y, pos);
      prox = new IntVector();
      g.setColor(Color.RED);
      for (Integer i : nearby) {
        int near = i.intValue();
        if (Distance.euclidean(new int[] {x, y},
                               new int[] {getXPixel(data[near][0]),
                                          getYPixel(data[near][1])}) < HIGHLIGHT_DISTANCE) {
          g.setColor(Color.RED);
          prox.add(near);
          if (excluded[near]) {
            g.fillOval(getXPixel(data[near][0]) - SIZE_FAILED / 2,
                       getYPixel(data[near][1]) - SIZE_FAILED / 2, SIZE_FAILED, SIZE_FAILED);
          } else {
            g.fillOval(getXPixel(data[near][0]) - SIZE / 2, getYPixel(data[near][1]) - SIZE / 2,
                       SIZE, SIZE);
          }
        }
      }
      if (prox != null && prox.size() > 0) {
        for (int i = 0; i < prox.size(); i++) {
          iv = sampLookup.get(samples[prox.elementAt(i)]);
          g.setColor(Color.YELLOW);
          for (int j = 0; j < Math.min(iv.size(), 10); j++) {
            if (excluded[iv.elementAt(j)]) {
              g.fillOval(getXPixel(data[iv.elementAt(j)][0]) - SIZE_FAILED / 2,
                         getYPixel(data[iv.elementAt(j)][1]) - SIZE_FAILED / 2, SIZE_FAILED,
                         SIZE_FAILED);
            } else {
              g.fillOval(getXPixel(data[iv.elementAt(j)][0]) - SIZE / 2,
                         getYPixel(data[iv.elementAt(j)][1]) - SIZE / 2, SIZE, SIZE);
            }
          }
        }
      }

      prevPos = pos;
    }
  }

  @Override
  public void mouseClicked(MouseEvent event) {
    JPopupMenu menu;

    if (prox != null && prox.size() > 0) {
      menu = new JPopupMenu();
      for (int i = 0; i < prox.size(); i++) {
        menu.add(new LaunchAction(proj, samples[prox.elementAt(i)], new String[] {"chr23", "chr24"},
                                  colorHash.containsKey(samples[prox.elementAt(i)] + "\t"
                                                        + "chr23") ? COLOR_SCHEME[Integer.parseInt(colorHash.get(samples[prox.elementAt(i)] + "\t" + "chr23"))] : COLOR_SCHEME[points[prox.elementAt(i)].getColor()]));
      }
      menu.show(this, event.getX(), event.getY());
    }
  }

  @Override
  public void mouseEntered(MouseEvent e) {}

  @Override
  public void mouseExited(MouseEvent e) {}

  @Override
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

  @Override
  public void highlightPoints() {}

  /**
   * Generate the array "protected PlotPoint[] points" from the following: samples - the array of
   * sample indices. One sample is represented as one data point in the graph type - Missing (4),
   * Not_A_Number (3), or Filled_Circle (1) data - double[][], with data[][0] being the X, and
   * data[][1] being the Y; size - the size of the points on the drawing; classcode - All (0),
   * GenotypeCode (1), Sex (2) or etc layer: - layer (0 or 1) of the drawing. For MosaicPanel, we
   * only use 0.
   */
  @Override
  public void generatePoints() {
    /*
     * byte color; String[] files; files = new File(proj.getDir(Project.IND_DIRECTORY)).list(new
     * FilenameFilter() { public boolean accept(File file, String filename) { return
     * filename.endsWith(""); } }); if (files==null) { files = new String[0]; }
     */

    // byte[] sexNew = estimatedSex();
    points = new PlotPoint[data.length];
    // setColorScheme(new Color[] {Color.GRAY, Color.BLUE, Color.PINK});
    setColorScheme(COLOR_SCHEME);
    // for (int i = 0; i<data.length&&getFlow(); i++) {
    for (int i = 0; i < data.length; i++) {
      /*
       * TODO would be nice to move this functionality (grey out those samples with no sample data)
       * to TwoDPlot if (colorHash.containsKey(samples[i][0]+"\t"+samples[i][1])) { //color =
       * colorScheme[Integer.parseInt(colorHash.get(samples[i][0]+"\t"+samples[i][1]))]; color =
       * (byte) Integer.parseInt(colorHash.get(samples[i][0]+"\t"+samples[i][1])); } else { color =
       * (byte) ((byte) ext.indexOfStr(samples[i][0]+".samp", files)>=0?0:1); // What is the color
       * code for Color.GRAY }
       */
      if (showExcluded || !excluded[i]) {
        points[i] = new PlotPoint("", PointType.FILLED_CIRCLE, (float) data[i][0],
                                  (float) data[i][1], (byte) SIZE, estimatedSexes[i], // sexes[i],
                                  (byte) (estimatedSexes[i] == sexes[i] ? 0 : 2));
      }
    }
    // Color color;
    setSwapable(false);
  }
}
