package org.genvisis.cnv.plots;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.Hashtable;

import javax.swing.JPopupMenu;

import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.CountVector;
import org.genvisis.common.Grafik;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.Matrix;
import org.genvisis.mining.Distance;

public class StratPanel extends AbstractPanel
                        implements MouseListener, MouseMotionListener, ComponentListener {
  public static final long serialVersionUID = 3L;
  public static final int HEAD_BUFFER = 25;
  // public static final int HEIGHT_X_AXIS = 55;
  public static final int HEIGHT_X_AXIS = 105;
  // public static final int WIDTH_Y_AXIS = 75;
  public static final int WIDTH_Y_AXIS = 125;
  public static final int WIDTH_BUFFER = 50;
  public static final byte SIZE = 12;
  public static final byte SIZE_TAGALONGS = 6;
  public static final int AXIS_THICKNESS = 4;
  public static final int TICK_THICKNESS = 3;
  public static final int TICK_LENGTH = 15;
  public static final int MARKER_SIZE = 6;
  public static final double X_STEP = 0.05;
  public static final double Y_STEP = 0.05;
  public static final int LOOKUP_RESOLUTION = 20;
  public static final double HIGHLIGHT_DISTANCE = SIZE * 2;// used to be Math.sqrt(SIZE*SIZE/2);
  public static final int DOUBLE_CLICK_INTERVAL = 500;

  // public static final Color[] DEFAULT_COLORS = {Color.BLACK,
  // new Color(240, 0, 0), // red
  // new Color(155, 128, 0), // orange
  // new Color(247, 235, 0), // yellow
  // new Color(0, 128, 0), // green
  // new Color(0, 0, 128), // blue
  // new Color(139, 0, 103), // purple
  // new Color(25, 64, 0), // orange red
  // new Color(255, 191, 0), // yellow orange
  // new Color(139, 188, 0), // yellow green
  // new Color(0, 74, 101), // blue green
  // new Color(71, 0, 113), // blue purple
  // new Color(174, 0, 68)}; // red purple

  // public static final Color[] DEFAULT_COLORS = {Color.BLACK,
  // new Color(203, 30, 44), // red
  // new Color(0, 89, 168), // blue
  // new Color(52, 168, 63), // green
  // new Color(94, 42, 132), // purple
  // new Color(237, 112, 35), // orange
  // new Color(255, 241, 0), // yellow
  // new Color(230, 73, 39), // orange red
  // new Color(249, 171, 22), // yellow orange
  // new Color(182, 209, 36), // yellow green
  // new Color(0, 157, 126), // blue green
  // new Color(62, 46, 133), // blue purple
  // new Color(126, 37, 131)}; // red purple

  public static final Color[] DEFAULT_COLORS = {Color.BLACK, new Color(201, 30, 10), // red
                                                new Color(55, 129, 252), // light blue
                                                new Color(140, 20, 180), // deep purple
                                                new Color(33, 87, 0), // dark green
                                                new Color(247, 150, 70), // orange
                                                new Color(94, 88, 214), // light purple
                                                new Color(217, 109, 194), // deep red/pink
                                                new Color(189, 243, 61), // light green
                                                new Color(230, 73, 39), // orange red
                                                new Color(255, 241, 0), // yellow
                                                new Color(0, 157, 126), // blue green
                                                new Color(62, 46, 133), // blue purple
                                                new Color(0, 0, 128), // blue
                                                new Color(102, 51, 0), // brown
                                                new Color(153, 102, 51), // light brown
  };

  public static final Color[] BLUES =
      {new Color(25, 25, 112), new Color(0, 0, 128), new Color(100, 149, 237),
       new Color(72, 61, 139), new Color(106, 90, 205), new Color(123, 104, 238),
       new Color(132, 112, 255), new Color(0, 0, 205), new Color(65, 105, 225),
       new Color(0, 0, 255), new Color(30, 144, 255), new Color(0, 191, 255),
       new Color(135, 206, 250), new Color(135, 206, 250), new Color(70, 130, 180),
       new Color(176, 196, 222), new Color(173, 216, 230), new Color(176, 224, 230),
       new Color(175, 238, 238), new Color(0, 206, 209), new Color(72, 209, 204),
       new Color(64, 224, 208), new Color(0, 255, 255), new Color(224, 255, 255),
       new Color(95, 158, 160)};

  private final String[][] names;
  private final Hashtable<String, float[][]> hash;
  private final String[] sampleList;
  private Hashtable<String, IntVector> locLookup;
  private IntVector prox;
  private final boolean finalImage;
  private final StratPlot sp;
  private SampleData sampleData;
  private final CountVector classCounts;
  private int[][] prevPair;

  public StratPanel(StratPlot sp, String[][] names, Hashtable<String, float[][]> hash) {
    super();

    this.sp = sp;
    this.names = names;
    this.hash = hash;
    finalImage = false;
    classCounts = new CountVector();

    setColorScheme(DEFAULT_COLORS);

    sampleData = sp.getSampleData();
    sampleList = HashVec.getKeys(hash);

    createLookup = false;

    setNullMessage("Select two factors to plot");

    // taken care of in AbstractPanel constructor
    // addMouseListener(this);
    // addMouseMotionListener(this);
    // addComponentListener(this);
    setZoomable(true, true);
  }

  public void pushSampleData() {
    sampleData = sp.getSampleData();
  }

  @Override
  public void assignAxisLabels() {
    int[][] currentPair;

    currentPair = sp.getCurrentPair();
    if (names.length == 0) {
      sp.setDescriptor("Error - no files with .mds extension were present in the project directory");
      displayXaxis = displayYaxis = false;
    } else if (currentPair[0][0] == -1 || currentPair[1][0] == -1) {
      sp.setDescriptor("Double click on one of the files within the box on the left-hand side and select two of its components");
      if (currentPair[0][0] == -1) {
        displayXaxis = false;
      } else {
        xAxisLabel =
            names[currentPair[0][0]][0] + "_" + names[currentPair[0][0]][currentPair[0][1] + 1];
        displayXaxis = true;
      }
      if (currentPair[1][0] == -1) {
        displayYaxis = false;
      } else {
        yAxisLabel =
            names[currentPair[1][0]][0] + "_" + names[currentPair[1][0]][currentPair[1][1] + 1];
        displayYaxis = true;
      }
    } else {
      xAxisLabel =
          names[currentPair[0][0]][0] + "_" + names[currentPair[0][0]][currentPair[0][1] + 1];
      yAxisLabel =
          names[currentPair[1][0]][0] + "_" + names[currentPair[1][0]][currentPair[1][1] + 1];
      sp.setDescriptor(xAxisLabel + " x " + yAxisLabel);
      displayXaxis = displayYaxis = true;
    }
  }

  @Override
  public void highlightPoints() {}

  @Override
  public void generatePoints() {
    int currentClass;
    byte colorCode;
    float[][] data;
    int[][] currentPair;
    String[] sampleID;
    boolean tagalong;
    int countFails;

    currentClass = sp.getCurrentVariable();
    currentPair = sp.getCurrentPair();
    // System.out.println(Array.toStr(currentPair[1]) +"\t"+ Array.toStr(currentPair[0])); // TODO
    // this method is being called twice, once before and after changes
    if (prevPair == null || !Matrix.equals(currentPair, prevPair)) {
      resetZoomProportions();
    }
    prevPair = Matrix.clone(currentPair);


    if (currentPair[0][0] == -1 || currentPair[1][0] == -1) {
      points = new PlotPoint[0];
    } else {
      if (sampleData.getActualClassName(currentClass).equals("Site")) {
        setColorScheme(BLUES);
      } else {
        setColorScheme(DEFAULT_COLORS);
      }

      countFails = 0;
      points = new PlotPoint[sampleList.length];
      classCounts.clear();
      for (int i = 0; i < sampleList.length; i++) {
        data = hash.get(sampleList[i]);
        if (data[currentPair[0][0]] != null && data[currentPair[1][0]] != null
            && !Float.isNaN(data[currentPair[0][0]][currentPair[0][1]])
            && !Float.isNaN(data[currentPair[1][0]][currentPair[1][1]])) {
          sampleID = sampleData.lookup(sampleList[i]);
          if (sampleID == null) {
            if (countFails < 10) {
              sp.getProject().getLog().reportError("Error - could not look up " + sampleList[i]); // looks
                                                                                                  // up
                                                                                                  // any
                                                                                                  // individual
                                                                                                  // present
                                                                                                  // in
                                                                                                  // any
                                                                                                  // .mds
                                                                                                  // file
                                                                                                  // that
                                                                                                  // was
                                                                                                  // loaded,
                                                                                                  // even
                                                                                                  // those
                                                                                                  // not
                                                                                                  // in
                                                                                                  // the
                                                                                                  // current
                                                                                                  // file
            } else if (countFails == 10) {
              sp.getProject().getLog().reportError("..."); // looks up any individual present in any
                                                           // .mds file that was loaded, even those
                                                           // not in the current file
              countFails++;
            }
            tagalong = true;
            colorCode = 0;
            countFails++;
          } else {
            tagalong = false;
            colorCode = sampleData.getClassForInd(sampleID[0], currentClass);
            if (colorCode == -1 && !sp.maskMissing()) {
              colorCode = 0;
            }
          }

          points[i] = new PlotPoint(sampleList[i], PlotPoint.FILLED_CIRCLE,
                                    data[currentPair[0][0]][currentPair[0][1]],
                                    data[currentPair[1][0]][currentPair[1][1]],
                                    tagalong ? SIZE_TAGALONGS : SIZE, colorCode,
                                    (byte) (colorCode > 0 ? 1 : 0));
          classCounts.add(colorCode + "");
        }
      }

      if (countFails > 10) {
        sp.getProject().getLog()
          .reportError("There were a total of " + (countFails - 1)
                       + " sample ID pairs found in any of the files that were not found in SampleData");
      }
    }
    sp.updateColorKey(getClassCounts().convertToHash());
  }

  @Override
  public void mouseMoved(MouseEvent event) {
    Graphics g = getGraphics();
    int[][] currentPair;
    float[][] data;
    IntVector iv;
    String pos;
    int x, y;

    if (finalImage) {
      x = event.getX();
      y = event.getY();

      pos = (int) Math.floor(x / LOOKUP_RESOLUTION) + "x" + (int) Math.floor(y / LOOKUP_RESOLUTION);
      // if (!pos.equals(prevPos)) {
      // repaint();
      // }
      iv = locLookup.get(pos);
      prox = new IntVector();
      g.setColor(Color.RED);
      currentPair = sp.getCurrentPair();
      for (int i = 0; iv != null && i < iv.size(); i++) {
        data = hash.get(sampleList[iv.elementAt(i)]);

        if (Distance.euclidean(new int[] {x, y},
                               new int[] {getXPixel(data[currentPair[0][0]][currentPair[0][1]]),
                                          getYPixel(data[currentPair[1][0]][currentPair[1][1]])}) < HIGHLIGHT_DISTANCE) {
          g.setColor(Color.RED);
          prox.add(iv.elementAt(i));
          Grafik.drawCircle(g, getXPixel(data[currentPair[0][0]][currentPair[0][1]]),
                            getYPixel(data[currentPair[1][0]][currentPair[1][1]]), SIZE, true);
        }
      }

      // prevPos = pos;
    }
  }

  @Override
  public void mouseClicked(MouseEvent event) {
    JPopupMenu menu;

    if (event.getButton() == MouseEvent.BUTTON1) { // left click
      // linkSamples = !linkSamples;
    } else if (event.getButton() == MouseEvent.BUTTON3) { // right click
    }

    if (prox != null && prox.size() > 0) {
      menu = new JPopupMenu();
      for (int i = 0; i < prox.size(); i++) {
        // menu.add(new LaunchAction(proj.getProjectDir()+StratPlot.DEFAULT_FILENAME, proj,
        // samples[prox.elementAt(i)],
        // colorHash.containsKey(samples[prox.elementAt(i)][0]+"\t"+samples[prox.elementAt(i)][1])?colorScheme[Integer.parseInt(colorHash.get(samples[prox.elementAt(i)][0]+"\t"+samples[prox.elementAt(i)][1]))]:Color.GRAY));
      }
      menu.show(this, event.getX(), event.getY());
    }
  }

  public CountVector getClassCounts() {
    return classCounts;
  }

  @Override
  public void mouseEntered(MouseEvent e) {}

  @Override
  public void mouseExited(MouseEvent e) {}

  // public void mousePressed(MouseEvent e) {}
  //
  // public void mouseReleased(MouseEvent e) {}
  //
  // public void mouseDragged(MouseEvent e) {}

  @Override
  public void componentHidden(ComponentEvent e) {}

  @Override
  public void componentMoved(ComponentEvent e) {}

  @Override
  public void componentResized(ComponentEvent e) {
    paintAgain();
  }

  @Override
  public void componentShown(ComponentEvent e) {}

  public static void main(String[] args) {
    StratPlot.main(new String[] {"-notJar"});
  }
}
