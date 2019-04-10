package org.genvisis.cnv.plots;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Map;
import javax.swing.JOptionPane;
import javax.swing.JPopupMenu;
import javax.swing.SwingUtilities;
import org.genvisis.cnv.filesys.ClusterFilter;
import org.genvisis.cnv.filesys.ClusterFilterCollection;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.cnv.gui.LaunchAction;
import org.genvisis.cnv.manage.SexOps;
import org.genvisis.cnv.manage.SexOps.SEX_LOAD_TYPE;
import org.genvisis.cnv.plots.PlotPoint.PointType;
import org.genvisis.cnv.plots.ScatterPlot.PLOT_TYPE;
import org.genvisis.cnv.qc.GcAdjustorParameter.GcAdjustorParameters;
import org.genvisis.cnv.qc.MendelErrors.MendelErrorCheck;
import org.genvisis.cnv.var.IndiPheno;
import org.genvisis.cnv.var.SampleData;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CountVector;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.IntVector;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF.Colors.BLUES;
import org.pankratzlab.common.PSF.Colors.GREENS;
import org.pankratzlab.common.PSF.Colors.ORANGES;
import org.pankratzlab.common.PSF.Colors.REDS;
import org.pankratzlab.common.PSF.Colors.VIOLETS;
import org.pankratzlab.common.PSF.Colors.YELLOWS;

// TODO Needs some cleanup, especially MouseMoved, MouseClicked, and generatePoints
public class ScatterPanel extends AbstractPanel implements MouseListener, MouseMotionListener {

  public static final long serialVersionUID = 3L;
  public static Color[] DEFAULT_COLORS = new Color[] {BLUES.MIDNIGHT_EXPRESS, // dark dark
                                                      BLUES.PERSIAN_BLUE, // dark blue
                                                      REDS.VENETIAN_RED, // deep red
                                                      VIOLETS.BLUE_VIOLET, // deep purple
                                                      GREENS.GREEN, // dark green
                                                      BLUES.DODGER_BLUE, // light blue
                                                      BLUES.SLATE_BLUE, // light purple
                                                      GREENS.GREEN_YELLOW, // light green
                                                      VIOLETS.ORCHID, // pink
                                                      BLUES.NAVY, BLUES.CORNFLOWER_BLUE,
                                                      BLUES.DARK_SLATE_BLUE, BLUES.SLATE_BLUE,
                                                      BLUES.MEDIUM_SLATE_BLUE,
                                                      BLUES.LIGHT_SLATE_BLUE, BLUES.MEDIUM_BLUE,
                                                      BLUES.ROYAL_BLUE, Color.BLUE,
                                                      BLUES.DODGER_BLUE, BLUES.DEEP_SKY_BLUE,
                                                      BLUES.LIGHT_SKY_BLUE, BLUES.LIGHT_SKY_BLUE,
                                                      BLUES.STEEL_BLUE, BLUES.LIGHT_STEEL_BLUE,
                                                      BLUES.LIGHT_BLUE, BLUES.POWDER_BLUE,
                                                      BLUES.PALE_TURQUOISE, BLUES.DARK_TURQUOISE,
                                                      BLUES.MEDIUM_TURQUOISE, BLUES.TURQUOISE,
                                                      BLUES.AQUA, BLUES.LIGHT_CYAN, YELLOWS.AMBER, // yellowy
                                                                                                   // orange
                                                      ORANGES.MANGO_TANGO, // halloween
                                                                           // orange
  };

  byte[] alleleCounts;
  protected ScatterPlot sp;
  protected String[] samples;
  protected SampleData sampleData;
  protected IntVector indicesOfNearbyPoints;
  private boolean updateQcPanel; // A control variable. Do not update QcPanel when resizing, or etc.
  private int mouseStartX;
  private int mouseStartY;
  private int mouseEndX;
  private int mouseEndY;
  private final int panelIndex;
  CountVector uniqueValueCounts;
  protected boolean shrunk = false;
  HashMap<String, Integer> indPtMap;
  HashMap<String, Integer> sampIndMap;

  public ScatterPanel(ScatterPlot sp, int index) {
    super();

    this.sp = sp;
    panelIndex = index;
    samples = sp.getSamples();
    sampIndMap = new HashMap<>();
    for (int i = 0; i < samples.length; i++) {
      sampIndMap.put(samples[i], i);
    }
    sampleData = sp.getSampleData();
    updateQcPanel = true;

    setColorScheme(DEFAULT_COLORS);

    // taken care of in AbstractPanel constructor
    // addComponentListener(this);
    setZoomable(true, true);
  }

  public Color[] getColorScheme() {
    return colorScheme;
  }

  @Override
  public void assignAxisLabels() {
    xAxisLabel = shrunk ? " " : sp.getPlotType(panelIndex).getAxis1();
    yAxisLabel = shrunk ? " " : sp.getPlotType(panelIndex).getAxis2();
  }

  public boolean invertX() {
    return false;
  }

  public boolean invertY() {
    return false;
  }

  public void toggleMasking() {

  }

  @Override
  public void highlightPoints() {
    byte defaultSize;

    defaultSize = sp.getPointSize();
    for (PlotPoint point : points) {
      if (point != null) {
        if (point.isHighlighted()) {
          point.setSize((byte) (defaultSize * 1.5));
        } else {
          point.setSize((defaultSize));
        }
      }
    }
  }

  @Override
  public void generatePoints() {
    int position, markerIndex;
    PLOT_TYPE plotType;
    byte chr, genotypeCode, classCode;
    PointType type;
    float[][] datapoints;
    byte layer;
    IndiPheno indi;
    byte size, xFontSize;
    boolean[] displayCents;
    float[][][][] cents;
    int numCents, count;
    byte centSize;
    float x, y;
    int[] genotype;
    String[] sex;
    String[] otherClass;
    MarkerData markerData;
    int currentClass;
    Hashtable<String, String> disabledClassValues;
    boolean shiftColorOfSexChromosomes;
    String newGenotypingFilename;
    boolean[] isNewGenotypingDifferent = null;
    Logger log;
    int countMissing;
    int numClassCodeLessThanZero = 0;

    log = sp.getProject().getLog();
    disabledClassValues = sp.getDisabledClassValues();// panelIndex);

    shiftColorOfSexChromosomes = sp.getProject().SHIFT_SEX_CHR_COLORS_YESNO.getValue();

    if (!sp.markerDataIsActive()) {
      return;
    }

    plotType = sp.getPlotType(panelIndex);
    currentClass = sp.getCurrentClass(panelIndex);
    if (currentClass < SampleData.BASIC_CLASSES.length
        && SampleData.BASIC_CLASSES[currentClass].equals(SampleData.HEATMAP)) {
      chartType = HEAT_MAP_TYPE;
    } else {
      chartType = SCATTER_PLOT_TYPE;
    }
    markerIndex = sp.getMarkerIndex();// this is the index from the markers loaded perspective
    // gcThreshold = sp.getGCthreshold();
    markerData = sp.getCurrentMarkerData();
    int markerProjectIndex = sp.getMarkerProjectIndices().get(markerData.getMarkerName()); // index
                                                                                           // of
                                                                                           // the
                                                                                           // marker
                                                                                           // in
                                                                                           // the
                                                                                           // project

    boolean[] toInclude = sp.hideExcludedSamples(panelIndex) ? sp.getProject()
                                                                 .getSamplesToInclude(null, false)
                                                             : ArrayUtils.booleanArray(samples.length,
                                                                                       true);
    out: if (plotType == PLOT_TYPE.BAF_LRR && sp.getDisplaygcAdjustor() != null
             && ArrayUtils.booleanArraySum(sp.getDisplaygcAdjustor()) == 1) {
      for (int i = 0; i < sp.getDisplaygcAdjustor().length; i++) {
        if (sp.getDisplaygcAdjustor()[i]) {
          if (sp.getGcAdjustorParameters()[i] == null) {
            sp.getGcAdjustorParameters()[i] = GcAdjustorParameters.readSerial(sp.getGcCList()[1][i],
                                                                              log);
            log.reportTimeInfo("Lazy loading " + sp.getGcCList()[1][i]);
          }
          datapoints = markerData.getGCCorrectedLRRBAF(sp.getGcAdjustorParameters()[i],
                                                       markerProjectIndex, log);
          break out;
        }
      }
      datapoints = markerData.getDatapoints(plotType.getLegacyIndex(),
                                            SexOps.getSampleSex(sp.getProject(),
                                                                SEX_LOAD_TYPE.MAPPED_SEX),
                                            toInclude, false, 1, sp.getGCthreshold(),
                                            sp.getClusterFilterCollection(), true, sp.getPcResids(),
                                            sp.getNumComponents(), 5, sp.getstdevFilter(),
                                            sp.getCorrectionRatio(),
                                            sp.getProject()
                                              .getProperty(sp.getProject().NUM_THREADS),
                                            sp.getCorrection(panelIndex), sp.getProject().getLog());

    } else {

      datapoints = markerData.getDatapoints(plotType.getLegacyIndex(),
                                            SexOps.getSampleSex(sp.getProject(),
                                                                SEX_LOAD_TYPE.MAPPED_SEX),
                                            toInclude, false, 1, sp.getGCthreshold(),
                                            sp.getClusterFilterCollection(), true, sp.getPcResids(),
                                            sp.getNumComponents(), 5, sp.getstdevFilter(),
                                            sp.getCorrectionRatio(),
                                            sp.getProject()
                                              .getProperty(sp.getProject().NUM_THREADS),
                                            sp.getCorrection(panelIndex), sp.getProject().getLog());
    }

    // alleleCounts = markerData[markerIndex].getAB_Genotypes();
    // alleleCounts = sp.getClusterFilterCollection().filterMarker(markerData[markerIndex],
    // sp.getGCthreshold());
    alleleCounts = markerData.getAbGenotypesAfterFilters(sp.getClusterFilterCollection(),
                                                         sp.getMarkerName(), sp.getGCthreshold(),
                                                         log);
    // Project r = sp.getProject();
    newGenotypingFilename = sp.getProject().DATA_DIRECTORY.getValue(false, true)
                            + sp.getMarkerName() + "_newGenotyping.xln";
    if (new File(newGenotypingFilename).exists()) {
      isNewGenotypingDifferent = loadNewGenotyping(newGenotypingFilename);
    }
    // sp.setCurrentClusterFilter(sp.getCurrentClusterFilter()); // what did this patch? this causes
    // a continuous loop
    sp.displayClusterFilterIndex();
    chr = markerData.getChr();
    position = markerData.getPosition();
    size = sp.getPointSize();
    xFontSize = (byte) (size * 2);
    displayCents = sp.getDisplayCentroids();
    cents = sp.getCentroids();
    centSize = 20;

    if (datapoints[0] == null || datapoints[1] == null) {
      errorMessage = "Data not available:";
      points = new PlotPoint[0];
      return;
    } else {
      errorMessage = null;
    }

    boolean checkCents = cents != null;
    for (int i = 0; i < displayCents.length && checkCents; i++) {
      checkCents = cents[i] != null && cents[i][markerIndex] != null
                   && cents[i][markerIndex].length == 3;
    }
    // if (plotType == 0 || plotType == 1 || plotType >= 4) {
    if (plotType != PLOT_TYPE.BAF_LRR && checkCents) {
      numCents = ArrayUtils.booleanArraySum(displayCents);
      points = new PlotPoint[samples.length + numCents * 3];

      count = 0;
      for (int i = 0; i < displayCents.length; i++) {
        if (displayCents[i] && checkCents) {
          for (int j = 0; j < 3; j++) {
            if (cents[i][markerIndex][j] == null) {
              x = 0;
              y = 0;
            } else if (plotType == PLOT_TYPE.X_Y) {
              x = (float) (cents[i][markerIndex][j][1]
                           / (1 + Math.sin(cents[i][markerIndex][j][0] * Math.PI / 2)
                                  / Math.cos(cents[i][markerIndex][j][0] * Math.PI / 2)));
              y = (float) (cents[i][markerIndex][j][1]
                           / (1 + Math.cos(cents[i][markerIndex][j][0] * Math.PI / 2)
                                  / Math.sin(cents[i][markerIndex][j][0] * Math.PI / 2)));
            } else {
              x = cents[i][markerIndex][j][0];
              y = cents[i][markerIndex][j][1];
            }
            if (x > 0 || y > 0) {
              points[count * 3 + j] = new PlotPoint("Centroids", PointType.FILLED_CIRCLE, x, y,
                                                    centSize, (byte) (5 + i), (byte) 10);
            } else {
              points[count * 3 + j] = new PlotPoint("Centroids", PointType.MISSING, x, y, centSize,
                                                    (byte) (5 + i), (byte) 10);
              points[count * 3 + j].setVisible(false);
            }

          }
          count++;
        }
      }
    } else {
      points = new PlotPoint[samples.length];
      numCents = 0;
    }

    // if (plotType < 1 || plotType >= 4) {
    if (plotType == PLOT_TYPE.X_Y) {
      forcePlotXmax = Float.NaN;
    } else {
      forcePlotXmax = 1;
    }

    // dataForQc = new int[3][samples.length];
    // genotype = markerData[markerIndex].getAB_GenotypesAfterFilters(null, sp.getGCthreshold());
    genotype = new int[samples.length];
    sex = new String[samples.length];
    otherClass = new String[samples.length];
    uniqueValueCounts = new CountVector();
    indPtMap = new HashMap<>();
    countMissing = 0;
    boolean[] incl = ArrayUtils.booleanArray(samples.length, true);
    for (int i = 0; i < samples.length; i++) {
      indi = sampleData.getIndiFromSampleHash(samples[i]);

      PlotPoint p = null;
      int index = (numCents * 3) + i;
      if (indi != null && (sp.hideExcludedSamples(panelIndex)
                           && sampleData.individualShouldBeExcluded(samples[i]))) {
        // if sample should be excluded then do nothing
        genotype[i] = -3;
        sex[i] = "e";
        otherClass[i] = "e";
        incl[i] = false;
      } else if (indi != null) {
        genotypeCode = (byte) (alleleCounts[i] + 1);

        // additional genotypeFilters
        if (currentClass == 1) {
          classCode = genotypeCode;
        } else if (sampleData.getClassName(currentClass)
                             .startsWith(SampleData.PLINK_CLASS_PREFIX)) {
          byte indiCode = sp.getPlinkGenotypeForIndi(samples[i], currentClass);// chr1:159,937,288-159,945,728
          classCode = (byte) (indiCode + 1);
        } else {
          classCode = sampleData.determineCodeFromClass(currentClass, alleleCounts[i], indi, chr,
                                                        position);
        }

        if (classCode <= -1 && !sp.maskMissing(panelIndex)) {
          classCode = 0;
        }
        if (Float.isNaN(datapoints[0][i]) || Float.isNaN(datapoints[1][i])) {
          type = PointType.NOT_A_NUMBER;
          // } else if (currentClass==1 && alleleCounts[i]==-1) {
        } else if (sp.getGCthreshold() > 0 && alleleCounts[i] == -1) {
          type = PointType.MISSING;
        } else if (isNewGenotypingDifferent != null && isNewGenotypingDifferent[i]) {
          type = PointType.OPEN_SQUARE;
        } else {
          type = PointType.FILLED_CIRCLE;
        }
        if (classCode == 100) {
          classCode = genotypeCode;
          if (classCode == 0) {
            type = PointType.FILLED_TRIANGLE;
            classCode = (byte) (colorScheme.length - 1);
          }
        }

        layer = (byte) ((sampleData.getClassCategoryAndIndex(currentClass)[0] == 2
                         && classCode > 0) ? 1 : 0);
        layer = classCode; // TODO temporary fix, since was always zero otherwise

        if (type == PointType.NOT_A_NUMBER || type == PointType.MISSING) {
          uniqueValueCounts.add(-1 + "");
          genotype[i] = 0;
        } else {
          uniqueValueCounts.add(classCode + "");
        }
        if (classCode < 0) {
          numClassCodeLessThanZero++;
        }
        if (currentClass < SampleData.BASIC_CLASSES.length
            && SampleData.BASIC_CLASSES[currentClass].equals(SampleData.GENOTYPE) && chr > 22
            && shiftColorOfSexChromosomes) {
          p = new PlotPoint(samples[i], type, datapoints[0][i], datapoints[1][i],
                            type == PointType.FILLED_CIRCLE ? size
                                                            : (type == PointType.FILLED_TRIANGLE ? (byte) (size
                                                                                                           + 5)
                                                                                                 : xFontSize),
                            classCode == 0 ? 0 : (byte) (classCode + 3), layer);
        } else {
          p = new PlotPoint(samples[i], type, datapoints[0][i], datapoints[1][i],
                            type == PointType.FILLED_CIRCLE ? size
                                                            : (type == PointType.FILLED_TRIANGLE ? (byte) (size
                                                                                                           + 5)
                                                                                                 : xFontSize),
                            classCode, layer);
        }
        genotype[i] = genotypeCode;
        if (sampleData.getSexClassIndex() >= 0) {
          int actSexClassIndex = sampleData.getSexClassIndex()
                                 + sampleData.getBasicClasses().length;
          sex[i] = sampleData.determineCodeFromClass(actSexClassIndex, alleleCounts[i], indi, chr,
                                                     position)
                   + "";
        } else {
          sex[i] = "0";
        }
        otherClass[i] = sampleData.determineCodeFromClass(currentClass, alleleCounts[i], indi, chr,
                                                          position)
                        + "";
      } else {
        if (countMissing < 10) {
          log.reportError("Error - no data pts for " + samples[i]);
        } else if (countMissing == 10) {
          log.reportError("...");
        }
        countMissing++;
        sex[i] = "missing";
        p = new PlotPoint(samples[i], PointType.MISSING, datapoints[0][i], datapoints[1][i],
                          (byte) (xFontSize * 2), (byte) 0, (byte) 99);
      }
      if (p != null) {
        points[index] = p;
        indPtMap.put(samples[i], index);
      }

      // create grid
    }
    if (countMissing >= 10) {
      log.reportError("Total of " + countMissing + " samples without data in SampleData");
    }
    if (numClassCodeLessThanZero >= 20) {
      log.reportError("Total of " + numClassCodeLessThanZero
                      + " samples with a missing class code (e.g., less than zero)");
    }

    // callRate=(samples.length-callRate)*100/samples.length;
    if (getUpdateQcPanel()) {
      sp.updateQcPanel(chr, ArrayUtils.subArray(genotype, incl), ArrayUtils.subArray(sex, incl),
                       otherClass, panelIndex);
      setUpdateQCPanel(false);
    }
    sp.updateColorKey(uniqueValueCounts.convertToHash(), panelIndex);

    Hashtable<String, String> hash = new Hashtable<>();
    for (PlotPoint point : points) { // only indi points? (i.e. not centroid points?)
      if (point != null
          && disabledClassValues.containsKey(currentClass + "\t" + point.getColor())) {
        point.setVisible(false);
      }
      if (point != null) {
        hash.put(point.getLayer() + "", "");
      }
    }
    // sort layers as bytes
    byte[] layers = ArrayUtils.toByteArray(HashVec.getKeys(hash, false));
    Arrays.sort(layers);
    setLayersInBase(layers);
    generateRectangles();
    setSwapable(false);
    if (sp.getCurrentClusterFilter() >= 0) {
      rectangles[sp.getCurrentClusterFilter()].setColor((byte) 0);
    }
    if (sp.displayMendelianErrors(panelIndex)) {
      generateLines(sex);
    } else {
      lines = new GenericLine[0];
    }
    // sp.setCurrentClusterFilter(sp.getCurrentClusterFilter()); // what did this patch? this causes
    // a continuous loop
  }

  private void generateLines(String[] sex) {
    ArrayList<GenericLine> linesList = new ArrayList<>();
    byte size = (byte) 3;
    byte momColor = (byte) 6;
    byte dadColor = (byte) 7;
    byte layer = (byte) 1;
    boolean swapAxes = false;

    if (sp.getPedigree() != null) {
      Map<String, MendelErrorCheck> mendelErrorChecks = Pedigree.checkMendelErrors(sp.getPedigree(),
                                                                                   sp.getCurrentMarkerData(),
                                                                                   sp.hideExcludedSamples(panelIndex) ? sp.getProject()
                                                                                                                          .getSamplesToInclude(null,
                                                                                                                                               false)
                                                                                                                      : null,
                                                                                   sex,
                                                                                   sp.getClusterFilterCollection(),
                                                                                   sp.getGCthreshold(),
                                                                                   sp.getProject()
                                                                                     .getLog());
      if (mendelErrorChecks != null) {
        for (int i = 0; i < samples.length; i++) {
          PlotPoint indiPoint = points[indPtMap.get(samples[i])];
          MendelErrorCheck mendelErrorCheck = mendelErrorChecks.get(samples[i]);
          if (mendelErrorCheck == null) {
            continue;
          }
          String[] lookup = sampleData.lookup(samples[i]);
          String fid = lookup[1].split("\t")[0];
          String iid = lookup[1].split("\t")[1];
          int indIndex = sp.getPedigree().getIndIndex(fid, iid);
          if (mendelErrorCheck.hasMoMendelError()) {
            int moIndex = sp.getPedigree().getMoDNAIndex(indIndex);
            PlotPoint momPoint = points[indPtMap.get(samples[moIndex])];
            GenericLine gl = new GenericLine(momPoint, indiPoint, size, momColor, layer, swapAxes,
                                             1, true);
            linesList.add(gl);
          }
          if (mendelErrorCheck.hasFaMendelError()) {
            int faIndex = sp.getPedigree().getFaDNAIndex(indIndex);
            PlotPoint dadPoint = points[indPtMap.get(samples[faIndex])];
            GenericLine gl = new GenericLine(dadPoint, indiPoint, size, dadColor, layer, swapAxes,
                                             1, true);
            linesList.add(gl);
          }

        }
      }
    }

    lines = linesList.toArray(new GenericLine[linesList.size()]);
  }

  private boolean[] loadNewGenotyping(String alternativeGenotypingFilename) {
    Hashtable<String, Byte> alleleCountsNew;
    boolean[] isNewGenotypingDifferent;
    BufferedReader reader;
    String[] line;

    try {
      alleleCountsNew = new Hashtable<>(alleleCounts.length);
      isNewGenotypingDifferent = new boolean[alleleCounts.length];
      reader = new BufferedReader(new FileReader(alternativeGenotypingFilename));
      reader.readLine();
      while (reader.ready()) {
        line = reader.readLine().split("\t");
        alleleCountsNew.put(line[0], Byte.parseByte(line[9]));
      }
      reader.close();

      for (int i = 0; i < samples.length; i++) {
        if (alleleCountsNew.containsKey(samples[i])
            && alleleCounts[i] != alleleCountsNew.get(samples[i])) {
          alleleCounts[i] = alleleCountsNew.get(samples[i]);
          isNewGenotypingDifferent[i] = true;
        }
      }

    } catch (FileNotFoundException fnfe) {
      // log.reportError("Error: file \"" + name + "\" not found in current directory");
      return null;
    } catch (IOException ioe) {
      // log.reportError("Error reading file \"" + name + "\"");
      return null;
    }
    return isNewGenotypingDifferent;
  }

  @Override
  public void mouseMoved(MouseEvent event) {
    Graphics g = getGraphics();
    String pos;
    int x, y;

    float[][] datapoints;
    IndiPheno indi;
    // float[] gcScores;
    int xWidth;
    PLOT_TYPE plotType;
    int currentClass;
    int i;
    byte chr;
    int position;
    byte size, xFontSize;
    MarkerData mData;

    plotType = sp.getPlotType(panelIndex);
    currentClass = sp.getCurrentClass(panelIndex);

    if (!sp.markerDataIsActive()) {
      return;
    }

    x = event.getX();
    y = event.getY();

    pos = (int) Math.floor(x / DEFAULT_LOOKUP_RESOLUTION) + "x"
          + (int) Math.floor(y / DEFAULT_LOOKUP_RESOLUTION);
    if (!pos.equals(prevPos)) {
      repaint();
    }
    indicesOfNearbyPoints = lookupNearbyPoints(x, y, pos);

    mData = sp.getCurrentMarkerData();
    datapoints = mData.getDatapoints(plotType.getLegacyIndex());
    chr = mData.getChr();
    position = mData.getPosition();

    size = sp.getPointSize();
    xFontSize = (byte) (size * 2);

    g.setFont(new Font("Arial", 0, (int) (xFontSize * 1.5)));
    xWidth = g.getFontMetrics(g.getFont()).stringWidth("X");

    for (int l = 0; indicesOfNearbyPoints != null && l < indicesOfNearbyPoints.size(); l++) {
      int i1 = indicesOfNearbyPoints.elementAt(l);
      if (sampIndMap.containsKey(points[i1].getId())) {
        i = sampIndMap.get(points[i1].getId());
        indi = sampleData.getIndiFromSampleHash(samples[i]);
        byte classCode;
        if (sampleData.getClassName(currentClass).startsWith(SampleData.PLINK_CLASS_PREFIX)) {
          classCode = (byte) (sp.getPlinkGenotypeForIndi(samples[i], currentClass) + 1);
        } else {
          classCode = sampleData.determineCodeFromClass(currentClass, alleleCounts[i], indi, chr,
                                                        position);
        }
        g.setColor(colorScheme[Math.max(0, Math.min(classCode, colorScheme.length - 1))]);
        if (sp.getGCthreshold() > 0 && alleleCounts[i] == -1) {
          g.drawString("X", getXPixel(datapoints[0][i]) - xWidth / 2,
                       getYPixel(datapoints[1][i]) + (int) (xFontSize / 2.0));
        } else {
          g.fillOval(getXPixel(datapoints[0][i]) - size * 2 / 2,
                     getYPixel(datapoints[1][i]) - size * 2 / 2, size * 2, size * 2);
        }
      }
    }
    prevPos = pos;
  }

  @Override
  public void mousePressed(MouseEvent e) {
    if (SwingUtilities.isLeftMouseButton(e) && !e.isControlDown()) {
      mouseStartX = e.getX();
      mouseStartY = e.getY();
    } else {
      super.mousePressed(e);
    }
  }

  @Override
  public void mouseReleased(MouseEvent e) {
    mouseEndX = e.getX();
    mouseEndY = e.getY();
    highlightRectangle = null;
    if (SwingUtilities.isLeftMouseButton(e) && !e.isControlDown()) {
      if (Math.abs(mouseEndX - mouseStartX) > (sp.getPointSize() / 2)
          || Math.abs(mouseEndY - mouseStartY) > (sp.getPointSize() / 2)) {
        // Check if custom GC is masking any markers that should be included
        float gcCurr = sp.getGCthreshold();
        float gcNonZeroMin = sp.findMinNonZeroGc();
        boolean drawFilter = true;
        if (gcCurr > gcNonZeroMin) {
          int select = JOptionPane.showConfirmDialog(sp.getFocusOwner(),
                                                     "A cluster filter was created while the GC slider was set to "
                                                                         + gcCurr
                                                                         + "; \nsince there were samples with genotypes with GC values less than this value (minimum was "
                                                                         + gcNonZeroMin
                                                                         + "), \nthe GC slider bar has been set to this minimum value, so that you can see all genotypes that may be exported. \nThis is to ensure that the user exports exactly what they saw when they manually reclustered the marker.",
                                                     "Custom GC Value",
                                                     JOptionPane.OK_CANCEL_OPTION);
          if (select == JOptionPane.CANCEL_OPTION) {
            drawFilter = false;
          } else {
            sp.updateGCSlider((int) (gcNonZeroMin * 100));
          }
        }

        if (drawFilter) {
          // Automatically predict the new genotype and assigns to the last filter.
          sp.getClusterFilterCollection()
            .addClusterFilter(sp.getMarkerName(),
                              new ClusterFilter((byte) sp.getPlotType(panelIndex).getLegacyIndex(),
                                                (float) Math.max(plotXmin,
                                                                 Math.min(getXValueFromXPixel(mouseStartX),
                                                                          getXValueFromXPixel(mouseEndX))),
                                                (float) Math.max(plotYmin,
                                                                 Math.min(getYValueFromYPixel(mouseStartY),
                                                                          getYValueFromYPixel(mouseEndY))),
                                                (float) Math.min(plotXmax,
                                                                 Math.max(getXValueFromXPixel(mouseStartX),
                                                                          getXValueFromXPixel(mouseEndX))),
                                                (float) Math.min(plotYmax,
                                                                 Math.max(getYValueFromYPixel(mouseStartY),
                                                                          getYValueFromYPixel(mouseEndY))),
                                                sp.getCurrentMarkerData()));
          // sp.startAutoSaveToTempFile();
          sp.setClusterFilterUpdated(true);
          setPointsGeneratable(true);
          setUpdateQCPanel(true);
          generateRectangles();
          sp.setCurrentClusterFilter((byte) (sp.getClusterFilterCollection()
                                               .getSize(sp.getMarkerName())
                                             - 1));
          sp.displayClusterFilterIndex();
          paintAgain();
        }
      }
    } else {
      super.mouseReleased(e);
    }

    // TODO Save the filters into files on hard drives;
  }

  @Override
  public void mouseDragged(MouseEvent e) {
    if (SwingUtilities.isLeftMouseButton(e) && !e.isControlDown()) {
      ClusterFilter clusterFilter;
      mouseEndX = e.getX();
      mouseEndY = e.getY();
      highlightRectangle = new GenericRectangle((float) getXValueFromXPixel(mouseStartX),
                                                (float) getYValueFromYPixel(mouseStartY),
                                                (float) getXValueFromXPixel(mouseEndX),
                                                (float) getYValueFromYPixel(mouseEndY), (byte) 1,
                                                false, false, (byte) 0, (byte) 99, true);

      clusterFilter = new ClusterFilter((byte) sp.getPlotType(panelIndex).getLegacyIndex(),
                                        (float) Math.min(getXValueFromXPixel(mouseStartX),
                                                         getXValueFromXPixel(mouseEndX)),
                                        (float) Math.min(getYValueFromYPixel(mouseStartY),
                                                         getYValueFromYPixel(mouseEndY)),
                                        (float) Math.max(getXValueFromXPixel(mouseStartX),
                                                         getXValueFromXPixel(mouseEndX)),
                                        (float) Math.max(getYValueFromYPixel(mouseStartY),
                                                         getYValueFromYPixel(mouseEndY)),
                                        (byte) 0);
      boolean[] highlight = new boolean[points.length];
      boolean[] markerHigh = sp.getCurrentMarkerData().getHighlightStatus(clusterFilter);
      for (int i = 0; i < samples.length; i++) {
        Integer ind = indPtMap.get(samples[i]);
        if (ind == null) {
          // hidden / excluded sample
        } else {
          highlight[ind.intValue()] = markerHigh[i];
        }
      }
      highlightPoints(highlight);
      setExtraLayersVisible(new byte[] {99});
      repaint();
    } else {
      super.mouseDragged(e);
    }
  }

  @Override
  public void mouseClicked(MouseEvent event) {
    JPopupMenu menu;
    MarkerData mData;
    String markerPosition;
    int window, position;
    int numberToInclude, currentClass;
    byte newClusterFilter, chr;

    if (SwingUtilities.isRightMouseButton(event)) {
      window = sp.getProject().WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER.getValue();
      mData = sp.getCurrentMarkerData();
      markerPosition = "chr" + mData.getChr() + ":" + (mData.getPosition() - window) + "-"
                       + (mData.getPosition() + window);
      currentClass = sp.getCurrentClass(panelIndex);
      chr = mData.getChr();
      position = mData.getPosition();
      if (indicesOfNearbyPoints != null && indicesOfNearbyPoints.size() > 0) {
        menu = new JPopupMenu();
        numberToInclude = Math.min(50, indicesOfNearbyPoints.size());
        for (int i = 0; i < numberToInclude; i++) {
          if (!sampIndMap.containsKey(points[indicesOfNearbyPoints.elementAt(i)].getId())) {
            continue;
          }
          int sampleIndex = sampIndMap.get(points[indicesOfNearbyPoints.elementAt(i)].getId());
          IndiPheno indi = sampleData.getIndiFromSampleHash(samples[sampleIndex]);
          byte classCode;
          if (sampleData.getClassName(currentClass).startsWith(SampleData.PLINK_CLASS_PREFIX)) {
            classCode = (byte) (sp.getPlinkGenotypeForIndi(samples[sampleIndex], currentClass) + 1);
          } else {
            classCode = sampleData.determineCodeFromClass(currentClass, alleleCounts[sampleIndex],
                                                          indi, chr, position);
          }
          menu.add(new LaunchAction(sp.getProject(), samples[sampleIndex], markerPosition,
                                    colorScheme[Math.max(0, Math.min(classCode,
                                                                     colorScheme.length - 1))]));
        }
        if (indicesOfNearbyPoints.size() > 50) {
          menu.add(new LaunchAction("Plus " + (indicesOfNearbyPoints.size() - 50)
                                    + " additional samples"));
        }

        menu.show(this, event.getX(), event.getY());
      }

    } else if (SwingUtilities.isLeftMouseButton(event)) {
      newClusterFilter = lookupNearbyRectangles(event.getX(), event.getY());
      if (newClusterFilter >= 0) {
        ClusterFilterCollection clusterFilterCollection;
        clusterFilterCollection = sp.getClusterFilterCollection();
        clusterFilterCollection.deleteClusterFilter(sp.getMarkerName(), newClusterFilter);
        sp.setCurrentClusterFilter((byte) Math.min(newClusterFilter,
                                                   clusterFilterCollection.getSize(sp.getMarkerName())
                                                                     - 1));
        sp.setClusterFilterUpdated(true);
        sp.displayClusterFilterIndex();
        setPointsGeneratable(true);
        setUpdateQCPanel(true);
        generateRectangles();
        sp.updateGUI();
      }

    }
  }

  public void setUpdateQCPanel(boolean updateQcPanel) {
    this.updateQcPanel = updateQcPanel;
  }

  public boolean getUpdateQcPanel() {
    return updateQcPanel;
  }

  public void generateRectangles() {
    rectangles = sp.getClusterFilterCollection()
                   .getRectangles(sp.getMarkerName(),
                                  (byte) sp.getPlotType(panelIndex).getLegacyIndex(), (byte) 1,
                                  false, false, (byte) 7, (byte) 99);
  }

  public GenericRectangle[] getRectangles() {
    return rectangles;
  }

  // public void setCurrentClass (byte newCurrentClass) {
  // currentClass = newCurrentClass;
  // }
  //

  // public boolean isCNV () {
  // DoubleVector x = new DoubleVector();
  // DoubleVector y = new DoubleVector();
  // float[][] datapoints;
  //
  // if (sp.getPlotType()==1) {
  // datapoints = markerData[sp.getMarkerIndex()].getDatapoints(sp.getPlotType());
  // for (int i=0; i<alleleCounts.length; i++) {
  // if (alleleCounts[i]==0 || alleleCounts[i]==1) {
  // x.add(datapoints[0][i]);
  // } else if (alleleCounts[i]==2 || alleleCounts[i]==1) {
  // y.add(datapoints[1][i]);
  // }
  // }
  // if (Array.isBimodal(x.toArray()) || Array.isBimodal(y.toArray())) {
  // return true;
  // } else {
  // return false;
  // }
  // }
  // return (Boolean) null;
  // }

}
