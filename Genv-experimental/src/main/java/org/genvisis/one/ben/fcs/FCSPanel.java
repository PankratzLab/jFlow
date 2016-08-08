package org.genvisis.one.ben.fcs;

import java.awt.Color;
import java.awt.Dialog.ModalityType;
import java.awt.MouseInfo;
import java.awt.PointerInfo;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.AffineTransform;
import java.awt.geom.Area;
import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.HashSet;

import javax.swing.JCheckBoxMenuItem;
import javax.swing.JColorChooser;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPopupMenu;
import javax.swing.SwingUtilities;
import javax.swing.colorchooser.AbstractColorChooserPanel;
import javax.swing.colorchooser.ColorSelectionModel;
import javax.swing.colorchooser.DefaultColorSelectionModel;
import javax.swing.event.ChangeListener;

import org.genvisis.cnv.plots.GenericLine;
import org.genvisis.cnv.plots.GenericPath;
import org.genvisis.cnv.plots.GenericRectangle;
import org.genvisis.cnv.plots.PlotPoint;
import org.genvisis.common.Array;
import org.genvisis.common.Numbers;
import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.gating.Gate;
import org.genvisis.one.ben.fcs.gating.Gate.PolygonGate;
import org.genvisis.one.ben.fcs.gating.Gate.RectangleGate;
import org.genvisis.one.ben.fcs.gating.GateDimension;
import org.genvisis.one.ben.fcs.gating.GateDimension.RectangleGateDimension;
import org.genvisis.one.ben.fcs.sub.PolygonGateEditor;
import org.genvisis.one.ben.fcs.sub.RectangleGateEditor;

public class FCSPanel extends AbstractPanel2 implements MouseListener, MouseMotionListener {
    public static final long serialVersionUID = 3L;
    public static final int LOOKUP_RESOLUTION = 20;
    public static final Color[] DEFAULT_COLORS = { new Color(33, 31, 53), // dark dark
                                                    new Color(201, 30, 10), // deep red
                                                    new Color(94, 88, 214), // light purple
                                                    new Color(189, 243, 61), // light green
                                                    new Color(217, 109, 194), // pink
                                                    new Color(33, 87, 0), // dark green
                                                    new Color(23, 58, 172), // dark blue
                                                    new Color(140, 20, 180), // deep purple
                                                    new Color(182, 182, 182), // light grey
                                                    new Color(220, 220, 220), // very light grey
                                                    new Color(0, 0, 128), // ALL KINDS OF BLUES
                                                    new Color(55, 129, 252), // light blue
                                                    new Color(100, 149, 237), 
                                                    new Color(72, 61, 139), 
                                                    new Color(106, 90, 205), 
                                                    new Color(123, 104, 238), 
                                                    new Color(132, 112, 255), 
                                                    new Color(0, 0, 205), 
                                                    new Color(65, 105, 225), 
                                                    new Color(0, 0, 255), 
                                                    new Color(30, 144, 255), 
                                                    new Color(0, 191, 255),
                                                    new Color(135, 206, 250), 
                                                    new Color(135, 206, 250), 
                                                    new Color(70, 130, 180), 
                                                    new Color(176, 196, 222), 
                                                    new Color(173, 216, 230), 
                                                    new Color(176, 224, 230), 
                                                    new Color(175, 238, 238), 
                                                    new Color(0, 206, 209), 
                                                    new Color(72, 209, 204),
                                                    new Color(64, 224, 208),
                                                    new Color(0, 255, 255), 
                                                    new Color(224, 255, 255),

    };
    private static final byte POINT_SIZE = 3;

    protected FCSPlot fcp;

    public FCSPanel(FCSPlot fcsPlot) {
        super();
        setDoubleBuffered(false);
        createLookup(false);

        this.fcp = fcsPlot;
        this.setAxisFontSize(24);
        this.setSymmetricAxes(false);
        setZoomable(true, true);
        setLayersInBase(new byte[] { 0 });

        setColorScheme(DEFAULT_COLORS);

        setNullMessage("Select two variables to plot");

    }

    int dataCount = -1;
    String xCol = null, yCol = null;
    double xMed = Double.NaN, xMin = Double.NaN, xMax = Double.NaN, yMed = Double.NaN, yMin = Double.NaN, yMax = Double.NaN, xSD = Double.NaN, ySD = Double.NaN;
    PLOT_TYPE prevPlotType;
    boolean[] showMedSD = { false, false, false, false }; // xMed, xSD, yMed, ySD
    double[] xData;
    double[] yData;
    ArrayList<GenericRectangle> rects = new ArrayList<GenericRectangle>();
    ArrayList<GenericPath> polys = new ArrayList<GenericPath>();
    ArrayList<boolean[]> gating = new ArrayList<boolean[]>();
    ArrayList<Gate> gates = new ArrayList<Gate>();

    ArrayList<GenericLine> histLines = new ArrayList<GenericLine>();
    ArrayList<double[]> tempPoly = new ArrayList<double[]>();
    volatile boolean forceGatesChanged = false;
    volatile boolean lackingData = true;

    public static enum GATING_TOOL {

        RECT_GATE("Rectangle"), POLY_GATE("Polygon");

        private String displayName;

        private GATING_TOOL(String disp) {
            this.displayName = disp;
        }

        public String getDisplayName() {
            return displayName;
        }

        static GATING_TOOL getGatingToolByDisplayName(String disp) {
            for (GATING_TOOL g : values()) {
                if (g.getDisplayName().equals(disp)) {
                    return g;
                }
            }
            return null;
        }
    };

    private volatile GATING_TOOL currentTool = GATING_TOOL.RECT_GATE;

    public void setGatingTool(GATING_TOOL tool) {
        tempPoly.clear();
        highlightPoly = null;
        highlightRectangle = null;
        currentTool = tool;
        paintAgain();
    }

    private boolean isHistogram() {
        return yCol != null && yCol.equals(FCSPlot.HISTOGRAM_COL);
    }

    private void setLines() {
        ArrayList<GenericLine> allLines = new ArrayList<GenericLine>();
        allLines.addAll(histLines);
        if (highlightRectangle != null) {
            allLines.add(new GenericLine(highlightRectangle.getStartXValue(), (float) plotYmin, highlightRectangle.getStartXValue(), (float) plotYmax, (byte) 1, (byte) 1, (byte) 99));
            allLines.add(new GenericLine(highlightRectangle.getStopXValue(), (float) plotYmin, highlightRectangle.getStopXValue(), (float) plotYmax, (byte) 1, (byte) 1, (byte) 99));
        }
        for (Gate g : fcp.getGatingForCurrentPlot()) {
            if (!(g instanceof RectangleGate))
                continue;
            RectangleGate rect = (RectangleGate) g;
            float xMin, xMax;
            xMin = rect.getDimension(xCol).getMin();
            xMax = rect.getDimension(xCol).getMax();
            if (!Numbers.isFinite(xMin)) {
                xMin = Integer.MIN_VALUE;
            }
            if (!Numbers.isFinite(xMax)) {
                xMax = Integer.MAX_VALUE;
            }
            allLines.add(new GenericLine(xMin, 
                                        (float) plotYmin, 
                                        xMin, 
                                        (float) plotYmax, 
                                        (byte) 1, (byte) 1, (byte) 99));
            allLines.add(new GenericLine(xMax, 
                                        (float) plotYmin, 
                                        xMax, 
                                        (float) plotYmax, 
                                        (byte) 1, (byte) 1, (byte) 99));
        }

        lines = allLines.toArray(new GenericLine[allLines.size()]);
    }

    private void resetAllForLackingData() {
        gates.clear();
        rects.clear();
        polys.clear();
        tempPoly.clear();
        histLines.clear();
        points = new PlotPoint[0];
        rectangles = new GenericRectangle[0];
        lines = new GenericLine[0];
        polygons = new GenericPath[0];
    }

    public void generatePointsRectanglesAndLines() {
        byte type;
        float xAxisValue, yAxisValue;
        byte size = POINT_SIZE;

        if (fcp.dataLoader == null) {
            lackingData = true;
            setNullMessage("Please load an FCS file..");
            resetAllForLackingData();
            return;
        }
        if (!fcp.isCurrentDataDisplayable()) {
            lackingData = true;
            setNullMessage("Please wait, data is loading...");
            resetAllForLackingData();
            return;
        }
        lackingData = false;

        boolean columnsChangedX = false;
        boolean columnsChangedY = false;
        boolean dataChanged = false;
        boolean optionsChanged = false;
        boolean gatesChanged = false;
        boolean typeChanged = false;
        {
            String newX = fcp.getXDataName();
            String newY = fcp.getYDataName();
            if (xCol == null || !fcp.getXDataName().equals(xCol)) {
                columnsChangedX = true;
            }
            if (yCol == null || !fcp.getYDataName().equals(yCol)) {
                columnsChangedY = true;
            }
            xCol = newX;
            yCol = newY;
        }
        {
            int count = fcp.getDataCount();
            if (count != dataCount) {
                dataChanged = true;
            }
            dataCount = count;
        }
        {
            PLOT_TYPE plotType = fcp.getPlotType();
            if (plotType != prevPlotType) {
                typeChanged = true;
            }
            prevPlotType = plotType;
        }
        {
            boolean mX = fcp.showMedian(false), mY = fcp.showMedian(true), sdX = fcp.showSD(false), sdY = fcp.showSD(true);
            if (mX != showMedSD[0] || sdX != showMedSD[1] || mY != showMedSD[2] || sdY != showMedSD[3]) {
                optionsChanged = true;
            }
            showMedSD[0] = mX;
            showMedSD[1] = sdX;
            showMedSD[2] = mY;
            showMedSD[3] = sdY;
        }

        {
            gatesChanged = rectangles != null && rectangles.length != rects.size() || polygons != null && polygons.length != polys.size() || forceGatesChanged;
            forceGatesChanged = false;
        }

        // TODO replace with loading gates from FCSPlot (and creating gates from mouse actions)
        rectangles = rects.toArray(new GenericRectangle[rects.size()]);
        polygons = polys.toArray(new GenericPath[polys.size()]);

        if (isHistogram()) {
            setLines();
        }

        // setForcePlotXMin(0);
        // setForcePlotYMin(0);

        boolean skip = !columnsChangedX && !columnsChangedY && !dataChanged && !optionsChanged && !gatesChanged/* && !typeChanged /* don't need to regen if only type has changed, for now */;
        if (skip)
            return;

        if (columnsChangedX || columnsChangedY || dataChanged || gatesChanged) {
            updateGating();
        }

        xData = columnsChangedX || dataChanged || xData == null ? fcp.getAxisData(false, true) : xData;
        yData = columnsChangedY || dataChanged || yData == null ? isHistogram() ? null : fcp.getAxisData(false, false) : yData;

        if (isHistogram()) {
            points = new PlotPoint[0];
            if (!columnsChangedX && !dataChanged && histLines != null && histLines.size() > 0)
                return;

            double[] minMax = Array.minMax(xData);
            int range = (int) Math.ceil(minMax[1]) - (int) Math.floor(minMax[0]) + 1;
//             int step = (int) Math.ceil(range / (2 * Array.iqrExclusive(xData) / Math.pow(xData.length, 1/3))); // bin size too small with this formula
            int step = Math.max(2, (int) (range / (2 * Math.sqrt(xData.length))));
            System.out.println("step: " + step);
            float[] histData = new float[range / step + 1];
            for (double x : xData) {
                histData[(int) (x - minMax[0]) / step]++;
            }

            histLines = new ArrayList<GenericLine>();
            for (int i = 0; i < histData.length - 1; i++) {
                histLines.add(new GenericLine((float) ((i * step) + minMax[0]), histData[i], (float) (((i + 1) * step) + minMax[0]), histData[i + 1], (byte) 2, (byte) 0, (byte) 0));
            }
            setLines();

            return;
        }

        ArrayList<GenericLine> lineList = new ArrayList<GenericLine>();
        if (showMedSD[0] || showMedSD[1]) {
            xMed = columnsChangedX || dataChanged || Double.isNaN(xMed) ? Array.median(xData) : xMed;
            xMin = columnsChangedX || dataChanged || Double.isNaN(xMin) ? Math.min(Math.min(0, plotXmin) - xMed, Array.min(xData) - xMed) : xMin;
            xMax = columnsChangedX || dataChanged || Double.isNaN(xMax) ? Math.max(this.plotXmax + xMed, Array.max(xData) + xMed) : xMax;
            if (showMedSD[0]) {
                lineList.add(new GenericLine((float) xMed, (float) xMin, (float) xMed, (float) xMax, (byte) 1, (byte) 8, (byte) 1, 0, false));
            }
            if (showMedSD[1]) {
                xSD = columnsChangedX || dataChanged || Double.isNaN(xSD) ? Array.stdev(xData, false) : xSD;
                lineList.add(new GenericLine((float) (xMed - xSD), (float) xMin, (float) (xMed - xSD), (float) xMax, (byte) 1, (byte) 9, (byte) 1, 0, false));
                lineList.add(new GenericLine((float) (xMed + xSD), (float) xMin, (float) (xMed + xSD), (float) xMax, (byte) 1, (byte) 9, (byte) 1, 0, false));
            }
        }
        if (showMedSD[2] || showMedSD[3]) {
            yMed = columnsChangedY || dataChanged || Double.isNaN(yMed) ? Array.median(yData) : yMed;
            yMin = columnsChangedY || dataChanged || Double.isNaN(yMax) ? Math.min(Math.min(0, plotYmin) - yMed, Array.min(yData)) : yMin;
            yMax = columnsChangedY || dataChanged || Double.isNaN(yMin) ? Math.max(this.plotYmax + yMed, Array.max(yData)) : yMax;
            if (showMedSD[2]) {
                lineList.add(new GenericLine((float) yMin, (float) yMed, (float) yMax, (float) yMed, (byte) 1, (byte) 8, (byte) 1, 0, false));
            }
            if (showMedSD[3]) {
                ySD = columnsChangedY || dataChanged || Double.isNaN(ySD) ? Array.stdev(yData, false) : ySD;
                lineList.add(new GenericLine((float) yMin, (float) (yMed - ySD), (float) yMax, (float) (yMed - ySD), (byte) 1, (byte) 9, (byte) 1, 0, false));
                lineList.add(new GenericLine((float) yMin, (float) (yMed + ySD), (float) yMax, (float) (yMed + ySD), (byte) 1, (byte) 9, (byte) 1, 0, false));
            }

        }

        lines = lineList.toArray(new GenericLine[lineList.size()]);
        lineList = null;

        byte color = 0;
        if (columnsChangedX || columnsChangedY || dataChanged) {
            points = new PlotPoint[dataCount];
            for (int i = 0; i < points.length; i++) {
                xAxisValue = (float) xData[i];
                yAxisValue = (float) yData[i];
                if (Float.isNaN(xAxisValue) || Float.isNaN(yAxisValue)) {
                    type = PlotPoint.NOT_A_NUMBER;
                } else {
                    type = PlotPoint.FILLED_CIRCLE;
                }

                color = (byte) 0;
                points[i] = new PlotPoint(i + "", type, xAxisValue, yAxisValue, size, color, (byte) 0);
            }
        }
        if (gatesChanged) {
            updateGateColor();
        }
    }

    private void refreshNonBaseLayers() {
        updateGating();
        if (!isHistogram()) {
            updateGateColor();
        }
        paintAgain();
        // repaint();
    }

    private void updateGating() {
        gates.clear();
        gating.clear();
        rects.clear();
        polys.clear();

        ArrayList<Gate> gatesForPlot = fcp.getGatingForCurrentPlot();

        for (Gate g : gatesForPlot) {
            boolean[] gt = g.gate(fcp.dataLoader);
            if (gt == null)
                continue;
            int sm = Array.booleanArraySum(gt);
            int sm1 = Array.booleanArraySum(g.getParentGating(fcp.dataLoader));
            float pctInt = 100 * ((float) sm / (float) sm1);
            String pct = ext.formDeci(pctInt, 2);
            String lbl = (g.getName() == null || "".equals(g.getName()) ? g.getID() : g.getName()) + "\n" + "(" + pct + "%)";
            gates.add(g);
            gating.add(gt);
            if (g instanceof RectangleGate) {
                RectangleGate rg = (RectangleGate) g;
                RectangleGateDimension rgdX = rg.getDimension(xCol);
                RectangleGateDimension rgdY = isHistogram() ? null : rg.getDimension(yCol);
                boolean editable = selectedGates.contains(g) || mouseGates.contains(g) || draggingVertexRects.contains(g);
                float xMin, xMax, yMin, yMax;
                xMin = !Numbers.isFinite(rgdX.getMin()) ? Integer.MIN_VALUE : rgdX.getMin();
                xMax = !Numbers.isFinite(rgdX.getMax()) ? Integer.MAX_VALUE : rgdX.getMax();
                yMin = isHistogram() || !Numbers.isFinite(rgdY.getMin()) ? Integer.MIN_VALUE : rgdY.getMin();
                yMax = isHistogram() || !Numbers.isFinite(rgdY.getMax()) ? Integer.MAX_VALUE : rgdY.getMax();
                rects.add(new GenericRectangle(lbl, xMin, yMin, xMax, yMax, (byte) 1, false, false, (byte) 0, (byte) 99, editable));
            } else if (g instanceof PolygonGate) {
                boolean editable = selectedGates.contains(g) || mouseGates.contains(g) || draggingPolys.contains(g);
                polys.add(new GenericPath(lbl, ((PolygonGate) g).getPath(), (byte) 0, (byte) 0, (byte) 99, false, editable));
            }
        }

        rectangles = rects.toArray(new GenericRectangle[rects.size()]);
        polygons = polys.toArray(new GenericPath[polys.size()]);
    }

    private void updateGateColor() {
        byte color = 0;
        boolean[] parentGating = fcp.getParentGate() == null ? null : fcp.getParentGate().gate(fcp.dataLoader);

        for (int i = 0, count = parentGating == null ? points.length : parentGating.length, index = 0; i < count; i++) {
            if (parentGating != null) {
                if (!parentGating[i]) {
                    continue;
                }
            }
            color = (byte) 0;
            for (int g = 0; g < gating.size(); g++) {
                if (gating.get(g)[i]) {
                    color = (byte) gates.get(g).getColorCode();
                    break;
                }
            }
            points[index++].setColor(color);
        }
    }

    private ArrayList<PolygonGate> getPolysWithVerticesNearClick(MouseEvent e) {
        int tempX = e.getX();
        int tempY = e.getY();

        ArrayList<Gate> gating = fcp.getGatingForCurrentPlot();
        ArrayList<PolygonGate> retPolys = new ArrayList<PolygonGate>();

        double[] coords = new double[6];
        for (int i = 0; i < gating.size(); i++) {
            if (gating.get(i) instanceof PolygonGate) {
                Path2D path = ((PolygonGate) gating.get(i)).getPath();
                PathIterator pi = path.getPathIterator(null);
                while (!pi.isDone()) {
                    pi.currentSegment(coords);
                    if (Math.abs(tempX - getXPixel(coords[0])) < DEFAULT_NEARBY_DIST && Math.abs(tempY - getYPixel(coords[1])) < DEFAULT_NEARBY_DIST) {
                        retPolys.add((PolygonGate) gating.get(i));
                        break;
                    }
                    pi.next();
                }
            }
        }

        return retPolys;
    }

    private ArrayList<PolygonGate> getPolysContainingClick(MouseEvent e) {
        int tempX = e.getX();
        int tempY = e.getY();

        ArrayList<Gate> gating = fcp.getGatingForCurrentPlot();
        ArrayList<PolygonGate> retPolys = new ArrayList<PolygonGate>();

        for (int i = 0; i < gating.size(); i++) {
            if (gating.get(i) instanceof PolygonGate) {
                Path2D path = ((PolygonGate) gating.get(i)).getPath();
                if (path.contains(getXValueFromXPixel(tempX), getYValueFromYPixel(tempY))) {
                    retPolys.add((PolygonGate) gating.get(i));
                }
            }
        }

        return retPolys;
    }

    private ArrayList<RectangleGate> getRectsWithVerticesNearClick(MouseEvent e) {
        int tempX = e.getX();
        int tempY = e.getY();

        ArrayList<Gate> gating = fcp.getGatingForCurrentPlot();
        ArrayList<RectangleGate> retRects = new ArrayList<RectangleGate>();

        for (int i = 0; i < gating.size(); i++) {
            if (gating.get(i) instanceof RectangleGate) {
                RectangleGate rect = (RectangleGate) gating.get(i);
                RectangleGateDimension gdX = rect.getDimension(xCol);
                RectangleGateDimension gdY = rect.getDimension(yCol);
                boolean aX = Numbers.isFinite(gdX.getMin()) && Math.abs(getXPixel(gdX.getMin()) - tempX) < DEFAULT_NEARBY_DIST;
                boolean aY = gdY == null || (Numbers.isFinite(gdY.getMin()) && Math.abs(getYPixel(gdY.getMin()) - tempY) < DEFAULT_NEARBY_DIST);
                boolean bX = Numbers.isFinite(gdX.getMax()) && Math.abs(getXPixel(gdX.getMax()) - tempX) < DEFAULT_NEARBY_DIST;
                boolean bY = gdY == null || (Numbers.isFinite(gdY.getMax()) && Math.abs(getYPixel(gdY.getMax()) - tempY) < DEFAULT_NEARBY_DIST);
                if (aX && aY || aX && bY || bX && aY || bX && bY) {
                    retRects.add(rect);
                }
            }
        }
        return retRects;
    }

    private ArrayList<RectangleGate> getRectsContainingClick(MouseEvent e) {
        int tempX = e.getX();
        int tempY = e.getY();

        ArrayList<Gate> gating = fcp.getGatingForCurrentPlot();
        ArrayList<RectangleGate> retRects = new ArrayList<RectangleGate>();

        for (int i = gating.size() - 1; i >= 0; i--) {
            if (gating.get(i) instanceof RectangleGate) {
                RectangleGate rect = (RectangleGate) gating.get(i);
                RectangleGateDimension gdX = rect.getDimension(xCol);
                double xLow, xHigh, yLow, yHigh;
                xLow = !Numbers.isFinite(gdX.getMin()) ? Integer.MIN_VALUE : getXPixel(gdX.getMin());
                xHigh = !Numbers.isFinite(gdX.getMax()) ? Integer.MAX_VALUE : getXPixel(gdX.getMax());
                if (xLow <= tempX && xHigh >= tempX) {
                    if (!isHistogram()) {
                        RectangleGateDimension gdY = rect.getDimension(yCol);
                        yHigh = !Numbers.isFinite(gdY.getMin()) ? Integer.MIN_VALUE : getYPixel(gdY.getMin());
                        yLow = !Numbers.isFinite(gdY.getMax()) ? Integer.MAX_VALUE : getYPixel(gdY.getMax());
                        if (yLow <= tempY && yHigh >= tempY) {
                            retRects.add(rect);
                        }
                    } else {
                        retRects.add(rect);
                    }
                }
            }
        }

        return retRects;
    }

    /*
     * 
     * - starting a new gate clears selections 
     * - clicking a selected gate deselects it 
     * - press and drag on a selected gate drags all selected gates 
     * - event order: press, drag, click(?), release(?) (Or release, click?) 
     * - drag flag? (InDrag) 
     * - if release and not in drag, and not creating new gate, 
     *          either select (actual logic in mousePressed, not mouseReleased) 
     *          or deselect (which has to be in mouseReleased cuz drag) single gate if inside, 
     *          or deselect all if outside 
     * - if release and inDrag, finalize drag and clear flag (either create rect gate

     *  mousePressed():
     * if drawing new poly, click == new point (i.e., set startX, startY, return) 
     * else
     * check mouse location vs all polys
     *     if not within any shapes, deselect all
     *     -   and start drawing next shape (i.e. set startX, startY) [more important in mouseReleased]
     *     else if within shape, select shape (begin drag)
     *     -   followed by either mouseDragged or mouseReleased
     *     -   mouseDragged == drag all selected polys/rects and mouse polys/rects
     *     -   mouseReleased == add poly to selected polys if still within poly and not dragging
     *     if within N pixels of vertex, start vertex drag (don't deselect everything; can drag/edit while preserving multiselect)
     *      
     */

    private static final int DEFAULT_NEARBY_DIST = 4;
    /** mouseGates is altered during the mousePressed event */
    private HashSet<Gate> mouseGates = new HashSet<Gate>();
    /** selectedGates is only removed-from during the mousePressed event */
    private HashSet<Gate> selectedGates = new HashSet<Gate>();

    private ArrayList<RectangleGate> draggingVertexRects = new ArrayList<RectangleGate>();
    private ArrayList<Integer> draggingVertexInds = new ArrayList<Integer>();
    private ArrayList<PolygonGate> draggingPolys = new ArrayList<PolygonGate>();
    private ArrayList<Integer> draggingPolyInds = new ArrayList<Integer>();

    private void setGateLevel(Gate g) {
        Area p = null;
        if (g instanceof PolygonGate) {
            p = new Area(((PolygonGate) g).getPath());
        } else if (g instanceof RectangleGate) {
            RectangleGateDimension rgd1 = ((RectangleGate) g).getDimension(xCol);
            double min, max, minY, maxY;
            min = !Numbers.isFinite(rgd1.getMin()) ? Integer.MIN_VALUE : Math.min(rgd1.getMin(), rgd1.getMax());
            max = !Numbers.isFinite(rgd1.getMax()) ? Integer.MAX_VALUE : Math.max(rgd1.getMin(), rgd1.getMax());
            if (g.getDimensions().size() == 2) {
                RectangleGateDimension rgd2 = ((RectangleGate) g).getDimension(yCol);
                minY = !Numbers.isFinite(rgd2.getMin()) ? Integer.MIN_VALUE : Math.min(rgd2.getMin(), rgd2.getMax());
                maxY = !Numbers.isFinite(rgd2.getMax()) ? Integer.MAX_VALUE : Math.max(rgd2.getMin(), rgd2.getMax());
            } else {
                minY = Integer.MIN_VALUE;
                maxY = Integer.MAX_VALUE;
            }
            p = new Area(new Rectangle2D.Double(min, minY, max - min, maxY - minY));
        }
        ArrayList<Gate> gating = fcp.getGatingForCurrentPlot();
        int lvl = 0;
        for (Gate gP : gating) {
            if (gP instanceof PolygonGate) {
                Area path = new Area(((PolygonGate) gP).getPath());
                path.intersect(p);
                if (!path.isEmpty()) {
                    lvl++;
                }
            } else if (gP instanceof RectangleGate) {
                RectangleGateDimension rgd1 = ((RectangleGate) gP).getDimension(xCol);
                double min, max, minY, maxY;
                min = !Numbers.isFinite(rgd1.getMin()) ? Integer.MIN_VALUE : Math.min(rgd1.getMin(), rgd1.getMax());
                max = !Numbers.isFinite(rgd1.getMax()) ? Integer.MAX_VALUE : Math.max(rgd1.getMin(), rgd1.getMax());
                if (gP.getDimensions().size() == 2) {
                    RectangleGateDimension rgd2 = ((RectangleGate) gP).getDimension(yCol);
                    minY = !Numbers.isFinite(rgd2.getMin()) ? Integer.MIN_VALUE : Math.min(rgd2.getMin(), rgd2.getMax());
                    maxY = !Numbers.isFinite(rgd2.getMax()) ? Integer.MAX_VALUE : Math.max(rgd2.getMin(), rgd2.getMax());
                } else {
                    minY = Integer.MIN_VALUE;
                    maxY = Integer.MAX_VALUE;
                }
                Rectangle2D r = new Rectangle2D.Double(min, minY, max - min, maxY - minY);
                if (p.intersects(r) || p.contains(r)) {
                    lvl++;
                }
            }
        }
        g.setLevel(lvl);
    }

    public void mouseReleased(MouseEvent e) {
        if (SwingUtilities.isLeftMouseButton(e) && !lackingData) {
            int mouseEndX;
            int mouseEndY;
            mouseEndX = e.getX();
            mouseEndY = e.getY();

            boolean didSelect = mouseGates.size() > 0;
            boolean didClear = !didSelect && selectedGates.size() > 0;
            boolean wasDrag = drag;
            if (!drag/* && (currentTool == GATING_TOOL.RECT_TOOL || tempPoly.isEmpty()) */) {
                if (!didSelect) {
                    if (!e.isShiftDown() && !e.isControlDown()) {
                        selectedGates.clear();
                    }
                } else {
                    for (Gate g : mouseGates) {
                        if (selectedGates.contains(g)) {
                            selectedGates.remove(g);
                        } else {
                            selectedGates.add(g);
                        }
                    }
                    mouseGates.clear();
                }
            } else {
                drag = false;
                for (Gate g : mouseGates) {
                    if (selectedGates.contains(g)) {
                    } else {
                        selectedGates.add(g);
                    }
                }
                mouseGates.clear();
                if (!draggingPolys.isEmpty()) {
                    draggingPolys.clear();
                    draggingPolyInds.clear();
                }
                if (draggingVertexRects.isEmpty()) {
                    if (currentTool == GATING_TOOL.RECT_GATE) {
                        highlightRectangle = null;

                        if (Math.abs(mouseEndX - startX) > DEFAULT_NEARBY_DIST) {
                            if (isHistogram() || (Math.abs(mouseEndY - startY) > DEFAULT_NEARBY_DIST)) {
                                String name = null;
                                int msgCode = 0;
                                do {
                                    String msg = "Gate Name:";
                                    if (msgCode == 1) {
                                        msg += " (must provide a name)";
                                    } else if (msgCode == 2) {
                                        msg += " (must be unique)";
                                    }
                                    name = JOptionPane.showInputDialog(fcp, msg, "Add Gate?", JOptionPane.QUESTION_MESSAGE);
                                } while ((msgCode = ("".equals(name) ? 1 : (fcp.duplicateGateName(name) ? 2 : 0))) > 0);
                                if (name != null) {
                                    RectangleGate rg = new RectangleGate(fcp.getParentGate(), name);
                                    rg.addDimension(new GateDimension.RectangleGateDimension(rg, xCol, fcp.getXScale(), (float) getXValueFromXPixel(startX), (float) getXValueFromXPixel(mouseEndX)));
                                    if (!isHistogram()) {
                                        rg.addDimension(new GateDimension.RectangleGateDimension(rg, yCol, fcp.getYScale(), (float) getYValueFromYPixel(startY), (float) getYValueFromYPixel(mouseEndY)));
                                    }
                                    rg.setColor(3);
                                    setGateLevel(rg);
                                    fcp.addGate(rg);
                                    selectedGates.add(rg);
                                }
                            }
                        }
                    }
                } else {
                    selectedGates.addAll(draggingVertexRects);
                    draggingVertexRects.clear();
                    draggingVertexInds.clear();
                }
            }
            if (!didSelect && !didClear && !wasDrag) {
                if (currentTool == GATING_TOOL.POLY_GATE) {
                    double tempValX = getXValueFromXPixel(mouseEndX);
                    double tempValY = getYValueFromYPixel(mouseEndY);
                    if (!tempPoly.isEmpty()) {
                        int initPtX = getXPixel(tempPoly.get(0)[0]);
                        int initPtY = getYPixel(tempPoly.get(0)[1]);
                        if (Math.abs(initPtX - mouseEndX) < 5 && Math.abs(initPtY - mouseEndY) < 5) {
                            if (tempPoly.size() > 2) {
                                String name = null;
                                int msgCode = 0;
                                do {
                                    String msg = "Gate Name:";
                                    if (msgCode == 1) {
                                        msg += " (must provide a name)";
                                    } else if (msgCode == 2) {
                                        msg += " (must be unique)";
                                    }
                                    name = JOptionPane.showInputDialog(fcp, msg, "Add Gate?", JOptionPane.QUESTION_MESSAGE);
                                } while ((msgCode = ("".equals(name) ? 1 : (fcp.duplicateGateName(name) ? 2 : 0))) > 0);
                                if (name != null) {
                                    Path2D path = new Path2D.Double(Path2D.WIND_EVEN_ODD);
                                    path.moveTo(tempPoly.get(0)[0], tempPoly.get(0)[1]);
                                    for (int i = 1; i < tempPoly.size(); ++i) {
                                        path.lineTo(tempPoly.get(i)[0], tempPoly.get(i)[1]);
                                    }
                                    path.closePath();
                                    PolygonGate pg = new PolygonGate(fcp.getParentGate(), name);
                                    pg.addDimension(new GateDimension(pg, xCol, fcp.getXScale()));
                                    pg.addDimension(new GateDimension(pg, yCol, fcp.getYScale()));
                                    pg.setPath(path);
                                    pg.setShouldMimicFlowJoGating(fcp.drawPolysAsFlowJo());
                                    pg.setColor(3);
                                    setGateLevel(pg);
                                    fcp.addGate(pg);
                                    selectedGates.add(pg);
                                }
                            }
                            tempPoly.clear();
                            highlightPoly = null;
                            setForceGatesChanged();
                            refreshNonBaseLayers();
                            return;
                        }
                    }
                    tempPoly.add(new double[] { tempValX, tempValY });
                    Path2D path = new Path2D.Double(Path2D.WIND_EVEN_ODD);
                    path.moveTo(tempPoly.get(0)[0], tempPoly.get(0)[1]);
                    for (int i = 1; i < tempPoly.size(); ++i) {
                        path.lineTo(tempPoly.get(i)[0], tempPoly.get(i)[1]);
                    }
                    highlightPoly = new GenericPath(null, path, (byte) 0, (byte) 0, (byte) 99, false, true);
                }

            }
            setForceGatesChanged();
            refreshNonBaseLayers();
        } else {
            super.mouseReleased(e);
        }
    }

    public void mousePressed(MouseEvent e) {
        if (SwingUtilities.isLeftMouseButton(e) && !lackingData) {
            int tempX = e.getX();
            int tempY = e.getY();
            if (tempPoly.isEmpty()) {
                ArrayList<PolygonGate> insidePolys = getPolysContainingClick(e);
                ArrayList<PolygonGate> closePolys = getPolysWithVerticesNearClick(e);
                ArrayList<RectangleGate> insideRects = getRectsContainingClick(e);
                ArrayList<RectangleGate> closeRects = getRectsWithVerticesNearClick(e);

                ArrayList<Gate> insideGates = new ArrayList<Gate>();
                insideGates.addAll(insidePolys);
                insideGates.addAll(insideRects);
                ArrayList<Gate> closeGates = new ArrayList<Gate>();
                closeGates.addAll(closePolys);
                closeGates.addAll(closeRects);

                if (closeGates.isEmpty()) {
                    if (!insideGates.isEmpty()) {
                        // if shift down, add lowest level gate to mouseRects
                        // if ctrl down, remove clicked rects from mouseRects
                        // if alt down, cycle to gate at next level up (wrap around?)
                        // if no modifiers, clear mouseRects and reselect new click (but only highest level)

                        if (e.isShiftDown() && !e.isControlDown() && !e.isAltDown()) {
                            Gate lowest = null;
                            for (Gate gr : insideGates) {
                                if (!selectedGates.contains(gr)) {
                                    if (lowest == null || lowest.getLevel() > gr.getLevel()) {
                                        lowest = gr;
                                    }
                                }
                            }
                            if (lowest != null) {
                                mouseGates.add(lowest);
                            }
                        } else if (!e.isShiftDown() && e.isControlDown() && !e.isAltDown()) {
                            boolean removed = false;
                            Gate highestSel = null;
                            for (Gate gr : insideGates) {
                                if ((selectedGates.contains(gr) || mouseGates.contains(gr)) && (highestSel == null || highestSel.getLevel() < gr.getLevel())) {
                                    highestSel = gr;
                                }
                            }
                            if (highestSel != null) {
                                selectedGates.remove(highestSel);
                                mouseGates.remove(highestSel);
                                removed = true;
                            }
                            if (!removed) {
                                Gate lowest = null;
                                for (Gate gr : insideGates) {
                                    if (lowest == null || lowest.getLevel() > gr.getLevel()) {
                                        lowest = gr;
                                    }
                                }
                                mouseGates.add(lowest);
                            }
                        } else if (!e.isShiftDown() && !e.isControlDown() && e.isAltDown()) {
                            int lowest = Integer.MAX_VALUE;
                            for (Gate g : mouseGates) {
                                lowest = Math.min(g.getLevel(), lowest);
                            }
                            Gate lowestSel = null;
                            for (Gate gr : insideGates) {
                                if (lowestSel == null || gr.getLevel() < lowest || (lowestSel != null && gr.getLevel() < lowestSel.getLevel())) {
                                    lowestSel = gr;
                                }
                            }
                            mouseGates.clear();
                            mouseGates.add(lowestSel);
                        } else {
                            mouseGates.clear();
                            Gate lowest = null;
                            for (Gate gr : insideGates) {
                                if (lowest == null || lowest.getLevel() > gr.getLevel()) {
                                    lowest = gr;
                                }
                            }
                            if (lowest != null) {
                                mouseGates.add(lowest);
                            }
                        }
                    } else {
                        mouseGates.clear();
                    }
                } else {
                    // RECTS
                    for (RectangleGate rect : closeRects) {
                        RectangleGateDimension gdX = rect.getDimension(xCol);
                        RectangleGateDimension gdY = rect.getDimension(yCol);
                        selectedGates.remove(rect);
                        mouseGates.remove(rect);

                        int ind = -1;
                        boolean closeToStartX = Numbers.isFinite(gdX.getMin()) && Math.abs(getXPixel(gdX.getMin()) - tempX) < 4;
                        boolean closeToStartY = isHistogram() || (Numbers.isFinite(gdY.getMin()) && Math.abs(getYPixel(gdY.getMin()) - tempY) < 4);
                        boolean closeToStopX = Numbers.isFinite(gdX.getMax()) && Math.abs(getXPixel(gdX.getMax()) - tempX) < 4;
                        boolean closeToStopY = isHistogram() || (Numbers.isFinite(gdY.getMax()) && Math.abs(getYPixel(gdY.getMax()) - tempY) < 4);
                        /*
                        1 ----- 2
                        |       |
                        |       |
                        |       |
                        0 ----- 3
                        */
                        if (closeToStartX && closeToStartY) {
                            ind = 0;
                        } else if (closeToStartX && closeToStopY) {
                            ind = 1;
                        } else if (closeToStopX && closeToStartY) {
                            ind = 3;
                        } else if (closeToStopX && closeToStopY) {
                            ind = 2;
                        }
                        if (ind != -1) {
                            draggingVertexRects.add(rect);
                            draggingVertexInds.add(ind);
                        }
                    }

                    // POLYS
                    selectedGates.addAll(closePolys);
                    mouseGates.removeAll(closePolys);

                    draggingPolys = closePolys;
                    draggingPolyInds.clear();

                    double[] coords = new double[6];
                    polyLoop: for (int i = 0; i < closePolys.size(); i++) {
                        Path2D path = closePolys.get(i).getPath();
                        PathIterator pi = path.getPathIterator(null);
                        int vInd = 0;
                        while (!pi.isDone()) {
                            pi.currentSegment(coords);
                            if (Math.abs(tempX - getXPixel(coords[0])) < 4 && Math.abs(tempY - getYPixel(coords[1])) < 4) {
                                draggingPolyInds.add(vInd);
                                continue polyLoop;
                            }
                            vInd++;
                            pi.next();
                        }
                        draggingPolyInds.add(-1); // this shouldn't technically happen, since closePolys only has polys with at least one vertex close to the mouse
                    }
                }
            } else {
                // do nothing, wait for mouseReleased to set new point or create poly
            }
            startX = e.getX();
            startY = e.getY();
            refreshNonBaseLayers();
        } else {
            super.mousePressed(e);
        }
    }

    // public void mouseMoved(MouseEvent e) {
    // super.mouseMoved(e);
    // int tempX = e.getX();
    // int tempY = e.getY();
    // for (int i = 0; i < rects.size(); i++) {
    // GenericRectangle rect = rects.get(i);
    // int xPix1, xPix2, yPix1, yPix2;
    // xPix1 = getXPixel(rect.getStartXValue());
    // xPix2 = getXPixel(rect.getStopXValue());
    // yPix1 = getYPixel(rect.getStartYValue());
    // yPix2 = getYPixel(rect.getStopYValue());
    // Rectangle myRect = new Rectangle(Math.min(xPix1, xPix2) - 10,
    // Math.min(yPix1, yPix2) - 10,
    // Math.max(xPix1, xPix2)
    // - Math.min(xPix1, xPix2) + 10,
    // Math.max(yPix1, yPix2)
    // - Math.min(yPix1, yPix2) + 10);
    // rect.setEditable(myRect.contains(tempX, tempY));
    // }
    // for (int i = 0; i < polys.size(); i++) {
    // GenericPath poly = polys.get(i);
    // Path2D path = ((Path2D)poly.myPath.clone());
    // poly.setEditable(path.contains(getXValueFromXPixel(tempX), getYValueFromYPixel(tempY)));
    // }
    // paintAgain();
    // }

    public void mouseClicked(MouseEvent e) {
        if (lackingData)
            return;
        if (SwingUtilities.isLeftMouseButton(e)) {
            if (e.getClickCount() == 2) {
                ArrayList<Gate> gates = new ArrayList<Gate>();
                gates.addAll(getPolysContainingClick(e));
                gates.addAll(getRectsContainingClick(e));

                Gate lowest = null;
                for (Gate g : gates) {
                    if (lowest == null || g.getLevel() < lowest.getLevel()) {
                        lowest = g;
                    }
                }
                fcp.gateSelected(lowest, true);
            }
        } else {
            tempPoly.clear();
            highlightPoly = null;
            highlightRectangle = null;

            ArrayList<Gate> allGates = new ArrayList<Gate>();
            allGates.addAll(getPolysContainingClick(e));
            allGates.addAll(getRectsContainingClick(e));
            if (allGates.size() == 0) {
                return;
            }
            Gate g = allGates.get(0);
            for (int i = 1; i < allGates.size(); i++) {
                if (allGates.get(i).getLevel() < g.getLevel()) {
                    g = allGates.get(i);
                }
            }
            selectedGates.clear();
            mouseGates.clear();
            selectedGates.add(g);
            refreshNonBaseLayers();
            JPopupMenu jpop = constructPopup(g);
            int x, y;
            x = MouseInfo.getPointerInfo().getLocation().x;
            y = MouseInfo.getPointerInfo().getLocation().y;
            jpop.setInvoker(fcp);
            jpop.setLocation(x, y);
            jpop.setVisible(true);

            // TODO delete right-clicked gate even if not selected? Delete clicked gate and NOT selected gates?
            // deleteSelectedGates();
        }
        paintAgain();
    }

    private JPopupMenu constructPopup(final Gate g) {
        JPopupMenu menu = new JPopupMenu("Gate Options:");

        ActionListener act = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                switch (e.getActionCommand().charAt(0)) {
                    case 'D':
                        int children = g.countAllChildren();
                        int opt = JOptionPane.showConfirmDialog(fcp, "Are you sure you wish to delete this gate" + (children > 0 ? " (and " + children + " child gates)" : "") + "?", "Confirm Delete Gate?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);
                        if (opt == JOptionPane.YES_OPTION) {
                            fcp.deleteGate(g);
                            clearGating();
                        }
                        break;
                    case 'R':
                        String newName = JOptionPane.showInputDialog(fcp, "Enter new name for gate:", "Rename Gate", JOptionPane.QUESTION_MESSAGE);
                        if (newName == null) return;
                        else {
                            if (fcp.duplicateGateName(newName)) {
                                do {
                                    newName = JOptionPane.showInputDialog(fcp, "Enter new name for gate: (must be unique)", "Rename Gate", JOptionPane.QUESTION_MESSAGE);
                                } while (fcp.duplicateGateName(newName));
                            }
                        }
                        if (newName != null) {
                            g.setName(newName);
                            setForceGatesChanged();
                            fcp.refreshGating();
                        }
                        break;
                    case 'C':
                        
                        // TODO colorPicker
                        break;
                    case 'F':
                        ((PolygonGate)g).setShouldMimicFlowJoGating(((JCheckBoxMenuItem)e.getSource()).isSelected());
                        refreshNonBaseLayers();
                        break;
                    case 'E':
                        if (g instanceof RectangleGate) {
                            RectangleGateEditor rge = new RectangleGateEditor();
                            rge.setGate(fcp, (RectangleGate) g);
                            rge.setModal(true);
                            rge.setModalityType(ModalityType.APPLICATION_MODAL);
                            rge.setVisible(true);
                            
                            if (!rge.isCancelled()) {
                                g.setName(rge.getName());
                                RectangleGateDimension xDim = (RectangleGateDimension) g.getDimensions().get(0);
                                xDim.setMin(rge.getXMin());
                                xDim.setMax(rge.getXMax());
                                if (g.getDimensions().size() == 2) {
                                    RectangleGateDimension yDim = (RectangleGateDimension) g.getDimensions().get(1);
                                    yDim.setMin(rge.getYMin());
                                    yDim.setMax(rge.getYMax());
                                }
                                setForceGatesChanged();
                                fcp.refreshGating();
                            }
                        } else if (g instanceof PolygonGate) {
                            PolygonGateEditor pge = new PolygonGateEditor(fcp, (PolygonGate) g);
                            pge.setModal(true);
                            pge.setModalityType(ModalityType.APPLICATION_MODAL);
                            pge.setVisible(true);
                            
                            if (!pge.isCancelled()) {
                                g.setName(pge.getName());
                                ((PolygonGate) g).setPath(pge.getNewPath());
                                ((PolygonGate) g).setShouldMimicFlowJoGating(pge.getMimicFlowJo());
                                setForceGatesChanged();
                                fcp.refreshGating();
                            }
                            
                        }
                        break;
                }
            }
        };

        JMenuItem delete = new JMenuItem("Delete");
        delete.addActionListener(act);
        delete.setActionCommand("D");
        menu.add(delete);

        JMenuItem rename = new JMenuItem("Rename");
        rename.addActionListener(act);
        rename.setActionCommand("R");
        menu.add(rename);

        JMenuItem color = new JMenuItem("Set Color");
        color.addActionListener(act);
        color.setActionCommand("C");
        menu.add(color);

        if (g instanceof PolygonGate) {
            JCheckBoxMenuItem flow = new JCheckBoxMenuItem("Bin Polygon (FlowJo)", ((PolygonGate) g).getMimicsFlowJo());
            flow.addActionListener(act);
            flow.setActionCommand("F");
            menu.add(flow);
        }

        JMenuItem edit = new JMenuItem("Manually Edit");
        edit.addActionListener(act);
        edit.setActionCommand("E");
        menu.add(edit);

        return menu;
    }

    public void deleteSelectedGates() {
        int cnt = selectedGates.size() + mouseGates.size();
        int childCnt = 0;
        for (Gate g : selectedGates) {
            childCnt += g.countAllChildren();
        }
        for (Gate g : mouseGates) {
            childCnt += g.countAllChildren();
        }
        if (cnt == 0) {
            return;
        }
        int opt = JOptionPane.showConfirmDialog(fcp, "Are you sure you wish to delete " + cnt + " gate(s)" + (childCnt > 0 ? " and " + childCnt + " downstream gate(s)?" : "?"), "Confirm Delete Gate?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);
        if (opt == JOptionPane.YES_OPTION) {
            for (Gate g : selectedGates) {
                fcp.deleteGate(g);
            }
            for (Gate g : mouseGates) {
                fcp.deleteGate(g);
            }
            clearGating();
        }
    }

    public void clearGating() {
        selectedGates.clear();
        mouseGates.clear();
        setForceGatesChanged();
        paintAgain();
    }

    public void mouseDragged(MouseEvent e) {
        if (SwingUtilities.isLeftMouseButton(e) && !e.isControlDown() && !lackingData) {
            drag = true;

            boolean updateStartCoords = false;
            boolean updateForceGating = false;

            int mouseEndX;
            int mouseEndY;
            mouseEndX = e.getX();
            mouseEndY = e.getY();

            if (draggingVertexRects.isEmpty() && draggingPolys.isEmpty()) {
                if (selectedGates.isEmpty() && mouseGates.isEmpty()) {
                    if (currentTool == GATING_TOOL.RECT_GATE) {
                        highlightRectangle = new GenericRectangle((float) getXValueFromXPixel(startX), (float) getYValueFromYPixel(startY), (float) getXValueFromXPixel(mouseEndX), isHistogram() ? (float) getYValueFromYPixel(startY) : (float) getYValueFromYPixel(mouseEndY), (byte) 1, false, false, (byte) 0, (byte) 99, true);
                    }
                } else {
                    float dx, dy;
                    dx = (float) (getXValueFromXPixel(mouseEndX) - getXValueFromXPixel(startX));
                    dy = (float) (getYValueFromYPixel(mouseEndY) - getYValueFromYPixel(startY));
                    AffineTransform at = AffineTransform.getTranslateInstance(dx, dy);
                    for (Gate gr : selectedGates) {
                        if (gr instanceof RectangleGate) {
                            RectangleGateDimension gdX = ((RectangleGate) gr).getDimension(xCol);
                            if (Numbers.isFinite(gdX.getMax())) {
                                gdX.setMax(gdX.getMax() + dx);
                            }
                            if (Numbers.isFinite(gdX.getMin())) {
                                gdX.setMin(gdX.getMin() + dx);
                            }
                            if (!isHistogram()) {
                                RectangleGateDimension gdY = ((RectangleGate) gr).getDimension(yCol);
                                if (Numbers.isFinite(gdY.getMin())) {
                                    gdY.setMin(gdY.getMin() + dy);
                                }
                                if (Numbers.isFinite(gdY.getMax())) {
                                    gdY.setMax(gdY.getMax() + dy);
                                }
                            }
                        } else if (gr instanceof PolygonGate) {
                            ((PolygonGate) gr).transform(at);
                        }
                    }
                    for (Gate gr : mouseGates) {
                        if (gr instanceof RectangleGate) {
                            if (!selectedGates.contains(gr)) {
                                RectangleGateDimension gdX = ((RectangleGate) gr).getDimension(xCol);
                                if (Numbers.isFinite(gdX.getMax())) {
                                    gdX.setMax(gdX.getMax() + dx);
                                }
                                if (Numbers.isFinite(gdX.getMin())) {
                                    gdX.setMin(gdX.getMin() + dx);
                                }
                                if (!isHistogram()) {
                                    RectangleGateDimension gdY = ((RectangleGate) gr).getDimension(yCol);
                                    if (Numbers.isFinite(gdY.getMin())) {
                                        gdY.setMin(gdY.getMin() + dy);
                                    }
                                    if (Numbers.isFinite(gdY.getMax())) {
                                        gdY.setMax(gdY.getMax() + dy);
                                    }
                                }
                            }
                        } else if (gr instanceof PolygonGate) {
                            if (!selectedGates.contains(gr)) {
                                ((PolygonGate) gr).transform(at);
                            }
                        }
                    }
                    for (Gate g : selectedGates) {
                        setGateLevel(g);
                    }
                    for (Gate g : mouseGates) {
                        if (g instanceof RectangleGate || !selectedGates.contains(g)) {
                            setGateLevel(g);
                        }
                    }
                    updateStartCoords = true;
                    updateForceGating = true;
                }
            }
            if (!draggingVertexRects.isEmpty()) {
                for (int i = 0; i < draggingVertexRects.size(); i++) {
                    float tempNewX, tempNewY;

                    RectangleGate rg = draggingVertexRects.get(i);
                    RectangleGateDimension rgdX = rg.getDimension(xCol);
                    /*
                    1 ----- 2
                    |       |
                    |       |
                    |       |
                    0 ----- 3
                    */
                    tempNewX = (float) getXValueFromXPixel(mouseEndX);

                    switch (draggingVertexInds.get(i)) {
                        case 0:
                        case 1:
                            rgdX.setMin(tempNewX);
                            break;
                        case 2:
                        case 3:
                            rgdX.setMax(tempNewX);
                            break;
                    }

                    if (!isHistogram()) {
                        RectangleGateDimension rgdY = rg.getDimension(yCol);
                        tempNewY = (float) getYValueFromYPixel(mouseEndY);

                        switch (draggingVertexInds.get(i)) {
                            case 0:
                            case 3:
                                rgdY.setMin(tempNewY);
                                break;
                            case 1:
                            case 2:
                                rgdY.setMax(tempNewY);
                                break;
                        }
                    }
                }
                updateForceGating = true;
            }
            if (!draggingPolys.isEmpty()) {
                for (int i = 0; i < draggingPolys.size(); i++) {
                    PolygonGate gp = draggingPolys.get(i);
                    int vertexInd = draggingPolyInds.get(i);
                    Path2D newPath = new Path2D.Double();
                    PathIterator pi = gp.getPath().getPathIterator(null);
                    double[] coords = new double[6];
                    int v = 0;
                    while (!pi.isDone()) {
                        int code = pi.currentSegment(coords);
                        if (v == vertexInd) {
                            coords[0] = getXValueFromXPixel(mouseEndX);
                            coords[1] = getYValueFromYPixel(mouseEndY);
                        }
                        switch (code) {
                            case PathIterator.SEG_MOVETO:
                                newPath.moveTo(coords[0], coords[1]);
                                break;
                            case PathIterator.SEG_LINETO:
                                newPath.lineTo(coords[0], coords[1]);
                                break;
                            case PathIterator.SEG_CLOSE:
                                newPath.closePath();
                                break;
                        }
                        v++;
                        pi.next();
                    }
                    gp.setPath(newPath);
                }
                updateStartCoords = true;
                updateForceGating = true;
            }

            if (updateStartCoords) {
                startX = e.getX();
                startY = e.getY();
            }
            if (updateForceGating) {
                setForceGatesChanged();
            }
            refreshNonBaseLayers();
        } else {
            super.mouseDragged(e);
        }
    }

    private void rightMouseClickedRect(MouseEvent e) {
        int tempX = e.getX();
        int tempY = e.getY();
        RectangleGate toRemove = null;

        ArrayList<Gate> gating = fcp.getGatingForCurrentPlot();
        for (int i = gating.size() - 1; i >= 0; i--) {
            if (!(gating.get(i) instanceof RectangleGate))
                continue;
            RectangleGate rect = (RectangleGate) gating.get(i);
            RectangleGateDimension gdX = rect.getDimension(xCol);
            RectangleGateDimension gdY = rect.getDimension(yCol);
            double xLow, xHigh, yLow, yHigh;
            xLow = !Numbers.isFinite(gdX.getMin()) ? Integer.MIN_VALUE : getXPixel(gdX.getMin());
            xHigh = !Numbers.isFinite(gdX.getMax()) ? Integer.MAX_VALUE : getXPixel(gdX.getMax());
            yHigh = gdY == null ? plotYmin : (!Numbers.isFinite(gdY.getMin()) ? Integer.MIN_VALUE : getYPixel(gdY.getMin()));
            yLow = gdY == null ? plotYmax : (!Numbers.isFinite(gdY.getMax()) ? Integer.MAX_VALUE : getYPixel(gdY.getMax()));
            if (isHistogram()) {
                if (xLow <= tempX && xHigh >= tempX) {
                    toRemove = rect;
                    break;
                }
            } else {
                if (xLow <= tempX && xHigh >= tempX && yLow <= tempY && yHigh >= tempY) {
                    toRemove = rect;
                    break;
                }
            }
        }
        if (toRemove == null) {
            for (Gate rect : selectedGates) {
                if (rect instanceof RectangleGate) {
                    RectangleGateDimension gdX = ((RectangleGate) rect).getDimension(xCol);
                    double xLow, xHigh, yLow, yHigh;
                    xLow = !Numbers.isFinite(gdX.getMin()) ? Integer.MIN_VALUE : getXPixel(gdX.getMin());
                    xHigh = !Numbers.isFinite(gdX.getMax()) ? Integer.MAX_VALUE : getXPixel(gdX.getMax());
                    if (isHistogram()) {
                        if (xLow <= tempX && xHigh >= tempX) {
                            toRemove = (RectangleGate) rect;
                            break;
                        }
                    } else {
                        RectangleGateDimension gdY = ((RectangleGate) rect).getDimension(yCol);
                        yHigh = gdY == null ? plotYmin : (!Numbers.isFinite(gdY.getMin()) ? Integer.MIN_VALUE : getYPixel(gdY.getMin()));
                        yLow = gdY == null ? plotYmax : (!Numbers.isFinite(gdY.getMax()) ? Integer.MAX_VALUE : getYPixel(gdY.getMax()));
                        if (xLow <= tempX && xHigh >= tempX && yLow <= tempY && yHigh >= tempY) {
                            toRemove = (RectangleGate) rect;
                            break;
                        }
                    }
                }
            }

        }
        if (toRemove != null) {
            selectedGates.remove(toRemove);
        }
        // TODO remove gate in gating data struct
        // TODO confirm gate removal
        super.mouseClicked(e);
        paintAgain();
    }

    private void rightMouseClickedPoly(MouseEvent e) {
        int tempX = e.getX();
        int tempY = e.getY();
        PolygonGate toRemove = null;

        if (tempPoly.size() > 0) {
            tempPoly.clear();
            highlightPoly = null;
            return;
        }

        ArrayList<Gate> gating = fcp.getGatingForCurrentPlot();
        double tempValX = getXValueFromXPixel(tempX);
        double tempValY = getYValueFromYPixel(tempY);
        for (int i = gating.size() - 1; i >= 0; i--) {
            if (!(gating.get(i) instanceof PolygonGate))
                continue;
            if (((PolygonGate) gating.get(i)).getPath().contains(tempValX, tempValY)) {
                toRemove = (PolygonGate) gating.get(i);
                break;
            }
        }

        if (toRemove != null) {
            // TODO remove gate in gating data struct
            // TODO confirm gate removal
            // polys.remove(toRemove);
            selectedGates.remove(toRemove);
        }
        super.mouseClicked(e);
        paintAgain();
    }

    boolean drag = false;

    private void setForceGatesChanged() {
        this.forceGatesChanged = true;
    }

    @Override
    public void highlightPoints() {
        return;
    }

    @Override
    public void assignAxisLabels() {
        setXAxisLabel(fcp.getXDataName());
        setYAxisLabel(fcp.getYDataName());
    }

}
