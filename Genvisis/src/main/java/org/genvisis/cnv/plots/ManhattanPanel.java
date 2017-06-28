package org.genvisis.cnv.plots;

import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.SwingUtilities;

import org.genvisis.cnv.plots.ManhattanPlot.ManhattanDataPoint;
import org.genvisis.common.Grafik;


public class ManhattanPanel extends AbstractPanel {

	ManhattanPlot mp;
	byte size = 2;
	byte layer = 1;
	int numPointColors = 3;

	int[] lineValuesToDraw = new int[] {5, 7};
	int[] sizeMult = new int[] {2, 4};
	int[] lineColors = new int[] {3, 3};
	int[] aboveLineColors = new int[] {4, 5};

	public ManhattanPanel(ManhattanPlot parent) {
		super();
		this.mp = parent;
		setZoomable(true, true);
		setSymmetricAxes(false);

		setColorScheme(new Color[] {
																Color.GRAY,
																Color.BLACK,
																Color.LIGHT_GRAY,
																Color.RED,
																Color.MAGENTA,
																Color.GREEN
		});
	}

	protected void drawXAxis(Graphics g, double[] plotMinMaxStep, FontMetrics fontMetrics) {
		String str;
		if (displayXAxisScale) {
			int prevEnd = -1;
			for (int i : linearizedChrBnds.keySet()) {
				int[] bnds = linearizedChrBnds.get(i);
				int x = (int) (bnds[0] + .5 * (bnds[1] - bnds[0]));
				if (x >= plotXmin || !truncate) {
					str = i + "";
					int xLoc = getXPixel(x) - fontMetrics.stringWidth(str) / 2;
					int len = TICK_LENGTH;
					if (xLoc <= prevEnd) {
						len -= TICK_LENGTH / 3;
					}
					Grafik.drawThickLine(g, getXPixel(x), getHeight() - canvasSectionMaximumY, getXPixel(x),
															 getHeight() - (canvasSectionMaximumY - len), TICK_THICKNESS,
															 Color.BLACK);
					if (xLoc > prevEnd) {
						g.drawString(str, xLoc,
												 getHeight() - (canvasSectionMaximumY - TICK_LENGTH - 30));
						prevEnd = xLoc + fontMetrics.stringWidth(str) + 1;
					}
				}
			}
		}
		if (displayXAxis) {
			Grafik.drawThickLine(g, canvasSectionMinimumX - (int) Math.ceil(AXIS_THICKNESS / 2.0),
													 getHeight() - canvasSectionMaximumY,
													 canvasSectionMaximumX + (int) Math.ceil(AXIS_THICKNESS / 2.0),
													 getHeight() - canvasSectionMaximumY, AXIS_THICKNESS, Color.BLACK);
		}
		if (xAxisLabel != null && !"".equals(xAxisLabel) && displayXLabel) {
			g.drawString(xAxisLabel,
									 (getWidth() - axisYWidth/* WIDTH_Y_AXIS */) / 2
											 - fontMetrics.stringWidth(xAxisLabel) / 2
											 + axisYWidth/* WIDTH_Y_AXIS */,
									 getHeight() - 20);
		}
	}

	HashMap<Integer, int[]> linearizedChrBnds = new HashMap<>();

	@Override
	public void generatePoints() {
		ArrayList<ManhattanDataPoint> dataPoints = mp.getData();
		if (dataPoints == null) {
			points = new PlotPoint[0];
			lines = new GenericLine[0];
			setNullMessage("No Data");
			return;
		} else if (dataPoints.isEmpty()) {
			points = new PlotPoint[0];
			lines = new GenericLine[0];
			setNullMessage("Data Loading, Please Wait...");
			return;
		} else {
			setNullMessage(null);
		}
		setForcePlotXmin(dataPoints.get(0).linearLoc - ManhattanPlot.LIN_CHR_BUFFER);
		setForcePlotXmax(dataPoints.get(dataPoints.size() - 1).linearLoc + ManhattanPlot.LIN_CHR_BUFFER);
		points = new PlotPoint[dataPoints.size()];
		int index = 0;
		for (ManhattanDataPoint mdp : dataPoints) {
			int[] bnds = linearizedChrBnds.get(mdp.chr);
			if (bnds == null) {
				bnds = new int[] {Integer.MAX_VALUE, Integer.MIN_VALUE};
				linearizedChrBnds.put(mdp.chr, bnds);
			}
			if (mdp.linearLoc < bnds[0]) {
				bnds[0] = mdp.linearLoc;
			}
			if (mdp.linearLoc > bnds[1]) {
				bnds[1] = mdp.linearLoc;
			}
			points[index++] = new PlotPoint(
																			mdp.mkr == null ? mdp.chr + ":" + mdp.pos : mdp.mkr,
																			PlotPoint.FILLED_CIRCLE,
																			(float) mdp.linearLoc,
																			(float) mdp.transformedPVal,
																			getSize(mdp.transformedPVal),
																			getPointColor(mdp.chr, mdp.transformedPVal), layer);
		}

		lines = new GenericLine[lineValuesToDraw.length];
		for (int i = 0; i < lineValuesToDraw.length; i++) {
			lines[i] = new GenericLine(Integer.MIN_VALUE, lineValuesToDraw[i], (float) Integer.MAX_VALUE,
																 lineValuesToDraw[i], (byte) 1, (byte) lineColors[i], (byte) 99);
		}
	}

	private byte getPointColor(int chr, double transP) {
		int c = (chr % numPointColors);
		for (int i = 0; i < lineValuesToDraw.length; i++) {
			if (transP < lineValuesToDraw[i]) {
				break;
			}
			c = aboveLineColors[i];
		}
		return (byte) c;
	}

	private byte getSize(double pVal) {
		int sz = size;
		for (int i = 0; i < lineValuesToDraw.length; i++) {
			if (pVal < lineValuesToDraw[i]) {
				break;
			}
			sz = size * sizeMult[i];
		}
		return (byte) sz;
	}

	@Override
	public void highlightPoints() {
		// TODO Auto-generated method stub

	}

	@Override
	public void mouseClicked(MouseEvent e) {
		if (SwingUtilities.isRightMouseButton(e)) {
			if (prox != null && prox.size() > 0) {
				for (int i = 0; i < prox.size(); i++) {
					System.out.println(points[prox.get(i)].getId());
				}
			}
		}
	}

	@Override
	public void assignAxisLabels() {
		xAxisLabel = "Chromosome";
		yAxisLabel = "-Log10(pVal)";
	}

}
