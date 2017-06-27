package org.genvisis.cnv.plots;

import java.awt.Color;
import java.awt.event.MouseEvent;
import java.util.ArrayList;

import javax.swing.SwingUtilities;

import org.genvisis.cnv.plots.ManhattanPlot.ManhattanDataPoint;


public class ManhattanPanel extends AbstractPanel {

	ManhattanPlot mp;
	byte size = 2;
	byte layer = 1;
	int numColors = 2;

	int[] lineValuesToDraw = new int[] {5, 7};
	int[] sizeMult = new int[] {3, 5};

	public ManhattanPanel(ManhattanPlot parent) {
		super();
		this.mp = parent;
		setZoomable(true, true);
		setSymmetricAxes(false);

		displayXAxisScale = false;
		displayXLabel = false;

		setColorScheme(new Color[] {
																Color.BLACK,
																Color.GRAY,
																Color.RED
		});
	}

	@Override
	public void generatePoints() {
		ArrayList<ManhattanDataPoint> dataPoints = mp.getData();
		if (dataPoints == null || dataPoints.isEmpty()) {
			points = new PlotPoint[0];
			errorMessage = "No Data";
			return;
		} else {
			errorMessage = null;
		}
		setForcePlotXmin(dataPoints.get(0).linearLoc - ManhattanPlot.LIN_CHR_BUFFER);
		setForcePlotXmax(dataPoints.get(dataPoints.size() - 1).linearLoc + ManhattanPlot.LIN_CHR_BUFFER);
		points = new PlotPoint[dataPoints.size()];
		int index = 0;
		for (ManhattanDataPoint mdp : dataPoints) {
			points[index++] = new PlotPoint(
																			mdp.mkr == null ? mdp.chr + ":" + mdp.pos : mdp.mkr,
																			PlotPoint.FILLED_CIRCLE,
																			(float) mdp.linearLoc,
																			(float) mdp.transformedPVal,
																			getSize(mdp.transformedPVal),
																			(byte) (mdp.chr % numColors), layer);
		}

		lines = new GenericLine[lineValuesToDraw.length];
		for (int i = 0; i < lineValuesToDraw.length; i++) {
			lines[i] = new GenericLine(Integer.MIN_VALUE, lineValuesToDraw[i], (float) Integer.MAX_VALUE,
																 lineValuesToDraw[i], (byte) 1, (byte) 2, (byte) 99);
		}
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
		xAxisLabel = "";
		yAxisLabel = "-Log10(pVal)";
	}

}
