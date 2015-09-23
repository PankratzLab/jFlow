package cnv.plots;

import common.*;

import java.awt.event.*;

public class QQPanel extends AbstractPanel implements ComponentListener {
	public static final long serialVersionUID = 1L;

	private double[][] pvals;
	private boolean log10;
	private boolean rotated;
	private float maxValue;
	
	public QQPanel(double[][] pvals, boolean log10, boolean rotated, float maxValue) {
		super();
		
		this.pvals = pvals;
		this.log10 = log10;
		this.rotated = rotated;
		this.maxValue = maxValue;
		
		createLookup(false);

		// taken care of in AbstractPanel constructor
		addComponentListener(this);
		setZoomable(true, true);
	}

	public void assignAxisLabels() {
		if (log10) {
			if (rotated) {
				xAxisLabel = "-log10(rank/n)";
				yAxisLabel = "-log10(p-value) - -log10(rank/n)";
			} else {
				xAxisLabel = "-log10(rank/n)";
				yAxisLabel = "-log10(p-value)";
			}
			plotXmin = 0;
			plotYmin = 0;
		} else {
			xAxisLabel = "Expected quantiles";
			yAxisLabel = "Observed quantiles";
		}
	}
	
	public boolean invertX() {
		return false;
	}

	public boolean invertY() {
		return false;
	}
	
	public void highlightPoints() {}

	public void generatePoints() {
		int[] keys;
		int count;
		int max;
		
		lines = new GenericLine[1];
		max = 0;
		for (int i = 0; i<pvals.length; i++) {
			max = Math.max(max, pvals[i].length);
        }
		max = (int)Math.ceil((-1*Math.log10(1.0/(double)max)));
		if (rotated) {
			lines[0] = new GenericLine(0, 0, max, 0, (byte)2, (byte)1, (byte)0);
		} else if (log10) { 
//			lines[0] = new PlotLine(0, 0, (float)(-1*Math.log10((1.0/pvals.length))), (float)(-1*Math.log10((1.0/pvals.length))), (byte)2, (byte)1);
			lines[0] = new GenericLine(0, 0, max, max, (byte)2, (byte)1, (byte)0);
		} else {
			lines[0] = new GenericLine(0, 0, 1, 1, (byte)2, (byte)1, (byte)0);
		}

		count = 0;
		for (int i = 0; i<pvals.length; i++) {
			count += pvals[i].length;
        }
		points = new PlotPoint[count];
		
		count = 0;
		for (int i = 0; i<pvals.length; i++) {
			keys = Sort.quicksort(pvals[i]);

			for (int j = 0; j<pvals[i].length; j++) {
				if (rotated) { 
					points[count] = new PlotPoint(keys[j]+"", PlotPoint.FILLED_CIRCLE, (float)(-1*Math.log10(((double)keys[j]+1)/pvals[i].length)), Math.min(maxValue, (float)(-1*Math.log10(pvals[i][j]))-(float)(-1*Math.log10(((double)keys[j]+1)/pvals[i].length))), (byte)6, (byte)(pvals.length==1?0:i+2), (byte)0);
				} else if (log10) { 
					points[count] = new PlotPoint(keys[j]+"", PlotPoint.FILLED_CIRCLE, (float)(-1*Math.log10(((double)keys[j]+1)/pvals[i].length)), Math.min(maxValue, (float)(-1*Math.log10(pvals[i][j]))), (byte)6, (byte)(pvals.length==1?0:i+2), (byte)0);
				} else {
					points[count] = new PlotPoint(keys[j]+"", PlotPoint.FILLED_CIRCLE, (float)(((double)keys[j]+1)/pvals[i].length), (float)pvals[i][j], (byte)6, (byte)(pvals.length==1?0:i+2), (byte)0);
				}
				count++;
			}
		}
	}

	public void componentHidden(ComponentEvent e) {}

	public void componentMoved(ComponentEvent e) {}

	public void componentResized(ComponentEvent e) {
		paintAgain();
	}

	public void componentShown(ComponentEvent e) {}

	public static void main(String[] args) {
		QQPlot.main(new String[] {});
	}
}