package cnv.plots;

import java.awt.event.*;

public class GenericPanel extends AbstractPanel implements ComponentListener {
	public static final long serialVersionUID = 1L;

	private float[][][] coords;
	
	public GenericPanel(float[][][] coords, String xLabel, String yLabel) {
		super();
		
		this.coords = coords;
		this.xAxisLabel = xLabel;
		this.yAxisLabel = yLabel;
		
		addComponentListener(this);
		setZoomable(true, true);
	}

	public void assignAxisLabels() {}
	
	public boolean invertX() {
		return false;
	}

	public boolean invertY() {
		return false;
	}
	
	public void highlightPoints() {}
	
	public void generatePoints() {
		int count;
		
		count = 0;
		for (int i = 0; i<coords.length; i++) {
			count += coords[i].length;
        }
		points = new PlotPoint[count];

		count = 0;
		for (int i = 0; i<coords.length; i++) {
			for (int j = 0; j<coords[i].length; j++) {
				points[count] = new PlotPoint(count+"", PlotPoint.FILLED_CIRCLE, coords[i][j][0], coords[i][j][1], (byte)6, (byte)(coords.length==1?0:i+2), (byte)0);
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
		GenericPlot.main(new String[] {});
	}
}