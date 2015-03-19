package cnv.plots;

import java.util.ArrayList;

import common.Logger;

class HistData {
	double binValue;
	int binCount;
	public HistData(double val, int cnt) {
		this.binValue = val;
		this.binCount = cnt;
	}
}
	
public class HistogramPanel extends AbstractPanel {
	
	private static final double DEFAULT_HALF_BIN_SIZE = 0.001;
	private Logger log;
	private HistogramPlot histPlot;
	
	public HistogramPanel(HistogramPlot histPlot, Logger log) {
		this.histPlot = histPlot;
		this.log = log;
		
	}
	
	@Override
	public void generatePoints() {
		ArrayList<HistData> data = this.histPlot.getData();
		points = new PlotPoint[data.size()];
		for (int i = 0; i < data.size(); i++) {
			points[i] = new PlotPoint("" + data.get(i).binValue, PlotPoint.FILLED_SQUARE, (float) data.get(i).binValue, (float) data.get(i).binCount, (byte) 3, (byte) 0, (byte) 0);
		}
//		points = new PlotPoint[0];
		generateRectangles();
	}
	
	private void generateRectangles() {
		ArrayList<HistData> data = this.histPlot.getData();
		rectangles = new GenericRectangle[data.size()];
		float binHalf = (float) (data.size() > 1 ? ((data.get(1).binValue - data.get(0).binValue) / 2) : DEFAULT_HALF_BIN_SIZE);
		for (int i = 0; i < data.size(); i++) {
			HistData currData = data.get(i);
			rectangles[i] = new GenericRectangle((float)(currData.binValue - binHalf), 0f, (float)(currData.binValue + binHalf), (float)currData.binCount, (byte) 5, true, false, (byte) 0, (byte) 0);
		}
	}

	@Override
	public void highlightPoints() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void assignAxisLabels() {
		displayXaxis = displayYaxis = true;
		xAxisLabel = " ";
		yAxisLabel = " ";
	}

}
