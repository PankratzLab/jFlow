package cnv.plots;

import java.io.Serializable;
import common.*;

public class PlotPointSet implements Serializable {
	public static final long serialVersionUID = 1L;
	private PlotPoint[] points;

	public PlotPointSet(PlotPoint[] points) {
		this.points = points;
	}
	
	public PlotPoint[] getPlotPoints() {
		return points;
	}
	
	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static PlotPointSet load(String filename, boolean jar) {
		return (PlotPointSet)Files.readSerial(filename, jar, true);
	}
}
