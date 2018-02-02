package org.genvisis.cnv.plots;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.plots.PlotPoint.PointType;
import org.genvisis.seq.manage.StrandOps;
import org.genvisis.seq.manage.StrandOps.CONFIG;

/**
 * {@link AbstractPanel} subtype that creates {@link PlotPoint}s from the data loaded by
 * {@link AFPlot}.
 * 
 */
public class AFPanel extends AbstractPanel {

	public static final long serialVersionUID = 1L;

	/**
	 * Point size
	 */
	private static final byte PTSIZE = 1;
	/**
	 * Point color
	 */
	private static final byte COLOR = 0;
	/**
	 * Layer of points
	 */
	private static final byte LAYER = 1;

	AFPlot plot;
	Map<CONFIG, Integer> configColorMap;

	private static final GenericLine[] CENTER_LINES = {
																										 new GenericLine(0f, 0.2f, 0.8f, 1f, (byte) 2,
																																		 (byte) 2, LAYER),
																										 new GenericLine(0.2f, 0f, 1f, 0.8f, (byte) 2,
																																		 (byte) 2, LAYER),
																										 new GenericLine(0f, 0f, 1f, 1f, (byte) 2,
																																		 (byte) 3, LAYER)
	};

	public AFPanel(AFPlot plot) {
		this.plot = plot;
		setColorScheme(new Color[] {
																Color.BLACK,
																Color.GRAY,
																Color.YELLOW,
																Color.RED,
																new Color(140, 20, 180), // deep purple
																new Color(33, 87, 0), // dark green
																new Color(55, 129, 252), // light blue
																new Color(94, 88, 214), // light purple
																new Color(189, 243, 61), // light green
																new Color(217, 109, 194), // pink
																new Color(0, 0, 128), // ALL KINDS OF BLUES
																new Color(100, 149, 237),
																new Color(72, 61, 139),
																new Color(106, 90, 205),
																new Color(123, 104, 238),
		});
		points = new PlotPoint[0];
		setNullMessage("No Data to Display.");

		setForcePlotXmin(0f);
		setForcePlotXmax(1f);
		setForcePlotYmin(0f);
		setForcePlotYmax(1f);
		setSymmetricAxes(true);
		setForceXAxisWholeNumbers(true);
		setForceYAxisWholeNumbers(true);

		configColorMap = new HashMap<>();
		for (int i = 0; i < CONFIG.values().length; i++) {
			configColorMap.put(CONFIG.values()[i], i);
		}
	}

	@Override
	public void highlightPoints() {

	}

	@Override
	public void generatePoints() {
		if (!plot.isForceRedraw() && (plot.getObservedData() == null || plot.isLoading()))
			return;
		plot.setForceRedraw(false);

		lines = plot.isMaskCenter() ? CENTER_LINES : null;

		List<PlotPoint> pointList = new ArrayList<>();
		float afObs, afExp;
		boolean add;
		String g1MkrNm;
		Marker g1Marker;
		CONFIG config;
		byte color = COLOR;
		Map<String, String[]> alleles = plot.getAlleleMap();
		Set<Entry<String, Double>> dataSet;
		plot.getAlleleInfo().clear();
		dataSet = plot.getObservedData().entrySet();
		for (Entry<String, Double> obs : dataSet) {
			g1MkrNm = plot.isChrPosLookup() ? plot.getG1KChrPosLookup().get(obs.getKey()) : obs.getKey();
			afObs = obs.getValue().floatValue();
			afExp = plot.getG1KData().get(g1MkrNm).get(plot.getSelectedPop()).floatValue();
			add = !plot.isMaskCenter() ? true : Math.abs(afObs - afExp) > 0.2;
			if (add) {
				g1Marker = plot.getG1KMarkers().get(g1MkrNm);
				color = COLOR;
				if (alleles.containsKey(obs.getKey()) && plot.isColorByConfig()) {
					config = StrandOps.determineStrandConfig(alleles.get(obs.getKey()),
																									 new String[] {g1Marker.getRef()
																																				 .getBaseString(),
																																 g1Marker.getAlt()
																																				 .getBaseString()});
					plot.getAlleleInfo().add(config);
					color = (byte) configColorMap.get(config).intValue();
				}
				pointList.add(new PlotPoint(obs.getKey() + " | " + g1MkrNm, PointType.FILLED_CIRCLE, afExp,
																		afObs, PTSIZE, color,
																		LAYER));
			}
		}
		points = pointList.toArray(new PlotPoint[pointList.size()]);
		plot.updateAlleleInfoPanel();
	}

	@Override
	public void assignAxisLabels() {
		displayXAxis = displayYAxis = true;
		yAxisLabel = "Observed";
		xAxisLabel = "Expected";
	}

}
