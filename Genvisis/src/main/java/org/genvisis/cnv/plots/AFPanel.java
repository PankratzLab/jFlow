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
import org.genvisis.common.PSF;
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
																PSF.Colors.BLUE_VIOLET,
																PSF.Colors.GREEN,
																PSF.Colors.DODGER_BLUE,
																PSF.Colors.SLATE_BLUE,
																PSF.Colors.GREEN_YELLOW,
																PSF.Colors.ORCHID,
																PSF.Colors.AMBER,
																PSF.Colors.MANGO_TANGO,
																PSF.Colors.NAVY,
																PSF.Colors.CORNFLOWER_BLUE,
																PSF.Colors.DARK_SLATE_BLUE,
																PSF.Colors.SLATE_BLUE,
																PSF.Colors.MEDIUM_SLATE_BLUE
		});
		points = new PlotPoint[0];
		resetMessage();

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

	public void resetMessage() {
		setNullMessage("No Data to Display.");
	}

	@Override
	public void highlightPoints() {

	}

	@Override
	public void generatePoints() {
		if (!plot.isForceRedraw() && (plot.getData().isEmpty() || plot.isLoading())) {
			resetMessage();
			points = new PlotPoint[0];
			return;
		}

		lines = plot.isMaskCenter() ? CENTER_LINES : null;

		List<PlotPoint> pointList = new ArrayList<>();
		float afObs, afExp;
		boolean add;
		Marker g1Marker;
		CONFIG config;
		byte color = COLOR;
		Map<Object, String[]> alleles = plot.getObservedAlleles();
		Set<Entry<Object, Double>> dataSet;
		plot.getAlleleInfo().clear();
		dataSet = plot.getData().entrySet();
		for (Entry<Object, Double> obs : dataSet) {
			afObs = obs.getValue().floatValue();
			afExp = plot.getG1KData().get(obs.getKey()).get(plot.getSelectedPop()).floatValue();
			add = !plot.isMaskCenter() ? true : Math.abs(afObs - afExp) > 0.2;
			if (add) {
				g1Marker = plot.getG1KMarkers().get(obs.getKey());
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
				pointList.add(new PlotPoint(obs.getKey().toString(), PointType.FILLED_CIRCLE, afExp,
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
