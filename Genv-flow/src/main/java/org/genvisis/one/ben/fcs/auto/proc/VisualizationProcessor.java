package org.genvisis.one.ben.fcs.auto.proc;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

import javax.swing.SwingConstants;

import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.AbstractPanel2.PLOT_TYPE;
import org.genvisis.one.ben.fcs.FCSDataLoader;
import org.genvisis.one.ben.fcs.FCSPlot;
import org.genvisis.one.ben.fcs.gating.Gate;
import org.genvisis.one.ben.fcs.gating.Workbench.SampleNode;

public class VisualizationProcessor implements SampleProcessor {

	final String autoDir;
	final String outDir;

	public VisualizationProcessor(String a, String o) {
		this.autoDir = a;
		this.outDir = o;
	}

	@Override
	public void processSample(SampleNode sn, Logger log) throws IOException {
		System.gc();
		final FCSPlot fcp = new FCSPlot();
		fcp.getPanel().setColorScheme(new Color[] {Color.BLACK, Color.RED,
																							 new Color(128, 128, 128, 64)});

		long time1 = System.nanoTime();
		fcp.loadFile(sn.fcsFile, true);
		FCSDataLoader loader = fcp.getDataLoader(sn.fcsFile);
		loader.waitForData();

		long time2 = System.nanoTime();
		int rowCnt = loader.getCount();

		fcp.loadWorkspaceFile(sn.wspFile);

		long time3 = System.nanoTime();

		// String fNum = fcp.discoverFNumFile(autoDir);
		// if (fNum == null)
		// return;
		// fcp.loadAutoValues(fNum);

		fcp.setSize(1000, 800);
		fcp.getPanel().setSize(800, 600);
		fcp.setSDVisible(false, false);
		fcp.setSDVisible(false, true);
		fcp.setMedianVisible(false, false);
		fcp.setMedianVisible(false, true);
		fcp.getPanel().setLayersInBase(new byte[] {0, 1, 99});
		fcp.setPlotType(PLOT_TYPE.DOT_PLOT);

		// for (String s : FCSProcessingPipeline.GATE_NAMES) {
		long time4 = System.nanoTime();
		ArrayList<String> gateNames = new ArrayList<>();
		ArrayList<long[]> gateTimes = new ArrayList<>();

		for (String s : fcp.getGatingStrategy().getAllGateNames()) {
			gateNames.add(s);
			long[] times = new long[4];

			Gate g = fcp.getGatingStrategy().gateMap.get(s);
			if (g.getParentGate() != null) {
				fcp.gateSelected(g.getParentGate(), false);
			}
			times[0] = System.nanoTime();

			if (g.getXDimension() == null && g.getYDimension() != null) {
				// correct for swapped histogram
				g.setXDimension(g.getYDimension());
				g.setYDimension(null);
			}
			fcp.setXDataName(g.getXDimension().getParam());
			if (g.getYDimension() != null) {
				fcp.setYDataName(g.getYDimension().getParam());
				fcp.setPlotType(PLOT_TYPE.HEATMAP);
			} else {
				fcp.setPlotType(PLOT_TYPE.HISTOGRAM);
			}
			g.setColor(1);
			g.setFillGate(true);
			g.setDisplayName(false);
			fcp.getPanel().setTitle(g.getName());
			fcp.getPanel().setTitleLocation(SwingConstants.NORTH);
			fcp.getPanel().setDisplayTitle(true);
			fcp.getPanel().setTitleFontSize(20f);

			// fcp.setClassifierGate(g.getName());

			String cleanedName = ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(sn.fcsFile));
			String outFile = outDir + cleanedName + "/" + cleanedName + "."
											 + ext.replaceWithLinuxSafeCharacters(g.getName());
			fcp.getPanel().classifierPrev = false;
			fcp.getPanel().setForceGatesChanged();
			times[1] = System.nanoTime();

			fcp.getPanel().createImage();
			times[2] = System.nanoTime();

			fcp.screencap(outFile + ".png");

			times[3] = System.nanoTime();

			g.setFillGate(false);

			// fcp.getPanel().classifierPrev = true;
			// fcp.getPanel().setForceGatesChanged();
			// fcp.getPanel().createImage();
			// fcp.screencap(outFile + "_prev.png");
			gateTimes.add(times);
		}

		loader.emptyAndReset();
		// fcp = null;
		System.gc();

		long time5 = System.nanoTime();

		StringBuilder sb1 = new StringBuilder("TIMING-HDR").append("\t");
		sb1.append("FCS_FILE").append("\t");
		sb1.append("ROWS").append("\t");
		sb1.append("INIT").append("\t");
		sb1.append("LOAD").append("\t");
		sb1.append("WSP_LOAD").append("\t");
		sb1.append("CONF").append("\t");
		for (String s : gateNames) {
			sb1.append(s).append("\t");
			sb1.append("CONF").append("\t");
			sb1.append("CREATE").append("\t");
			sb1.append("SCREEN").append("\t");
		}
		sb1.append("CLEANUP");

		StringBuilder sb = new StringBuilder("TIMING").append("\t");
		sb.append(sn.fcsFile).append("\t");
		sb.append(rowCnt).append("\t");
		sb.append(TimeUnit.NANOSECONDS.toSeconds(time1)).append("\t");
		sb.append(TimeUnit.NANOSECONDS.toSeconds(time2)).append("\t");
		sb.append(TimeUnit.NANOSECONDS.toSeconds(time3)).append("\t");
		sb.append(TimeUnit.NANOSECONDS.toSeconds(time4)).append("\t");
		for (long[] timing : gateTimes) {
			for (long l : timing) {
				sb.append(TimeUnit.NANOSECONDS.toSeconds(l)).append("\t");
			}
		}
		sb.append(TimeUnit.NANOSECONDS.toSeconds(time5));

		System.out.println(sb1.toString());
		System.out.println(sb.toString());
	}
}
