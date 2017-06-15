package org.genvisis.one.ben.fcs.auto.proc;

import java.awt.Color;
import java.io.IOException;

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

		fcp.loadFile(sn.fcsFile, true);
		FCSDataLoader loader = fcp.getDataLoader(sn.fcsFile);
		loader.waitForData();

		fcp.loadWorkspaceFile(sn.wspFile);

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
		for (String s : fcp.getGatingStrategy().getAllGateNames()) {
			Gate g = fcp.getGatingStrategy().gateMap.get(s);
			if (g.getParentGate() != null) {
				fcp.gateSelected(g.getParentGate(), false);
			}

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
			fcp.getPanel().createImage();

			fcp.screencap(outFile + ".png");

			g.setFillGate(false);

			// fcp.getPanel().classifierPrev = true;
			// fcp.getPanel().setForceGatesChanged();
			// fcp.getPanel().createImage();
			// fcp.screencap(outFile + "_prev.png");

		}

		loader.emptyAndReset();
		// fcp = null;
		System.gc();
	}
}
