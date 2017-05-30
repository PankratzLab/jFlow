package org.genvisis.one.ben.fcs.auto.proc;

import java.io.IOException;

import org.genvisis.common.Logger;
import org.genvisis.one.ben.fcs.AbstractPanel2.PLOT_TYPE;
import org.genvisis.one.ben.fcs.FCSDataLoader;
import org.genvisis.one.ben.fcs.FCSPlot;
import org.genvisis.one.ben.fcs.auto.FCSProcessingPipeline;
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
		FCSPlot fcp = new FCSPlot();

		fcp.loadFile(sn.fcsFile, true);
		FCSDataLoader loader = fcp.getDataLoader(sn.fcsFile);
		loader.waitForData();

		fcp.loadWorkspaceFile(sn.wspFile);

		String fNum = fcp.discoverFNumFile(autoDir);
		if (fNum == null)
			return;
		fcp.loadAutoValues(fNum);

		fcp.setSize(1000, 800);
		fcp.getPanel().setSize(800, 600);
		fcp.setSDVisible(false, false);
		fcp.setSDVisible(false, true);
		fcp.setMedianVisible(false, false);
		fcp.setMedianVisible(false, true);
		fcp.getPanel().setLayersInBase(new byte[] {0, 1, 99});
		fcp.setPlotType(PLOT_TYPE.DOT_PLOT);

		for (String s : FCSProcessingPipeline.GATE_NAMES) {
			Gate g = fcp.getGatingStrategy().gateMap.get(s);
			fcp.setXDataName(g.getXDimension().getParam());
			fcp.setYDataName(g.getYDimension().getParam());

			fcp.setClassifierGate(g.getName());

			// String name = g.getName();
			// String fNumD = FCSProcessingPipeline.getFNum(sn.fcsFile);
			String outFile = outDir + sn.fcsFile + "/"
											 + sn.fcsFile + "." + g.getXDimension().getParam() + "v"
											 + g.getYDimension().getParam();
			fcp.getPanel().classifierPrev = false;
			fcp.getPanel().setForceGatesChanged();
			fcp.getPanel().createImage();
			fcp.screencap(outFile + ".png");

			fcp.getPanel().classifierPrev = true;
			fcp.getPanel().setForceGatesChanged();
			fcp.getPanel().createImage();
			fcp.screencap(outFile + "_prev.png");

			fcp.gateSelected(g, true);
		}

		loader.emptyAndReset();
		fcp = null;
		System.gc();
	}

}
