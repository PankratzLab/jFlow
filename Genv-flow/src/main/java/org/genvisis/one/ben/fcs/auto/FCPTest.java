package org.genvisis.one.ben.fcs.auto;

import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.AbstractPanel2.PLOT_TYPE;
import org.genvisis.one.ben.fcs.FCSDataLoader;
import org.genvisis.one.ben.fcs.FCSPlot;
import org.genvisis.one.ben.fcs.gating.Gate;

public class FCPTest {


	static String DIR = "F:/Flow/concordance/";
	static String FCS = DIR + "2016-08-01_PANEL 1_DHS_Group one_F1636850_009.fcs";
	static String WSP = DIR + "801_Panel1_DHS.wsp";
	static String CONC = DIR
											 + "2016-08-01_PANEL 1_DHS_Group one_F1636850_009.fcs_panel1_gate_def.txt.gz";

	static String[] GATE_NAMES = {
																"Lymphocytes (SSC-A v FSC-A)",
																"Single Cells (FSC-H v FSC-W)",
																"Live cells (PE-)",
																"Tcells (CD3+ CD19-)",
	};

	public static void main(String[] args) {

		FCSPlot fcp = new FCSPlot();
		fcp.loadFile(FCS, true);
		FCSDataLoader loader = fcp.getDataLoader(FCS);

		fcp.loadWorkspaceFile(WSP);

		fcp.setSDVisible(false, false);
		fcp.setMedianVisible(false, false);
		fcp.setSDVisible(false, true);
		fcp.setMedianVisible(false, true);
		fcp.getPanel().setSize(800, 600);
		fcp.getPanel().setLayersInBase(new byte[] {0, 1, 99});
		fcp.setPlotType(PLOT_TYPE.DOT_PLOT);

		fcp.loadAutoValues(fcp.discoverFNumFile(DIR));

		for (String s : GATE_NAMES) {
			Gate g = fcp.getGatingStrategy().gateMap.get(s);
			fcp.setXDataName(g.getXDimension().getParam());
			fcp.setYDataName(g.getYDimension().getParam());
			fcp.setClassifierGate(g.getName());
			fcp.getPanel().setForceGatesChanged();
			fcp.updateGUI();
			fcp.getPanel().createImage();
			String name = g.getName();
			try {
				Thread.sleep(200);
			} catch (InterruptedException e) {
			}
			fcp.screencap(DIR + ext.replaceWithLinuxSafeCharacters(name + ".png"));
			fcp.gateSelected(g, true);
		}
	}


}
