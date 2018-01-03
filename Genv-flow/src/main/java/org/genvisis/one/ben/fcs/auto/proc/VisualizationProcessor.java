package org.genvisis.one.ben.fcs.auto.proc;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;

import javax.swing.SwingConstants;

import org.genvisis.common.Files;
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

	private static String[][] hardcodedAddlImages = {
																									 {"effector helper Tcells (CD95/CD28)",
																										"effector helper Tcells (CCR7- CD45RA+)"},
																									 {"naive helper Tcells (CD95/CD28)",
																										"naive helper Tcells (CCR7+ CD45RA+)"},
																									 {"effector memory helper Tcells (CD95/CD28)",
																										"effector memory helper Tcells (CCR7- CD45RA-)"},
																									 {"central memory helper Tcells (CD95/CD28)",
																										"central memory helper Tcells (CCR7+ CD45RA-)"},
																									 {"effector cytotoxic Tcells (CD95/CD28)",
																										"effector cytotoxic Tcells  (CCR7-  CD45RA+)"},
																									 {"naive cytotoxic Tcells (CD95/CD28)",
																										"naive cytotoxic Tcells (CCR7+ , CD45RA+)"},
																									 {
																										"effector memory cytotoxic Tcells (CD95/CD28)",
																										"effector memory cytotoxic Tcells (CCR7- , CD45RA-)"},
																									 {
																										"central memory cytotoxic Tcells (CD95/CD28)",
																										"central memory cytotoxic Tcells (CCR7+ , CD45RA-)"},
	};

	static class AddlImage {
		String xDim;
		String yDim;
		String parentName;
		String name;

		public AddlImage(String x, String y, String p, String n) {
			this.xDim = x;
			this.yDim = y;
			this.parentName = p;
			this.name = n;
		}
	}

	static final Map<String, List<AddlImage>> addlImgs = new HashMap<>();

	{
		for (String[] addlImg : hardcodedAddlImages) {
			String parent;
			parent = addlImg[1];
			addlImgs.put(parent, new ArrayList<AddlImage>());
			addlImgs.get(parent).add(new AddlImage("Comp-BV 605-A (CD95)", "Comp-BV 510-A (CD28)",
																						 parent, addlImg[0]));
		}

		//
		// String key;
		// key = "effector helper Tcells (CCR7- CD45RA+)";
		// addlImgs.put(key, new ArrayList<AddlImage>());
		// addlImgs.get(key).add(new AddlImage("Comp-BV 605-A (CD95)", "Comp-BV 510-A (CD28)", key,
		// "effector helper Tcells (CD95/CD28)"));
		//
		// key = "naive helper Tcells (CCR7+ CD45RA+)";
		// addlImgs.put(key, new ArrayList<AddlImage>());
		// addlImgs.get(key).add(new AddlImage("Comp-BV 605-A (CD95)", "Comp-BV 510-A (CD28)", key,
		// "naive helper Tcells (CD95/CD28)"));
		//
		// key = "effector memory helper Tcells (CCR7- CD45RA-)";
		// addlImgs.put(key, new ArrayList<AddlImage>());
		// addlImgs.get(key).add(new AddlImage("Comp-BV 605-A (CD95)", "Comp-BV 510-A (CD28)", key,
		// "effector memory helper Tcells (CD95/CD28)"));
		//
		// key = "central memory helper Tcells (CCR7+ CD45RA-)";
		// addlImgs.put(key, new ArrayList<AddlImage>());
		// addlImgs.get(key).add(new AddlImage("Comp-BV 605-A (CD95)", "Comp-BV 510-A (CD28)", key,
		// "central memory helper Tcells (CD95/CD28)"));
		//
		// key = "effector cytotoxic Tcells (CCR7- CD45RA+)";
		// addlImgs.put(key, new ArrayList<AddlImage>());
		// addlImgs.get(key).add(new AddlImage("Comp-BV 605-A (CD95)", "Comp-BV 510-A (CD28)", key,
		// "effector cytotoxic Tcells (CD95/CD28)"));
		//
		// key = "naive cytotoxic Tcells (CCR7+ , CD45RA+)";
		// addlImgs.put(key, new ArrayList<AddlImage>());
		// addlImgs.get(key).add(new AddlImage("Comp-BV 605-A (CD95)", "Comp-BV 510-A (CD28)", key,
		// "naive cytotoxic Tcells (CD95/CD28)"));
		//
		// key = "effector memory cytotoxic Tcells (CCR7- , CD45RA-)";
		// addlImgs.put(key, new ArrayList<AddlImage>());
		// addlImgs.get(key).add(new AddlImage("Comp-BV 605-A (CD95)", "Comp-BV 510-A (CD28)", key,
		// "effector memory cytotoxic Tcells (CD95/CD28)"));
		//
		// key = "central memory cytotoxic Tcells (CCR7+ , CD45RA-)";
		// addlImgs.put(key, new ArrayList<AddlImage>());
		// addlImgs.get(key).add(new AddlImage("Comp-BV 605-A (CD95)", "Comp-BV 510-A (CD28)", key,
		// "central memory cytotoxic Tcells (CD95/CD28)"));
	}

	@Override
	public void processSample(SampleNode sn, Logger log) throws IOException {
		System.gc();
		final FCSPlot fcp = new FCSPlot();
		fcp.getPanel().setColorScheme(new Color[] {Color.BLACK, Color.RED,
																							 new Color(128, 128, 128, 64)});
		fcp.getPanel().allowSkip = false;

		long time1 = System.nanoTime();
		fcp.loadWorkspaceFile(sn.wspFile);
		Set<String> sampeIds = fcp.getWorkbench().getAllSamples();
		String id = sampeIds.toArray(new String[1])[0];
		if (fcp.getWorkbench().getSample(id).fcsFile.equals(sn.fcsFile)) {
			fcp.setCurrentSampleInWSP(sn.fcsFile);
		} else if (fcp.getWorkbench().getSample(id).fcsFile.equals(ext.rootOf(sn.fcsFile, true))) {
			fcp.setCurrentSampleInWSP(ext.rootOf(sn.fcsFile, true));
		}
		fcp.setCurrentSampleInWSP(id);

		boolean hasAll = true;

		for (String s : fcp.getGatingStrategy().getAllGateNames()) {
			Gate g = fcp.getGatingStrategy().gateMap.get(s);

			String cleanedName = ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(sn.fcsFile));
			String outFile = outDir + cleanedName + "/" + cleanedName + "."
											 + ext.replaceWithLinuxSafeCharacters(g.getName());


			if (Files.exists(outFile + ".png")) {
				if (addlImgs.containsKey(g.getName())) {
					boolean all = true;
					for (AddlImage addl : addlImgs.get(g.getName())) {
						String outFile2 = outDir + cleanedName + "/" + cleanedName + "."
															+ ext.replaceWithLinuxSafeCharacters(addl.name);
						if (!Files.exists(outFile2 + ".png")) {
							all = false;
							break;
						}
					}
					if (all) {
						continue;
					}
				} else {
					continue;
				}
			}

			hasAll = false;
			break;
		}
		if (hasAll) {
			log.report("All screenshots found for " + fcp.getGatingStrategy().getAllGateNames().size()
								 + " gates; Skipping FCS file: " + sn.fcsFile);
			return;
		}

		long time2 = System.nanoTime();
		fcp.loadFile(sn.fcsFile, true);
		FCSDataLoader loader = fcp.getDataLoader(sn.fcsFile);
		int rowCnt = loader.getCount();

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

		long time4 = System.nanoTime();
		ArrayList<String> gateNames = new ArrayList<>();
		ArrayList<long[]> gateTimes = new ArrayList<>();

		for (String s : fcp.getGatingStrategy().getAllGateNames()) {
			long[] times = new long[4];

			Gate g = fcp.getGatingStrategy().gateMap.get(s);

			String cleanedName = ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(sn.fcsFile));
			String outFile = outDir + cleanedName + "/" + cleanedName + "."
											 + ext.replaceWithLinuxSafeCharacters(g.getName());

			if (Files.exists(outFile + ".png")) {
				if (addlImgs.containsKey(g.getName())) {
					boolean all = true;
					for (AddlImage addl : addlImgs.get(g.getName())) {
						String outFile2 = outDir + cleanedName + "/" + cleanedName + "."
															+ ext.replaceWithLinuxSafeCharacters(addl.name);
						if (!Files.exists(outFile2 + ".png")) {
							all = false;
							break;
						}
					}
					if (all) {
						continue;
					}
				} else {
					continue;
				}
			}

			gateNames.add(s);

			fcp.gateSelected(g.getParentGate(), false);
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

			if (addlImgs.containsKey(g.getName())) {
				fcp.gateSelected(g, false);
				for (AddlImage addl : addlImgs.get(g.getName())) {
					String outFile2 = outDir + cleanedName + "/" + cleanedName + "."
														+ ext.replaceWithLinuxSafeCharacters(addl.name);
					if (Files.exists(outFile2 + ".png")) {
						continue;
					}
					fcp.setXDataName(addl.xDim);
					fcp.setYDataName(addl.yDim);
					fcp.getPanel().setTitle(addl.name);
					fcp.getPanel().setForceGatesChanged();
					fcp.getPanel().createImage();
					fcp.screencap(outFile2 + ".png");
				}
			}

		}

		loader.emptyAndReset();
		// fcp = null;
		System.gc();

		long time5 = System.nanoTime();

		StringBuilder sb1 = new StringBuilder("TIMING-HDR").append("\t");
		sb1.append("FCS_FILE").append("\t");
		sb1.append("ROWS").append("\t");
		sb1.append("INIT").append("\t");
		sb1.append("WSP_LOAD").append("\t");
		sb1.append("LOAD").append("\t");
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
