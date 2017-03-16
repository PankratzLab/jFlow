package org.genvisis.one.ben.fcs.auto;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.HashSet;
import java.util.concurrent.ConcurrentHashMap;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.AbstractPanel2.PLOT_TYPE;
import org.genvisis.one.ben.fcs.FCSDataLoader;
import org.genvisis.one.ben.fcs.FCSFileDuplicator;
import org.genvisis.one.ben.fcs.FCSPlot;
import org.genvisis.one.ben.fcs.gating.Gate;
import org.genvisis.one.ben.fcs.gating.Workbench.SampleNode;

public class GateExportTest {

	private static final String[][] MATCH = {
                                    {"lymph", "Lymphocytes (SSC-A v FSC-A)"},
                                    {"Singlets", "Single Cells (FSC-H v FSC-W)"},
                                    {"PE.A","Live cells (PE-)"},
                                    {"CD3.","Tcells (CD3+ CD19-)"},                                   
		};

	static String[] GATE_NAMES = {
    	"Lymphocytes (SSC-A v FSC-A)",
			"Single Cells (FSC-H v FSC-W)",
			"Live cells (PE-)",
			"Tcells (CD3+ CD19-)",        
	};
	
	public GateExportTest(String fcs, String wsp, String auto, String out) {
		fcsDir = fcs;
		wspDir = wsp;
		autoDir = auto;
		outDir = out;
	}
	
	private String fcsDir, wspDir, autoDir, outDir;
	
	private void run() throws IOException {

		ProcessorFactory<? extends SampleProcessor> pf;

		// pf = new GateAssignmentFactory();
		// pf = new PercentageWriterFactory();
		// pf = new LeafDataSamplerFactory("/scratch.global/cole0482/FCS/", new Logger());
//		pf = new ProcessorFactory<ConcordanceProcessor>() {
//
//			ConcurrentHashMap<String, ConcurrentHashMap<String, String>> resultMap = new ConcurrentHashMap<>();
//			ConcurrentHashMap<String, ConcurrentHashMap<String, String>> resultMapTree = new ConcurrentHashMap<>();
//			public void cleanup(Object owner) {
//				PrintWriter writer = Files.getAppropriateWriter(OUT + "concordance.xln");
//				
//				Set<String> files = resultMap.keySet();
//				for (String s : files) {
//					for (int i = 0; i < MATCH.length; i++) {
//						writer.println(s + "\t" + MATCH[i][1] + "\t" + resultMap.get(s).get(MATCH[i][1])); 
//					}
//				}
//				writer.flush();
//				writer.close();
//
//				writer = Files.getAppropriateWriter(OUT + "concordance_tree.xln");
//				
//				files = resultMapTree.keySet();
//				for (String s : files) {
//					for (int i = 0; i < MATCH.length; i++) {
//						writer.println(s + "\t" + MATCH[i][1] + "\t" + resultMapTree.get(s).get(MATCH[i][1])); 
//					}
//				}
//				writer.flush();
//				writer.close();
//			};
//
//			@Override
//			public ConcordanceProcessor createProcessor(Object owner, int index) {
//				return new ConcordanceProcessor(resultMap, resultMapTree);
//			}
//		};
//
//		SamplingPipeline sp = new SamplingPipeline(1, null, WSP, FCS, null, OUT, pf);
		
		pf = new ProcessorFactory<SampleProcessor>() {
			
			@Override
			public void cleanup(Object owner) {}
			@Override
			public SampleProcessor createProcessor(Object owner, int index) {
				return new VisualizationProcessor();
			}
		};
		
		SamplingPipeline sp = new SamplingPipeline(1, null, wspDir, fcsDir, null, outDir, pf);
		
		sp.run();

		while (!sp.checkFinished()) {
			Thread.yield();
		}

	}
	
	private static String getFNum(String file) {
		String[] pts = ext.removeDirectoryInfo(file).split("_");
		String fNum = null;
		for (String p : pts) {
			if (p.startsWith("F")) {
				fNum = p;
				break;
			}
		}
		return fNum;
	}
	
	
	class VisualizationProcessor implements SampleProcessor {
		
		@Override
		public void processSample(SampleNode sn, Logger log) throws IOException {
			FCSPlot fcp = new FCSPlot();
			
			fcp.loadFile(sn.fcsFile, true);
			FCSDataLoader loader = fcp.getDataLoader(sn.fcsFile);
			loader.waitForData();
			
			fcp.loadWorkspaceFile(sn.wspFile);
			
			String fNum = fcp.discoverFNumFile(autoDir);
			if (fNum == null) return;
			fcp.loadAutoValues(fNum);
			
			fcp.setSize(1000, 800);
			fcp.getPanel().setSize(800, 600);
			fcp.setSDVisible(false, false);
			fcp.setSDVisible(false, true);
			fcp.setMedianVisible(false, false);
			fcp.setMedianVisible(false, true);
			fcp.getPanel().setLayersInBase(new byte[]{0, 1, 99});
			fcp.setPlotType(PLOT_TYPE.DOT_PLOT);
			
			for (String s : GATE_NAMES) {
				Gate g = fcp.getGatingStrategy().gateMap.get(s);
				fcp.setXDataName(g.getXDimension().getParam());
				fcp.setYDataName(g.getYDimension().getParam());
				
				fcp.setClassifierGate(g.getName());
				
				String name = g.getName();
				String fNumD = getFNum(sn.fcsFile);

				fcp.getPanel().classifierPrev = false;
				fcp.getPanel().setForceGatesChanged();
				fcp.getPanel().createImage();
				fcp.screencap(outDir + fNumD + "/" + ext.replaceWithLinuxSafeCharacters(fNumD + "_" + name + ".png"));
				
				fcp.getPanel().classifierPrev = true;
				fcp.getPanel().setForceGatesChanged();
				fcp.getPanel().createImage();
				fcp.screencap(outDir + fNumD + "/" + ext.replaceWithLinuxSafeCharacters(fNumD + "_" + name + "_prev.png"));
				
				fcp.gateSelected(g, true);
			}
			
			loader.emptyAndReset();
			fcp = null;
			System.gc();
		}
		
	}
	
	
	class ConcordanceProcessor extends AbstractSampleProcessor {
		private static final String AUTO_FILE_SUFF = "def.txt.gz";

		private String autoDir;
		private String[] filesInAutoDir;
		String[][] autoData;
		private ConcurrentHashMap<String, ConcurrentHashMap<String, String>> results;
		private ConcurrentHashMap<String, ConcurrentHashMap<String, String>> resultsTree;
		
		public ConcordanceProcessor(String autoDir, ConcurrentHashMap<String, ConcurrentHashMap<String, String>> resultMap, ConcurrentHashMap<String, ConcurrentHashMap<String, String>> resultMap2) {
			this.autoDir = autoDir;
			filesInAutoDir = (new File(autoDir)).list(new FilenameFilter() {
				@Override
				public boolean accept(File dir, String name) {
					return name.endsWith(AUTO_FILE_SUFF);
				}
			});
			this.results = resultMap;
			this.resultsTree = resultMap2;
		}
		
		@Override
		public void processSample(SampleNode sn, Logger log) throws IOException {
			loadPopsAndGates(sn);
			loadData(sn);
			
			loadAutoValues(sn);
			if (autoData == null) {
				System.err.println("Error - problem occured with auto-data.  Check log and try again for sample: " + sn.fcsFile);
				return;
			}
			
			boolean[][] hand = new boolean[autoData[0].length][d.getCount()];
			String[] gateNames = new String[hand.length];
			for (int i = 0; i < hand.length; i++) {
				Gate gate = sn.gating.gateMap.get(MATCH[i][1]);
				hand[i] = gate.gate(d);
				gateNames[i] = gate.getName();
			}
			
			boolean[][] auto = new boolean[autoData[0].length][];
			for (int i = 0; i < auto.length; i++) {
				String[] data = Matrix.extractColumn(autoData, i);
				auto[i] = new boolean[data.length];
				for (int k = 0; k < data.length; k++) {
					auto[i][k] = Boolean.valueOf(data[k]); 
				}
			}
			
			ConcurrentHashMap<String, String> fileResults = new ConcurrentHashMap<String, String>();
			ConcurrentHashMap<String, String> treeResults = new ConcurrentHashMap<String, String>();
			for (int i = 0; i < auto.length; i++) {
				boolean[] prevHand = i == 0 ? ArrayUtils.booleanArray(hand[i].length, true) : hand[i - 1];
				boolean[] prevAuto = i == 0 ? ArrayUtils.booleanArray(auto[i].length, true) : auto[i - 1];
				boolean[] handV = hand[i];
				boolean[] autoV = auto[i];
				int[] conc = concordance(handV, autoV);
				int[] prevConc = prevConcordance(prevHand, handV, prevAuto, autoV);
				fileResults.put(gateNames[i], ArrayUtils.toStr(conc, "\t"));
				treeResults.put(gateNames[i], ArrayUtils.toStr(prevConc, "\t"));
			}
			results.put(sn.fcsFile, fileResults);
			resultsTree.put(sn.fcsFile, treeResults);
			
			System.out.println("Done with " + sn.fcsFile);
			
			d = null;
			System.gc();
		}
		
		private int[] prevConcordance(boolean[] prevH, boolean[] h, boolean[] prevA, boolean[] a) {
			int tp = 0, tn = 0, fp = 0, fn = 0;
			for (int i = 0; i < h.length; i++) {
				if (prevH[i] && prevA[i]) {
					if (h[i] && a[i]) {
						tp++;
					} else if (h[i] && !a[i]) {
						fn++;
					} else if (!h[i] && a[i]) {
						fp++;
					} else if (!h[i] && !a[i]) {
						tn++;
					} 
				}
			}
			return new int[]{tp, tn, fp, fn};
		}
		
		/**
		 * @param hand
		 * @param auto
		 * @return TP, TN, FP, FN
		 */
		private int[] concordance(boolean[] hand, boolean[] auto) {
			int tp = 0, tn = 0, fp = 0, fn = 0;
			for (int i = 0; i < hand.length; i++) {
				if (hand[i] && auto[i]) {
					tp++;
				} else if (hand[i] && !auto[i]) {
					fn++;
				} else if (!hand[i] && auto[i]) {
					fp++;
				} else if (!hand[i] && !auto[i]) {
					tn++;
				} 
			}
			return new int[]{tp, tn, fp, fn};
		}
		
		private void loadAutoValues(SampleNode sn) {
			String[] pts = sn.fcsFile.split("_");
			String fNum = null;
			for (String p : pts) {
				if (p.startsWith("F")) {
					fNum = p;
					break;
				}
			}
			if (fNum == null) {
				System.err.println("Error - couldn't parse F-number for fcs file: " + sn.fcsFile);
				return;
			}
			String file = null;
			for (String s : filesInAutoDir) {
				pts = s.split("_");
				for (String p : pts) {
					if (fNum.equals(p)) {
						file = s;
						break;
					}
				}
			}
			if (file == null) {
				System.err.println("Error - Couldn't match F-number of fcs file with auto-gated files: " + sn.fcsFile + " | " + fNum);
				return;
			}
			autoData = HashVec.loadFileToStringMatrix(this.autoDir + file, true, null, false);
		}
		
	}

	class DuplicatorProcessor extends AbstractSampleProcessor {
		@Override
		public void processSample(SampleNode sn, Logger log) throws IOException {
			loadPopsAndGates(sn);
			loadData(sn);

			String[] gatings = ArrayUtils.stringArray(d.getCount(), "");

			HashSet<Gate> leaves = sn.gating.getAllLeafGates();
			for (Gate g : leaves) {
				boolean[] b = g.gate(d);
				for (int i = 0; i < b.length; i++) {
					if (b[i]) {
						if ("".equals(gatings[i])) {
							gatings[i] = g.getName();
						} else {
							int g1 = countDepth(g);
							int g2 = countDepth(sn.gating.gateMap.get(gatings[i]));
							if (g1 > g2) {
								gatings[i] = g.getName();
							}
						}
					}
				}
			}

			String file = sn.fcsFile;
			String newFile = ext.rootOf(file, false) + "_gated.fcs";
			FCSFileDuplicator.createFrom(file, newFile, log, gatings);

		}

		private int countDepth(Gate g) {
			int cnt = 0;
			Gate g1 = g;
			while (g1.getParentGate() != null) {
				cnt++;
				g1 = g1.getParentGate();
			}
			return cnt;
		}
	}



	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		
//	 private static final String FCS = "/home/pankrat2/shared/flow/fcs2/";
//	private static final String FCS = "/scratch.global/cole0482/FCS/src/";
//	private static final String FCS = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD FCS HRS CTL/";
//	private static final String WSP = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD WSP/";
//	private static final String OUT = "/scratch.global/cole0482/FCS/testConcordance/";
		
		String fcs = "/home/pankrat2/shared/flow/fcs2/";
		String wsp = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD WSP/";
		String auto = "/home/pankrat2/shared/flow/testAutoGate/test2/gates2/";
		String out = "/scratch.global/cole0482/FCS/testConcordance/";
		
		for (String arg : args) {
			if (arg.startsWith("fcs=")) {
				fcs = ext.parseStringArg(arg);
				numArgs--;
			} else if (arg.startsWith("wsp=")) {
				wsp = ext.parseStringArg(arg);
				numArgs--;
			} else if (arg.startsWith("auto=")) {
				auto = ext.parseStringArg(arg);
				numArgs--;
			} else if (arg.startsWith("out=")) {
				out = ext.parseStringArg(arg);
				numArgs--;
			}
		}
		
		new GateExportTest(fcs, wsp, auto, out).run();
	}

}
