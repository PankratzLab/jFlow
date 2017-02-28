package org.genvisis.one.ben.fcs.auto;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.FCSDataLoader;
import org.genvisis.one.ben.fcs.FCSFileDuplicator;
import org.genvisis.one.ben.fcs.auto.LeafDataSampler.GateAssignmentFactory;
import org.genvisis.one.ben.fcs.gating.Gate;
import org.genvisis.one.ben.fcs.gating.Workbench.SampleNode;

public class GateExportTest {

	// private static final String WSP = "F:/Flow/gatingTest/";
	// private static final String FCS = "F:/Flow/gatingTest/2016-07-13_PANEL 1_ZF_Group
	// one_F1631196_002.fcs";

	 private static final String FCS = "/home/pankrat2/shared/flow/fcs2/";
//	private static final String FCS = "/scratch.global/cole0482/FCS/src/";
	private static final String WSP = FCS;// "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD WSP/";
	private static final String OUT = "/scratch.global/cole0482/FCS/testConcordance/";

	private void run() throws IOException {

		ProcessorFactory<? extends SampleProcessor> pf;

		// pf = new GateAssignmentFactory();
		// pf = new PercentageWriterFactory();
		// pf = new LeafDataSamplerFactory("/scratch.global/cole0482/FCS/", new Logger());
		pf = new ProcessorFactory<ConcordanceProcessor>() {

			ConcurrentHashMap<String, String> resultMap = new ConcurrentHashMap<>();
			public void cleanup(Object owner) {
				PrintWriter writer = Files.getAppropriateWriter(OUT + "concordance.xln");
				for (Entry<String, String> ent : resultMap.entrySet()) {
					writer.println(ent.getKey() + "\t" + ent.getValue()); 
				}
				writer.flush();
				writer.close();
			};

			@Override
			public ConcordanceProcessor createProcessor(Object owner, int index) {
				return new ConcordanceProcessor(resultMap);
			}
		};


		SamplingPipeline sp = new SamplingPipeline(1, null, WSP, FCS, null, OUT, pf);

		sp.run();

		while (!sp.checkFinished()) {
			Thread.yield();
		}

	}

	
	class ConcordanceProcessor extends AbstractSampleProcessor {
		private static final String AUTO_DATA_DIR = "/home/pankrat2/shared/flow/testAutoGate/test2/gates2/";
		private static final String AUTO_FILE_SUFF = "def.txt.gz";
		
		private String[] filesInAutoDir = (new File(AUTO_DATA_DIR)).list(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return name.endsWith(AUTO_FILE_SUFF);
			}
		});
		
		String[][] autoData;
		private ConcurrentHashMap<String, String> results;
		
		public ConcordanceProcessor(ConcurrentHashMap<String, String> resultMap) {
			this.results = resultMap;
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
			// TODO refactor for gate header matching
			String[] gateNames = new String[hand.length];
			Gate gate = sn.gating.getRootGates().get(0);
			for (int i = 0; i < hand.length; i++) {
				hand[i] = gate.gate(d);
				gateNames[i] = gate.getName();
				gate = gate.getChildGates().get(0);
			}
			
			boolean[][] auto = new boolean[autoData[0].length][];
			for (int i = 0; i < auto.length; i++) {
				String[] data = Matrix.extractColumn(autoData, i);
				auto[i] = new boolean[data.length];
				for (int k = 0; k < data.length; k++) {
					auto[i][k] = Boolean.valueOf(data[k]); 
				}
			}
			
			for (int i = 0; i < auto.length; i++) {
				boolean[] handV = hand[i];
				boolean[] autoV = auto[i];
				int[] conc = concordance(handV, autoV);
				results.put(sn.fcsFile + "\t" + gateNames[i], ArrayUtils.toStr(conc, "\t"));
			}
			
			System.out.println("Done with " + sn.fcsFile);
			
			d = null;
			System.gc();
		}
		
		/**
		 * @param std
		 * @param chk
		 * @return TP, TN, FP, FN
		 */
		private int[] concordance(boolean[] std, boolean[] chk) {
			int tp = 0, tn = 0, fp = 0, fn = 0;
			for (int i = 0; i < std.length; i++) {
				if (std[i] && chk[i]) {
					tp++;
				} else if (std[i] && !chk[i]) {
					fn++;
				} else if (!std[i] && chk[i]) {
					fp++;
				} else if (!std[i] && !chk[i]) {
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
			autoData = HashVec.loadFileToStringMatrix(AUTO_DATA_DIR + file, true, null, false);
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
		new GateExportTest().run();
	}

}
