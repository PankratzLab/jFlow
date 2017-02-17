package org.genvisis.one.ben.fcs.auto;

import java.io.IOException;
import java.util.HashSet;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Logger;
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

	// private static final String FCS = "/home/pankrat2/shared/flow/fcs/";
	private static final String FCS = "/scratch.global/cole0482/FCS/src/";
	private static final String WSP = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD WSP/";
	private static final String OUT = "/scratch.global/cole0482/FCS/withGates/";

	private void run() throws IOException {

		ProcessorFactory<? extends SampleProcessor> pf;

		// pf = new GateAssignmentFactory();
		// pf = new PercentageWriterFactory();
		// pf = new LeafDataSamplerFactory("/scratch.global/cole0482/FCS/", new Logger());
		pf = new ProcessorFactory<DuplicatorProcessor>() {
			public void cleanup(Object owner) {};

			@Override
			public DuplicatorProcessor createProcessor(Object owner, int index) {
				return new DuplicatorProcessor();
			}
		};


		SamplingPipeline sp = new SamplingPipeline(1, null, WSP, FCS, null, OUT, pf);

		sp.run();

		while (!sp.checkFinished()) {
			Thread.yield();
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
