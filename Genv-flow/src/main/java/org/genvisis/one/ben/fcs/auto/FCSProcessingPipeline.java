package org.genvisis.one.ben.fcs.auto;

import java.io.IOException;

import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.auto.proc.ProcessorFactory;
import org.genvisis.one.ben.fcs.auto.proc.SampleProcessor;
import org.genvisis.one.ben.fcs.auto.proc.VisualizationProcessor;

public class FCSProcessingPipeline {

	public static final String[][] MATCH = {
																					{"lymph", "Lymphocytes (SSC-A v FSC-A)"},
																					{"Singlets", "Single Cells (FSC-H v FSC-W)"},
																					{"PE.A", "Live cells (PE-)"},
																					{"CD3.", "Tcells (CD3+ CD19-)"},
	};

	public static String[] GATE_NAMES = {
																			 "Lymphocytes (SSC-A v FSC-A)",
																			 "Single Cells (FSC-H v FSC-W)",
																			 "Live cells (PE-)",
																			 "Tcells (CD3+ CD19-)",
	};

	public FCSProcessingPipeline(String fcs, String wsp, String auto, String out) {
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
		// pf = new ProcessorFactory<ConcordanceProcessor>() {
		//
		// ConcurrentHashMap<String, ConcurrentHashMap<String, String>> resultMap = new
		// ConcurrentHashMap<>();
		// ConcurrentHashMap<String, ConcurrentHashMap<String, String>> resultMapTree = new
		// ConcurrentHashMap<>();
		// public void cleanup(Object owner) {
		// PrintWriter writer = Files.getAppropriateWriter(OUT + "concordance.xln");
		//
		// Set<String> files = resultMap.keySet();
		// for (String s : files) {
		// for (int i = 0; i < MATCH.length; i++) {
		// writer.println(s + "\t" + MATCH[i][1] + "\t" + resultMap.get(s).get(MATCH[i][1]));
		// }
		// }
		// writer.flush();
		// writer.close();
		//
		// writer = Files.getAppropriateWriter(OUT + "concordance_tree.xln");
		//
		// files = resultMapTree.keySet();
		// for (String s : files) {
		// for (int i = 0; i < MATCH.length; i++) {
		// writer.println(s + "\t" + MATCH[i][1] + "\t" + resultMapTree.get(s).get(MATCH[i][1]));
		// }
		// }
		// writer.flush();
		// writer.close();
		// };
		//
		// @Override
		// public ConcordanceProcessor createProcessor(Object owner, int index) {
		// return new ConcordanceProcessor(resultMap, resultMapTree);
		// }
		// };
		//
		// SamplingPipeline sp = new SamplingPipeline(1, null, WSP, FCS, null, OUT, pf);

		pf = new ProcessorFactory<SampleProcessor>() {

			@Override
			public void cleanup(Object owner) {}

			@Override
			public SampleProcessor createProcessor(Object owner, int index) {
				return new VisualizationProcessor(autoDir, outDir);
			}
		};

		SamplingPipeline sp = new SamplingPipeline(1, null, wspDir, fcsDir, null, outDir, pf);

		sp.run();

		while (!sp.checkFinished()) {
			Thread.yield();
		}

	}

	public static String getFNum(String file) {
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


	public static void main(String[] args) throws IOException {
		int numArgs = args.length;

		// private static final String FCS = "/home/pankrat2/shared/flow/fcs2/";
		// private static final String FCS = "/scratch.global/cole0482/FCS/src/";
		// private static final String FCS =
		// "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD FCS HRS CTL/";
		// private static final String WSP = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD WSP/";
		// private static final String OUT = "/scratch.global/cole0482/FCS/testConcordance/";

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

		new FCSProcessingPipeline(fcs, wsp, auto, out).run();
	}

}