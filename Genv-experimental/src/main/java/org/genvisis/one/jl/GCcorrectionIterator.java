package org.genvisis.one.jl;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.concurrent.Callable;

import javax.jms.IllegalStateException;

import org.genvisis.cnv.analysis.CentroidCompute;
import org.genvisis.cnv.analysis.CentroidCompute.CentroidBuilder;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.GcAdjustor;
import org.genvisis.cnv.qc.GcAdjustorParameter;
import org.genvisis.cnv.qc.GcAdjustor.GCAdjustorBuilder;
import org.genvisis.cnv.qc.GcAdjustor.GC_CORRECTION_METHOD;
import org.genvisis.cnv.qc.GcAdjustor.GcModel;
import org.genvisis.cnv.qc.GcAdjustorParameter.GcAdjustorParameters;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.ext;
import org.genvisis.common.PSF.Ext;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.stats.Rscript;
import org.genvisis.stats.Rscript.COLUMNS_MULTIPLOT;
import org.genvisis.stats.Rscript.PLOT_DEVICE;
import org.genvisis.stats.Rscript.RScatter;
import org.genvisis.stats.Rscript.RScatters;
import org.genvisis.stats.Rscript.SCATTER_TYPE;

public class GCcorrectionIterator {
	private static final String CENT_TAG = "_CENT";

	private static void batch(Project proj, String outputRootDir, int[] bpModels, int[] regressDistance, int[] snpMAD, int numThreads) {
		String batchRoot = proj.PROJECT_DIRECTORY.getValue() + outputRootDir;
		ArrayList<String> pbs = new ArrayList<String>();
		for (int i = 0; i < bpModels.length; i++) {
			String currentSub = batchRoot + "gcmodel_bp_" + bpModels[i] + ".pbs";
			pbs.add(currentSub);
			ArrayList<String> command = new ArrayList<String>();
			command.add("java -cp ~/parkGC.jar one.JL.GCcorrectionIterator ");
			command.add("proj=" + proj.PROJECT_PROPERTIES_FILENAME.getValue());
			command.add("numthreads=" + numThreads);
			command.add("bpGcModel=" + bpModels[i]);
			Files.qsub(currentSub, Array.toStr(Array.toStringArray(command), " "), 55000, 48.00, numThreads);
		}
		String batchMaster = batchRoot + "master.pbs";
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(batchMaster));
			for (int i = 0; i < pbs.size(); i++) {
				writer.println("qsub -q small " + pbs.get(i));
			}
			writer.close();
		} catch (Exception e) {
			proj.getLog().reportError("Error writing to " + batchMaster);
			proj.getLog().reportException(e);
		}
		Files.chmod(batchMaster);
	}

	private static void iterate(Project proj, String outputRootDir, int[] bpModels, int[] regressDistance, int[] snpMAD, int numThreads) {
		new File(proj.PROJECT_DIRECTORY.getValue() + outputRootDir).mkdirs();
		String freshCents = proj.PROJECT_DIRECTORY.getValue() + outputRootDir + "freshCents.cent";

		proj.getLog().reportTimeInfo("total iterations currently at (2X) " + (bpModels.length * regressDistance.length * snpMAD.length));
		if (!Files.exists(freshCents)) {
			CentroidCompute.computeAndDumpCentroids(proj, freshCents, new CentroidBuilder(), numThreads, 2);
		}
		ArrayList<IterationParameters> finals = new ArrayList<IterationParameters>();
		GCAdjustorBuilder builder = new GCAdjustorBuilder();
		builder.verbose(true);
		String outputGz = proj.PROJECT_DIRECTORY.getValue() + outputRootDir + "finalSummaryRaw.gz";

		if (!Files.exists(outputGz)) {
			for (int i = 0; i < bpModels.length; i++) {
				String model = proj.PROJECT_DIRECTORY.getValue() + outputRootDir + "gcmodel_bp_" + bpModels[i] + ".ser";
				GcModel gcModel = null;
				if (bpModels[i] < 0) {
					gcModel = GcAdjustor.GcModel.populateFromFile(proj.GC_MODEL_FILENAME.getValue(), true, proj.getLog());
					gcModel.Serialize(model);
				}
				// if (Files.exists(model)) {
				// proj.getLog().reportTimeError("JOHN remember to remove this gcmodel skipper");

				if (!Files.exists(model)) {
					gcModel = GcModel.generateSnpWindowModel(proj, bpModels[i]);
					gcModel.Serialize(model);
				} else {
					gcModel = GcModel.loadSerial(model);
				}
				// }
				ArrayList<GCAdjustorBuilder> builders = new ArrayList<GCAdjustorBuilder>();
				ArrayList<String> outs = new ArrayList<String>();
				for (int j = 0; j < regressDistance.length; j++) {
					for (int j2 = 0; j2 < snpMAD.length; j2++) {
						String root = outputRootDir + "gcmodel_bp_" + bpModels[i] + "_regress_" + regressDistance[j] + "_snpMad_" + snpMAD[j2];
						outs.add(root);
						builder.regressionDistance(regressDistance[j]);
						builder.numSnpMAD(snpMAD[j2]);
						builders.add(new GCAdjustorBuilder(builder));
						proj.getLog().reportTimeError("JOHN remember to remove this");
					}
				}
				proj.getLog().reportTimeInfo("Beginnning iteration group for gc model " + bpModels[i] + " (" + builders.size() + " iterations");
				String[][][] generated = null;

				generated = GcAdjustorParameter.generateAdjustmentParameters(proj, builders.toArray(new GCAdjustorBuilder[builders.size()]), new String[] { freshCents }, new GC_CORRECTION_METHOD[] { GC_CORRECTION_METHOD.GENVISIS_GC }, gcModel, Array.toStringArray(outs), numThreads, true);
				proj.getLog().reportTimeInfo("Adding gc Content");
				double[] gcs = gcModel.getGCsFor(proj.getMarkerNames());
//				for (int j = 0; j < generated.length; j++) {
//					for (int j2 = 0; j2 < generated[j].length; j2++) {
//						for (int k = 0; k < generated[j][j2].length; k++) {
//							GcAdjustorParameters tmp = GcAdjustorParameters.readSerial(generated[j][j2][k], proj.getLog());
//							tmp.setGcContent(gcs);
//							tmp.writeSerial(generated[j][j2][k]);
//						}
//					}
//				}
				IterationParameters[] tmp = getParameters(generated, bpModels[i], builders);
				for (int j = 0; j < tmp.length; j++) {
					finals.add(tmp[j]);
				}

			}
		}

		try {
			summarizeSampleQC(proj, outputGz, numThreads,finals);
		} catch (IllegalStateException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	// private static void GenerateMarkerQC(Project proj, ArrayList<IterationParameters> finals, String finalOut) {
	// String[] header = new String[] { "MarkerName", "CENT", "gcmodel_bp", "regress_bp", "snpMAD", "LRR_MEAN_PRIOR", "LRR_MEAN_POST", "LRR_SD_PRIOR", "LRR_SD_POST" };
	//
	// }

	private static void format(String base, Hashtable<String, String[]> plotCombos) {
		plotCombos.put(base + "_PRIOR", new String[] { base + "_POST", base + "_PRIOR" + CENT_TAG, base + "_POST" + CENT_TAG });

	}

	private static class IPloadProducer extends AbstractProducer<IterationParameters> {
		private ArrayList<IterationParameters> finals;
		private Logger log;
		private int index;

		public IPloadProducer(ArrayList<IterationParameters> finals, Logger log) {
			super();
			this.finals = finals;
			this.log = log;
		}

		@Override
		public boolean hasNext() {
			return index < finals.size();
		}

		@Override
		public Callable<IterationParameters> next() {
			IPLoadWorker lw = new IPLoadWorker(finals.get(index), log);
			index++;
			// TODO Auto-generated method stub
			return lw;
		}
	}

	private static class IPLoadWorker implements Callable<IterationParameters> {
		private IterationParameters ip;
		private Logger log;

		public IPLoadWorker(IterationParameters ip, Logger log) {
			super();
			this.ip = ip;
			this.log = log;
		}

		@Override
		public IterationParameters call() throws Exception {
			ip.load(log);
			// TODO Auto-generated method stub
			return ip;
		}

	}

	private static void summarizeSampleQC(Project proj, String outputGZ, int numthreads, ArrayList<IterationParameters> finals) throws IllegalStateException {
		String[] commonHeader = new String[] { "SampleName", "gcmodel_bp", "regress_bp", "snpMAD" };

		String[] specificHeader = new String[] { "BETA_0", "BETA_1", "WF_PRIOR", "WF_POST", "GCWF_PRIOR", "GCWF_POST", "LRR_MEAN_PRIOR", "LRR_MEAN_POST", "LRR_SD_PRIOR", "LRR_SD_POST" };
		Hashtable<String, String[]> plotCombos = new Hashtable<String, String[]>();
		format("WF", plotCombos);
		format("GCWF", plotCombos);
		format("LRR_MEAN", plotCombos);
		format("LRR_SD", plotCombos);
		IPloadProducer ipp = new IPloadProducer(finals, proj.getLog());
		WorkerTrain<IterationParameters> train = new WorkerTrain<GCcorrectionIterator.IterationParameters>(ipp, numthreads, 2, proj.getLog());
		if (!Files.exists(outputGZ)) {
			String[] withoutCent = Array.tagOn(specificHeader, null, "");
			String[] withCent = Array.tagOn(specificHeader, null, CENT_TAG);
			PrintWriter writer = Files.getAppropriateWriter(outputGZ);
			writer.println(Array.toStr(commonHeader) + "\t" + Array.toStr(withoutCent) + "\t" + Array.toStr(withCent));
			int index = 0;
			while (train.hasNext()) {
				IterationParameters cur = train.next();
				if (index % 50 == 0) {
					proj.getLog().reportTimeInfo("Summarized " + index);
				}
				if (cur.getSerFiles().length != 2) {
					throw new IllegalStateException("Ser replicates must be in two-fers");

				} else {
					GcAdjustorParameters noCents = cur.getNoCents();
					GcAdjustorParameters cents = cur.getCents();
					String[] allSamples = proj.getSamples();
					GcAdjustorParameter[] noCentParams = noCents.getGcAdjustorParameters();
					GcAdjustorParameter[] centParams = cents.getGcAdjustorParameters();
					for (int j = 0; j < noCentParams.length; j++) {
						GcAdjustorParameter noC = noCentParams[j];
						GcAdjustorParameter c = centParams[j];
						if (!allSamples[j].equals(noC.getSample()) || !allSamples[j].equals(c.getSample())) {
							throw new IllegalStateException("MisMatched sample order");
						} else {
							writer.println(noC.getSample() + "\t" + Array.toStr(cur.getParams()) + "\t" + Array.toStr(noC.getQCString()) + "\t" + Array.toStr(c.getQCString()));
						}
					}
				}
				finals.set(index, null);
				index++;
				proj.getLog().memoryPercentTotalFree();
			}
			writer.close();
		}
		ArrayList<RScatter> allLooks = new ArrayList<RScatter>();

		for (String base : plotCombos.keySet()) {
			
//			System.out.println(base);
//			System.out.println(Array.toStr(plotCombos.get(base)));
//			// String root = ext.rootOf(outputGZ, false) + base;
//			// RScatter rScatter = new RScatter(outputGZ, root + ".rscript", ext.removeDirectoryInfo(root), root + ".pdf", base, plotCombos.get(base), null, SCATTER_TYPE.POINT, proj.getLog());
//			// rScatter.setyLabel(base.replaceAll("_PRIOR", ""));
//			// rScatter.setxLabel(base);
//			// rScatter.execute();
//			// allLooks.add(rScatter);
//			// //
//			String rootBox = ext.rootOf(outputGZ, false) + base + ".box";
//			RScatter rScatterBox = new RScatter(outputGZ, rootBox + ".rscript", ext.removeDirectoryInfo(rootBox), rootBox + ".pdf", "gcmodel_bp", Array.concatAll(new String[] { base }, plotCombos.get(base)), null, SCATTER_TYPE.BOX, proj.getLog());
//			rScatterBox.setyLabel(rootBox.replaceAll("_PRIOR", ""));
//			rScatterBox.setFontsize(3);
//			rScatterBox.setyLabel(base.replaceAll("_PRIOR", ""));
//			rScatterBox.setTitle("Original AND re-computed LRR");
//			// rScatterBox.setxLabel(rootBox);
//			rScatterBox.execute();
//			allLooks.add(rScatterBox);
//
//			String rootBoxTrim = ext.rootOf(outputGZ, false) + base + ".trim";
//			RScatter rScatterBoxTrim = new RScatter(outputGZ, rootBoxTrim + ".rscript", ext.removeDirectoryInfo(rootBoxTrim), rootBoxTrim + ".pdf", "gcmodel_bp", new String[] { plotCombos.get(base)[0], plotCombos.get(base)[2] }, null, SCATTER_TYPE.BOX, proj.getLog());
//			rScatterBoxTrim.setyLabel(rootBoxTrim.replaceAll("_PRIOR", ""));
//			rScatterBoxTrim.setFontsize(3);
//			rScatterBoxTrim.setyLabel(base.replaceAll("_PRIOR", ""));
//			rScatterBoxTrim.setTitle("re-computed LRR Only");
//
//			// rScatterBox.setxLabel(rootBox);
//			rScatterBoxTrim.execute();
//			allLooks.add(rScatterBoxTrim);
//
//			String rootBoxTrimTrim = ext.rootOf(outputGZ, false) + base + ".trimTrim";
//			RScatter rScatterBoxTrimTrim = new RScatter(outputGZ, rootBoxTrimTrim + ".rscript", ext.removeDirectoryInfo(rootBoxTrimTrim), rootBoxTrimTrim + ".pdf", "gcmodel_bp", new String[] { plotCombos.get(base)[0], plotCombos.get(base)[2] }, null, SCATTER_TYPE.BOX, proj.getLog());
//			rScatterBoxTrimTrim.setyLabel(rootBoxTrimTrim.replaceAll("_PRIOR", ""));
//			rScatterBoxTrimTrim.setFontsize(3);
//			rScatterBoxTrimTrim.setyLabel(base.replaceAll("_PRIOR", ""));
//			rScatterBoxTrimTrim.setyRange(new double[] { -.5, .5 });
//			rScatterBoxTrimTrim.setTitle("re-computed LRR Only");
//
//			// rScatterBox.setxLabel(rootBox);
//			rScatterBoxTrimTrim.execute();
//			allLooks.add(rScatterBoxTrimTrim);
//
//			// String rootVsNoCent = ext.rootOf(outputGZ, false) + base + ".vs_noCent";
//			// RScatter rScattervs = new RScatter(outputGZ, rootVsNoCent + ".rscript", ext.removeDirectoryInfo(rootVsNoCent), rootVsNoCent + ".pdf", base, new String[] { plotCombos.get(base)[0] }, "gcmodel_bp", SCATTER_TYPE.POINT, proj.getLog());
//			// rScattervs.setxLabel(base);
//			// rScattervs.setyLabel(plotCombos.get(base)[0]);
//			// rScattervs.setTitle("NO Recompute");
//			// rScattervs.setFontsize(3);
//			// rScattervs.execute();
//			// allLooks.add(rScattervs);
//
//			String rootVsCent = ext.rootOf(outputGZ, false) + base + ".vs_Cent";
//			RScatter rScattervsCent = new RScatter(outputGZ, rootVsCent + ".rscript", ext.removeDirectoryInfo(rootVsCent), rootVsCent + ".pdf", plotCombos.get(base)[1], new String[] { plotCombos.get(base)[2] }, "gcmodel_bp", SCATTER_TYPE.POINT, proj.getLog());
//			rScattervsCent.setFontsize(3);
//			rScattervsCent.setxLabel("gcmodel_bp");
//			rScattervsCent.setyLabel(plotCombos.get(base)[2]);
//			rScattervsCent.setTitle("re-computed LRR Only");
//
//			// rScatterBox.setxLabel(rootBox);
//			rScattervsCent.execute();
//
//			String rootVsCentBP = ext.rootOf(outputGZ, false) + base + ".vs_Cent";
//			RScatter rScattervsCentBP = new RScatter(outputGZ, rootVsCentBP + ".rscript", ext.removeDirectoryInfo(rootVsCentBP), rootVsCentBP + ".pdf", "gcmodel_bp", new String[] { plotCombos.get(base)[2] }, "regress_bp", SCATTER_TYPE.POINT, proj.getLog());
//			rScattervsCentBP.setFontsize(3);
//			rScattervsCentBP.setxLabel("gcmodel_bp");
//			rScattervsCentBP.setyLabel(plotCombos.get(base)[2]);
//			rScattervsCentBP.setTitle("re-computed LRR Only");
//
//			allLooks.add(rScattervsCentBP);

			String rootVsCentDensity = ext.rootOf(outputGZ, false) + base + ".vs_CentDensity";
			RScatter rScattervsCentDensity = new RScatter(outputGZ, rootVsCentDensity + ".rscript", ext.removeDirectoryInfo(rootVsCentDensity), rootVsCentDensity + ".pdf", "gcmodel_bp", new String[] { "regress_bp" }, plotCombos.get(base)[2], SCATTER_TYPE.POINT, proj.getLog());
			rScattervsCentDensity.setFontsize(3);
			rScattervsCentDensity.setxLabel("regress_bp");
			rScattervsCentDensity.setyLabel(plotCombos.get(base)[2]);
			rScattervsCentDensity.setTitle("re-computed LRR Only");
			rScattervsCentDensity.setRestrictions(new Rscript.Restrictions[] { new Rscript.Restrictions(new String[] { "snpMAD" }, new double[] { 0 }, new String[] { "==" }, null) });
			rScattervsCentDensity.setScaleDensity(true);
			
			allLooks.add(rScattervsCentDensity);

		}

		String finalSummaryRoot = ext.rootOf(outputGZ, false) + "finalSummary";
		RScatters rscScatters = new RScatters(allLooks.toArray(new RScatter[allLooks.size()]), finalSummaryRoot + ".rscript", finalSummaryRoot + ".pdf", COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1, PLOT_DEVICE.PDF, proj.getLog());
		rscScatters.execute();

	}

	private static IterationParameters[] getParameters(String[][][] generated, int bpModel, ArrayList<GCAdjustorBuilder> builders) {
		IterationParameters[] params = new IterationParameters[builders.size()];
		for (int i = 0; i < generated.length; i++) {
			ArrayList<String> sers = new ArrayList<String>();
			for (int j = 0; j < generated[i].length; j++) {
				for (int j2 = 0; j2 < generated[i][j].length; j2++) {
					sers.add(generated[i][j][j2]);
				}
			}
			params[i] = new IterationParameters(bpModel, builders.get(i).getRegressionDistance(), builders.get(i).getNumSnpMAD(), Array.toStringArray(sers));

		}
		return params;
	}

	private static class IterationParameters {
		private int bpModel;
		private int regressDistance;
		private int numSnpMad;
		private String[] serFiles;
		private GcAdjustorParameters noCents;
		private GcAdjustorParameters cents;

		private IterationParameters(int bpModel, int regressDistance, int numSnpMad, String[] serFiles) {
			super();
			this.bpModel = bpModel;
			this.regressDistance = regressDistance;
			this.numSnpMad = numSnpMad;
			this.serFiles = serFiles;
		}

		private String[] getSerFiles() {
			return serFiles;
		}

		private void load(Logger log) {
			this.noCents = GcAdjustorParameters.readSerial(serFiles[0], log);
			this.cents = GcAdjustorParameters.readSerial(serFiles[1], log);
		}

		public GcAdjustorParameters getNoCents() {
			return noCents;
		}

		public GcAdjustorParameters getCents() {
			return cents;
		}

		private String[] getParams() {
			ArrayList<String> params = new ArrayList<String>();
			params.add(bpModel + "");
			params.add(regressDistance + "");
			params.add(numSnpMad + "");
			return Array.toStringArray(params);

		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String rootDir = "gcCorrectionIterations/";
		int numThreads = 24;
		String filename = null;
		int[] bpModels = new int[] { -1, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000 };

		System.out.println("JOHN add 1000000 and remove final.gz");
		bpModels = new int[] { -1,50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000, 100000, 250000, 500000 };
		int[] regressDistance = new int[] { 1000000, 10, 100, 1000, 2000, 4000, 8000, 10000, 20000, 40000, 80000, 100000, 500000 };// eq 13
		int[] snpMAD = new int[] { 0, 1, 2, 5, 10, 15 };
		boolean batch = false;
		String usage = "\n" + "one.JL.GCcorrectionIterator requires 0-1 arguments\n";
		usage += "   (1) project filename (i.e. proj= (no default))\n" + "";
		usage += "   (2) root directory under project directory (i.e. root=" + rootDir + " (default))\n" + "";
		usage += "   (3) gcModel bp, comma delimited (i.e. bpGcModel=" + Array.toStr(Array.toStringArray(bpModels), ",") + " (default))\n" + "";
		usage += "   (4) batch by gc model (i.e -batch, not the default)\n" + "";
		usage += "   (5) regressDistance, comma delimited (i.e. regress=" + Array.toStr(Array.toStringArray(regressDistance), ",") + " (default))\n" + "";
		usage += "   (6) snpMAD, comma delimited (i.e. mad=" + Array.toStr(Array.toStringArray(snpMAD), ",") + " (default))\n" + "";

		usage += Ext.getNumThreadsCommand(5, numThreads);
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("root=")) {
				rootDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("bpGcModel=")) {
				bpModels = Array.toIntArray(ext.parseStringArg(args[i], "").split(","));
				numArgs--;
			} else if (args[i].startsWith("regress=")) {
				regressDistance = Array.toIntArray(ext.parseStringArg(args[i], "").split(","));
				numArgs--;
			} else if (args[i].startsWith("mad=")) {
				snpMAD = Array.toIntArray(ext.parseStringArg(args[i], "").split(","));
				numArgs--;
			} else if (args[i].startsWith("-batch")) {
				batch = true;
				numArgs--;
			} else if (args[i].startsWith(Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			Project proj = new Project(filename, false);
			if (!batch) {
				iterate(proj, rootDir, bpModels, regressDistance, snpMAD, numThreads);
			} else {
				proj.PROJECT_PROPERTIES_FILENAME.setValue(filename);
				batch(proj, rootDir, bpModels, regressDistance, snpMAD, numThreads);
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	//
	// public static void main(String[] args) {
	//
	//
	// Project proj = new Project(args[0], false);
	// iterate(proj, "gcCorrectionIterations/", 24);
	// // generateAdjustmentParameters(proj, GC_CORRECTION_METHOD.values(), 4);
	// }

	// String meanRoot = ext.rootOf(outputGZ, false) + "means";
	//
	// RScatter rScattermean = new RScatter(outputGZ, meanRoot + ".rscript", ext.removeDirectoryInfo(meanRoot), meanRoot + ".pdf", "LRR_MEAN_PRIOR", new String[] { "LRR_MEAN_POST", "LRR_MEAN_PRIOR" + CENT_TAG, "LRR_MEAN_POST" + CENT_TAG }, null, SCATTER_TYPE.POINT, proj.getLog());
	// rScattermean.setxLabel("LRR_MEAN_PRIOR (no recompute)");
	// rScattermean.setyLabel("LRR_MEAN");
	// rScattermean.execute();
	// allLooks.add(rScattermean);
	//
	// String meanBoxRoot = ext.rootOf(outputGZ, false) + "meansBox";
	//
	// RScatter rScattermeanBox = new RScatter(outputGZ, meanBoxRoot + ".rscript", ext.removeDirectoryInfo(meanBoxRoot), meanBoxRoot + ".pdf", "gcmodel_bp", new String[] { "LRR_MEAN_PRIOR", "LRR_MEAN_POST", "LRR_MEAN_PRIOR" + CENT_TAG, "LRR_MEAN_POST" + CENT_TAG }, null, SCATTER_TYPE.BOX, proj.getLog());
	// rScattermeanBox.setyRange(new double[] { 0, .5 });
	//
	// rScattermeanBox.execute();
	//
	// allLooks.add(rScattermeanBox);
	//
	// String sdRoot = ext.rootOf(outputGZ, false) + "SD";
	//
	// RScatter rScatterSD = new RScatter(outputGZ, sdRoot + ".rscript", ext.removeDirectoryInfo(sdRoot), sdRoot + ".pdf", "LRR_SD_PRIOR", new String[] { "LRR_SD_POST", "LRR_SD_PRIOR" + CENT_TAG, "LRR_SD_POST" + CENT_TAG }, null, SCATTER_TYPE.POINT, proj.getLog());
	// rScatterSD.setxLabel("LRR_SD_PRIOR (no recompute)");
	// rScatterSD.setyLabel("LRR_SD");
	// rScatterSD.execute();
	// allLooks.add(rScatterSD);
	//
	// String sdBoxRoot = ext.rootOf(outputGZ, false) + "SDBox";
	//
	// RScatter rScattersdBox = new RScatter(outputGZ, sdBoxRoot + ".rscript", ext.removeDirectoryInfo(sdBoxRoot), sdBoxRoot + ".pdf", "gcmodel_bp", new String[] { "LRR_SD_PRIOR", "LRR_SD_POST", "LRR_SD_PRIOR" + CENT_TAG, "LRR_SD_POST" + CENT_TAG }, null, SCATTER_TYPE.BOX, proj.getLog());
	// rScattersdBox.setyRange(new double[] { 0, .5 });
	//
	// rScattersdBox.execute();
	//
	// allLooks.add(rScattersdBox);
	//

}
