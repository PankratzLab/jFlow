package one.JL;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;

import javax.jms.IllegalStateException;

import stats.Rscript.COLUMNS_MULTIPLOT;
import stats.Rscript.PLOT_DEVICE;
import stats.Rscript.RScatter;
import stats.Rscript.RScatters;
import stats.Rscript.SCATTER_TYPE;
import common.Array;
import common.Files;
import common.ext;
import common.PSF.Ext;
import cnv.analysis.CentroidCompute;
import cnv.analysis.CentroidCompute.Builder;
import cnv.filesys.Project;
import cnv.qc.GcAdjustor.GCAdjustorBuilder;
import cnv.qc.GcAdjustor.GC_CORRECTION_METHOD;
import cnv.qc.GcAdjustor.GcModel;
import cnv.qc.GcAdjustorParameter;
import cnv.qc.GcAdjustorParameter.GcAdjustorParameters;

public class GCcorrectionIterator {

	private static void batch(Project proj, String outputRootDir, int[] bpModels, int numThreads) {
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
			Files.qsub(currentSub, Array.toStr(Array.toStringArray(command), " "), 245000, 48.00, numThreads);
		}
		String batchMaster = batchRoot + "master.pbs";
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(batchMaster));
			for (int i = 0; i < pbs.size(); i++) {
				writer.println("qsub -q ram256g " + pbs.get(i));
			}
			writer.close();
		} catch (Exception e) {
			proj.getLog().reportError("Error writing to " + batchMaster);
			proj.getLog().reportException(e);
		}
		Files.chmod(batchMaster);
	}

	private static void iterate(Project proj, String outputRootDir, int[] bpModels, int numThreads) {
		new File(proj.PROJECT_DIRECTORY.getValue() + outputRootDir).mkdirs();
		int[] regressDistance = new int[] { 10, 100, 1000, 2000, 4000, 8000, 10000, 20000, 40000, 80000, 100000, 500000, 1000000 };// eq
		int[] snpMAD = new int[] { 0, 1, 2, 5, 10, 15 };// gt
		String freshCents = proj.PROJECT_DIRECTORY.getValue() + outputRootDir + "freshCents.cent";

		proj.getLog().reportTimeInfo("total iterations currently at (2X) " + (bpModels.length * regressDistance.length * snpMAD.length));
		if (!Files.exists(freshCents)) {
			CentroidCompute.computeAndDumpCentroids(proj, null, freshCents, new Builder(), numThreads, 2);
		}
		ArrayList<IterationParameters> finals = new ArrayList<IterationParameters>();
		GCAdjustorBuilder builder = new GCAdjustorBuilder();
		builder.verbose(true);
		for (int i = 0; i < bpModels.length; i++) {
			String model = proj.PROJECT_DIRECTORY.getValue() + outputRootDir + "gcmodel_bp_" + bpModels[i] + ".ser";
			GcModel gcModel = null;
			if (!Files.exists(model)) {
				gcModel = GcModel.generateSnpWindowModel(proj, bpModels[i]);
				gcModel.Serialize(model);
			} else {
				gcModel = GcModel.loadSerial(model);
			}
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
					break;
				}
				break;
			}
			proj.getLog().reportTimeInfo("Beginnning iteration group for gc model " + bpModels[i] + " (" + builders.size() + " iterations");
			String[][][] generated = GcAdjustorParameter.generateAdjustmentParameters(proj, builders.toArray(new GCAdjustorBuilder[builders.size()]), new String[] { freshCents }, new GC_CORRECTION_METHOD[] { GC_CORRECTION_METHOD.GENVISIS_GC }, gcModel, Array.toStringArray(outs), numThreads);
			IterationParameters[] tmp = getParameters(generated, bpModels[i], builders);
			for (int j = 0; j < tmp.length; j++) {
				finals.add(tmp[j]);
			}
		}
		String outputGz = proj.PROJECT_DIRECTORY.getValue() + outputRootDir + "finalSummaryRaw.gz";

		try {
			summarize(proj, outputGz, finals);
		} catch (IllegalStateException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static void summarize(Project proj, String outputGZ, ArrayList<IterationParameters> finals) throws IllegalStateException {
		String[] commonHeader = new String[] { "SampleName", "gcmodel_bp", "regress_bp", "snpMAD" };

		String[] specificHeader = new String[] { "BETA_0", "BETA_1", "WF_PRIOR", "WF_POST", "GCWF_PRIOR", "GCWF_POST", "LRR_MEAN_PRIOR", "LRR_MEAN_POST", "LRR_SD_PRIOR", "LRR_SD_POST" };
		String centTag = "_Cent";
		String[] withoutCent = Array.tagOn(specificHeader, null, "");
		String[] withCent = Array.tagOn(specificHeader, null, centTag);
		PrintWriter writer = Files.getAppropriateWriter(outputGZ);
		writer.println(Array.toStr(commonHeader) + "\t" + Array.toStr(withoutCent) + "\t" + Array.toStr(withCent));

		for (int i = 0; i < finals.size(); i++) {
			IterationParameters cur = finals.get(i);
			if (cur.getSerFiles().length != 2) {
				throw new IllegalStateException("Ser replicates must be in two-fers");

			} else {
				GcAdjustorParameters noCents = GcAdjustorParameters.readSerial(cur.getSerFiles()[0], proj.getLog());
				GcAdjustorParameters cents = GcAdjustorParameters.readSerial(cur.getSerFiles()[1], proj.getLog());
				String[] allSamples = proj.getSamples();
				GcAdjustorParameter[] noCentParams = noCents.getGcAdjustorParameters();
				GcAdjustorParameter[] centParams = cents.getGcAdjustorParameters();
				for (int j = 0; j < noCentParams.length; j++) {
					GcAdjustorParameter noC = noCentParams[j];
					GcAdjustorParameter c = centParams[j];
					if (!allSamples[i].equals(noC.getSample()) || !allSamples[i].equals(c.getSample())) {
						throw new IllegalStateException("MisMatched sample order");
					} else {
						writer.print(noC.getSample() + "\t" + Array.toStr(cur.getParams()) + "\t" + Array.toStr(noC.getQCString()) + "\t" + Array.toStr(c.getQCString()));
					}
				}
			}
		}
		writer.close();
		ArrayList<RScatter> allLooks = new ArrayList<RScatter>();
		String meanRoot = ext.rootOf(outputGZ, false) + "means";

		RScatter rScatter = new RScatter(outputGZ, meanRoot + ".rscript", ext.removeDirectoryInfo(meanRoot), meanRoot + ".pdf", "LRR_MEAN_PRIOR", new String[] { "LRR_MEAN_POST", "LRR_MEAN_PRIOR" + centTag, "LRR_MEAN_POST" + centTag }, null, SCATTER_TYPE.POINT, proj.getLog());
		rScatter.execute();
		allLooks.add(rScatter);
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

		public IterationParameters(int bpModel, int regressDistance, int numSnpMad, String[] serFiles) {
			super();
			this.bpModel = bpModel;
			this.regressDistance = regressDistance;
			this.numSnpMad = numSnpMad;
			this.serFiles = serFiles;
		}

		public String[] getSerFiles() {
			return serFiles;
		}

		public String[] getParams() {
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
		int[] bpModels = new int[] { 50, 100, 250, 500, 1000, 2500, 5000, 10000, 1000000 };
		boolean batch = false;
		String usage = "\n" + "one.JL.GCcorrectionIterator requires 0-1 arguments\n";
		usage += "   (1) project filename (i.e. proj= (no default))\n" + "";
		usage += "   (2) root directory under project directory (i.e. root=" + rootDir + " (default))\n" + "";
		usage += "   (3) gcModel bp, comma delimited (i.e. bpGcModel=" + Array.toStr(Array.toStringArray(bpModels), ",") + " (default))\n" + "";
		usage += "   (4) batch by gc model (i.e -batch, not the default)\n" + "";

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
				iterate(proj, rootDir, bpModels, numThreads);
			} else {
				proj.PROJECT_PROPERTIES_FILENAME.setValue(filename);
				batch(proj, rootDir, bpModels, numThreads);
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

}
