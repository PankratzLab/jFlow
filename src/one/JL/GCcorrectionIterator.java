package one.JL;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;

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
		int[] regressDistance = new int[] { 10, 100, 1000, 2000, 4000, 8000, 10000, 20000, 40000, 80000, 100000, 500000, 1000000 };
		int[] snpMAD = new int[] { 0, 1, 2, 5, 10, 15 };
		String freshCents = proj.PROJECT_DIRECTORY.getValue() + outputRootDir + "freshCents.cent";
		proj.getLog().reportTimeInfo("total iterations currently at (2X) " + (bpModels.length * regressDistance.length * snpMAD.length));
		if (!Files.exists(freshCents)) {
			CentroidCompute.computeAndDumpCentroids(proj, null, freshCents, new Builder(), numThreads, 2);
		}
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
				}
			}
			proj.getLog().reportTimeInfo("Beginnning iteration group for gc model " + bpModels[i] + " (" + builders.size() + " iterations");
			String[][][] generated = GcAdjustorParameter.generateAdjustmentParameters(proj, builders.toArray(new GCAdjustorBuilder[builders.size()]), new String[] { freshCents }, new GC_CORRECTION_METHOD[] { GC_CORRECTION_METHOD.GENVISIS_GC }, gcModel, Array.toStringArray(outs), numThreads);
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
