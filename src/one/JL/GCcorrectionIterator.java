package one.JL;

import java.io.File;
import java.util.ArrayList;

import common.Array;
import common.Files;
import cnv.analysis.CentroidCompute;
import cnv.analysis.CentroidCompute.Builder;
import cnv.filesys.Project;
import cnv.qc.GcAdjustor.GCAdjustorBuilder;
import cnv.qc.GcAdjustor.GC_CORRECTION_METHOD;
import cnv.qc.GcAdjustor.GcModel;
import cnv.qc.GcAdjustorParameter;

public class GCcorrectionIterator {

	private static void iterate(Project proj, String outputRootDir, int numThreads) {

		new File(proj.PROJECT_DIRECTORY.getValue() + outputRootDir).mkdirs();
		int[] bpModels = new int[] { 50, 100, 250, 500, 1000, 2500, 5000, 10000, 1000000 };
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
		Project proj = new Project(args[0], false);
		iterate(proj, "gcCorrectionIterations/", 24);
		// generateAdjustmentParameters(proj, GC_CORRECTION_METHOD.values(), 4);
	}

}
