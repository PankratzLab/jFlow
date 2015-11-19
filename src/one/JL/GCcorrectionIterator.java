package one.JL;

import java.io.File;

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
		int[] bpModels = new int[] { 50, 100, 250, 500, 1000, 2500, 5000, 10000 };
		int[] regressDistance = new int[] { 1000, 10000, 100000, 1000000 };
		int[] snpMAD = new int[] { 0, 1, 2, 5, 10, 15 };
		String freshCents = proj.PROJECT_DIRECTORY.getValue() + outputRootDir + "freshCents.cent";
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
			for (int j = 0; j < regressDistance.length; j++) {
				for (int j2 = 0; j2 < snpMAD.length; j2++) {
					String root = outputRootDir + "gcmodel_bp_" + bpModels[i] + "_regress_" + regressDistance[j] + "_snpMad_" + snpMAD[j2];
					builder.regressionDistance(regressDistance[j]);
					builder.numSnpMAD(snpMAD[j2]);
					String[][] generated = GcAdjustorParameter.generateAdjustmentParameters(proj, builder, new String[] { freshCents }, new GC_CORRECTION_METHOD[] { GC_CORRECTION_METHOD.GENVISIS_GC }, gcModel, root, numThreads);

				}
			}
		}
	}

	public static void main(String[] args) {
		Project proj = new Project(args[0], false);
		iterate(proj, "gcCorrectionIterations/", 1);
		// generateAdjustmentParameters(proj, GC_CORRECTION_METHOD.values(), 4);
	}

}
