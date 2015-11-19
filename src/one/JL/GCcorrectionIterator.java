package one.JL;

import java.io.File;

import cnv.filesys.Project;

public class GCcorrectionIterator {

	private static void iterate(Project proj, String outputRootDir) {
		new File(outputRootDir).mkdirs();
		int[] bpModels = new int[] { 50, 100, 250, 500, 1000, 2500, 5000, 10000 };
		
		
		for (int i = 0; i < bpModels.length; i++) {

		}
	}
	
	public static void main(String[] args) {
		Project proj = new Project("C:/workspace/Genvisis/projects/OSv2_hg19.properties", false);
		iterate(proj, proj.PROJECT_DIRECTORY.getValue()+"gcCorrectionIterations/");
		
		//generateAdjustmentParameters(proj, GC_CORRECTION_METHOD.values(), 4);
	}

}
