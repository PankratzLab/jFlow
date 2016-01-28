package seq.analysis;


import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Set;






import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;
import seq.manage.BamOps;
import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;

/**
 * Wrapper and parser/assembler for Somatic Sniper (https://github.com/genome/somatic-sniper)
 *
 */
public class SomaticSniper {

	private static void run(String bams, String vcfPop,String outputDir, Logger log) {
		VcfPopulation vpop = VcfPopulation.load(vcfPop, POPULATION_TYPE.TUMOR_NORMAL, log);
		String[] bamFiles = null;
		if (Files.isDirectory(bams)) {
			bamFiles = Files.list(bams, ".bam", false);
		} else {
			bamFiles = HashVec.loadFileToStringArray(bams, false, new int[] { 0 }, true);
		}
		log.reportTimeInfo("Detected " + bamFiles.length + " total bam files");
		TNSample[] tnSamples = matchSamples(bamFiles, outputDir, vpop, log);
		
	}

	private static TNSample[] matchSamples(String[] bamFiles, String outputDir, VcfPopulation vpop, Logger log) {
		Hashtable<String, String> all = new Hashtable<String, String>();
		for (int i = 0; i < bamFiles.length; i++) {
			all.put(BamOps.getSampleName(bamFiles[i]), bamFiles[i]);
		}
		ArrayList<String> analysisBams =new ArrayList<String>();
		ArrayList<TNSample> tnSamples = new ArrayList<TNSample>();
		for (String tnPair : vpop.getSubPop().keySet()) {
			Set<String> samps = vpop.getSubPop().get(tnPair);
			String tumor = null;
			String normal = null;
			for (String samp : samps) {
				if (vpop.getPopulationForInd(samp, RETRIEVE_TYPE.SUPER)[0].equals(VcfPopulation.TUMOR)) {
					tumor = samp;
				} else if (vpop.getPopulationForInd(samp, RETRIEVE_TYPE.SUPER)[0].equals(VcfPopulation.NORMAL)) {
					normal = samp;
				} else {
					throw new IllegalArgumentException("Unknown types");
				}
			}
			
			if (!all.containsKey(tumor) || !all.containsKey(normal)) {
				throw new IllegalArgumentException("Could not find bam file for Tumor " + tumor + " or for Normal " + normal);
			} else {
				TNSample tSample = new TNSample(all.get(normal), normal, all.get(tumor), tumor);
				tnSamples.add(tSample);
				analysisBams.add(all.get(normal));
				analysisBams.add(all.get(tumor));

			}
		}
		// Files.write(Array.toStringArray(analysisBams), outputDir + ext.rootOf(vpop.getFileName() + ".analysis.bams.txt"));
		log.reportTimeInfo("Matched " + tnSamples.size() + " tumor normals with bams ");
		return tnSamples.toArray(new TNSample[tnSamples.size()]);
	}

	private static class TNSample {
		private String normalBam;
		private String normalSample;
		private String tumorBam;
		private String tumorSample;

		public TNSample(String normalBam, String normalSample, String tumorBam, String tumorSample) {
			super();
			this.normalBam = normalBam;
			this.normalSample = normalSample;
			this.tumorBam = tumorBam;
			this.tumorSample = tumorSample;
		}

	}

	public static void main(String[] args) {
		String bams = "/scratch.global/lane0212/Project_Tsai_Project_028/151120_D00635_0090_BC81JNANXX/somaticSniper/allBams.txt";
		String vpop = "/scratch.global/lane0212/Project_Tsai_Project_028/151120_D00635_0090_BC81JNANXX/somaticSniper/TN.vpop";
		String outputDir = "/scratch.global/lane0212/Project_Tsai_Project_028/151120_D00635_0090_BC81JNANXX/somaticSniper/results/";
		new File(outputDir).mkdirs();
		Logger log = new Logger(outputDir + "somaticSniper.log");
		run(bams, vpop, outputDir, log);
	}

}
