package seq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Set;
import java.util.concurrent.Callable;

import common.Array;
import common.CmdLine;
import common.Files;
import common.HashVec;
import common.Logger;
import common.WorkerTrain;
import common.ext;
import common.WorkerTrain.Producer;
import seq.manage.BamOps;
import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;

/**
 * Wrapper and parser/assembler for Somatic Sniper (https://github.com/genome/somatic-sniper)
 *
 */
public class SomaticSniper {

	private static void run(String bams, String vcfPop, String outputDir, String referenceGenome, String somaticSnipLoc, int numThreads, Logger log) {
		VcfPopulation vpop = VcfPopulation.load(vcfPop, POPULATION_TYPE.TUMOR_NORMAL, log);
		String[] bamFiles = null;
		if (Files.isDirectory(bams)) {
			bamFiles = Files.list(bams, ".bam", false);
		} else {
			bamFiles = HashVec.loadFileToStringArray(bams, false, new int[] { 0 }, true);
		}
		log.reportTimeInfo("Detected " + bamFiles.length + " total bam files");

		SomaticParams params = new SomaticParams("", referenceGenome, outputDir, true);
		TNSample[] tnSamples = matchSamples(bamFiles, outputDir, vpop, params, log);

		TNProducer producer = new TNProducer(tnSamples);
		WorkerTrain<Boolean> train = new WorkerTrain<Boolean>(producer, numThreads, 2, log);
		while (train.hasNext()) {
			train.next();
		}

	}

	private static class TNProducer implements Producer<Boolean> {
		private TNSample[] tnSamples;
		private int index;

		public TNProducer(TNSample[] tnSamples) {
			super();
			this.tnSamples = tnSamples;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < tnSamples.length;
		}

		@Override
		public Callable<Boolean> next() {
			TNSample current = tnSamples[index];
			index++;

			return current;
		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}

	}

	private static TNSample[] matchSamples(String[] bamFiles, String outputDir, VcfPopulation vpop, SomaticParams somaticParams, Logger log) {
		Hashtable<String, String> all = new Hashtable<String, String>();
		for (int i = 0; i < bamFiles.length; i++) {
			all.put(BamOps.getSampleName(bamFiles[i]), bamFiles[i]);
		}
		ArrayList<String> analysisBams = new ArrayList<String>();
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
				TNSample tSample = new TNSample(all.get(normal), normal, all.get(tumor), tumor, somaticParams, log);
				tnSamples.add(tSample);
				analysisBams.add(all.get(normal));
				analysisBams.add(all.get(tumor));

			}
		}
		if (analysisBams.size() < bamFiles.length) {
			Files.writeList(Array.toStringArray(analysisBams), outputDir + ext.rootOf(vpop.getFileName() + ".analysis.bams.txt"));
		}
		log.reportTimeInfo("Matched " + tnSamples.size() + " tumor normals with bams ");
		return tnSamples.toArray(new TNSample[tnSamples.size()]);
	}

	private static class SomaticParams {
		private String somaticSnipeLoc;
		private String refGenome;
		private String outputDir;
		private boolean overwriteExisting;

		public SomaticParams(String somaticSnipeLoc, String refGenome, String outputDir, boolean overwriteExisting) {
			super();
			this.somaticSnipeLoc = somaticSnipeLoc;
			this.refGenome = refGenome;
			this.outputDir = outputDir;
			this.overwriteExisting = overwriteExisting;

		}

		public String getOutputDir() {
			return outputDir;
		}

		public boolean isOverwriteExisting() {
			return overwriteExisting;
		}

		public String getRefGenome() {
			return refGenome;
		}

		public String getSnipe() {
			return somaticSnipeLoc + "bam-somaticsniper";
		}

		private String[] generateCommand(String normalBam, String tumorBam, String output) {
			ArrayList<String> command = new ArrayList<String>();
			command.add(getSnipe());
			command.add("-Q");
			command.add("40");
			command.add("-G");
			command.add("-L");
			command.add("-f");
			command.add(refGenome);
			command.add(tumorBam);
			command.add(normalBam);
			command.add(output);
			return Array.toStringArray(command);
		}
	}

	private static class TNSample implements Callable<Boolean> {
		private String normalBam;
		private String normalSample;
		private String tumorBam;
		private String tumorSample;
		private SomaticParams somaticParams;
		private Logger log;

		public TNSample(String normalBam, String normalSample, String tumorBam, String tumorSample, SomaticParams somaticParams, Logger log) {
			super();
			this.normalBam = normalBam;
			this.normalSample = normalSample;
			this.tumorBam = tumorBam;
			this.tumorSample = tumorSample;
			this.somaticParams = somaticParams;
			this.log = log;
		}

		@Override
		public Boolean call() throws Exception {
			String[] inputs = new String[] { somaticParams.getRefGenome(), somaticParams.getSnipe(), tumorBam, normalBam };
			String output = somaticParams.getOutputDir() + "T_" + tumorSample + "_N_" + normalSample + ".vcf";
			String[] outputs = new String[] { output };
			String[] command = somaticParams.generateCommand(normalBam, tumorBam, output);
			return CmdLine.runCommandWithFileChecks(command, "", inputs, outputs, true, somaticParams.isOverwriteExisting(), false, log);
		}
	}

	public static void main(String[] args) {
		String bams = "/scratch.global/lane0212/Project_Tsai_Project_028/151120_D00635_0090_BC81JNANXX/somaticSniper/results/TN.vpop.analysis.bams";
		String vpop = "/scratch.global/lane0212/Project_Tsai_Project_028/151120_D00635_0090_BC81JNANXX/somaticSniper/TN.vpop";
		String outputDir = "/scratch.global/lane0212/Project_Tsai_Project_028/151120_D00635_0090_BC81JNANXX/somaticSniper/results/";
		String refGenome = "/home/pankrat2/public/bin/ref/hg19_canonical.fa";
		String somaticSnipLoc = "/home/pankrat2/public/bin/somaticSniper/somatic-sniper/build/bin/";
		int numThreads = 1;
		new File(outputDir).mkdirs();
		Logger log = new Logger(outputDir + "somaticSniper.log");
		run(bams, vpop, outputDir, refGenome, somaticSnipLoc, numThreads, log);

	}

}
