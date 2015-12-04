package seq.analysis;

import java.util.Hashtable;
import java.util.concurrent.Callable;

import common.Logger;
import common.WorkerTrain.Producer;
import seq.analysis.GATK.Mutect2Normal;
import seq.manage.BamOps;

/**
 * @author lane0212 For using the native Mutect in GATK 3.5+
 */
public class Mutect2 {

	private static class PopulationOfNormals {
		private GATK gatk;
		private String[] bamFilesFullPath;
		private String[] samples;
		private Logger log;

		private NormalSample[] normalSamples;

		private PopulationOfNormals(String[] bamFilesFullPath, String[] samples, Logger log) {
			super();
			this.bamFilesFullPath = bamFilesFullPath;
			this.samples = samples;
			if (samples.length != bamFilesFullPath.length) {
				throw new IllegalArgumentException("Sample array must be same length as bam file array");
			}
			this.log = log;
			log.reportTimeInfo("Mapping normals to bam files ");
			mapNormals();
		}

		private void mapNormals() {
			this.normalSamples = new NormalSample[samples.length];
			Hashtable<String, String> map = new Hashtable<String, String>();
			for (int i = 0; i < bamFilesFullPath.length; i++) {
				map.put(BamOps.getSampleName(bamFilesFullPath[i]), bamFilesFullPath[i]);
			}

			for (int i = 0; i < samples.length; i++) {
				if (!map.containsKey(samples[i])) {
					throw new IllegalArgumentException("Could not determine bam map for " + samples[i]);
				} else {
					normalSamples[i] = new NormalSample(samples[i], map.get(samples[i]));
				}
			}
		}

		private void generatePON(String outputDir, String rootponVCF, int numThreads) {
			String finalVcf = outputDir + rootponVCF;
			for (int i = 0; i < normalSamples.length; i++) {
				String normalPon = outputDir + samples[i] + ".normal.vcf";
			}
		}

	}

	private static class NormalProducer implements Producer<Mutect2Normal> {

		@Override
		public boolean hasNext() {
			// TODO Auto-generated method stub
			return false;
		}

		@Override
		public Callable<Mutect2Normal> next() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}

	}

	private static class NormalSample {
		private String sample;
		private String bamFile;

		private NormalSample(String sample, String bamFile) {
			super();
			this.sample = sample;
			this.bamFile = bamFile;
			if (!BamOps.getSampleName(bamFile).equals(sample)) {
				throw new IllegalArgumentException("Sample detected in bam file " + bamFile + " shoul have been " + sample + " but saw " + BamOps.getSampleName(bamFile) + "+ instead");

			}
		}

	}

}
