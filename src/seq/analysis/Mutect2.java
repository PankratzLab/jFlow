package seq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.concurrent.Callable;

import javax.jms.IllegalStateException;

import common.Files;
import common.HashVec;
import common.Logger;
import common.WorkerTrain;
import common.WorkerTrain.Producer;
import common.ext;
import seq.analysis.GATK.Mutect2Normal;
import seq.manage.BamOps;
import seq.manage.VCFOps;
import seq.manage.VCOps;

/**
 * @author lane0212 For using the native Mutect in GATK 3.5+
 */
public class Mutect2 {

	private static class PopulationOfNormals {
		private GATK gatk;
		private String[] samples;
		private Logger log;

		private NormalSample[] normalSamples;

		private PopulationOfNormals(GATK gatk, String[] bamFilesFullPath, String[] samples, Logger log) {
			super();
			this.samples = samples;
			this.gatk = gatk;
			this.log = log;
			log.reportTimeInfo("Mapping normals to bam files ");
			mapNormals(bamFilesFullPath);
		}

		private void mapNormals(String[] bamFilesFullPath) {
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

		private void generatePON(String outputNormalDir, String ponVCF, int numThreads) throws IllegalStateException {
			NormalProducer producer = new NormalProducer(gatk, normalSamples, outputNormalDir, log);
			WorkerTrain<Mutect2Normal> train = new WorkerTrain<GATK.Mutect2Normal>(producer, numThreads, 2, log);
			ArrayList<String> vcfsToCombine = new ArrayList<String>();

			while (train.hasNext()) {
				Mutect2Normal validate = train.next();
				if (validate.isFail() || !Files.exists(validate.getOutputVCF())) {
					throw new IllegalStateException("PON generation failed");
				} else {
					vcfsToCombine.add(validate.getOutputVCF());
				}
			}

		}

	}

	private static class NormalProducer implements Producer<Mutect2Normal> {
		private GATK gatk;
		private NormalSample[] normalSamples;
		private String outputDir;
		private Logger log;
		private int index;

		public NormalProducer(GATK gatk, NormalSample[] normalSamples, String outputDir, Logger log) {
			super();
			this.gatk = gatk;
			this.normalSamples = normalSamples;
			this.outputDir = outputDir;

			this.log = log;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < normalSamples.length;
		}

		@Override
		public Callable<Mutect2Normal> next() {
			final NormalSample current = normalSamples[index];
			Mutect2Worker worker = new Mutect2Worker(current, outputDir, gatk, log);
			index++;
			return worker;
		}

		@Override
		public void shutdown() {
		}
	}

	private static class Mutect2Worker implements Callable<Mutect2Normal> {
		private NormalSample current;
		private String outputDir;
		private GATK gatk;
		private Logger log;

		public Mutect2Worker(NormalSample current, String outputDir, GATK gatk, Logger log) {
			super();
			this.current = current;
			this.outputDir = outputDir;
			this.gatk = gatk;
			this.log = log;
		}

		@Override
		public Mutect2Normal call() throws Exception {
			String outputVCF = outputDir + current.getSample() + ".normal.vcf";
			return gatk.generateMutect2Normal(current.getBamFile(), outputVCF, 8, log);

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

		public String getBamFile() {
			return bamFile;
		}

		public String getSample() {
			return sample;
		}
	}

	public static void run(String fileOfBams, String outputDir, GATK gatk, int numThreads, Logger log) {
		String[] bams = HashVec.loadFileToStringArray(fileOfBams, false, new int[] { 0 }, true);
		String[] samples = new String[bams.length];
		for (int i = 0; i < samples.length; i++) {
			try {
				samples[i] = BamOps.getSampleName(bams[i]);
			} catch (Exception e) {
				log.reportTimeError("Could not find sample for " + bams[i]);
				return;
			}
		}
		PopulationOfNormals normals = new PopulationOfNormals(gatk, bams, samples, log);
		String ponDir = outputDir + "pon/";
		new File(ponDir).mkdirs();
		String ponVCF = ponDir + "pon.vcf";
		try {
			normals.generatePON(ponDir, ponVCF, numThreads);
		} catch (IllegalStateException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String fileOfBams = "fileOfPonBams.txt";
		String referenceGenomeFasta = "/home/pankrat2/public/bin/ref/hg19_canonical.fa";
		String knownSnps = "/home/tsaim/lane0212/bin/ref/dbsnp_138.hg19.vcf";
		String gatkLocation = "/home/tsaim/lane0212/bin/GATK_3_5/";
		String bams = "/scratch.global/lane0212/Project_Tsai_Project_028/151120_D00635_0090_BC81JNANXX/mutect/ponBams.txt";
		String cosmic = "/home/pankrat2/public/bin/ref/COSMIC/b37_cosmic_v54_120711.hg19_chr.vcf";
		String regions = "/home/pankrat2/public/bin/ref/AgilentCaptureRegions.bpBuffer100.bed";
		int numthreads = 1;
		String ouptputDir = "/scratch.global/lane0212/Project_Tsai_Project_028/151120_D00635_0090_BC81JNANXX/mutect/";
		new File(ouptputDir).mkdirs();
		Logger log = new Logger(ouptputDir + "TN.log");

		String usage = "\n" + "seq.analysis.Mutect2 requires 0-1 arguments\n";
		usage += "   (1) full path to a file of bams (i.e. bams=" + fileOfBams + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("bams=")) {
				fileOfBams = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		GATK gatk = new GATK(gatkLocation, referenceGenomeFasta, knownSnps, regions, cosmic, true, false, true, log);
		run(bams, ouptputDir, gatk, numthreads, log);

	}
}
// VCFOps.b_ToHg(cosmic, ext.addToRoot(cosmic, ".hg19_chr"), referenceGenomeFasta, log);
// System.exit(1);

// int buffer = 100;
// BEDFileReader reader = new BEDFileReader(regions, false);
// LocusSet<BEDFeatureSeg> segs = reader.loadAll(new Logger());
// reader.close();
// String buffered = ext.addToRoot(regions, ".bpBuffer" + buffer);
// // segs.writeRegions(buffered, TO_STRING_TYPE.REGULAR, false, new Logger());
// try {
// PrintWriter writer = new PrintWriter(new FileWriter(buffered));
// for (int i = 0; i < segs.getLoci().length; i++) {
// Segment buff = segs.getLoci()[i].getBufferedSegment(buffer);
// writer.println(buff.getChromosomeUCSC() + "\t" + buff.getStart() + "\t" + buff.getStop());
// }
// writer.close();
// } catch (Exception e) {
// new Logger().reportError("Error writing to " + buffered);
// new Logger().reportException(e);
// }
// System.exit(1);
