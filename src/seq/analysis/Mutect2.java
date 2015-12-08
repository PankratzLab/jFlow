package seq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.concurrent.Callable;

import javax.jms.IllegalStateException;

import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.WorkerTrain;
import common.WorkerTrain.Producer;
import common.ext;
import seq.analysis.GATK.Mutect2Normal;
import seq.manage.BamOps;

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

		private void generatePON(String outputNormalDir, int numThreads) throws IllegalStateException {
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
			return gatk.generateMutect2Normal(current.getBamFile(), outputVCF, 1, log);
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

		try {
			normals.generatePON(outputDir, numThreads);
		} catch (IllegalStateException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static void batchPON(int numNormalBatches, GATK gatk, String bamFilesFullPath, String outputDir, int numthreads, Logger log) {
		ArrayList<String[]> splits = Array.splitUpArray(HashVec.loadFileToStringArray(bamFilesFullPath, false, new int[] { 0 }, true), numNormalBatches, log);
		ArrayList<String> command = new ArrayList<String>();
		String[][] batches = new String[splits.size()][1];
		String baseOut = "[%0]";
		for (int i = 0; i < batches.length; i++) {
			batches[i][0] = outputDir + "batch_" + i + "_" + "pon.txt";
			Files.writeList(splits.get(i), batches[i][0]);
		}
		getJava(command);
		command.addAll(getBaseArgs(gatk, outputDir, numthreads));
		command.add("normalBams=" + baseOut);

		Files.qsub(outputDir + "Contam.txt", Array.toStr(Array.toStringArray(command), " "), batches, 255000, 120, numthreads, "pankratz");
	}

	private static void getJava(ArrayList<String> command) {
		command.add("java");
		command.add("-Xmx256000m");

		command.add("-cp");
		command.add("~/parkMutect.jar");
		command.add("seq.analysis.Mutect2");

	}

	private static ArrayList<String> getBaseArgs(GATK gatk, String outputDir, int numthreads) {
		ArrayList<String> base = new ArrayList<String>();
		base.add("gatk=" + gatk.getGATKLocation());
		base.add("ref=" + gatk.getReferenceGenomeFasta());
		base.add("knownSnps=" + gatk.getDbSnpKnownSites());
		base.add("cosmic=" + gatk.getCosmicKnownSites());
		base.add("regions=" + gatk.getRegionsFile());
		base.add("outputDir=" + outputDir);
		base.add("numthreads=" + numthreads);
		return base;
	}

	public enum MUTECT_RUN_TYPES {
		/**
		 * Generate individual normal calls
		 */
		GEN_NORMALS, /**
		 * Combine individual normal calls
		 */
		COMBINE_NORMALS, /**
		 * Call somatic variants
		 */
		CALL_SOMATIC;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String fileOfNormalBams = "";
		String referenceGenomeFasta = "hg19_canonical.fa";
		String knownSnps = "dbsnp_138.hg19.vcf";
		String gatkLocation = "GATK_3_5/";
		String cosmic = "b37_cosmic_v54_120711.hg19_chr.vcf";
		String regions = "AgilentCaptureRegions.bed";
		int numthreads = 1;
		int numbatches = 0;
		String outputDir = "mutect/";
		MUTECT_RUN_TYPES run = MUTECT_RUN_TYPES.GEN_NORMALS;

		String usage = "\n" + "seq.analysis.Mutect2 requires 0-1 arguments\n";
		usage += "   (1) full path to a file of normal sample bams (i.e. normalBams=" + fileOfNormalBams + " (default))\n" + "";
		usage += "   (2) full path to a reference genome (i.e. ref=" + referenceGenomeFasta + " (default))\n" + "";
		usage += "   (3) known dbsnp snps (i.e. knownSnps=" + knownSnps + " (default))\n" + "";
		usage += "   (4) cosmic snps (i.e. cosmic=" + cosmic + " (default))\n" + "";
		usage += "   (5) regions for calling (i.e. regions=" + regions + " (default))\n" + "";
		usage += "   (6) output root directory (i.e. outputDir=" + outputDir + " (default))\n" + "";
		usage += "   (7) number of batches- triggers qsub generation (i.e. numNormalBatches=" + numbatches + " (default))\n" + "";
		usage += "   (8) number of threads (i.e. numthreads=" + numthreads + " (default))\n" + "";
		usage += "   (9) gatk directory (i.e. gatk=" + gatkLocation + " (default))\n" + "";
		usage += "   (10) type of analysis (i.e. run=" + run + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("normalBams=")) {
				fileOfNormalBams = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("ref=")) {
				referenceGenomeFasta = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("run=")) {
				run = MUTECT_RUN_TYPES.valueOf(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("gatk=")) {
				gatkLocation = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("knownSnps=")) {
				knownSnps = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("cosmic=")) {
				cosmic = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("regions=")) {
				regions = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("outputDir=")) {
				outputDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("numNormalBatches=")) {
				numbatches = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("numthreads=")) {
				numthreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		new File(outputDir).mkdirs();
		Logger log = new Logger(outputDir + "TN.log");
		GATK gatk = new GATK(gatkLocation, referenceGenomeFasta, knownSnps, regions, cosmic, true, false, true, log);
		switch (run) {
		case CALL_SOMATIC:
			break;
		case COMBINE_NORMALS:
			break;
		case GEN_NORMALS:

			if (numbatches > 0) {
				String outputNormal = outputDir + "pon/";
				new File(outputNormal).mkdirs();
				log.reportTimeInfo("Generating normal sample vcfs batches in " + outputNormal);
				batchPON(numbatches, gatk, fileOfNormalBams, outputNormal, numthreads, log);
			} else if (fileOfNormalBams != "") {
				log.reportTimeInfo("Generating normal sample vcfs");
				run(fileOfNormalBams, outputDir, gatk, numthreads, log);
			}
			break;
		default:

			break;

		}

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
