package seq.analysis;

import filesys.Segment;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
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
import seq.analysis.GATK.MutectTumorNormal;
import seq.analysis.GATK_Genotyper.ANNOVCF;
import seq.manage.BamOps;
import seq.manage.VCFOps;
import seq.manage.VCFTumorNormalOps;
import seq.manage.VCOps;

/**
 * @author lane0212 For using the native Mutect in GATK 3.5+
 */
public class Mutect2 implements Producer<MutectTumorNormal> {

	private GATK gatk;
	private String[][] tumorNormalMatchedBams;
	private String pon;

	private String outputDir;
	private Logger log;
	private int index;

	/**
	 * @param tumorNormalMatchedBams
	 *            String[sampleCombo][0]=normal and String[sampleCombo][1]=tumor
	 */
	public Mutect2(GATK gatk, String[][] tumorNormalMatchedBams, String pon, String outputDir, int numSampleThreads, Logger log) {
		super();
		this.gatk = gatk;
		this.tumorNormalMatchedBams = tumorNormalMatchedBams;
		this.pon = pon;
		this.outputDir = outputDir;
		new File(outputDir).mkdirs();
		this.log = log;
		verify();
		this.index = 0;
	}

	private void verify() {
		for (int i = 0; i < tumorNormalMatchedBams.length; i++) {
			if (tumorNormalMatchedBams[i] == null || tumorNormalMatchedBams[i].length != 2) {
				throw new IllegalArgumentException("All bams must have a single normal and a single tumor");
			} else if (!Files.exists("", tumorNormalMatchedBams[i])) {
				throw new IllegalArgumentException("Missing one of " + Array.toStr(tumorNormalMatchedBams[i], "\n"));
			}
		}
	}

	@Override
	public boolean hasNext() {
		return index < tumorNormalMatchedBams.length;
	}

	@Override
	public Callable<MutectTumorNormal> next() {
		String normalBam = tumorNormalMatchedBams[index][0];
		String normalSample = BamOps.getSampleName(normalBam);
		String tumorBam = tumorNormalMatchedBams[index][1];
		String tumorSample = BamOps.getSampleName(tumorBam);
		index++;
		return new MutectTumorNormalWorker(gatk, normalBam, tumorBam, outputDir + normalSample + "_normal_" + tumorSample + "_tumor" + ".vcf.gz", pon, log);
	}

	@Override
	public void shutdown() {
	}

	private static class MutectTumorNormalWorker implements Callable<MutectTumorNormal> {
		private GATK gatk;
		private String normalBam;
		private String tumorBam;
		private String outputVCF;
		private String pon;
		private Logger log;

		public MutectTumorNormalWorker(GATK gatk, String normalBam, String tumorBam, String outputVCF, String pon, Logger log) {
			super();
			this.gatk = gatk;
			this.normalBam = normalBam;
			this.tumorBam = tumorBam;
			this.outputVCF = outputVCF;
			this.pon = pon;
			this.log = log;
		}

		@Override
		public MutectTumorNormal call() throws Exception {
			return gatk.callTumor(normalBam, tumorBam, outputVCF, pon, true, log);
		}
	}

	private static class PopulationOfNormals {
		private GATK gatk;
		private String[] samples;
		private Logger log;
		private String ponDir;

		private NormalSample[] normalSamples;

		private PopulationOfNormals(GATK gatk, String[] bamFilesFullPath, String ponDir, String[] samples, Logger log) {
			super();
			this.samples = samples;
			this.gatk = gatk;
			this.log = log;
			this.ponDir = ponDir;
			log.reportTimeInfo("Mapping normals to bam files ");
			mapNormals(bamFilesFullPath);
		}

		public NormalSample[] getNormalSamples() {
			return normalSamples;
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
					normalSamples[i] = new NormalSample(samples[i], map.get(samples[i]), ponDir + samples[i] + ".normal.vcf");
				}
			}
		}

		private void generatePON(int numThreads, int numSampleThreads) throws IllegalStateException {
			NormalProducer producer = new NormalProducer(gatk, normalSamples, ponDir, numSampleThreads, log);
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
		private Logger log;
		private int index;
		private int numSampleThreads;

		private NormalProducer(GATK gatk, NormalSample[] normalSamples, String outputDir, int numSampleThreads, Logger log) {
			super();
			this.gatk = gatk;
			this.normalSamples = normalSamples;
			this.numSampleThreads = numSampleThreads;
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
			Mutect2Worker worker = new Mutect2Worker(current, gatk, numSampleThreads, log);
			index++;
			return worker;
		}

		@Override
		public void shutdown() {
		}
	}

	private static class Mutect2Worker implements Callable<Mutect2Normal> {
		private NormalSample current;
		private GATK gatk;
		private int numSampleThreads;
		private Logger log;

		private Mutect2Worker(NormalSample current, GATK gatk, int numSampleThreads, Logger log) {
			super();
			this.current = current;
			this.gatk = gatk;
			this.numSampleThreads = numSampleThreads;
			this.log = log;
		}

		@Override
		public Mutect2Normal call() throws Exception {
			log.reportTimeInfo("Calling pon for sample " + current.getSample());
			return gatk.generateMutect2Normal(current.getBamFile(), current.getPonVCF(), numSampleThreads, log);
		}
	}

	private static class NormalSample {
		private String sample;
		private String bamFile;
		private String ponVCF;

		private NormalSample(String sample, String bamFile, String ponVCF) {
			super();
			this.sample = sample;
			this.bamFile = bamFile;
			this.ponVCF = ponVCF;
			if (!BamOps.getSampleName(bamFile).equals(sample)) {
				throw new IllegalArgumentException("Sample detected in bam file " + bamFile + " should have been " + sample + " but saw " + BamOps.getSampleName(bamFile) + "+ instead");
			}
		}

		private String getSample() {
			return sample;
		}

		private String getBamFile() {
			return bamFile;
		}

		private String getPonVCF() {
			return ponVCF;
		}

	}

	/**
	 * @param fileOfBams
	 * @param fileOftumorNormalMatchedBams
	 *            first column is normal, second is tumor
	 * @param outputDir
	 * @param ponVcf
	 * @param gatk
	 * @param type
	 * @param numThreads
	 * @param numSampleThreads
	 * @param log
	 * @throws IllegalStateException
	 */
	public static void run(String fileOfBams, String fileOftumorNormalMatchedBams, String outputDir, String ponVcf, GATK gatk, MUTECT_RUN_TYPES type, int numThreads, int numSampleThreads, Logger log) throws IllegalStateException {
		if (type == MUTECT_RUN_TYPES.GEN_NORMALS || type == MUTECT_RUN_TYPES.COMBINE_NORMALS) {
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
			PopulationOfNormals normals = new PopulationOfNormals(gatk, bams, outputDir, samples, log);

			try {
				normals.generatePON(numThreads, numSampleThreads);
			} catch (IllegalStateException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			if (type == MUTECT_RUN_TYPES.COMBINE_NORMALS) {
				ArrayList<String> ponVcfs = new ArrayList<String>();
				for (int i = 0; i < normals.getNormalSamples().length; i++) {
					if (Files.exists(normals.getNormalSamples()[i].getPonVCF())) {
						ponVcfs.add(normals.getNormalSamples()[i].getPonVCF());
					} else {
						throw new IllegalStateException("Missing pon vcf " + normals.getNormalSamples()[i].getPonVCF());
					}
				}
				gatk.combinePonVcfs(Array.toStringArray(ponVcfs), ponVcf, 2, log);
			}
		} else if (type == MUTECT_RUN_TYPES.CALL_SOMATIC) {
			callSomatic(fileOftumorNormalMatchedBams, outputDir, ponVcf, gatk, null, null, null, numThreads, numSampleThreads, true, log);
		}
	}

	public static MutectTumorNormal[] callSomatic(String fileOftumorNormalMatchedBams, String outputDir, String ponVcf, GATK gatk, ANNOVCF annoVCF, String finalMergeWithVCF, String tparams, int numThreads, int numSampleThreads, boolean merge, Logger log) {
		if (fileOftumorNormalMatchedBams == null || !Files.exists(fileOftumorNormalMatchedBams)) {
			throw new IllegalArgumentException("Missing file " + fileOftumorNormalMatchedBams);
		}
		if (!merge) {
			log.reportTimeInfo("Will not merge results as strict tumor normal comparison");
		}
		String[][] tumorNormalMatchedBams = HashVec.loadFileToStringMatrix(fileOftumorNormalMatchedBams, false, new int[] { 0, 1 }, false);
		Mutect2 mutect2 = new Mutect2(gatk, tumorNormalMatchedBams, ponVcf, outputDir, numSampleThreads, log);
		WorkerTrain<MutectTumorNormal> train = new WorkerTrain<GATK.MutectTumorNormal>(mutect2, numThreads, 2, log);
		ArrayList<MutectTumorNormal> results = new ArrayList<GATK.MutectTumorNormal>();
		ArrayList<String> finalTNnoFiltVCfs = new ArrayList<String>();
		ArrayList<String> finalTNFiltVCFS = new ArrayList<String>();

		while (train.hasNext()) {
			MutectTumorNormal tmp = train.next();
			results.add(tmp);
			finalTNnoFiltVCfs.add(tmp.getReNamedOutputVCF());
			finalTNFiltVCFS.add(tmp.getReNamedFilteredVCF());

		}

		if (merge) {
			String rootNoFilter = outputDir + ext.rootOf(fileOftumorNormalMatchedBams) + ".merged";
			mergeAndAnnotate(outputDir, gatk, annoVCF, finalMergeWithVCF, numThreads, log, tumorNormalMatchedBams, finalTNnoFiltVCfs, false, rootNoFilter);

			String rootFilter = outputDir + ext.rootOf(fileOftumorNormalMatchedBams) + ".merged.filtered";
			String filtAnno = mergeAndAnnotate(outputDir, gatk, annoVCF, finalMergeWithVCF, numThreads, log, tumorNormalMatchedBams, finalTNFiltVCFS, true, rootFilter);

			if (tparams != null) {
				runTally(tparams, log, filtAnno);
			}
		}
		return results.toArray(new MutectTumorNormal[results.size()]);
	}

	private static void runTally(String tparams, Logger log, String vcf) {
		log.reportTimeInfo("Loading vcf tally params from " + tparams);
		log.reportTimeWarning("Relying on specific file format... you will likely need to see the code to get this to work");
		String[] params = HashVec.loadFileToStringArray(tparams, false, new int[] { 0 }, true);
		String vpops = null;
		String popDir = null;
		String omim = null;
		String[] extras = null;
		boolean controlSpecifiComp = true;
		double[] mafs = null;
		String genesetDir = null;

		for (int i = 0; i < params.length; i++) {
			if (params[i].startsWith("vpops=")) {
				vpops = ext.parseStringArg(params[i], null);
			} else if (params[i].startsWith("popDir=")) {
				popDir = ext.parseStringArg(params[i], null);
			} else if (params[i].startsWith("omim=")) {
				omim = ext.parseStringArg(params[i], null);
			} else if (params[i].startsWith("extra=")) {
				extras = ext.parseStringArg(params[i], null).split(",");
			} else if (params[i].startsWith("controlSpecifiComp=")) {
				controlSpecifiComp = ext.parseBooleanArg(params[i]);
			} else if (params[i].startsWith("genesetDir=")) {
				genesetDir = ext.parseStringArg(params[i], null);
			} else if (params[i].startsWith("maf=")) {
				mafs = Array.toDoubleArray(ext.parseStringArg(params[i], null).split(","));
			}
		}

		String annoCp = popDir + VCFOps.getAppropriateRoot(vcf, true) + ".vcf.gz";
		Files.copyFileUsingFileChannels(vcf, annoCp, log);
		Files.copyFileUsingFileChannels(vcf + ".tbi", annoCp + ".tbi", log);
		popDir = popDir + "freq_" + VCFOps.getAppropriateRoot(vcf, true);
		new File(popDir).mkdirs();
		Files.copyFile(vpops, popDir + ext.removeDirectoryInfo(vpops));

		System.exit(1);
		for (int i = 0; i < mafs.length; i++) {
			VCFSimpleTally.test(vcf, popDir, new String[] { ext.removeDirectoryInfo(vpops) }, omim, extras, genesetDir, mafs[i], controlSpecifiComp);
		}
	}

	private static String mergeAndAnnotate(String outputDir, GATK gatk, ANNOVCF annoVCF, String finalMergeWithVCF, int numThreads, Logger log, String[][] tumorNormalMatchedBams, ArrayList<String> finalTNVCfs, boolean extract, String root) {
		String outMergeVCF = root + ".vcf.gz";
		String outMergeRenameVCF = root + ".renamed.vcf.gz";
		// String outMergeRenameAnnoVCF = root + ".renamed.anno.vcf.gz";
		String vcfToReturn = null;
		if (!Files.exists(outMergeVCF)) {
			gatk.mergeVCFs(Array.toStringArray(finalTNVCfs), outMergeVCF, numThreads, false, log);
		}

		if (!Files.exists(outMergeRenameVCF)) {
			VCFTumorNormalOps.renameMergeVCF(outMergeVCF, outMergeRenameVCF);
		}
		log.reportTimeWarning("Attempting to annotate using default locations...");
		String finalAnno = GATK_Genotyper.annotateOnlyWithDefualtLocations(outMergeRenameVCF, annoVCF, true, false, log);

		if (finalMergeWithVCF != null) {
			String finalMerge = VCFOps.getAppropriateRoot(finalAnno, false) + ".finalMerge.vcf.gz";
			log.reportTimeInfo("Merging with " + finalMergeWithVCF);
			gatk.mergeVCFs(new String[] { outMergeRenameVCF, finalMergeWithVCF }, finalMerge, numThreads, false, log);
			vcfToReturn = finalMerge;
		} else {
			vcfToReturn = finalAnno;
		}

		if (extract) {
			ArrayList<String> bamsToExtract = new ArrayList<String>();
			for (int i = 0; i < tumorNormalMatchedBams.length; i++) {
				for (int j = 0; j < tumorNormalMatchedBams[i].length; j++) {
					bamsToExtract.add(tumorNormalMatchedBams[i][j]);
				}
			}
			String extractFiltDir = outputDir + "extract_" + VCFOps.getAppropriateRoot(finalAnno, true);
			extractBamsTo(extractFiltDir, Array.toStringArray(bamsToExtract), finalAnno, numThreads, log);
		}

		return vcfToReturn;

	}

	private static void extractBamsTo(String extractDir, String[] bams, String vcf, int numThreads, Logger log) {
		new File(extractDir).mkdirs();
		String segFile = extractDir + "segments.txt";
		String bamFile = extractDir + "bams.txt";
		Files.writeList(bams, bamFile);
		if (!Files.exists(segFile)) {
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(segFile));
				VCFFileReader reader = new VCFFileReader(vcf, true);
				for (VariantContext vc : reader) {
					Segment seg = VCOps.getSegment(vc);
					writer.println(seg.getChr() + "\t" + seg.getStart() + "\t" + seg.getStop());
				}
				writer.close();
				reader.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + segFile);
				log.reportException(e);
			}

		}
		VCFOps.extractSegments(vcf, segFile, 300, bamFile, extractDir, false, true, false, false, null, numThreads, log);
	}

	private static void batchPON(int numNormalBatches, GATK gatk, String bamFilesFullPath, String outputDir, int numthreads, int numSampleThreads, Logger log) {
		ArrayList<String[]> splits = Array.splitUpArray(HashVec.loadFileToStringArray(bamFilesFullPath, false, new int[] { 0 }, true), numNormalBatches, log);
		ArrayList<String> command = new ArrayList<String>();
		String[][] batches = new String[splits.size()][1];
		String baseOut = "[%0]";
		for (int i = 0; i < batches.length; i++) {
			batches[i][0] = outputDir + "batch_" + i + "_" + "pon.txt";
			Files.writeList(splits.get(i), batches[i][0]);
		}
		getJava(command);
		command.addAll(getBaseArgs(gatk, outputDir, numthreads, numSampleThreads));
		command.add("normalBams=" + baseOut);

		Files.qsub(outputDir + "Contam.txt", Array.toStr(Array.toStringArray(command), " "), batches, 62000, 40, numthreads * numSampleThreads, "small");
	}

	private static void getJava(ArrayList<String> command) {
		command.add("java");
		command.add("-Xmx62000m");
		command.add("-cp");
		command.add("~/parkMutect.jar");
		command.add("seq.analysis.Mutect2");

	}

	private static ArrayList<String> getBaseArgs(GATK gatk, String outputDir, int numthreads, int numSampleThreads) {
		ArrayList<String> base = new ArrayList<String>();
		base.add("gatk=" + gatk.getGATKLocation());
		base.add("ref=" + gatk.getReferenceGenomeFasta());
		base.add("knownSnps=" + gatk.getDbSnpKnownSites());
		base.add("cosmic=" + gatk.getCosmicKnownSites());
		base.add("regions=" + gatk.getRegionsFile());
		base.add("outputDir=" + outputDir);
		base.add("numthreads=" + numthreads);
		base.add("numSampleThreads=" + numSampleThreads);
		return base;
	}

	public enum MUTECT_RUN_TYPES {
		/**
		 * Generate individual normal calls
		 */
		GEN_NORMALS, /**
		 * Combine individual normal calls, will generate normals if needed
		 */
		COMBINE_NORMALS, /**
		 * Call somatic variants, will combine variants and generate normals if needed
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
		String ponVCF = null;
		int numthreads = 1;
		int numSampleThreads = 4;
		int numbatches = 0;
		String outputDir = "mutect/";
		String tumorNormalBams = null;
		MUTECT_RUN_TYPES run = MUTECT_RUN_TYPES.GEN_NORMALS;
		String finalMergeVCF = null;
		ANNOVCF annoVCF = null;
		String tparams = null;
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
		usage += "   (11) pon vcf (i.e. ponVcf= (no default))\n" + "";
		usage += "   (12) number of threads per sample (i.e. numSampleThreads=" + numSampleThreads + " (default))\n" + "";
		usage += "   (13) full path to a file of tumor-normal matched (tab-delimited) bam files, normal in first column, tumor in second (i.e. tumorNormalBams= (no default))\n" + "";
		usage += "   (14) full path to a final merge vcf to merge somatic calls with (i.e. finalVCF= (no default))\n" + "";
		usage += "   (15) full path to a vcf tallyparamFile (i.e. tparams= (no default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("normalBams=")) {
				fileOfNormalBams = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("finalVCF=")) {
				finalMergeVCF = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("tparams=")) {
				tparams = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("ref=")) {
				referenceGenomeFasta = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("tumorNormalBams=")) {
				tumorNormalBams = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("run=")) {
				run = MUTECT_RUN_TYPES.valueOf(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("gatk=")) {
				gatkLocation = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("ponVCF=")) {
				ponVCF = args[i].split("=")[1];
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
			} else if (args[i].startsWith("numSampleThreads=")) {
				numSampleThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith(GATK_Genotyper.EXTRA_VCF_ANNOTATIONS)) {
				annoVCF = ANNOVCF.fromArg(args[i]);
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
		case COMBINE_NORMALS:
		case CALL_SOMATIC:
			callSomatic(tumorNormalBams, outputDir, ponVCF, gatk, annoVCF, finalMergeVCF, tparams, numthreads, numSampleThreads, true, log);

			break;
		case GEN_NORMALS:

			if (numbatches > 0) {
				String outputNormal = outputDir + "pon/";
				new File(outputNormal).mkdirs();
				log.reportTimeInfo("Generating normal sample vcfs batches in " + outputNormal);
				batchPON(numbatches, gatk, fileOfNormalBams, outputNormal, numthreads, numSampleThreads, log);
			} else if (fileOfNormalBams != "") {
				log.reportTimeInfo("Generating normal sample vcfs");
				try {
					run(fileOfNormalBams, tumorNormalBams, outputDir, ponVCF, gatk, run, numthreads, numSampleThreads, log);
				} catch (IllegalStateException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
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
