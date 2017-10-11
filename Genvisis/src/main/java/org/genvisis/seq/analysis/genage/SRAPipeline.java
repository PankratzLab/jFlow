package org.genvisis.seq.analysis.genage;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

import org.genvisis.CLI;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.qsub.Qsub;
import org.genvisis.seq.SeqVariables.ASSAY_TYPE;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.seq.SeqVariables.PLATFORM;
import org.genvisis.seq.analysis.genage.Pipeline.PIPELINE_PARTS;
import org.genvisis.seq.analysis.genage.Pipeline.PipelinePart;
import org.genvisis.seq.manage.BamImport;
import org.genvisis.seq.manage.BamSample.NORMALIZATON_METHOD;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.sra.SRARunTable;
import org.genvisis.sra.SRASample;
import org.genvisis.sra.SRAUtils;
import org.genvisis.sra.SRAUtils.SRABamWorker;
import org.genvisis.sra.SRAUtils.SRAConversionResult;

/**
 * more specific version of {@link Pipeline} that starts with a single SRA file
 *
 */
public class SRAPipeline implements Callable<List<PipelinePart>> {
	private static final String SRA_INPUT = "sraInput";
	private static final String SRA_RUN_TABLE = "sraRunTable";
	private static final String NUM_THREADS = "threads";
	private static final String NUM_THREADS_PIPELINE = "threadsPipe";

	private static final String CAPTURE_BED = "bed";
	private static final String BIN_BED = "bin";
	private static final String VCF = "vcf";
	private static final String NUM_BATCHES = "batch";
	private static final String COMPILE = "compile";
	private static final String COMPUTEL_LOCATION = "computelLocation";
	private static final String CLEANUP = "clean";
	private static final String FULL_PIPELINE = "full";

	private static final String PLATFORM_TYPE = "platform";

	private static final String GENVISIS_PART = "genvisis";
	private static final String MTDNACN_PART = "mtDNACN";
	private static final String TELSEQ_PART = "telseq";
	private static final String COMPUTEL_PART = "computel";
	private static final String ALL_PART = "all";
	private static final String SKIP_PRELIM = "skipPrelim";

	private final SRASample sraSample;
	private final String inputSRA;
	private final String rootOutDir;
	private final String captureBed;
	private final String refGenome;
	private final String binBed;
	private final String vcfFile;
	private final String computelLocation;
	private final List<PIPELINE_PARTS> partsToRun;
	private boolean cleanupTmpFiles;
	private final int numThreads;
	private final Logger log;

	/**
	 * @param sraSample an {@link SRASample} to analyze
	 * @param inputSRA the input sra file in appropriate sra-toolkit directory
	 * @param rootOutDir the output directory for the analysis
	 * @param genomeBuild reference genome build
	 * @param captureBed the capture bed, only utilized with {@link ASSAY_TYPE#WXS}
	 * @param atype {@link ASSAY_TYPE} of the sample
	 * @param aName {@link ASSEMBLY_NAME} for the sample
	 * @param numThreads number of threads for the pipeline branches
	 * @param cleanup clean up by deleting the .sra and .bam file after completion
	 * @param log
	 */
	public SRAPipeline(SRASample sraSample, String inputSRA, String rootOutDir,
										 String refGenome, String captureBed, String binBed, String vcf,
										 String computelLocation, int numThreads, boolean cleanup,
										 List<PIPELINE_PARTS> parts, Logger log) {
		super();
		this.sraSample = sraSample;
		this.inputSRA = inputSRA;
		this.rootOutDir = rootOutDir;
		this.captureBed = captureBed;
		this.refGenome = refGenome;
		this.binBed = binBed;
		vcfFile = vcf;
		this.computelLocation = computelLocation;
		this.numThreads = numThreads;
		this.cleanupTmpFiles = cleanup;
		partsToRun = parts;
		this.log = log;
	}

	@Override
	public List<PipelinePart> call() throws Exception {
		String bamDir = getBamDirectory(rootOutDir);
		if (!Files.exists(getCompleteFile(rootOutDir, inputSRA))) {
			new File(bamDir).mkdirs();
			String bam = bamDir + ext.rootOf(inputSRA) + ".bam";
			if (!Files.exists(bam)) {
				WorkerHive<SRAConversionResult> hive = new WorkerHive<SRAUtils.SRAConversionResult>(1, 10,
																																														log);
				hive.addCallable(new SRABamWorker(inputSRA, bam, log));
				// String vdbcache = inputSRA + ".vdbcache";
				// if (!Files.exists(vdbcache) && sraSample.getPlatform() != PLATFORM.ABI_SOLID) {
				// log.reportError(vdbcache
				// + " was was not seen, this sample might fail conversion to bam, skipping");
				// // String bamFailDirectory = rootOutDir + "bamFail/";
				// // new File(bamFailDirectory).mkdirs();
				// // Files.write(inputSRA, bamFailDirectory + ext.rootOf(inputSRA) + ".fail");
				// // return new ArrayList<Pipeline.PipelinePart>();
				// }
				hive.execute(true);
			}
			if (cleanupTmpFiles) {
				cleanupTmpFiles = new File(inputSRA).delete();
			}

			List<PipelinePart> parts = Pipeline.pipeline(bam, rootOutDir, captureBed, binBed,
																									 vcfFile, sraSample, refGenome, partsToRun,
																									 computelLocation,
																									 numThreads, log);

			if (cleanupTmpFiles) {
				cleanupTmpFiles = new File(bam).delete();
			}
			makeComplete(rootOutDir, inputSRA);
			return parts;

		} else {

			return new ArrayList<Pipeline.PipelinePart>();
		}
	}

	private static void makeComplete(String rootOutDir, String sra) {
		String file = getCompleteFile(rootOutDir, sra);
		new File(ext.parseDirectoryOfFile(file)).mkdirs();
		Files.write("complete", file);
	}

	private static String getCompleteFile(String rootOutDir, String sraFile) {
		return rootOutDir + "completes/" + ext.rootOf(sraFile) + ".complete";
	}

	private static String getBamDirectory(String rootOutDir) {
		return rootOutDir + "bams/";
	}

	private static List<SRASample> loadSraSamples(String sraInput, String sraRunTableFile,
																								PLATFORM platform, Logger log) {
		ArrayList<SRASample> samples = new ArrayList<SRASample>();
		SRARunTable srRunTable = SRARunTable.load(sraRunTableFile, platform, log);

		String[] sraFiles;
		if (Files.isDirectory(sraInput)) {
			log.reportTimeInfo("Gathering sra files from " + sraInput);
			sraFiles = Files.listFullPaths(sraInput, ".sra");
		} else {
			log.reportTimeInfo("Reading sra files from " + sraInput);
			sraFiles = HashVec.loadFileToStringArray(sraInput, false, new int[] {0}, true);
		}
		log.reportTimeInfo("Found " + sraFiles.length + " sra files in " + sraInput);
		for (String sraFile : sraFiles) {
			SRASample sample = srRunTable.get(ext.rootOf(sraFile));
			sample.setSraFile(sraFile);
			samples.add(sample);
		}

		return samples;
	}

	private static void runCompile(List<SRASample> samples, String rootOutDir,
																 String captureBed, String binBed,
																 String vcf, String referenceGenomeFasta, int numThreads,
																 Logger log) {


		ASSAY_TYPE atType = samples.get(0).getaType();
		ASSEMBLY_NAME aName = samples.get(0).getaName();

		for (SRASample sample : samples) {
			if (sample.getaType() != atType) {
				throw new IllegalArgumentException("Mismatched assay types");
			}
			if (sample.getaName() != aName) {
				throw new IllegalArgumentException("Mismatched assembly name");
			}

		}
		Project proj = Pipeline.getProjectFor(atType, rootOutDir);
		String serDir = proj.PROJECT_DIRECTORY.getValue() + "tmpBamSer/";
		String[] serFiles = Files.listFullPaths(serDir, ".ser");
		log.reportTimeInfo("Allocating fake bam files associated with " + serFiles.length
											 + " .ser files in " + serDir);

		String[] fakeBams = new String[serFiles.length];
		for (int i = 0; i < fakeBams.length; i++) {
			fakeBams[i] = ext.rootOf(serFiles[i], false) + ".bam";
		}

		BamImport.importTheWholeBamProject(proj, binBed, captureBed, vcf, BamImport.CAPTURE_BUFFER, 4,
																			 true, atType, aName, NORMALIZATON_METHOD.GENOME, fakeBams,
																			 referenceGenomeFasta, false, true,
																			 numThreads);
	}

	private static void compilePrep(String sraInput, String sraRunTableFile, String rootOutDir,
																	String refGenome, String captureBed, String binBed,
																	String vcf, PLATFORM platform, int numThreads) {
		Logger log = new Logger(rootOutDir + "compile.log");

		List<SRASample> samples = loadSraSamples(sraInput, sraRunTableFile, platform, log);
		ArrayList<SRASample> wgsSamples = new ArrayList<SRASample>();
		ArrayList<SRASample> wxsSamples = new ArrayList<SRASample>();
		for (SRASample sample : samples) {
			switch (sample.getaType()) {
				case WGS:
					wgsSamples.add(sample);
					break;
				case WXS:
					wxsSamples.add(sample);
					break;
				default:
					throw new IllegalArgumentException("Invalid assay type " + sample.getaType());
			}
		}

		log.reportTimeInfo("Found " + wgsSamples.size() + " " + ASSAY_TYPE.WGS + " samples and "
											 + wxsSamples.size() + " " + ASSAY_TYPE.WXS + " samples");
		runCompile(wxsSamples, rootOutDir, captureBed, binBed, vcf, refGenome, numThreads, log);
		runCompile(wgsSamples, rootOutDir, captureBed, binBed, vcf, refGenome, numThreads, log);

	}

	private static void runAll(String sraInput, String sraRunTableFile, String rootOutDir,
														 String refGenome, String captureBed, String binBed, String vcf,
														 String computelLocation, int numThreads, int numThreadsPipeline,
														 int numBatches, CLI c, boolean cleanup) {
		Logger log = new Logger();

		WorkerHive<List<PipelinePart>> hive = new WorkerHive<List<PipelinePart>>(numThreads, 10, log);
		boolean prelimGenvisisWGS = false;
		boolean prelimGenvisisWXS = false;
		List<SRASample> samples = loadSraSamples(sraInput, sraRunTableFile,
																						 PLATFORM.valueOf(c.get(PLATFORM_TYPE)), log);
		ArrayList<String> sampleSummary = new ArrayList<String>();
		ArrayList<String> sraFiles = new ArrayList<String>();

		ArrayList<PIPELINE_PARTS> partsToRun = new ArrayList<Pipeline.PIPELINE_PARTS>();
		if (c.has(ALL_PART)) {
			for (PIPELINE_PARTS part : PIPELINE_PARTS.values()) {
				partsToRun.add(part);
			}
		} else {
			if (c.has(COMPUTEL_PART)) {
				partsToRun.add(PIPELINE_PARTS.COMPUTEL);
			}
			if (c.has(TELSEQ_PART)) {
				partsToRun.add(PIPELINE_PARTS.TELSEQ);
			}
			if (c.has(MTDNACN_PART)) {
				partsToRun.add(PIPELINE_PARTS.MTDNACN);
			}
			if (c.has(GENVISIS_PART)) {
				partsToRun.add(PIPELINE_PARTS.GENVISIS);
			}

		}

		for (SRASample sample : samples) {
			sraFiles.add(sample.getSraFile());
			sampleSummary.add(sample.getSraFile() + "\t" + sample.toString());
			SRAPipeline pipeline = new SRAPipeline(sample, sample.getSraFile(), rootOutDir, refGenome,
																						 captureBed, binBed, vcf, computelLocation,
																						 numThreadsPipeline, cleanup, partsToRun, log);
			switch (sample.getaType()) {// create the required markerSets for
				// import...prior to threading
				case WGS:
					if (!prelimGenvisisWGS && !c.has(SKIP_PRELIM)) {
						generatePrelim(rootOutDir, refGenome, null, null, vcf, log, sample.getaType());
						prelimGenvisisWGS = true;
					}

					break;
				case WXS:
					if (!prelimGenvisisWXS && !c.has(SKIP_PRELIM)) {
						generatePrelim(rootOutDir, refGenome, captureBed, binBed, vcf, log,
													 sample.getaType());
						prelimGenvisisWXS = true;
					}
					break;
				default:
					throw new IllegalArgumentException("Invalid assay type " + sample.getaType());

			}
			hive.addCallable(pipeline);
		}
		Files.writeIterable(sampleSummary, rootOutDir + "sampleAnalysis.summary.txt");
		if (numBatches > 0) {
			batch(ArrayUtils.toStringArray(sraFiles), rootOutDir, c, log);
		} else {
			hive.execute(true);
		}
	}

	private static void generatePrelim(String rootOutDir, String refGenome, String captureBed,
																		 String binBed, String vcf, Logger log, ASSAY_TYPE aType) {
		Project proj = Pipeline.getProjectFor(aType, rootOutDir);
		if (!Files.exists(proj.MARKERSET_FILENAME.getValue())) {

			BamImport.generateAnalysisSet(proj, binBed, captureBed, vcf, BamImport.CAPTURE_BUFFER, aType,
																		true,
																		log, new ReferenceGenome(refGenome,
																														 log));
		}
	}

	private static void fullPipeline(CLI c) {
		Logger log = new Logger();
		PLATFORM platform = PLATFORM.valueOf(c.get(PLATFORM_TYPE));
		SRARunTable srRunTable = SRARunTable.load(c.get(SRA_RUN_TABLE), platform, log);
		String processDir = c.get(CLI.ARG_OUTDIR) + "process/";
		new File(processDir).mkdirs();

		String processFile = processDir + "process.sh";

		ArrayList<String> sraFilesToAnalyze = new ArrayList<String>();
		String[] sraFiles = ArrayUtils.tagOn(srRunTable.getAllRunSFiles(), c.get(SRA_INPUT), ".sra");
		for (int i = 0; i < sraFiles.length; i++) {
			if (!Files.exists(getCompleteFile(c.get(CLI.ARG_OUTDIR), sraFiles[i]))) {
				sraFilesToAnalyze.add(sraFiles[i]);
			}
		}
		log.reportTimeInfo("Detected " + sraFilesToAnalyze.size() + " samples to analyze for platform "
											 + platform + "," + (sraFiles.length - sraFilesToAnalyze.size())
											 + " samples are already complete");

		String[][] batches = batch(ArrayUtils.toStringArray(sraFilesToAnalyze), c.get(CLI.ARG_OUTDIR),
															 c, log);

		ArrayList<String> process = new ArrayList<String>();
		int num = 0;
		int processBatch = 0;
		// FIXME

		String plug = c.get(CLI.ARG_OUTDIR) + "sraPipline.plug";
		Files.write("Delete this file to stop disowned jobs", plug);

		log.reportTimeWarning("set up for lanej aspera");
		for (int i = 0; i < batches.length; i++) {
			for (int j = 0; j < batches[i].length; j++) {
				process.add("cd " + ext.parseDirectoryOfFile(batches[i][j]));
				process.add("echo \"start " + ext.rootOf(batches[i][j]) + "\" `date` >>" + processDir
										+ "sraDL.times");

				process.add("prefetch.2.8.2 -a \"/home/pankrat2/lane0212/.aspera/connect/bin/ascp|/home/pankrat2/lane0212/.aspera/connect/etc/asperaweb_id_dsa.openssh\" --max-size 100000000000 -Tr -Q -l 600M -L "
										+ ext.rootOf(batches[i][j]));
				process.add("echo \"end " + ext.rootOf(batches[i][j]) + "\" `date` >>"
										+ ext.parseDirectoryOfFile(batches[i][j]) + ".times");
				process.add("FILE=\"" + plug + "\"");
				process.add("if [ -f $FILE ]; then");
				process.add("echo \"File $FILE does exist.\"");
				process.add("else");
				process.add("echo \"File $FILE does not exist.\"");
				process.add("exit 1");
				process.add("fi");
				num++;
			}
			process.add("cd " + processDir);// so the qsubs get placed there
			process.add("qsub -q small " + getBatch(getBatchDirectory(c.get(CLI.ARG_OUTDIR)), i)
									+ ".qsub");
			if (num >= 1500) {
				Files.writeIterable(process, ext.addToRoot(processFile, "_" + processBatch));
				num = 0;
				processBatch++;
				process.clear();
			}
		}
		if (num > 0) {
			Files.writeIterable(process, ext.addToRoot(processFile, "_" + processBatch));
		}
		if (c.has(GENVISIS_PART) || c.has(ALL_PART)) {
			generatePrelim(c.get(CLI.ARG_OUTDIR), c.get(CLI.ARG_REFERENCE_GENOME), c.get(CAPTURE_BED),
										 c.get(BIN_BED),
										 c.get(VCF), log, ASSAY_TYPE.WXS);
			generatePrelim(c.get(CLI.ARG_OUTDIR), c.get(CLI.ARG_REFERENCE_GENOME), null, null, c.get(VCF),
										 log,
										 ASSAY_TYPE.WGS);
		}

	}

	private static String[][] batch(String[] sraFiles, String rootOutDir, CLI c, Logger log) {
		String[][] splits = ArrayUtils.splitUpStringArray(sraFiles, c.getI(NUM_BATCHES), log);
		String runningJar = SRAPipeline.class.getProtectionDomain().getCodeSource().getLocation()
																				 .getFile();
		PLATFORM platform = PLATFORM.valueOf(c.get(PLATFORM_TYPE));

		String jarRun = ext.parseDirectoryOfFile(runningJar) + "sraJars/" + ext.rootOf(runningJar)
										+ ext.getTimestampForFilename() + ".jar";
		new File(ext.parseDirectoryOfFile(jarRun)).mkdirs();
		Files.copyFileUsingFileChannels(runningJar, jarRun, log);
		ArrayList<String> baseCommand = new ArrayList<String>();
		baseCommand.add("module load gcc/4.8.1\n");
		baseCommand.add("java -Xmx60g -jar " + jarRun + " seq.analysis.genage.SRAPipeline");
		baseCommand.add(CLI.ARG_OUTDIR + "=" + c.get(CLI.ARG_OUTDIR));
		baseCommand.add(SRA_RUN_TABLE + "=" + c.get(SRA_RUN_TABLE));
		baseCommand.add(NUM_THREADS + "=" + c.get(NUM_THREADS));
		baseCommand.add(NUM_THREADS_PIPELINE + "=" + c.get(NUM_THREADS_PIPELINE));
		baseCommand.add(CLI.ARG_REFERENCE_GENOME + "=" + c.get(CLI.ARG_REFERENCE_GENOME));
		baseCommand.add(CAPTURE_BED + "=" + c.get(CAPTURE_BED));
		baseCommand.add(BIN_BED + "=" + c.get(BIN_BED));
		baseCommand.add(VCF + "=" + c.get(VCF));
		baseCommand.add(COMPUTEL_LOCATION + "=" + c.get(COMPUTEL_LOCATION));
		baseCommand.add(PLATFORM_TYPE + "=" + c.get(PLATFORM_TYPE));

		if (c.has(CLEANUP)) {
			baseCommand.add("-" + CLEANUP);
		}
		if (c.has(ALL_PART)) {
			baseCommand.add("-" + ALL_PART);
			if (platform != PLATFORM.ILLUMINA) {
				throw new IllegalArgumentException("Pipeline " + ALL_PART
																					 + " not completely implemented for " + platform);
			}
		}
		if (c.has(COMPUTEL_PART)) {
			baseCommand.add("-" + COMPUTEL_PART);
			if (platform != PLATFORM.ILLUMINA) {
				throw new IllegalArgumentException("Pipeline " + COMPUTEL_PART
																					 + " not completely implemented for " + platform);
			}
		}
		if (c.has(TELSEQ_PART)) {
			baseCommand.add("-" + TELSEQ_PART);
			if (platform != PLATFORM.ILLUMINA) {
				throw new IllegalArgumentException("Pipeline " + TELSEQ_PART
																					 + " not completely implemented for " + platform);
			}
		}
		if (c.has(MTDNACN_PART)) {
			baseCommand.add("-" + MTDNACN_PART);
		}
		if (c.has(GENVISIS_PART)) {
			baseCommand.add("-" + GENVISIS_PART);
			if (platform != PLATFORM.ILLUMINA) {
				throw new IllegalArgumentException("Pipeline " + GENVISIS_PART
																					 + " not completely implemented for " + platform);
			}
		}
		if (c.has(SKIP_PRELIM)) {
			baseCommand.add("-" + SKIP_PRELIM);
		}
		String batchDir = getBatchDirectory(rootOutDir);
		new File(batchDir).mkdirs();
		for (int i = 0; i < splits.length; i++) {
			String batch = getBatch(batchDir, i) + ".txt";
			String qsub = getBatch(batchDir, i) + ".qsub";
			Files.writeArray(splits[i], batch);
			ArrayList<String> currentCommand = new ArrayList<String>();
			currentCommand.addAll(baseCommand);
			currentCommand.add(SRA_INPUT + "=" + batch);
			Qsub.qsub(qsub, ArrayUtils.toStr(ArrayUtils.toStringArray(currentCommand), " "), 55000, 55,
								c.getI(NUM_THREADS) * c.getI(NUM_THREADS_PIPELINE));
		}
		return splits;
	}

	private static String getBatch(String batchDir, int index) {
		return batchDir + "batch_" + index;
	}

	private static String getBatchDirectory(String rootOutDir) {
		return rootOutDir + "batches/";
	}

	public static void main(String[] args) {

		CLI c = new CLI(SRAPipeline.class);

		String sraDirDefault = "sra/";
		c.addArgWithDefault(SRA_INPUT, "directory or filename with .sra files", sraDirDefault);

		String outDir = "out/";

		c.addArgWithDefault(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, outDir);

		String sraRunTableDefault = "sraRuntable.txt";

		c.addArgWithDefault(SRA_RUN_TABLE, "a sra run table providing sample information",
												sraRunTableDefault);

		int numThreads = 24;

		c.addArgWithDefault(NUM_THREADS, "number of threads across samples",
												Integer.toString(numThreads));

		int numThreadsPipe = 1;
		c.addArgWithDefault(NUM_THREADS_PIPELINE, "number of threads within samples",
												Integer.toString(numThreadsPipe));


		c.addArgWithDefault(CLI.ARG_REFERENCE_GENOME, CLI.DESC_REFERENCE_GENOME, "hg19.fa");

		String captureBedFile = "VCRome_2_1_hg19_capture_targets.bed";

		c.addArgWithDefault(CAPTURE_BED, "bed file of targeted capture", captureBedFile);

		String binBed = "targetsOfInterest.bed";

		c.addArgWithDefault(BIN_BED, "bed file of targets of interests", binBed);

		String vcf = "vcf.vcf";
		c.addArgWithDefault(VCF, "vcf file of variants", vcf);

		String computelLocation = null;
		c.addArgWithDefault(COMPUTEL_LOCATION, "directory of computel", computelLocation);

		int batch = -1;
		c.addArgWithDefault(NUM_BATCHES, "number of batches", Integer.toString(batch));

		c.addFlag(COMPILE, "Compile the genvisis portion of the pipeline");
		c.addFlag(CLEANUP, "Cleanup by deleting .sra and .bam after completion");
		c.addFlag(FULL_PIPELINE,
							"prepare for a run of the full pipeline, which will include downloading .sra files");

		c.addFlag(GENVISIS_PART, "run the genvisis portion of the pipeline");
		c.addFlag(MTDNACN_PART, "run the mtDNA CN portion of the pipeline");
		c.addFlag(TELSEQ_PART, "run the telseq portion of the pipeline");
		c.addFlag(COMPUTEL_PART, "run the computel portion of the pipeline");
		c.addFlag(ALL_PART, "run the entire pipeline");
		c.addFlag(SKIP_PRELIM, "skip preliminary genvisis setup");

		c.addArgWithDefault(PLATFORM_TYPE, "platform to analyze", PLATFORM.ILLUMINA.toString());

		c.parseWithExit(args);


		if (c.has(COMPILE)) {
			compilePrep(c.get(SRA_INPUT), c.get(SRA_RUN_TABLE), c.get(CLI.ARG_OUTDIR),
									c.get(CLI.ARG_REFERENCE_GENOME),
									c.get(CAPTURE_BED), c.get(BIN_BED), c.get(VCF),
									PLATFORM.valueOf(c.get(PLATFORM_TYPE)), c.getI(NUM_THREADS));
		} else if (c.has(FULL_PIPELINE)) {
			fullPipeline(c);
		} else {
			runAll(c.get(SRA_INPUT), c.get(SRA_RUN_TABLE), c.get(CLI.ARG_OUTDIR),
						 c.get(CLI.ARG_REFERENCE_GENOME),
						 c.get(CAPTURE_BED), c.get(BIN_BED), c.get(VCF), c.get(COMPUTEL_LOCATION),
						 c.getI(NUM_THREADS), c.getI(NUM_THREADS_PIPELINE), c.getI(NUM_BATCHES), c,
						 c.has(CLEANUP));
		}
	}
}
