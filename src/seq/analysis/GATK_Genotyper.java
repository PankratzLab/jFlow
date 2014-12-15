package seq.analysis;

import java.io.File;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;

public class GATK_Genotyper {
	public static final String SPACE = " ";
	private GATK gatk;
	private GATK.SingleSampleHaplotypeCaller[] siSampleHaplotypeCallers;
	private boolean fail, verbose;
	private int numBetweenSampleThreads;
	private int numWithinSampleThreads;
	private Logger log;

	public GATK_Genotyper(GATK gatk, int numBetweenSampleThreads, int numWithinSampleThreads, boolean verbose, Logger log) {
		super();
		this.gatk = gatk;
		this.numBetweenSampleThreads = numBetweenSampleThreads;
		this.numWithinSampleThreads = numWithinSampleThreads;
		this.log = log;
	}

	public boolean isFail() {
		return fail;
	}

	public void setFail(boolean fail) {
		this.fail = fail;
	}

	public boolean isVerbose() {
		return verbose;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	public void runJointGenotyping(JointGATKGenotyper jGatkGenotyper) {
		boolean progress = true;
		System.out.println("HIHD");
		if (!jGatkGenotyper.isFail()) {
			progress = gatk.jointGenotypeGVCFs(jGatkGenotyper.getInputGVCFs(), jGatkGenotyper.getOutputVCF(), numWithinSampleThreads, jGatkGenotyper.getLog());
			jGatkGenotyper.setFail(!progress);
		}
	}

	public void runRecalibration(JointGATKGenotyper jGatkGenotyper) {
		if (!jGatkGenotyper.isFail()) {
			jGatkGenotyper = gatk.recalibrateAVCF(jGatkGenotyper, numWithinSampleThreads, log);
		}
	}

	public void batch(JointGATKGenotyper jointGATKGenotyper, String rootOutputDir, int memoryInMB, int wallTimeInHours, String baseName) {
		// TODO, change classpath
		String command = "module load R\nmodule load java\njava -cp parkGATK.jar -Xmx" + memoryInMB + "m seq.analysis.GATK_Genotyper ";
		command += GATK_LanePrep.ROOT_INPUT_COMMAND + jointGATKGenotyper.getRootInputDir() + SPACE;
		command += GATK_LanePrep.ROOT_OUTPUT_COMMAND + rootOutputDir + SPACE;
		command += GATK_LanePrep.REFERENCE_GENOME_COMMAND + gatk.getReferenceGenomeFasta() + SPACE;
		command += GATK.GATK_LOCATION_COMMAND + gatk.getGATKLocation() + SPACE;
		command += NUM_THREADS + numWithinSampleThreads + SPACE;
		command += OUTPUT + jointGATKGenotyper.getOutput() + SPACE;
		command += GATK_LanePrep.LOG_FILE_COMMAND + jointGATKGenotyper.getLog().getFilename() + SPACE;
		if (jointGATKGenotyper.getFileOfGVCFs() != null) {
			command += FILE_OF_GVCFS + jointGATKGenotyper.getFileOfGVCFs() + SPACE;
		}
		command += HAPMAP_TRAIN + gatk.getHapMapTraining() + SPACE;
		command += OMNI_TRAIN + gatk.getOmniTraining() + SPACE;
		command += G1000 + gatk.getThousandGTraining() + SPACE;
		command += DBSNP + gatk.getDbSnpTraining() + SPACE;
		command += MILLS + gatk.getMillsIndelTraining() + SPACE;
		Files.qsub("GATK_Genotype_" + baseName, command, memoryInMB, wallTimeInHours, numWithinSampleThreads);
	}

	public void runSingleSampleAllSites(String[] inputBams) {

		if (!isFail() && Files.checkAllFiles("", inputBams, verbose, log)) {
			if (inputBams != null) {
				this.siSampleHaplotypeCallers = new GATK.SingleSampleHaplotypeCaller[inputBams.length];
				int[] actualWithinSampleThreads = optimizeThreads(siSampleHaplotypeCallers.length, numBetweenSampleThreads, numWithinSampleThreads, log);

				ExecutorService executor = Executors.newFixedThreadPool(numBetweenSampleThreads);
				Hashtable<String, Future<GATK.SingleSampleHaplotypeCaller>> tmpResults = new Hashtable<String, Future<GATK.SingleSampleHaplotypeCaller>>();
				for (int i = 0; i < inputBams.length; i++) {
					Logger altLog = new Logger(ext.rootOf(inputBams[i], false) + ".HC_ERC.log");
					tmpResults.put(i + "", executor.submit(new WorkerSingleSampleAllSites(gatk, inputBams[i], ext.rootOf(inputBams[i]), actualWithinSampleThreads[i], altLog)));
				}
				for (int i = 0; i < siSampleHaplotypeCallers.length; i++) {
					try {
						siSampleHaplotypeCallers[i] = tmpResults.get(i + "").get();
						if (siSampleHaplotypeCallers[i].isFail() && !isFail()) {
							log.reportError("Error - failed single sample haplotype calling for " + siSampleHaplotypeCallers[i].getInputBam());
							setFail(true);
						}
					} catch (InterruptedException e) {
						log.reportError("Error - when running GATK single sample haplotype calling on internal index " + i);
						log.reportException(e);
						setFail(true);
					} catch (ExecutionException e) {
						log.reportError("Error - when running GATK single sample haplotype calling on internal index " + i);
						log.reportException(e);
						setFail(true);
					}
				}
				executor.shutdown();
				try {
					executor.awaitTermination(10, TimeUnit.DAYS);
				} catch (InterruptedException e) {
					log.reportException(e);
				}
			} else {
				// TODO better check
			}
		}
	}

	private static int[] optimizeThreads(int numInputs, int numBetweenSampleThreads, int numWithinSampleThreads, Logger log) {
		int[] optimizedWithin = new int[numInputs];
		Arrays.fill(optimizedWithin, numWithinSampleThreads);
		if (numInputs < numBetweenSampleThreads) {
			int numExtra = numBetweenSampleThreads - numInputs;
			int index = 0;
			while (numExtra > 0) {
				numExtra--;
				optimizedWithin[index]++;
				index++;
				if (index == numInputs) {
					index = 0;
				}
			}
			log.report(ext.getTime() + " Info - since we have extra between sample threads, we will allocate some to within sample(s) threads");
			log.report(ext.getTime() + " Info - allocated " + (numBetweenSampleThreads - numInputs) + " extra thread(s) within sample(s) as follows: " + Array.toStr(optimizedWithin));
		}
		return optimizedWithin;
	}

	private static class WorkerSingleSampleAllSites implements Callable<GATK.SingleSampleHaplotypeCaller> {
		private GATK GATK;
		private String inputBam, baseId;
		private int numWithinSampleThreads;
		private Logger altLog;

		public WorkerSingleSampleAllSites(seq.analysis.GATK gATK, String inputBam, String baseId, int numWithinSampleThreads, Logger altLog) {
			super();
			this.GATK = gATK;
			this.inputBam = inputBam;
			this.baseId = baseId;
			this.numWithinSampleThreads = numWithinSampleThreads;
			this.altLog = altLog;
		}

		@Override
		public GATK.SingleSampleHaplotypeCaller call() {
			return GATK.haplotypeCallABam(baseId, inputBam, numWithinSampleThreads, altLog);
		}
	}

	public static class JointGATKGenotyper {
		public static final String RECAL_EXT = ".recal";
		public static final String TRANCHES_EXT = ".tranches";
		public static final String RScript_EXT = ".R";

		private String rootInputDir, rootOutputDir, output;
		private String rawVCF;
		private String fileOfGVCFs, recalSNP_VCF_File, recalSNP_Indel_VCF_File;
		private String recalSNPFile, tranchesSNPFile, rscriptSNPFile;
		private String recalINDELFile, tranchesINDELFile, rscriptINDELFile;

		private String[] inputGVCFs;
		private boolean fail;
		private Logger log;

		public JointGATKGenotyper(String rootInputDir, String rootOutputDir, String output, Logger log) {
			super();
			this.rootInputDir = rootInputDir;
			this.rootOutputDir = rootOutputDir;
			this.output = output;
			this.fileOfGVCFs = null;
			this.log = log;
			this.fail = false;
		}

		public void init(String fileOfGVCFs) {
			String currentRoot = rootOutputDir + ext.rootOf(output);
			this.fileOfGVCFs = fileOfGVCFs;
			this.rawVCF = currentRoot + GATK.VCF;

			this.recalSNPFile = currentRoot + GATK.SNP + RECAL_EXT;
			this.tranchesSNPFile = currentRoot + GATK.SNP + TRANCHES_EXT;
			this.rscriptSNPFile = currentRoot + GATK.SNP + RScript_EXT;

			this.recalINDELFile = currentRoot + GATK.INDEL + RECAL_EXT;
			this.tranchesINDELFile = currentRoot + GATK.INDEL + TRANCHES_EXT;
			this.rscriptINDELFile = currentRoot + GATK.INDEL + RScript_EXT;

			this.recalSNP_VCF_File = ext.addToRoot(rawVCF, RECAL_EXT + GATK.SNP);
			this.recalSNP_Indel_VCF_File = ext.addToRoot(rawVCF, RECAL_EXT + GATK.INDEL);

			if (fileOfGVCFs != null) {
				log.report(ext.getTime() + " Info - using GVCF files listed in the first column of" + fileOfGVCFs);
				this.inputGVCFs = HashVec.loadFileToStringArray(fileOfGVCFs, false, new int[] { 0 }, true);
			} else if (rootInputDir == null) {
				log.reportError("Error - a file listing GVCF files was not provided and the root input directory was not provided, halting...");
				fail = true;
			} else {
				log.report(ext.getTime() + " Info - finding files with extension " + GATK.GVCF + " in " + rootInputDir);
				this.inputGVCFs = Files.toFullPaths(Files.list(rootInputDir, GATK.GVCF, false), rootInputDir);
			}
			if (inputGVCFs == null || inputGVCFs.length < 1) {
				log.reportError("Error - could not find any GVCF files to joint genotype");
				fail = true;
			} else {
				log.report(ext.getTime() + " Info - using " + inputGVCFs.length + " file(s) for joint Genotyping");
			}

		}

		public String getOutput() {
			return output;
		}

		public String getOutputVCF() {
			return rawVCF;
		}

		public String getFileOfGVCFs() {
			return fileOfGVCFs;
		}

		public boolean isFail() {
			return fail;
		}

		public void setFail(boolean fail) {
			this.fail = fail;
		}

		public String getRootInputDir() {
			return rootInputDir;
		}

		public String[] getInputGVCFs() {
			return inputGVCFs;
		}

		public Logger getLog() {
			return log;
		}

		public String getRootOutputDir() {
			return rootOutputDir;
		}

		public static String getRecalExt() {
			return RECAL_EXT;
		}

		public static String getTranchesExt() {
			return TRANCHES_EXT;
		}

		public static String getRscriptExt() {
			return RScript_EXT;
		}

		public String getRawVCF() {
			return rawVCF;
		}

		public String getRecalSNP_VCF_File() {
			return recalSNP_VCF_File;
		}

		public String getRecalSNP_Indel_VCF_File() {
			return recalSNP_Indel_VCF_File;
		}

		public String getRecalSNPFile() {
			return recalSNPFile;
		}

		public String getTranchesSNPFile() {
			return tranchesSNPFile;
		}

		public String getRscriptSNPFile() {
			return rscriptSNPFile;
		}

		public String getRecalINDELFile() {
			return recalINDELFile;
		}

		public String getTranchesINDELFile() {
			return tranchesINDELFile;
		}

		public String getRscriptINDELFile() {
			return rscriptINDELFile;
		}

	}

	public static void jointGenotype(String rootInputDir, String rootOutputDir, String output, String gATKLocation, String referenceGenomeFasta, String fileOfGVCFs, String hapMapTraining, String omniTraining, String thousandGTraining, String dbSnpTraining, String millsIndelTraining, boolean verbose, boolean overwriteExisting, boolean batch, int numThreads, int memoryInMB, int wallTimeInHours, Logger log) {
		GATK gatk = new GATK(gATKLocation, referenceGenomeFasta, null, null, null, verbose, overwriteExisting, log);
		gatk.setDbSnpTraining(dbSnpTraining);
		gatk.setHapMapTraining(hapMapTraining);
		gatk.setOmniTraining(omniTraining);
		gatk.setThousandGTraining(thousandGTraining);
		gatk.setMillsIndelTraining(millsIndelTraining);
		GATK_Genotyper genotyper = new GATK_Genotyper(gatk, 0, numThreads, verbose, log);
		JointGATKGenotyper jGatkGenotyper = new JointGATKGenotyper(rootInputDir, rootOutputDir, output, log);
		jGatkGenotyper.init(fileOfGVCFs);
		new File(rootOutputDir).mkdirs();
		if (batch) {
			genotyper.batch(jGatkGenotyper, rootOutputDir, memoryInMB, wallTimeInHours, output);
		} else {
			genotyper.runJointGenotyping(jGatkGenotyper);
			genotyper.runRecalibration(jGatkGenotyper);
		}
	}

	public static final String NUM_THREADS = "numThreads=";
	public static final String FILE_OF_GVCFS = "gvcfs=";
	public static final String OUTPUT = "output=";
	public static final String HAPMAP_TRAIN = "hapmap=";
	public static final String OMNI_TRAIN = "omni=";
	public static final String G1000 = "1000G=";
	public static final String DBSNP = "dbSNP=";
	public static final String MILLS = "mills=";

	public static void main(String[] args) {

		int numArgs = args.length;
		String rootInputDir = null;
		String rootOutputDir = null;
		String referenceGenomeFasta = "";
		String gATKLocation = "";
		String output = "joint_genotypes";
		String fileOfGVCFs = null;
		boolean verbose = true;
		boolean batch = false;
		int memoryInMB = 23000;
		int wallTimeInHours = 24;
		int numThreads = 1;
		boolean overwriteExisting = false;
		String hapMapTraining = null;
		String omniTraining = null;
		String thousandGTraining = null;
		String dbSnpTraining = null;
		String millsIndelTraining = null;

		String logFile = "GATK_GENOTYPE.log";

		String usage = "\n" + "seq.GATK_Genotyper requires 2 argument\n";
		usage += "   (1) root input directory (i.e. " + GATK_LanePrep.ROOT_INPUT_COMMAND + rootInputDir + " (no default))\n" + "";
		usage += "   (2) root output directory (i.e. " + GATK_LanePrep.ROOT_OUTPUT_COMMAND + rootOutputDir + " (no default))\n" + "";
		usage += "   (3) tab-delimited file with no header of (i.e. " + FILE_OF_GVCFS + fileOfGVCFs + " (optional, no default))\n" + "";
		usage += "   (4) the full path to a  reference genome in fasta format (i.e." + GATK_LanePrep.REFERENCE_GENOME_COMMAND + referenceGenomeFasta + " (no default))\n" + "";
		usage += "   (5) the full path to the GATK executable (i.e. " + GATK.GATK_LOCATION_COMMAND + gATKLocation + " (defualts to systems path))\n" + "";

		usage += "   (6) run in quiet mode (i.e. " + GATK_LanePrep.QUIET_COMMAND + " (not tbe default))\n" + "";
		usage += "   (7) number of  threads for analysis(i.e." + NUM_THREADS + numThreads + " (default))\n" + "";

		usage += "   (8) filename for a log (i.e. " + GATK_LanePrep.LOG_FILE_COMMAND + logFile + " (default))\n" + "";
		usage += "   (9) over-write exsiting files (i.e. " + GATK_LanePrep.OVERWRITE_EXISTING_COMMAND + " (not the default))\n" + "";
		usage += "   (10) set up a batch analysis for the root input directory for a log (i.e. " + GATK_LanePrep.BATCH_COMMAND + " (not the default))\n" + "";
		usage += "   (11) root output for analysis (i.e. " + OUTPUT + output + " ( default))\n" + "";
		usage += "   (12) HapMap SNP Training Referenece (i.e. " + HAPMAP_TRAIN + " ( no default))\n" + "";
		usage += "   (13) Omni SNP Training Referenece (i.e. " + OMNI_TRAIN + " ( no default))\n" + "";
		usage += "   (14) 1000 SNP genomes Training Referenece (i.e. " + G1000 + " ( no default))\n" + "";
		usage += "   (15) dbSNP SNP Training Referenece (i.e. " + DBSNP + " ( no default))\n" + "";
		usage += "   (16) mills INDEL Training Referenece (i.e. " + MILLS + " ( no default))\n" + "";
		usage += "   (17) log file name (i.e. " + MILLS + " ( no default))\n" + "";

		for (int i = 0; i < args.length; i++) {

			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith(HAPMAP_TRAIN)) {
				hapMapTraining = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(OMNI_TRAIN)) {
				omniTraining = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(G1000)) {
				thousandGTraining = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(DBSNP)) {
				dbSnpTraining = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(MILLS)) {
				millsIndelTraining = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(GATK_LanePrep.ROOT_INPUT_COMMAND)) {
				rootInputDir = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(GATK_LanePrep.ROOT_OUTPUT_COMMAND)) {
				rootOutputDir = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(GATK_LanePrep.FILE_OF_SAMPLE_PAIRS_COMMAND)) {
				fileOfGVCFs = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(GATK_LanePrep.REFERENCE_GENOME_COMMAND)) {
				referenceGenomeFasta = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(GATK_LanePrep.LOG_FILE_COMMAND)) {
				logFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(NUM_THREADS)) {
				numThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("memoryInMB=")) {
				memoryInMB = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("wallTimeInHours=")) {
				wallTimeInHours = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith(GATK_LanePrep.QUIET_COMMAND)) {
				verbose = false;
				numArgs--;
			} else if (args[i].startsWith(GATK_LanePrep.BATCH_COMMAND)) {
				batch = true;
				numArgs--;
			} else if (args[i].startsWith(GATK_LanePrep.OVERWRITE_EXISTING_COMMAND)) {
				batch = true;
				numArgs--;
			} else if (args[i].startsWith(GATK.GATK_LOCATION_COMMAND)) {
				gATKLocation = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(OUTPUT)) {
				output = ext.parseStringArg(args[i], "");
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.exit(1);
		}
		Logger log = new Logger(rootOutputDir + logFile);
		jointGenotype(rootInputDir, rootOutputDir, output, gATKLocation, referenceGenomeFasta, fileOfGVCFs, hapMapTraining, omniTraining, thousandGTraining, dbSnpTraining, millsIndelTraining, verbose, overwriteExisting, batch, numThreads, memoryInMB, wallTimeInHours, log);
	}

}
