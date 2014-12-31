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

import seq.analysis.SNPEFF.SnpEffResult;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;

public class GATK_Genotyper {
	public static final String SPACE = " ";
	private GATK gatk;
	private SNPEFF snpeff;
	private GATK.SingleSampleHaplotypeCaller[] siSampleHaplotypeCallers;
	private boolean fail, verbose;
	private int numBetweenSampleThreads;
	private int numWithinSampleThreads;
	private Logger log;

	public GATK_Genotyper(GATK gatk, SNPEFF snpeff, int numBetweenSampleThreads, int numWithinSampleThreads, boolean verbose, Logger log) {
		super();
		this.gatk = gatk;
		this.snpeff = snpeff;
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

	public SNPEFF getSnpeff() {
		return snpeff;
	}

	public void runJointGenotyping(JointGATKGenotyper jGatkGenotyper) {
		boolean progress = true;
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

	public void annotateVCF(String inputVCF, String build) {
		SnpEffResult snpEffResult =snpeff.annotateAVCF(inputVCF, build);
		gatk.annotateAVcfWithSnpEFF(snpEffResult);
	}

	public void batch(JointGATKGenotyper jointGATKGenotyper, String rootOutputDir, boolean annotate, int memoryInMB, int wallTimeInHours, String baseName) {
		// TODO, change classpath
		String command = "module load R\nmodule load java\njava -cp parkGATK.jar -Xmx" + memoryInMB + "m seq.analysis.GATK_Genotyper ";
		command += GATK_LanePrep.ROOT_INPUT_COMMAND + jointGATKGenotyper.getRootInputDir() + SPACE;
		command += GATK_LanePrep.ROOT_OUTPUT_COMMAND + rootOutputDir + SPACE;
		command += GATK_LanePrep.REFERENCE_GENOME_COMMAND + gatk.getReferenceGenomeFasta() + SPACE;
		command += GATK.GATK_LOCATION_COMMAND + gatk.getGATKLocation() + SPACE;
		command += NUM_THREADS + numWithinSampleThreads + SPACE;
		command += OUTPUT_COMMAND + jointGATKGenotyper.getOutput() + SPACE;
		command += GATK_LanePrep.LOG_FILE_COMMAND + ext.removeDirectoryInfo(jointGATKGenotyper.getLog().getFilename()) + SPACE;
		if (jointGATKGenotyper.getFileOfGVCFs() != null) {
			command += FILE_OF_GVCFS + jointGATKGenotyper.getFileOfGVCFs() + SPACE;
		}
		command += HAPMAP_TRAIN_COMMAND + gatk.getHapMapTraining() + SPACE;
		command += OMNI_TRAIN_COMMAND + gatk.getOmniTraining() + SPACE;
		command += G1000_COMMAND + gatk.getThousandGTraining() + SPACE;
		command += DBSNP_COMMMAND + gatk.getDbSnpTraining() + SPACE;
		command += MILLS + gatk.getMillsIndelTraining() + SPACE;
		if (annotate) {
			command += SNPEFF.SNP_EFF_COMMAND + snpeff.getSnpEffLocation() + SPACE;
		} else {
			command += SNPEFF.SNP_EFF_NO_ANNO_COMMAND;
		}
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

			this.recalSNPFile = currentRoot + "." + GATK.SNP + RECAL_EXT;
			this.tranchesSNPFile = currentRoot + "." + GATK.SNP + TRANCHES_EXT;
			this.rscriptSNPFile = currentRoot + "." + GATK.SNP + RScript_EXT;

			this.recalINDELFile = currentRoot + "." + GATK.INDEL + RECAL_EXT;
			this.tranchesINDELFile = currentRoot + "." + GATK.INDEL + TRANCHES_EXT;
			this.rscriptINDELFile = currentRoot + "." + GATK.INDEL + RScript_EXT;

			this.recalSNP_VCF_File = ext.addToRoot(rawVCF, "." + GATK.SNP + RECAL_EXT);
			this.recalSNP_Indel_VCF_File = ext.addToRoot(recalSNP_VCF_File, "." + GATK.INDEL + RECAL_EXT);

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

	public static void jointGenotype(String rootInputDir, String rootOutputDir, String output, String gATKLocation, String referenceGenomeFasta, String fileOfGVCFs, String hapMapTraining, String omniTraining, String thousandGTraining, String dbSnpTraining, String millsIndelTraining, String snpEffLocation, String annoBuild, boolean verbose, boolean overwriteExisting, boolean batch, boolean annotate, int numThreads, int memoryInMB, int wallTimeInHours, Logger log) {
		GATK gatk = new GATK(gATKLocation, referenceGenomeFasta, null, null, null, verbose, overwriteExisting, log);
		gatk.setDbSnpTraining(dbSnpTraining);
		gatk.setHapMapTraining(hapMapTraining);
		gatk.setOmniTraining(omniTraining);
		gatk.setThousandGTraining(thousandGTraining);
		gatk.setMillsIndelTraining(millsIndelTraining);
		SNPEFF snpeff = new SNPEFF(snpEffLocation, verbose, overwriteExisting, log);
		GATK_Genotyper genotyper = new GATK_Genotyper(gatk, snpeff, 0, numThreads, verbose, log);
		JointGATKGenotyper jGatkGenotyper = new JointGATKGenotyper(rootInputDir, rootOutputDir, output, log);
		jGatkGenotyper.init(fileOfGVCFs);
		new File(rootOutputDir).mkdirs();
		if (batch) {
			genotyper.batch(jGatkGenotyper, rootOutputDir, annotate, memoryInMB, wallTimeInHours, output);
		} else {
			genotyper.runJointGenotyping(jGatkGenotyper);
			genotyper.runRecalibration(jGatkGenotyper);
			if (annotate) {
				// /
				genotyper.annotateVCF(jGatkGenotyper.getRecalSNP_Indel_VCF_File(), annoBuild);
			}
		}
	}

	public static final String NUM_THREADS = "numThreads=";
	public static final String FILE_OF_GVCFS = "gvcfs=";
	public static final String OUTPUT_COMMAND = "output=";
	public static final String HAPMAP_TRAIN_COMMAND = "hapmap=";
	public static final String OMNI_TRAIN_COMMAND = "omni=";
	public static final String G1000_COMMAND = "1000G=";
	public static final String DBSNP_COMMMAND = "dbSNP=";
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
		String snpEffLocation = "";
		String annoBuild = SNPEFF.BUILDS[0];
		boolean annotate = true;

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
		usage += "   (11) root output for analysis (i.e. " + OUTPUT_COMMAND + output + " ( default))\n" + "";
		usage += "   (12) HapMap SNP Training Referenece (i.e. " + HAPMAP_TRAIN_COMMAND + " ( no default))\n" + "";
		usage += "   (13) Omni SNP Training Referenece (i.e. " + OMNI_TRAIN_COMMAND + " ( no default))\n" + "";
		usage += "   (14) 1000 SNP genomes Training Referenece (i.e. " + G1000_COMMAND + " ( no default))\n" + "";
		usage += "   (15) dbSNP SNP Training Referenece (i.e. " + DBSNP_COMMMAND + " ( no default))\n" + "";
		usage += "   (16) mills INDEL Training Referenece (i.e. " + MILLS + " ( no default))\n" + "";
		usage += "   (17) full path to the SNP EFF directory (i.e. " + SNPEFF.SNP_EFF_COMMAND + " ( no default))\n" + "";
		usage += "   (18) the build version for SNP EFF annotation (options are " + Array.toStr(SNPEFF.BUILDS, ", ") + " (i.e. " + SNPEFF.SNP_EFF_BUILD_COMMAND + annoBuild + " ( default))\n" + "";

		usage += "   (19) do not annotate with SNPEFF (i.e. " + SNPEFF.SNP_EFF_NO_ANNO_COMMAND + " ( not the default))\n" + "";

		usage += "   (18) log file name (i.e. " + MILLS + " ( no default))\n" + "";

		for (int i = 0; i < args.length; i++) {

			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith(HAPMAP_TRAIN_COMMAND)) {
				hapMapTraining = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(OMNI_TRAIN_COMMAND)) {
				omniTraining = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(G1000_COMMAND)) {
				thousandGTraining = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(DBSNP_COMMMAND)) {
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
			} else if (args[i].startsWith(SNPEFF.SNP_EFF_COMMAND)) {
				snpEffLocation = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(SNPEFF.SNP_EFF_NO_ANNO_COMMAND)) {
				annotate = false;
				numArgs--;
			} else if (args[i].startsWith(OUTPUT_COMMAND)) {
				output = ext.parseStringArg(args[i], "");
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.exit(1);
		}
		new File(rootOutputDir).mkdirs();
		Logger log = new Logger(rootOutputDir + logFile);
		jointGenotype(rootInputDir, rootOutputDir, output, gATKLocation, referenceGenomeFasta, fileOfGVCFs, hapMapTraining, omniTraining, thousandGTraining, dbSnpTraining, millsIndelTraining, snpEffLocation, annoBuild, verbose, overwriteExisting, batch, annotate, numThreads, memoryInMB, wallTimeInHours, log);
	}

}
