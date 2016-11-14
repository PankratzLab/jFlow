package org.genvisis.seq.analysis;

import java.io.File;
import java.util.Arrays;
import java.util.concurrent.Callable;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.seq.analysis.ANNOVAR.AnnovarResults;
import org.genvisis.seq.analysis.SNPEFF.SnpEffResult;
import org.genvisis.seq.manage.VCFOps;

public class GATK_Genotyper {
	public static final String SPACE = " ";
	private final GATK gatk;
	private final SNPEFF snpeff;
	private final ANNOVAR annovar;
	private final SNPSIFT snpsift;
	private GATK.SingleSampleHaplotypeCaller[] siSampleHaplotypeCallers;
	private boolean fail, verbose;
	private final int numBetweenSampleThreads;
	private final int numWithinSampleThreads;
	private final Logger log;

	public GATK_Genotyper(GATK gatk, SNPEFF snpeff, SNPSIFT snpsift, ANNOVAR annovar,
												int numBetweenSampleThreads, int numWithinSampleThreads, boolean verbose,
												Logger log) {
		super();
		this.gatk = gatk;
		this.snpeff = snpeff;
		this.snpsift = snpsift;
		this.annovar = annovar;
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

	public ANNOVAR getAnnovar() {
		return annovar;
	}

	public boolean runJointGenotyping(JointGATKGenotyper jGatkGenotyper) {
		boolean progress = !jGatkGenotyper.isFail();
		if (progress) {
			progress = gatk.jointGenotypeGVCFs(	jGatkGenotyper.getInputGVCFs(),
																					jGatkGenotyper.getOutputVCF(),
																					jGatkGenotyper.getRestrictionContig(),
																					numWithinSampleThreads, jGatkGenotyper.getLog());
			jGatkGenotyper.setFail(!progress);
		}
		return progress;
	}

	public boolean runRecalibration(final JointGATKGenotyper jGatkGenotyper) {
		boolean progress = !jGatkGenotyper.isFail();
		if (progress) {
			gatk.recalibrateAVCF(jGatkGenotyper, numWithinSampleThreads, log);
		}
		return progress;
	}

	public String annotateVCF(String inputVCF, String build, MergeVCF mergeVCF, ANNOVCF annoVCF) {
		if (!fail) {
			String in = inputVCF;
			String out = "";
			if (mergeVCF != null) {
				log.reportTimeInfo("Applying merge " + mergeVCF.getArg() + " prior to annotation");
				in = VCFOps.getAppropriateRoot(inputVCF, false) + ".merge_" + mergeVCF.getTag() + ".vcf";
				log.reportTimeInfo("Output merge: " + in);
				gatk.mergeVCFs(	Array.concatAll(new String[] {inputVCF}, mergeVCF.getVcfsToMergeWith()), in,
												numWithinSampleThreads, false, log);
				out = in;
			}

			if (!annovar.isFail()) {
				AnnovarResults annovarResults = annovar.AnnovarAVCF(in, build, numWithinSampleThreads, log);
				in = annovarResults.getOutputVCF();
				out = annovarResults.getOutputVCF();
			}
			if (!snpeff.isFail()) {
				SnpEffResult snpEffResult = snpeff.annotateAVCF(in, build);
				gatk.annotateAVcfWithSnpEFF(snpEffResult, false);
				out = snpEffResult.getOutputGatkSnpEffVCF();
				if (!snpEffResult.isFail()) {
					if (annoVCF != null) {
						log.reportTimeInfo("Applying annotations from " + annoVCF.getArg());
						out = VCFOps.getAppropriateRoot(snpEffResult.getOutputGatkSnpEffVCF(), false)+ ".anno_"
									+ annoVCF.getTag() + ".vcf";
						log.reportTimeInfo("Output anno: " + out);
						gatk.annotateWithAnotherVCF(snpEffResult.getOutputGatkSnpEffVCF(), annoVCF.getVcf(),
																				out, annoVCF.getAnnos(), annoVCF.getTag(), null,
																				numWithinSampleThreads);
						in = out;
					}
					log.reportTimeInfo("Since Annovar uses invalid char sequence for 1000g2014oct_* and 1000g2015aug_* (starts with number), replacing with g10002014oct_* or g10002015aug_*");
					log.reportTimeInfo("Note this is a hard-coded sed, so...");
					in = sed1000g(out, log);
					out = VCFOps.gzipAndIndex(in, log);
					// SnpSiftResult ssr = snpsift.annotateDbnsfp(in, log);
					// out = ssr.getOutputVCF();
				}
			}
			return out;
		}
		return null;
	}

	private String sed1000g(String in, Logger log) {
		String out = ext.addToRoot(in, ".sed1000g");
		String command = "cat "+ in
											+ "|sed 's/1000g2014oct_/g10002014oct_/g'|sed 's/1000g2015aug_/g10002015aug_/g'>"
											+ out;
		String[] bat = CmdLine.prepareBatchForCommandLine(new String[] {command}, out + ".bat", true,
																											log);
		if (CmdLine.runCommandWithFileChecks(	bat, "", new String[] {in}, new String[] {out}, true,
																					false, false, log)) {
			return out;
		} else {
			return null;
		}
	}

	public boolean determineTsTV(String inputVCF) {
		return !snpsift.tsTv(inputVCF, log).isFail();
	}

	public void batch(JointGATKGenotyper jointGATKGenotyper, String rootOutputDir, MergeVCF mergeVCF,
										ANNOVCF annoVCF, boolean annotate, int memoryInMB, int wallTimeInHours,
										String baseName) {
		// TODO, change classpath
		String command = Array.toStr(PSF.Load.getAllModules(), "\n");
		command += "\njava -Xmx" + memoryInMB + "m -jar parkGATK.jar seq.analysis.GATK_Genotyper ";
		command += GATK_LanePrep.ROOT_INPUT_COMMAND + jointGATKGenotyper.getRootInputDir() + SPACE;
		command += GATK_LanePrep.ROOT_OUTPUT_COMMAND + rootOutputDir + SPACE;
		command += GATK_LanePrep.REFERENCE_GENOME_COMMAND + gatk.getReferenceGenomeFasta() + SPACE;
		command += GATK.GATK_LOCATION_COMMAND + gatk.getGATKLocation() + SPACE;
		command += NUM_THREADS + numWithinSampleThreads + SPACE;
		command += OUTPUT_COMMAND + jointGATKGenotyper.getOutput() + SPACE;
		command += GATK_LanePrep.LOG_FILE_COMMAND
								+ ext.removeDirectoryInfo(jointGATKGenotyper.getLog().getFilename()) + SPACE;
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
			command += SNPSIFT.SNP_SIFT_LOCATION_COMMAND + snpsift.getSnpSiftLocation() + SPACE;
			command += ANNOVAR.ANNOVAR_COMMAND + annovar.getAnnovarLocation() + SPACE;
		} else {
			command += SNPEFF.SNP_EFF_NO_ANNO_COMMAND + SPACE;
		}
		if (mergeVCF != null) {
			command += MERGE_WITH + mergeVCF.getArg() + SPACE;
		}
		if (annoVCF != null) {
			command += EXTRA_VCF_ANNOTATIONS + annoVCF.getArg() + SPACE;
		}
		if (gatk.getRegionsFile() != null) {
			command += REGIONS_FILE + gatk.getRegionsFile();
		}
		Files.qsub("GATK_Genotype_"+ baseName, command, memoryInMB, wallTimeInHours,
								numWithinSampleThreads);
	}

	public void runSingleSampleAllSites(String[] inputBams) {

		if (!isFail() && Files.checkAllFiles("", inputBams, verbose, log)) {
			if (inputBams != null) {
				siSampleHaplotypeCallers = new GATK.SingleSampleHaplotypeCaller[inputBams.length];
				int[] actualWithinSampleThreads = optimizeThreads(inputBams.length, numBetweenSampleThreads,
																													numWithinSampleThreads, log);
				WorkerHive<GATK.SingleSampleHaplotypeCaller> hive =
																													new WorkerHive<GATK.SingleSampleHaplotypeCaller>(	numBetweenSampleThreads,
																																																						10,
																																																						log);
				WorkerSingleSampleAllSites[] workers = new WorkerSingleSampleAllSites[inputBams.length];
				for (int i = 0; i < inputBams.length; i++) {
					Logger altLog = new Logger(ext.rootOf(inputBams[i], false) + ".HC_ERC.log");
					workers[i] = new WorkerSingleSampleAllSites(gatk, inputBams[i], ext.rootOf(inputBams[i]),
																											actualWithinSampleThreads[i], altLog);
				}
				hive.addCallables(workers);
				hive.execute(true);
				siSampleHaplotypeCallers = hive	.getResults()
																				.toArray(new GATK.SingleSampleHaplotypeCaller[hive.getResults()
																																													.size()]);
				for (int i = 0; i < workers.length; i++) {
					if (siSampleHaplotypeCallers[i].isFail()) {
						log.reportTimeError("Failed single sample haplotype calling for "
																+ siSampleHaplotypeCallers[i].getInputBam());
						fail = true;
					}
				}
			} else {
				// TODO better check
			}
		}
	}

	private static int[] optimizeThreads(	int numInputs, int numBetweenSampleThreads,
																				int numWithinSampleThreads, Logger log) {
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
			log.report(ext.getTime()
									+ " Info - since we have extra between sample threads, we will allocate some to within sample(s) threads");
			log.report(ext.getTime()+ " Info - allocated " + (numBetweenSampleThreads - numInputs)
									+ " extra thread(s) within sample(s) as follows: "
									+ Array.toStr(optimizedWithin));
		}
		return optimizedWithin;
	}

	private static class WorkerSingleSampleAllSites	implements
																									Callable<GATK.SingleSampleHaplotypeCaller> {
		private final GATK GATK;
		private final String inputBam, baseId;
		private final int numWithinSampleThreads;
		private final Logger altLog;

		public WorkerSingleSampleAllSites(org.genvisis.seq.analysis.GATK gATK, String inputBam,
																			String baseId, int numWithinSampleThreads, Logger altLog) {
			super();
			GATK = gATK;
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

		private final String rootInputDir, rootOutputDir;
		private String output;
		private String rawVCF;
		private String fileOfGVCFs, recalSNP_VCF_File, recalSNP_Indel_VCF_File;
		private String recalSNPFile, tranchesSNPFile, rscriptSNPFile;
		private String recalINDELFile, tranchesINDELFile, rscriptINDELFile;
		private String restrictionContig;
		private String[] inputGVCFs;
		private boolean fail;
		private final Logger log;

		/**
		 * @param rootInputDir
		 * @param rootOutputDir
		 * @param output
		 * @param restrictionContig can be null
		 * @param log
		 */
		public JointGATKGenotyper(String rootInputDir, String rootOutputDir, String output,
															String restrictionContig, Logger log) {
			super();
			this.rootInputDir = rootInputDir;
			this.rootOutputDir = rootOutputDir;
			this.output = output;
			this.restrictionContig = restrictionContig;
			fileOfGVCFs = null;
			this.log = log;
			fail = false;
		}

		public String getRestrictionContig() {
			return restrictionContig;
		}

		public void setOutput(String output) {
			this.output = output;
		}

		public void init(String fileOfGVCFs) {
			this.fileOfGVCFs = fileOfGVCFs;
			initOutputs();
			if (fileOfGVCFs != null) {
				log.report(ext.getTime()+ " Info - using GVCF files listed in the first column of"
										+ fileOfGVCFs);
				inputGVCFs = HashVec.loadFileToStringArray(fileOfGVCFs, false, new int[] {0}, true);
			} else if (rootInputDir == null) {
				log.reportError("Error - a file listing GVCF files was not provided and the root input directory was not provided, halting...");
				fail = true;
			} else {
				log.report(ext.getTime()+ " Info - finding files with extension " + GATK.GVCF + " in "
										+ rootInputDir);
				inputGVCFs = Files.toFullPaths(Files.list(rootInputDir, GATK.GVCF, false), rootInputDir);

			}
			if (inputGVCFs == null || inputGVCFs.length < 1) {
				log.reportError("Error - could not find any GVCF files to joint genotype");
				fail = true;
			} else {
				log.report(ext.getTime()+ " Info - using " + inputGVCFs.length
										+ " file(s) for joint Genotyping");
			}

		}

		public void initOutputs() {
			String currentRoot = rootOutputDir + ext.rootOf(output);
			rawVCF = currentRoot + GATK.VCF;

			recalSNPFile = currentRoot + "." + GATK.SNP + RECAL_EXT;
			tranchesSNPFile = currentRoot + "." + GATK.SNP + TRANCHES_EXT;
			rscriptSNPFile = currentRoot + "." + GATK.SNP + RScript_EXT;

			recalINDELFile = currentRoot + "." + GATK.INDEL + RECAL_EXT;
			tranchesINDELFile = currentRoot + "." + GATK.INDEL + TRANCHES_EXT;
			rscriptINDELFile = currentRoot + "." + GATK.INDEL + RScript_EXT;

			recalSNP_VCF_File = ext.addToRoot(rawVCF, "." + GATK.SNP + RECAL_EXT);
			recalSNP_Indel_VCF_File = ext.addToRoot(recalSNP_VCF_File, "." + GATK.INDEL + RECAL_EXT);
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

	public static void jointGenotype(	String rootInputDir, String rootOutputDir, String output,
																		String gATKLocation, String referenceGenomeFasta,
																		String fileOfGVCFs, String hapMapTraining, String omniTraining,
																		String thousandGTraining, String dbSnpTraining,
																		String millsIndelTraining, String snpEffLocation,
																		String snpSiftLocation, String annovarLocation,
																		String annoBuild, String regionsFile, MergeVCF mergeVCF,
																		ANNOVCF annoVCF, String restrictionContig, boolean verbose,
																		boolean overwriteExisting, boolean batch, boolean annotate,
																		boolean skipRecalibration, int numThreads, int memoryInMB,
																		int wallTimeInHours, Logger log) {
		GATK gatk = new GATK(	gATKLocation, referenceGenomeFasta, null, null, null, verbose,
													overwriteExisting, log);
		gatk.setDbSnpTraining(dbSnpTraining);
		gatk.setHapMapTraining(hapMapTraining);
		gatk.setOmniTraining(omniTraining);
		gatk.setThousandGTraining(thousandGTraining);
		gatk.setMillsIndelTraining(millsIndelTraining);
		gatk.setRegionsFile(regionsFile);
		SNPEFF snpeff = new SNPEFF(snpEffLocation, verbose, overwriteExisting, log);
		SNPSIFT snpsift = new SNPSIFT(snpSiftLocation, verbose, overwriteExisting, log);
		ANNOVAR annovar = new ANNOVAR(annovarLocation, verbose, overwriteExisting, log);
		GATK_Genotyper genotyper = new GATK_Genotyper(gatk, snpeff, snpsift, annovar, 0, numThreads,
																									verbose, log);
		JointGATKGenotyper jGatkGenotyper = new JointGATKGenotyper(	rootInputDir, rootOutputDir, output,
																																restrictionContig, log);
		jGatkGenotyper.init(fileOfGVCFs);
		new File(rootOutputDir).mkdirs();
		if (batch) {
			genotyper.batch(jGatkGenotyper, rootOutputDir, mergeVCF, annoVCF, annotate, memoryInMB,
											wallTimeInHours, output);
		} else {
			boolean progress = genotyper.runJointGenotyping(jGatkGenotyper);
			if (progress) {
				if (gatk.getRegionsFile() != null) {
					log.reportTimeInfo("Subsetting the vcf to regions defined in "+ regionsFile
															+ " with a 100 bp buffer");
					String subsetVcf = VCFOps.extractSegments(jGatkGenotyper.getRawVCF(), regionsFile, 100,
																										null, rootOutputDir, false, false, 1, log);
					String newOutput = ext.removeDirectoryInfo(subsetVcf);
					jGatkGenotyper.setOutput(newOutput);
					jGatkGenotyper.initOutputs();// now starting with the subset vcf
				}
				if (skipRecalibration) {
					if (annotate) {
						genotyper.annotateVCF(jGatkGenotyper.output, annoBuild, mergeVCF, annoVCF);
					}
				} else {
					progress = genotyper.runRecalibration(jGatkGenotyper);
					if (progress) {
						if (progress && annotate) {
							genotyper.annotateVCF(jGatkGenotyper.getRecalSNP_Indel_VCF_File(), annoBuild,
																		mergeVCF, annoVCF);
							genotyper.determineTsTV(jGatkGenotyper.getRecalSNP_Indel_VCF_File());
						}
					}
				}
			}
		}
	}

	public static String annotateOnlyWithDefualtLocations(String vcf, ANNOVCF annoVCF,
																												boolean verbose, boolean overwriteExisting,
																												Logger log) {
		return annotateOnly(vcf, AnnotationDefaults.GATK_LOC, AnnotationDefaults.REF, null, null, null,
												null, null, null, AnnotationDefaults.SNP_EFF, null,
												AnnotationDefaults.ANNOVAR, SNPEFF.BUILDS[0], annoVCF, verbose,
												overwriteExisting, log);
	}

	public static String annotateOnly(String vcf, String gATKLocation, String referenceGenomeFasta,
																		String fileOfGVCFs, String hapMapTraining, String omniTraining,
																		String thousandGTraining, String dbSnpTraining,
																		String millsIndelTraining, String snpEffLocation,
																		String snpSiftLocation, String annovarLocation,
																		String annoBuild, ANNOVCF annoVCF, boolean verbose,
																		boolean overwriteExisting, Logger log) {
		String snpSiftLoc = snpSiftLocation;
		if (snpSiftLoc == null || snpSiftLoc.equals(PSF.Ext.BLANK)) {
			snpSiftLoc = snpEffLocation;
		}
		GATK gatk = new GATK(	gATKLocation, referenceGenomeFasta, null, null, null, verbose,
													overwriteExisting, log);
		SNPEFF snpeff = new SNPEFF(snpEffLocation, verbose, overwriteExisting, log);
		SNPSIFT snpsift = new SNPSIFT(snpSiftLoc, verbose, overwriteExisting, log);
		ANNOVAR annovar = new ANNOVAR(annovarLocation, verbose, overwriteExisting, log);
		GATK_Genotyper genotyper = new GATK_Genotyper(gatk, snpeff, snpsift, annovar, 0, 1, verbose,
																									log);
		return genotyper.annotateVCF(vcf, annoBuild, null, annoVCF);
	}

	public static String annotateOnly(String vcf, String gATKLocation, String referenceGenomeFasta,
																		String snpEffLocation, String snpSiftLocation,
																		String annovarLocation, String annoBuild, boolean verbose,
																		boolean overwriteExisting, Logger log) {
		return annotateOnly(vcf, gATKLocation, referenceGenomeFasta, null, null, null, null, null, null,
												snpEffLocation, snpSiftLocation, annovarLocation, annoBuild, null, verbose,
												overwriteExisting, log);
	}

	private static class MergeVCF {
		private final String tag;
		private final String[] vcfsToMergeWith;

		private MergeVCF(String tag, String[] vcfsToMergeWith) {
			super();
			this.tag = tag;
			this.vcfsToMergeWith = vcfsToMergeWith;
		}

		private String getArg() {
			return tag + ":" + Array.toStr(vcfsToMergeWith, ",");

		}

		public String getTag() {
			return tag;
		}

		public String[] getVcfsToMergeWith() {
			return vcfsToMergeWith;
		}

		private static MergeVCF fromArg(String arg) {
			try {
				String[] tmp = ext.parseStringArg(arg, "").split(":");
				String tag = tmp[0];
				String[] vcfsToMergeWith = tmp[1].split(",");

				return new MergeVCF(tag, vcfsToMergeWith);
			} catch (Exception e) {
				System.err.println("Could not process merge arg " + arg + " format is tag:vcf1,vcf2");
				return null;
			}
		}

	}

	public static class ANNOVCF {
		private final String[] annos;
		private final String vcf;
		private final String tag;

		private ANNOVCF(String tag, String[] annos, String vcf) {
			super();
			this.tag = tag;
			this.annos = annos;
			this.vcf = vcf;
		}

		public String getArg() {
			return tag + ":" + vcf + ":" + Array.toStr(annos, ",");

		}

		public String[] getAnnos() {
			return annos;
		}

		public String getVcf() {
			return vcf;
		}

		public String getTag() {
			return tag;
		}

		public static ANNOVCF fromArg(String arg) {
			try {
				String[] tmp = ext.parseStringArg(arg, "").split(":");
				String tag = tmp[0];
				String vcf = tmp[1];
				String[] annos = tmp[2].split(",");

				return new ANNOVCF(tag, annos, vcf);
			} catch (Exception e) {
				System.err.println("Could not process merge arg "+ arg
														+ " format is tag:vcf.vcf:anno1,anno2");
				return null;
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
	public static final String ANNOTATE_VCF = "annotateThisVCFOnly=";
	public static final String REGIONS_FILE = "regionsFile=";
	public static final String MERGE_WITH = "mergeWith=";
	public static final String EXTRA_VCF_ANNOTATIONS = "extraVCFAnno=";

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
		String snpEffLocation = PSF.Ext.BLANK;
		String snpSiftLocation = PSF.Ext.BLANK;
		String annovarLocation = PSF.Ext.BLANK;
		String annoBuild = SNPEFF.BUILDS[0];
		String vcfToAnnotate = null;
		String restrictionContig = null;
		boolean annotate = true;
		String regionsFile = null;
		String logFile = "GATK_GENOTYPE.log";
		boolean skipRecalibration = false;
		MergeVCF mergeVCF = null;
		ANNOVCF annoVCF = null;
		String usage = "\n" + "seq.GATK_Genotyper requires 2 argument\n";
		usage += "   (1) root input directory (i.e. "+ GATK_LanePrep.ROOT_INPUT_COMMAND + rootInputDir
							+ " (no default))\n" + "";
		usage += "   (2) root output directory (i.e. "+ GATK_LanePrep.ROOT_OUTPUT_COMMAND
							+ rootOutputDir + " (no default))\n" + "";
		usage += "   (3) tab-delimited file with no header of (i.e. "+ FILE_OF_GVCFS + fileOfGVCFs
							+ " (optional, no default))\n" + "";
		usage += "   (4) the full path to a  reference genome in fasta format (i.e."
							+ GATK_LanePrep.REFERENCE_GENOME_COMMAND + referenceGenomeFasta + " (no default))\n"
							+ "";
		usage += "   (5) the full path to the GATK executable (i.e. "+ GATK.GATK_LOCATION_COMMAND
							+ gATKLocation + " (defaults to systems path))\n" + "";

		usage += "   (6) run in quiet mode (i.e. "+ GATK_LanePrep.QUIET_COMMAND
							+ " (not tbe default))\n" + "";
		usage += "   (7) number of  threads for analysis(i.e."+ NUM_THREADS + numThreads
							+ " (default))\n" + "";

		usage += "   (8) filename for a log (i.e. "+ GATK_LanePrep.LOG_FILE_COMMAND + logFile
							+ " (default))\n" + "";
		usage += "   (9) over-write exsiting files (i.e. "+ GATK_LanePrep.OVERWRITE_EXISTING_COMMAND
							+ " (not the default))\n" + "";
		usage += "   (10) set up a batch analysis for the root input directory for a log (i.e. "
							+ GATK_LanePrep.BATCH_COMMAND + " (not the default))\n" + "";
		usage += "   (11) root output for analysis (i.e. "+ OUTPUT_COMMAND + output + " ( default))\n"
							+ "";
		usage += "   (12) HapMap SNP Training Referenece (i.e. "+ HAPMAP_TRAIN_COMMAND
							+ " ( no default))\n" + "";
		usage += "   (13) Omni SNP Training Referenece (i.e. "+ OMNI_TRAIN_COMMAND
							+ " ( no default))\n" + "";
		usage += "   (14) 1000 SNP genomes Training Referenece (i.e. "+ G1000_COMMAND
							+ " ( no default))\n" + "";
		usage += "   (15) dbSNP SNP Training Referenece (i.e. "+ DBSNP_COMMMAND + " ( no default))\n"
							+ "";
		usage += "   (16) mills INDEL Training Referenece (i.e. " + MILLS + " ( no default))\n" + "";
		usage += "   (17) full path to the SNP EFF directory (i.e. "+ SNPEFF.SNP_EFF_COMMAND
							+ " ( no default))\n" + "";
		usage += "   (18) the build version for SNP EFF annotation (options are "
							+ Array.toStr(SNPEFF.BUILDS, ", ") + " (i.e. " + SNPEFF.SNP_EFF_BUILD_COMMAND
							+ annoBuild + " ( default))\n" + "";
		usage +=
					"   (19) full path to the SNP SIFT directory (only if different from the SNP EFF directory) (i.e. "
							+ SNPSIFT.SNP_SIFT_LOCATION_COMMAND + " ( no default))\n" + "";
		usage += "   (20) do not annotate with SNP EFF/SIFT/ANNOVAR (i.e. "
							+ SNPEFF.SNP_EFF_NO_ANNO_COMMAND + " ( not the default))\n" + "";
		usage += "   (21) full path to the ANNOVAR directory (i.e. "+ ANNOVAR.ANNOVAR_COMMAND
							+ " ( no default))\n" + "";
		usage +=
					"   (23) full path to a file for restricting the vcf...and subsequent recalibrations (i.e. "
							+ REGIONS_FILE + " ( no default))\n" + "";

		usage += "   (24) annotate this vcf only (skipping all previous steps) (i.e. "+ ANNOTATE_VCF
							+ " ( no default))\n" + "";
		usage +=
					"   (25) merge in these vcfs prior to annotating( like ARIC) (ex. MERGE_WITH=ARIC:aric1.vcf,aric2.vcf) ( (i.e. "
							+ MERGE_WITH + " ( no default))\n" + "";
		usage +=
					"   (26) use another vcf to add more annotations (ex. EXTRA_VCF_ANNOTATIONS=charge.vcf:charge.maf1,charge.maf2) ( (i.e. "
							+ EXTRA_VCF_ANNOTATIONS + " ( no default))\n" + "";

		usage += "   (27) restrict genotyping to a specific contig ( (i.e. "+ "chrM"
							+ " ( no default))\n" + "";
		usage +=
					"   (28) since targeted sequencing, or mtDNA genotyping will generally fail variant recalibration, skip it and see caveats here https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php ( (i.e. "
							+ "-skipRecal" + " ( not the default))\n" + "";

		// usage += " (18) log file name (i.e. " + MILLS + " ( no default))\n" + "";

		for (String arg : args) {

			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith(HAPMAP_TRAIN_COMMAND)) {
				hapMapTraining = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(OMNI_TRAIN_COMMAND)) {
				omniTraining = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(G1000_COMMAND)) {
				thousandGTraining = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(DBSNP_COMMMAND)) {
				dbSnpTraining = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(MILLS)) {
				millsIndelTraining = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(GATK_LanePrep.ROOT_INPUT_COMMAND)) {
				rootInputDir = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(GATK_LanePrep.ROOT_OUTPUT_COMMAND)) {
				rootOutputDir = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(FILE_OF_GVCFS)) {
				fileOfGVCFs = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(GATK_LanePrep.REFERENCE_GENOME_COMMAND)) {
				referenceGenomeFasta = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(GATK_LanePrep.LOG_FILE_COMMAND)) {
				logFile = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(NUM_THREADS)) {
				numThreads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("memoryInMB=")) {
				memoryInMB = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("wallTimeInHours=")) {
				wallTimeInHours = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(GATK_LanePrep.QUIET_COMMAND)) {
				verbose = false;
				numArgs--;
			} else if (arg.startsWith(GATK_LanePrep.BATCH_COMMAND)) {
				batch = true;
				numArgs--;
			} else if (arg.startsWith(GATK_LanePrep.OVERWRITE_EXISTING_COMMAND)) {
				batch = true;
				numArgs--;
			} else if (arg.startsWith(GATK.GATK_LOCATION_COMMAND)) {
				gATKLocation = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(SNPEFF.SNP_EFF_COMMAND)) {
				snpEffLocation = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(SNPSIFT.SNP_SIFT_LOCATION_COMMAND)) {
				snpSiftLocation = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(SNPEFF.SNP_EFF_NO_ANNO_COMMAND)) {
				annotate = false;
				numArgs--;
			} else if (arg.startsWith("-skipRecal")) {
				skipRecalibration = true;
				numArgs--;
			} else if (arg.startsWith(OUTPUT_COMMAND)) {
				output = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(ANNOVAR.ANNOVAR_COMMAND)) {
				annovarLocation = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(ANNOTATE_VCF)) {
				vcfToAnnotate = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(GATK.L)) {
				restrictionContig = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(REGIONS_FILE)) {
				regionsFile = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(MERGE_WITH)) {
				mergeVCF = MergeVCF.fromArg(arg);
				numArgs--;
			} else if (arg.startsWith(EXTRA_VCF_ANNOTATIONS)) {
				annoVCF = ANNOVCF.fromArg(arg);
				numArgs--;
			}

			else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.exit(1);
		}
		if (rootOutputDir != null) {
			new File(rootOutputDir).mkdirs();
		}
		Logger log = new Logger(rootOutputDir + logFile);
		if (snpSiftLocation.equals(PSF.Ext.BLANK)) {
			snpSiftLocation = snpEffLocation;
		}
		if (vcfToAnnotate != null) {
			log.reportTimeInfo("Attempting to annotate " + vcfToAnnotate);
			annotateOnly(	vcfToAnnotate, gATKLocation, referenceGenomeFasta, fileOfGVCFs, hapMapTraining,
										omniTraining, thousandGTraining, dbSnpTraining, millsIndelTraining,
										snpEffLocation, snpSiftLocation, annovarLocation, annoBuild, null, verbose,
										overwriteExisting, log);
		} else {
			jointGenotype(rootInputDir, rootOutputDir, output, gATKLocation, referenceGenomeFasta,
										fileOfGVCFs, hapMapTraining, omniTraining, thousandGTraining, dbSnpTraining,
										millsIndelTraining, snpEffLocation, snpSiftLocation, annovarLocation, annoBuild,
										regionsFile, mergeVCF, annoVCF, restrictionContig, verbose, overwriteExisting,
										batch, annotate, skipRecalibration, numThreads, memoryInMB, wallTimeInHours,
										log);
		}
	}
}
