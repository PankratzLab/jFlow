package org.genvisis.gwas;

import java.io.File;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.List;
import java.util.Vector;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

public class Qc {

	public static final String QC_DIR = "quality_control/";

	public static final String MARKER_QC_DIR = "marker_qc/";
	public static final String SAMPLE_QC_DIR = "sample_qc/";
	public static final String LD_PRUNING_DIR = "ld_pruning/";
	public static final String GENOME_DIR = "genome/";
	public static final String ANCESTRY_DIR = "ancestry/";
	public static final String FURTHER_ANALYSIS_DIR = "further_analysis_QC/";

	public static final String FURTHER_ANALYSIS_QC_PLINK_SUFFIX = "_QCd";

	/** A rough listing of the Folders created by fullGamut */
	public static String[] FOLDERS_CREATED = {QC_DIR + MARKER_QC_DIR, QC_DIR + SAMPLE_QC_DIR,
																						QC_DIR + LD_PRUNING_DIR,
																						QC_DIR + GENOME_DIR, QC_DIR + ANCESTRY_DIR,
																						QC_DIR + FURTHER_ANALYSIS_DIR};
	/** A rough listing of the files created, by folder, by fullGamut */
	// TODO: This does not accommodate cases where the plinkroot is something other than "plink"
	// Also ought to be automated...
	public static String[][] FILES_CREATED = {{"plink.bed", "freq.frq", "missing.imiss",
																						 /* "test.missing.missing", *//*
																																					 * not actually necessary
																																					 */ "hardy.hwe",
																						 "mishap.missing.hap", "gender.assoc", "gender.missing",
																						 "miss_drops.dat"},
																						{"plink.bed", "missing.imiss"},
																						{"plink.bed", "plink.prune.in"},
																						{"plink.bed", "plink.genome", "plink.genome_keep.dat"},
																						{"plink.bed", "unrelateds.txt"},
																						PSF.Plink.getPlinkBedBimFam("plink"
																																				+ FURTHER_ANALYSIS_QC_PLINK_SUFFIX)};

	private final String dir;
	private final String plink;
	private final Logger log;



	/**
	 * @param dir Directory with plink files to run from
	 * @param plinkPrefix prefix of plink binaries
	 * @param log
	 */
	public Qc(String dir, String plinkPrefix, Logger log) {
		super();
		dir = ext.verifyDirFormat(dir) + QC_DIR;
		if (!dir.startsWith("/") && !dir.contains(":")) {
			dir = ext.verifyDirFormat((new File("./" + dir)).getAbsolutePath());
		}
		this.dir = dir;
		this.plink = plinkPrefix == null ? "plink" : plinkPrefix;
		this.log = log;
	}

	private boolean markerQc(String subDir, final Collection<String> markerQcParams) {
		subDir = ext.verifyDirFormat(subDir);
		new File(dir + subDir).mkdirs();
		if (!Files.exists(dir + subDir + plink + ".bed")) {
			log.report(ext.getTime() + "]\tRunning initial --make-bed");
			CmdLine.runDefaults("plink2 --bfile ../../" + plink + " --make-bed --noweb --out ./" + plink,
													dir + subDir, log);
		}
		if (!Files.exists(dir + subDir + plink + ".bed")) {
			System.err.println("Error - CRITICAL ERROR, creating initial PLINK files failed.");
			return false;
		}
		String geno20 = plink + "_geno20";
		if (!Files.exists(dir + subDir + geno20 + ".bed")) {
			log.report(ext.getTime() + "]\tRunning --geno 0.2");
			CmdLine.runDefaults("plink2 --bfile " + plink + " --geno 0.2 --make-bed --noweb --out ./"
													+ geno20, dir + subDir, log);
		}
		PSF.checkInterrupted();
		String geno20mind10 = geno20 + "_mind10";
		if (!Files.exists(dir + subDir + geno20mind10 + ".bed")) {
			log.report(ext.getTime() + "]\tRunning --mind 0.1");
			CmdLine.runDefaults("plink2 --bfile " + geno20 + " --mind 0.1 --make-bed --noweb --out ./"
													+ geno20mind10, dir + subDir, log);
		}
		PSF.checkInterrupted();
		String mind10 = plink + "_mind10";
		if (!Files.exists(dir + subDir + mind10 + ".bed")) {
			log.report(ext.getTime() + "]\tRemoving trimmed samples from --mind 0.1");
			CmdLine.runDefaults("plink2 --bfile " + plink + " --keep " + geno20mind10
													+ ".fam --make-bed --noweb --out ./" + mind10, dir + subDir, log);
		}
		PSF.checkInterrupted();
		if (!Files.exists(dir + subDir + "freq.frq")) {
			log.report(ext.getTime() + "]\tRunning --freq");
			CmdLine.runDefaults("plink2 --bfile " + mind10
													+ " --geno 1 --mind 1 --freq --out freq --noweb", dir + subDir, log);
		}
		PSF.checkInterrupted();
		if (!Files.exists(dir + subDir + "missing.imiss")) {
			log.report(ext.getTime() + "]\tRunning --missing");
			CmdLine.runDefaults("plink2 --bfile " + mind10
													+ " --geno 1 --mind 1 --missing --out missing --noweb", dir + subDir,
													log);
		}
		PSF.checkInterrupted();
		if (!Files.exists(dir + subDir + "test.missing.missing")) {
			log.report(ext.getTime() + "]\tRunning --test-missing");
			CmdLine.runDefaults("plink2 --bfile " + mind10
													+ " --geno 1 --mind 1 --test-missing --out test.missing --noweb",
													dir + subDir, log);
		}
		PSF.checkInterrupted();
		if (!Files.exists(dir + subDir + "hardy.hwe")) {
			log.report(ext.getTime() + "]\tRunning --hardy");
			CmdLine.runDefaults("plink2 --bfile " + mind10
													+ " --geno 1 --mind 1 --hardy --out hardy --noweb", dir + subDir, log);
		}
		PSF.checkInterrupted();
		if (!Files.exists(dir + subDir + "mishap.missing.hap")) {
			log.report(ext.getTime() + "]\tRunning --test-mishap");
			CmdLine.runDefaults("plink2 --bfile " + geno20mind10
													+ " --geno 1 --mind 1 --test-mishap --out mishap --noweb", dir + subDir,
													log);
		}
		PSF.checkInterrupted();
		if (!Files.exists(dir + subDir + "gender.assoc")) {
			log.report(ext.getTime() + "]\tRunning --assoc gender");
			CmdLine.runDefaults("plink2 --bfile " + mind10 + " --geno 1 --mind 1 --pheno " + mind10
													+ ".fam --mpheno 3 --assoc --out gender --noweb", dir + subDir, log);
		}
		PSF.checkInterrupted();
		if (!Files.exists(dir + subDir + "gender.missing")) {
			log.report(ext.getTime() + "]\tRunning --test-missing gender");
			CmdLine.runDefaults("plink2 --bfile " + mind10 + " --geno 1 --mind 1 --pheno " + mind10
													+ ".fam --mpheno 3 --test-missing --out gender --noweb", dir + subDir,
													log);
		}
		PSF.checkInterrupted();
		if (!Files.exists(dir + subDir + "miss_drops.dat")) {
			List<String> writeToFile = Lists.newArrayList();
			writeToFile.add(MarkerQC.COMMAND);
			writeToFile.add("dir=" + dir + subDir);
			writeToFile.addAll(markerQcParams);
			Files.writeIterable(writeToFile, dir + subDir + "miss.crf");
			int runCode = MarkerQC.parseParameters(dir + subDir + "miss.crf", log, false);
			if (runCode != 0) {
				log.reportError("Failed to perform marker QC with " + dir + subDir + "miss.crf");
				return false;
			}
		}
		PSF.checkInterrupted();
		return true;
	}

	private void runRelationAncestryQC(boolean keepGenomeInfoForRelatedsOnly) {
		if (!markerQc(MARKER_QC_DIR, MarkerQC.DEFAULT_STRICT_QC_THRESHOLDS)) {
			return;
		}

		new File(dir + SAMPLE_QC_DIR).mkdirs();
		if (!Files.exists(dir + SAMPLE_QC_DIR + plink + ".bed")) {
			log.report(ext.getTime() + "]\tRunning --exclude miss_drops.dat");
			CmdLine.runDefaults("plink2 --bfile ../" + MARKER_QC_DIR + plink
													+ " --exclude ../" + MARKER_QC_DIR
													+ "miss_drops.dat --make-bed --noweb --out "
													+ plink, dir + SAMPLE_QC_DIR, log);
		}
		PSF.checkInterrupted();
		if (!Files.exists(dir + SAMPLE_QC_DIR + "missing.imiss")) {
			log.report(ext.getTime() + "]\tRunning --missing");
			CmdLine.runDefaults("plink2 --bfile " + plink
													+ " --geno 1 --mind 1 --missing --out missing --noweb",
													dir + SAMPLE_QC_DIR, log);
		}
		PSF.checkInterrupted();

		new File(dir + LD_PRUNING_DIR).mkdirs();
		if (!Files.exists(dir + LD_PRUNING_DIR + plink + ".bed")) {
			log.report(ext.getTime()
								 + "]\tRunning --mind 0.05 (removes samples with callrate <95% for the markers that did pass QC)");
			CmdLine.runDefaults("plink2 --bfile ../" + SAMPLE_QC_DIR + plink
													+ " --mind 0.05 --make-bed --noweb --out " + plink, dir + LD_PRUNING_DIR,
													log);
		}
		PSF.checkInterrupted();
		if (!Files.exists(dir + LD_PRUNING_DIR + plink + ".prune.in")) {
			log.report(ext.getTime() + "]\tRunning --indep-pairwise 50 5 0.3");
			CmdLine.runDefaults("plink2 --noweb --bfile " + plink + " --indep-pairwise 50 5 0.3 --out "
													+ plink, dir + LD_PRUNING_DIR, log);
		}
		PSF.checkInterrupted();

		new File(dir + GENOME_DIR).mkdirs();
		if (!Files.exists(dir + GENOME_DIR + plink + ".bed")) {
			log.report(ext.getTime() + "]\tRunning --extract " + plink + ".prune.in");
			CmdLine.runDefaults("plink2 --bfile ../" + LD_PRUNING_DIR + plink + " --extract ../"
													+ LD_PRUNING_DIR
													+ plink + ".prune.in --make-bed --noweb --out " + plink, dir + GENOME_DIR,
													log);
		}
		PSF.checkInterrupted();
		if (!Files.exists(dir + GENOME_DIR + plink + ".genome")) {
			log.report(ext.getTime() + "]\tRunning --genome"
								 + (keepGenomeInfoForRelatedsOnly ? " --min 0.1" : ""));
			CmdLine.runDefaults("plink2 --noweb --bfile " + plink + " --genome"
													+ (keepGenomeInfoForRelatedsOnly ? " --min 0.1" : "") + " --out " + plink,
													dir + GENOME_DIR, log);
		}
		PSF.checkInterrupted();
		if (!keepGenomeInfoForRelatedsOnly && !Files.exists(dir + GENOME_DIR + "mds20.mds")) {
			log.report(ext.getTime() + "]\tRunning --mds-plot 20");
			CmdLine.runDefaults("plink2 --bfile " + plink + " --read-genome " + plink
													+ ".genome --cluster --mds-plot 20 --out mds20 --noweb", dir + GENOME_DIR,
													log);
		}
		PSF.checkInterrupted();
		if (!Files.exists(dir + GENOME_DIR + plink + ".genome_keep.dat")) {
			log.report(ext.getTime() + "]\tRunning flagRelateds");
			String lrrFile = dir + GENOME_DIR + "lrr_sd.xln";
			if (!Files.exists(lrrFile)) {
				lrrFile = dir + "../../lrr_sd.xln";
			}
			Plink.flagRelateds(dir + GENOME_DIR + plink + ".genome", dir + GENOME_DIR + plink + ".fam",
												 dir + MARKER_QC_DIR + "missing.imiss", lrrFile, Plink.FLAGS,
												 Plink.THRESHOLDS, 4,
												 false);
		}
		PSF.checkInterrupted();

		new File(dir + ANCESTRY_DIR).mkdirs();
		if (!Files.exists(dir + ANCESTRY_DIR + "unrelateds.txt")) {
			log.report(ext.getTime() + "]\tCopying " + GENOME_DIR + plink
								 + ".genome_keep.dat to " + ANCESTRY_DIR + "unrelateds.txt");
			Files.copyFile(dir + GENOME_DIR + plink + ".genome_keep.dat",
										 dir + ANCESTRY_DIR + "unrelateds.txt");
		}
		PSF.checkInterrupted();
		if (!Files.exists(dir + ANCESTRY_DIR + plink + ".bed")) {
			log.report(ext.getTime() + "]\tRunning --extract " + plink
								 + ".prune.in (again, this time to " + ANCESTRY_DIR + ")");
			CmdLine.runDefaults("plink2 --bfile ../" + GENOME_DIR + plink + " --make-bed --noweb --out "
													+ plink, dir + ANCESTRY_DIR, log);
		}
		PSF.checkInterrupted();

	}

	private void runFurtherAnalysisQC() {
		final String subDir = FURTHER_ANALYSIS_DIR;
		final String plinkQCd = plink + FURTHER_ANALYSIS_QC_PLINK_SUFFIX;
		if (!markerQc(subDir, MarkerQC.DEFAULT_GWAS_QC_THRESHOLDS))
			return;
		List<String> applyQCCommand = ImmutableList.of("plink2", "--noweb", "--bfile", plink,
																									 "--exclude", "miss_drops.dat", "--mind", "0.05",
																									 "--make-bed", "--out", plinkQCd);
		List<String> requiredOutputs = Lists.newArrayList(PSF.Plink.getPlinkBedBimFam(plinkQCd));
		requiredOutputs.add("miss_drops.dat");
		requiredOutputs = Collections.unmodifiableList(requiredOutputs);
		List<String> requiredInputs = ImmutableList.copyOf(PSF.Plink.getPlinkBedBimFam(plink));
		CmdLine.runCommandWithFileChecks(applyQCCommand, dir + subDir, requiredInputs, requiredOutputs,
																		 true, false, true, log);
	}

	public void run() {
		run(false);
	}

	public void run(boolean keepGenomeInfoForRelatedsOnly) {

		long time = new Date().getTime();

		runRelationAncestryQC(keepGenomeInfoForRelatedsOnly);
		runFurtherAnalysisQC();



		System.out.println("Finished this round in " + ext.getTimeElapsed(time));
	}

	/**
	 * 
	 * @param keepGenomeInfoForRelatedsOnly true to save disk usage if unrelated genome info is not
	 *        required
	 * 
	 * @return full path to plinkroot of QC'd plink dataset
	 */
	public static void fullGamut(String dir, String plinkPrefix,
															 boolean keepGenomeInfoForRelatedsOnly, Logger log) {
		new Qc(dir, plinkPrefix, log).run(keepGenomeInfoForRelatedsOnly);
	}

	public static void fromParameters(String filename, Logger log) {
		Vector<String> params;

		params = Files.parseControlFile(filename, "gwas.Qc",
																		new String[] {"dir=./",
																									"# Make keepGenomeInfoForRelatedsOnly=false if the sample size is small and you want to run MDS plot as well",
																									"keepGenomeInfoForRelatedsOnly=true"},
																		log);

		if (params != null) {
			params.add("log=" + log.getFilename());
			main(ArrayUtils.toStringArray(params));
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = "./";
		boolean keepGenomeInfoForRelatedsOnly = true;
		String logfile = null;
		Logger log;

		String usage = "\n" + "gwas.Qc requires 0-1 arguments\n"
									 + "   (1) directory with plink.* files (i.e. dir=" + dir + " (default))\n"
									 + "   (2) if no MDS will be run, smaller file (i.e. keepGenomeInfoForRelatedsOnly="
									 + keepGenomeInfoForRelatedsOnly + " (default))\n";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("dir=")) {
				dir = ext.parseStringArg(arg, "./");
				numArgs--;
			} else if (arg.startsWith("log=")) {
				logfile = ext.parseStringArg(arg, null);
				numArgs--;
			} else if (arg.startsWith("keepGenomeInfoForRelatedsOnly=")) {
				keepGenomeInfoForRelatedsOnly = ext.parseBooleanArg(arg);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		if (logfile == null) {
			log = new Logger(dir + "fullGamutOfMarkerAndSampleQC.log");
		} else {
			log = new Logger(logfile);
		}
		try {
			fullGamut(dir, null, keepGenomeInfoForRelatedsOnly, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
