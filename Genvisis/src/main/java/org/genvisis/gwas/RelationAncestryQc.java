package org.genvisis.gwas;

import java.io.File;
import java.util.Collections;
import java.util.Date;
import java.util.EnumSet;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import org.genvisis.CLI;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.gwas.MarkerQC.QC_METRIC;

import com.google.common.collect.Maps;

public class RelationAncestryQc extends Qc {

	public static final String MARKER_QC_DIR = "marker_qc/";
	public static final String SAMPLE_QC_DIR = "sample_qc/";
	public static final String LD_PRUNING_DIR = "ld_pruning/";
	public static final String GENOME_DIR = "genome/";
	public static final String ANCESTRY_DIR = "ancestry/";
	/** A rough listing of the Folders created by fullGamut */
	public static String[] FOLDERS_CREATED = {Qc.QC_DIR + MARKER_QC_DIR, Qc.QC_DIR + SAMPLE_QC_DIR,
																						Qc.QC_DIR + LD_PRUNING_DIR,
																						Qc.QC_DIR + GENOME_DIR, Qc.QC_DIR + ANCESTRY_DIR};
	/** A rough listing of the files created, by folder, by fullGamut */
	// TODO: This does not accommodate cases where the plinkroot is something other than
	// Qc.DEFAULT_PLINKROOT
	// Also ought to be automated...
	public static String[][] FILES_CREATED = {{Qc.DEFAULT_PLINKROOT + ".bed", "freq.frq",
																						 "missing.imiss",
																						 /* "test.missing.missing", *//*
																																					 * not actually necessary
																																					 */ "hardy.hwe",
																						 "mishap.missing.hap", "gender.assoc", "gender.missing",
																						 "miss_drops.dat"},
																						{Qc.DEFAULT_PLINKROOT + ".bed", "missing.imiss"},
																						{Qc.DEFAULT_PLINKROOT + ".bed",
																						 Qc.DEFAULT_PLINKROOT + ".prune.in"},
																						{Qc.DEFAULT_PLINKROOT + ".bed",
																						 Qc.DEFAULT_PLINKROOT + ".genome",
																						 Qc.DEFAULT_PLINKROOT + ".genome_keep.dat"},
																						{Qc.DEFAULT_PLINKROOT + ".bed", "unrelateds.txt"}};
	public static final String ARGS_KEEPGENOME = "keepGenomeInfoForRelatedsOnly";

	private static final Set<QC_METRIC> CUSTOMIZABLE_QC_METRICS = Collections.unmodifiableSet(EnumSet.of(QC_METRIC.CALLRATE));

	/**
	 * @see Qc#Qc(String, String, Map, Logger)
	 */
	public RelationAncestryQc(String dir, String plinkPrefix,
														Map<QC_METRIC, String> markerQCThresholds, Logger log) {
		super(dir, plinkPrefix, markerQCThresholds, log);
	}

	public void run(boolean keepGenomeInfoForRelatedsOnly) {
		long time = new Date().getTime();

		if (!markerQc(RelationAncestryQc.MARKER_QC_DIR)) {
			return;
		}

		new File(dir + RelationAncestryQc.SAMPLE_QC_DIR).mkdirs();
		if (!Files.exists(dir + RelationAncestryQc.SAMPLE_QC_DIR + plink + ".bed")) {
			log.report(ext.getTime() + "]\tRunning --exclude miss_drops.dat");
			CmdLine.runDefaults("plink2 --bfile ../" + RelationAncestryQc.MARKER_QC_DIR + plink
													+ " --exclude ../" + RelationAncestryQc.MARKER_QC_DIR
													+ "miss_drops.dat --make-bed --noweb --out "
													+ plink, dir + RelationAncestryQc.SAMPLE_QC_DIR, log);
		}
		PSF.checkInterrupted();
		if (!Files.exists(dir + RelationAncestryQc.SAMPLE_QC_DIR + "missing.imiss")) {
			log.report(ext.getTime() + "]\tRunning --missing");
			CmdLine.runDefaults("plink2 --bfile " + plink
													+ " --geno 1 --mind 1 --missing --out missing --noweb",
													dir + RelationAncestryQc.SAMPLE_QC_DIR, log);
		}
		PSF.checkInterrupted();

		new File(dir + RelationAncestryQc.LD_PRUNING_DIR).mkdirs();
		if (!Files.exists(dir + RelationAncestryQc.LD_PRUNING_DIR + plink + ".bed")) {
			log.report(ext.getTime()
								 + "]\tRunning --mind 0.05 (removes samples with callrate <95% for the markers that did pass QC)");
			CmdLine.runDefaults("plink2 --bfile ../" + RelationAncestryQc.SAMPLE_QC_DIR + plink
													+ " --mind 0.05 --make-bed --noweb --out " + plink,
													dir + RelationAncestryQc.LD_PRUNING_DIR,
													log);
		}
		PSF.checkInterrupted();
		if (!Files.exists(dir + RelationAncestryQc.LD_PRUNING_DIR + plink + ".prune.in")) {
			log.report(ext.getTime() + "]\tRunning --indep-pairwise 50 5 0.3");
			CmdLine.runDefaults("plink2 --noweb --bfile " + plink + " --indep-pairwise 50 5 0.3 --out "
													+ plink, dir + RelationAncestryQc.LD_PRUNING_DIR, log);
		}
		PSF.checkInterrupted();

		new File(dir + RelationAncestryQc.GENOME_DIR).mkdirs();
		if (!Files.exists(dir + RelationAncestryQc.GENOME_DIR + plink + ".bed")) {
			log.report(ext.getTime() + "]\tRunning --extract " + plink + ".prune.in");
			CmdLine.runDefaults("plink2 --bfile ../" + RelationAncestryQc.LD_PRUNING_DIR + plink
													+ " --extract ../"
													+ RelationAncestryQc.LD_PRUNING_DIR
													+ plink + ".prune.in --make-bed --noweb --out " + plink,
													dir + RelationAncestryQc.GENOME_DIR,
													log);
		}
		PSF.checkInterrupted();
		if (!Files.exists(dir + RelationAncestryQc.GENOME_DIR + plink + ".genome")) {
			log.report(ext.getTime() + "]\tRunning --genome"
								 + (keepGenomeInfoForRelatedsOnly ? " --min 0.1" : ""));
			CmdLine.runDefaults("plink2 --noweb --bfile " + plink + " --genome"
													+ (keepGenomeInfoForRelatedsOnly ? " --min 0.1" : "") + " --out " + plink,
													dir + RelationAncestryQc.GENOME_DIR, log);
		}
		PSF.checkInterrupted();
		if (!keepGenomeInfoForRelatedsOnly
				&& !Files.exists(dir + RelationAncestryQc.GENOME_DIR + "mds20.mds")) {
			log.report(ext.getTime() + "]\tRunning --mds-plot 20");
			CmdLine.runDefaults("plink2 --bfile " + plink + " --read-genome " + plink
													+ ".genome --cluster --mds-plot 20 --out mds20 --noweb",
													dir + RelationAncestryQc.GENOME_DIR,
													log);
		}
		PSF.checkInterrupted();
		if (!Files.exists(dir + RelationAncestryQc.GENOME_DIR + plink + ".genome_keep.dat")) {
			log.report(ext.getTime() + "]\tRunning flagRelateds");
			String lrrFile = dir + RelationAncestryQc.GENOME_DIR + "lrr_sd.xln";
			if (!Files.exists(lrrFile)) {
				lrrFile = dir + "../../lrr_sd.xln";
			}
			Plink.flagRelateds(dir + RelationAncestryQc.GENOME_DIR + plink + ".genome",
												 dir + RelationAncestryQc.GENOME_DIR + plink + ".fam",
												 dir + RelationAncestryQc.MARKER_QC_DIR + "missing.imiss", lrrFile,
												 Plink.FLAGS,
												 Plink.THRESHOLDS, 4,
												 false);
		}
		PSF.checkInterrupted();

		new File(dir + RelationAncestryQc.ANCESTRY_DIR).mkdirs();
		if (!Files.exists(dir + RelationAncestryQc.ANCESTRY_DIR + "unrelateds.txt")) {
			log.report(ext.getTime() + "]\tCopying " + RelationAncestryQc.GENOME_DIR + plink
								 + ".genome_keep.dat to " + RelationAncestryQc.ANCESTRY_DIR + "unrelateds.txt");
			Files.copyFile(dir + RelationAncestryQc.GENOME_DIR + plink + ".genome_keep.dat",
										 dir + RelationAncestryQc.ANCESTRY_DIR + "unrelateds.txt");
		}
		PSF.checkInterrupted();
		if (!Files.exists(dir + RelationAncestryQc.ANCESTRY_DIR + plink + ".bed")) {
			log.report(ext.getTime() + "]\tRunning --extract " + plink
								 + ".prune.in (again, this time to " + RelationAncestryQc.ANCESTRY_DIR + ")");
			CmdLine.runDefaults("plink2 --bfile ../" + RelationAncestryQc.GENOME_DIR + plink
													+ " --make-bed --noweb --out "
													+ plink, dir + RelationAncestryQc.ANCESTRY_DIR, log);
		}
		PSF.checkInterrupted();

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
		fullGamut(dir, plinkPrefix, keepGenomeInfoForRelatedsOnly, log,
							MarkerQC.DEFAULT_METRIC_THRESHOLDS);
	}

	/**
	 * 
	 * @param keepGenomeInfoForRelatedsOnly true to save disk usage if unrelated genome info is not
	 *        required
	 * @param markerQCThresholds TODO
	 * @return full path to plinkroot of QC'd plink dataset
	 */
	public static void fullGamut(String dir, String plinkPrefix,
															 boolean keepGenomeInfoForRelatedsOnly, Logger log,
															 Map<QC_METRIC, String> markerQCThresholds) {
		new RelationAncestryQc(dir, plinkPrefix, markerQCThresholds,
													 log).run(keepGenomeInfoForRelatedsOnly);
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
		CLI c = new CLI(RelationAncestryQc.class);
		c.addArgWithDefault(CLI.ARG_INDIR, "directory with binary plink dataset", "./");
		c.addArgWithDefault(CLI.ARG_PLINKROOT, CLI.DESC_PLINKROOT, Qc.DEFAULT_PLINKROOT);
		c.addArgWithDefault(RelationAncestryQc.ARGS_KEEPGENOME, "if no MDS will be run, smaller file",
												String.valueOf(true));
		for (QC_METRIC metric : CUSTOMIZABLE_QC_METRICS) {
			String defaultThreshold = MarkerQC.DEFAULT_METRIC_THRESHOLDS.get(metric);
			c.addArgWithDefault(metric.getKey(), metric.getCLIDescription(), defaultThreshold);
		}
		c.addArgWithDefault(CLI.ARG_LOG, CLI.DESC_LOG, "fullGamutOfMarkerAndSampleQC.log");

		c.parseWithExit(args);

		String dir = c.get(CLI.ARG_INDIR);
		String inputPlinkroot = c.get(CLI.ARG_PLINKROOT);
		boolean keepGenomeInfoForRelatedsOnly = Boolean.parseBoolean(c.get(RelationAncestryQc.ARGS_KEEPGENOME));
		Map<QC_METRIC, String> markerQCThresholds = Maps.newEnumMap(QC_METRIC.class);
		for (QC_METRIC metric : CUSTOMIZABLE_QC_METRICS) {
			String defaultThreshold = MarkerQC.DEFAULT_METRIC_THRESHOLDS.get(metric);
			c.addArgWithDefault(metric.getKey(), metric.getCLIDescription(), defaultThreshold);
		}
		Logger log = new Logger(dir + c.get(CLI.ARG_LOG));

		try {
			RelationAncestryQc.fullGamut(dir, inputPlinkroot, keepGenomeInfoForRelatedsOnly, log,
																	 markerQCThresholds);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
