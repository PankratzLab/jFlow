package org.genvisis.gwas;

import java.io.File;
import java.util.Map;

import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.gwas.MarkerQC.QC_METRIC;

public abstract class Qc {

	protected static final String DEFAULT_PLINKROOT = "plink";
	public static final String QC_DIR = "quality_control/";

	protected final String dir;
	protected final String plink;
	protected final Logger log;



	/**
	 * @param dir Directory with plink files to run from
	 * @param plinkPrefix prefix of plink binaries
	 * @param log
	 */
	protected Qc(String dir, String plinkPrefix, Logger log) {
		super();
		dir = ext.verifyDirFormat(dir) + QC_DIR;
		if (!dir.startsWith("/") && !dir.contains(":")) {
			dir = ext.verifyDirFormat((new File("./" + dir)).getAbsolutePath());
		}
		this.dir = dir;
		this.plink = plinkPrefix == null ? DEFAULT_PLINKROOT : plinkPrefix;
		this.log = log;
	}



	protected boolean markerQc(String subDir, final Map<QC_METRIC, String> markerQcThresholds) {
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
			MarkerQC.generateCRF(dir + subDir, dir + subDir + "miss.crf", markerQcThresholds);
			int runCode = MarkerQC.parseParameters(dir + subDir + "miss.crf", log, false);
			if (runCode != 0) {
				log.reportError("Failed to perform marker QC with " + dir + subDir + "miss.crf");
				return false;
			}
		}
		PSF.checkInterrupted();
		return true;
	}

	/**
	 * @deprecated Use {@link RelationAncestryQc#fromParameters(String,Logger)} instead
	 */
	@Deprecated
	public static void fromParameters(String filename, Logger log) {
		RelationAncestryQc.fromParameters(filename, log);
	}

	/**
	 * @deprecated Use {@link RelationAncestryQc#main(String[])} instead
	 */
	@Deprecated
	public static void main(String[] args) {
		RelationAncestryQc.main(args);
	}
}
