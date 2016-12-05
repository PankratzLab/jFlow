package org.genvisis.one.JL.mtDNA;

import java.io.File;

import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCOps;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * @author Kitty Quick for processing and analyzing mtDNA variants with skat
 */
public class SkatMtDNA {
	private static final String MIT_IMPACT = "MitImpact_id";

	private static void filter(String vcf, String anno, String annoFilt, String outDir) {
		VCFFileReader reader = new VCFFileReader(new File(vcf), true);
		String outputVCF = outDir + VCFOps.getAppropriateRoot(vcf, true) + anno + ".vcf";
		VariantContextWriter writer = VCFOps.initWriterWithHeader(reader, outputVCF, VCFOps.DEFUALT_WRITER_OPTIONS,
				new Logger());
		for (VariantContext vc : reader) {
			String annoVc = VCOps.getAnnotationsFor(new String[] { anno }, vc, ".")[0];
			if (!annoVc.equals(annoFilt)) {
				writer.add(vc);
			}
		}

		reader.close();
		writer.close();
	}

	private static void filterCR(String inputVCF, double cr, String outDir, String pseqDir, Logger log) {

		// cat i-statsats.mtdna.ALL.txt | grep -v NALT|awk '$6<=.95'>
		// excludeSamples.txt
		// vcftools --vcf
		// ARIC_CUSHING_EPP_OSTEO_FP_MITO.chrM.rcrs.poly.disease.conv.hg19_multianno.eff.gatk.sed1000g.vcf
		// --remove excludeSamples.txt --recode --recode-INFO-all --out
		// ARIC_CUSHING_EPP_OSTEO_FP_MITO.chrM.rcrs.poly.disease.conv.hg19_multianno.eff.gatk.sed1000g.HQ.SAMPS
		// vcftools --vcf
		// ARIC_CUSHING_EPP_OSTEO_FP_MITO.chrM.rcrs.poly.disease.conv.hg19_multianno.eff.gatk.sed1000g.HQ.SAMPS.recode.vcf
		// --max-missing $callRate --recode --recode-INFO-all --out
		// ARIC_CUSHING_EPP_OSTEO_FP_MITO.chrM.rcrs.poly.disease.conv.hg19_multianno.eff.gatk.sed1000g.HQ.SAMPS.$callRate.CR

	}

	private static String runIstats(String vcf, String pseqDir, String outDir, Logger log) {
		String out = outDir + ext.rootOf(vcf) + ".istats";
		String bat = out + ".bat";
		String proj = "proj_" + VCFOps.getAppropriateRoot(vcf, true);
		StringBuilder istatGetter = new StringBuilder(pseqDir + "pseq " + proj + " new-project\n");
		istatGetter.append(pseqDir + "pseq " + proj + " load-vcf --vcf " + vcf + "\n");
		istatGetter.append(pseqDir + "pseq " + proj + " i-stats >" + out);

		// "$plinkseqDir"pseq proj new-project
		// #load vcf
		// "$plinkseqDir"pseq proj load-vcf --vcf
		// ARIC_CUSHING_EPP_OSTEO_FP_MITO.chrM.rcrs.poly.disease.conv.hg19_multianno.eff.gatk.sed1000g.vcf
		// "$plinkseqDir"pseq proj i-stats >i-statsats.mtdna.ALL.txt
		if (!Files.exists(out)) {
			Files.write(istatGetter.toString(), bat);
			Files.chmod(bat);
			log.reportTimeInfo(istatGetter.toString());
			CmdLine.run(bat, outDir);
		}
		return out;
	}

	public static void main(String[] args) {
		String pseqDir = "/Users/Kitty/bin/plinkseq-0.10/";
		String outDir = "/Volumes/Beta/data/mtDNA-dev/vcf/analysis/skat/";
		new File(outDir).mkdirs();
		String inputVCF = "/Volumes/Beta/data/mtDNA-dev/vcf/analysis/ARIC_CUSHING_EPP_OSTEO_FP_MITO.chrM.rcrs.poly.disease.conv.hg19_multianno.eff.gatk.sed1000g.vcf";
		Logger log = new Logger(outDir + "log.log");
		runIstats(inputVCF, pseqDir, outDir, log);
	}

}
