package org.genvisis.one.JL.mtDNA;

import java.io.File;

import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCOps;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
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
				VariantContextBuilder builder = new VariantContextBuilder(vc);
//				builder.id(ID)
				writer.add(vc);
			}
		}
		reader.close();
		writer.close();
	}

	private static String filterCR(String inputVCF, double cr, String outDir, String pseqDir, Logger log) {
		String istats = runIstats(inputVCF, pseqDir, outDir, log);
		StringBuilder filter = new StringBuilder();
		String root = outDir + VCFOps.getAppropriateRoot(inputVCF, true) + "cr" + cr;
		String excludes = root + ".excludedSamps.txt";
		filter.append("cat " + istats + " | grep -v NALT|awk '$6<=" + cr + "'>" + excludes + "\n");
		String out1 = root + ".HQ.vcf";
		filter.append("vcftools --vcf " + inputVCF + " --remove excludeSamples.txt --max-missing " + cr
				+ " --recode --recode-INFO-all --out " + out1);
		if (!Files.exists(out1)) {
			String bat = out1 + ".bat";
			Files.write(filter.toString(), bat);
			Files.chmod(bat);
			log.reportTimeInfo(filter.toString());
			CmdLine.run(bat, outDir);
		}
		return out1;
	}

	private static String runIstats(String vcf, String pseqDir, String outDir, Logger log) {
		String out = outDir + ext.rootOf(vcf) + ".istats";
		String bat = out + ".bat";
		String proj = "proj_" + VCFOps.getAppropriateRoot(vcf, true);
		StringBuilder istatGetter = new StringBuilder(pseqDir + "pseq " + proj + " new-project\n");
		istatGetter.append(pseqDir + "pseq " + proj + " load-vcf --vcf " + vcf + "\n");
		istatGetter.append(pseqDir + "pseq " + proj + " i-stats >" + out);

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
//		runIstats(inputVCF, pseqDir, outDir, log);
	}

}
