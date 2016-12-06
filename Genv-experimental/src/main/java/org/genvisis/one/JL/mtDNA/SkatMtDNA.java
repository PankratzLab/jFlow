package org.genvisis.one.JL.mtDNA;

import java.io.File;


import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;
import org.genvisis.seq.manage.VCOps;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * @author Kitty Quick for processing and analyzing mtDNA variants with skat
 * 
 *         https://cran.r-project.org/web/packages/SKAT/SKAT.pdf
 * 
 *         https://github.com/ttimbers/SKAT_NGS-2015/blob/master/
 *         NGS_GWAS_via_SKAT.md
 */
public class SkatMtDNA {

	// TODO, generate fam with pheno /covar
	private static final String MIT_IMPACT = "MitImpact_id";

	private static String filter(String vcf, String anno, String annoFilt, String outDir) {
		VCFFileReader reader = new VCFFileReader(new File(vcf), true);
		String outputVCF = outDir + VCFOps.getAppropriateRoot(vcf, true) + anno + ".vcf";
		VariantContextWriter writer = VCFOps.initWriterWithHeader(reader, outputVCF, VCFOps.DEFUALT_WRITER_OPTIONS,
				new Logger());
		for (VariantContext vc : reader) {
			String annoVc = VCOps.getAnnotationsFor(new String[] { anno }, vc, ".")[0];
			if (!annoVc.equals(annoFilt)) {
				VariantContextBuilder builder = new VariantContextBuilder(vc);
				builder.id(new VCOps.LocusID(vc).getId());
				writer.add(vc);
			}
		}
		reader.close();
		writer.close();
		return outputVCF;
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

	private static void addPhenoToFam(String famFile, VcfPopulation vpop, String cases, String controls) {
		String[][] fam = HashVec.loadFileToStringMatrix(famFile, false, null, false);

		for (String[] famline : fam) {
			if (vpop.getPopulationForInd(famline[0], RETRIEVE_TYPE.SUPER).equals(cases)) {
				famline[famline.length - 1] = "2";

			} else if (vpop.getPopulationForInd(famline[0], RETRIEVE_TYPE.SUPER).equals(cases)) {
				famline[famline.length - 1] = "1";

			} else {
				famline[famline.length - 1] = "-9";

			}
		}
	}

	private static class SNPInfo {
		private String fileName;
		private String setHeader;

		private SNPInfo(String fileName, String setHeader) {
			super();
			this.fileName = fileName;
			this.setHeader = setHeader;
		}

	}

	private static SNPInfo generateSNPInfo(String vcf, String outputDir, String anno, boolean force) {
		String fileName = outputDir + VCFOps.getAppropriateRoot(vcf, true) + "." + anno + ".snpInfo";

		StringBuilder snpInfo = new StringBuilder();
		snpInfo.append("Name\t" + anno);
		VCFFileReader reader = new VCFFileReader(new File(vcf), true);
		for (VariantContext vc : reader) {
			if (force) {
				snpInfo.append(vc.getID() + "\t" + anno);
			} else {
				String ofInterest = VCOps.getAnnotationsFor(new String[] { anno }, vc, ".")[0];
				if (!ofInterest.equals(".")) {
					snpInfo.append(vc.getID() + "\t" + ofInterest);
				}
			}
		}
		return new SNPInfo(fileName, anno);

	}

	public static void main(String[] args) {
		String pseqDir = "/Users/Kitty/bin/plinkseq-0.10/";
		String outDir = "/Volumes/Beta/data/mtDNA-dev/vcf/analysis/skat/";
		new File(outDir).mkdirs();
		String inputVCF = "/Volumes/Beta/data/mtDNA-dev/vcf/analysis/ARIC_CUSHING_EPP_OSTEO_FP_MITO.chrM.rcrs.poly.disease.conv.hg19_multianno.eff.gatk.sed1000g.vcf";
		String vpopFile = "";
		String cases = "Cushing";
		String controls = "Aric";
		String covarFile = "";

		Logger log = new Logger(outDir + "log.log");
		VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.ANY, log);
		String hqVcf = filterCR(inputVCF, 0.95, outDir, pseqDir, log);
		String nsVcf = filter(hqVcf, MIT_IMPACT, ".", outDir);
		SNPInfo nsSNPInfo = generateSNPInfo(nsVcf, outDir, "mtDNAFull", true);

		// runIstats(inputVCF, pseqDir, outDir, log);
	}

}
