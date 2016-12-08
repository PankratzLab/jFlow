package org.genvisis.one.JL.mtDNA;

import java.io.File;
import java.util.ArrayList;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.analysis.mtdna.HaplogroupSelector;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCOps;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

public class MtDNABurden {

	private enum VARIANT_TYPES {
		LHON_ASSOC("Associated_disease", new String[] { "Lebers_optic_atrophy", "Optic_neuropathy", "LHON" },
				true), DISEASE_ASSOC("Associated_disease", new String[] { "-", "." },
						false), MIT_IMPACT("MitImpact_id", new String[] { "." }, false);

		private String flag;
		private String[] filt;
		private boolean contains;

		private VARIANT_TYPES(String flag, String[] filt, boolean contains) {
			this.flag = flag;
			this.filt = filt;
			this.contains = contains;
		}
	}

	private enum VARIANT_SETS {
		UNIPROT_NAME("Uniprot_name", false), MTDNA_FULL("mtDNAFull", true);
		private String flag;
		private boolean force;

		private VARIANT_SETS(String flag, boolean force) {
			this.flag = flag;
			this.force = force;
		}
	}

	private static String filter(String vcf, VARIANT_TYPES vTypes, String outDir, double genBankMaf, Logger log) {
		String outputVCF = outDir + VCFOps.getAppropriateRoot(vcf, true) + vTypes.toString() + ".gb." + genBankMaf
				+ ".vcf";
		String finalOut = ext.addToRoot(outputVCF, ".sed");

		if (!Files.exists(finalOut)) {
			VCFFileReader reader = new VCFFileReader(new File(vcf), false);

			VariantContextWriter writer = VCFOps.initWriterWithHeader(reader, outputVCF, VCFOps.DEFUALT_WRITER_OPTIONS,
					new Logger());

			int keep = 0;
			for (VariantContext vc : reader) {
				String annoVc = VCOps.getAnnotationsFor(new String[] { vTypes.flag }, vc, ".")[0];
				boolean use = vTypes.contains ? false : true;
				if (!vc.isMonomorphicInSamples()) {
					String gb = VCOps.getAnnotationsFor(new String[] { "mtPoly.AF" }, vc, ".")[0];
					boolean gbPass = true;
					try {
						double gf = Double.parseDouble(gb) / 100;
						gbPass = gf < genBankMaf;
					} catch (NumberFormatException nfe) {

					}
					if (gbPass) {
						for (String filt : vTypes.filt) {
							if (!vTypes.contains) {
								if (annoVc.equals(filt)) {
									use = false;
								}
							} else {
								if (annoVc.contains(filt)) {
									use = true;
								}

							}
						}
						if (use) {
							keep++;
							VariantContextBuilder builder = new VariantContextBuilder(vc);
							builder.id(new VCOps.LocusID(vc).getId().replaceAll(":", "_"));
							writer.add(builder.make());
						}
					}
				}
			}
			log.reportTimeInfo("Kept " + keep + " variants for " + vTypes.toString());
			reader.close();
			writer.close();
			StringBuilder sed = new StringBuilder();
			sed.append("sed 's/##fileformat=VCFv4.2/##fileformat=VCFv4.1/g' " + outputVCF + " > " + finalOut);
			String bat = outputVCF + ".bat";
			Files.write(sed.toString(), bat);
			Files.chmod(bat);
			CmdLine.runCommandWithFileChecks(new String[] { bat }, "", null, null, true, true, false, log);
		}
		return finalOut;
	}

	private static String filterVpop(String inputVCF, String outDir, VcfPopulation vpop, String casePop,
			String[] controlPop, Logger log) {
		StringBuilder filter = new StringBuilder();
		String root = outDir + VCFOps.getAppropriateRoot(inputVCF, true) + "vp" + ext.rootOf(vpop.getFileName()) + "_"
				+ casePop + "_" + Array.toStr(controlPop, "_");
		String out = root;

		if (!Files.exists(out + ".recode.vcf")) {
			String excludes = root + ".excludedSamps.txt";
			ArrayList<String> excluded = new ArrayList<>();
			String[] samps = VCFOps.getSamplesInFile(inputVCF);
			for (String sample : samps) {
				String[] pop = vpop.getPopulationForInd(sample, RETRIEVE_TYPE.SUPER);
				if (pop.length == 0 || (!pop[0].equals(casePop) && ext.indexOfStr(pop[0], controlPop) < 0)) {
					excluded.add(sample);
				}
			}

			Files.writeIterable(excluded, excludes);
			filter.append("vcftools --vcf " + inputVCF + " --remove " + excludes + " --recode --recode-INFO-all --out "
					+ out);
			String bat = out + ".bat";
			log.reportTimeInfo("Removing " + excluded.size() + " samples for non-case/control status");
			Files.write(filter.toString(), bat);
			Files.chmod(bat);
			log.reportTimeInfo(filter.toString());
			CmdLine.runCommandWithFileChecks(new String[] { bat }, "", null, null, true, true, false, log);
		}
		return out + ".recode.vcf";

	}

	private static String runKeeper(String inputVcf, String keepFile, String outputDir, Logger log) {
		String out = outputDir + ext.rootOf(inputVcf) + "_" + ext.rootOf(keepFile);
		StringBuilder filter = new StringBuilder(
				"vcftools --vcf " + inputVcf + " --keep " + keepFile + " --recode --recode-INFO-all --out " + out);
		if (!Files.exists(out + ".recode.vcf")) {
			String bat = out + ".bat";
			Files.write(filter.toString(), bat);
			Files.chmod(bat);
			log.reportTimeInfo(filter.toString());
			CmdLine.runCommandWithFileChecks(new String[] { bat }, "", null, null, true, true, false, log);
		}
		return out + ".recode.vcf";

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
			CmdLine.runCommandWithFileChecks(new String[] { bat }, "", null, null, true, true, false, log);
		}
		return out;
	}

	private static class FilterResults {
		private String vcf;
		private String excludeFile;

		public FilterResults(String vcf, String excludeFile) {
			super();
			this.vcf = vcf;
			this.excludeFile = excludeFile;
		}

	}

	private static FilterResults filterCR(String inputVCF, double cr, String outDir, String pseqDir, Logger log) {
		String istats = runIstats(inputVCF, pseqDir, outDir, log);
		StringBuilder filter = new StringBuilder();
		String root = outDir + VCFOps.getAppropriateRoot(inputVCF, true) + "cr" + cr;
		String excludes = root + ".excludedSamps.txt";
		filter.append("cat " + istats + " | grep -v NALT|awk '$6<=" + cr + "'>" + excludes + "\n");
		String out1 = root + ".HQ";
		filter.append("vcftools --vcf " + inputVCF + " --remove " + excludes + " --max-missing " + cr
				+ " --recode --recode-INFO-all --out " + out1 + "\n");

		// filter.append("sed -i 's/##fileformat=VCFv4.2/##fileformat=VCFv4.1/g'
		// " + out1 + ".recode.vcf");
		if (!Files.exists(out1 + ".recode.vcf")) {
			String bat = out1 + ".bat";
			Files.write(filter.toString(), bat);
			Files.chmod(bat);
			log.reportTimeInfo(filter.toString());
			CmdLine.runCommandWithFileChecks(new String[] { bat }, "", null, null, true, true, false, log);
		}
		return new FilterResults(out1 + ".recode.vcf", excludes);
	}

	private static String filterMaf(String inputVCF, double maf, String outDir, Logger log) {
		String root = outDir + VCFOps.getAppropriateRoot(inputVCF, true) + "maf" + maf;
		StringBuilder filter = new StringBuilder();
		filter.append("vcftools --vcf " + inputVCF + " --max-maf  " + maf + " --recode --recode-INFO-all --out " + root
				+ "\n");
		if (!Files.exists(root + ".recode.vcf")) {
			String bat = root + ".bat";
			Files.write(filter.toString(), bat);
			Files.chmod(bat);
			log.reportTimeInfo(filter.toString());
			CmdLine.run(bat, outDir);
		}
		return root + ".recode.vcf";
	}

	private static void runPseq(String pseqDir, String outDir, String vcfIn, VARIANT_TYPES vTypes, String mtGeneFile,
			String phe, double maf, Logger log) {
		String vcf = filter(vcfIn, vTypes, outDir, maf, log);

		StringBuilder builder = new StringBuilder();
		String root = outDir + "proj_" + VCFOps.getAppropriateRoot(vcf, true);
		String pse = "pseq ";
		builder.append(pseqDir + pse + root + " new-project\n");
		builder.append(pseqDir + pse + root + " load-vcf --vcf " + vcf + "\n");
		builder.append(pseqDir + pse + root + " load-pheno --file " + phe + "\n");
		builder.append(
				pseqDir + pse + root + " loc-load --file " + mtGeneFile + " --group " + ext.rootOf(mtGeneFile) + "\n");
		String istatsFile = root + ".istats";
		if (!Files.exists(istatsFile)) {

			builder.append(pseqDir + pse + root + " i-stats >" + istatsFile + "\n");
		}
		String results = root + ".maf" + maf + ".results";
		if (!Files.exists(results)) {
			builder.append(pseqDir + pse + root + "  assoc --phenotype phe1 --mask loc.group=" + ext.rootOf(mtGeneFile)
					+ " maf=0-0.05 --tests burden fw --perm 50000 >" + results + " \n");

		}
		String bat = root + ".bat";
		Files.write(builder.toString(), bat);
		Files.chmod(bat);
		log.reportTimeInfo(builder.toString());
		CmdLine.runCommandWithFileChecks(new String[] { bat }, "", null, null, true, true, false, log);

		// #run T1 and T5 tests
		// # "$plinkseqDir"pseq projNS assoc --phenotype phe1 --mask
		// loc.group=mtDNA maf=0-0.05 --tests burden fw --perm 50000
		// >results.T5.mtdna.NS.txt
		// # "$plinkseqDir"pseq projNS assoc --phenotype phe1 --mask
		// loc.group=mtDNA maf=0-0.01 --tests burden fw --perm 50000
		// >results.T1.mtdna.NS.txt
		//

	}

	public static void main(String[] args) {
		String pseqDir = "/Users/Kitty/bin/plinkseq-0.10/";
		String outDir = "/Volumes/Beta/data/mtDNA-dev/vcf/analysis/burden/";
		new File(outDir).mkdirs();
		String inputVCF = "/Volumes/Beta/data/mtDNA-dev/vcf/analysis/mtvar.final.vcf";
		String vpopFile = "/Volumes/Beta/data/mtDNA-dev/vcf/analysis/skat.vpop";
		String cases = "CUSHING_FREQ_V2";
		String[] controls = new String[] { "ARIC", "CUSHINGS" };
		String haps = "/Volumes/Beta/data/mtDNA-dev/vcf/analysis/ARIC_CUSHING_EPP_OSTEO_FP_MITO.chrM.rcrs.poly.disease.conv.hg19_multianno.eff.gatk.sed1000g.haplotypes";

		Logger log = new Logger(outDir + "log.log");

		VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.ANY, log);
		vpop.report();

		String caseControlVcf = filterVpop(inputVCF, outDir, vpop, cases, controls, log);
		FilterResults filterResults = filterCR(caseControlVcf, 0.95, outDir, pseqDir, log);
		String hqVcf = filterResults.vcf;
		String phe = HaplogroupSelector.run(haps, "CUSHING_FREQ_V2", controls, vpopFile, filterResults.excludeFile,
				outDir, 5, 1);
		// System.exit(1);

		double[] mafs = new double[] { 0.01 };
		String analysisVCF = runKeeper(hqVcf, ext.rootOf(phe, false) + HaplogroupSelector.KEEP_EXT, outDir, log);
		String mtDNADefs = "/Volumes/Beta/data/mtDNA-dev/vcf/analysis/mtUniport.reg";

		for (VARIANT_TYPES vType : VARIANT_TYPES.values()) {
			for (int i = 0; i < mafs.length; i++) {
				runPseq(pseqDir, outDir, analysisVCF, vType, mtDNADefs, phe, mafs[i], log);
			}
		}
	}

}
