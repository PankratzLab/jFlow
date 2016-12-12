package org.genvisis.one.JL.mtDNA;

import java.io.File;
import java.util.ArrayList;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.analysis.mtdna.HaplogroupSelector;
import org.genvisis.seq.manage.GenotypeOps;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCOps;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

public class MtDNABurden {

	private enum FILTER_TYPE {
		NOT_EQUAL, CONTAINS, LESS_THAN;
	}

	private enum VARIANT_TYPES {
		APOGEE("APOGEE", new String[] { "P" }, FILTER_TYPE.CONTAINS), LHON_ASSOC("Associated_disease",
				new String[] { "Lebers_optic_atrophy", "Optic_neuropathy", "LHON" },
				FILTER_TYPE.CONTAINS), DISEASE_ASSOC("Associated_disease", new String[] { "-", "." },
						FILTER_TYPE.NOT_EQUAL), MIT_IMPACT("MitImpact_id", new String[] { "." },
								FILTER_TYPE.NOT_EQUAL), GB("mtPoly.AF", null, FILTER_TYPE.LESS_THAN);

		private String flag;
		private String[] filt;
		private FILTER_TYPE fType;

		private VARIANT_TYPES(String flag, String[] filt, FILTER_TYPE fType) {
			this.flag = flag;
			this.filt = filt;
			this.fType = fType;
		}
	}

	private static String filter(String vcf, VARIANT_TYPES vTypes, String outDir, double genBankMaf, Logger log) {
		String outputVCF = outDir + VCFOps.getAppropriateRoot(vcf, true) + vTypes.toString()
				+ (vTypes == VARIANT_TYPES.GB ? ".gb." + genBankMaf : "") + ".vcf";
		String finalOut = ext.addToRoot(outputVCF, ".sed");

		if (!Files.exists(finalOut)) {
			VCFFileReader reader = new VCFFileReader(new File(vcf), false);

			VariantContextWriter writer = VCFOps.initWriterWithHeader(reader, outputVCF, VCFOps.DEFUALT_WRITER_OPTIONS,
					new Logger());

			int keep = 0;
			for (VariantContext vc : reader) {
				String annoVc = VCOps.getAnnotationsFor(new String[] { vTypes.flag }, vc, ".")[0];
				if (!vc.isMonomorphicInSamples()) {
					boolean use = false;
					switch (vTypes.fType) {
					case CONTAINS:
						use = false;
						for (String filt : vTypes.filt) {
							if (annoVc.contains(filt)) {
								use = true;
							}
						}
						break;
					case LESS_THAN:
						use = true;
						if (vTypes == VARIANT_TYPES.GB) {
							try {
								double gf = Double.parseDouble(annoVc) / 100;
								use = gf < genBankMaf;
							} catch (NumberFormatException nfe) {

							}
						} else {
							writer.close();
							throw new IllegalArgumentException("Must use GB");
						}
						break;
					case NOT_EQUAL:
						use = true;
						for (String filt : vTypes.filt) {
							if (annoVc.equals(filt)) {
								use = false;
							}
						}
						break;
					default:
						break;

					}

					if (use) {
						keep++;
						VariantContextBuilder builder = new VariantContextBuilder(vc);
						builder.id(new VCOps.LocusID(vc).getId().replaceAll(":", "_"));
						writer.add(builder.make());

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
		String proj = outDir + "proj_" + VCFOps.getAppropriateRoot(vcf, true);
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

		if (!Files.exists(out1 + ".recode.vcf")) {
			String bat = out1 + ".bat";
			Files.write(filter.toString(), bat);
			Files.chmod(bat);
			log.reportTimeInfo(filter.toString());
			CmdLine.runCommandWithFileChecks(new String[] { bat }, "", null, null, true, true, false, log);
		}
		return new FilterResults(out1 + ".recode.vcf", excludes);
	}

	private static String makeHetVcf(String inputVcf, String outDir, Logger log) {
		VCFFileReader reader = new VCFFileReader(new File(inputVcf), false);
		String root = outDir + VCFOps.getAppropriateRoot(inputVcf, true) + "_HETONLY";
		String filtVcf = root + ".vcf";
		VariantContextWriter writer = VCFOps.initWriterWithHeader(reader, filtVcf, VCFOps.DEFUALT_WRITER_OPTIONS,
				new Logger());

		for (VariantContext vc : reader) {

			VariantContextBuilder builder = new VariantContextBuilder(vc);
			ArrayList<Genotype> genotypes = new ArrayList<>();
			for (Genotype g : vc.getGenotypes()) {
				GenotypeBuilder gb = new GenotypeBuilder(g);
				if (!g.isHet()) {
					gb.alleles(GenotypeOps.getNoCall());
				}
				genotypes.add(gb.make());
			}

			GenotypesContext bc = GenotypesContext.create(genotypes);
			builder.genotypes(bc);
			writer.add(builder.make());
		}
		reader.close();
		writer.close();
		String finalOut = root + ".sed.vcf";
		StringBuilder sed = new StringBuilder();
		sed.append("sed 's/##fileformat=VCFv4.2/##fileformat=VCFv4.1/g' " + filtVcf + " > " + finalOut);
		String bat = finalOut + ".bat";
		Files.write(sed.toString(), bat);
		Files.chmod(bat);
		CmdLine.runCommandWithFileChecks(new String[] { bat }, "", null, null, true, true, false, log);
		return finalOut;
	}

	private static void runPseq(String pseqDir, String outDir, String vcf, String mtGeneFile, String phe, double maf,
			Logger log) {

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
					+ " maf=0-" + maf + " --tests burden fw --perm 50000 >" + results + " \n");

		}
		String singleVar = root + ".maf" + maf + ".SingleVarResults";
		if (!Files.exists(singleVar)) {
			builder.append(pseqDir + pse + root + "  v-assoc --phenotype phe1 " + "--mask maf=0-" + maf + " >"
					+ singleVar + " \n");
		}

		builder.append("head -n1 " + results + " >" + outDir + "results.txt\n");
		builder.append("grep -v \"=\\|LOCUS\" " + outDir + "*.results >>" + outDir + "results.txt\n");
		String bat = root + ".bat";
		Files.write(builder.toString(), bat);
		Files.chmod(bat);
		log.reportTimeInfo(builder.toString());

		CmdLine.runCommandWithFileChecks(new String[] { bat }, "", null, null, true, true, false, log);
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

		for (int i = 4; i < 20; i++) {
			String phe = HaplogroupSelector.run(haps, "CUSHING_FREQ_V2", controls, vpopFile, filterResults.excludeFile,
					outDir, i, 1);
			double[] mafs = new double[] { 0.001, 0.01 };
			String analysisVCF = runKeeper(hqVcf, ext.rootOf(phe, false) + HaplogroupSelector.KEEP_EXT, outDir, log);
			runIstats(analysisVCF, pseqDir, outDir, log);

			String mtDNADefs = "/Volumes/Beta/data/mtDNA-dev/vcf/analysis/mtUniport.reg";

			for (VARIANT_TYPES vType : VARIANT_TYPES.values()) {
				if (vType != VARIANT_TYPES.GB) {
					for (double maf : mafs) {
						String vcf = filter(analysisVCF, vType, outDir, maf, log);
						String hetVcf = makeHetVcf(vcf, outDir, log);
						runPseq(pseqDir, outDir, hetVcf, mtDNADefs, phe, 1, log);
						runPseq(pseqDir, outDir, vcf, mtDNADefs, phe, maf, log);
						String gbVCF = filter(vcf, VARIANT_TYPES.GB, outDir, maf, log);
						runPseq(pseqDir, outDir, gbVCF, mtDNADefs, phe, 1, log);
					}
				}
			}
		}
	}
}
