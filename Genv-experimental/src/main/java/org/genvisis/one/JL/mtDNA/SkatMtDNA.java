package org.genvisis.one.JL.mtDNA;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
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

	private static final String MIT_IMPACT = "MitImpact_id";
	private static final String DISEASE_ASSOC = "Associated_disease";
	private static final String UNIPROT_NAME = "Uniprot_name";

	private enum VARIANT_TYPES {
		MIT_IMPACT("MitImpact_id", "."), DISEASE_ASSOC("Associated_disease", "-,.");

		private String flag;
		private String filt;

		private VARIANT_TYPES(String flag, String filt) {
			this.flag = flag;
			this.filt = filt;
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

	private static String filter(String vcf, String anno, String annoFilt, String outDir, Logger log) {
		String outputVCF = outDir + VCFOps.getAppropriateRoot(vcf, true) + anno + ".vcf";
		String finalOut = ext.addToRoot(outputVCF, ".sed");
		String[] filterAnnos = annoFilt.split(",");
		if (!Files.exists(finalOut)) {
			VCFFileReader reader = new VCFFileReader(new File(vcf), false);

			VariantContextWriter writer = VCFOps.initWriterWithHeader(reader, outputVCF, VCFOps.DEFUALT_WRITER_OPTIONS,
					new Logger());

			int keep = 0;
			for (VariantContext vc : reader) {
				String annoVc = VCOps.getAnnotationsFor(new String[] { anno }, vc, ".")[0];
				boolean use = true;
				for (String filt : filterAnnos) {
					if (annoVc.equals(filt)) {
						use = false;
					}
				}
				if (use) {
					keep++;
					VariantContextBuilder builder = new VariantContextBuilder(vc);
					builder.id(new VCOps.LocusID(vc).getId().replaceAll(":", "_"));
					writer.add(builder.make());
				}
			}
			log.reportTimeInfo("Kept " + keep + " variants for " + anno);
			reader.close();
			writer.close();
			StringBuilder sed = new StringBuilder();
			sed.append("sed 's/##fileformat=VCFv4.2/##fileformat=VCFv4.1/g' " + outputVCF + " > " + finalOut);
			String bat = outputVCF + ".bat";
			Files.write(sed.toString(), bat);
			Files.chmod(bat);
			CmdLine.run(bat, outDir);
		}
		return finalOut;
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

	private static String filterCR(String inputVCF, double cr, String outDir, String pseqDir, Logger log) {
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
			CmdLine.run(bat, outDir);
		}
		return out1 + ".recode.vcf";
	}

	private static String filterVpop(String inputVCF, String outDir, VcfPopulation vpop, String casePop,
			String controlPop, Logger log) {
		StringBuilder filter = new StringBuilder();
		String root = outDir + VCFOps.getAppropriateRoot(inputVCF, true) + "vp" + ext.rootOf(vpop.getFileName()) + "_"
				+ casePop + "_" + controlPop;
		String out = root;

		// if (!Files.exists(out)) {

		if (!Files.exists(out + ".recode.vcf")) {
			String excludes = root + ".excludedSamps.txt";
			ArrayList<String> excluded = new ArrayList<>();
			String[] samps = VCFOps.getSamplesInFile(inputVCF);
			for (String sample : samps) {
				String[] pop = vpop.getPopulationForInd(sample, RETRIEVE_TYPE.SUPER);
				if (pop.length == 0 || (!pop[0].equals(casePop) && !pop[0].equals(controlPop))) {
					excluded.add(sample);
				}
			}

			// }
			Files.writeIterable(excluded, excludes);
			filter.append("vcftools --vcf " + inputVCF + " --remove " + excludes + " --recode --recode-INFO-all --out "
					+ out);
			String bat = out + ".bat";
			log.reportTimeInfo("Removing " + excluded.size() + " samples for non-case/control status");
			Files.write(filter.toString(), bat);
			Files.chmod(bat);
			log.reportTimeInfo(filter.toString());
			CmdLine.run(bat, outDir);
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
			CmdLine.run(bat, outDir);
		}
		return out;
	}

	private static void addPhenoToFam(String famFile, VcfPopulation vpop, String cases, String controls) {
		String[][] fam = HashVec.loadFileToStringMatrix(famFile, false, null, false);
		for (int i = 0; i < fam.length; i++) {
			String[] famline = fam[i];
			String[] sPop = vpop.getPopulationForInd(famline[0], RETRIEVE_TYPE.SUPER);
			if (sPop.length > 0 && sPop[0].equals(cases)) {
				famline[famline.length - 1] = "2";
			} else if (sPop.length > 0l && sPop[0].equals(controls)) {
				famline[famline.length - 1] = "1";

			} else {
				famline[famline.length - 1] = "-9";

			}
			fam[i] = famline;
		}
		Files.writeMatrix(fam, famFile, "\t");
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

		if (!Files.exists(fileName)) {
			StringBuilder snpInfo = new StringBuilder();
			// snpInfo.append("Name\t" + anno);
			VCFFileReader reader = new VCFFileReader(new File(vcf), false);
			boolean first = true;
			for (VariantContext vc : reader) {
				if (force) {
					snpInfo.append((first ? "" : "\n") + anno + "\t" + vc.getID());
				} else {
					String ofInterest = VCOps.getAnnotationsFor(new String[] { anno }, vc, ".")[0];
					if (!ofInterest.equals(".")) {
						snpInfo.append((first ? "" : "\n") + ofInterest + "\t" + vc.getID());
					}
				}
				first = false;
			}
			reader.close();
			Files.write(snpInfo.toString(), fileName);
		}
		return new SNPInfo(fileName, anno);
	}

	private static String[] convertToPlink(String inputVCF, String outDir, Logger log) {
		String rootOut = outDir + VCFOps.getAppropriateRoot(inputVCF, true);
		String[] outFiles = PSF.Plink.getPlinkBedBimFam(rootOut);
		if (!Files.exists("", outFiles)) {
			String[] plinkCommand = new String[] { "plink2", "--double-id", "--vcf", inputVCF, "--no-fid",
					"--no-parents", "--no-sex", "--no-pheno", "--out", rootOut };
			if (CmdLine.runCommandWithFileChecks(plinkCommand, "", new String[] { inputVCF }, outFiles, true, true,
					false, log)) {

			}
		}
		return outFiles;
	}

	private static void runSkat(String[] plinks, SNPInfo snInfo, String hapFile, Logger log) {
		StringBuilder skatBuilder = new StringBuilder();
		String root = ext.rootOf(snInfo.fileName, false);
		String ssd = root + ".SSD";
		String info = root + ".info";
		skatBuilder.append("library(SKAT)\n");
		skatBuilder.append("Generate_SSD_SetID('" + plinks[0] + "','" + plinks[1] + "','" + plinks[2] + "','"
				+ snInfo.fileName + "','" + ssd + "','" + info + "')\n");

		skatBuilder.append("fam =Read_Plink_FAM_Cov('" + plinks[2] + "','" + hapFile + "', Is.binary=TRUE)\n");
		skatBuilder.append("SSD.info <- Open_SSD('" + ssd + "', '" + info + "')\n");
		skatBuilder.append(
				"Null_Model <- SKAT_Null_Model(formula = fam$Phenotype ~  fam$HAPH+fam$HAPL, out_type=\"D\")\n");
		skatBuilder.append("All_SKAT_Data  <- SKATBinary.SSD.All(SSD.INFO = SSD.info, obj = Null_Model) \n");
		skatBuilder.append(
				"All_SKAT_Data$results$P.value.bonf =p.adjust(All_SKAT_Data$results$P.value, method = \"bonferroni\", n = length(All_SKAT_Data$results$P.value))\n");
		skatBuilder.append("write.table(x = All_SKAT_Data$results, file = \"" + root
				+ "covar.pvalues\", row.names = FALSE, col.names = TRUE, quote = FALSE, append = FALSE,sep=\"\\t\")\n");
		
		skatBuilder.append(
				"Null_Model <- SKAT_Null_Model(formula = fam$Phenotype ~ 1, out_type=\"D\")\n");
		skatBuilder.append("All_SKAT_Data  <- SKATBinary.SSD.All(SSD.INFO = SSD.info, obj = Null_Model) \n");
		skatBuilder.append(
				"All_SKAT_Data$results$P.value.bonf =p.adjust(All_SKAT_Data$results$P.value, method = \"bonferroni\", n = length(All_SKAT_Data$results$P.value))\n");
		skatBuilder.append("write.table(x = All_SKAT_Data$results, file = \"" + root
				+ "No.covar.pvalues\", row.names = FALSE, col.names = TRUE, quote = FALSE, append = FALSE,sep=\"\\t\")");
		
		String script = root + ".rscript";
		Files.write(skatBuilder.toString(), script);
		CmdLine.run("Rscript " + script, ext.parseDirectoryOfFile(root));

	}

	private static String parseHaplogroupFile(String haplogroups, String outputDir, Logger log) {
		String out = outputDir + ext.rootOf(haplogroups) + ".parsed.txt";
		int hapIndex = ext.indexOfStr("Haplogroup", Files.getHeaderOfFile(haplogroups, log));

		try {
			BufferedReader reader = Files.getAppropriateReader(haplogroups);
			StringBuilder haps = new StringBuilder();
			haps.append("FID\tIID\tHAPH\tHAPL");

			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("\t");
				String haplogroup = line[hapIndex].substring(0, 1);
				if (haplogroup.equals("H")) {
					haps.append("\n" + line[0] + "\t" + line[0] + "\t1\t0");
				} else if (haplogroup.equals("L")) {
					haps.append("\n" + line[0] + "\t" + line[0] + "\t0\t1");
				} else {
					haps.append("\n" + line[0] + "\t" + line[0] + "\t0\t0");
				}
			}
			reader.close();
			Files.write(haps.toString(), out);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return out;

	}

	public static void main(String[] args) {
		String pseqDir = "/Users/Kitty/bin/plinkseq-0.10/";
		String outDir = "/Volumes/Beta/data/mtDNA-dev/vcf/analysis/skat/";
		new File(outDir).mkdirs();
		String inputVCF = "/Volumes/Beta/data/mtDNA-dev/vcf/analysis/mtvar.final.vcf";
		String vpopFile = "/Volumes/Beta/data/mtDNA-dev/vcf/analysis/skat.vpop";
		String cases = "CUSHING_FREQ_V2";
		String controls = "ARIC";
		String haps = "/Volumes/Beta/data/mtDNA-dev/vcf/analysis/ARIC_CUSHING_EPP_OSTEO_FP_MITO.chrM.rcrs.poly.disease.conv.hg19_multianno.eff.gatk.sed1000g.haplotypes";
		Logger log = new Logger(outDir + "log.log");
		String parsedHaps = parseHaplogroupFile(haps, outDir, log);

		VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.ANY, log);
		vpop.report();

		String caseControlVcf = filterVpop(inputVCF, outDir, vpop, cases, controls, log);
		String hqVcf = filterCR(caseControlVcf, 0.95, outDir, pseqDir, log);
		double[] mafs = new double[] { 0.01 };
		for (VARIANT_TYPES vType : VARIANT_TYPES.values()) {
			for (VARIANT_SETS vSets : VARIANT_SETS.values()) {
				for (double d : mafs) {
					String mafVcf = filterMaf(hqVcf, d, outDir, log);
					runType(pseqDir, vType, vSets, outDir, cases, controls, log, vpop, mafVcf, parsedHaps);
				}
			}
		}
	}

	private static void runType(String pseqDir, VARIANT_TYPES vtTypes, VARIANT_SETS vSets, String outDir, String cases,
			String controls, Logger log, VcfPopulation vpop, String hqVcf, String haps) {
		String setVCF = filter(hqVcf, vtTypes.flag, vtTypes.filt, outDir, log);
		runIstats(setVCF, pseqDir, outDir, log);

		SNPInfo nsUniprot = generateSNPInfo(setVCF, outDir, vSets.flag, vSets.force);

		String[] plinks = convertToPlink(setVCF, outDir, log);
		addPhenoToFam(plinks[2], vpop, cases, controls);
		runSkat(plinks, nsUniprot, haps, log);
	}

}
