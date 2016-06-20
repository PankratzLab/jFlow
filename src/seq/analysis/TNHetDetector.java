package seq.analysis;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;

import seq.manage.GenotypeOps;
import seq.manage.VCFOps;
import seq.manage.VCOps;
import seq.manage.VCFOps.HEADER_COPY_TYPE;
import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;
import seq.manage.VCOps.VC_SUBSET_TYPE;
import common.Logger;
import common.ext;

/**
 * @author lane0212 fairly specific class that uses the results of {@link VCFSimpleTally} to find genes with potential compound hets in tumor normal subsets<br>
 *         note this assumes pretty small vcf sizes<br>
 *         Kinda a hack
 * 
 *
 */
public class TNHetDetector {

	private static void filter(String vcfTumor, String vcfNormals, String vcfPop, String outputDir, String omimDir, String[] otherGenesOfInterest, double maf, Logger log) {
		log.reportTimeInfo("Building normal hash");
		VcfPopulation vpop = VcfPopulation.load(vcfPop, POPULATION_TYPE.TUMOR_NORMAL, log);
		Hashtable<String, Genotype> normalHash = new Hashtable<String, Genotype>();
		Hashtable<String, String> match = matchSamples(vpop, log);

		VcfPopulation.splitVcfByPopulation(vcfNormals, vcfPop, true, true, false, log);
		VCFFileReader readerNormal = new VCFFileReader(vcfNormals, true);

		for (VariantContext vc : readerNormal) {

			VariantContext vcNormal = VCOps.getSubset(vc, vpop.getSuperPop().get(VcfPopulation.NORMAL), VC_SUBSET_TYPE.SUBSET_STRICT);
			for (String normalSamp : vcNormal.getSampleNames()) {
				Genotype g = vcNormal.getGenotype(normalSamp);

				if (!g.isHomRef() && g.isCalled() && ext.indexOfStr(VCOps.getSNP_EFFImpact(vc), VCFSimpleTally.EFF) >= 0) {
					normalHash.put(normalSamp + "_" + VCOps.getSNP_EFFGeneName(vc), g);
				}
			}
		}
		readerNormal.close();

		VCFFileReader readerTumor = new VCFFileReader(vcfTumor, true);
		String output = outputDir + VCFOps.getAppropriateRoot(vcfTumor, true) + ".HET.vcf.gz";
		VariantContextWriter writer = VCFOps.initWriter(output, VCFOps.DEFUALT_WRITER_OPTIONS, readerTumor.getFileHeader().getSequenceDictionary());
		VCFOps.copyHeader(readerTumor, writer, vpop.getSuperPop().get(VcfPopulation.TUMOR), HEADER_COPY_TYPE.FULL_COPY, log);
		int num = 0;
		for (VariantContext vc : readerTumor) {
			VariantContext vcTumor = VCOps.getSubset(vc, vpop.getSuperPop().get(VcfPopulation.TUMOR), VC_SUBSET_TYPE.SUBSET_STRICT);
			VariantContextBuilder builder = new VariantContextBuilder(vcTumor);
			ArrayList<Genotype> genotypes = new ArrayList<Genotype>();
			for (String tumorSamp : vcTumor.getSampleNames()) {
				Genotype g = vcTumor.getGenotype(tumorSamp);
				boolean add = false;
				GenotypeBuilder bGenotypeBuilder = new GenotypeBuilder(g);
				bGenotypeBuilder.alleles(GenotypeOps.getNoCall());
				if (!g.isHomRef() && g.isCalled()) {
					String key = match.get(tumorSamp) + "_" + VCOps.getSNP_EFFGeneName(vc);
					if (normalHash.containsKey(key)) {
						if (!g.sameGenotype(normalHash.get(key))) {
							add = true;
							genotypes.add(g);
						} else {
							genotypes.add(bGenotypeBuilder.make());
							log.reportTimeInfo("Removing identical geno");
						}
					}
				} else {
					genotypes.add(bGenotypeBuilder.make());
				}
				if (add) {
					builder.genotypes(genotypes);
					writer.add(builder.make());
					num++;
				}
			}
		}
		log.reportTimeInfo("Found " + num + " potential candidates");
		readerTumor.close();
		writer.close();

		String outCase = outputDir + VcfPopulation.TUMOR + ".vpop";
		VcfPopulation tumorCase = developTNVpopCase(vpop);
		tumorCase.dump(outCase);

		VCFSimpleTally.test(output, new String[] { outCase }, omimDir, otherGenesOfInterest, null, maf, true, true,null);
	}

	private static VcfPopulation developTNVpopCase(VcfPopulation vpop) {
		Set<String> tmpSet = new HashSet<String>();
		Set<String> tmpSet2 = new HashSet<String>();

		Hashtable<String, Set<String>> casePop = new Hashtable<String, Set<String>>();

		casePop.put(VcfPopulation.TUMOR, tmpSet);
		casePop.put(VcfPopulation.NORMAL, tmpSet2);

		casePop.get(VcfPopulation.TUMOR).addAll(vpop.getSuperPop().get(VcfPopulation.TUMOR));
		casePop.get(VcfPopulation.NORMAL).addAll(vpop.getSuperPop().get(VcfPopulation.NORMAL));

		VcfPopulation caseTumor = new VcfPopulation(casePop, casePop, POPULATION_TYPE.CASE_CONTROL, new Logger());
		return caseTumor;
	}

	private static Hashtable<String, String> matchSamples(VcfPopulation vpop, Logger log) {
		Hashtable<String, String> matched = new Hashtable<String, String>();

		for (String tnPair : vpop.getSubPop().keySet()) {
			Set<String> samps = vpop.getSubPop().get(tnPair);
			String tumor = null;
			String normal = null;
			for (String samp : samps) {
				if (vpop.getPopulationForInd(samp, RETRIEVE_TYPE.SUPER)[0].equals(VcfPopulation.TUMOR)) {
					tumor = samp;
				} else if (vpop.getPopulationForInd(samp, RETRIEVE_TYPE.SUPER)[0].equals(VcfPopulation.NORMAL)) {
					normal = samp;
				} else {
					throw new IllegalArgumentException("Unknown types");
				}
			}
			matched.put(tumor, normal);
		}
		return matched;

	}

	private static void test() {
		// String dir = "/home/tsaim/shared/Project_Tsai_21_25_26_28_Spector_Joint/vcf/Freq/";
		String dir = "D:/data/Project_Tsai_21_25_26_28_spector/TumorNormal/hets/";
		String omimDir = "C:/bin/ref/OMIM/";
		String cosmic = "C:/bin/ref/cancer_gene_census.txt";
		String gdi = "C:/bin/ref/GDI.txt";

		String vcfTumor = dir + "CUSHINGS_TUMOR.maf_0.001.final.vcf.gz";
		String vcfNormal = dir + "CUSHINGS_TUMOR_CONTROL.maf_0.001.final.CUSHING_FREQ.vcf.gz";

		String outDir = dir + "TN_HET/";
		String vcfPop = outDir + "TN.vpop";
		Logger log = new Logger(outDir + "TNHet.log");
		filter(vcfTumor, vcfNormal, vcfPop, outDir, omimDir, new String[] { cosmic, gdi }, .0001, log);

	}

	public static void main(String[] args) {
		test();
	}

}
