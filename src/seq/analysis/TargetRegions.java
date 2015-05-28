package seq.analysis;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Set;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import common.Array;
import common.Logger;
import common.Positions;
import common.ext;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import seq.manage.VCFOps;
import seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;
import seq.manage.VCOps;
import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCOps.ALT_ALLELE_CONTEXT_TYPE;
import seq.manage.VCOps.COMMON_INFO;
import seq.manage.VCOps.GENOTYPE_INFO;
import seq.manage.VCOps.VC_SUBSET_TYPE;
import cnv.var.LocusSet;
import filesys.Segment;

public class TargetRegions<T extends Segment> {
	private static final String[] TO_REPORT = new String[] { "SNPEFF_GENE_NAME", "SNPEFF_EFFECT", "SNPEFF_IMPACT", "AAChange.refGene", "SNPEFF_EXON_ID", "culprit", "snp138", "esp6500si_all", "g10002014oct_all" };

	private String vcfFile;
	private LocusSet<T> targetRegions;
	private VcfPopulation vpop;
	private Logger log;

	public TargetRegions(String vcfFile, LocusSet<T> targetRegions, VcfPopulation vpop, Logger log) {
		super();
		this.vcfFile = vcfFile;
		this.targetRegions = targetRegions;
		this.vpop = vpop;
		this.log = log;
	}

	public void summarizeRegions(String fullPathToOutput, String toMatchVCF, String[] toMatchAnnotations) {
		VCFFileReader reader = new VCFFileReader(vcfFile, true);
		T[] regions = targetRegions.getLoci();

		String[] subpop = vpop.getSubPop().keySet().toArray(new String[vpop.getSubPop().keySet().size()]);
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(fullPathToOutput));
			writer.print("CHR\tStart\tStop\tRef\tAlt\tFILTER");
			for (int i = 0; i < subpop.length; i++) {
				writer.print("\t" + subpop[i] + "_AAC");
				writer.print("\t" + subpop[i] + "_SamplesWithAlt");
				writer.print("\t" + subpop[i] + "_NumExcluded");

			}
			writer.print("\t" + Array.toStr(TO_REPORT));
			if (toMatchVCF != null && toMatchAnnotations != null) {
				for (int i = 0; i < toMatchAnnotations.length; i++) {
					writer.print("\t" + ext.rootOf(toMatchVCF) + "_" + toMatchAnnotations[i]);
				}
			}
			for (int i = 0; i < subpop.length; i++) {
				writer.print("\t" + subpop[i] + "_AVG_DP");
				writer.print("\t" + subpop[i] + "_AVG_GQ");
			}
			writer.println();
			for (int i = 0; i < regions.length; i++) {
				CloseableIterator<VariantContext> cIterator = reader.query(Positions.getChromosomeUCSC(regions[i].getChr(), true), regions[i].getStart(), regions[i].getStop());

				while (cIterator.hasNext()) {
					VariantContext vc = cIterator.next();
					writer.print(Positions.getChromosomeUCSC(VCOps.getSegment(vc).getChr(), true));
					writer.print("\t" + VCOps.getSegment(vc).getStart());
					writer.print("\t" + VCOps.getSegment(vc).getStop());
					writer.print("\t" + vc.getReference().getBaseString());
					writer.print("\t" + vc.getAlternateAlleles().toString());
					writer.print("\t" + vc.getFilters().toString());

					for (int j = 0; j < subpop.length; j++) {
						writer.print("\t" + VCOps.getAAC(vc, vpop.getSubPop().get(subpop[j])));
						VariantContext vcAlt = VCOps.getAltAlleleContext(VCOps.getSubset(vc, vpop.getSubPop().get(subpop[j])), null, null, ALT_ALLELE_CONTEXT_TYPE.ALL, log);
						Set<String> tmpSamps = vcAlt.getSampleNames();
						writer.print("\t");
						int index = 0;
						for (String altSamp : tmpSamps) {
							VariantContext curContext = VCOps.getSubset(vcAlt, altSamp, VC_SUBSET_TYPE.SUBSET_STRICT);
							writer.print((index > 0 ? "|" : "") + altSamp + ":GQ=" + curContext.getGenotype(0).getGQ() + ":AD=" + Array.toStr(curContext.getGenotype(0).getAD(), ","));
							index++;
						}
						int numExcluded = 0;
						for (String altSamp : vcAlt.getSampleNames()) {
							if (vpop.getPopulationForInd(altSamp, RETRIEVE_TYPE.SUPER)[0].equals(VcfPopulation.EXCLUDE)) {
								numExcluded++;
							}
						}
						writer.print("\t" + numExcluded);
					}
					writer.print("\t" + Array.toStr(VCOps.getAnnotationsFor(TO_REPORT, vc, ".")));
					if (toMatchVCF != null && toMatchAnnotations != null) {
						VariantContext vcMatch = VCFOps.lookupExactVariant(toMatchVCF, vc, log);
						if (vcMatch == null) {
							for (int j = 0; j < toMatchAnnotations.length; j++) {
								writer.print("\tNA");
							}
						} else {
							writer.print("\t" + Array.toStr(VCOps.getAnnotationsFor(toMatchAnnotations, vcMatch, ".")));
						}
					}
					for (int j = 0; j < subpop.length; j++) {
						writer.print("\t" + VCOps.getAvgGenotypeInfo(vc, vpop.getSubPop().get(subpop[j]), GENOTYPE_INFO.DP, log));
						writer.print("\t" + VCOps.getAvgGenotypeInfo(vc, vpop.getSubPop().get(subpop[j]), GENOTYPE_INFO.GQ, log));
					}

					writer.println();
				}
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + fullPathToOutput);
			log.reportException(e);
		}
		reader.close();
	}

	public static void test() {
		String vcf = "D:/data/Project_Tsai_21_25_26_spector/joint_genotypes_tsai_21_25_26_spector.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.vcf.gz";
		String vpopFile = "D:/data/Project_Tsai_21_25_26_spector/candidateGenes/USP8/usp8.vpop";
		String toMatchVCF = "D:/data/CHARGE/CHARGE_MAFS/charge_fibrinogen_mafs_and_counts.xln.hg19_multianno.eff.gatk.sed.vcf";
		String[] toMatchAnnotations = new String[] { "MAF_blacks", "MAF_whites" };
		String output = ext.parseDirectoryOfFile(vpopFile) + "targetRegions.txt";
		Logger log = new Logger(ext.parseDirectoryOfFile(vpopFile) + "target.log");
		VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.ANY, log);
		vpop.report();
		Segment[] segs = new Segment[] { new Segment("chr15:50714579-50795277") };
		LocusSet<Segment> set = new LocusSet<Segment>(segs, true, log) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;
		};
		TargetRegions<Segment> targetRegions = new TargetRegions<Segment>(vcf, set, vpop, log);
		targetRegions.summarizeRegions(output, toMatchVCF, toMatchAnnotations);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "TargetRegions.dat";
		String logfile = null;
		Logger log;

		String usage = "\n" + "seq.analysis.TargetRegions requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			log = new Logger(logfile);
			test();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
