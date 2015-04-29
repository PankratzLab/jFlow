package seq.analysis;

import filesys.Segment;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Set;

import common.Array;
import common.Files;
import common.Logger;
import common.ext;
import seq.manage.VCFOps;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;
import seq.manage.VCOps;
import seq.manage.VCFOps.HEADER_COPY_TYPE;
import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCOps.VC_SUBSET_TYPE;
import seq.qc.FilterNGS;
import seq.qc.FilterNGS.VARIANT_FILTER_BOOLEAN;
import seq.qc.FilterNGS.VARIANT_FILTER_DOUBLE;
import seq.qc.FilterNGS.VariantContextFilter;

public class RegNovo {
	private static final String[] HEADER = new String[] { "Sample", "ID", "CHR", "Start", "Stop", "PassesPop", "RefDepth", "AltDepth", "GQ", "Pop_AAC", "AllParents_AAC", "AllOffspring_AAC", "P1RefDepth", "P1AltDepth", "P1GQ", "P1Filter", "P1_AAC", "P2RefDepth", "P2AltDepth", "P2GQ", "P2Filter", "P2_AAC", "GATK_FILTER" };
	private static final String[] TO_REPORT = new String[] { "SNPEFF_GENE_NAME", "SNPEFF_EFFECT", "SNPEFF_IMPACT", "AAChange.refGene", "SNPEFF_EXON_ID", "esp6500si_all", "g10002014oct_all", "culprit" };
	public static final String OFFSPRING = "OFFSPRING";
	private String vcf;
	private VcfPopulation vpop;
	private VariantContextFilter filterControl;
	private VariantContextFilter filterOffspring;
	private FilterNGS readDepths;
	private String outputDir;
	private Logger log;

	public RegNovo(String vcf, VcfPopulation vpop, VariantContextFilter filterControl, VariantContextFilter filterOffspring, FilterNGS readDepths, String outputDir, Logger log) {
		super();
		this.vcf = vcf;
		this.vpop = vpop;
		this.filterControl = filterControl;
		this.filterOffspring = filterOffspring;
		this.outputDir = outputDir;
		new File(outputDir).mkdirs();
		this.readDepths = readDepths;
		this.log = log;
	}

	public void scanForDenovo() {
		String outputVCF = outputDir + ext.rootOf(vcf) + ".denovo.vcf.gz";
		String outputSummary = outputDir + ext.rootOf(vcf) + ".denovoSummary.txt";
		if (!vpop.getSuperPop().containsKey(OFFSPRING) || !vpop.getSuperPop().containsKey(VcfPopulation.CONTROL)) {
			log.reportTimeError(vpop.getFileName() + " must contain both " + OFFSPRING + " and " + VcfPopulation.CONTROL);
			return;
		}

		VCFFileReader reader = new VCFFileReader(vcf, true);
		VariantContextWriter writer = VCFOps.initWriter(outputVCF, null, reader.getFileHeader().getSequenceDictionary());
		VCFOps.copyHeader(reader, writer, null, HEADER_COPY_TYPE.FULL_COPY, log);
		Set<String> offspring = vpop.getSuperPop().get(OFFSPRING);
		Set<String> controls = vpop.getSuperPop().get(VcfPopulation.CONTROL);
		HashSet<String> all = new HashSet<String>();
		all.addAll(controls);
		all.addAll(offspring);
		PrintWriter summaryWriter = Files.getAppropriateWriter(outputSummary);
		summaryWriter.println(Array.toStr(HEADER) + "\t" + Array.toStr(TO_REPORT));
		int total = 0;
		int denovo = 0;
		int denovoExtra = 0;
		for (VariantContext vc : reader) {
			total++;
			VariantContext vcControls = VCOps.getSubset(vc, controls, VC_SUBSET_TYPE.SUBSET_STRICT);
			for (String off : offspring) {
				VariantContext vcOff = VCOps.getSubset(vc, off, VC_SUBSET_TYPE.SUBSET_STRICT);
				if (filterOffspring.filter(vcOff).passed()) {
					if (VCOps.getAltAlleleContext(vcOff, readDepths.getReadDepthFilter()[1],log).getSampleNames().size() > 0) {
						String[] fam = vpop.getPopulationForInd(off, RETRIEVE_TYPE.SUB);
						for (int i = 0; i < fam.length; i++) {
							boolean passesAll = true;
							Set<String> curFam = vpop.getSubPop().get(fam[i]);
							if (passesAll) {
								passesAll = curFam.size() == 3;
							}
							String pString = "";
							HashSet<String> parents = new HashSet<String>();
							for (String famInd : curFam) {
								if (!famInd.equals(off)) {
									parents.add(famInd);
									VariantContext vcFam = VCOps.getSubset(vc, famInd, VC_SUBSET_TYPE.SUBSET_STRICT);
									pString += "\t" + vcFam.getGenotypes().get(0).getAD()[0];
									pString += "\t" + vcFam.getGenotypes().get(0).getAD()[1];
									pString += "\t" + vcFam.getGenotypes().get(0).getGQ() + "";
									pString += "\t" + filterControl.filter(vcFam).getTestPerformed() + "\t" + VCOps.getAAC(vcFam, null);

								}
							}
							if (VCOps.getAAC(vcOff, null) > VCOps.getAAC(vc, parents)) {
								String passesAllControls = filterControl.filter(vcControls).getTestPerformed();
								denovo++;

								if (filterControl.filter(vcControls).passed()) {
									denovoExtra++;
								}
								writer.add(vc);
								String[] anno = VCOps.getAnnotationsFor(TO_REPORT, vc, ".");

								String sum = "\t" + vcOff.getGenotypes().get(0).getAD()[0];
								sum += "\t" + vcOff.getGenotypes().get(0).getAD()[1];
								sum += "\t" + vcOff.getGenotypes().get(0).getGQ() + "";
								sum += "\t" + VCOps.getAAC(VCOps.getAltAlleleContext(vc, readDepths.getReadDepthFilter()[0],log), all);
								sum += "\t" + VCOps.getAAC(VCOps.getAltAlleleContext(vc, readDepths.getReadDepthFilter()[0],log), controls);
								sum += "\t" + VCOps.getAAC(VCOps.getAltAlleleContext(vc, readDepths.getReadDepthFilter()[0],log), offspring);

								summaryWriter.println(off + "\t" + vcOff.getID() + "\t" + vc.getChr() + "\t" + vc.getStart() + "\t" + vc.getEnd() + "\t" + passesAllControls + sum + pString + "\t" + vc.getFilters().toString() + "\t" + Array.toStr(anno));
								summaryWriter.flush();
							}
						}
					}
				}
			}
			if (total % 10000 == 0) {
				log.reportTimeInfo("Scanned " + total + " variants and detected" + denovo + " denovos (denovoextra=" + denovoExtra + ")");
			}
		}
		reader.close();
		writer.close();
		summaryWriter.close();
	}

	public void scanForSegs(Segment[] segs, VariantContextFilter variantContextFilter, String root) {
		String output = outputDir + root + ".extracted.txt";
		String outputVCF = outputDir + root + ".extracted.vcf";
		Set<String> offspring = vpop.getSuperPop().get(OFFSPRING);
		Set<String> controls = vpop.getSuperPop().get(VcfPopulation.CONTROL);
		Set<String> all = new HashSet<String>();
		all.addAll(offspring);
		all.addAll(controls);

		PrintWriter writer = Files.getAppropriateWriter(output);
		writer.println("CHR\tSTART\tSTOP\tNUM_OFFSPRING\tAAC_OFFSPRING\tNUM_RENTS\tAAC_RENTS\t" + Array.toStr(TO_REPORT) + "\tFILTER\tHQ_Pop\tHQ_ALT\t\tINFO");
		VCFFileReader reader = new VCFFileReader(vcf, true);
		int count = 0;
		VariantContextWriter writerVC = VCFOps.initWriter(outputVCF, null, reader.getFileHeader().getSequenceDictionary());
		VCFOps.copyHeader(reader, writerVC, null, HEADER_COPY_TYPE.FULL_COPY, log);
		for (VariantContext vc : reader) {
			count++;
			if (count % 100000 == 0) {
				log.reportTimeInfo("count " + count);
			}
			if (VCOps.isInTheseSegments(vc, segs)) {
				VariantContext vcSub = VCOps.getSubset(vc, all);

				if (VCOps.getAltAlleleContext(vcSub, 0,log).getSampleNames().size() > 0) {
					Segment vcSeg = VCOps.getSegment(vc);
					String out = "";
					out += vcSeg.getChr();
					out += "\t" + vcSeg.getStart();
					out += "\t" + vcSeg.getStop();
					out += "\t" + VCOps.getAltAlleleContext(VCOps.getSubset(vc, offspring), 0,log).getSampleNames().size();
					out += "\t" + VCOps.getAAC(vcSub, offspring);
					out += "\t" + VCOps.getAltAlleleContext(VCOps.getSubset(vc, controls), 0,log).getSampleNames().size();
					out += "\t" + VCOps.getAAC(vcSub, controls);
					out += "\t" + Array.toStr(VCOps.getAnnotationsFor(TO_REPORT, vcSub, "."));
					out += "\t" + vc.getFilters().toString();
					out += "\t" + variantContextFilter.filter(vcSub).getTestPerformed();
					out += "\t" + variantContextFilter.filter(VCOps.getAltAlleleContext(vcSub, 0,log)).getTestPerformed();
					writerVC.add(vcSub);
					writer.println(out + "\t" + vcSub.toStringWithoutGenotypes());
				}
			}
		}
		reader.close();
		writer.close();
	}

	private static VariantContextFilter getQualityFilter(Logger log) {
		VARIANT_FILTER_DOUBLE callRate = VARIANT_FILTER_DOUBLE.CALL_RATE;
		VARIANT_FILTER_DOUBLE gq = VARIANT_FILTER_DOUBLE.GQ;
		gq.setDFilter(80);
		VARIANT_FILTER_DOUBLE dp = VARIANT_FILTER_DOUBLE.DP;
		dp.setDFilter(8);
		VARIANT_FILTER_DOUBLE vqslod = VARIANT_FILTER_DOUBLE.VQSLOD_LOOSE;
		VARIANT_FILTER_BOOLEAN fail = VARIANT_FILTER_BOOLEAN.FAILURE_FILTER;
		VARIANT_FILTER_DOUBLE[] qualFilts = new VARIANT_FILTER_DOUBLE[] { callRate, gq, vqslod };
		VariantContextFilter vContextFilter = new VariantContextFilter(qualFilts, new VARIANT_FILTER_BOOLEAN[] { fail }, null, null, log);
		return vContextFilter;
	}

	public static void detectDenovo(String vcf, String vpopFile, Logger log) {
		VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.ANY, log);
		vpop.report();
		FilterNGS readDepths = new FilterNGS(0, 0, new int[] { 0, 5 });
		RegNovo regNovo = new RegNovo(vcf, vpop, getQualityFilter(log), getQualityFilter(log), readDepths, ext.parseDirectoryOfFile(vpopFile), log);

		Segment[] segs = new Segment[] { new Segment("chr17:7571520-7590968") };
		// regNovo.scanForSegs(segs, getQualityFilter(log), "TP53");
		regNovo.scanForDenovo();
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "RegNovo.dat";
		String vcf = "D:/data/Project_Tsai_Spector_Joint/joint_genotypes_tsai_21_25_spector_mt.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.vcf.gz";
		String vpopFile = "D:/data/logan/OSv2_seq/RegNovo/OsSamps.vcfPop_outliersRemoved.txt";
		String logfile = null;
		Logger log;

		String usage = "\n" + "seq.analysis.RegNovo requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

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
			detectDenovo(vcf, vpopFile, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
