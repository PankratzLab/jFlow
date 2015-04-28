package seq.manage;

import javax.swing.text.Position;

import cnv.var.LocusSet;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import seq.manage.BamPileUp.PILE_TYPE;
import seq.manage.VCOps.VC_SUBSET_TYPE;
import seq.qc.FilterNGS;
import seq.qc.FilterNGS.SAM_FILTER_TYPE;
import seq.qc.FilterNGS.VARIANT_FILTER_BOOLEAN;
import seq.qc.FilterNGS.VARIANT_FILTER_DOUBLE;
import common.Logger;
import common.Positions;
import common.ext;
import filesys.Segment;

public class BamImport {

	private BamPileUp bamPileUp;
	private String vcf;
	private String bam;
	private FilterNGS readDepthFilter;
	private LocusSet<Segment> intervals;
	private ReferenceGenome referenceGenome;
	private Logger log;

	public BamImport(String vcf, String bam, FilterNGS readDepthFilter, LocusSet<Segment> intervals, ReferenceGenome referenceGenome, Logger log) {
		super();
		this.vcf = vcf;
		this.bam = bam;
		this.readDepthFilter = readDepthFilter;
		this.referenceGenome = referenceGenome;
		this.bamPileUp = new BamPileUp(bam, referenceGenome, 1, readDepthFilter, intervals.getLoci(), PILE_TYPE.REGULAR, SAM_FILTER_TYPE.COPY_NUMBER, log);
		this.intervals = new LocusSet<Segment>(BamOps.converQItoSegs(bamPileUp.getQueryIntervals(), BamOps.getHeader(bam), log), true, log) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;
		};

		this.log = log;
	}

	public void importBam() {
		VCFFileReader reader = new VCFFileReader(vcf, true);
		String bamSample = BamOps.getSampleName(bam);

		if (ext.indexOfStr(bamSample, VCFOps.getSamplesInFile(reader)) < 0) {
			log.reportTimeError("Could not find sample " + bamSample + " in the vcf " + vcf);
			return;
		} else {
			log.reportTimeInfo("Detected sample " + bamSample + " in vcf " + vcf);
		}
		log.reportTimeWarning("Only un-ambigous and biallelic variants will be imported from " + vcf);
		FilterNGS.VariantContextFilter niceAllele = new FilterNGS.VariantContextFilter(new VARIANT_FILTER_DOUBLE[] {}, new VARIANT_FILTER_BOOLEAN[] { VARIANT_FILTER_BOOLEAN.BIALLELIC_FILTER, VARIANT_FILTER_BOOLEAN.AMBIGUOUS_FILTER }, null, null, log);
		SampleNGS ngsSample = new SampleNGS(bamSample);
		Segment curSegment = null;
		while (bamPileUp.hasNext()) {

			BamPile bamPile = bamPileUp.next();
			Segment curSeg = bamPile.getBin();
			Segment[] overlaps = intervals.getOverLappingLoci(curSeg);
			if (overlaps == null || overlaps.length == 0) {
				log.reportTimeInfo("Un-matched segments");
				return;
			}
			if (overlaps.length > 1) {
				log.reportTimeInfo("Non-unique (overlapping) segments were supplied, halting");
				return;
			}
			CloseableIterator<VariantContext> reg = reader.query(Positions.getChromosomeUCSC(overlaps[0].getChr(), true), overlaps[0].getStart(), overlaps[0].getStop());
			while (reg.hasNext()) {
				VariantContext vc = reg.next();
				if (VCOps.getSegment(vc).equals(curSeg) && niceAllele.filter(vc).passed()) {
					ngsSample.addGeno(VCOps.getGenotypeFor(vc, bamSample, VC_SUBSET_TYPE.SUBSET_STRICT), log);
				}
			}
		}
	}

	private static class tmpGenotype {

	}

	public static void test() {
		String bam = "D:/data/Project_Tsai_Project_021/testPileUp/rrd_lane_CONTROL_4_CTCTCTAC-CTCTCTAT.merge.sorted.dedup.realigned.bam";
		String ref = "C:/bin/ref/hg19_canonical.fa";
		String segs = "C:/bin/Agilent/captureLibraries/SureSelectHumanAllExonV5UTRs/AgilentCaptureRegions_chr1.txt";
		String vcf = "D:/data/Project_Tsai_Spector_Joint/joint_genotypes_tsai_21_25_spector_mt.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.vcf.gz";
		Logger log = new Logger();
		Segment[] q = segs == null ? null : Segment.loadRegions(segs, 0, 1, 2, 0, true, true, true, 0);
		FilterNGS filterNGS = new FilterNGS(0, 0, new int[] { 0, 0 });
		ReferenceGenome referenceGenome = ref == null ? null : new ReferenceGenome(ref, log);
		LocusSet<Segment> locusSet = new LocusSet<Segment>(q, true, log) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};
		BamImport bamImport = new BamImport(vcf, bam, filterNGS, locusSet, referenceGenome, log);
		bamImport.importBam();

	}

	public static void main(String[] args) {
		String bam = "D:/data/Project_Tsai_Project_021/testPileUp/rrd_lane_CONTROL_4_CTCTCTAC-CTCTCTAT.merge.sorted.dedup.realigned.bam";
		String ref = "C:/bin/ref/hg19_canonical.fa";
		String segs = "C:/bin/Agilent/captureLibraries/SureSelectHumanAllExonV5UTRs/AgilentCaptureRegions_chr1.txt";
		String vcf = "D:/data/Project_Tsai_Spector_Joint/joint_genotypes_tsai_21_25_spector_mt.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.vcf.gz";
		test();
	}
}
