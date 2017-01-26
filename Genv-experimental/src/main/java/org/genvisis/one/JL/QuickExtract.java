package org.genvisis.one.JL;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.manage.VCOps;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * @author lane0212 Quickly looking at few variants in specific regions and pulling out some info of
 *         interest
 */
public class QuickExtract {

	public static void main(String[] args) {
		String vcf = "/home/tsaim/shared/Project_Tsai_Project_026/vcf/joint_genotypes_tsai_21_25_26_spector.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.sed.EPP.vcf.gz";
		String vpop = "/home/tsaim/shared/Project_Tsai_21_25_26_Spector_Joint/EPP/EPP.vpop";
		String thing = "/home/tsaim/shared/Project_Tsai_Project_026/vcf/HM_Pos.txt";
		String out = VCFOps.getAppropriateRoot(vcf, false) + ".alittleSummary.txt";

		// String[] segsS = HashVec.loadFileToStringArray(thing, true, new int[] { 2 }, true);
		// Segment[] segs = new Segment[segsS.length];
		// for (int i = 0; i < segs.length; i++) {
		// segs[i] = new Segment(segsS[i]);
		// }
		Logger log = new Logger(ext.rootOf(out) + ".log");
		VcfPopulation vpoppeer = VcfPopulation.load(vpop, POPULATION_TYPE.CASE_CONTROL, log);

		LocusSet<Segment> set = LocusSet.loadSegmentSetFromFile(thing, 0, 1, 2, 0, true, true, 0, log);

		VCFFileReader reader = new VCFFileReader(new File(vcf), true);
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(out));
			writer.println("CHROM\tPOS\tREF\tALT\tNumAlt\tNumAltMinus3\t1000G\tSNPEFF_IMPACT\tSNPEFF_GENE_NAME");
			int num = 0;
			int found = 0;
			for (VariantContext vc : reader) {
				num++;
				if (num % 10000 == 0) {
					log.reportTimeInfo(num + " scanned");
					log.reportTimeInfo(found + " found");
				}

				Segment vSegment = VCOps.getSegment(vc);
				int[] indices = set.getOverlappingIndices(vSegment);
				if (indices != null) {
					VariantContext vcSub = VCOps.getSubset(vc, vpoppeer.getSubPop().get(VcfPopulation.CASE));
					int numAlt = vc.getHetCount() + vc.getHomVarCount();
					int numAltMinux = vcSub.getHetCount() + vcSub.getHomVarCount();
					writer.println(vc.getContig()	+ "\t" + vc.getStart() + "\t"
													+ vc.getReference().getBaseString() + "\t"
													+ vc.getAlternateAlleles().toString() + "\t" + numAlt + "\t" + numAltMinux
													+ "\t"
													+ ArrayUtils.toStr(VCOps.getAnnotationsFor(new String[] {"g10002014oct_all",
																																							"SNPEFF_IMPACT",
																																							"SNPEFF_GENE_NAME"},
																																vc, ".")));
					found++;
				}
			}

			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + out);
			log.reportException(e);
		}
		reader.close();

	}
}
