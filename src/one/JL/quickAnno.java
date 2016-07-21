package one.JL;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import common.Logger;
import common.ext;
import seq.analysis.GATK_Genotyper;
import seq.manage.VCFOps;

public class quickAnno {

	public static void main(String[] args) {
		String file = "/home/tsaim/shared/Project_Tsai_21_25_26_28_Spector_Joint/mitochondrial/joint_genotypes_tsai_21_25_26_28_spector.chrM.vcf";
		Logger log = new Logger(ext.parseDirectoryOfFile(file) + "mitoLog");
		log.reportTimeInfo("Converting positions of " + file);

		String oneMinus = adjustPositions(file, -1, log);

		String anno = GATK_Genotyper.annotateOnlyWithDefualtLocations(oneMinus, null, true, false, new Logger());
		String finalVcf = adjustPositions(anno, 1, log);

	}

	private static String adjustPositions(String in, int adjust, Logger log) {
		VCFFileReader reader = new VCFFileReader(in, false);
		String output = VCFOps.getAppropriateRoot(in, false)+ ".posAdjust_" + adjust+".vcf";
		VariantContextWriter writer = VCFOps.initWriterWithHeader(reader, output, VCFOps.DEFUALT_WRITER_OPTIONS, log);

		for (VariantContext vc : reader) {
			VariantContextBuilder builder = new VariantContextBuilder(vc);
			builder.start(vc.getStart() + adjust);
			builder.stop(vc.getEnd() + adjust);
			if(!vc.getContig().equals("chrM")){
				throw new IllegalArgumentException("why are you using this method");
			}
			writer.add(builder.make());
		}
		return output;
	}

}
