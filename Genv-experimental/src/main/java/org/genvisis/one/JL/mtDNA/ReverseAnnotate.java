package org.genvisis.one.JL.mtDNA;

import java.io.File;

import org.genvisis.CLI;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.VCOps;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class ReverseAnnotate {

	private static void run(String vcf, String regions) {
		VCFFileReader reader = new VCFFileReader(new File(vcf), false);
		String out = ext.parseDirectoryOfFile(vcf) + ext.rootOf(regions) + "annot.txt";
		Segment[] segs = Segment.loadRegions(regions, 0, 1, 2, true);
		Logger log = new Logger(out + ".log");
		StringBuilder builder = new StringBuilder("#CHR\tBP1\tBP2\tID");
		for (VariantContext vc : reader) {
			for (Segment seg : segs) {
				if (VCOps.getSegment(vc).overlaps(seg)) {
					builder.append(seg.getChromosomeUCSC() + "\t" + seg.getStart() + "\t" + seg.getStop()
							+ VCOps.getAnnotationsFor(new String[] { "" }, vc, "NOTHING"));
				}
			}
		}
		Files.write(builder.toString(), out);

	}

	public static void main(String[] args) {
		CLI c = new CLI(ReverseAnnotate.class);
		String vcf = "a.vcf";
		String regions = "a.regions";

		c.addArgWithDefault("vcf", "vcf", vcf);
		c.addArgWithDefault("regions", "regions", regions);

		c.parseWithExit(args);
		run(c.get("vcf"), c.get("regions"));
	}

}
