package org.genvisis.one.JL.mica;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;

import org.genvisis.CLI;
import org.genvisis.cnv.filesys.ABLookup;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.analysis.Blast;
import org.genvisis.seq.analysis.Blast.FastaEntry;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCOps;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class MStr {

	private static void run(String[] vcfs, String outDir) {

		for (String vcf : vcfs) {
			new File(outDir).mkdirs();
			Segment mSeg = new Segment("chr6:31380150-31380165");
			String outfile = outDir + ext.rootOf(vcf) + "_"
					+ ext.replaceWithLinuxSafeCharacters(mSeg.getUCSClocation(), true) + ".matrix.txt";
			PrintWriter writer = Files.getAppropriateWriter(outfile);
			Logger log = new Logger(outDir + "log.log");
			VCFFileReader reader = new VCFFileReader(new File(vcf), false);
			String[] samples = VCFOps.getSamplesInFile(reader);
			boolean[] samplesWithOneGenotype = Array.booleanArray(samples.length, false);
			writer.println("CHR\tPOS\tREF\tFULL_ALT\tSAMPLE\tA1\tA2\tGQ\tAD\tHOM_HET\tFullGeno");
			for (VariantContext vc : reader) {
				if (VCOps.getSegment(vc).overlaps(mSeg)) {

					if (vc.isIndel()) {
						int index = 0;
						StringBuilder builder = new StringBuilder(vc.getContig() + "\t" + vc.getStart() + "\t"
								+ vc.getReference().getBaseString() + "\t" + vc.getAlternateAlleles().toString());
						for (Genotype g : vc.getGenotypes()) {

							if (!g.getSampleName().equals(samples[index])) {
								writer.close();

								throw new IllegalStateException("Sample mismatch");
							}
							StringBuilder builder2 = new StringBuilder(builder.toString());
							builder2.append("\t" + g.getSampleName());

							if (g.isCalled()) {
								samplesWithOneGenotype[index] = true;

								ArrayList<String> alts = new ArrayList<String>();
								String ref = null;
								for (Allele a : g.getAlleles()) {
									if (a.isNonReference()) {
										alts.add(a.getBaseString());
									} else {
										ref = a.getBaseString();
									}

								}
								if (ref != null) {
									if (alts.isEmpty()) {
										builder2.append("\t" + ref + "\t" + ref);
									} else {
										builder2.append("\t" + ref + "\t" + alts.get(0));
										if (alts.size() != 1) {
											writer.close();

											throw new IllegalStateException("geno mismatch");
										}
									}
								} else {
									if (alts.size() != 2) {
										writer.close();

										throw new IllegalStateException("geno mismatch");

									}
									for (String al : alts) {
										builder2.append("\t" + al);
									}
								}
							} else {
								builder2.append("\t.\t.");
							}
							builder2.append("\t" + g.getGQ() + "\t" + (g.getAD() == null ? g.getAnyAttribute("DPR")
									: Array.toStr(Array.toStringArray((g.getAD())), ",")));

							builder2.append("\t" + (g.isCalled() ? (g.isHom() ? "HOM" : "HET") : "NA"));
							builder2.append("\t" + g.toString());
							if (!g.getSampleName().contains("H20")) {
								writer.println(builder2.toString());
							}
							index++;
						}
					}
				}
			}
			writer.close();
			reader.close();
			log.reportTimeInfo(
					Array.booleanArraySum(samplesWithOneGenotype) + " of " + samples.length + " samples had genotypes");
		}

	}

	public static void test() {
		String fastaDb = "/Volumes/Beta/ref/hg19_canonical.fa";
		Blast blast = new Blast(fastaDb, 60, 100, new Logger(), true, true);
		FastaEntry fastaEntry = new FastaEntry("HDSIF",
				"GAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCTCATTCTATTATTTATCGCACCTACGTTCAATATTACAGGCGAACATACCTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAAT");
		blast.blastSequence(new FastaEntry[] { fastaEntry }, null);
	}

	public static void main(String[] args) {

		CLI c = new CLI(ABLookup.class);
		c.addArgWithDefault("vcfs", "vcf files to explore,comma-delimited list", null);
		c.addArgWithDefault(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, null);
		c.parseWithExit(args);
		run(c.get("vcfs").split(","), c.get(CLI.ARG_OUTDIR));
	}

}
