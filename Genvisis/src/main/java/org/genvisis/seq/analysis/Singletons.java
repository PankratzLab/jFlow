package org.genvisis.seq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import org.genvisis.CLI;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCOps;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * Getting singleton estimates when lacking large number of samples
 *
 */
public class Singletons {


	private static void detectSingletons(	String[] vcfs, String outDir, int gqThreshold,
																				List<String> dotAnnos, HashSet<String> sampsToExclude) {
		new File(outDir).mkdirs();
		Logger log = new Logger(outDir + "singleton.log");
		for (String vcf : vcfs) {
			VCFFileReader reader = new VCFFileReader(new File(vcf), true);
			HashMap<String, Integer> singletonCounts = new HashMap<String, Integer>();
			HashMap<String, Integer> altCounts = new HashMap<String, Integer>();

			String[] samples = VCFOps.getSamplesInFile(vcf);
			for (int i = 0; i < samples.length; i++) {
				if (!sampsToExclude.contains(samples[i])) {
					singletonCounts.put(samples[i], 0);
					altCounts.put(samples[i], 0);
				}
			}
			int num = 0;
			for (VariantContext vcRaw : reader) {
				VariantContext vc = VCOps.getSubset(vcRaw, singletonCounts.keySet());
				num++;
				if (num % 10000 == 0) {
					log.reportTimeInfo("Scanned " + num + " variants");
				}
				if (!vc.isIndel() && !vc.isFiltered()) {
					for (Genotype g : vc.getGenotypes()) {
						if (g.isHet()) {
							altCounts.put(g.getSampleName(), altCounts.get(g.getSampleName()) + 1);
						} else if (g.isHomVar()) {
							altCounts.put(g.getSampleName(), altCounts.get(g.getSampleName()) + 2);
						}
						if (vc.getHetCount() + vc.getHomVarCount() == 1) {// singleton
																															// in
																															// VCF
							boolean found = false;
							if (g.isHet() || g.isHomVar()) {
								if (found) {
									reader.close();
									throw new IllegalStateException("Invalid singleton");
								} else {
									found = true;
									if (g.getGQ() > gqThreshold) {
										boolean externalSingleton = true;
										for (String dot : dotAnnos) {
											String anno = VCOps.getAnnotationsFor(new String[] {dot}, vc, ".")[0];
											if (!".".equals(anno)
													&& (anno.startsWith("rs") || Double.parseDouble(anno) > 0)) {
												externalSingleton = false;
											}
										}
										if (externalSingleton) {
											singletonCounts.put(g.getSampleName(),
																					singletonCounts.get(g.getSampleName()) + 1);
										}
									}
								}
							}
						}
					}
				}
			}
			reader.close();

			String outFile = outDir	+ ext.rootOf(vcf) + "_singleton.counts_GQ_" + gqThreshold
												+ "_sampleRemoved_" + sampsToExclude.size() + ".txt";
			for (String anno : dotAnnos) {
				outFile = ext.addToRoot(outFile, anno);
			}
			StringBuilder buildr = new StringBuilder();
			buildr.append("SAMPLE\tSINGLETONS\tALT_ALLELES\n");
			for (String sample : singletonCounts.keySet()) {
				buildr.append(sample	+ "\t" + singletonCounts.get(sample) + "\t" + altCounts.get(sample)
											+ "\n");
			}
			for (String sample : sampsToExclude) {
				buildr.append(sample + "\tNA\tNA\n");
			}
			Files.write(buildr.toString(), outFile);
		}
	}

	public static void main(String[] args) {
		CLI c = new CLI(Singletons.class);
		ArrayList<String> dotAnnos = new ArrayList<String>();
		dotAnnos.add("esp6500si_all");
		dotAnnos.add("ExAC_ALL");
		dotAnnos.add("charge.n_blacks");
		dotAnnos.add("charge.n_whites");


		HashSet<String> sampsToExclude = new HashSet<String>();
		sampsToExclude.add("PT434_T030");
		sampsToExclude.add("PT431_T029");
		sampsToExclude.add("PT427_T028");


		c.addArgWithDefault("vcf", "vcf to annotate with default methods", "a.vcf");
		c.addArgWithDefault(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, "out/");
		c.addArgWithDefault("gq", "GQ threshold", "20");


		c.parseWithExit(args);
		detectSingletons(c.get("vcf").split(","), c.get(CLI.ARG_OUTDIR), 50, dotAnnos, sampsToExclude);
		detectSingletons(	c.get("vcf").split(","), c.get(CLI.ARG_OUTDIR), 50, dotAnnos,
											new HashSet<String>());

		dotAnnos.add("PopFreqMax");

		detectSingletons(c.get("vcf").split(","), c.get(CLI.ARG_OUTDIR), 50, dotAnnos, sampsToExclude);
		detectSingletons(	c.get("vcf").split(","), c.get(CLI.ARG_OUTDIR), 50, dotAnnos,
											new HashSet<String>());
		dotAnnos.add("snp138");
		detectSingletons(c.get("vcf").split(","), c.get(CLI.ARG_OUTDIR), 50, dotAnnos, sampsToExclude);
		// detectSingletons(c.get("vcf").split(","), c.get("outDir"), 50, dotAnnos, new
		// HashSet<String>());


	}
}
