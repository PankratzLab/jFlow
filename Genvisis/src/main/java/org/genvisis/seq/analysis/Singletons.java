package org.genvisis.seq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
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


	private static void detectSingletons(	String vcf, String outDir, int gqThreshold,
																				List<String> dotAnnos) {
		new File(outDir).mkdirs();
		Logger log = new Logger(outDir + "singleton.log");
		VCFFileReader reader = new VCFFileReader(new File(vcf), true);
		HashMap<String, Integer> singletonCounts = new HashMap<String, Integer>();
		String[] samples = VCFOps.getSamplesInFile(vcf);
		for (int i = 0; i < samples.length; i++) {
			singletonCounts.put(samples[i], 0);
		}
		int num = 0;
		for (VariantContext vc : reader) {
			num++;
			if (num % 10000 == 0) {
				log.reportTimeInfo("Scanned " + num + " variants");
			}
			if (!vc.isIndel() && !vc.isFiltered() && vc.getHetCount() + vc.getHomVarCount() == 1) {// singleton
																																															// in
																																															// VCF
				boolean found = false;
				for (Genotype g : vc.getGenotypes()) {
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
									if (!".".equals(anno)) {
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
		reader.close();

		String outFile = outDir + "Singleton.counts_GQ_" + gqThreshold + ".txt";
		for (String anno : dotAnnos) {
			outFile = ext.addToRoot(outFile, anno);
		}
		StringBuilder buildr = new StringBuilder();
		buildr.append("SAMPLE\tSINGLETONS\n");
		for (String sample : singletonCounts.keySet()) {
			buildr.append(sample + "\t" + singletonCounts.get(sample) + "\n");
		}
		Files.write(buildr.toString(), outFile);
	}

	public static void main(String[] args) {
		CLI c = new CLI(Singletons.class);
		ArrayList<String> dotAnnos = new ArrayList<String>();
		dotAnnos.add("PopFreqMax");
		c.addArgWithDefault("vcf", "vcf to annotate with default methods", "a.vcf");
		c.addArgWithDefault("outDir", "output directory", "out/");
		c.addArgWithDefault("gq", "GQ threshold", "20");


		c.parseWithExit(args);
		for (int i = 20; i < 80; i += 40) {
			detectSingletons(c.get("vcf"), c.get("outDir"), i, dotAnnos);
		}
	}
}
