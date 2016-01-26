package seq.analysis;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Set;

import javax.jms.IllegalStateException;

import common.Files;
import common.HashVec;
import common.Logger;
import seq.manage.BamOps;
import seq.manage.VCOps;
import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;

/**
 * Wrapper and parser/assembler for Somatic Sniper (https://github.com/genome/somatic-sniper)
 *
 */
public class SomaticSniper {

	private static void run(String bams, String vcfPop, Logger log) {
		VcfPopulation vpop = VcfPopulation.load(vcfPop, POPULATION_TYPE.TUMOR_NORMAL, log);
		String[] bamFiles = null;
		if (Files.isDirectory(bams)) {
			bamFiles = Files.list(bams, ".bam", false);
		} else {
			bamFiles = HashVec.loadFileToStringArray(bams, false, new int[] { 0 }, true);
		}
	}

	private static TNSample[] matchSamples(String[] bamFiles, VcfPopulation vpop, Logger log) {
		Hashtable<String, String> all = new Hashtable<String, String>();
		for (int i = 0; i < bamFiles.length; i++) {
			all.put(BamOps.getSampleName(bamFiles[i]), bamFiles[i]);
		}
		ArrayList<TNSample> tnSamples = new ArrayList<TNSample>();
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
			if (!all.containsKey(tumor) || !all.containsKey(normal)) {
				throw new IllegalArgumentException("Could not find bam file for Tumor " + tumor + " or for Normal " + normal);
			}
		}
		return tnSamples.toArray(new TNSample[tnSamples.size()]);
	}

	private static class TNSample {
		private String normalBam;
		private String normalSam;
		private String tumorBam;
		private String tumorSample;

	}

}
