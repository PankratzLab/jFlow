package org.genvisis.one.JL;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.genvisis.common.AlleleFreq;
import org.genvisis.common.Logger;
import org.genvisis.seq.analysis.VCFSourceReader;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class CMAF {

	private static class TrackCmaf extends HashMap<String, Double> {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		private TrackCmaf() {

		}

	}

	private static String getKey(String group, double freq) {
		return group + "_" + freq;

	}

	private static void compute(String vcf, double[] freqs, List<String> groups, int numSamples, Logger log) {
		VCFFileReader reader = new VCFSourceReader(vcf, false);
		HashMap<String, TrackCmaf> tracker = new HashMap<String, CMAF.TrackCmaf>();
		for (String group : groups) {
			for (int i = 0; i < freqs.length; i++) {
				TrackCmaf current = new TrackCmaf();
				tracker.put(getKey(group, freqs[i]), current);
			}

		}

		for (VariantContext vc : reader) {
			for (String group : groups) {
				if (vc.hasAttribute(group)) {
					for (int i = 0; i < freqs.length; i++) {
						double freq = freqs[i];
						try {
							int ac = Integer.parseInt(vc.getAttributeAsString("AC", "NaN"));
							double maf = AlleleFreq.calcMAF(numSamples - ac, 0, ac);
							if (maf <= freq) {
								log.reportTimeInfo(maf + "");
								tracker.get("GD");
							}
						} catch (NumberFormatException nfe) {

						}
					}
				}
			}
		}
		reader.close();

	}

	public static void main(String[] args) {
		String vcf = "/Volumes/Beta/data/Cushings/mito/processDir/polymorphismsMTVariants.pos_1.posAdjust_-1.hg19_multianno.eff.gatk.sed1000g.posAdjust_1.vcf";
		ArrayList<String> groups = new ArrayList<String>();
		groups.add("Uniprot_name");
		groups.add("OXPHOS_complex");
		double[] freqs = new double[] { 1, .1, .01 };
		int numSamples = 32059;
		compute(vcf, freqs, groups, numSamples, new Logger());
	}

}
