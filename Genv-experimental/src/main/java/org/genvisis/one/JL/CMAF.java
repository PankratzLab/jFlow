package org.genvisis.one.JL;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.genvisis.common.AlleleFreq;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.seq.analysis.VCFSourceReader;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCOps;

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
		return group + "\t" + freq;

	}

	private static void compute(String vcf, double[] freqs, List<String> groups, int numSamples, boolean useAC,
			String outDir, String requiredAnno, Logger log) {
		VCFFileReader reader = new VCFSourceReader(vcf, false);
		HashMap<String, TrackCmaf> tracker = new HashMap<String, CMAF.TrackCmaf>();
		for (String group : groups) {
			for (int i = 0; i < freqs.length; i++) {
				TrackCmaf current = new TrackCmaf();
				tracker.put(getKey(group, freqs[i]), current);
			}

		}

		int[] numVars = new int[freqs.length];
		int[] numMuts = new int[freqs.length];

		for (VariantContext vc : reader) {

			for (int i = 0; i < freqs.length; i++) {
				int groupI = 0;
				for (String group : groups) {

					if (vc.hasAttribute(group)
							&& (requiredAnno == null || !vc.getAttributeAsString(requiredAnno, ".").equals("."))) {
						double freq = freqs[i];
						int ac = -1;
						try {
							double maf = Double.NaN;
							if (useAC) {
								ac = Array
										.sum(Array.toIntArray(vc.getAttributeAsString("AC", "NaN").replaceAll("\\[", "")
												.replaceAll("\\]", "").replaceAll(" ", "").trim().split(",")));
								maf = AlleleFreq.calcMAF(numSamples - ac, 0, ac);
							} else {
								maf = VCOps.getMAF(vc, null);

							}
							String anno = vc.getAttributeAsString(group, ".");
							if (maf <= freq) {
								if (groupI == 0 && ac > 0) {
									numVars[i]++;
									numMuts[i] += ac;
								}
								if (tracker.get(getKey(group, freqs[i])).containsKey(anno)) {
									double current = tracker.get(getKey(group, freqs[i])).get(anno);
									tracker.get(getKey(group, freqs[i])).put(anno, maf + current);
								} else {
									tracker.get(getKey(group, freqs[i])).put(anno, maf);
								}
							}
						} catch (NumberFormatException nfe) {

						}
					}
					groupI++;
				}
			}
		}
		reader.close();

		ArrayList<String> out = new ArrayList<String>();
		for (String group : tracker.keySet()) {
			if (!group.startsWith("AC")) {
				for (String subGroup : tracker.get(group).keySet()) {
					out.add(group + "\t" + subGroup + "\t" + tracker.get(group).get(subGroup));
				}
			}
		}
		log.reportTimeInfo(Array.toStr(numVars));
		log.reportTimeInfo(Array.toStr(numMuts));

		Files.writeIterable(out, outDir + VCFOps.getAppropriateRoot(vcf, true)
				+ (requiredAnno == null ? "" : requiredAnno) + ".cmaf.txt");

	}

	public static void main(String[] args) {
		String vcf = "/Volumes/Beta/data/Cushings/mito/processDir/polymorphismsMTVariants.pos_1.posAdjust_-1.hg19_multianno.eff.gatk.sed1000g.posAdjust_1.vcf";
		ArrayList<String> groups = new ArrayList<String>();
		groups.add("AC");
		groups.add("Uniprot_name");

		groups.add("OXPHOS_complex");
		groups.add("Associated_disease");
		double[] freqs = new double[] { 1000, .01, .001, .0001 };
		int numSamples = 32059;
		String outDir = "/Volumes/Beta/data/Cushings/mito/processDir/";
		String second = "/Volumes/Beta/data/Cushings/mito/CUSHING_FREQ_V2_Mito/CUSHING_FREQ_V2.maf_0.0.final.CUSHING_FREQ_V2.vcf.gz";

		compute(vcf, freqs, groups, numSamples, true, outDir, null, new Logger());
		compute(second, freqs, groups, numSamples, false, outDir, null, new Logger());

		String required = "MitImpact_id";
		compute(vcf, freqs, groups, numSamples, true, outDir, required, new Logger());
		compute(second, freqs, groups, numSamples, false, outDir, required, new Logger());

	}

}
